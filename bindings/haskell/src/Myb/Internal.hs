-- SPDX-License-Identifier: BSD-3-Clause
{-# LANGUAGE RankNTypes #-}

-- | What the modules of the upper layer share: outcomes from statuses, strings for C, the
-- ownership of objects of the library. Not exported from the package.
module Myb.Internal
  ( -- * Outcomes
    check
  , failure
  , withDiag
  , create
    -- * Strings
  , withText
  , withTextLen
  , peekText
  , writeText
    -- * Ownership
  , Kept
  , newKept
  , keepAlive
  , keepRelease
  , keepFunPtr
  , keepValue
  , owned
  , ownedLeaf
  , rawPtr
    -- * Numbers
  , toSize
  , fromSize
    -- * Calls
  , andThen
  , fill
  , fillZero
  , fill2
  , fillArray
  , fillDoubles
  , fillVec
  , withDoubles
  , withDoublesLen
  , withVec
  , withVecs
  , peekDoubles
  , pureCall
  , arrayForm
  , requireLength
  , withMany
  , newOwned
  , withOwned
  , dependent
  , copyBytes'
    -- * Providers
  , StateSource (..)
  , withStateSource
  , keepStateSource
  , OrientationSource (..)
  , withOrientationSource
  , keepOrientationSource
  ) where

import Control.Monad (join)
import Data.Coerce (coerce)
import Data.IORef (IORef, modifyIORef', newIORef, readIORef)
import Data.Word (Word8)
import Foreign.C.String (CString)
import Foreign.C.Types (CDouble (..), CSize)
import qualified Foreign.Concurrent as Concurrent
import Foreign.ForeignPtr (FinalizerPtr, ForeignPtr, mallocForeignPtr, newForeignPtr, touchForeignPtr, withForeignPtr)
import Foreign.ForeignPtr.Unsafe (unsafeForeignPtrToPtr)
import Foreign.Marshal.Alloc (alloca, allocaBytes)
import Foreign.Marshal.Array (allocaArray, peekArray, withArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (FunPtr, Ptr, castPtr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (Storable, peek, poke)
import System.IO.Unsafe (unsafeDupablePerformIO)
import qualified GHC.Foreign as Foreign
import GHC.IO.Encoding.Failure (CodingFailureMode (RoundtripFailure))
import GHC.IO.Encoding.Types (TextEncoding)
import GHC.IO.Encoding.UTF8 (mkUTF8)
import qualified Myb.C as C
import Myb.Error

-- ---- Outcomes ------------------------------------------------------------------------------

-- | The outcome of a status of the C interface.
check :: C.Myb_status -> Result ()
check code = maybe (Right ()) (\status -> Left (Error status "")) (statusFromCode code)

failure :: Status -> Result a
failure status = Left (Error status "")

-- | A call that fills a diagnostic of the library: its outcome, with the message (and the
-- line) the library wrote.
withDiag :: (Ptr C.Myb_diag -> IO C.Myb_status) -> IO (Result ())
withDiag call = with (C.zero :: C.Myb_diag) $ \diag -> do
  code <- call diag
  case statusFromCode code of
    Nothing -> pure (Right ())
    Just status -> do
      written <- peek diag
      let bytes = takeWhile (/= 0) (C.myb_diag'message written)
          line = C.myb_diag'line written
      text <- withArrayLen bytes $ \n p -> Foreign.peekCStringLen encoding (p, n)
      pure (Left (Error status (if line > 0 then "line " ++ show line ++ ": " ++ text else text)))

-- | An object the library creates into a pointer, wrapped once it is there.
create :: (Ptr (Ptr a) -> IO (Result ())) -> (Ptr a -> IO b) -> IO (Result b)
create call wrap = alloca $ \out -> do
  poke out nullPtr
  made <- call out
  case made of
    Left e -> pure (Left e)
    Right () -> do
      p <- peek out
      if p == nullPtr then pure (failure NoMemory) else Right <$> wrap p

-- ---- Strings -------------------------------------------------------------------------------

-- | UTF-8, as the library reads paths and names. Bytes that are not UTF-8 come back as the
-- escapes of GHC and go out again as the same bytes: decoding never fails.
encoding :: TextEncoding
encoding = mkUTF8 RoundtripFailure

-- | Whether C can take the character: not NUL, not a surrogate other than the escape of a byte.
acceptable :: Char -> Bool
acceptable c = c /= '\0' && (c < '\xD800' || c > '\xDFFF' || (c >= '\xDC80' && c <= '\xDCFF'))

-- | A string as C takes it, NUL-terminated; one it cannot take is a failure of the domain.
withText :: String -> (CString -> IO (Result a)) -> IO (Result a)
withText s k
  | all acceptable s = Foreign.withCString encoding s k
  | otherwise = pure (failure Domain)

-- | A text with its length in bytes.
withTextLen :: String -> (CString -> CSize -> IO (Result a)) -> IO (Result a)
withTextLen s k
  | all (\c -> c == '\0' || acceptable c) s = Foreign.withCStringLen encoding s (\(p, n) -> k p (toSize n))
  | otherwise = pure (failure Domain)

peekText :: CString -> IO String
peekText = Foreign.peekCString encoding

-- | A text of the library's writers, which tell the length when the buffer is too small: a
-- buffer of 256 bytes first (some writers take no NULL buffer), then one of the length told.
writeText :: (CString -> CSize -> Ptr CSize -> IO C.Myb_status) -> IO (Result String)
writeText write = alloca $ \lengthPtr -> do
  poke lengthPtr 0
  first <- allocaBytes 256 $ \buffer -> do
    code <- write buffer 256 lengthPtr
    if code == C.MYB_OK
      then Right . Right <$> (peek lengthPtr >>= \n -> Foreign.peekCStringLen encoding (buffer, fromSize n))
      else pure (if code == C.MYB_E_OVERFLOW then Left () else Right (check code >> failure Domain))
  case first of
    Right done -> pure done
    Left () -> do
      size <- peek lengthPtr
      allocaBytes (fromSize size + 1) $ \buffer -> do
        written <- check <$> write buffer (size + 1) lengthPtr
        case written of
          Left e -> pure (Left e)
          Right () -> do
            n <- peek lengthPtr
            Right <$> Foreign.peekCStringLen encoding (buffer, fromSize n)

-- ---- Ownership -----------------------------------------------------------------------------

-- | What an object of the library points to and must outlive (kept alive until the object is
-- destroyed), and what is released after it (Haskell functions C held, work spaces).
data Kept = Kept (IORef [IO ()]) (IORef [IO ()])

newKept :: IO Kept
newKept = Kept <$> newIORef [] <*> newIORef []

-- | The object keeps another alive: its finalizer touches it after destroying the object, so
-- that the other is not finalized first.
keepAlive :: Kept -> ForeignPtr a -> IO ()
keepAlive (Kept alive _) fp = modifyIORef' alive (touchForeignPtr fp :)

-- | Something to do once the object is destroyed.
keepRelease :: Kept -> IO () -> IO ()
keepRelease (Kept _ after) action = modifyIORef' after (action :)

-- | A Haskell function C holds: freed once the object is destroyed.
keepFunPtr :: Kept -> FunPtr a -> IO ()
keepFunPtr kept f = keepRelease kept (freeHaskellFunPtr f)

-- | A value in memory the object points to; the pointer is valid as long as the object.
keepValue :: Storable a => Kept -> a -> IO (Ptr a)
keepValue kept value = do
  fp <- mallocForeignPtr
  withForeignPtr fp (`poke` value)
  keepAlive kept fp
  pure (unsafeForeignPtrToPtr fp)

-- | An object of the library that keeps others: once unreachable it is destroyed, then what
-- it held is released and what it kept may go.
owned :: Ptr a -> (Ptr a -> IO ()) -> Kept -> IO (ForeignPtr a)
owned p destroy (Kept alive after) = Concurrent.newForeignPtr p $ do
  destroy p
  join (sequence_ <$> readIORef after)
  join (sequence_ <$> readIORef alive)

-- | An object of the library that keeps nothing: its C function destroys it.
ownedLeaf :: FinalizerPtr a -> Ptr a -> IO (ForeignPtr a)
ownedLeaf destroy p = newForeignPtr destroy p

-- | The pointer of an object, to be stored in a structure of C: valid only while the object
-- is kept alive (keepAlive).
rawPtr :: ForeignPtr a -> Ptr a
rawPtr = unsafeForeignPtrToPtr

-- ---- Numbers -------------------------------------------------------------------------------

toSize :: Int -> CSize
toSize = fromIntegral

fromSize :: CSize -> Int
fromSize = fromIntegral

-- ---- Calls ---------------------------------------------------------------------------------

-- | The second action when the status is MYB_OK; the failure otherwise.
andThen :: C.Myb_status -> IO (Result a) -> IO (Result a)
andThen code next = either (pure . Left) (const next) (check code)

-- | A call that writes one value through a pointer.
fill :: Storable a => (Ptr a -> IO C.Myb_status) -> IO (Result a)
fill f = alloca $ \p -> f p >>= \code -> andThen code (Right <$> peek p)

-- | The same with the place zeroed first: for a structure the library fills in part.
fillZero :: (Storable a, C.Zero a) => (Ptr a -> IO C.Myb_status) -> IO (Result a)
fillZero f = with C.zero $ \p -> f p >>= \code -> andThen code (Right <$> peek p)

-- | A call that writes two values.
fill2 :: (Storable a, Storable b) => (Ptr a -> Ptr b -> IO C.Myb_status) -> IO (Result (a, b))
fill2 f = alloca $ \p -> alloca $ \q -> f p q >>= \code -> andThen code (Right <$> ((,) <$> peek p <*> peek q))

-- | A call that writes n values in a row.
fillArray :: Storable a => Int -> (Ptr a -> IO C.Myb_status) -> IO (Result [a])
fillArray n f = allocaArray (max 1 n) $ \p -> f p >>= \code -> andThen code (Right <$> peekArray n p)

-- | A call that writes n doubles.
fillDoubles :: Int -> (Ptr CDouble -> IO C.Myb_status) -> IO (Result [Double])
fillDoubles n f = fmap (coerce :: [CDouble] -> [Double]) <$> fillArray n f

-- | A call that writes a value of doubles (a Vec3, a State, ...) through a pointer to doubles.
fillVec :: Storable v => (Ptr CDouble -> IO C.Myb_status) -> IO (Result v)
fillVec f = fill (f . castPtr)

withDoubles :: [Double] -> (Ptr CDouble -> IO a) -> IO a
withDoubles xs = withArray (coerce xs :: [CDouble])

withDoublesLen :: [Double] -> (Int -> Ptr CDouble -> IO a) -> IO a
withDoublesLen xs = withArrayLen (coerce xs :: [CDouble])

-- | A value of doubles (Vec3, State) as the pointer to its doubles.
withVec :: Storable v => v -> (Ptr CDouble -> IO a) -> IO a
withVec v k = with v (k . castPtr)

-- | Values of doubles one after another (positions of bodies: 3 doubles each).
withVecs :: Storable v => [v] -> (Ptr CDouble -> IO a) -> IO a
withVecs vs k = withArray vs (k . castPtr)

peekDoubles :: Int -> Ptr CDouble -> IO [Double]
peekDoubles n p = (coerce :: [CDouble] -> [Double]) <$> peekArray n p

-- | A function of numbers alone, pure as it is in C.
pureCall :: IO a -> a
pureCall = unsafeDupablePerformIO

-- | An array form: on a failure the error names the element that failed ("element 3").
arrayForm :: (Ptr CSize -> IO C.Myb_status) -> IO (Result ())
arrayForm f = alloca $ \failed -> do
  poke failed 0
  code <- f failed
  case statusFromCode code of
    Nothing -> pure (Right ())
    Just status -> do
      k <- peek failed
      pure (Left (Error status ("element " ++ show k)))

-- | Goes on only when the list has the length the call needs; a failure of the domain otherwise.
requireLength :: Int -> [a] -> IO (Result b) -> IO (Result b)
requireLength n xs k = if length xs == n then k else pure (failure Domain)

-- | Continuations over a list: each item marshalled for C by the same function.
withMany :: (a -> (b -> IO r) -> IO r) -> [a] -> ([b] -> IO r) -> IO r
withMany _ [] k = k []
withMany f (x : rest) k = f x $ \y -> withMany f rest (k . (y :))

-- | An object of the library with what it keeps: created into a pointer, destroyed by its
-- function once unreachable.
newOwned :: (Ptr (Ptr a) -> IO (Result ())) -> (Ptr a -> IO ()) -> IO (Result (ForeignPtr a, Kept))
newOwned call destroy = do
  kept <- newKept
  create call (\p -> (\fp -> (fp, kept)) <$> owned p destroy kept)

withOwned :: ForeignPtr a -> (Ptr a -> IO b) -> IO b
withOwned = withForeignPtr

-- | A pointer that another object owns and keeps (the point masses of a system): no
-- finalizer of its own, and the owner lives as long as it does.
dependent :: ForeignPtr owner -> Ptr a -> IO (ForeignPtr a)
dependent owner p = Concurrent.newForeignPtr p (touchForeignPtr owner)

-- | Bytes from a buffer of C.
copyBytes' :: Ptr () -> Int -> IO [Word8]
copyBytes' p n = peekArray n (castPtr p)

-- ---- Providers -----------------------------------------------------------------------------

-- | Where the state of a point comes from (a myb_state_function and its context): for the
-- length of one call, or for the life of an object that keeps what it needs.
data StateSource
  = StateSource
      (forall r. (C.Myb_state_function -> Ptr () -> IO r) -> IO r)
      (Kept -> IO (C.Myb_state_function, Ptr ()))

withStateSource :: StateSource -> (C.Myb_state_function -> Ptr () -> IO r) -> IO r
withStateSource (StateSource w _) = w

keepStateSource :: Kept -> StateSource -> IO (C.Myb_state_function, Ptr ())
keepStateSource kept (StateSource _ k) = k kept

-- | Where the orientation of a frame comes from (a myb_orientation_function and its context).
data OrientationSource
  = OrientationSource
      (forall r. (C.Myb_orientation_function -> Ptr () -> IO r) -> IO r)
      (Kept -> IO (C.Myb_orientation_function, Ptr ()))

withOrientationSource :: OrientationSource -> (C.Myb_orientation_function -> Ptr () -> IO r) -> IO r
withOrientationSource (OrientationSource w _) = w

keepOrientationSource :: Kept -> OrientationSource -> IO (C.Myb_orientation_function, Ptr ())
keepOrientationSource kept (OrientationSource _ k) = k kept
