-- SPDX-License-Identifier: BSD-3-Clause
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
  ) where

import Control.Monad (join)
import Data.IORef (IORef, modifyIORef', newIORef, readIORef)
import Foreign.C.String (CString)
import Foreign.C.Types (CSize)
import qualified Foreign.Concurrent as Concurrent
import Foreign.ForeignPtr (FinalizerPtr, ForeignPtr, mallocForeignPtr, newForeignPtr, touchForeignPtr, withForeignPtr)
import Foreign.ForeignPtr.Unsafe (unsafeForeignPtrToPtr)
import Foreign.Marshal.Alloc (alloca, allocaBytes)
import Foreign.Marshal.Array (withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (FunPtr, Ptr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (Storable, peek, poke)
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

-- | A text of the library's writers, which tell the length when the buffer is too small.
writeText :: (CString -> CSize -> Ptr CSize -> IO C.Myb_status) -> IO (Result String)
writeText write = alloca $ \lengthPtr -> do
  poke lengthPtr 0
  code <- write nullPtr 0 lengthPtr
  if code /= C.MYB_OK && code /= C.MYB_E_OVERFLOW
    then pure (check code >> failure Domain)
    else do
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
