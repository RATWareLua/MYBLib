-- SPDX-License-Identifier: BSD-3-Clause
-- | Kernels and time scales (myb/io.h, myb/time.h). Scales are named: "TT", "TAI", "TDB", ...;
-- UTC is a calendar label.
module Myb.Time
  ( -- * A pool of text kernels
    KernelPool
  , newKernelPool
  , loadKernelFile
  , loadKernelText
  , KernelType (..)
  , describeKernel
  , kernelNumbers
  , kernelString
  , kernelPoolPtr
    -- * Time scales
  , TimeSystem
  , newTimeSystem
  , convertTime
  , fromUtc
  , toUtc
  , timeSystemPtr
    -- * Instants
  , Time (..)
  , Calendar (..)
  , timeFromSeconds
  , seconds
    -- * Array forms
  , timeFromSecondsMany
  , secondsMany
  , calendarToTimeMany
  , timeToCalendarMany
  , convertTimeMany
  , fromUtcMany
  , toUtcMany
  ) where

import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, nullPtr)
import Foreign.Storable (peek)
import qualified GHC.Foreign as Foreign
import GHC.IO.Encoding (char8)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | A pool of text kernels.
newtype KernelPool = KernelPool (ForeignPtr C.Myb_kpool)

newKernelPool :: IO (Result KernelPool)
newKernelPool =
  create (\out -> check <$> C.myb_kpool_create nullPtr out) (fmap KernelPool . ownedLeaf C.p'myb_kpool_destroy)

loadKernelFile :: KernelPool -> FilePath -> IO (Result ())
loadKernelFile (KernelPool pool) path =
  withText path $ \cpath -> withForeignPtr pool $ \p -> withDiag (C.myb_kpool_load_file p cpath)

-- | The object of the library, to mix with "Myb.C".
kernelPoolPtr :: KernelPool -> ForeignPtr C.Myb_kpool
kernelPoolPtr (KernelPool pool) = pool

-- | The time scales of a pool. It keeps nothing of the pool.
newtype TimeSystem = TimeSystem (ForeignPtr C.Myb_timesys)

newTimeSystem :: KernelPool -> IO (Result TimeSystem)
newTimeSystem (KernelPool pool) =
  withForeignPtr pool $ \p ->
    create
      (\out -> withDiag (C.myb_timesys_create p nullPtr out))
      (fmap TimeSystem . ownedLeaf C.p'myb_timesys_destroy)

scale :: TimeSystem -> String -> IO (Result C.Myb_scale)
scale (TimeSystem system) name =
  withText name $ \cname -> withForeignPtr system $ \s -> alloca $ \out -> do
    code <- C.myb_timesys_scale s cname out
    case check code of
      Left e -> pure (Left e)
      Right () -> Right <$> peek out

-- | The instant on another scale: convertTime system t "TT" "TAI".
convertTime :: TimeSystem -> Time -> String -> String -> IO (Result Time)
convertTime ts@(TimeSystem system) t from to = do
  scales <- (,) <$> scale ts from <*> scale ts to
  case scales of
    (Right a, Right b) -> withForeignPtr system $ \s -> alloca $ \out -> do
      code <- C.myb_timesys_convert s (timeToC t) a b out
      either (pure . Left) (const (Right . timeFromC <$> peek out)) (check code)
    (Left e, _) -> pure (Left e)
    (_, Left e) -> pure (Left e)

-- | The instant of a UTC label on a scale.
fromUtc :: TimeSystem -> Calendar -> String -> IO (Result Time)
fromUtc ts@(TimeSystem system) label to = do
  target <- scale ts to
  case target of
    Left e -> pure (Left e)
    Right b -> withForeignPtr system $ \s -> with (calendarToC label) $ \l -> alloca $ \out -> do
      code <- C.myb_timesys_from_utc s l b out
      either (pure . Left) (const (Right . timeFromC <$> peek out)) (check code)

-- | The UTC label of an instant of a scale.
toUtc :: TimeSystem -> Time -> String -> IO (Result Calendar)
toUtc ts@(TimeSystem system) t from = do
  source <- scale ts from
  case source of
    Left e -> pure (Left e)
    Right a -> withForeignPtr system $ \s -> alloca $ \out -> do
      code <- C.myb_timesys_to_utc s (timeToC t) a out
      either (pure . Left) (const (Right . calendarFromC <$> peek out)) (check code)

timeSystemPtr :: TimeSystem -> ForeignPtr C.Myb_timesys
timeSystemPtr (TimeSystem system) = system

-- ---- The pool's variables (myb/io.h) ---------------------------------------------------------

-- | Loads kernel text into the pool (atomic: on an error the pool is unchanged).
loadKernelText :: KernelPool -> String -> IO (Result ())
loadKernelText (KernelPool pool) text =
  withTextLen text $ \ctext n -> withForeignPtr pool $ \p -> withDiag (C.myb_kpool_load_text p ctext n)

-- | What a variable holds.
data KernelType = KernelNumbers | KernelStrings
  deriving (Eq, Show)

-- | The type and the number of values of a variable.
describeKernel :: KernelPool -> String -> IO (Result (KernelType, Int))
describeKernel (KernelPool pool) name =
  withText name $ \cname -> withForeignPtr pool $ \p -> alloca $ \kind -> alloca $ \count -> do
    code <- C.myb_kpool_describe p cname kind count
    andThen code $ do
      k <- peek kind
      n <- peek count
      pure (Right (if k == C.MYB_KVAR_STRINGS then KernelStrings else KernelNumbers, fromSize n))

-- | The numeric values of a variable.
kernelNumbers :: KernelPool -> String -> IO (Result [Double])
kernelNumbers (KernelPool pool) name =
  withText name $ \cname -> withForeignPtr pool $ \p -> alloca $ \values -> alloca $ \count -> do
    code <- C.myb_kpool_numbers p cname values count
    andThen code $ do
      n <- fromSize <$> peek count
      v <- peek values
      Right <$> peekDoubles n v

-- | One string value of a variable.
kernelString :: KernelPool -> String -> Int -> IO (Result String)
kernelString (KernelPool pool) name index =
  withText name $ \cname -> withForeignPtr pool $ \p -> alloca $ \value -> alloca $ \len -> do
    code <- C.myb_kpool_string p cname (toSize index) value len
    andThen code $ do
      v <- peek value
      n <- fromSize <$> peek len
      Right <$> Foreign.peekCStringLen char8 (v, n)

-- ---- Array forms ---------------------------------------------------------------------------

withTimes :: [Time] -> (Int -> Ptr C.Myb_time -> IO a) -> IO a
withTimes ts = withArrayLen (map timeToC ts)

-- | myb_time_from_seconds of every value.
timeFromSecondsMany :: [Double] -> [Time]
timeFromSecondsMany xs = pureCall $ withDoublesLen xs $ \n px -> allocaArray (max 1 n) $ \o -> do
  C.myb_time_from_seconds_many px (toSize n) o
  map timeFromC <$> peekArray n o

-- | myb_time_seconds of every instant.
secondsMany :: [Time] -> [Double]
secondsMany ts = pureCall $ withTimes ts $ \n pt -> allocaArray (max 1 n) $ \o -> do
  C.myb_time_seconds_many pt (toSize n) o
  peekDoubles n o

-- | The instants of calendar labels on their own scale (myb_calendar_to_time).
calendarToTimeMany :: [Calendar] -> Result [Time]
calendarToTimeMany labels = pureCall $ withArrayLen (map calendarToC labels) $ \n pc -> allocaArray (max 1 n) $ \o -> do
  done <- arrayForm (C.myb_calendar_to_time_many pc (toSize n) o)
  either (pure . Left) (const (Right . map timeFromC <$> peekArray n o)) done

-- | The calendar labels of instants (myb_time_to_calendar).
timeToCalendarMany :: [Time] -> [Calendar]
timeToCalendarMany ts = pureCall $ withTimes ts $ \n pt -> allocaArray (max 1 n) $ \o -> do
  C.myb_time_to_calendar_many pt (toSize n) o
  map calendarFromC <$> peekArray n o

-- | 'convertTime' of every instant.
convertTimeMany :: TimeSystem -> [Time] -> String -> String -> IO (Result [Time])
convertTimeMany ts@(TimeSystem system) instants from to = do
  scales <- (,) <$> scale ts from <*> scale ts to
  case scales of
    (Right a, Right b) -> withForeignPtr system $ \s -> withTimes instants $ \n pt -> allocaArray (max 1 n) $ \o -> do
      done <- arrayForm (C.myb_timesys_convert_many s pt (toSize n) a b o)
      either (pure . Left) (const (Right . map timeFromC <$> peekArray n o)) done
    (Left e, _) -> pure (Left e)
    (_, Left e) -> pure (Left e)

-- | 'fromUtc' of every label.
fromUtcMany :: TimeSystem -> [Calendar] -> String -> IO (Result [Time])
fromUtcMany ts@(TimeSystem system) labels to = do
  target <- scale ts to
  case target of
    Left e -> pure (Left e)
    Right b -> withForeignPtr system $ \s -> withArrayLen (map calendarToC labels) $ \n pc -> allocaArray (max 1 n) $ \o -> do
      done <- arrayForm (C.myb_timesys_from_utc_many s pc (toSize n) b o)
      either (pure . Left) (const (Right . map timeFromC <$> peekArray n o)) done

-- | 'toUtc' of every instant.
toUtcMany :: TimeSystem -> [Time] -> String -> IO (Result [Calendar])
toUtcMany ts@(TimeSystem system) instants from = do
  source <- scale ts from
  case source of
    Left e -> pure (Left e)
    Right a -> withForeignPtr system $ \s -> withTimes instants $ \n pt -> allocaArray (max 1 n) $ \o -> do
      done <- arrayForm (C.myb_timesys_to_utc_many s pt (toSize n) a o)
      either (pure . Left) (const (Right . map calendarFromC <$> peekArray n o)) done
