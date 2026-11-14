-- SPDX-License-Identifier: BSD-3-Clause
-- | Kernels and time scales (myb/io.h, myb/time.h). Scales are named: "TT", "TAI", "TDB", ...;
-- UTC is a calendar label.
module Myb.Time
  ( -- * A pool of text kernels
    KernelPool
  , newKernelPool
  , loadKernelFile
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
  ) where

import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (nullPtr)
import Foreign.Storable (peek)
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
