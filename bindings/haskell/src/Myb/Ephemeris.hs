-- SPDX-License-Identifier: BSD-3-Clause
-- | A set of SPK files (myb/kin.h): the states of their bodies at instants of the files' time
-- scale.
module Myb.Ephemeris
  ( Ephemeris
  , newEphemeris
  , loadEphemerisFile
  , ephemerisState
  , Segment (..)
  , segments
  , ephemerisPtr
  ) where

import Data.Int (Int32)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, nullPtr)
import Foreign.Storable (peek)
import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import qualified GHC.Foreign as Foreign
import GHC.IO.Encoding (char8)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

newtype Ephemeris = Ephemeris (ForeignPtr C.Myb_spk)

newEphemeris :: IO (Result Ephemeris)
newEphemeris = create (\out -> check <$> C.myb_spk_create nullPtr out) (fmap Ephemeris . ownedLeaf C.p'myb_spk_destroy)

loadEphemerisFile :: Ephemeris -> FilePath -> IO (Result ())
loadEphemerisFile (Ephemeris spk) path =
  withText path $ \cpath -> withForeignPtr spk $ \s -> withDiag (C.myb_spk_load_file s cpath)

-- | The state of the target relative to the center at an instant of the files' time scale.
ephemerisState :: Ephemeris -> Int32 -> Int32 -> Time -> IO (Result State)
ephemerisState (Ephemeris spk) target center t =
  withForeignPtr spk $ \s -> alloca $ \out -> do
    code <- C.myb_spk_state s target center (timeToC t) (castPtr out) nullPtr
    either (pure . Left) (const (Right <$> peek out)) (check code)

-- | A segment of the set, as the file describes it.
data Segment = Segment
  { segmentTarget :: !Int32
  , segmentCenter :: !Int32
  , segmentFrame :: !Int32
  , segmentType :: !Int32
  , -- | The interval of the segment, in seconds past J2000 of the files' time scale.
    segmentStart :: !Double
  , segmentStop :: !Double
  , -- | The index of the file in the order of loading.
    segmentFile :: !Int
  , segmentName :: !String
  }
  deriving (Eq, Show)

-- | The segments of the set, in the order of loading.
segments :: Ephemeris -> IO [Segment]
segments (Ephemeris spk) = withForeignPtr spk $ \s -> do
  count <- C.myb_spk_segment_count s
  if count == 0 then pure [] else mapM (segment s) [0 .. count - 1]
  where
    segment s index = with (C.zero :: C.Myb_spk_segment_info) $ \out -> do
      C.myb_spk_segment s index out
      info <- peek out
      name <- withArrayLen (takeWhile (/= 0) (C.myb_spk_segment_info'name info)) $ \n p -> Foreign.peekCStringLen char8 (p, n)
      pure
        Segment
          { segmentTarget = C.myb_spk_segment_info'target info
          , segmentCenter = C.myb_spk_segment_info'center info
          , segmentFrame = C.myb_spk_segment_info'frame info
          , segmentType = C.myb_spk_segment_info'type info
          , segmentStart = coerce (C.myb_spk_segment_info'start info)
          , segmentStop = coerce (C.myb_spk_segment_info'stop info)
          , segmentFile = fromSize (C.myb_spk_segment_info'file info)
          , segmentName = name
          }

ephemerisPtr :: Ephemeris -> ForeignPtr C.Myb_spk
ephemerisPtr (Ephemeris spk) = spk
