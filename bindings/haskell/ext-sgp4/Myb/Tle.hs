-- SPDX-License-Identifier: BSD-3-Clause
-- | Two-line elements and SGP4 (ext/sgp4, the library myb:tle over myb:sgp4): records read,
-- written and converted to and from OMM, the propagation of the theory, the fit of a set to
-- positions, and the SPK segments of type 10. Functions of numbers alone are pure; the records
-- are those of "Myb.C.Ext.Sgp4".
module Myb.Tle
  ( S.Myb_tle (..)
  , tleParse
  , tleEpoch
  , tleFormat
  , tleFromOmm
  , ommFromTle
    -- * SGP4
  , Sgp4Gravity (..)
  , Sgp4Mode (..)
  , S.Myb_sgp4_constants (..)
  , S.Myb_sgp4 (..)
  , sgp4Constants
  , sgp4Init
  , sgp4InitElements
  , sgp4Propagate
  , sgp4PropagateKm
  , sgp4Gmst
    -- * Fits
  , S.Myb_tle_fit_options (..)
  , tleFit
  , tleFromState
    -- * SPK segments of type 10
  , spk10State
  , spk10Constants
  , spk10Packet
  , spk10Words
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CInt, CSize)
import Foreign.C.String (withCString)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca, allocaBytes)
import Foreign.Marshal.Array (allocaArray, peekArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import Foreign.Storable (peek, poke)
import Foreign.C.String (peekCString)
import qualified Myb.C as C
import qualified Myb.C.Ext.Sgp4 as S
import Myb.Ephemeris (Ephemeris, ephemerisPtr)
import Myb.Error
import Myb.Time (TimeSystem, timeSystemPtr)
import Myb.Types
import System.IO.Unsafe (unsafeDupablePerformIO)

pureCall :: IO a -> a
pureCall = unsafeDupablePerformIO

outcome :: C.Myb_status -> Result ()
outcome code = maybe (Right ()) (\s -> Left (Error s "")) (statusFromCode code)

after :: C.Myb_status -> IO (Result a) -> IO (Result a)
after code next = either (pure . Left) (const next) (outcome code)

-- | The diagnostic of a parser as an error: the status, the line and the message.
diagnosed :: C.Myb_status -> Ptr C.Myb_diag -> IO (Result ())
diagnosed code diag = case statusFromCode code of
  Nothing -> pure (Right ())
  Just status -> do
    d <- peek diag
    let text = map (toEnum . fromIntegral) (takeWhile (/= 0) (C.myb_diag'message d))
    pure (Left (Error status text))

-- | A record from its two lines of 69 characters (no checksum: the checksums are not checked).
tleParse :: String -> String -> Bool -> Result S.Myb_tle
tleParse line1 line2 noChecksum = pureCall $
  withCString line1 $ \p1 -> withCString line2 $ \p2 -> with C.zero $ \tle -> with C.zero $ \diag -> do
    code <- S.myb_tle_parse p1 p2 (if noChecksum then S.MYB_TLE_NO_CHECKSUM else 0) tle diag
    done <- diagnosed code diag
    either (pure . Left) (const (Right <$> peek tle)) done

-- | The epoch of the record as a UTC label.
tleEpoch :: S.Myb_tle -> Calendar
tleEpoch tle = pureCall $ with tle $ \pt -> with C.zero $ \label -> S.myb_tle_epoch pt label >> calendarFromC <$> peek label

-- | The two lines of the record.
tleFormat :: S.Myb_tle -> Result (String, String)
tleFormat tle = pureCall $ with tle $ \pt -> allocaBytes 70 $ \l1 -> allocaBytes 70 $ \l2 -> do
  code <- S.myb_tle_format pt l1 l2
  after code (Right <$> ((,) <$> peekCString l1 <*> peekCString l2))

-- | The record of an OMM with the TLE parameters.
tleFromOmm :: C.Myb_omm -> Result S.Myb_tle
tleFromOmm omm = pureCall $ with omm $ \po -> with C.zero $ \pt -> S.myb_tle_from_omm po pt >>= \code -> after code (Right <$> peek pt)

-- | The OMM of a record.
ommFromTle :: S.Myb_tle -> Result C.Myb_omm
ommFromTle tle = pureCall $ with tle $ \pt -> with C.zero $ \po -> S.myb_omm_from_tle pt po >>= \code -> after code (Right <$> peek po)

data Sgp4Gravity = Wgs72Old | Wgs72 | Wgs84
  deriving (Eq, Show)

data Sgp4Mode = Improved | Afspc
  deriving (Eq, Show)

gravityCode :: Sgp4Gravity -> S.Myb_sgp4_gravity
gravityCode g = case g of
  Wgs72Old -> S.MYB_SGP4_WGS72_OLD
  Wgs72 -> S.MYB_SGP4_WGS72
  Wgs84 -> S.MYB_SGP4_WGS84

modeCode :: Sgp4Mode -> S.Myb_sgp4_mode
modeCode m = case m of
  Improved -> S.MYB_SGP4_IMPROVED
  Afspc -> S.MYB_SGP4_AFSPC

-- | The constants of a model of gravity.
sgp4Constants :: Sgp4Gravity -> S.Myb_sgp4_constants
sgp4Constants g = pureCall $ with C.zero $ \p -> S.myb_sgp4_constants_of (gravityCode g) p >> peek p

-- | The propagator of a record.
sgp4Init :: S.Myb_tle -> Sgp4Gravity -> Sgp4Mode -> Result S.Myb_sgp4
sgp4Init tle g m = pureCall $ with tle $ \pt -> with C.zero $ \ps -> S.myb_sgp4_init pt (gravityCode g) (modeCode m) ps >>= \code -> after code (Right <$> peek ps)

-- | The propagator from elements in the units of the theory: epoch, B*, eccentricity, argument
-- of perigee, inclination, mean anomaly, mean motion (rad/min), node.
sgp4InitElements :: S.Myb_sgp4_constants -> Sgp4Mode -> (Double, Double, Double, Double, Double, Double, Double, Double) -> Result S.Myb_sgp4
sgp4InitElements constants m (epoch, bstar, e, argp, i, ma, n, node) = pureCall $ with constants $ \pc -> with C.zero $ \ps -> do
  code <- S.myb_sgp4_init_elements pc (modeCode m) (coerce epoch) (coerce bstar) (coerce e) (coerce argp) (coerce i) (coerce ma) (coerce n) (coerce node) ps
  after code (Right <$> peek ps)

propagateWith :: (Ptr S.Myb_sgp4 -> CDouble -> Ptr CDouble -> Ptr CDouble -> Ptr S.Myb_sgp4_error -> IO C.Myb_status) -> S.Myb_sgp4 -> Double -> Result (Vec3, Vec3)
propagateWith f sgp4 minutes = pureCall $ with sgp4 $ \ps -> alloca $ \r -> alloca $ \v -> alloca $ \err -> do
  poke err (S.Myb_sgp4_error 0)
  code <- f ps (coerce minutes) (castPtr r) (castPtr v) err
  case statusFromCode code of
    Nothing -> Right <$> ((,) <$> peek r <*> peek v)
    Just status -> do
      S.Myb_sgp4_error e <- peek err
      pure (Left (Error status ("SGP4 error " ++ show e)))

-- | The TEME position (m) and velocity (m/s) at minutes since the epoch; a failure names the
-- error code of the theory.
sgp4Propagate :: S.Myb_sgp4 -> Double -> Result (Vec3, Vec3)
sgp4Propagate = propagateWith S.myb_sgp4_propagate

-- | The same in km and km/s.
sgp4PropagateKm :: S.Myb_sgp4 -> Double -> Result (Vec3, Vec3)
sgp4PropagateKm = propagateWith S.myb_sgp4_propagate_km

-- | Greenwich mean sidereal time of IAU 1982 at the UT1 Julian date jd1 + jd2.
sgp4Gmst :: Double -> Double -> Double
sgp4Gmst jd1 jd2 = coerce (pureCall (S.myb_sgp4_gmst (coerce jd1) (coerce jd2)))

-- | The elements of the record (the first guess) fitted to TEME positions at minutes from its
-- epoch: the record, the RMS of the residuals (m) and the steps.
tleFit :: [Double] -> [Vec3] -> S.Myb_tle_fit_options -> S.Myb_tle -> Result (S.Myb_tle, Double, Int)
tleFit minutes positions options tle
  | length minutes /= length positions = Left (Error Domain "")
  | otherwise = pureCall $ withArrayLen (coerce minutes :: [CDouble]) $ \n pm -> withArrayLen positions $ \_ pp -> with options $ \po -> with tle $ \pt -> alloca $ \rms -> alloca $ \its -> do
      poke rms 0
      poke its (0 :: CInt)
      code <- S.myb_tle_fit pm (castPtr pp) (fromIntegral n :: CSize) po pt rms its
      after code (Right <$> ((,,) <$> peek pt <*> (coerce <$> peek rms) <*> (fromIntegral <$> peek its)))

-- | A first guess for the fit from a TEME state at the epoch of the record (whose epoch, B* and
-- identification are kept).
tleFromState :: Vec3 -> Vec3 -> Sgp4Gravity -> S.Myb_tle -> Result S.Myb_tle
tleFromState r v g tle = pureCall $ with r $ \pr -> with v $ \pv -> with tle $ \pt -> do
  code <- S.myb_tle_from_state (castPtr pr) (castPtr pv) (gravityCode g) pt
  after code (Right <$> peek pt)

-- | The state of segment `index` of an SPK set (type 10) at the TDB instant, in J2000; the
-- time system converts TDB to UTC, the equinox model (an object of "Myb.C") gives TEME.
spk10State :: Ephemeris -> Int -> TimeSystem -> ForeignPtr C.Myb_equinox -> Time -> IO (Result State)
spk10State eph index system equinox t =
  withForeignPtr (ephemerisPtr eph) $ \ps -> withForeignPtr (timeSystemPtr system) $ \pt -> withForeignPtr equinox $ \pe -> alloca $ \y -> do
    code <- S.myb_sgp4_spk10_state ps (fromIntegral index) pt pe (timeToC t) (castPtr y)
    after code (Right <$> peek y)

-- | The eight constants of a segment of type 10.
spk10Constants :: S.Myb_sgp4_constants -> [Double]
spk10Constants model = pureCall $ with model $ \pm -> allocaArray 8 $ \o -> S.myb_sgp4_spk10_constants pm o >> (coerce <$> peekArray 8 o)

-- | The packet (14) of a record, with the nutation from the equinox model ('Nothing': zeros).
spk10Packet :: S.Myb_tle -> TimeSystem -> Maybe (ForeignPtr C.Myb_equinox) -> IO (Result [Double])
spk10Packet tle system equinox =
  with tle $ \pt -> withForeignPtr (timeSystemPtr system) $ \ps -> withEquinox $ \pe -> allocaArray 14 $ \o -> do
    code <- S.myb_sgp4_spk10_packet pt ps pe o
    after code (Right . (coerce :: [CDouble] -> [Double]) <$> peekArray 14 o)
  where
    withEquinox k = maybe (k nullPtr) (`withForeignPtr` k) equinox

-- | The words of a segment of type 10 from the constants (8) and the packets (14 each).
spk10Words :: [Double] -> [Double] -> Result [Double]
spk10Words constants packets
  | length constants /= 8 || length packets `mod` 14 /= 0 = Left (Error Domain "")
  | otherwise = pureCall $ withArrayLen (coerce constants :: [CDouble]) $ \_ pc -> withArrayLen (coerce packets :: [CDouble]) $ \_ pp -> do
      let count = length packets `div` 14
      size <- S.myb_sgp4_spk10_size (fromIntegral count)
      allocaArray (max 1 (fromIntegral size)) $ \o -> do
        code <- S.myb_sgp4_spk10_words pc pp (fromIntegral count) o size
        after code (Right . (coerce :: [CDouble] -> [Double]) <$> peekArray (fromIntegral size) o)
