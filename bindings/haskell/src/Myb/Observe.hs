-- SPDX-License-Identifier: BSD-3-Clause
-- | Observables (myb/observe.h): light time with the relativistic delay, round trips, the
-- Shapiro delay, aberration, light deflection, directions, and the dispersive delay in a plasma.
module Myb.Observe
  ( LightModel (..)
  , LightDirection (..)
  , C.Myb_light_path (..)
  , lightTime
  , roundTrip
  , shapiroDelay
  , aberration
  , lightDeflection
  , directionRadec
  , plasmaGroupDelay
  , C.Myb_plasma_profile (..)
  , plasmaColumn
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (peek)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | The speed of light (m/s), the PPN gamma, and the bodies that delay light (GM, where they
-- are; none: the Newtonian light time).
data LightModel = LightModel {modelSpeed :: Double, modelGamma :: Double, modelBodies :: [(Double, StateSource)]}

data LightDirection
  = -- | the observer receives at t what the target emitted earlier
    Reception
  | -- | the observer emits at t, the target receives later
    Transmission
  deriving (Eq, Show)

withModel :: LightModel -> (Ptr C.Myb_light_model -> IO a) -> IO a
withModel model k =
  withMany (\(gm, source) go -> withStateSource source (\f c -> go (C.zero {C.myb_light_body'gm = coerce gm, C.myb_light_body'state = f, C.myb_light_body'context = c})))
    (modelBodies model) $ \bodies -> withArrayLen bodies $ \n pb ->
    with
      C.zero
        { C.myb_light_model'c = coerce (modelSpeed model)
        , C.myb_light_model'gamma = coerce (modelGamma model)
        , C.myb_light_model'bodies = pb
        , C.myb_light_model'count = toSize n
        }
      k

-- | The one-way light time between the observer (its state at t) and the target.
lightTime :: LightModel -> StateSource -> Time -> State -> LightDirection -> IO (Result C.Myb_light_path)
lightTime model target t observer direction =
  withModel model $ \pm -> withStateSource target $ \f c -> withVec observer $ \po ->
    fillZero (C.myb_light_time pm f c (timeToC t) po side)
  where
    side = case direction of
      Reception -> C.MYB_LIGHT_RECEPTION
      Transmission -> C.MYB_LIGHT_TRANSMISSION

-- | A round trip received at t by a station: the total time, the instants of the bounce and of
-- the transmission.
roundTrip :: LightModel -> StateSource -> StateSource -> Time -> IO (Result (Double, Time, Time))
roundTrip model station target t =
  withModel model $ \pm -> withStateSource station $ \fs cs -> withStateSource target $ \ft ct -> alloca $ \total -> alloca $ \bounce -> alloca $ \transmit -> do
    code <- C.myb_round_trip pm fs cs ft ct (timeToC t) total bounce transmit
    andThen code (Right <$> ((,,) <$> (coerce <$> peek total) <*> (timeFromC <$> peek bounce) <*> (timeFromC <$> peek transmit)))

-- | The relativistic delay (s) of a ray from x1 to x2 by a body of the GM at xb.
shapiroDelay :: Double -> Double -> Double -> Vec3 -> Vec3 -> Vec3 -> Double
shapiroDelay gm c gamma x1 x2 xb = coerce $ pureCall $ withVec x1 $ \p1 -> withVec x2 $ \p2 -> withVec xb $ \pb ->
  C.myb_shapiro_delay (coerce gm) (coerce c) (coerce gamma) p1 p2 pb

-- | The direction p as seen by an observer moving with the velocity v.
aberration :: Vec3 -> Vec3 -> Double -> Vec3
aberration p v c = pureCall $ withVec p $ \pp -> withVec v $ \pv -> alloca $ \o -> C.myb_aberration pp pv (coerce c) (castPtr o) >> peek o

-- | The direction p deflected by a body (first order in gm, PPN): q body -> source, e body ->
-- observer, em their distance, dlim the limiter near the limb.
lightDeflection :: Double -> Double -> Double -> Vec3 -> Vec3 -> Vec3 -> Double -> Double -> Vec3
lightDeflection gm c gamma p q e em dlim = pureCall $ withVec p $ \pp -> withVec q $ \pq -> withVec e $ \pe -> alloca $ \o ->
  C.myb_light_deflection (coerce gm) (coerce c) (coerce gamma) pp pq pe (coerce em) (coerce dlim) (castPtr o) >> peek o

-- | Right ascension in [0, 2 pi) and declination of a direction.
directionRadec :: Vec3 -> (Double, Double)
directionRadec p = pureCall $ withVec p $ \pp -> alloca $ \ra -> alloca $ \dec -> do
  C.myb_direction_radec pp ra dec
  (\a b -> (coerce a, coerce b)) <$> peek ra <*> peek dec

-- | K N / f^2 (m) for the electron column (m^-2), the frequency (Hz) and the constants e, m_e,
-- eps0.
plasmaGroupDelay :: Double -> Double -> Double -> Double -> Double -> Double
plasmaGroupDelay column frequency e me eps0 = coerce (pureCall (C.myb_plasma_group_delay (coerce column) (coerce frequency) (coerce e) (coerce me) (coerce eps0)))

-- | The electron column (m^-2) of a profile along the segment from p1 to p2.
plasmaColumn :: C.Myb_plasma_profile -> Vec3 -> Vec3 -> Result Double
plasmaColumn profile p1 p2 = pureCall $ with profile $ \pp -> withVec p1 $ \a -> withVec p2 $ \b -> fmap coerce <$> fill (C.myb_plasma_column pp a b)
