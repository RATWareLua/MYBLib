-- SPDX-License-Identifier: BSD-3-Clause
-- | Thrust (myb/thrust.h): arcs of an engine as a term of a force model, the Gauss equations in
-- equinoctial elements, Edelbaum's transfer, and the Q-law as a direction, a term and a flight
-- with the mass in the state. The descriptions are the records of "Myb.C".
module Myb.Thrust
  ( -- * Arcs of thrust
    C.Myb_thrust_arc (..)
  , C.Myb_thrust (..)
  , thrustInit
  , thrustMass
  , thrustMassAt
  , thrustAcceleration
  , thrustRtn
  , addThrust
    -- * Equinoctial elements and Edelbaum
  , C.Myb_equinoctial (..)
  , equinoctialRates
  , edelbaum
  , edelbaumState
    -- * The Q-law
  , C.Myb_qlaw (..)
  , C.Myb_qlaw_term (..)
  , C.Myb_qlaw_arrival (..)
  , qlawDirection
  , qlawArrived
  , qlawAcceleration
  , addQlaw
  , qlawFlightRhs
  , qlawFly
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (mallocForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullFunPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Force (ForceModel, addRaw, forceModelKeep, forceModelPtr)
import Myb.Internal
import Myb.Types

-- | The arcs of an engine from the mass at the start of the first.
thrustInit :: Double -> [C.Myb_thrust_arc] -> Result C.Myb_thrust
thrustInit mass0 arcs = pureCall $ withArrayLen arcs $ \n pa -> fillZero (\pt -> C.myb_thrust_init pt (coerce mass0) pa (toSize n))

-- | The mass at t (kg).
thrustMass :: C.Myb_thrust -> Double -> Double
thrustMass thrust t = coerce (pureCall (with thrust $ \p -> C.myb_thrust_mass p (coerce t)))

-- | The same through the form of a mass function of an object (myb_thrust_mass_at).
thrustMassAt :: C.Myb_thrust -> Double -> Double
thrustMassAt thrust t = coerce (pureCall (with thrust $ \p -> C.myb_thrust_mass_at (castPtr p) (coerce t)))

-- | The acceleration of the thrust at t, x and v relative to the central body.
thrustAcceleration :: C.Myb_thrust -> Double -> Vec3 -> Vec3 -> Result Vec3
thrustAcceleration thrust t x v = pureCall $ with thrust $ \p -> withVec x $ \px -> withVec v $ \pv ->
  fillVec (C.myb_thrust_acceleration (castPtr p) (coerce t) px pv)

-- | The same in the radial, transverse, normal frame of (x, v).
thrustRtn :: C.Myb_thrust -> Double -> Vec3 -> Vec3 -> Result Vec3
thrustRtn thrust t x v = pureCall $ with thrust $ \p -> withVec x $ \px -> withVec v $ \pv ->
  fillVec (C.myb_thrust_rtn p (coerce t) px pv)

-- | The engine as a term of a force model (myb_thrust_acceleration); the model keeps it.
addThrust :: ForceModel -> String -> C.Myb_thrust -> IO (Result Int)
addThrust model name thrust = do
  memory <- mallocForeignPtr
  withForeignPtr memory (`poke` thrust)
  forceModelKeep model memory
  addRaw model name C.p'myb_thrust_acceleration nullFunPtr (castPtr (rawPtr memory))

-- | The rates of the modified equinoctial elements under the perturbing acceleration (R, T, N).
equinoctialRates :: Double -> C.Myb_equinoctial -> Vec3 -> [Double]
equinoctialRates mu elements a = pureCall $ with elements $ \pe -> withVec a $ \pa -> allocaArray 6 $ \rates ->
  C.myb_equinoctial_rates (coerce mu) pe pa rates >> peekDoubles 6 rates

-- | Edelbaum's transfer between circles of speeds v0, v1 with the plane change di: dV and the
-- initial yaw beta0.
edelbaum :: Double -> Double -> Double -> Result (Double, Double)
edelbaum v0 v1 di = pureCall $ fmap (\(a, b) -> (coerce a, coerce b)) <$> fill2 (C.myb_edelbaum (coerce v0) (coerce v1) (coerce di))

-- | The speed and the inclination change after dv along Edelbaum's law.
edelbaumState :: Double -> Double -> Double -> (Double, Double)
edelbaumState v0 beta0 dv = pureCall $ alloca $ \v -> alloca $ \di -> do
  C.myb_edelbaum_state (coerce v0) (coerce beta0) (coerce dv) v di
  (\a b -> (coerce a, coerce b)) <$> peek v <*> peek di

-- | The unit direction of the thrust (R, T, N) for the orbit, Q, its rate and the effectivity.
qlawDirection :: C.Myb_qlaw -> Elements -> Result (Vec3, Double, Double, Double)
qlawDirection law now = pureCall $ with law $ \pl -> with (elementsToC now) $ \pe -> alloca $ \d -> alloca $ \q -> alloca $ \rate -> alloca $ \eff -> do
  code <- C.myb_qlaw_direction pl pe (castPtr d) q rate eff
  andThen code (Right <$> ((,,,) <$> peek d <*> (coerce <$> peek q) <*> (coerce <$> peek rate) <*> (coerce <$> peek eff)))

-- | Whether every targeted element is within its tolerance (five values).
qlawArrived :: C.Myb_qlaw -> Elements -> [Double] -> Result Bool
qlawArrived law now tolerance
  | length tolerance /= 5 = failure Domain
  | otherwise = pureCall $ with law $ \pl -> with (elementsToC now) $ \pe -> withDoubles tolerance $ \pt ->
      Right . (/= 0) <$> C.myb_qlaw_arrived pl pe pt

withTerm :: C.Myb_qlaw -> C.Myb_qlaw_term -> (Ptr C.Myb_qlaw_term -> IO a) -> IO a
withTerm law term k = with law $ \pl -> with term {C.myb_qlaw_term'law = pl} k

-- | The acceleration of the engine of the Q-law at t (myb_term_qlaw).
qlawAcceleration :: C.Myb_qlaw -> C.Myb_qlaw_term -> Double -> Vec3 -> Vec3 -> Result Vec3
qlawAcceleration law term t x v = pureCall $ withTerm law term $ \pt -> withVec x $ \px -> withVec v $ \pv ->
  fillVec (C.myb_term_qlaw (castPtr pt) (coerce t) px pv)

-- | The engine of the Q-law as a term of a force model; the model keeps the law and the term.
addQlaw :: ForceModel -> String -> C.Myb_qlaw -> C.Myb_qlaw_term -> IO (Result Int)
addQlaw model name law term = do
  pl <- mallocForeignPtr
  withForeignPtr pl (`poke` law)
  pt <- mallocForeignPtr
  withForeignPtr pt (`poke` term {C.myb_qlaw_term'law = rawPtr pl})
  forceModelKeep model pl
  forceModelKeep model pt
  addRaw model name C.p'myb_term_qlaw nullFunPtr (castPtr (rawPtr pt))

withFlight :: C.Myb_qlaw -> C.Myb_qlaw_term -> Maybe ForceModel -> (Double, Double) -> (Ptr C.Myb_qlaw_flight -> IO a) -> IO a
withFlight law term field (rtol, atol) k = withTerm law term $ \pt -> withField $ \pf ->
  with C.zero {C.myb_qlaw_flight'engine = pt, C.myb_qlaw_flight'field = pf, C.myb_qlaw_flight'rtol = coerce rtol, C.myb_qlaw_flight'atol = coerce atol} k
  where
    withField k' = maybe (k' nullPtr) (\m -> withForeignPtr (forceModelPtr m) k') field

-- | The right-hand side of the flight (7 components: x, v, the mass) at t.
qlawFlightRhs :: C.Myb_qlaw -> C.Myb_qlaw_term -> Maybe ForceModel -> Double -> [Double] -> IO (Result [Double])
qlawFlightRhs law term field t y =
  requireLength 7 y $ withFlight law term field (0, 0) $ \pf -> withDoubles y $ \py -> fillDoubles 7 (C.myb_qlaw_flight_rhs (castPtr pf) (coerce t) py)

-- | Flies the law by DOP853 (its tolerances rtol, atol; 0: the defaults) from y0 (7) at t0
-- until the orbit is within the tolerances, looked at every step seconds, or until tLimit.
qlawFly :: C.Myb_qlaw -> C.Myb_qlaw_term -> Maybe ForceModel -> (Double, Double) -> Double -> [Double] -> Double -> Double -> IO (Result C.Myb_qlaw_arrival)
qlawFly law term field tolerances t0 y0 tLimit step =
  requireLength 7 y0 $ withFlight law term field tolerances $ \pf -> withDoubles y0 $ \py ->
    fillZero (C.myb_qlaw_fly pf (coerce t0) py (coerce tLimit) (coerce step) nullPtr)
