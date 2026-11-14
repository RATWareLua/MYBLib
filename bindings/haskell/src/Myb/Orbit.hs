-- SPDX-License-Identifier: BSD-3-Clause
-- | One object under a force model (myb/orbit.h): states at the times asked for, events,
-- impulses, the transition matrix. An orbit keeps its model alive.
module Myb.Orbit
  ( Orbit
  , Integrator (..)
  , OrbitOptions (..)
  , defaultOrbitOptions
  , Direction (..)
  , newOrbit
  , start
  , stateAt
  , statesAt
  , nextEvent
  , kick
  , partials
  , orbitTime
  , orbitSteps
  , orbitPtr
  ) where

import Data.Coerce (coerce)
import Data.Int (Int64)
import Foreign.C.Types (CDouble (..), CInt)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Force (ForceModel, forceModelPtr)
import Myb.Internal
import Myb.Types

-- | An orbit: its object and what it keeps (its model).
data Orbit = Orbit (ForeignPtr C.Myb_orbit) Kept

data Integrator = Ias15 | Dop853 | Rkn12
  deriving (Eq, Show)

-- | How an orbit is integrated; zeros stand for the defaults of the library.
data OrbitOptions = OrbitOptions
  { -- | The precision of IAS15 (0: 1e-10).
    epsilon :: !Double
  , -- | The largest step in s (0: unlimited).
    hMax :: !Double
  , -- | The budget of steps of one call (0: 10^9).
    maxSteps :: !Int64
  , -- | Integrate the transition matrix too, for 'partials'.
    variational :: !Bool
  , integrator :: !Integrator
  , -- | DOP853: the relative tolerance (0: 1e-12) and the absolute one (0: 1e-9).
    rtol :: !Double
  , atol :: !Double
  }
  deriving (Eq, Show)

defaultOrbitOptions :: OrbitOptions
defaultOrbitOptions = OrbitOptions 0 0 0 False Ias15 0 0

-- | Which zeros of an event function count.
data Direction = Any | Rising | Falling
  deriving (Eq, Show)

newOrbit :: ForceModel -> OrbitOptions -> IO (Result Orbit)
newOrbit model options = do
  kept <- newKept
  let fm = forceModelPtr model
  keepAlive kept fm
  with made $ \o -> withForeignPtr fm $ \m ->
    create
      (\out -> check <$> C.myb_orbit_create m o nullPtr out)
      (\p -> (`Orbit` kept) <$> owned p C.myb_orbit_destroy kept)
  where
    made =
      C.zero
        { C.myb_orbit_options'epsilon = coerce (epsilon options)
        , C.myb_orbit_options'h_max = coerce (hMax options)
        , C.myb_orbit_options'max_steps = maxSteps options
        , C.myb_orbit_options'variational = if variational options then 1 else 0
        , C.myb_orbit_options'integrator = case integrator options of
            Ias15 -> C.MYB_ORBIT_IAS15
            Dop853 -> C.MYB_ORBIT_DOP853
            Rkn12 -> C.MYB_ORBIT_RKN12
        , C.myb_orbit_options'rtol = coerce (rtol options)
        , C.myb_orbit_options'atol = coerce (atol options)
        }

-- | The object is at the state at t0 and will be asked about up to tEnd (before t0: backwards).
start :: Orbit -> Double -> State -> Double -> IO (Result ())
start (Orbit orbit _) t0 y0 tEnd =
  withForeignPtr orbit $ \o -> with y0 $ \y -> check <$> C.myb_orbit_start o (coerce t0) (castPtr y) (coerce tEnd)

-- | The state at t.
stateAt :: Orbit -> Double -> IO (Result State)
stateAt (Orbit orbit _) t =
  withForeignPtr orbit $ \o -> alloca $ \out -> do
    code <- C.myb_orbit_state o (coerce t) (castPtr out)
    either (pure . Left) (const (Right <$> peek out)) (check code)

-- | The states at times in the order of the propagation.
statesAt :: Orbit -> [Double] -> IO (Result [State])
statesAt (Orbit orbit _) times =
  withForeignPtr orbit $ \o -> withArrayLen (coerce times :: [CDouble]) $ \n ts -> allocaArray n $ \out -> do
    code <- C.myb_orbit_states o ts (toSize n) (castPtr out)
    either (pure . Left) (const (Right <$> peekArray n out)) (check code)

-- | The next zero of g(t, state) up to tLimit: its time and the state there, or 'Nothing' if
-- there is none before tLimit. The function is called from C during the call alone.
nextEvent :: Orbit -> (Double -> State -> IO Double) -> Double -> Direction -> IO (Result (Maybe (Double, State)))
nextEvent (Orbit orbit _) g tLimit direction = do
  function <- C.mk'myb_event_function call
  found <-
    withForeignPtr orbit $ \o -> alloca $ \hit -> alloca $ \tHit -> alloca $ \yHit -> do
      poke hit (0 :: CInt)
      code <- C.myb_orbit_event o function nullPtr (coerce tLimit) side hit tHit (castPtr yHit)
      case check code of
        Left e -> pure (Left e)
        Right () -> do
          any' <- peek hit
          if any' == 0
            then pure (Right Nothing)
            else (\t y -> Right (Just (coerce (t :: CDouble), y))) <$> peek tHit <*> peek yHit
  freeHaskellFunPtr function
  pure found
  where
    side = case direction of
      Any -> 0
      Rising -> 1
      Falling -> -1
    call _ t y = do
      now <- peek (castPtr y)
      coerce <$> g (coerce t) now

-- | An impulse: the velocity changes by dv where the object is now.
kick :: Orbit -> Vec3 -> IO (Result ())
kick (Orbit orbit _) dv = withForeignPtr orbit $ \o -> with dv $ \d -> check <$> C.myb_orbit_kick o (castPtr d)

-- | The state at t and the transition matrix from the start, 36 numbers by rows
-- ('variational' must be set).
partials :: Orbit -> Double -> IO (Result (State, [Double]))
partials (Orbit orbit _) t =
  withForeignPtr orbit $ \o -> alloca $ \y -> allocaArray 36 $ \phi -> do
    code <- C.myb_orbit_partials o (coerce t) (castPtr y) phi nullPtr
    case check code of
      Left e -> pure (Left e)
      Right () -> (\s m -> Right (s, coerce (m :: [CDouble]))) <$> peek y <*> peekArray 36 phi

-- | The time the propagation has reached.
orbitTime :: Orbit -> IO Double
orbitTime (Orbit orbit _) = coerce <$> withForeignPtr orbit C.myb_orbit_time

orbitSteps :: Orbit -> IO Int64
orbitSteps (Orbit orbit _) = withForeignPtr orbit C.myb_orbit_steps

orbitPtr :: Orbit -> ForeignPtr C.Myb_orbit
orbitPtr (Orbit orbit _) = orbit
