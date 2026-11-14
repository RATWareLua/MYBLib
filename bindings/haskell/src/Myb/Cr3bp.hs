-- SPDX-License-Identifier: BSD-3-Clause
-- | The circular restricted three-body problem (myb/cr3bp.h), in its units: the distance of the
-- primaries, their total mass and their mean motion are 1, mu = m2 / (m1 + m2). Lagrange points,
-- the Jacobi constant, periodic orbits and their families, monodromy and stability, manifolds,
-- quasi-periodic tori with their stability and manifolds, propagation to a section, the frames
-- of real bodies, and the problem as a term of a force model. Functions of numbers alone are
-- pure; the orbits, tori and spectra are the records of "Myb.C".
module Myb.Cr3bp
  ( cr3bpAcceleration
  , cr3bpJacobian
  , jacobi
  , lagrange
  , addCr3bp
    -- * Periodic orbits symmetric about a plane or an axis
  , Symmetry (..)
  , Correction (..)
  , correction
  , correct
  , monodromy
  , stability
  , manifold
  , richardson
    -- * Frames
  , cr3bpFromEphemeris
  , cr3bpToEphemeris
    -- * Periodic orbits of any shape and their families
  , C.Myb_cr3bp_orbit (..)
  , orbitFromState
  , orbitCorrect
  , orbitContinue
    -- * Tori
  , C.Myb_cr3bp_torus (..)
  , C.Myb_cr3bp_torus_spectrum (..)
  , TorusFamily (..)
  , torusFromOrbit
  , torusContinue
  , torusState
  , torusStability
  , torusManifold
    -- * Propagation
  , Section (..)
  , propagate
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CInt)
import Foreign.ForeignPtr (mallocForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Force (ForceModel, addRaw, forceModelKeep)
import Myb.Internal
import Myb.Types

-- | The acceleration in the rotating frame.
cr3bpAcceleration :: Double -> Vec3 -> Vec3 -> Vec3
cr3bpAcceleration mu x v = pureCall $ withVec x $ \px -> withVec v $ \pv -> alloca $ \a ->
  C.myb_cr3bp_acceleration (coerce mu) px pv (castPtr a) >> peek a

-- | d a / d x and d a / d v.
cr3bpJacobian :: Double -> Vec3 -> (Mat3, Mat3)
cr3bpJacobian mu x = pureCall $ withVec x $ \px -> alloca $ \dx -> alloca $ \dv ->
  C.myb_cr3bp_jacobian (coerce mu) px (castPtr dx) (castPtr dv) >> ((,) <$> peek dx <*> peek dv)

-- | The Jacobi constant C = 2 Omega - |v|^2.
jacobi :: Double -> State -> Double
jacobi mu y = coerce (pureCall (withVec y (C.myb_cr3bp_jacobi (coerce mu))))

-- | The five Lagrange points.
lagrange :: Double -> Result [Vec3]
lagrange mu = pureCall (fillArray 5 (C.myb_cr3bp_lagrange (coerce mu) . castPtr))

-- | The problem as a term of a force model (myb_cr3bp_term with its Jacobian); the model keeps mu.
addCr3bp :: ForceModel -> String -> Double -> IO (Result Int)
addCr3bp model name mu = do
  context <- mallocForeignPtr
  withForeignPtr context (`poke` (coerce mu :: CDouble))
  forceModelKeep model context
  addRaw model name C.p'myb_cr3bp_term C.p'myb_cr3bp_term_jacobian (castPtr (rawPtr context))

data Symmetry
  = -- | y0 = (x, 0, z, 0, vy, 0): halo, Lyapunov, distant retrograde
    SymmetryXZ
  | -- | y0 = (x, 0, 0, 0, vy, vz): vertical, axial
    SymmetryX
  deriving (Eq, Show)

-- | How a symmetric orbit is corrected: the component of y0 held fixed, the crossing that ends
-- the half period (1: the first), the tolerance (0: 1e-12) and the iterations (0: 30).
data Correction = Correction
  { correctionSymmetry :: !Symmetry
  , correctionFixed :: !Int
  , correctionCrossing :: !Int
  , correctionTolerance :: !Double
  , correctionMaxIterations :: !Int
  }
  deriving (Eq, Show)

-- | The first crossing, the defaults of the library, the component given held fixed.
correction :: Symmetry -> Int -> Correction
correction symmetry fixed = Correction symmetry fixed 1 0 0

-- | The periodic orbit through y0 (corrected) and its period.
correct :: Double -> Correction -> State -> Result (State, Double)
correct mu c y0 = pureCall $ with made $ \pc -> with y0 $ \py -> alloca $ \period -> do
  code <- C.myb_cr3bp_correct (coerce mu) pc (castPtr py) period
  andThen code (Right <$> ((,) <$> peek py <*> (coerce <$> peek period)))
  where
    made =
      C.zero
        { C.myb_cr3bp_correction'symmetry = case correctionSymmetry c of
            SymmetryXZ -> C.MYB_CR3BP_SYMMETRY_XZ
            SymmetryX -> C.MYB_CR3BP_SYMMETRY_X
        , C.myb_cr3bp_correction'fixed = fromIntegral (correctionFixed c)
        , C.myb_cr3bp_correction'crossing = fromIntegral (correctionCrossing c)
        , C.myb_cr3bp_correction'tolerance = coerce (correctionTolerance c)
        , C.myb_cr3bp_correction'max_iterations = fromIntegral (correctionMaxIterations c)
        }

-- | The monodromy matrix (36 numbers by rows).
monodromy :: Double -> State -> Double -> Result [Double]
monodromy mu y0 period = pureCall $ withVec y0 $ \py -> fillDoubles 36 (C.myb_cr3bp_monodromy (coerce mu) py (coerce period))

-- | The eigenvalues of a monodromy matrix (real and imaginary parts, 6 each) and the stability
-- index.
stability :: [Double] -> Result ([Double], [Double], Double)
stability m
  | length m /= 36 = failure Domain
  | otherwise = pureCall $ withDoubles m $ \pm -> alloca $ \re -> alloca $ \im -> alloca $ \index -> do
      code <- C.myb_cr3bp_stability pm (castPtr re) (castPtr im) index
      andThen code (Right <$> ((,,) <$> (stateList <$> peek re) <*> (stateList <$> peek im) <*> (coerce <$> peek index)))

-- | A state on the unstable (True) or stable manifold at a phase (0 .. 1) of the period,
-- displaced by epsilon on the side sign (+1, -1).
manifold :: Double -> State -> Double -> Double -> Bool -> Int -> Double -> Result State
manifold mu y0 period phase unstable sign epsilon = pureCall $ withVec y0 $ \py ->
  fillVec (C.myb_cr3bp_manifold (coerce mu) py (coerce period) (coerce phase) (if unstable then 1 else 0) (fromIntegral sign) (coerce epsilon))

-- | The third-order halo orbit of Richardson about L1 (1) or L2 (2) with the amplitude az and
-- the branch (+1 northern, -1): the state at the crossing of the xz-plane and the period.
richardson :: Double -> Int -> Double -> Int -> Result (State, Double)
richardson mu point az branch = pureCall $ alloca $ \py -> alloca $ \period -> do
  code <- C.myb_cr3bp_richardson (coerce mu) (fromIntegral point) (coerce az) (fromIntegral branch) (castPtr py) period
  andThen code (Right <$> ((,) <$> peek py <*> (coerce <$> peek period)))

-- | A state (m, m/s) in the rotating frame of two real bodies (their states and GM).
cr3bpFromEphemeris :: State -> State -> Double -> Double -> State -> Result State
cr3bpFromEphemeris b1 b2 gm1 gm2 y = pureCall $ withVec b1 $ \p1 -> withVec b2 $ \p2 -> withVec y $ \py ->
  fillVec (C.myb_cr3bp_from_ephemeris p1 p2 (coerce gm1) (coerce gm2) py)

cr3bpToEphemeris :: State -> State -> Double -> Double -> State -> Result State
cr3bpToEphemeris b1 b2 gm1 gm2 y = pureCall $ withVec b1 $ \p1 -> withVec b2 $ \p2 -> withVec y $ \py ->
  fillVec (C.myb_cr3bp_to_ephemeris p1 p2 (coerce gm1) (coerce gm2) py)

-- | The cuts of the orbit through y0 with the period, by propagation.
orbitFromState :: Double -> State -> Double -> Int -> Result C.Myb_cr3bp_orbit
orbitFromState mu y0 period patches = pureCall $ withVec y0 $ \py ->
  fillZero (C.myb_cr3bp_orbit_from_state (coerce mu) py (coerce period) (toSize patches))

-- | The orbit corrected (tolerance 0: 1e-11; iterations 0: 30), with its tangent.
orbitCorrect :: Double -> C.Myb_cr3bp_orbit -> Double -> Int -> Result C.Myb_cr3bp_orbit
orbitCorrect mu orbit tolerance iterations = pureCall $ with orbit $ \po -> do
  code <- C.myb_cr3bp_orbit_correct (coerce mu) po (coerce tolerance) (fromIntegral iterations)
  andThen code (Right <$> peek po)

-- | The next orbit of the family a step along the tangent.
orbitContinue :: Double -> C.Myb_cr3bp_orbit -> Double -> Double -> Int -> Result C.Myb_cr3bp_orbit
orbitContinue mu orbit step tolerance iterations = pureCall $ with orbit $ \po -> do
  code <- C.myb_cr3bp_orbit_continue (coerce mu) po (coerce step) (coerce tolerance) (fromIntegral iterations)
  andThen code (Right <$> peek po)

-- | What a family of tori keeps.
data TorusFamily = KeepJacobi | KeepPeriod | KeepRotation
  deriving (Eq, Show)

-- | The smallest torus about the periodic orbit through y0: the points, the pair of
-- eigenvalues (0: the first) and the size.
torusFromOrbit :: Double -> State -> Double -> Int -> Int -> Double -> Result C.Myb_cr3bp_torus
torusFromOrbit mu y0 period points pair size = pureCall $ withVec y0 $ \py ->
  fillZero (C.myb_cr3bp_torus_from_orbit (coerce mu) py (coerce period) (toSize points) (fromIntegral pair) (coerce size))

-- | The torus a step along the tangent of the family, corrected (a step of 0 corrects it).
torusContinue :: Double -> C.Myb_cr3bp_torus -> TorusFamily -> Double -> Double -> Int -> Result C.Myb_cr3bp_torus
torusContinue mu torus family step tolerance iterations = pureCall $ with torus $ \pt -> do
  code <- C.myb_cr3bp_torus_continue (coerce mu) pt kind (coerce step) (coerce tolerance) (fromIntegral iterations)
  andThen code (Right <$> peek pt)
  where
    kind = case family of
      KeepJacobi -> C.MYB_CR3BP_TORUS_JACOBI
      KeepPeriod -> C.MYB_CR3BP_TORUS_PERIOD
      KeepRotation -> C.MYB_CR3BP_TORUS_ROTATION

-- | The point of the circle at an angle (rad).
torusState :: C.Myb_cr3bp_torus -> Double -> State
torusState torus angle = pureCall $ with torus $ \pt -> alloca $ \py ->
  C.myb_cr3bp_torus_state pt (coerce angle) (castPtr py) >> peek py

-- | The spectrum of the torus (tolerance 0: 1e-4): its eigenvalues, their classes, the
-- hyperbolic directions and the normal pair of an elliptic torus.
torusStability :: Double -> C.Myb_cr3bp_torus -> Double -> Result C.Myb_cr3bp_torus_spectrum
torusStability mu torus tolerance = pureCall $ with torus $ \pt ->
  fillZero (C.myb_cr3bp_torus_stability (coerce mu) pt (coerce tolerance))

-- | A state on the unstable (True) or stable manifold of the torus at an angle.
torusManifold :: C.Myb_cr3bp_torus -> C.Myb_cr3bp_torus_spectrum -> Double -> Bool -> Int -> Double -> Result State
torusManifold torus spectrum angle unstable sign epsilon = pureCall $ with torus $ \pt -> with spectrum $ \ps ->
  fillVec (C.myb_cr3bp_torus_manifold pt ps (coerce angle) (if unstable then 1 else 0) (fromIntegral sign) (coerce epsilon))

-- | A plane y[axis] = value crossed growing (+1), falling (-1) or either (0).
data Section = Section {sectionAxis :: !Int, sectionValue :: !Double, sectionDirection :: !Int}
  deriving (Eq, Show)

-- | The flow from y0 for a duration (negative: backward), stopping at the first crossing of
-- the section if one is given: the state there (or at the end), the time elapsed and whether
-- the plane was reached.
propagate :: Double -> State -> Double -> Maybe Section -> Result (State, Double, Bool)
propagate mu y0 duration section = pureCall $ withVec y0 $ \py0 -> withSection $ \ps -> alloca $ \py -> alloca $ \elapsed -> alloca $ \crossed -> do
  poke crossed (0 :: CInt)
  code <- C.myb_cr3bp_propagate (coerce mu) py0 (coerce duration) ps (castPtr py) elapsed crossed
  andThen code (Right <$> ((,,) <$> peek py <*> (coerce <$> peek elapsed) <*> ((/= 0) <$> peek crossed)))
  where
    withSection k = case section of
      Nothing -> k nullPtr
      Just s ->
        with
          C.zero
            { C.myb_cr3bp_section'axis = fromIntegral (sectionAxis s)
            , C.myb_cr3bp_section'value = coerce (sectionValue s)
            , C.myb_cr3bp_section'direction = fromIntegral (sectionDirection s)
            }
          k
