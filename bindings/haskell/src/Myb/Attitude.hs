-- SPDX-License-Identifier: BSD-3-Clause
-- | Attitude (myb/attitude.h): quaternions, Euler angles, principal axes, the equations of a
-- rigid body and of a body with a fluid core, gravitational torques, the exact torque-free
-- rotation and the splitting integrator. Functions of numbers alone are pure.
module Myb.Attitude
  ( Quaternion (..)
  , Mat3 (..)
  , mat3List
    -- * Quaternions
  , quatMultiply
  , quatConjugate
  , quatNormalize
  , quatToMatrix
  , quatFromMatrix
  , quatRotate
  , quatFromRotationVector
  , quatToRotationVector
  , quatSlerp
  , quatDerivative
    -- * Euler angles and inertia
  , eulerToMatrix
  , eulerFromMatrix
  , principalAxes
    -- * Equations of motion
  , rigidBodyRates
  , coreMantleRates
    -- * Torques
  , torquePointMass
  , torqueField
  , torqueFigureFigure
    -- * Torque-free rotation
  , C.Myb_free_rotation (..)
  , freeRotationInit
  , freeRotationState
  , rotationSplit
  ) where

import Data.Coerce (coerce)
import Data.Int (Int64)
import Foreign.C.Types (CDouble (..), CInt)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (Storable, peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | A value written by a function of C that cannot fail.
value :: Storable a => (Ptr a -> IO ()) -> a
value f = pureCall (alloca $ \p -> f p >> peek p)

q' :: Quaternion -> (Ptr CDouble -> IO a) -> IO a
q' = withVec

quatMultiply :: Quaternion -> Quaternion -> Quaternion
quatMultiply a b = value $ \o -> q' a $ \pa -> q' b $ \pb -> C.myb_quat_multiply pa pb (castPtr o)

quatConjugate :: Quaternion -> Quaternion
quatConjugate q = value $ \o -> q' q $ \pq -> C.myb_quat_conjugate pq (castPtr o)

-- | The quaternion made unit, and its former norm (unchanged if it was 0).
quatNormalize :: Quaternion -> (Quaternion, Double)
quatNormalize q = pureCall $ with q $ \pq -> do
  norm <- C.myb_quat_normalize (castPtr pq)
  (\u -> (u, coerce norm)) <$> peek pq

quatToMatrix :: Quaternion -> Mat3
quatToMatrix q = value $ \o -> q' q $ \pq -> C.myb_quat_to_matrix pq (castPtr o)

-- | The unit quaternion of a rotation matrix, w >= 0.
quatFromMatrix :: Mat3 -> Quaternion
quatFromMatrix m = value $ \o -> withVec m $ \pm -> C.myb_quat_from_matrix pm (castPtr o)

-- | M(q) v.
quatRotate :: Quaternion -> Vec3 -> Vec3
quatRotate q v = value $ \o -> q' q $ \pq -> withVec v $ \pv -> C.myb_quat_rotate pq pv (castPtr o)

quatFromRotationVector :: Vec3 -> Quaternion
quatFromRotationVector phi = value $ \o -> withVec phi $ \p -> C.myb_quat_from_rotation_vector p (castPtr o)

quatToRotationVector :: Quaternion -> Vec3
quatToRotationVector q = value $ \o -> q' q $ \pq -> C.myb_quat_to_rotation_vector pq (castPtr o)

-- | Spherical linear interpolation along the shorter arc, t in [0, 1].
quatSlerp :: Quaternion -> Quaternion -> Double -> Quaternion
quatSlerp a b t = value $ \o -> q' a $ \pa -> q' b $ \pb -> C.myb_quat_slerp pa pb (coerce t) (castPtr o)

-- | dq/dt for the body angular velocity.
quatDerivative :: Quaternion -> Vec3 -> Quaternion
quatDerivative q omega = value $ \o -> q' q $ \pq -> withVec omega $ \pw -> C.myb_quat_derivative pq pw (castPtr o)

axesOf :: (Int, Int, Int) -> [CInt]
axesOf (i, j, k) = map fromIntegral [i, j, k]

-- | M = R_k(c) R_j(b) R_i(a) for the axes (i, j, k) and the angles (a, b, c).
eulerToMatrix :: (Int, Int, Int) -> Vec3 -> Result Mat3
eulerToMatrix axes angles = pureCall $ withArray (axesOf axes) $ \pa -> withVec angles $ \pg -> fillVec (C.myb_euler_to_matrix pa pg)

eulerFromMatrix :: (Int, Int, Int) -> Mat3 -> Result Vec3
eulerFromMatrix axes m = pureCall $ withArray (axesOf axes) $ \pa -> withVec m $ \pm -> fillVec (C.myb_euler_from_matrix pa pm)

-- | The principal moments (ascending) and the rotation to the principal frame.
principalAxes :: Mat3 -> Result (Vec3, Mat3)
principalAxes inertia = pureCall $ withVec inertia $ \pi' -> fill2 (\pm pa -> C.myb_principal_axes pi' (castPtr pm) (castPtr pa))

-- | dq/dt and dL/dt for the principal moments, the attitude, the body angular momentum and
-- the body torque ('Nothing': free).
rigidBodyRates :: Vec3 -> Quaternion -> Vec3 -> Maybe Vec3 -> (Quaternion, Vec3)
rigidBodyRates moments q l torque = pureCall $
  withVec moments $ \pm -> q' q $ \pq -> withVec l $ \pl -> withMaybe torque $ \pt -> alloca $ \dq -> alloca $ \dl -> do
    C.myb_rigid_body_rates pm pq pl pt (castPtr dq) (castPtr dl)
    (,) <$> peek dq <*> peek dl

withMaybe :: Maybe Vec3 -> (Ptr CDouble -> IO a) -> IO a
withMaybe Nothing k = k nullPtr
withMaybe (Just v) k = withVec v k

-- | The rates of the angular velocities of the mantle and of the core: the inertia of the
-- mantle and its rate ('Nothing': constant), the core moments (A_c, C_c), k_v, the two
-- angular velocities and the external torque.
coreMantleRates :: Mat3 -> Maybe Mat3 -> (Double, Double) -> Double -> Vec3 -> Vec3 -> Vec3 -> Result (Vec3, Vec3)
coreMantleRates inertia rate (ac, cc) kv omegaMantle omegaCore torque = pureCall $
  withVec inertia $ \pi' -> withMaybeM rate $ \pr -> withDoubles [ac, cc] $ \pc -> withVec omegaMantle $ \pm -> withVec omegaCore $ \po -> withVec torque $ \pt ->
    fill2 (\rm rc -> C.myb_core_mantle_rates pi' pr pc (coerce kv) pm po pt (castPtr rm) (castPtr rc))
  where
    withMaybeM Nothing k = k nullPtr
    withMaybeM (Just m) k = withVec m k

-- | 3 gm / r^5 (r x I r) on a body of inertia diag(moments).
torquePointMass :: Double -> Vec3 -> Vec3 -> Vec3
torquePointMass gm r moments = value $ \o -> withVec r $ \pr -> withVec moments $ \pm -> C.myb_torque_point_mass (coerce gm) pr pm (castPtr o)

-- | -(m_other / m_body) r x a, per unit mass of the body.
torqueField :: Double -> Vec3 -> Vec3 -> Vec3
torqueField ratio r a = value $ \o -> withVec r $ \pr -> withVec a $ \pa -> C.myb_torque_field (coerce ratio) pr pa (castPtr o)

-- | The torque per unit mass of A from the quadrupole-quadrupole potential.
torqueFigureFigure :: Double -> Mat3 -> Mat3 -> Vec3 -> Result Vec3
torqueFigureFigure gmA gqA gqB r = pureCall $
  withVec gqA $ \pa -> withVec gqB $ \pb -> withVec r $ \pr -> fillVec (C.myb_torque_figure_figure (coerce gmA) pa pb pr)

-- | The exact torque-free rotation from the principal moments, the attitude and the body
-- angular momentum at t = 0.
freeRotationInit :: Vec3 -> Quaternion -> Vec3 -> Result C.Myb_free_rotation
freeRotationInit moments q l = pureCall $
  withVec moments $ \pm -> q' q $ \pq -> withVec l $ \pl -> fill (C.myb_free_rotation_init pm pq pl)

-- | The attitude and the body angular momentum at t.
freeRotationState :: C.Myb_free_rotation -> Double -> (Quaternion, Vec3)
freeRotationState rotation t = pureCall $
  with rotation $ \pr -> alloca $ \pq -> alloca $ \pl -> do
    C.myb_free_rotation_state pr (coerce t) (castPtr pq) (castPtr pl)
    (,) <$> peek pq <*> peek pl

-- | Advances (q, l) from t0 by steps of h with the splitting integrator of the order given
-- (2, 4, 6, 8), the torque a Haskell function of the time and the attitude.
rotationSplit :: Vec3 -> (Double -> Quaternion -> IO (Either Status Vec3)) -> Int -> Double -> Double -> Int64 -> Quaternion -> Vec3 -> IO (Result (Quaternion, Vec3))
rotationSplit moments torque order t0 h steps q l = do
  function <- C.mk'myb_torque_function act
  done <- withVec moments $ \pm -> with q $ \pq -> with l $ \pl -> do
    code <- C.myb_rotation_split pm function nullPtr (fromIntegral order) (coerce t0) (coerce h) steps (castPtr pq) (castPtr pl)
    andThen code (Right <$> ((,) <$> peek pq <*> peek pl))
  freeHaskellFunPtr function
  pure done
  where
    act _ t pq out' = do
      attitude <- peek (castPtr pq)
      given <- torque (coerce t) attitude
      case given of
        Right n -> poke (castPtr out') n >> pure C.MYB_OK
        Left status -> pure (statusCode status)
