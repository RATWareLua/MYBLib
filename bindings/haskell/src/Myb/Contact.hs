-- SPDX-License-Identifier: BSD-3-Clause
-- | Contact (myb/contact.h): the surface nearest to a point, the force and the torque of a
-- surface on a body through its contact points — with Coulomb's friction made smooth, or the
-- friction that sticks, a spring of its own whose deflections are integrated with the motion —
-- and the damping that gives a restitution. The surface, the law and the points are the records
-- of "Myb.C"; a mesh a surface points to is an object of "Myb.C" the caller keeps alive.
-- Functions of numbers, so pure. In a system: 'Myb.System.addContact', 'Myb.System.packContact'.
module Myb.Contact
  ( C.Myb_contact_surface (..)
  , C.Myb_contact_law (..)
  , C.Myb_contact_point (..)
  , C.Myb_contact_result (..)
  , Friction (..)
  , frictionCode
  , contactSurfaceNearest
  , contactForce
  , contactForceStick
  , contactRestitution
  , contactDamping
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr)
import Foreign.Storable (peek)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | The friction of a law.
data Friction
  = -- | Coulomb's made smooth: a function of the state (the default)
    FrictionSmooth
  | -- | a spring with a slip limit: a state of its own
    FrictionStick
  deriving (Eq, Show)

-- | The value of the field `kind` of a law: @zero {myb_contact_law'kind = frictionCode FrictionStick, ...}@.
frictionCode :: Friction -> C.Myb_contact_friction
frictionCode FrictionSmooth = C.MYB_CONTACT_FRICTION_SMOOTH
frictionCode FrictionStick = C.MYB_CONTACT_FRICTION_STICK

-- | The point of the surface nearest to a position, the outward unit normal there and the height
-- of the position over it (negative inside).
contactSurfaceNearest :: C.Myb_contact_surface -> Vec3 -> Result (Vec3, Vec3, Double)
contactSurfaceNearest surface position' = pureCall $ with surface $ \ps -> withVec position' $ \pp -> alloca $ \point -> alloca $ \normal -> alloca $ \height -> do
  code <- C.myb_contact_surface_nearest ps pp (castPtr point) (castPtr normal) height
  andThen code (Right <$> ((,,) <$> peek point <*> peek normal <*> (coerce <$> peek height)))

-- | The force (per unit mass, axes of the surface) and the torque (per unit mass, body axes) of
-- the surface on a body: its contact points, the place and velocity of its centre of mass, its
-- attitude (surface -> body) and angular velocity, all in the frame of the surface.
contactForce :: C.Myb_contact_surface -> C.Myb_contact_law -> [C.Myb_contact_point] -> Vec3 -> Vec3 -> Mat3 -> Vec3 -> Result C.Myb_contact_result
contactForce surface law points x v attitude' omega = pureCall $
  with surface $ \ps -> with law $ \pl -> withArrayLen points $ \n pp -> withVec x $ \px -> withVec v $ \pv -> withVec attitude' $ \pa -> withVec omega $ \pw ->
    fillZero (C.myb_contact_force ps pl pp (toSize n) px pv pa pw)

-- | The same with the state of a friction that sticks: the deflection of every point (m, axes
-- of the surface); with the result, the rates of the deflections, to be integrated with the
-- motion. With the smooth friction the deflections are not read.
contactForceStick :: C.Myb_contact_surface -> C.Myb_contact_law -> [C.Myb_contact_point] -> Vec3 -> Vec3 -> Mat3 -> Vec3 -> [Vec3] -> Result (C.Myb_contact_result, [Vec3])
contactForceStick surface law points x v attitude' omega deflections
  | length deflections /= length points = failure Domain
  | otherwise = pureCall $
      with surface $ \ps -> with law $ \pl -> withArrayLen points $ \n pp -> withVec x $ \px -> withVec v $ \pv -> withVec attitude' $ \pa -> withVec omega $ \pw ->
        withVecs deflections $ \pd -> allocaArray (max 1 n) $ \rates -> with C.zero $ \result -> do
          code <- C.myb_contact_force_stick ps pl pp (toSize n) px pv pa pw pd (castPtr rates) result
          andThen code (Right <$> ((,) <$> peek result <*> peekArray n rates))

-- | The restitution e of an impact at the approach speed v (m/s) for the damping chi (s/m).
contactRestitution :: Double -> Double -> Result Double
contactRestitution damping speed = pureCall $ fmap coerce <$> fill (C.myb_contact_restitution (coerce damping) (coerce speed))

-- | The damping chi (s/m) that gives the restitution e at the approach speed v (m/s).
contactDamping :: Double -> Double -> Result Double
contactDamping restitution speed = pureCall $ fmap coerce <$> fill (C.myb_contact_damping (coerce restitution) (coerce speed))
