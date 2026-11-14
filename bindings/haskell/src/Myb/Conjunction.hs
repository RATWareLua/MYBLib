-- SPDX-License-Identifier: BSD-3-Clause
-- | Close approaches and the probability of collision (myb/conjunction.h): approaches of two
-- providers or of two orbits being propagated, the plane of the encounter, the probability by
-- integration, by Chan's series and its largest value, approaches with their covariances
-- carried, and long encounters.
module Myb.Conjunction
  ( C.Myb_approach (..)
  , approachFind
  , approachOrbits
  , C.Myb_encounter (..)
  , encounterPlane
  , collisionProbability
  , collisionProbabilityChan
  , collisionProbabilityMax
  , C.Myb_conjunction (..)
  , conjunctionFind
  , collisionRate
  , collisionProbabilityLong
  , conjunctionLong
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CSize)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, pokeArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (Storable, peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Orbit (Orbit, orbitPtr)
import Myb.Types

doubles :: (CDouble, CDouble) -> (Double, Double)
doubles (a, b) = (coerce a, coerce b)

-- | Up to `capacity` results written, and how many were found (it may be more).
found :: Storable a => Int -> (Ptr a -> Ptr CSize -> IO C.Myb_status) -> IO (Result ([a], Int))
found capacity f = allocaArray (max 1 capacity) $ \items -> alloca $ \count -> do
  poke count 0
  code <- f items count
  andThen code $ do
    k <- fromSize <$> peek count
    xs <- peekArray (min k capacity) items
    pure (Right (xs, k))

-- | The closest approaches of two providers over [t0, t1], looked for in steps of `step`.
approachFind :: StateSource -> StateSource -> Double -> Double -> Double -> Int -> IO (Result ([C.Myb_approach], Int))
approachFind first other t0 t1 step capacity =
  withStateSource first $ \f1 c1 -> withStateSource other $ \f2 c2 ->
    found capacity (\items count -> C.myb_approach_find f1 c1 f2 c2 (coerce t0) (coerce t1) (coerce step) items (toSize capacity) count)

-- | The closest approaches of two orbits being propagated, standing at one time, up to t1;
-- both are left at t1.
approachOrbits :: Orbit -> Orbit -> Double -> Double -> Int -> IO (Result ([C.Myb_approach], Int))
approachOrbits first other t1 step capacity =
  withForeignPtr (orbitPtr first) $ \o1 -> withForeignPtr (orbitPtr other) $ \o2 ->
    found capacity (\items count -> C.myb_approach_orbits o1 o2 (coerce t1) (coerce step) items (toSize capacity) count)

-- | The plane of the encounter from the relative position and velocity at the closest
-- approach, the covariance of the relative position and the sum of the radii.
encounterPlane :: Vec3 -> Vec3 -> Mat3 -> Double -> Result C.Myb_encounter
encounterPlane position' velocity' covariance radius = pureCall $
  withVec position' $ \pp -> withVec velocity' $ \pv -> withVec covariance $ \pc -> fill (C.myb_encounter_plane pp pv pc (coerce radius))

-- | The integral of the Gaussian over the disk.
collisionProbability :: C.Myb_encounter -> Result Double
collisionProbability e = pureCall $ with e $ \pe -> fmap coerce <$> fill (C.myb_collision_probability pe)

-- | The series of Chan (NaN for an invalid encounter).
collisionProbabilityChan :: C.Myb_encounter -> Double
collisionProbabilityChan e = coerce (pureCall (with e C.myb_collision_probability_chan))

-- | The largest probability for a covariance known in shape only, and the factor of it.
collisionProbabilityMax :: C.Myb_encounter -> Result (Double, Double)
collisionProbabilityMax e = pureCall $ with e $ \pe -> fmap (\(a, b) -> (coerce a, coerce b)) <$> fill2 (C.myb_collision_probability_max pe)

-- | The approaches of two orbits (created with variational) with the covariances of their
-- states at the time they stand (36 each, by rows) carried to every approach.
conjunctionFind :: Orbit -> [Double] -> Orbit -> [Double] -> Double -> Double -> Double -> Int -> IO (Result ([C.Myb_conjunction], Int))
conjunctionFind first p1 other p2 radius t1 step capacity =
  requireLength 36 p1 $ requireLength 36 p2 $ withForeignPtr (orbitPtr first) $ \o1 -> withForeignPtr (orbitPtr other) $ \o2 -> withDoubles p1 $ \c1 -> withDoubles p2 $ \c2 ->
    found capacity (\items count -> C.myb_conjunction_find o1 c1 o2 c2 (coerce radius) (coerce t1) (coerce step) items (toSize capacity) count)

-- | The rate of entry (1/s) into the sphere of the radius for the relative state with its
-- mean and covariance (36, by rows).
collisionRate :: State -> [Double] -> Double -> Result Double
collisionRate mean covariance radius
  | length covariance /= 36 = failure Domain
  | otherwise = pureCall $ withVec mean $ \pm -> withDoubles covariance $ \pc -> fmap coerce <$> fill (C.myb_collision_rate pm pc (coerce radius))

-- | The expected number of entries over [t0, t1] for the relative state given as a Gaussian by
-- a Haskell function of the time (its mean and covariance), and the error estimate.
collisionProbabilityLong :: (Double -> IO (Either Status (State, [Double]))) -> Double -> Double -> Double -> Double -> IO (Result (Double, Double))
collisionProbabilityLong gaussian t0 t1 step radius = do
  function <- C.mk'myb_relative_gaussian_function act
  done <- fmap doubles <$> fill2 (C.myb_collision_probability_long function nullPtr (coerce t0) (coerce t1) (coerce step) (coerce radius))
  freeHaskellFunPtr function
  pure done
  where
    act _ t mean covariance = do
      given <- gaussian (coerce t)
      case given of
        Left status -> pure (statusCode status)
        Right (m, c)
          | length c /= 36 -> pure C.MYB_E_DOMAIN
          | otherwise -> poke (castPtr mean) m >> pokeArray covariance (coerce c) >> pure C.MYB_OK

-- | The expected number of entries for two orbits (with variational) from the time they stand
-- to t1, and the error estimate; both are left at t1.
conjunctionLong :: Orbit -> [Double] -> Orbit -> [Double] -> Double -> Double -> Double -> IO (Result (Double, Double))
conjunctionLong first p1 other p2 radius t1 step =
  requireLength 36 p1 $ requireLength 36 p2 $ withForeignPtr (orbitPtr first) $ \o1 -> withForeignPtr (orbitPtr other) $ \o2 -> withDoubles p1 $ \c1 -> withDoubles p2 $ \c2 ->
    fmap doubles <$> fill2 (C.myb_conjunction_long o1 c1 o2 c2 (coerce radius) (coerce t1) (coerce step))
