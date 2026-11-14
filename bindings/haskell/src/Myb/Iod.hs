-- SPDX-License-Identifier: BSD-3-Clause
-- | Initial orbit determination (myb/iod.h): from three positions (Gibbs, Herrick-Gibbs) and
-- from three directions (Gauss, Laplace, the double-r iteration, Gooding). Pure.
module Myb.Iod
  ( gibbs
  , herrickGibbs
  , Sighting (..)
  , C.Myb_iod_solution (..)
  , gauss
  , laplace
  , doubleR
  , gooding
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CSize)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | The velocity at the second of three positions of one orbit, and the angle of the first
-- out of the plane of the other two.
gibbs :: Double -> Vec3 -> Vec3 -> Vec3 -> Result (Vec3, Double)
gibbs mu r1 r2 r3 = pureCall $ withVec r1 $ \p1 -> withVec r2 $ \p2 -> withVec r3 $ \p3 ->
  fmap (fmap coerce) <$> fill2 (\v2 angle -> C.myb_iod_gibbs (coerce mu) p1 p2 p3 (castPtr v2) (angle :: Ptr CDouble))

-- | The velocity at the second of three close positions at the times t1 < t2 < t3.
herrickGibbs :: Double -> (Vec3, Double) -> (Vec3, Double) -> (Vec3, Double) -> Result Vec3
herrickGibbs mu (r1, t1) (r2, t2) (r3, t3) = pureCall $ withVec r1 $ \p1 -> withVec r2 $ \p2 -> withVec r3 $ \p3 ->
  fillVec (C.myb_iod_herrick_gibbs (coerce mu) p1 (coerce t1) p2 (coerce t2) p3 (coerce t3))

-- | Three observations: the times (increasing), the positions of the observers relative to the
-- centre and the unit directions to the object.
data Sighting = Sighting {sightTimes :: !(Double, Double, Double), sightObservers :: !Mat3, sightDirections :: !Mat3}
  deriving (Eq, Show)

withSighting :: Sighting -> (Ptr CDouble -> Ptr CDouble -> Ptr CDouble -> IO a) -> IO a
withSighting (Sighting (t1, t2, t3) observers sights) k =
  withDoubles [t1, t2, t3] $ \pt -> withVec observers $ \po -> withVec sights $ \ps -> k pt po ps

many' :: (Ptr C.Myb_iod_solution -> Ptr CSize -> IO C.Myb_status) -> IO (Result [C.Myb_iod_solution])
many' f = allocaArray C.MYB_IOD_MAX_SOLUTIONS $ \solutions -> alloca $ \count -> do
  poke count 0
  code <- f solutions count
  andThen code (peek count >>= \k -> Right <$> peekArray (fromSize k) solutions)

-- | Gauss's method, every candidate refined; up to three orbits.
gauss :: Double -> Sighting -> Result [C.Myb_iod_solution]
gauss mu s = pureCall $ withSighting s $ \pt po ps -> many' (C.myb_iod_gauss (coerce mu) pt po ps)

-- | Laplace's method, every candidate refined.
laplace :: Double -> Sighting -> Result [C.Myb_iod_solution]
laplace mu s = pureCall $ withSighting s $ \pt po ps -> many' (C.myb_iod_laplace (coerce mu) pt po ps)

-- | The double-r iteration from guesses of the distances from the centre at the first two times.
doubleR :: Double -> Sighting -> (Double, Double) -> Result C.Myb_iod_solution
doubleR mu s (d1, d2) = pureCall $ withSighting s $ \pt po ps -> withDoubles [d1, d2] $ \pd ->
  fillZero (C.myb_iod_double_r (coerce mu) pt po ps pd)

-- | Gooding's method from guesses of the ranges along the first and the third sight, with the
-- whole revolutions and the branch (-1, +1).
gooding :: Double -> Sighting -> (Double, Double) -> Int -> Int -> Result C.Myb_iod_solution
gooding mu s (r1, r3) revolutions branch = pureCall $ withSighting s $ \pt po ps -> withDoubles [r1, r3] $ \pr ->
  fillZero (C.myb_iod_gooding (coerce mu) pt po ps pr (fromIntegral revolutions) (fromIntegral branch))
