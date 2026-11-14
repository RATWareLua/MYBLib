-- SPDX-License-Identifier: BSD-3-Clause
-- | Targeting (myb/transfer.h, myb/lambert.h): the B-plane, patched conics, Newton's method
-- for the shooting of trajectories, and Lambert's problem.
module Myb.Transfer
  ( C.Myb_bplane (..)
  , bplane
  , sphereOfInfluence
  , hillRadius
  , flybyTurn
  , departureImpulse
    -- * Shooting
  , ShootOptions (..)
  , Shot (..)
  , shoot
  , shootWorkSize
    -- * Lambert's problem
  , C.Myb_lambert_solution (..)
  , lambert
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CInt)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, pokeArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (freeHaskellFunPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | The B-plane of the two-body hyperbola through (r, v) about a target of the GM given, with
-- its pole.
bplane :: Double -> Vec3 -> Vec3 -> Vec3 -> Result C.Myb_bplane
bplane gm r v pole = pureCall $ withVec r $ \pr -> withVec v $ \pv -> withVec pole $ \pp -> fill (C.myb_bplane_compute (coerce gm) pr pv pp)

-- | a (m/M)^(2/5).
sphereOfInfluence :: Double -> Double -> Double
sphereOfInfluence a ratio = coerce (pureCall (C.myb_sphere_of_influence (coerce a) (coerce ratio)))

-- | a (1 - e) (m / (3M))^(1/3).
hillRadius :: Double -> Double -> Double -> Double
hillRadius a e ratio = coerce (pureCall (C.myb_hill_radius (coerce a) (coerce e) (coerce ratio)))

-- | The turn angle of a flyby: gm, the excess speed, the periapsis radius.
flybyTurn :: Double -> Double -> Double -> Double
flybyTurn gm v periapsis = coerce (pureCall (C.myb_flyby_turn (coerce gm) (coerce v) (coerce periapsis)))

-- | The impulse from a circular orbit of a radius to a hyperbola of an excess speed.
departureImpulse :: Double -> Double -> Double -> Double
departureImpulse gm radius v = coerce (pureCall (C.myb_departure_impulse (coerce gm) (coerce radius) (coerce v)))

-- | The tolerance on max |f_i| (> 0) and the iterations (0: 20).
data ShootOptions = ShootOptions {shootTolerance :: !Double, shootMaxIterations :: !Int}
  deriving (Eq, Show)

-- | The unknowns found, the iterations and the residual left (max |f_i|).
data Shot = Shot {shotUnknowns :: ![Double], shotIterations :: !Int, shotResidual :: !Double}
  deriving (Eq, Show)

-- | Doubles of work memory the method needs for n unknowns and m equations.
shootWorkSize :: Int -> Int -> Int
shootWorkSize n m = fromSize (pureCall (C.myb_shoot_work_size (toSize n) (toSize m)))

-- | Newton's method for f(u) = 0 with the unknowns u (their number n) and m equations: the
-- residual function gives the m values and, when asked (True), the m x n Jacobian by rows.
shoot :: ([Double] -> Bool -> IO (Either Status ([Double], [Double]))) -> Int -> [Double] -> ShootOptions -> IO (Result Shot)
shoot residual m u options = do
  function <- C.mk'myb_shoot_function act
  done <- withDoubles u $ \pu -> allocaArray (max 1 (shootWorkSize n m)) $ \work -> alloca $ \iterations -> alloca $ \left -> do
    poke iterations 0
    poke left 0
    code <- with made $ \po -> C.myb_shoot function nullPtr (toSize n) (toSize m) pu po work iterations left
    andThen code (Right <$> (Shot <$> peekDoubles n pu <*> (fromIntegral <$> peek iterations) <*> (coerce <$> peek left)))
  freeHaskellFunPtr function
  pure done
  where
    n = length u
    made = C.zero {C.myb_shoot_options'tolerance = coerce (shootTolerance options), C.myb_shoot_options'max_iterations = fromIntegral (shootMaxIterations options)}
    act _ pu f jacobian = do
      unknowns <- peekDoubles n pu
      given <- residual unknowns (jacobian /= nullPtr)
      case given of
        Left status -> pure (statusCode status)
        Right (values, derivatives)
          | length values /= m || (jacobian /= nullPtr && length derivatives /= m * n) -> pure C.MYB_E_DOMAIN
          | otherwise -> do
              pokeArray f (coerce values)
              if jacobian /= nullPtr then pokeArray jacobian (coerce derivatives) else pure ()
              pure C.MYB_OK

-- | Every solution with 0 .. maxRevolutions complete revolutions from r1 to r2 in the time of
-- flight; retrograde: the motion clockwise about +z.
lambert :: Double -> Vec3 -> Vec3 -> Double -> Bool -> Int -> Result [C.Myb_lambert_solution]
lambert mu r1 r2 tof retrograde revolutions
  | revolutions < 0 = failure Domain
  | otherwise = pureCall $ withVec r1 $ \p1 -> withVec r2 $ \p2 -> allocaArray capacity $ \solutions -> alloca $ \count -> do
      poke count 0
      code <- C.myb_lambert (coerce mu) p1 p2 (coerce tof) (if retrograde then 1 else 0 :: CInt) (fromIntegral revolutions) solutions (toSize capacity) count
      andThen code (peek count >>= \k -> Right <$> peekArray (fromSize k) solutions)
  where
    capacity = 1 + 2 * revolutions
