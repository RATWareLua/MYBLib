-- SPDX-License-Identifier: BSD-3-Clause
-- | Systems of point masses (myb/dyn.h): the Newtonian N-body problem with a graph of
-- interactions and the first post-Newtonian terms, its invariants, the fast multipole method,
-- and the Wisdom-Holman integrators. Bodies are given by GM; positions and velocities are one
-- 'Vec3' per body, the massive bodies first, then the test bodies.
module Myb.NBody
  ( -- * The system
    NBody
  , newNBody
  , nbodyCount
  , nbodyMassive
  , nbodyGm
  , Interaction (..)
  , setInteraction
  , setRelativity
  , accelerations
  , accelerationsXV
  , invariants
  , pnInvariants
  , C.Myb_nbody_invariants (..)
  , C.Myb_nbody_pn_invariants (..)
  , C.Myb_dd (..)
  , C.Myb_ppn (..)
    -- * Threads
  , Pool
  , newPool
  , poolThreads
  , setPool
    -- * The fast multipole method
  , setApproximation
  , approximation
    -- * Wisdom-Holman integrators
  , WhScheme (..)
  , WhOptions (..)
  , Wh
  , newWh
  , whStart
  , whStep
  , whIntegrate
  , whTime
  , whState
  , whDense
  , whStatistics
  , C.Myb_wh_stats (..)
    -- * The objects of C
  , nbodyPtr
  , poolPtr
  , whPtr
  , wrapNBody
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, nullPtr)
import Foreign.Storable (peek)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | A system of point masses and what it keeps (its pool of threads).
data NBody = NBody (ForeignPtr C.Myb_nbody) Kept

-- | massive >= 1 bodies with the GM given (m^3 s^-2, one each) and test bodies after them.
newNBody :: Int -> Int -> [Double] -> IO (Result NBody)
newNBody massive test gm
  | massive < 0 || test < 0 || length gm /= massive = pure (failure Domain)
  | otherwise =
      withDoubles gm $ \g ->
        fmap (uncurry NBody) <$> newOwned (\o -> check <$> C.myb_nbody_create (toSize massive) (toSize test) g nullPtr o) C.myb_nbody_destroy

-- | The point masses of a system of myb/system.h: the pointer the system owns, kept with the
-- system (for "Myb.System").
wrapNBody :: ForeignPtr C.Myb_nbody -> Kept -> NBody
wrapNBody = NBody

nbodyPtr :: NBody -> ForeignPtr C.Myb_nbody
nbodyPtr (NBody p _) = p

-- | Massive and test bodies.
nbodyCount :: NBody -> IO Int
nbodyCount (NBody p _) = fromSize <$> withForeignPtr p C.myb_nbody_count

nbodyMassive :: NBody -> IO Int
nbodyMassive (NBody p _) = fromSize <$> withForeignPtr p C.myb_nbody_massive

-- | The GM of a body (0 for a test body); 'Nothing' for no such body.
nbodyGm :: NBody -> Int -> IO (Maybe Double)
nbodyGm nb@(NBody p _) body = do
  count <- nbodyCount nb
  if body < 0 || body >= count then pure Nothing else Just . coerce <$> withForeignPtr p (\n -> C.myb_nbody_gm n (toSize body))

-- | How a source acts on a receiver.
data Interaction = NoInteraction | Newton | Eih
  deriving (Eq, Show)

-- | The level of the action of the source on the receiver.
setInteraction :: NBody -> Int -> Int -> Interaction -> IO (Result ())
setInteraction (NBody p _) receiver source level
  | receiver < 0 || source < 0 = pure (failure Domain)
  | otherwise = withForeignPtr p $ \n -> check <$> C.myb_nbody_set_interaction n (toSize receiver) (toSize source) code
  where
    code = case level of
      NoInteraction -> C.MYB_INTERACTION_NONE
      Newton -> C.MYB_INTERACTION_NEWTON
      Eih -> C.MYB_INTERACTION_EIH

-- | The 1PN terms with the parameters given (c, beta, gamma); 'Nothing': back to Newtonian.
setRelativity :: NBody -> Maybe C.Myb_ppn -> IO (Result ())
setRelativity (NBody p _) ppn = withForeignPtr p $ \n -> case ppn of
  Nothing -> check <$> C.myb_nbody_set_relativity n nullPtr
  Just value -> with value $ \v -> check <$> C.myb_nbody_set_relativity n v

-- | The accelerations at the positions (one per body), at the time t.
accelerations :: NBody -> Double -> [Vec3] -> IO (Result [Vec3])
accelerations nb@(NBody p _) t x = do
  count <- nbodyCount nb
  requireLength count x $ withForeignPtr p $ \n -> withVecs x $ \px ->
    fillArray count (\a -> C.myb_nbody_acceleration (castPtr n) (coerce t) px (castPtr a))

-- | The same with the velocities, for the 1PN terms.
accelerationsXV :: NBody -> Double -> [Vec3] -> [Vec3] -> IO (Result [Vec3])
accelerationsXV nb@(NBody p _) t x v = do
  count <- nbodyCount nb
  requireLength count x $ requireLength count v $ withForeignPtr p $ \n -> withVecs x $ \px -> withVecs v $ \pv ->
    fillArray count (\a -> C.myb_nbody_acceleration_xv (castPtr n) (coerce t) px pv (castPtr a))

-- | The conserved quantities of the massive bodies, in double-word arithmetic.
invariants :: NBody -> [Vec3] -> [Vec3] -> IO (Result C.Myb_nbody_invariants)
invariants nb@(NBody p _) x v = do
  count <- nbodyCount nb
  requireLength count x $ requireLength count v $ withForeignPtr p $ \n -> withVecs x $ \px -> withVecs v $ \pv ->
    fill (C.myb_nbody_invariants_of n px pv)

-- | The first post-Newtonian conserved quantities (the parameters of 'setRelativity').
pnInvariants :: NBody -> [Vec3] -> [Vec3] -> IO (Result C.Myb_nbody_pn_invariants)
pnInvariants nb@(NBody p _) x v = do
  count <- nbodyCount nb
  requireLength count x $ requireLength count v $ withForeignPtr p $ \n -> withVecs x $ \px -> withVecs v $ \pv ->
    fill (C.myb_nbody_pn_invariants_of n px pv)

-- | Threads for the accelerations (myb/pool.h); the results do not depend on them.
newtype Pool = Pool (ForeignPtr C.Myb_pool)

-- | A pool of the threads given (0: as many as the processors).
newPool :: Int -> IO (Result Pool)
newPool threads = fmap (Pool . fst) <$> newOwned (\o -> check <$> C.myb_pool_create (fromIntegral threads) nullPtr o) C.myb_pool_destroy

poolThreads :: Pool -> IO Int
poolThreads (Pool p) = fromIntegral <$> withForeignPtr p C.myb_pool_threads

poolPtr :: Pool -> ForeignPtr C.Myb_pool
poolPtr (Pool p) = p

-- | The threads of the system ('Nothing': the calling thread only). The system keeps the pool.
setPool :: NBody -> Maybe Pool -> IO ()
setPool (NBody p kept) pool = withForeignPtr p $ \n -> case pool of
  Nothing -> C.myb_nbody_set_pool n nullPtr
  Just (Pool q) -> keepAlive kept q >> withForeignPtr q (C.myb_nbody_set_pool n)

-- | The fast multipole method of the order given (0: exact sums again) and opening angle.
setApproximation :: NBody -> Int -> Double -> IO (Result ())
setApproximation (NBody p _) order theta = withForeignPtr p $ \n -> check <$> C.myb_nbody_set_approximation n (fromIntegral order) (coerce theta)

-- | The order and the opening angle of the approximation; (0, 0) when the sums are exact.
approximation :: NBody -> IO (Int, Double)
approximation (NBody p _) = withForeignPtr p $ \n -> alloca $ \theta -> do
  order <- C.myb_nbody_approximation n theta
  (,) (fromIntegral order) . coerce <$> peek theta

-- | The schemes of the splitting.
data WhScheme = Dkd | Saba2 | Saba3 | Saba4 | Sbab1 | Sbab2 | Sbab3 | Sbab4
  deriving (Eq, Show)

-- | The step (s; its sign the direction), the scheme, the order of the corrector (0, or odd
-- 3 .. 17 with 'Dkd').
data WhOptions = WhOptions {whStepSize :: !Double, whScheme :: !WhScheme, whCorrector :: !Int}
  deriving (Eq, Show)

-- | An integrator, what it keeps and its system.
data Wh = Wh (ForeignPtr C.Myb_wh) Kept NBody

newWh :: NBody -> WhOptions -> IO (Result Wh)
newWh nb@(NBody p _) options =
  with made $ \o -> withForeignPtr p $ \n -> do
    created <- newOwned (\w -> check <$> C.myb_wh_create n o nullPtr w) C.myb_wh_destroy
    case created of
      Left e -> pure (Left e)
      Right (fp, kept) -> keepAlive kept p >> pure (Right (Wh fp kept nb))
  where
    made =
      C.zero
        { C.myb_wh_options'h = coerce (whStepSize options)
        , C.myb_wh_options'scheme = case whScheme options of
            Dkd -> C.MYB_WH_DKD
            Saba2 -> C.MYB_WH_SABA2
            Saba3 -> C.MYB_WH_SABA3
            Saba4 -> C.MYB_WH_SABA4
            Sbab1 -> C.MYB_WH_SBAB1
            Sbab2 -> C.MYB_WH_SBAB2
            Sbab3 -> C.MYB_WH_SBAB3
            Sbab4 -> C.MYB_WH_SBAB4
        , C.myb_wh_options'corrector = fromIntegral (whCorrector options)
        }

whPtr :: Wh -> ForeignPtr C.Myb_wh
whPtr (Wh p _ _) = p

-- | Starts at t0 from the inertial positions and velocities.
whStart :: Wh -> Double -> [Vec3] -> [Vec3] -> IO (Result ())
whStart (Wh p _ nb) t0 x v = do
  count <- nbodyCount nb
  requireLength count x $ requireLength count v $ withForeignPtr p $ \w -> withVecs x $ \px -> withVecs v $ \pv ->
    check <$> C.myb_wh_start w (coerce t0) px pv

whStep :: Wh -> IO (Result ())
whStep (Wh p _ _) = withForeignPtr p $ \w -> check <$> C.myb_wh_step w

-- | Steps to the node nearest to tEnd.
whIntegrate :: Wh -> Double -> IO (Result ())
whIntegrate (Wh p _ _) tEnd = withForeignPtr p $ \w -> check <$> C.myb_wh_integrate w (coerce tEnd)

whTime :: Wh -> IO Double
whTime (Wh p _ _) = coerce <$> withForeignPtr p C.myb_wh_time

-- | The corrected state at the current node: positions and velocities.
whState :: Wh -> IO (Result ([Vec3], [Vec3]))
whState (Wh p _ nb) = do
  count <- nbodyCount nb
  withForeignPtr p $ \w -> allocaArray (max 1 count) $ \px -> allocaArray (max 1 count) $ \pv -> do
    code <- C.myb_wh_state w (castPtr px) (castPtr pv)
    andThen code (Right <$> ((,) <$> peekArray count px <*> peekArray count pv))

-- | The state at t within a step of the node.
whDense :: Wh -> Double -> IO (Result ([Vec3], [Vec3]))
whDense (Wh p _ nb) t = do
  count <- nbodyCount nb
  withForeignPtr p $ \w -> allocaArray (max 1 count) $ \px -> allocaArray (max 1 count) $ \pv -> do
    code <- C.myb_wh_dense w (coerce t) (castPtr px) (castPtr pv)
    andThen code (Right <$> ((,) <$> peekArray count px <*> peekArray count pv))

whStatistics :: Wh -> IO C.Myb_wh_stats
whStatistics (Wh p _ _) = withForeignPtr p C.myb_wh_statistics
