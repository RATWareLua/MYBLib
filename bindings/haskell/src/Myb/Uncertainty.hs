-- SPDX-License-Identifier: BSD-3-Clause
-- | Uncertainty (myb/uncert.h): the variational equations of a force model, transformations of
-- covariances, and deterministic random numbers (the same sequence on every platform).
module Myb.Uncertainty
  ( -- * Variational equations
    Variational
  , C.Myb_variational_options (..)
  , newVariational
  , variationalStart
  , variationalStep
  , variationalIntegrate
  , variationalTime
  , variationalGet
  , variationalDense
  , variationalPtr
    -- * Covariances
  , covarianceTransform
  , covarianceRtn
    -- * Random numbers
  , Rng
  , newRng
  , rngNext
  , rngUniform
  , rngNormal
  , rngGaussian
  , rngState
  , sampleStatistics
  ) where

import Data.Coerce (coerce)
import Data.Word (Word64)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (ForeignPtr, mallocForeignPtr, withForeignPtr)
import Foreign.Marshal.Array (allocaArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Force (ForceModel, forceModelPtr)
import Myb.Internal
import Myb.Types

-- | The variational equations of a model: the object, what it keeps (the model) and the number
-- of the parameters differentiated.
data Variational = Variational (ForeignPtr C.Myb_variational) Kept Int

variationalPtr :: Variational -> ForeignPtr C.Myb_variational
variationalPtr (Variational p _ _) = p

newVariational :: ForceModel -> C.Myb_variational_options -> IO (Result Variational)
newVariational model options = withForeignPtr (forceModelPtr model) $ \m -> do
  np <- fromSize <$> C.myb_force_model_parameter_count m
  made <- with options $ \o -> newOwned (\v -> check <$> C.myb_variational_create m o nullPtr v) C.myb_variational_destroy
  case made of
    Left e -> pure (Left e)
    Right (fp, kept) -> keepAlive kept (forceModelPtr model) >> pure (Right (Variational fp kept np))

-- | Starts at t0 from the state with Phi = I and S = 0.
variationalStart :: Variational -> Double -> State -> IO (Result ())
variationalStart (Variational p _ _) t0 y0 = withForeignPtr p $ \v -> withVec y0 $ \py -> check <$> C.myb_variational_start v (coerce t0) py

-- | One step towards tLimit.
variationalStep :: Variational -> Double -> IO (Result ())
variationalStep (Variational p _ _) t = withForeignPtr p $ \v -> check <$> C.myb_variational_step v (coerce t)

variationalIntegrate :: Variational -> Double -> IO (Result ())
variationalIntegrate (Variational p _ _) t = withForeignPtr p $ \v -> check <$> C.myb_variational_integrate v (coerce t)

variationalTime :: Variational -> IO Double
variationalTime (Variational p _ _) = coerce <$> withForeignPtr p C.myb_variational_time

-- | The state, Phi (36 by rows) and S (6 x parameters by rows) now.
variationalGet :: Variational -> IO (State, [Double], [Double])
variationalGet (Variational p _ np) = withForeignPtr p $ \v -> with (State (Vec3 0 0 0) (Vec3 0 0 0)) $ \py -> allocaArray 36 $ \phi -> allocaArray (max 1 (6 * np)) $ \s -> do
  C.myb_variational_get v (castPtr py) phi s
  (,,) <$> peek py <*> peekDoubles 36 phi <*> peekDoubles (6 * np) s

-- | The same at t within the last step.
variationalDense :: Variational -> Double -> IO (Result (State, [Double], [Double]))
variationalDense (Variational p _ np) t = withForeignPtr p $ \v -> with (State (Vec3 0 0 0) (Vec3 0 0 0)) $ \py -> allocaArray 36 $ \phi -> allocaArray (max 1 (6 * np)) $ \s -> do
  code <- C.myb_variational_dense v (coerce t) (castPtr py) phi s
  andThen code (Right <$> ((,,) <$> peek py <*> peekDoubles 36 phi <*> peekDoubles (6 * np) s))

-- | M P M^T for n x n matrices (by rows).
covarianceTransform :: Int -> [Double] -> [Double] -> Result [Double]
covarianceTransform n m p
  | n < 1 || length m /= n * n || length p /= n * n = failure Domain
  | otherwise = pureCall $ withDoubles m $ \pm -> withDoubles p $ \pp -> allocaArray (n * n) $ \o -> do
      C.myb_covariance_transform (toSize n) pm pp o
      Right <$> peekDoubles (n * n) o

-- | A state covariance (36) in the radial, transverse, normal axes of (r, v).
covarianceRtn :: Vec3 -> Vec3 -> [Double] -> Result [Double]
covarianceRtn r v p
  | length p /= 36 = failure Domain
  | otherwise = pureCall $ withVec r $ \pr -> withVec v $ \pv -> withDoubles p $ \pp -> fillDoubles 36 (C.myb_covariance_rtn pr pv pp)

-- | A generator of xoshiro256** (mutable).
newtype Rng = Rng (ForeignPtr C.Myb_rng)

-- | A generator seeded by SplitMix64.
newRng :: Word64 -> IO Rng
newRng seed = do
  memory <- mallocForeignPtr
  withForeignPtr memory $ \p -> poke p C.zero >> C.myb_rng_seed p seed
  pure (Rng memory)

rngNext :: Rng -> IO Word64
rngNext (Rng p) = withForeignPtr p C.myb_rng_next

-- | In [0, 1), 53 random bits.
rngUniform :: Rng -> IO Double
rngUniform (Rng p) = coerce <$> withForeignPtr p C.myb_rng_uniform

-- | A standard normal deviate.
rngNormal :: Rng -> IO Double
rngNormal (Rng p) = coerce <$> withForeignPtr p C.myb_rng_normal

-- | A sample of N(mean, L L^T), L lower triangular (n x n by rows).
rngGaussian :: Rng -> [Double] -> [Double] -> IO (Result [Double])
rngGaussian (Rng p) mean l
  | length l /= n * n || n < 1 = pure (failure Domain)
  | otherwise = withForeignPtr p $ \r -> withDoubles mean $ \pm -> withDoubles l $ \pl -> allocaArray n $ \o -> do
      C.myb_rng_gaussian r (toSize n) pm pl o
      Right <$> peekDoubles n o
  where
    n = length mean

-- | The structure of C of the generator (its state).
rngState :: Rng -> IO C.Myb_rng
rngState (Rng p) = withForeignPtr p peek

-- | The mean (n) and the covariance (n x n, divisor count - 1) of samples of dimension n
-- stored one after the other.
sampleStatistics :: Int -> [Double] -> Result ([Double], [Double])
sampleStatistics n samples
  | n < 1 || length samples `mod` n /= 0 = failure Domain
  | otherwise = pureCall $ withDoubles samples $ \ps -> allocaArray n $ \pm -> allocaArray (n * n) $ \pc -> do
      code <- C.myb_sample_statistics (toSize (length samples `div` n)) (toSize n) ps pm pc
      andThen code (Right <$> ((,) <$> peekDoubles n pm <*> peekDoubles (n * n) pc))
