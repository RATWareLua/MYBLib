-- SPDX-License-Identifier: BSD-3-Clause
-- | The forces on one object (myb/model.h) with the ready-made terms of myb/env.h, and terms
-- written in Haskell. The model owns the structures of its terms and keeps alive what they
-- point to: environments, fields, Haskell functions.
module Myb.Force
  ( ForceModel
  , newForceModel
  , Object (..)
  , ExponentialAtmosphere (..)
  , addPointMass
  , addThirdBody
  , addField
  , addRadiation
  , addDrag
  , addSchwarzschild
  , addSchwarzschildPpn
  , addTerm
  , addRaw
  , acceleration
  , termCount
  , termName
  , termEnabled
  , setEnabled
  , BreakScheme (..)
  , setBreakScheme
  , breakScheme
  , forceModelKeep
  , forceModelPtr
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullFunPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Environment (Environment, bodyCount, environmentPtr)
import Myb.Error
import Myb.Gravity (GravityField, fieldDegree, gravityFieldPtr)
import Myb.Internal
import Myb.Types

-- | A force model: its object and what it keeps.
data ForceModel = ForceModel (ForeignPtr C.Myb_force_model) Kept

newForceModel :: IO (Result ForceModel)
newForceModel = do
  kept <- newKept
  create
    (\out -> check <$> C.myb_force_model_create nullPtr out)
    (\p -> (`ForceModel` kept) <$> owned p C.myb_force_model_destroy kept)

-- | The object the forces act on.
data Object = Object
  { -- | kg
    objectMass :: !Double
  , -- | m^2: the cross-section
    objectArea :: !Double
  , -- | the radiation pressure coefficient, 1 + reflectivity
    objectCr :: !Double
  , -- | the drag coefficient
    objectCd :: !Double
  }
  deriving (Eq, Show)

objectToC :: Object -> C.Myb_object
objectToC o =
  C.zero
    { C.myb_object'mass = coerce (objectMass o)
    , C.myb_object'area = coerce (objectArea o)
    , C.myb_object'cr = coerce (objectCr o)
    , C.myb_object'cd = coerce (objectCd o)
    }

-- | A layered exponential atmosphere over an ellipsoid: (base height m, density kg/m^3 at the
-- base, scale height m) of each layer, by height.
data ExponentialAtmosphere = ExponentialAtmosphere
  { equatorialRadius :: !Double
  , flattening :: !Double
  , layers :: ![(Double, Double, Double)]
  }
  deriving (Eq, Show)

-- | Any term of the C interface. Returns its index. The function and the context must stay
-- valid as long as the model: keep what the context points to with 'forceModelKeep'.
addRaw :: ForceModel -> String -> C.Myb_term_function -> C.Myb_term_jacobian -> Ptr () -> IO (Result Int)
addRaw (ForceModel model _) name function jacobian context =
  withText name $ \cname -> withForeignPtr model $ \m -> alloca $ \index -> do
    code <- C.myb_force_model_add m cname function jacobian context index
    either (pure . Left) (const (Right . fromSize <$> peek index)) (check code)

-- | A ready-made term whose structure the model keeps.
addKept :: ForceModel -> String -> C.Myb_term_function -> C.Myb_term_jacobian -> Ptr a -> IO (Result Int)
addKept model name function jacobian term = addRaw model name function jacobian (castPtr term)

-- | The index of a body of the environment, or a failure of the domain if it has none such.
bodyOf :: Environment -> Int -> IO (Result Int)
bodyOf environment body = do
  count <- bodyCount environment
  pure (if body >= 0 && body < count then Right body else failure Domain)

-- | The attraction of a point mass at the origin (gm in m^3/s^2).
addPointMass :: ForceModel -> String -> Double -> IO (Result Int)
addPointMass model@(ForceModel _ kept) name gm = do
  term <- keepValue kept (C.zero {C.myb_point_mass_term'gm = coerce gm})
  addKept model name C.p'myb_term_point_mass C.p'myb_term_point_mass_jacobian term

-- | The attraction of a body of the environment away from the origin, with its indirect part.
addThirdBody :: ForceModel -> String -> Environment -> Int -> IO (Result Int)
addThirdBody model@(ForceModel _ kept) name environment body = do
  index <- bodyOf environment body
  case index of
    Left e -> pure (Left e)
    Right b -> do
      let env = environmentPtr environment
      term <- keepValue kept (C.zero {C.myb_third_body_term'env = rawPtr env, C.myb_third_body_term'body = toSize b})
      keepAlive kept env
      addKept model name C.p'myb_term_third_body C.p'myb_term_third_body_jacobian term

-- | The field of a body of the environment in harmonics of the degrees nMin .. nMax
-- ('Nothing': all the field has).
addField :: ForceModel -> String -> Environment -> Int -> GravityField -> Int -> Maybe Int -> IO (Result Int)
addField model@(ForceModel _ kept) name environment body field nMin nMax = do
  index <- bodyOf environment body
  degree <- fieldDegree field
  case index of
    Left e -> pure (Left e)
    Right b -> do
      made <- create (\out -> check <$> C.myb_gravity_work_create (fromIntegral degree) nullPtr out) pure
      case made of
        Left e -> pure (Left e)
        Right work -> do
          keepRelease kept (C.myb_gravity_work_destroy work)
          let env = environmentPtr environment
              harmonics = gravityFieldPtr field
          term <-
            keepValue kept $
              C.zero
                { C.myb_field_term'env = rawPtr env
                , C.myb_field_term'body = toSize b
                , C.myb_field_term'field = rawPtr harmonics
                , C.myb_field_term'work = work
                , C.myb_field_term'n_min = fromIntegral (max 0 nMin)
                , C.myb_field_term'n_max = maybe (-1) (fromIntegral . max 0) nMax
                }
          keepAlive kept env
          keepAlive kept harmonics
          addKept model name C.p'myb_term_field C.p'myb_term_field_jacobian term

-- | The pressure of the light of a body of the environment on the object, in the shadows of
-- the occulting bodies; c: the speed of light (m/s).
addRadiation :: ForceModel -> String -> Environment -> Int -> Object -> Double -> [Int] -> IO (Result Int)
addRadiation model@(ForceModel _ kept) name environment source object c occulters = do
  count <- bodyCount environment
  let valid b = b >= 0 && b < count
  if not (valid source) || not (all valid occulters) || length occulters > C.MYB_ENV_MAX_OCCULTERS
    then pure (failure Domain)
    else do
      let env = environmentPtr environment
      o <- keepValue kept (objectToC object)
      term <-
        keepValue kept $
          C.zero
            { C.myb_radiation_term'env = rawPtr env
            , C.myb_radiation_term'source = toSize source
            , C.myb_radiation_term'object = o
            , C.myb_radiation_term'c = coerce c
            , C.myb_radiation_term'occulters = map toSize occulters
            , C.myb_radiation_term'occulter_count = toSize (length occulters)
            }
      keepAlive kept env
      addKept model name C.p'myb_term_radiation C.p'myb_term_radiation_jacobian term

-- | The drag of an exponential atmosphere that turns with a body of the environment.
addDrag :: ForceModel -> String -> Environment -> Int -> ExponentialAtmosphere -> Object -> IO (Result Int)
addDrag model@(ForceModel _ kept) name environment body air object = do
  index <- bodyOf environment body
  let count = length (layers air)
  case index of
    Left e -> pure (Left e)
    Right _ | count == 0 || count > C.MYB_ATMOSPHERE_LAYERS -> pure (failure Domain)
    Right b -> do
      let env = environmentPtr environment
          (bases, densities, scales) = unzip3 (layers air)
      atmosphere <-
        keepValue kept $
          C.zero
            { C.myb_exponential_atmosphere'ellipsoid =
                C.zero {C.myb_ellipsoid'a = coerce (equatorialRadius air), C.myb_ellipsoid'f = coerce (flattening air)}
            , C.myb_exponential_atmosphere'count = toSize count
            , C.myb_exponential_atmosphere'base = coerce bases
            , C.myb_exponential_atmosphere'density = coerce densities
            , C.myb_exponential_atmosphere'scale = coerce scales
            }
      o <- keepValue kept (objectToC object)
      term <-
        keepValue kept $
          C.zero
            { C.myb_drag_term'env = rawPtr env
            , C.myb_drag_term'body = toSize b
            , C.myb_drag_term'density = C.p'myb_density_exponential
            , C.myb_drag_term'density_context = castPtr atmosphere
            , C.myb_drag_term'object = o
            }
      keepAlive kept env
      addKept model name C.p'myb_term_drag C.p'myb_term_drag_jacobian term

-- | The post-Newtonian correction of a central body of the environment (Schwarzschild), in
-- general relativity (beta = gamma = 1); c: the speed of light (m/s).
addSchwarzschild :: ForceModel -> String -> Environment -> Int -> Double -> IO (Result Int)
addSchwarzschild model name environment central c = addSchwarzschildPpn model name environment central c 1 1

-- | The same with the parameters beta and gamma of the parametrized post-Newtonian theory.
addSchwarzschildPpn :: ForceModel -> String -> Environment -> Int -> Double -> Double -> Double -> IO (Result Int)
addSchwarzschildPpn model@(ForceModel _ kept) name environment central c beta gamma = do
  index <- bodyOf environment central
  case index of
    Left e -> pure (Left e)
    Right b -> do
      let env = environmentPtr environment
      term <-
        keepValue kept $
          C.zero
            { C.myb_relativity_term'env = rawPtr env
            , C.myb_relativity_term'ppn = C.zero {C.myb_ppn'c = coerce c, C.myb_ppn'beta = coerce beta, C.myb_ppn'gamma = coerce gamma}
            , C.myb_relativity_term'central = toSize b
            , C.myb_relativity_term'external = toSize b
            }
      keepAlive kept env
      addKept model name C.p'myb_term_schwarzschild C.p'myb_term_schwarzschild_jacobian term

-- | A term of the caller's own: the acceleration at the model time t and the state (x, v),
-- or the status that stops the evaluation. Its derivatives are taken by differences of the
-- model. The model keeps the function.
addTerm :: ForceModel -> String -> (Double -> Vec3 -> Vec3 -> IO (Either Status Vec3)) -> IO (Result Int)
addTerm model@(ForceModel _ kept) name term = do
  function <- C.mk'myb_term_function call
  keepFunPtr kept function
  addRaw model name function nullFunPtr nullPtr
  where
    call _ t x v a = do
      position' <- peek (castPtr x)
      velocity' <- peek (castPtr v)
      given <- term (coerce t) position' velocity'
      case given of
        Right acc -> poke (castPtr a) acc >> pure C.MYB_OK
        Left status -> pure (statusCode status)

-- | The acceleration of the object at the model time t and the state (x, v).
acceleration :: ForceModel -> Double -> Vec3 -> Vec3 -> IO (Result Vec3)
acceleration (ForceModel model _) t x v =
  withForeignPtr model $ \m -> with x $ \px -> with v $ \pv -> alloca $ \out -> do
    code <- C.myb_force_model_acceleration m (coerce t) (castPtr px) (castPtr pv) (castPtr out)
    either (pure . Left) (const (Right <$> peek out)) (check code)

termCount :: ForceModel -> IO Int
termCount (ForceModel model _) = fromSize <$> withForeignPtr model C.myb_force_model_count

-- | The name of a term; 'Nothing' if there is no such term.
termName :: ForceModel -> Int -> IO (Maybe String)
termName fm@(ForceModel model _) index = do
  count <- termCount fm
  if index < 0 || index >= count
    then pure Nothing
    else withForeignPtr model $ \m -> Just <$> (C.myb_force_model_name m (toSize index) >>= peekText)

-- | Whether a term is switched on; 'Nothing' if there is no such term.
termEnabled :: ForceModel -> Int -> IO (Maybe Bool)
termEnabled fm@(ForceModel model _) index = do
  count <- termCount fm
  if index < 0 || index >= count
    then pure Nothing
    else Just . (/= 0) <$> withForeignPtr model (\m -> C.myb_force_model_enabled m (toSize index))

-- | Switches a term on or off; False if there is no such term.
setEnabled :: ForceModel -> Int -> Bool -> IO Bool
setEnabled fm@(ForceModel model _) index on = do
  count <- termCount fm
  if index < 0 || index >= count
    then pure False
    else True <$ withForeignPtr model (\m -> C.myb_force_model_set_enabled m (toSize index) (if on then 1 else 0))

-- | How an orbit integrated by IAS15 ends its steps at the breaks of the model.
data BreakScheme
  = -- | a step found to hold a break is taken again up to it (the default)
    BreaksRetake
  | -- | the steps are aimed at the breaks predicted before them
    BreaksPredict
  deriving (Eq, Show)

-- | The scheme of the breaks, from the next start of an orbit.
setBreakScheme :: ForceModel -> BreakScheme -> IO ()
setBreakScheme (ForceModel model _) scheme = withForeignPtr model $ \m ->
  C.myb_force_model_set_break_scheme m (case scheme of BreaksRetake -> C.MYB_BREAKS_RETAKE; BreaksPredict -> C.MYB_BREAKS_PREDICT)

breakScheme :: ForceModel -> IO BreakScheme
breakScheme (ForceModel model _) = do
  scheme <- withForeignPtr model C.myb_force_model_break_scheme
  pure (if scheme == C.MYB_BREAKS_PREDICT then BreaksPredict else BreaksRetake)

-- | Ties the life of an object to the model: for the context of a term added with 'addRaw'.
forceModelKeep :: ForceModel -> ForeignPtr a -> IO ()
forceModelKeep (ForceModel _ kept) = keepAlive kept

forceModelPtr :: ForceModel -> ForeignPtr C.Myb_force_model
forceModelPtr (ForceModel model _) = model
