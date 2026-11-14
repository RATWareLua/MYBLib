-- SPDX-License-Identifier: BSD-3-Clause
-- | Estimation of an orbit from observations (myb/estim.h): stations, measurements with their
-- partial derivatives, the weighted least squares over an arc and the sequential filter with
-- its smoother. The unknowns are the state at the epoch, the parameters of the force model
-- (myb_force_model_add_parameter) and the systematic shifts of the measurements, in that order.
module Myb.Estimation
  ( -- * Stations
    Station
  , newStation
  , stationAt
  , setStationMotion
  , stationState
  , stationRecord
  , C.Myb_station_motion (..)
    -- * Measurements
  , Measure (..)
  , Light (..)
  , Offset (..)
  , AttitudeLaw (..)
  , evaluateMeasure
    -- * Observations
  , Observation (..)
  , observation
    -- * Least squares
  , FitOptions (..)
  , defaultFitOptions
  , Fit (..)
  , C.Myb_fit_result (..)
  , fitOrbit
    -- * The filter
  , FilterOptions (..)
  , defaultFilterOptions
  , Filter
  , newFilter
  , filterStart
  , filterUpdate
  , filterPredict
  , filterGet
  , filterSetBias
  , filterBias
  , filterSteps
  , filterSmooth
  , filterStep
  , filterPtr
  ) where

import Data.Coerce (coerce)
import Data.IORef (IORef, modifyIORef', newIORef, readIORef)
import Foreign.C.Types (CDouble (..), CUChar)
import Foreign.ForeignPtr (ForeignPtr, castForeignPtr, mallocForeignPtr, touchForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, pokeArray, withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (Storable, peek, poke)
import qualified Myb.C as C
import Myb.Environment (Environment, environmentPtr)
import Myb.Error
import Myb.Force (ForceModel, forceModelPtr)
import Myb.Internal
import Myb.Types

-- ---- Stations ------------------------------------------------------------------------------

-- | An observer on a body of an environment: its structure of C, the environment and what its
-- motion points to.
data Station = Station (ForeignPtr C.Myb_station) Environment (IORef [ForeignPtr ()])

-- | The station for a call: the environment and what the motion needs stay alive meanwhile.
withStation :: Station -> (Ptr C.Myb_station -> IO a) -> IO a
withStation (Station memory env extras) k =
  withForeignPtr (environmentPtr env) $ \_ -> withForeignPtr memory $ \p -> do
    done <- k p
    readIORef extras >>= mapM_ touchForeignPtr
    pure done

-- | A station at geodetic longitude, latitude (rad) and height (m) on an ellipsoid of the
-- body, fixed on it.
newStation :: Environment -> Int -> C.Myb_ellipsoid -> Double -> Double -> Double -> IO (Result Station)
newStation env body ellipsoid longitude latitude height = do
  memory <- mallocForeignPtr
  extras <- newIORef []
  withForeignPtr (environmentPtr env) $ \e -> withForeignPtr memory $ \p -> with ellipsoid $ \pe -> do
    poke p C.zero
    code <- C.myb_station_init p e (toSize body) pe (coerce longitude) (coerce latitude) (coerce height)
    andThen code (pure (Right (Station memory env extras)))

-- | A station at a body-fixed position with its axes east, north, up (rows, body frame).
stationAt :: Environment -> Int -> Vec3 -> Mat3 -> IO Station
stationAt env body place enu = do
  memory <- mallocForeignPtr
  extras <- newIORef []
  withForeignPtr memory $ \p ->
    poke p $
      C.zero
        { C.myb_station'env = rawPtr (environmentPtr env)
        , C.myb_station'body = toSize body
        , C.myb_station'position = coerce (vecList place)
        , C.myb_station'enu = coerce (mat3List enu)
        }
  pure (Station memory env extras)
  where
    vecList (Vec3 a b c) = [a, b, c]

-- | What moves the station on its body (its position is then the place without its tides).
-- The tides and the loads the motion points to are objects of "Myb.C": pass their owners,
-- which the station keeps alive.
setStationMotion :: Station -> C.Myb_station_motion -> [ForeignPtr ()] -> IO ()
setStationMotion (Station memory _ extras) motion owners = do
  kept <- mallocForeignPtr
  withForeignPtr kept (`poke` motion)
  modifyIORef' extras ((castForeignPtr kept : owners) ++)
  withForeignPtr memory $ \p -> do
    s <- peek p
    poke p s {C.myb_station'motion = rawPtr kept}

-- | The state of the station at the model time t relative to the origin of the environment.
stationState :: Station -> Double -> IO (Result State)
stationState station t = withStation station $ \p -> fillVec (C.myb_station_state p (coerce t))

-- | The structure of C of the station.
stationRecord :: Station -> IO C.Myb_station
stationRecord (Station memory _ _) = withForeignPtr memory peek

-- ---- Measurements --------------------------------------------------------------------------

-- | A measurement: a function of the state of the object at the time of the observation, with
-- its derivatives.
data Measure
  = -- | the position of the object (3 values)
    MeasurePosition
  | -- | its state (6 values)
    MeasureState
  | -- | |d| from a station (1 value)
    MeasureRange Station
  | -- | d |d| / dt (1 value)
    MeasureRangeRate Station
  | -- | right ascension and declination (2 values; periodic = 1)
    MeasureRaDec Station
  | -- | azimuth and elevation (2 values; periodic = 1)
    MeasureAzEl Station
  | -- | with the light time
    MeasureLightRange Light
  | MeasureLightRangeRate Light
  | MeasureLightRaDec Light
  | MeasureLightAzEl Light
  | -- | a measurement of a point of the object off its centre of mass
    MeasureOffset Offset
  | -- | a Haskell function of the time and the state: the values and their partials by the
    -- state (values x 6, by rows)
    MeasureFunction (Double -> State -> IO (Either Status ([Double], [Double])))
  | -- | a measurement of C and its context, which the caller keeps alive
    MeasureC C.Myb_measure_function (Ptr ())

-- | The station of a measurement with the light time, the speed of light (m/s), the GM of the
-- body at the origin (0: a straight line over the light time), two-way, the aberration.
data Light = Light
  { lightStation :: Station
  , lightSpeed :: Double
  , lightGm :: Double
  , lightTwoWay :: Bool
  , lightAberration :: Bool
  }

-- | The measurement taken at a point offset (m, body frame) from the centre of mass, turned by
-- the attitude law; gm of the body at the origin (0: the velocity kept over the second).
data Offset = Offset {offsetMeasure :: Measure, offsetAttitude :: AttitudeLaw, offsetVector :: Vec3, offsetGm :: Double}

-- | The attitude of the object: the rotation inertial -> body at the time and the state.
data AttitudeLaw
  = AttitudeFunction (Double -> Vec3 -> Vec3 -> IO (Either Status Mat3))
  | -- | a law of myb/surface.h and its context, which the caller keeps alive
    AttitudeC C.Myb_attitude_function (Ptr ())

withAttitude :: AttitudeLaw -> (C.Myb_attitude_function -> Ptr () -> IO r) -> IO r
withAttitude law k = case law of
  AttitudeC f c -> k f c
  AttitudeFunction f -> do
    function <- C.mk'myb_attitude_function $ \_ t x v m -> do
      position' <- peek (castPtr x)
      velocity' <- peek (castPtr v)
      given <- f (coerce t) position' velocity'
      case given of
        Right matrix -> poke (castPtr m) matrix >> pure C.MYB_OK
        Left status -> pure (statusCode status)
    done <- k function nullPtr
    freeHaskellFunPtr function
    pure done

withLight :: Light -> (Ptr C.Myb_light_measure -> IO r) -> IO r
withLight light k = withStation (lightStation light) $ \s ->
  with
    ( C.zero
        { C.myb_light_measure'station = s
        , C.myb_light_measure'c = coerce (lightSpeed light)
        , C.myb_light_measure'gm = coerce (lightGm light)
        , C.myb_light_measure'two_way = if lightTwoWay light then 1 else 0
        , C.myb_light_measure'aberration = if lightAberration light then 1 else 0
        }
    )
    k

-- | The function of C and its context for a measurement of `count` values, for one call.
withMeasure :: Int -> Measure -> (C.Myb_measure_function -> Ptr () -> IO r) -> IO r
withMeasure count measure k = case measure of
  MeasurePosition -> k C.p'myb_measure_position nullPtr
  MeasureState -> k C.p'myb_measure_state nullPtr
  MeasureRange s -> withStation s (k C.p'myb_measure_range . castPtr)
  MeasureRangeRate s -> withStation s (k C.p'myb_measure_range_rate . castPtr)
  MeasureRaDec s -> withStation s (k C.p'myb_measure_radec . castPtr)
  MeasureAzEl s -> withStation s (k C.p'myb_measure_azel . castPtr)
  MeasureLightRange l -> withLight l (k C.p'myb_measure_light_range . castPtr)
  MeasureLightRangeRate l -> withLight l (k C.p'myb_measure_light_range_rate . castPtr)
  MeasureLightRaDec l -> withLight l (k C.p'myb_measure_light_radec . castPtr)
  MeasureLightAzEl l -> withLight l (k C.p'myb_measure_light_azel . castPtr)
  MeasureOffset o ->
    withMeasure count (offsetMeasure o) $ \f c -> withAttitude (offsetAttitude o) $ \af ac ->
      with
        ( C.zero
            { C.myb_offset_measure'measure = f
            , C.myb_offset_measure'context = c
            , C.myb_offset_measure'attitude = af
            , C.myb_offset_measure'attitude_context = ac
            , C.myb_offset_measure'offset = let Vec3 a b d = offsetVector o in coerce [a, b, d]
            , C.myb_offset_measure'gm = coerce (offsetGm o)
            }
        )
        (k C.p'myb_measure_offset . castPtr)
  MeasureC f c -> k f c
  MeasureFunction f -> do
    function <- C.mk'myb_measure_function $ \_ t y values partials -> do
      state <- peek (castPtr y)
      given <- f (coerce t) state
      case given of
        Left status -> pure (statusCode status)
        Right (vs, ps)
          | length vs /= count || (partials /= nullPtr && length ps /= 6 * count) -> pure C.MYB_E_DOMAIN
          | otherwise -> do
              pokeArray values (coerce vs)
              if partials /= nullPtr then pokeArray partials (coerce ps) else pure ()
              pure C.MYB_OK
    done <- k function nullPtr
    freeHaskellFunPtr function
    pure done

-- | The `count` values of a measurement at the time and the state, and their partials by the
-- state (count x 6, by rows).
evaluateMeasure :: Measure -> Int -> Double -> State -> IO (Result ([Double], [Double]))
evaluateMeasure measure count t y
  | count < 1 || count > C.MYB_MEASURE_MAX = pure (failure Domain)
  | otherwise =
      withMeasure count measure $ \f c -> withVec y $ \py -> allocaArray count $ \values -> allocaArray (6 * count) $ \partials -> do
        code <- C.call'myb_measure_function f c (coerce t) py values partials
        andThen code (Right <$> ((,) <$> peekDoubles count values <*> peekDoubles (6 * count) partials))

-- ---- Observations --------------------------------------------------------------------------

-- | An observation at a model time: the measurement, the values observed and their standard
-- deviations; periodic: bit i for a value that is an angle; bias: 0, or k > 0 for the shifts
-- number k - 1 + i of the estimator added to what is computed for value i.
data Observation = Observation
  { obsTime :: Double
  , obsMeasure :: Measure
  , obsValue :: [Double]
  , obsSigma :: [Double]
  , obsPeriodic :: Word
  , obsBias :: Int
  }

-- | An observation without angles and without shifts.
observation :: Double -> Measure -> [Double] -> [Double] -> Observation
observation t m values sigmas = Observation t m values sigmas 0 0

valid :: Observation -> Bool
valid o = let n = length (obsValue o) in n >= 1 && n <= C.MYB_MEASURE_MAX && length (obsSigma o) == n

withObservation :: Observation -> (C.Myb_observation -> IO r) -> IO r
withObservation o k = withMeasure count (obsMeasure o) $ \f c ->
  k
    C.zero
      { C.myb_observation't = coerce (obsTime o)
      , C.myb_observation'measure = f
      , C.myb_observation'context = c
      , C.myb_observation'count = toSize count
      , C.myb_observation'value = coerce (obsValue o)
      , C.myb_observation'sigma = coerce (obsSigma o)
      , C.myb_observation'periodic = fromIntegral (obsPeriodic o)
      , C.myb_observation'bias = fromIntegral (obsBias o)
      }
  where
    count = length (obsValue o)

-- ---- Least squares -------------------------------------------------------------------------

-- | Zeros stand for the defaults of the library; fixed and the prior sigmas, when given, have
-- an entry for every unknown (6 + parameters + shifts).
data FitOptions = FitOptions
  { fitMaxIterations :: Int
  , fitTolerance :: Double
  , fitEpsilon :: Double
  , fitFixed :: Maybe [Bool]
  , fitPriorSigma :: Maybe [Double]
  , -- | the first guess of the systematic shifts (their number is that of the list)
    fitBias :: [Double]
  }

defaultFitOptions :: FitOptions
defaultFitOptions = FitOptions 0 0 0 Nothing Nothing []

-- | The estimate: the state at the epoch, the parameters of the model (also set in it), the
-- shifts, what the iterations did, the covariance of the unknowns (by rows) and the residuals
-- (observed - computed) of every value.
data Fit = Fit
  { fitState :: State
  , fitParameters :: [Double]
  , fitBiases :: [Double]
  , fitSummary :: C.Myb_fit_result
  , fitCovariance :: [Double]
  , fitResiduals :: [Double]
  }

parameters :: Ptr C.Myb_force_model -> IO [Double]
parameters m = do
  n <- C.myb_force_model_parameter_count m
  if n == 0 then pure [] else mapM (\i -> C.myb_force_model_parameter m i >>= fmap coerce . peek) [0 .. n - 1]

-- | Fits the state at t0 (the first guess given) and the parameters of the model to the
-- observations, ordered in time away from t0.
fitOrbit :: ForceModel -> Double -> State -> [Observation] -> FitOptions -> IO (Result Fit)
fitOrbit model t0 y0 observations options
  | not (all valid observations) = pure (failure Domain)
  | otherwise = withForeignPtr (forceModelPtr model) $ \m -> do
      np <- fromSize <$> C.myb_force_model_parameter_count m
      let nb = length (fitBias options)
          unknowns = 6 + np + nb
          values = sum (map (length . obsValue) observations)
          sized xs = maybe True ((== unknowns) . length) xs
      if not (sized (fitFixed options)) || not (sized (fitPriorSigma options)) || nb > C.MYB_FIT_MAX_BIASES
        then pure (failure Domain)
        else withMany withObservation observations $ \made -> withArray made $ \obs ->
          withMaybeArray (fmap (map flag) (fitFixed options)) $ \fixed ->
            withMaybeArray (fmap (coerce :: [Double] -> [CDouble]) (fitPriorSigma options)) $ \prior ->
              withDoubles (fitBias options ++ [0]) $ \bias -> with y0 $ \py -> alloca $ \result ->
                allocaArray (unknowns * unknowns) $ \covariance -> allocaArray (max 1 values) $ \residuals -> do
                  let fitting =
                        C.zero
                          { C.myb_fit_options'max_iterations = fromIntegral (fitMaxIterations options)
                          , C.myb_fit_options'tolerance = coerce (fitTolerance options)
                          , C.myb_fit_options'epsilon = coerce (fitEpsilon options)
                          , C.myb_fit_options'fixed = fixed
                          , C.myb_fit_options'prior_sigma = prior
                          , C.myb_fit_options'bias_count = toSize nb
                          , C.myb_fit_options'bias = bias
                          }
                  code <- with fitting $ \po ->
                    C.myb_fit_orbit m (coerce t0) (castPtr py) obs (toSize (length made)) po nullPtr result covariance residuals
                  andThen code $
                    fmap Right $
                      Fit
                        <$> peek py
                        <*> parameters m
                        <*> peekDoubles nb bias
                        <*> peek result
                        <*> peekDoubles (unknowns * unknowns) covariance
                        <*> peekDoubles values residuals
  where
    flag b = if b then 1 else 0 :: CUChar

withMaybeArray :: Storable a => Maybe [a] -> (Ptr a -> IO r) -> IO r
withMaybeArray Nothing k = k nullPtr
withMaybeArray (Just xs) k = withArray xs k

-- ---- The filter ----------------------------------------------------------------------------

-- | Zeros stand for the defaults: the IAS15 precision, the process noise q (m^2 s^-3), the
-- number of systematic shifts, the instants kept for the smoother.
data FilterOptions = FilterOptions {filterEpsilon :: Double, filterProcessNoise :: Double, filterBiasCount :: Int, filterHistory :: Int}

defaultFilterOptions :: FilterOptions
defaultFilterOptions = FilterOptions 0 0 0 0

-- | A filter, what it keeps (its model) and the numbers of the parameters and the shifts.
data Filter = Filter (ForeignPtr C.Myb_filter) Kept Int Int

filterPtr :: Filter -> ForeignPtr C.Myb_filter
filterPtr (Filter p _ _ _) = p

newFilter :: ForceModel -> FilterOptions -> IO (Result Filter)
newFilter model options = withForeignPtr (forceModelPtr model) $ \m -> do
  np <- fromSize <$> C.myb_force_model_parameter_count m
  let made =
        C.zero
          { C.myb_filter_options'epsilon = coerce (filterEpsilon options)
          , C.myb_filter_options'process_noise = coerce (filterProcessNoise options)
          , C.myb_filter_options'bias_count = toSize (filterBiasCount options)
          , C.myb_filter_options'history = toSize (filterHistory options)
          }
  created <- with made $ \o -> newOwned (\f -> check <$> C.myb_filter_create m o nullPtr f) C.myb_filter_destroy
  case created of
    Left e -> pure (Left e)
    Right (fp, kept) -> do
      keepAlive kept (forceModelPtr model)
      pure (Right (Filter fp kept np (filterBiasCount options)))

-- | Starts at t0 with the state, the present parameters and the covariance of the 6 + np
-- unknowns (by rows); observations follow towards tEnd.
filterStart :: Filter -> Double -> State -> [Double] -> Double -> IO (Result ())
filterStart (Filter p _ np _) t0 y0 covariance tEnd =
  requireLength ((6 + np) * (6 + np)) covariance $ withForeignPtr p $ \f -> withVec y0 $ \py -> withDoubles covariance $ \pc ->
    check <$> C.myb_filter_start f (coerce t0) py pc (coerce tEnd)

-- | Propagates to the observation and applies it; returns the innovation (observed - predicted).
filterUpdate :: Filter -> Observation -> IO (Result [Double])
filterUpdate (Filter p _ _ _) o
  | not (valid o) = pure (failure Domain)
  | otherwise = withForeignPtr p $ \f -> withObservation o $ \made -> with made $ \po ->
      fillDoubles (length (obsValue o)) (C.myb_filter_update f po)

-- | Propagates to t without an observation.
filterPredict :: Filter -> Double -> IO (Result ())
filterPredict (Filter p _ _ _) t = withForeignPtr p $ \f -> check <$> C.myb_filter_predict f (coerce t)

-- | The time, the state and the covariance of all the unknowns (by rows).
filterGet :: Filter -> IO (Double, State, [Double])
filterGet (Filter p _ np nb) = withForeignPtr p $ \f -> alloca $ \t -> alloca $ \y -> allocaArray (n * n) $ \c -> do
  C.myb_filter_get f t (castPtr y) c
  (,,) <$> (coerce <$> peek t) <*> peek y <*> peekDoubles (n * n) c
  where
    n = 6 + np + nb

-- | The first values of the shifts (after the start, before the first observation).
filterSetBias :: Filter -> [Double] -> IO (Result ())
filterSetBias (Filter p _ _ nb) bias = requireLength nb bias $ withForeignPtr p $ \f -> withDoubles (bias ++ [0]) $ \pb -> check <$> C.myb_filter_set_bias f pb

-- | The present estimates of the shifts.
filterBias :: Filter -> IO [Double]
filterBias (Filter p _ _ nb) = withForeignPtr p $ \f -> allocaArray (max 1 nb) $ \pb -> C.myb_filter_bias f pb >> peekDoubles nb pb

-- | The instants kept so far.
filterSteps :: Filter -> IO Int
filterSteps (Filter p _ _ _) = fromSize <$> withForeignPtr p C.myb_filter_steps

-- | The backward pass of the smoother over the instants kept.
filterSmooth :: Filter -> IO (Result ())
filterSmooth (Filter p _ _ _) = withForeignPtr p $ \f -> check <$> C.myb_filter_smooth f

-- | The estimate at an instant kept: its time, the unknowns and their covariance; 'Nothing'
-- for no such instant.
filterStep :: Filter -> Int -> IO (Maybe (Double, [Double], [Double]))
filterStep fl@(Filter p _ np nb) step = do
  steps <- filterSteps fl
  if step < 0 || step >= steps
    then pure Nothing
    else withForeignPtr p $ \f -> alloca $ \t -> allocaArray n $ \u -> allocaArray (n * n) $ \c -> do
      C.myb_filter_step f (toSize step) t u c
      (\a b d -> Just (a, b, d)) <$> (coerce <$> peek t) <*> peekDoubles n u <*> peekDoubles (n * n) c
  where
    n = 6 + np + nb
