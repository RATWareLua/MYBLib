-- SPDX-License-Identifier: BSD-3-Clause
-- | The upper layer against the C interface it wraps:
--
--  * the low orbit of examples/low_orbit.c built from the upper layer gives the nodes, the
--    final state and the steps of the same calls made through "Myb.C" (the C functions) bit
--    for bit;
--  * a failed call is a Left with the status of the library and the message of a parser;
--  * terms, providers and events as Haskell functions; one that fails stops the call with its
--    status;
--  * a model keeps alive the environment and the field it was given, an environment its
--    provider and its ephemeris, an orbit its model, across collections that run finalizers;
--  * kernels and time scales, an SPK set and a body of it in an environment, a message
--    written and read in both notations, a library of ext/;
--  * the lower layer: the constants, a function that passes a structure by value, a callback
--    that takes one (the trampolines of generated/cbits);
--  * an exception in a Haskell function that C calls ends the process: it is not caught (the
--    test runs itself again as the process that dies).
--
-- No exception is thrown or caught by the checks: each one is counted, and the program ends
-- with the number of those that failed.
module Main (main) where

import Control.Concurrent (threadDelay)
import Control.Monad (forM_, unless, when)
import Data.Coerce (coerce)
import Data.IORef
import Data.Int (Int64)
import Data.List (isInfixOf, isPrefixOf)
import Foreign.C.String (withCString)
import Foreign.C.Types (CDouble (..), CInt)
import qualified Foreign.Concurrent as Concurrent
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr
import Foreign.Storable
import GHC.Float (castDoubleToWord64)
import Myb (State (..), Vec3 (..))
import qualified Myb
import qualified Myb.C as C
import qualified Myb.C.Ext.Sgp4 as Sgp4
import System.Environment (getEnvironment, getExecutablePath, lookupEnv)
import System.Exit (ExitCode (..), exitWith)
import System.FilePath ((</>))
import System.IO (hPutStrLn, stderr)
import System.Mem (performMajorGC)
import qualified System.Process as Process

-- ---- Checks --------------------------------------------------------------------------------

data Tally = Tally (IORef Int) (IORef Int)

expect :: Tally -> String -> Bool -> IO ()
expect (Tally count failed) what ok = do
  modifyIORef' count (+ 1)
  unless ok $ do
    modifyIORef' failed (+ 1)
    hPutStrLn stderr ("FAILED: " ++ what)

-- | Goes on with the value of a successful outcome; a failure is a failed check.
onRight :: Tally -> String -> Myb.Result a -> (a -> IO ()) -> IO ()
onRight tally what outcome k = case outcome of
  Right value -> expect tally what True >> k value
  Left e -> expect tally (what ++ ": " ++ Myb.errorText e) False

-- | Whether the outcome is a failure with the status.
failsWith :: Myb.Status -> Myb.Result a -> Bool
failsWith status outcome = either ((== status) . Myb.errorStatus) (const False) outcome

bits :: State -> [Word]
bits = map (fromIntegral . castDoubleToWord64) . Myb.stateList

-- | Lets the finalizers of what is no longer reachable run.
collect :: IO ()
collect = performMajorGC >> threadDelay 20000 >> performMajorGC >> threadDelay 20000

-- ---- The example ---------------------------------------------------------------------------

gm, radius, j2, day :: Double
gm = 3.986004418e14
radius = 6378137.0
j2 = 1.08263e-3
day = 86400.0

elements :: Myb.Elements
elements =
  Myb.Elements
    { Myb.semiLatusRectum = (radius + 400000.0) * (1.0 - e * e)
    , Myb.eccentricity = e
    , Myb.inclination = 51.6 * pi / 180.0
    , Myb.raan = 0.5
    , Myb.argumentOfPeriapsis = 1.0
    , Myb.trueAnomaly = 0.0
    }
  where
    e = 0.001

coefficients :: ([Double], [Double])
coefficients = ([1.0, 0.0, 0.0, -j2 / sqrt 5.0, 0.0, 0.0], replicate 6 0.0)

air :: Myb.ExponentialAtmosphere
air = Myb.ExponentialAtmosphere radius 0.0 [(400000.0, 3.0e-12, 58000.0)]

satellite :: Myb.Object
satellite = Myb.Object {Myb.objectMass = 500.0, Myb.objectArea = 4.0, Myb.objectCr = 0.0, Myb.objectCd = 2.2}

-- | The example through the C functions: the nodes, the state at the end of the day, the steps.
lowOrbitRaw :: Tally -> IO ([Double], [Double], Int64)
lowOrbitRaw tally = do
  let (c, s) = coefficients
      ok what code = expect tally ("raw: " ++ what) (code == C.MYB_OK)
  withArray (coerce c :: [CDouble]) $ \cs -> withArray (coerce s :: [CDouble]) $ \ss ->
    alloca $ \fieldOut -> alloca $ \workOut -> alloca $ \envOut -> alloca $ \planetOut -> alloca $ \modelOut -> alloca $ \orbitOut -> do
      C.myb_gravity_field_create (coerce gm) (coerce radius) 2 cs ss nullPtr fieldOut >>= ok "field"
      C.myb_gravity_work_create 2 nullPtr workOut >>= ok "work"
      field <- peek fieldOut
      work <- peek workOut
      epoch <- C.myb_time_from_seconds 0.0
      with (C.zero {C.myb_env_options'epoch = epoch}) (\options -> C.myb_env_create options nullPtr envOut) >>= ok "environment"
      env <- peek envOut
      with (C.zero {C.myb_env_body'id = 1, C.myb_env_body'gm = coerce gm}) (\body -> C.myb_env_add_body env body planetOut) >>= ok "planet"
      planet <- peek planetOut
      let atmosphere =
            C.zero
              { C.myb_exponential_atmosphere'ellipsoid = C.zero {C.myb_ellipsoid'a = coerce radius}
              , C.myb_exponential_atmosphere'count = 1
              , C.myb_exponential_atmosphere'base = [400000.0]
              , C.myb_exponential_atmosphere'density = [3.0e-12]
              , C.myb_exponential_atmosphere'scale = [58000.0]
              }
          object = C.zero {C.myb_object'mass = 500.0, C.myb_object'area = 4.0, C.myb_object'cr = 0.0, C.myb_object'cd = 2.2}
          gravity =
            C.zero
              { C.myb_field_term'env = env
              , C.myb_field_term'body = planet
              , C.myb_field_term'field = field
              , C.myb_field_term'work = work
              , C.myb_field_term'n_max = -1
              }
      with atmosphere $ \airP -> with object $ \objectP -> with gravity $ \gravityP -> do
        let drag =
              C.zero
                { C.myb_drag_term'env = env
                , C.myb_drag_term'body = planet
                , C.myb_drag_term'density = C.p'myb_density_exponential
                , C.myb_drag_term'density_context = castPtr airP
                , C.myb_drag_term'object = objectP
                }
        with drag $ \dragP -> do
          C.myb_force_model_create nullPtr modelOut >>= ok "model"
          model <- peek modelOut
          withCString "gravity" (\n -> C.myb_force_model_add model n C.p'myb_term_field C.p'myb_term_field_jacobian (castPtr gravityP) nullPtr) >>= ok "gravity"
          withCString "drag" (\n -> C.myb_force_model_add model n C.p'myb_term_drag C.p'myb_term_drag_jacobian (castPtr dragP) nullPtr) >>= ok "drag"
          pv <- with (Myb.elementsToC elements) $ \e -> alloca $ \out -> do
            C.myb_pv_from_elements (coerce gm) e out >>= ok "elements"
            peek out
          let y0 = C.myb_pv'r pv ++ C.myb_pv'v pv
          C.myb_orbit_create model nullPtr nullPtr orbitOut >>= ok "orbit"
          orbit <- peek orbitOut
          withArray y0 (\y -> C.myb_orbit_start orbit 0.0 y (coerce day)) >>= ok "start"
          -- the event of the example, y[2], read where C passes it
          equator <- C.mk'myb_event_function (\_ _ y -> peekElemOff y 2)
          let nodes k acc
                | k >= (32 :: Int) = pure (reverse acc)
                | otherwise = alloca $ \found -> alloca $ \tNode -> do
                    poke found (0 :: CInt)
                    C.myb_orbit_event orbit equator nullPtr (coerce day) 1 found tNode nullPtr >>= ok "event"
                    hit <- peek found
                    if hit == 0 then pure (reverse acc) else peek tNode >>= \t -> nodes (k + 1) (coerce t : acc)
          times <- nodes 0 []
          final <- allocaArray 6 $ \out -> do
            C.myb_orbit_state orbit (coerce day) out >>= ok "state"
            coerce <$> peekArray 6 out
          steps <- C.myb_orbit_steps orbit
          C.myb_orbit_destroy orbit
          C.myb_force_model_destroy model
          C.myb_env_destroy env
          C.myb_gravity_work_destroy work
          C.myb_gravity_field_destroy field
          freeHaskellFunPtr equator
          pure (times, final, steps)

-- | The model of the example from the upper layer; what it is built from is dropped here.
lowOrbitModel :: IO (Myb.Result Myb.ForceModel)
lowOrbitModel = do
  let (c, s) = coefficients
  made <- Myb.newGravityField gm radius 2 c s
  environment <- Myb.newEnvironment (Myb.timeFromSeconds 0.0)
  model <- Myb.newForceModel
  case (made, environment, model) of
    (Right field, Right env, Right m) -> do
      planet <- Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 1, Myb.bodyGm = gm}
      case planet of
        Left e -> pure (Left e)
        Right p -> do
          g <- Myb.addField m "gravity" env p field 0 Nothing
          d <- Myb.addDrag m "drag" env p air satellite
          pure (g >> d >> Right m)
    (Left e, _, _) -> pure (Left e)
    (_, Left e, _) -> pure (Left e)
    (_, _, Left e) -> pure (Left e)

-- | Collects the rising zeros of z up to the day.
ascendingNodes :: Tally -> Myb.Orbit -> IO [Double]
ascendingNodes tally orbit = go (0 :: Int) []
  where
    go k acc
      | k > 32 = expect tally "no more than 32 nodes" False >> pure (reverse acc)
      | otherwise = do
          next <- Myb.nextEvent orbit (\_ y -> pure (let Vec3 _ _ z = position y in z)) day Myb.Rising
          case next of
            Left e -> expect tally ("an event: " ++ Myb.errorText e) False >> pure (reverse acc)
            Right Nothing -> pure (reverse acc)
            Right (Just (t, y)) -> do
              let Vec3 _ _ z = position y
                  Vec3 _ _ vz = velocity y
              expect tally "a node is on the equator, rising" (abs z < 1e-6 && vz > 0)
              go (k + 1) (t : acc)

theUpperLayerGivesTheBitsOfTheCCalls :: Tally -> IO ()
theUpperLayerGivesTheBitsOfTheCCalls tally = do
  (rawNodes, rawFinal, rawSteps) <- lowOrbitRaw tally
  made <- lowOrbitModel
  onRight tally "the model" made $ \model -> do
    collect -- the field and the environment live in the model alone
    count <- Myb.termCount model
    expect tally "two terms" (count == 2)
    name <- Myb.termName model 1
    expect tally "the name of a term" (name == Just "drag")
    none <- Myb.termName model 2
    expect tally "no third term" (none == Nothing)
    orbit <- Myb.newOrbit model Myb.defaultOrbitOptions
    onRight tally "the orbit" orbit $ \o -> do
      collect -- and the model in the orbit
      onRight tally "the initial state" (Myb.stateFromElements gm elements) $ \y0 -> do
        started <- Myb.start o 0.0 y0 day
        onRight tally "start" started $ \() -> do
          nodes <- ascendingNodes tally o
          final <- Myb.stateAt o day
          onRight tally "the final state" final $ \y -> do
            steps <- Myb.orbitSteps o
            reached <- Myb.orbitTime o
            expect tally "at least 14 nodes" (length nodes >= 14)
            expect tally "the nodes have the bits of the C calls" (map castDoubleToWord64 nodes == map castDoubleToWord64 rawNodes)
            expect tally "the final state has the bits of the C calls" (bits y == map (fromIntegral . castDoubleToWord64) rawFinal)
            expect tally "the steps of the C calls" (steps == rawSteps)
            expect tally "the time reached" (reached == day)
            putStrLn ("the low orbit: " ++ show (length nodes) ++ " ascending nodes, " ++ show steps ++ " steps, bit for bit with the C calls")
            -- the node drifts as the first-order theory says
            let a0 = radius + 400000.0
                theory = -1.5 * sqrt (gm / a0 ^ (3 :: Int)) * j2 * (radius / Myb.semiLatusRectum elements) ^ (2 :: Int) * cos (Myb.inclination elements)
            onRight tally "the elements at the end" (Myb.elementsFromState gm y) $ \after ->
              expect tally "the drift of the node" (abs ((Myb.raan after - Myb.raan elements) / day / theory - 1.0) < 0.02)

-- ---- Failures ------------------------------------------------------------------------------

dataDir, testDataDir :: IO FilePath
dataDir = maybe ("../.." </> "data") id <$> lookupEnv "MYB_DATA_DIR"
testDataDir = maybe ("../.." </> "tests" </> "data") id <$> lookupEnv "MYB_TEST_DATA_DIR"

aFailureIsAValueWithTheStatus :: Tally -> IO ()
aFailureIsAValueWithTheStatus tally = do
  let wrong = Myb.stateFromElements gm elements {Myb.semiLatusRectum = -1.0}
  expect tally "a domain error" (failsWith Myb.Domain wrong)
  expect tally "the name of the status" (Myb.statusName Myb.Domain == "MYB_E_DOMAIN")
  expect tally "the text of an error" (either Myb.errorText (const "") wrong == "MYB_E_DOMAIN")
  -- a parser names the line
  parsed <- Myb.parseOem "CCSDS_OEM_VERS = 2.0\nMETA_START\n"
  expect tally "a format error" (failsWith Myb.Format parsed)
  expect tally "the line of the error" (either (("line 2" `isPrefixOf`) . Myb.errorMessage) (const False) parsed)
  pool <- Myb.newKernelPool
  onRight tally "a pool" pool $ \p -> do
    dir <- dataDir
    missing <- Myb.loadKernelFile p (dir </> "no-such-kernel.tk")
    expect tally "a file that is not there" (failsWith Myb.Io missing)
    nul <- Myb.loadKernelFile p "a\0b"
    expect tally "a NUL in a path" (failsWith Myb.Domain nul)
  short <- Myb.newGravityField gm radius 4 [1.0] [0.0]
  expect tally "too few coefficients" (failsWith Myb.Domain short)
  -- a body the environment does not have
  environment <- Myb.newEnvironment (Myb.timeFromSeconds 0.0)
  model <- Myb.newForceModel
  case (environment, model) of
    (Right env, Right m) -> do
      nobody <- Myb.addThirdBody m "nobody" env 3
      expect tally "a third body that is not there" (failsWith Myb.Domain nobody)
      state <- Myb.bodyState env 0.0 0
      expect tally "the state of a body that is not there" (failsWith Myb.Domain state)
      switched <- Myb.setEnabled m 0 False
      expect tally "no term to switch" (not switched)
    _ -> expect tally "an environment and a model" False

-- ---- Haskell functions ---------------------------------------------------------------------

termsAndEventsAreHaskellFunctions :: Tally -> IO ()
termsAndEventsAreHaskellFunctions tally = do
  let x = Vec3 7.0e6 1.0e5 (-2.0e5)
      v = Vec3 0.0 7.5e3 1.0e3
      push = Vec3 1.0e-6 (-2.0e-6) 3.0e-6
  calls <- newIORef (0 :: Int)
  made <- Myb.newForceModel
  onRight tally "a model" made $ \model -> do
    _ <- Myb.addPointMass model "planet" gm
    added <- Myb.addTerm model "push" (\_ _ _ -> modifyIORef' calls (+ 1) >> pure (Right push))
    onRight tally "a term of Haskell" added $ \index -> do
      both <- Myb.acceleration model 0.0 x v
      off <- Myb.setEnabled model index False
      isOff <- Myb.termEnabled model index
      planet <- Myb.acceleration model 0.0 x v
      on <- Myb.setEnabled model index True
      isOn <- Myb.termEnabled model index
      noSuch <- Myb.termEnabled model 5
      expect tally "the term switched off and on" (off && on && isOff == Just False && isOn == Just True && noSuch == Nothing)
      case (both, planet) of
        (Right (Vec3 a b c), Right (Vec3 a' b' c')) -> do
          let Vec3 p q r = push
          expect tally "the push is in the sum" (abs (a - a' - p) < 1e-14 && abs (b - b' - q) < 1e-14 && abs (c - c' - r) < 1e-14)
        _ -> expect tally "two accelerations" False
      once <- readIORef calls
      expect tally "the term was called once" (once == 1)
      -- an orbit with the push: an event, many states, the transition matrix, an impulse
      orbit <- Myb.newOrbit model Myb.defaultOrbitOptions {Myb.variational = True}
      onRight tally "an orbit with the transition matrix" orbit $ \o ->
        onRight tally "the initial state" (Myb.stateFromElements gm elements) $ \y0 -> do
          _ <- Myb.start o 0.0 y0 6000.0
          crossing <- Myb.nextEvent o (\_ y -> pure (let Vec3 _ _ z = position y in z)) 6000.0 Myb.Falling
          onRight tally "an event" crossing $ \found -> case found of
            Just (_, at) -> do
              let Vec3 _ _ z = position at
                  Vec3 _ _ vz = velocity at
              expect tally "the descending node" (abs z < 1e-6 && vz < 0)
            Nothing -> expect tally "the descending node of the first revolution" False
          withPartials <- Myb.partials o 5000.0
          plain <- Myb.stateAt o 5000.0
          case (withPartials, plain) of
            (Right (y, phi), Right y') -> do
              expect tally "the state with the matrix is the state" (bits y == bits y')
              expect tally "a matrix of 36 numbers" (length phi == 36 && take 1 phi /= [0] && not (any isNaN phi))
            _ -> expect tally "the partials and the state" False
          many <- Myb.statesAt o [5100.0, 5200.0, 5300.0]
          expect tally "three states" (either (const False) ((== 3) . length) many)
          kicked <- Myb.kick o (Vec3 0.0 1.0 0.0)
          expect tally "an impulse" (kicked == Right ())
          more <- readIORef calls
          expect tally "the term was called by the propagation" (more > 100)
  -- DOP853 agrees with IAS15 to its tolerance, and both with Kepler
  plain <- Myb.newForceModel
  onRight tally "a plain model" plain $ \model -> do
    _ <- Myb.addPointMass model "planet" gm
    ias15 <- Myb.newOrbit model Myb.defaultOrbitOptions
    dop853 <- Myb.newOrbit model Myb.defaultOrbitOptions {Myb.integrator = Myb.Dop853, Myb.rtol = 1e-11, Myb.atol = 1e-6}
    case (ias15, dop853, Myb.stateFromElements gm elements) of
      (Right a, Right b, Right y0) -> do
        _ <- Myb.start a 0.0 y0 3000.0
        _ <- Myb.start b 0.0 y0 3000.0
        ya <- Myb.stateAt a 3000.0
        yb <- Myb.stateAt b 3000.0
        case (ya, yb, Myb.keplerPropagate gm y0 3000.0) of
          (Right p, Right q, Right k) -> do
            let close tol s s' = and (zipWith (\u w -> abs (u - w) < tol) (take 3 (Myb.stateList s)) (take 3 (Myb.stateList s')))
            expect tally "DOP853 and IAS15" (close 1e-2 p q)
            expect tally "IAS15 and Kepler" (close 1e-6 p k)
          _ -> expect tally "three states" False
      _ -> expect tally "two orbits" False
  -- a term that fails: its status comes back from the evaluation
  failing <- Myb.newForceModel
  onRight tally "a model" failing $ \model -> do
    _ <- Myb.addTerm model "refuses" (\_ _ _ -> pure (Left Myb.Range))
    refused <- Myb.acceleration model 0.0 x v
    expect tally "the status of the term" (failsWith Myb.Range refused)

aProviderIsAHaskellFunction :: Tally -> IO ()
aProviderIsAHaskellFunction tally = do
  let rate = 2.0e-7
      distance = 1.5e11
      sun t =
        let angle = rate * t
         in State (Vec3 (distance * cos angle) (distance * sin angle) 0.0) (Vec3 (-distance * rate * sin angle) (distance * rate * cos angle) 0.0)
      x = Vec3 7.0e6 1.0e5 (-2.0e5)
      v = Vec3 0.0 7.5e3 1.0e3
  model <- Myb.newForceModel
  index <- do
    made <- Myb.newEnvironment (Myb.timeFromSeconds 0.0)
    case (made, model) of
      (Right env, Right m) -> do
        _ <- Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 399, Myb.bodyGm = gm}
        added <- Myb.addWith env Myb.defaultBody {Myb.bodyId = 10, Myb.bodyGm = 1.32712440018e20} (pure . Right . sun)
        count <- Myb.bodyCount env
        expect tally "two bodies" (count == 2)
        found <- Myb.findBody env 10
        expect tally "found by its identifier" (found == added)
        lost <- Myb.findBody env 11
        expect tally "not found" (failsWith Myb.NotFound lost)
        case added of
          Right i -> do
            at <- Myb.bodyState env 1000.0 i
            expect tally "the state of the provider" (fmap bits at == Right (bits (sun 1000.0)))
            _ <- Myb.addPointMass m "planet" gm
            third <- Myb.addThirdBody m "sun" env i
            expect tally "the sun as a third body" (third == Right 1)
          Left e -> expect tally ("a provider: " ++ Myb.errorText e) False
        pure added
      _ -> pure (Left (Myb.Error Myb.NoMemory ""))
  onRight tally "the provider" index $ \_ -> onRight tally "the model" model $ \m -> do
    collect -- the model keeps the environment, and the environment its function
    withSun <- Myb.acceleration m 0.0 x v
    _ <- Myb.setEnabled m 1 False
    alone <- Myb.acceleration m 0.0 x v
    case (withSun, alone) of
      (Right (Vec3 a b c), Right (Vec3 a' b' c')) -> do
        let tide = sqrt ((a - a') ^ (2 :: Int) + (b - b') ^ (2 :: Int) + (c - c') ^ (2 :: Int))
        expect tally "the tide of the sun" (tide > 1e-7 && tide < 1e-5) -- 2 GM r / d^3 at 7000 km is 5e-7 m/s^2
      _ -> expect tally "two accelerations" False
  -- a provider that fails stops the evaluation with its status
  broken <- Myb.newEnvironment (Myb.timeFromSeconds 0.0)
  onRight tally "an environment" broken $ \env -> do
    _ <- Myb.addOrigin env Myb.defaultBody
    lost <- Myb.addWith env Myb.defaultBody {Myb.bodyId = 5} (\_ -> pure (Left Myb.Range))
    onRight tally "a provider that fails" lost $ \i -> do
      state <- Myb.bodyState env 0.0 i
      expect tally "the status of the provider" (failsWith Myb.Range state)

-- | A model keeps its environment: the environment is not finalized while the model lives,
-- and both are once neither is reachable (the finalizers are watched through IORefs).
ownershipFollowsTheModel :: Tally -> IO ()
ownershipFollowsTheModel tally = do
  envGone <- newIORef False
  modelGone <- newIORef False
  made <- build envGone modelGone
  onRight tally "a model over an environment" made $ \model -> do
    collect
    early <- readIORef envGone
    expect tally "the environment lives with the model" (not early)
    a <- Myb.acceleration model 0.0 (Vec3 7.0e6 0 0) (Vec3 0 7.5e3 0)
    expect tally "the model still works" (either (const False) (const True) a)
  collect
  collect
  gone <- (,) <$> readIORef modelGone <*> readIORef envGone
  expect tally "the model, then the environment, finalized" (gone == (True, True))
  where
    build envGone modelGone = do
      environment <- Myb.newEnvironment (Myb.timeFromSeconds 0.0)
      model <- Myb.newForceModel
      case (environment, model) of
        (Right env, Right m) -> do
          Concurrent.addForeignPtrFinalizer (Myb.environmentPtr env) (writeIORef envGone True)
          Concurrent.addForeignPtrFinalizer (Myb.forceModelPtr m) (writeIORef modelGone True)
          _ <- Myb.addOrigin env Myb.defaultBody {Myb.bodyGm = gm}
          moon <- Myb.addWith env Myb.defaultBody {Myb.bodyId = 301, Myb.bodyGm = 4.9e12} (\_ -> pure (Right (State (Vec3 3.8e8 0 0) (Vec3 0 1.0e3 0))))
          case moon of
            Left e -> pure (Left e)
            Right i -> fmap (const m) <$> Myb.addThirdBody m "moon" env i
        (Left e, _) -> pure (Left e)
        (_, Left e) -> pure (Left e)

-- | The pressure of sunlight and the correction of relativity: terms whose structures hold an
-- array (the occulters) and a nested structure (the parameters of the theory).
lightAndRelativity :: Tally -> IO ()
lightAndRelativity tally = do
  let sun = State (Vec3 1.5e11 0.0 0.0) (Vec3 0.0 0.0 0.0)
      c = 299792458.0
  made <- Myb.newEnvironment (Myb.timeFromSeconds 0.0)
  onRight tally "an environment" made $ \env -> do
    earth <- Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 399, Myb.bodyGm = gm, Myb.bodyRadius = radius}
    star <- Myb.addWith env Myb.defaultBody {Myb.bodyId = 10, Myb.bodyGm = 1.32712440018e20, Myb.bodyRadius = 6.957e8, Myb.bodyLuminosity = 3.828e26} (\_ -> pure (Right sun))
    case (earth, star) of
      (Right e, Right s) -> do
        tooMany <- Myb.newForceModel
        onRight tally "a model" tooMany $ \m -> do
          wrong <- Myb.addRadiation m "light" env s satellite c (replicate 99 e)
          expect tally "too many occulters" (failsWith Myb.Domain wrong)
        light <- magnitude (\m -> Myb.addRadiation m "light" env s satellite {Myb.objectCr = 1.3} c [e])
        -- L / (4 pi d^2 c) cr A / m at 1.5e11 m: 4.7e-8 m/s^2
        expect tally "the pressure of sunlight" (maybe False (\a -> a > 3e-8 && a < 7e-8) light)
        relativity <- magnitude (\m -> Myb.addSchwarzschild m "relativity" env e c)
        -- GM / (c^2 r) GM / r^2 at 7000 km: 5e-9 m/s^2
        expect tally "the correction of relativity" (maybe False (\a -> a > 1e-9 && a < 2e-8) relativity)
        ppn <- magnitude (\m -> Myb.addSchwarzschildPpn m "relativity" env e c 1.0 1.0)
        expect tally "the same with beta = gamma = 1" (ppn == relativity)
      _ -> expect tally "two bodies" False
  where
    x0 = Vec3 7.0e6 1.0e5 (-2.0e5)
    v0 = Vec3 0.0 7.5e3 1.0e3
    -- the acceleration of a model of the one term
    magnitude add = do
      made <- Myb.newForceModel
      case made of
        Left _ -> pure Nothing
        Right m -> do
          added <- add m
          a <- Myb.acceleration m 0.0 x0 v0
          pure $ case (added, a) of
            (Right _, Right (Vec3 p q r)) -> Just (sqrt (p * p + q * q + r * r))
            _ -> Nothing

-- ---- Files ---------------------------------------------------------------------------------

timeScalesAndAnSpkSet :: Tally -> IO ()
timeScalesAndAnSpkSet tally = do
  dir <- dataDir
  tests <- testDataDir
  made <- Myb.newKernelPool
  system <- case made of
    Left e -> pure (Left e)
    Right pool -> do
      loaded <- mapM (Myb.loadKernelFile pool . (dir </>)) ["time" </> "time-scales.tk", "time" </> "tdb-tt-fb1990.tk", "time" </> "utc-tai.tk"]
      case sequence loaded of
        Left e -> pure (Left e)
        Right _ -> Myb.newTimeSystem pool -- the system keeps nothing of the pool
  onRight tally "a time system" system $ \ts -> do
    collect
    let label = Myb.Calendar {Myb.year = 2026, Myb.month = 11, Myb.day = 11, Myb.hour = 12, Myb.minute = 30, Myb.second = 15.25}
    instant <- Myb.fromUtc ts label "TT"
    onRight tally "an instant of TT" instant $ \t -> do
      back <- Myb.toUtc ts t "TT"
      onRight tally "its label" back $ \l -> do
        expect tally "the label again" ((Myb.year l, Myb.month l, Myb.day l, Myb.hour l, Myb.minute l) == (2026, 11, 11, 12, 30))
        expect tally "the second again" (abs (Myb.second l - 15.25) < 1e-9)
      tai <- Myb.convertTime ts t "TT" "TAI"
      onRight tally "TAI" tai $ \t' -> expect tally "TT - TAI" (abs (Myb.seconds t - Myb.seconds t' - 32.184) < 1e-6)
      nowhere <- Myb.convertTime ts t "TT" "NO SUCH SCALE"
      expect tally "a scale that is not there" (failsWith Myb.NotFound nowhere)
  set <- Myb.newEphemeris
  onRight tally "an SPK set" set $ \eph -> do
    loaded <- Myb.loadEphemerisFile eph (tests </> "spk" </> "t02_03.bsp")
    expect tally "an SPK file" (loaded == Right ())
    missing <- Myb.loadEphemerisFile eph (tests </> "spk" </> "missing.bsp")
    expect tally "an SPK file that is not there" (either (const True) (const False) missing)
    segments <- Myb.segments eph
    case segments of
      [] -> expect tally "the segments of the set" False
      (first : _) -> do
        let middle = Myb.timeFromSeconds (0.5 * (Myb.segmentStart first + Myb.segmentStop first))
            target = Myb.segmentTarget first
            center = Myb.segmentCenter first
        there <- Myb.ephemerisState eph target center middle
        raw <- withForeignPtr (Myb.ephemerisPtr eph) $ \p -> allocaArray 6 $ \out -> do
          code <- C.myb_spk_state p target center (Myb.timeToC middle) out nullPtr
          values <- peekArray 6 out
          pure (code, coerce values :: [Double])
        expect tally "the state of the set has the bits of the C call" (fmap Myb.stateList there == Right (snd raw) && fst raw == C.MYB_OK)
        late <- Myb.ephemerisState eph target center (Myb.timeFromSeconds (Myb.segmentStop first + 1.0e9))
        expect tally "outside the segments" (failsWith Myb.Range late)
        env <- Myb.newEnvironment middle
        onRight tally "an environment at the middle" env $ \e -> do
          _ <- Myb.addOrigin e Myb.defaultBody {Myb.bodyId = fromIntegral center, Myb.bodyGm = gm}
          index <- Myb.addFromEphemeris e Myb.defaultBody {Myb.bodyId = fromIntegral target, Myb.bodyGm = 1.0} eph target center 0
          onRight tally "a body of the set" index $ \i -> do
            state <- Myb.bodyState e 0.0 i
            expect tally "the body has the state of the set" (fmap bits state == fmap bits there)

aMessageInBothNotations :: Tally -> IO ()
aMessageInBothNotations tally = do
  let segment =
        Myb.OemSegment
          { Myb.objectName = "A & B"
          , Myb.objectId = "2026-001A"
          , Myb.centerName = "EARTH"
          , Myb.refFrame = "EME2000"
          , Myb.timeSystem = "UTC"
          , Myb.interpolation = Just ("HERMITE", 7)
          , Myb.epochs =
              [ Myb.Calendar 2026 3 1 12 0 0.0
              , Myb.Calendar 2026 3 1 12 1 0.5
              ]
          , Myb.states =
              [ State (Vec3 7.0e6 1.0 (-2.5)) (Vec3 0.25 7.5e3 (-1.0e3))
              , State (Vec3 (7.0e6 + 15.0) 4.5e5 (-6.0e4)) (Vec3 (-480.0) 7.4e3 (-990.0))
              ]
          }
      created = Myb.Calendar 2026 3 2 0 0 0.0
  forM_ [Myb.Kvn, Myb.Xml] $ \notation -> do
    text <- Myb.writeOem "myb" created [segment] notation
    onRight tally ("a message in " ++ show notation) text $ \t -> do
      expect tally "the notation" (("<" `isPrefixOf` t) == (notation == Myb.Xml))
      parsed <- Myb.parseOem t
      onRight tally "the message read" parsed $ \message -> do
        counts <- (,,) <$> Myb.oemSegmentCount message <*> Myb.oemStateCount message 0 <*> Myb.oemStateCount message 1
        expect tally "one segment of two states" (counts == (1, 2, 0))
        name <- Myb.oemMeta message 0 "OBJECT_NAME"
        expect tally "the name of the object" (name == Just "A & B")
        degree <- Myb.oemMeta message 0 "INTERPOLATION_DEGREE"
        expect tally "the degree" (degree == Just "7")
        nothing <- Myb.oemMeta message 0 "NO_SUCH_KEY"
        expect tally "a keyword that is not there" (nothing == Nothing)
        line <- Myb.oemState message 0 1
        case line of
          Just (epoch, y) -> do
            expect tally "the epoch of the second line" ((Myb.minute epoch, Myb.second epoch) == (1, 0.5))
            expect tally "the state of the second line, bit for bit" (bits y == bits (Myb.states segment !! 1))
          Nothing -> expect tally "the second line" False
        past <- Myb.oemState message 0 2
        expect tally "no third line" (past == Nothing)
  empty <- Myb.writeOem "myb" created [] Myb.Kvn
  expect tally "a message without segments" (failsWith Myb.Domain empty)
  uneven <- Myb.writeOem "myb" created [segment {Myb.states = [State (Vec3 0 0 0) (Vec3 0 0 0)]}] Myb.Kvn
  expect tally "epochs and states of different numbers" (failsWith Myb.Domain uneven)

-- ---- The lower layer -----------------------------------------------------------------------

theLowerLayer :: Tally -> IO ()
theLowerLayer tally = do
  version <- C.myb_version_number
  expect tally "the library is the one declared" (version == C.MYB_VERSION_NUMBER)
  t <- C.myb_time_from_seconds 12.5 -- a structure returned by value, through a shim
  expect tally "a structure by value" (C.myb_time's t == 12 && C.myb_time'f t == 0.5)
  back <- C.myb_time_seconds t -- and passed by value
  expect tally "a structure passed by value" (back == 12.5)
  -- a Haskell function as a callback that takes a structure by value: the trampoline, and
  -- the C function that calls one (call')
  (function, context) <- C.mk'myb_state_function $ \_ instant out -> do
    s <- C.myb_time_seconds instant
    pokeElemOff out 0 s
    pure C.MYB_OK
  status <- allocaArray 6 $ \out -> do
    code <- C.call'myb_state_function function context (C.zero {C.myb_time's = 7, C.myb_time'f = 0.25}) out
    first <- peekElemOff out 0
    pure (code == C.MYB_OK && first == 7.25)
  expect tally "a callback that takes a structure by value" status
  freeHaskellFunPtr (castPtrToFunPtr context)
  -- a structure goes to C and back unchanged
  let body = C.zero {C.myb_env_body'id = 42, C.myb_env_body'gm = 1.5, C.myb_env_body'limb = [0.25, 0.5]}
  again <- with body peek
  expect tally "a structure through its Storable instance" (again == body)
  expect tally "the constants of an enumeration" (C.MYB_E_FORMAT /= C.MYB_OK && C.MYB_ODM_XML /= C.MYB_ODM_KVN)
  expect tally "a constant of an anonymous enumeration" ((C.MYB_ENV_MAX_OCCULTERS :: Int) > 0)
  -- a library of ext/
  let line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
      line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
  parsed <- withCString line1 $ \l1 -> withCString line2 $ \l2 -> with (C.zero :: Sgp4.Myb_tle) $ \tle -> with (C.zero :: C.Myb_diag) $ \diag -> do
    code <- Sgp4.myb_tle_parse l1 l2 0 tle diag
    (,) code <$> peek tle
  expect tally "a TLE through myb:sgp4" (fst parsed == C.MYB_OK && Sgp4.myb_tle'catalog_number (snd parsed) == 5)
  expect tally "its inclination" (abs (coerce (Sgp4.myb_tle'inclination (snd parsed)) - 34.2682 :: Double) < 1e-12)

-- ---- An exception in a Haskell function called from C ---------------------------------------

-- | The variable that makes this program the process that throws.
dies :: String
dies = "MYB_TEST_THROW_IN"

throwIn :: String -> IO ()
throwIn what = do
  made <- Myb.newForceModel
  case made of
    Left _ -> pure ()
    Right model -> do
      _ <- Myb.addPointMass model "centre" gm
      case what of
        "term" -> do
          _ <- Myb.addTerm model "broken" (\_ _ _ -> ioError (userError "in a term"))
          _ <- Myb.acceleration model 0.0 (Vec3 7.0e6 0 0) (Vec3 0 7.5e3 0)
          pure ()
        _ -> do
          orbit <- Myb.newOrbit model Myb.defaultOrbitOptions
          case orbit of
            Left _ -> pure ()
            Right o -> do
              _ <- Myb.start o 0.0 (State (Vec3 7.0e6 0 0) (Vec3 0 7.5e3 1.0e3)) 6000.0
              _ <- Myb.nextEvent o (\_ _ -> ioError (userError "in an event")) 6000.0 Myb.Any
              pure ()
  putStrLn "the process went on"

anExceptionInAFunctionCalledFromCEndsTheProcess :: Tally -> IO ()
anExceptionInAFunctionCalledFromCEndsTheProcess tally = do
  this <- getExecutablePath
  inherited <- getEnvironment
  forM_ ["term", "event"] $ \what -> do
    let child = (Process.proc this []) {Process.env = Just ((dies, what) : filter ((/= dies) . fst) inherited)}
    (code, out, err) <- Process.readCreateProcessWithExitCode child ""
    let message = "in " ++ (if what == "term" then "a " else "an ") ++ what
    expect tally (what ++ ": the message of the exception (" ++ err ++ ")") (message `isInfixOf` err)
    expect tally (what ++ ": the process ended") (code /= ExitSuccess && not ("went on" `isInfixOf` out))

-- ---- Main ----------------------------------------------------------------------------------

main :: IO ()
main = do
  throwing <- lookupEnv dies
  case throwing of
    Just what -> throwIn what
    Nothing -> do
      tally@(Tally count failed) <- Tally <$> newIORef 0 <*> newIORef 0
      theLowerLayer tally
      theUpperLayerGivesTheBitsOfTheCCalls tally
      aFailureIsAValueWithTheStatus tally
      termsAndEventsAreHaskellFunctions tally
      aProviderIsAHaskellFunction tally
      lightAndRelativity tally
      ownershipFollowsTheModel tally
      timeScalesAndAnSpkSet tally
      aMessageInBothNotations tally
      anExceptionInAFunctionCalledFromCEndsTheProcess tally
      collect
      n <- readIORef count
      bad <- readIORef failed
      putStrLn (show n ++ " checks, " ++ show bad ++ " failed")
      when (bad > 0) (exitWith (ExitFailure 1))
