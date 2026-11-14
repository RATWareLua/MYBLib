-- SPDX-License-Identifier: BSD-3-Clause
-- | Optimal transfers with low thrust (myb/optimal.h): the indirect method, the solver of
-- nonlinear programs, the direct method of Sims and Flanagan and chains of its legs joined by
-- flybys, powered or not. The problems and results are the records of "Myb.C"; a force model
-- for the field is passed beside them ('Nothing': the point mass mu).
module Myb.Optimal
  ( -- * The indirect method
    C.Myb_indirect_problem (..)
  , C.Myb_indirect_solution (..)
  , C.Myb_indirect_point (..)
  , indirectSolve
  , indirectSample
    -- * Nonlinear programs
  , NlpValues (..)
  , NlpProblem (..)
  , C.Myb_nlp_options (..)
  , C.Myb_nlp_result (..)
  , nlpSolve
    -- * Sims and Flanagan
  , C.Myb_sims_flanagan (..)
  , C.Myb_sims_flanagan_result (..)
  , simsFlanaganMismatch
  , simsFlanaganSolve
    -- * Chains with flybys
  , chainFlyby
  , ChainBody (..)
  , Chain (..)
  , C.Myb_chain_result (..)
  , chainSize
  , chainGuess
  , chainEvaluate
  , chainSolve
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, pokeArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, freeHaskellFunPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Force (ForceModel, forceModelPtr)
import Myb.Internal
import Myb.Types

withField :: Maybe ForceModel -> (Ptr C.Myb_force_model -> IO a) -> IO a
withField Nothing k = k nullPtr
withField (Just m) k = withForeignPtr (forceModelPtr m) k

withNlpOptions :: Maybe C.Myb_nlp_options -> (Ptr C.Myb_nlp_options -> IO a) -> IO a
withNlpOptions Nothing k = k nullPtr
withNlpOptions (Just o) k = with o k

-- | Solves the problem of the indirect method; guess: the eight costates of a solution nearby.
indirectSolve :: C.Myb_indirect_problem -> Maybe ForceModel -> Maybe [Double] -> IO (Result C.Myb_indirect_solution)
indirectSolve problem field guess
  | maybe False ((/= 8) . length) guess = pure (failure Domain)
  | otherwise = withField field $ \f -> with problem {C.myb_indirect_problem'field = f} $ \pp -> withGuess $ \pg ->
      fillZero (C.myb_indirect_solve pp pg nullPtr)
  where
    withGuess k = maybe (k nullPtr) (`withDoubles` k) guess

-- | The solution at times in increasing order between t0 and the arrival.
indirectSample :: C.Myb_indirect_problem -> Maybe ForceModel -> C.Myb_indirect_solution -> [Double] -> IO (Result [C.Myb_indirect_point])
indirectSample problem field solution times =
  withField field $ \f -> with problem {C.myb_indirect_problem'field = f} $ \pp -> with solution $ \ps -> withDoublesLen times $ \n pt ->
    fillArray n (C.myb_indirect_sample pp ps pt (toSize n) nullPtr)

-- | What the function of a program gives at a point: the objective, the constraints
-- (equalities first) and, when asked, the gradient and the Jacobian of the constraints (by rows).
data NlpValues = NlpValues {nlpObjective :: !Double, nlpConstraints :: ![Double], nlpGradient :: ![Double], nlpJacobian :: ![Double]}
  deriving (Eq, Show)

-- | minimise f(x) with c_i(x) = 0 for the equalities, c_i(x) <= 0 for the inequalities after
-- them, lower <= x <= upper; the function is told (True) when the derivatives are wanted.
data NlpProblem = NlpProblem
  { nlpUnknowns :: Int
  , nlpEqualities :: Int
  , nlpInequalities :: Int
  , nlpLower :: Maybe [Double]
  , nlpUpper :: Maybe [Double]
  , nlpFunction :: [Double] -> Bool -> IO (Either Status NlpValues)
  }

-- | The solution from the first guess: x, the multipliers of the constraints, what the solver
-- did.
nlpSolve :: NlpProblem -> Maybe C.Myb_nlp_options -> [Double] -> IO (Result ([Double], [Double], C.Myb_nlp_result))
nlpSolve problem options x0
  | length x0 /= n || maybe False ((/= n) . length) (nlpLower problem) || maybe False ((/= n) . length) (nlpUpper problem) = pure (failure Domain)
  | otherwise = do
      function <- C.mk'myb_nlp_function act
      done <- withBounds (nlpLower problem) $ \lower -> withBounds (nlpUpper problem) $ \upper -> do
        let made =
              C.zero
                { C.myb_nlp_problem'n = toSize n
                , C.myb_nlp_problem'equalities = toSize (nlpEqualities problem)
                , C.myb_nlp_problem'inequalities = toSize (nlpInequalities problem)
                , C.myb_nlp_problem'lower = lower
                , C.myb_nlp_problem'upper = upper
                , C.myb_nlp_problem'function = function
                }
        with made $ \pp -> withNlpOptions options $ \po -> withDoubles x0 $ \px -> allocaArray (max 1 m) $ \multipliers -> alloca $ \result -> do
          code <- C.myb_nlp_solve pp po px multipliers nullPtr result
          andThen code (Right <$> ((,,) <$> peekDoubles n px <*> peekDoubles m multipliers <*> peek result))
      freeHaskellFunPtr function
      pure done
  where
    n = nlpUnknowns problem
    m = nlpEqualities problem + nlpInequalities problem
    withBounds Nothing k = k nullPtr
    withBounds (Just b) k = withDoubles b k
    act _ px f c gradient jacobian = do
      x <- peekDoubles n px
      given <- nlpFunction problem x (gradient /= nullPtr)
      case given of
        Left status -> pure (statusCode status)
        Right v
          | length (nlpConstraints v) /= m
              || (gradient /= nullPtr && (length (nlpGradient v) /= n || length (nlpJacobian v) /= m * n)) ->
              pure C.MYB_E_DOMAIN
          | otherwise -> do
              poke f (coerce (nlpObjective v))
              pokeArray c (coerce (nlpConstraints v))
              if gradient /= nullPtr
                then pokeArray gradient (coerce (nlpGradient v)) >> pokeArray jacobian (coerce (nlpJacobian v))
                else pure ()
              pure C.MYB_OK

-- | The mismatch at the meeting point (6: m, m/s) for the throttles (3 per segment), the mass
-- at the arrival and the Jacobian of the mismatch by the throttles (6 x 3 segments, by rows).
simsFlanaganMismatch :: C.Myb_sims_flanagan -> Maybe ForceModel -> [Double] -> IO (Result ([Double], Double, [Double]))
simsFlanaganMismatch leg field throttles =
  requireLength (3 * segments) throttles $ withField field $ \f -> with leg {C.myb_sims_flanagan'field = f} $ \pl -> withDoubles throttles $ \pt ->
    allocaArray 6 $ \mismatch -> alloca $ \mass -> allocaArray (18 * segments) $ \jacobian -> do
      code <- C.myb_sims_flanagan_mismatch pl pt mismatch mass jacobian
      andThen code (Right <$> ((,,) <$> peekDoubles 6 mismatch <*> (coerce <$> peek mass) <*> peekDoubles (18 * segments) jacobian))
  where
    segments = fromSize (C.myb_sims_flanagan'segments leg)

-- | The throttles that join the two states with the least propellant, from the first guess
-- (all zeros: the one built in).
simsFlanaganSolve :: C.Myb_sims_flanagan -> Maybe ForceModel -> Maybe C.Myb_nlp_options -> [Double] -> IO (Result ([Double], C.Myb_sims_flanagan_result))
simsFlanaganSolve leg field options throttles =
  requireLength (3 * segments) throttles $ withField field $ \f -> with leg {C.myb_sims_flanagan'field = f} $ \pl -> withNlpOptions options $ \po -> withDoubles throttles $ \pt -> alloca $ \result -> do
    code <- C.myb_sims_flanagan_solve pl po pt nullPtr result
    andThen code (Right <$> ((,) <$> peekDoubles (3 * segments) pt <*> peek result))
  where
    segments = fromSize (C.myb_sims_flanagan'segments leg)

-- | A flyby with an impulse at the periapsis: the impulse dv and the turn of the excess
-- velocity for the excess speeds coming in and going out and the periapsis radius, with their
-- partial derivatives (d dv, then d turn, by (v_in, v_out, periapsis)).
chainFlyby :: Double -> Double -> Double -> Double -> Result (Double, Double, [Double])
chainFlyby gm vIn vOut periapsis = pureCall $ alloca $ \dv -> alloca $ \turn -> allocaArray 6 $ \partials -> do
  code <- C.myb_chain_flyby (coerce gm) (coerce vIn) (coerce vOut) (coerce periapsis) dv turn partials
  andThen code (Right <$> ((,,) <$> (coerce <$> peek dv) <*> (coerce <$> peek turn) <*> peekDoubles 6 partials))

-- | A body of a chain: where it is, its GM and the least distance of a flyby, and the
-- effective exhaust velocity of an impulse at its periapsis (0: an unpowered flyby).
data ChainBody = ChainBody {chainBodySource :: StateSource, chainBodyGm :: Double, chainBodyPeriapsis :: Double, chainBodyImpulseExhaust :: Double}

-- | Legs of Sims and Flanagan from body to body: the centre (and the units of the solver), the
-- field of the legs, the epoch of the time 0, the engine, the first mass, the segments of each
-- leg, the bodies (legs + 1), the bounds of the departure and of each flight (s), the largest
-- excess speeds at the departure and at the arrival (m/s).
data Chain = Chain
  { chainMu :: Double
  , chainField :: Maybe ForceModel
  , chainEpoch :: Time
  , chainThrust :: Double
  , chainExhaust :: Double
  , chainMass0 :: Double
  , chainSegments :: [Int]
  , chainBodies :: [ChainBody]
  , chainDeparture :: (Double, Double)
  , chainFlight :: [(Double, Double)]
  , chainLaunch :: Double
  , chainArrival :: Double
  }

-- | The chain as C reads it, for one call.
withChain :: Chain -> (Ptr C.Myb_chain -> IO (Result a)) -> IO (Result a)
withChain chain k
  | legs < 1 || legs > C.MYB_CHAIN_MAX_LEGS || length (chainBodies chain) /= legs + 1 || length (chainFlight chain) /= legs = pure (failure Domain)
  | otherwise =
      withField (chainField chain) $ \f ->
        withMany (\b go -> withStateSource (chainBodySource b) (\fn ctx -> go (b, fn, ctx))) (chainBodies chain) $ \bodies ->
          with
            C.zero
              { C.myb_chain'mu = coerce (chainMu chain)
              , C.myb_chain'field = f
              , C.myb_chain'epoch = timeToC (chainEpoch chain)
              , C.myb_chain'thrust = coerce (chainThrust chain)
              , C.myb_chain'exhaust = coerce (chainExhaust chain)
              , C.myb_chain'mass0 = coerce (chainMass0 chain)
              , C.myb_chain'legs = toSize legs
              , C.myb_chain'segments = map toSize (chainSegments chain)
              , C.myb_chain'bodies =
                  [ C.zero {C.myb_chain_body'state = fn, C.myb_chain_body'context = ctx, C.myb_chain_body'gm = coerce (chainBodyGm b), C.myb_chain_body'periapsis = coerce (chainBodyPeriapsis b)}
                  | (b, fn, ctx) <- bodies
                  ]
              , C.myb_chain'departure = let (a, b) = chainDeparture chain in coerce [a, b]
              , C.myb_chain'flight = coerce (concat [[a, b] | (a, b) <- chainFlight chain])
              , C.myb_chain'launch = coerce (chainLaunch chain)
              , C.myb_chain'arrival = coerce (chainArrival chain)
              , C.myb_chain'impulse_exhaust = coerce (map chainBodyImpulseExhaust (chainBodies chain))
              }
            k
  where
    legs = length (chainSegments chain)

-- | The numbers of the unknowns, of the equalities and of the inequalities.
chainSize :: Chain -> IO (Result (Int, Int, Int))
chainSize chain = withChain chain $ \pc -> alloca $ \u -> alloca $ \e -> alloca $ \i -> do
  code <- C.myb_chain_size pc u e i
  andThen code (Right <$> ((,,) <$> (fromSize <$> peek u) <*> (fromSize <$> peek e) <*> (fromSize <$> peek i)))

-- | The first guess built in.
chainGuess :: Chain -> IO (Result [Double])
chainGuess chain = do
  sized <- chainSize chain
  case sized of
    Left e -> pure (Left e)
    Right (n, _, _) -> withChain chain $ \pc -> fillDoubles n (C.myb_chain_guess pc)

-- | The propellant (kg), the constraints, the gradient of the propellant and the Jacobian of
-- the constraints (by rows) at x.
chainEvaluate :: Chain -> [Double] -> IO (Result (Double, [Double], [Double], [Double]))
chainEvaluate chain x = do
  sized <- chainSize chain
  case sized of
    Left e -> pure (Left e)
    Right (n, eq, ineq) ->
      let m = eq + ineq
       in requireLength n x $ withChain chain $ \pc -> withDoubles x $ \px -> alloca $ \propellant -> allocaArray (max 1 m) $ \c -> allocaArray n $ \g -> allocaArray (max 1 (m * n)) $ \j -> do
            code <- C.myb_chain_evaluate pc px propellant c g j nullPtr
            andThen code (Right <$> ((,,,) <$> (coerce <$> peek propellant) <*> peekDoubles m c <*> peekDoubles n g <*> peekDoubles (m * n) j))

-- | The chain with the least propellant from the first guess (all zeros: that of
-- 'chainGuess'): the unknowns found and what the solver gives of the chain.
chainSolve :: Chain -> Maybe C.Myb_nlp_options -> [Double] -> IO (Result ([Double], C.Myb_chain_result))
chainSolve chain options x = do
  sized <- chainSize chain
  case sized of
    Left e -> pure (Left e)
    Right (n, _, _) ->
      requireLength n x $ withChain chain $ \pc -> withNlpOptions options $ \po -> withDoubles x $ \px -> with C.zero $ \result -> do
        code <- C.myb_chain_solve pc po px nullPtr result
        andThen code (Right <$> ((,) <$> peekDoubles n px <*> peek result))
