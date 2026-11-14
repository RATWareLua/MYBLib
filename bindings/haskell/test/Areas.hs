-- SPDX-License-Identifier: BSD-3-Clause
-- | The upper layer of the rest of the interface against the same calls made through "Myb.C",
-- one area after another: numbers by their bits, structures by their bytes; and failures as
-- Left with the status of the library.
--
--   A  N bodies and the system of bodies       E  attitude
--   B  estimation                               F  observables
--   C  the restricted three-body problem        G  files: numbers, tables, SPK, SP3, OPM, OMM,
--   D  targeting, optimal transfers, thrust,       OCM, CK and SCLK, TLE
--      IOD, relative motion, conjunctions,      H  the array forms
--      uncertainty
module Areas (theAreas) where

import Checks
import Control.Monad (forM)
import qualified Data.ByteString as B
import Data.Coerce (coerce)
import Data.List (isInfixOf, isPrefixOf)
import Data.Word (Word8)
import Foreign.C.String (castCharToCChar, withCString)
import Foreign.C.Types (CChar, CDouble (..), CInt)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Alloc (alloca, allocaBytes)
import Foreign.Marshal.Array (allocaArray, peekArray, withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, nullPtr)
import Foreign.Storable (Storable, peek, poke, sizeOf)
import Myb (Mat3 (..), Quaternion (..), State (..), Vec3 (..))
import qualified Myb
import qualified Myb.C as C
import qualified Myb.C.Ext.Sgp4 as S
import qualified Myb.Tle as Tle
import System.Directory (getTemporaryDirectory)
import System.Environment (lookupEnv)
import System.FilePath ((</>))

theAreas :: Tally -> IO ()
theAreas tally = do
  nBodies tally
  aSystem tally
  estimation tally
  restricted tally
  targeting tally
  attitude tally
  observables tally
  files tally
  arrayForms tally

-- ---- Helpers -------------------------------------------------------------------------------

gm, radius :: Double
gm = 3.986004418e14
radius = 6378137.0

ok :: C.Myb_status -> Bool
ok = (== C.MYB_OK)

v3 :: Vec3 -> [Double]
v3 (Vec3 a b c) = [a, b, c]

vecsBits :: [Vec3] -> [Word]
vecsBits = bitsOf . concatMap v3

-- | The bytes of a structure as its Storable instance writes it (padding zeroed).
bytesOf :: Storable a => a -> IO [Word8]
bytesOf x = with x $ \p -> peekArray (sizeOf x) (castPtr p)

sameBytes :: Storable a => a -> a -> IO Bool
sameBytes a b = (==) <$> bytesOf a <*> bytesOf b

cd :: [Double] -> [CDouble]
cd = coerce

dd :: [CDouble] -> [Double]
dd = coerce

y0 :: State
y0 = State (Vec3 7.0e6 1.0e5 2.0e5) (Vec3 (-100.0) 7.5e3 1.0e3)

leo :: Myb.Elements
leo = Myb.Elements ((radius + 500.0e3) * (1 - 0.01 * 0.01)) 0.01 0.9 0.5 1.0 0.2

-- ---- A: N bodies ---------------------------------------------------------------------------

nBodies :: Tally -> IO ()
nBodies tally = do
  let gms = [1.32712440018e20, 3.986004418e14, 4.9028e12]
      x = [Vec3 0 0 0, Vec3 1.496e11 0 0, Vec3 1.49984e11 0 0]
      v = [Vec3 0 0 0, Vec3 0 29780 0, Vec3 0 30802 0]
  empty <- Myb.newNBody 0 0 []
  expect tally "A: no massive body" (failsWith Myb.Domain empty)
  made <- Myb.newNBody 3 0 gms
  onRight tally "A: a system of point masses" made $ \nb -> do
    upper <- Myb.accelerations nb 0 x
    wrong <- Myb.accelerations nb 0 (take 2 x)
    expect tally "A: positions of another number" (failsWith Myb.Domain wrong)
    (rawOk, rawAccs, rawState, rawInv) <- withArray (cd gms) $ \pg -> alloca $ \out -> do
      code <- C.myb_nbody_create 3 0 pg nullPtr out
      p <- peek out
      a <- withArray x $ \px -> allocaArray 3 $ \pa -> do
        c <- C.myb_nbody_acceleration (castPtr p) 0 (castPtr px) (castPtr pa)
        (,) c <$> peekArray 3 pa
      inv <- withArray x $ \px -> withArray v $ \pv -> alloca $ \pi' -> C.myb_nbody_invariants_of p (castPtr px) (castPtr pv) pi' >> peek pi'
      st <- with C.zero {C.myb_wh_options'h = 86400, C.myb_wh_options'scheme = C.MYB_WH_DKD} $ \po -> alloca $ \wout -> do
        c1 <- C.myb_wh_create p po nullPtr wout
        w <- peek wout
        c2 <- withArray x $ \px -> withArray v $ \pv -> C.myb_wh_start w 0 (castPtr px) (castPtr pv)
        c3 <- C.myb_wh_integrate w (10 * 86400)
        (c4, sx, sv) <- allocaArray 3 $ \px -> allocaArray 3 $ \pv -> do
          c <- C.myb_wh_state w (castPtr px) (castPtr pv)
          (,,) c <$> peekArray 3 px <*> peekArray 3 pv
        C.myb_wh_destroy w
        pure (all ok [c1, c2, c3, c4], sx, sv)
      C.myb_nbody_destroy p
      pure (ok code && ok (fst a), snd a, st, inv)
    expect tally "A: accelerations bit for bit" (either (const False) ((== vecsBits rawAccs) . vecsBits) upper && rawOk)
    inv <- Myb.invariants nb x v
    onRight tally "A: invariants" inv $ \i -> sameBytes i rawInv >>= expect tally "A: invariants byte for byte"
    wh <- Myb.newWh nb (Myb.WhOptions 86400 Myb.Dkd 0)
    onRight tally "A: a Wisdom-Holman integrator" wh $ \w -> do
      _ <- Myb.whStart w 0 x v
      _ <- Myb.whIntegrate w (10 * 86400)
      state <- Myb.whState w
      let (good, sx, sv) = rawState
      expect tally "A: the integrator bit for bit" (good && either (const False) (\(a, b) -> vecsBits a == vecsBits sx && vecsBits b == vecsBits sv) state)
      t <- Myb.whTime w
      expect tally "A: the time of the node" (t == 10 * 86400)
    approximated <- Myb.setApproximation nb 8 0.5
    order <- Myb.approximation nb
    expect tally "A: the fast multipole method set" (approximated == Right () && order == (8, 0.5))
    tooHigh <- Myb.setApproximation nb 99 0.5
    expect tally "A: an order too high" (failsWith Myb.Domain tooHigh)
    _ <- Myb.setApproximation nb 0 0
    self <- Myb.setInteraction nb 1 1 Myb.NoInteraction
    expect tally "A: a body on itself" (failsWith Myb.Domain self)
    _ <- Myb.setRelativity nb (Just C.zero {C.myb_ppn'c = 299792458, C.myb_ppn'beta = 1, C.myb_ppn'gamma = 1})
    newtonian <- Myb.accelerations nb 0 x
    withVelocities <- Myb.accelerationsXV nb 0 x v
    expect tally "A: 1PN needs the velocities" (failsWith Myb.Domain newtonian && either (const False) ((== 3) . length) withVelocities)

-- ---- A: the system of bodies -----------------------------------------------------------------

aSystem :: Tally -> IO ()
aSystem tally = do
  let gms = [3.986004418e14, 4.9028e12]
      x = [0, 0, 0, 3.844e8, 0, 0]
      v = [0, 0, 0, 0, 1022, 0]
  made <- Myb.newSystem gms
  onRight tally "A: a system of bodies" made $ \s -> do
    n <- Myb.systemCoordinates s
    expect tally "A: the coordinates" (n == 6)
    upper <- Myb.systemAcceleration s 0 x v
    (rawA, rawY) <- withArray (cd gms) $ \pg -> alloca $ \out -> do
      _ <- C.myb_system_create 2 pg nullPtr out
      p <- peek out
      a <- withArray (cd x) $ \px -> withArray (cd v) $ \pv -> allocaArray 6 $ \pa -> C.myb_system_acceleration (castPtr p) 0 px pv pa >> (dd <$> peekArray 6 pa)
      y <- alloca $ \iout -> do
        _ <- C.myb_system_integrator p nullPtr nullPtr iout
        i <- peek iout
        _ <- withArray (cd (x ++ v)) $ \py -> C.myb_ias15_start i C.p'myb_system_acceleration (castPtr p) 0 py
        _ <- C.myb_ias15_integrate i 3600
        y <- C.myb_ias15_state i >>= peekArray 12
        C.myb_ias15_destroy i
        pure (dd y)
      C.myb_system_destroy p
      pure (a, y)
    expect tally "A: the system's accelerations bit for bit" (fmap bitsOf upper == Right (bitsOf rawA))
    integrator <- Myb.newIntegrator s Nothing
    onRight tally "A: an integrator of the system" integrator $ \i -> do
      _ <- Myb.startIntegrator i 0 (x ++ v)
      _ <- Myb.integratorIntegrate i 3600
      y <- Myb.integratorState i
      expect tally "A: the system integrated bit for bit" (bitsOf y == bitsOf rawY)
      _ <- Myb.setRadius s 0 radius
      _ <- Myb.setRadius s 1 1737400
      event <- Myb.advance s i 7200 0
      expect tally "A: an advance that arrives" (fmap Myb.eventKind event == Right Myb.Arrived)
    badRotation <- Myb.setRotation s 7 C.zero
    expect tally "A: a rotation of no body" (failsWith Myb.Domain badRotation)
    same <- Myb.addFieldPair s 0 0
    expect tally "A: a field on itself" (failsWith Myb.Domain same)
    pushed <- Myb.newSystem gms
    onRight tally "A: a system with an action" pushed $ \sp -> do
      _ <- Myb.addAction sp 1 (\_ _ -> pure (Right (Vec3 1.0e-6 0 0, Vec3 0 0 0)))
      with' <- Myb.systemAcceleration sp 0 x v
      case (with', upper) of
        (Right a, Right b) -> expect tally "A: the action pushes" (abs ((a !! 3 - b !! 3) - 1.0e-6) < 1.0e-12)
        _ -> expect tally "A: two accelerations" False
    nb <- Myb.systemNBody s
    count <- Myb.nbodyCount nb
    expect tally "A: the point masses of the system" (count == 2)

-- ---- B: estimation -------------------------------------------------------------------------

estimation :: Tally -> IO ()
estimation tally = do
  modelMade <- Myb.newForceModel
  onRight tally "B: a model" modelMade $ \model -> do
    _ <- Myb.addPointMass model "earth" gm
    let times = [60 * fromIntegral k | k <- [1 .. 12 :: Int]]
        truth t = either (const y0) id (Myb.keplerPropagate gm y0 t)
        observations = [Myb.observation t Myb.MeasurePosition (v3 (Myb.position (truth t))) [1, 1, 1] | t <- times]
        guess = State (Vec3 (7.0e6 + 100) 1.0e5 2.0e5) (Vec3 (-100.0) 7.5e3 1.0e3)
    fitted <- Myb.fitOrbit model 0 guess observations Myb.defaultFitOptions
    raw <- withForeignPtr (Myb.forceModelPtr model) $ \m ->
      withArray [C.zero {C.myb_observation't = coerce t, C.myb_observation'measure = C.p'myb_measure_position, C.myb_observation'count = 3, C.myb_observation'value = cd (v3 (Myb.position (truth t))), C.myb_observation'sigma = [1, 1, 1]} | t <- times] $ \po ->
        with guess $ \py -> alloca $ \result -> do
          code <- C.myb_fit_orbit m 0 (castPtr py) po (fromIntegral (length times)) nullPtr nullPtr result nullPtr nullPtr
          (,,) code <$> peek py <*> peek result
    let (code, rawState, rawResult) = raw
    onRight tally "B: a fit" fitted $ \f -> do
      expect tally "B: the fit bit for bit" (ok code && bits (Myb.fitState f) == bits rawState)
      same <- sameBytes (Myb.fitSummary f) rawResult
      expect tally "B: what the fit did, byte for byte" same
      expect tally ("B: the fit near the truth (" ++ show (C.myb_fit_result'converged rawResult, C.myb_fit_result'iterations rawResult, C.myb_fit_result'rms rawResult, Myb.position (Myb.fitState f)) ++ ")") (case Myb.position (Myb.fitState f) of Vec3 a _ _ -> abs (a - 7.0e6) < 1.0e-2)
      expect tally "B: the residuals and the covariance" (length (Myb.fitResiduals f) == 36 && length (Myb.fitCovariance f) == 36)
    none <- Myb.fitOrbit model 0 guess [] Myb.defaultFitOptions
    expect tally "B: no observations" (failsWith Myb.Domain none)
    -- the filter against the same calls
    let covariance = [if i == j then (if i < 3 then 1.0e4 else 1.0) else 0 | i <- [0 .. 5 :: Int], j <- [0 .. 5 :: Int]]
    filt <- Myb.newFilter model Myb.defaultFilterOptions
    onRight tally "B: a filter" filt $ \fl -> do
      _ <- Myb.filterStart fl 0 guess covariance 1000
      innovation <- Myb.filterUpdate fl (Myb.observation 60 Myb.MeasurePosition (v3 (Myb.position (truth 60))) [1, 1, 1])
      (t, state, _) <- Myb.filterGet fl
      rawFilter <- withForeignPtr (Myb.forceModelPtr model) $ \m -> with C.zero $ \po -> alloca $ \fout -> do
        _ <- C.myb_filter_create m po nullPtr fout
        f <- peek fout
        _ <- with guess $ \py -> withArray (cd covariance) $ \pc -> C.myb_filter_start f 0 (castPtr py) pc 1000
        _ <- with C.zero {C.myb_observation't = 60, C.myb_observation'measure = C.p'myb_measure_position, C.myb_observation'count = 3, C.myb_observation'value = cd (v3 (Myb.position (truth 60))), C.myb_observation'sigma = [1, 1, 1]} $ \pob ->
          C.myb_filter_update f pob nullPtr
        s <- alloca $ \py -> C.myb_filter_get f nullPtr (castPtr py) nullPtr >> peek py
        C.myb_filter_destroy f
        pure s
      expect tally "B: the filter bit for bit" (t == 60 && bits state == bits rawFilter && fmap length innovation == Right 3)
  -- a station's measurement against the same call
  envMade <- Myb.newEnvironment (Myb.timeFromSeconds 0)
  onRight tally "B: an environment" envMade $ \env -> do
    _ <- Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 399, Myb.bodyGm = gm, Myb.bodyRadius = radius}
    station <- Myb.newStation env 0 C.zero {C.myb_ellipsoid'a = coerce radius, C.myb_ellipsoid'f = 1 / 298.257223563} 0.5 0.8 100
    onRight tally "B: a station" station $ \st -> do
      measured <- Myb.evaluateMeasure (Myb.MeasureRange st) 1 0 y0
      record <- Myb.stationRecord st
      raw <- with record $ \ps -> with y0 $ \py -> alloca $ \value -> allocaArray 6 $ \partials -> do
        code <- C.myb_measure_range (castPtr ps) 0 (castPtr py) value partials
        (,,) code <$> peek value <*> peekArray 6 partials
      let (code, value, partials) = raw
      expect tally "B: a range bit for bit" (ok code && fmap (\(a, b) -> bitsOf (a ++ b)) measured == Right (bitsOf (dd (value : partials))))
      zero <- Myb.evaluateMeasure (Myb.MeasureRange st) 0 0 y0
      expect tally "B: a measurement of no values" (failsWith Myb.Domain zero)
      haskell <- Myb.evaluateMeasure (Myb.MeasureFunction (\_ y -> pure (Right (take 1 (Myb.stateList y), [1, 0, 0, 0, 0, 0])))) 1 0 y0
      expect tally "B: a measurement of Haskell" (haskell == Right ([7.0e6], [1, 0, 0, 0, 0, 0]))

-- ---- C: the restricted three-body problem ----------------------------------------------------

restricted :: Tally -> IO ()
restricted tally = do
  let mu = 0.012150585609624
  rawPoints <- allocaArray 15 $ \p -> C.myb_cr3bp_lagrange (coerce mu) p >> (dd <$> peekArray 15 p)
  expect tally "C: the Lagrange points bit for bit" (fmap (bitsOf . concatMap v3) (Myb.lagrange mu) == Right (bitsOf rawPoints))
  expect tally "C: mu outside" (failsWith Myb.Domain (Myb.lagrange 0.7))
  onRight tally "C: Richardson's halo" (Myb.richardson mu 1 0.01 1) $ \(guess, _) -> do
    let corrected = Myb.correct mu (Myb.correction Myb.SymmetryXZ 2) guess
    raw <- with C.zero {C.myb_cr3bp_correction'symmetry = C.MYB_CR3BP_SYMMETRY_XZ, C.myb_cr3bp_correction'fixed = 2, C.myb_cr3bp_correction'crossing = 1} $ \pc -> with guess $ \py -> alloca $ \period -> do
      code <- C.myb_cr3bp_correct (coerce mu) pc (castPtr py) period
      (,,) code <$> peek py <*> peek period
    let (code, rawY, rawPeriod) = raw
    onRight tally "C: a halo orbit" corrected $ \(y, period) -> do
      expect tally "C: the correction bit for bit" (ok code && bits y == bits rawY && bitsOf [period] == bitsOf (dd [rawPeriod]))
      rawM <- with y $ \py -> allocaArray 36 $ \pm -> C.myb_cr3bp_monodromy (coerce mu) (castPtr py) (coerce period) pm >> (dd <$> peekArray 36 pm)
      expect tally "C: the monodromy bit for bit" (fmap bitsOf (Myb.monodromy mu y period) == Right (bitsOf rawM))
      onRight tally "C: the stability" (Myb.monodromy mu y period >>= Myb.stability) $ \(_, _, index) -> expect tally "C: an unstable halo" (index > 1)
      rawJ <- with y $ \py -> C.myb_cr3bp_jacobi (coerce mu) (castPtr py)
      expect tally "C: the Jacobi constant bit for bit" (bitsOf [Myb.jacobi mu y] == bitsOf (dd [rawJ]))
      let half = Myb.propagate mu y period (Just (Myb.Section 1 0 0))
      rawHalf <- with y $ \py -> with C.zero {C.myb_cr3bp_section'axis = 1} $ \ps -> alloca $ \out -> alloca $ \el -> alloca $ \cr -> do
        poke cr (0 :: CInt)
        c <- C.myb_cr3bp_propagate (coerce mu) (castPtr py) (coerce period) ps (castPtr out) el cr
        (,) c <$> peek out
      expect tally "C: a section bit for bit" (ok (fst rawHalf) && fmap (\(s, _, crossed) -> (bits s, crossed)) half == Right (bits (snd rawHalf), True))
      onRight tally "C: an orbit cut in patches" (Myb.orbitFromState mu y period 2) $ \orbit -> do
        rawOrbit <- with y $ \py -> with C.zero $ \po -> C.myb_cr3bp_orbit_from_state (coerce mu) (castPtr py) (coerce period) 2 po >> peek po
        sameBytes orbit rawOrbit >>= expect tally "C: the patches byte for byte"
      -- a quasi-halo torus about it: the circle of the linear motion, its points by the
      -- interpolant, its spectrum (with the normal pair of an elliptic torus)
      let torus = Myb.torusFromOrbit mu y period 5 0 1.0e-4
      rawTorus <- with y $ \py -> with C.zero $ \pt -> C.myb_cr3bp_torus_from_orbit (coerce mu) (castPtr py) (coerce period) 5 0 1.0e-4 pt >>= \c -> (,) c <$> peek pt
      onRight tally "C: a torus" torus $ \t -> do
        same <- sameBytes t (snd rawTorus)
        expect tally "C: the torus byte for byte" (ok (fst rawTorus) && same)
        rawPoint <- with t $ \pt -> alloca $ \py -> C.myb_cr3bp_torus_state pt 0.3 (castPtr py) >> peek py
        expect tally "C: a point of the circle bit for bit" (bits (Myb.torusState t 0.3) == bits rawPoint)
        case Myb.torusStability mu t 0 of
          Right spectrum -> do
            rawSpectrum <- with t $ \pt -> with C.zero $ \ps -> C.myb_cr3bp_torus_stability (coerce mu) pt 0 ps >> peek ps
            sameBytes spectrum rawSpectrum >>= expect tally "C: the spectrum of the torus byte for byte"
          Left e -> expect tally ("C: the spectrum of the torus (" ++ Myb.errorText e ++ ")") (Myb.errorStatus e == Myb.NoConvergence)
      expect tally "C: a torus of no size" (failsWith Myb.Domain (Myb.torusFromOrbit mu y period 5 0 0))

-- ---- D: targeting, optimal transfers, thrust, IOD, relative motion, conjunctions ------------

targeting :: Tally -> IO ()
targeting tally = do
  let r1 = Vec3 7.0e6 0 0
      r2 = Vec3 0 8.0e6 1.0e5
  rawLambert <- withArray (cd (v3 r1)) $ \p1 -> withArray (cd (v3 r2)) $ \p2 -> allocaArray 1 $ \ps -> alloca $ \count -> do
    poke count 0
    code <- C.myb_lambert (coerce gm) p1 p2 3000 0 0 ps 1 count
    (,) code <$> peekArray 1 ps
  onRight tally "D: Lambert" (Myb.lambert gm r1 r2 3000 False 0) $ \solutions -> do
    same <- and <$> sequence (zipWith sameBytes solutions (snd rawLambert))
    expect tally "D: Lambert byte for byte" (ok (fst rawLambert) && length solutions == 1 && same)
  expect tally "D: a B-plane of an ellipse" (failsWith Myb.Domain (Myb.bplane gm (Vec3 7.0e6 0 0) (Vec3 0 7.5e3 0) (Vec3 0 0 1)))
  shot <- Myb.shoot (\u _ -> let x = sum u in pure (Right ([x * x - 2], [2 * x]))) 1 [1] (Myb.ShootOptions 1e-14 0)
  expect tally "D: shooting with a Haskell function" (either (const False) (\s -> abs (sum (Myb.shotUnknowns s) - sqrt 2) < 1e-12) shot)
  -- IOD and relative motion
  let positions = [either (const r1) Myb.position (Myb.keplerPropagate gm y0 t) | t <- [0, 600, 1200]]
  case positions of
    [p1, p2, p3] -> do
      rawGibbs <- withArray (cd (v3 p1)) $ \a -> withArray (cd (v3 p2)) $ \b -> withArray (cd (v3 p3)) $ \c -> allocaArray 3 $ \v -> alloca $ \angle -> do
        code <- C.myb_iod_gibbs (coerce gm) a b c v angle
        (,) code . dd <$> peekArray 3 v
      expect tally "D: Gibbs bit for bit" (ok (fst rawGibbs) && fmap (bitsOf . v3 . fst) (Myb.gibbs gm p1 p2 p3) == Right (bitsOf (snd rawGibbs)))
    _ -> expect tally "D: three positions" False
  rawCw <- allocaArray 36 $ \p -> C.myb_cw_stm 0.0011 300 p >> (dd <$> peekArray 36 p)
  expect tally "D: Clohessy-Wiltshire bit for bit" (bitsOf (Myb.cwStm 0.0011 300) == bitsOf rawCw)
  let deputy = State (Vec3 7.0001e6 1.0e5 2.0e5) (Vec3 (-100.0) 7.5001e3 1.0e3)
  rawRel <- with y0 $ \pc -> with deputy $ \pd -> allocaArray 6 $ \o -> do
    code <- C.myb_relative_from_inertial (castPtr pc) (castPtr pd) o
    (,) code . dd <$> peekArray 6 o
  expect tally "D: relative motion bit for bit" (ok (fst rawRel) && fmap bits (Myb.relativeFromInertial y0 deputy) == Right (bitsOf (snd rawRel)))
  -- thrust
  rawEdelbaum <- alloca $ \dv -> alloca $ \beta -> C.myb_edelbaum 7.7e3 3.1e3 0.5 dv beta >> ((\a b -> dd [a, b]) <$> peek dv <*> peek beta)
  expect tally "D: Edelbaum bit for bit" (fmap (\(a, b) -> bitsOf [a, b]) (Myb.edelbaum 7.7e3 3.1e3 0.5) == Right (bitsOf rawEdelbaum))
  expect tally "D: Edelbaum without speed" (failsWith Myb.Domain (Myb.edelbaum 0 3.1e3 0.5))
  let law = C.zero {C.myb_qlaw'mu = coerce gm, C.myb_qlaw'target = [8.0e6, 0.001, 1.0, 0, 0], C.myb_qlaw'weight = [1, 1, 1, 0, 0]}
  rawQ <- with law $ \pl -> with (Myb.elementsToC leo) $ \pe -> allocaArray 3 $ \d -> alloca $ \q -> do
    code <- C.myb_qlaw_direction pl pe d q nullPtr nullPtr
    (,,) code <$> (dd <$> peekArray 3 d) <*> peek q
  let (qcode, qd, qq) = rawQ
  expect tally "D: the Q-law bit for bit" (ok qcode && fmap (\(d, q, _, _) -> bitsOf (v3 d ++ [q])) (Myb.qlawDirection law leo) == Right (bitsOf (qd ++ dd [qq])))
  -- optimal: the flyby of a chain (new), a nonlinear program, Sims and Flanagan
  rawFlyby <- alloca $ \dv -> alloca $ \turn -> allocaArray 6 $ \p -> do
    code <- C.myb_chain_flyby 3.2486e14 5000 5200 7.0e6 dv turn p
    (,) code . dd <$> ((\a b c -> a : b : c) <$> peek dv <*> peek turn <*> peekArray 6 p)
  expect tally "D: a powered flyby bit for bit" (ok (fst rawFlyby) && fmap (\(a, b, c) -> bitsOf (a : b : c)) (Myb.chainFlyby 3.2486e14 5000 5200 7.0e6) == Right (bitsOf (snd rawFlyby)))
  expect tally "D: a flyby without speed" (failsWith Myb.Domain (Myb.chainFlyby 3.2486e14 0 5200 7.0e6))
  let program =
        Myb.NlpProblem 2 1 0 Nothing Nothing $ \x _ -> case x of
          [a, b] -> pure (Right (Myb.NlpValues ((a - 1) ^ (2 :: Int) + (b - 2) ^ (2 :: Int)) [a + b - 1] [2 * (a - 1), 2 * (b - 2)] [1, 1]))
          _ -> pure (Left Myb.Domain)
  solved <- Myb.nlpSolve program Nothing [0.5, 0.5]
  expect tally "D: a nonlinear program" (either (const False) (\(x, _, _) -> case x of [a, b] -> abs a < 1e-5 && abs (b - 1) < 1e-5; _ -> False) solved)
  badProgram <- Myb.nlpSolve program Nothing [0]
  expect tally "D: a first guess of another length" (failsWith Myb.Domain badProgram)
  let target = either (const y0) id (Myb.keplerPropagate gm y0 1500)
      leg = C.zero {C.myb_sims_flanagan'mu = coerce gm, C.myb_sims_flanagan'thrust = 0.1, C.myb_sims_flanagan'exhaust = 30000, C.myb_sims_flanagan'start = cd (Myb.stateList y0), C.myb_sims_flanagan'mass0 = 1000, C.myb_sims_flanagan't1 = 1500, C.myb_sims_flanagan'target = cd (Myb.stateList target), C.myb_sims_flanagan'segments = 2}
  mismatch <- Myb.simsFlanaganMismatch leg Nothing (replicate 6 0)
  rawMismatch <- with leg $ \pl -> withArray (replicate 6 0) $ \pt -> allocaArray 6 $ \m -> alloca $ \mass -> do
    code <- C.myb_sims_flanagan_mismatch pl pt m mass nullPtr
    (,) code . dd <$> peekArray 6 m
  expect tally "D: Sims and Flanagan bit for bit" (ok (fst rawMismatch) && fmap (\(m, _, _) -> bitsOf m) mismatch == Right (bitsOf (snd rawMismatch)))
  -- conjunctions
  let rel = Vec3 100 50 0
      relV = Vec3 0 0 10
      cov = Mat3 (Vec3 400 0 0) (Vec3 0 900 0) (Vec3 0 0 100)
  rawEnc <- withArray (cd (v3 rel)) $ \pr -> withArray (cd (v3 relV)) $ \pv -> withArray (cd (Myb.mat3List cov)) $ \pc -> with C.zero $ \pe -> do
    code <- C.myb_encounter_plane pr pv pc 20 pe
    e <- peek pe
    pr' <- alloca $ \p -> C.myb_collision_probability pe p >> peek p
    pure (code, e, pr')
  let (ecode, rawE, rawP) = rawEnc
  onRight tally "D: an encounter" (Myb.encounterPlane rel relV cov 20) $ \e -> do
    same <- sameBytes e rawE
    expect tally "D: the encounter byte for byte" (ok ecode && same)
    expect tally "D: the probability bit for bit" (fmap (\p -> bitsOf [p]) (Myb.collisionProbability e) == Right (bitsOf (dd [rawP])))
  expect tally "D: an encounter without relative velocity" (failsWith Myb.Domain (Myb.encounterPlane rel (Vec3 0 0 0) cov 20))
  approaches <- Myb.approachFind (Myb.fromFunction (\_ -> pure (Right (State (Vec3 0 0 0) (Vec3 0 0 0))))) (Myb.fromFunction (\t -> pure (Right (State (Vec3 (-1000 + 10 * t) 1000 0) (Vec3 10 0 0))))) 0 300 10 4
  expect tally "D: the closest approach of two Haskell providers" (either (const False) (\(as, n) -> n == 1 && length as == 1 && all (\a -> abs (coerce (C.myb_approach't a) - (100 :: Double)) < 1e-6) as) approaches)
  -- uncertainty
  rng <- Myb.newRng 42
  draws <- forM [1 .. 3 :: Int] (const (Myb.rngUniform rng))
  rawDraws <- with C.zero $ \p -> C.myb_rng_seed p 42 >> forM [1 .. 3 :: Int] (const (C.myb_rng_uniform p))
  expect tally "D: random numbers bit for bit" (bitsOf draws == bitsOf (dd rawDraws))
  let samples = [1, 2, 2, 3, 4, 1, 0, 5]
  rawStats <- withArray (cd samples) $ \ps -> allocaArray 2 $ \pm -> allocaArray 4 $ \pc -> C.myb_sample_statistics 4 2 ps pm pc >> ((++) <$> (dd <$> peekArray 2 pm) <*> (dd <$> peekArray 4 pc))
  expect tally "D: sample statistics bit for bit" (fmap (\(m, c) -> bitsOf (m ++ c)) (Myb.sampleStatistics 2 samples) == Right (bitsOf rawStats))
  expect tally "D: one sample" (failsWith Myb.Domain (Myb.sampleStatistics 2 [1, 2]))

-- ---- E: attitude ---------------------------------------------------------------------------

attitude :: Tally -> IO ()
attitude tally = do
  let qa = Quaternion 0.9 0.1 (-0.3) 0.2
      qb = Quaternion 0.5 0.5 0.5 0.5
  rawProduct <- with qa $ \pa -> with qb $ \pb -> alloca $ \o -> C.myb_quat_multiply (castPtr pa) (castPtr pb) (castPtr o) >> peek o
  expect tally "E: a product of quaternions bit for bit" (qBits (Myb.quatMultiply qa qb) == qBits rawProduct)
  rawMatrix <- with qa $ \pa -> alloca $ \o -> C.myb_quat_to_matrix (castPtr pa) (castPtr o) >> peek o
  expect tally "E: a rotation matrix bit for bit" (bitsOf (Myb.mat3List (Myb.quatToMatrix qa)) == bitsOf (Myb.mat3List rawMatrix))
  rawEuler <- withArray [3, 1, 3 :: CInt] $ \pa -> withArray [0.1, 0.2, 0.3] $ \pg -> alloca $ \o -> C.myb_euler_to_matrix pa pg (castPtr o) >>= \c -> (,) c <$> peek o
  expect tally "E: Euler angles bit for bit" (ok (fst rawEuler) && fmap (bitsOf . Myb.mat3List) (Myb.eulerToMatrix (3, 1, 3) (Vec3 0.1 0.2 0.3)) == Right (bitsOf (Myb.mat3List (snd rawEuler))))
  expect tally "E: an invalid sequence" (failsWith Myb.Domain (Myb.eulerToMatrix (1, 1, 2) (Vec3 0.1 0.2 0.3)))
  let moments = Vec3 1 2 3
      l = Vec3 0.1 0.2 1.5
  onRight tally "E: a free rotation" (Myb.freeRotationInit moments qa l) $ \rotation -> do
    raw <- with moments $ \pm -> with qa $ \pq -> with l $ \pl -> with C.zero $ \pr -> do
      _ <- C.myb_free_rotation_init (castPtr pm) (castPtr pq) (castPtr pl) pr
      alloca $ \q -> alloca $ \ll -> C.myb_free_rotation_state pr 10 (castPtr q) (castPtr ll) >> ((,) <$> peek q <*> peek ll)
    let (q, ll) = Myb.freeRotationState rotation 10
    expect tally "E: the free rotation bit for bit" (qBits q == qBits (fst raw) && bitsOf (v3 ll) == bitsOf (v3 (snd raw)))
  split <- Myb.rotationSplit moments (\_ _ -> pure (Right (Vec3 0 0 0))) 4 0 0.1 100 qa l
  rawSplit <- do
    f <- C.mk'myb_torque_function (\_ _ _ o -> poke (castPtr o) (Vec3 0 0 0) >> pure C.MYB_OK)
    r <- with moments $ \pm -> with qa $ \pq -> with l $ \pl -> do
      c <- C.myb_rotation_split (castPtr pm) f nullPtr 4 0 0.1 100 (castPtr pq) (castPtr pl)
      (,,) c <$> peek pq <*> peek pl
    pure r
  let (scode, sq, sl) = rawSplit
  expect tally "E: the splitting integrator bit for bit" (ok scode && fmap (\(q, ll) -> (qBits q, bitsOf (v3 ll))) split == Right (qBits sq, bitsOf (v3 sl)))
  where
    qBits (Quaternion a b c d) = bitsOf [a, b, c, d]

-- ---- F: observables ------------------------------------------------------------------------

observables :: Tally -> IO ()
observables tally = do
  let p = Vec3 0.6 0.8 0
      vel = Vec3 3.0e4 0 0
  rawAb <- withArray (cd (v3 p)) $ \pp -> withArray (cd (v3 vel)) $ \pv -> allocaArray 3 $ \o -> C.myb_aberration pp pv 299792458 o >> (dd <$> peekArray 3 o)
  expect tally "F: aberration bit for bit" (bitsOf (v3 (Myb.aberration p vel 299792458)) == bitsOf rawAb)
  rawRadec <- withArray (cd (v3 p)) $ \pp -> alloca $ \ra -> alloca $ \de -> C.myb_direction_radec pp ra de >> ((\a b -> dd [a, b]) <$> peek ra <*> peek de)
  expect tally "F: a direction bit for bit" ((\(a, b) -> bitsOf [a, b]) (Myb.directionRadec p) == bitsOf rawRadec)
  expect tally "F: an invalid profile" (failsWith Myb.Domain (Myb.plasmaColumn C.zero (Vec3 1 0 0) (Vec3 0 1 0)))
  tests <- testDataDir
  loaded <- Myb.newEphemeris
  onRight tally "F: an SPK set" loaded $ \eph -> do
    _ <- Myb.loadEphemerisFile eph (tests </> "spk" </> "t02_03.bsp")
    segs <- Myb.segments eph
    case segs of
      (first : _) -> do
        let middle = Myb.timeFromSeconds (0.5 * (Myb.segmentStart first + Myb.segmentStop first))
            model = Myb.LightModel 299792458 1 []
            source = Myb.fromEphemeris eph (Myb.segmentTarget first) (Myb.segmentCenter first) 0
            observer = State (Vec3 1.0e7 0 0) (Vec3 0 0 0)
        path <- Myb.lightTime model source middle observer Myb.Reception
        raw <- withForeignPtr (Myb.ephemerisPtr eph) $ \spk ->
          with C.zero {C.myb_spk_body'spk = spk, C.myb_spk_body'target = Myb.segmentTarget first, C.myb_spk_body'center = Myb.segmentCenter first} $ \body ->
            with C.zero {C.myb_light_model'c = 299792458, C.myb_light_model'gamma = 1} $ \pm -> with observer $ \po -> with C.zero $ \pp -> do
              code <- C.myb_light_time pm C.p'myb_spk_body_state (castPtr body) (Myb.timeToC middle) (castPtr po) C.MYB_LIGHT_RECEPTION pp
              (,) code <$> peek pp
        onRight tally "F: a light time" path $ \lp -> do
          same <- sameBytes lp (snd raw)
          expect tally "F: the light time byte for byte" (ok (fst raw) && same)
      [] -> expect tally "F: the segments" False

-- ---- G: files ------------------------------------------------------------------------------

cchars :: String -> [CChar]
cchars = map castCharToCChar

files :: Tally -> IO ()
files tally = do
  expect tally "G: a number read" (Myb.parseDouble "1.25e3x" == Right (1250, 6))
  expect tally "G: no number" (failsWith Myb.Format (Myb.parseDouble "abc"))
  expect tally "G: the shortest decimal" (Myb.formatDouble 0.1 == "0.1" && Myb.formatFixed 2.5 1 == Right "2.5")
  -- a table
  written <- Myb.writeTable (Just ["a", "b"]) 2 [1, 2, 3, 4]
  onRight tally "G: a table written" written $ \text -> do
    table <- Myb.loadTable text
    onRight tally "G: the table read" table $ \t -> do
      rows <- Myb.tableRows t
      row <- Myb.tableRow t 1
      name <- Myb.tableName t 1
      expect tally "G: the table again" (rows == 2 && row == Just [3, 4] && name == Just "b")
  broken <- Myb.loadTable "1,2\n3\n"
  expect tally "G: rows of other lengths" (either (\e -> Myb.errorStatus e == Myb.Format && "line 2" `isPrefixOf` Myb.errorMessage e) (const False) broken)
  -- an SPK file written, against the same calls, and read back
  let epochs = [60 * fromIntegral k | k <- [0 .. 10 :: Int]]
      states = [either (const y0) id (Myb.keplerPropagate gm y0 t) | t <- epochs]
      header = Myb.SegmentHeader (-1000) 399 1 "TEST"
  writer <- Myb.newSpkWriter Nothing
  onRight tally "G: an SPK writer" writer $ \w -> do
    added <- Myb.spkAddStates w header 13 7 epochs states
    upper <- Myb.spkWriterBytes w
    raw <- alloca $ \out -> do
      _ <- C.myb_spk_writer_create nullPtr nullPtr out
      pw <- peek out
      _ <- withCString "TEST" $ \name -> with C.zero {C.myb_spk_segment_header'target = -1000, C.myb_spk_segment_header'center = 399, C.myb_spk_segment_header'frame = 1, C.myb_spk_segment_header'name = name} $ \ph ->
        withArray (cd epochs) $ \pe -> withArray states $ \ps -> C.myb_spk_writer_add_states pw ph 13 7 pe (castPtr ps) 11
      size <- C.myb_spk_writer_size pw
      bytes <- allocaBytes (fromIntegral size) $ \buffer -> C.myb_spk_writer_write pw buffer size >> B.packCStringLen (castPtr buffer, fromIntegral size)
      C.myb_spk_writer_destroy pw
      pure bytes
    expect tally "G: the SPK file byte for byte" (added == Right () && upper == Right raw)
    directory <- maybe getTemporaryDirectory pure =<< lookupEnv "MYB_TEST_OUTPUT_DIR"
    let path = directory </> "myb-haskell-test.bsp"
    _ <- Myb.spkWriterFile w path
    readBack <- Myb.newEphemeris
    onRight tally "G: the SPK file read back" readBack $ \eph -> do
      _ <- Myb.loadEphemerisFile eph path
      at <- Myb.ephemerisState eph (-1000) 399 (Myb.timeFromSeconds 120)
      let near a b = abs (a - b) <= 1e-9 * max 1 (abs b)
      expect tally "G: a state read back" (either (const False) (\s -> and (zipWith near (Myb.stateList s) (Myb.stateList (states !! 2)))) at)
  badHeader <- Myb.newSpkWriter Nothing >>= either (pure . Left) (\w -> Myb.spkAddStates w header {Myb.headerName = ""} 13 7 epochs states)
  expect tally "G: a segment without a name" (failsWith Myb.Domain badHeader)
  -- SP3 written, against the same call, and read back
  let sp3 = Myb.Sp3Data 'c' 'L' "UTC" "SLR20" "FIT" "TEST" "SLR" ["L51"] [Myb.Calendar 2026 3 1 0 0 0, Myb.Calendar 2026 3 1 0 2 0] 120 [Vec3 7.0e6 1.0e5 2.0e5, Vec3 7.1e6 1.0e5 2.0e5] Nothing Nothing Nothing Nothing Nothing Nothing
  text <- Myb.writeSp3 sp3
  rawSp3 <- withCString "UTC" $ \ts -> withCString "SLR20" $ \fr -> withCString "FIT" $ \ot -> withCString "TEST" $ \ag -> withCString "SLR" $ \du -> withCString "L51" $ \s1 -> withArray [s1] $ \sats ->
    withArray [Myb.calendarToC (Myb.Calendar 2026 3 1 0 0 0), Myb.calendarToC (Myb.Calendar 2026 3 1 0 2 0)] $ \pe -> withArray [7.0e6, 1.0e5, 2.0e5, 7.1e6, 1.0e5, 2.0e5] $ \pp ->
      with C.zero {C.myb_sp3_data'version = castCharToCChar 'c', C.myb_sp3_data'file_type = castCharToCChar 'L', C.myb_sp3_data'time_system = ts, C.myb_sp3_data'frame = fr, C.myb_sp3_data'orbit_type = ot, C.myb_sp3_data'agency = ag, C.myb_sp3_data'data_used = du, C.myb_sp3_data'satellite_count = 1, C.myb_sp3_data'satellites = sats, C.myb_sp3_data'epoch_count = 2, C.myb_sp3_data'epochs = pe, C.myb_sp3_data'interval = 120, C.myb_sp3_data'positions = pp} $ \pd ->
        alloca $ \len -> do
          poke len 0
          _ <- C.myb_sp3_write pd nullPtr 0 len
          n <- peek len
          allocaBytes (fromIntegral n + 1) $ \buffer -> C.myb_sp3_write pd buffer (n + 1) len >> map (toEnum . fromIntegral) <$> (peekArray (fromIntegral n) (castPtr buffer) :: IO [Word8])
  onRight tally "G: an SP3 file written" text $ \t -> do
    expect tally "G: the SP3 file character for character" (t == rawSp3)
    parsed <- Myb.loadSp3 t
    onRight tally "G: the SP3 file read" parsed $ \s -> do
      sats <- Myb.sp3Satellites s
      st <- Myb.sp3State s 0 1
      expect tally "G: a position read back" (sats == ["L51"] && fmap (fmap v3 . fst) st == Just (Just [7.1e6, 1.0e5, 2.0e5]))
  -- OPM through the writer and the parser
  let opm =
        C.zero
          { C.myb_opm'originator = cchars "MYB"
          , C.myb_opm'object_name = cchars "SAT"
          , C.myb_opm'object_id = cchars "2026-001A"
          , C.myb_opm'center_name = cchars "EARTH"
          , C.myb_opm'ref_frame = cchars "EME2000"
          , C.myb_opm'time_system = cchars "UTC"
          , C.myb_opm'creation_date = Myb.calendarToC (Myb.Calendar 2026 3 2 0 0 0)
          , C.myb_opm'epoch = Myb.calendarToC (Myb.Calendar 2026 3 1 12 0 0)
          , C.myb_opm'state = cd (Myb.stateList y0)
          }
  opmWritten <- Myb.opmText opm Myb.Kvn
  rawOpm <- with opm $ \po -> alloca $ \len -> do
    poke len 0
    _ <- C.myb_opm_write po nullPtr 0 len
    n <- peek len
    allocaBytes (fromIntegral n + 1) $ \buffer -> C.myb_opm_write po buffer (n + 1) len >> map (toEnum . fromIntegral) <$> (peekArray (fromIntegral n) (castPtr buffer) :: IO [Word8])
  onRight tally "G: an OPM written" opmWritten $ \t -> do
    expect tally "G: the OPM character for character" (t == rawOpm)
    back <- Myb.parseOpm t
    expect tally "G: the state of the OPM bit for bit" (fmap (bitsOf . dd . C.myb_opm'state) back == Right (bits y0))
  garbage <- Myb.parseOpm "NOT AN OPM\n"
  expect tally "G: not an OPM" (failsWith Myb.Format garbage)
  -- an OCM written and read
  let ocm =
        Myb.OcmData
          "MYB"
          (Myb.Calendar 2026 3 2 0 0 0)
          [("OBJECT_NAME", "SAT"), ("TIME_SYSTEM", "UTC"), ("EPOCH_TZERO", "2026-03-01T00:00:00")]
          [Myb.OcmBlock Myb.OcmTraj [("CENTER_NAME", "EARTH"), ("TRAJ_REF_FRAME", "EME2000"), ("TRAJ_TYPE", "CARTPV")] ["0.0 7000.0 100.0 200.0 -0.1 7.5 1.0"]]
  ocmText <- Myb.writeOcm ocm Myb.Kvn
  onRight tally "G: an OCM written" ocmText $ \t -> do
    loadedOcm <- Myb.loadOcm t
    onRight tally "G: the OCM read" loadedOcm $ \o -> do
      blocks <- Myb.ocmBlockCount o Myb.OcmTraj
      st <- Myb.ocmState o 0 0
      name <- Myb.ocmMeta o "OBJECT_NAME"
      expect tally "G: the OCM again" (blocks == 1 && name == Just "SAT" && fmap snd st == Right [7.0e6, 1.0e5, 2.0e5, -100, 7.5e3, 1.0e3])
  -- an OEM with a covariance
  let segment = Myb.OemSegment "SAT" "2026-001A" "EARTH" "EME2000" "UTC" Nothing [Myb.Calendar 2026 3 1 12 0 0] [y0] [(Myb.Calendar 2026 3 1 12 0 0, [if i == j then 2 else 0 | i <- [0 .. 5 :: Int], j <- [0 .. 5 :: Int]])] Nothing
  oem <- Myb.writeOem "MYB" (Myb.Calendar 2026 3 2 0 0 0) [segment] Myb.Kvn
  onRight tally "G: an OEM with a covariance" oem $ \t -> do
    parsedOem <- Myb.parseOem t
    onRight tally "G: the OEM read" parsedOem $ \o -> do
      c <- Myb.oemCovarianceCount o 0
      first <- Myb.oemCovariance o 0 0
      version <- Myb.oemHeader o "CCSDS_OEM_VERS"
      expect tally "G: the covariance read back" (c == 1 && fmap (\(_, _, m) -> take 1 m) first == Just [2] && version == Just "2.0")
  -- CK and SCLK
  tests <- testDataDir
  ck <- Myb.newCk
  onRight tally "G: a CK set" ck $ \set -> do
    loadedCk <- Myb.loadCkFile set (tests </> "ck" </> "test.bc")
    pointing <- Myb.ckPointing set (-9099001) 6360436.339305353 0 False
    raw <- withForeignPtr (Myb.ckPtr set) $ \p -> with C.zero $ \pp -> C.myb_ck_pointing_at p (-9099001) 6360436.339305353 0 0 pp >>= \code -> (,) code <$> peek pp
    onRight tally "G: a pointing" (loadedCk >> pointing) $ \pt -> do
      same <- sameBytes pt (snd raw)
      expect tally "G: the pointing byte for byte" (ok (fst raw) && same)
    far <- Myb.ckPointing set (-9099001) 1.0e15 0 False
    expect tally "G: no pointing there" (failsWith Myb.Range far)
  pool <- Myb.newKernelPool
  onRight tally "G: a kernel pool" pool $ \kp -> do
    _ <- Myb.loadKernelFile kp (tests </> "ck" </> "test.tsc")
    described <- Myb.describeKernel kp "SCLK01_MODULI_9099"
    numbers <- Myb.kernelNumbers kp "SCLK01_MODULI_9099"
    expect tally "G: the variables of the pool" (described == Right (Myb.KernelNumbers, 3) && numbers == Right [10000000, 256, 8])
    clock <- Myb.newSclk kp 9099
    onRight tally "G: a clock" clock $ \s -> do
      t <- Myb.sclkToTime s 0
      raw <- withForeignPtr (Myb.sclkPtr s) $ \p -> alloca $ \pt -> C.myb_sclk_to_time p 0 pt >> peek pt
      expect tally "G: ticks to time bit for bit" (fmap Myb.timeToC t == Right raw && fmap Myb.seconds t == Right 599996400.0)
      decoded <- Myb.sclkDecode s 0
      expect tally "G: a clock string" (decoded == Right "1/0000000:000:0")
    missing <- Myb.newSclk kp 1234
    expect tally "G: a clock that is not there" (failsWith Myb.NotFound missing)
  ckWriter <- Myb.newCkWriter Nothing
  onRight tally "G: a CK writer" ckWriter $ \w -> do
    added <- Myb.ckAddPointings w (Myb.CkHeader (-9099001) 1 "TEST") 1 [0, 100] [Quaternion 1 0 0 0, Quaternion 0.6 0.8 0 0] Nothing Nothing
    bytes <- Myb.ckWriterBytes w
    raw <- alloca $ \out -> do
      _ <- C.myb_ck_writer_create nullPtr nullPtr out
      pw <- peek out
      _ <- withCString "TEST" $ \name -> with C.zero {C.myb_ck_segment_header'instrument = -9099001, C.myb_ck_segment_header'frame = 1, C.myb_ck_segment_header'name = name} $ \ph ->
        withArray [0, 100] $ \pt -> withArray [1, 0, 0, 0, 0.6, 0.8, 0, 0] $ \pq -> C.myb_ck_writer_add_pointings pw ph 1 pt pq nullPtr 2 nullPtr 0
      size <- C.myb_ck_writer_size pw
      b <- allocaBytes (fromIntegral size) $ \buffer -> C.myb_ck_writer_write pw buffer size >> B.packCStringLen (castPtr buffer, fromIntegral size)
      C.myb_ck_writer_destroy pw
      pure b
    expect tally "G: the CK file byte for byte" (added == Right () && bytes == Right raw)
    readBack <- Myb.newCk
    case (readBack, bytes) of
      (Right set, Right b) -> do
        _ <- Myb.loadCkBytes set b
        segs <- Myb.ckSegments set
        expect tally "G: the CK file read back" (length segs == 1)
      _ -> expect tally "G: a CK set from memory" False
  -- two-line elements and SGP4 (myb:tle)
  let line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
      line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
  onRight tally "G: a TLE" (Tle.tleParse line1 line2 False) $ \tle -> do
    onRight tally "G: SGP4" (Tle.sgp4Init tle Tle.Wgs72 Tle.Improved) $ \sgp4 -> do
      raw <- with tle $ \pt -> with C.zero $ \ps -> do
        _ <- S.myb_sgp4_init pt S.MYB_SGP4_WGS72 S.MYB_SGP4_IMPROVED ps
        allocaArray 3 $ \r -> allocaArray 3 $ \v -> S.myb_sgp4_propagate ps 100 r v nullPtr >>= \c -> (,) c . dd <$> ((++) <$> peekArray 3 r <*> peekArray 3 v)
      expect tally "G: SGP4 bit for bit" (ok (fst raw) && fmap (\(r, v) -> bitsOf (v3 r ++ v3 v)) (Tle.sgp4Propagate sgp4 100) == Right (bitsOf (snd raw)))
    expect tally "G: the two lines again" (fmap fst (Tle.tleFormat tle) == Right line1)
    onRight tally "G: an OMM of the TLE" (Tle.ommFromTle tle) $ \omm -> do
      ommWritten <- Myb.ommText omm Myb.Kvn
      back <- either (pure . Left) Myb.parseOmm ommWritten
      expect tally "G: the TLE through an OMM" (fmap S.myb_tle'catalog_number (back >>= Tle.tleFromOmm) == Right 5)
  expect tally "G: not a TLE" (failsWith Myb.Format (Tle.tleParse "1 x" "2 y" False))

-- ---- H: the array forms --------------------------------------------------------------------

arrayForms :: Tally -> IO ()
arrayForms tally = do
  let xs = [0.5, 1.0e9 + 0.25, -3.75]
  expect tally "H: instants from seconds as single calls" (Myb.timeFromSecondsMany xs == map Myb.timeFromSeconds xs)
  expect tally "H: seconds of instants as single calls" (bitsOf (Myb.secondsMany (map Myb.timeFromSeconds xs)) == bitsOf (map (Myb.seconds . Myb.timeFromSeconds) xs))
  let label = Myb.Calendar 2026 3 1 12 30 15.25
  rawT <- with (Myb.calendarToC label) $ \pl -> alloca $ \pt -> C.myb_calendar_to_time pl pt >>= \c -> (,) c <$> peek pt
  expect tally "H: calendar labels as single calls" (ok (fst rawT) && fmap (map Myb.timeToC) (Myb.calendarToTimeMany [label]) == Right [snd rawT])
  expect tally "H: instants to labels and back" (Myb.timeToCalendarMany (either (const []) id (Myb.calendarToTimeMany [label])) == [label])
  let dts = [60, 120, 3600]
  expect tally "H: Kepler as single calls" (fmap (map bits) (Myb.keplerPropagateMany gm y0 dts) == mapM (fmap bits . Myb.keplerPropagate gm y0) dts)
  stms <- pure (Myb.keplerStmMany gm y0 dts)
  expect tally "H: transition matrices" (either (const False) (all ((== 36) . length . snd)) stms)
  let bad = leo {Myb.semiLatusRectum = -1}
  expect tally "H: the element that failed" (either (\e -> Myb.errorStatus e == Myb.Domain && "element 1" `isInfixOf` Myb.errorMessage e) (const False) (Myb.stateFromElementsMany gm [leo, bad]))
  expect tally "H: elements as single calls" (fmap (map Myb.elementsToC) (Myb.stateFromElementsMany gm [leo] >>= Myb.elementsFromStateMany gm) == fmap (map Myb.elementsToC) (Myb.stateFromElements gm leo >>= \s -> (: []) <$> Myb.elementsFromState gm s))
  tests <- testDataDir
  loaded <- Myb.newEphemeris
  onRight tally "H: an SPK set" loaded $ \eph -> do
    _ <- Myb.loadEphemerisFile eph (tests </> "spk" </> "t02_03.bsp")
    segs <- Myb.segments eph
    case segs of
      (first : _) -> do
        let ts = [Myb.timeFromSeconds (Myb.segmentStart first + f * (Myb.segmentStop first - Myb.segmentStart first)) | f <- [0.25, 0.5, 0.75]]
            target = Myb.segmentTarget first
            center = Myb.segmentCenter first
        many <- Myb.ephemerisStateMany eph target center ts
        single <- mapM (Myb.ephemerisState eph target center) ts
        expect tally "H: ephemerides as single calls" (fmap (map (bits . fst)) many == fmap (map bits) (sequence single))
      [] -> expect tally "H: the segments" False
  let rotation = C.zero {C.myb_rotation'm = [0, 1, 0, -1, 0, 0, 0, 0, 1]}
  rawRot <- with rotation $ \pr -> withArray [1, 2, 3] $ \pin -> allocaArray 3 $ \o -> C.myb_rotation_apply_position pr pin o >> (dd <$> peekArray 3 o)
  expect tally "H: rotations as single calls" (fmap (concatMap v3) (Myb.rotationApplyPositionMany [rotation] [Vec3 1 2 3]) == Right rawRot)
  expect tally "H: rotations of another number" (failsWith Myb.Domain (Myb.rotationApplyPositionMany [rotation] []))
