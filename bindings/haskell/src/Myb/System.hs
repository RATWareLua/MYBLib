-- SPDX-License-Identifier: BSD-3-Clause
-- | A system of bodies with their models (myb/system.h): the point masses of "Myb.NBody" joined
-- by figures, tides, distortions, pairs beyond 1PN, actions, frame dragging, guided bodies and
-- contacts, with the rotation of chosen bodies integrated together with the motion; an IAS15
-- integrator of the system, its encounters, mergers and bounces, its partial derivatives.
--
-- A state of the system is 'systemCoordinates' coordinates and then their rates: a list of
-- 2 n doubles (3 per body, 9 per integrated rotation, 6 per pair with spins).
module Myb.System
  ( System
  , newSystem
  , systemBodies
  , systemNBody
  , systemKeep
    -- * The model of a body
  , setOrientation
  , setField
  , setRotation
  , setTide
  , C.Myb_system_rotation (..)
  , C.Myb_system_tide (..)
    -- * The levels of the graph
  , addFieldPair
  , addFigures
  , addTide
  , addDistortion
  , addPnPair
  , C.Myb_system_pn_pair (..)
  , PnOrder (..)
  , pnOrders
  , SystemView (..)
  , addAction
  , addLightAction
  , C.Myb_system_light (..)
  , addFrameDragging
  , setMassRate
  , addGuided
  , systemGuided
  , addContact
    -- * The state
  , systemCoordinates
  , packRotation
  , unpackRotation
  , packSpins
  , unpackSpins
  , packContact
  , unpackContact
  , systemAcceleration
  , systemInertia
  , SystemPartials (..)
  , systemPartials
    -- * The integrator
  , SystemIntegrator
  , newIntegrator
  , startIntegrator
  , integratorStep
  , integratorIntegrate
  , integratorTime
  , integratorState
  , integratorDense
  , integratorStatistics
  , C.Myb_ias15_options (..)
  , C.Myb_ias15_stats (..)
    -- * Encounters
  , setRadius
  , systemRadius
  , systemPresent
  , EventKind (..)
  , SystemEvent (..)
  , advance
  , merge
  , C.Myb_system_merger (..)
  , bounce
    -- * The objects of C
  , systemPtr
  , integratorPtr
  ) where

import Data.Bits ((.|.))
import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CSize, CUInt)
import Foreign.ForeignPtr (ForeignPtr, mallocForeignPtrArray, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, pokeArray, withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Gravity (GravityField, gravityFieldPtr)
import Myb.Internal
import Myb.NBody (NBody, wrapNBody)
import Myb.Types

-- | A system and what it keeps.
data System = System (ForeignPtr C.Myb_system) Kept

-- | Massive bodies with the GM given (m^3 s^-2, one each).
newSystem :: [Double] -> IO (Result System)
newSystem gm =
  withDoublesLen gm $ \n g ->
    fmap (uncurry System) <$> newOwned (\o -> check <$> C.myb_system_create (toSize n) g nullPtr o) C.myb_system_destroy

systemPtr :: System -> ForeignPtr C.Myb_system
systemPtr (System p _) = p

-- | Ties the life of an object to the system: what a provider or a level added through
-- "Myb.C" points to (a mesh, the tides of a station).
systemKeep :: System -> ForeignPtr a -> IO ()
systemKeep (System _ kept) = keepAlive kept

systemBodies :: System -> IO Int
systemBodies (System p _) = fromSize <$> withForeignPtr p C.myb_system_bodies

-- | The point masses inside: their graph, relativity and threads are set there.
systemNBody :: System -> IO NBody
systemNBody (System p kept) = do
  raw <- withForeignPtr p C.myb_system_nbody
  inner <- dependent p raw
  pure (wrapNBody inner kept)

call :: System -> (Ptr C.Myb_system -> IO C.Myb_status) -> IO (Result ())
call (System p _) f = withForeignPtr p $ \s -> check <$> f s

-- | An index for C: a negative one becomes a size no system has, which the library refuses.
index :: Int -> CSize
index = toSize

-- | The body-fixed frame of a body from outside. The system keeps the provider.
setOrientation :: System -> Int -> OrientationSource -> IO (Result ())
setOrientation system@(System _ kept) body source = do
  (function, context) <- keepOrientationSource kept source
  call system $ \s -> C.myb_system_set_orientation s (index body) function context

-- | The figure of a body as harmonics in its frame up to a degree (0: all of the field). The
-- system keeps the field.
setField :: System -> Int -> GravityField -> Int -> IO (Result ())
setField system@(System _ kept) body field degree = do
  keepAlive kept (gravityFieldPtr field)
  call system $ \s -> withForeignPtr (gravityFieldPtr field) $ \f -> C.myb_system_set_field s (index body) f (fromIntegral degree)

-- | An integrated rotation of a body.
setRotation :: System -> Int -> C.Myb_system_rotation -> IO (Result ())
setRotation system body rotation = call system $ \s -> with rotation (C.myb_system_set_rotation s (index body))

-- | The answer of a body to tides.
setTide :: System -> Int -> C.Myb_system_tide -> IO (Result ())
setTide system body tide = call system $ \s -> with tide (C.myb_system_set_tide s (index body))

-- | The figure of the body acts on the point.
addFieldPair :: System -> Int -> Int -> IO (Result ())
addFieldPair system body point = call system $ \s -> C.myb_system_add_field s (index body) (index point)

-- | The quadrupoles of two bodies on each other.
addFigures :: System -> Int -> Int -> IO (Result ())
addFigures system a b = call system $ \s -> C.myb_system_add_figures s (index a) (index b)

-- | The tide the raiser raises on the body acts on the point.
addTide :: System -> Int -> Int -> Int -> IO (Result ())
addTide system body raiser point = call system $ \s -> C.myb_system_add_tide s (index body) (index raiser) (index point)

-- | The raiser distorts the body, whose rotation is integrated.
addDistortion :: System -> Int -> Int -> IO (Result ())
addDistortion system body raiser = call system $ \s -> C.myb_system_add_distortion s (index body) (index raiser)

-- | The terms of a PN pair (MYB_PN_* of myb/rel.h), by name.
data PnOrder
  = Pn1
  | Pn2
  | Pn25
  | Pn3
  | Pn35
  | PnSo
  | PnSs
  | Pn4
  | Pn4rr
  | PnSo25
  | PnSo35
  | PnSorr
  | -- | spin-spin at next-to-leading order (1/c^6), with its precession
    PnSs3
  | -- | spin-orbit radiation reaction at next-to-leading order (1/c^10)
    PnSorr5
  | -- | the spin-orbit part of the dissipative tail, quasi-circular orbits
    PnSotail
  | -- | the tails of Pn4 and Pn4rr on any osculating conic
    PnTailx
  deriving (Eq, Show, Enum, Bounded)

-- | The field `orders` of a PN pair: @zero {myb_system_pn_pair'orders = pnOrders [Pn2, Pn25]}@.
pnOrders :: [PnOrder] -> CUInt
pnOrders = foldr ((.|.) . bit') 0
  where
    bit' o = case o of
      Pn1 -> C.MYB_PN_1
      Pn2 -> C.MYB_PN_2
      Pn25 -> C.MYB_PN_25
      Pn3 -> C.MYB_PN_3
      Pn35 -> C.MYB_PN_35
      PnSo -> C.MYB_PN_SO
      PnSs -> C.MYB_PN_SS
      Pn4 -> C.MYB_PN_4
      Pn4rr -> C.MYB_PN_4RR
      PnSo25 -> C.MYB_PN_SO25
      PnSo35 -> C.MYB_PN_SO35
      PnSorr -> C.MYB_PN_SORR
      PnSs3 -> C.MYB_PN_SS3
      PnSorr5 -> C.MYB_PN_SORR5
      PnSotail -> C.MYB_PN_SOTAIL
      PnTailx -> C.MYB_PN_TAILX

-- | The terms beyond 1PN of two bodies.
addPnPair :: System -> Int -> Int -> C.Myb_system_pn_pair -> IO (Result ())
addPnPair system a b pair = call system $ \s -> with pair (C.myb_system_add_pn_pair s (index a) (index b))

-- | What an action sees: the body acted on, the positions, velocities and GM of all the bodies
-- as they are now, the rotation inertial -> body (9 values by rows) if it has a frame.
data SystemView = SystemView
  { viewBody :: !Int
  , viewPositions :: ![Vec3]
  , viewVelocities :: ![Vec3]
  , viewGm :: ![Double]
  , viewFrame :: !(Maybe [Double])
  }
  deriving (Eq, Show)

peekView :: Ptr C.Myb_system_view -> IO SystemView
peekView p = do
  v <- peek p
  let n = fromSize (C.myb_system_view'bodies v)
      frame = C.myb_system_view'frame v
  SystemView (fromSize (C.myb_system_view'body v))
    <$> peekArray n (castPtr (C.myb_system_view'x v))
    <*> peekArray n (castPtr (C.myb_system_view'v v))
    <*> peekDoubles n (C.myb_system_view'gm v)
    <*> (if frame == nullPtr then pure Nothing else Just <$> peekDoubles 9 frame)

-- | Anything else that pushes a body or turns it: a Haskell function of the time and the view
-- giving the acceleration (m/s^2, inertial axes) and the torque per unit mass (m^2/s^2, in
-- the frame of the body). The system keeps it.
addAction :: System -> Int -> (Double -> SystemView -> IO (Either Status (Vec3, Vec3))) -> IO (Result ())
addAction system@(System _ kept) body action = do
  function <- C.mk'myb_system_action act
  keepFunPtr kept function
  call system $ \s -> C.myb_system_add_action s (index body) function nullPtr
  where
    act _ t view acceleration torque = do
      seen <- peekView view
      given <- action (coerce t) seen
      case given of
        Right (a, n) -> poke (castPtr acceleration) a >> poke (castPtr torque) n >> pure C.MYB_OK
        Left status -> pure (statusCode status)

-- | Light on the mesh of a body (myb_system_light_action). The system keeps the description;
-- the mesh it points to (myb/surface.h, through "Myb.C") is the caller's: tie it with
-- 'systemKeep'.
addLightAction :: System -> Int -> C.Myb_system_light -> IO (Result ())
addLightAction system@(System _ kept) body light = do
  context <- keepValue kept light
  call system $ \s -> C.myb_system_add_action s (index body) C.p'myb_system_light_action (castPtr context)

-- | The spin of the body (its angular momentum per unit mass, m^2/s) drags the point.
addFrameDragging :: System -> Int -> Int -> Vec3 -> IO (Result ())
addFrameDragging system body point spin = call system $ \s -> withVec spin (C.myb_system_add_frame_dragging s (index body) (index point))

-- | A mass that changes with the time: rate in 1/s from the epoch.
setMassRate :: System -> Int -> Double -> Double -> IO (Result ())
setMassRate system body rate epoch = call system $ \s -> C.myb_system_set_mass_rate s (index body) (coerce rate) (coerce epoch)

-- | A body with a given motion that attracts the bodies of the system. Returns its number among
-- the guided bodies; the system keeps the provider.
addGuided :: System -> Double -> StateSource -> IO (Result Int)
addGuided (System p kept) gm source = do
  (function, context) <- keepStateSource kept source
  withForeignPtr p $ \s -> fmap fromSize <$> fill (C.myb_system_add_guided s (coerce gm) function context)

systemGuided :: System -> IO Int
systemGuided (System p _) = fromSize <$> withForeignPtr p C.myb_system_guided

-- | The body touches the surface of the ground with its contact points. The system keeps the
-- points.
addContact :: System -> Int -> Int -> C.Myb_contact_surface -> C.Myb_contact_law -> [C.Myb_contact_point] -> IO (Result ())
addContact system@(System _ kept) body ground surface law points = do
  memory <- mallocForeignPtrArray (max 1 (length points))
  withForeignPtr memory (`pokeArray` points)
  keepAlive kept memory
  let contact =
        C.zero
          { C.myb_system_contact'surface = surface
          , C.myb_system_contact'law = law
          , C.myb_system_contact'points = rawPtr memory
          , C.myb_system_contact'point_count = toSize (length points)
          }
  call system $ \s -> with contact (C.myb_system_add_contact s (index body) (index ground))

-- | The coordinates of a state (a state holds them and then their rates).
systemCoordinates :: System -> IO Int
systemCoordinates (System p _) = fromSize <$> withForeignPtr p C.myb_system_coordinates

-- | The rotation of a body into a state: the attitude, the angular velocity of the mantle and,
-- if any, of the core. Returns the state changed.
packRotation :: System -> Int -> Quaternion -> Vec3 -> Maybe Vec3 -> [Double] -> IO (Result [Double])
packRotation system@(System p _) body q omega core y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> with q $ \pq -> withVec omega $ \po -> withDoubles y $ \py -> do
    code <- case core of
      Nothing -> C.myb_system_pack_rotation s (index body) (castPtr pq) po nullPtr py
      Just c -> withVec c $ \pc -> C.myb_system_pack_rotation s (index body) (castPtr pq) po pc py
    andThen code (Right <$> peekDoubles (2 * n) py)

-- | The rotation of a body out of a state: the attitude, the angular velocities of the mantle
-- and of the core.
unpackRotation :: System -> Int -> [Double] -> IO (Result (Quaternion, Vec3, Vec3))
unpackRotation system@(System p _) body y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> withDoubles y $ \py -> alloca $ \pq -> alloca $ \po -> alloca $ \pc -> do
    code <- C.myb_system_unpack_rotation s (index body) py (castPtr pq) (castPtr po) (castPtr pc)
    andThen code (Right <$> ((,,) <$> peek pq <*> peek po <*> peek pc))

-- | The spins of a PN pair into a state. Returns the state changed.
packSpins :: System -> Int -> Vec3 -> Vec3 -> [Double] -> IO (Result [Double])
packSpins system@(System p _) pair chi1 chi2 y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> withVec chi1 $ \p1 -> withVec chi2 $ \p2 -> withDoubles y $ \py -> do
    code <- C.myb_system_pack_spins s (index pair) p1 p2 py
    andThen code (Right <$> peekDoubles (2 * n) py)

-- | The spins of a PN pair out of a state.
unpackSpins :: System -> Int -> [Double] -> IO (Result (Vec3, Vec3))
unpackSpins system@(System p _) pair y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> withDoubles y $ \py ->
    fill2 (\p1 p2 -> C.myb_system_unpack_spins s (index pair) py (castPtr p1) (castPtr p2))

-- | The deflections of the springs of a contact with the friction that sticks (in the order of
-- 'addContact'; one per point, axes of its surface) into a state. Returns the state changed.
packContact :: System -> Int -> [Vec3] -> [Double] -> IO (Result [Double])
packContact system@(System p _) contact deflections y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> withVecs deflections $ \pd -> withDoubles y $ \py -> do
    code <- C.myb_system_pack_contact s (index contact) pd py
    andThen code (Right <$> peekDoubles (2 * n) py)

-- | The deflections of the springs of a contact (of the number of points given) out of a state.
unpackContact :: System -> Int -> Int -> [Double] -> IO (Result [Vec3])
unpackContact system@(System p _) contact points y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> withDoubles y $ \py ->
    fillArray points (C.myb_system_unpack_contact s (index contact) py . castPtr)

-- | The second derivatives of the coordinates at t for the coordinates and their rates.
systemAcceleration :: System -> Double -> [Double] -> [Double] -> IO (Result [Double])
systemAcceleration system@(System p _) t x v = do
  n <- systemCoordinates system
  requireLength n x $ requireLength n v $ withForeignPtr p $ \s -> withDoubles x $ \px -> withDoubles v $ \pv ->
    fillDoubles n (C.myb_system_acceleration (castPtr s) (coerce t) px pv)

-- | The inertia tensor (9 values) of the whole body of an integrated rotation at the state (t, y).
systemInertia :: System -> Int -> Double -> [Double] -> IO (Result [Double])
systemInertia system@(System p _) body t y = do
  n <- systemCoordinates system
  requireLength (2 * n) y $ withForeignPtr p $ \s -> withDoubles y $ \py ->
    fillDoubles 9 (C.myb_system_inertia s (index body) (coerce t) py)

-- | The state at t1, d y1 / d y0 (2 n x 2 n by rows) and d y1 / d gm of the bodies asked
-- (2 n x count by rows).
data SystemPartials = SystemPartials {partialsState :: ![Double], partialsPhi :: ![Double], partialsSensitivity :: ![Double]}
  deriving (Eq, Show)

-- | The derivatives of the state at t1 with respect to the state at t0 and to the GM of the
-- bodies given, by central differences of whole propagations.
systemPartials :: System -> Maybe C.Myb_ias15_options -> Double -> [Double] -> Double -> [Int] -> IO (Result SystemPartials)
systemPartials system@(System p _) options t0 y0 t1 bodies = do
  n <- systemCoordinates system
  let m = 2 * n
      k = length bodies
  requireLength m y0 $ withForeignPtr p $ \s -> withOptions options $ \o -> withDoubles y0 $ \py0 -> withArray (map index bodies) $ \pb ->
    allocaArray (max 1 m) $ \y1 -> allocaArray (max 1 (m * m)) $ \phi -> allocaArray (max 1 (m * k)) $ \sens -> do
      code <- C.myb_system_partials s o (coerce t0) py0 (coerce t1) pb (toSize k) y1 phi sens
      andThen code (Right <$> (SystemPartials <$> peekDoubles m y1 <*> peekDoubles (m * m) phi <*> peekDoubles (m * k) sens))

withOptions :: Maybe C.Myb_ias15_options -> (Ptr C.Myb_ias15_options -> IO a) -> IO a
withOptions Nothing k = k nullPtr
withOptions (Just o) k = with o k

-- | An IAS15 integrator of a system, what it keeps (the system) and its coordinates.
data SystemIntegrator = SystemIntegrator (ForeignPtr C.Myb_ias15) Kept System Int

integratorPtr :: SystemIntegrator -> ForeignPtr C.Myb_ias15
integratorPtr (SystemIntegrator p _ _ _) = p

-- | An integrator of the system as it stands (myb_system_integrator: the positions, each
-- rotation and each pair of spins as groups of their own).
newIntegrator :: System -> Maybe C.Myb_ias15_options -> IO (Result SystemIntegrator)
newIntegrator system@(System p _) options = do
  n <- systemCoordinates system
  withForeignPtr p $ \s -> withOptions options $ \o -> do
    made <- newOwned (\i -> check <$> C.myb_system_integrator s o nullPtr i) C.myb_ias15_destroy
    case made of
      Left e -> pure (Left e)
      Right (fp, kept) -> keepAlive kept p >> pure (Right (SystemIntegrator fp kept system n))

-- | Starts the integrator at t0 on the accelerations of its system from the state y (2 n).
startIntegrator :: SystemIntegrator -> Double -> [Double] -> IO (Result ())
startIntegrator (SystemIntegrator p _ (System s _) n) t0 y =
  requireLength (2 * n) y $ withForeignPtr p $ \i -> withForeignPtr s $ \sp -> withDoubles y $ \py ->
    check <$> C.myb_ias15_start i C.p'myb_system_acceleration (castPtr sp) (coerce t0) py

-- | One step towards tLimit.
integratorStep :: SystemIntegrator -> Double -> IO (Result ())
integratorStep (SystemIntegrator p _ _ _) tLimit = withForeignPtr p $ \i -> check <$> C.myb_ias15_step i (coerce tLimit)

integratorIntegrate :: SystemIntegrator -> Double -> IO (Result ())
integratorIntegrate (SystemIntegrator p _ _ _) tEnd = withForeignPtr p $ \i -> check <$> C.myb_ias15_integrate i (coerce tEnd)

integratorTime :: SystemIntegrator -> IO Double
integratorTime (SystemIntegrator p _ _ _) = coerce <$> withForeignPtr p C.myb_ias15_time

-- | The state (2 n values: the coordinates, then their rates).
integratorState :: SystemIntegrator -> IO [Double]
integratorState (SystemIntegrator p _ _ n) = withForeignPtr p $ \i -> C.myb_ias15_state i >>= peekDoubles (2 * n)

-- | The state at t within the last step.
integratorDense :: SystemIntegrator -> Double -> IO (Result [Double])
integratorDense (SystemIntegrator p _ _ n) t = withForeignPtr p $ \i -> fillDoubles (2 * n) (C.myb_ias15_dense i (coerce t))

integratorStatistics :: SystemIntegrator -> IO C.Myb_ias15_stats
integratorStatistics (SystemIntegrator p _ _ _) = withForeignPtr p C.myb_ias15_statistics

-- | The radius of a body (a sphere that may touch).
setRadius :: System -> Int -> Double -> IO (Result ())
setRadius system body radius = call system $ \s -> C.myb_system_set_radius s (index body) (coerce radius)

systemRadius :: System -> Int -> IO Double
systemRadius (System p _) body = coerce <$> withForeignPtr p (\s -> C.myb_system_radius s (index body))

-- | False for a body merged into another.
systemPresent :: System -> Int -> IO Bool
systemPresent (System p _) body = (/= 0) <$> withForeignPtr p (\s -> C.myb_system_present s (index body))

data EventKind = Arrived | Contact | Approach
  deriving (Eq, Show)

-- | Where an advance stopped: the kind, the two bodies, the time and the distance of the centres.
data SystemEvent = SystemEvent
  { eventKind :: !EventKind
  , eventA :: !Int
  , eventB :: !Int
  , eventTime :: !Double
  , eventDistance :: !Double
  }
  deriving (Eq, Show)

-- | Steps the integrator towards tEnd to the first contact, or (watch > 0) the closest
-- approach of two bodies within watch of each other.
advance :: System -> SystemIntegrator -> Double -> Double -> IO (Result SystemEvent)
advance (System p _) (SystemIntegrator i _ _ _) tEnd watch =
  withForeignPtr p $ \s -> withForeignPtr i $ \ip -> fmap event <$> fillZero (C.myb_system_advance s ip (coerce tEnd) (coerce watch))
  where
    event e =
      SystemEvent
        { eventKind = case C.myb_system_event'kind e of
            C.MYB_SYSTEM_CONTACT -> Contact
            C.MYB_SYSTEM_APPROACH -> Approach
            _ -> Arrived
        , eventA = fromSize (C.myb_system_event'a e)
        , eventB = fromSize (C.myb_system_event'b e)
        , eventTime = coerce (C.myb_system_event't e)
        , eventDistance = coerce (C.myb_system_event'distance e)
        }

-- | Merges body b into body a at the state of the integrator, which starts again.
merge :: System -> SystemIntegrator -> Int -> Int -> IO (Result C.Myb_system_merger)
merge (System p _) (SystemIntegrator i _ _ _) a b =
  withForeignPtr p $ \s -> withForeignPtr i $ \ip -> fill (C.myb_system_merge s ip (index a) (index b))

-- | Two hard spheres rebound with the restitution given.
bounce :: System -> SystemIntegrator -> Int -> Int -> Double -> IO (Result ())
bounce (System p _) (SystemIntegrator i _ _ _) a b restitution =
  withForeignPtr p $ \s -> withForeignPtr i $ \ip -> check <$> C.myb_system_bounce s ip (index a) (index b) (coerce restitution)
