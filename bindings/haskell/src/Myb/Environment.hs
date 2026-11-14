-- SPDX-License-Identifier: BSD-3-Clause
-- | The bodies of a scenario (myb/env.h): what the terms of a force model ask for. A body is
-- at the origin, or its state relative to the origin comes from an SPK set or from a Haskell
-- function. The environment keeps alive what its bodies point to.
module Myb.Environment
  ( Environment
  , newEnvironment
  , Body (..)
  , defaultBody
  , addOrigin
  , addFromEphemeris
  , addWith
  , findBody
  , bodyState
  , bodyCount
  , environmentKeep
  , environmentPtr
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Data.Int (Int32, Int64)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, castPtrToFunPtr, nullFunPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Ephemeris (Ephemeris, ephemerisPtr)
import Myb.Error
import Myb.Internal
import Myb.Types

-- | An environment: its object and what it keeps.
data Environment = Environment (ForeignPtr C.Myb_env) Kept

-- | An environment whose model time 0 is the instant given (one time scale: that of the
-- instants of its bodies).
newEnvironment :: Time -> IO (Result Environment)
newEnvironment epoch = do
  kept <- newKept
  with (C.zero {C.myb_env_options'epoch = timeToC epoch}) $ \options ->
    create
      (\out -> check <$> C.myb_env_create options nullPtr out)
      (\p -> (`Environment` kept) <$> owned p C.myb_env_destroy kept)

-- | A body: where it is comes with the function that adds it; this is what it is to the terms
-- of a model.
data Body = Body
  { -- | The caller's identifier (a NAIF code, say).
    bodyId :: !Int64
  , -- | m^3 s^-2
    bodyGm :: !Double
  , -- | m: the sphere that casts shadows and bounds the disk of a source.
    bodyRadius :: !Double
  , -- | W: the body as a source of radiation.
    bodyLuminosity :: !Double
  , -- | The limb darkening u1, u2 of the disk of a source.
    bodyLimb :: !(Double, Double)
  }
  deriving (Eq, Show)

-- | All zeros: a record update gives what matters, @defaultBody {bodyId = 399, bodyGm = gm}@.
defaultBody :: Body
defaultBody = Body 0 0 0 0 (0, 0)

add :: Environment -> Body -> C.Myb_state_function -> Ptr () -> IO (Result Int)
add (Environment env _) body function context =
  with made $ \b -> withForeignPtr env $ \e -> alloca $ \index -> do
    code <- C.myb_env_add_body e b index
    either (pure . Left) (const (Right . fromSize <$> peek index)) (check code)
  where
    (u1, u2) = bodyLimb body
    made =
      C.zero
        { C.myb_env_body'id = bodyId body
        , C.myb_env_body'gm = coerce (bodyGm body)
        , C.myb_env_body'state = function
        , C.myb_env_body'state_context = context
        , C.myb_env_body'radius = coerce (bodyRadius body)
        , C.myb_env_body'luminosity = coerce (bodyLuminosity body)
        , C.myb_env_body'limb = coerce [u1, u2]
        }

-- | The body at the origin of the environment. Returns its index.
addOrigin :: Environment -> Body -> IO (Result Int)
addOrigin environment body = add environment body nullFunPtr nullPtr

-- | A body whose state relative to the origin is that of the target relative to the center in
-- an SPK set (frame 0: whatever frame the segments are in). The environment keeps the set.
addFromEphemeris :: Environment -> Body -> Ephemeris -> Int32 -> Int32 -> Int32 -> IO (Result Int)
addFromEphemeris environment@(Environment _ kept) body ephemeris target center frame = do
  let spk = ephemerisPtr ephemeris
  context <-
    keepValue kept $
      C.zero
        { C.myb_spk_body'spk = rawPtr spk
        , C.myb_spk_body'target = target
        , C.myb_spk_body'center = center
        , C.myb_spk_body'frame = frame
        }
  keepAlive kept spk
  add environment body C.p'myb_spk_body_state (castPtr context)

-- | A body whose state relative to the origin a Haskell function gives, at an instant in
-- seconds past J2000 of the model's scale. The function reports a failure by its status; the
-- environment keeps it.
addWith :: Environment -> Body -> (Double -> IO (Either Status State)) -> IO (Result Int)
addWith environment@(Environment _ kept) body provider = do
  (function, context) <- C.mk'myb_state_function call
  keepFunPtr kept (castPtrToFunPtr context)
  add environment body function context
  where
    call _ t out = do
      at <- C.myb_time_seconds t
      given <- provider (coerce at)
      case given of
        Right y -> poke (castPtr out) y >> pure C.MYB_OK
        Left status -> pure (statusCode status)

-- | The index of the body with the identifier.
findBody :: Environment -> Int64 -> IO (Result Int)
findBody (Environment env _) identifier =
  withForeignPtr env $ \e -> alloca $ \index -> do
    code <- C.myb_env_find e identifier index
    either (pure . Left) (const (Right . fromSize <$> peek index)) (check code)

-- | The state of a body relative to the origin at the model time t.
bodyState :: Environment -> Double -> Int -> IO (Result State)
bodyState environment@(Environment env _) t body = do
  count <- bodyCount environment
  if body < 0 || body >= count
    then pure (failure Domain)
    else withForeignPtr env $ \e -> alloca $ \out -> do
      code <- C.myb_env_state e (coerce t) (toSize body) (castPtr out)
      either (pure . Left) (const (Right <$> peek out)) (check code)

bodyCount :: Environment -> IO Int
bodyCount (Environment env _) = fromSize <$> withForeignPtr env C.myb_env_body_count

-- | Ties the life of an object to the environment: for what a body added through "Myb.C"
-- points to.
environmentKeep :: Environment -> ForeignPtr a -> IO ()
environmentKeep (Environment _ kept) = keepAlive kept

environmentPtr :: Environment -> ForeignPtr C.Myb_env
environmentPtr (Environment env _) = env
