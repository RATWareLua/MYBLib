-- SPDX-License-Identifier: BSD-3-Clause
-- | Providers: where the state of a point (myb_state_function of myb/kin.h) or the orientation
-- of a frame (myb_orientation_function of myb/frames.h) comes from, for the functions that
-- ask for one — the guided bodies of a system, the bodies of a chain of legs, light time,
-- conjunctions, the Chebyshev segments of an SPK writer. A body of an SPK set is passed to the
-- library as its own provider (myb_spk_body_state), so that nothing is converted on the way;
-- a Haskell function is called with the seconds past J2000 of the time scale of the call.
module Myb.Source
  ( StateSource
  , fromEphemeris
  , fromFunction
  , fromC
  , OrientationSource
  , orientationFunction
  , orientationC
  ) where

import Data.Coerce (coerce)
import Data.Int (Int32)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, castPtrToFunPtr, freeHaskellFunPtr)
import Foreign.Storable (poke)
import qualified Myb.C as C
import Myb.Ephemeris (Ephemeris, ephemerisPtr)
import Myb.Error
import Myb.Internal
import Myb.Types

-- | Target relative to center in an SPK set (frame 0: whatever frame the segments are in).
fromEphemeris :: Ephemeris -> Int32 -> Int32 -> Int32 -> StateSource
fromEphemeris ephemeris target center frame = StateSource for keep
  where
    spk = ephemerisPtr ephemeris
    body p = C.zero {C.myb_spk_body'spk = p, C.myb_spk_body'target = target, C.myb_spk_body'center = center, C.myb_spk_body'frame = frame}
    for k = withForeignPtr spk $ \p -> with (body p) $ \b -> k C.p'myb_spk_body_state (castPtr b)
    keep kept = do
      context <- keepValue kept (body (rawPtr spk))
      keepAlive kept spk
      pure (C.p'myb_spk_body_state, castPtr context)

-- | A Haskell function of the seconds past J2000; it reports a failure by its status.
fromFunction :: (Double -> IO (Either Status State)) -> StateSource
fromFunction f = StateSource for keep
  where
    for k = do
      (function, context) <- C.mk'myb_state_function call
      done <- k function context
      freeHaskellFunPtr (castPtrToFunPtr context)
      pure done
    keep kept = do
      (function, context) <- C.mk'myb_state_function call
      keepFunPtr kept (castPtrToFunPtr context)
      pure (function, context)
    call _ t y = do
      at <- C.myb_time_seconds t
      given <- f (coerce at)
      case given of
        Right state -> poke (castPtr y) state >> pure C.MYB_OK
        Left status -> pure (statusCode status)

-- | A provider of C and its context, which the caller keeps alive as long as it is used.
fromC :: C.Myb_state_function -> Ptr () -> StateSource
fromC function context = StateSource (\k -> k function context) (\_ -> pure (function, context))

-- | A Haskell function of the seconds past J2000 giving the orientation (myb/frames.h).
orientationFunction :: (Double -> IO (Either Status C.Myb_orientation)) -> OrientationSource
orientationFunction f = OrientationSource for keep
  where
    for k = do
      (function, context) <- C.mk'myb_orientation_function call
      done <- k function context
      freeHaskellFunPtr (castPtrToFunPtr context)
      pure done
    keep kept = do
      (function, context) <- C.mk'myb_orientation_function call
      keepFunPtr kept (castPtrToFunPtr context)
      pure (function, context)
    call _ t o = do
      at <- C.myb_time_seconds t
      given <- f (coerce at)
      case given of
        Right orientation -> poke o orientation >> pure C.MYB_OK
        Left status -> pure (statusCode status)

-- | An orientation provider of C (myb_orientation_cip, ...) and its context, which the caller
-- keeps alive as long as it is used.
orientationC :: C.Myb_orientation_function -> Ptr () -> OrientationSource
orientationC function context = OrientationSource (\k -> k function context) (\_ -> pure (function, context))
