-- SPDX-License-Identifier: BSD-3-Clause
-- | Relative motion (myb/relative.h): the state of a deputy relative to a chief in the turning
-- frame of the chief, its transition matrix, a rendezvous of two impulses and the relative
-- orbital elements. Pure.
module Myb.Relative
  ( relativeFromInertial
  , relativeToInertial
  , relativeStm
  , cwStm
  , rendezvous
  , roeFromStates
  , roeToPosition
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import Foreign.Storable (peek)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | A state as myb_pv (the same six doubles in a row).
withPv :: State -> (Ptr C.Myb_pv -> IO a) -> IO a
withPv y k = with y (k . castPtr)

-- | The relative state of the deputy from the inertial states of both.
relativeFromInertial :: State -> State -> Result State
relativeFromInertial chief deputy = pureCall $ withPv chief $ \pc -> withPv deputy $ \pd -> fillVec (C.myb_relative_from_inertial pc pd)

relativeToInertial :: State -> State -> Result State
relativeToInertial chief relative = pureCall $ withPv chief $ \pc -> withVec relative $ \pr -> alloca $ \pd -> do
  code <- C.myb_relative_to_inertial pc pr (castPtr pd)
  andThen code (Right <$> peek pd)

-- | relative(t0 + dt) = phi relative(t0) for a chief on a two-body orbit of mu (36 by rows).
relativeStm :: Double -> State -> Double -> Result [Double]
relativeStm mu chief dt = pureCall $ withPv chief $ \pc -> fillDoubles 36 (C.myb_relative_stm (coerce mu) pc (coerce dt))

-- | The same for a chief on a circle of mean motion n (Clohessy-Wiltshire).
cwStm :: Double -> Double -> [Double]
cwStm n dt = pureCall $ allocaArray 36 $ \phi -> C.myb_cw_stm (coerce n) (coerce dt) phi >> peekDoubles 36 phi

-- | The two impulses of a rendezvous in the linear motion: from the relative state to the
-- relative position target after dt, arriving with the relative velocity given ('Nothing':
-- at rest).
rendezvous :: Double -> State -> Double -> State -> Vec3 -> Maybe Vec3 -> Result (Vec3, Vec3)
rendezvous mu chief dt relative target arrival = pureCall $
  withPv chief $ \pc -> withVec relative $ \pr -> withVec target $ \pt -> withMaybe arrival $ \pv ->
    fill2 (\d0 d1 -> C.myb_relative_rendezvous (coerce mu) pc (coerce dt) pr pt pv (castPtr d0) (castPtr d1))
  where
    withMaybe Nothing k = k nullPtr
    withMaybe (Just v) k = withVec v k

-- | The relative orbital elements (da, dlambda, dex, dey, dix, diy) of the deputy.
roeFromStates :: Double -> State -> State -> Result [Double]
roeFromStates mu chief deputy = pureCall $ withPv chief $ \pc -> withPv deputy $ \pd -> fillDoubles 6 (C.myb_roe_from_states (coerce mu) pc pd)

-- | The relative position (R, T, N) the elements give at the mean argument of latitude u of a
-- chief of semi-major axis a.
roeToPosition :: Double -> Double -> [Double] -> Result Vec3
roeToPosition a u roe
  | length roe /= 6 = failure Domain
  | otherwise = pureCall $ withDoubles roe $ \pr -> alloca $ \pos -> do
      C.myb_roe_to_position (coerce a) (coerce u) pr (castPtr pos)
      Right <$> peek pos
