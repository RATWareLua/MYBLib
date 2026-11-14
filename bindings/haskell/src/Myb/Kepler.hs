-- SPDX-License-Identifier: BSD-3-Clause
-- | Two bodies (myb/kepler.h): states and elements, the state a time later. Functions of
-- numbers alone, so they are pure.
module Myb.Kepler
  ( Elements (..)
  , stateFromElements
  , elementsFromState
  , keplerPropagate
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Utils (with)
import Foreign.Storable (peek)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types
import System.IO.Unsafe (unsafeDupablePerformIO)

pvOf :: State -> C.Myb_pv
pvOf (State (Vec3 x y z) (Vec3 vx vy vz)) = C.zero {C.myb_pv'r = coerce [x, y, z], C.myb_pv'v = coerce [vx, vy, vz]}

-- | The state of a myb_pv: its arrays are of three numbers, as its Storable instance reads them.
stateOf :: C.Myb_pv -> State
stateOf pv = State (vec (C.myb_pv'r pv)) (vec (C.myb_pv'v pv))
  where
    vec xs = case coerce xs of
      [x, y, z] -> Vec3 x y z
      _ -> error "myb_pv: an array of other than three numbers"

-- | The state of the elements about a body of the gravitational parameter mu (m^3/s^2).
stateFromElements :: Double -> Elements -> Result State
stateFromElements mu elements = unsafeDupablePerformIO $
  with (elementsToC elements) $ \e -> alloca $ \out -> do
    code <- C.myb_pv_from_elements (coerce mu) e out
    either (pure . Left) (const (Right . stateOf <$> peek out)) (check code)

elementsFromState :: Double -> State -> Result Elements
elementsFromState mu y = unsafeDupablePerformIO $
  with (pvOf y) $ \pv -> alloca $ \out -> do
    code <- C.myb_elements_from_pv (coerce mu) pv out
    either (pure . Left) (const (Right . elementsFromC <$> peek out)) (check code)

-- | The state a time dt later on the orbit of two bodies.
keplerPropagate :: Double -> State -> Double -> Result State
keplerPropagate mu y dt = unsafeDupablePerformIO $
  with (pvOf y) $ \pv -> alloca $ \out -> do
    code <- C.myb_kepler_propagate (coerce mu) pv (coerce dt) out
    either (pure . Left) (const (Right . stateOf <$> peek out)) (check code)
