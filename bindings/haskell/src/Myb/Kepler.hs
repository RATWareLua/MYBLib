-- SPDX-License-Identifier: BSD-3-Clause
-- | Two bodies (myb/kepler.h): states and elements, the state a time later. Functions of
-- numbers alone, so they are pure.
module Myb.Kepler
  ( Elements (..)
  , stateFromElements
  , elementsFromState
  , keplerPropagate
    -- * Array forms
  , C.Myb_equinoctial (..)
  , elementsFromStateMany
  , stateFromElementsMany
  , equinoctialFromStateMany
  , stateFromEquinoctialMany
  , keplerPropagateMany
  , keplerStmMany
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CSize)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, peekArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable, peek)
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

-- | An array form over a list: the outputs, or the failure of the element that failed.
manyOf :: (Storable a, Storable b) => [a] -> (Ptr a -> CSize -> Ptr b -> Ptr CSize -> IO C.Myb_status) -> Result [b]
manyOf xs f = pureCall $ withArrayLen xs $ \n px -> allocaArray (max 1 n) $ \o -> do
  done <- arrayForm (f px (toSize n) o)
  either (pure . Left) (const (Right <$> peekArray n o)) done

-- | 'elementsFromState' of every state.
elementsFromStateMany :: Double -> [State] -> Result [Elements]
elementsFromStateMany mu ys = map elementsFromC <$> manyOf (map pvOf ys) (C.myb_elements_from_pv_many (coerce mu))

-- | 'stateFromElements' of every set of elements.
stateFromElementsMany :: Double -> [Elements] -> Result [State]
stateFromElementsMany mu es = map stateOf <$> manyOf (map elementsToC es) (C.myb_pv_from_elements_many (coerce mu))

-- | The modified equinoctial elements of every state.
equinoctialFromStateMany :: Double -> [State] -> Result [C.Myb_equinoctial]
equinoctialFromStateMany mu ys = manyOf (map pvOf ys) (C.myb_equinoctial_from_pv_many (coerce mu))

-- | The state of every set of modified equinoctial elements.
stateFromEquinoctialMany :: Double -> [C.Myb_equinoctial] -> Result [State]
stateFromEquinoctialMany mu es = map stateOf <$> manyOf es (C.myb_pv_from_equinoctial_many (coerce mu))

-- | 'keplerPropagate' of one state to every time offset.
keplerPropagateMany :: Double -> State -> [Double] -> Result [State]
keplerPropagateMany mu y dts = pureCall $ with (pvOf y) $ \pv -> withDoublesLen dts $ \n pd -> allocaArray (max 1 n) $ \o -> do
  done <- arrayForm (C.myb_kepler_propagate_many (coerce mu) pv pd (toSize n) o)
  either (pure . Left) (const (Right . map stateOf <$> peekArray n o)) done

-- | The states and the transition matrices (36 by rows) of one state to every time offset.
keplerStmMany :: Double -> State -> [Double] -> Result [(State, [Double])]
keplerStmMany mu y dts = pureCall $ with (pvOf y) $ \pv -> withDoublesLen dts $ \n pd -> allocaArray (max 1 n) $ \o -> allocaArray (max 1 (36 * n)) $ \stm -> do
  done <- arrayForm (C.myb_kepler_stm_many (coerce mu) pv pd (toSize n) o stm)
  either (pure . Left) (const (Right <$> (zip <$> (map stateOf <$> peekArray n o) <*> (chunks <$> peekDoubles (36 * n) stm)))) done
  where
    chunks [] = []
    chunks xs = let (a, b) = splitAt 36 xs in a : chunks b

-- | The state a time dt later on the orbit of two bodies.
keplerPropagate :: Double -> State -> Double -> Result State
keplerPropagate mu y dt = unsafeDupablePerformIO $
  with (pvOf y) $ \pv -> alloca $ \out -> do
    code <- C.myb_kepler_propagate (coerce mu) pv (coerce dt) out
    either (pure . Left) (const (Right . stateOf <$> peek out)) (check code)
