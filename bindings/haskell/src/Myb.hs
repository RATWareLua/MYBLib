-- SPDX-License-Identifier: BSD-3-Clause
-- | myb for Haskell: the motion of objects in space.
--
-- Two layers over the C library:
--
-- * "Myb.C" declares every function, structure and constant of the C interface by the names
--   of the C headers, generated from them; it is in IO and unchecked, as C is.
-- * The modules re-exported here wrap the path most programs take: kernels and time,
--   ephemerides, gravity fields, the environment and the force model with its ready-made
--   terms, propagation, events, impulses, the transition matrix, the text messages. An object
--   of the library is released by its finalizer ('Foreign.ForeignPtr.ForeignPtr'); a call that
--   can fail returns @'Result' a@, that is @Either 'Error' a@ with the status of the library;
--   terms, providers and events are Haskell functions. A force model keeps alive the
--   environments and fields its terms point to, an environment its ephemerides and providers,
--   an orbit its model.
--
-- Nothing is thrown or caught here. A Haskell function the library calls reports a failure
-- through its status; an exception that leaves it is not caught: it ends the program, as a
-- broken contract stops at once in the library (MYB_CHECK).
--
-- No arithmetic is done here: a number obtained through this package has the bits of the
-- same calls made in C. States are position (m) and velocity (m/s); the times of a model are
-- seconds from the epoch of its environment.
--
-- The rest of the interface has its modules too: N bodies ("Myb.NBody", "Myb.System"),
-- estimation ("Myb.Estimation"), the restricted three-body problem ("Myb.Cr3bp"), targeting and
-- optimal transfers ("Myb.Transfer", "Myb.Optimal", "Myb.Thrust", "Myb.Iod", "Myb.Relative"),
-- conjunctions and uncertainty ("Myb.Conjunction", "Myb.Uncertainty"), attitude
-- ("Myb.Attitude"), observables ("Myb.Observe"), files ("Myb.Files", "Myb.Ck", "Myb.Message"),
-- providers ("Myb.Source") and the array forms of frames ("Myb.Frames"); two-line elements are
-- in the library myb:tle ("Myb.Tle"). What has no upper layer is taken through "Myb.C"; the two
-- mix freely: each object gives its 'Foreign.ForeignPtr.ForeignPtr' (forceModelPtr, ...).
--
-- > import Control.Monad.Trans.Except (ExceptT (..), runExceptT)
-- > import qualified Myb
-- >
-- > main :: IO ()
-- > main = do
-- >   later <- runExceptT $ do
-- >     model <- ExceptT Myb.newForceModel
-- >     _ <- ExceptT (Myb.addPointMass model "planet" 3.986004418e14)
-- >     orbit <- ExceptT (Myb.newOrbit model Myb.defaultOrbitOptions)
-- >     ExceptT (Myb.start orbit 0 (Myb.State (Myb.Vec3 7.0e6 0 0) (Myb.Vec3 0 7.5e3 0)) 6000)
-- >     ExceptT (Myb.stateAt orbit 6000)
-- >   either (putStrLn . Myb.errorText) print later
--
-- (bindings/haskell/examples/LowOrbit.hs is the example of the low orbit.)
module Myb
  ( module Myb.Error
  , module Myb.Types
  , module Myb.Time
  , module Myb.Kepler
  , module Myb.Ephemeris
  , module Myb.Gravity
  , module Myb.Environment
  , module Myb.Force
  , module Myb.Orbit
  , module Myb.Message
  , module Myb.Source
  , module Myb.NBody
  , module Myb.System
  , module Myb.Estimation
  , module Myb.Cr3bp
  , module Myb.Transfer
  , module Myb.Optimal
  , module Myb.Thrust
  , module Myb.Iod
  , module Myb.Relative
  , module Myb.Conjunction
  , module Myb.Uncertainty
  , module Myb.Attitude
  , module Myb.Contact
  , module Myb.Observe
  , module Myb.Files
  , module Myb.Ck
  , module Myb.Frames
  ) where

import Myb.Attitude
import Myb.Contact
import Myb.Ck
import Myb.Conjunction
import Myb.Cr3bp
import Myb.Environment
import Myb.Ephemeris
import Myb.Error
import Myb.Estimation
import Myb.Files
import Myb.Force
import Myb.Frames
import Myb.Gravity
import Myb.Iod
import Myb.Kepler
import Myb.Message
import Myb.NBody
import Myb.Observe
import Myb.Optimal
import Myb.Orbit
import Myb.Relative
import Myb.Source
import Myb.System
import Myb.Thrust
import Myb.Time
import Myb.Transfer
import Myb.Types
import Myb.Uncertainty
