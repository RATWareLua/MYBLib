-- SPDX-License-Identifier: BSD-3-Clause
-- | The example of the README: a satellite in a low orbit about an oblate planet with an
-- atmosphere (examples/low_orbit.c), its first ascending node and its state a day later.
-- Every call returns its outcome; ExceptT stops at the first failure and the program says it.
module Main (main) where

import Control.Monad.IO.Class (liftIO)
import Control.Monad.Trans.Except (ExceptT (..), runExceptT)
import Myb (State (..), Vec3 (..))
import qualified Myb
import System.Exit (exitFailure)

main :: IO ()
main = do
  let gm = 3.986004418e14
      radius = 6378137.0
  done <- runExceptT $ do
    field <- ExceptT (Myb.newGravityField gm radius 2 [1, 0, 0, -1.08263e-3 / sqrt 5, 0, 0] (replicate 6 0))
    env <- ExceptT (Myb.newEnvironment (Myb.timeFromSeconds 0))
    planet <- ExceptT (Myb.addOrigin env Myb.defaultBody {Myb.bodyId = 1, Myb.bodyGm = gm})
    model <- ExceptT Myb.newForceModel
    _ <- ExceptT (Myb.addField model "gravity" env planet field 0 Nothing)
    _ <-
      ExceptT
        ( Myb.addDrag model "drag" env planet (Myb.ExponentialAtmosphere radius 0 [(400e3, 3e-12, 58e3)])
            (Myb.Object 500 4 1.3 2.2)
        )
    orbit <- ExceptT (Myb.newOrbit model Myb.defaultOrbitOptions)
    y0 <- ExceptT (pure (Myb.stateFromElements gm (Myb.Elements 6.778e6 0.001 0.9 0.5 1.0 0)))
    ExceptT (Myb.start orbit 0 y0 86400)
    node <- ExceptT (Myb.nextEvent orbit (\_ y -> pure (let Vec3 _ _ z = position y in z)) 86400 Myb.Rising)
    liftIO (print node) -- Just (time, state), or Nothing
    ExceptT (Myb.stateAt orbit 86400)
  either (\e -> putStrLn (Myb.errorText e) >> exitFailure) print done
