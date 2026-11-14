-- SPDX-License-Identifier: BSD-3-Clause
-- | What the checks share: a tally of checks and failures (no exception is thrown or caught),
-- and the comparisons of bits.
module Checks
  ( Tally (..)
  , newTally
  , expect
  , onRight
  , failsWith
  , bits
  , bitsOf
  , collect
  , dataDir
  , testDataDir
  ) where

import Control.Concurrent (threadDelay)
import Control.Monad (unless)
import Data.IORef
import GHC.Float (castDoubleToWord64)
import qualified Myb
import System.Environment (lookupEnv)
import System.FilePath ((</>))
import System.IO (hPutStrLn, stderr)
import System.Mem (performMajorGC)

data Tally = Tally (IORef Int) (IORef Int)

newTally :: IO Tally
newTally = Tally <$> newIORef 0 <*> newIORef 0

expect :: Tally -> String -> Bool -> IO ()
expect (Tally count failed) what ok = do
  modifyIORef' count (+ 1)
  unless ok $ do
    modifyIORef' failed (+ 1)
    hPutStrLn stderr ("FAILED: " ++ what)

-- | Goes on with the value of a successful outcome; a failure is a failed check.
onRight :: Tally -> String -> Myb.Result a -> (a -> IO ()) -> IO ()
onRight tally what outcome k = case outcome of
  Right value -> expect tally what True >> k value
  Left e -> expect tally (what ++ ": " ++ Myb.errorText e) False

-- | Whether the outcome is a failure with the status.
failsWith :: Myb.Status -> Myb.Result a -> Bool
failsWith status outcome = either ((== status) . Myb.errorStatus) (const False) outcome

bits :: Myb.State -> [Word]
bits = bitsOf . Myb.stateList

bitsOf :: [Double] -> [Word]
bitsOf = map (fromIntegral . castDoubleToWord64)

-- | Lets the finalizers of what is no longer reachable run.
collect :: IO ()
collect = performMajorGC >> threadDelay 20000 >> performMajorGC >> threadDelay 20000

dataDir, testDataDir :: IO FilePath
dataDir = maybe ("../.." </> "data") id <$> lookupEnv "MYB_DATA_DIR"
testDataDir = maybe ("../.." </> "tests" </> "data") id <$> lookupEnv "MYB_TEST_DATA_DIR"
