-- SPDX-License-Identifier: BSD-3-Clause
-- | The values the upper layer passes: vectors and states of doubles, orbital elements,
-- calendar labels and instants. A conversion to or from a structure of C moves the same bits
-- (a CDouble is a Double): nothing is computed.
module Myb.Types
  ( Vec3 (..)
  , State (..)
  , stateList
  , stateFromList
  , Quaternion (..)
  , Mat3 (..)
  , mat3List
  , Elements (..)
  , Calendar (..)
  , Time (..)
  , timeFromSeconds
  , seconds
    -- * Structures of C
  , elementsToC
  , elementsFromC
  , calendarToC
  , calendarFromC
  , timeToC
  , timeFromC
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..))
import Data.Int (Int64)
import Foreign.Ptr (castPtr)
import Foreign.Storable (Storable (..))
import qualified Myb.C as C
import System.IO.Unsafe (unsafeDupablePerformIO)

-- | Three numbers: a position (m), a velocity (m/s), an acceleration (m/s^2).
data Vec3 = Vec3 !Double !Double !Double
  deriving (Eq, Show)

-- | Position (m) and velocity (m/s), relative to the origin of a model.
data State = State {position :: !Vec3, velocity :: !Vec3}
  deriving (Eq, Show)

-- | The six numbers of a state in the order of C: x, y, z, vx, vy, vz.
stateList :: State -> [Double]
stateList (State (Vec3 x y z) (Vec3 vx vy vz)) = [x, y, z, vx, vy, vz]

-- | Three doubles in a row, as C writes double[3].
instance Storable Vec3 where
  sizeOf _ = 3 * sizeOf (0 :: Double)
  alignment _ = alignment (0 :: Double)
  peek p = Vec3 <$> peekElemOff q 0 <*> peekElemOff q 1 <*> peekElemOff q 2
    where
      q = castPtr p
  poke p (Vec3 x y z) = pokeElemOff q 0 x >> pokeElemOff q 1 y >> pokeElemOff q 2 z
    where
      q = castPtr p

-- | Six doubles in a row, as C writes double[6].
instance Storable State where
  sizeOf _ = 2 * sizeOf (Vec3 0 0 0)
  alignment _ = alignment (0 :: Double)
  peek p = State <$> peek (castPtr p) <*> peekElemOff (castPtr p) 1
  poke p (State r v) = poke (castPtr p) r >> pokeElemOff (castPtr p) 1 v

-- | A state from six numbers in the order of C; 'Nothing' for another count.
stateFromList :: [Double] -> Maybe State
stateFromList [x, y, z, vx, vy, vz] = Just (State (Vec3 x y z) (Vec3 vx vy vz))
stateFromList _ = Nothing

-- | A quaternion in the order of myb/attitude.h: the scalar first, (w, x, y, z).
data Quaternion = Quaternion !Double !Double !Double !Double
  deriving (Eq, Show)

-- | Four doubles in a row, as C writes double[4].
instance Storable Quaternion where
  sizeOf _ = 4 * sizeOf (0 :: Double)
  alignment _ = alignment (0 :: Double)
  peek p = Quaternion <$> peekElemOff q 0 <*> peekElemOff q 1 <*> peekElemOff q 2 <*> peekElemOff q 3
    where
      q = castPtr p
  poke p (Quaternion w x y z) = pokeElemOff q 0 w >> pokeElemOff q 1 x >> pokeElemOff q 2 y >> pokeElemOff q 3 z
    where
      q = castPtr p

-- | A 3 x 3 matrix by its rows, as C writes double[9] (row-major).
data Mat3 = Mat3 !Vec3 !Vec3 !Vec3
  deriving (Eq, Show)

-- | The nine numbers by rows.
mat3List :: Mat3 -> [Double]
mat3List (Mat3 (Vec3 a b c) (Vec3 d e f) (Vec3 g h i)) = [a, b, c, d, e, f, g, h, i]

instance Storable Mat3 where
  sizeOf _ = 3 * sizeOf (Vec3 0 0 0)
  alignment _ = alignment (0 :: Double)
  peek p = Mat3 <$> peekElemOff q 0 <*> peekElemOff q 1 <*> peekElemOff q 2
    where
      q = castPtr p
  poke p (Mat3 a b c) = pokeElemOff q 0 a >> pokeElemOff q 1 b >> pokeElemOff q 2 c
    where
      q = castPtr p

-- | Classical elements.
data Elements = Elements
  { -- | p (m)
    semiLatusRectum :: !Double
  , eccentricity :: !Double
  , -- | (rad)
    inclination :: !Double
  , -- | the right ascension of the ascending node (rad)
    raan :: !Double
  , -- | the argument of periapsis (rad)
    argumentOfPeriapsis :: !Double
  , -- | (rad)
    trueAnomaly :: !Double
  }
  deriving (Eq, Show)

-- | A calendar label: year, month, day, hour, minute, second.
data Calendar = Calendar
  { year :: !Int64
  , month :: !Int
  , day :: !Int
  , hour :: !Int
  , minute :: !Int
  , second :: !Double
  }
  deriving (Eq, Show)

-- | An instant: whole seconds and a fraction past J2000 of its time scale.
data Time = Time {timeSeconds :: !Int64, timeFraction :: !Double}
  deriving (Eq, Show)

-- | An instant from seconds past J2000 (myb_time_from_seconds).
timeFromSeconds :: Double -> Time
timeFromSeconds s = timeFromC (unsafeDupablePerformIO (C.myb_time_from_seconds (coerce s)))

-- | The seconds past J2000 of an instant, as one number (myb_time_seconds).
seconds :: Time -> Double
seconds t = coerce (unsafeDupablePerformIO (C.myb_time_seconds (timeToC t)))

elementsToC :: Elements -> C.Myb_elements
elementsToC e =
  C.zero
    { C.myb_elements'p = coerce (semiLatusRectum e)
    , C.myb_elements'e = coerce (eccentricity e)
    , C.myb_elements'i = coerce (inclination e)
    , C.myb_elements'raan = coerce (raan e)
    , C.myb_elements'argp = coerce (argumentOfPeriapsis e)
    , C.myb_elements'nu = coerce (trueAnomaly e)
    }

elementsFromC :: C.Myb_elements -> Elements
elementsFromC e =
  Elements
    { semiLatusRectum = coerce (C.myb_elements'p e)
    , eccentricity = coerce (C.myb_elements'e e)
    , inclination = coerce (C.myb_elements'i e)
    , raan = coerce (C.myb_elements'raan e)
    , argumentOfPeriapsis = coerce (C.myb_elements'argp e)
    , trueAnomaly = coerce (C.myb_elements'nu e)
    }

calendarToC :: Calendar -> C.Myb_calendar
calendarToC c =
  C.zero
    { C.myb_calendar'year = year c
    , C.myb_calendar'month = fromIntegral (month c)
    , C.myb_calendar'day = fromIntegral (day c)
    , C.myb_calendar'hour = fromIntegral (hour c)
    , C.myb_calendar'minute = fromIntegral (minute c)
    , C.myb_calendar'second = coerce (second c)
    }

calendarFromC :: C.Myb_calendar -> Calendar
calendarFromC c =
  Calendar
    { year = C.myb_calendar'year c
    , month = fromIntegral (C.myb_calendar'month c)
    , day = fromIntegral (C.myb_calendar'day c)
    , hour = fromIntegral (C.myb_calendar'hour c)
    , minute = fromIntegral (C.myb_calendar'minute c)
    , second = coerce (C.myb_calendar'second c)
    }

timeToC :: Time -> C.Myb_time
timeToC t = C.zero {C.myb_time's = timeSeconds t, C.myb_time'f = coerce (timeFraction t)}

timeFromC :: C.Myb_time -> Time
timeFromC t = Time {timeSeconds = C.myb_time's t, timeFraction = coerce (C.myb_time'f t)}
