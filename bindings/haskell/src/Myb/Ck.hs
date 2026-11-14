-- SPDX-License-Identifier: BSD-3-Clause
-- | Orientation files and spacecraft clocks (myb/ck.h): SCLK of type 1, CK files read (types 1
-- to 6) and their pointings, an instrument frame as an orientation provider, and the writer of
-- CK files.
module Myb.Ck
  ( -- * Spacecraft clocks
    Sclk
  , newSclk
  , sclkTimeSystem
  , sclkPartitions
  , sclkToTime
  , sclkFromTime
  , sclkEncode
  , sclkDecode
    -- * Reading
  , Ck
  , newCk
  , loadCkFile
  , loadCkBytes
  , C.Myb_ck_segment_info (..)
  , ckSegments
  , C.Myb_ck_pointing (..)
  , ckPointing
  , ckOrientation
    -- * Writing
  , CkWriter
  , newCkWriter
  , CkHeader (..)
  , CkPackets (..)
  , ckAddPointings
  , ckAddIntervals
  , ckAddChebyshev
  , ckAddPackets
  , ckAddMinisegments
  , ckWriterSize
  , ckWriterBytes
  , ckWriterFile
  , sclkPtr
  , ckPtr
  ) where

import qualified Data.ByteString as B
import qualified Data.ByteString.Unsafe as BU
import Data.Coerce (coerce)
import Data.Int (Int32, Int64)
import Foreign.C.Types (CDouble (..), CSize)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (allocaBytes)
import Foreign.Marshal.Array (withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Time (KernelPool, kernelPoolPtr)
import Myb.Types

-- ---- Clocks --------------------------------------------------------------------------------

newtype Sclk = Sclk (ForeignPtr C.Myb_sclk)

sclkPtr :: Sclk -> ForeignPtr C.Myb_sclk
sclkPtr (Sclk p) = p

-- | The clock whose kernel variables end in _<id> in the pool (what is needed is copied).
newSclk :: KernelPool -> Int64 -> IO (Result Sclk)
newSclk pool identifier = withForeignPtr (kernelPoolPtr pool) $ \p ->
  create (\o -> withDiag (C.myb_sclk_create p identifier nullPtr o)) (fmap Sclk . ownedLeaf C.p'myb_sclk_destroy)

-- | The scale of the times of the table: 1 = TDB, 2 = TT.
sclkTimeSystem :: Sclk -> IO Int
sclkTimeSystem (Sclk s) = fromIntegral <$> withForeignPtr s C.myb_sclk_time_system

sclkPartitions :: Sclk -> IO Int
sclkPartitions (Sclk s) = fromSize <$> withForeignPtr s C.myb_sclk_partitions

-- | Encoded ticks -> the instant on the clock's time scale.
sclkToTime :: Sclk -> Double -> IO (Result Time)
sclkToTime (Sclk s) ticks = withForeignPtr s $ \p -> fmap timeFromC <$> fill (C.myb_sclk_to_time p (coerce ticks))

-- | The instant -> encoded ticks.
sclkFromTime :: Sclk -> Time -> IO (Result Double)
sclkFromTime (Sclk s) t = withForeignPtr s $ \p -> fmap coerce <$> fill (C.myb_sclk_from_time p (timeToC t))

-- | A clock string ("partition/field field ...") -> encoded ticks.
sclkEncode :: Sclk -> String -> IO (Result Double)
sclkEncode (Sclk s) text = withTextLen text $ \pt n -> withForeignPtr s $ \p -> fmap coerce <$> fill (C.myb_sclk_encode p pt n)

-- | Encoded ticks -> the clock string with the partition.
sclkDecode :: Sclk -> Double -> IO (Result String)
sclkDecode (Sclk s) ticks = withForeignPtr s $ \p -> writeText (C.myb_sclk_decode p (coerce ticks))

-- ---- Reading -------------------------------------------------------------------------------

newtype Ck = Ck (ForeignPtr C.Myb_ck)

ckPtr :: Ck -> ForeignPtr C.Myb_ck
ckPtr (Ck p) = p

newCk :: IO (Result Ck)
newCk = create (\o -> check <$> C.myb_ck_create nullPtr o) (fmap Ck . ownedLeaf C.p'myb_ck_destroy)

loadCkFile :: Ck -> FilePath -> IO (Result ())
loadCkFile (Ck c) path = withText path $ \p -> withForeignPtr c $ \pc -> withDiag (C.myb_ck_load_file pc p)

-- | A file in memory.
loadCkBytes :: Ck -> B.ByteString -> IO (Result ())
loadCkBytes (Ck c) bytes = BU.unsafeUseAsCStringLen bytes $ \(p, n) -> withForeignPtr c $ \pc -> withDiag (C.myb_ck_load_memory pc (castPtr p) (toSize n))

-- | The segments of the files, in the order of loading.
ckSegments :: Ck -> IO [C.Myb_ck_segment_info]
ckSegments (Ck c) = withForeignPtr c $ \p -> do
  n <- C.myb_ck_segment_count p
  if n == 0 then pure [] else mapM (\i -> either (const C.zero) id <$> fillZero (\o -> C.myb_ck_segment p i o >> pure C.MYB_OK)) [0 .. n - 1]

-- | The orientation of an instrument at ticks of its clock, or the nearest within the
-- tolerance (ticks); need the rate: pass over segments without angular velocities.
ckPointing :: Ck -> Int32 -> Double -> Double -> Bool -> IO (Result C.Myb_ck_pointing)
ckPointing (Ck c) instrument ticks tolerance needRate = withForeignPtr c $ \p ->
  fillZero (C.myb_ck_pointing_at p instrument (coerce ticks) (coerce tolerance) (if needRate then 1 else 0))

-- | An instrument frame as an orientation provider (myb_orientation_ck): the files and the
-- clock, the instrument, the base frame (0: any) and the tolerance (ticks). The provider keeps
-- the files and the clock alive.
ckOrientation :: Ck -> Sclk -> Int32 -> Int32 -> Double -> OrientationSource
ckOrientation (Ck c) (Sclk s) instrument base tolerance = OrientationSource for keep
  where
    record pc ps =
      C.zero
        { C.myb_ck_orientation'ck = pc
        , C.myb_ck_orientation'sclk = ps
        , C.myb_ck_orientation'instrument = instrument
        , C.myb_ck_orientation'base = base
        , C.myb_ck_orientation'tolerance = coerce tolerance
        }
    for k = withForeignPtr c $ \pc -> withForeignPtr s $ \ps -> with (record pc ps) $ \r -> k C.p'myb_orientation_ck (castPtr r)
    keep kept = do
      context <- keepValue kept (record (rawPtr c) (rawPtr s))
      keepAlive kept c
      keepAlive kept s
      pure (C.p'myb_orientation_ck, castPtr context)

-- ---- Writing -------------------------------------------------------------------------------

newtype CkWriter = CkWriter (ForeignPtr C.Myb_ck_writer)

-- | A writer with the label of the file ('Nothing': "myb").
newCkWriter :: Maybe String -> IO (Result CkWriter)
newCkWriter name = withName $ \p -> create (\o -> check <$> C.myb_ck_writer_create p nullPtr o) (fmap CkWriter . ownedLeaf C.p'myb_ck_writer_destroy)
  where
    withName k = maybe (k nullPtr) (`withText` k) name

-- | The instrument, the base frame (1 = J2000) and the name of a segment.
data CkHeader = CkHeader {ckInstrument :: !Int32, ckFrame :: !Int32, ckName :: !String}
  deriving (Eq, Show)

withHeader :: CkHeader -> (Ptr C.Myb_ck_segment_header -> IO (Result a)) -> IO (Result a)
withHeader h k = withText (ckName h) $ \name ->
  with C.zero {C.myb_ck_segment_header'instrument = ckInstrument h, C.myb_ck_segment_header'frame = ckFrame h, C.myb_ck_segment_header'name = name} k

withOptional :: (a -> (Ptr b -> IO r) -> IO r) -> Maybe a -> (Ptr b -> IO r) -> IO r
withOptional _ Nothing k = k nullPtr
withOptional f (Just x) k = f x k

-- | Pointings at increasing ticks: type 1 or 3, quaternions, rates ('Nothing': none), the
-- starts of the intervals of interpolation ('Nothing': one interval).
ckAddPointings :: CkWriter -> CkHeader -> Int -> [Double] -> [Quaternion] -> Maybe [Vec3] -> Maybe [Double] -> IO (Result ())
ckAddPointings (CkWriter w) header kind ticks quaternions rates starts
  | length quaternions /= n || maybe False ((/= n) . length) rates = pure (failure Domain)
  | otherwise = withForeignPtr w $ \pw -> withHeader header $ \ph -> withDoubles ticks $ \pt -> withVecs quaternions $ \pq ->
      withOptional withVecs rates $ \pr -> withOptional withDoubles starts $ \ps ->
        check <$> C.myb_ck_writer_add_pointings pw ph (fromIntegral kind) pt pq pr (toSize n) ps (toSize (maybe 0 length starts))
  where
    n = length ticks

-- | Intervals of type 2: starts, stops, the quaternion at each start, the angular velocity and
-- the seconds a tick lasts.
ckAddIntervals :: CkWriter -> CkHeader -> [Double] -> [Double] -> [Quaternion] -> [Vec3] -> [Double] -> IO (Result ())
ckAddIntervals (CkWriter w) header starts stops quaternions rates tickSeconds
  | any (/= n) [length stops, length quaternions, length rates, length tickSeconds] = pure (failure Domain)
  | otherwise = withForeignPtr w $ \pw -> withHeader header $ \ph -> withDoubles starts $ \pa -> withDoubles stops $ \pb -> withVecs quaternions $ \pq -> withVecs rates $ \pr -> withDoubles tickSeconds $ \ps ->
      check <$> C.myb_ck_writer_add_intervals pw ph pa pb pq pr ps (toSize n)
  where
    n = length starts

-- | Intervals of type 4: Chebyshev series of the quaternion (and the rates): 7 counts per
-- interval and the coefficients one after another.
ckAddChebyshev :: CkWriter -> CkHeader -> [Double] -> [Double] -> [Int] -> [Double] -> Bool -> IO (Result ())
ckAddChebyshev (CkWriter w) header starts stops counts coefficients hasRates
  | length stops /= n || length counts /= 7 * n = pure (failure Domain)
  | otherwise = withForeignPtr w $ \pw -> withHeader header $ \ph -> withDoubles starts $ \pa -> withDoubles stops $ \pb -> withArray (map fromIntegral counts) $ \pc -> withDoubles coefficients $ \pk ->
      check <$> C.myb_ck_writer_add_chebyshev pw ph pa pb pc pk (toSize n) (if hasRates then 1 else 0)
  where
    n = length starts

-- | Packets for interpolation in sliding windows (types 5 and 6): the subtype (0 .. 3), the
-- window, the seconds a tick lasts, the ticks and the packets (of the subtype's size each).
data CkPackets = CkPackets {packetsSubtype :: !Int, packetsWindow :: !Int, packetsSecondsPerTick :: !Double, packetsTicks :: ![Double], packetsData :: ![Double]}
  deriving (Eq, Show)

withPackets :: CkPackets -> (C.Myb_ck_packets -> IO r) -> IO r
withPackets p k = withDoubles (packetsTicks p) $ \pt -> withDoubles (packetsData p) $ \pd ->
  k
    C.zero
      { C.myb_ck_packets'subtype = fromIntegral (packetsSubtype p)
      , C.myb_ck_packets'window = fromIntegral (packetsWindow p)
      , C.myb_ck_packets'seconds_per_tick = coerce (packetsSecondsPerTick p)
      , C.myb_ck_packets'ticks = pt
      , C.myb_ck_packets'packets = pd
      , C.myb_ck_packets'count = toSize (length (packetsTicks p))
      }

-- | A segment of type 5 over [start, stop] with the starts of its intervals ('Nothing': one).
ckAddPackets :: CkWriter -> CkHeader -> CkPackets -> Bool -> Double -> Double -> Maybe [Double] -> IO (Result ())
ckAddPackets (CkWriter w) header packets hasRates start stop starts =
  withForeignPtr w $ \pw -> withHeader header $ \ph -> withPackets packets $ \made -> with made $ \pp -> withOptional withDoubles starts $ \ps ->
    check <$> C.myb_ck_writer_add_packets pw ph pp (if hasRates then 1 else 0) (coerce start) (coerce stop) ps (toSize (maybe 0 length starts))

-- | A segment of type 6 over [start, stop]: mini-segments with their boundaries (one more).
ckAddMinisegments :: CkWriter -> CkHeader -> [CkPackets] -> [Double] -> Bool -> Bool -> Double -> Double -> IO (Result ())
ckAddMinisegments (CkWriter w) header minisegments boundaries hasRates selectLast start stop
  | length boundaries /= length minisegments + 1 = pure (failure Domain)
  | otherwise = withForeignPtr w $ \pw -> withHeader header $ \ph -> withMany withPackets minisegments $ \made -> withArray made $ \pm -> withDoubles boundaries $ \pb ->
      check <$> C.myb_ck_writer_add_minisegments pw ph pm pb (toSize (length minisegments)) (if hasRates then 1 else 0) (if selectLast then 1 else 0) (coerce start) (coerce stop)

ckWriterSize :: CkWriter -> IO Int
ckWriterSize (CkWriter w) = fromSize <$> withForeignPtr w C.myb_ck_writer_size

-- | The file in memory.
ckWriterBytes :: CkWriter -> IO (Result B.ByteString)
ckWriterBytes writer@(CkWriter w) = do
  size <- ckWriterSize writer
  withForeignPtr w $ \pw -> allocaBytes (max 1 size) $ \buffer -> do
    code <- C.myb_ck_writer_write pw buffer (toSize size :: CSize)
    andThen code (Right <$> B.packCStringLen (castPtr buffer, size))

-- | The file at a path, created or replaced.
ckWriterFile :: CkWriter -> FilePath -> IO (Result ())
ckWriterFile (CkWriter w) path = withText path $ \p -> withForeignPtr w $ \pw -> check <$> C.myb_ck_writer_write_file pw p
