-- SPDX-License-Identifier: BSD-3-Clause
-- | Files (myb/io.h, myb/kin.h): numbers read and printed exactly, tables of numbers (CSV and
-- the like), the writer of SPK files (types 2, 3, 9, 13 and raw segments) and SP3 precise
-- orbits, read and written.
module Myb.Files
  ( -- * Numbers
    parseDouble
  , formatDouble
  , formatFixed
  , formatScientific
    -- * Tables
  , Table
  , loadTable
  , loadTableFile
  , tableRows
  , tableColumns
  , tableName
  , tableRow
  , writeTable
    -- * Writing SPK files
  , SpkWriter
  , newSpkWriter
  , SegmentHeader (..)
  , spkAddStates
  , spkAddChebyshev
  , spkAddRaw
  , spkWriterSize
  , spkWriterBytes
  , spkWriterFile
    -- * SP3 precise orbits
  , Sp3
  , loadSp3
  , loadSp3File
  , sp3Satellites
  , sp3EpochCount
  , sp3TimeSystem
  , sp3Frame
  , sp3Epoch
  , sp3State
  , sp3Clock
  , sp3Sigma
  , sp3Accuracy
  , sp3Flags
  , Sp3Data (..)
  , writeSp3
  ) where

import qualified Data.ByteString as B
import Data.Coerce (coerce)
import Data.Int (Int32)
import Data.Word (Word8)
import Foreign.C.String (CString)
import Foreign.C.Types (CDouble (..), CSize)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca, allocaBytes)
import Foreign.Marshal.Array (allocaArray, peekArray, withArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- ---- Numbers -------------------------------------------------------------------------------

-- | The number at the start of the text, correctly rounded, and the characters it used.
parseDouble :: String -> Result (Double, Int)
parseDouble text = pureCall $ withTextLen text $ \p n -> alloca $ \value -> alloca $ \used -> do
  code <- C.myb_parse_double p n value used
  andThen code (Right <$> ((,) <$> (coerce <$> peek value) <*> (fromSize <$> peek used)))

-- | The shortest decimal that reads back to the same value.
formatDouble :: Double -> String
formatDouble x = pureCall $ allocaBytes C.MYB_DOUBLE_CHARS $ \text -> do
  n <- C.myb_format_double (coerce x) text
  peekTextLen text (fromSize n)

peekTextLen :: CString -> Int -> IO String
peekTextLen p n = map (toEnum . fromIntegral) <$> (peekArray n (castPtr p) :: IO [Word8])

-- | `decimals` (0 .. 60) digits after the point, as "%.Nf".
formatFixed :: Double -> Int -> Result String
formatFixed x decimals = pureCall $ writeText (C.myb_format_fixed (coerce x) (fromIntegral decimals))

-- | `digits` (1 .. 60) significant digits with an exponent.
formatScientific :: Double -> Int -> Result String
formatScientific x digits = pureCall $ writeText (C.myb_format_scientific (coerce x) (fromIntegral digits))

-- ---- Tables --------------------------------------------------------------------------------

newtype Table = Table (ForeignPtr C.Myb_table)

-- | A table from text; a malformed one fails with the line of the error.
loadTable :: String -> IO (Result Table)
loadTable text = withTextLen text $ \p n ->
  create (\o -> withDiag (C.myb_table_load p n nullPtr o)) (fmap Table . ownedLeaf C.p'myb_table_destroy)

loadTableFile :: FilePath -> IO (Result Table)
loadTableFile path = withText path $ \p ->
  create (\o -> withDiag (C.myb_table_load_file p nullPtr o)) (fmap Table . ownedLeaf C.p'myb_table_destroy)

tableRows :: Table -> IO Int
tableRows (Table t) = fromSize <$> withForeignPtr t C.myb_table_rows

tableColumns :: Table -> IO Int
tableColumns (Table t) = fromSize <$> withForeignPtr t C.myb_table_columns

-- | The name of a column ("" without a row of names); 'Nothing' for no such column.
tableName :: Table -> Int -> IO (Maybe String)
tableName table@(Table t) column = do
  n <- tableColumns table
  if column < 0 || column >= n then pure Nothing else withForeignPtr t $ \p -> Just <$> (C.myb_table_name p (toSize column) >>= peekText)

-- | The values of a row; 'Nothing' for no such row.
tableRow :: Table -> Int -> IO (Maybe [Double])
tableRow table@(Table t) row = do
  rows <- tableRows table
  columns <- tableColumns table
  if row < 0 || row >= rows then pure Nothing else withForeignPtr t $ \p -> Just <$> (C.myb_table_row p (toSize row) >>= peekDoubles columns)

-- | Values of the columns given, row by row, as comma-separated text, with a first row of
-- names if given.
writeTable :: Maybe [String] -> Int -> [Double] -> IO (Result String)
writeTable names columns values
  | columns < 1 || length values `mod` columns /= 0 || maybe False ((/= columns) . length) names = pure (failure Domain)
  | otherwise = withNames names $ \pn -> withDoubles values $ \pv ->
      writeText (C.myb_table_write pn (toSize columns) pv (toSize (length values `div` columns)))
  where
    withNames Nothing k = k nullPtr
    withNames (Just ns) k = withMany withText ns (\ps -> withArray ps k)

-- ---- Writing SPK files ---------------------------------------------------------------------

newtype SpkWriter = SpkWriter (ForeignPtr C.Myb_spk_writer)

-- | A writer with the label of the file ('Nothing': "myb").
newSpkWriter :: Maybe String -> IO (Result SpkWriter)
newSpkWriter name = withName $ \p -> create (\o -> check <$> C.myb_spk_writer_create p nullPtr o) (fmap SpkWriter . ownedLeaf C.p'myb_spk_writer_destroy)
  where
    withName k = maybe (k nullPtr) (`withText` k) name

-- | What a segment describes: target relative to center in the frame (NAIF codes), its name.
data SegmentHeader = SegmentHeader {headerTarget :: !Int32, headerCenter :: !Int32, headerFrame :: !Int32, headerName :: !String}
  deriving (Eq, Show)

withHeader :: SegmentHeader -> (Ptr C.Myb_spk_segment_header -> IO (Result a)) -> IO (Result a)
withHeader h k = withText (headerName h) $ \name ->
  with C.zero {C.myb_spk_segment_header'target = headerTarget h, C.myb_spk_segment_header'center = headerCenter h, C.myb_spk_segment_header'frame = headerFrame h, C.myb_spk_segment_header'name = name} k

-- | A segment of states at increasing epochs (s past J2000, TDB): type 9 (Lagrange) or 13
-- (Hermite) of the degree given.
spkAddStates :: SpkWriter -> SegmentHeader -> Int -> Int -> [Double] -> [State] -> IO (Result ())
spkAddStates (SpkWriter w) header kind degree epochs states
  | length epochs /= length states = pure (failure Domain)
  | otherwise = withForeignPtr w $ \pw -> withHeader header $ \ph -> withDoubles epochs $ \pe -> withVecs states $ \ps ->
      check <$> C.myb_spk_writer_add_states pw ph (fromIntegral kind) (fromIntegral degree) pe ps (toSize (length epochs))

-- | A segment of Chebyshev polynomials (type 2 or 3) of the degree on intervals of [start,
-- stop], fitted to a provider; with True also the largest distance (m) between the provider
-- and the polynomials half-way between the nodes.
spkAddChebyshev :: SpkWriter -> SegmentHeader -> Int -> StateSource -> Double -> Double -> Int -> Int -> Bool -> IO (Result (Maybe Double))
spkAddChebyshev (SpkWriter w) header kind source start stop intervals degree askError =
  withForeignPtr w $ \pw -> withHeader header $ \ph -> withStateSource source $ \f c -> alloca $ \err -> do
    poke err 0
    code <- C.myb_spk_writer_add_chebyshev pw ph (fromIntegral kind) f c (coerce start) (coerce stop) (toSize intervals) (fromIntegral degree) (if askError then err else nullPtr)
    andThen code (if askError then Right . Just . coerce <$> peek err else pure (Right Nothing))

-- | A segment of any type from its words, written as given.
spkAddRaw :: SpkWriter -> SegmentHeader -> Int -> Double -> Double -> [Double] -> IO (Result ())
spkAddRaw (SpkWriter w) header kind start stop words' =
  withForeignPtr w $ \pw -> withHeader header $ \ph -> withDoublesLen words' $ \n pd ->
    check <$> C.myb_spk_writer_add_raw pw ph (fromIntegral kind) (coerce start) (coerce stop) pd (toSize n)

-- | Bytes of the file with the segments added so far.
spkWriterSize :: SpkWriter -> IO Int
spkWriterSize (SpkWriter w) = fromSize <$> withForeignPtr w C.myb_spk_writer_size

-- | The file in memory.
spkWriterBytes :: SpkWriter -> IO (Result B.ByteString)
spkWriterBytes writer@(SpkWriter w) = do
  size <- spkWriterSize writer
  withForeignPtr w $ \pw -> fileImage size (C.myb_spk_writer_write pw)

-- | The file at a path, created or replaced.
spkWriterFile :: SpkWriter -> FilePath -> IO (Result ())
spkWriterFile (SpkWriter w) path = withText path $ \p -> withForeignPtr w $ \pw -> check <$> C.myb_spk_writer_write_file pw p

-- | The bytes a writer of C lays out into a buffer of the size given.
fileImage :: Int -> (Ptr () -> CSize -> IO C.Myb_status) -> IO (Result B.ByteString)
fileImage size write = allocaBytes (max 1 size) $ \buffer -> do
  code <- write buffer (toSize size)
  andThen code (Right <$> B.packCStringLen (castPtr buffer, size))

-- ---- SP3 -----------------------------------------------------------------------------------

newtype Sp3 = Sp3 (ForeignPtr C.Myb_sp3)

loadSp3 :: String -> IO (Result Sp3)
loadSp3 text = withTextLen text $ \p n -> create (\o -> withDiag (C.myb_sp3_load p n nullPtr o)) (fmap Sp3 . ownedLeaf C.p'myb_sp3_destroy)

loadSp3File :: FilePath -> IO (Result Sp3)
loadSp3File path = withText path $ \p -> create (\o -> withDiag (C.myb_sp3_load_file p nullPtr o)) (fmap Sp3 . ownedLeaf C.p'myb_sp3_destroy)

-- | The identifiers of the satellites ("L51", ...).
sp3Satellites :: Sp3 -> IO [String]
sp3Satellites (Sp3 s) = withForeignPtr s $ \p -> do
  n <- C.myb_sp3_satellite_count p
  if n == 0 then pure [] else mapM (\i -> C.myb_sp3_satellite p i >>= peekText) [0 .. n - 1]

sp3EpochCount :: Sp3 -> IO Int
sp3EpochCount (Sp3 s) = fromSize <$> withForeignPtr s C.myb_sp3_epoch_count

sp3TimeSystem :: Sp3 -> IO String
sp3TimeSystem (Sp3 s) = withForeignPtr s (\p -> C.myb_sp3_time_system p >>= peekText)

sp3Frame :: Sp3 -> IO String
sp3Frame (Sp3 s) = withForeignPtr s (\p -> C.myb_sp3_frame p >>= peekText)

-- | Whether the satellite and the epoch exist.
inRange :: Sp3 -> Int -> Int -> IO Bool
inRange sp3@(Sp3 s) satellite epoch = do
  n <- fromSize <$> withForeignPtr s C.myb_sp3_satellite_count
  m <- sp3EpochCount sp3
  pure (satellite >= 0 && satellite < n && epoch >= 0 && epoch < m)

-- | The calendar label of an epoch in the file's time system.
sp3Epoch :: Sp3 -> Int -> IO (Maybe Calendar)
sp3Epoch sp3@(Sp3 s) epoch = do
  m <- sp3EpochCount sp3
  if epoch < 0 || epoch >= m then pure Nothing else withForeignPtr s $ \p -> with C.zero $ \l -> C.myb_sp3_epoch p (toSize epoch) l >> Just . calendarFromC <$> peek l

-- | The position and the velocity of a satellite at an epoch, each if the file gives it.
sp3State :: Sp3 -> Int -> Int -> IO (Maybe (Maybe Vec3, Maybe Vec3))
sp3State sp3@(Sp3 s) satellite epoch = do
  ok <- inRange sp3 satellite epoch
  if not ok then pure Nothing else withForeignPtr s $ \p -> alloca $ \pos -> alloca $ \vel -> do
    flags <- C.myb_sp3_state p (toSize satellite) (toSize epoch) (castPtr pos) (castPtr vel)
    r <- peek pos
    v <- peek vel
    pure (Just (if odd flags then Just r else Nothing, if flags `div` 2 `mod` 2 == 1 then Just v else Nothing))

-- | The clock correction (s) and its rate (s/s) of a satellite at an epoch, each if given.
sp3Clock :: Sp3 -> Int -> Int -> IO (Maybe (Maybe Double, Maybe Double))
sp3Clock sp3@(Sp3 s) satellite epoch = do
  ok <- inRange sp3 satellite epoch
  if not ok then pure Nothing else withForeignPtr s $ \p -> alloca $ \clock -> alloca $ \rate -> do
    poke clock 0
    poke rate 0
    flags <- C.myb_sp3_clock p (toSize satellite) (toSize epoch) clock rate
    c <- coerce <$> peek clock
    r <- coerce <$> peek rate
    pure (Just (if odd flags then Just c else Nothing, if flags `div` 2 `mod` 2 == 1 then Just r else Nothing))

-- | The eight standard deviations the records carry (0 where none).
sp3Sigma :: Sp3 -> Int -> Int -> IO (Maybe [Double])
sp3Sigma sp3@(Sp3 s) satellite epoch = do
  ok <- inRange sp3 satellite epoch
  if not ok then pure Nothing else withForeignPtr s $ \p -> allocaArray 8 $ \o -> C.myb_sp3_sigma p (toSize satellite) (toSize epoch) o >> Just <$> peekDoubles 8 o

-- | The accuracy of the orbit of a satellite the header states (m), 0 if unknown.
sp3Accuracy :: Sp3 -> Int -> IO (Maybe Double)
sp3Accuracy sp3@(Sp3 s) satellite = do
  ok <- inRange sp3 satellite 0
  if not ok then pure Nothing else Just . coerce <$> withForeignPtr s (\p -> C.myb_sp3_accuracy p (toSize satellite))

-- | The flags of the position record (MYB_SP3_*).
sp3Flags :: Sp3 -> Int -> Int -> IO (Maybe Word)
sp3Flags sp3@(Sp3 s) satellite epoch = do
  ok <- inRange sp3 satellite epoch
  if not ok then pure Nothing else Just . fromIntegral <$> withForeignPtr s (\p -> C.myb_sp3_flags p (toSize satellite) (toSize epoch))

-- | What an SP3 file is written from; the records are by epoch, then satellite; the optional
-- parts may be absent.
data Sp3Data = Sp3Data
  { sp3Version :: Char
  , sp3FileType :: Char
  , sp3DataTimeSystem :: String
  , sp3DataFrame :: String
  , sp3OrbitType :: String
  , sp3Agency :: String
  , sp3DataUsed :: String
  , sp3DataSatellites :: [String]
  , sp3DataEpochs :: [Calendar]
  , sp3Interval :: Double
  , sp3Positions :: [Vec3]
  , sp3Velocities :: Maybe [Vec3]
  , sp3Clocks :: Maybe [Double]
  , sp3ClockRates :: Maybe [Double]
  , sp3Accuracies :: Maybe [Double]
  , sp3Sigmas :: Maybe [Double]
  , sp3DataFlags :: Maybe [Word8]
  }

-- | The file as text.
writeSp3 :: Sp3Data -> IO (Result String)
writeSp3 d
  | length (sp3Positions d) /= records = pure (failure Domain)
  | otherwise =
      withText (sp3DataTimeSystem d) $ \ts -> withText (sp3DataFrame d) $ \fr -> withText (sp3OrbitType d) $ \ot -> withText (sp3Agency d) $ \ag -> withText (sp3DataUsed d) $ \du ->
        withMany withText (sp3DataSatellites d) $ \sats -> withArray (sats ++ [nullPtr]) $ \psats ->
          withArrayLen (map calendarToC (sp3DataEpochs d)) $ \ne pe -> withVecs (sp3Positions d) $ \pp ->
            withOptional withVecs (sp3Velocities d) $ \pv -> withOptional withDoubles (sp3Clocks d) $ \pc -> withOptional withDoubles (sp3ClockRates d) $ \pr ->
              withOptional withDoubles (sp3Accuracies d) $ \pa -> withOptional withDoubles (sp3Sigmas d) $ \psg -> withOptional withArray (sp3DataFlags d) $ \pf ->
                with
                  C.zero
                    { C.myb_sp3_data'version = toEnum (fromEnum (sp3Version d))
                    , C.myb_sp3_data'file_type = toEnum (fromEnum (sp3FileType d))
                    , C.myb_sp3_data'time_system = ts
                    , C.myb_sp3_data'frame = fr
                    , C.myb_sp3_data'orbit_type = ot
                    , C.myb_sp3_data'agency = ag
                    , C.myb_sp3_data'data_used = du
                    , C.myb_sp3_data'satellite_count = toSize (length sats)
                    , C.myb_sp3_data'satellites = psats
                    , C.myb_sp3_data'epoch_count = toSize ne
                    , C.myb_sp3_data'epochs = pe
                    , C.myb_sp3_data'interval = coerce (sp3Interval d)
                    , C.myb_sp3_data'positions = pp
                    , C.myb_sp3_data'velocities = pv
                    , C.myb_sp3_data'clocks = pc
                    , C.myb_sp3_data'clock_rates = pr
                    , C.myb_sp3_data'accuracy = pa
                    , C.myb_sp3_data'sigmas = psg
                    , C.myb_sp3_data'flags = castPtr pf
                    }
                  (writeText . C.myb_sp3_write)
  where
    records = length (sp3DataSatellites d) * length (sp3DataEpochs d)
    withOptional :: (a -> (Ptr b -> IO r) -> IO r) -> Maybe a -> (Ptr c -> IO r) -> IO r
    withOptional _ Nothing k = k nullPtr
    withOptional f (Just x) k = f x (k . castPtr)
