-- SPDX-License-Identifier: BSD-3-Clause
-- | The orbit data messages of the CCSDS (myb/ccsds.h): ephemeris (OEM), parameter (OPM),
-- mean-elements (OMM) and comprehensive (OCM) messages, written and read in either notation.
module Myb.Message
  ( Notation (..)
  , OemSegment (..)
  , writeOem
  , Oem
  , parseOem
  , oemSegmentCount
  , oemStateCount
  , oemMeta
  , oemState
  , oemPtr
  , loadOemFile
  , oemHeader
  , oemAcceleration
  , oemCovarianceCount
  , oemCovariance
    -- * OPM and OMM
  , C.Myb_opm (..)
  , C.Myb_opm_maneuver (..)
  , parseOpm
  , opmText
  , C.Myb_omm (..)
  , parseOmm
  , ommText
    -- * OCM
  , Ocm
  , OcmKind (..)
  , loadOcm
  , loadOcmFile
  , ocmHeader
  , ocmMeta
  , ocmBlockCount
  , ocmValue
  , ocmLineCount
  , ocmLineText
  , ocmLine
  , ocmState
  , OcmBlock (..)
  , OcmData (..)
  , writeOcm
  , ocmFormatLine
  ) where

import Data.Coerce (coerce)
import Foreign.C.String (CString)
import Foreign.C.Types (CDouble (..))
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (allocaArray, withArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import Foreign.Storable (peek, poke)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Types

-- | The notation of a message.
data Notation
  = -- | KEYWORD = value lines.
    Kvn
  | Xml
  deriving (Eq, Show)

-- | A segment of an orbit ephemeris message: states in m and m/s at calendar labels.
data OemSegment = OemSegment
  { objectName :: String
  , objectId :: String
  , centerName :: String
  , refFrame :: String
  , timeSystem :: String
  , -- | The interpolation and its degree, if the message is to name them.
    interpolation :: Maybe (String, Int)
  , epochs :: [Calendar]
  , states :: [State]
  , -- | the covariance matrices (36 by rows: m^2, m^2/s, m^2/s^2) at their epochs
    covariances :: [(Calendar, [Double])]
  , -- | COV_REF_FRAME of each, if the message is to name one
    covarianceFrame :: Maybe String
  }
  deriving (Eq, Show)

withMaybeText :: Maybe String -> (CString -> IO (Result a)) -> IO (Result a)
withMaybeText Nothing k = k nullPtr
withMaybeText (Just s) k = withText s k

-- | The structures of C of the segments, valid inside the continuation.
withSegments :: [OemSegment] -> ([C.Myb_oem_segment_data] -> IO (Result a)) -> IO (Result a)
withSegments [] k = k []
withSegments (s : rest) k =
  withText (objectName s) $ \name ->
    withText (objectId s) $ \identifier ->
      withText (centerName s) $ \center ->
        withText (refFrame s) $ \frame ->
          withText (timeSystem s) $ \system ->
            withMaybeText (fst <$> interpolation s) $ \method ->
              withArray (map calendarToC (epochs s)) $ \labels ->
                withArray (states s) $ \ys ->
                 withArray (map (calendarToC . fst) (covariances s)) $ \covLabels ->
                  withArray (coerce (concatMap snd (covariances s)) :: [CDouble]) $ \covValues ->
                   withMaybeText (covarianceFrame s) $ \covFrame ->
                    withSegments rest $ \others ->
                    k
                      ( C.zero
                          { C.myb_oem_segment_data'object_name = name
                          , C.myb_oem_segment_data'object_id = identifier
                          , C.myb_oem_segment_data'center_name = center
                          , C.myb_oem_segment_data'ref_frame = frame
                          , C.myb_oem_segment_data'time_system = system
                          , C.myb_oem_segment_data'interpolation = method
                          , C.myb_oem_segment_data'interpolation_degree = maybe 0 (fromIntegral . snd) (interpolation s)
                          , C.myb_oem_segment_data'count = toSize (length (epochs s))
                          , C.myb_oem_segment_data'epochs = labels
                          , C.myb_oem_segment_data'states = castPtr ys
                          , C.myb_oem_segment_data'covariance_count = toSize (length (covariances s))
                          , C.myb_oem_segment_data'covariance_epochs = covLabels
                          , C.myb_oem_segment_data'covariance_frame = covFrame
                          , C.myb_oem_segment_data'covariances = covValues
                          }
                          : others
                      )

-- | An orbit ephemeris message as text: the originator, the creation date, the segments.
writeOem :: String -> Calendar -> [OemSegment] -> Notation -> IO (Result String)
writeOem originator created segments notation
  | any (\s -> length (epochs s) /= length (states s) || any ((/= 36) . length . snd) (covariances s)) segments = pure (failure Domain)
  | otherwise =
      withText originator $ \by -> withSegments segments $ \made -> withArray made $ \list ->
        with
          ( C.zero
              { C.myb_oem_data'originator = by
              , C.myb_oem_data'creation_date = calendarToC created
              , C.myb_oem_data'segment_count = toSize (length made)
              , C.myb_oem_data'segments = list
              }
          )
          $ \d -> writeText (case notation of Kvn -> C.myb_oem_write d; Xml -> C.myb_oem_write_as d C.MYB_ODM_XML)

-- | An orbit ephemeris message read from text.
newtype Oem = Oem (ForeignPtr C.Myb_oem)

-- | Reads a message in either notation; a malformed one fails with the line of the error.
parseOem :: String -> IO (Result Oem)
parseOem text =
  withTextLen text $ \p n ->
    create (\out -> withDiag (C.myb_oem_load p n nullPtr out)) (fmap Oem . ownedLeaf C.p'myb_oem_destroy)

oemSegmentCount :: Oem -> IO Int
oemSegmentCount (Oem oem) = fromSize <$> withForeignPtr oem C.myb_oem_segment_count

-- | The number of states of a segment (0 if there is no such segment).
oemStateCount :: Oem -> Int -> IO Int
oemStateCount m@(Oem oem) segment = do
  count <- oemSegmentCount m
  if segment < 0 || segment >= count
    then pure 0
    else fromSize <$> withForeignPtr oem (\o -> C.myb_oem_state_count o (toSize segment))

-- | A keyword of the metadata of a segment, as written; 'Nothing' if the message lacks it.
oemMeta :: Oem -> Int -> String -> IO (Maybe String)
oemMeta m@(Oem oem) segment key = do
  count <- oemSegmentCount m
  if segment < 0 || segment >= count
    then pure Nothing
    else do
      found <- withText key $ \ckey -> withForeignPtr oem $ \o -> do
        value <- C.myb_oem_meta o (toSize segment) ckey
        if value == nullPtr then pure (Right Nothing) else Right . Just <$> peekText value
      pure (either (const Nothing) id found)

-- | The epoch and the state of a line of a segment; 'Nothing' if there is no such line.
oemState :: Oem -> Int -> Int -> IO (Maybe (Calendar, State))
oemState m@(Oem oem) segment index = do
  count <- oemStateCount m segment
  if index < 0 || index >= count
    then pure Nothing
    else withForeignPtr oem $ \o -> with (C.zero :: C.Myb_calendar) $ \label -> with (State (Vec3 0 0 0) (Vec3 0 0 0)) $ \y -> do
      _ <- C.myb_oem_state o (toSize segment) (toSize index) label (castPtr y) nullPtr
      (\l s -> Just (calendarFromC l, s)) <$> peek label <*> peek y

oemPtr :: Oem -> ForeignPtr C.Myb_oem
oemPtr (Oem oem) = oem

-- | A message read from a file.
loadOemFile :: FilePath -> IO (Result Oem)
loadOemFile path = withText path $ \p ->
  create (\out -> withDiag (C.myb_oem_load_file p nullPtr out)) (fmap Oem . ownedLeaf C.p'myb_oem_destroy)

-- | A keyword of the header, as written.
oemHeader :: Oem -> String -> IO (Maybe String)
oemHeader (Oem oem) key = do
  found <- withText key $ \ckey -> withForeignPtr oem $ \o -> do
    value <- C.myb_oem_header o ckey
    if value == nullPtr then pure (Right Nothing) else Right . Just <$> peekText value
  pure (either (const Nothing) id found)

-- | The acceleration of a line of a segment, if the line gives one.
oemAcceleration :: Oem -> Int -> Int -> IO (Maybe Vec3)
oemAcceleration m@(Oem oem) segment index = do
  count <- oemStateCount m segment
  if index < 0 || index >= count
    then pure Nothing
    else withForeignPtr oem $ \o -> with (Vec3 0 0 0) $ \a -> do
      flags <- C.myb_oem_state o (toSize segment) (toSize index) nullPtr nullPtr (castPtr a)
      if flags == 3 then Just <$> peek a else pure Nothing

-- | The number of covariance matrices of a segment (0 if there is no such segment).
oemCovarianceCount :: Oem -> Int -> IO Int
oemCovarianceCount m@(Oem oem) segment = do
  count <- oemSegmentCount m
  if segment < 0 || segment >= count then pure 0 else fromSize <$> withForeignPtr oem (\o -> C.myb_oem_covariance_count o (toSize segment))

-- | A covariance matrix of a segment: its epoch, its frame (if the message names one) and the
-- 36 numbers by rows.
oemCovariance :: Oem -> Int -> Int -> IO (Maybe (Calendar, Maybe String, [Double]))
oemCovariance m@(Oem oem) segment index = do
  count <- oemCovarianceCount m segment
  if index < 0 || index >= count
    then pure Nothing
    else withForeignPtr oem $ \o -> with (C.zero :: C.Myb_calendar) $ \label -> with nullPtr $ \frame -> allocaArray 36 $ \c -> do
      C.myb_oem_covariance o (toSize segment) (toSize index) label frame c
      f <- peek frame
      name <- if f == nullPtr then pure Nothing else Just <$> peekText f
      (\l v -> Just (calendarFromC l, name, v)) <$> peek label <*> peekDoubles 36 c

-- | An orbit parameter message read from text in either notation.
parseOpm :: String -> IO (Result C.Myb_opm)
parseOpm text = withTextLen text $ \p n -> with C.zero $ \o -> do
  done <- withDiag (C.myb_opm_parse p n o)
  either (pure . Left) (const (Right <$> peek o)) done

-- | An orbit parameter message as text.
opmText :: C.Myb_opm -> Notation -> IO (Result String)
opmText opm notation = with opm $ \o -> writeText (case notation of Kvn -> C.myb_opm_write o; Xml -> C.myb_opm_write_as o C.MYB_ODM_XML)

-- | An orbit mean-elements message read from text in either notation.
parseOmm :: String -> IO (Result C.Myb_omm)
parseOmm text = withTextLen text $ \p n -> with C.zero $ \o -> do
  done <- withDiag (C.myb_omm_parse p n o)
  either (pure . Left) (const (Right <$> peek o)) done

-- | An orbit mean-elements message as text.
ommText :: C.Myb_omm -> Notation -> IO (Result String)
ommText omm notation = with omm $ \o -> writeText (case notation of Kvn -> C.myb_omm_write o; Xml -> C.myb_omm_write_as o C.MYB_ODM_XML)

-- ---- OCM -----------------------------------------------------------------------------------

-- | An orbit comprehensive message read from text.
newtype Ocm = Ocm (ForeignPtr C.Myb_ocm)

-- | The kinds of the blocks of an OCM.
data OcmKind = OcmTraj | OcmPhys | OcmCov | OcmMan | OcmPert | OcmOd | OcmUser
  deriving (Eq, Show)

kindCode :: OcmKind -> C.Myb_ocm_kind
kindCode kind = case kind of
  OcmTraj -> C.MYB_OCM_TRAJ
  OcmPhys -> C.MYB_OCM_PHYS
  OcmCov -> C.MYB_OCM_COV
  OcmMan -> C.MYB_OCM_MAN
  OcmPert -> C.MYB_OCM_PERT
  OcmOd -> C.MYB_OCM_OD
  OcmUser -> C.MYB_OCM_USER

loadOcm :: String -> IO (Result Ocm)
loadOcm text = withTextLen text $ \p n -> create (\out -> withDiag (C.myb_ocm_load p n nullPtr out)) (fmap Ocm . ownedLeaf C.p'myb_ocm_destroy)

loadOcmFile :: FilePath -> IO (Result Ocm)
loadOcmFile path = withText path $ \p -> create (\out -> withDiag (C.myb_ocm_load_file p nullPtr out)) (fmap Ocm . ownedLeaf C.p'myb_ocm_destroy)

optionalText :: IO CString -> IO (Maybe String)
optionalText get = get >>= \p -> if p == nullPtr then pure Nothing else Just <$> peekText p

-- | A keyword of the header.
ocmHeader :: Ocm -> String -> IO (Maybe String)
ocmHeader (Ocm m) key = either (const Nothing) id <$> withText key (\k -> withForeignPtr m $ \o -> Right <$> optionalText (C.myb_ocm_header o k))

-- | A keyword of the metadata.
ocmMeta :: Ocm -> String -> IO (Maybe String)
ocmMeta (Ocm m) key = either (const Nothing) id <$> withText key (\k -> withForeignPtr m $ \o -> Right <$> optionalText (C.myb_ocm_meta o k))

-- | The blocks of a kind.
ocmBlockCount :: Ocm -> OcmKind -> IO Int
ocmBlockCount (Ocm m) kind = fromSize <$> withForeignPtr m (\o -> C.myb_ocm_block_count o (kindCode kind))

blockOk :: Ocm -> OcmKind -> Int -> IO Bool
blockOk ocm kind index = (\n -> index >= 0 && index < n) <$> ocmBlockCount ocm kind

-- | A keyword of a block of a kind.
ocmValue :: Ocm -> OcmKind -> Int -> String -> IO (Maybe String)
ocmValue ocm@(Ocm m) kind index key = do
  ok <- blockOk ocm kind index
  if not ok then pure Nothing else either (const Nothing) id <$> withText key (\k -> withForeignPtr m $ \o -> Right <$> optionalText (C.myb_ocm_value o (kindCode kind) (toSize index) k))

-- | The lines of data of a block.
ocmLineCount :: Ocm -> OcmKind -> Int -> IO Int
ocmLineCount ocm@(Ocm m) kind index = do
  ok <- blockOk ocm kind index
  if not ok then pure 0 else fromSize <$> withForeignPtr m (\o -> C.myb_ocm_line_count o (kindCode kind) (toSize index))

-- | A line of data as text.
ocmLineText :: Ocm -> OcmKind -> Int -> Int -> IO (Maybe String)
ocmLineText ocm@(Ocm m) kind index line = do
  n <- ocmLineCount ocm kind index
  if line < 0 || line >= n then pure Nothing else withForeignPtr m $ \o -> optionalText (C.myb_ocm_line_text o (kindCode kind) (toSize index) (toSize line))

-- | A line of data as numbers: its time (a calendar label, or seconds from EPOCH_TZERO) and its
-- values in the units of the message.
ocmLine :: Ocm -> OcmKind -> Int -> Int -> IO (Result (Either Double Calendar, [Double]))
ocmLine ocm@(Ocm m) kind index line = do
  n <- ocmLineCount ocm kind index
  if line < 0 || line >= n
    then pure (failure Domain)
    else withForeignPtr m $ \o -> do
      counted <- lineOf o 0 nullPtr
      case counted of
        Left e -> pure (Left e)
        Right (_, count) -> allocaArray (max 1 count) $ \values -> do
          full <- lineOf o count values
          case full of
            Left e -> pure (Left e)
            Right (time, _) -> Right . (,) time <$> peekDoubles count values
  where
    lineOf :: Ptr C.Myb_ocm -> Int -> Ptr CDouble -> IO (Result (Either Double Calendar, Int))
    lineOf o capacity values = alloca $ \absolute -> with (C.zero :: C.Myb_calendar) $ \epoch -> alloca $ \offset -> alloca $ \count -> do
      poke absolute 0
      poke offset 0
      poke count 0
      code <- C.myb_ocm_line o (kindCode kind) (toSize index) (toSize line) absolute epoch offset values (toSize capacity) count
      andThen code $ do
        a <- peek absolute
        time <- if a /= 0 then Right . calendarFromC <$> peek epoch else Left . coerce <$> peek offset
        k <- fromSize <$> peek count
        pure (Right (time, k))

-- | A line of a trajectory in m, m/s, m/s^2: its time and 3, 6 or 9 values.
ocmState :: Ocm -> Int -> Int -> IO (Result (Either Double Calendar, [Double]))
ocmState (Ocm m) index line = withForeignPtr m $ \o -> alloca $ \absolute -> with (C.zero :: C.Myb_calendar) $ \epoch -> alloca $ \offset -> allocaArray 9 $ \state -> alloca $ \count -> do
  poke absolute 0
  poke offset 0
  poke count 0
  code <- C.myb_ocm_state o (toSize index) (toSize line) absolute epoch offset state count
  andThen code $ do
    a <- peek absolute
    time <- if a /= 0 then Right . calendarFromC <$> peek epoch else Left . coerce <$> peek offset
    k <- fromIntegral <$> peek count
    Right . (,) time <$> peekDoubles k state

-- | A block of an OCM to write: its kind, its keywords and values, its lines of data.
data OcmBlock = OcmBlock {blockKind :: OcmKind, blockPairs :: [(String, String)], blockLines :: [String]}
  deriving (Eq, Show)

-- | What an OCM is written from: the originator, the creation date, the metadata (with
-- TIME_SYSTEM and EPOCH_TZERO) and the blocks.
data OcmData = OcmData {ocmOriginator :: String, ocmCreated :: Calendar, ocmMetadata :: [(String, String)], ocmBlocks :: [OcmBlock]}
  deriving (Eq, Show)

withPairs :: [(String, String)] -> (Ptr C.Myb_ocm_pair -> Int -> IO (Result a)) -> IO (Result a)
withPairs pairs k =
  withMany (\(key, value) go -> withText key $ \pk -> withText value $ \pv -> go (C.zero {C.myb_ocm_pair'key = pk, C.myb_ocm_pair'value = pv})) pairs $ \made ->
    withArrayLen made $ \n p -> k p n

-- | An OCM (version 3.0) as text.
writeOcm :: OcmData -> Notation -> IO (Result String)
writeOcm d notation =
  withText (ocmOriginator d) $ \by -> withPairs (ocmMetadata d) $ \pm nm ->
    withMany block (ocmBlocks d) $ \blocks -> withArrayLen blocks $ \nb pb ->
      with
        C.zero
          { C.myb_ocm_data'originator = by
          , C.myb_ocm_data'creation_date = calendarToC (ocmCreated d)
          , C.myb_ocm_data'meta_count = toSize nm
          , C.myb_ocm_data'meta = pm
          , C.myb_ocm_data'block_count = toSize nb
          , C.myb_ocm_data'blocks = pb
          }
        (\pd -> writeText (C.myb_ocm_write pd (case notation of Kvn -> C.MYB_ODM_KVN; Xml -> C.MYB_ODM_XML)))
  where
    block b go = withPairs (blockPairs b) $ \pp np -> withMany withText (blockLines b) $ \ls -> withArrayLen ls $ \nl pl ->
      go
        C.zero
          { C.myb_ocm_block_data'kind = kindCode (blockKind b)
          , C.myb_ocm_block_data'pair_count = toSize np
          , C.myb_ocm_block_data'pairs = pp
          , C.myb_ocm_block_data'line_count = toSize nl
          , C.myb_ocm_block_data'lines = pl
          }

-- | A line of data: the time (a calendar label, or seconds) and the values, each the shortest
-- decimal that reads back.
ocmFormatLine :: Either Double Calendar -> [Double] -> IO (Result String)
ocmFormatLine time values = withDoublesLen values $ \n pv -> case time of
  Left offset -> writeText (C.myb_ocm_format_line nullPtr (coerce offset) pv (toSize n))
  Right label -> with (calendarToC label) $ \pl -> writeText (C.myb_ocm_format_line pl 0 pv (toSize n))
