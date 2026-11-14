-- SPDX-License-Identifier: BSD-3-Clause
-- | Orbit ephemeris messages of the CCSDS (myb/ccsds.h), written and read in either notation.
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
  ) where

import Foreign.C.String (CString)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Array (withArray)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (castPtr, nullPtr)
import Foreign.Storable (peek)
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
                          }
                          : others
                      )

-- | An orbit ephemeris message as text: the originator, the creation date, the segments.
writeOem :: String -> Calendar -> [OemSegment] -> Notation -> IO (Result String)
writeOem originator created segments notation
  | any (\s -> length (epochs s) /= length (states s)) segments = pure (failure Domain)
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
          $ \d -> writeText (C.myb_oem_write_as d side)
  where
    side = case notation of
      Kvn -> C.MYB_ODM_KVN
      Xml -> C.MYB_ODM_XML

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
