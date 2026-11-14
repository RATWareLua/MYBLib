-- SPDX-License-Identifier: BSD-3-Clause
-- | The array forms of myb/frames.h as list functions: a series of rotations or orientations in
-- one call, the same bits as the single calls in a loop. The objects of the frames (the CIP
-- model, the EOP, the equinox model, the WGCCRE model, the PCK set, the graph of frames) have no
-- upper layer yet: these take the object of "Myb.C" in a ForeignPtr (made with its destroy
-- function, newForeignPtr p'myb_cip_destroy ...). A failure names the element that failed.
module Myb.Frames
  ( C.Myb_rotation (..)
  , C.Myb_orientation (..)
  , C.Myb_eop_values (..)
  , rotationApplyPositionMany
  , rotationApplyStateMany
  , cipXysMany
  , eopAtMany
  , cipTerrestrialMany
  , equinoxRotationMany
  , equinoxTemePefMany
  , wgccreRotationMany
  , pckRotationMany
  , orientationMany
  , frameGraphRotationMany
  ) where

import Data.Int (Int32)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Array (allocaArray, peekArray, withArray, withArrayLen)
import Foreign.Ptr (Ptr, castPtr, nullPtr)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal
import Myb.Time (TimeSystem, timeSystemPtr)
import Myb.Types

-- | Each position carried by its rotation (one rotation per position).
rotationApplyPositionMany :: [C.Myb_rotation] -> [Vec3] -> Result [Vec3]
rotationApplyPositionMany rotations positions
  | length rotations /= length positions = failure Domain
  | otherwise = pureCall $ withArray rotations $ \pr -> withVecs positions $ \pin -> allocaArray (max 1 n) $ \o -> do
      C.myb_rotation_apply_position_many pr pin (toSize n) (castPtr o)
      Right <$> peekArray n o
  where
    n = length positions

-- | Each state carried by its rotation, with the rate of the rotation.
rotationApplyStateMany :: [C.Myb_rotation] -> [State] -> Result [State]
rotationApplyStateMany rotations states
  | length rotations /= length states = failure Domain
  | otherwise = pureCall $ withArray rotations $ \pr -> withVecs states $ \pin -> allocaArray (max 1 n) $ \o -> do
      C.myb_rotation_apply_state_many pr pin (toSize n) (castPtr o)
      Right <$> peekArray n o
  where
    n = length states

times :: [Time] -> (Int -> Ptr C.Myb_time -> IO a) -> IO a
times ts k = withArrayLen (map timeToC ts) k

-- | X, Y, s (rad) and their rates at the TT instants.
cipXysMany :: ForeignPtr C.Myb_cip -> [Time] -> IO ([Vec3], [Vec3])
cipXysMany cip tt = withForeignPtr cip $ \c -> times tt $ \n pt -> allocaArray (max 1 n) $ \xys -> allocaArray (max 1 n) $ \rates -> do
  C.myb_cip_xys_many c pt (toSize n) (castPtr xys) (castPtr rates)
  (,) <$> peekArray n xys <*> peekArray n rates

-- | The values of the EOP and UT1 at instants of a scale (by its name).
eopAtMany :: ForeignPtr C.Myb_eop -> TimeSystem -> [Time] -> String -> IO (Result [(C.Myb_eop_values, Time)])
eopAtMany eop system ts scale =
  withText scale $ \name -> withForeignPtr (timeSystemPtr system) $ \s -> do
    found <- fill (C.myb_timesys_scale s name)
    case found of
      Left e -> pure (Left e)
      Right sc -> withForeignPtr eop $ \pe -> times ts $ \n pt -> allocaArray (max 1 n) $ \values -> allocaArray (max 1 n) $ \ut1 -> do
        done <- arrayForm (C.myb_eop_at_many pe s pt (toSize n) sc values ut1)
        either (pure . Left) (const (Right <$> (zip <$> peekArray n values <*> (map timeFromC <$> peekArray n ut1)))) done

-- | The celestial -> terrestrial rotations at the TT and UT1 instants with the values of the
-- EOP ('Nothing': all zero).
cipTerrestrialMany :: ForeignPtr C.Myb_cip -> [Time] -> [Time] -> Maybe [C.Myb_eop_values] -> IO (Result [C.Myb_rotation])
cipTerrestrialMany cip tt ut1 eop
  | length ut1 /= n || maybe False ((/= n) . length) eop = pure (failure Domain)
  | otherwise = withForeignPtr cip $ \c -> times tt $ \_ ptt -> times ut1 $ \_ put -> withEop $ \pe -> allocaArray (max 1 n) $ \o -> do
      C.myb_cip_terrestrial_many c ptt put pe (toSize n) o
      Right <$> peekArray n o
  where
    n = length tt
    withEop k = maybe (k nullPtr) (`withArray` k) eop

-- | The rotations J2000 -> the frame of the equinox model at the TT instants.
equinoxRotationMany :: ForeignPtr C.Myb_equinox -> C.Myb_equinox_frame -> [Time] -> IO [C.Myb_rotation]
equinoxRotationMany model frame tt = withForeignPtr model $ \m -> times tt $ \n pt -> allocaArray (max 1 n) $ \o ->
  C.myb_equinox_rotation_many m frame pt (toSize n) o >> peekArray n o

-- | The rotations TEME -> PEF at the UT1 instants.
equinoxTemePefMany :: ForeignPtr C.Myb_equinox -> [Time] -> IO [C.Myb_rotation]
equinoxTemePefMany model ut1 = withForeignPtr model $ \m -> times ut1 $ \n pt -> allocaArray (max 1 n) $ \o ->
  C.myb_equinox_teme_pef_many m pt (toSize n) o >> peekArray n o

-- | The rotations reference frame -> body-fixed frame of a WGCCRE model at the TDB instants.
wgccreRotationMany :: ForeignPtr C.Myb_wgccre -> [Time] -> IO [C.Myb_rotation]
wgccreRotationMany model ts = withForeignPtr model $ \m -> times ts $ \n pt -> allocaArray (max 1 n) $ \o ->
  C.myb_wgccre_rotation_many m pt (toSize n) o >> peekArray n o

-- | The rotations from the base frame to the frame class of a PCK set, and the bases.
pckRotationMany :: ForeignPtr C.Myb_pck -> Int32 -> [Time] -> IO (Result [(C.Myb_rotation, Int32)])
pckRotationMany pck class' ts = withForeignPtr pck $ \p -> times ts $ \n pt -> allocaArray (max 1 n) $ \o -> allocaArray (max 1 n) $ \bases -> do
  done <- arrayForm (C.myb_pck_rotation_many p class' pt (toSize n) o bases)
  either (pure . Left) (const (Right <$> (zip <$> peekArray n o <*> peekArray n bases))) done

-- | What an orientation provider gives at the instants.
orientationMany :: OrientationSource -> [Time] -> IO (Result [C.Myb_orientation])
orientationMany source ts = withOrientationSource source $ \f c -> times ts $ \n pt -> allocaArray (max 1 n) $ \o -> do
  done <- arrayForm (C.myb_orientation_many f c pt (toSize n) o)
  either (pure . Left) (const (Right <$> peekArray n o)) done

-- | The rotations from one frame of a graph to another at the instants.
frameGraphRotationMany :: ForeignPtr C.Myb_frame_graph -> Int32 -> Int32 -> [Time] -> IO (Result [C.Myb_rotation])
frameGraphRotationMany graph from to ts = withForeignPtr graph $ \g -> times ts $ \n pt -> allocaArray (max 1 n) $ \o -> do
  done <- arrayForm (C.myb_frame_graph_rotation_many g from to pt (toSize n) o)
  either (pure . Left) (const (Right <$> peekArray n o)) done
