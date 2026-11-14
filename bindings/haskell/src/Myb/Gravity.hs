-- SPDX-License-Identifier: BSD-3-Clause
-- | The field of a body in spherical harmonics (myb/grav.h).
module Myb.Gravity
  ( GravityField
  , newGravityField
  , loadIcgemFile
  , fieldDegree
  , fieldGm
  , fieldRadius
  , gravityFieldPtr
  ) where

import Data.Coerce (coerce)
import Foreign.C.Types (CDouble (..), CInt)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Marshal.Array (withArray)
import Foreign.Ptr (nullPtr)
import qualified Myb.C as C
import Myb.Error
import Myb.Internal

newtype GravityField = GravityField (ForeignPtr C.Myb_gravity_field)

-- | A field from its coefficients: gm (m^3/s^2), the reference radius (m), the degree, and the
-- fully normalised C and S at the index n (n + 1) / 2 + m, to the degree.
newGravityField :: Double -> Double -> Int -> [Double] -> [Double] -> IO (Result GravityField)
newGravityField gm radius degree c s
  | degree < 0 || degree > fromIntegral (maxBound :: CInt) = pure (failure Domain)
  | length c < needed || length s < needed = pure (failure Domain)
  | otherwise =
      withArray (coerce c :: [CDouble]) $ \cs -> withArray (coerce s :: [CDouble]) $ \ss ->
        create
          (\out -> check <$> C.myb_gravity_field_create (coerce gm) (coerce radius) (fromIntegral degree) cs ss nullPtr out)
          (fmap GravityField . ownedLeaf C.p'myb_gravity_field_destroy)
  where
    needed = (degree + 1) * (degree + 2) `div` 2

-- | A field of an ICGEM file (gfc), to a degree; 'Nothing': all of it.
loadIcgemFile :: FilePath -> Maybe Int -> IO (Result GravityField)
loadIcgemFile path maxDegree =
  withText path $ \cpath ->
    create
      (\out -> withDiag (C.myb_gravity_field_load_icgem_file cpath degree nullPtr out))
      (fmap GravityField . ownedLeaf C.p'myb_gravity_field_destroy)
  where
    degree = maybe (-1) (fromIntegral . min (fromIntegral (maxBound :: CInt)) . max 0) maxDegree

fieldDegree :: GravityField -> IO Int
fieldDegree (GravityField field) = fromIntegral <$> withForeignPtr field C.myb_gravity_field_degree

-- | The gravitational parameter of the field (m^3/s^2).
fieldGm :: GravityField -> IO Double
fieldGm (GravityField field) = coerce <$> withForeignPtr field C.myb_gravity_field_gm

-- | The reference radius of the field (m).
fieldRadius :: GravityField -> IO Double
fieldRadius (GravityField field) = coerce <$> withForeignPtr field C.myb_gravity_field_radius

gravityFieldPtr :: GravityField -> ForeignPtr C.Myb_gravity_field
gravityFieldPtr (GravityField field) = field
