-- SPDX-License-Identifier: BSD-3-Clause
-- | Outcomes: a call of the library that can fail returns @'Result' a@, the value or an
-- 'Error' with the status of the library and, from a parser, what was wrong and on which line.
-- Nothing is thrown: a failure is a value, as in C.
module Myb.Error
  ( Status (..)
  , Error (..)
  , Result
  , statusName
  , statusCode
  , statusFromCode
  , errorText
  ) where

import qualified Myb.C as C
import System.IO.Unsafe (unsafeDupablePerformIO)
import qualified GHC.Foreign as Foreign
import GHC.IO.Encoding (utf8)

-- | What the library answers when a call fails.
data Status
  = -- | An argument outside the mathematical domain of the model.
    Domain
  | -- | An epoch or an argument outside the coverage of the data.
    Range
  | -- | The iteration budget was exhausted without convergence.
    NoConvergence
  | -- | The step of an integrator fell below the allowed minimum.
    StepMin
  | -- | An allocation failed.
    NoMemory
  | -- | A size or an index computation overflowed.
    Overflow
  | -- | A file could not be read.
    Io
  | -- | Malformed input data.
    Format
  | -- | A named item does not exist.
    NotFound
  | -- | A status this binding does not know (a newer library), by its code.
    Other Int
  deriving (Eq, Show)

-- | A failed call: the status and, from parsers, the message ("line 2: ...").
data Error = Error {errorStatus :: !Status, errorMessage :: !String}
  deriving (Eq, Show)

type Result a = Either Error a

-- | The status of a code of the C interface; 'Nothing' for MYB_OK.
statusFromCode :: C.Myb_status -> Maybe Status
statusFromCode code = case code of
  C.MYB_OK -> Nothing
  C.MYB_E_DOMAIN -> Just Domain
  C.MYB_E_RANGE -> Just Range
  C.MYB_E_NOCONV -> Just NoConvergence
  C.MYB_E_STEP_MIN -> Just StepMin
  C.MYB_E_NOMEM -> Just NoMemory
  C.MYB_E_OVERFLOW -> Just Overflow
  C.MYB_E_IO -> Just Io
  C.MYB_E_FORMAT -> Just Format
  C.MYB_E_NOT_FOUND -> Just NotFound
  C.Myb_status other -> Just (Other (fromIntegral other))

-- | The code of the C interface.
statusCode :: Status -> C.Myb_status
statusCode status = case status of
  Domain -> C.MYB_E_DOMAIN
  Range -> C.MYB_E_RANGE
  NoConvergence -> C.MYB_E_NOCONV
  StepMin -> C.MYB_E_STEP_MIN
  NoMemory -> C.MYB_E_NOMEM
  Overflow -> C.MYB_E_OVERFLOW
  Io -> C.MYB_E_IO
  Format -> C.MYB_E_FORMAT
  NotFound -> C.MYB_E_NOT_FOUND
  Other code -> C.Myb_status (fromIntegral code)

-- | The name the library gives the status: "MYB_E_DOMAIN", ...
statusName :: Status -> String
statusName status =
  -- myb_status_name returns a static string for any code
  unsafeDupablePerformIO (C.myb_status_name (statusCode status) >>= Foreign.peekCString utf8)

-- | The status by its name, then the message if there is one: "MYB_E_FORMAT: line 2: ...".
errorText :: Error -> String
errorText (Error status message)
  | null message = statusName status
  | otherwise = statusName status ++ ": " ++ message
