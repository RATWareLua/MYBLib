/* SPDX-License-Identifier: BSD-3-Clause */
#include "core/base.h"
#include "myb/core.h"

const char* myb_status_name(myb_status status) {
  switch (status) {
    case MYB_OK:
      return "MYB_OK";
    case MYB_E_DOMAIN:
      return "MYB_E_DOMAIN";
    case MYB_E_RANGE:
      return "MYB_E_RANGE";
    case MYB_E_NOCONV:
      return "MYB_E_NOCONV";
    case MYB_E_STEP_MIN:
      return "MYB_E_STEP_MIN";
    case MYB_E_NOMEM:
      return "MYB_E_NOMEM";
    case MYB_E_OVERFLOW:
      return "MYB_E_OVERFLOW";
    case MYB_E_IO:
      return "MYB_E_IO";
    case MYB_E_FORMAT:
      return "MYB_E_FORMAT";
    case MYB_E_NOT_FOUND:
      return "MYB_E_NOT_FOUND";
  }
  /* A value outside the enumeration (cast from an integer) is not a contract violation of
   * this diagnostic helper: it reports it instead of terminating. */
  return "MYB_E_UNKNOWN";
}

long myb_version_number(void) {
  return MYB_VERSION_NUMBER;
}
