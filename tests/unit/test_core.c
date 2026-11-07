/* SPDX-License-Identifier: BSD-3-Clause */
#include <string.h>

#include "core/base.h"
#include "mtest.h"
#include "myb/core.h"

static void test_version(void) {
  MT_CHECK(myb_version_number() == MYB_VERSION_NUMBER);
}

static void test_status_names(void) {
  static const struct {
    myb_status status;
    const char* name;
  } expected[] = {
      {MYB_OK, "MYB_OK"},
      {MYB_E_DOMAIN, "MYB_E_DOMAIN"},
      {MYB_E_RANGE, "MYB_E_RANGE"},
      {MYB_E_NOCONV, "MYB_E_NOCONV"},
      {MYB_E_STEP_MIN, "MYB_E_STEP_MIN"},
      {MYB_E_NOMEM, "MYB_E_NOMEM"},
      {MYB_E_OVERFLOW, "MYB_E_OVERFLOW"},
      {MYB_E_IO, "MYB_E_IO"},
      {MYB_E_FORMAT, "MYB_E_FORMAT"},
  };

  for (size_t i = 0; i < sizeof expected / sizeof expected[0]; ++i) {
    MT_CHECK(strcmp(myb_status_name(expected[i].status), expected[i].name) == 0);
    MT_CHECK((int)expected[i].status == (int)i); /* values are stable and contiguous */
  }
  MT_CHECK(strcmp(myb_status_name((myb_status)42), "MYB_E_UNKNOWN") == 0);
}

int main(void) {
  test_version();
  test_status_names();
  return mt_finish("test_core");
}
