/* SPDX-License-Identifier: BSD-3-Clause */
/* Coefficient tables of the IGRF form (myb_magnetic_field_load_table) at an epoch drawn from
 * the input, the field evaluated at one place. */
#include "fuzz.h"
#include "myb/force.h"
#include "myb/grav.h"

MYB_FUZZ_NO_INITIALIZE

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  const uint64_t hash = fuzz_hash(data, size);
  const double year = 1890.0 + (double)(hash % 1600u) / 10.0;
  myb_magnetic_field* field = NULL;
  myb_diag diag = {0};
  if (myb_magnetic_field_load_table((const char*)data, size, 6371200.0, year, NULL, &field,
                                    &diag) != MYB_OK) {
    MYB_CHECK(field == NULL);
    return 0;
  }
  const int degree = myb_magnetic_field_degree(field);
  MYB_CHECK(degree >= 1);
  myb_gravity_work* work = NULL;
  if (degree <= 200 && myb_gravity_work_create(degree, NULL, &work) == MYB_OK) {
    const double r[3] = {7.0e6, -1.0e6, 2.0e6};
    double b[3];
    if (myb_magnetic_field_eval(field, work, r, b) == MYB_OK) {
      MYB_CHECK(fuzz_finite(b, 3));
    }
    myb_gravity_work_destroy(work);
  }
  myb_magnetic_field_destroy(field);
  return 0;
}
