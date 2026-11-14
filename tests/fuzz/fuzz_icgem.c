/* SPDX-License-Identifier: BSD-3-Clause */
/* ICGEM gravity fields (myb_gravity_field_load_icgem), truncated to a degree drawn from the
 * input: the coefficients asked for, the time-variable terms set to an instant, the potential,
 * the acceleration and its gradient evaluated outside the body. */
#include "fuzz.h"
#include "myb/grav.h"

MYB_FUZZ_NO_INITIALIZE

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  const uint64_t hash = fuzz_hash(data, size);
  const int max_degree = 1 + (int)(hash % 40u);
  myb_gravity_field* field = NULL;
  myb_diag diag = {0};
  if (myb_gravity_field_load_icgem((const char*)data, size, max_degree, NULL, &field, &diag) !=
      MYB_OK) {
    MYB_CHECK(field == NULL);
    return 0;
  }
  const int degree = myb_gravity_field_degree(field);
  MYB_CHECK(degree >= 0 && degree <= max_degree);
  const double gm = myb_gravity_field_gm(field);
  const double radius = myb_gravity_field_radius(field);
  MYB_CHECK(gm > 0.0 && radius > 0.0);
  (void)myb_gravity_field_tide_system(field);
  for (int n = 0; n <= degree + 1; ++n) {
    double c = 0.0;
    double s = 0.0;
    (void)myb_gravity_field_coefficient(field, n, n / 2, &c, &s);
  }
  myb_gravity_field_set_time(field, myb_time_from_seconds((double)(int64_t)(hash >> 34)));
  myb_gravity_work* work = NULL;
  if (myb_gravity_work_create(degree, NULL, &work) == MYB_OK) {
    const double r[3] = {1.1 * radius, 0.2 * radius, -0.3 * radius};
    double potential = 0.0;
    double a[3];
    (void)myb_gravity_field_eval(field, work, r, 0, -1, &potential, a);
    double gradient[9];
    (void)myb_gravity_field_gradient(field, work, r, 2, -1, gradient);
    myb_gravity_work_destroy(work);
  }
  myb_gravity_field_destroy(field);
  return 0;
}
