/* SPDX-License-Identifier: BSD-3-Clause */
#include "myb/body.h"

#include <math.h> /* sqrt only: correctly rounded by IEEE-754 */

#include "core/base.h"

static int finite(double x) {
  return x - x == 0.0;
}

static int valid(myb_param p) {
  return finite(p.value) && finite(p.sigma) && p.sigma >= 0.0;
}

/* sqrt(a^2 + b^2) without overflow or harmful underflow (a, b >= 0). */
static double norm2(double a, double b) {
  const double big = a > b ? a : b;
  const double small = a > b ? b : a;
  if (big == 0.0) {
    return 0.0;
  }
  const double r = small / big;
  return big * sqrt(1.0 + r * r);
}

myb_status myb_param_product(myb_param a, myb_param b, myb_param* product) {
  MYB_CHECK(product != NULL);
  if (!valid(a) || !valid(b)) {
    return MYB_E_DOMAIN;
  }
  const double value = a.value * b.value;
  const double sigma = norm2(fabs(b.value) * a.sigma, fabs(a.value) * b.sigma);
  if (!finite(value) || !finite(sigma)) {
    return MYB_E_DOMAIN;
  }
  *product = (myb_param){.value = value, .sigma = sigma};
  return MYB_OK;
}

myb_status myb_gm_from_mass(myb_param mass, myb_param g, myb_param* gm) {
  MYB_CHECK(gm != NULL);
  if (!(mass.value >= 0.0) || !(g.value > 0.0)) {
    return MYB_E_DOMAIN;
  }
  return myb_param_product(mass, g, gm);
}

static myb_status read_one(const myb_kpool* pool, const char* name, double* value) {
  const double* values = NULL;
  size_t count = 0;
  const myb_status status = myb_kpool_numbers(pool, name, &values, &count);
  if (status != MYB_OK) {
    return status;
  }
  if (count != 1) {
    return MYB_E_FORMAT;
  }
  *value = values[0];
  return MYB_OK;
}

myb_status myb_param_from_kpool(const myb_kpool* pool, const char* name, const char* sigma_name,
                                double scale, myb_param* param) {
  MYB_CHECK(pool != NULL && name != NULL && param != NULL);
  double value = 0.0;
  double sigma = 0.0;
  myb_status status = read_one(pool, name, &value);
  if (status == MYB_OK && sigma_name != NULL) {
    status = read_one(pool, sigma_name, &sigma);
  }
  if (status != MYB_OK) {
    return status;
  }
  const myb_param p = {.value = value * scale, .sigma = sigma * fabs(scale)};
  if (!valid(p)) {
    return MYB_E_DOMAIN;
  }
  *param = p;
  return MYB_OK;
}

myb_status myb_body_validate(const myb_body* body) {
  MYB_CHECK(body != NULL);
  if (!valid(body->gm) || !(body->gm.value >= 0.0)) {
    return MYB_E_DOMAIN;
  }
  if (body->motion != MYB_MOTION_DYNAMIC && body->motion != MYB_MOTION_KINEMATIC) {
    return MYB_E_DOMAIN;
  }
  return MYB_OK;
}
