/* SPDX-License-Identifier: BSD-3-Clause */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/frames.h"
#include "myb/num.h"

static const double degree = 0x1.1df46a2529d39p-6; /* pi / 180 */
static const double seconds_per_day = 86400.0;
static const double seconds_per_century = 3155760000.0;

enum { MAX_PHASE_DEGREE = 4 };

struct myb_wgccre {
  myb_allocator allocator;
  int32_t frame;
  myb_time epoch;      /* TDB instant of the model's epoch */
  double ra[3];        /* deg, deg/century, deg/century^2 */
  double dec[3];
  double pm[3];        /* deg, deg/day, deg/day^2 */
  size_t angles;       /* phase angles */
  int degree;          /* of their polynomials */
  double* phase;       /* angles x (degree + 1) coefficients, deg and centuries */
  double* nut_ra;      /* angles coefficients each (missing ones zero) */
  double* nut_dec;
  double* nut_pm;
};

void myb_wgccre_destroy(myb_wgccre* model) {
  if (model == NULL) {
    return;
  }
  const myb_allocator a = model->allocator;
  myb_internal_release(&a, model->phase);
  myb_internal_release(&a, model->nut_ra);
  myb_internal_release(&a, model->nut_dec);
  myb_internal_release(&a, model->nut_pm);
  myb_internal_release(&a, model);
}

/* Up to `max` numbers of a variable; *count = 0 if absent. */
static myb_status optional(const myb_kpool* pool, const char* name, const double** values,
                           size_t* count, myb_diag* diag) {
  const myb_status status = myb_kpool_numbers(pool, name, values, count);
  if (status == MYB_E_NOT_FOUND) {
    *values = NULL;
    *count = 0;
    return MYB_OK;
  }
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: not numeric", name);
    return status;
  }
  for (size_t i = 0; i < *count; ++i) {
    if (!isfinite((*values)[i])) {
      myb_internal_diag(diag, 0, "%s: non-finite value", name);
      return MYB_E_FORMAT;
    }
  }
  return MYB_OK;
}

static myb_status polynomial3(const myb_kpool* pool, int64_t body, const char* what, int required,
                              double out[3], myb_diag* diag) {
  char name[64];
  snprintf(name, sizeof name, "BODY%" PRId64 "_%s", body, what);
  const double* v = NULL;
  size_t n = 0;
  const myb_status status = optional(pool, name, &v, &n, diag);
  if (status != MYB_OK) {
    return status;
  }
  if (n == 0 && required) {
    myb_internal_diag(diag, 0, "%s: not defined", name);
    return MYB_E_NOT_FOUND;
  }
  if (n > 3) {
    myb_internal_diag(diag, 0, "%s: more than 3 coefficients", name);
    return MYB_E_FORMAT;
  }
  out[0] = out[1] = out[2] = 0.0;
  for (size_t i = 0; i < n; ++i) {
    out[i] = v[i];
  }
  return MYB_OK;
}

myb_status myb_wgccre_create(const myb_kpool* pool, int64_t body, const myb_allocator* allocator,
                             myb_wgccre** model, myb_diag* diag) {
  MYB_CHECK(pool != NULL && model != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_wgccre), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_wgccre* m = block;
  memset(m, 0, sizeof *m);
  m->allocator = a;
  m->frame = 1;
  m->degree = 1;
  status = polynomial3(pool, body, "POLE_RA", 1, m->ra, diag);
  if (status == MYB_OK) {
    status = polynomial3(pool, body, "POLE_DEC", 1, m->dec, diag);
  }
  if (status == MYB_OK) {
    status = polynomial3(pool, body, "PM", 1, m->pm, diag);
  }
  char name[64];
  const double* v = NULL;
  size_t n = 0;
  if (status == MYB_OK) {
    snprintf(name, sizeof name, "BODY%" PRId64 "_CONSTANTS_REF_FRAME", body);
    status = optional(pool, name, &v, &n, diag);
    if (status == MYB_OK && n > 0) {
      if (n != 1 || !(v[0] >= 1.0 && v[0] <= 2147483647.0) || v[0] != (double)(int32_t)v[0]) {
        myb_internal_diag(diag, 0, "%s: expected one frame code", name);
        status = MYB_E_FORMAT;
      } else {
        m->frame = (int32_t)v[0];
      }
    }
  }
  if (status == MYB_OK) {
    snprintf(name, sizeof name, "BODY%" PRId64 "_CONSTANTS_JED_EPOCH", body);
    status = optional(pool, name, &v, &n, diag);
    if (status == MYB_OK && n > 0) {
      if (n != 1) {
        myb_internal_diag(diag, 0, "%s: expected one Julian date", name);
        status = MYB_E_FORMAT;
      } else {
        m->epoch = myb_time_from_jd(v[0], 0.0);
      }
    }
  }
  /* nutation-precession terms and their phase angles */
  const double* terms[3] = {NULL, NULL, NULL};
  size_t counts[3] = {0, 0, 0};
  const char* suffixes[3] = {"NUT_PREC_RA", "NUT_PREC_DEC", "NUT_PREC_PM"};
  size_t needed = 0;
  for (int k = 0; k < 3 && status == MYB_OK; ++k) {
    snprintf(name, sizeof name, "BODY%" PRId64 "_%s", body, suffixes[k]);
    status = optional(pool, name, &terms[k], &counts[k], diag);
    if (counts[k] > needed) {
      needed = counts[k];
    }
  }
  const int64_t barycenter = (body > 0 && body < 1000) ? body / 100 : body;
  const double* phase = NULL;
  size_t phase_count = 0;
  if (status == MYB_OK && needed > 0) {
    snprintf(name, sizeof name, "BODY%" PRId64 "_MAX_PHASE_DEGREE", barycenter);
    status = optional(pool, name, &v, &n, diag);
    if (status == MYB_OK && n > 0) {
      if (n != 1 || !(v[0] >= 1.0 && v[0] <= (double)MAX_PHASE_DEGREE) ||
          v[0] != (double)(int)v[0]) {
        myb_internal_diag(diag, 0, "%s: expected an integer in 1 .. 4", name);
        status = MYB_E_FORMAT;
      } else {
        m->degree = (int)v[0];
      }
    }
    if (status == MYB_OK) {
      snprintf(name, sizeof name, "BODY%" PRId64 "_NUT_PREC_ANGLES", barycenter);
      status = optional(pool, name, &phase, &phase_count, diag);
      if (status == MYB_OK && phase_count == 0) {
        myb_internal_diag(diag, 0, "%s: not defined", name);
        status = MYB_E_NOT_FOUND;
      } else if (status == MYB_OK && phase_count % (size_t)(m->degree + 1) != 0) {
        myb_internal_diag(diag, 0, "%s: not a multiple of %d coefficients", name,
                          m->degree + 1);
        status = MYB_E_FORMAT;
      } else if (status == MYB_OK && phase_count / (size_t)(m->degree + 1) < needed) {
        myb_internal_diag(diag, 0, "BODY%" PRId64 ": more nutation terms than phase angles",
                          body);
        status = MYB_E_FORMAT;
      }
    }
  }
  if (status == MYB_OK && needed > 0) {
    m->angles = needed;
    void* blocks[4] = {NULL, NULL, NULL, NULL};
    const size_t sizes[4] = {needed * (size_t)(m->degree + 1), needed, needed, needed};
    for (int k = 0; k < 4 && status == MYB_OK; ++k) {
      status = myb_internal_alloc_array(&a, sizes[k], sizeof(double), &blocks[k]);
    }
    m->phase = blocks[0];
    m->nut_ra = blocks[1];
    m->nut_dec = blocks[2];
    m->nut_pm = blocks[3];
    if (status == MYB_OK) {
      memcpy(m->phase, phase, sizes[0] * sizeof(double));
      double* targets[3] = {m->nut_ra, m->nut_dec, m->nut_pm};
      for (int k = 0; k < 3; ++k) {
        for (size_t i = 0; i < needed; ++i) {
          targets[k][i] = i < counts[k] ? terms[k][i] : 0.0;
        }
      }
    }
  }
  if (status != MYB_OK) {
    myb_wgccre_destroy(m);
    return status;
  }
  *model = m;
  return MYB_OK;
}

int32_t myb_wgccre_frame(const myb_wgccre* model) {
  MYB_CHECK(model != NULL);
  return model->frame;
}

/* x mod 360 for x = high + low (degrees), as a double in [0, 360). */
static double reduce_degrees(double high, double low) {
  const double turns = floor(high / 360.0);
  myb_dd r = myb_dd_add_d((myb_dd){high, low}, -360.0 * turns);
  double x = myb_dd_to_double(r);
  if (x < 0.0) {
    x += 360.0;
  } else if (x >= 360.0) {
    x -= 360.0;
  }
  return x;
}

void myb_wgccre_euler(const myb_wgccre* model, myb_time t, double angles[3], double rates[3]) {
  MYB_CHECK(model != NULL && angles != NULL && rates != NULL);
  const myb_dd elapsed = myb_time_diff(t, model->epoch); /* s */
  const double centuries = elapsed.hi / seconds_per_century;
  const myb_dd days = myb_dd_div_d(elapsed, seconds_per_day);
  double ra = model->ra[0] + centuries * (model->ra[1] + centuries * model->ra[2]);
  double dec = model->dec[0] + centuries * (model->dec[1] + centuries * model->dec[2]);
  double ra_rate = (model->ra[1] + 2.0 * centuries * model->ra[2]) / seconds_per_century;
  double dec_rate = (model->dec[1] + 2.0 * centuries * model->dec[2]) / seconds_per_century;
  /* W = w0 + w1 d + w2 d^2 with w1 d in double-word arithmetic, reduced mod 360 */
  const myb_dd linear = myb_dd_add_d(myb_dd_mul_d(days, model->pm[1]), model->pm[0]);
  double w = reduce_degrees(linear.hi, linear.lo) + model->pm[2] * days.hi * days.hi;
  double w_rate = (model->pm[1] + 2.0 * model->pm[2] * days.hi) / seconds_per_day;
  const int k = model->degree;
  for (size_t i = 0; i < model->angles; ++i) {
    const double* c = model->phase + i * (size_t)(k + 1);
    double theta = 0.0;
    double theta_rate = 0.0;
    for (int j = k; j >= 0; --j) {
      theta = theta * centuries + c[j];
      if (j >= 1) {
        theta_rate = theta_rate * centuries + (double)j * c[j];
      }
    }
    theta = reduce_degrees(theta, 0.0) * degree;
    theta_rate *= degree / seconds_per_century; /* rad/s */
    const double s = myb_sin(theta);
    const double co = myb_cos(theta);
    ra += model->nut_ra[i] * s;
    dec += model->nut_dec[i] * co;
    w += model->nut_pm[i] * s;
    ra_rate += model->nut_ra[i] * co * theta_rate; /* deg/s */
    dec_rate -= model->nut_dec[i] * s * theta_rate;
    w_rate += model->nut_pm[i] * co * theta_rate;
  }
  angles[0] = (90.0 + ra) * degree;
  angles[1] = (90.0 - dec) * degree;
  angles[2] = w * degree;
  rates[0] = ra_rate * degree;
  rates[1] = -dec_rate * degree;
  rates[2] = w_rate * degree;
}

void myb_wgccre_rotation(const myb_wgccre* model, myb_time t, myb_rotation* rotation) {
  MYB_CHECK(rotation != NULL);
  double angles[3];
  double rates[3];
  myb_wgccre_euler(model, t, angles, rates);
  myb_rotation_euler313(angles, rates, rotation);
}
