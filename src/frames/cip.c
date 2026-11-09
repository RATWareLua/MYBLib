/* SPDX-License-Identifier: BSD-3-Clause */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "frames/rotation.h"
#include "myb/frames.h"
#include "myb/num.h"
#include "num/poisson.h"

static const double two_pi = 0x1.921fb54442d18p+2;
static const double seconds_per_day = 86400.0;
static const double seconds_per_century = 3155760000.0; /* 36525 days */

struct myb_cip {
  myb_allocator allocator;
  myb_internal_poisson_group series; /* X, Y, s + XY/2 */
  double era[2];                  /* ERA = 2 pi (era[0] + Tu + era[1] Tu) */
  double sprime_rate;             /* s' = sprime_rate * t * scale of X */
};

static myb_status read_numbers(const myb_kpool* pool, const char* name, size_t expected,
                               double* out, myb_diag* diag) {
  const double* v = NULL;
  size_t n = 0;
  const myb_status status = myb_kpool_numbers(pool, name, &v, &n);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "%s: %s", name,
                      status == MYB_E_NOT_FOUND ? "not defined" : "not numeric");
    return status;
  }
  if (n != expected) {
    myb_internal_diag(diag, 0, "%s: expected %zu numbers", name, expected);
    return MYB_E_FORMAT;
  }
  for (size_t i = 0; i < n; ++i) {
    if (!isfinite(v[i])) {
      myb_internal_diag(diag, 0, "%s: non-finite value", name);
      return MYB_E_FORMAT;
    }
    out[i] = v[i];
  }
  return MYB_OK;
}

myb_status myb_cip_create(const myb_kpool* pool, const char* prefix,
                          const myb_allocator* allocator, myb_cip** cip, myb_diag* diag) {
  MYB_CHECK(pool != NULL && prefix != NULL && cip != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_cip), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_cip* c = block;
  memset(c, 0, sizeof *c);
  c->allocator = a;
  char names[3][128];
  const char* list[3] = {names[0], names[1], names[2]};
  const char* suffixes[3] = {"X", "Y", "S"};
  for (int i = 0; i < 3; ++i) {
    snprintf(names[i], sizeof names[i], "%s_%s", prefix, suffixes[i]);
  }
  status = myb_internal_poisson_load(pool, prefix, list, 3, &a, &c->series, diag);
  char name[128];
  if (status == MYB_OK) {
    snprintf(name, sizeof name, "%s_ERA", prefix);
    status = read_numbers(pool, name, 2, c->era, diag);
  }
  if (status == MYB_OK) {
    snprintf(name, sizeof name, "%s_SPRIME_RATE", prefix);
    status = read_numbers(pool, name, 1, &c->sprime_rate, diag);
  }
  if (status != MYB_OK) {
    myb_cip_destroy(c);
    return status;
  }
  *cip = c;
  return MYB_OK;
}

void myb_cip_destroy(myb_cip* cip) {
  if (cip == NULL) {
    return;
  }
  const myb_allocator a = cip->allocator;
  myb_internal_poisson_free(&cip->series, &a);
  myb_internal_release(&a, cip);
}

static double centuries(myb_time tt) {
  return myb_time_seconds(tt) / seconds_per_century;
}

void myb_cip_xys(const myb_cip* cip, myb_time tt, double xys[3], double rates[3]) {
  MYB_CHECK(cip != NULL && xys != NULL);
  const double t = centuries(tt);
  double value[3];
  double rate[3];
  myb_internal_poisson_eval(&cip->series, t, value, rate);
  /* s = (s + XY/2) - XY/2 */
  xys[0] = value[0];
  xys[1] = value[1];
  xys[2] = value[2] - 0.5 * value[0] * value[1];
  if (rates != NULL) {
    const double sdot = rate[2] - 0.5 * (rate[0] * value[1] + value[0] * rate[1]);
    rates[0] = rate[0] / seconds_per_century;
    rates[1] = rate[1] / seconds_per_century;
    rates[2] = sdot / seconds_per_century;
  }
}

double myb_cip_era(const myb_cip* cip, myb_time ut1, double* rate) {
  MYB_CHECK(cip != NULL);
  /* Tu = days + fraction, the fraction of the day since 12h UT1 kept exact. */
  const int64_t day = ut1.s >= 0 ? ut1.s / 86400 : -((-ut1.s + 86399) / 86400);
  const double fraction = ((double)(ut1.s - day * 86400) + ut1.f) / seconds_per_day;
  /* turns = era[0] + fraction + era[1] day + era[1] fraction in double-word arithmetic (the
   * whole days are multiplied exactly), then reduced mod 1 */
  double low;
  const double high = myb_two_prod(cip->era[1], (double)day, &low);
  myb_dd turns = myb_dd_add_d(myb_dd_from_double(cip->era[0]), fraction);
  turns = myb_dd_add(turns, (myb_dd){high, low});
  turns = myb_dd_add_d(turns, cip->era[1] * fraction);
  const double whole = floor(turns.hi);
  turns = myb_dd_add_d(turns, -whole);
  double angle = two_pi * myb_dd_to_double(turns);
  if (angle < 0.0) {
    angle += two_pi;
  } else if (angle >= two_pi) {
    angle -= two_pi;
  }
  if (rate != NULL) {
    *rate = two_pi * (1.0 + cip->era[1]) / seconds_per_day;
  }
  return angle;
}

double myb_cip_sprime(const myb_cip* cip, myb_time tt) {
  MYB_CHECK(cip != NULL);
  return cip->sprime_rate * centuries(tt) * cip->series.scale[0];
}

void myb_cip_c2i(const double xys[3], const double rates[3], myb_rotation* rotation) {
  MYB_CHECK(xys != NULL && rotation != NULL);
  const double zero[3] = {0.0, 0.0, 0.0};
  const double* d = rates != NULL ? rates : zero;
  const double x = xys[0];
  const double y = xys[1];
  const double r2 = x * x + y * y;
  const double z = sqrt(1.0 - r2);
  const double a = 1.0 / (1.0 + z);
  const double zdot = -(x * d[0] + y * d[1]) / z;
  const double adot = -zdot * a * a;
  /* B^T, the transpose of eq. 5.10's matrix, and its derivative */
  const double bt[9] = {1.0 - a * x * x, -a * x * y, -x, -a * x * y, 1.0 - a * y * y, -y,
                        x, y, z};
  const double cross = -(adot * x * y + a * d[0] * y + a * x * d[1]);
  const double dbt[9] = {-(adot * x * x + 2.0 * a * x * d[0]), cross, -d[0],
                         cross, -(adot * y * y + 2.0 * a * y * d[1]), -d[1],
                         d[0], d[1], zdot};
  double r3[9];
  double dr3[9];
  myb_internal_rot_axis(2, -xys[2], r3, dr3);
  double first[9];
  double second[9];
  myb_internal_mat3_mul(r3, bt, rotation->m);
  myb_internal_mat3_mul(dr3, bt, first);
  myb_internal_mat3_mul(r3, dbt, second);
  for (int i = 0; i < 9; ++i) {
    rotation->dm[i] = -d[2] * first[i] + second[i];
  }
}

void myb_cip_terrestrial(const myb_cip* cip, myb_time tt, myb_time ut1, const myb_eop_values* eop,
                         myb_rotation* rotation) {
  MYB_CHECK(cip != NULL && rotation != NULL);
  const myb_eop_values none = {0};
  const myb_eop_values* e = eop != NULL ? eop : &none;
  double xys[3];
  double rates[3];
  myb_cip_xys(cip, tt, xys, rates);
  xys[0] += e->dx;
  xys[1] += e->dy;
  myb_rotation c2i;
  myb_cip_c2i(xys, rates, &c2i);
  double era_rate;
  const double era = myb_cip_era(cip, ut1, &era_rate);
  era_rate *= 1.0 - e->lod / seconds_per_day;
  myb_rotation spin;
  double dr[9];
  myb_internal_rot_axis(2, era, spin.m, dr);
  for (int i = 0; i < 9; ++i) {
    spin.dm[i] = dr[i] * era_rate;
  }
  /* W = R1(-yp) R2(-xp) R3(s') */
  double r1[9];
  double r2[9];
  double r3[9];
  double unused[9];
  myb_internal_rot_axis(0, -e->yp, r1, unused);
  myb_internal_rot_axis(1, -e->xp, r2, unused);
  myb_internal_rot_axis(2, myb_cip_sprime(cip, tt), r3, unused);
  myb_rotation w;
  myb_internal_mat3_mul(r2, r3, w.m);
  myb_internal_mat3_mul(r1, w.m, w.m);
  memset(w.dm, 0, sizeof w.dm);
  myb_rotation r;
  myb_rotation_compose(&spin, &c2i, &r);
  myb_rotation_compose(&w, &r, rotation);
}
