/* SPDX-License-Identifier: BSD-3-Clause */
/* Equinox-based systems: IAU 1976 precession, IAU 1980 obliquity and nutation, TEME, GMST82. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "frames/rotation.h"
#include "myb/frames.h"
#include "num/poisson.h"
#include "num/vec3.h"

static const double two_pi = 0x1.921fb54442d18p+2;
static const double arcsec = 4.848136811095359935899141e-6; /* rad */
static const double seconds_per_day = 86400.0;
static const double seconds_per_century = 3155760000.0; /* 36525 days */

enum { CUBIC = 4 };

struct myb_equinox {
  myb_allocator allocator;
  myb_internal_poisson_group nutation; /* Delta psi, Delta eps */
  double obliquity[CUBIC];             /* arcsec */
  double precession[3][CUBIC];         /* zeta, z, theta (arcsec) */
  double gmst[CUBIC];                  /* s */
};

static myb_status read_numbers(const myb_kpool* pool, const char* prefix, const char* suffix,
                               size_t expected, double* out, myb_diag* diag) {
  char name[128];
  snprintf(name, sizeof name, "%s_%s", prefix, suffix);
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

myb_status myb_equinox_create(const myb_kpool* pool, const char* prefix,
                              const myb_allocator* allocator, myb_equinox** equinox,
                              myb_diag* diag) {
  MYB_CHECK(pool != NULL && prefix != NULL && equinox != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_equinox), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_equinox* e = block;
  memset(e, 0, sizeof *e);
  e->allocator = a;
  char names[2][128];
  const char* list[2] = {names[0], names[1]};
  snprintf(names[0], sizeof names[0], "%s_DPSI", prefix);
  snprintf(names[1], sizeof names[1], "%s_DEPS", prefix);
  status = myb_internal_poisson_load(pool, prefix, list, 2, &a, &e->nutation, diag);
  const char* angles[3] = {"PRECESSION_ZETA", "PRECESSION_Z", "PRECESSION_THETA"};
  for (int i = 0; i < 3 && status == MYB_OK; ++i) {
    status = read_numbers(pool, prefix, angles[i], CUBIC, e->precession[i], diag);
  }
  if (status == MYB_OK) {
    status = read_numbers(pool, prefix, "OBLIQUITY", CUBIC, e->obliquity, diag);
  }
  if (status == MYB_OK) {
    status = read_numbers(pool, prefix, "GMST", CUBIC, e->gmst, diag);
  }
  if (status != MYB_OK) {
    myb_equinox_destroy(e);
    return status;
  }
  *equinox = e;
  return MYB_OK;
}

void myb_equinox_destroy(myb_equinox* equinox) {
  if (equinox == NULL) {
    return;
  }
  const myb_allocator a = equinox->allocator;
  myb_internal_poisson_free(&equinox->nutation, &a);
  myb_internal_release(&a, equinox);
}

static double centuries(myb_time tt) {
  return myb_time_seconds(tt) / seconds_per_century;
}

/* c0 + c1 t + c2 t^2 + c3 t^3 (Horner) and its derivative in t. */
static double cubic(const double c[CUBIC], double t, double* derivative) {
  *derivative = (3.0 * c[3] * t + 2.0 * c[2]) * t + c[1];
  return ((c[3] * t + c[2]) * t + c[1]) * t + c[0];
}

void myb_equinox_precession(const myb_equinox* equinox, myb_time tt, double angles[3],
                            double rates[3]) {
  MYB_CHECK(equinox != NULL && angles != NULL);
  const double t = centuries(tt);
  for (int i = 0; i < 3; ++i) {
    double d;
    angles[i] = cubic(equinox->precession[i], t, &d) * arcsec;
    if (rates != NULL) {
      rates[i] = d * arcsec / seconds_per_century;
    }
  }
}

double myb_equinox_obliquity(const myb_equinox* equinox, myb_time tt, double* rate) {
  MYB_CHECK(equinox != NULL);
  double d;
  const double eps = cubic(equinox->obliquity, centuries(tt), &d) * arcsec;
  if (rate != NULL) {
    *rate = d * arcsec / seconds_per_century;
  }
  return eps;
}

void myb_equinox_nutation(const myb_equinox* equinox, myb_time tt, double nutation[2],
                          double rates[2]) {
  MYB_CHECK(equinox != NULL && nutation != NULL);
  double value[2];
  double rate[2];
  myb_internal_poisson_eval(&equinox->nutation, centuries(tt), value, rate);
  nutation[0] = value[0];
  nutation[1] = value[1];
  if (rates != NULL) {
    rates[0] = rate[0] / seconds_per_century;
    rates[1] = rate[1] / seconds_per_century;
  }
}

/* r = R_axis(angle) with dr = dR/dt for the angle's rate. */
static void axis_rotation(int axis, double angle, double rate, myb_rotation* r) {
  double dm[9];
  myb_internal_rot_axis(axis, angle, r->m, dm);
  for (int i = 0; i < 9; ++i) {
    r->dm[i] = dm[i] * rate;
  }
}

static void precession_matrix(const myb_equinox* e, myb_time tt, myb_rotation* p) {
  double a[3];
  double ra[3];
  myb_equinox_precession(e, tt, a, ra);
  myb_rotation r1;
  myb_rotation r2;
  myb_rotation r3;
  axis_rotation(2, -a[0], -ra[0], &r1); /* R3(-zeta) */
  axis_rotation(1, a[2], ra[2], &r2);   /* R2(theta) */
  axis_rotation(2, -a[1], -ra[1], &r3); /* R3(-z) */
  myb_rotation w;
  myb_rotation_compose(&r2, &r1, &w);
  myb_rotation_compose(&r3, &w, p);
}

static void nutation_matrix(const myb_equinox* e, myb_time tt, myb_rotation* n) {
  double rate_eps;
  const double eps = myb_equinox_obliquity(e, tt, &rate_eps);
  double nut[2];
  double rnut[2];
  myb_equinox_nutation(e, tt, nut, rnut);
  myb_rotation r1;
  myb_rotation r2;
  myb_rotation r3;
  axis_rotation(0, eps, rate_eps, &r1);                              /* R1(eps_A) */
  axis_rotation(2, -nut[0], -rnut[0], &r2);                          /* R3(-Delta psi) */
  axis_rotation(0, -(eps + nut[1]), -(rate_eps + rnut[1]), &r3);     /* R1(-eps) */
  myb_rotation w;
  myb_rotation_compose(&r2, &r1, &w);
  myb_rotation_compose(&r3, &w, n);
}

void myb_equinox_rotation(const myb_equinox* equinox, myb_equinox_frame frame, myb_time tt,
                          myb_rotation* rotation) {
  MYB_CHECK(equinox != NULL && rotation != NULL &&
            (frame == MYB_EQUINOX_MOD || frame == MYB_EQUINOX_TOD || frame == MYB_EQUINOX_TEME));
  myb_rotation p;
  precession_matrix(equinox, tt, &p);
  if (frame == MYB_EQUINOX_MOD) {
    *rotation = p;
    return;
  }
  myb_rotation n;
  nutation_matrix(equinox, tt, &n);
  myb_rotation np;
  myb_rotation_compose(&n, &p, &np);
  if (frame == MYB_EQUINOX_TOD) {
    *rotation = np;
    return;
  }
  /* TEME: rows are the axes in J2000. z = the TOD pole, x = the mean equinox (first row of P)
   * without its component along z, normalised; y = z x x. */
  const double* z = np.m + 6;
  const double* zd = np.dm + 6;
  const double* xm = p.m;
  const double* xmd = p.dm;
  const double c = myb_vec3_dot(xm, z);
  const double cd = myb_vec3_dot(xmd, z) + myb_vec3_dot(xm, zd);
  double a[3];
  double ad[3];
  for (int i = 0; i < 3; ++i) {
    a[i] = xm[i] - c * z[i];
    ad[i] = xmd[i] - cd * z[i] - c * zd[i];
  }
  const double na = myb_vec3_norm(a);
  double x[3];
  double xd[3];
  for (int i = 0; i < 3; ++i) {
    x[i] = a[i] / na;
  }
  const double xad = myb_vec3_dot(x, ad);
  for (int i = 0; i < 3; ++i) {
    xd[i] = (ad[i] - x[i] * xad) / na;
  }
  double y[3];
  double y1[3];
  double y2[3];
  myb_vec3_cross(z, x, y);
  myb_vec3_cross(zd, x, y1);
  myb_vec3_cross(z, xd, y2);
  for (int i = 0; i < 3; ++i) {
    rotation->m[i] = x[i];
    rotation->m[3 + i] = y[i];
    rotation->m[6 + i] = z[i];
    rotation->dm[i] = xd[i];
    rotation->dm[3 + i] = y1[i] + y2[i];
    rotation->dm[6 + i] = zd[i];
  }
}

double myb_equinox_gmst(const myb_equinox* equinox, myb_time ut1, double* rate) {
  MYB_CHECK(equinox != NULL);
  const double* g = equinox->gmst;
  const double tu = myb_time_seconds(ut1) / seconds_per_century;
  /* seconds of UT1 since the last 0h: J2000.0 is 12h, whole days drop out */
  const int64_t day = ut1.s >= 0 ? ut1.s / 86400 : -((-ut1.s + 86399) / 86400);
  const double noon = (double)(ut1.s - day * 86400) + ut1.f; /* [0, 86400] since 12h */
  const double d = (3.0 * g[3] * tu + 2.0 * g[2]) * tu + g[1];
  const double drift = ((g[3] * tu + g[2]) * tu + g[1]) * tu;
  double seconds = (g[0] - 0.5 * seconds_per_day) + drift + noon;
  seconds -= seconds_per_day * floor(seconds / seconds_per_day);
  double angle = seconds * (two_pi / seconds_per_day);
  if (angle >= two_pi) {
    angle -= two_pi;
  }
  if (rate != NULL) {
    *rate = (two_pi / seconds_per_day) * (1.0 + d / seconds_per_century);
  }
  return angle;
}

void myb_equinox_teme_pef(const myb_equinox* equinox, myb_time ut1, myb_rotation* rotation) {
  MYB_CHECK(equinox != NULL && rotation != NULL);
  double rate;
  const double gmst = myb_equinox_gmst(equinox, ut1, &rate);
  axis_rotation(2, gmst, rate, rotation);
}
