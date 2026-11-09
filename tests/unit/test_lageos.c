/* SPDX-License-Identifier: BSD-3-Clause */
/* LAGEOS-1 over one week against the combined orbit of the International Laser Ranging Service
 * (ILRS, ilrsa.orb.lageos1.250607.v80.sp3: 2025-06-01 .. 06-07, SLRF2020; the excerpt
 * tests/data/sp3/ilrsa.lageos1.250607.excerpt.sp3 keeps every 3rd epoch, 6 min, 1680 positions;
 * the file of EDC/DGFI, ftp://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/lageos1/2025/250607/,
 * SHA-256 of the .gz ba1d9933c870c44b4e0533f423c343340e09947bdbb6096f6410de86a8bb7ef5).
 *
 * The model is assembled from myb's pieces with every number from data:
 *  - Earth: GOCO06s to degree 30 (zero-tide, with its trends and annual terms), solid tides of
 *    IERS 2010 (steps 1-3), FES2004 ocean tides to degree 20, solid and ocean pole tides; the
 *    orientation of IAU 2006/2000A with EOP 20 C04 and its diurnal corrections;
 *  - Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn from DE440 (GM of header.440);
 *  - relativity of IERS 2010 eq. 10.12 (Schwarzschild, Lense-Thirring, de Sitter);
 *  - solar radiation pressure on a sphere (cr fitted) with the penumbra of the Earth and the
 *    Moon and the limb darkening of the Sun; Earth's albedo and infrared after Knocke et al.;
 *  - empirical accelerations as estimated by the ILRS analysis centres: a constant along-track
 *    term and once-per-revolution along- and cross-track terms.
 * The initial state, cr and the five empirical parameters are fitted to the positions by
 * Gauss-Newton with partials by finite differences; the residual is the misfit of our dynamics
 * to the ILRS orbit (itself fitted to laser ranges at the cm level by models of the same kind).
 * Needs data/download (scripts/fetch-data.cmd: de440.bsp, header.440, eopc04.1962-now,
 * GOCO06s.gfc, fes2004_Cnm-Snm.dat); skipped otherwise. Measured contributions (RMS over the
 * week after the fit, the full model 2.0 cm): without the solid tides 21 cm, without the ocean
 * tides 6.3 cm, without relativity 2.06 cm (radial 5.3 mm against 2.6), without shadows 2.1 cm;
 * MYB_LAGEOS_OFF selects parts to leave out, MYB_LAGEOS_PROFILE times an evaluation. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mtest.h"
#include "myb/myb.h"
#include "core/pool.h"
#include "num/vec3.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

enum { DEGREE = 30, OCEAN_DEGREE = 20, PARAMETERS = 12, STRIDE = 1, BODIES = 7 };

static const double c_light = 299792458.0;
static const double area_over_mass = 3.141592653589793 * 0.30 * 0.30 / 406.965; /* LAGEOS-1 */
static const double solar_flux_1au = 1361.0;  /* W m^-2 (Kopp & Lean 2011) */
static const double au_m = 149597870700.0;
static const double earth_j = 9.8e8;           /* angular momentum per unit mass, m^2/s */
static const double sun_radius = 695700000.0;  /* IAU 2015 B3 */
static const double moon_radius = 1737400.0;
static const double earth_radius = 6378137.0;

typedef struct model {
  myb_timesys* sys;
  myb_scale tt;
  myb_scale tdb;
  myb_cip* cip;
  myb_eop* eop;
  myb_spk* spk;
  myb_gravity_field* field;
  myb_gravity_work* work;
  myb_tides* tides;
  myb_time t0; /* TT */
  int body[BODIES];
  double gm_body[BODIES];
  double gm_earth;
  double cr;
  double empirical[5]; /* along: constant, cos u, sin u; cross: cos u, sin u */
  /* what depends on t only, kept for the next call at the same t (IAS15 evaluates every
   * substep of a step several times while it iterates) */
  int off; /* MYB_LAGEOS_OFF: 1 solid tides, 2 ocean, 4 pole, 8 relativity, 16 albedo,
              32 planets, 128 shadows (for the table of contributions) */
  int cached;
  double cache_t;
  myb_rotation rot;
  myb_rotation inv;
  myb_eop_values eop_values;
  /* TDB - TT (s) every hour from t = -1 h: the series of 787 terms costs more than all the
   * rest of an evaluation; 4-point interpolation of the table is exact to 1e-13 s */
  double tdb_tt[200];
  double s[BODIES][6];
} model;

static int read_constant(const char* header, const char* name, double* value) {
  const char* g1040 = strstr(header, "GROUP   1040");
  const char* g1041 = strstr(header, "GROUP   1041");
  if (g1040 == NULL || g1041 == NULL) {
    return 0;
  }
  const size_t len = strlen(name);
  const char* p = g1040 + 12;
  double count = 0.0;
  size_t used = 0;
  /* bound: white space */
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  int index = -1;
  for (int k = 0; k < (int)count; ++k) {
    /* bound: white space */
    while (*p == ' ' || *p == '\n' || *p == '\r') {
      ++p;
    }
    size_t n = 0;
    /* bound: a name */
    while (p[n] != ' ' && p[n] != '\n' && p[n] != '\r' && p[n] != '\0') {
      ++n;
    }
    if (n == len && memcmp(p, name, len) == 0) {
      index = k;
    }
    p += n;
  }
  if (index < 0) {
    return 0;
  }
  p = g1041 + 12;
  /* bound: white space */
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  for (int k = 0; k <= index; ++k) {
    /* bound: white space */
    while (*p == ' ' || *p == '\n' || *p == '\r') {
      ++p;
    }
    if (myb_parse_double(p, strlen(p), value, &used) != MYB_OK) {
      return 0;
    }
    p += used;
  }
  return 1;
}

static int read_file(const char* path, char* buffer, size_t size) {
  FILE* f = fopen(path, "rb");
  if (f == NULL) {
    return 0;
  }
  const size_t n = fread(buffer, 1, size - 1, f);
  fclose(f);
  buffer[n] = '\0';
  return n > 0;
}

static void scale_add3(double* a, double k, const double* b) {
  a[0] += k * b[0];
  a[1] += k * b[1];
  a[2] += k * b[2];
}

/* Earth orientation at TT instant t: celestial -> terrestrial, EOP values, UT1. */
static myb_status orientation(const model* m, myb_time t, myb_rotation* rot,
                              myb_eop_values* eop, myb_time* ut1) {
  const myb_status st = myb_eop_at(m->eop, m->sys, t, m->tt, eop, ut1);
  if (st != MYB_OK) {
    return st;
  }
  myb_cip_terrestrial(m->cip, t, *ut1, eop, rot);
  return MYB_OK;
}

static myb_status acceleration(void* context, double t, const double* x, const double* v,
                               double* a) {
  model* m = context;
  const myb_time tt = myb_time_add(m->t0, t);
  myb_status st = MYB_OK;
  if (!m->cached || m->cache_t != t) {
    m->cached = 0;
    const double u = t / 3600.0 + 1.0;
    const int k0 = (int)floor(u) - 1;
    if (k0 < 0 || k0 + 3 >= 200) {
      return MYB_E_RANGE;
    }
    const double w = u - (double)k0; /* in [1, 2) */
    const double* f = m->tdb_tt + k0;
    const double dtdb = -f[0] * (w - 1.0) * (w - 2.0) * (w - 3.0) / 6.0 +
                        f[1] * w * (w - 2.0) * (w - 3.0) / 2.0 -
                        f[2] * w * (w - 1.0) * (w - 3.0) / 2.0 +
                        f[3] * w * (w - 1.0) * (w - 2.0) / 6.0;
    const myb_time tdb = myb_time_add(tt, dtdb);
    myb_time ut1;
    st = orientation(m, tt, &m->rot, &m->eop_values, &ut1);
    if (st != MYB_OK) {
      return st;
    }
    myb_rotation_invert(&m->rot, &m->inv);
    /* bodies, geocentric (GCRS axes = ICRF axes) */
    for (int b = 0; b < BODIES; ++b) {
      st = myb_spk_state(m->spk, m->body[b], 399, tdb, m->s[b], NULL);
      if (st != MYB_OK) {
        return st;
      }
    }
    /* the field with its time-variable terms and tides */
    myb_gravity_field_set_time(m->field, tt);
    myb_tide_raiser raisers[2];
    raisers[0].gm = m->gm_body[0];
    myb_rotation_apply_position(&m->rot, m->s[0], raisers[0].r);
    raisers[1].gm = m->gm_body[1];
    myb_rotation_apply_position(&m->rot, m->s[1], raisers[1].r);
    const myb_tides_input tin = {.raisers = raisers, .count = 2, .t = ut1, .xp = m->eop_values.xp,
                                 .yp = m->eop_values.yp,
                                 .flags = MYB_TIDES_ALL &
                                          ~((m->off & 1 ? MYB_TIDES_SOLID | MYB_TIDES_FREQUENCY : 0) |
                                            (m->off & 256 ? MYB_TIDES_FREQUENCY : 0) |
                                            (m->off & 2 ? MYB_TIDES_OCEAN : 0) |
                                            (m->off & 4 ? MYB_TIDES_POLE_SOLID |
                                                              MYB_TIDES_POLE_OCEAN
                                                        : 0))};
    st = myb_tides_apply(m->tides, m->field, &tin);
    if (st != MYB_OK) {
      return st;
    }
    m->cached = 1;
    m->cache_t = t;
  }
  const myb_rotation rot = m->rot;
  const myb_rotation inv = m->inv;
  const double (*s)[6] = (const double (*)[6])m->s;
  const double* sun = s[0];
  const double* moon = s[1];
  double r_itrs[3];
  myb_rotation_apply_position(&rot, x, r_itrs);
  double g_itrs[3];
  st = myb_gravity_field_eval(m->field, m->work, r_itrs, 0, -1, NULL, g_itrs);
  if (st != MYB_OK) {
    return st;
  }
  myb_rotation_apply_position(&inv, g_itrs, a);
  /* third bodies */
  for (int b = 0; b < ((m->off & 32) ? 2 : BODIES); ++b) {
    myb_force_third_body(m->gm_body[b], s[b], x, a);
  }
  /* relativity */
  const myb_ppn ppn = {c_light, 1.0, 1.0};
  double tmp[3];
  const double rel_on = (m->off & 8) ? 0.0 : 1.0;
  st = myb_rel_schwarzschild(&ppn, m->gm_earth, x, v, tmp);
  if (st != MYB_OK) {
    return st;
  }
  scale_add3(a, rel_on, tmp);
  const double j_itrs[3] = {0.0, 0.0, earth_j};
  double j[3];
  myb_rotation_apply_position(&inv, j_itrs, j);
  st = myb_rel_lense_thirring(&ppn, m->gm_earth, j, x, v, tmp);
  if (st != MYB_OK) {
    return st;
  }
  scale_add3(a, rel_on, tmp);
  const double R[3] = {-sun[0], -sun[1], -sun[2]};
  const double V[3] = {-sun[3], -sun[4], -sun[5]};
  st = myb_rel_de_sitter(&ppn, m->gm_body[0], R, V, v, tmp);
  if (st != MYB_OK) {
    return st;
  }
  scale_add3(a, rel_on, tmp);
  /* solar radiation pressure with shadows */
  const double d[3] = {x[0] - sun[0], x[1] - sun[1], x[2] - sun[2]};
  const double dn = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double u[3] = {d[0] / dn, d[1] / dn, d[2] / dn};
  const double source[3] = {-d[0], -d[1], -d[2]};
  const double earth[3] = {-x[0], -x[1], -x[2]};
  const double luna[3] = {moon[0] - x[0], moon[1] - x[1], moon[2] - x[2]};
  const double limb[2] = {0.6, 0.0};
  const double nu = (m->off & 128) ? 1.0 : myb_shadow_fraction(source, sun_radius, earth, earth_radius, limb) *
                    myb_shadow_fraction(source, sun_radius, luna, moon_radius, limb);
  const double pressure = solar_flux_1au * (au_m / dn) * (au_m / dn) / c_light;
  myb_force_radiation_sphere(nu * pressure, m->cr, area_over_mass, u, a);
  /* albedo and infrared of the Earth (Knocke et al. 1988), seasonal phase from 1981-12-22 */
  static const myb_planet_radiation knocke = {
      .radius = 6378137.0,
      .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
      .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
      .omega = 2.0 * 3.141592653589793 / (365.25 * 86400.0)};
  const double rs = sqrt(sun[0] * sun[0] + sun[1] * sun[1] + sun[2] * sun[2]);
  const double phase = knocke.omega * (myb_time_seconds(tt) + 568915200.0);
  double sun_itrs[3];
  myb_rotation_apply_position(&rot, sun, sun_itrs);
  double albedo[3] = {0.0, 0.0, 0.0};
  st = myb_force_planet_radiation(&knocke, r_itrs, sun_itrs,
                                  solar_flux_1au * (au_m / rs) * (au_m / rs), phase,
                                  m->cr * area_over_mass, 10, albedo);
  if (st != MYB_OK) {
    return st;
  }
  myb_rotation_apply_position(&inv, albedo, tmp);
  scale_add3(a, (m->off & 16) ? 0.0 : 1.0, tmp);
  /* empirical: argument of latitude u from the node */
  double h[3];
  myb_vec3_cross(x, v, h);
  const double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  const double nodal[3] = {-h[1], h[0], 0.0}; /* z x h */
  const double nn = hypot(nodal[0], nodal[1]);
  const double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double hxn[3];
  const double hu[3] = {h[0] / hn, h[1] / hn, h[2] / hn};
  const double nu_[3] = {nodal[0] / nn, nodal[1] / nn, 0.0};
  myb_vec3_cross(hu, nu_, hxn);
  const double cu = (nu_[0] * x[0] + nu_[1] * x[1]) / r;
  const double su = (hxn[0] * x[0] + hxn[1] * x[1] + hxn[2] * x[2]) / r;
  const double* e = m->empirical;
  const double rtn[3] = {0.0, e[0] + e[1] * cu + e[2] * su, e[3] * cu + e[4] * su};
  return myb_force_empirical_rtn(x, v, rtn, a);
}

/* MYB_LAGEOS_PROFILE: time of the pieces of one evaluation of the acceleration. */
static void profile(model* m, const double* x, const double* v) {
  enum { REPS = 300 };
  double t_rot = 0.0;
  double t_tdb = 0.0;
  double t_spk = 0.0;
  double t_set = 0.0;
  double t_tides = 0.0;
  double t_eval = 0.0;
  double t_alb = 0.0;
  double t_all = 0.0;
  myb_rotation rot;
  myb_eop_values eop;
  myb_time ut1;
  myb_time tdb;
  double s[6];
  double a[3];
  for (int k = 0; k < REPS; ++k) {
    const myb_time tt = myb_time_add(m->t0, 100.0 * k);
    clock_t c0 = clock();
    (void)orientation(m, tt, &rot, &eop, &ut1);
    clock_t c1 = clock();
    (void)myb_timesys_convert(m->sys, tt, m->tt, m->tdb, &tdb);
    clock_t c2 = clock();
    for (int b = 0; b < BODIES; ++b) {
      (void)myb_spk_state(m->spk, m->body[b], 399, tdb, s, NULL);
    }
    clock_t c3 = clock();
    myb_gravity_field_set_time(m->field, tt);
    clock_t c4 = clock();
    myb_tide_raiser raisers[2] = {{m->gm_body[0], {1.4e11, 2e10, 1e10}},
                                  {m->gm_body[1], {3e8, 1e8, 1e8}}};
    const myb_tides_input tin = {.raisers = raisers, .count = 2, .t = ut1, .xp = eop.xp,
                                 .yp = eop.yp, .flags = MYB_TIDES_ALL};
    (void)myb_tides_apply(m->tides, m->field, &tin);
    clock_t c5 = clock();
    (void)myb_gravity_field_eval(m->field, m->work, x, 0, -1, NULL, a);
    clock_t c6 = clock();
    static const myb_planet_radiation knocke = {
        .radius = 6378137.0, .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
        .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18}, .omega = 1.99e-7};
    (void)myb_force_planet_radiation(&knocke, x, (const double[3]){1.4e11, 2e10, 1e10}, 1361.0,
                                     0.3, 1e-3, 10, a);
    clock_t c7 = clock();
    (void)acceleration(m, 100.0 * k, x, v, a);
    clock_t c8 = clock();
    t_rot += (double)(c1 - c0);
    t_tdb += (double)(c2 - c1);
    t_spk += (double)(c3 - c2);
    t_set += (double)(c4 - c3);
    t_tides += (double)(c5 - c4);
    t_eval += (double)(c6 - c5);
    t_alb += (double)(c7 - c6);
    t_all += (double)(c8 - c7);
  }
  const double us = 1e6 / (double)CLOCKS_PER_SEC / (double)REPS;
  printf("  profile (us per call): orientation %.1f, TT->TDB %.1f, 7 bodies %.1f, set_time %.1f, "
         "tides %.1f, field %.1f, albedo %.1f; whole acceleration %.1f\n",
         t_rot * us, t_tdb * us, t_spk * us, t_set * us, t_tides * us, t_eval * us, t_alb * us,
         t_all * us);
}

typedef struct observations {
  size_t count;
  double t[2000];       /* TT seconds from t0 */
  double r[2000][3];    /* GCRS positions, m */
  double rtn[2000][9];  /* rows R, T, N of the orbital frame at the observation */
} observations;

/* Propagates from the parameters p and fills the residuals obs - model (3 per epoch). */
static myb_status residuals(model* m, const observations* obs, const double* p, double* res,
                            long* steps) {
  m->cr = p[6];
  memcpy(m->empirical, p + 7, 5 * sizeof(double));
  myb_ias15* ias = NULL;
  const myb_ias15_options opt = {.epsilon = 1e-10};
  myb_status st = myb_ias15_create(3, &opt, NULL, &ias);
  if (st != MYB_OK) {
    return st;
  }
  st = myb_ias15_start(ias, acceleration, m, 0.0, p);
  size_t k = 0;
  /* bound: the observations, one accepted step per iteration */
  while (st == MYB_OK && k < obs->count) {
    const double before = myb_ias15_time(ias);
    if (obs->t[k] <= before) {
      double y[6];
      if (k == 0) {
        memcpy(y, p, 6 * sizeof(double));
      } else {
        st = myb_ias15_dense(ias, obs->t[k], y);
      }
      for (int i = 0; i < 3 && st == MYB_OK; ++i) {
        res[3 * k + (size_t)i] = obs->r[k][i] - y[i];
      }
      ++k;
      continue;
    }
    st = myb_ias15_step(ias, obs->t[obs->count - 1]);
  }
  if (steps != NULL) {
    *steps = (long)myb_ias15_statistics(ias).accepted;
  }
  myb_ias15_destroy(ias);
  return st;
}

/* Solves the symmetric positive definite system a x = b (n <= 12) by Cholesky. */
static int solve(int n, double a[PARAMETERS][PARAMETERS], double* b) {
  for (int j = 0; j < n; ++j) {
    double d = a[j][j];
    for (int k = 0; k < j; ++k) {
      d -= a[j][k] * a[j][k];
    }
    if (!(d > 0.0)) {
      return 0;
    }
    a[j][j] = sqrt(d);
    for (int i = j + 1; i < n; ++i) {
      double s = a[i][j];
      for (int k = 0; k < j; ++k) {
        s -= a[i][k] * a[j][k];
      }
      a[i][j] = s / a[j][j];
    }
  }
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    for (int k = 0; k < i; ++k) {
      s -= a[i][k] * b[k];
    }
    b[i] = s / a[i][i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    for (int k = i + 1; k < n; ++k) {
      s -= a[k][i] * b[k];
    }
    b[i] = s / a[i][i];
  }
  return 1;
}

typedef struct partial_job {
  model models[PARAMETERS];
  const observations* obs;
  const double* p;
  const double* step;
  const double* res;
  double (*jac)[6000];
  size_t nres;
  myb_status status[PARAMETERS];
  double scratch[PARAMETERS][6000];
} partial_job;

static void partial_range(void* context, size_t begin, size_t end) {
  partial_job* job = context;
  for (size_t j = begin; j < end; ++j) {
    double q[PARAMETERS];
    memcpy(q, job->p, sizeof q);
    q[j] += job->step[j];
    job->status[j] = residuals(&job->models[j], job->obs, q, job->scratch[j], NULL);
    for (size_t i = 0; i < job->nres; ++i) {
      job->jac[j][i] = -(job->scratch[j][i] - job->res[i]) / job->step[j]; /* d(model)/dp */
    }
  }
}

static void statistics(const observations* obs, const double* res, double* rms, double rtn[3]) {
  double s = 0.0;
  double c[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k < obs->count; ++k) {
    for (int i = 0; i < 3; ++i) {
      const double* row = obs->rtn[k] + 3 * i;
      const double q = row[0] * res[3 * k] + row[1] * res[3 * k + 1] + row[2] * res[3 * k + 2];
      c[i] += q * q;
      s += res[3 * k + (size_t)i] * res[3 * k + (size_t)i];
    }
  }
  *rms = sqrt(s / (double)obs->count);
  for (int i = 0; i < 3; ++i) {
    rtn[i] = sqrt(c[i] / (double)obs->count);
  }
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static char header[32768];
  static model m;
  static observations obs;
  myb_diag diag = {0};
  myb_sp3* sp3 = NULL;
  if (!read_file(DOWNLOAD "header.440", header, sizeof header) ||
      myb_sp3_load_file(MYB_TEST_DATA_DIR "/sp3/ilrsa.lageos1.250607.excerpt.sp3", NULL, &sp3, &diag) !=
          MYB_OK) {
    printf("header.440 not downloaded (or the SP3 excerpt missing): skipped\n");
    return MTEST_SKIP;
  }
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* kernels[] = {MYB_DATA_DIR "/time/time-scales.tk",
                           MYB_DATA_DIR "/time/tdb-tt-fb1990.tk", MYB_DATA_DIR "/time/utc-tai.tk",
                           MYB_DATA_DIR "/frames/iau2006-cip.tk",
                           MYB_DATA_DIR "/frames/eop-tides.tk",
                           MYB_DATA_DIR "/grav/iers2010-tides.tk"};
  for (size_t i = 0; i < sizeof kernels / sizeof kernels[0]; ++i) {
    MT_CHECK(myb_kpool_load_file(pool, kernels[i], &diag) == MYB_OK);
  }
  MT_CHECK(myb_timesys_create(pool, NULL, &m.sys, &diag) == MYB_OK);
  MT_CHECK(myb_timesys_scale(m.sys, "TT", &m.tt) == MYB_OK);
  MT_CHECK(myb_timesys_scale(m.sys, "TDB", &m.tdb) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &m.cip, &diag) == MYB_OK);
  MT_CHECK(myb_eop_create(NULL, &m.eop) == MYB_OK);
  MT_CHECK(myb_eop_set_tides(m.eop, pool, "IERS10_EOPT", &diag) == MYB_OK);
  MT_CHECK(myb_tides_create(pool, "IERS10_TIDE", NULL, &m.tides, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  MT_CHECK(myb_spk_create(NULL, &m.spk) == MYB_OK);
  const int ok =
      myb_eop_load_c04_file(m.eop, DOWNLOAD "eopc04.1962-now", &diag) == MYB_OK &&
      myb_spk_load_file(m.spk, DOWNLOAD "de440.bsp", &diag) == MYB_OK &&
      myb_gravity_field_load_icgem_file(DOWNLOAD "GOCO06s.gfc", DEGREE, NULL, &m.field, &diag) ==
          MYB_OK &&
      myb_tides_load_ocean_file(m.tides, DOWNLOAD "fes2004_Cnm-Snm.dat", OCEAN_DEGREE,
                                &diag) == MYB_OK;
  if (!ok) {
    printf("data/download incomplete (%s): skipped\n", diag.message);
    myb_sp3_destroy(sp3);
    return MTEST_SKIP;
  }
  MT_CHECK(myb_gravity_work_create(DEGREE, NULL, &m.work) == MYB_OK);
  m.gm_earth = myb_gravity_field_gm(m.field);
  /* GM of the perturbing bodies from DE440 (au^3/day^2) */
  double au = 0.0;
  double emrat = 0.0;
  MT_CHECK(read_constant(header, "AU", &au) && read_constant(header, "EMRAT", &emrat));
  const double unit = (au * 1e3) * (au * 1e3) * (au * 1e3) / (86400.0 * 86400.0);
  const char* names[BODIES] = {"GMS", "GMB", "GM1", "GM2", "GM4", "GM5", "GM6"};
  const int ids[BODIES] = {10, 301, 1, 2, 4, 5, 6};
  for (int b = 0; b < BODIES; ++b) {
    double g = 0.0;
    MT_CHECK(read_constant(header, names[b], &g));
    m.gm_body[b] = g * unit;
    m.body[b] = ids[b];
  }
  m.gm_body[1] /= 1.0 + emrat; /* Moon from the Earth-Moon barycentre */
  /* observations: every STRIDE-th epoch, ITRS (SLRF2020) -> GCRS */
  MT_CHECK(strcmp(myb_sp3_time_system(sp3), "UTC") == 0);
  const size_t epochs = myb_sp3_epoch_count(sp3);
  double p[PARAMETERS] = {0.0};
  for (size_t e = 0; e < epochs && obs.count < 2000; e += STRIDE) {
    myb_calendar label;
    myb_sp3_epoch(sp3, e, &label);
    myb_time t;
    MT_CHECK(myb_timesys_from_utc(m.sys, &label, m.tt, &t) == MYB_OK);
    if (e == 0) {
      m.t0 = t;
    }
    double itrs[6];
    const int flags = myb_sp3_state(sp3, 0, e, itrs, itrs + 3);
    if (!(flags & 1)) {
      continue;
    }
    myb_rotation rot;
    myb_eop_values eop;
    myb_time ut1;
    MT_CHECK(orientation(&m, t, &rot, &eop, &ut1) == MYB_OK);
    myb_rotation inv;
    myb_rotation_invert(&rot, &inv);
    double gcrs[6];
    myb_rotation_apply_state(&inv, itrs, gcrs);
    const size_t k = obs.count++;
    obs.t[k] = myb_dd_to_double(myb_time_diff(t, m.t0));
    memcpy(obs.r[k], gcrs, 3 * sizeof(double));
    double rot_rtn[9];
    double n[3];
    myb_vec3_cross(gcrs, gcrs + 3, n);
    const double rr = myb_vec3_norm(gcrs);
    const double nn = myb_vec3_norm(n);
    for (int i = 0; i < 3; ++i) {
      rot_rtn[i] = gcrs[i] / rr;
      rot_rtn[6 + i] = n[i] / nn;
    }
    myb_vec3_cross(rot_rtn + 6, rot_rtn, rot_rtn + 3);
    memcpy(obs.rtn[k], rot_rtn, sizeof rot_rtn);
    if (e == 0) {
      MT_CHECK(flags == 3);
      memcpy(p, gcrs, 6 * sizeof(double));
    }
  }
  myb_sp3_destroy(sp3);
  for (int k = 0; k < 200; ++k) {
    const myb_time tt = myb_time_add(m.t0, 3600.0 * (k - 1));
    myb_time tdb;
    MT_CHECK(myb_timesys_convert(m.sys, tt, m.tt, m.tdb, &tdb) == MYB_OK);
    m.tdb_tt[k] = myb_dd_to_double(myb_time_diff(tdb, tt));
  }
  printf("LAGEOS-1, ILRS combined orbit 2025-06-01 .. 06-07: %zu positions (every %d min)\n",
         obs.count, 6 * STRIDE);
  p[6] = 1.13;
  const int estimate_cr = getenv("MYB_LAGEOS_CR") != NULL;
  if (getenv("MYB_LAGEOS_OFF") != NULL) {
    m.off = atoi(getenv("MYB_LAGEOS_OFF"));
    printf("  model without the parts %d\n", m.off);
  }
  if (getenv("MYB_LAGEOS_PROFILE") != NULL) {
    profile(&m, p, p + 3);
  }
  static double res[6000];
  static double jac[PARAMETERS][6000];
  const double step[PARAMETERS] = {1.0,   1.0,   1.0,   1e-3,  1e-3,  1e-3,
                                   0.01,  1e-11, 1e-11, 1e-11, 1e-11, 1e-11};
  long steps = 0;
  MT_CHECK(residuals(&m, &obs, p, res, &steps) == MYB_OK);
  double rms = 0.0;
  double rtn[3];
  statistics(&obs, res, &rms, rtn);
  printf("  from the ILRS state at the first epoch, cr 1.13, no empirical: RMS %.3f m "
         "(R %.3f, T %.3f, N %.3f), %ld steps\n",
         rms, rtn[0], rtn[1], rtn[2], steps);
  const size_t nres = 3 * obs.count;
  for (int iteration = 0; iteration < 3; ++iteration) {
    if (iteration == 0) {
      /* the 12 propagations of the partials in parallel, each with its own field */
      static partial_job job;
      job.obs = &obs;
      job.p = p;
      job.step = step;
      job.res = res;
      job.jac = jac;
      job.nres = nres;
      int fields_ok = 1;
      for (int j = 0; j < PARAMETERS; ++j) {
        job.models[j] = m;
        job.models[j].cached = 0;
        fields_ok = fields_ok &&
                    myb_gravity_field_clone(m.field, NULL, &job.models[j].field) == MYB_OK &&
                    myb_gravity_work_create(DEGREE, NULL, &job.models[j].work) == MYB_OK;
      }
      MT_CHECK(fields_ok);
      myb_pool* pool_threads = NULL;
      MT_CHECK(myb_pool_create(0, NULL, &pool_threads) == MYB_OK);
      myb_internal_pool_for(pool_threads, PARAMETERS, 1, partial_range, &job);
      myb_pool_destroy(pool_threads);
      for (int j = 0; j < PARAMETERS; ++j) {
        MT_CHECK(job.status[j] == MYB_OK);
        myb_gravity_work_destroy(job.models[j].work);
        myb_gravity_field_destroy(job.models[j].field);
      }
    }
    /* normal equations in scaled parameters */
    /* cr is held at 1.13 as by the ILRS analysis centres: over a week a force of fixed
     * inertial direction is nearly a once-per-revolution signal in R, T, N, so cr and the
     * once-per-revolution terms cannot be separated */
    int free_index[PARAMETERS];
    int nfree = 0;
    for (int i = 0; i < PARAMETERS; ++i) {
      if (i != 6 || estimate_cr) {
        free_index[nfree++] = i;
      }
    }
    double n[PARAMETERS][PARAMETERS];
    double b[PARAMETERS];
    for (int i = 0; i < nfree; ++i) {
      const int fi = free_index[i];
      b[i] = myb_dot2(jac[fi], res, nres) * step[fi];
      for (int j = 0; j <= i; ++j) {
        const int fj = free_index[j];
        n[i][j] = myb_dot2(jac[fi], jac[fj], nres) * step[fi] * step[fj];
        n[j][i] = n[i][j];
      }
    }
    MT_CHECK(solve(nfree, n, b));
    for (int i = 0; i < nfree; ++i) {
      p[free_index[i]] += b[i] * step[free_index[i]];
    }
    MT_CHECK(residuals(&m, &obs, p, res, NULL) == MYB_OK);
    statistics(&obs, res, &rms, rtn);
    printf("  iteration %d: RMS %.4f m (R %.4f, T %.4f, N %.4f)\n", iteration + 1, rms, rtn[0],
           rtn[1], rtn[2]);
  }
  printf("  fitted: cr %.4f; along-track %.2f, %.2f cos u, %.2f sin u; cross-track %.2f cos u, "
         "%.2f sin u (pm/s^2)\n",
         p[6], p[7] * 1e12, p[8] * 1e12, p[9] * 1e12, p[10] * 1e12, p[11] * 1e12);
  /* measured: 2.0 cm (R 2.6 mm, T 1.8 cm, N 8 mm), constant along-track -1.4 pm/s^2 (the
   * thermal drag of LAGEOS); the once-per-revolution cross-track terms (~600 pm/s^2) absorb
   * what differs from the ILRS models (ocean tides EOT11a and their admittance, S1/S2
   * atmospheric tides) — the terms the analysis centres estimate for the same reason */
  if (m.off == 0 && !estimate_cr) {
    MT_CHECK(rms < 0.03);
    MT_CHECK(fabs(p[7]) < 5e-12);
  }
  myb_gravity_work_destroy(m.work);
  myb_gravity_field_destroy(m.field);
  myb_tides_destroy(m.tides);
  myb_spk_destroy(m.spk);
  myb_eop_destroy(m.eop);
  myb_cip_destroy(m.cip);
  myb_timesys_destroy(m.sys);
  return mt_finish("test_lageos");
}
