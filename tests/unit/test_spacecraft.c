/* SPDX-License-Identifier: BSD-3-Clause */
/* Phase 8 criterion: real spacecraft of JPL Horizons (tests/data/horizons), each fitted with
 * myb's force models and variational equations (Gauss-Newton on the initial state and the
 * parameters of non-gravitational forces, partials from myb/uncert.h):
 *  - the Tesla Roadster, heliocentric, 2018-03-20 .. 2025-03-20 (7 years, a Mars encounter at
 *    0.05 au in 2020): Sun with its Schwarzschild term, planets and Moon of DE440, the radial
 *    non-gravitational acceleration A1 (1 au / r)^2 of JPL's solution #11;
 *  - OSIRIS-REx through its Earth gravity assist of 2017-09-22 (perigee 17 000 km altitude):
 *    the Earth with J2..J4 about the CIP, Sun, Moon, planets, solar pressure with the Earth's
 *    penumbra; the perigee and the B-plane against Horizons;
 *  - JWST in its orbit about Sun-Earth L2: in the pulsating frame of the Sun and the Earth-Moon
 *    barycentre (myb/cr3bp.h) it stays near L2 at the Jacobi constant of a halo orbit; a week
 *    between station-keeping manoeuvres fitted with a solar pressure of fixed direction relative
 *    to the Sun line.
 * Needs data/download (de440.bsp, header.440); skipped otherwise. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"
#include "num/vec3.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

enum { MAX_OBS = 1000, MAX_BODIES = 12 };

static const double au_m = 149597870700.0;
static const double c_light = 299792458.0;

typedef struct track {
  size_t count;
  double t[MAX_OBS];  /* TDB seconds from J2000 */
  double y[MAX_OBS][6];
} track;

static int read_track(const char* path, track* tr) {
  FILE* f = fopen(path, "r");
  if (f == NULL) {
    return 0;
  }
  char line[512];
  tr->count = 0;
  while (fgets(line, sizeof line, f) != NULL && tr->count < MAX_OBS) {
    if (line[0] == '#') {
      continue;
    }
    double jd = 0.0;
    double v[6];
    char date[64];
    if (sscanf(line, "%lf, A.D. %63[^,], %lf, %lf, %lf, %lf, %lf, %lf", &jd, date, &v[0], &v[1],
               &v[2], &v[3], &v[4], &v[5]) == 8) {
      tr->t[tr->count] = (jd - 2451545.0) * 86400.0;
      for (int i = 0; i < 6; ++i) {
        tr->y[tr->count][i] = v[i] * 1e3;
      }
      tr->count++;
    }
  }
  fclose(f);
  return tr->count > 0;
}

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
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  int index = -1;
  for (int k = 0; k < (int)count; ++k) {
    while (*p == ' ' || *p == '\n' || *p == '\r') {
      ++p;
    }
    size_t n = 0;
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
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  for (int k = 0; k <= index; ++k) {
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

/* ---- Terms ----------------------------------------------------------------------------------- */

typedef struct point_mass {
  double gm;
} point_mass;

static myb_status central_accel(void* c, double t, const double x[3], const double v[3],
                                double a[3]) {
  (void)t;
  (void)v;
  const double r = myb_vec3_norm(x);
  const double k = -((const point_mass*)c)->gm / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * x[i];
  }
  return MYB_OK;
}

static myb_status central_jac(void* c, double t, const double x[3], const double v[3],
                              double dadx[9], double dadv[9]) {
  (void)t;
  (void)v;
  myb_point_mass_jacobian(((const point_mass*)c)->gm, x, dadx);
  memset(dadv, 0, 9 * sizeof(double));
  return MYB_OK;
}

static myb_status schwarzschild_accel(void* c, double t, const double x[3], const double v[3],
                                      double a[3]) {
  (void)t;
  const myb_ppn ppn = {c_light, 1.0, 1.0};
  return myb_rel_schwarzschild(&ppn, ((const point_mass*)c)->gm, x, v, a);
}

/* Perturbing bodies from the SPK, positions relative to `center` (0: barycentric, direct
 * terms only; otherwise direct minus indirect). */
typedef struct bodies {
  const myb_spk* spk;
  int center;
  int count;
  int id[MAX_BODIES];
  double gm[MAX_BODIES];
} bodies;

static myb_status bodies_positions(const bodies* b, double t, double s[MAX_BODIES][3]) {
  for (int k = 0; k < b->count; ++k) {
    double st[6];
    const myb_status r = myb_spk_state(b->spk, b->id[k], b->center, myb_time_from_seconds(t), st, NULL);
    if (r != MYB_OK) {
      return r;
    }
    memcpy(s[k], st, 3 * sizeof(double));
  }
  return MYB_OK;
}

static myb_status bodies_accel(void* c, double t, const double x[3], const double v[3],
                               double a[3]) {
  (void)v;
  const bodies* b = c;
  double s[MAX_BODIES][3];
  const myb_status st = bodies_positions(b, t, s);
  if (st != MYB_OK) {
    return st;
  }
  a[0] = a[1] = a[2] = 0.0;
  for (int k = 0; k < b->count; ++k) {
    if (b->center == 0) {
      const double d[3] = {s[k][0] - x[0], s[k][1] - x[1], s[k][2] - x[2]};
      const double dd = myb_vec3_norm(d);
      const double q = b->gm[k] / (dd * dd * dd);
      for (int i = 0; i < 3; ++i) {
        a[i] += q * d[i];
      }
    } else {
      myb_force_third_body(b->gm[k], s[k], x, a);
    }
  }
  return MYB_OK;
}

static myb_status bodies_jac(void* c, double t, const double x[3], const double v[3],
                             double dadx[9], double dadv[9]) {
  (void)v;
  const bodies* b = c;
  double s[MAX_BODIES][3];
  const myb_status st = bodies_positions(b, t, s);
  if (st != MYB_OK) {
    return st;
  }
  memset(dadx, 0, 9 * sizeof(double));
  memset(dadv, 0, 9 * sizeof(double));
  for (int k = 0; k < b->count; ++k) {
    const double r[3] = {x[0] - s[k][0], x[1] - s[k][1], x[2] - s[k][2]};
    double j[9];
    myb_point_mass_jacobian(b->gm[k], r, j);
    for (int i = 0; i < 9; ++i) {
      dadx[i] += j[i];
    }
  }
  return MYB_OK;
}

/* Radial acceleration A1 (1 au / r)^2 from the Sun at `sun` (NULL: the origin). */
typedef struct radial {
  double a1;
  const bodies* sun; /* body 0 of this set is the Sun (barycentric runs) */
  double shadow;     /* 1: Earth's penumbra (geocentric runs) */
  const bodies* geo; /* geocentric set whose body 0 is the Sun */
} radial;

static myb_status radial_accel(void* c, double t, const double x[3], const double v[3],
                               double a[3]) {
  (void)v;
  const radial* r = c;
  double d[3] = {x[0], x[1], x[2]};
  double nu = 1.0;
  if (r->geo != NULL) {
    double s[MAX_BODIES][3];
    const myb_status st = bodies_positions(r->geo, t, s);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      d[i] = x[i] - s[0][i];
    }
    const double src[3] = {-d[0], -d[1], -d[2]};
    const double earth[3] = {-x[0], -x[1], -x[2]};
    nu = myb_shadow_fraction(src, 6.957e8, earth, 6.378137e6, (const double[2]){0.6, 0.0});
  }
  const double rr = myb_vec3_norm(d);
  const double k = nu * r->a1 * (au_m / rr) * (au_m / rr) / rr;
  for (int i = 0; i < 3; ++i) {
    a[i] = k * d[i];
  }
  return MYB_OK;
}

/* Earth's zonal harmonics J2..J4 about the CIP pole. */
typedef struct zonal {
  const myb_cip* cip;
  myb_gravity_field* field;
  myb_gravity_work* work;
} zonal;

static myb_status zonal_accel(void* c, double t, const double x[3], const double v[3], double a[3]) {
  (void)v;
  const zonal* z = c;
  double xys[3];
  myb_cip_xys(z->cip, myb_time_from_seconds(t), xys, NULL);
  const double p[3] = {xys[0], xys[1], sqrt(1.0 - xys[0] * xys[0] - xys[1] * xys[1])};
  /* any frame with the pole as third axis */
  double e1[3] = {p[2], 0.0, -p[0]};
  const double n1 = myb_vec3_norm(e1);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= n1;
  }
  double e2[3];
  myb_vec3_cross(p, e1, e2);
  const double xb[3] = {myb_vec3_dot(e1, x), myb_vec3_dot(e2, x), myb_vec3_dot(p, x)};
  double ab[3];
  const myb_status st = myb_gravity_field_eval(z->field, z->work, xb, 2, -1, NULL, ab);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    a[i] = e1[i] * ab[0] + e2[i] * ab[1] + p[i] * ab[2];
  }
  return MYB_OK;
}

/* Solar pressure of fixed direction in the frame of the Sun line: R from the Sun, T = z x R,
 * N = R x T, components c[0..2] at 1 au, scaled by (1 au / r)^2. */
typedef struct sail {
  double c[3];
  const bodies* sun;
} sail;

static myb_status sail_accel(void* ctx, double t, const double x[3], const double v[3], double a[3]) {
  (void)v;
  const sail* s = ctx;
  double sp[MAX_BODIES][3];
  const myb_status st = bodies_positions(s->sun, t, sp);
  if (st != MYB_OK) {
    return st;
  }
  const double d[3] = {x[0] - sp[0][0], x[1] - sp[0][1], x[2] - sp[0][2]};
  const double rr = myb_vec3_norm(d);
  const double r[3] = {d[0] / rr, d[1] / rr, d[2] / rr};
  double tt[3] = {-r[1], r[0], 0.0};
  const double tn = myb_vec3_norm(tt);
  for (int i = 0; i < 3; ++i) {
    tt[i] /= tn;
  }
  double nn[3];
  myb_vec3_cross(r, tt, nn);
  const double k = (au_m / rr) * (au_m / rr);
  for (int i = 0; i < 3; ++i) {
    a[i] = k * (s->c[0] * r[i] + s->c[1] * tt[i] + s->c[2] * nn[i]);
  }
  return MYB_OK;
}

/* ---- Fitting ----------------------------------------------------------------------------------- */

typedef struct fit_result {
  double rms;
  double max;
  double y0[6];
} fit_result;

/* Gauss-Newton on (y0, parameters) over observations first .. last of the track. */
static myb_status fit(myb_force_model* m, const track* tr, size_t first, size_t last, int iterations,
                      fit_result* out) {
  const size_t np = myb_force_model_parameter_count(m);
  const size_t nu = 6 + np;
  double y0[6];
  memcpy(y0, tr->y[first], sizeof y0);
  const double t0 = tr->t[first];
  myb_status st = MYB_OK;
  for (int it = 0; it <= iterations && st == MYB_OK; ++it) {
    myb_variational* var = NULL;
    const myb_variational_options o = {.epsilon = 1e-11};
    st = myb_variational_create(m, &o, NULL, &var);
    if (st == MYB_OK) {
      st = myb_variational_start(var, t0, y0);
    }
    double n[22 * 22];
    double b[22];
    memset(n, 0, sizeof n);
    memset(b, 0, sizeof b);
    double s2 = 0.0;
    double mx = 0.0;
    size_t k = first;
    /* bound: the observations; one integrator step per iteration */
    while (st == MYB_OK && k <= last) {
      if (tr->t[k] <= myb_variational_time(var)) {
        double y[6];
        double phi[36];
        double sens[6 * 16];
        st = k == first ? MYB_OK : myb_variational_dense(var, tr->t[k], y, phi, sens);
        if (k == first) {
          memcpy(y, y0, sizeof y);
          memset(phi, 0, sizeof phi);
          for (int i = 0; i < 6; ++i) {
            phi[i * 6 + i] = 1.0;
          }
          memset(sens, 0, sizeof sens);
        }
        double h[3][22];
        for (int i = 0; i < 3; ++i) {
          const double res = tr->y[k][i] - y[i];
          s2 += res * res;
          for (size_t j = 0; j < 6; ++j) {
            h[i][j] = phi[(size_t)i * 6 + j];
          }
          for (size_t j = 0; j < np; ++j) {
            h[i][6 + j] = sens[(size_t)i * np + j];
          }
          for (size_t p = 0; p < nu; ++p) {
            b[p] += h[i][p] * res;
            for (size_t q = 0; q < nu; ++q) {
              n[p * nu + q] += h[i][p] * h[i][q];
            }
          }
        }
        const double dr = sqrt((tr->y[k][0] - y[0]) * (tr->y[k][0] - y[0]) +
                               (tr->y[k][1] - y[1]) * (tr->y[k][1] - y[1]) +
                               (tr->y[k][2] - y[2]) * (tr->y[k][2] - y[2]));
        mx = fmax(mx, dr);
        ++k;
        continue;
      }
      st = myb_variational_step(var, tr->t[last]);
    }
    myb_variational_destroy(var);
    if (st != MYB_OK) {
      break;
    }
    out->rms = sqrt(s2 / (double)(last - first + 1));
    out->max = mx;
    memcpy(out->y0, y0, sizeof y0);
    if (it == iterations) {
      break;
    }
    /* scaled normal equations (diagonal equilibration) */
    double d[22];
    for (size_t p = 0; p < nu; ++p) {
      d[p] = sqrt(n[p * nu + p]);
      if (!(d[p] > 0.0)) {
        return MYB_E_DOMAIN;
      }
    }
    for (size_t p = 0; p < nu; ++p) {
      b[p] /= d[p];
      for (size_t q = 0; q < nu; ++q) {
        n[p * nu + q] /= d[p] * d[q];
      }
    }
    st = myb_cholesky(nu, n);
    if (st != MYB_OK) {
      break;
    }
    myb_cholesky_solve(nu, n, b);
    for (size_t p = 0; p < 6; ++p) {
      y0[p] += b[p] / d[p];
    }
    for (size_t p = 0; p < np; ++p) {
      *myb_force_model_parameter(m, p) += b[6 + p] / d[6 + p];
    }
  }
  return st;
}

/* ---- Scenarios ---------------------------------------------------------------------------------- */

typedef struct ephemeris {
  myb_spk* spk;
  char header[32768];
  double au_km;
  double emrat;
  double unit;
} ephemeris;

static double gm_of(const ephemeris* e, const char* name) {
  double g = 0.0;
  MT_CHECK(read_constant(e->header, name, &g));
  return g * e->unit;
}

static void scenario_roadster(const ephemeris* e) {
  track tr;
  MT_CHECK(read_track(MYB_TEST_DATA_DIR "/horizons/roadster.txt", &tr));
  point_mass sun = {gm_of(e, "GMS")};
  bodies planets = {.spk = e->spk, .center = 10};
  static const char* names[9] = {"GM1", "GM2", "GMB", "GMB", "GM4", "GM5", "GM6", "GM7", "GM8"};
  static const int ids[9] = {1, 2, 399, 301, 4, 5, 6, 7, 8};
  for (int k = 0; k < 9; ++k) {
    planets.id[k] = ids[k];
    planets.gm[k] = gm_of(e, names[k]);
  }
  planets.gm[2] *= e->emrat / (1.0 + e->emrat);
  planets.gm[3] /= 1.0 + e->emrat;
  planets.count = 9;
  radial nongrav = {.a1 = 0.0};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "sun", central_accel, central_jac, &sun, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "sun 1PN", schwarzschild_accel, NULL, &sun, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "planets", bodies_accel, bodies_jac, &planets, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "A1", radial_accel, NULL, &nongrav, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "A1", &nongrav.a1, 1e-12, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_set_scales(m, au_m, 3e4) == MYB_OK);
  fit_result r;
  const myb_status st = fit(m, &tr, 0, tr.count - 1, 4, &r);
  MT_CHECK(st == MYB_OK);
  printf("Tesla Roadster, 2018-03 .. 2025-03 (512 positions, heliocentric): RMS %.3f km, max "
         "%.3f km; A1 = %.3e m/s^2 at 1 au\n",
         r.rms / 1e3, r.max / 1e3, nongrav.a1);
  MT_CHECK(r.rms < 600.0 && r.max < 2500.0); /* measured 306 m, 1.34 km */
  MT_CHECK(nongrav.a1 > 1e-8 && nongrav.a1 < 2e-7);
  myb_force_model_destroy(m);
}

static void scenario_osiris_rex(const ephemeris* e, const myb_cip* cip) {
  track tr;
  MT_CHECK(read_track(MYB_TEST_DATA_DIR "/horizons/osiris-rex-earth-flyby.txt", &tr));
  const double gm_e = gm_of(e, "GMB") * e->emrat / (1.0 + e->emrat);
  point_mass earth = {gm_e};
  bodies third = {.spk = e->spk, .center = 399};
  static const char* names[5] = {"GMS", "GMB", "GM2", "GM4", "GM5"};
  static const int ids[5] = {10, 301, 2, 4, 5};
  for (int k = 0; k < 5; ++k) {
    third.id[k] = ids[k];
    third.gm[k] = gm_of(e, names[k]);
  }
  third.gm[1] /= 1.0 + e->emrat;
  third.count = 5;
  double re = 0.0;
  double j[3];
  MT_CHECK(read_constant(e->header, "RE", &re) && read_constant(e->header, "J2E", &j[0]) &&
           read_constant(e->header, "J3E", &j[1]) && read_constant(e->header, "J4E", &j[2]));
  double c[15] = {1.0};
  double s[15] = {0.0};
  for (int n = 2; n <= 4; ++n) {
    c[n * (n + 1) / 2] = -j[n - 2] / sqrt(2.0 * n + 1.0);
  }
  zonal z = {.cip = cip};
  MT_CHECK(myb_gravity_field_create(gm_e, re * 1e3, 4, c, s, NULL, &z.field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(4, NULL, &z.work) == MYB_OK);
  radial srp = {.a1 = 5e-8, .geo = &third};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "earth", central_accel, central_jac, &earth, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "J2..J4", zonal_accel, NULL, &z, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "earth 1PN", schwarzschild_accel, NULL, &earth, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "third bodies", bodies_accel, bodies_jac, &third, NULL) ==
           MYB_OK);
  MT_CHECK(myb_force_model_add(m, "pressure", radial_accel, NULL, &srp, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "pressure", &srp.a1, 1e-10, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_set_scales(m, 1e7, 1e3) == MYB_OK);
  fit_result r;
  MT_CHECK(fit(m, &tr, 0, tr.count - 1, 4, &r) == MYB_OK);
  printf("OSIRIS-REx, Earth gravity assist 2017-09-19 .. 25 (865 positions, geocentric): RMS %.1f m, "
         "max %.1f m; pressure %.3e m/s^2 at 1 au\n",
         r.rms, r.max, srp.a1);
  /* perigee: Horizons (bisection on r.v between samples, Hermite) vs the fitted trajectory */
  size_t k = 0;
  for (size_t i = 0; i + 1 < tr.count; ++i) {
    const double* a = tr.y[i];
    const double* b = tr.y[i + 1];
    if (a[0] * a[3] + a[1] * a[4] + a[2] * a[5] < 0.0 && b[0] * b[3] + b[1] * b[4] + b[2] * b[5] >= 0.0) {
      k = i;
    }
  }
  /* the fitted trajectory near the perigee */
  /* steps of at most the sample interval: from the start of the bracket one step covers it,
   * and the bisection runs on the dense output of that step */
  myb_variational* var = NULL;
  const myb_variational_options vo = {.epsilon = 1e-11, .h_max = tr.t[k + 1] - tr.t[k]};
  MT_CHECK(myb_variational_create(m, &vo, NULL, &var) == MYB_OK);
  MT_CHECK(myb_variational_start(var, tr.t[0], r.y0) == MYB_OK);
  MT_CHECK(myb_variational_integrate(var, tr.t[k]) == MYB_OK);
  /* step until r.v changes sign; the bisection then runs within that last step */
  double lo = tr.t[k];
  double hi = tr.t[k + 1];
  for (int guard = 0; guard < 10000; ++guard) {
    const double t_before = myb_variational_time(var);
    MT_CHECK(myb_variational_step(var, tr.t[k + 1] + 3600.0) == MYB_OK);
    double ys[6];
    myb_variational_get(var, ys, NULL, NULL);
    if (ys[0] * ys[3] + ys[1] * ys[4] + ys[2] * ys[5] >= 0.0) {
      lo = t_before;
      hi = myb_variational_time(var);
      break;
    }
  }
  double ym[6];
  for (int it = 0; it < 80; ++it) {
    const double mid = 0.5 * (lo + hi);
    MT_CHECK(myb_variational_dense(var, mid, ym, NULL, NULL) == MYB_OK);
    if (ym[0] * ym[3] + ym[1] * ym[4] + ym[2] * ym[5] < 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  const double t_ours = 0.5 * (lo + hi);
  MT_CHECK(myb_variational_dense(var, t_ours, ym, NULL, NULL) == MYB_OK);
  myb_variational_destroy(var);
  /* Horizons perigee by cubic Hermite interpolation between the bracketing samples */
  const double h = tr.t[k + 1] - tr.t[k];
  double best_t = tr.t[k];
  double best_r = 1e300;
  for (int q = 0; q <= 6000; ++q) {
    const double u = q / 6000.0;
    const double h00 = 2 * u * u * u - 3 * u * u + 1;
    const double h10 = u * u * u - 2 * u * u + u;
    const double h01 = -2 * u * u * u + 3 * u * u;
    const double h11 = u * u * u - u * u;
    double p[3];
    for (int i = 0; i < 3; ++i) {
      p[i] = h00 * tr.y[k][i] + h10 * h * tr.y[k][3 + i] + h01 * tr.y[k + 1][i] +
             h11 * h * tr.y[k + 1][3 + i];
    }
    const double rr = myb_vec3_norm(p);
    if (rr < best_r) {
      best_r = rr;
      best_t = tr.t[k] + u * h;
    }
  }
  const double r_ours = myb_vec3_norm(ym);
  printf("  perigee: radius %.3f km (Horizons %.3f km, altitude %.0f km), time differs by %.2f s\n",
         r_ours / 1e3, best_r / 1e3, (best_r - 6378137.0) / 1e3, t_ours - best_t);
  MT_CHECK(fabs(r_ours - best_r) < 150.0 && fabs(t_ours - best_t) < 0.5); /* 54 m, 0.17 s */
  /* B-plane one day before the perigee: Horizons state vs ours */
  size_t d1 = 0;
  while (d1 + 1 < tr.count && tr.t[d1] < best_t - 86400.0) {
    ++d1;
  }
  myb_variational* v2 = NULL;
  double yd[6];
  MT_CHECK(myb_variational_create(m, NULL, NULL, &v2) == MYB_OK);
  MT_CHECK(myb_variational_start(v2, tr.t[0], r.y0) == MYB_OK);
  MT_CHECK(myb_variational_integrate(v2, tr.t[d1]) == MYB_OK);
  myb_variational_get(v2, yd, NULL, NULL);
  myb_variational_destroy(v2);
  const double pole[3] = {0.0, -0.3977771559319137, 0.9174820620691818}; /* ecliptic J2000 */
  myb_bplane bh;
  myb_bplane bo;
  MT_CHECK(myb_bplane_compute(gm_e, tr.y[d1], tr.y[d1] + 3, pole, &bh) == MYB_OK);
  MT_CHECK(myb_bplane_compute(gm_e, yd, yd + 3, pole, &bo) == MYB_OK);
  printf("  B-plane a day before (ecliptic): B.T %.3f / %.3f km, B.R %.3f / %.3f km (ours / "
         "Horizons), v_inf %.4f km/s\n",
         bo.b_t / 1e3, bh.b_t / 1e3, bo.b_r / 1e3, bh.b_r / 1e3, bh.v_infinity / 1e3);
  MT_CHECK(fabs(bo.b_t - bh.b_t) < 100.0 && fabs(bo.b_r - bh.b_r) < 100.0); /* 32, 30 m */
  MT_CHECK(r.rms < 30.0); /* measured 12.6 m */
  myb_force_model_destroy(m);
  myb_gravity_work_destroy(z.work);
  myb_gravity_field_destroy(z.field);
}

static void scenario_jwst(const ephemeris* e) {
  track tr;
  MT_CHECK(read_track(MYB_TEST_DATA_DIR "/horizons/jwst.txt", &tr));
  /* the pulsating frame of the Sun and the Earth-Moon barycentre */
  const double gm_s = gm_of(e, "GMS");
  const double gm_b = gm_of(e, "GMB");
  const double mu = gm_b / (gm_s + gm_b);
  double l[5][3];
  MT_CHECK(myb_cr3bp_lagrange(mu, l) == MYB_OK);
  double xmin = 1e9;
  double xmax = -1e9;
  double ymax = 0.0;
  double zmax = 0.0;
  double cmin = 1e9;
  double cmax = -1e9;
  for (size_t k = 0; k < tr.count; ++k) {
    double sun[6];
    double emb[6];
    MT_CHECK(myb_spk_state(e->spk, 10, 0, myb_time_from_seconds(tr.t[k]), sun, NULL) == MYB_OK);
    MT_CHECK(myb_spk_state(e->spk, 3, 0, myb_time_from_seconds(tr.t[k]), emb, NULL) == MYB_OK);
    double rot[6];
    MT_CHECK(myb_cr3bp_from_ephemeris(sun, emb, gm_s, gm_b, tr.y[k], rot) == MYB_OK);
    xmin = fmin(xmin, rot[0]);
    xmax = fmax(xmax, rot[0]);
    ymax = fmax(ymax, fabs(rot[1]));
    zmax = fmax(zmax, fabs(rot[2]));
    const double cj = myb_cr3bp_jacobi(mu, rot);
    cmin = fmin(cmin, cj);
    cmax = fmax(cmax, cj);
  }
  const double c_l2 = myb_cr3bp_jacobi(mu, (const double[6]){l[1][0], 0, 0, 0, 0, 0});
  printf("JWST 2023-01 .. 03 in the Sun-EMB rotating frame: x %.5f .. %.5f (L2 %.5f), |y| < %.5f, "
         "|z| < %.5f; Jacobi %.6f .. %.6f (L2 %.6f)\n",
         xmin, xmax, l[1][0], ymax, zmax, cmin, cmax, c_l2);
  MT_CHECK(xmin > l[1][0] - 0.004 && xmax < l[1][0] + 0.004 && ymax < 0.01 && zmax < 0.01);
  MT_CHECK(cmax < c_l2 && cmin > c_l2 - 0.002); /* energy just above L2: libration orbits */
  /* a week fitted in the full model, 2023-02-01 .. 02-08 */
  bodies all = {.spk = e->spk, .center = 0};
  static const char* names[11] = {"GMS", "GM1", "GM2", "GMB", "GMB", "GM4", "GM5", "GM6", "GM7", "GM8", "GM9"};
  static const int ids[11] = {10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9};
  for (int k = 0; k < 11; ++k) {
    all.id[k] = ids[k];
    all.gm[k] = gm_of(e, names[k]);
  }
  all.gm[3] *= e->emrat / (1.0 + e->emrat);
  all.gm[4] /= 1.0 + e->emrat;
  all.count = 11;
  sail press = {.c = {1e-7, 0.0, 0.0}, .sun = &all};
  myb_force_model* m = NULL;
  MT_CHECK(myb_force_model_create(NULL, &m) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "bodies", bodies_accel, bodies_jac, &all, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(m, "pressure", sail_accel, NULL, &press, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "radial", &press.c[0], 1e-10, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "transverse", &press.c[1], 1e-10, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add_parameter(m, "normal", &press.c[2], 1e-10, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_set_scales(m, au_m, 3e4) == MYB_OK);
  size_t first = 0;
  while (first < tr.count && tr.t[first] < tr.t[0] + 31.0 * 86400.0) {
    ++first;
  }
  size_t last = first;
  while (last + 1 < tr.count && tr.t[last + 1] <= tr.t[first] + 7.0 * 86400.0) {
    ++last;
  }
  fit_result r;
  MT_CHECK(fit(m, &tr, first, last, 4, &r) == MYB_OK);
  printf("  a week (2023-02-01 .. 08, %zu positions) in the ephemeris model: RMS %.1f m, max %.1f m; "
         "pressure %.3e, %.3e, %.3e m/s^2 at 1 au\n",
         last - first + 1, r.rms, r.max, press.c[0], press.c[1], press.c[2]);
  MT_CHECK(r.rms < 500.0);
  MT_CHECK(press.c[0] > 0.0);
  myb_force_model_destroy(m);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static ephemeris e;
  MT_CHECK(myb_spk_create(NULL, &e.spk) == MYB_OK);
  FILE* f = fopen(DOWNLOAD "header.440", "rb");
  size_t n = 0;
  if (f != NULL) {
    n = fread(e.header, 1, sizeof e.header - 1, f);
    fclose(f);
  }
  e.header[n] = '\0';
  if (n == 0 || myb_spk_load_file(e.spk, DOWNLOAD "de440.bsp", NULL) != MYB_OK) {
    printf("de440.bsp or header.440 not downloaded: skipped\n");
    myb_spk_destroy(e.spk);
    return MTEST_SKIP;
  }
  MT_CHECK(read_constant(e.header, "AU", &e.au_km) && read_constant(e.header, "EMRAT", &e.emrat));
  e.unit = (e.au_km * 1e3) * (e.au_km * 1e3) * (e.au_km * 1e3) / (86400.0 * 86400.0);
  myb_kpool* pool = NULL;
  myb_cip* cip = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau2006-cip.tk", NULL) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &cip, NULL) == MYB_OK);
  myb_kpool_destroy(pool);
  scenario_jwst(&e);
  scenario_osiris_rex(&e, cip);
  scenario_roadster(&e);
  myb_cip_destroy(cip);
  myb_spk_destroy(e.spk);
  return mt_finish("test_spacecraft");
}
