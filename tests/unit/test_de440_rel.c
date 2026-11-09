/* SPDX-License-Identifier: BSD-3-Clause */
/* Phase 6 criterion: the Sun, the planets (barycentres of the planetary systems beyond the
 * Earth), the Earth and the Moon as point masses, integrated with IAS15 from the DE440 states
 * at J2000 with the GM and c of DE440 (header.440), compared with DE440 itself after 30, 365
 * and 3652 days — with the EIH equations (as DE440, beta = gamma = 1) and with Newton's; and
 * (phase 7) with EIH and the figures: the Earth's J2..J5 about the CIP, the Moon's harmonics to
 * degree 4 in the principal-axis frame of DE440 (moon_pa_de440_200625.bpc), the Sun's J2 about
 * its IAU pole, the Earth-Moon quadrupole-quadrupole coupling, all from header.440. DE440 also
 * models the tides of the Earth and the Moon and 373 small bodies; their absence bounds the
 * agreement. Needs data/download (fetch-data); skipped otherwise. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/dyn.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/integ.h"
#include "myb/io.h"
#include "myb/kin.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

enum { BODIES = 11 };

static const char* names[BODIES] = {"Sun",    "Mercury", "Venus",  "Earth",   "Moon",  "Mars",
                                    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
static const int32_t naif[BODIES] = {10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9};

/* Reads the constants of header.440 (GROUP 1040 names, GROUP 1041 values). */
static int read_constant(const char* text, const char* name, double* value) {
  const char* g1040 = strstr(text, "GROUP   1040");
  const char* g1041 = strstr(text, "GROUP   1041");
  if (g1040 == NULL || g1041 == NULL) {
    return 0;
  }
  /* index of the name */
  const char* p = g1040 + 12;
  double count = 0.0;
  size_t used = 0;
  /* bound: skips white space before the count */
  while (*p == ' ' || *p == '\n' || *p == '\r') {
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  int index = -1;
  for (int k = 0; k < (int)count && p < g1041; ++k) {
    /* bound: white space */
    while (*p == ' ' || *p == '\n' || *p == '\r') {
      ++p;
    }
    size_t n = 0;
    /* bound: a name */
    while (p[n] != ' ' && p[n] != '\n' && p[n] != '\r' && p[n] != '\0') {
      ++n;
    }
    if (n == strlen(name) && strncmp(p, name, n) == 0) {
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

static void state_at(const myb_spk* spk, double t, double* x, double* v) {
  for (int b = 0; b < BODIES; ++b) {
    double s[6];
    MT_CHECK(myb_spk_state(spk, naif[b], 0, myb_time_from_seconds(t), s, NULL) == MYB_OK);
    memcpy(x + 3 * b, s, 3 * sizeof(double));
    memcpy(v + 3 * b, s + 3, 3 * sizeof(double));
  }
}

/* ---- Figures (phase 7): the Earth's zonal harmonics J2..J5 about the CIP, the Moon's degree
 * 2..4 in its principal-axis frame of DE440 (moon_pa_de440_200625.bpc), the Sun's J2 about the
 * IAU pole, each acting on the point masses of the others with the reaction on its own centre,
 * and the quadrupole-quadrupole coupling of the Earth and the Moon. ------------------------ */

enum { EARTH = 3, MOON = 4, MOON_DEGREE = 4 };

typedef struct figures {
  myb_nbody* nbody;
  myb_cip* cip;
  myb_pck* pck;
  myb_gravity_field* earth;
  myb_gravity_field* moon;
  myb_gravity_field* sun;
  myb_gravity_work* work;
  double gm[BODIES];
  double gq_earth[9]; /* body frame */
  double gq_moon[9];
  double sun_pole[3];
} figures;

/* Rows of a rotation whose third row is the unit vector p (zonal fields need no more). */
static void pole_frame(const double p[3], double m[9]) {
  const double h[3] = {fabs(p[0]) < 0.9 ? 1.0 : 0.0, fabs(p[0]) < 0.9 ? 0.0 : 1.0, 0.0};
  double e1[3] = {h[1] * p[2] - h[2] * p[1], h[2] * p[0] - h[0] * p[2], h[0] * p[1] - h[1] * p[0]};
  const double n = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= n;
  }
  const double e2[3] = {p[1] * e1[2] - p[2] * e1[1], p[2] * e1[0] - p[0] * e1[2],
                        p[0] * e1[1] - p[1] * e1[0]};
  for (int i = 0; i < 3; ++i) {
    m[i] = e1[i];
    m[3 + i] = e2[i];
    m[6 + i] = p[i];
  }
}

/* Field of body `owner` (rotation m: inertial -> body) on every body in `targets`, with the
 * reaction on the owner. */
static myb_status apply_field(figures* f, myb_gravity_field* field, const double m[9], int owner,
                              const int* targets, int count, const double* x, double* a) {
  for (int k = 0; k < count; ++k) {
    const int b = targets[k];
    const double d[3] = {x[3 * b] - x[3 * owner], x[3 * b + 1] - x[3 * owner + 1],
                         x[3 * b + 2] - x[3 * owner + 2]};
    double db[3];
    for (int i = 0; i < 3; ++i) {
      db[i] = m[3 * i] * d[0] + m[3 * i + 1] * d[1] + m[3 * i + 2] * d[2];
    }
    double ab[3];
    const myb_status st = myb_gravity_field_eval(field, f->work, db, 2, -1, NULL, ab);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 3; ++i) {
      const double ai = m[i] * ab[0] + m[3 + i] * ab[1] + m[6 + i] * ab[2];
      a[3 * b + i] += ai;
      a[3 * owner + i] -= f->gm[b] / f->gm[owner] * ai;
    }
  }
  return MYB_OK;
}

static myb_status figure_acceleration(void* context, double t, const double* x, const double* v,
                                      double* a) {
  figures* f = context;
  myb_status st = myb_nbody_acceleration_xv(f->nbody, t, x, v, a);
  if (st != MYB_OK) {
    return st;
  }
  const myb_time tdb = myb_time_from_seconds(t);
  /* Earth: pole of the CIP (TT ~ TDB to 2 ms) */
  double xys[3];
  myb_cip_xys(f->cip, tdb, xys, NULL);
  const double pole[3] = {xys[0], xys[1], sqrt(1.0 - xys[0] * xys[0] - xys[1] * xys[1])};
  double me[9];
  pole_frame(pole, me);
  static const int earth_targets[2] = {0, MOON};
  st = apply_field(f, f->earth, me, EARTH, earth_targets, 2, x, a);
  if (st != MYB_OK) {
    return st;
  }
  /* Moon: principal axes of DE440 */
  myb_rotation pa;
  st = myb_pck_rotation(f->pck, 31008, tdb, &pa, NULL);
  if (st != MYB_OK) {
    return st;
  }
  static const int moon_targets[2] = {0, EARTH};
  st = apply_field(f, f->moon, pa.m, MOON, moon_targets, 2, x, a);
  if (st != MYB_OK) {
    return st;
  }
  /* Sun: J2 on every other body */
  double ms[9];
  pole_frame(f->sun_pole, ms);
  static const int sun_targets[BODIES - 1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  st = apply_field(f, f->sun, ms, 0, sun_targets, BODIES - 1, x, a);
  if (st != MYB_OK) {
    return st;
  }
  /* Earth-Moon quadrupole-quadrupole */
  double qe[9];
  double qm[9];
  myb_quadrupole_rotate(me, f->gq_earth, qe);
  myb_quadrupole_rotate(pa.m, f->gq_moon, qm);
  const double r[3] = {x[3 * MOON] - x[3 * EARTH], x[3 * MOON + 1] - x[3 * EARTH + 1],
                       x[3 * MOON + 2] - x[3 * EARTH + 2]};
  return myb_figure_figure(f->gm[EARTH], qe, f->gm[MOON], qm, r, a + 3 * EARTH, a + 3 * MOON);
}

static double factorial_ratio(int n, int m) { /* (n - m)! / (n + m)! */
  double r = 1.0;
  for (int k = n - m + 1; k <= n + m; ++k) {
    r /= k;
  }
  return r;
}

/* Normalised coefficients from the unnormalised ones of header.440 (J_n = -C_n0). */
static int setup_figures(figures* f, const char* header, double au) {
  (void)au;
  double re = 0.0;
  double am = 0.0;
  double asun = 0.0;
  double j2sun = 0.0;
  if (!read_constant(header, "RE", &re) || !read_constant(header, "AM", &am) ||
      !read_constant(header, "ASUN", &asun) || !read_constant(header, "J2SUN", &j2sun)) {
    return 0;
  }
  double ce[21] = {0.0};
  double se[21] = {0.0};
  ce[0] = 1.0;
  static const char* jn[4] = {"J2E", "J3E", "J4E", "J5E"};
  for (int n = 2; n <= 5; ++n) {
    double j = 0.0;
    if (!read_constant(header, jn[n - 2], &j)) {
      return 0;
    }
    ce[n * (n + 1) / 2] = -j / sqrt(2.0 * n + 1.0);
  }
  double cm[15] = {0.0};
  double sm[15] = {0.0};
  cm[0] = 1.0;
  for (int n = 2; n <= MOON_DEGREE; ++n) {
    for (int m = 0; m <= n; ++m) {
      char name[16];
      double c = 0.0;
      double s = 0.0;
      if (m == 0) {
        snprintf(name, sizeof name, "J%dM", n);
        if (!read_constant(header, name, &c)) {
          return 0;
        }
        c = -c;
      } else {
        snprintf(name, sizeof name, "C%d%dM", n, m);
        if (!read_constant(header, name, &c)) {
          return 0;
        }
        snprintf(name, sizeof name, "S%d%dM", n, m);
        if (!read_constant(header, name, &s)) {
          return 0;
        }
      }
      const double norm = sqrt((m == 0 ? 1.0 : 2.0) * (2.0 * n + 1.0) * factorial_ratio(n, m));
      cm[n * (n + 1) / 2 + m] = c / norm;
      sm[n * (n + 1) / 2 + m] = s / norm;
    }
  }
  double cs[6] = {1.0, 0.0, 0.0, -j2sun / sqrt(5.0), 0.0, 0.0};
  double ss[6] = {0.0};
  if (myb_gravity_field_create(f->gm[EARTH], re * 1e3, 5, ce, se, NULL, &f->earth) != MYB_OK ||
      myb_gravity_field_create(f->gm[MOON], am * 1e3, MOON_DEGREE, cm, sm, NULL, &f->moon) !=
          MYB_OK ||
      myb_gravity_field_create(f->gm[0], asun * 1e3, 2, cs, ss, NULL, &f->sun) != MYB_OK ||
      myb_gravity_work_create(5, NULL, &f->work) != MYB_OK) {
    return 0;
  }
  const double c2e[3] = {ce[3], 0.0, 0.0};
  const double s20[2] = {0.0, 0.0};
  myb_quadrupole_from_coefficients(f->gm[EARTH], re * 1e3, c2e, s20, f->gq_earth);
  const double c2m[3] = {cm[3], cm[4], cm[5]};
  const double s2m[2] = {sm[4], sm[5]};
  myb_quadrupole_from_coefficients(f->gm[MOON], am * 1e3, c2m, s2m, f->gq_moon);
  /* the Sun's pole, IAU: alpha 286.13 deg, delta 63.87 deg (ICRF) */
  const double ra = 286.13 * 3.141592653589793 / 180.0;
  const double de = 63.87 * 3.141592653589793 / 180.0;
  f->sun_pole[0] = cos(de) * cos(ra);
  f->sun_pole[1] = cos(de) * sin(ra);
  f->sun_pole[2] = sin(de);
  return 1;
}

static void release_figures(figures* f) {
  myb_gravity_work_destroy(f->work);
  myb_gravity_field_destroy(f->earth);
  myb_gravity_field_destroy(f->moon);
  myb_gravity_field_destroy(f->sun);
  myb_pck_destroy(f->pck);
  myb_cip_destroy(f->cip);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static char header[32768];
  myb_spk* spk = NULL;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  if (!read_file(DOWNLOAD "header.440", header, sizeof header) ||
      myb_spk_load_file(spk, DOWNLOAD "de440.bsp", NULL) != MYB_OK) {
    printf("de440.bsp or header.440 not downloaded: skipped\n");
    myb_spk_destroy(spk);
    return MTEST_SKIP;
  }
  double au = 0.0;
  double clight = 0.0;
  double emrat = 0.0;
  double g[11];
  const char* gm_names[] = {"GMS", "GM1", "GM2", "GMB", "GM4", "GM5", "GM6", "GM7", "GM8",
                            "GM9"};
  int ok = read_constant(header, "AU", &au) && read_constant(header, "CLIGHT", &clight) &&
           read_constant(header, "EMRAT", &emrat);
  for (int k = 0; k < 10; ++k) {
    ok = ok && read_constant(header, gm_names[k], &g[k]);
  }
  MT_CHECK(ok && fabs(au - 149597870.7) < 1e-6 && clight == 299792.458);
  if (!ok) {
    myb_spk_destroy(spk);
    return mt_finish("test_de440_rel");
  }
  /* au^3/day^2 -> m^3/s^2 */
  const double unit = (au * 1e3) * (au * 1e3) * (au * 1e3) / (86400.0 * 86400.0);
  double gm[BODIES];
  gm[0] = g[0] * unit;
  gm[1] = g[1] * unit;
  gm[2] = g[2] * unit;
  gm[3] = g[3] * emrat / (1.0 + emrat) * unit;
  gm[4] = g[3] / (1.0 + emrat) * unit;
  for (int k = 5; k < BODIES; ++k) {
    gm[k] = g[k - 1] * unit;
  }
  const double t0 = 0.0; /* J2000 TDB */
  double x0[3 * BODIES];
  double v0[3 * BODIES];
  state_at(spk, t0, x0, v0);
  const double days[3] = {30.0, 365.25, 3652.5};
  double worst[3][3][BODIES] = {{{0}}};
  /* model 2: EIH and the figures */
  static figures fig;
  memcpy(fig.gm, gm, sizeof fig.gm);
  int have_figures = setup_figures(&fig, header, au);
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  have_figures = have_figures &&
                 myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau2006-cip.tk", NULL) == MYB_OK &&
                 myb_cip_create(pool, "IAU2006A", NULL, &fig.cip, NULL) == MYB_OK &&
                 myb_pck_create(NULL, &fig.pck) == MYB_OK &&
                 myb_pck_load_file(fig.pck, DOWNLOAD "moon_pa_de440_200625.bpc", NULL) == MYB_OK;
  myb_kpool_destroy(pool);
  MT_CHECK(have_figures);
  for (int model = 0; model < (have_figures ? 3 : 2); ++model) {
    myb_nbody* s = NULL;
    MT_CHECK(myb_nbody_create(BODIES, 0, gm, NULL, &s) == MYB_OK);
    if (model >= 1) {
      const myb_ppn ppn = {clight * 1e3, 1.0, 1.0};
      MT_CHECK(myb_nbody_set_relativity(s, &ppn) == MYB_OK);
    }
    const myb_ias15_options o = {.epsilon = 1e-10, .block = 3};
    myb_ias15* in = NULL;
    MT_CHECK(myb_ias15_create(3 * BODIES, &o, NULL, &in) == MYB_OK);
    double y[6 * BODIES];
    memcpy(y, x0, sizeof x0);
    memcpy(y + 3 * BODIES, v0, sizeof v0);
    fig.nbody = s;
    MT_CHECK(myb_ias15_start(in, model == 2 ? figure_acceleration : myb_nbody_acceleration_xv,
                             model == 2 ? (void*)&fig : (void*)s, t0, y) == MYB_OK);
    for (int k = 0; k < 3; ++k) {
      const double t = t0 + days[k] * 86400.0;
      MT_CHECK(myb_ias15_integrate(in, t) == MYB_OK);
      const double* z = myb_ias15_state(in);
      double xr[3 * BODIES];
      double vr[3 * BODIES];
      state_at(spk, t, xr, vr);
      /* barycentric for the Sun, heliocentric for the planets and the Earth, geocentric for
       * the Moon: the barycentre itself drifts by the momentum of the missing small bodies */
      for (int b = 0; b < BODIES; ++b) {
        const int origin = b == 0 ? -1 : (b == 4 ? 3 : 0);
        double d[3];
        for (int c = 0; c < 3; ++c) {
          d[c] = z[3 * b + c] - xr[3 * b + c];
          if (origin >= 0) {
            d[c] -= z[3 * origin + c] - xr[3 * origin + c];
          }
        }
        worst[model][k][b] = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      }
    }
    myb_ias15_destroy(in);
    myb_nbody_destroy(s);
  }
  myb_spk_destroy(spk);
  release_figures(&fig);
  printf("position differences from DE440 (m; Sun barycentric, planets heliocentric, Moon "
         "geocentric)\n%-8s %11s %11s %11s  %11s %11s %11s  %11s %11s %11s\n",
         "", "EIH+fig 30d", "1 yr", "10 yr", "EIH 30 d", "1 yr", "10 yr", "Newton 30d", "1 yr",
         "10 yr");
  for (int b = 0; b < BODIES; ++b) {
    printf("%-8s %11.3f %11.3f %11.3f  %11.3f %11.3f %11.3f  %11.3f %11.3f %11.3f\n", names[b],
           worst[2][0][b], worst[2][1][b], worst[2][2][b], worst[1][0][b], worst[1][1][b],
           worst[1][2][b], worst[0][0][b], worst[0][1][b], worst[0][2][b]);
  }
  /* Measured (Windows and Gentoo alike): after 1 year EIH leaves Mercury 106 m, Venus 10 m,
   * Mars 68 m off DE440 against 58, 99 and 40 km with Newton. What remains: the Moon (23 km:
   * the figures of the Earth and the Moon and the tides, phase 7), the Earth through it (1/82
   * of the Moon's error), the barycentre (~270 m/yr: the momentum of DE440's 373 small bodies)
   * seen in the heliocentric outer planets, the small bodies themselves on Mars. */
  const double limit_1yr[BODIES] = {400, 200, 30, 400, 30000, 150, 400, 400, 400, 400, 400};
  for (int b = 0; b < BODIES; ++b) {
    MT_CHECK(worst[1][1][b] < limit_1yr[b]);
  }
  MT_CHECK(worst[1][0][1] < 10.0 && worst[1][2][1] < 2000.0); /* Mercury, 30 days and 10 years */
  MT_CHECK(worst[1][2][2] < 500.0);                          /* Venus, 10 years */
  for (int b = 1; b < BODIES; ++b) {
    if (b != 4 && b != 3) { /* GR improves every planet over Newton at 10 years */
      MT_CHECK(worst[1][2][b] <= worst[0][2][b]);
    }
  }
  MT_CHECK(worst[1][1][1] < 0.005 * worst[0][1][1]); /* GR accounts for Mercury */
  /* Phase 7, the figures (measured: Moon 0.65 m after 1 year, 239 m after 10 years — the
   * secular tidal acceleration of DE440, 25.8"/cy^2, gives 240 m in 10 years; Mercury 14 m
   * after 1 year with the solar J2, against 106 m without) */
  MT_CHECK(worst[2][1][MOON] < 5.0 && worst[2][2][MOON] < 400.0);
  MT_CHECK(worst[2][1][MOON] < 1e-3 * worst[1][1][MOON]);
  MT_CHECK(worst[2][1][1] < 30.0 && worst[2][2][1] < 300.0);
  MT_CHECK(worst[2][1][EARTH] < worst[1][1][EARTH] && worst[2][2][EARTH] < worst[1][2][EARTH]);
  return mt_finish("test_de440_rel");
}
