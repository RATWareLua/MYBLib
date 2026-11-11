/* SPDX-License-Identifier: BSD-3-Clause */
/* Phase 18 criterion: DE440 integrated again by the library's own system of bodies
 * (myb/system.h) — the motion of the Sun, the planets and the Moon together with the rotation
 * of the Moon, nothing taken from the ephemeris but the state at the start:
 *   - point masses with the 1PN terms (EIH);
 *   - the figures: J2 of the Sun on every body; the zonals J2..J5 of the Earth and the field of
 *     the Moon to degree 6 on the Earth or the Moon, the Sun, Mercury, Venus, Mars and
 *     Jupiter; the quadrupoles of the Earth and the Moon on each other;
 *   - the tides raised on the Earth by the Moon and the Sun acting on the Moon, by band, with
 *     the Love numbers and the five time lags of DE440;
 *   - the rotation of the Moon: an anelastic mantle distorted by the tide of the Earth and by
 *     its own spin, with a fluid core; its degree-2 field follows the inertia tensor.
 * Constants from header.440; initial conditions of the DE440 integration at JD 2440400.5 TDB
 * (positions from de440.bsp, the lunar angles and angular velocities from the header). The
 * pole of the Earth is the CIP (IAU 2006/2000A); DE440 itself uses a simpler model of it.
 * Not in the model: the 343 asteroids and 30 trans-Neptunian objects of DE440 (they move the
 * barycentre and Mars), the Lense-Thirring effect of the Sun, the loss of mass of the Sun.
 * Compared: the positions with de440.bsp, the attitude of the Moon with
 * moon_pa_de440_200625.bpc. Needs data/download (fetch-data); skipped otherwise. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/attitude.h"
#include "myb/frames.h"
#include "myb/grav.h"
#include "myb/integ.h"
#include "myb/io.h"
#include "myb/kin.h"
#include "myb/system.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

static const double DAY = 86400.0;
static const double ARCSEC = 4.848136811095359935899141e-6;

enum { BODIES = 11, SUN = 0, EARTH = 3, MOON = 4, MOON_DEGREE = 6 };

static const char* names[BODIES] = {"Sun",    "Mercury", "Venus",  "Earth",   "Moon",  "Mars",
                                    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
static const int32_t naif[BODIES] = {10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9};

/* The constants of header.440 (GROUP 1040 names, GROUP 1041 values). */
static int read_constant(const char* text, const char* name, double* value) {
  const char* g1040 = strstr(text, "GROUP   1040");
  const char* g1041 = strstr(text, "GROUP   1041");
  if (g1040 == NULL || g1041 == NULL) {
    return 0;
  }
  const char* p = g1040 + 12;
  double count = 0.0;
  size_t used = 0;
  while (*p == ' ' || *p == '\n' || *p == '\r') { /* bound: white space */
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  int index = -1;
  for (int k = 0; k < (int)count && p < g1041; ++k) {
    while (*p == ' ' || *p == '\n' || *p == '\r') { /* bound: white space */
      ++p;
    }
    size_t n = 0;
    while (p[n] != ' ' && p[n] != '\n' && p[n] != '\r' && p[n] != '\0') { /* bound: a name */
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
  while (*p == ' ' || *p == '\n' || *p == '\r') { /* bound: white space */
    ++p;
  }
  if (myb_parse_double(p, strlen(p), &count, &used) != MYB_OK) {
    return 0;
  }
  p += used;
  for (int k = 0; k <= index; ++k) {
    while (*p == ' ' || *p == '\n' || *p == '\r') { /* bound: white space */
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

static double factorial_ratio(int n, int m) { /* (n - m)! / (n + m)! */
  double r = 1.0;
  for (int k = n - m + 1; k <= n + m; ++k) {
    r /= k;
  }
  return r;
}

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

/* The frame of the Earth: its pole is the CIP (TT taken as TDB, 2 ms apart). */
static myb_status earth_frame(void* context, myb_time t, myb_orientation* o) {
  double xys[3];
  myb_cip_xys(context, t, xys, NULL);
  const double pole[3] = {xys[0], xys[1], sqrt(1.0 - xys[0] * xys[0] - xys[1] * xys[1])};
  memset(o, 0, sizeof *o);
  pole_frame(pole, o->rotation.m);
  return MYB_OK;
}

/* The frame of the Sun: its IAU pole, alpha 286.13 deg, delta 63.87 deg. */
static myb_status sun_frame(void* context, myb_time t, myb_orientation* o) {
  (void)context;
  (void)t;
  const double ra = 286.13 * 3.141592653589793 / 180.0;
  const double de = 63.87 * 3.141592653589793 / 180.0;
  const double pole[3] = {cos(de) * cos(ra), cos(de) * sin(ra), sin(de)};
  memset(o, 0, sizeof *o);
  pole_frame(pole, o->rotation.m);
  return MYB_OK;
}

static void state_at(const myb_spk* spk, double t, double* x, double* v) {
  for (int b = 0; b < BODIES; ++b) {
    double s[6];
    MT_CHECK(myb_spk_state(spk, naif[b], 0, myb_time_from_seconds(t), s, NULL) == MYB_OK);
    memcpy(x + 3 * b, s, 3 * sizeof(double));
    memcpy(v + 3 * b, s + 3, 3 * sizeof(double));
  }
}

static double angle_between(const double a[9], const double b[9]) {
  double c[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      c[3 * i + j] = a[3 * i] * b[3 * j] + a[3 * i + 1] * b[3 * j + 1] + a[3 * i + 2] * b[3 * j + 2];
    }
  }
  const double s[3] = {0.5 * (c[5] - c[7]), 0.5 * (c[6] - c[2]), 0.5 * (c[1] - c[3])};
  return asin(fmin(1.0, sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2])));
}

typedef struct constants {
  double gm[BODIES];
  double clight;
  double re, am, asun, j2sun;
  double jn_earth[4];
  double k2e[3], taue[3], tauer[2];
  double k2m, taum, lbet, lgam, ifac, coblat, kvc, j2m;
  double phi, tht, psi, omega[3], omega_core[3];
  double jdepoc;
} constants;

static int read_constants(const char* header, constants* k) {
  double au = 0, emrat = 0, gmb = 0;
  const char* gm_names[10] = {"GMS", "GM1", "GM2", "GMB", "GM4", "GM5", "GM6", "GM7", "GM8", "GM9"};
  double g[10];
  int ok = read_constant(header, "AU", &au) && read_constant(header, "EMRAT", &emrat) &&
           read_constant(header, "CLIGHT", &k->clight);
  for (int i = 0; i < 10; ++i) {
    ok = ok && read_constant(header, gm_names[i], &g[i]);
  }
  const struct {
    const char* name;
    double* value;
  } list[] = {{"RE", &k->re},          {"AM", &k->am},          {"ASUN", &k->asun},
              {"J2SUN", &k->j2sun},    {"J2E", &k->jn_earth[0]}, {"J3E", &k->jn_earth[1]},
              {"J4E", &k->jn_earth[2]}, {"J5E", &k->jn_earth[3]}, {"K2E0", &k->k2e[0]},
              {"K2E1", &k->k2e[1]},    {"K2E2", &k->k2e[2]},    {"TAUE0", &k->taue[0]},
              {"TAUE1", &k->taue[1]},  {"TAUE2", &k->taue[2]},  {"TAUER1", &k->tauer[0]},
              {"TAUER2", &k->tauer[1]}, {"K2M", &k->k2m},        {"TAUM", &k->taum},
              {"LBET", &k->lbet},      {"LGAM", &k->lgam},      {"IFAC", &k->ifac},
              {"COBLAT", &k->coblat},  {"KVC", &k->kvc},        {"J2M", &k->j2m},
              {"PHI", &k->phi},        {"THT", &k->tht},        {"PSI", &k->psi},
              {"OMEGAX", &k->omega[0]}, {"OMEGAY", &k->omega[1]}, {"OMEGAZ", &k->omega[2]},
              {"OMGCX", &k->omega_core[0]}, {"OMGCY", &k->omega_core[1]},
              {"OMGCZ", &k->omega_core[2]}, {"JDEPOC", &k->jdepoc}};
  for (size_t i = 0; i < sizeof list / sizeof list[0]; ++i) {
    ok = ok && read_constant(header, list[i].name, list[i].value);
  }
  if (!ok) {
    return 0;
  }
  gmb = g[3];
  const double unit = (au * 1e3) * (au * 1e3) * (au * 1e3) / (DAY * DAY);
  k->gm[0] = g[0] * unit;
  k->gm[1] = g[1] * unit;
  k->gm[2] = g[2] * unit;
  k->gm[EARTH] = gmb * emrat / (1.0 + emrat) * unit;
  k->gm[MOON] = gmb / (1.0 + emrat) * unit;
  for (int b = 5; b < BODIES; ++b) {
    k->gm[b] = g[b - 1] * unit;
  }
  return 1;
}

typedef struct model {
  myb_system* system;
  myb_gravity_field* earth;
  myb_gravity_field* moon;
  myb_gravity_field* sun;
  myb_cip* cip;
} model;

static void release(model* m) {
  myb_system_destroy(m->system);
  myb_gravity_field_destroy(m->earth);
  myb_gravity_field_destroy(m->moon);
  myb_gravity_field_destroy(m->sun);
  myb_cip_destroy(m->cip);
}

/* tides: the tides of the Earth on the Moon; rotation: the Moon's rotation integrated. */
static int build(model* m, const constants* k, const char* header, int tides) {
  memset(m, 0, sizeof *m);
  myb_kpool* pool = NULL;
  int ok = myb_kpool_create(NULL, &pool) == MYB_OK &&
           myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau2006-cip.tk", NULL) == MYB_OK &&
           myb_cip_create(pool, "IAU2006A", NULL, &m->cip, NULL) == MYB_OK;
  myb_kpool_destroy(pool);
  /* the fields, normalised from the unnormalised constants (J_n = -C_n0) */
  double ce[21] = {0.0};
  double se[21] = {0.0};
  ce[0] = 1.0;
  for (int n = 2; n <= 5; ++n) {
    ce[n * (n + 1) / 2] = -k->jn_earth[n - 2] / sqrt(2.0 * n + 1.0);
  }
  double cm[28] = {0.0};
  double sm[28] = {0.0};
  cm[0] = 1.0;
  for (int n = 3; n <= MOON_DEGREE && ok; ++n) {
    for (int mm = 0; mm <= n; ++mm) {
      char name[16];
      double c = 0.0;
      double s = 0.0;
      if (mm == 0) {
        snprintf(name, sizeof name, "J%dM", n);
        ok = ok && read_constant(header, name, &c);
        c = -c;
      } else {
        snprintf(name, sizeof name, "C%d%dM", n, mm);
        ok = ok && read_constant(header, name, &c);
        snprintf(name, sizeof name, "S%d%dM", n, mm);
        ok = ok && read_constant(header, name, &s);
      }
      const double norm = sqrt((mm == 0 ? 1.0 : 2.0) * (2.0 * n + 1.0) * factorial_ratio(n, mm));
      cm[n * (n + 1) / 2 + mm] = c / norm;
      sm[n * (n + 1) / 2 + mm] = s / norm;
    }
  }
  const double cs[6] = {1.0, 0.0, 0.0, -k->j2sun / sqrt(5.0), 0.0, 0.0};
  const double ss[6] = {0.0};
  ok = ok &&
       myb_gravity_field_create(k->gm[EARTH], k->re * 1e3, 5, ce, se, NULL, &m->earth) == MYB_OK &&
       myb_gravity_field_create(k->gm[MOON], k->am * 1e3, MOON_DEGREE, cm, sm, NULL, &m->moon) ==
           MYB_OK &&
       myb_gravity_field_create(k->gm[SUN], k->asun * 1e3, 2, cs, ss, NULL, &m->sun) == MYB_OK &&
       myb_system_create(BODIES, k->gm, NULL, &m->system) == MYB_OK;
  if (!ok) {
    return 0;
  }
  myb_system* s = m->system;
  const myb_ppn ppn = {k->clight * 1e3, 1.0, 1.0};
  ok = myb_nbody_set_relativity(myb_system_nbody(s), &ppn) == MYB_OK;
  /* the Sun */
  ok = ok && myb_system_set_orientation(s, SUN, sun_frame, NULL) == MYB_OK &&
       myb_system_set_field(s, SUN, m->sun, 0) == MYB_OK;
  for (size_t b = 1; b < BODIES && ok; ++b) {
    ok = myb_system_add_field(s, SUN, b) == MYB_OK;
  }
  /* the Earth */
  static const size_t near[5] = {0, 1, 2, 5, 6}; /* Sun, Mercury, Venus, Mars, Jupiter */
  ok = ok && myb_system_set_orientation(s, EARTH, earth_frame, m->cip) == MYB_OK &&
       myb_system_set_field(s, EARTH, m->earth, 0) == MYB_OK &&
       myb_system_add_field(s, EARTH, MOON) == MYB_OK;
  for (int i = 0; i < 5 && ok; ++i) {
    ok = myb_system_add_field(s, EARTH, near[i]) == MYB_OK;
  }
  if (tides) {
    const myb_system_tide tide = {
        .radius = k->re * 1e3,
        .k2 = {k->k2e[0], k->k2e[1], k->k2e[2]},
        .orbital_lag = {k->taue[0] * DAY, k->taue[1] * DAY, k->taue[2] * DAY},
        .rotation_lag = {0.0, k->tauer[0] * DAY, k->tauer[1] * DAY},
        .rotation_rate = 7.292115e-5};
    ok = ok && myb_system_set_tide(s, EARTH, &tide) == MYB_OK &&
         myb_system_add_tide(s, EARTH, MOON, MOON) == MYB_OK &&
         myb_system_add_tide(s, EARTH, SUN, MOON) == MYB_OK;
  }
  /* the Moon: undistorted moments A, B, C = 2 (1 - b g, 1 + g, 1 + b) / (2 b - g + b g) R^2 J2
   * (Folkner et al. 2014, eq. 36) */
  const double rm = k->am * 1e3;
  const double den = 2.0 * k->lbet - k->lgam + k->lbet * k->lgam;
  const double r2j2 = rm * rm * k->j2m;
  myb_system_rotation moon = {
      .moments = {2.0 * (1.0 - k->lbet * k->lgam) / den * r2j2, 2.0 * (1.0 + k->lgam) / den * r2j2,
                  2.0 * (1.0 + k->lbet) / den * r2j2},
      .radius = rm,
      .k2 = k->k2m,
      .lag = k->taum * DAY,
      .mean_motion = 0.2299708 / DAY};
  moon.core[1] = k->ifac * moon.moments[2];
  moon.core[0] = (1.0 - k->coblat) * moon.core[1];
  moon.core_friction = k->kvc * moon.moments[2] / DAY;
  ok = ok && myb_system_set_rotation(s, MOON, &moon) == MYB_OK &&
       myb_system_set_field(s, MOON, m->moon, 0) == MYB_OK &&
       myb_system_add_field(s, MOON, EARTH) == MYB_OK &&
       myb_system_add_distortion(s, MOON, EARTH) == MYB_OK &&
       myb_system_add_figures(s, EARTH, MOON) == MYB_OK;
  for (int i = 0; i < 5 && ok; ++i) {
    ok = myb_system_add_field(s, MOON, near[i]) == MYB_OK;
  }
  return ok;
}

typedef struct outcome {
  double position[3][BODIES]; /* m, at 30 days, 1 year and the last epoch */
  double attitude[3];         /* rad: the worst up to each epoch */
  double rate;                /* rad/s: the worst angular velocity difference */
  long long evaluations;
  long long steps;
} outcome;

static int run(const constants* k, const char* header, const myb_spk* spk, const myb_pck* pck,
               int tides, double years, outcome* out) {
  memset(out, 0, sizeof *out);
  model m;
  if (!build(&m, k, header, tides)) {
    release(&m);
    return 0;
  }
  const size_t n = myb_system_coordinates(m.system);
  MT_CHECK(n == 3 * BODIES + 9);
  static double y[2 * (3 * BODIES + 9)];
  memset(y, 0, sizeof y);
  const double t0 = (k->jdepoc - 2451545.0) * DAY;
  state_at(spk, t0, y, y + n);
  const double angles[3] = {k->phi, k->tht, k->psi};
  const int sequence[3] = {3, 1, 3};
  double m0[9];
  double q0[4];
  MT_CHECK(myb_euler_to_matrix(sequence, angles, m0) == MYB_OK);
  myb_quat_from_matrix(m0, q0);
  const double w[3] = {k->omega[0] / DAY, k->omega[1] / DAY, k->omega[2] / DAY};
  const double wc[3] = {k->omega_core[0] / DAY, k->omega_core[1] / DAY, k->omega_core[2] / DAY};
  MT_CHECK(myb_system_pack_rotation(m.system, MOON, q0, w, wc, y) == MYB_OK);
  const myb_ias15_options o = {.epsilon = 1e-10};
  myb_ias15* in = NULL;
  MT_CHECK(myb_system_integrator(m.system, &o, NULL, &in) == MYB_OK);
  MT_CHECK(myb_ias15_start(in, myb_system_acceleration, m.system, t0, y) == MYB_OK);
  const double epochs[3] = {30.0, 365.25, 365.25 * years};
  double worst_attitude = 0.0;
  int epoch = 0;
  const int samples = (int)(365.25 * years / 5.0) + 1; /* every 5 days, and the epochs */
  int ok = 1;
  for (int step = 1; step <= samples && epoch < 3 && ok; ++step) {
    double t = t0 + 5.0 * DAY * step;
    const double t_epoch = t0 + epochs[epoch] * DAY;
    const int at_epoch = t >= t_epoch;
    if (at_epoch) {
      t = t_epoch;
      --step; /* the sample is taken after the epoch */
    }
    ok = myb_ias15_integrate(in, t) == MYB_OK;
    MT_CHECK(ok);
    const double* z = myb_ias15_state(in);
    double q[4];
    double wm[3];
    MT_CHECK(myb_system_unpack_rotation(m.system, MOON, z, q, wm, NULL) == MYB_OK);
    double mq[9];
    myb_quat_to_matrix(q, mq);
    myb_rotation pa;
    MT_CHECK(myb_pck_rotation(pck, 31008, myb_time_from_seconds(t), &pa, NULL) == MYB_OK);
    worst_attitude = fmax(worst_attitude, angle_between(mq, pa.m));
    /* the angular velocity of the PCK: dM M^T = -[w x] */
    double wx[9];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        wx[3 * i + j] = pa.dm[3 * i] * pa.m[3 * j] + pa.dm[3 * i + 1] * pa.m[3 * j + 1] +
                        pa.dm[3 * i + 2] * pa.m[3 * j + 2];
      }
    }
    const double wpck[3] = {wx[5], -wx[2], wx[1]};
    for (int i = 0; i < 3; ++i) {
      out->rate = fmax(out->rate, fabs(wm[i] - wpck[i]));
    }
    if (at_epoch) {
      static double xr[3 * BODIES];
      static double vr[3 * BODIES];
      state_at(spk, t, xr, vr);
      /* barycentric for the Sun, heliocentric for the planets and the Earth, geocentric for
       * the Moon: the barycentre itself drifts by the momentum of the missing small bodies */
      for (int b = 0; b < BODIES; ++b) {
        const int origin = b == 0 ? -1 : (b == MOON ? EARTH : 0);
        double d[3];
        for (int c = 0; c < 3; ++c) {
          d[c] = z[3 * b + c] - xr[3 * b + c];
          if (origin >= 0) {
            d[c] -= z[3 * origin + c] - xr[3 * origin + c];
          }
        }
        out->position[epoch][b] = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      }
      out->attitude[epoch] = worst_attitude;
      ++epoch;
    }
  }
  out->evaluations = (long long)myb_ias15_statistics(in).evaluations;
  out->steps = (long long)myb_ias15_statistics(in).accepted;
  myb_ias15_destroy(in);
  release(&m);
  return ok && epoch == 3;
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static char header[32768];
  myb_spk* spk = NULL;
  myb_pck* pck = NULL;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  if (!read_file(DOWNLOAD "header.440", header, sizeof header) ||
      myb_spk_load_file(spk, DOWNLOAD "de440.bsp", NULL) != MYB_OK ||
      myb_pck_load_file(pck, DOWNLOAD "moon_pa_de440_200625.bpc", NULL) != MYB_OK) {
    printf("de440.bsp, header.440 or moon_pa_de440_200625.bpc not downloaded: skipped\n");
    myb_pck_destroy(pck);
    myb_spk_destroy(spk);
    return MTEST_SKIP;
  }
  static constants k;
  MT_CHECK(read_constants(header, &k));
  const double years = 5.0;
  static outcome full;
  static outcome no_tides;
  MT_CHECK(run(&k, header, spk, pck, 1, years, &full));
  MT_CHECK(run(&k, header, spk, pck, 0, years, &no_tides));
  myb_pck_destroy(pck);
  myb_spk_destroy(spk);
  printf("DE440 integrated again from JD %.1f by the system of the library (%lld steps, %lld "
         "evaluations over %.0f years)\n",
         k.jdepoc, full.steps, full.evaluations, years);
  printf("position differences from DE440 (m; Sun barycentric, planets heliocentric, Moon "
         "geocentric)\n%-8s %12s %12s %12s   %s\n",
         "", "30 days", "1 year", "5 years", "without the tides of the Earth, 5 years");
  for (int b = 0; b < BODIES; ++b) {
    printf("%-8s %12.4f %12.4f %12.4f   %12.4f\n", names[b], full.position[0][b],
           full.position[1][b], full.position[2][b], no_tides.position[2][b]);
  }
  printf("attitude of the Moon against the PCK of DE440 (arcsec, the worst up to the epoch): "
         "%.4f, %.4f, %.4f (%.2f m at the surface); angular velocity %.2e rad/day\n",
         full.attitude[0] / ARCSEC, full.attitude[1] / ARCSEC, full.attitude[2] / ARCSEC,
         full.attitude[2] * k.am * 1e3, full.rate * DAY);
  /* the Moon: with the tides of the Earth it follows DE440; without them it runs ahead by
   * the secular acceleration, 25.8"/cy^2 (58 m in 5 years) */
  MT_CHECK(full.position[0][MOON] < 0.05);
  MT_CHECK(full.position[1][MOON] < 1.0);
  MT_CHECK(full.position[2][MOON] < 10.0);
  MT_CHECK(no_tides.position[2][MOON] > 30.0 && no_tides.position[2][MOON] < 120.0);
  MT_CHECK(full.position[2][MOON] < 0.1 * no_tides.position[2][MOON]);
  /* the librations, integrated with the motion */
  MT_CHECK(full.attitude[1] < 0.06 * ARCSEC);
  MT_CHECK(full.attitude[2] < 0.5 * ARCSEC);
  MT_CHECK(full.rate * DAY < 1e-6);
  /* the planets: what the missing small bodies allow */
  MT_CHECK(full.position[1][1] < 30.0); /* Mercury, a year */
  MT_CHECK(full.position[1][2] < 30.0); /* Venus */
  MT_CHECK(full.position[1][EARTH] < 30.0);
  return mt_finish("test_de440_system");
}
