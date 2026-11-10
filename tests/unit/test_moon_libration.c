/* SPDX-License-Identifier: BSD-3-Clause */
/* The Moon's physical librations integrated with the library's rotational dynamics against
 * DE440 (moon_pa_de440_200625.bpc). The model is that of the JPL ephemerides (W. M. Folkner,
 * J. G. Williams, D. H. Boggs, R. S. Park, P. Kuchynka, IPN Progress Report 42-196, 2014, §IV),
 * every constant from header.440:
 *   - mantle with the undistorted moments from beta = LBET, gamma = LGAM and J2 = J2M,
 *     distorted by the Earth's tide (Love number K2M, time lag TAUM: the Earth's position and
 *     the mantle frame at t - tau) and by the spin (k2 R^5 / 3G (w w^T - (w^2 - n^2)/3 1 -
 *     n^2 z z^T)), with the rate of the tidal part;
 *   - fluid core of polar moment IFAC C_T, oblateness COBLAT, friction KVC
 *     (myb_core_mantle_rates);
 *   - torques of the Earth, Sun, Mercury, Venus, Mars, Jupiter on the Moon's figure: degree 2
 *     from the distorted inertia (myb_torque_point_mass), degrees 3..6 from the harmonics
 *     (myb_torque_field), and the Earth's J2 (pole of the CIP) on the Moon's quadrupole
 *     (myb_torque_figure_figure);
 *   - initial conditions of the DE440 integration (PHI, THT, PSI, OMEGA*, OMGC*) at
 *     JD 2440400.5 TDB; positions from de440.bsp.
 * The attitude is integrated as a quaternion with DOP853. */
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

#define DOWNLOAD MYB_DATA_DIR "/download/"

static const double DAY = 86400.0;
static const double ARCSEC = 4.848136811095359935899141e-6;

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

enum { BODIES = 6, EARTH = 0, DEGREE = 6 };
static const int32_t naif[BODIES] = {399, 10, 1, 2, 4, 5};

typedef struct lunar {
  myb_spk* spk;
  myb_cip* cip;
  myb_gravity_field* field; /* degrees 3..6 of the Moon */
  myb_gravity_work* work;
  double gm_moon;
  double gm[BODIES];
  double undistorted[3]; /* A, B, C of the whole Moon per unit mass (m^2) */
  double core[2];        /* A_c, C_c per unit mass */
  double kv;             /* m^2/s */
  double k2;
  double tau;            /* s */
  double rm;             /* m */
  double n;              /* rad/s */
  double re;
  double j2e;
  int64_t evaluations;
} lunar;

static void mat_vec(const double m[9], const double v[3], double out[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = m[3 * i] * v[0] + m[3 * i + 1] * v[1] + m[3 * i + 2] * v[2];
  }
}

static void cross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static myb_status rhs(void* context, double t, const double* y, double* dy) {
  lunar* L = context;
  ++L->evaluations;
  double q[4] = {y[0], y[1], y[2], y[3]};
  myb_quat_normalize(q);
  const double* wm = y + 4;
  const double* wc = y + 7;
  double m[9];
  myb_quat_to_matrix(q, m);

  /* delayed tide: the Earth at t - tau in the mantle frame of t - tau */
  const double rot[3] = {-wm[0] * L->tau, -wm[1] * L->tau, -wm[2] * L->tau};
  double qe[4];
  double qd[4];
  myb_quat_from_rotation_vector(rot, qe);
  myb_quat_multiply(q, qe, qd);
  double md[9];
  myb_quat_to_matrix(qd, md);
  double sd[6];
  myb_status st = myb_spk_state(L->spk, naif[EARTH], 301, myb_time_from_seconds(t - L->tau), sd,
                                NULL);
  if (st != MYB_OK) {
    return st;
  }
  double r[3];
  double rv[3];
  mat_vec(md, sd, r);
  mat_vec(md, sd + 3, rv);
  double wxr[3];
  cross(wm, r, wxr);
  for (int i = 0; i < 3; ++i) {
    rv[i] -= wxr[i]; /* rate of the components in the rotating frame */
  }
  const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  const double rr = sqrt(r2);
  const double r5 = r2 * r2 * rr;
  const double rdot = r[0] * rv[0] + r[1] * rv[1] + r[2] * rv[2];
  const double r5r = L->rm * L->rm * L->rm * L->rm * L->rm;
  const double ft = -L->k2 * (L->gm[EARTH] / L->gm_moon) * r5r;
  const double fs = L->k2 * r5r / (3.0 * L->gm_moon);
  const double w2 = wm[0] * wm[0] + wm[1] * wm[1] + wm[2] * wm[2];
  const double n2 = L->n * L->n;
  double im[9];
  double dim[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double d = i == j ? 1.0 : 0.0;
      const double tide = (r[i] * r[j] - d * r2 / 3.0) / r5;
      const double tide_rate = (rv[i] * r[j] + r[i] * rv[j]) / r5 -
                               5.0 * rdot * r[i] * r[j] / (r5 * r2) + d * rdot / r5;
      const double spin = wm[i] * wm[j] - d * (w2 - n2) / 3.0 - (i == 2 && j == 2 ? n2 : 0.0);
      const double undistorted_mantle =
          d * (L->undistorted[i] - (i == 2 ? L->core[1] : L->core[0]));
      im[3 * i + j] = undistorted_mantle + ft * tide + fs * spin;
      dim[3 * i + j] = ft * tide_rate;
    }
  }
  double it[9];
  memcpy(it, im, sizeof it);
  it[0] += L->core[0];
  it[4] += L->core[0];
  it[8] += L->core[1];

  /* torques per unit lunar mass */
  double torque[3] = {0.0, 0.0, 0.0};
  double re_m[3] = {0.0, 0.0, 0.0};
  for (int b = 0; b < BODIES; ++b) {
    double s[6];
    st = myb_spk_state(L->spk, naif[b], 301, myb_time_from_seconds(t), s, NULL);
    if (st != MYB_OK) {
      return st;
    }
    double rb[3];
    mat_vec(m, s, rb);
    if (b == EARTH) {
      memcpy(re_m, rb, sizeof re_m);
    }
    const double b2 = rb[0] * rb[0] + rb[1] * rb[1] + rb[2] * rb[2];
    const double f = 3.0 * L->gm[b] / (b2 * b2 * sqrt(b2));
    double irb[3];
    mat_vec(it, rb, irb);
    double c[3];
    cross(rb, irb, c);
    double a[3];
    st = myb_gravity_field_eval(L->field, L->work, rb, 3, DEGREE, NULL, a);
    if (st != MYB_OK) {
      return st;
    }
    double tf[3];
    myb_torque_field(L->gm[b] / L->gm_moon, rb, a, tf);
    for (int i = 0; i < 3; ++i) {
      torque[i] += f * c[i] + tf[i];
    }
  }
  /* the Earth's J2 on the Moon's quadrupole */
  double xys[3];
  myb_cip_xys(L->cip, myb_time_from_seconds(t), xys, NULL);
  const double pole_i[3] = {xys[0], xys[1], sqrt(1.0 - xys[0] * xys[0] - xys[1] * xys[1])};
  double p[3];
  mat_vec(m, pole_i, p);
  double gq_e[9];
  double gq_m[9];
  const double tr = it[0] + it[4] + it[8];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double d = i == j ? 1.0 : 0.0;
      gq_e[3 * i + j] = -L->gm[EARTH] * L->re * L->re * L->j2e * (3.0 * p[i] * p[j] - d);
      gq_m[3 * i + j] = L->gm_moon * (tr * d - 3.0 * it[3 * i + j]);
    }
  }
  double tff[3];
  st = myb_torque_figure_figure(L->gm_moon, gq_m, gq_e, re_m, tff);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    torque[i] += tff[i];
  }

  st = myb_core_mantle_rates(im, dim, L->core, L->kv, wm, wc, torque, dy + 4, dy + 7);
  if (st != MYB_OK) {
    return st;
  }
  myb_quat_derivative(y, wm, dy);
  return MYB_OK;
}

static double factorial_ratio(int n, int m) { /* (n - m)! / (n + m)! */
  double r = 1.0;
  for (int k = n - m + 1; k <= n + m; ++k) {
    r /= k;
  }
  return r;
}

static double angle_between(const double a[9], const double b[9]) {
  /* rotation angle of a b^T from its trace, via the antisymmetric part for small angles */
  double c[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      c[3 * i + j] = a[3 * i] * b[3 * j] + a[3 * i + 1] * b[3 * j + 1] + a[3 * i + 2] * b[3 * j + 2];
    }
  }
  const double s[3] = {0.5 * (c[5] - c[7]), 0.5 * (c[6] - c[2]), 0.5 * (c[1] - c[3])};
  return asin(fmin(1.0, sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2])));
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  static char header[32768];
  lunar L;
  memset(&L, 0, sizeof L);
  myb_pck* pck = NULL;
  MT_CHECK(myb_spk_create(NULL, &L.spk) == MYB_OK);
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  if (!read_file(DOWNLOAD "header.440", header, sizeof header) ||
      myb_spk_load_file(L.spk, DOWNLOAD "de440.bsp", NULL) != MYB_OK ||
      myb_pck_load_file(pck, DOWNLOAD "moon_pa_de440_200625.bpc", NULL) != MYB_OK) {
    printf("de440.bsp, header.440 or moon_pa_de440_200625.bpc not downloaded: skipped\n");
    myb_pck_destroy(pck);
    myb_spk_destroy(L.spk);
    return 77;
  }
  /* constants */
  double au = 0, gmb = 0, emrat = 0, gms = 0, gm1 = 0, gm2 = 0, gm4 = 0, gm5 = 0, am = 0, re = 0;
  double j2e = 0, k2m = 0, taum = 0, lbet = 0, lgam = 0, ifac = 0, coblat = 0, kvc = 0, j2m = 0;
  double phi = 0, tht = 0, psi = 0, ox = 0, oy = 0, oz = 0, cx = 0, cy = 0, cz = 0, jdepoc = 0;
  const char* names[] = {"AU",   "GMB",  "EMRAT",  "GMS",    "GM1",    "GM2",    "GM4",
                         "GM5",  "AM",   "RE",     "J2E",    "K2M",    "TAUM",   "LBET",
                         "LGAM", "IFAC", "COBLAT", "KVC",    "J2M",    "PHI",    "THT",
                         "PSI",  "OMEGAX", "OMEGAY", "OMEGAZ", "OMGCX", "OMGCY", "OMGCZ",
                         "JDEPOC"};
  double* values[] = {&au,   &gmb, &emrat, &gms, &gm1, &gm2, &gm4, &gm5, &am,  &re,
                      &j2e,  &k2m, &taum,  &lbet, &lgam, &ifac, &coblat, &kvc, &j2m, &phi,
                      &tht,  &psi, &ox,    &oy,   &oz,   &cx,  &cy,  &cz,  &jdepoc};
  int ok = 1;
  for (size_t k = 0; k < sizeof names / sizeof names[0]; ++k) {
    ok = ok && read_constant(header, names[k], values[k]);
  }
  MT_CHECK(ok);
  const double gm_unit = (au * 1e3) * (au * 1e3) * (au * 1e3) / (DAY * DAY);
  L.gm_moon = gmb / (1.0 + emrat) * gm_unit;
  L.gm[EARTH] = gmb * emrat / (1.0 + emrat) * gm_unit;
  L.gm[1] = gms * gm_unit;
  L.gm[2] = gm1 * gm_unit;
  L.gm[3] = gm2 * gm_unit;
  L.gm[4] = gm4 * gm_unit;
  L.gm[5] = gm5 * gm_unit;
  L.rm = am * 1e3;
  L.re = re * 1e3;
  L.j2e = j2e;
  L.k2 = k2m;
  L.tau = taum * DAY;
  L.n = 0.2299708 / DAY; /* the lunar mean motion (only the constant part of the spin term) */
  /* undistorted moments (eq. 36): A, B, C = 2 (1 - b g, 1 + g, 1 + b) / (2 b - g + b g) R^2 J2 */
  const double den = 2.0 * lbet - lgam + lbet * lgam;
  const double r2j2 = L.rm * L.rm * j2m;
  L.undistorted[0] = 2.0 * (1.0 - lbet * lgam) / den * r2j2;
  L.undistorted[1] = 2.0 * (1.0 + lgam) / den * r2j2;
  L.undistorted[2] = 2.0 * (1.0 + lbet) / den * r2j2;
  L.core[1] = ifac * L.undistorted[2];
  L.core[0] = (1.0 - coblat) * L.core[1];
  L.kv = kvc * L.undistorted[2] / DAY;
  printf("C/MR^2 = %.6f, beta %.6e, gamma %.6e, core C_c/C_T %.4g, f_c %.4g, k_v/C %.4e /day\n",
         L.undistorted[2] / (L.rm * L.rm), lbet, lgam, ifac, coblat, kvc);

  /* the lunar field, degrees 3..6, normalised from header.440 */
  double c[28] = {0.0};
  double s[28] = {0.0};
  c[0] = 1.0;
  for (int n = 3; n <= DEGREE && ok; ++n) {
    for (int mm = 0; mm <= n; ++mm) {
      char name[16];
      double cv = 0.0;
      double sv = 0.0;
      if (mm == 0) {
        snprintf(name, sizeof name, "J%dM", n);
        ok = ok && read_constant(header, name, &cv);
        cv = -cv;
      } else {
        snprintf(name, sizeof name, "C%d%dM", n, mm);
        ok = ok && read_constant(header, name, &cv);
        snprintf(name, sizeof name, "S%d%dM", n, mm);
        ok = ok && read_constant(header, name, &sv);
      }
      const double norm = sqrt((mm == 0 ? 1.0 : 2.0) * (2.0 * n + 1.0) * factorial_ratio(n, mm));
      c[n * (n + 1) / 2 + mm] = cv / norm;
      s[n * (n + 1) / 2 + mm] = sv / norm;
    }
  }
  MT_CHECK(ok);
  MT_CHECK(myb_gravity_field_create(L.gm_moon, L.rm, DEGREE, c, s, NULL, &L.field) == MYB_OK);
  MT_CHECK(myb_gravity_work_create(DEGREE, NULL, &L.work) == MYB_OK);
  myb_kpool* pool = NULL;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  MT_CHECK(myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau2006-cip.tk", NULL) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &L.cip, NULL) == MYB_OK);
  myb_kpool_destroy(pool);

  /* initial state of the DE440 integration */
  const double t0 = (jdepoc - 2451545.0) * DAY;
  const double ang[3] = {phi, tht, psi};
  const int seq[3] = {3, 1, 3};
  double m0[9];
  MT_CHECK(myb_euler_to_matrix(seq, ang, m0) == MYB_OK);
  myb_rotation pa;
  MT_CHECK(myb_pck_rotation(pck, 31008, myb_time_from_seconds(t0), &pa, NULL) == MYB_OK);
  const double start_diff = angle_between(m0, pa.m);
  printf("initial attitude: header vs PCK %.3e arcsec\n", start_diff / ARCSEC);
  MT_CHECK(start_diff < 1e-3 * ARCSEC);
  double y[10];
  myb_quat_from_matrix(m0, y);
  y[4] = ox / DAY;
  y[5] = oy / DAY;
  y[6] = oz / DAY;
  y[7] = cx / DAY;
  y[8] = cy / DAY;
  y[9] = cz / DAY;

  const myb_dop853_options o = {.rtol = 1e-13, .atol = 1e-18, .h_max = 0.25 * DAY};
  myb_dop853* d = NULL;
  MT_CHECK(myb_dop853_create(10, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_dop853_start(d, rhs, &L, t0, y) == MYB_OK);
  const int years = 2;
  double worst = 0.0;
  double worst_year1 = 0.0;
  double worst_rate = 0.0;
  for (int k = 1; k <= years * 73; ++k) { /* every 5 days */
    const double t = t0 + 5.0 * DAY * k;
    const myb_status st = myb_dop853_integrate(d, t);
    MT_CHECK(st == MYB_OK);
    if (st != MYB_OK) {
      break;
    }
    const double* ys = myb_dop853_state(d);
    double q[4] = {ys[0], ys[1], ys[2], ys[3]};
    double m[9];
    myb_quat_to_matrix(q, m);
    MT_CHECK(myb_pck_rotation(pck, 31008, myb_time_from_seconds(t), &pa, NULL) == MYB_OK);
    const double e = angle_between(m, pa.m);
    worst = fmax(worst, e);
    if (k <= 73) {
      worst_year1 = fmax(worst_year1, e);
    }
    /* the angular velocity of the PCK: dM M^T = -[w x] */
    double w[9];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        w[3 * i + j] = pa.dm[3 * i] * pa.m[3 * j] + pa.dm[3 * i + 1] * pa.m[3 * j + 1] +
                       pa.dm[3 * i + 2] * pa.m[3 * j + 2];
      }
    }
    const double wpck[3] = {w[5], -w[2], w[1]};
    for (int i = 0; i < 3; ++i) {
      worst_rate = fmax(worst_rate, fabs(ys[4 + i] - wpck[i]));
    }
    if (k % 73 == 0) {
      printf("  %d year(s): attitude vs DE440 %.4f arcsec (max so far %.4f), angular velocity "
             "%.2e rad/day\n",
             k / 73, e / ARCSEC, worst / ARCSEC, worst_rate * DAY);
    }
  }
  printf("lunar librations over %d years: max %.4f arcsec (%.2f m at the surface), first year "
         "%.4f arcsec; %lld evaluations\n",
         years, worst / ARCSEC, worst * L.rm, worst_year1 / ARCSEC, (long long)L.evaluations);
  /* measured: 0.040 arcsec in the first year, 0.071 over two (0.6 m at the surface), the
   * angular velocity to 6.6e-8 rad/day (3e-7). Leaving out a part of the model: the tidal
   * distortion 0.64 arcsec, its rate 0.97, the core friction 0.11, the figure-figure torque
   * 0.085, the spin distortion 0.077; the tide in the mantle frame of t instead of t - tau 29. */
  MT_CHECK(worst_year1 < 0.06 * ARCSEC);
  MT_CHECK(worst < 0.1 * ARCSEC);
  MT_CHECK(worst_rate * DAY < 1e-7);
  myb_dop853_destroy(d);
  myb_cip_destroy(L.cip);
  myb_gravity_work_destroy(L.work);
  myb_gravity_field_destroy(L.field);
  myb_pck_destroy(pck);
  myb_spk_destroy(L.spk);
  return mt_finish("test_moon_libration");
}
