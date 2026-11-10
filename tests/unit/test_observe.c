/* SPDX-License-Identifier: BSD-3-Clause */
/* Observables: light time (reception and transmission) against CSPICE N0067 on DE440 (spkezr
 * CN, XCN and their stellar aberration), aberration and deflection against ERFA (eraAb, eraLd),
 * the relativistic delay against the quadrature of the 1PN null line, the range rate against
 * differences of light-time solutions, round trips. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "erfa.h"
#include "mtest.h"
#include "myb/kin.h"
#include "myb/observe.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

static const double C = 299792458.0;
static const double GM_SUN = 1.32712440041e20;

typedef struct body_ref {
  const myb_spk* spk;
  int32_t id;
} body_ref;

static myb_status spk_state(void* context, myb_time t, double state[6]) {
  const body_ref* b = context;
  return myb_spk_state(b->spk, b->id, 0, t, state, NULL);
}

static double norm(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double angle(const double a[3], const double b[3]) {
  const double c[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                       a[0] * b[1] - a[1] * b[0]};
  return atan2(norm(c), a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static void test_spice(void) {
  myb_spk* spk = NULL;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  if (myb_spk_load_file(spk, DOWNLOAD "de440.bsp", NULL) != MYB_OK) {
    printf("de440.bsp not downloaded: light time against CSPICE skipped\n");
    myb_spk_destroy(spk);
    return;
  }
  FILE* f = fopen(MYB_TEST_DATA_DIR "/observe/spice_lighttime.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    myb_spk_destroy(spk);
    return;
  }
  const myb_light_model model = {.c = C, .gamma = 1.0};
  body_ref earth = {spk, 399};
  static char line[1024];
  double worst_pos = 0.0;
  double worst_vel = 0.0;
  double worst_lt = 0.0;
  double worst_aberr = 0.0;
  int rows = 0;
  while (fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    char corr[8];
    int target = 0;
    int used = 0;
    if (sscanf(line, "%7s %d%n", corr, &target, &used) != 2) {
      continue;
    }
    double v[8];
    const char* p = line + used;
    int n = 0;
    while (n < 8) {
      char* end = NULL;
      v[n] = strtod(p, &end);
      if (end == p) {
        break;
      }
      p = end;
      ++n;
    }
    if (n != 8) {
      continue;
    }
    const myb_time t = myb_time_from_seconds(v[0]);
    double obs[6];
    MT_CHECK(spk_state(&earth, t, obs) == MYB_OK);
    body_ref tb = {spk, target};
    const int transmit = corr[0] == 'X';
    myb_light_path path;
    MT_CHECK(myb_light_time(&model, spk_state, &tb, t, obs,
                            transmit ? MYB_LIGHT_TRANSMISSION : MYB_LIGHT_RECEPTION,
                            &path) == MYB_OK);
    double rel[3];
    double relv[3];
    const double sign = transmit ? 1.0 : -1.0;
    for (int i = 0; i < 3; ++i) {
      rel[i] = path.target[i] - obs[i];
      /* d/dt of x_T(t + sign lt) - x_o(t) */
      relv[i] = path.target[3 + i] * (1.0 + sign * path.range_rate / C) - obs[3 + i];
    }
    const double spice_pos[3] = {v[2] * 1e3, v[3] * 1e3, v[4] * 1e3};
    const double spice_vel[3] = {v[5] * 1e3, v[6] * 1e3, v[7] * 1e3};
    if (strchr(corr, 'S') == NULL) {
      for (int i = 0; i < 3; ++i) {
        worst_pos = fmax(worst_pos, fabs(rel[i] - spice_pos[i]));
        worst_vel = fmax(worst_vel, fabs(relv[i] - spice_vel[i]));
      }
      worst_lt = fmax(worst_lt, fabs(path.light_time - v[1]));
    } else {
      /* stellar aberration: the observer's velocity (reversed for transmission) */
      double u[3] = {rel[0] / path.range, rel[1] / path.range, rel[2] / path.range};
      const double vo[3] = {-sign * obs[3], -sign * obs[4], -sign * obs[5]};
      double ab[3];
      myb_aberration(u, vo, C, ab);
      worst_aberr = fmax(worst_aberr, angle(ab, spice_pos));
    }
    ++rows;
  }
  fclose(f);
  printf("light time vs CSPICE (%d rows, Moon, Sun, Mars, Jupiter from the Earth, 1974-2026): "
         "position %.2e m, velocity %.2e m/s, light time %.2e s; aberration (Lorentz) vs "
         "SPICE's stellar aberration %.2e rad\n",
         rows, worst_pos, worst_vel, worst_lt, worst_aberr);
  MT_CHECK(rows == 208);
  MT_CHECK(worst_pos < 5e-3); /* 1.7 mm at Jupiter, 1e12 m: 2e-15 */
  MT_CHECK(worst_vel < 1e-9);
  MT_CHECK(worst_lt < 1e-11);
  /* SPICE corrects to first order in v/c: the difference is (v/c)^2 ~ 1e-8 */
  MT_CHECK(worst_aberr < 2e-8);

  /* round trip Earth -> Moon -> Earth (geocentres), with the Sun's delay */
  body_ref moon = {spk, 301};
  body_ref sun = {spk, 10};
  const myb_light_body sun_body = {.gm = GM_SUN, .state = spk_state, .context = &sun};
  const myb_light_model relativistic = {.c = C, .gamma = 1.0, .bodies = &sun_body, .count = 1};
  const myb_time t = myb_time_from_seconds(7.0e8);
  double total = 0.0;
  myb_time tb;
  myb_time tt;
  MT_CHECK(myb_round_trip(&relativistic, spk_state, &earth, spk_state, &moon, t, &total, &tb,
                          &tt) == MYB_OK);
  double obs[6];
  MT_CHECK(spk_state(&earth, t, obs) == MYB_OK);
  myb_light_path down;
  MT_CHECK(myb_light_time(&relativistic, spk_state, &moon, t, obs, MYB_LIGHT_RECEPTION, &down) ==
           MYB_OK);
  /* the legs differ by the Earth's barycentric motion along the line of sight (30 km/s x 2.7 s),
   * their sum is twice the geocentric distance at the bounce up to (v/c)^2 and the delays */
  double eb[6];
  double mb[6];
  MT_CHECK(spk_state(&earth, tb, eb) == MYB_OK);
  MT_CHECK(spk_state(&moon, tb, mb) == MYB_OK);
  const double geo[3] = {mb[0] - eb[0], mb[1] - eb[1], mb[2] - eb[2]};
  const double geocentric = 2.0 * norm(geo) / C;
  printf("round trip Earth-Moon: %.12f s (down-link %.12f s, Sun's delay %.3e s); "
         "2 |x_M - x_E| / c at the bounce %.12f s\n",
         total, down.light_time, down.delay, geocentric);
  MT_CHECK(fabs(total - geocentric - 2.0 * down.delay) < 1e-7);
  MT_CHECK(down.delay > 1e-8 && down.delay < 1e-7);
  myb_spk_destroy(spk);
}

/* ---- ERFA ---------------------------------------------------------------------------------- */

static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

static void test_erfa(void) {
  mt_rng rng = {31};
  double worst_ab = 0.0;
  double worst_ld = 0.0;
  const double au = 149597870700.0;
  for (int k = 0; k < 2000; ++k) {
    double p[3];
    double v[3];
    for (int i = 0; i < 3; ++i) {
      p[i] = 2.0 * uniform(&rng) - 1.0;
      v[i] = (2.0 * uniform(&rng) - 1.0) * 5e4;
    }
    const double np = norm(p);
    for (int i = 0; i < 3; ++i) {
      p[i] /= np;
    }
    double ours[3];
    myb_aberration(p, v, C, ours);
    double vc[3] = {v[0] / C, v[1] / C, v[2] / C};
    const double bm1 = sqrt(1.0 - (vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]));
    double erfa[3];
    eraAb(p, vc, 1e300, bm1, erfa);
    for (int i = 0; i < 3; ++i) {
      worst_ab = fmax(worst_ab, fabs(ours[i] - erfa[i]));
    }
    /* deflection by the Sun, observer at ~1 au, source at infinity in direction p */
    double e[3];
    for (int i = 0; i < 3; ++i) {
      e[i] = 2.0 * uniform(&rng) - 1.0;
    }
    const double ne = norm(e);
    for (int i = 0; i < 3; ++i) {
      e[i] /= ne;
    }
    const double em = (0.9 + 0.2 * uniform(&rng)) * au;
    double ld_ours[3];
    double ld_erfa[3];
    myb_light_deflection(GM_SUN, C, 1.0, p, p, e, em, 1e-6, ld_ours);
    eraLd(1.0, p, p, e, em / au, 1e-6, ld_erfa);
    const double defl = norm((double[3]){ld_erfa[0] - p[0], ld_erfa[1] - p[1], ld_erfa[2] - p[2]});
    for (int i = 0; i < 3; ++i) {
      /* ERFA's SRS 1.97412574336e-8 au against 2 GM_sun / c^2 / au: 9e-9 relative */
      worst_ld = fmax(worst_ld, fabs(ld_ours[i] - ld_erfa[i]) / fmax(defl, 1e-300));
    }
  }
  printf("vs ERFA: aberration %.2e, deflection %.2e (relative to the deflection)\n", worst_ab,
         worst_ld);
  MT_CHECK(worst_ab < 5e-16);
  MT_CHECK(worst_ld < 2e-8);
}

/* ---- The relativistic delay against quadrature --------------------------------------------- */

static void test_shapiro(void) {
  /* (1 + gamma) gm / c^3 int dl / r along the straight segment, Gauss-Legendre on 4000 panels */
  static const double xg[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                               0.8611363115940526};
  static const double wg[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                               0.3478548451374538};
  mt_rng rng = {77};
  double worst = 0.0;
  for (int k = 0; k < 50; ++k) {
    double x1[3];
    double x2[3];
    const double xb[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
      x1[i] = (2.0 * uniform(&rng) - 1.0) * 2e11;
      x2[i] = (2.0 * uniform(&rng) - 1.0) * 2e11;
    }
    const double d[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
    const double len = norm(d);
    double integral = 0.0;
    const int panels = 4000;
    for (int pnl = 0; pnl < panels; ++pnl) {
      for (int g = 0; g < 4; ++g) {
        const double s = (pnl + 0.5 * (1.0 + xg[g])) / panels;
        const double x[3] = {x1[0] + s * d[0], x1[1] + s * d[1], x1[2] + s * d[2]};
        integral += wg[g] * 0.5 / panels / norm(x);
      }
    }
    const double quad = 2.0 * GM_SUN / (C * C * C) * len * integral;
    /* myb's form with the (1 + gamma) gm / c^2 terms removed (gm -> 0 inside the logarithm) */
    const double r1 = norm(x1);
    const double r2 = norm(x2);
    const double closed = 2.0 * GM_SUN / (C * C * C) * log((r1 + r2 + len) / (r1 + r2 - len));
    const double ours = myb_shapiro_delay(GM_SUN, C, 1.0, x1, x2, xb);
    worst = fmax(worst, fabs(closed - quad) / quad);
    /* the second-order terms change it by k / (r1 + r2 - r12) relative at most */
    MT_CHECK(fabs(ours - closed) <= 2.0 * 2.0 * GM_SUN / (C * C) / (r1 + r2 - len) * closed + 1e-20);
  }
  printf("relativistic delay: closed form vs quadrature %.2e (relative)\n", worst);
  MT_CHECK(worst < 1e-10);
  /* grazing the solar limb between 1 au and 9.5 au (superior conjunction of Saturn) */
  const double x1[3] = {-149597870700.0, 6.96e8, 0.0};
  const double x2[3] = {9.5 * 149597870700.0, 6.96e8, 0.0};
  const double xb[3] = {0.0, 0.0, 0.0};
  const double grazing = myb_shapiro_delay(GM_SUN, C, 1.0, x1, x2, xb);
  /* 2 GM / c^3 ln(4 r1 r2 / b^2) for b << r */
  const double approx = 2.0 * GM_SUN / (C * C * C) *
                        log(4.0 * norm(x1) * norm(x2) / (6.96e8 * 6.96e8));
  printf("one-way delay grazing the Sun, 1 -> 9.5 au: %.3f us (2 GM/c^3 ln(4 r1 r2 / b^2) = "
         "%.3f us)\n",
         grazing * 1e6, approx * 1e6);
  MT_CHECK(fabs(grazing - approx) < 1e-3 * approx);
}

/* ---- Range rate against differences ------------------------------------------------------- */

typedef struct moving {
  double x0[3];
  double v[3];
  double a[3];
} moving;

static myb_status moving_state(void* context, myb_time t, double state[6]) {
  const moving* m = context;
  const double s = myb_time_seconds(t);
  for (int i = 0; i < 3; ++i) {
    state[i] = m->x0[i] + m->v[i] * s + 0.5 * m->a[i] * s * s;
    state[3 + i] = m->v[i] + m->a[i] * s;
  }
  return MYB_OK;
}

static void test_range_rate(void) {
  moving target = {{3e11, -1e11, 5e10}, {2.5e4, 1.8e4, -3e3}, {1e-3, -2e-3, 5e-4}};
  moving observer = {{1e11, 5e10, 0.0}, {-1e4, 2.9e4, 1e3}, {-5e-3, 1e-3, 0.0}};
  const myb_light_model model = {.c = C, .gamma = 1.0};
  double worst = 0.0;
  for (int dir = 0; dir < 2; ++dir) {
    const myb_time t = myb_time_from_seconds(1e5);
    double obs[6];
    MT_CHECK(moving_state(&observer, t, obs) == MYB_OK);
    myb_light_path path;
    MT_CHECK(myb_light_time(&model, moving_state, &target, t, obs, (myb_light_direction)dir,
                            &path) == MYB_OK);
    const double h = 10.0;
    double r[2];
    for (int s = 0; s < 2; ++s) {
      const myb_time ts = myb_time_add(t, s == 0 ? h : -h);
      double o[6];
      MT_CHECK(moving_state(&observer, ts, o) == MYB_OK);
      myb_light_path p2;
      MT_CHECK(myb_light_time(&model, moving_state, &target, ts, o, (myb_light_direction)dir,
                              &p2) == MYB_OK);
      r[s] = p2.range;
    }
    const double diff = (r[0] - r[1]) / (2.0 * h);
    worst = fmax(worst, fabs(path.range_rate - diff) / fabs(diff));
    printf("range rate (%s): %.9f m/s, difference %.9f m/s, %d iterations\n",
           dir == 0 ? "reception" : "transmission", path.range_rate, diff, path.iterations);
  }
  MT_CHECK(worst < 1e-9);
  double ra;
  double dec;
  myb_direction_radec((double[3]){0.0, -1.0, 1.0}, &ra, &dec);
  MT_CHECK(fabs(ra - 1.5 * 3.141592653589793) < 1e-15 && fabs(dec - 0.25 * 3.141592653589793) < 1e-15);
  const myb_light_model bad = {.c = 0.0};
  myb_light_path path;
  double obs[6] = {0};
  MT_CHECK(myb_light_time(&bad, moving_state, &target, myb_time_from_seconds(0.0), obs,
                          MYB_LIGHT_RECEPTION, &path) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_erfa();
  test_shapiro();
  test_range_rate();
  test_spice();
  return mt_finish("test_observe");
}
