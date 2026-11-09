/* SPDX-License-Identifier: BSD-3-Clause */
/* JPL's acceptance test of an ephemeris reader: testpo.440 lists positions and velocities (au,
 * au/day) of the planets, the Moon, the Sun, the barycentres, and the lunar libration angles
 * (rad, rad/day) of DE440 at 13 000 dates from 1550 to 2650, computed by JPL from the same
 * Chebyshev coefficients. myb reads them from NAIF's SPK conversion de440.bsp and the lunar
 * orientation from moon_pa_de440_200625.bpc (data/download, see data/README.md; the test is
 * skipped without them). Nutations (target 14) are not part of the SPK conversion. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/kin.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

static const double au_km = 149597870.7; /* DE440 header constant AU */

/* NAIF identifiers of the testpo bodies 1 .. 13. */
static const int32_t naif[14] = {0, 1, 2, 399, 4, 5, 6, 7, 8, 9, 301, 10, 0, 3};

static int next_number(const char** p, double* value) {
  /* bound: skips at most the rest of the line */
  while (**p == ' ') {
    ++*p;
  }
  size_t used = 0;
  if (myb_parse_double(*p, strlen(*p), value, &used) != MYB_OK) {
    return 0;
  }
  *p += used;
  return 1;
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  FILE* f = fopen(DOWNLOAD "testpo.440", "r");
  if (f == NULL) {
    printf("testpo.440 not downloaded: skipped\n");
    return MTEST_SKIP;
  }
  myb_spk* spk = NULL;
  myb_pck* pck = NULL;
  myb_diag diag;
  MT_CHECK(myb_spk_create(NULL, &spk) == MYB_OK);
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  if (myb_spk_load_file(spk, DOWNLOAD "de440.bsp", &diag) != MYB_OK) {
    printf("de440.bsp not available (%s): skipped\n", diag.message);
    fclose(f);
    myb_spk_destroy(spk);
    myb_pck_destroy(pck);
    return MTEST_SKIP;
  }
  const int have_moon = myb_pck_load_file(pck, DOWNLOAD "moon_pa_de440_200625.bpc", &diag) ==
                        MYB_OK;
  if (!have_moon) {
    printf("moon_pa_de440_200625.bpc not available: librations not checked\n");
  }
  char line[256];
  int in_data = 0;
  long checked = 0;
  long librations = 0;
  long nutations = 0;
  double worst_position = 0.0; /* au */
  double worst_velocity = 0.0; /* au/day */
  double worst_angle = 0.0;    /* rad */
  double worst_rate = 0.0;     /* rad/day */
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    if (!in_data) {
      in_data = strncmp(line, "EOT", 3) == 0;
      continue;
    }
    int de = 0;
    char date[16];
    int consumed = 0;
    if (sscanf(line, "%d %15s%n", &de, date, &consumed) != 2 || de != 440) {
      continue;
    }
    const char* p = line + consumed;
    double jd = 0.0;
    double target = 0.0;
    double center = 0.0;
    double coordinate = 0.0;
    double expected = 0.0;
    MT_CHECK(next_number(&p, &jd) && next_number(&p, &target) && next_number(&p, &center) &&
             next_number(&p, &coordinate) && next_number(&p, &expected));
    const int t = (int)target;
    const int c = (int)center;
    const int k = (int)coordinate - 1;
    const myb_time when = myb_time_from_jd(jd, 0.0);
    if (t == 14) {
      ++nutations;
      continue;
    }
    if (t == 15) {
      if (!have_moon) {
        continue;
      }
      double angles[3];
      double rates[3];
      int32_t base = -1;
      MT_CHECK(myb_pck_euler(pck, 31008, when, angles, rates, &base) == MYB_OK && base == 1);
      if (k < 3) {
        /* the angles may be stored reduced: compare modulo 2 pi */
        const double two_pi = 6.283185307179586;
        double d = fmod(angles[k] - expected, two_pi);
        if (d > 0.5 * two_pi) {
          d -= two_pi;
        } else if (d < -0.5 * two_pi) {
          d += two_pi;
        }
        worst_angle = fmax(worst_angle, fabs(d) / fmax(1.0, fabs(expected)));
      } else {
        worst_rate = fmax(worst_rate, fabs(rates[k - 3] * 86400.0 - expected));
      }
      ++librations;
      continue;
    }
    MT_CHECK(t >= 1 && t <= 13 && c >= 0 && c <= 13 && k >= 0 && k < 6);
    double state[6];
    int32_t frame = -1;
    const int32_t center_id = c == 0 ? 0 : naif[c];
    const myb_status status = myb_spk_state(spk, naif[t], center_id, when, state, &frame);
    MT_CHECK(status == MYB_OK && (frame == 1 || naif[t] == center_id));
    if (k < 3) {
      worst_position = fmax(worst_position, fabs(state[k] / 1000.0 / au_km - expected));
    } else {
      worst_velocity = fmax(worst_velocity, fabs(state[k] / 1000.0 / au_km * 86400.0 - expected));
    }
    ++checked;
  }
  fclose(f);
  printf("states: %ld values, max |error| %.2e au (position), %.2e au/day (velocity)\n", checked,
         worst_position, worst_velocity);
  printf("librations: %ld values, max |error| %.2e (angles, relative to max(1 rad, |angle|)), "
         "%.2e rad/day (rates)\n",
         librations, worst_angle, worst_rate);
  printf("nutations: %ld values not in the SPK conversion\n", nutations);
  MT_CHECK(checked > 11000);
  /* JPL's testeph accepts differences below 1e-13 (au, au/day, rad). The libration angle psi
   * accumulates (|psi| ~ 3.5e4 rad, one ulp = 7.3e-12 rad), so angles are compared relative to
   * their magnitude: measured 4.3e-16, about one ulp. */
  MT_CHECK(worst_position < 1e-13 && worst_velocity < 1e-13);
  if (have_moon) {
    MT_CHECK(librations > 900);
    MT_CHECK(worst_angle < 1e-13 && worst_rate < 1e-13);
  }
  myb_spk_destroy(spk);
  myb_pck_destroy(pck);
  return mt_finish("test_testpo");
}
