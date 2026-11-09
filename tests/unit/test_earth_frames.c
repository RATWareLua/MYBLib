/* SPDX-License-Identifier: BSD-3-Clause */
/* The whole terrestrial chain against JPL Horizons: geodetic site on WGS-84 -> ITRS ->
 * (EOP 20 C04 with tidal corrections, IAU 2006/2000A CIO-based rotation, TT and UT1 from the
 * time scales) -> GCRS, compared with the site vectors Horizons implies (the Moon relative to
 * the Earth's centre minus the Moon relative to the site, tests/data/frames). And the orbital
 * frames: orthonormality, known axes, derivatives against finite differences along a conic. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"
#include "myb/kepler.h"
#include "myb/time.h"

static const double deg = 0.017453292519943295;

typedef struct site_row {
  double jd;
  double s[6];
} site_row;

static int read_rows(const char* tag, site_row* rows, int max) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/frames/horizons_moon_sites.txt", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return 0;
  }
  char line[512];
  int n = 0;
  /* bound: one iteration per line of the file */
  while (fgets(line, sizeof line, f) != NULL && n < max) {
    if (line[0] == '#' || strncmp(line, tag, strlen(tag)) != 0 || line[strlen(tag)] != ',') {
      continue;
    }
    const char* p = line + strlen(tag) + 1;
    site_row r;
    int used = 0;
    if (sscanf(p, "%lf,%n", &r.jd, &used) != 1) {
      continue;
    }
    p = strchr(p + used, ','); /* skip the calendar label */
    for (int i = 0; i < 6 && p != NULL; ++i) {
      MT_CHECK(sscanf(p + 1, "%lf%n", &r.s[i], &used) == 1);
      p = strchr(p + 1, ',');
    }
    rows[n++] = r;
  }
  fclose(f);
  return n;
}

static void test_horizons_sites(void) {
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const char* files[] = {MYB_DATA_DIR "/time/time-scales.tk", MYB_DATA_DIR "/time/tdb-tt-fb1990.tk",
                         MYB_DATA_DIR "/time/utc-tai.tk", MYB_DATA_DIR "/frames/iau2006-cip.tk",
                         MYB_DATA_DIR "/frames/eop-tides.tk"};
  for (size_t i = 0; i < sizeof files / sizeof files[0]; ++i) {
    MT_CHECK(myb_kpool_load_file(pool, files[i], &diag) == MYB_OK);
  }
  myb_timesys* sys = NULL;
  myb_cip* cip = NULL;
  myb_eop* eop = NULL;
  MT_CHECK(myb_timesys_create(pool, NULL, &sys, &diag) == MYB_OK);
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &cip, &diag) == MYB_OK);
  MT_CHECK(myb_eop_create(NULL, &eop) == MYB_OK);
  MT_CHECK(myb_eop_load_c04_file(eop, MYB_TEST_DATA_DIR "/eop/eopc04-excerpt.txt", &diag) ==
           MYB_OK);
  MT_CHECK(myb_eop_set_tides(eop, pool, "IERS10_EOPT", &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  myb_scale tdb = {0};
  myb_scale tt = {0};
  MT_CHECK(myb_timesys_scale(sys, "TDB", &tdb) == MYB_OK);
  MT_CHECK(myb_timesys_scale(sys, "TT", &tt) == MYB_OK);

  site_row geo[16];
  site_row s1[16];
  site_row s2[16];
  const int n = read_rows("GEO", geo, 16);
  MT_CHECK(n == 10 && read_rows("S1", s1, 16) == n && read_rows("S2", s2, 16) == n);
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  const double sites[2][3] = {{30.3261 * deg, 59.7717 * deg, 75.0},
                              {-116.8895 * deg, -35.4 * deg, 1200.0}};
  double worst_r = 0.0;
  double worst_v = 0.0;
  for (int k = 0; k < n; ++k) {
    const myb_time t_tdb = myb_time_from_jd(geo[k].jd, 0.0);
    myb_time t_tt;
    MT_CHECK(myb_timesys_convert(sys, t_tdb, tdb, tt, &t_tt) == MYB_OK);
    myb_eop_values v;
    myb_time ut1;
    MT_CHECK(myb_eop_at(eop, sys, t_tt, tt, &v, &ut1) == MYB_OK);
    myb_rotation c2t;
    myb_cip_terrestrial(cip, t_tt, ut1, &v, &c2t);
    myb_rotation t2c;
    myb_rotation_invert(&c2t, &t2c);
    for (int s = 0; s < 2; ++s) {
      const site_row* topo = s == 0 ? &s1[k] : &s2[k];
      MT_CHECK(topo->jd == geo[k].jd);
      double itrs[6] = {0};
      MT_CHECK(myb_geodetic_to_cartesian(&wgs84, sites[s][0], sites[s][1], sites[s][2], itrs) ==
               MYB_OK);
      double gcrs[6];
      myb_rotation_apply_state(&t2c, itrs, gcrs);
      double dr = 0.0;
      double dv = 0.0;
      for (int i = 0; i < 3; ++i) {
        const double hr = (geo[k].s[i] - topo->s[i]) * 1000.0;
        const double hv = (geo[k].s[i + 3] - topo->s[i + 3]) * 1000.0;
        dr += (gcrs[i] - hr) * (gcrs[i] - hr);
        dv += (gcrs[i + 3] - hv) * (gcrs[i + 3] - hv);
      }
      worst_r = fmax(worst_r, sqrt(dr));
      worst_v = fmax(worst_v, sqrt(dv));
    }
  }
  printf("site vectors vs Horizons (2 sites x %d dates of 2020), IERS 2010 chain: position "
         "%.3f m, velocity %.2e m/s\n",
         n, worst_r, worst_v);
  /* Horizons orients the Earth with NAIF's ITRF93 model (IAU 1976/1980 precession-nutation
   * with the EOP's nutation corrections, J2000 taken as the ICRF); myb follows IERS 2010
   * (IAU 2006/2000A, CIO-based, frame bias included). The difference is a rotation of 3 to
   * 11 mas (no radial part): 0.19 m at most. The ITRF93 PCK itself reproduces Horizons to
   * 0.08 m (below). */
  MT_CHECK(worst_r < 0.3 && worst_v < 2e-5);

  /* The same sites through NAIF's binary PCK of ITRF93 relative to ECLIPJ2000 (frame 17 of
   * data/frames/naif-inertial.tk), when downloaded. */
  myb_pck* pck = NULL;
  MT_CHECK(myb_pck_create(NULL, &pck) == MYB_OK);
  if (myb_pck_load_file(pck, MYB_DATA_DIR "/download/earth_latest_high_prec.bpc", &diag) !=
      MYB_OK) {
    printf("earth_latest_high_prec.bpc not downloaded: ITRF93 comparison skipped\n");
  } else {
    myb_kpool* frames = NULL;
    MT_CHECK(myb_kpool_create(NULL, &frames) == MYB_OK);
    MT_CHECK(myb_kpool_load_file(frames, MYB_DATA_DIR "/frames/naif-inertial.tk", &diag) ==
             MYB_OK);
    myb_rotation ecliptic;
    char relative[33];
    MT_CHECK(myb_tkframe(frames, 17, &ecliptic, relative, &diag) == MYB_OK);
    myb_kpool_destroy(frames);
    double pck_r = 0.0;
    double pck_v = 0.0;
    for (int k = 0; k < n; ++k) {
      myb_rotation itrf93;
      int32_t base = 0;
      MT_CHECK(myb_pck_rotation(pck, 3000, myb_time_from_jd(geo[k].jd, 0.0), &itrf93, &base) ==
                   MYB_OK &&
               base == 17);
      myb_rotation j2000_to_itrf93;
      myb_rotation_compose(&itrf93, &ecliptic, &j2000_to_itrf93);
      myb_rotation back;
      myb_rotation_invert(&j2000_to_itrf93, &back);
      for (int s = 0; s < 2; ++s) {
        const site_row* topo = s == 0 ? &s1[k] : &s2[k];
        double itrs[6] = {0};
        MT_CHECK(myb_geodetic_to_cartesian(&wgs84, sites[s][0], sites[s][1], sites[s][2],
                                           itrs) == MYB_OK);
        double g[6];
        myb_rotation_apply_state(&back, itrs, g);
        double dr = 0.0;
        double dv = 0.0;
        for (int i = 0; i < 3; ++i) {
          const double hr = (geo[k].s[i] - topo->s[i]) * 1000.0;
          const double hv = (geo[k].s[i + 3] - topo->s[i + 3]) * 1000.0;
          dr += (g[i] - hr) * (g[i] - hr);
          dv += (g[i + 3] - hv) * (g[i + 3] - hv);
        }
        pck_r = fmax(pck_r, sqrt(dr));
        pck_v = fmax(pck_v, sqrt(dv));
      }
    }
    printf("site vectors vs Horizons through the ITRF93 binary PCK: position %.3f m, velocity "
           "%.2e m/s\n",
           pck_r, pck_v);
    /* Horizons' EOP file (eop.261106) is not the one the PCK was made from, and the PCK is a
     * Chebyshev fit of the EOP model. */
    MT_CHECK(pck_r < 0.1 && pck_v < 1e-5);
  }
  myb_pck_destroy(pck);
  myb_eop_destroy(eop);
  myb_cip_destroy(cip);
  myb_timesys_destroy(sys);
}

static void test_orbital(void) {
  const double mu = 3.986004418e14;
  const myb_pv start = {{7.0e6, -1.2e6, 0.4e6}, {1.1e3, 7.2e3, 1.9e3}};
  const myb_orbital_frame kinds[3] = {MYB_ORBITAL_RTN, MYB_ORBITAL_VNC, MYB_ORBITAL_LVLH};
  double worst_orth = 0.0;
  double worst_rate = 0.0;
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 20; ++j) {
      const double t = 300.0 * j;
      myb_pv pv;
      myb_pv pp;
      myb_pv pm;
      const double h = 0.5;
      MT_CHECK(myb_kepler_propagate(mu, &start, t, &pv) == MYB_OK);
      MT_CHECK(myb_kepler_propagate(mu, &start, t + h, &pp) == MYB_OK);
      MT_CHECK(myb_kepler_propagate(mu, &start, t - h, &pm) == MYB_OK);
      double state[6];
      double acc[3];
      const double r = sqrt(pv.r[0] * pv.r[0] + pv.r[1] * pv.r[1] + pv.r[2] * pv.r[2]);
      for (int i = 0; i < 3; ++i) {
        state[i] = pv.r[i];
        state[i + 3] = pv.v[i];
        acc[i] = -mu * pv.r[i] / (r * r * r);
      }
      myb_rotation m;
      MT_CHECK(myb_orbital_rotation(kinds[k], state, acc, &m) == MYB_OK);
      /* M M^T = I */
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          double dot = 0.0;
          for (int c = 0; c < 3; ++c) {
            dot += m.m[3 * a + c] * m.m[3 * b + c];
          }
          worst_orth = fmax(worst_orth, fabs(dot - (a == b ? 1.0 : 0.0)));
        }
      }
      /* the position in RTN is (r, 0, 0); in LVLH (0, 0, -r) */
      double local[3];
      myb_rotation_apply_position(&m, state, local);
      if (kinds[k] == MYB_ORBITAL_RTN) {
        MT_CHECK(fabs(local[0] - r) < 1e-8 && fabs(local[1]) < 1e-8 && fabs(local[2]) < 1e-8);
      } else if (kinds[k] == MYB_ORBITAL_LVLH) {
        MT_CHECK(fabs(local[2] + r) < 1e-8 && fabs(local[0]) < 1e-8 && fabs(local[1]) < 1e-8);
      }
      double sp[6];
      double sm[6];
      double ap[3];
      double am[3];
      const double rp = sqrt(pp.r[0] * pp.r[0] + pp.r[1] * pp.r[1] + pp.r[2] * pp.r[2]);
      const double rm = sqrt(pm.r[0] * pm.r[0] + pm.r[1] * pm.r[1] + pm.r[2] * pm.r[2]);
      for (int i = 0; i < 3; ++i) {
        sp[i] = pp.r[i];
        sp[i + 3] = pp.v[i];
        sm[i] = pm.r[i];
        sm[i + 3] = pm.v[i];
        ap[i] = -mu * pp.r[i] / (rp * rp * rp);
        am[i] = -mu * pm.r[i] / (rm * rm * rm);
      }
      myb_rotation mp;
      myb_rotation mm;
      MT_CHECK(myb_orbital_rotation(kinds[k], sp, ap, &mp) == MYB_OK);
      MT_CHECK(myb_orbital_rotation(kinds[k], sm, am, &mm) == MYB_OK);
      for (int i = 0; i < 9; ++i) {
        worst_rate = fmax(worst_rate, fabs(m.dm[i] - (mp.m[i] - mm.m[i]) / (2.0 * h)));
      }
    }
  }
  printf("orbital frames: orthonormality %.2e, derivative vs central differences %.2e /s\n",
         worst_orth, worst_rate);
  MT_CHECK(worst_orth < 1e-15 && worst_rate < 1e-10);
  /* degenerate states */
  const double radial[6] = {7e6, 0, 0, 1e3, 0, 0};
  const double zero_acc[3] = {0, 0, 0};
  myb_rotation m;
  MT_CHECK(myb_orbital_rotation(MYB_ORBITAL_RTN, radial, zero_acc, &m) == MYB_E_DOMAIN);
}

int main(void) {
  test_horizons_sites();
  test_orbital();
  return mt_finish("test_earth_frames");
}
