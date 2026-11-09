/* SPDX-License-Identifier: BSD-3-Clause */
/* The CIO-based Earth orientation of myb, driven by data/frames/iau2006-cip.tk (IERS
 * Conventions tables), against the IAU SOFA implementation through ERFA: X, Y (eraXy06), s
 * (eraS06), ERA (eraEra00), s' (eraSp00), the celestial-to-terrestrial matrix (eraC2ixys,
 * eraPom00, eraC2tcio); and the analytic time derivatives against finite differences. */
#include <math.h>
#include <stdio.h>

#include "erfa.h"
#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"

static myb_cip* load_cip(void) {
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const myb_status status =
      myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau2006-cip.tk", &diag);
  if (status != MYB_OK) {
    fprintf(stderr, "iau2006-cip.tk: %s (line %ld: %s)\n", myb_status_name(status), diag.line,
            diag.message);
  }
  MT_CHECK(status == MYB_OK);
  myb_cip* cip = NULL;
  MT_CHECK(myb_cip_create(pool, "IAU2006A", NULL, &cip, &diag) == MYB_OK);
  myb_kpool_destroy(pool);
  return cip;
}

/* An instant as myb_time and as an ERFA two-part Julian date (whole days, fraction). */
typedef struct instant {
  myb_time t;
  double jd1;
  double jd2;
} instant;

static instant make(int64_t seconds, double fraction) {
  instant i;
  i.t = myb_time_add((myb_time){seconds, 0.0}, fraction);
  const int64_t day = i.t.s >= 0 ? i.t.s / 86400 : -((-i.t.s + 86399) / 86400);
  i.jd1 = 2451545.0 + (double)day;
  i.jd2 = ((double)(i.t.s - day * 86400) + i.t.f) / 86400.0;
  return i;
}

static void test_against_erfa(const myb_cip* cip) {
  mt_rng rng = {20261109};
  double worst_xy = 0.0;
  double worst_s = 0.0;
  double worst_era = 0.0;
  double worst_sp = 0.0;
  double worst_matrix = 0.0;
  const int samples = 3000;
  for (int k = 0; k < samples; ++k) {
    /* 1800 .. 2200 */
    const int64_t seconds = (int64_t)mt_int(&rng, -1000000, 1000000) * 6311;
    const double fraction = (double)mt_int(&rng, 0, 1 << 30) / (double)(1 << 30);
    const instant tt = make(seconds, fraction);
    const instant ut1 = make(seconds - 69, fraction + 0.184);
    double xys[3];
    myb_cip_xys(cip, tt.t, xys, NULL);
    double x;
    double y;
    eraXy06(tt.jd1, tt.jd2, &x, &y);
    const double s = eraS06(tt.jd1, tt.jd2, x, y);
    worst_xy = fmax(worst_xy, fmax(fabs(xys[0] - x), fabs(xys[1] - y)));
    worst_s = fmax(worst_s, fabs(xys[2] - s));
    const double era = myb_cip_era(cip, ut1.t, NULL);
    double d_era = fabs(era - eraEra00(ut1.jd1, ut1.jd2));
    d_era = fmin(d_era, fabs(d_era - 6.283185307179586));
    worst_era = fmax(worst_era, d_era);
    worst_sp = fmax(worst_sp, fabs(myb_cip_sprime(cip, tt.t) - eraSp00(tt.jd1, tt.jd2)));

    const myb_eop_values eop = {.xp = 0.2e-6 * (k % 7), .yp = 1.5e-6, .dx = 1e-9, .dy = -2e-9};
    myb_rotation m;
    myb_cip_terrestrial(cip, tt.t, ut1.t, &eop, &m);
    double rc2i[3][3];
    double rpom[3][3];
    double rc2t[3][3];
    eraC2ixys(x + eop.dx, y + eop.dy, s, rc2i);
    eraPom00(eop.xp, eop.yp, eraSp00(tt.jd1, tt.jd2), rpom);
    eraC2tcio(rc2i, eraEra00(ut1.jd1, ut1.jd2), rpom, rc2t);
    for (int i = 0; i < 9; ++i) {
      worst_matrix = fmax(worst_matrix, fabs(m.m[i] - rc2t[i / 3][i % 3]));
    }
  }
  printf("vs ERFA over 1800-2200 (%d dates): X,Y %.2e rad, s %.2e rad, ERA %.2e rad, "
         "s' %.2e rad, matrix %.2e\n",
         samples, worst_xy, worst_s, worst_era, worst_sp, worst_matrix);
  /* Measured: X, Y 1.7e-17 rad (3.6e-6 µas), s 1.5e-19, s' 1.6e-25. ERA and the matrix differ
   * by ERFA's rounding: eraEra00 sums up to ~200 turns in one double (one ulp = 1.8e-13 rad),
   * while myb_cip_era is exact to 3.5e-16 rad (test_mpfr_frames). */
  MT_CHECK(worst_xy < 5e-17);
  MT_CHECK(worst_s < 1e-18);
  MT_CHECK(worst_era < 6e-13);
  MT_CHECK(worst_sp < 1e-20);
  MT_CHECK(worst_matrix < 6e-13);
}

static void test_derivatives(const myb_cip* cip) {
  mt_rng rng = {7};
  double worst_rate = 0.0;
  double worst_dm = 0.0;
  for (int k = 0; k < 200; ++k) {
    const int64_t seconds = (int64_t)mt_int(&rng, -1000000, 1000000) * 3000;
    const myb_time t = {seconds, 0.25};
    const double h = 50.0;
    double xys[3];
    double rates[3];
    double plus[3];
    double minus[3];
    myb_cip_xys(cip, t, xys, rates);
    myb_cip_xys(cip, myb_time_add(t, h), plus, NULL);
    myb_cip_xys(cip, myb_time_add(t, -h), minus, NULL);
    for (int i = 0; i < 3; ++i) {
      worst_rate = fmax(worst_rate, fabs(rates[i] - (plus[i] - minus[i]) / (2.0 * h)));
    }
    const myb_eop_values eop = {.xp = 1e-6, .yp = 2e-6, .lod = 1.5e-3};
    myb_rotation m;
    myb_rotation mp;
    myb_rotation mm;
    const double g = 0.5;
    /* UT1 advances by (1 - LOD/86400) per TT second */
    const double k_ut1 = 1.0 - eop.lod / 86400.0;
    myb_cip_terrestrial(cip, t, t, &eop, &m);
    myb_cip_terrestrial(cip, myb_time_add(t, g), myb_time_add(t, g * k_ut1), &eop, &mp);
    myb_cip_terrestrial(cip, myb_time_add(t, -g), myb_time_add(t, -g * k_ut1), &eop, &mm);
    for (int i = 0; i < 9; ++i) {
      worst_dm = fmax(worst_dm, fabs(m.dm[i] - (mp.m[i] - mm.m[i]) / (2.0 * g)));
    }
  }
  printf("derivatives vs central differences: X,Y,s %.2e rad/s, matrix %.2e /s\n", worst_rate,
         worst_dm);
  MT_CHECK(worst_rate < 1e-18);
  MT_CHECK(worst_dm < 1e-13);
}

/* Vermeille's closed form on WGS-84 from deep inside the Earth to beyond the Moon: the
 * geodetic coordinates that generated each point (through ERFA's eraGd2gce) come back, and
 * near the surface they agree with ERFA's eraGc2gde (Fukushima 2006). Far from the surface
 * eraGc2gde loses accuracy (2.6e-11 rad at 1.3e7 .. 1.9e7 m), so it is not the reference
 * there. */
static void test_geodetic(void) {
  const myb_ellipsoid wgs84 = {6378137.0, 1.0 / 298.257223563};
  mt_rng rng = {99};
  double worst_angle = 0.0;
  double worst_height = 0.0;
  double worst_erfa = 0.0;
  double worst_cartesian = 0.0;
  double worst_inner = 0.0;
  for (int k = 0; k < 20000; ++k) {
    const double lon = (mt_int(&rng, 0, 1 << 30) / (double)(1 << 30) - 0.5) * 6.283185307179586;
    double lat = (mt_int(&rng, 0, 1 << 30) / (double)(1 << 30) - 0.5) * 3.141592653589793;
    if (k % 50 == 0) {
      lat = (k % 100 == 0 ? 1.0 : -1.0) * 1.5707963267948966; /* poles */
    } else if (k % 50 == 1) {
      lat = 0.0;
    }
    const int band = k % 4;
    const double u = mt_int(&rng, 0, 1 << 30) / (double)(1 << 30);
    const double h = band == 0 ? -6.0e6 * u : band == 1 ? 2.0e4 * (u - 0.5)
                                         : band == 2 ? 4.0e7 * u : 4.0e8 * u;
    double xyz[3];
    MT_CHECK(myb_geodetic_to_cartesian(&wgs84, lon, lat, h, xyz) == MYB_OK);
    double e_xyz[3];
    MT_CHECK(eraGd2gce(wgs84.a, wgs84.f, lon, lat, h, e_xyz) == 0);
    for (int i = 0; i < 3; ++i) {
      worst_cartesian = fmax(worst_cartesian, fabs(xyz[i] - e_xyz[i]) / (wgs84.a + fabs(h)));
    }
    double lo;
    double la;
    double he;
    MT_CHECK(myb_geodetic_from_cartesian(&wgs84, e_xyz, &lo, &la, &he) == MYB_OK);
    if (band == 0) {
      /* inside, the nearest surface point may differ from the generating one: round trip */
      double back[3];
      MT_CHECK(myb_geodetic_to_cartesian(&wgs84, lo, la, he, back) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        worst_inner = fmax(worst_inner, fabs(back[i] - e_xyz[i]));
      }
      continue;
    }
    double dlon = fabs(lo - lon);
    dlon = fmin(dlon, fabs(dlon - 6.283185307179586));
    if (fabs(lat) > 1.5707963) {
      dlon = 0.0; /* the longitude of a pole is arbitrary */
    }
    worst_angle = fmax(worst_angle, fmax(dlon * cos(lat), fabs(la - lat)));
    worst_height = fmax(worst_height, fabs(he - h) / (wgs84.a + fabs(h)));
    if (band == 1) {
      double elon;
      double phi;
      double height;
      MT_CHECK(eraGc2gde(wgs84.a, wgs84.f, e_xyz, &elon, &phi, &height) == 0);
      worst_erfa = fmax(worst_erfa, fmax(fabs(la - phi), fabs(he - height) / wgs84.a));
    }
  }
  printf("geodetic: angles %.2e rad, height %.2e (relative to a + |h|), vs eraGc2gde near "
         "the surface %.2e, to_cartesian vs eraGd2gce %.2e, interior round trip %.2e m\n",
         worst_angle, worst_height, worst_erfa, worst_cartesian, worst_inner);
  MT_CHECK(worst_angle < 2e-15 && worst_height < 1e-15 && worst_erfa < 2e-15 &&
           worst_cartesian < 1e-15 && worst_inner < 1e-8);
}

int main(void) {
  test_geodetic();
  myb_cip* cip = load_cip();
  if (cip == NULL) {
    return mt_finish("test_erfa_frames");
  }
  test_against_erfa(cip);
  test_derivatives(cip);
  myb_cip_destroy(cip);
  return mt_finish("test_erfa_frames");
}
