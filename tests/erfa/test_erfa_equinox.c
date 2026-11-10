/* SPDX-License-Identifier: BSD-3-Clause */
/* Equinox-based systems (myb_equinox, data/frames/iau1980.tk) against ERFA: IAU 1976 precession
 * (eraPrec76, eraPmat76), IAU 1980 obliquity and nutation (eraObl80, eraNut80, eraNutm80),
 * GMST 1982 (eraGmst82) over 1800-2200; the TEME axes against their definition; derivatives
 * against central differences. */
#include <math.h>
#include <stdio.h>

#include "erfa.h"
#include "mtest.h"
#include "myb/frames.h"
#include "myb/io.h"

static myb_equinox* load(void) {
  myb_kpool* pool = NULL;
  myb_diag diag;
  MT_CHECK(myb_kpool_create(NULL, &pool) == MYB_OK);
  const myb_status status = myb_kpool_load_file(pool, MYB_DATA_DIR "/frames/iau1980.tk", &diag);
  if (status != MYB_OK) {
    fprintf(stderr, "iau1980.tk: %s (line %ld: %s)\n", myb_status_name(status), diag.line,
            diag.message);
  }
  MT_CHECK(status == MYB_OK);
  myb_equinox* e = NULL;
  const myb_status created = myb_equinox_create(pool, "IAU1980", NULL, &e, &diag);
  if (created != MYB_OK) {
    fprintf(stderr, "myb_equinox_create: %s\n", diag.message);
  }
  MT_CHECK(created == MYB_OK);
  myb_kpool_destroy(pool);
  return e;
}

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

static double angle_difference(double a, double b) {
  double d = fabs(a - b);
  return fmin(d, fabs(d - 6.283185307179586));
}

static void test_against_erfa(const myb_equinox* e) {
  mt_rng rng = {20261110};
  double worst_prec = 0.0;
  double worst_pmat = 0.0;
  double worst_obl = 0.0;
  double worst_nut = 0.0;
  double worst_nmat = 0.0;
  double worst_gmst = 0.0;
  double worst_tod = 0.0;
  const int samples = 3000;
  for (int k = 0; k < samples; ++k) {
    const int64_t seconds = (int64_t)mt_int(&rng, -1000000, 1000000) * 6311;
    const double fraction = (double)mt_int(&rng, 0, 1 << 30) / (double)(1 << 30);
    const instant tt = make(seconds, fraction);
    const instant ut1 = make(seconds - 64, fraction + 0.3);

    double angles[3];
    myb_equinox_precession(e, tt.t, angles, NULL);
    double zeta;
    double z;
    double theta;
    eraPrec76(2451545.0, 0.0, tt.jd1, tt.jd2, &zeta, &z, &theta);
    worst_prec = fmax(worst_prec, fmax(fabs(angles[0] - zeta),
                                       fmax(fabs(angles[1] - z), fabs(angles[2] - theta))));
    myb_rotation p;
    myb_equinox_rotation(e, MYB_EQUINOX_MOD, tt.t, &p);
    double rp[3][3];
    eraPmat76(tt.jd1, tt.jd2, rp);
    for (int i = 0; i < 9; ++i) {
      worst_pmat = fmax(worst_pmat, fabs(p.m[i] - rp[i / 3][i % 3]));
    }

    worst_obl = fmax(worst_obl, fabs(myb_equinox_obliquity(e, tt.t, NULL) -
                                     eraObl80(tt.jd1, tt.jd2)));
    double nut[2];
    myb_equinox_nutation(e, tt.t, nut, NULL);
    double dpsi;
    double deps;
    eraNut80(tt.jd1, tt.jd2, &dpsi, &deps);
    worst_nut = fmax(worst_nut, fmax(fabs(nut[0] - dpsi), fabs(nut[1] - deps)));

    myb_rotation np;
    myb_equinox_rotation(e, MYB_EQUINOX_TOD, tt.t, &np);
    double rn[3][3];
    eraNutm80(tt.jd1, tt.jd2, rn);
    double rnp[3][3];
    eraRxr(rn, rp, rnp);
    for (int i = 0; i < 9; ++i) {
      worst_tod = fmax(worst_tod, fabs(np.m[i] - rnp[i / 3][i % 3]));
    }
    /* N alone: TOD P^T */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        double s = 0.0;
        for (int m = 0; m < 3; ++m) {
          s += np.m[3 * i + m] * p.m[3 * j + m];
        }
        worst_nmat = fmax(worst_nmat, fabs(s - rn[i][j]));
      }
    }

    worst_gmst = fmax(worst_gmst, angle_difference(myb_equinox_gmst(e, ut1.t, NULL),
                                                   eraGmst82(ut1.jd1, ut1.jd2)));
  }
  printf("vs ERFA over 1800-2200 (%d dates): precession angles %.2e rad, P %.2e, obliquity "
         "%.2e rad, nutation %.2e rad, N %.2e, NP %.2e, GMST82 %.2e rad\n",
         samples, worst_prec, worst_pmat, worst_obl, worst_nut, worst_nmat, worst_tod,
         worst_gmst);
  MT_CHECK(worst_prec < 1e-16);
  MT_CHECK(worst_pmat < 1e-15);
  MT_CHECK(worst_obl < 1e-16);
  /* the arguments: ERFA reduces the whole revolutions separately (fmod(1325 t, 1)), the kernel
   * folds them into the linear coefficient: 1e-7 arcsec rounding of a coefficient over 2 centuries */
  MT_CHECK(worst_nut < 2e-15);
  MT_CHECK(worst_nmat < 2e-15);
  MT_CHECK(worst_tod < 2e-15);
  /* both sum seconds of up to 1.7e7 in one double (half an ulp: 1.9e-9 s = 1.4e-13 rad) */
  MT_CHECK(worst_gmst < 6e-13);
}

static void test_teme(const myb_equinox* e) {
  /* TEME: orthonormal, z the TOD pole, x in the plane of the TOD pole and the mean equinox, and
   * the rotation TOD -> TEME about z by almost Delta psi cos eps (the equation of the equinoxes
   * of 1982) */
  double worst_ortho = 0.0;
  double worst_z = 0.0;
  double worst_eqeq = 0.0;
  for (int k = 0; k < 200; ++k) {
    const myb_time tt = {(int64_t)(k - 100) * 31557600, 0.123};
    myb_rotation teme;
    myb_rotation tod;
    myb_rotation mod;
    myb_equinox_rotation(e, MYB_EQUINOX_TEME, tt, &teme);
    myb_equinox_rotation(e, MYB_EQUINOX_TOD, tt, &tod);
    myb_equinox_rotation(e, MYB_EQUINOX_MOD, tt, &mod);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        double s = 0.0;
        for (int m = 0; m < 3; ++m) {
          s += teme.m[3 * i + m] * teme.m[3 * j + m];
        }
        worst_ortho = fmax(worst_ortho, fabs(s - (i == j ? 1.0 : 0.0)));
      }
      worst_z = fmax(worst_z, fabs(teme.m[6 + i] - tod.m[6 + i]));
    }
    /* angle from the true equinox (TOD x) to the TEME x, about z */
    double cx = 0.0;
    double sx = 0.0;
    for (int m = 0; m < 3; ++m) {
      cx += teme.m[m] * tod.m[m];
      sx += teme.m[m] * tod.m[3 + m];
    }
    const double angle = atan2(sx, cx);
    double nut[2];
    myb_equinox_nutation(e, tt, nut, NULL);
    const double eps = myb_equinox_obliquity(e, tt, NULL);
    /* exact for this definition: atan(tan(Delta psi) cos(eps + Delta eps)) */
    const double expected = atan(tan(nut[0]) * cos(eps + nut[1]));
    worst_eqeq = fmax(worst_eqeq, fabs(angle - expected));
  }
  printf("TEME: orthonormal %.2e, pole %.2e, angle from the true equinox %.2e rad\n",
         worst_ortho, worst_z, worst_eqeq);
  MT_CHECK(worst_ortho < 1e-15);
  MT_CHECK(worst_z == 0.0);
  MT_CHECK(worst_eqeq < 1e-15);
}

static void test_derivatives(const myb_equinox* e) {
  double worst = 0.0;
  double worst_gmst = 0.0;
  for (int k = 0; k < 100; ++k) {
    const myb_time t = {(int64_t)(k - 50) * 63113851, 0.5};
    const double h = 1000.0;
    for (int frame = 0; frame < 3; ++frame) {
      myb_rotation r;
      myb_rotation rp;
      myb_rotation rm;
      myb_equinox_rotation(e, (myb_equinox_frame)frame, t, &r);
      myb_equinox_rotation(e, (myb_equinox_frame)frame, myb_time_add(t, h), &rp);
      myb_equinox_rotation(e, (myb_equinox_frame)frame, myb_time_add(t, -h), &rm);
      for (int i = 0; i < 9; ++i) {
        worst = fmax(worst, fabs(r.dm[i] - (rp.m[i] - rm.m[i]) / (2.0 * h)));
      }
    }
    /* the rate of GMST; the angle sums seconds up to 1e7 (rounding 2e-9 s), so a wide step */
    double rate;
    myb_equinox_gmst(e, t, &rate);
    const double gp = myb_equinox_gmst(e, myb_time_add(t, 100.0), NULL);
    const double gm = myb_equinox_gmst(e, myb_time_add(t, -100.0), NULL);
    double dg = gp - gm;
    if (dg < 0.0) {
      dg += 6.283185307179586;
    }
    worst_gmst = fmax(worst_gmst, fabs(rate - dg / 200.0));
    myb_rotation g;
    myb_equinox_teme_pef(e, t, &g);
    MT_CHECK(fabs(g.dm[1] - rate * g.m[0]) <= 1e-20 && fabs(g.dm[0] + rate * g.m[1]) <= 1e-20);
  }
  printf("derivatives against central differences: frames %.2e /s, GMST %.2e rad/s\n", worst,
         worst_gmst);
  /* the 13.7-day nutation term (0.23 arcsec) leaves A w^3 h^2 / 6 = 3e-17 /s with h = 1000 s */
  MT_CHECK(worst < 2e-16);
  MT_CHECK(worst_gmst < 1e-14);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  myb_equinox* e = load();
  if (e != NULL) {
    test_against_erfa(e);
    test_teme(e);
    test_derivatives(e);
  }
  myb_equinox_destroy(e);
  return mt_finish("test_erfa_equinox");
}
