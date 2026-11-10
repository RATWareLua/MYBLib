/* SPDX-License-Identifier: BSD-3-Clause */
/* Orbit-motion-limited charging (myb_charging_current, myb_charging_equilibrium):
 *   - the current of one species against the direct integration of the OML cross-section over
 *     the drifting Maxwellian (mpmath, 30 digits) for attracted and repelled particles, at rest
 *     and drifting up to Mach 8 (ram), including the series branch of small drifts;
 *   - the floating potential of a sphere in an electron-proton plasma of equal temperatures,
 *     e phi / k T = 1 - W(e sqrt(m_p / m_e)) with Lambert's W (Spitzer's -2.5);
 *   - a dust grain in the solar wind at 1 AU charged positive by photoemission (a few volts),
 *     and the current balance at the root. */
#include <math.h>
#include <stdio.h>

#include "mtest.h"
#include "myb/force.h"
#include "num/erf.h"

#define QE 1.602176634e-19
#define ME 9.1093837139e-31
#define MP 1.67262192595e-27

static void test_oml_currents(void) {
  static const double cases[10][3] = {
      {0.0, -2.0, 3.385137501286537721688},    {0.0, 1.5, 0.2517754242621735826393},
      {0.2, 2.5, 0.100101512430734264985},     {0.2, -1.0, 2.256877330473115657649},
      {8.0, -0.5, 8.125},                      {8.0, 0.3, 8.025000000000000001388},
      {8.0, 70.0, 0.2968663591530668273424},   {1e-4, 2.0, 0.1527095167223228876045},
      {3.0, 5.0, 1.564695348199868916652},     {0.05, 0.7, 0.5614569697211371630823}};
  const double kt = 10.0 * QE; /* 10 eV protons */
  const double vth = sqrt(2.0 * kt / MP);
  const double radius = 1e-6;
  const double area = 3.14159265358979323846 * radius * radius;
  double worst = 0.0;
  for (int k = 0; k < 10; ++k) {
    const myb_plasma_species s = {.charge = QE, .mass = MP, .density = 5e6, .temperature = kt,
                                  .drift = cases[k][0] * vth};
    double i = 0.0;
    MT_CHECK(myb_charging_current(&s, 1, NULL, radius, cases[k][1] * kt / QE, &i) == MYB_OK);
    const double normalised = i / (QE * 5e6 * area * vth);
    const double err = fabs(normalised / cases[k][2] - 1.0);
    worst = fmax(worst, err);
    if (err > 1e-13) {
      printf("M %.4g, chi %.4g: %.17g (mpmath %.17g)\n", cases[k][0], cases[k][1], normalised,
             cases[k][2]);
    }
  }
  printf("OML currents vs the integration over the drifting Maxwellian: %.1e\n", worst);
  MT_CHECK(worst < 1e-13);
}

static void test_floating_potential(void) {
  const double kt = 1.0 * QE;
  const myb_plasma_species plasma[2] = {
      {.charge = -QE, .mass = ME, .density = 1e10, .temperature = kt},
      {.charge = QE, .mass = MP, .density = 1e10, .temperature = kt}};
  double phi = 0.0;
  MT_CHECK(myb_charging_equilibrium(plasma, 2, NULL, 1e-6, &phi) == MYB_OK);
  const double y = phi * QE / kt;
  printf("electron-proton plasma, equal temperatures: e phi / k T = %.16f (Lambert W "
         "-2.503851165125069269)\n",
         y);
  MT_CHECK(fabs(y + 2.503851165125069269137) < 1e-13);
  /* electrons only: no balance */
  MT_CHECK(myb_charging_equilibrium(plasma, 1, NULL, 1e-6, &phi) == MYB_E_DOMAIN);
}

static void test_solar_wind_grain(void) {
  /* 1 AU: n = 5 cm^-3, T_e = 12 eV, T_p = 8 eV, 400 km/s; photoemission 2.5e14 electrons per
   * m^2 and s from the cross-section (yield of a dielectric), photoelectrons of 2 eV */
  const myb_plasma_species wind[2] = {
      {.charge = -QE, .mass = ME, .density = 5e6, .temperature = 12.0 * QE, .drift = 4e5},
      {.charge = QE, .mass = MP, .density = 5e6, .temperature = 8.0 * QE, .drift = 4e5}};
  const myb_photoemission photo = {.current_density = 2.5e14 * QE, .temperature = 2.0 * QE,
                                   .charge = QE};
  double phi = 0.0;
  MT_CHECK(myb_charging_equilibrium(wind, 2, &photo, 1e-6, &phi) == MYB_OK);
  double i = 0.0;
  double ie = 0.0;
  MT_CHECK(myb_charging_current(wind, 2, &photo, 1e-6, phi, &i) == MYB_OK);
  MT_CHECK(myb_charging_current(wind, 1, NULL, 1e-6, phi, &ie) == MYB_OK);
  printf("dust grain in the solar wind at 1 AU: +%.4f V, residual current %.1e of %.1e A\n", phi,
         i, fabs(ie));
  /* photoelectrons held back by a few k T_p / e balance the collected electrons */
  MT_CHECK(phi > 2.0 * 2.0 && phi < 6.0 * 2.0);
  MT_CHECK(fabs(i) < 1e-12 * fabs(ie));
  /* without sunlight the grain charges negative, a few k T_e / e */
  MT_CHECK(myb_charging_equilibrium(wind, 2, NULL, 1e-6, &phi) == MYB_OK);
  printf("in the shadow: %.4f V\n", phi);
  MT_CHECK(phi < -12.0 && phi > -60.0);
  /* domain */
  const myb_plasma_species bad = {.charge = QE, .mass = 0.0, .density = 1.0, .temperature = 1.0};
  MT_CHECK(myb_charging_current(&bad, 1, NULL, 1e-6, 0.0, &i) == MYB_E_DOMAIN);
  MT_CHECK(myb_charging_current(wind, 2, NULL, 0.0, 0.0, &i) == MYB_E_DOMAIN);
}

/* the internal error function (fdlibm algorithm) against mpmath: within 2 ulp */
static void test_erf(void) {
  static const double v[11][3] = {
      {1e-10, 1.1283791670955125739e-10, 0.99999999988716208329},
      {0.1, 0.1124629160182848922, 0.8875370839817151078},
      {0.5, 0.52049987781304653768, 0.47950012218695346232},
      {0.84, 0.76514271145499453466, 0.23485728854500546534},
      {1.0, 0.84270079294971486934, 0.15729920705028513066},
      {1.2, 0.91031397822963538024, 0.089686021770364619762},
      {2.0, 0.99532226501895273416, 0.0046777349810472658379},
      {2.9, 0.99995890212190054116, 0.000041097878099458835684},
      {4.0, 0.99999998458274209972, 1.5417257900280018852e-8},
      {5.5, 0.99999999999999264215, 7.3578479179743980631e-15},
      {-1.3, -0.9340079449406524366, 1.9340079449406524366}};
  double worst = 0.0;
  for (int k = 0; k < 11; ++k) {
    worst = fmax(worst, fabs(myb_internal_erf(v[k][0]) / v[k][1] - 1.0));
    worst = fmax(worst, fabs(myb_internal_erfc(v[k][0]) / v[k][2] - 1.0));
  }
  worst = fmax(worst, fabs(myb_internal_erfc(10.0) / 2.088487583762544757e-45 - 1.0));
  printf("erf, erfc vs mpmath: %.1e\n", worst);
  MT_CHECK(worst < 5e-16);
  MT_CHECK(myb_internal_erf(0.0) == 0.0 && myb_internal_erfc(0.0) == 1.0);
  MT_CHECK(myb_internal_erf(-30.0) == -1.0 && myb_internal_erfc(-30.0) == 2.0);
  MT_CHECK(myb_internal_erfc(30.0) == 0.0);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_erf();
  test_oml_currents();
  test_floating_potential();
  test_solar_wind_grain();
  return mt_finish("test_charging");
}
