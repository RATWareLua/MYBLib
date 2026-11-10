/* SPDX-License-Identifier: BSD-3-Clause */
/* Orbit-motion-limited charging of a sphere in a plasma (myb/force.h). */
#include <math.h>

#include "core/check.h"
#include "myb/cmath.h"
#include "myb/force.h"
#include "num/erf.h"

#define SQRT_PI 1.7724538509055160273

/* Current of one species divided by q n pi a^2 v_th, for the Mach number M = u / v_th and
 * chi = q phi / (k T). */
static double oml_normalised(double mach, double chi) {
  if (chi <= 0.0) {
    /* attracted: <v> - chi v_th^2 <1/v> for the drifting Maxwellian, <1/v> = erf(M) / u */
    if (mach == 0.0) {
      return 2.0 / SQRT_PI * (1.0 - chi);
    }
    const double erf_m = myb_internal_erf(mach);
    return (mach + 0.5 / mach) * erf_m + myb_exp(-mach * mach) / SQRT_PI - chi * erf_m / mach;
  }
  const double x0 = sqrt(chi);
  if (mach < 1e-3) {
    /* the series in M (the closed form below cancels as M -> 0):
     * (2 / sqrt(pi)) e^(-chi) [1 + M^2 (2 chi / 3 + 1/3) + O(M^4)] */
    return 2.0 / SQRT_PI * myb_exp(-chi) * (1.0 + mach * mach * (2.0 * chi / 3.0 + 1.0 / 3.0));
  }
  /* (1 / (sqrt(pi) M)) [J(x0 - M, +M) - J(x0 + M, -M)],
   * J(y0, s) = int_y0^inf ((y + s)^2 - chi) e^(-y^2) dy
   *          = (y0 / 2 + s) e^(-y0^2) + (sqrt(pi) / 4) erfc(y0) (1 + 2 (s^2 - chi)) */
  const double ya = x0 - mach;
  const double yb = x0 + mach;
  const double s2 = mach * mach - chi;
  const double ja = (0.5 * ya + mach) * myb_exp(-ya * ya) +
                    SQRT_PI / 4.0 * myb_internal_erfc(ya) * (1.0 + 2.0 * s2);
  const double jb = (0.5 * yb - mach) * myb_exp(-yb * yb) +
                    SQRT_PI / 4.0 * myb_internal_erfc(yb) * (1.0 + 2.0 * s2);
  return (ja - jb) / (SQRT_PI * mach);
}

static int species_ok(const myb_plasma_species* s) {
  return s->charge != 0.0 && isfinite(s->charge) && s->mass > 0.0 && isfinite(s->mass) &&
         s->density >= 0.0 && isfinite(s->density) && s->temperature > 0.0 &&
         isfinite(s->temperature) && s->drift >= 0.0 && isfinite(s->drift);
}

static int photo_ok(const myb_photoemission* p) {
  return p == NULL || (p->current_density >= 0.0 && isfinite(p->current_density) &&
                       p->temperature > 0.0 && isfinite(p->temperature) && p->charge > 0.0 &&
                       isfinite(p->charge));
}

static double total_current(const myb_plasma_species* species, size_t count,
                            const myb_photoemission* photo, double radius, double potential) {
  const double area = 3.14159265358979323846 * radius * radius;
  double sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    const myb_plasma_species* s = &species[k];
    const double vth = sqrt(2.0 * s->temperature / s->mass);
    const double chi = s->charge * potential / s->temperature;
    sum += s->charge * s->density * area * vth * oml_normalised(s->drift / vth, chi);
  }
  if (photo != NULL && photo->current_density > 0.0) {
    double f = 1.0;
    if (potential > 0.0) {
      const double chi = photo->charge * potential / photo->temperature;
      f = (1.0 + chi) * myb_exp(-chi);
    }
    sum += area * photo->current_density * f;
  }
  return sum;
}

myb_status myb_charging_current(const myb_plasma_species* species, size_t count,
                                const myb_photoemission* photo, double radius, double potential,
                                double* current) {
  MYB_CHECK(current != NULL && (species != NULL || count == 0));
  if (!(radius > 0.0) || !isfinite(radius) || !isfinite(potential) || !photo_ok(photo)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < count; ++k) {
    if (!species_ok(&species[k])) {
      return MYB_E_DOMAIN;
    }
  }
  *current = total_current(species, count, photo, radius, potential);
  return MYB_OK;
}

myb_status myb_charging_equilibrium(const myb_plasma_species* species, size_t count,
                                    const myb_photoemission* photo, double radius,
                                    double* potential) {
  MYB_CHECK(potential != NULL && (species != NULL || count == 0));
  double i0 = 0.0;
  myb_status st = myb_charging_current(species, count, photo, radius, 0.0, &i0);
  if (st != MYB_OK) {
    return st;
  }
  if (i0 == 0.0) {
    *potential = 0.0;
    return MYB_OK;
  }
  /* the scale of the potential: the largest k T / |q| */
  double scale = 0.0;
  for (size_t k = 0; k < count; ++k) {
    scale = fmax(scale, species[k].temperature / fabs(species[k].charge));
  }
  if (photo != NULL && photo->current_density > 0.0) {
    scale = fmax(scale, photo->temperature / photo->charge);
  }
  if (!(scale > 0.0)) {
    return MYB_E_DOMAIN;
  }
  /* bracket: the current decreases with phi; step outwards from 0 */
  const double dir = i0 > 0.0 ? 1.0 : -1.0;
  double lo = 0.0;
  double hi = 0.0;
  double step = scale;
  int bracketed = 0;
  for (int k = 0; k < 200; ++k) { /* bound: 200 doublings (potentials up to 2^200 kT/q) */
    const double x = dir * step;
    const double i = total_current(species, count, photo, radius, x);
    if (!isfinite(i)) {
      return MYB_E_NOCONV;
    }
    /* a current that underflows to 0 (all particles repelled) is not a balance */
    if (i != 0.0 && (i > 0.0) != (i0 > 0.0)) {
      lo = dir > 0.0 ? step / 2.0 * (k > 0) : x;
      hi = dir > 0.0 ? x : -step / 2.0 * (k > 0);
      bracketed = 1;
      break;
    }
    step *= 2.0;
  }
  if (!bracketed) {
    return MYB_E_DOMAIN; /* currents of one sign at every potential */
  }
  /* bisection: I(lo) > 0 > I(hi) */
  for (int k = 0; k < 2200; ++k) { /* bound: halvings to adjacent doubles */
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) {
      break;
    }
    const double i = total_current(species, count, photo, radius, mid);
    if (i == 0.0) {
      *potential = mid;
      return MYB_OK;
    }
    if (i > 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  *potential = 0.5 * (lo + hi);
  return MYB_OK;
}
