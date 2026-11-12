/* SPDX-License-Identifier: BSD-3-Clause */
/* Charging of a sphere in a plasma (myb/force.h): orbit-motion-limited currents, the thin
 * sheath, photoelectrons and secondary electrons. */
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

/* ---- Secondary electrons and the thin sheath ------------------------------------------------ */

enum {
  YIELD_PANELS = 24,  /* of the integral over the incident energies */
  SCAN_STEPS = 3200   /* of the search for the equilibria: +-200 scales in eighths */
};

static const double gl_node[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                  0.9602898564975363};
static const double gl_weight[4] = {0.3626837833783620, 0.3137066458778873,
                                    0.2223810344533745, 0.1012285362903763};

double myb_secondary_yield(const myb_secondary_emission* emission, double energy) {
  MYB_CHECK(emission != NULL);
  if (!(energy > 0.0)) {
    return 0.0;
  }
  const double x = energy / emission->max_energy;
  return 7.38905609893065 * emission->max_yield * x * myb_exp(-2.0 * sqrt(x));
}

static int secondary_ok(const myb_secondary_emission* s) {
  return s == NULL || (s->max_yield >= 0.0 && isfinite(s->max_yield) && s->max_energy > 0.0 &&
                       isfinite(s->max_energy) && s->temperature > 0.0 &&
                       isfinite(s->temperature) && s->charge > 0.0 && isfinite(s->charge));
}

/* The mean yield of a Maxwellian of electrons (k T) hitting with the energies E0 + w, w >= 0
 * weighted by (w + offset) e^(-w / k T) / (k T)^2 (offset: the energy the weight counts from
 * below the impact energy), integrated in w / k T on panels to 48 temperatures. */
static double mean_yield(const myb_secondary_emission* s, double temperature, double e0,
                         double offset) {
  double sum = 0.0;
  if (e0 == 0.0) {
    /* the yield goes as sqrt near zero energy: in u = sqrt(x) the integrand is smooth */
    const double du = 7.0 / (double)YIELD_PANELS;
    for (int p = 0; p < YIELD_PANELS; ++p) {
      const double centre = ((double)p + 0.5) * du;
      for (int k = 0; k < 8; ++k) {
        const double u = centre + (k < 4 ? -1.0 : 1.0) * 0.5 * du * gl_node[k % 4];
        const double x = u * u;
        sum += 0.5 * du * gl_weight[k % 4] * myb_secondary_yield(s, x * temperature) *
               (x + offset / temperature) * myb_exp(-x) * 2.0 * u;
      }
    }
    return sum;
  }
  const double width = 2.0;
  for (int p = 0; p < YIELD_PANELS; ++p) {
    const double centre = ((double)p + 0.5) * width;
    for (int k = 0; k < 8; ++k) {
      const double x = centre + (k < 4 ? -1.0 : 1.0) * 0.5 * width * gl_node[k % 4];
      const double w = x * temperature;
      sum += 0.5 * width * gl_weight[k % 4] * myb_secondary_yield(s, e0 + w) *
             (x + offset / temperature) * myb_exp(-x);
    }
  }
  return sum;
}

/* Is the species an electron (negative, light)? */
static int electron(const myb_plasma_species* s) {
  return s->charge < 0.0 && s->mass < 1e-29;
}

static double model_current(const myb_plasma_species* species, size_t count,
                            const myb_charging_model* model, double radius, double potential) {
  const double pi = 3.14159265358979323846;
  const double area = pi * radius * radius;
  const int thin = model->sheath == MYB_SHEATH_THIN;
  double sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    const myb_plasma_species* s = &species[k];
    const double vth = sqrt(2.0 * s->temperature / s->mass);
    const double chi = s->charge * potential / s->temperature;
    double normalised = oml_normalised(s->drift / vth, thin && chi < 0.0 ? 0.0 : chi);
    double current = s->charge * s->density * area * vth * normalised;
    if (thin && chi < 0.0 && s->charge > 0.0) {
      /* attracted ions: at least the flux of Bohm through the whole surface */
      double hottest = 0.0;
      for (size_t j = 0; j < count; ++j) {
        if (species[j].charge < 0.0) {
          hottest = fmax(hottest, species[j].temperature);
        }
      }
      const double bohm = s->charge * s->density * 0.6065306597126334 *
                          sqrt(hottest / s->mass) * 4.0 * area;
      current = fmax(current, bohm);
    }
    sum += current;
    if (model->secondary != NULL && electron(s) && s->density > 0.0) {
      /* the secondaries the electrons release: the electron flux at zero potential, the
       * incident energies as the sheath and the potential make them */
      const myb_secondary_emission* se = model->secondary;
      const double e_phi = -s->charge * potential; /* e phi for the electrons */
      const double random = -s->charge * s->density * area * vth * 2.0 / SQRT_PI;
      double yield = 0.0;
      if (e_phi < 0.0) {
        /* repelled: those that pass, with E - e|phi| at the surface */
        yield = myb_exp(e_phi / s->temperature) * mean_yield(se, s->temperature, 0.0, 0.0);
      } else if (thin) {
        /* the sheath accelerates the random flux: E + e phi */
        yield = mean_yield(se, s->temperature, e_phi, 0.0);
      } else {
        /* the orbits gather more of them, weighted by the energy at the surface */
        yield = mean_yield(se, s->temperature, e_phi, e_phi);
      }
      double escape = 1.0;
      if (potential > 0.0) {
        const double chis = se->charge * potential / se->temperature;
        escape = (1.0 + chis) * myb_exp(-chis);
      }
      sum += random * yield * escape;
    }
  }
  const myb_photoemission* photo = model->photo;
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

double myb_debye_length(const myb_plasma_species* species, size_t count) {
  MYB_CHECK(species != NULL || count == 0);
  double sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    sum += species[k].charge * species[k].charge * species[k].density / species[k].temperature;
  }
  return sum > 0.0 ? sqrt(8.8541878128e-12 / sum) : 0.0;
}

static myb_status model_ok(const myb_plasma_species* species, size_t count,
                           const myb_charging_model* model, double radius) {
  if (!(radius > 0.0) || !isfinite(radius) || !photo_ok(model->photo) ||
      !secondary_ok(model->secondary) ||
      (model->sheath != MYB_SHEATH_THICK && model->sheath != MYB_SHEATH_THIN)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < count; ++k) {
    if (!species_ok(&species[k])) {
      return MYB_E_DOMAIN;
    }
  }
  return MYB_OK;
}

myb_status myb_charging_current_model(const myb_plasma_species* species, size_t count,
                                      const myb_charging_model* model, double radius,
                                      double potential, double* current) {
  MYB_CHECK(current != NULL && model != NULL && (species != NULL || count == 0));
  const myb_status st = model_ok(species, count, model, radius);
  if (st != MYB_OK || !isfinite(potential)) {
    return st != MYB_OK ? st : MYB_E_DOMAIN;
  }
  *current = model_current(species, count, model, radius, potential);
  return MYB_OK;
}

myb_status myb_charging_equilibria(const myb_plasma_species* species, size_t count,
                                   const myb_charging_model* model, double radius,
                                   double potentials[MYB_CHARGING_MAX_EQUILIBRIA],
                                   size_t* found) {
  MYB_CHECK(model != NULL && potentials != NULL && found != NULL &&
            (species != NULL || count == 0));
  *found = 0;
  const myb_status st = model_ok(species, count, model, radius);
  if (st != MYB_OK) {
    return st;
  }
  double scale = 0.0;
  for (size_t k = 0; k < count; ++k) {
    scale = fmax(scale, species[k].temperature / fabs(species[k].charge));
  }
  if (model->photo != NULL && model->photo->current_density > 0.0) {
    scale = fmax(scale, model->photo->temperature / model->photo->charge);
  }
  if (model->secondary != NULL) {
    scale = fmax(scale, model->secondary->temperature / model->secondary->charge);
  }
  if (!(scale > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double step = scale / 8.0;
  double before = -200.0 * scale;
  double i_before = model_current(species, count, model, radius, before);
  for (int k = 1; k <= SCAN_STEPS && *found < MYB_CHARGING_MAX_EQUILIBRIA; ++k) {
    const double after = -200.0 * scale + (double)k * step;
    const double i_after = model_current(species, count, model, radius, after);
    if (!isfinite(i_after)) {
      return MYB_E_DOMAIN;
    }
    /* a stable balance: positive below, negative above */
    if (i_before > 0.0 && i_after <= 0.0) {
      double lo = before;
      double hi = after;
      for (int h = 0; h < 2200; ++h) { /* bound: halvings to adjacent doubles */
        const double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) {
          break;
        }
        if (model_current(species, count, model, radius, mid) > 0.0) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      potentials[(*found)++] = lo; /* the last potential with a positive current */
    }
    before = after;
    i_before = i_after;
  }
  return *found > 0 ? MYB_OK : MYB_E_DOMAIN;
}
