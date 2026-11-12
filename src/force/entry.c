/* SPDX-License-Identifier: BSD-3-Clause */
/* Atmospheric entry (myb/surface.h): the hypersonic continuum of the modified Newtonian theory
 * on plates and on a mesh, the heating of the stagnation point and the ablation it drives. The
 * term that bridges the continuum to the free molecules is with the other surface terms, in
 * force/surface.c. */
#include <math.h>

#include "core/base.h"
#include "force/surface_internal.h"
#include "myb/cmath.h"
#include "myb/surface.h"
#include "num/vec3.h"

static const double gas_constant = 8.314462618; /* J mol^-1 K^-1, exact since 2019 */

/* ---- Continuum flow -------------------------------------------------------------------------- */

int myb_internal_air_valid(const myb_air* air) {
  return air->temperature > 0.0 && air->molar_mass > 0.0 && air->heat_ratio > 1.0 &&
         isfinite(air->temperature) && isfinite(air->molar_mass) && isfinite(air->heat_ratio);
}

double myb_air_heat_ratio(double monatomic, double diatomic) {
  if (!(monatomic >= 0.0) || !(diatomic >= 0.0) || !isfinite(monatomic) || !isfinite(diatomic) ||
      !(monatomic + diatomic > 0.0)) {
    return 0.0;
  }
  /* c_v = 3/2 R per monatomic molecule, 5/2 R per diatomic one; c_p = c_v + R */
  const double total = monatomic + diatomic;
  const double one = monatomic / total;
  const double two = diatomic / total;
  return (5.0 * one + 7.0 * two) / (3.0 * one + 5.0 * two);
}

myb_status myb_newtonian_cp_max(double mach, double heat_ratio, double* cp_max) {
  MYB_CHECK(cp_max != NULL);
  const double g = heat_ratio;
  if (!(mach >= 0.0) || !isfinite(mach) || !(g > 1.0) || !isfinite(g)) {
    return MYB_E_DOMAIN;
  }
  const double power = g / (g - 1.0);
  if (mach >= 1.0) {
    /* Rayleigh's pitot formula divided through by M^2, so that nothing grows with M:
     *   p_02 / (p M^2) = [(g + 1)^2 / (4 g - 2 (g - 1) / M^2)]^(g / (g - 1))
     *                    (2 g - (g - 1) / M^2) / (g + 1),
     * C_p,max = (2 / g) (p_02 / (p M^2) - 1 / M^2) */
    const double inverse = 1.0 / (mach * mach);
    const double ratio = (g + 1.0) * (g + 1.0) / (4.0 * g - 2.0 * (g - 1.0) * inverse);
    const double behind =
        myb_exp(power * myb_log(ratio)) * (2.0 * g - (g - 1.0) * inverse) / (g + 1.0);
    *cp_max = 2.0 / g * (behind - inverse);
    return MYB_OK;
  }
  const double m2 = mach * mach;
  if (mach < 1.0e-8) {
    /* 1 + M^2 / 4 + (2 - g) M^4 / 24 + ...: the second term is the last one a double sees */
    *cp_max = 1.0 + 0.25 * m2;
    return MYB_OK;
  }
  /* no shock: the isentropic stagnation pressure, p_0 / p - 1 without the cancellation */
  *cp_max = 2.0 / (g * m2) * myb_expm1(power * myb_log1p(0.5 * (g - 1.0) * m2));
  return MYB_OK;
}

/* The speed of the object relative to the air and C_p,max at its Mach number. */
static myb_status flow_of(double density, const myb_air* air, const double v[3], double* speed,
                          double* cp_max) {
  *speed = myb_vec3_norm(v);
  if (!myb_internal_air_valid(air) || !(density >= 0.0) || !isfinite(density) ||
      !(*speed > 0.0) || !isfinite(*speed)) {
    return MYB_E_DOMAIN;
  }
  const double sound =
      sqrt(air->heat_ratio * gas_constant * air->temperature / air->molar_mass);
  return myb_newtonian_cp_max(*speed / sound, air->heat_ratio, cp_max);
}

myb_status myb_force_newtonian_plates(double density, const myb_air* air, double mass,
                                      const myb_plate* plates, size_t count, const double v[3],
                                      double a[3]) {
  MYB_CHECK(air != NULL && (plates != NULL || count == 0) && v != NULL && a != NULL);
  if (!(mass > 0.0) || !isfinite(mass)) {
    return MYB_E_DOMAIN;
  }
  double speed = 0.0;
  double cp_max = 0.0;
  const myb_status status = flow_of(density, air, v, &speed, &cp_max);
  if (status != MYB_OK) {
    return status;
  }
  /* the direction of the motion: the flow comes from it */
  const double ahead[3] = {v[0] / speed, v[1] / speed, v[2] / speed};
  double sum[3] = {0.0, 0.0, 0.0};
  for (size_t k = 0; k < count; ++k) {
    const myb_plate* p = &plates[k];
    const double n2 = myb_vec3_dot(p->normal, p->normal);
    if (!(p->area >= 0.0) || !isfinite(p->area) || !(fabs(n2 - 1.0) <= 1e-12)) {
      return MYB_E_DOMAIN;
    }
    const double sine = myb_vec3_dot(p->normal, ahead);
    if (!(sine > 0.0)) {
      continue; /* leeward or edge-on */
    }
    const double push = p->area * sine * sine;
    for (int i = 0; i < 3; ++i) {
      sum[i] -= push * p->normal[i];
    }
  }
  const double scale = 0.5 * density * speed * speed * cp_max / mass;
  for (int i = 0; i < 3; ++i) {
    a[i] += scale * sum[i];
  }
  return MYB_OK;
}

myb_status myb_force_newtonian_mesh(const myb_polyhedron* shape, double density,
                                    const myb_air* air, const double v[3], int shadows,
                                    double force[3], double torque[3]) {
  MYB_CHECK(shape != NULL && air != NULL && v != NULL && shadows >= 0);
  double speed = 0.0;
  double cp_max = 0.0;
  const myb_status status = flow_of(density, air, v, &speed, &cp_max);
  if (status != MYB_OK) {
    return status;
  }
  /* a mirror under light of pressure p takes 2 p A sin^2 theta along -n on every face turned
   * to the light: the Newtonian pressure with p = (1/4) rho V^2 C_p,max, the flow going along
   * -v as the light goes along u */
  const myb_plate mirror = {.specular = 1.0};
  const double u[3] = {-v[0] / speed, -v[1] / speed, -v[2] / speed};
  return myb_force_radiation_mesh(shape, &mirror, 0.25 * density * speed * speed * cp_max, u,
                                  shadows, force, torque);
}

/* ---- Heat and ablation ----------------------------------------------------------------------- */

myb_status myb_heating_sutton_graves(double k, double density, double nose_radius, double speed,
                                     double* q) {
  MYB_CHECK(q != NULL);
  if (!(k >= 0.0) || !(density >= 0.0) || !(nose_radius > 0.0) || !(speed >= 0.0) ||
      !isfinite(k) || !isfinite(density) || !isfinite(nose_radius) || !isfinite(speed)) {
    return MYB_E_DOMAIN;
  }
  *q = k * sqrt(density / nose_radius) * speed * speed * speed;
  return MYB_OK;
}

static int ablation_valid(const myb_ablation* a) {
  return a->k >= 0.0 && a->nose_radius > 0.0 && a->area >= 0.0 && a->heat > 0.0 &&
         isfinite(a->k) && isfinite(a->nose_radius) && isfinite(a->area) && isfinite(a->heat);
}

myb_status myb_ablation_rate(const myb_ablation* ablation, double t, const double x[3],
                             const double v[3], double* heating, double* rate) {
  MYB_CHECK(ablation != NULL && ablation->env != NULL && ablation->density != NULL && x != NULL &&
            v != NULL);
  if (!ablation_valid(ablation)) {
    return MYB_E_DOMAIN;
  }
  myb_time when;
  double fixed[3];
  double density = 0.0;
  double w[3];
  myb_status status = myb_internal_gas_at(ablation->env, t, ablation->body, ablation->density,
                                          ablation->density_context, ablation->wind,
                                          ablation->wind_context, x, v, &when, fixed, &density, w);
  if (status != MYB_OK) {
    return status;
  }
  double q = 0.0;
  status = myb_heating_sutton_graves(ablation->k, density, ablation->nose_radius,
                                     myb_vec3_norm(w), &q);
  if (status != MYB_OK) {
    return status;
  }
  if (heating != NULL) {
    *heating = q;
  }
  if (rate != NULL) {
    *rate = -q * ablation->area / ablation->heat;
  }
  return MYB_OK;
}

myb_status myb_ablation_rhs(void* context, double t, const double* y, double* dydt) {
  MYB_CHECK(context != NULL && y != NULL && dydt != NULL);
  myb_ablation* ablation = context;
  if (!(y[6] > 0.0) || !isfinite(y[6])) {
    return MYB_E_DOMAIN; /* the object is consumed */
  }
  /* the mass that the terms of the field read through myb_ablation_mass */
  ablation->mass = y[6];
  double a[3] = {0.0, 0.0, 0.0};
  myb_status status = MYB_OK;
  if (ablation->field != NULL) {
    status = myb_force_model_acceleration(ablation->field, t, y, y + 3, a);
    if (status != MYB_OK) {
      return status;
    }
  }
  double rate = 0.0;
  status = myb_ablation_rate(ablation, t, y, y + 3, NULL, &rate);
  if (status != MYB_OK) {
    return status;
  }
  for (int i = 0; i < 3; ++i) {
    dydt[i] = y[3 + i];
    dydt[3 + i] = a[i];
  }
  dydt[6] = rate;
  return MYB_OK;
}

double myb_ablation_mass(void* context, double t) {
  MYB_CHECK(context != NULL);
  (void)t;
  const myb_ablation* ablation = context;
  return ablation->mass;
}
