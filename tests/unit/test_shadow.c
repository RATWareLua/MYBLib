/* SPDX-License-Identifier: BSD-3-Clause */
/* Shadows beyond one airless occulter (myb/force.h):
 *  - several occulters at once: one in front is myb_shadow_fraction to the bit; two apart on
 *    the disk add up; two that overlap against a count of the points of the disk they
 *    leave, with and without limb darkening;
 *  - the atmosphere of the occulter: with next to no air the integration over the true disk
 *    gives the airless fraction; the bending of a ray that grazes the Earth; without
 *    extinction the light only moves — the flux missing over a plane behind the body is that
 *    of the rays that meet its surface; extinction only takes away; an atmosphere dense
 *    enough to trap light;
 *  - the gradient of the radiation of a planet against differences of the acceleration;
 *  - the terms of a model (myb_term_radiation, myb_term_surface_radiation): with two bodies
 *    on the disk of the Sun the light is that of the union of what they cover, with one it
 *    is the product as before, to the bit; with an atmosphere given for an occulter, that of
 *    myb_shadow_fraction_refracted. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double pi = 3.141592653589793;
static const double au = 1.495978707e11;
static const double sun_radius = 6.957e8;

static double norm3(const double v[3]) {
  return hypot(hypot(v[0], v[1]), v[2]);
}

/* The fraction of the disk of the source left by two occulters, by counting: circles about
 * the centre of the source, each cut into equal arcs. */
static double counted(const double source[3], double rs, const double o1[3], double r1,
                      const double o2[3], double r2, const double limb[2]) {
  const double ds = norm3(source);
  const double a = asin(rs / ds);
  const double s[3] = {source[0] / ds, source[1] / ds, source[2] / ds};
  const double helper[3] = {0.0, 0.0, 1.0};
  double e1[3] = {helper[1] * s[2] - helper[2] * s[1], helper[2] * s[0] - helper[0] * s[2],
                  helper[0] * s[1] - helper[1] * s[0]};
  const double n1 = norm3(e1);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= n1;
  }
  const double e2[3] = {s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2],
                        s[0] * e1[1] - s[1] * e1[0]};
  const double d1 = norm3(o1);
  const double d2 = norm3(o2);
  const double cb1 = cos(asin(r1 / d1));
  const double cb2 = cos(asin(r2 / d2));
  const int rings = 3000;
  const int arcs = 1500;
  double total = 0.0;
  double seen = 0.0;
  for (int i = 0; i < rings; ++i) {
    const double rho = a * (i + 0.5) / rings;
    const double mu = sqrt(fmax(0.0, 1.0 - sin(rho) * sin(rho) / (sin(a) * sin(a))));
    const double light = limb != NULL
                             ? 1.0 - limb[0] * (1.0 - mu) - limb[1] * (1.0 - mu) * (1.0 - mu)
                             : 1.0;
    const double weight = light * sin(rho);
    int open = 0;
    for (int j = 0; j < arcs; ++j) {
      const double phi = 2.0 * pi * (j + 0.5) / arcs;
      double p[3];
      for (int k = 0; k < 3; ++k) {
        p[k] = cos(rho) * s[k] + sin(rho) * (cos(phi) * e1[k] + sin(phi) * e2[k]);
      }
      const double c1 = (p[0] * o1[0] + p[1] * o1[1] + p[2] * o1[2]) / d1;
      const double c2 = (p[0] * o2[0] + p[1] * o2[1] + p[2] * o2[2]) / d2;
      open += !(c1 > cb1 || c2 > cb2);
    }
    total += weight;
    seen += weight * open / arcs;
  }
  return seen / total;
}

static void test_many(void) {
  /* seen from the origin: the Sun along x; occulters a little off that direction */
  const double source[3] = {au, 0.0, 0.0};
  const double a = asin(sun_radius / au);
  const double far[3] = {4.0e8, 4.0e8 * 0.02, 0.0}; /* 1.1 deg off: not in front */
  const double radii[2] = {1.7374e6, 1.7374e6};
  /* one in front */
  double pair[6] = {3.8e8, 3.8e8 * 0.004, 0.0, far[0], far[1], far[2]};
  const double single = myb_shadow_fraction(source, sun_radius, pair, radii[0], NULL);
  MT_CHECK(single > 0.05 && single < 0.95);
  MT_CHECK(myb_shadow_fraction_many(source, sun_radius, pair, radii, 2, NULL) == single);
  MT_CHECK(myb_shadow_fraction_many(source, sun_radius, pair, radii, 1, NULL) == single);
  MT_CHECK(myb_shadow_fraction_many(source, sun_radius, pair, radii, 0, NULL) == 1.0);
  /* the same occulter twice covers what it covers once */
  memcpy(pair + 3, pair, 3 * sizeof(double));
  MT_CHECK(fabs(myb_shadow_fraction_many(source, sun_radius, pair, radii, 2, NULL) - single) <
           1e-13);
  /* two small ones on opposite sides of the disk: what each takes adds up */
  const double small[2] = {5.0e5, 6.0e5};
  double apart[6] = {3.8e8, 3.8e8 * 0.6 * a, 0.0, 3.8e8, -3.8e8 * 0.5 * a, 3.8e8 * 0.2 * a};
  const double limb[2] = {0.4, 0.25};
  for (int dark = 0; dark < 2; ++dark) {
    const double* l = dark ? limb : NULL;
    const double f1 = myb_shadow_fraction(source, sun_radius, apart, small[0], l);
    const double f2 = myb_shadow_fraction(source, sun_radius, apart + 3, small[1], l);
    const double both = myb_shadow_fraction_many(source, sun_radius, apart, small, 2, l);
    MT_CHECK(f1 < 1.0 && f2 < 1.0);
    if (fabs(both - (f1 + f2 - 1.0)) >= 1e-11) {
      printf("two occulters apart: %.3e from the sum\n", both - (f1 + f2 - 1.0));
    }
    /* the two integrations are cut into different pieces */
    MT_CHECK(fabs(both - (f1 + f2 - 1.0)) < 1e-11);
    MT_CHECK(fabs(both - f1 * f2) > 1e-4); /* the product is not it */
  }
  /* two that overlap on the disk */
  const double big[2] = {1.2e6, 1.5e6};
  double over[6] = {3.8e8, 3.8e8 * 0.35 * a, 0.0, 3.8e8, 3.8e8 * 0.1 * a, 3.8e8 * 0.45 * a};
  double worst = 0.0;
  double product_error = 0.0;
  for (int dark = 0; dark < 2; ++dark) {
    const double* l = dark ? limb : NULL;
    const double both = myb_shadow_fraction_many(source, sun_radius, over, big, 2, l);
    const double count = counted(source, sun_radius, over, big[0], over + 3, big[1], l);
    const double f1 = myb_shadow_fraction(source, sun_radius, over, big[0], l);
    const double f2 = myb_shadow_fraction(source, sun_radius, over + 3, big[1], l);
    worst = fmax(worst, fabs(both - count));
    product_error = fmax(product_error, fabs(f1 * f2 - count));
    MT_CHECK(both > 0.0 && both < f1 && both < f2);
    MT_CHECK(both > f1 + f2 - 1.0); /* the overlap is counted once */
  }
  printf("two occulters overlapping on the disk: the union of arcs within %.1e of a count of "
         "points; the product of the two fractions is off by %.1e\n",
         worst, product_error);
  MT_CHECK(worst < 3e-4);
  MT_CHECK(product_error > 20.0 * worst);
  /* three, one of them over the whole disk */
  double three[9];
  memcpy(three, over, sizeof over);
  three[6] = 2.0e7;
  three[7] = 0.0;
  three[8] = 0.0;
  const double radii3[3] = {big[0], big[1], 6.4e6};
  MT_CHECK(myb_shadow_fraction_many(source, sun_radius, three, radii3, 3, NULL) == 0.0);
}

static void test_refraction(void) {
  const double radius = 6378137.0;
  const myb_refraction air = {.refractivity = 2.77e-4, .scale_height = 8430.0};
  const double grazing = myb_refraction_bending(&air, radius, 0.0);
  printf("a ray that grazes the Earth is bent by %.1f' (twice the horizontal refraction, "
         "2 x 34.5')\n",
         grazing * 180.0 * 60.0 / pi);
  MT_CHECK(grazing * 180.0 * 60.0 / pi > 60.0 && grazing * 180.0 * 60.0 / pi < 72.0);
  MT_CHECK(fabs(myb_refraction_bending(&air, radius, 8430.0) / grazing - exp(-1.0)) < 1e-3);

  /* next to no air: the integration over the true disk against the airless fraction */
  const myb_refraction thin = {.refractivity = 1.0e-13, .scale_height = 8430.0};
  const double limb[2] = {0.4, 0.25};
  double worst = 0.0;
  for (int dark = 0; dark < 2; ++dark) {
    for (int k = 0; k < 12; ++k) {
      /* across the penumbra, 12 270 km from the centre of the Earth */
      const double off = 6.30e6 + 2.0e4 * k;
      const double along = sqrt(1.227e7 * 1.227e7 - off * off);
      const double to_earth[3] = {along, -off, 0.0};
      const double to_sun[3] = {au + along, -off, 0.0};
      const double plain =
          myb_shadow_fraction(to_sun, sun_radius, to_earth, radius, dark ? limb : NULL);
      const double through =
          myb_shadow_fraction_refracted(to_sun, sun_radius, to_earth, radius, &thin,
                                        dark ? limb : NULL);
      worst = fmax(worst, fabs(plain - through));
    }
  }
  printf("with next to no air the fraction is the airless one within %.1e\n", worst);
  MT_CHECK(worst < 1e-7);
  const double to_earth0[3] = {1.0e7, -6.4e6, 0.0};
  const double to_sun0[3] = {au + 1.0e7, -6.4e6, 0.0};
  MT_CHECK(myb_shadow_fraction_refracted(to_sun0, sun_radius, to_earth0, radius, NULL, NULL) ==
           myb_shadow_fraction(to_sun0, sun_radius, to_earth0, radius, NULL));

  /* the flux missing over a plane 40 000 km behind the Earth: without extinction, that of
   * the rays that meet the surface, pi p_min^2. The Sun is put far away so that its light
   * is parallel, as the model takes it. */
  const double far = 1.0e15;
  const double far_radius = far * sun_radius / au;
  const double behind = 4.0e7;
  const myb_refraction hazy = {.refractivity = 2.77e-4, .scale_height = 8430.0,
                               .optical_depth = 0.3};
  double missing[3] = {0.0, 0.0, 0.0}; /* airless, clear air, hazy air */
  double inside[3] = {1.0, 1.0, 1.0}; /* 6000 km from the axis: within the airless umbra */
  const int steps = 3500;
  const double edge = 7.0e6;
  int ordered = 1;
  for (int k = 0; k <= steps; ++k) {
    const double rho = edge * k / steps;
    const double to_earth[3] = {behind, -rho, 0.0};
    const double to_sun[3] = {far + behind, -rho, 0.0};
    const double f[3] = {
        myb_shadow_fraction(to_sun, far_radius, to_earth, radius, NULL),
        myb_shadow_fraction_refracted(to_sun, far_radius, to_earth, radius, &air, NULL),
        myb_shadow_fraction_refracted(to_sun, far_radius, to_earth, radius, &hazy, NULL)};
    const double weight = (k == 0 || k == steps ? 0.5 : 1.0) * 2.0 * pi * rho * edge / steps;
    for (int i = 0; i < 3; ++i) {
      missing[i] += weight * (1.0 - f[i]);
      if (k == 3000) {
        inside[i] = f[i];
      }
    }
    ordered = ordered && f[2] <= f[1] + 1e-12 && f[1] >= 0.0 && f[1] <= 1.0;
  }
  const double scale = (far + behind) / far;
  const double solid = pi * radius * radius * scale * scale;
  const double p_min = radius * (1.0 + air.refractivity);
  const double grazed = pi * p_min * p_min * scale * scale;
  printf("the flux missing behind the Earth, in its cross-section: airless %.6f, clear air "
         "%.6f (the rays that meet the ground: %.6f), hazy air %.4f; the light 6000 km from "
         "the axis, inside the airless umbra: %.1e, %.1e, %.1e\n",
         missing[0] / solid, missing[1] / solid, grazed / solid, missing[2] / solid, inside[0],
         inside[1], inside[2]);
  MT_CHECK(fabs(missing[0] / solid - 1.0) < 2e-4);
  MT_CHECK(fabs(missing[1] / grazed - 1.0) < 1e-3);
  MT_CHECK(missing[2] > 1.005 * missing[1]);
  MT_CHECK(ordered);
  /* the airless shadow is black inside; air lets light in, haze takes most of it back */
  MT_CHECK(inside[0] == 0.0 && inside[1] > 1e-4 && inside[2] < 0.5 * inside[1]);

  /* an atmosphere that traps light below a level over the surface: the limb is that level */
  const myb_refraction dense = {.refractivity = 1.4e-2, .scale_height = 15900.0};
  const double venus = 6.0518e6;
  const double level = venus + 60000.0; /* somewhere above the trapping level */
  int sane = 1;
  double previous = -1.0;
  for (int k = 0; k <= 60; ++k) {
    const double rho = 5.5e6 + 3.0e4 * k;
    const double to_body[3] = {3.0e7, -rho, 0.0};
    const double to_sun[3] = {far + 3.0e7, -rho, 0.0};
    const double f = myb_shadow_fraction_refracted(to_sun, far_radius, to_body, venus, &dense,
                                                   NULL);
    sane = sane && f >= 0.0 && f <= 1.0 && f >= previous - 1e-9;
    previous = f;
  }
  MT_CHECK(sane && previous > 0.99);
  (void)level;
  /* inside the body: dark; the body behind the Sun or behind the object: nothing */
  const double within[3] = {1.0e6, 0.0, 0.0};
  MT_CHECK(myb_shadow_fraction_refracted(to_sun0, sun_radius, within, radius, &air, NULL) ==
           0.0);
  const double beyond[3] = {2.0 * au, 0.0, 0.0};
  MT_CHECK(myb_shadow_fraction_refracted(to_sun0, sun_radius, beyond, radius, &air, NULL) ==
           1.0);
  const double backward[3] = {-1.0e7, 0.0, 0.0};
  MT_CHECK(myb_shadow_fraction_refracted(to_sun0, sun_radius, backward, radius, &air, NULL) ==
           1.0);
}

static void test_planet_gradient(void) {
  const myb_planet_radiation knocke = {.radius = 6378136.6,
                                       .albedo = {0.34, 0.0, 0.10, 0.0, 0.29},
                                       .emissivity = {0.68, 0.0, -0.07, 0.0, -0.18},
                                       .omega = 2.0 * pi / (365.25 * 86400.0)};
  const double r[3] = {7.0e6, -3.0e6, 9.5e6};
  const double sun[3] = {0.8 * au, 0.55 * au, 0.2 * au};
  const double k = 1.13 * 0.28 / 407.0;
  double gradient[9];
  MT_CHECK(myb_force_planet_radiation_gradient(&knocke, r, sun, 1361.0, 0.7, k, 50, gradient) ==
           MYB_OK);
  double worst = 0.0;
  double largest = 0.0;
  for (int j = 0; j < 3; ++j) {
    const double h = 2.0e4;
    double q[3];
    double plus[3] = {0.0, 0.0, 0.0};
    double minus[3] = {0.0, 0.0, 0.0};
    memcpy(q, r, sizeof q);
    q[j] = r[j] + h;
    MT_CHECK(myb_force_planet_radiation(&knocke, q, sun, 1361.0, 0.7, k, 50, plus) == MYB_OK);
    q[j] = r[j] - h;
    MT_CHECK(myb_force_planet_radiation(&knocke, q, sun, 1361.0, 0.7, k, 50, minus) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst = fmax(worst, fabs((plus[i] - minus[i]) / (2.0 * h) - gradient[3 * i + j]));
      largest = fmax(largest, fabs(gradient[3 * i + j]));
    }
  }
  printf("the gradient of the radiation of a planet against differences of the acceleration: "
         "%.1e of its largest element\n",
         worst / largest);
  MT_CHECK(worst < 5e-3 * largest);
  /* with fewer rings it is the gradient of the same field, a little coarser */
  double coarse[9];
  MT_CHECK(myb_force_planet_radiation_gradient(&knocke, r, sun, 1361.0, 0.7, k, 10, coarse) ==
           MYB_OK);
  for (int i = 0; i < 9; ++i) {
    MT_CHECK(fabs(coarse[i] - gradient[i]) < 2e-2 * largest);
  }
  const double below[3] = {1.0e6, 0.0, 0.0};
  MT_CHECK(myb_force_planet_radiation_gradient(&knocke, below, sun, 1361.0, 0.7, k, 10,
                                               coarse) == MYB_E_DOMAIN);
}

static myb_status fixed_state(void* context, myb_time t, double y[6]) {
  (void)t;
  memcpy(y, context, 6 * sizeof(double));
  return MYB_OK;
}

/* The Sun at the origin; a planet and its moon near the line to an object behind them. */
static void test_terms(void) {
  const double c_light = 299792458.0;
  const double luminosity = 3.828e26;
  const double planet_radius = 6.4e6;
  const double moon_radius = 1.7e6;
  double planet_state[6] = {au, 0.0, 0.0, 0.0, 0.0, 0.0};
  double moon_state[6] = {au + 3.0e8, 0.0, 0.0, 0.0, 0.0, 0.0};
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body sun = {.id = 10, .gm = 1.327e20, .radius = sun_radius,
                            .luminosity = luminosity, .limb = {0.4, 0.2}};
  const myb_env_body planet = {.id = 399, .gm = 3.986e14, .state = fixed_state,
                               .state_context = planet_state, .radius = planet_radius};
  const myb_env_body moon = {.id = 301, .gm = 4.9e12, .state = fixed_state,
                             .state_context = moon_state, .radius = moon_radius};
  size_t i_sun = 0;
  size_t i_planet = 0;
  size_t i_moon = 0;
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &planet, &i_planet) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &moon, &i_moon) == MYB_OK);
  const myb_object object = {.mass = 400.0, .area = 3.0, .cr = 1.4};
  myb_radiation_term light = {.env = env, .source = i_sun, .object = &object, .c = c_light,
                              .occulters = {i_planet, i_moon}, .occulter_count = 2};
  myb_radiation_term open = light;
  open.occulter_count = 0;
  const double zero[3] = {0.0, 0.0, 0.0};
  const double limb[2] = {0.4, 0.2};
  /* the object 4e8 m behind the planet and 1e8 m behind the moon: across the line the planet
   * (0.016 rad) comes on the disk of the Sun (0.0047 rad) from one side and the moon
   * (0.017 rad) from the other — apart, overlapping, covering it between them — and the
   * planet alone. The object's offsets, then the moon's from the object: */
  static const double offsets[5][4] = {{7.3e6, 0.4e6, 1.9e6, 0.0},
                                       {6.9e6, 0.0, 2.0e6, 0.3e6},
                                       {6.0e6, -0.2e6, 1.6e6, 0.0},
                                       {5.6e6, 0.1e6, 1.5e6, -0.2e6},
                                       {8.0e6, 0.0, 4.0e7, 0.0}};
  double worst_union = 0.0;
  double worst_product = 0.0;
  int partial_both = 0;
  int exact_single = 0;
  for (int k = 0; k < 5; ++k) {
    const double x[3] = {au + 4.0e8, offsets[k][0], offsets[k][1]};
    moon_state[1] = offsets[k][0] + offsets[k][2];
    moon_state[2] = offsets[k][1] + offsets[k][3];
    double full[3];
    double shaded[3];
    /* the environment keeps the states of an instant: another instant for another place
     * of the moon */
    const double t = 100.0 * (double)(k + 1);
    MT_CHECK(myb_term_radiation(&open, t, x, zero, full) == MYB_OK);
    MT_CHECK(myb_term_radiation(&light, t, x, zero, shaded) == MYB_OK);
    /* what the object sees: the bodies relative to it */
    const double to_sun[3] = {-x[0], -x[1], -x[2]};
    const double bodies[6] = {planet_state[0] - x[0], planet_state[1] - x[1],
                              planet_state[2] - x[2], moon_state[0] - x[0],
                              moon_state[1] - x[1],   moon_state[2] - x[2]};
    const double radii[2] = {planet_radius, moon_radius};
    const double together = myb_shadow_fraction_many(to_sun, sun_radius, bodies, radii, 2, limb);
    const double first = myb_shadow_fraction(to_sun, sun_radius, bodies, planet_radius, limb);
    const double second = myb_shadow_fraction(to_sun, sun_radius, bodies + 3, moon_radius, limb);
    const double seen = norm3(shaded) / norm3(full);
    if (first < 1.0 && second < 1.0) {
      ++partial_both;
      worst_union = fmax(worst_union, fabs(seen - together));
      worst_product = fmax(worst_product, fabs(first * second - together));
    } else {
      /* one body on the disk: the product, as before the union was there */
      exact_single += first < 1.0 && first > 0.0 && fabs(seen - first * second) < 1e-14;
    }
    printf("  the planet leaves %.4f, the moon %.4f, both %.4f (their product %.4f)\n", first,
           second, together, first * second);
  }
  printf("the radiation term with two bodies on the Sun's disk: %d cases within %.1e of the "
         "union of what they cover (their product is up to %.1e off); one body: %d as the "
         "product\n",
         partial_both, worst_union, worst_product, exact_single);
  MT_CHECK(partial_both == 4 && exact_single == 1);
  MT_CHECK(worst_union < 1e-14 && worst_product > 1e-3);

  /* an atmosphere for the planet: the light that the air bends into the shadow */
  const myb_refraction air = {.refractivity = 2.9e-4, .scale_height = 8.0e3,
                              .optical_depth = 0.3};
  const myb_refraction atmospheres[2] = {air, {0}};
  myb_radiation_term through_air = light;
  through_air.occulter_count = 1;
  through_air.atmospheres = atmospheres;
  myb_radiation_term airless = through_air;
  airless.atmospheres = NULL;
  moon_state[1] = 4.0e8; /* out of the way */
  double worst_air = 0.0;
  double most_gain = 0.0;
  for (int k = 0; k < 6; ++k) {
    /* from the umbra to the full light, 1.4e9 m behind the planet */
    const double x[3] = {au + 1.4e9, 0.5e6 + 2.5e6 * k, 0.0};
    double full[3];
    double bare[3];
    double bent[3];
    MT_CHECK(myb_term_radiation(&open, 0.0, x, zero, full) == MYB_OK);
    MT_CHECK(myb_term_radiation(&airless, 0.0, x, zero, bare) == MYB_OK);
    MT_CHECK(myb_term_radiation(&through_air, 0.0, x, zero, bent) == MYB_OK);
    const double to_sun[3] = {-x[0], -x[1], -x[2]};
    const double to_planet[3] = {planet_state[0] - x[0], planet_state[1] - x[1],
                                 planet_state[2] - x[2]};
    const double expected =
        myb_shadow_fraction_refracted(to_sun, sun_radius, to_planet, planet_radius, &air, limb);
    worst_air = fmax(worst_air, fabs(norm3(bent) / norm3(full) - expected));
    most_gain = fmax(most_gain, (norm3(bent) - norm3(bare)) / norm3(full));
  }
  printf("the radiation term with an atmosphere of the occulter: within %.1e of "
         "myb_shadow_fraction_refracted; the air lets up to %.1e of the full light into the "
         "shadow\n",
         worst_air, most_gain);
  MT_CHECK(worst_air < 1e-14 && most_gain > 1e-5);

  /* the same light on the plates of a spacecraft */
  const myb_plate plate = {.area = 3.0, .normal = {-1.0, 0.0, 0.0}, .specular = 0.2,
                           .diffuse = 0.3}; /* it faces the Sun */
  double still[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  myb_surface_radiation_term plates = {.env = env, .source = i_sun, .object = &object,
                                       .c = c_light, .plates = &plate, .plate_count = 1,
                                       .attitude = myb_attitude_fixed, .attitude_context = still,
                                       .occulters = {i_planet}, .occulter_count = 1};
  myb_surface_radiation_term lit = plates;
  lit.occulter_count = 0;
  myb_surface_radiation_term in_air = plates;
  in_air.atmospheres = atmospheres;
  const double x[3] = {au + 1.4e9, 5.5e6, 0.0};
  double full[3];
  double bare[3];
  double bent[3];
  MT_CHECK(myb_term_surface_radiation(&lit, 0.0, x, zero, full) == MYB_OK);
  MT_CHECK(myb_term_surface_radiation(&plates, 0.0, x, zero, bare) == MYB_OK);
  MT_CHECK(myb_term_surface_radiation(&in_air, 0.0, x, zero, bent) == MYB_OK);
  const double to_sun[3] = {-x[0], -x[1], -x[2]};
  const double to_planet[3] = {planet_state[0] - x[0], planet_state[1] - x[1],
                               planet_state[2] - x[2]};
  const double with_air =
      myb_shadow_fraction_refracted(to_sun, sun_radius, to_planet, planet_radius, &air, limb);
  const double without = myb_shadow_fraction(to_sun, sun_radius, to_planet, planet_radius, limb);
  printf("plates behind the planet: %.6f of the full light without air, %.6f with\n",
         norm3(bare) / norm3(full), norm3(bent) / norm3(full));
  MT_CHECK(norm3(full) > 0.0);
  MT_CHECK(fabs(norm3(bare) / norm3(full) - without) < 1e-14);
  MT_CHECK(fabs(norm3(bent) / norm3(full) - with_air) < 1e-14);
  myb_env_destroy(env);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_many();
  test_refraction();
  test_planet_gradient();
  test_terms();
  return mt_finish("test_shadow");
}
