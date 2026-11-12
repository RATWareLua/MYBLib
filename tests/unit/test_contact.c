/* SPDX-License-Identifier: BSD-3-Clause */
/* Contact and landing (myb/contact.h, myb_system_add_contact, myb_polyhedron_nearest):
 *  - the restitution of the law of Hunt and Crossley in closed form and its inverse, against
 *    the relation itself;
 *  - the nearest point of a mesh: a cube against its closed form, a star-shaped mesh against a
 *    search over all faces done here and the inside told by the solid angle;
 *  - a ball thrown at a plane without gravity: the restitution against the closed form; dropped
 *    under gravity: the difference is the impulse of the gravity during the contact, and falls
 *    with g in proportion; the bounces end at rest at the depth of the closed form, the energy
 *    never growing;
 *  - a block sliding with friction: the times against the closed form of the regularised law,
 *    and the stop against v0 / (mu g);
 *  - on a mesh: the rebound from a face against that from a plane, from an edge and from a
 *    corner along the line to the ball;
 *  - a ball on the equator of a turning asteroid (a term with an environment): carried along by
 *    the friction to the rotation of the ground, pressed by gravity less the centrifugal force;
 *  - a lander on four legs (a system with an integrated rotation): settling on a slope below
 *    the angle of friction at the creep speed of the closed form, sliding down a steeper one
 *    at g (sin a - mu cos a); the energy never growing;
 *  - momentum and angular momentum of a spinning ball bouncing off a turning body;
 *  - what is refused. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/attitude.h"
#include "myb/contact.h"
#include "myb/env.h"
#include "myb/geometry.h"
#include "myb/integ.h"
#include "myb/model.h"
#include "myb/orbit.h"
#include "myb/system.h"

static const double PI = 3.141592653589793;

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double norm3(const double a[3]) {
  return sqrt(dot3(a, a));
}

static void cross3(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

/* ---- The restitution of the law -------------------------------------------------------------- */

static void test_restitution(void) {
  const double es[] = {0.02, 0.3, 0.7, 0.95, 0.999, 0.999999};
  const double v = 1.3;
  double worst = 0.0;
  double residual = 0.0;
  for (size_t k = 0; k < sizeof es / sizeof es[0]; ++k) {
    double chi = 0.0;
    double back = 0.0;
    MT_CHECK(myb_contact_damping(es[k], v, &chi) == MYB_OK);
    MT_CHECK(myb_contact_restitution(chi, v, &back) == MYB_OK);
    worst = fmax(worst, fabs(back - es[k]) / es[k]);
    /* the relation itself, where the closed forms keep their digits (not where e chi v rounds
     * to 1) */
    const double y = chi * v;
    if (y > 0.1 && es[k] * y < 0.9) {
      const double up = y - log1p(y);
      const double down = -es[k] * y - log1p(-es[k] * y);
      residual = fmax(residual, fabs(up - down) / up);
    }
  }
  /* small and large y */
  double e_small = 0.0;
  double e_large = 0.0;
  MT_CHECK(myb_contact_restitution(1.0e-6, 1.0, &e_small) == MYB_OK);
  MT_CHECK(myb_contact_restitution(1.0e4, 1.0, &e_large) == MYB_OK);
  const double slope = (1.0 - e_small) / 1.0e-6;
  printf(
      "restitution: e -> chi -> e within %.1e, the relation %.1e; (1 - e) / (chi v) at 1e-6 "
      "%.9f (2/3), e chi v at 1e4 %.15f\n",
      worst, residual, slope, e_large * 1.0e4);
  MT_CHECK(worst < 1e-14 && residual < 1e-14);
  MT_CHECK(fabs(slope - 2.0 / 3.0) < 1e-6);
  MT_CHECK(fabs(e_large * 1.0e4 - 1.0) < 1e-12);
  double e = 0.0;
  double chi = -1.0;
  MT_CHECK(myb_contact_restitution(0.0, 3.0, &e) == MYB_OK && e == 1.0);
  MT_CHECK(myb_contact_restitution(2.0, 0.0, &e) == MYB_OK && e == 1.0);
  MT_CHECK(myb_contact_damping(1.0, 3.0, &chi) == MYB_OK && chi == 0.0);
  MT_CHECK(myb_contact_restitution(-1.0, 1.0, &e) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_restitution(1.0, (double)NAN, &e) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_restitution((double)INFINITY, 1.0, &e) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_damping(0.0, 1.0, &chi) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_damping(1.5, 1.0, &chi) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_damping(0.5, 0.0, &chi) == MYB_E_DOMAIN);
}

/* ---- The nearest point of a mesh ------------------------------------------------------------- */

static const double cube_vertices[24] = {-1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1,
                                         -1, -1, 1,  1, -1, 1,  -1, 1, 1,  1, 1, 1};
static const size_t cube_faces[36] = {0, 2, 1, 1, 2, 3, 4, 5, 6, 5, 7, 6, 0, 1, 4, 1, 5, 4,
                                      2, 6, 3, 3, 6, 7, 0, 4, 2, 2, 4, 6, 1, 3, 5, 3, 7, 5};

/* The signed distance to the cube |x|, |y|, |z| <= 1. */
static double cube_distance(const double p[3]) {
  double outside = 0.0;
  double inside = 1.0;
  for (int i = 0; i < 3; ++i) {
    const double d = fabs(p[i]) - 1.0;
    outside += d > 0.0 ? d * d : 0.0;
    inside = fmin(inside, -d);
  }
  return outside > 0.0 ? sqrt(outside) : -inside;
}

/* The distance from r to a triangle, done here in another way: the projection on its plane if
 * it falls inside, the nearest of the three edges otherwise. */
static double triangle_distance(const double* a, const double* b, const double* c,
                                const double r[3]) {
  double ab[3];
  double ac[3];
  double ar[3];
  double n[3];
  for (int i = 0; i < 3; ++i) {
    ab[i] = b[i] - a[i];
    ac[i] = c[i] - a[i];
    ar[i] = r[i] - a[i];
  }
  cross3(ab, ac, n);
  const double nn = norm3(n);
  for (int i = 0; i < 3; ++i) {
    n[i] /= nn;
  }
  const double h = dot3(ar, n);
  double p[3];
  for (int i = 0; i < 3; ++i) {
    p[i] = r[i] - h * n[i];
  }
  const double* corners[3] = {a, b, c};
  int inside = 1;
  for (int k = 0; k < 3; ++k) {
    const double* u = corners[k];
    const double* w = corners[(k + 1) % 3];
    double e[3];
    double q[3];
    double s[3];
    for (int i = 0; i < 3; ++i) {
      e[i] = w[i] - u[i];
      q[i] = p[i] - u[i];
    }
    cross3(e, q, s);
    inside = inside && dot3(s, n) >= 0.0;
  }
  if (inside) {
    return fabs(h);
  }
  double best = (double)INFINITY;
  for (int k = 0; k < 3; ++k) {
    const double* u = corners[k];
    const double* w = corners[(k + 1) % 3];
    double e[3];
    double q[3];
    for (int i = 0; i < 3; ++i) {
      e[i] = w[i] - u[i];
      q[i] = r[i] - u[i];
    }
    double s = dot3(q, e) / dot3(e, e);
    s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
    double d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = q[i] - s * e[i];
    }
    best = fmin(best, norm3(d));
  }
  return best;
}

enum { LATITUDES = 12, LONGITUDES = 24 };

enum { STAR_VERTICES = 2 + (LATITUDES - 1) * LONGITUDES };

enum { STAR_FACES = 2 * LONGITUDES + 2 * (LATITUDES - 2) * LONGITUDES };

/* A star-shaped, non-convex mesh: a sphere of latitudes and longitudes with a wavy radius. */
static void star_mesh(double* vertices, size_t* faces) {
  const size_t south = 1 + (LATITUDES - 1) * LONGITUDES;
  vertices[0] = 0.0;
  vertices[1] = 0.0;
  vertices[2] = 1.1;
  vertices[3 * south] = 0.0;
  vertices[3 * south + 1] = 0.0;
  vertices[3 * south + 2] = -1.1;
  for (int j = 1; j < LATITUDES; ++j) {
    const double theta = PI * j / (double)LATITUDES;
    for (int i = 0; i < LONGITUDES; ++i) {
      const double phi = 2.0 * PI * i / (double)LONGITUDES;
      const double r =
          1.0 + 0.3 * sin(theta) * sin(theta) * cos(3.0 * phi) + 0.1 * cos(2.0 * theta);
      const size_t at = 3 * (1 + (size_t)(j - 1) * LONGITUDES + (size_t)i);
      vertices[at] = r * sin(theta) * cos(phi);
      vertices[at + 1] = r * sin(theta) * sin(phi);
      vertices[at + 2] = r * cos(theta);
    }
  }
  size_t f = 0;
  for (int i = 0; i < LONGITUDES; ++i) {
    const size_t a = 1 + (size_t)i;
    const size_t b = 1 + (size_t)((i + 1) % LONGITUDES);
    faces[f++] = 0;
    faces[f++] = a;
    faces[f++] = b;
  }
  for (int j = 1; j < LATITUDES - 1; ++j) {
    for (int i = 0; i < LONGITUDES; ++i) {
      const size_t a = 1 + (size_t)(j - 1) * LONGITUDES + (size_t)i;
      const size_t b = 1 + (size_t)(j - 1) * LONGITUDES + (size_t)((i + 1) % LONGITUDES);
      const size_t c = a + LONGITUDES;
      const size_t d = b + LONGITUDES;
      faces[f++] = a;
      faces[f++] = c;
      faces[f++] = d;
      faces[f++] = a;
      faces[f++] = d;
      faces[f++] = b;
    }
  }
  for (int i = 0; i < LONGITUDES; ++i) {
    const size_t a = 1 + (size_t)(LATITUDES - 2) * LONGITUDES + (size_t)i;
    const size_t b = 1 + (size_t)(LATITUDES - 2) * LONGITUDES + (size_t)((i + 1) % LONGITUDES);
    faces[f++] = a;
    faces[f++] = south;
    faces[f++] = b;
  }
}

static void test_nearest(void) {
  /* a plane and a sphere */
  myb_contact_surface plane = {.shape = MYB_CONTACT_PLANE, .point = {0, 0, 1}, .normal = {0, 0, 2}};
  double point[3];
  double normal[3];
  double h = 0.0;
  MT_CHECK(myb_contact_surface_nearest(&plane, (const double[3]){3, 4, 5}, point, normal, &h) ==
           MYB_OK);
  MT_CHECK(h == 4.0 && point[0] == 3.0 && point[1] == 4.0 && point[2] == 1.0 && normal[2] == 1.0);
  myb_contact_surface ball = {.shape = MYB_CONTACT_ELLIPSOID, .radii = {2, 2, 2}};
  MT_CHECK(myb_contact_surface_nearest(&ball, (const double[3]){0, 0, 5}, point, normal, &h) ==
           MYB_OK);
  MT_CHECK(fabs(h - 3.0) < 1e-15 && fabs(point[2] - 2.0) < 1e-15 && fabs(normal[2] - 1.0) < 1e-15);

  /* the cube against its closed form */
  myb_polyhedron* cube = NULL;
  MT_CHECK(myb_polyhedron_create(cube_vertices, 8, cube_faces, 12, 1.0, NULL, &cube) == MYB_OK);
  mt_rng rng = {12345};
  double worst = 0.0;
  for (int k = 0; k < 4000; ++k) {
    const double p[3] = {6.0 * uniform(&rng) - 3.0, 6.0 * uniform(&rng) - 3.0,
                         6.0 * uniform(&rng) - 3.0};
    double d = 0.0;
    MT_CHECK(myb_polyhedron_nearest(cube, p, NULL, NULL, &d, NULL) == MYB_OK);
    worst = fmax(worst, fabs(d - cube_distance(p)));
  }
  /* the normal on a face, an edge, a corner, from inside */
  const double s2 = sqrt(0.5);
  const double s3 = sqrt(1.0 / 3.0);
  const double probes[4][3] = {{0.2, 0.3, 1.7}, {1.5, 1.5, 0.1}, {2.0, 2.0, 2.0}, {0.2, 0.3, 0.9}};
  const double normals[4][3] = {{0, 0, 1}, {s2, s2, 0}, {s3, s3, s3}, {0, 0, 1}};
  const double distances[4] = {0.7, sqrt(0.5), sqrt(3.0), -0.1};
  double worst_normal = 0.0;
  for (int k = 0; k < 4; ++k) {
    double d = 0.0;
    size_t face = 99;
    MT_CHECK(myb_polyhedron_nearest(cube, probes[k], point, normal, &d, &face) == MYB_OK);
    MT_CHECK(fabs(d - distances[k]) < 1e-15 && face < 12);
    for (int i = 0; i < 3; ++i) {
      worst_normal = fmax(worst_normal, fabs(normal[i] - normals[k][i]));
    }
  }
  /* a point on an edge itself: the normal of the edge */
  MT_CHECK(myb_polyhedron_nearest(cube, (const double[3]){0.4, 1.0, 1.0}, point, normal, &h,
                                  NULL) == MYB_OK);
  MT_CHECK(h == 0.0 && fabs(normal[0]) < 1e-15 && fabs(normal[1] - s2) < 1e-15 &&
           fabs(normal[2] - s2) < 1e-15);
  MT_CHECK(myb_polyhedron_nearest(cube, (const double[3]){(double)NAN, 0, 0}, NULL, NULL, &h,
                                  NULL) == MYB_E_DOMAIN);
  printf(
      "nearest point of a cube: distance within %.1e of the closed form (4000 points), "
      "normals of a face, an edge, a corner within %.1e\n",
      worst, worst_normal);
  MT_CHECK(worst < 1e-15 && worst_normal < 1e-15);
  myb_polyhedron_destroy(cube);

  /* a star-shaped mesh against all its faces, and inside against the solid angle */
  static double vertices[3 * STAR_VERTICES];
  static size_t faces[3 * STAR_FACES];
  star_mesh(vertices, faces);
  myb_polyhedron* star = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, STAR_VERTICES, faces, STAR_FACES, 1.0, NULL, &star) ==
           MYB_OK);
  if (star == NULL) {
    return;
  }
  double worst_distance = 0.0;
  int wrong_side = 0;
  int inside_count = 0;
  for (int k = 0; k < 3000; ++k) {
    const double p[3] = {3.2 * uniform(&rng) - 1.6, 3.2 * uniform(&rng) - 1.6,
                         3.0 * uniform(&rng) - 1.5};
    double d = 0.0;
    MT_CHECK(myb_polyhedron_nearest(star, p, point, normal, &d, NULL) == MYB_OK);
    double best = (double)INFINITY;
    for (size_t f = 0; f < STAR_FACES; ++f) {
      best =
          fmin(best, triangle_distance(vertices + 3 * faces[3 * f], vertices + 3 * faces[3 * f + 1],
                                       vertices + 3 * faces[3 * f + 2], p));
    }
    worst_distance = fmax(worst_distance, fabs(fabs(d) - best));
    const int inside = myb_polyhedron_solid_angle(star, p) > 2.0 * PI;
    inside_count += inside;
    wrong_side += (d < 0.0) != inside && fabs(d) > 1e-12;
  }
  printf(
      "nearest point of a star-shaped mesh of %d faces: within %.1e of the search over all "
      "faces, %d of 3000 points (%d inside) on the wrong side\n",
      STAR_FACES, worst_distance, wrong_side, inside_count);
  MT_CHECK(worst_distance < 1e-14 && wrong_side == 0 && inside_count > 300);
  myb_polyhedron_destroy(star);
}

/* ---- A ball on a fixed surface under constant gravity: a force model ------------------------- */

typedef struct gravity {
  double g[3];
} gravity;

static myb_status gravity_term(void* context, double t, const double x[3], const double v[3],
                               double a[3]) {
  (void)t;
  (void)x;
  (void)v;
  memcpy(a, ((const gravity*)context)->g, 3 * sizeof *a);
  return MYB_OK;
}

typedef struct scene {
  gravity g;
  myb_contact_term contact;
  myb_force_model* model;
  myb_orbit* orbit;
} scene;

static void scene_create(scene* s, const myb_contact_surface* surface, const myb_contact_law* law,
                         double radius, const double g[3], double rtol) {
  memset(s, 0, sizeof *s);
  memcpy(s->g.g, g, sizeof s->g.g);
  s->contact = (myb_contact_term){.env = NULL, .surface = *surface, .law = *law, .radius = radius};
  MT_CHECK(myb_force_model_create(NULL, &s->model) == MYB_OK);
  MT_CHECK(myb_force_model_add(s->model, "gravity", gravity_term, myb_term_no_jacobian, &s->g,
                               NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(s->model, "contact", myb_term_contact, NULL, &s->contact, NULL) ==
           MYB_OK);
  MT_CHECK(myb_force_model_add_break(s->model, myb_event_contact, &s->contact) == MYB_OK);
  const myb_orbit_options options = {
      .integrator = MYB_ORBIT_DOP853, .rtol = rtol, .atol = 1e-3 * rtol, .max_steps = 0};
  MT_CHECK(myb_orbit_create(s->model, &options, NULL, &s->orbit) == MYB_OK);
}

static void scene_destroy(scene* s) {
  myb_orbit_destroy(s->orbit);
  myb_force_model_destroy(s->model);
}

/* E = |v|^2 / 2 - g . x + the energy of the springs. */
static double scene_energy(const scene* s, const double y[6]) {
  const myb_contact_point sphere = {.offset = {0, 0, 0}, .radius = s->contact.radius};
  myb_contact_result r;
  if (myb_contact_force(&s->contact.surface, &s->contact.law, &sphere, 1, y, y + 3, NULL, NULL,
                        &r) != MYB_OK) {
    return (double)NAN;
  }
  return 0.5 * dot3(y + 3, y + 3) - dot3(s->g.g, y) + r.energy;
}

/* One contact from y0: the speeds along the normal n at the touch and at the parting. */
static void one_contact(scene* s, const double y0[6], const double n[3], double t_end, double* v_in,
                        double* v_out, double y_out[6]) {
  int found = 0;
  double t = 0.0;
  double y[6];
  MT_CHECK(myb_orbit_start(s->orbit, 0.0, y0, t_end) == MYB_OK);
  MT_CHECK(myb_orbit_event(s->orbit, myb_event_contact, &s->contact, t_end, -1, &found, &t, y) ==
           MYB_OK);
  MT_CHECK(found);
  *v_in = -dot3(y + 3, n);
  MT_CHECK(myb_orbit_event(s->orbit, myb_event_contact, &s->contact, t_end, 1, &found, &t, y) ==
           MYB_OK);
  MT_CHECK(found);
  *v_out = dot3(y + 3, n);
  if (y_out != NULL) {
    memcpy(y_out, y, sizeof y);
  }
}

static const myb_contact_surface ground = {
    .shape = MYB_CONTACT_PLANE, .point = {0, 0, 0}, .normal = {0, 0, 1}};

static void test_drop(void) {
  const double r = 0.1;
  const myb_contact_law law = {.stiffness = 1.0e9, .exponent = 1.5, .damping = 0.5};
  const double up[3] = {0, 0, 1};
  /* thrown at the plane without gravity */
  scene s;
  scene_create(&s, &ground, &law, r, (const double[3]){0, 0, 0}, 1e-13);
  const double y0[6] = {0.0, 0.0, r + 0.01, 0.3, 0.0, -1.0};
  double v_in = 0.0;
  double v_out = 0.0;
  one_contact(&s, y0, up, 1.0, &v_in, &v_out, NULL);
  double e_closed = 0.0;
  MT_CHECK(myb_contact_restitution(law.damping, v_in, &e_closed) == MYB_OK);
  const double free_error = fabs(v_out / v_in - e_closed);
  scene_destroy(&s);

  /* dropped from 0.5 m under the gravity of the Moon, and from 5 m under a tenth of it: the
   * same speed at the touch */
  double deviation[2] = {0.0, 0.0};
  double e_measured[2] = {0.0, 0.0};
  for (int k = 0; k < 2; ++k) {
    const double g = k == 0 ? 1.62 : 0.162;
    scene_create(&s, &ground, &law, r, (const double[3]){0, 0, -g}, 1e-13);
    const double drop[6] = {0.0, 0.0, r + 0.5 * 1.62 / g, 0.0, 0.0, 0.0};
    one_contact(&s, drop, up, 10.0, &v_in, &v_out, NULL);
    MT_CHECK(myb_contact_restitution(law.damping, v_in, &e_closed) == MYB_OK);
    e_measured[k] = v_out / v_in;
    deviation[k] = e_measured[k] - e_closed;
    scene_destroy(&s);
  }
  const double ratio = deviation[0] / deviation[1];
  printf(
      "restitution, thrown at a plane: against the closed form %.1e; dropped under 1.62 and "
      "0.162 m/s^2 to the same speed: e = %.12f, %.12f, off the closed form by %.2e and %.2e "
      "(ratio %.3f: the impulse of the gravity during the contact)\n",
      free_error, e_measured[0], e_measured[1], deviation[0], deviation[1], ratio);
  MT_CHECK(free_error < 1e-10);
  MT_CHECK(fabs(deviation[0]) < 1e-2 && ratio > 9.5 && ratio < 10.5);
}

/* The bounces of a ball dropped on the Moon until it rests; the energy at every step. */
static void test_rest(void) {
  const double r = 0.1;
  const double g = 1.62;
  const myb_contact_law law = {.stiffness = 1.0e6, .exponent = 1.5, .damping = 2.0};
  scene s;
  scene_create(&s, &ground, &law, r, (const double[3]){0, 0, -g}, 1e-12);
  const double y0[6] = {0.0, 0.0, r + 0.2, 0.0, 0.0, 0.0};
  const double t_end = 30.0;
  MT_CHECK(myb_orbit_start(s.orbit, 0.0, y0, t_end) == MYB_OK);
  double energy = scene_energy(&s, y0);
  const double energy0 = energy;
  double rise = 0.0;
  int bounces = 0;
  double y[6];
  /* bound: the samples of 30 s every millisecond */
  for (int k = 1; k <= 30000; ++k) {
    const double t = t_end * k / 30000.0;
    MT_CHECK(myb_orbit_state(s.orbit, t, y) == MYB_OK);
    const double e = scene_energy(&s, y);
    rise = fmax(rise, e - energy);
    energy = e;
  }
  /* the touches, counted again from the start */
  MT_CHECK(myb_orbit_start(s.orbit, 0.0, y0, t_end) == MYB_OK);
  /* bound: the bounces of 30 s */
  for (int k = 0; k < 100000; ++k) {
    int found = 0;
    double t = 0.0;
    MT_CHECK(myb_orbit_event(s.orbit, myb_event_contact, &s.contact, t_end, -1, &found, &t, y) ==
             MYB_OK);
    if (!found) {
      break;
    }
    ++bounces;
  }
  MT_CHECK(myb_orbit_state(s.orbit, t_end, y) == MYB_OK);
  const double depth = pow(g / law.stiffness, 1.0 / law.exponent);
  const double off = fabs(y[2] - (r - depth));
  printf(
      "dropped 0.2 m on the Moon: %d touches, at rest after %.0f s %.2e m from the depth "
      "(g / k)^(2/3) = %.4e m, speed %.1e m/s; energy %.4f -> %.6f m^2/s^2, never up by more "
      "than %.1e\n",
      bounces, t_end, off, depth, norm3(y + 3), energy0, energy, rise);
  MT_CHECK(bounces > 10);
  MT_CHECK(off < 1e-12 && norm3(y + 3) < 1e-9);
  MT_CHECK(rise < 1e-11);
  scene_destroy(&s);
}

/* ---- Sliding with friction ------------------------------------------------------------------- */

typedef struct speed_event {
  double target;
} speed_event;

static double speed_crossing(void* context, double t, const double* y) {
  (void)t;
  return hypot(y[3], y[4]) - ((const speed_event*)context)->target;
}

/* The time of the regularised law to slow from v0 to v: the integral of
 * sqrt(v^2 + a^2) / (mu g v) dv. */
static double slide_time(double v0, double v, double a, double mu_g) {
  const double s0 = sqrt(v0 * v0 + a * a);
  const double s = sqrt(v * v + a * a);
  return (s0 - s + a * log((a + s) * v0 / ((a + s0) * v))) / mu_g;
}

static void test_slide(void) {
  const double r = 0.05;
  const double g = 9.81;
  const double mu = 0.4;
  const double slip = 1.0e-3;
  const myb_contact_law law = {
      .stiffness = 1.0e5, .exponent = 1.0, .damping = 1.0, .friction = mu, .slip = slip};
  scene s;
  scene_create(&s, &ground, &law, r, (const double[3]){0, 0, -g}, 1e-13);
  /* at rest on the plane, pressed by its weight, and pushed along it */
  const double v0 = 2.0;
  const double y0[6] = {0.0, 0.0, r - g / law.stiffness, v0, 0.0, 0.0};
  const double stop = v0 / (mu * g);
  MT_CHECK(myb_orbit_start(s.orbit, 0.0, y0, 2.0 * stop) == MYB_OK);
  const double targets[3] = {0.5 * v0, 0.01 * v0, 10.0 * slip};
  double worst = 0.0;
  double times[3] = {0.0, 0.0, 0.0};
  double y[6];
  for (int k = 0; k < 3; ++k) {
    speed_event e = {targets[k]};
    int found = 0;
    MT_CHECK(myb_orbit_event(s.orbit, speed_crossing, &e, 2.0 * stop, -1, &found, &times[k], y) ==
             MYB_OK);
    MT_CHECK(found);
    worst = fmax(worst, fabs(times[k] - slide_time(v0, targets[k], slip, mu * g)));
  }
  MT_CHECK(myb_orbit_state(s.orbit, 2.0 * stop, y) == MYB_OK);
  const double travelled = y[0];
  const double coulomb = v0 * v0 / (2.0 * mu * g);
  printf(
      "sliding at %.1f m/s with mu = %.1f: half the speed at %.9f s (Coulomb %.9f), %.0e m/s "
      "at %.6f s, the regularised law to %.1e s; Coulomb stops at %.6f s; slid %.6f m "
      "(Coulomb %.6f)\n",
      v0, mu, times[0], 0.5 * stop, 10.0 * slip, times[2], worst, stop, travelled, coulomb);
  MT_CHECK(worst < 1e-9);
  MT_CHECK(fabs(times[0] - 0.5 * stop) < 1e-6 * stop);
  MT_CHECK(fabs(times[2] - stop) < 0.01 * stop && times[2] < stop);
  MT_CHECK(fabs(travelled - coulomb) < 1e-3 * coulomb);
  scene_destroy(&s);
}

/* ---- On a mesh: a face, an edge, a corner ---------------------------------------------------- */

static void test_mesh(void) {
  myb_polyhedron* cube = NULL;
  MT_CHECK(myb_polyhedron_create(cube_vertices, 8, cube_faces, 12, 1.0, NULL, &cube) == MYB_OK);
  const myb_contact_surface surface = {.shape = MYB_CONTACT_MESH, .mesh = cube};
  const myb_contact_law law = {.stiffness = 1.0e8, .exponent = 1.5, .damping = 0.3};
  const double r = 0.1;
  const double zero[3] = {0, 0, 0};
  const double s2 = sqrt(0.5);
  const double s3 = sqrt(1.0 / 3.0);
  /* aimed along -n at a point of the face, of the edge y = z = 1, of the corner (1, 1, 1) */
  const double targets[3][3] = {{0.2, 0.3, 1.0}, {0.3, 1.0, 1.0}, {1.0, 1.0, 1.0}};
  const double normals[3][3] = {{0, 0, 1}, {0, s2, s2}, {s3, s3, s3}};
  double outgoing[3][3];
  double e_mesh[3];
  scene s;
  for (int k = 0; k < 3; ++k) {
    scene_create(&s, &surface, &law, r, zero, 1e-13);
    double y0[6];
    for (int i = 0; i < 3; ++i) {
      y0[i] = targets[k][i] + (r + 0.02) * normals[k][i];
      y0[3 + i] = -normals[k][i];
    }
    double v_in = 0.0;
    double v_out = 0.0;
    double y[6];
    one_contact(&s, y0, normals[k], 1.0, &v_in, &v_out, y);
    e_mesh[k] = v_out / v_in;
    memcpy(outgoing[k], y + 3, sizeof outgoing[k]);
    scene_destroy(&s);
  }
  /* the plane, from the same height at the same speed */
  scene_create(&s, &ground, &law, r, zero, 1e-13);
  const double y0[6] = {0.2, 0.3, r + 0.02, 0.0, 0.0, -1.0};
  double v_in = 0.0;
  double v_out = 0.0;
  one_contact(&s, y0, (const double[3]){0, 0, 1}, 1.0, &v_in, &v_out, NULL);
  scene_destroy(&s);
  double e_closed = 0.0;
  MT_CHECK(myb_contact_restitution(law.damping, 1.0, &e_closed) == MYB_OK);
  double across = 0.0;
  for (int k = 0; k < 3; ++k) {
    double c[3];
    cross3(outgoing[k], normals[k], c);
    across = fmax(across, norm3(c) / norm3(outgoing[k]));
  }
  printf(
      "on a cube: the face %.15f, the plane %.15f, the edge %.15f, the corner %.15f "
      "(closed form %.15f); off the line to the ball %.1e\n",
      e_mesh[0], v_out / v_in, e_mesh[1], e_mesh[2], e_closed, across);
  MT_CHECK(fabs(e_mesh[0] - v_out / v_in) < 1e-12);
  MT_CHECK(fabs(e_mesh[1] - e_closed) < 1e-10 && fabs(e_mesh[2] - e_closed) < 1e-10);
  MT_CHECK(across < 1e-13);
  myb_polyhedron_destroy(cube);
}

/* ---- A ball carried by a turning asteroid (a term with an environment) ----------------------- */

typedef struct spinner {
  double rate;
} spinner;

static myb_status spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  const double w = ((const spinner*)context)->rate;
  const double angle = w * myb_time_seconds(t);
  const double c = cos(angle);
  const double sn = sin(angle);
  memset(o, 0, sizeof *o);
  const double m[9] = {c, sn, 0.0, -sn, c, 0.0, 0.0, 0.0, 1.0};
  const double dm[9] = {-sn * w, c * w, 0.0, -c * w, -sn * w, 0.0, 0.0, 0.0, 0.0};
  memcpy(o->rotation.m, m, sizeof m);
  memcpy(o->rotation.dm, dm, sizeof dm);
  o->ut = t;
  return MYB_OK;
}

/* A point just off a sphere where the Newton iteration of myb_ellipsoid_nearest once crept by
 * steps below the resolution of t + a^2 until its budget ran out (found by the ball below). */
static void test_sphere_nearest(void) {
  const double radii[3] = {1000.0, 1000.0, 1000.0};
  const double p[3] = {0x1.f437b5d0a9fa1p+9, -0x1.6af8dfcf22181p+3, 0.0};
  double point[3];
  double h = 0.0;
  MT_CHECK(myb_ellipsoid_nearest(radii, p, point, &h) == MYB_OK);
  const double exact = sqrt(p[0] * p[0] + p[1] * p[1]) - 1000.0;
  MT_CHECK(fabs(h - exact) < 1e-12 && fabs(norm3(point) - 1000.0) < 1e-12);
}

static void test_turning_ground(void) {
  const double radius = 1000.0;
  const double rate = 1.0e-3; /* rad/s: the ground moves at 1 m/s on the equator */
  const double gm = 1.0e4;    /* 0.01 m/s^2 on the surface */
  const double r = 0.5;
  spinner spin = {rate};
  myb_env* env = NULL;
  const myb_env_options eo = {.epoch = myb_time_from_seconds(0.0)};
  MT_CHECK(myb_env_create(&eo, NULL, &env) == MYB_OK);
  const myb_env_body asteroid = {
      .gm = gm, .orientation = spinner_orientation, .orientation_context = &spin};
  size_t body = 0;
  MT_CHECK(myb_env_add_body(env, &asteroid, &body) == MYB_OK);
  const myb_contact_law law = {
      .stiffness = 1.0e3, .exponent = 1.5, .damping = 5.0, .friction = 0.5, .slip = 1.0e-3};
  myb_contact_term term = {
      .env = env,
      .body = body,
      .surface = {.shape = MYB_CONTACT_ELLIPSOID, .radii = {radius, radius, radius}},
      .law = law,
      .radius = r};
  myb_point_mass_term pull = {gm};
  myb_force_model* model = NULL;
  MT_CHECK(myb_force_model_create(NULL, &model) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "gravity", myb_term_point_mass, NULL, &pull, NULL) == MYB_OK);
  MT_CHECK(myb_force_model_add(model, "contact", myb_term_contact, NULL, &term, NULL) == MYB_OK);
  const myb_orbit_options options = {.integrator = MYB_ORBIT_DOP853, .rtol = 1e-12, .atol = 1e-12};
  myb_orbit* orbit = NULL;
  MT_CHECK(myb_orbit_create(model, &options, NULL, &orbit) == MYB_OK);
  /* at rest in space on the equator: the ground slides under it */
  const double rest = pow(gm / (radius * radius) / law.stiffness, 1.0 / 1.5);
  const double y0[6] = {radius + r - rest, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double t_end = 1500.0;
  MT_CHECK(myb_orbit_start(orbit, 0.0, y0, t_end) == MYB_OK);
  double y[6];
  MT_CHECK(myb_orbit_state(orbit, t_end, y) == MYB_OK);
  /* carried along: the velocity of the ground under it, w x x */
  const double carried[3] = {-rate * y[1], rate * y[0], 0.0};
  const double relative[3] = {y[3] - carried[0], y[4] - carried[1], y[5] - carried[2]};
  const double centre = norm3(y);
  const double pressed = gm / (centre * centre) - rate * rate * centre;
  const double depth = pow(pressed / law.stiffness, 1.0 / 1.5);
  const double off = fabs(centre - (radius + r - depth));
  printf(
      "a ball on the equator of a turning asteroid: after %.0f s it moves with the ground to "
      "%.1e m/s, %.1e m from the depth of gravity less the centrifugal force (%.4e m)\n",
      t_end, norm3(relative), off, depth);
  MT_CHECK(norm3(relative) < 1e-9 && off < 1e-9);
  myb_orbit_destroy(orbit);
  myb_force_model_destroy(model);
  myb_env_destroy(env);
}

/* ---- A lander on its legs: a system with an integrated rotation ------------------------------ */

typedef struct first_order {
  myb_system* system;
  size_t n;
} first_order;

static myb_status system_ode(void* context, double t, const double* y, double* dydt) {
  const first_order* f = context;
  memcpy(dydt, y + f->n, f->n * sizeof *dydt);
  return myb_system_acceleration(f->system, t, y, y + f->n, dydt + f->n);
}

static myb_status weight(void* context, double t, const myb_system_view* view,
                         double acceleration[3], double torque[3]) {
  (void)t;
  (void)view;
  (void)torque;
  memcpy(acceleration, ((const gravity*)context)->g, 3 * sizeof *acceleration);
  return MYB_OK;
}

static const myb_contact_point legs[4] = {{{1.0, 1.0, -1.0}, 0.05},
                                          {{1.0, -1.0, -1.0}, 0.05},
                                          {{-1.0, 1.0, -1.0}, 0.05},
                                          {{-1.0, -1.0, -1.0}, 0.05}};
static const double lander_moments[3] = {0.7, 0.7, 0.9};

typedef struct landing {
  double slope;  /* rad */
  double along;  /* m/s down the slope at the end */
  double across; /* m/s off the slope's line */
  double spin;   /* rad/s at the end */
  double rise;   /* the most the energy rose between steps */
  double energy0;
  double energy1;
  double support; /* |contact force - weight| / g at the end */
  double accel;   /* down the slope, between the last two samples */
  size_t touching;
  long steps;
} landing;

/* The lander (body 1) dropped 5 cm onto a plane through the origin of body 0 tilted by the
 * slope, its legs square to it, with a small turn to rock it. */
static void land(double slope, double t_end, double sample, landing* out) {
  const double g = 1.62;
  const double mu = 0.6;
  /* per leg: k d^1.5 = g / 4 at 3 mm */
  const myb_contact_law law = {.stiffness = 0.25 * g / pow(3.0e-3, 1.5),
                               .exponent = 1.5,
                               .damping = 10.0,
                               .friction = mu,
                               .slip = 1.0e-3};
  const double n[3] = {sin(slope), 0.0, cos(slope)};
  const double down[3] = {cos(slope), 0.0, -sin(slope)};
  const double gm[2] = {0.0, 1.0e-20};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
  const myb_system_rotation rotation = {.moments = {0.7, 0.7, 0.9}};
  MT_CHECK(myb_system_set_rotation(s, 1, &rotation) == MYB_OK);
  gravity weight_of = {{0.0, 0.0, -g}};
  MT_CHECK(myb_system_add_action(s, 1, weight, &weight_of) == MYB_OK);
  const myb_system_contact contact = {
      .surface = {.shape = MYB_CONTACT_PLANE, .point = {0, 0, 0}, .normal = {n[0], n[1], n[2]}},
      .law = law,
      .points = legs,
      .point_count = 4};
  MT_CHECK(myb_system_add_contact(s, 1, 0, &contact) == MYB_OK);
  const size_t nc = myb_system_coordinates(s);
  MT_CHECK(nc == 15);
  double y[30];
  memset(y, 0, sizeof y);
  const double height = 1.05 + 0.05;
  for (int i = 0; i < 3; ++i) {
    y[3 + i] = height * n[i];
  }
  /* body axes: x down the slope, y level, z along the normal */
  const double m[9] = {down[0], down[1], down[2], 0.0, 1.0, 0.0, n[0], n[1], n[2]};
  double q[4];
  myb_quat_from_matrix(m, q);
  const double omega[3] = {0.05, -0.03, 0.02};
  MT_CHECK(myb_system_pack_rotation(s, 1, q, omega, NULL, y) == MYB_OK);
  first_order f = {s, nc};
  const myb_dop853_options options = {.rtol = 1e-10, .atol = 1e-12};
  myb_dop853* in = NULL;
  MT_CHECK(myb_dop853_create(2 * nc, &options, NULL, &in) == MYB_OK);
  MT_CHECK(myb_dop853_start(in, system_ode, &f, 0.0, y) == MYB_OK);
  /* E = |v|^2 / 2 + w . I w / 2 + g z + the springs */
  double energy = 0.0;
  double rise = 0.0;
  double last_along = 0.0;
  double last_t = 0.0;
  double accel = 0.0;
  long steps = 0;
  /* bound: the steps of t_end, within the budget of the integrator */
  for (long k = 0; k < 10000000; ++k) {
    const double* z = myb_dop853_state(in);
    double qz[4];
    double w[3];
    MT_CHECK(myb_system_unpack_rotation(s, 1, z, qz, w, NULL) == MYB_OK);
    double mz[9];
    myb_quat_to_matrix(qz, mz);
    myb_contact_result r;
    MT_CHECK(myb_contact_force(&contact.surface, &law, legs, 4, z + 3, z + nc + 3, mz, w, &r) ==
             MYB_OK);
    const double* v = z + nc + 3;
    const double e = 0.5 * dot3(v, v) +
                     0.5 * (lander_moments[0] * w[0] * w[0] + lander_moments[1] * w[1] * w[1] +
                            lander_moments[2] * w[2] * w[2]) +
                     g * z[5] + r.energy;
    if (k == 0) {
      out->energy0 = e;
    } else {
      rise = fmax(rise, e - energy);
    }
    energy = e;
    const double t = myb_dop853_time(in);
    if (t >= t_end) {
      out->along = dot3(v, down);
      double off[3];
      for (int i = 0; i < 3; ++i) {
        off[i] = v[i] - out->along * down[i];
      }
      out->across = norm3(off);
      out->spin = norm3(w);
      const double weight_back[3] = {r.force[0], r.force[1], r.force[2] - g};
      out->support = norm3(weight_back) / g;
      out->touching = r.touching;
      break;
    }
    if (t >= t_end - sample && last_t == 0.0) {
      last_t = t;
      last_along = dot3(v, down);
    }
    if (last_t > 0.0 && t - last_t > 0.5 * sample) {
      accel = (dot3(v, down) - last_along) / (t - last_t);
    }
    MT_CHECK(myb_dop853_step(in, t_end) == MYB_OK);
    ++steps;
  }
  out->slope = slope;
  out->rise = rise;
  out->energy1 = energy;
  out->accel = accel;
  out->steps = steps;
  myb_dop853_destroy(in);
  myb_system_destroy(s);
}

static void test_lander(void) {
  const double g = 1.62;
  const double mu = 0.6;
  const double slip = 1.0e-3;
  landing settle;
  memset(&settle, 0, sizeof settle);
  land(15.0 * PI / 180.0, 20.0, 2.0, &settle);
  const double t = tan(settle.slope) / mu;
  const double creep = slip * t / sqrt(1.0 - t * t);
  printf(
      "lander on a slope of 15 deg (mu %.1f): after 20 s it creeps at %.9e m/s, the closed "
      "form to %.1e; %.1e m/s across, turning at %.1e rad/s, %zu legs down, the legs bear its "
      "weight to %.1e; energy %.4f -> %.4f m^2/s^2, never up by more than %.1e (%ld steps)\n",
      mu, settle.along, fabs(settle.along / creep - 1.0), settle.across, settle.spin,
      settle.touching, settle.support, settle.energy0, settle.energy1, settle.rise, settle.steps);
  MT_CHECK(fabs(settle.along - creep) < 1e-10 * creep);
  MT_CHECK(settle.across < 1e-9 && settle.spin < 1e-9 && settle.touching == 4);
  MT_CHECK(settle.support < 1e-8);
  MT_CHECK(settle.rise < 1e-9);

  landing slide;
  memset(&slide, 0, sizeof slide);
  land(40.0 * PI / 180.0, 16.0, 4.0, &slide);
  const double expected = g * (sin(slide.slope) - mu * cos(slide.slope));
  printf(
      "lander on a slope of 40 deg: slides at %.3f m/s, accelerating at %.9f m/s^2, "
      "g (sin a - mu cos a) to %.1e; %zu legs down; energy never up by more than %.1e\n",
      slide.along, slide.accel, fabs(slide.accel / expected - 1.0), slide.touching, slide.rise);
  MT_CHECK(fabs(slide.accel - expected) < 1e-6 * expected);
  MT_CHECK(slide.touching == 4 && slide.rise < 1e-9);
}

/* ---- Momentum and angular momentum through a contact of two turning bodies ------------------- */

static void momenta(myb_system* s, const double* y, size_t nc, const double* gm,
                    const double moments[2][3], double p[3], double l[3]) {
  memset(p, 0, 3 * sizeof *p);
  memset(l, 0, 3 * sizeof *l);
  for (size_t b = 0; b < 2; ++b) {
    const double* x = y + 3 * b;
    const double* v = y + nc + 3 * b;
    double xv[3];
    cross3(x, v, xv);
    double q[4];
    double w[3];
    MT_CHECK(myb_system_unpack_rotation(s, b, y, q, w, NULL) == MYB_OK);
    double m[9];
    myb_quat_to_matrix(q, m);
    const double body[3] = {moments[b][0] * w[0], moments[b][1] * w[1], moments[b][2] * w[2]};
    for (int i = 0; i < 3; ++i) {
      const double spin = m[i] * body[0] + m[3 + i] * body[1] + m[6 + i] * body[2];
      p[i] += gm[b] * v[i];
      l[i] += gm[b] * (xv[i] + spin);
    }
  }
}

static void test_momentum(void) {
  const double G = 6.6743e-11;
  const double gm[2] = {G * 2.0e4, G * 10.0};
  /* a homogeneous ellipsoid of semi-axes 5, 4, 3 m and a ball of 0.1 m, per unit mass */
  const double moments[2][3] = {{5.0, 6.8, 8.2}, {0.004, 0.004, 0.004}};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(2, gm, NULL, &s) == MYB_OK);
  for (size_t b = 0; b < 2; ++b) {
    const myb_system_rotation rotation = {.moments = {moments[b][0], moments[b][1], moments[b][2]}};
    MT_CHECK(myb_system_set_rotation(s, b, &rotation) == MYB_OK);
  }
  const myb_contact_point ball[1] = {{{0.0, 0.0, 0.0}, 0.1}};
  const myb_system_contact contact = {
      .surface = {.shape = MYB_CONTACT_ELLIPSOID, .radii = {5.0, 4.0, 3.0}},
      .law = {.stiffness = 1.0e7, .exponent = 1.5, .damping = 0.5, .friction = 0.3, .slip = 1.0e-3},
      .points = ball,
      .point_count = 1};
  MT_CHECK(myb_system_add_contact(s, 1, 0, &contact) == MYB_OK);
  const size_t nc = myb_system_coordinates(s);
  double y[2 * 24];
  memset(y, 0, sizeof y);
  /* the ball comes in obliquely, spinning, at the turning body from outside its sphere */
  const double x1[3] = {6.0, 5.0, 4.8};
  const double v1[3] = {-0.7, -0.4, -0.5};
  memcpy(y + 3, x1, sizeof x1);
  memcpy(y + nc + 3, v1, sizeof v1);
  const double q0[4] = {0.9, 0.1, 0.3, -0.2};
  const double q1[4] = {1.0, 0.0, 0.0, 0.0};
  MT_CHECK(myb_system_pack_rotation(s, 0, q0, (const double[3]){0.01, 0.02, 0.2}, NULL, y) ==
           MYB_OK);
  MT_CHECK(myb_system_pack_rotation(s, 1, q1, (const double[3]){3.0, -2.0, 5.0}, NULL, y) ==
           MYB_OK);
  double p0[3];
  double l0[3];
  momenta(s, y, nc, gm, moments, p0, l0);
  first_order f = {s, nc};
  const myb_dop853_options options = {.rtol = 1e-13, .atol = 1e-15};
  myb_dop853* in = NULL;
  MT_CHECK(myb_dop853_create(2 * nc, &options, NULL, &in) == MYB_OK);
  MT_CHECK(myb_dop853_start(in, system_ode, &f, 0.0, y) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(in, 12.0) == MYB_OK);
  double p1[3];
  double l1[3];
  const double* z = myb_dop853_state(in);
  momenta(s, z, nc, gm, moments, p1, l1);
  double dp[3];
  double dl[3];
  for (int i = 0; i < 3; ++i) {
    dp[i] = p1[i] - p0[i];
    dl[i] = l1[i] - l0[i];
  }
  /* what the ball lost of its own momentum: the contact took place */
  double ball_change[3];
  for (int i = 0; i < 3; ++i) {
    ball_change[i] = gm[1] * (z[nc + 3 + (size_t)i] - v1[i]);
  }
  const double scale_p = norm3(ball_change);
  const double scale_l = gm[1] * norm3(x1) * norm3(v1);
  printf(
      "a spinning ball off a turning ellipsoid: the ball's momentum changed by %.3e of "
      "itself, the total by %.1e, the angular momentum by %.1e\n",
      scale_p / (gm[1] * norm3(v1)), norm3(dp) / scale_p, norm3(dl) / scale_l);
  MT_CHECK(scale_p > 0.5 * gm[1] * norm3(v1));
  MT_CHECK(norm3(dp) < 1e-11 * scale_p && norm3(dl) < 1e-10 * scale_l);
  myb_dop853_destroy(in);
  myb_system_destroy(s);
}

/* ---- What is refused ------------------------------------------------------------------------- */

static void test_refused(void) {
  const myb_contact_point point = {{0, 0, 0}, 0.1};
  const myb_contact_law good = {.stiffness = 1.0, .exponent = 1.5};
  const double x[3] = {0, 0, 0.05};
  const double v[3] = {0, 0, 0};
  myb_contact_result r;
  MT_CHECK(myb_contact_force(&ground, &good, &point, 1, x, v, NULL, NULL, &r) == MYB_OK);
  MT_CHECK(r.touching == 1 && r.gap < 0.0 && r.force[2] > 0.0);
  const myb_contact_law bad[5] = {{.stiffness = 0.0, .exponent = 1.5},
                                  {.stiffness = 1.0, .exponent = 0.5},
                                  {.stiffness = 1.0, .exponent = 1.0, .damping = -1.0},
                                  {.stiffness = 1.0, .exponent = 1.0, .friction = 0.5},
                                  {.stiffness = (double)NAN, .exponent = 1.0}};
  for (int k = 0; k < 5; ++k) {
    MT_CHECK(myb_contact_force(&ground, &bad[k], &point, 1, x, v, NULL, NULL, &r) == MYB_E_DOMAIN);
  }
  const myb_contact_surface flat = {.shape = MYB_CONTACT_PLANE, .normal = {0, 0, 0}};
  const myb_contact_surface egg = {.shape = MYB_CONTACT_ELLIPSOID, .radii = {1, 0, 1}};
  const myb_contact_surface none = {.shape = MYB_CONTACT_MESH, .mesh = NULL};
  MT_CHECK(myb_contact_force(&flat, &good, &point, 1, x, v, NULL, NULL, &r) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_force(&egg, &good, &point, 1, x, v, NULL, NULL, &r) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_force(&none, &good, &point, 1, x, v, NULL, NULL, &r) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_force(&ground, &good, &point, 0, x, v, NULL, NULL, &r) == MYB_E_DOMAIN);
  const myb_contact_point negative = {{0, 0, 0}, -0.1};
  MT_CHECK(myb_contact_force(&ground, &good, &negative, 1, x, v, NULL, NULL, &r) == MYB_E_DOMAIN);
  MT_CHECK(myb_contact_force(&ground, &good, &point, 1, (const double[3]){0, (double)NAN, 0}, v,
                             NULL, NULL, &r) == MYB_E_DOMAIN);
  /* the term and its event */
  myb_contact_term term = {.surface = flat, .law = good, .radius = 0.1};
  double a[3];
  MT_CHECK(myb_term_contact(&term, 0.0, x, v, a) == MYB_E_DOMAIN);
  const double y[6] = {0, 0, 0.05, 0, 0, 0};
  MT_CHECK(isnan(myb_event_contact(&term, 0.0, y)));
  term.surface = ground;
  MT_CHECK(myb_event_contact(&term, 0.0, y) == 0.05 - 0.1);
  /* the system */
  const double gm[3] = {1.0, 1.0, 1.0};
  myb_system* s = NULL;
  MT_CHECK(myb_system_create(3, gm, NULL, &s) == MYB_OK);
  myb_system_contact c = {.surface = ground, .law = good, .points = &point, .point_count = 1};
  MT_CHECK(myb_system_add_contact(s, 1, 1, &c) == MYB_E_DOMAIN);
  MT_CHECK(myb_system_add_contact(s, 1, 3, &c) == MYB_E_DOMAIN);
  c.point_count = 0;
  MT_CHECK(myb_system_add_contact(s, 1, 0, &c) == MYB_E_DOMAIN);
  c.point_count = 1;
  c.law.exponent = 0.9;
  MT_CHECK(myb_system_add_contact(s, 1, 0, &c) == MYB_E_DOMAIN);
  c.law = good;
  for (int k = 0; k < MYB_SYSTEM_MAX_CONTACTS; ++k) {
    MT_CHECK(myb_system_add_contact(s, 1 + (size_t)(k % 2), 0, &c) == MYB_OK);
  }
  MT_CHECK(myb_system_add_contact(s, 1, 0, &c) == MYB_E_OVERFLOW);
  myb_system_destroy(s);
}

int main(void) {
  test_restitution();
  test_nearest();
  test_drop();
  test_rest();
  test_slide();
  test_mesh();
  test_sphere_nearest();
  test_turning_ground();
  test_lander();
  test_momentum();
  test_refused();
  return mt_finish("test_contact");
}
