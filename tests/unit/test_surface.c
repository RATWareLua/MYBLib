/* SPDX-License-Identifier: BSD-3-Clause */
/* The surface of the object (myb/surface.h):
 *  - free-molecular flow (Sentman): plates laid over a sphere by a quadrature integrate to the
 *    closed form for the sphere and to no lift; a plate across a fast flow; the direction of
 *    the lift;
 *  - light on a mesh: a cube as six plates; two boxes, one in the shadow of the other — the
 *    lit area is the silhouette; a faceted sphere against (1 + 4/9) pi R^2 for a Lambertian
 *    surface; the torque of Rubincam's model averaged over a rotation: none for a symmetric
 *    body, opposite for a body and its mirror image (the YORP effect);
 *  - attitude laws;
 *  - the terms: a plate against the sphere of myb_term_radiation, in the light and in the
 *    shadow; plates over a sphere in the gas against myb_term_drag with the coefficient of the
 *    sphere; a wind. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double pi = 3.141592653589793;
static const double two_pi = 6.283185307179586;
static const double c_light = 299792458.0;

static double norm3(const double v[3]) {
  return hypot(hypot(v[0], v[1]), v[2]);
}

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* ---- Plates over a sphere: Gauss-Legendre in the cosine of the angle from an axis, on each
 * hemisphere apart, equal steps round the axis; the areas add up to 4 pi R^2. With the axis
 * along the flow the force of the gas depends on that cosine alone, and steeply only at the
 * edge between the hemispheres, where the nodes are densest ----------------------------------- */

enum { RINGS = 40, SECTORS = 48, SPHERE_PLATES = 2 * RINGS * SECTORS };

static const double rings = 40.0;   /* RINGS */
static const double sectors = 48.0; /* SECTORS */

static myb_plate sphere_plates[SPHERE_PLATES];

static void make_sphere(double radius, const double axis[3]) {
  /* two directions across the axis */
  const double an = norm3(axis);
  const double z[3] = {axis[0] / an, axis[1] / an, axis[2] / an};
  const double helper[3] = {fabs(z[0]) < 0.6 ? 1.0 : 0.0, fabs(z[0]) < 0.6 ? 0.0 : 1.0, 0.0};
  double e1[3] = {helper[1] * z[2] - helper[2] * z[1], helper[2] * z[0] - helper[0] * z[2],
                  helper[0] * z[1] - helper[1] * z[0]};
  const double e1n = norm3(e1);
  for (int i = 0; i < 3; ++i) {
    e1[i] /= e1n;
  }
  const double e2[3] = {z[1] * e1[2] - z[2] * e1[1], z[2] * e1[0] - z[0] * e1[2],
                        z[0] * e1[1] - z[1] * e1[0]};
  /* nodes and weights of Gauss-Legendre on [-1, 1] by Newton's iteration */
  for (int i = 0; i < RINGS; ++i) {
    double x = cos(pi * (i + 0.75) / (rings + 0.5));
    double dp = 1.0;
    for (int iteration = 0; iteration < 100; ++iteration) {
      double p0 = 1.0;
      double p1 = x;
      for (int k = 2; k <= RINGS; ++k) {
        const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = rings * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / dp;
      x -= step;
      if (fabs(step) < 1e-16) {
        break;
      }
    }
    const double weight = 2.0 / ((1.0 - x * x) * dp * dp);
    for (int half = 0; half < 2; ++half) {
      /* the node on [0, 1] and on [-1, 0] */
      const double cosine = half == 0 ? 0.5 * (x + 1.0) : 0.5 * (x - 1.0);
      const double sine = sqrt(1.0 - cosine * cosine);
      for (int j = 0; j < SECTORS; ++j) {
        const double lon = two_pi * (j + 0.5) / sectors;
        myb_plate* p = &sphere_plates[(2 * i + half) * SECTORS + j];
        p->area = radius * radius * 0.5 * weight * two_pi / sectors;
        for (int q = 0; q < 3; ++q) {
          p->normal[q] = sine * cos(lon) * e1[q] + sine * sin(lon) * e2[q] + cosine * z[q];
        }
        const double n = norm3(p->normal);
        for (int q = 0; q < 3; ++q) {
          p->normal[q] /= n; /* unit to the last bit */
        }
      }
    }
  }
}

static void test_gas(void) {
  const double radius = 0.7;
  const double flow[3] = {0.36, -0.48, 0.8};
  make_sphere(radius, flow);
  double area = 0.0;
  for (int k = 0; k < SPHERE_PLATES; ++k) {
    area += sphere_plates[k].area;
  }
  MT_CHECK(fabs(area / (4.0 * pi * radius * radius) - 1.0) < 1e-13);
  double worst = 0.0;
  double worst_lift = 0.0;
  static const double speeds[4] = {300.0, 1500.0, 7600.0, 30000.0};
  static const double accommodation[3] = {0.0, 0.8, 1.0};
  for (int s = 0; s < 4; ++s) {
    for (int k = 0; k < 3; ++k) {
      const myb_gas gas = {.temperature = 900.0, .molar_mass = 0.016, .wall_temperature = 300.0,
                           .accommodation = accommodation[k]};
      const double v[3] = {0.36 * speeds[s], -0.48 * speeds[s], 0.8 * speeds[s]};
      double a[3] = {0.0, 0.0, 0.0};
      const double density = 3.0e-12;
      const double mass = 50.0;
      MT_CHECK(myb_force_gas_plates(density, &gas, mass, sphere_plates, SPHERE_PLATES, v, a) ==
               MYB_OK);
      double cd = 0.0;
      MT_CHECK(myb_gas_sphere_cd(&gas, speeds[s], &cd) == MYB_OK);
      const double expected = 0.5 * density * speeds[s] * speeds[s] * cd * pi * radius * radius /
                              mass;
      const double along = -dot3(a, v) / speeds[s];
      const double across[3] = {a[0] + along * v[0] / speeds[s], a[1] + along * v[1] / speeds[s],
                                a[2] + along * v[2] / speeds[s]};
      worst = fmax(worst, fabs(along / expected - 1.0));
      worst_lift = fmax(worst_lift, norm3(across) / expected);
      if (s == 2 && k == 1) {
        printf("a sphere at 7.6 km/s in atomic oxygen at 900 K, accommodation 0.8: C_D = %.4f\n",
               cd);
        MT_CHECK(cd > 2.0 && cd < 2.6);
      }
    }
  }
  printf("Sentman's plates over a sphere against the closed form: drag within %.1e, lift %.1e "
         "of it\n",
         worst, worst_lift);
  MT_CHECK(worst < 1e-12 && worst_lift < 1e-13);

  /* one plate across a flow much faster than the molecules: C_D = 2 + sqrt(pi) v_out / v_in */
  const myb_gas cold = {.temperature = 1.0, .molar_mass = 0.016, .wall_temperature = 300.0,
                        .accommodation = 1.0};
  const myb_plate front = {.area = 2.0, .normal = {1.0, 0.0, 0.0}};
  const double fast[3] = {7600.0, 0.0, 0.0};
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_gas_plates(1.0e-11, &cold, 10.0, &front, 1, fast, a) == MYB_OK);
  const double leaving = sqrt(0.5 * 4.0 * 8.314462618 * 300.0 / (0.016 * 7600.0 * 7600.0));
  const double cd = -a[0] / (0.5 * 1.0e-11 * 7600.0 * 7600.0 * 2.0 / 10.0);
  MT_CHECK(fabs(cd - (2.0 + sqrt(pi) * leaving)) < 1e-4);
  MT_CHECK(a[1] == 0.0 && a[2] == 0.0);
  /* the same plate turned away feels next to nothing; tilted, it is pushed along -normal */
  const double back[3] = {-7600.0, 0.0, 0.0};
  double b[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_gas_plates(1.0e-11, &cold, 10.0, &front, 1, back, b) == MYB_OK);
  MT_CHECK(fabs(b[0]) < 1e-30 * fabs(a[0]) + 1e-300);
  const double tilted[3] = {7600.0 * 0.8, 7600.0 * 0.6, 0.0};
  double t[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_gas_plates(1.0e-11, &cold, 10.0, &front, 1, tilted, t) == MYB_OK);
  /* drag against the motion; and besides it a lift that pushes the plate back along -x */
  MT_CHECK(dot3(t, tilted) < 0.0);
  const double lift_x = t[0] - dot3(t, tilted) * tilted[0] / (7600.0 * 7600.0);
  MT_CHECK(lift_x < 0.0 && t[2] == 0.0);
  /* errors */
  const double still[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_gas_plates(1.0e-11, &cold, 10.0, &front, 1, still, t) == MYB_E_DOMAIN);
  myb_gas bad = cold;
  bad.accommodation = 1.5;
  MT_CHECK(myb_force_gas_plates(1.0e-11, &bad, 10.0, &front, 1, fast, t) == MYB_E_DOMAIN);
  double unused = 0.0;
  MT_CHECK(myb_gas_sphere_cd(&bad, 7600.0, &unused) == MYB_E_DOMAIN);
  const myb_plate skew = {.area = 1.0, .normal = {1.0, 0.1, 0.0}};
  MT_CHECK(myb_force_gas_plates(1.0e-11, &cold, 10.0, &skew, 1, fast, t) == MYB_E_DOMAIN);
}

/* ---- Meshes ------------------------------------------------------------------------------------ */

/* A box [lo, hi] as 8 vertices and 12 triangles, appended to a mesh. */
static void add_box(double* vertices, size_t* vertex_count, size_t* faces, size_t* face_count,
                    const double lo[3], const double hi[3]) {
  static const size_t quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4},
                                     {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
  const size_t base = *vertex_count;
  for (int k = 0; k < 8; ++k) {
    double* v = vertices + 3 * (base + (size_t)k);
    v[0] = (k & 1) ? hi[0] : lo[0];
    v[1] = (k & 2) ? hi[1] : lo[1];
    v[2] = (k & 4) ? hi[2] : lo[2];
  }
  for (int q = 0; q < 6; ++q) {
    size_t* f = faces + 3 * *face_count;
    f[0] = base + quads[q][0];
    f[1] = base + quads[q][1];
    f[2] = base + quads[q][2];
    f[3] = base + quads[q][0];
    f[4] = base + quads[q][2];
    f[5] = base + quads[q][3];
    *face_count += 2;
  }
  *vertex_count += 8;
}

/* An icosahedron with every face split in four `levels` times, its vertices on the sphere. */
enum { MAX_VERTICES = 700, MAX_FACES = 1300 };

static double ico_vertices[3 * MAX_VERTICES];
static size_t ico_faces[3 * MAX_FACES];
static size_t ico_vertex_count;
static size_t ico_face_count;

static size_t midpoint(size_t a, size_t b) {
  double m[3];
  for (int i = 0; i < 3; ++i) {
    m[i] = 0.5 * (ico_vertices[3 * a + (size_t)i] + ico_vertices[3 * b + (size_t)i]);
  }
  const double n = norm3(m);
  for (int i = 0; i < 3; ++i) {
    m[i] /= n;
  }
  for (size_t k = 0; k < ico_vertex_count; ++k) {
    if (fabs(ico_vertices[3 * k] - m[0]) + fabs(ico_vertices[3 * k + 1] - m[1]) +
            fabs(ico_vertices[3 * k + 2] - m[2]) < 1e-9) {
      return k;
    }
  }
  memcpy(ico_vertices + 3 * ico_vertex_count, m, sizeof m);
  return ico_vertex_count++;
}

static void make_icosphere(int levels) {
  const double g = (1.0 + sqrt(5.0)) / 2.0;
  const double v[12][3] = {{-1, g, 0}, {1, g, 0},  {-1, -g, 0}, {1, -g, 0},
                           {0, -1, g}, {0, 1, g},  {0, -1, -g}, {0, 1, -g},
                           {g, 0, -1}, {g, 0, 1},  {-g, 0, -1}, {-g, 0, 1}};
  static const size_t faces[60] = {
      0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
      3, 9,  4, 3, 4, 2, 3, 2, 6, 3, 6, 8,  3, 8,  9,  4, 9, 5, 2, 4,  11, 6, 2,  10, 8, 6, 7, 9, 8, 1};
  ico_vertex_count = 12;
  ico_face_count = 20;
  for (int k = 0; k < 12; ++k) {
    const double n = norm3(v[k]);
    for (int i = 0; i < 3; ++i) {
      ico_vertices[3 * k + i] = v[k][i] / n;
    }
  }
  memcpy(ico_faces, faces, sizeof faces);
  static size_t next[3 * MAX_FACES];
  for (int level = 0; level < levels; ++level) {
    size_t count = 0;
    for (size_t f = 0; f < ico_face_count; ++f) {
      const size_t a = ico_faces[3 * f];
      const size_t b = ico_faces[3 * f + 1];
      const size_t c = ico_faces[3 * f + 2];
      const size_t ab = midpoint(a, b);
      const size_t bc = midpoint(b, c);
      const size_t ca = midpoint(c, a);
      const size_t four[12] = {a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca};
      memcpy(next + 3 * count, four, sizeof four);
      count += 4;
    }
    memcpy(ico_faces, next, 3 * count * sizeof(size_t));
    ico_face_count = count;
  }
}

/* The torque about z of Rubincam's model averaged over a rotation about z, the Sun in the
 * plane of the equator. */
static double spin_torque(const myb_polyhedron* shape, int phases) {
  const myb_plate lambert = {.thermal = 1};
  double sum = 0.0;
  for (int k = 0; k < phases; ++k) {
    const double angle = two_pi * (k + 0.5) / phases;
    const double u[3] = {-cos(angle), -sin(angle), 0.0}; /* the light goes along u */
    double torque[3];
    const double un = norm3(u);
    const double unit[3] = {u[0] / un, u[1] / un, 0.0};
    MT_CHECK(myb_force_radiation_mesh(shape, &lambert, 1.0, unit, 0, NULL, torque) == MYB_OK);
    sum += torque[2];
  }
  return sum / phases;
}

static void test_mesh_light(void) {
  static double vertices[3 * 16];
  static size_t faces[3 * 24];
  size_t vertex_count = 0;
  size_t face_count = 0;
  const double lo[3] = {-0.5, -0.5, -0.5};
  const double hi[3] = {0.5, 0.5, 0.5};
  add_box(vertices, &vertex_count, faces, &face_count, lo, hi);
  myb_polyhedron* cube = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, vertex_count, faces, face_count, 1.0, NULL, &cube) ==
           MYB_OK);
  /* a cube is six plates */
  const myb_plate surface = {.specular = 0.3, .diffuse = 0.2, .thermal = 1};
  myb_plate plates[6];
  for (int k = 0; k < 6; ++k) {
    plates[k] = surface;
    plates[k].area = 1.0;
    plates[k].normal[k / 2] = k % 2 ? 1.0 : -1.0;
  }
  const double raw[3] = {-0.5, 0.3, -0.81};
  const double rn = norm3(raw);
  const double u[3] = {raw[0] / rn, raw[1] / rn, raw[2] / rn};
  double force[3];
  double torque[3];
  double ref[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_radiation_mesh(cube, &surface, 4.5e-6, u, 0, force, torque) == MYB_OK);
  MT_CHECK(myb_force_radiation_plates(4.5e-6, 1.0, plates, 6, u, ref) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(force[i] - ref[i]) / norm3(ref));
  }
  printf("light on a cube as a mesh and as six plates: %.1e apart; torque %.1e of force x "
         "size\n",
         worst, norm3(torque) / norm3(ref));
  MT_CHECK(worst < 1e-15 && norm3(torque) < 1e-15 * norm3(ref));
  /* a convex body casts no shadow on itself */
  double shadowed[3];
  MT_CHECK(myb_force_radiation_mesh(cube, &surface, 4.5e-6, u, 3, shadowed, NULL) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(shadowed[i] - force[i]) < 1e-13 * norm3(force));
  }
  myb_polyhedron_destroy(cube);

  /* the cube in front of a slab of 3 x 3: the light along -x reaches 8 of its 9 square metres */
  const double slab_lo[3] = {-3.2, -1.5, -1.5};
  const double slab_hi[3] = {-3.0, 1.5, 1.5};
  add_box(vertices, &vertex_count, faces, &face_count, slab_lo, slab_hi);
  myb_polyhedron* two = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, vertex_count, faces, face_count, 1.0, NULL, &two) ==
           MYB_OK);
  const myb_plate black = {0};
  const double along[3] = {-1.0, 0.0, 0.0};
  MT_CHECK(myb_force_radiation_mesh(two, &black, 1.0, along, 0, force, NULL) == MYB_OK);
  MT_CHECK(fabs(force[0] + 10.0) < 1e-13); /* without shadows: both faces whole */
  MT_CHECK(myb_force_radiation_mesh(two, &black, 1.0, along, 24, force, torque) == MYB_OK);
  printf("a cube before a slab of 9 m^2, shadows by 24^2 parts a face: lit area %.12f m^2 (the "
         "silhouette: 9)\n",
         -force[0]);
  MT_CHECK(fabs(force[0] + 9.0) < 1e-12 && fabs(force[1]) < 1e-15 && fabs(force[2]) < 1e-15);
  MT_CHECK(norm3(torque) < 1e-12);
  /* from the other side the slab hides the cube: the silhouette again */
  const double back[3] = {1.0, 0.0, 0.0};
  MT_CHECK(myb_force_radiation_mesh(two, &black, 1.0, back, 24, force, NULL) == MYB_OK);
  MT_CHECK(fabs(force[0] - 9.0) < 1e-12);
  /* a coarse division is off by the parts cut by the edge of the shadow */
  MT_CHECK(myb_force_radiation_mesh(two, &black, 1.0, along, 5, force, NULL) == MYB_OK);
  MT_CHECK(fabs(force[0] + 9.0) < 0.5 && fabs(force[0] + 9.0) > 1e-6);
  const double skew[3] = {-1.0, 0.1, 0.0};
  MT_CHECK(myb_force_radiation_mesh(two, &black, 1.0, skew, 0, force, NULL) == MYB_E_DOMAIN);
  myb_polyhedron_destroy(two);

  /* a sphere of facets, absorbing and re-emitting as a Lambertian surface:
   * (1 + 4/9) pi R^2 times the pressure */
  make_icosphere(3);
  myb_polyhedron* ball = NULL;
  MT_CHECK(myb_polyhedron_create(ico_vertices, ico_vertex_count, ico_faces, ico_face_count, 1.0,
                                 NULL, &ball) == MYB_OK);
  const myb_plate lambert = {.thermal = 1};
  MT_CHECK(myb_force_radiation_mesh(ball, &lambert, 1.0, u, 0, force, torque) == MYB_OK);
  const double expected = (1.0 + 4.0 / 9.0) * pi;
  printf("a sphere of %zu facets, Lambertian: force %.4f of (1 + 4/9) pi R^2 p, off the "
         "direction of the light by %.1e\n",
         ico_face_count, dot3(force, u) / expected,
         sqrt(fmax(0.0, dot3(force, force) - dot3(force, u) * dot3(force, u))) / expected);
  MT_CHECK(fabs(dot3(force, u) / expected - 1.0) < 0.01);
  MT_CHECK(norm3(torque) < 1e-3 * expected);
  /* symmetric about planes through the axis: no torque along it on the average */
  const double symmetric = spin_torque(ball, 72);
  myb_polyhedron_destroy(ball);

  /* a body with a twist: two wedges on the equator leaning the same way round (Rubincam's
   * windmill) — every vertex of the sphere near the equator is pushed along the direction of
   * rotation in proportion to its height over the mean radius of a two-lobed bump */
  static double twisted[3 * MAX_VERTICES];
  static double mirrored[3 * MAX_VERTICES];
  static size_t mirrored_faces[3 * MAX_FACES];
  for (size_t k = 0; k < ico_vertex_count; ++k) {
    const double* v = ico_vertices + 3 * k;
    const double lon = atan2(v[1], v[0]);
    const double equator = 1.0 - v[2] * v[2];
    /* a saw-tooth in longitude with two teeth: slow rise, steep fall */
    const double phase = fmod(lon / pi + 2.0, 1.0);
    const double height = 0.25 * equator * equator * phase;
    for (int i = 0; i < 3; ++i) {
      twisted[3 * k + (size_t)i] = v[i] * (1.0 + height);
    }
    mirrored[3 * k] = twisted[3 * k];
    mirrored[3 * k + 1] = -twisted[3 * k + 1];
    mirrored[3 * k + 2] = twisted[3 * k + 2];
  }
  for (size_t f = 0; f < ico_face_count; ++f) { /* the mirror turns the faces over */
    mirrored_faces[3 * f] = ico_faces[3 * f];
    mirrored_faces[3 * f + 1] = ico_faces[3 * f + 2];
    mirrored_faces[3 * f + 2] = ico_faces[3 * f + 1];
  }
  myb_polyhedron* mill = NULL;
  myb_polyhedron* mirror = NULL;
  MT_CHECK(myb_polyhedron_create(twisted, ico_vertex_count, ico_faces, ico_face_count, 1.0, NULL,
                                 &mill) == MYB_OK);
  MT_CHECK(myb_polyhedron_create(mirrored, ico_vertex_count, mirrored_faces, ico_face_count, 1.0,
                                 NULL, &mirror) == MYB_OK);
  const double yorp = spin_torque(mill, 72);
  const double yorp_mirror = spin_torque(mirror, 72);
  printf("the torque about the axis averaged over a rotation (in p R^3): sphere %.1e, a body "
         "with a twist %.4e, its mirror image %.4e\n",
         symmetric, yorp, yorp_mirror);
  MT_CHECK(fabs(symmetric) < 1e-12);
  MT_CHECK(fabs(yorp) > 1e-4);
  MT_CHECK(fabs(yorp + yorp_mirror) < 1e-12 * fabs(yorp) + 1e-15);
  myb_polyhedron_destroy(mill);
  myb_polyhedron_destroy(mirror);
}

/* ---- Attitude laws and the terms -------------------------------------------------------------- */

typedef struct mover {
  double radius;
  double rate;
  double phase;
} mover;

static myb_status mover_state(void* context, myb_time t, double y[6]) {
  const mover* m = context;
  const double u = m->phase + m->rate * myb_time_seconds(t);
  y[0] = m->radius * myb_cos(u);
  y[1] = m->radius * myb_sin(u);
  y[2] = 0.0;
  y[3] = -m->radius * m->rate * myb_sin(u);
  y[4] = m->radius * m->rate * myb_cos(u);
  y[5] = 0.0;
  return MYB_OK;
}

typedef struct spinner {
  double rate;
} spinner;

static myb_status spinner_orientation(void* context, myb_time t, myb_orientation* o) {
  const spinner* s = context;
  const double w = s->rate * myb_time_seconds(t);
  const double c = myb_cos(w);
  const double sn = myb_sin(w);
  const double m[9] = {c, sn, 0.0, -sn, c, 0.0, 0.0, 0.0, 1.0};
  const double dm[9] = {-sn * s->rate, c * s->rate, 0.0, -c * s->rate, -sn * s->rate, 0.0,
                        0.0,           0.0,         0.0};
  memcpy(o->rotation.m, m, sizeof m);
  memcpy(o->rotation.dm, dm, sizeof dm);
  o->ut = t;
  o->xp = 0.0;
  o->yp = 0.0;
  return MYB_OK;
}

static myb_status constant_density(void* context, myb_time t, const double r[3],
                                   double* density) {
  (void)t;
  (void)r;
  *density = *(const double*)context;
  return MYB_OK;
}

/* A wind to the east, in proportion to the distance from the axis. */
static myb_status east_wind(void* context, myb_time t, const double r[3], double w[3]) {
  (void)t;
  const double rate = *(const double*)context;
  w[0] = -rate * r[1];
  w[1] = rate * r[0];
  w[2] = 0.0;
  return MYB_OK;
}

static int orthonormal(const double m[9]) {
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double d = dot3(m + 3 * i, m + 3 * j) - (i == j ? 1.0 : 0.0);
      worst = fmax(worst, fabs(d));
    }
  }
  /* and right-handed */
  const double* x = m;
  const double* y = m + 3;
  const double* z = m + 6;
  const double det = x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0]) +
                     x[2] * (y[0] * z[1] - y[1] * z[0]);
  return worst < 1e-15 && fabs(det - 1.0) < 1e-15;
}

static void test_terms(void) {
  const double mu = 3.986004418e14;
  mover sun_motion = {.radius = 1.495978707e11, .rate = 1.99e-7, .phase = 0.3};
  spinner spin = {.rate = 7.292115e-5};
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = mu, .radius = 6378136.6,
                              .orientation = spinner_orientation, .orientation_context = &spin};
  const myb_env_body sun = {.id = 10, .gm = 1.32712440018e20, .state = mover_state,
                            .state_context = &sun_motion, .radius = 6.957e8,
                            .luminosity = 3.828e26};
  size_t i_earth = 0;
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &earth, &i_earth) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  const double t = 1234.5;
  const double x[3] = {-2.0e6, 6.5e6, 1.5e6};
  const double v[3] = {-7.0e3, -2.0e3, 1.5e3};
  double sun_state[6];
  MT_CHECK(myb_env_state(env, t, i_sun, sun_state) == MYB_OK);

  /* attitude laws */
  double m[9];
  myb_orbital_attitude orbital = {.env = env, .body = i_earth};
  MT_CHECK(myb_attitude_orbital(&orbital, t, x, v, m) == MYB_OK);
  MT_CHECK(orthonormal(m));
  MT_CHECK(fabs(dot3(m, x) - norm3(x)) < 1e-9);            /* x from the body */
  MT_CHECK(fabs(dot3(m + 6, x)) < 1e-9 && fabs(dot3(m + 6, v)) < 1e-12); /* z across the orbit */
  MT_CHECK(dot3(m + 3, v) > 0.0);                          /* y towards the motion */
  myb_sun_attitude yaw = {.env = env, .body = i_earth, .source = i_sun};
  MT_CHECK(myb_attitude_sun_nadir(&yaw, t, x, v, m) == MYB_OK);
  MT_CHECK(orthonormal(m));
  const double to_sun[3] = {sun_state[0] - x[0], sun_state[1] - x[1], sun_state[2] - x[2]};
  MT_CHECK(fabs(dot3(m + 6, x) + norm3(x)) < 1e-9);        /* z to the body */
  MT_CHECK(fabs(dot3(m + 3, to_sun)) < 1e-4);              /* y across the Sun */
  MT_CHECK(dot3(m, to_sun) > 0.0);                         /* x on the side of the Sun */
  const double in_line[3] = {-sun_state[0] * 1e-4, -sun_state[1] * 1e-4, -sun_state[2] * 1e-4};
  MT_CHECK(myb_attitude_sun_nadir(&yaw, t, in_line, v, m) == MYB_E_DOMAIN);
  myb_spin_attitude turning = {.m0 = {0, 1, 0, 0, 0, 1, 1, 0, 0}, .rate = 0.01, .phase = 0.2};
  MT_CHECK(myb_attitude_spin(&turning, 30.0, x, v, m) == MYB_OK);
  MT_CHECK(orthonormal(m));
  MT_CHECK(m[6] == 1.0 && m[7] == 0.0 && m[8] == 0.0);     /* the axis stays */
  MT_CHECK(fabs(m[1] - cos(0.5)) < 1e-16 && fabs(m[2] - sin(0.5)) < 1e-16);
  double fixed[9] = {0, 0, 1, 1, 0, 0, 0, 1, 0};
  MT_CHECK(myb_attitude_fixed(fixed, t, x, v, m) == MYB_OK && memcmp(m, fixed, sizeof m) == 0);

  /* a black plate facing the Sun is the sphere with cr = 1 */
  myb_object object = {.mass = 400.0, .area = 6.0, .cr = 1.0, .cd = 2.2};
  myb_radiation_term sphere = {.env = env, .source = i_sun, .object = &object, .c = c_light,
                               .occulters = {i_earth}, .occulter_count = 1};
  myb_plate panel = {.area = 6.0, .normal = {1.0, 0.0, 0.0}};
  /* the body x axis towards the Sun: rows x, y, z */
  const double sn = norm3(to_sun);
  double facing[9];
  for (int i = 0; i < 3; ++i) {
    facing[i] = to_sun[i] / sn;
  }
  const double helper[3] = {0.0, 0.0, 1.0};
  double yv[3] = {helper[1] * facing[2] - helper[2] * facing[1],
                  helper[2] * facing[0] - helper[0] * facing[2],
                  helper[0] * facing[1] - helper[1] * facing[0]};
  const double yn = norm3(yv);
  for (int i = 0; i < 3; ++i) {
    facing[3 + i] = yv[i] / yn;
  }
  facing[6] = facing[1] * facing[5] - facing[2] * facing[4];
  facing[7] = facing[2] * facing[3] - facing[0] * facing[5];
  facing[8] = facing[0] * facing[4] - facing[1] * facing[3];
  myb_surface_radiation_term plate_term = {.env = env, .source = i_sun, .object = &object,
                                           .c = c_light, .occulters = {i_earth},
                                           .occulter_count = 1, .attitude = myb_attitude_fixed,
                                           .attitude_context = facing, .plates = &panel,
                                           .plate_count = 1};
  double a_plate[3];
  double a_sphere[3];
  MT_CHECK(myb_term_surface_radiation(&plate_term, t, x, v, a_plate) == MYB_OK);
  MT_CHECK(myb_term_radiation(&sphere, t, x, v, a_sphere) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(a_plate[i] - a_sphere[i]) / norm3(a_sphere));
  }
  printf("a black plate facing the Sun against the sphere: %.1e\n", worst);
  MT_CHECK(worst < 1e-15 && norm3(a_sphere) > 0.0);
  /* edge-on: nothing; in the shadow of the Earth: nothing */
  double edge[9];
  memcpy(edge, facing + 3, 3 * sizeof(double));
  memcpy(edge + 3, facing + 6, 3 * sizeof(double));
  memcpy(edge + 6, facing, 3 * sizeof(double));
  plate_term.attitude_context = edge;
  MT_CHECK(myb_term_surface_radiation(&plate_term, t, x, v, a_plate) == MYB_OK);
  MT_CHECK(norm3(a_plate) < 1e-20);
  plate_term.attitude_context = facing;
  const double dark[3] = {-sun_state[0] * 6.0e-5, -sun_state[1] * 6.0e-5, 0.0};
  MT_CHECK(myb_term_surface_radiation(&plate_term, t, dark, v, a_plate) == MYB_OK);
  MT_CHECK(norm3(a_plate) == 0.0);
  /* a mesh in place of the plates: the cube, black: its silhouette is its face */
  static double vertices[3 * 8];
  static size_t faces[3 * 12];
  size_t vertex_count = 0;
  size_t face_count = 0;
  const double lo[3] = {-0.5, -0.5, -0.5};
  const double hi[3] = {0.5, 0.5, 0.5};
  add_box(vertices, &vertex_count, faces, &face_count, lo, hi);
  myb_polyhedron* cube = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, vertex_count, faces, face_count, 1.0, NULL, &cube) ==
           MYB_OK);
  const myb_plate black = {0};
  myb_surface_radiation_term mesh_term = plate_term;
  mesh_term.plates = NULL;
  mesh_term.plate_count = 0;
  mesh_term.shape = cube;
  mesh_term.surface = &black;
  mesh_term.shadows = 2;
  object.area = 1.0;
  double a_mesh[3];
  MT_CHECK(myb_term_surface_radiation(&mesh_term, t, x, v, a_mesh) == MYB_OK);
  MT_CHECK(myb_term_radiation(&sphere, t, x, v, a_sphere) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(a_mesh[i] - a_sphere[i]) < 1e-14 * norm3(a_sphere));
  }
  myb_polyhedron_destroy(cube);

  /* plates over a sphere in a gas against the drag of a sphere with its coefficient */
  const double radius = 0.9;
  const double pole[3] = {0.0, 0.0, 1.0};
  make_sphere(radius, pole); /* not along the flow: the quadrature is the limit of this check */
  double density = 2.5e-12;
  const myb_gas gas = {.temperature = 1000.0, .molar_mass = 0.017, .wall_temperature = 300.0,
                       .accommodation = 0.9};
  myb_surface_drag_term gas_term = {.env = env, .body = i_earth, .density = constant_density,
                                    .density_context = &density, .object = &object, .gas = gas,
                                    .attitude = myb_attitude_orbital,
                                    .attitude_context = &orbital, .plates = sphere_plates,
                                    .plate_count = SPHERE_PLATES};
  myb_drag_term ball = {.env = env, .body = i_earth, .density = constant_density,
                        .density_context = &density, .object = &object};
  /* the speed relative to the turning air */
  const double omega[3] = {0.0, 0.0, spin.rate};
  const double relative[3] = {v[0] - (omega[1] * x[2] - omega[2] * x[1]),
                              v[1] - (omega[2] * x[0] - omega[0] * x[2]),
                              v[2] - (omega[0] * x[1] - omega[1] * x[0])};
  double cd = 0.0;
  MT_CHECK(myb_gas_sphere_cd(&gas, norm3(relative), &cd) == MYB_OK);
  object.cd = cd;
  object.area = pi * radius * radius;
  double a_gas[3];
  double a_ball[3];
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_gas) == MYB_OK);
  MT_CHECK(myb_term_drag(&ball, t, x, v, a_ball) == MYB_OK);
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(a_gas[i] - a_ball[i]) / norm3(a_ball));
  }
  printf("plates over a sphere in the gas against the drag of the sphere (C_D %.4f): %.1e\n", cd,
         worst);
  MT_CHECK(worst < 1e-6);
  /* a wind that carries the air round a tenth faster: the same as a faster rotation */
  double rate = 0.1 * spin.rate;
  ball.wind = east_wind;
  ball.wind_context = &rate;
  gas_term.wind = east_wind;
  gas_term.wind_context = &rate;
  double a_wind[3];
  MT_CHECK(myb_term_drag(&ball, t, x, v, a_wind) == MYB_OK);
  const double faster[3] = {v[0] + 1.1 * spin.rate * x[1], v[1] - 1.1 * spin.rate * x[0], v[2]};
  const double speed = norm3(faster);
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double expected = -0.5 * density * cd * (object.area / object.mass) * speed * faster[i];
    worst = fmax(worst, fabs(a_wind[i] - expected) / norm3(a_ball));
  }
  MT_CHECK(worst < 1e-14);
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_gas) == MYB_OK);
  MT_CHECK(myb_gas_sphere_cd(&gas, speed, &cd) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    const double expected = -0.5 * density * cd * (object.area / object.mass) * speed * faster[i];
    MT_CHECK(fabs(a_gas[i] - expected) < 1e-6 * norm3(a_ball));
  }
  /* the Jacobian of the drag with the wind against differences in the velocity */
  double dadx[9];
  double dadv[9];
  MT_CHECK(myb_term_drag_jacobian(&ball, t, x, v, dadx, dadv) == MYB_OK);
  worst = 0.0;
  double largest = 0.0;
  for (int j = 0; j < 3; ++j) {
    double q[3];
    double plus[3];
    double minus[3];
    memcpy(q, v, sizeof q);
    q[j] = v[j] + 1.0;
    MT_CHECK(myb_term_drag(&ball, t, x, q, plus) == MYB_OK);
    q[j] = v[j] - 1.0;
    MT_CHECK(myb_term_drag(&ball, t, x, q, minus) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst = fmax(worst, fabs(0.5 * (plus[i] - minus[i]) - dadv[3 * i + j]));
      largest = fmax(largest, fabs(dadv[3 * i + j]));
    }
  }
  MT_CHECK(worst < 1e-7 * largest);
  /* the transition to the continuum: in air so dense that the free path is a millimetre the
   * term is the drag given for the continuum; in thin air it is the free molecules */
  gas_term.wind = NULL;
  gas_term.length = 2.0;
  gas_term.molecule = 3.7e-10;
  gas_term.continuum_cd = 1.0;
  gas_term.continuum_area = 3.0;
  double a_thin[3];
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_thin) == MYB_OK);
  gas_term.length = 0.0;
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_gas) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(a_thin[i] == a_gas[i]); /* Knudsen 1e7: free molecules */
  }
  gas_term.length = 2.0;
  density = 1.0e-3;
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_thin) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    const double expected = -0.5 * density * 1.0 * (3.0 / object.mass) * norm3(relative) *
                            relative[i];
    MT_CHECK(fabs(a_thin[i] - expected) <= 1e-14 * fabs(expected));
  }
  gas_term.molecule = 0.0;
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_thin) == MYB_E_DOMAIN);
  gas_term.length = 0.0;
  density = 0.0; /* no gas, no force */
  MT_CHECK(myb_term_surface_drag(&gas_term, t, x, v, a_gas) == MYB_OK && norm3(a_gas) == 0.0);
  myb_env_destroy(env);
}

static void test_bridge_and_yorp(void) {
  /* the bridge between the free molecules and the continuum */
  MT_CHECK(myb_gas_bridge(1.0e-4) == 0.0 && myb_gas_bridge(1.0e-3) == 0.0);
  MT_CHECK(myb_gas_bridge(10.0) == 1.0 && myb_gas_bridge(1.0e3) == 1.0);
  MT_CHECK(fabs(myb_gas_bridge(0.1) - 0.5) < 1e-15); /* the middle of the four decades */
  MT_CHECK(myb_gas_bridge(0.01) < myb_gas_bridge(0.1) && myb_gas_bridge(1.0) > 0.5);
  /* air at 100 km (5.6e-7 kg/m^3, 28.4 g/mol, 3.7e-10 m): a free path of 14 cm */
  const double path = myb_gas_knudsen(5.6e-7, 0.0284, 3.7e-10, 1.0);
  printf("the mean free path of air at 100 km: %.3f m; a body of 2 m there: Knudsen %.3f, "
         "the free-molecular part %.3f\n",
         path, path / 2.0, myb_gas_bridge(path / 2.0));
  MT_CHECK(path > 0.10 && path < 0.20);
  MT_CHECK(myb_gas_knudsen(0.0, 0.0284, 3.7e-10, 1.0) == (double)INFINITY);

  /* YORP with the heat conducted into the surface */
  make_icosphere(3);
  static double twisted[3 * MAX_VERTICES];
  for (size_t k = 0; k < ico_vertex_count; ++k) {
    const double* v = ico_vertices + 3 * k;
    const double lon = atan2(v[1], v[0]);
    const double equator = 1.0 - v[2] * v[2];
    const double phase = fmod(lon / pi + 2.0, 1.0);
    const double height = 0.25 * equator * equator * phase * (1.0 + 0.3 * v[2]);
    for (int i = 0; i < 3; ++i) {
      twisted[3 * k + (size_t)i] = v[i] * (1.0 + height);
    }
  }
  myb_polyhedron* mill = NULL;
  MT_CHECK(myb_polyhedron_create(twisted, ico_vertex_count, ico_faces, ico_face_count, 1.0, NULL,
                                 &mill) == MYB_OK);
  const double colatitude = 1.1;
  /* without thermal inertia: the average of the torque of the re-emission, turned into the
   * axes that do not rotate — the mesh with and without re-emission, subtracted */
  const myb_plate lambert = {.thermal = 1};
  const myb_plate black = {0};
  const int phases = 96;
  double direct[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < phases; ++k) {
    const double phi = two_pi * (k + 0.5) / phases;
    const double u[3] = {-sin(colatitude) * cos(phi), sin(colatitude) * sin(phi),
                         -cos(colatitude)};
    double with[3];
    double without[3];
    MT_CHECK(myb_force_radiation_mesh(mill, &lambert, 1.0, u, 0, NULL, with) == MYB_OK);
    MT_CHECK(myb_force_radiation_mesh(mill, &black, 1.0, u, 0, NULL, without) == MYB_OK);
    const double m[3] = {with[0] - without[0], with[1] - without[1], with[2] - without[2]};
    direct[0] += (m[0] * cos(phi) - m[1] * sin(phi)) / phases;
    direct[1] += (m[0] * sin(phi) + m[1] * cos(phi)) / phases;
    direct[2] += m[2] / phases;
  }
  double cold[3];
  double warm[3];
  double hot[3];
  MT_CHECK(myb_force_yorp(mill, 1.0, colatitude, 0.0, 0, phases, cold) == MYB_OK);
  MT_CHECK(myb_force_yorp(mill, 1.0, colatitude, 1.0, 0, phases, warm) == MYB_OK);
  MT_CHECK(myb_force_yorp(mill, 1.0, colatitude, 30.0, 0, phases, hot) == MYB_OK);
  double worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(cold[i] - direct[i]));
  }
  printf("YORP (p R^3): without thermal inertia (%.4e, %.4e, %.4e), within %.1e of the mesh "
         "averaged; thermal parameter 1: (%.4e, %.4e, %.4e); 30: (%.4e, %.4e, %.4e)\n",
         cold[0], cold[1], cold[2], worst, warm[0], warm[1], warm[2], hot[0], hot[1], hot[2]);
  MT_CHECK(worst < 1e-14);
  /* along the axis: the same whatever the conduction; across it: turned and weakened */
  MT_CHECK(warm[2] == cold[2] && hot[2] == cold[2] && fabs(cold[2]) > 1e-4);
  const double across_cold = hypot(cold[0], cold[1]);
  const double across_warm = hypot(warm[0], warm[1]);
  const double across_hot = hypot(hot[0], hot[1]);
  MT_CHECK(across_cold > 1e-5);
  MT_CHECK(fabs(across_warm / across_cold - 1.0 / sqrt(2.5)) < 1e-12); /* |1 / (1.5 + 0.5 i)| */
  MT_CHECK(across_hot < 0.1 * across_cold);
  /* the lag: a quarter of the phase of 1 / (1 + (1 + i) / 2) */
  const double turned = atan2(warm[1], warm[0]) - atan2(cold[1], cold[0]);
  MT_CHECK(fabs(remainder(turned - atan2(0.5, 1.5), two_pi)) < 1e-12);
  /* a convex body casts no shadows on itself: the same with them; a dented one differs */
  make_icosphere(2);
  myb_polyhedron* ball = NULL;
  MT_CHECK(myb_polyhedron_create(ico_vertices, ico_vertex_count, ico_faces, ico_face_count, 1.0,
                                 NULL, &ball) == MYB_OK);
  double open[3];
  double shaded[3];
  MT_CHECK(myb_force_yorp(ball, 1.0, colatitude, 0.5, 0, 16, open) == MYB_OK);
  MT_CHECK(myb_force_yorp(ball, 1.0, colatitude, 0.5, 2, 16, shaded) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(open[i] - shaded[i]) < 1e-13);
  }
  MT_CHECK(myb_force_yorp(ball, 1.0, colatitude, -1.0, 0, 16, open) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_yorp(ball, 1.0, colatitude, 0.5, 0, 4, open) == MYB_E_DOMAIN);
  myb_polyhedron_destroy(ball);
  myb_polyhedron_destroy(mill);
}

/* ---- Jacobians ------------------------------------------------------------------------------- */

/* The Jacobian of a term by central differences at two steps, the larger error taken out
 * (Richardson): the reference the formulas are held against. */
static void differences(myb_term_function term, void* context, double t, const double x[3],
                        const double v[3], double hx, double hv, double dadx[9],
                        double dadv[9]) {
  for (int j = 0; j < 6; ++j) {
    double slope[2][3];
    for (int level = 0; level < 2; ++level) {
      const double h = (j < 3 ? hx : hv) / (level ? 2.0 : 1.0);
      double q[6] = {x[0], x[1], x[2], v[0], v[1], v[2]};
      double plus[3];
      double minus[3];
      q[j] += h;
      MT_CHECK(term(context, t, q, q + 3, plus) == MYB_OK);
      q[j] -= 2.0 * h;
      MT_CHECK(term(context, t, q, q + 3, minus) == MYB_OK);
      for (int i = 0; i < 3; ++i) {
        slope[level][i] = (plus[i] - minus[i]) / (2.0 * h);
      }
    }
    double* out = j < 3 ? dadx : dadv;
    for (int i = 0; i < 3; ++i) {
      out[3 * i + j % 3] = (4.0 * slope[1][i] - slope[0][i]) / 3.0;
    }
  }
}

/* The largest difference of two Jacobians over the largest element of the reference, the
 * part in the position and the part in the velocity each by its own. */
static double compare(myb_term_function term, myb_term_jacobian jacobian, void* context,
                      double t, const double x[3], const double v[3], double hx, double hv,
                      double* largest_dadv) {
  double dadx[9];
  double dadv[9];
  double ref_x[9];
  double ref_v[9];
  MT_CHECK(jacobian(context, t, x, v, dadx, dadv) == MYB_OK);
  differences(term, context, t, x, v, hx, hv, ref_x, ref_v);
  double scale_x = 0.0;
  double scale_v = 0.0;
  double worst_x = 0.0;
  double worst_v = 0.0;
  for (int i = 0; i < 9; ++i) {
    scale_x = fmax(scale_x, fabs(ref_x[i]));
    scale_v = fmax(scale_v, fabs(ref_v[i]));
    worst_x = fmax(worst_x, fabs(dadx[i] - ref_x[i]));
    worst_v = fmax(worst_v, fabs(dadv[i] - ref_v[i]));
  }
  double a[3];
  MT_CHECK(term(context, t, x, v, a) == MYB_OK);
  /* a part that is zero is zero in both to the rounding of the differences */
  const double floor_v = 1e-9 * norm3(a) / norm3(v);
  if (largest_dadv != NULL) {
    *largest_dadv = scale_v;
  }
  MT_CHECK(scale_x > 0.0);
  return fmax(worst_x / scale_x, scale_v > floor_v ? worst_v / scale_v : worst_v / floor_v * 1e-9);
}

/* A law of the caller's own: x along the velocity, z across the orbit. */
static myb_status along_track(void* context, double t, const double x[3], const double v[3],
                              double m[9]) {
  (void)context;
  (void)t;
  const double vn = norm3(v);
  const double h[3] = {x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2],
                       x[0] * v[1] - x[1] * v[0]};
  const double hn = norm3(h);
  for (int i = 0; i < 3; ++i) {
    m[i] = v[i] / vn;
    m[6 + i] = h[i] / hn;
  }
  m[3] = m[7] * m[2] - m[8] * m[1];
  m[4] = m[8] * m[0] - m[6] * m[2];
  m[5] = m[6] * m[1] - m[7] * m[0];
  return MYB_OK;
}

static void test_jacobians(void) {
  const double mu = 3.986004418e14;
  mover sun_motion = {.radius = 1.495978707e11, .rate = 1.99e-7, .phase = 0.3};
  spinner spin = {.rate = 7.292115e-5};
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = mu, .radius = 6378136.6,
                              .orientation = spinner_orientation, .orientation_context = &spin};
  const myb_env_body sun = {.id = 10, .gm = 1.32712440018e20, .state = mover_state,
                            .state_context = &sun_motion, .radius = 6.957e8,
                            .luminosity = 3.828e26};
  size_t i_earth = 0;
  size_t i_sun = 0;
  MT_CHECK(myb_env_add_body(env, &earth, &i_earth) == MYB_OK);
  MT_CHECK(myb_env_add_body(env, &sun, &i_sun) == MYB_OK);
  const double t = 1234.5;
  const double x[3] = {5.1e6, 4.2e6, 1.9e6}; /* 6.87e6 m from the centre, in the light */
  const double v[3] = {-4.4e3, 4.9e3, 3.7e3};

  /* a box of six plates of different surfaces and a panel */
  myb_plate plates[7];
  memset(plates, 0, sizeof plates);
  for (int k = 0; k < 6; ++k) {
    plates[k].area = 1.0 + 0.5 * k;
    plates[k].normal[k / 2] = k % 2 ? 1.0 : -1.0;
    plates[k].specular = 0.1 * k;
    plates[k].diffuse = 0.05 * k;
    plates[k].thermal = k % 2;
  }
  plates[6] = (myb_plate){.area = 8.0, .normal = {0.6, 0.0, 0.8}, .specular = 0.2, .diffuse = 0.1};
  myb_object object = {.mass = 400.0, .area = 6.0, .cr = 1.3, .cd = 2.2};

  /* the attitude laws: the four of the library and one of the caller's */
  double tilted[9] = {0.36, 0.48, 0.8, -0.8, 0.6, 0.0, -0.48, -0.64, 0.6};
  myb_spin_attitude turning = {.m0 = {0.36, 0.48, 0.8, -0.8, 0.6, 0.0, -0.48, -0.64, 0.6},
                               .rate = 0.01, .phase = 0.2};
  myb_orbital_attitude orbital = {.env = env, .body = i_earth};
  myb_sun_attitude yaw = {.env = env, .body = i_earth, .source = i_sun};
  const struct {
    const char* name;
    myb_attitude_function law;
    void* context;
    int by_velocity;
    double step; /* of the reference for the light: the attitude changes over the distance
                    from the planet, the light alone over that from the Sun */
  } laws[5] = {{"fixed", myb_attitude_fixed, tilted, 0, 1.0e5},
               {"spin", myb_attitude_spin, &turning, 0, 1.0e5},
               {"orbital", myb_attitude_orbital, &orbital, 1, 100.0},
               {"Sun and nadir", myb_attitude_sun_nadir, &yaw, 0, 100.0},
               {"a law of the caller", along_track, NULL, 1, 100.0}};

  /* light on the plates */
  myb_surface_radiation_term light = {.env = env, .source = i_sun, .object = &object,
                                      .c = c_light, .occulters = {i_earth}, .occulter_count = 1,
                                      .plates = plates, .plate_count = 7};
  printf("Jacobians of the surface terms against differences of the terms:\n  light on plates:");
  for (int k = 0; k < 5; ++k) {
    light.attitude = laws[k].law;
    light.attitude_context = laws[k].context;
    double by_velocity = 0.0;
    const double worst = compare(myb_term_surface_radiation, myb_term_surface_radiation_jacobian,
                                 &light, t, x, v, laws[k].step, 0.5, &by_velocity);
    printf(" %s %.1e;", laws[k].name, worst);
    MT_CHECK(worst < 2e-8);
    MT_CHECK((by_velocity > 0.0) == (laws[k].by_velocity != 0) || by_velocity < 1e-20);
  }
  printf("\n");
  /* in the full shadow nothing; in the penumbra the differences of the term itself */
  double sun_state[6];
  MT_CHECK(myb_env_state(env, t, i_sun, sun_state) == MYB_OK);
  light.attitude = myb_attitude_orbital;
  light.attitude_context = &orbital;
  const double sn = norm3(sun_state);
  const double dark[3] = {-sun_state[0] / sn * 7.0e6, -sun_state[1] / sn * 7.0e6, 1.0e5};
  double dadx[9];
  double dadv[9];
  MT_CHECK(myb_term_surface_radiation_jacobian(&light, t, dark, v, dadx, dadv) == MYB_OK);
  double sum = 0.0;
  for (int i = 0; i < 9; ++i) {
    sum += fabs(dadx[i]) + fabs(dadv[i]);
  }
  MT_CHECK(sum == 0.0);
  /* across the axis of the shadow: its edge is 6.378e6 m from the axis there, 30 km wide */
  double across[3] = {-sun_state[1] / sn, sun_state[0] / sn, 0.0};
  double edge[3];
  int found = 0;
  for (int k = 0; k < 400 && !found; ++k) {
    const double off = 6.34e6 + 200.0 * k;
    for (int i = 0; i < 3; ++i) {
      edge[i] = -sun_state[i] / sn * 7.0e6 + off * across[i];
    }
    double lit[3];
    double full[3];
    myb_surface_radiation_term open = light;
    open.occulter_count = 0;
    MT_CHECK(myb_term_surface_radiation(&light, t, edge, v, lit) == MYB_OK);
    MT_CHECK(myb_term_surface_radiation(&open, t, edge, v, full) == MYB_OK);
    const double part = norm3(lit) / norm3(full);
    found = part > 0.2 && part < 0.8;
  }
  MT_CHECK(found);
  if (found) {
    const double worst = compare(myb_term_surface_radiation, myb_term_surface_radiation_jacobian,
                                 &light, t, edge, v, 20.0, 0.5, NULL);
    printf("  in a penumbra (differences of the term): %.1e\n", worst);
    MT_CHECK(worst < 1e-5);
  }

  /* light on a mesh: a box of three unequal sides, with and without the ray casting */
  static double vertices[3 * 8];
  static size_t faces[3 * 12];
  size_t vertex_count = 0;
  size_t face_count = 0;
  const double lo[3] = {-0.5, -1.0, -1.5};
  const double hi[3] = {0.5, 1.0, 1.5};
  add_box(vertices, &vertex_count, faces, &face_count, lo, hi);
  myb_polyhedron* box = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, vertex_count, faces, face_count, 1.0, NULL, &box) ==
           MYB_OK);
  const myb_plate skin = {.specular = 0.3, .diffuse = 0.2, .thermal = 1};
  myb_surface_radiation_term mesh = light;
  mesh.plates = NULL;
  mesh.plate_count = 0;
  mesh.shape = box;
  mesh.surface = &skin;
  printf("  light on a mesh:");
  for (int shadows = 0; shadows <= 3; shadows += 3) {
    mesh.shadows = shadows;
    for (int k = 2; k < 5; ++k) {
      mesh.attitude = laws[k].law;
      mesh.attitude_context = laws[k].context;
      const double worst = compare(myb_term_surface_radiation,
                                   myb_term_surface_radiation_jacobian, &mesh, t, x, v, 100.0,
                                   0.5, NULL);
      printf(" %s%s %.1e;", laws[k].name, shadows ? " with rays" : "", worst);
      MT_CHECK(worst < 2e-8);
    }
  }
  printf("\n");

  /* a gas on the plates: an exponential atmosphere that turns with the planet, a wind */
  myb_exponential_atmosphere air = {.ellipsoid = {.a = 6378136.6, .f = 1.0 / 298.257},
                                    .count = 1,
                                    .base = {400000.0},
                                    .density = {3.0e-12},
                                    .scale = {58000.0}};
  double rate = 2.0e-5;
  myb_surface_drag_term drag = {.env = env, .body = i_earth, .density = myb_density_exponential,
                                .density_context = &air, .object = &object,
                                .gas = {.temperature = 900.0, .molar_mass = 0.016,
                                        .wall_temperature = 300.0, .accommodation = 0.85},
                                .plates = plates, .plate_count = 7};
  printf("  a gas on plates:");
  for (int k = 0; k < 5; ++k) {
    drag.attitude = laws[k].law;
    drag.attitude_context = laws[k].context;
    const double worst = compare(myb_term_surface_drag, myb_term_surface_drag_jacobian, &drag, t,
                                 x, v, 50.0, 0.5, NULL);
    printf(" %s %.1e;", laws[k].name, worst);
    MT_CHECK(worst < 2e-7);
  }
  drag.attitude = myb_attitude_orbital;
  drag.attitude_context = &orbital;
  drag.wind = east_wind;
  drag.wind_context = &rate;
  const double windy = compare(myb_term_surface_drag, myb_term_surface_drag_jacobian, &drag, t, x,
                               v, 50.0, 0.5, NULL);
  /* a gas dense enough for the transition to the continuum (Knudsen 0.2) */
  air.density[0] = 1.0e-7;
  drag.length = 2.0;
  drag.molecule = 3.7e-10;
  drag.continuum_cd = 1.1;
  drag.continuum_area = 6.0;
  double a[3];
  MT_CHECK(myb_term_surface_drag(&drag, t, x, v, a) == MYB_OK);
  double rho = 0.0;
  double fixed_position[3] = {x[0], x[1], x[2]};
  MT_CHECK(myb_density_exponential(&air, myb_time_from_seconds(t), fixed_position, &rho) ==
           MYB_OK);
  const double weight = myb_gas_bridge(myb_gas_knudsen(rho, 0.016, 3.7e-10, 2.0));
  const double bridged = compare(myb_term_surface_drag, myb_term_surface_drag_jacobian, &drag, t,
                                 x, v, 50.0, 0.5, NULL);
  printf(" with a wind %.1e; in the transition (weight %.2f) %.1e\n", windy, weight, bridged);
  MT_CHECK(windy < 2e-7 && bridged < 2e-7 && weight > 0.05 && weight < 0.95);
  /* the Jacobian of the sphere's drag with the same wind: its gradient is in */
  air.density[0] = 3.0e-12;
  myb_drag_term ball = {.env = env, .body = i_earth, .density = myb_density_exponential,
                        .density_context = &air, .object = &object, .wind = east_wind,
                        .wind_context = &rate};
  const double sphere = compare(myb_term_drag, myb_term_drag_jacobian, &ball, t, x, v, 50.0, 0.5,
                                NULL);
  printf("  the drag of a sphere in a wind that changes from place to place: %.1e\n", sphere);
  MT_CHECK(sphere < 2e-7);
  /* above the gas: nothing */
  const double high[3] = {4.0e7, 1.0e6, 2.0e6};
  air.count = 1;
  air.density[0] = 0.0;
  MT_CHECK(myb_term_surface_drag_jacobian(&drag, t, high, v, dadx, dadv) == MYB_OK);
  sum = 0.0;
  for (int i = 0; i < 9; ++i) {
    sum += fabs(dadx[i]) + fabs(dadv[i]);
  }
  MT_CHECK(sum == 0.0);

  /* the recoil of an antenna fixed in the body */
  myb_recoil_term recoil = {.object = &object, .c = c_light, .power = 60.0, .efficiency = 0.9,
                            .direction = {0.6, 0.0, 0.8}, .epoch = 0.0, .decay = 5000.0};
  printf("  the recoil of an antenna:");
  for (int k = 0; k < 5; ++k) {
    recoil.attitude = laws[k].law;
    recoil.attitude_context = laws[k].context;
    double by_velocity = 0.0;
    double ref_x[9];
    double ref_v[9];
    MT_CHECK(myb_term_recoil_jacobian(&recoil, t, x, v, dadx, dadv) == MYB_OK);
    differences(myb_term_recoil, &recoil, t, x, v, 100.0, 0.5, ref_x, ref_v);
    double scale = 0.0;
    double worst = 0.0;
    MT_CHECK(myb_term_recoil(&recoil, t, x, v, a) == MYB_OK);
    for (int i = 0; i < 9; ++i) {
      scale = fmax(scale, fmax(fabs(ref_x[i]) * norm3(x), fabs(ref_v[i]) * norm3(v)));
      by_velocity = fmax(by_velocity, fabs(dadv[i]));
    }
    scale = fmax(scale, 1e-6 * norm3(a)); /* a fixed attitude: zero against the rounding */
    for (int i = 0; i < 9; ++i) {
      worst = fmax(worst, fabs(dadx[i] - ref_x[i]) * norm3(x) / scale);
      worst = fmax(worst, fabs(dadv[i] - ref_v[i]) * norm3(v) / scale);
    }
    printf(" %s %.1e;", laws[k].name, worst);
    MT_CHECK(worst < 2e-8);
    MT_CHECK((by_velocity > 0.0) == (laws[k].by_velocity != 0));
  }
  printf("\n");
  myb_polyhedron_destroy(box);
  myb_env_destroy(env);
}

/* The light of the moment with the heat conducted into the surface: with theta = 0 the light
 * of the moment; over a turn the torque of the late re-emission averages to that of
 * myb_force_yorp; the late answer turns the force forward in the sense of the rotation. */
static void test_lagged_light(void) {
  myb_polyhedron* shape = NULL;
  /* a lumpy body: an octahedron stretched unevenly */
  const double vertices[18] = {1.3, 0.1, 0.0, -0.9, 0.0, 0.2, 0.0, 1.1, 0.1,
                               0.1, -0.8, 0.0, 0.0, 0.2, 1.2, -0.1, 0.0, -0.7};
  const size_t faces[24] = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
  MT_CHECK(myb_polyhedron_create(vertices, 6, faces, 8, 1.0, NULL, &shape) == MYB_OK);
  if (shape == NULL) {
    return;
  }
  const myb_plate surface = {.specular = 0.1, .diffuse = 0.2, .thermal = 1};
  const double axis[3] = {0.0, 0.0, 1.0};
  const double pressure = 4.56e-6;
  double worst_same = 0.0;
  for (int shadows = 0; shadows <= 3; shadows += 3) {
    const double u0[3] = {-0.6, 0.48, -0.64};
    double force_a[3];
    double torque_a[3];
    double force_b[3];
    double torque_b[3];
    MT_CHECK(myb_force_radiation_mesh(shape, &surface, pressure, u0, shadows, force_a,
                                      torque_a) == MYB_OK);
    MT_CHECK(myb_force_radiation_mesh_lagged(shape, &surface, pressure, u0, axis, 0.0, 64,
                                             shadows, force_b, torque_b) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst_same = fmax(worst_same, fabs(force_a[i] - force_b[i]) / pressure);
      worst_same = fmax(worst_same, fabs(torque_a[i] - torque_b[i]) / pressure);
    }
  }
  /* the mean over a turn of the late re-emission against myb_force_yorp */
  const int phases = 48;
  const double colatitude = 1.1;
  double worst_mean = 0.0;
  for (int k = 0; k < 3; ++k) {
    const double theta = k == 0 ? 0.0 : (k == 1 ? 0.7 : 3.0);
    const myb_plate emitting = {.thermal = 1};
    const myb_plate dark = {.thermal = 0};
    double mean[3] = {0.0, 0.0, 0.0};
    for (int step = 0; step < phases; ++step) {
      const double phi = 6.283185307179586 * ((double)step + 0.5) / phases;
      const double to_source[3] = {sin(colatitude) * cos(phi), -sin(colatitude) * sin(phi),
                                   cos(colatitude)};
      const double u[3] = {-to_source[0], -to_source[1], -to_source[2]};
      double f_all[3];
      double t_all[3];
      double f_dark[3];
      double t_dark[3];
      MT_CHECK(myb_force_radiation_mesh_lagged(shape, &emitting, pressure, u, axis, theta,
                                               phases, 0, f_all, t_all) == MYB_OK);
      MT_CHECK(myb_force_radiation_mesh(shape, &dark, pressure, u, 0, f_dark, t_dark) == MYB_OK);
      /* the re-emitted part, into the inertial axes (the body has turned by phi) */
      const double tb[3] = {t_all[0] - t_dark[0], t_all[1] - t_dark[1], t_all[2] - t_dark[2]};
      mean[0] += (tb[0] * cos(phi) - tb[1] * sin(phi)) / phases;
      mean[1] += (tb[0] * sin(phi) + tb[1] * cos(phi)) / phases;
      mean[2] += tb[2] / phases;
    }
    double yorp[3];
    MT_CHECK(myb_force_yorp(shape, pressure, colatitude, theta, 0, phases, yorp) == MYB_OK);
    double scale = 0.0;
    for (int i = 0; i < 3; ++i) {
      scale = fmax(scale, fabs(yorp[i]));
    }
    for (int i = 0; i < 3; ++i) {
      worst_mean = fmax(worst_mean, fabs(mean[i] - yorp[i]) / scale);
    }
  }
  printf("late re-emission: with theta = 0 the light of the moment to %.1e of the pressure; its "
         "mean over a turn against myb_force_yorp (theta 0, 0.7, 3) %.1e\n",
         worst_same, worst_mean);
  MT_CHECK(worst_same < 4e-15 && worst_mean < 1e-12);
  /* refused */
  const double bent[3] = {0.0, 0.6, 0.6};
  const double u1[3] = {0.0, 0.0, -1.0};
  double f[3];
  double tq[3];
  MT_CHECK(myb_force_radiation_mesh_lagged(shape, &surface, pressure, u1, bent, 1.0, 64, 0, f,
                                           tq) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_radiation_mesh_lagged(shape, &surface, pressure, u1, axis, 1.0, 257, 0, f,
                                           tq) == MYB_E_DOMAIN);
  myb_polyhedron_destroy(shape);
}

int main(void) {
  test_lagged_light();
  setvbuf(stdout, NULL, _IONBF, 0);
  test_bridge_and_yorp();
  test_gas();
  test_mesh_light();
  test_terms();
  test_jacobians();
  return mt_finish("test_surface");
}
