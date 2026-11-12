/* SPDX-License-Identifier: BSD-3-Clause */
/* Atmospheric entry (myb/surface.h):
 *  - C_p,max against Rayleigh's pitot formula written out and against the tables of NACA
 *    Report 1135; its limit 1.839 for gamma = 1.4; the isentropic branch below M = 1 meeting
 *    the other at M = 1;
 *  - the modified Newtonian flow against closed forms: a flat plate at incidence (C_N =
 *    C_p,max sin^2 alpha); plates over a sphere by a quadrature exact for cos^3 (C_D =
 *    C_p,max / 2); a sphere of facets; a pyramid, exactly C_p,max sin^2 of the inclination of
 *    its faces, and the cone it approaches; two cubes in a line, the one behind hidden from
 *    the flow, with the torque about the origin;
 *  - the bridge: the Newtonian continuum on a sphere of plates against the old continuum_cd =
 *    C_p,max / 2 on its cross-section at every Knudsen number; the air of a function serving
 *    the whole term; a mesh as the continuum; the Jacobian;
 *  - Sutton and Graves against a hand computation; a body slowed by its drag in air of
 *    constant density and burning away, integrated with its mass in the state, against the
 *    closed form m = m0 exp(-(B / 2C) (v0^2 - v^2)). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/myb.h"

static const double pi = 3.141592653589793;
static const double two_pi = 6.283185307179586;
static const double gas_constant = 8.314462618;
static const double air_molar_mass = 0.0289644; /* kg/mol, the US Standard Atmosphere 1976 */

static double norm3(const double v[3]) {
  return hypot(hypot(v[0], v[1]), v[2]);
}

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* C_p,max of Rayleigh's pitot formula as the textbooks write it. */
static double rayleigh(double mach, double g) {
  const double m2 = mach * mach;
  const double behind = pow((g + 1.0) * (g + 1.0) * m2 / (4.0 * g * m2 - 2.0 * (g - 1.0)),
                            g / (g - 1.0)) *
                        (1.0 - g + 2.0 * g * m2) / (g + 1.0);
  return 2.0 / (g * m2) * (behind - 1.0);
}

/* The Mach number of a speed in air of a temperature, as the library forms it. */
static double mach_of(double speed, const myb_air* air) {
  return speed / sqrt(air->heat_ratio * gas_constant * air->temperature / air->molar_mass);
}

static double cp_max_of(double speed, const myb_air* air) {
  double cp = 0.0;
  MT_CHECK(myb_newtonian_cp_max(mach_of(speed, air), air->heat_ratio, &cp) == MYB_OK);
  return cp;
}

static void test_cp_max(void) {
  const double g = 1.4;
  /* the formula of the library (divided through by M^2) against the formula written out */
  static const double machs[9] = {1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 100.0};
  double worst = 0.0;
  for (int k = 0; k < 9; ++k) {
    double cp = 0.0;
    MT_CHECK(myb_newtonian_cp_max(machs[k], g, &cp) == MYB_OK);
    worst = fmax(worst, fabs(cp / rayleigh(machs[k], g) - 1.0));
  }
  printf("C_p,max against Rayleigh's formula written out, M = 1 .. 100: %.1e\n", worst);
  MT_CHECK(worst < 1e-14);
  /* NACA Report 1135, table II (normal shock, gamma = 1.4): p_02 / p_1 */
  static const double table_mach[4] = {2.0, 3.0, 5.0, 10.0};
  static const double table_pitot[4] = {5.6404, 12.061, 32.653, 129.22};
  double table_worst = 0.0;
  printf("p_02 / p_inf from C_p,max (NACA 1135):");
  for (int k = 0; k < 4; ++k) {
    double cp = 0.0;
    MT_CHECK(myb_newtonian_cp_max(table_mach[k], g, &cp) == MYB_OK);
    const double pitot = 1.0 + 0.5 * g * table_mach[k] * table_mach[k] * cp;
    table_worst = fmax(table_worst, fabs(pitot / table_pitot[k] - 1.0));
    printf(" M %.0f: %.4f (%.4f), C_p,max %.4f;", table_mach[k], pitot, table_pitot[k], cp);
  }
  printf("\n");
  MT_CHECK(table_worst < 5e-5);
  /* M -> infinity: 4 / (g + 1) ((g + 1)^2 / (4 g))^(g / (g - 1)) */
  const double limit = 4.0 / (g + 1.0) * pow((g + 1.0) * (g + 1.0) / (4.0 * g), g / (g - 1.0));
  double far = 0.0;
  MT_CHECK(myb_newtonian_cp_max(1.0e8, g, &far) == MYB_OK);
  const double g3 = 5.0 / 3.0;
  const double limit3 = 4.0 / (g3 + 1.0) * pow((g3 + 1.0) * (g3 + 1.0) / (4.0 * g3), 2.5);
  double far3 = 0.0;
  MT_CHECK(myb_newtonian_cp_max(1.0e8, g3, &far3) == MYB_OK);
  printf("C_p,max as M -> infinity: %.6f for gamma 1.4 (closed form %.6f), %.6f for 5/3 (%.6f)\n",
         far, limit, far3, limit3);
  MT_CHECK(fabs(far / limit - 1.0) < 1e-14 && fabs(far3 / limit3 - 1.0) < 1e-14);
  MT_CHECK(fabs(far - 1.839) < 5e-4);
  /* below the speed of sound: the isentropic stagnation pressure; NACA 1135 table I gives
   * p / p_0 = 0.8430 at M = 0.5 */
  double half = 0.0;
  MT_CHECK(myb_newtonian_cp_max(0.5, g, &half) == MYB_OK);
  const double isentropic = 2.0 / (g * 0.25) * (pow(1.0 + 0.2 * 0.25, 3.5) - 1.0);
  MT_CHECK(fabs(half / isentropic - 1.0) < 1e-14);
  MT_CHECK(fabs(1.0 / (1.0 + 0.5 * g * 0.25 * half) - 0.8430) < 5e-5);
  /* the two meet at M = 1 */
  double below = 0.0;
  double at = 0.0;
  MT_CHECK(myb_newtonian_cp_max(nextafter(1.0, 0.0), g, &below) == MYB_OK);
  MT_CHECK(myb_newtonian_cp_max(1.0, g, &at) == MYB_OK);
  const double sonic = 2.0 / g * (pow((g + 1.0) / 2.0, g / (g - 1.0)) - 1.0);
  printf("at M = 1: %.15f from the shock, %.15f without it (closed form %.15f); M = 0.5: %.6f\n",
         at, below, sonic, half);
  MT_CHECK(fabs(at - sonic) < 1e-15 && fabs(below - sonic) < 1e-14);
  /* M -> 0: the incompressible stagnation pressure, 1 + M^2 / 4 + ... */
  double still = 0.0;
  double slow = 0.0;
  MT_CHECK(myb_newtonian_cp_max(0.0, g, &still) == MYB_OK && still == 1.0);
  MT_CHECK(myb_newtonian_cp_max(1.0e-3, g, &slow) == MYB_OK);
  MT_CHECK(fabs(slow - (1.0 + 0.25e-6 + (2.0 - g) / 24.0 * 1e-12)) < 1e-15);
  /* errors */
  double unused = 0.0;
  MT_CHECK(myb_newtonian_cp_max(5.0, 1.0, &unused) == MYB_E_DOMAIN);
  MT_CHECK(myb_newtonian_cp_max(-1.0, g, &unused) == MYB_E_DOMAIN);
  MT_CHECK(myb_newtonian_cp_max((double)NAN, g, &unused) == MYB_E_DOMAIN);
  MT_CHECK(myb_newtonian_cp_max((double)INFINITY, g, &unused) == MYB_E_DOMAIN);
  /* the heat ratio of a mixture */
  MT_CHECK(fabs(myb_air_heat_ratio(1.0, 0.0) - 5.0 / 3.0) < 1e-15);
  MT_CHECK(fabs(myb_air_heat_ratio(0.0, 3.0e19) - 1.4) < 1e-15);
  MT_CHECK(fabs(myb_air_heat_ratio(0.2, 0.8) - 6.6 / 4.6) < 1e-15);
  MT_CHECK(myb_air_heat_ratio(0.0, 0.0) == 0.0 && myb_air_heat_ratio(-1.0, 1.0) == 0.0);
}

/* ---- Plates over a sphere: Gauss-Legendre in the cosine of the angle from an axis, on each
 * hemisphere apart, equal steps round the axis (as in test_surface.c). With the axis along the
 * flow the Newtonian pressure is cos^2 on [0, 1] and 0 on [-1, 0]: the quadrature of the force,
 * cos^3, is exact ---------------------------------------------------------------------------- */

enum { RINGS = 24, SECTORS = 32, SPHERE_PLATES = 2 * RINGS * SECTORS };

static myb_plate sphere_plates[SPHERE_PLATES];

static void make_sphere(double radius, const double axis[3]) {
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
  const double rings = (double)RINGS;
  const double sectors = (double)SECTORS;
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
      const double cosine = half == 0 ? 0.5 * (x + 1.0) : 0.5 * (x - 1.0);
      const double sine = sqrt(1.0 - cosine * cosine);
      for (int j = 0; j < SECTORS; ++j) {
        const double lon = two_pi * (j + 0.5) / sectors;
        myb_plate* p = &sphere_plates[(2 * i + half) * SECTORS + j];
        memset(p, 0, sizeof *p);
        p->area = radius * radius * 0.5 * weight * two_pi / sectors;
        for (int q = 0; q < 3; ++q) {
          p->normal[q] = sine * cos(lon) * e1[q] + sine * sin(lon) * e2[q] + cosine * z[q];
        }
        const double n = norm3(p->normal);
        for (int q = 0; q < 3; ++q) {
          p->normal[q] /= n;
        }
      }
    }
  }
}

static void test_plates(void) {
  const myb_air air = {.temperature = 250.0, .molar_mass = air_molar_mass, .heat_ratio = 1.4};
  const double speed = 6000.0;
  const double density = 3.0e-4;
  const double mass = 200.0;
  const double q = 0.5 * density * speed * speed;
  const double cp = cp_max_of(speed, &air);
  /* a flat plate of 2 m^2 facing down, at the incidence alpha: C_N = C_p,max sin^2 alpha,
   * C_D = C_p,max sin^3 alpha, C_L = C_p,max sin^2 alpha cos alpha */
  const myb_plate flat = {.area = 2.0, .normal = {0.0, 0.0, -1.0}};
  static const double degrees[4] = {5.0, 20.0, 45.0, 90.0};
  double worst = 0.0;
  for (int k = 0; k < 4; ++k) {
    const double alpha = degrees[k] * pi / 180.0;
    const double v[3] = {speed * cos(alpha), 0.0, -speed * sin(alpha)};
    double a[3] = {0.0, 0.0, 0.0};
    MT_CHECK(myb_force_newtonian_plates(density, &air, mass, &flat, 1, v, a) == MYB_OK);
    const double s = sin(alpha);
    const double scale = q * 2.0 / mass;
    /* the drag against v, the lift across it in the plane of v and the normal */
    const double drag = -dot3(a, v) / speed / scale;
    const double lift = (a[2] * cos(alpha) + a[0] * sin(alpha)) / scale;
    worst = fmax(worst, fabs(a[2] / scale - cp * s * s) / cp);
    worst = fmax(worst, fabs(drag - cp * s * s * s) / cp);
    worst = fmax(worst, fabs(lift - cp * s * s * cos(alpha)) / cp);
    MT_CHECK(a[0] == 0.0 && a[1] == 0.0);
  }
  printf("modified Newtonian, M %.2f, C_p,max %.6f: a flat plate at 5 .. 90 deg against "
         "C_p,max sin^2 alpha: %.1e\n",
         mach_of(speed, &air), cp, worst);
  MT_CHECK(worst < 1e-15);
  /* turned away from the flow, or along it: nothing */
  const double back[3] = {speed, 0.0, speed};
  const double along[3] = {speed, 0.0, 0.0};
  double nothing[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_newtonian_plates(density, &air, mass, &flat, 1, back, nothing) == MYB_OK);
  MT_CHECK(myb_force_newtonian_plates(density, &air, mass, &flat, 1, along, nothing) == MYB_OK);
  MT_CHECK(nothing[0] == 0.0 && nothing[1] == 0.0 && nothing[2] == 0.0);

  /* plates over a sphere with the axis along the flow: C_D = C_p,max / 2, no lift */
  const double radius = 0.8;
  const double v[3] = {0.48 * speed, -0.6 * speed, 0.64 * speed};
  make_sphere(radius, v);
  double a[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_newtonian_plates(density, &air, mass, sphere_plates, SPHERE_PLATES, v, a) ==
           MYB_OK);
  const double expected = q * 0.5 * cp * pi * radius * radius / mass;
  const double along_flow = -dot3(a, v) / speed;
  double across[3];
  for (int i = 0; i < 3; ++i) {
    across[i] = a[i] + along_flow * v[i] / speed;
  }
  printf("plates over a sphere: C_D %.15f against C_p,max / 2 = %.15f, %.1e apart, lift %.1e\n",
         along_flow / (q * pi * radius * radius / mass), 0.5 * cp,
         fabs(along_flow / expected - 1.0), norm3(across) / expected);
  MT_CHECK(fabs(along_flow / expected - 1.0) < 1e-14 && norm3(across) < 1e-14 * expected);
  /* errors */
  myb_air bad = air;
  bad.heat_ratio = 1.0;
  MT_CHECK(myb_force_newtonian_plates(density, &bad, mass, &flat, 1, v, a) == MYB_E_DOMAIN);
  bad = air;
  bad.temperature = 0.0;
  MT_CHECK(myb_force_newtonian_plates(density, &bad, mass, &flat, 1, v, a) == MYB_E_DOMAIN);
  const double still[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_newtonian_plates(density, &air, mass, &flat, 1, still, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_newtonian_plates(-1.0, &air, mass, &flat, 1, v, a) == MYB_E_DOMAIN);
  MT_CHECK(myb_force_newtonian_plates(density, &air, 0.0, &flat, 1, v, a) == MYB_E_DOMAIN);
  const myb_plate skew = {.area = 1.0, .normal = {0.0, 0.1, -1.0}};
  MT_CHECK(myb_force_newtonian_plates(density, &air, mass, &skew, 1, v, a) == MYB_E_DOMAIN);
}

/* ---- Meshes ---------------------------------------------------------------------------------- */

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

/* An icosahedron with every face split in four `levels` times, its vertices on the unit
 * sphere. */
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
  static const size_t faces[60] = {0, 11, 5, 0, 5,  1,  0,  1,  7,  0,  7, 10, 0, 10, 11,
                                   1, 5,  9, 5, 11, 4,  11, 10, 2,  10, 7, 6,  7, 1,  8,
                                   3, 9,  4, 3, 4,  2,  3,  2,  6,  3,  6, 8,  3, 8,  9,
                                   4, 9,  5, 2, 4,  11, 6,  2,  10, 8,  6, 7,  9, 8,  1};
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

/* A pyramid of `sides` faces round the x axis, its apex at x = height and its base a regular
 * polygon inscribed in the circle of the radius in the plane x = 0, closed by a fan. */
enum { MAX_SIDES = 128 };

static myb_polyhedron* make_pyramid(int sides, double radius, double height) {
  static double vertices[3 * (MAX_SIDES + 2)];
  static size_t faces[3 * 2 * MAX_SIDES];
  const size_t n = (size_t)sides;
  vertices[0] = height; /* the apex */
  vertices[1] = 0.0;
  vertices[2] = 0.0;
  vertices[3] = 0.0; /* the centre of the base */
  vertices[4] = 0.0;
  vertices[5] = 0.0;
  for (int j = 0; j < sides; ++j) {
    const double phi = two_pi * j / (double)sides;
    double* p = vertices + 3 * (2 + (size_t)j);
    p[0] = 0.0;
    p[1] = radius * cos(phi);
    p[2] = radius * sin(phi);
  }
  for (size_t j = 0; j < n; ++j) {
    const size_t here = 2 + j;
    const size_t next = 2 + (j + 1) % n;
    size_t* side = faces + 3 * j;
    side[0] = 0;
    side[1] = here;
    side[2] = next;
    size_t* base = faces + 3 * (n + j);
    base[0] = 1;
    base[1] = next;
    base[2] = here;
  }
  myb_polyhedron* shape = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, n + 2, faces, 2 * n, 1.0, NULL, &shape) == MYB_OK);
  return shape;
}

static void test_mesh(void) {
  const myb_air air = {.temperature = 230.0, .molar_mass = air_molar_mass, .heat_ratio = 1.4};
  const double speed = 7000.0;
  const double density = 1.0e-4;
  const double q = 0.5 * density * speed * speed;
  const double cp = cp_max_of(speed, &air);
  double force[3];
  double torque[3];

  /* a sphere of 1280 facets: C_D -> C_p,max / 2 as the facets shrink */
  make_icosphere(3);
  myb_polyhedron* ball = NULL;
  MT_CHECK(myb_polyhedron_create(ico_vertices, ico_vertex_count, ico_faces, ico_face_count, 1.0,
                                 NULL, &ball) == MYB_OK);
  const double v[3] = {0.36 * speed, -0.48 * speed, 0.8 * speed};
  MT_CHECK(myb_force_newtonian_mesh(ball, density, &air, v, 0, force, torque) == MYB_OK);
  const double cd = -dot3(force, v) / speed / (q * pi);
  printf("a sphere of %zu facets: C_D %.5f, %.5f of C_p,max / 2; torque %.1e of the force\n",
         ico_face_count, cd, cd / (0.5 * cp), norm3(torque) / norm3(force));
  MT_CHECK(fabs(cd / (0.5 * cp) - 1.0) < 0.01);
  MT_CHECK(norm3(torque) < 1e-3 * norm3(force));
  /* convex: no part hides another */
  double shaded[3];
  MT_CHECK(myb_force_newtonian_mesh(ball, density, &air, v, 2, shaded, NULL) == MYB_OK);
  for (int i = 0; i < 3; ++i) {
    MT_CHECK(fabs(shaded[i] - force[i]) < 1e-13 * norm3(force));
  }
  myb_polyhedron_destroy(ball);

  /* a pyramid of 64 faces apex first: exactly C_p,max sin^2 delta_f on the area of its base,
   * tan delta_f = R cos(pi / 64) / h; the cone of half-angle delta = 20 deg it approaches has
   * C_p,max sin^2 delta */
  const double delta = 20.0 * pi / 180.0;
  const double radius = 1.0;
  const double height = radius / tan(delta);
  const int sides = 64;
  myb_polyhedron* cone = make_pyramid(sides, radius, height);
  const double forward[3] = {speed, 0.0, 0.0};
  MT_CHECK(myb_force_newtonian_mesh(cone, density, &air, forward, 0, force, torque) == MYB_OK);
  const double base = 0.5 * sides * radius * radius * sin(two_pi / sides);
  const double inclination = atan(radius * cos(pi / sides) / height);
  const double pyramid = cp * sin(inclination) * sin(inclination);
  const double measured = -force[0] / (q * base);
  printf("a pyramid of 64 faces, 20 deg: C_D %.15f against C_p,max sin^2 delta_f %.15f (%.1e); "
         "the cone C_p,max sin^2 delta = %.6f\n",
         measured, pyramid, fabs(measured / pyramid - 1.0), cp * sin(delta) * sin(delta));
  MT_CHECK(fabs(measured / pyramid - 1.0) < 1e-13);
  MT_CHECK(fabs(force[1]) < 1e-13 * fabs(force[0]) && fabs(force[2]) < 1e-13 * fabs(force[0]));
  MT_CHECK(norm3(torque) < 1e-13 * fabs(force[0]) * height);
  /* the faces lean in by (1 - cos(pi / 64)) tan delta cos^2 delta: sin^2 falls by 2.1e-3 */
  const double lean = 2.0 * (1.0 - cos(pi / sides)) * cos(delta) * cos(delta);
  MT_CHECK(fabs(1.0 - measured / (cp * sin(delta) * sin(delta)) - lean) < 1e-5);
  /* base first: the base alone faces the flow, square to it */
  const double behind[3] = {-speed, 0.0, 0.0};
  MT_CHECK(myb_force_newtonian_mesh(cone, density, &air, behind, 0, force, NULL) == MYB_OK);
  MT_CHECK(fabs(force[0] - q * cp * base) < 1e-13 * q * cp * base);
  myb_polyhedron_destroy(cone);

  /* two cubes in a line along the flow, 1 m apart, off the axis by 2 m: the one behind is
   * hidden from the flow — its front face takes nothing once the rays are cast */
  static double vertices[3 * 16];
  static size_t faces[3 * 24];
  size_t vertex_count = 0;
  size_t face_count = 0;
  const double front_lo[3] = {2.0, 1.5, -0.5};
  const double front_hi[3] = {3.0, 2.5, 0.5};
  const double rear_lo[3] = {0.0, 1.5, -0.5};
  const double rear_hi[3] = {1.0, 2.5, 0.5};
  add_box(vertices, &vertex_count, faces, &face_count, front_lo, front_hi);
  add_box(vertices, &vertex_count, faces, &face_count, rear_lo, rear_hi);
  myb_polyhedron* two = NULL;
  MT_CHECK(myb_polyhedron_create(vertices, vertex_count, faces, face_count, 1.0, NULL, &two) ==
           MYB_OK);
  const double p = q * cp; /* on a face square to the flow, of 1 m^2 */
  MT_CHECK(myb_force_newtonian_mesh(two, density, &air, forward, 0, force, torque) == MYB_OK);
  MT_CHECK(fabs(force[0] + 2.0 * p) < 1e-13 * p && fabs(force[1]) < 1e-13 * p &&
           fabs(force[2]) < 1e-13 * p);
  MT_CHECK(fabs(torque[2] - 4.0 * p) < 1e-13 * p);
  MT_CHECK(myb_force_newtonian_mesh(two, density, &air, forward, 6, force, torque) == MYB_OK);
  printf("two cubes in a line, the rays cast: force %.12f, torque about z %.12f (of q C_p,max "
         "m^2, m: -1 and 2, the front face alone)\n",
         force[0] / p, torque[2] / p);
  MT_CHECK(fabs(force[0] + p) < 1e-13 * p && fabs(force[1]) < 1e-13 * p &&
           fabs(force[2]) < 1e-13 * p);
  MT_CHECK(fabs(torque[0]) < 1e-13 * p && fabs(torque[1]) < 1e-13 * p &&
           fabs(torque[2] - 2.0 * p) < 1e-13 * p);
  /* errors */
  myb_air bad = air;
  bad.molar_mass = -1.0;
  MT_CHECK(myb_force_newtonian_mesh(two, density, &bad, forward, 0, force, NULL) == MYB_E_DOMAIN);
  const double still[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_newtonian_mesh(two, density, &air, still, 0, force, NULL) == MYB_E_DOMAIN);
  myb_polyhedron_destroy(two);
}

/* ---- The bridge in the term ------------------------------------------------------------------ */

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

static myb_status constant_air(void* context, myb_time t, const double r[3], myb_air* air) {
  (void)t;
  (void)r;
  *air = *(const myb_air*)context;
  return MYB_OK;
}

/* The Jacobian of a term by central differences at two steps, the larger error taken out
 * (Richardson), as in test_surface.c. */
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

static void test_bridge(void) {
  spinner spin = {.rate = 7.292115e-5};
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = 3.986004418e14, .radius = 6378136.6,
                              .orientation = spinner_orientation, .orientation_context = &spin};
  size_t i_earth = 0;
  MT_CHECK(myb_env_add_body(env, &earth, &i_earth) == MYB_OK);
  const double t = 321.0;
  const double x[3] = {-2.0e6, 6.0e6, 1.0e6};
  const double v[3] = {-6.9e3, -2.4e3, 1.1e3};
  const double omega[3] = {0.0, 0.0, spin.rate};
  const double relative[3] = {v[0] - (omega[1] * x[2] - omega[2] * x[1]),
                              v[1] - (omega[2] * x[0] - omega[0] * x[2]),
                              v[2] - (omega[0] * x[1] - omega[1] * x[0])};
  const double speed = norm3(relative);
  /* plates over a sphere about the flow, the body axes those of the model */
  const double radius = 0.9;
  make_sphere(radius, relative);
  static double identity[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  myb_object object = {.mass = 300.0, .area = 1.0, .cr = 1.0, .cd = 2.2};
  double density = 0.0;
  const myb_gas gas = {.temperature = 220.0, .molar_mass = air_molar_mass,
                       .wall_temperature = 300.0, .accommodation = 1.0};
  const myb_air air = {.temperature = 220.0, .molar_mass = air_molar_mass, .heat_ratio = 1.4};
  const double cp = cp_max_of(speed, &air);
  myb_surface_drag_term old = {.env = env, .body = i_earth, .density = constant_density,
                               .density_context = &density, .object = &object, .gas = gas,
                               .attitude = myb_attitude_fixed, .attitude_context = identity,
                               .plates = sphere_plates, .plate_count = SPHERE_PLATES,
                               .length = 1.8, .molecule = 3.7e-10, .continuum_cd = 0.5 * cp,
                               .continuum_area = pi * radius * radius};
  const myb_continuum newton = {.heat_ratio = 1.4};
  myb_surface_drag_term term = old;
  term.continuum = &newton;
  /* the Newtonian sphere is the old number C_p,max / 2 at every Knudsen number */
  double worst = 0.0;
  double lowest = 1.0;
  double highest = 0.0;
  for (int k = 0; k <= 20; ++k) {
    density = pow(10.0, -12.0 + 0.5 * k);
    double a_old[3];
    double a_new[3];
    MT_CHECK(myb_term_surface_drag(&old, t, x, v, a_old) == MYB_OK);
    MT_CHECK(myb_term_surface_drag(&term, t, x, v, a_new) == MYB_OK);
    for (int i = 0; i < 3; ++i) {
      worst = fmax(worst, fabs(a_new[i] - a_old[i]) / norm3(a_old));
    }
    const double weight = myb_gas_bridge(myb_gas_knudsen(density, air_molar_mass, 3.7e-10, 1.8));
    lowest = fmin(lowest, weight);
    highest = fmax(highest, weight);
  }
  printf("the bridge with the Newtonian sphere against continuum_cd = C_p,max / 2 = %.6f, "
         "density 1e-12 .. 1e-2 (free-molecular weight %.0f .. %.0f): %.1e\n",
         0.5 * cp, lowest, highest, worst);
  MT_CHECK(worst < 1e-14 && lowest == 0.0 && highest == 1.0);
  /* without a length the term is the free molecules, to the bit, with or without a continuum */
  density = 1.0e-6;
  old.length = 0.0;
  term.length = 0.0;
  double a_old[3];
  double a_new[3];
  MT_CHECK(myb_term_surface_drag(&old, t, x, v, a_old) == MYB_OK);
  MT_CHECK(myb_term_surface_drag(&term, t, x, v, a_new) == MYB_OK);
  MT_CHECK(memcmp(a_old, a_new, sizeof a_old) == 0);
  old.length = 1.8;
  term.length = 1.8;

  /* the air of a function serves the whole term: the free molecules at its temperature and
   * molar mass, the Knudsen number of its molar mass, the continuum at its Mach number */
  myb_air warm = {.temperature = 300.0, .molar_mass = 0.025, .heat_ratio = 1.3};
  const myb_continuum local = {.heat_ratio = 1.4, .air = constant_air, .air_context = &warm};
  term.continuum = &local;
  density = 3.0e-7;
  MT_CHECK(myb_term_surface_drag(&term, t, x, v, a_new) == MYB_OK);
  myb_gas warm_gas = gas;
  warm_gas.temperature = 300.0;
  warm_gas.molar_mass = 0.025;
  const double weight = myb_gas_bridge(myb_gas_knudsen(density, 0.025, 3.7e-10, 1.8));
  double rare[3] = {0.0, 0.0, 0.0};
  double dense[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_force_gas_plates(density, &warm_gas, object.mass, sphere_plates, SPHERE_PLATES,
                                relative, rare) == MYB_OK);
  MT_CHECK(myb_force_newtonian_plates(density, &warm, object.mass, sphere_plates, SPHERE_PLATES,
                                      relative, dense) == MYB_OK);
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(a_new[i] - (weight * rare[i] + (1.0 - weight) * dense[i])) /
                            norm3(a_new));
  }
  printf("the air of a function (300 K, 25 g/mol, gamma 1.3; weight %.3f): %.1e from the parts "
         "computed apart\n",
         weight, worst);
  MT_CHECK(worst < 1e-14 && weight > 0.05 && weight < 0.95);
  myb_air cold = warm;
  cold.temperature = 0.0;
  const myb_continuum frozen = {.heat_ratio = 1.4, .air = constant_air, .air_context = &cold};
  term.continuum = &frozen;
  MT_CHECK(myb_term_surface_drag(&term, t, x, v, a_new) == MYB_E_DOMAIN);
  const myb_continuum monatomic_one = {.heat_ratio = 1.0};
  term.continuum = &monatomic_one;
  MT_CHECK(myb_term_surface_drag(&term, t, x, v, a_new) == MYB_E_DOMAIN);

  /* a mesh as the continuum, in air so dense that the free molecules weigh nothing */
  make_icosphere(2);
  myb_polyhedron* ball = NULL;
  MT_CHECK(myb_polyhedron_create(ico_vertices, ico_vertex_count, ico_faces, ico_face_count, 1.0,
                                 NULL, &ball) == MYB_OK);
  const myb_continuum meshed = {.heat_ratio = 1.4, .shape = ball, .shadows = 1};
  term.continuum = &meshed;
  density = 1.0e-2;
  MT_CHECK(myb_term_surface_drag(&term, t, x, v, a_new) == MYB_OK);
  double force[3];
  MT_CHECK(myb_force_newtonian_mesh(ball, density, &air, relative, 1, force, NULL) == MYB_OK);
  worst = 0.0;
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(a_new[i] - force[i] / object.mass) / norm3(a_new));
  }
  printf("a mesh as the continuum of the term: %.1e from myb_force_newtonian_mesh\n", worst);
  MT_CHECK(worst < 1e-14);
  myb_polyhedron_destroy(ball);

  /* the Jacobian with a continuum: differences of the term, in an exponential atmosphere at
   * 90 km in the transition, with a body turned along the orbit */
  myb_exponential_atmosphere layers = {.ellipsoid = {.a = 6378136.6, .f = 1.0 / 298.257},
                                       .count = 1,
                                       .base = {90000.0},
                                       .density = {3.4e-6},
                                       .scale = {5500.0}};
  myb_orbital_attitude orbital = {.env = env, .body = i_earth};
  myb_plate box[6];
  memset(box, 0, sizeof box);
  for (int k = 0; k < 6; ++k) {
    box[k].area = 1.0 + 0.5 * k;
    box[k].normal[k / 2] = k % 2 ? 1.0 : -1.0;
  }
  myb_surface_drag_term entry = {.env = env, .body = i_earth,
                                 .density = myb_density_exponential, .density_context = &layers,
                                 .object = &object, .gas = gas,
                                 .attitude = myb_attitude_orbital, .attitude_context = &orbital,
                                 .plates = box, .plate_count = 6, .length = 1.8,
                                 .molecule = 3.7e-10, .continuum = &newton};
  const double low[3] = {6468136.6 * 0.6, 6468136.6 * 0.8, 0.0}; /* 90 km over the equator */
  const double descent[3] = {-7400.0 * 0.8 - 150.0 * 0.6, 7400.0 * 0.6 - 150.0 * 0.8, 400.0};
  double dadx[9];
  double dadv[9];
  double ref_x[9];
  double ref_v[9];
  MT_CHECK(myb_term_surface_drag_jacobian(&entry, t, low, descent, dadx, dadv) == MYB_OK);
  differences(myb_term_surface_drag, &entry, t, low, descent, 50.0, 0.5, ref_x, ref_v);
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
  double rho = 0.0;
  MT_CHECK(myb_density_exponential(&layers, myb_time_from_seconds(t), low, &rho) == MYB_OK);
  printf("the Jacobian of the term with a continuum (free-molecular weight %.3f) against "
         "Richardson's differences: %.1e in x, %.1e in v\n",
         myb_gas_bridge(myb_gas_knudsen(rho, air_molar_mass, 3.7e-10, 1.8)), worst_x / scale_x,
         worst_v / scale_v);
  MT_CHECK(worst_x < 1e-4 * scale_x && worst_v < 1e-6 * scale_v && scale_x > 0.0);
  myb_env_destroy(env);
}

/* ---- Heat and ablation ----------------------------------------------------------------------- */

static void test_ablation(void) {
  /* Sutton and Graves by hand: 1.7415e-4 sqrt(3e-4 / 1) 7500^3 = 1.7415e-4 x 0.0173205 x
   * 4.21875e11 = 1.27253e6 W m^-2 */
  double q = 0.0;
  MT_CHECK(myb_heating_sutton_graves(1.7415e-4, 3.0e-4, 1.0, 7500.0, &q) == MYB_OK);
  printf("Sutton and Graves at 7.5 km/s, 3e-4 kg/m^3, a nose of 1 m: %.6e W/m^2 (by hand "
         "1.27253e6)\n",
         q);
  MT_CHECK(fabs(q / 1.27253e6 - 1.0) < 5e-6);
  /* q grows as V^3 and as 1 / sqrt(R_n) */
  double q2 = 0.0;
  MT_CHECK(myb_heating_sutton_graves(1.7415e-4, 3.0e-4, 4.0, 15000.0, &q2) == MYB_OK);
  MT_CHECK(fabs(q2 / q - 4.0) < 1e-14);
  MT_CHECK(myb_heating_sutton_graves(1.7415e-4, 3.0e-4, 0.0, 7500.0, &q) == MYB_E_DOMAIN);
  MT_CHECK(myb_heating_sutton_graves(1.7415e-4, -1.0, 1.0, 7500.0, &q) == MYB_E_DOMAIN);

  /* a body slowed by its drag in air of constant density, nothing else on it: m v' = -C v^2,
   * m' = -B v^3 with C = rho c_D A / 2, B = k sqrt(rho / R_n) A_q / Q*, so that
   * d ln m / dv = (B / C) v and m = m0 exp(-(B / 2C) (v0^2 - v^2)) */
  spinner still = {.rate = 0.0}; /* the air at rest */
  myb_env* env = NULL;
  const myb_env_options options = {0};
  MT_CHECK(myb_env_create(&options, NULL, &env) == MYB_OK);
  const myb_env_body earth = {.id = 399, .gm = 3.986004418e14, .radius = 6378136.6,
                              .orientation = spinner_orientation, .orientation_context = &still};
  size_t i_earth = 0;
  MT_CHECK(myb_env_add_body(env, &earth, &i_earth) == MYB_OK);
  double density = 1.0e-4;
  const double k = 1.7415e-4;
  const double nose = 0.13;
  const double heated = 0.05;
  const double heat = 1.0e7;
  myb_ablation ablation = {.env = env, .body = i_earth, .density = constant_density,
                           .density_context = &density, .k = k, .nose_radius = nose,
                           .area = heated, .heat = heat, .mass = 10.0};
  myb_object body = {.mass = 10.0, .area = 0.05, .cr = 1.0, .cd = 1.0,
                     .mass_function = myb_ablation_mass, .mass_context = &ablation};
  myb_drag_term drag = {.env = env, .body = i_earth, .density = constant_density,
                        .density_context = &density, .object = &body};
  myb_force_model* field = NULL;
  MT_CHECK(myb_force_model_create(NULL, &field) == MYB_OK);
  MT_CHECK(myb_force_model_add(field, "drag", myb_term_drag, myb_term_drag_jacobian, &drag,
                               NULL) == MYB_OK);
  ablation.field = field;
  const double y0[7] = {7.0e6, 0.0, 0.0, 0.0, 12000.0, 0.0, 10.0};
  /* the rate at the start by hand */
  double heating = 0.0;
  double rate = 0.0;
  MT_CHECK(myb_ablation_rate(&ablation, 0.0, y0, y0 + 3, &heating, &rate) == MYB_OK);
  const double b = k * sqrt(density / nose) * heated / heat;
  MT_CHECK(fabs(rate / (-b * 12000.0 * 12000.0 * 12000.0) - 1.0) < 1e-15);
  MT_CHECK(fabs(rate + heating * heated / heat) <= 1e-16 * fabs(rate));
  /* integrated with the mass in the state */
  const myb_dop853_options o = {.rtol = 1e-12, .atol = 1e-12};
  myb_dop853* integrator = NULL;
  MT_CHECK(myb_dop853_create(7, &o, NULL, &integrator) == MYB_OK);
  MT_CHECK(myb_dop853_start(integrator, myb_ablation_rhs, &ablation, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_dop853_integrate(integrator, 200.0) == MYB_OK);
  const double* y = myb_dop853_state(integrator);
  const double v_end = norm3(y + 3);
  const double c = 0.5 * density * 1.0 * 0.05;
  const double closed = 10.0 * exp(-0.5 * b / c * (12000.0 * 12000.0 - v_end * v_end));
  printf("a body of 10 kg at 12 km/s in air of 1e-4 kg/m^3 after 200 s: %.6f km/s, %.12f kg "
         "against the closed form %.12f kg (%.1e); %lld steps\n",
         v_end / 1000.0, y[6], closed, fabs(y[6] / closed - 1.0),
         (long long)myb_dop853_statistics(integrator).accepted);
  MT_CHECK(fabs(y[6] / closed - 1.0) < 1e-10);
  MT_CHECK(y[6] < 7.0 && y[6] > 3.0); /* a good part of it is gone */
  MT_CHECK(y[0] == 7.0e6 && y[2] == 0.0 && y[3] == 0.0 && y[5] == 0.0); /* a straight line */
  myb_dop853_destroy(integrator);
  /* without the ablation the drag of the same body would be that of 10 kg */
  double dydt[7];
  double gone[7];
  memcpy(gone, y0, sizeof gone);
  gone[6] = 0.0;
  MT_CHECK(myb_ablation_rhs(&ablation, 0.0, gone, dydt) == MYB_E_DOMAIN);
  gone[6] = 5.0;
  MT_CHECK(myb_ablation_rhs(&ablation, 0.0, gone, dydt) == MYB_OK);
  MT_CHECK(fabs(dydt[4] / (-c * 12000.0 * 12000.0 / 5.0) - 1.0) < 1e-15);
  MT_CHECK(myb_ablation_mass(&ablation, 123.0) == 5.0);
  ablation.heat = 0.0;
  MT_CHECK(myb_ablation_rate(&ablation, 0.0, y0, y0 + 3, &heating, &rate) == MYB_E_DOMAIN);
  myb_force_model_destroy(field);
  myb_env_destroy(env);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_cp_max();
  test_plates();
  test_mesh();
  test_bridge();
  test_ablation();
  return mt_finish("test_entry");
}
