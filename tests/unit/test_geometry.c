/* SPDX-License-Identifier: BSD-3-Clause */
/* The geometry of surfaces (myb/geometry.h):
 *  - a ray and a triaxial ellipsoid, the nearest point and the altitude against CSPICE N0067
 *    (surfpt_c, nearpt_c; tests/data/geometry/ellipsoid_cspice.ref), from outside and from
 *    inside, on five shapes;
 *  - what the results must satisfy: the point lies on the surface and on the ray, the offset
 *    of the nearest point is along the normal;
 *  - a ray and a triangular mesh against an independent solution (the planes of the faces and
 *    a point-in-triangle test by areas) on a convex and on a dented icosahedron, from outside
 *    and from inside; rays through edges and vertices;
 *  - the solid angle of a mesh: 4 pi inside, 0 outside;
 *  - the angles of illumination on simple configurations; errors. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mtest.h"
#include "myb/myb.h"

static const double pi = 3.141592653589793;

static double norm3(const double v[3]) {
  return hypot(hypot(v[0], v[1]), v[2]);
}

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static void test_ellipsoid(void) {
  FILE* f = fopen(MYB_TEST_DATA_DIR "/geometry/ellipsoid_cspice.ref", "r");
  MT_CHECK(f != NULL);
  if (f == NULL) {
    return;
  }
  char line[1024];
  int rays = 0;
  int hits = 0;
  int nears = 0;
  double worst_ray = 0.0;
  double worst_near = 0.0;
  double worst_altitude = 0.0;
  double worst_surface = 0.0;
  double worst_normal = 0.0;
  /* bound: the lines of the file */
  while (fgets(line, sizeof line, f) != NULL) {
    double v[16];
    int found = 0;
    if (sscanf(line, "RAY %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf", &v[0], &v[1],
               &v[2], &v[3], &v[4], &v[5], &v[6], &v[7], &v[8], &found, &v[9], &v[10],
               &v[11]) == 13) {
      const double size = fmax(v[0], fmax(v[1], v[2]));
      double point[3];
      double distance = 0.0;
      const int hit = myb_ray_ellipsoid(v, v + 3, v + 6, &distance, point);
      MT_CHECK(hit == found);
      ++rays;
      if (hit && found) {
        ++hits;
        /* the direction is given to 1e-16: the point moves by that times the distance */
        const double reach = fmax(size, norm3(v + 3));
        double on = 0.0;
        for (int i = 0; i < 3; ++i) {
          worst_ray = fmax(worst_ray, fabs(point[i] - v[9 + i]) / reach);
          on += (point[i] / v[i]) * (point[i] / v[i]);
        }
        worst_surface = fmax(worst_surface, fabs(on - 1.0) * size / reach);
        MT_CHECK(distance >= 0.0);
      }
    } else if (sscanf(line, "NEAR %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2],
                      &v[3], &v[4], &v[5], &v[6], &v[7], &v[8], &v[9]) == 10) {
      const double size = fmax(v[0], fmax(v[1], v[2]));
      double point[3];
      double height = 0.0;
      const myb_status status = myb_ellipsoid_nearest(v, v + 3, point, &height);
      MT_CHECK(status == MYB_OK);
      if (status != MYB_OK) {
        printf("no nearest point: %s", line);
        continue;
      }
      ++nears;
      const double reach = fmax(size, norm3(v + 3));
      double on = 0.0;
      double offset[3];
      for (int i = 0; i < 3; ++i) {
        worst_near = fmax(worst_near, fabs(point[i] - v[6 + i]) / size);
        on += (point[i] / v[i]) * (point[i] / v[i]);
        offset[i] = v[3 + i] - point[i];
      }
      worst_altitude = fmax(worst_altitude, fabs(height - v[9]) / reach);
      worst_surface = fmax(worst_surface, fabs(on - 1.0));
      /* the offset is along the normal, outwards for a positive height */
      double normal[3];
      double across[3];
      myb_ellipsoid_normal(v, point, normal);
      cross3(offset, normal, across);
      worst_normal = fmax(worst_normal, norm3(across) / reach);
      MT_CHECK(fabs(dot3(offset, normal) - height) <= 1e-9 * reach);
    }
  }
  fclose(f);
  printf("ellipsoids against CSPICE: %d rays (%d meet the surface) within %.1e of the distance, "
         "%d nearest points within %.1e of the size, altitudes %.1e; on the surface to %.1e, "
         "offsets along the normal to %.1e\n",
         rays, hits, worst_ray, nears, worst_near, worst_altitude, worst_surface, worst_normal);
  MT_CHECK(rays == 300 && hits > 100 && nears == 300);
  MT_CHECK(worst_ray < 5e-15 && worst_altitude < 5e-15);
  MT_CHECK(worst_near < 1e-9); /* the nearest point is ill-determined near the evolute */
  MT_CHECK(worst_surface < 1e-14 && worst_normal < 1e-9);

  /* special rays */
  const double radii[3] = {3.0, 2.0, 1.0};
  double point[3];
  double distance = 0.0;
  const double along_x[3] = {-1.0, 0.0, 0.0};
  const double from_x[3] = {10.0, 0.0, 0.0};
  MT_CHECK(myb_ray_ellipsoid(radii, from_x, along_x, &distance, point) == 1);
  MT_CHECK(fabs(distance - 7.0) < 1e-14 && fabs(point[0] - 3.0) < 1e-14 && point[1] == 0.0);
  const double away[3] = {1.0, 0.0, 0.0};
  MT_CHECK(myb_ray_ellipsoid(radii, from_x, away, &distance, point) == 0); /* behind the ray */
  const double centre[3] = {0.0, 0.0, 0.0};
  const double up[3] = {0.0, 0.0, 5.0};
  MT_CHECK(myb_ray_ellipsoid(radii, centre, up, &distance, point) == 1);
  MT_CHECK(distance == 1.0 && point[2] == 1.0);
  const double grazing[3] = {10.0, 2.0, 0.0}; /* tangent at (0, 2, 0) */
  MT_CHECK(myb_ray_ellipsoid(radii, grazing, along_x, &distance, point) == 1);
  MT_CHECK(fabs(distance - 10.0) < 1e-7 && fabs(point[1] - 2.0) < 1e-15);
  const double above[3] = {10.0, 2.0 + 1e-9, 0.0};
  MT_CHECK(myb_ray_ellipsoid(radii, above, along_x, NULL, NULL) == 0);
  const double zero[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_ray_ellipsoid(radii, from_x, zero, &distance, point) == 0);
  const double bad[3] = {(double)NAN, 0.0, 0.0};
  MT_CHECK(myb_ray_ellipsoid(radii, bad, along_x, &distance, point) == 0);
  /* nearest: on an axis, on the surface, the centre */
  double height = 0.0;
  MT_CHECK(myb_ellipsoid_nearest(radii, from_x, point, &height) == MYB_OK);
  MT_CHECK(point[0] == 3.0 && point[1] == 0.0 && point[2] == 0.0 && height == 7.0);
  const double on_surface[3] = {0.0, 2.0, 0.0};
  MT_CHECK(myb_ellipsoid_nearest(radii, on_surface, point, &height) == MYB_OK);
  MT_CHECK(point[1] == 2.0 && height == 0.0);
  const double below_pole[3] = {0.0, 0.0, 0.5};
  MT_CHECK(myb_ellipsoid_nearest(radii, below_pole, point, &height) == MYB_OK);
  MT_CHECK(fabs(point[2] - 1.0) < 1e-15 && fabs(height + 0.5) < 1e-15);
  /* inside the evolute the nearest points leave the plane of the long axes: from the
   * centre they are the ends of the shortest axis; the positive one is given */
  MT_CHECK(myb_ellipsoid_nearest(radii, centre, point, &height) == MYB_OK);
  MT_CHECK(point[0] == 0.0 && point[1] == 0.0 && point[2] == 1.0 && height == -1.0);
  const double deep[3] = {0.5, 0.0, 0.0};
  MT_CHECK(myb_ellipsoid_nearest(radii, deep, point, &height) == MYB_OK);
  /* x = a^2 y / (a^2 - c^2) = 9/16, z = sqrt(1 - (x/3)^2) */
  MT_CHECK(fabs(point[0] - 0.5625) < 1e-15 && point[1] == 0.0 &&
           fabs(point[2] - sqrt(1.0 - 0.1875 * 0.1875)) < 1e-15);
  MT_CHECK(fabs(height + hypot(point[0] - 0.5, point[2])) < 1e-15);
  /* nearer than the end of the axis in the plane */
  MT_CHECK(-height < 2.5);
  /* near the long end of the plane the point stays in it */
  const double shallow[3] = {2.9, 0.0, 0.0};
  MT_CHECK(myb_ellipsoid_nearest(radii, shallow, point, &height) == MYB_OK);
  MT_CHECK(fabs(point[0] - 3.0) < 1e-15 && point[2] == 0.0 && fabs(height + 0.1) < 1e-15);
  /* a prolate body, a position on its long axis: a circle of nearest points */
  const double prolate[3] = {1.0, 1.0, 4.0};
  const double on_axis[3] = {0.0, 0.0, 1.5};
  MT_CHECK(myb_ellipsoid_nearest(prolate, on_axis, point, &height) == MYB_OK);
  MT_CHECK(fabs(point[2] - 1.6) < 1e-15 && point[1] == 0.0 &&
           fabs(point[0] - sqrt(1.0 - 0.16)) < 1e-15);
  MT_CHECK(myb_ellipsoid_nearest(radii, bad, point, &height) == MYB_E_DOMAIN);
}

/* ---- Meshes --------------------------------------------------------------------------------- */

enum { VERTICES = 12, FACES = 20 };

static double mesh_vertices[3 * VERTICES];
static size_t mesh_faces[3 * FACES];

/* An icosahedron stretched along the axes; with `dent` one vertex is pulled towards the
 * centre, which makes the five faces around it a hollow. */
static myb_polyhedron* make_mesh(int dent) {
  const double g = (1.0 + sqrt(5.0)) / 2.0;
  const double v[VERTICES][3] = {{-1, g, 0}, {1, g, 0},  {-1, -g, 0}, {1, -g, 0},
                                 {0, -1, g}, {0, 1, g},  {0, -1, -g}, {0, 1, -g},
                                 {g, 0, -1}, {g, 0, 1},  {-g, 0, -1}, {-g, 0, 1}};
  static const size_t faces[3 * FACES] = {
      0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
      3, 9,  4, 3, 4, 2, 3, 2, 6, 3, 6, 8,  3, 8,  9,  4, 9, 5, 2, 4,  11, 6, 2,  10, 8, 6, 7, 9, 8, 1};
  for (int k = 0; k < VERTICES; ++k) {
    mesh_vertices[3 * k] = 1000.0 * v[k][0];
    mesh_vertices[3 * k + 1] = 700.0 * v[k][1];
    mesh_vertices[3 * k + 2] = 500.0 * v[k][2];
  }
  if (dent) {
    for (int i = 0; i < 3; ++i) {
      mesh_vertices[i] *= 0.15; /* vertex 0 */
    }
  }
  memcpy(mesh_faces, faces, sizeof faces);
  myb_polyhedron* mesh = NULL;
  MT_CHECK(myb_polyhedron_create(mesh_vertices, VERTICES, mesh_faces, FACES, 1.0, NULL, &mesh) ==
           MYB_OK);
  return mesh;
}

/* The nearest face met by a ray, by the plane of each face and the areas of the three
 * triangles the point makes with its edges. */
static int reference_ray(const double origin[3], const double d[3], double* distance) {
  int hit = 0;
  for (int f = 0; f < FACES; ++f) {
    const double* a = mesh_vertices + 3 * mesh_faces[3 * f];
    const double* b = mesh_vertices + 3 * mesh_faces[3 * f + 1];
    const double* c = mesh_vertices + 3 * mesh_faces[3 * f + 2];
    const double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    double n[3];
    cross3(ab, ac, n);
    const double along = dot3(n, d);
    if (fabs(along) < 1e-12 * norm3(n)) {
      continue;
    }
    const double oa[3] = {a[0] - origin[0], a[1] - origin[1], a[2] - origin[2]};
    const double s = dot3(n, oa) / along;
    if (s < 0.0) {
      continue;
    }
    const double p[3] = {origin[0] + s * d[0], origin[1] + s * d[1], origin[2] + s * d[2]};
    const double pa[3] = {a[0] - p[0], a[1] - p[1], a[2] - p[2]};
    const double pb[3] = {b[0] - p[0], b[1] - p[1], b[2] - p[2]};
    const double pc[3] = {c[0] - p[0], c[1] - p[1], c[2] - p[2]};
    double x1[3];
    double x2[3];
    double x3[3];
    cross3(pa, pb, x1);
    cross3(pb, pc, x2);
    cross3(pc, pa, x3);
    /* inside: the three sub-triangles turn the same way as the face */
    if (dot3(x1, n) >= 0.0 && dot3(x2, n) >= 0.0 && dot3(x3, n) >= 0.0 &&
        (!hit || s < *distance)) {
      hit = 1;
      *distance = s;
    }
  }
  return hit;
}

static uint64_t next_random(uint64_t* state) {
  *state = *state * 6364136223846793005ull + 1442695040888963407ull;
  uint64_t z = *state;
  z ^= z >> 33;
  z *= 0xff51afd7ed558ccdull;
  z ^= z >> 33;
  return z;
}

static double uniform(uint64_t* state) {
  return (double)(next_random(state) >> 11) * 0x1p-53;
}

static void test_mesh(void) {
  uint64_t state = 16;
  for (int dent = 0; dent < 2; ++dent) {
    myb_polyhedron* mesh = make_mesh(dent);
    if (mesh == NULL) {
      return;
    }
    MT_CHECK(myb_polyhedron_face_count(mesh) == FACES);
    int hits = 0;
    int misses = 0;
    int disagreements = 0;
    double worst = 0.0;
    double worst_inside = 0.0;
    double worst_outside = 0.0;
    for (int k = 0; k < 4000; ++k) {
      double origin[3];
      double d[3];
      const int inside = k % 4 == 3;
      for (int i = 0; i < 3; ++i) {
        origin[i] = (inside ? 150.0 : 6000.0) * (2.0 * uniform(&state) - 1.0);
        d[i] = (inside ? 1.0 : 0.0) * (2.0 * uniform(&state) - 1.0);
      }
      if (!inside) { /* towards a point near the body */
        for (int i = 0; i < 3; ++i) {
          d[i] = 2200.0 * (2.0 * uniform(&state) - 1.0) - origin[i];
        }
      }
      const double length = norm3(d);
      const double unit[3] = {d[0] / length, d[1] / length, d[2] / length};
      double distance = 0.0;
      double point[3];
      size_t face = 0;
      const int hit = myb_polyhedron_ray(mesh, origin, d, &distance, point, &face);
      double expected = 0.0;
      const int reference = reference_ray(origin, unit, &expected);
      if (hit != reference) {
        ++disagreements;
        continue;
      }
      if (!hit) {
        ++misses;
        continue;
      }
      ++hits;
      worst = fmax(worst, fabs(distance - expected) / (1.0 + expected));
      /* the point is on the plane of the face reported, and the surface is crossed there */
      double n[3];
      myb_polyhedron_face_normal(mesh, face, n);
      const double* v0 = mesh_vertices + 3 * mesh_faces[3 * face];
      const double off[3] = {point[0] - v0[0], point[1] - v0[1], point[2] - v0[2]};
      MT_CHECK(fabs(dot3(off, n)) < 1e-9);
      const double before[3] = {point[0] - 1e-3 * unit[0], point[1] - 1e-3 * unit[1],
                                point[2] - 1e-3 * unit[2]};
      const double after[3] = {point[0] + 1e-3 * unit[0], point[1] + 1e-3 * unit[1],
                               point[2] + 1e-3 * unit[2]};
      const double w_before = myb_polyhedron_solid_angle(mesh, before);
      const double w_after = myb_polyhedron_solid_angle(mesh, after);
      if (distance > 2e-3) {
        const double in = dot3(unit, n) < 0.0 ? w_after : w_before;
        const double out = dot3(unit, n) < 0.0 ? w_before : w_after;
        worst_inside = fmax(worst_inside, fabs(in - 4.0 * pi));
        worst_outside = fmax(worst_outside, fabs(out));
      }
    }
    printf("%s icosahedron: %d rays meet it, %d miss, %d disagree with the reference; distances "
           "within %.1e; solid angle 4 pi inside to %.1e, 0 outside to %.1e\n",
           dent ? "dented" : "convex", hits, misses, disagreements, worst, worst_inside,
           worst_outside);
    MT_CHECK(hits > 1500 && misses > 300 && disagreements == 0);
    MT_CHECK(worst < 1e-13 && worst_inside < 1e-11 && worst_outside < 1e-11);

    /* through a vertex and through the middle of an edge: no gap between the faces */
    const double far[3] = {4000.0, 3000.0, 2500.0};
    for (int v = 1; v < VERTICES; ++v) {
      const double* p = mesh_vertices + 3 * v;
      const double* q = mesh_vertices + 3 * ((v % (VERTICES - 1)) + 1);
      const double to_vertex[3] = {p[0] - far[0], p[1] - far[1], p[2] - far[2]};
      const double to_middle[3] = {0.5 * (p[0] + q[0]) - far[0], 0.5 * (p[1] + q[1]) - far[1],
                                   0.5 * (p[2] + q[2]) - far[2]};
      MT_CHECK(myb_polyhedron_ray(mesh, far, to_vertex, NULL, NULL, NULL) == 1);
      MT_CHECK(myb_polyhedron_ray(mesh, far, to_middle, NULL, NULL, NULL) == 1);
    }
    /* on a face: between the two */
    const double* a = mesh_vertices + 3 * mesh_faces[3 * 12];
    const double* b = mesh_vertices + 3 * mesh_faces[3 * 12 + 1];
    const double* c = mesh_vertices + 3 * mesh_faces[3 * 12 + 2];
    const double middle[3] = {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0,
                              (a[2] + b[2] + c[2]) / 3.0};
    const double on_face = myb_polyhedron_solid_angle(mesh, middle);
    MT_CHECK(on_face > -1e-9 && on_face < 4.0 * pi + 1e-9);
    const double zero[3] = {0.0, 0.0, 0.0};
    MT_CHECK(myb_polyhedron_ray(mesh, far, zero, NULL, NULL, NULL) == 0);
    myb_polyhedron_destroy(mesh);
  }
}

/* ---- A large mesh: the tree of boxes ------------------------------------------------------- */

enum { RINGS = 64, SECTORS = 128 };
enum { BIG_VERTICES = 2 + (RINGS - 1) * SECTORS, BIG_FACES = 2 * SECTORS * (RINGS - 1) };

static double big_vertices[3 * BIG_VERTICES];
static size_t big_faces[3 * BIG_FACES];

/* A lumpy sphere of 16 128 faces: the poles and RINGS - 1 rings of SECTORS vertices. */
static myb_polyhedron* make_big_mesh(void) {
  size_t nv = 0;
  size_t nf = 0;
  big_vertices[0] = 0.0;
  big_vertices[1] = 0.0;
  big_vertices[2] = 1000.0;
  nv = 1;
  for (int i = 1; i < RINGS; ++i) {
    const double theta = pi * (double)i / (double)RINGS;
    for (int j = 0; j < SECTORS; ++j) {
      const double phi = 2.0 * pi * (double)j / (double)SECTORS;
      const double r = 1000.0 * (1.0 + 0.25 * sin(3.0 * theta) * cos(5.0 * phi));
      big_vertices[3 * nv] = r * sin(theta) * cos(phi);
      big_vertices[3 * nv + 1] = r * sin(theta) * sin(phi);
      big_vertices[3 * nv + 2] = r * cos(theta);
      ++nv;
    }
  }
  big_vertices[3 * nv] = 0.0;
  big_vertices[3 * nv + 1] = 0.0;
  big_vertices[3 * nv + 2] = -1000.0;
  const size_t south = nv++;
  for (size_t j = 0; j < SECTORS; ++j) {
    const size_t next = (j + 1) % SECTORS;
    big_faces[3 * nf] = 0;
    big_faces[3 * nf + 1] = 1 + j;
    big_faces[3 * nf + 2] = 1 + next;
    ++nf;
    const size_t last = 1 + (size_t)(RINGS - 2) * SECTORS;
    big_faces[3 * nf] = south;
    big_faces[3 * nf + 1] = last + next;
    big_faces[3 * nf + 2] = last + j;
    ++nf;
  }
  for (size_t i = 0; i + 2 < RINGS; ++i) {
    for (size_t j = 0; j < SECTORS; ++j) {
      const size_t next = (j + 1) % SECTORS;
      const size_t a = 1 + i * SECTORS + j;
      const size_t b = 1 + i * SECTORS + next;
      const size_t c = 1 + (i + 1) * SECTORS + j;
      const size_t d = 1 + (i + 1) * SECTORS + next;
      big_faces[3 * nf] = a;
      big_faces[3 * nf + 1] = c;
      big_faces[3 * nf + 2] = d;
      ++nf;
      big_faces[3 * nf] = a;
      big_faces[3 * nf + 1] = d;
      big_faces[3 * nf + 2] = b;
      ++nf;
    }
  }
  MT_CHECK(nv == BIG_VERTICES && nf == BIG_FACES);
  myb_polyhedron* mesh = NULL;
  MT_CHECK(myb_polyhedron_create(big_vertices, BIG_VERTICES, big_faces, BIG_FACES, 1.0, NULL,
                                 &mesh) == MYB_OK);
  return mesh;
}

/* Every face in turn: the nearest met, by the plane of the face and the signs of the three
 * sub-triangles. */
static int every_face(const double origin[3], const double d[3], double* distance, size_t* face) {
  int hit = 0;
  for (size_t f = 0; f < BIG_FACES; ++f) {
    const double* a = big_vertices + 3 * big_faces[3 * f];
    const double* b = big_vertices + 3 * big_faces[3 * f + 1];
    const double* c = big_vertices + 3 * big_faces[3 * f + 2];
    const double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    double n[3];
    cross3(ab, ac, n);
    const double along = dot3(n, d);
    if (fabs(along) < 1e-12 * norm3(n)) {
      continue;
    }
    const double oa[3] = {a[0] - origin[0], a[1] - origin[1], a[2] - origin[2]};
    const double s = dot3(n, oa) / along;
    if (s < 0.0 || (hit && s >= *distance)) {
      continue;
    }
    const double p[3] = {origin[0] + s * d[0], origin[1] + s * d[1], origin[2] + s * d[2]};
    const double pa[3] = {a[0] - p[0], a[1] - p[1], a[2] - p[2]};
    const double pb[3] = {b[0] - p[0], b[1] - p[1], b[2] - p[2]};
    const double pc[3] = {c[0] - p[0], c[1] - p[1], c[2] - p[2]};
    double x1[3];
    double x2[3];
    double x3[3];
    cross3(pa, pb, x1);
    cross3(pb, pc, x2);
    cross3(pc, pa, x3);
    if (dot3(x1, n) >= 0.0 && dot3(x2, n) >= 0.0 && dot3(x3, n) >= 0.0) {
      hit = 1;
      *distance = s;
      *face = f;
    }
  }
  return hit;
}

static void test_big_mesh(void) {
  myb_polyhedron* mesh = make_big_mesh();
  if (mesh == NULL) {
    return;
  }
  enum { RAYS = 3000 };
  static double origins[3 * RAYS];
  static double directions[3 * RAYS];
  static double found[RAYS];
  static int met[RAYS];
  uint64_t state = 2026;
  for (int k = 0; k < RAYS; ++k) {
    const int inside = k % 5 == 4;
    double d[3];
    for (int i = 0; i < 3; ++i) {
      origins[3 * k + i] = (inside ? 300.0 : 5000.0) * (2.0 * uniform(&state) - 1.0);
      d[i] = inside ? 2.0 * uniform(&state) - 1.0
                    : 1600.0 * (2.0 * uniform(&state) - 1.0) - origins[3 * k + i];
    }
    const double length = norm3(d);
    for (int i = 0; i < 3; ++i) {
      directions[3 * k + i] = d[i] / length;
    }
  }
  /* along the axes too: the boxes divide by zero there */
  for (int k = 0; k < 6; ++k) {
    for (int i = 0; i < 3; ++i) {
      origins[3 * k + i] = i == k / 2 ? (k % 2 ? 4000.0 : -4000.0) : 37.5;
      directions[3 * k + i] = i == k / 2 ? (k % 2 ? -1.0 : 1.0) : 0.0;
    }
  }
  const clock_t t0 = clock();
  for (int k = 0; k < RAYS; ++k) {
    met[k] = myb_polyhedron_ray(mesh, origins + 3 * k, directions + 3 * k, &found[k], NULL, NULL);
  }
  const clock_t t1 = clock();
  int hits = 0;
  int disagreements = 0;
  int other_faces = 0;
  double worst = 0.0;
  for (int k = 0; k < RAYS; ++k) {
    double expected = 0.0;
    size_t expected_face = 0;
    const int reference = every_face(origins + 3 * k, directions + 3 * k, &expected, &expected_face);
    if (reference != met[k]) {
      ++disagreements;
      continue;
    }
    if (!reference) {
      continue;
    }
    ++hits;
    worst = fmax(worst, fabs(found[k] - expected) / (1.0 + expected));
    size_t face = 0;
    MT_CHECK(myb_polyhedron_ray(mesh, origins + 3 * k, directions + 3 * k, NULL, NULL, &face) ==
             1);
    other_faces += face != expected_face;
  }
  const clock_t t2 = clock();
  const double tree = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
  const double all = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
  printf("mesh of %d faces, %d rays: %d meet it, %d disagree with a pass over every face, %d "
         "name another face; distances within %.1e; %.2f us a ray through the tree, %.0f us "
         "over every face\n",
         (int)BIG_FACES, (int)RAYS, hits, disagreements, other_faces, worst,
         1e6 * tree / (double)RAYS, 1e6 * all / (double)RAYS);
  MT_CHECK(hits > 1000 && hits < RAYS - 300 && disagreements == 0 && other_faces == 0);
  MT_CHECK(worst < 1e-13);
  MT_CHECK(met[0] && met[1] && met[2] && met[3] && met[4] && met[5]);
  MT_CHECK(tree < all);
  myb_polyhedron_destroy(mesh);
}

static void test_illumination(void) {
  const double point[3] = {1.0, 0.0, 0.0};
  const double normal[3] = {3.0, 0.0, 0.0};
  const double source[3] = {1.0 + 5.0, 5.0, 0.0};   /* 45 degrees from the normal */
  const double observer[3] = {1.0 + 2.0, 0.0, -2.0}; /* 45 degrees, in another plane */
  double phase = 0.0;
  double incidence = 0.0;
  double emission = 0.0;
  MT_CHECK(myb_illumination_angles(point, normal, source, observer, &phase, &incidence,
                                   &emission) == MYB_OK);
  MT_CHECK(fabs(incidence - pi / 4.0) < 1e-15 && fabs(emission - pi / 4.0) < 1e-15);
  MT_CHECK(fabs(phase - pi / 3.0) < 1e-15); /* cos = 1/2 */
  const double below[3] = {1.0 - 1.0, 1.0, 0.0};
  MT_CHECK(myb_illumination_angles(point, normal, below, observer, NULL, &incidence, NULL) ==
           MYB_OK);
  MT_CHECK(fabs(incidence - 3.0 * pi / 4.0) < 1e-15); /* the source below the horizon */
  MT_CHECK(myb_illumination_angles(point, normal, point, observer, &phase, &incidence,
                                   &emission) == MYB_E_DOMAIN);
  const double none[3] = {0.0, 0.0, 0.0};
  MT_CHECK(myb_illumination_angles(point, none, source, observer, &phase, &incidence,
                                   &emission) == MYB_E_DOMAIN);
  /* the margin of a shadow: the disks touch where the fraction leaves 1 and reaches 0 */
  const double sun[3] = {1.5e11, 0.0, 0.0};
  const double rs = 6.957e8;
  const double re = 6.371e6;
  for (int inner = 0; inner < 2; ++inner) {
    /* an object behind the Earth, moved across the shadow: bisect the margin for its zero */
    double lo = inner ? 0.0 : 6.0e6;
    double hi = inner ? 6.3e6 : 7.0e6;
    /* bound: 100 halvings */
    for (int k = 0; k < 100; ++k) {
      const double y = 0.5 * (lo + hi);
      const double to_sun[3] = {sun[0] + 4.0e7, -y, 0.0};
      const double to_earth[3] = {4.0e7, -y, 0.0};
      if (myb_shadow_margin(to_sun, rs, to_earth, re, inner) < 0.0) {
        lo = y;
      } else {
        hi = y;
      }
    }
    const double eps = 1.0;
    const double inside_sun[3] = {sun[0] + 4.0e7, -(lo - eps), 0.0};
    const double inside_earth[3] = {4.0e7, -(lo - eps), 0.0};
    const double outside_sun[3] = {sun[0] + 4.0e7, -(lo + eps), 0.0};
    const double outside_earth[3] = {4.0e7, -(lo + eps), 0.0};
    const double f_in = myb_shadow_fraction(inside_sun, rs, inside_earth, re, NULL);
    const double f_out = myb_shadow_fraction(outside_sun, rs, outside_earth, re, NULL);
    printf("shadow 40 000 km behind the Earth, %s edge %.3f km off the axis: visible %.3e a "
           "metre inside, %.3e a metre outside\n",
           inner ? "umbra" : "penumbra", lo / 1e3, f_in, f_out);
    if (inner) {
      MT_CHECK(f_in == 0.0 && f_out > 0.0 && f_out < 1e-6);
    } else {
      MT_CHECK(f_out == 1.0 && f_in < 1.0 && f_in > 1.0 - 1e-6);
    }
  }
  const double within[3] = {1.0e6, 0.0, 0.0};
  MT_CHECK(myb_shadow_margin(sun, rs, within, re, 0) == -1.0);
  const double behind[3] = {2.0e11, 0.0, 0.0};
  MT_CHECK(myb_shadow_margin(sun, rs, behind, re, 0) == 1.0);
}

int main(void) {
  test_ellipsoid();
  test_mesh();
  test_big_mesh();
  test_illumination();
  return mt_finish("test_geometry");
}
