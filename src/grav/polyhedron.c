/* SPDX-License-Identifier: BSD-3-Clause */
/* Exact gravity of a homogeneous polyhedron (R. A. Werner, D. J. Scheeres, Celest. Mech. Dyn.
 * Astron. 65, 313 (1997), eq. 13-16):
 *   U       =  G rho / 2 sum_e r_e.E_e.r_e L_e  -  G rho / 2 sum_f r_f.F_f.r_f w_f
 *   grad U  = -G rho sum_e E_e r_e L_e          +  G rho sum_f F_f r_f w_f
 *   lap U   = -G rho sum_f w_f
 * with r_e, r_f from the field point to a point of the edge or face, F_f = n_f n_f^T,
 * E_e = n_A n_12^A^T + n_B n_21^B^T (face normals and outward edge normals of the two faces of
 * edge e), L_e = ln((a + b + e) / (a + b - e)) and the signed solid angle w_f of face f. */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/grav.h"
#include "num/vec3.h"

typedef struct poly_edge {
  size_t v0;
  size_t v1;
  double e[9]; /* the edge dyad */
} poly_edge;

struct myb_polyhedron {
  myb_allocator allocator;
  size_t vertex_count;
  size_t face_count;
  size_t edge_count;
  double* vertices;
  size_t* faces;
  double* normals; /* 3 per face */
  poly_edge* edges;
  double g_rho;
  double volume;
  double center[3];
};

typedef struct directed {
  size_t lo;
  size_t hi;
  size_t from; /* the edge goes from -> to in its face */
  size_t face;
} directed;

static int directed_before(const directed* a, const directed* b) {
  return a->lo < b->lo || (a->lo == b->lo && a->hi < b->hi);
}

/* Bottom-up merge sort of the directed edges by their undirected key. */
static void sort_directed(directed* d, directed* scratch, size_t n) {
  /* bound: the run width doubles, at most 64 passes */
  for (size_t width = 1; width < n; width *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * width) {
      const size_t mid = lo + width < n ? lo + width : n;
      const size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
      size_t i = lo;
      size_t j = mid;
      for (size_t k = lo; k < hi; ++k) {
        if (i < mid && (j >= hi || !directed_before(&d[j], &d[i]))) {
          scratch[k] = d[i++];
        } else {
          scratch[k] = d[j++];
        }
      }
    }
    memcpy(d, scratch, n * sizeof *d);
  }
}

static void sub3(const double* a, const double* b, double* out) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

/* Outward normal of the edge from -> to of a face with normal n: (to - from) x n, unit. */
static void edge_normal(const double* from, const double* to, const double* n, double* out) {
  double d[3];
  sub3(to, from, d);
  myb_vec3_cross(d, n, out);
  const double k = myb_vec3_norm(out);
  for (int i = 0; i < 3; ++i) {
    out[i] /= k;
  }
}

void myb_polyhedron_destroy(myb_polyhedron* polyhedron) {
  if (polyhedron == NULL) {
    return;
  }
  const myb_allocator a = polyhedron->allocator;
  myb_internal_release(&a, polyhedron->vertices);
  myb_internal_release(&a, polyhedron->faces);
  myb_internal_release(&a, polyhedron->normals);
  myb_internal_release(&a, polyhedron->edges);
  myb_internal_release(&a, polyhedron);
}

myb_status myb_polyhedron_create(const double* vertices, size_t vertex_count, const size_t* faces,
                                 size_t face_count, double g_rho, const myb_allocator* allocator,
                                 myb_polyhedron** polyhedron) {
  MYB_CHECK(vertices != NULL && faces != NULL && polyhedron != NULL);
  if (vertex_count < 4 || face_count < 4 || !(g_rho >= 0.0) || !isfinite(g_rho) ||
      face_count > SIZE_MAX / 6) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < 3 * vertex_count; ++i) {
    if (!isfinite(vertices[i])) {
      return MYB_E_DOMAIN;
    }
  }
  for (size_t i = 0; i < 3 * face_count; ++i) {
    if (faces[i] >= vertex_count) {
      return MYB_E_DOMAIN;
    }
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_polyhedron), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_polyhedron* p = block;
  memset(p, 0, sizeof *p);
  p->allocator = a;
  p->vertex_count = vertex_count;
  p->face_count = face_count;
  p->g_rho = g_rho;
  directed* d = NULL;
  const size_t directed_count = 3 * face_count;
  status = myb_internal_alloc_array(&a, 3 * vertex_count, sizeof(double), &block);
  if (status == MYB_OK) {
    p->vertices = block;
    memcpy(p->vertices, vertices, 3 * vertex_count * sizeof(double));
    status = myb_internal_alloc_array(&a, 3 * face_count, sizeof(size_t), &block);
  }
  if (status == MYB_OK) {
    p->faces = block;
    memcpy(p->faces, faces, 3 * face_count * sizeof(size_t));
    status = myb_internal_alloc_array(&a, 3 * face_count, sizeof(double), &block);
  }
  if (status == MYB_OK) {
    p->normals = block;
    status = myb_internal_alloc_array(&a, 2 * directed_count, sizeof(directed), &block);
  }
  if (status != MYB_OK) {
    myb_polyhedron_destroy(p);
    return status;
  }
  d = block;
  /* face normals, volume and centre (tetrahedra with the origin) */
  double volume6 = 0.0;
  double moment[3] = {0.0, 0.0, 0.0};
  for (size_t f = 0; f < face_count; ++f) {
    const double* v0 = p->vertices + 3 * faces[3 * f];
    const double* v1 = p->vertices + 3 * faces[3 * f + 1];
    const double* v2 = p->vertices + 3 * faces[3 * f + 2];
    double e1[3];
    double e2[3];
    double n[3];
    sub3(v1, v0, e1);
    sub3(v2, v0, e2);
    myb_vec3_cross(e1, e2, n);
    const double k = myb_vec3_norm(n);
    if (!(k > 0.0)) {
      status = MYB_E_DOMAIN; /* degenerate face */
      break;
    }
    for (int i = 0; i < 3; ++i) {
      p->normals[3 * f + (size_t)i] = n[i] / k;
    }
    double c[3];
    myb_vec3_cross(v1, v2, c);
    const double t = myb_vec3_dot(v0, c);
    volume6 += t;
    for (int i = 0; i < 3; ++i) {
      moment[i] += t * (v0[i] + v1[i] + v2[i]) / 4.0;
    }
    for (int j = 0; j < 3; ++j) {
      const size_t from = faces[3 * f + (size_t)j];
      const size_t to = faces[3 * f + (size_t)((j + 1) % 3)];
      d[3 * f + (size_t)j] = (directed){.lo = from < to ? from : to,
                                        .hi = from < to ? to : from,
                                        .from = from,
                                        .face = f};
    }
  }
  if (status == MYB_OK && !(volume6 > 0.0)) {
    status = MYB_E_DOMAIN; /* inward orientation */
  }
  if (status == MYB_OK) {
    p->volume = volume6 / 6.0;
    for (int i = 0; i < 3; ++i) {
      p->center[i] = moment[i] / volume6;
    }
    sort_directed(d, d + directed_count, directed_count);
    /* every undirected edge exactly twice, in opposite directions */
    if (directed_count % 2 != 0) {
      status = MYB_E_DOMAIN;
    }
    for (size_t k = 0; k + 1 < directed_count && status == MYB_OK; k += 2) {
      if (d[k].lo != d[k + 1].lo || d[k].hi != d[k + 1].hi || d[k].from == d[k + 1].from ||
          (k + 2 < directed_count && d[k + 2].lo == d[k].lo && d[k + 2].hi == d[k].hi) ||
          d[k].lo == d[k].hi) {
        status = MYB_E_DOMAIN;
      }
    }
  }
  if (status == MYB_OK) {
    p->edge_count = directed_count / 2;
    status = myb_internal_alloc_array(&a, p->edge_count, sizeof(poly_edge), &block);
    if (status == MYB_OK) {
      p->edges = block;
    }
  }
  if (status == MYB_OK) {
    for (size_t k = 0; k < p->edge_count; ++k) {
      const directed* x = &d[2 * k];
      const directed* y = &d[2 * k + 1];
      poly_edge* e = &p->edges[k];
      e->v0 = x->lo;
      e->v1 = x->hi;
      const double* na = p->normals + 3 * x->face;
      const double* nb = p->normals + 3 * y->face;
      const size_t xto = x->from == x->lo ? x->hi : x->lo;
      const size_t yto = y->from == y->lo ? y->hi : y->lo;
      double ea[3];
      double eb[3];
      edge_normal(p->vertices + 3 * x->from, p->vertices + 3 * xto, na, ea);
      edge_normal(p->vertices + 3 * y->from, p->vertices + 3 * yto, nb, eb);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          e->e[3 * i + j] = na[i] * ea[j] + nb[i] * eb[j];
        }
      }
    }
  }
  myb_internal_release(&a, d);
  if (status != MYB_OK) {
    myb_polyhedron_destroy(p);
    return status;
  }
  *polyhedron = p;
  return MYB_OK;
}

double myb_polyhedron_volume(const myb_polyhedron* polyhedron, double center[3]) {
  MYB_CHECK(polyhedron != NULL);
  if (center != NULL) {
    memcpy(center, polyhedron->center, sizeof polyhedron->center);
  }
  return polyhedron->volume;
}

myb_status myb_polyhedron_eval(const myb_polyhedron* polyhedron, const double r[3],
                               double* potential, double a[3], double* laplacian) {
  MYB_CHECK(polyhedron != NULL && r != NULL);
  const myb_polyhedron* p = polyhedron;
  if (!myb_vec3_finite(r)) {
    return MYB_E_DOMAIN;
  }
  double u = 0.0;
  double g[3] = {0.0, 0.0, 0.0};
  double omega_sum = 0.0;
  /* faces */
  for (size_t f = 0; f < p->face_count; ++f) {
    double r1[3];
    double r2[3];
    double r3[3];
    sub3(p->vertices + 3 * p->faces[3 * f], r, r1);
    sub3(p->vertices + 3 * p->faces[3 * f + 1], r, r2);
    sub3(p->vertices + 3 * p->faces[3 * f + 2], r, r3);
    const double l1 = myb_vec3_norm(r1);
    const double l2 = myb_vec3_norm(r2);
    const double l3 = myb_vec3_norm(r3);
    double c[3];
    myb_vec3_cross(r2, r3, c);
    const double num = myb_vec3_dot(r1, c);
    const double den = l1 * l2 * l3 + l1 * myb_vec3_dot(r2, r3) + l2 * myb_vec3_dot(r3, r1) +
                       l3 * myb_vec3_dot(r1, r2);
    const double w = 2.0 * myb_atan2(num, den);
    const double* n = p->normals + 3 * f;
    const double nr = myb_vec3_dot(n, r1); /* n.r_f */
    u -= nr * nr * w;
    for (int i = 0; i < 3; ++i) {
      g[i] += n[i] * nr * w;
    }
    omega_sum += w;
  }
  /* edges */
  for (size_t k = 0; k < p->edge_count; ++k) {
    const poly_edge* e = &p->edges[k];
    double r1[3];
    double r2[3];
    double d[3];
    sub3(p->vertices + 3 * e->v0, r, r1);
    sub3(p->vertices + 3 * e->v1, r, r2);
    sub3(p->vertices + 3 * e->v1, p->vertices + 3 * e->v0, d);
    const double a1 = myb_vec3_norm(r1);
    const double a2 = myb_vec3_norm(r2);
    const double len = myb_vec3_norm(d);
    const double den = a1 + a2 - len;
    if (!(den > 0.0)) {
      /* on the edge: E_e r_e = 0 there (r_e runs along the edge, which E_e annihilates), and
       * the limit of the term is 0 */
      continue;
    }
    const double le = myb_log1p(2.0 * len / den);
    double er[3];
    for (int i = 0; i < 3; ++i) {
      er[i] = e->e[3 * i] * r1[0] + e->e[3 * i + 1] * r1[1] + e->e[3 * i + 2] * r1[2];
    }
    u += myb_vec3_dot(r1, er) * le;
    for (int i = 0; i < 3; ++i) {
      g[i] -= er[i] * le;
    }
  }
  if (potential != NULL) {
    *potential = 0.5 * p->g_rho * u;
  }
  if (a != NULL) {
    for (int i = 0; i < 3; ++i) {
      a[i] = p->g_rho * g[i];
    }
  }
  if (laplacian != NULL) {
    *laplacian = -p->g_rho * omega_sum;
  }
  return MYB_OK;
}
