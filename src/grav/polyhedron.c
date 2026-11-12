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
#include "myb/geometry.h"
#include "myb/grav.h"
#include "force/surface_internal.h"
#include "myb/surface.h"
#include "num/vec3.h"

typedef struct poly_edge {
  size_t v0;
  size_t v1;
  double e[9]; /* the edge dyad */
} poly_edge;

/* A node of the tree of boxes that the rays descend: the box of its faces and either the
 * faces themselves (a leaf) or two nodes that share them. */
typedef struct tree_node {
  double lo[3];
  double hi[3];
  size_t first; /* in the order of the tree */
  size_t count; /* the faces of a leaf; 0: an inner node */
  size_t left;
  size_t right;
} tree_node;

enum { LEAF_FACES = 4, TREE_STACK = 128 };

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
  tree_node* nodes; /* at most face_count */
  size_t node_count;
  size_t* order; /* the faces, leaf by leaf */
  /* what tells the inside from the outside at the nearest point (myb_polyhedron_nearest) */
  size_t* across;         /* 3 per face: the face across its edge from corner j to corner j + 1 */
  double* corner_normals; /* 3 per vertex: the normals of its faces weighted by their angles */
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
  myb_internal_release(&a, polyhedron->nodes);
  myb_internal_release(&a, polyhedron->order);
  myb_internal_release(&a, polyhedron->across);
  myb_internal_release(&a, polyhedron->corner_normals);
  myb_internal_release(&a, polyhedron);
}

/* ---- The tree of boxes ------------------------------------------------------------------------- */

static double face_centre(const myb_polyhedron* p, size_t f, int axis) {
  const size_t* v = p->faces + 3 * f;
  return (p->vertices[3 * v[0] + (size_t)axis] + p->vertices[3 * v[1] + (size_t)axis] +
          p->vertices[3 * v[2] + (size_t)axis]) / 3.0;
}

static void sift_down(const myb_polyhedron* p, size_t* order, size_t start, size_t end,
                      int axis) {
  size_t parent = start;
  /* bound: the depth of the heap */
  while (2 * parent + 1 < end) {
    size_t child = 2 * parent + 1;
    if (child + 1 < end &&
        face_centre(p, order[child], axis) < face_centre(p, order[child + 1], axis)) {
      ++child;
    }
    if (!(face_centre(p, order[parent], axis) < face_centre(p, order[child], axis))) {
      return;
    }
    const size_t swap = order[parent];
    order[parent] = order[child];
    order[child] = swap;
    parent = child;
  }
}

/* The faces order[0 .. count) by their centres along an axis: a heap sort. */
static void sort_faces(const myb_polyhedron* p, size_t* order, size_t count, int axis) {
  for (size_t start = count / 2; start-- > 0;) {
    sift_down(p, order, start, count, axis);
  }
  for (size_t end = count; end-- > 1;) {
    const size_t swap = order[0];
    order[0] = order[end];
    order[end] = swap;
    sift_down(p, order, 0, end, axis);
  }
}

/* Halves the faces along the longest side of their box until four or fewer are left: at most
 * face_count nodes, a depth of log2(face_count). */
static myb_status build_tree(myb_polyhedron* p) {
  const myb_allocator a = p->allocator;
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, p->face_count, sizeof(size_t), &block);
  if (status != MYB_OK) {
    return status;
  }
  p->order = block;
  status = myb_internal_alloc_array(&a, p->face_count, sizeof(tree_node), &block);
  if (status != MYB_OK) {
    return status;
  }
  p->nodes = block;
  for (size_t f = 0; f < p->face_count; ++f) {
    p->order[f] = f;
  }
  /* faces count a part in 1e12 larger than they are: so do their boxes, and more */
  double reach = 0.0;
  for (size_t i = 0; i < 3 * p->vertex_count; ++i) {
    reach = fmax(reach, fabs(p->vertices[i]));
  }
  const double pad = 1.0e-9 * reach;
  size_t stack[TREE_STACK];
  size_t top = 0;
  p->nodes[0] = (tree_node){.first = 0, .count = p->face_count};
  p->node_count = 1;
  stack[top++] = 0;
  /* bound: the nodes of the tree */
  while (top > 0) {
    tree_node* node = &p->nodes[stack[--top]];
    for (int i = 0; i < 3; ++i) {
      node->lo[i] = (double)INFINITY;
      node->hi[i] = -(double)INFINITY;
    }
    for (size_t k = 0; k < node->count; ++k) {
      const size_t* v = p->faces + 3 * p->order[node->first + k];
      for (int corner = 0; corner < 3; ++corner) {
        for (int i = 0; i < 3; ++i) {
          const double x = p->vertices[3 * v[corner] + (size_t)i];
          node->lo[i] = fmin(node->lo[i], x);
          node->hi[i] = fmax(node->hi[i], x);
        }
      }
    }
    int axis = 0;
    for (int i = 0; i < 3; ++i) {
      if (node->hi[i] - node->lo[i] > node->hi[axis] - node->lo[axis]) {
        axis = i;
      }
      node->lo[i] -= pad;
      node->hi[i] += pad;
    }
    if (node->count <= LEAF_FACES) {
      continue;
    }
    sort_faces(p, p->order + node->first, node->count, axis);
    const size_t half = node->count / 2;
    MYB_CHECK(p->node_count + 2 <= p->face_count && top + 2 <= TREE_STACK);
    node->left = p->node_count;
    node->right = p->node_count + 1;
    p->nodes[node->left] = (tree_node){.first = node->first, .count = half};
    p->nodes[node->right] = (tree_node){.first = node->first + half, .count = node->count - half};
    p->node_count += 2;
    node->count = 0;
    stack[top++] = node->left;
    stack[top++] = node->right;
  }
  return MYB_OK;
}

/* The corner of face f at the vertex `vertex` (one of its three). */
static size_t corner_of(const myb_polyhedron* p, size_t f, size_t vertex) {
  const size_t* v = p->faces + 3 * f;
  return v[0] == vertex ? 0u : (v[1] == vertex ? 1u : 2u);
}

/* The normal of every vertex as the sum of the normals of its faces weighted by the angle of
 * each face at the vertex (G. Thurmer, C. A. Wuthrich, J. Graphics Tools 3, 43, 1998). */
static void corner_normals(myb_polyhedron* p) {
  memset(p->corner_normals, 0, 3 * p->vertex_count * sizeof(double));
  for (size_t f = 0; f < p->face_count; ++f) {
    const double* n = p->normals + 3 * f;
    for (size_t j = 0; j < 3; ++j) {
      const size_t at = p->faces[3 * f + j];
      const double* v0 = p->vertices + 3 * at;
      double e1[3];
      double e2[3];
      double c[3];
      sub3(p->vertices + 3 * p->faces[3 * f + (j + 1) % 3], v0, e1);
      sub3(p->vertices + 3 * p->faces[3 * f + (j + 2) % 3], v0, e2);
      myb_vec3_cross(e1, e2, c);
      const double angle = myb_atan2(myb_vec3_norm(c), myb_vec3_dot(e1, e2));
      for (int i = 0; i < 3; ++i) {
        p->corner_normals[3 * at + (size_t)i] += angle * n[i];
      }
    }
  }
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
      status = myb_internal_alloc_array(&a, 3 * face_count, sizeof(size_t), &block);
    }
    if (status == MYB_OK) {
      p->across = block;
      status = myb_internal_alloc_array(&a, 3 * vertex_count, sizeof(double), &block);
    }
    if (status == MYB_OK) {
      p->corner_normals = block;
    }
  }
  if (status == MYB_OK) {
    for (size_t k = 0; k < p->edge_count; ++k) {
      const directed* x = &d[2 * k];
      const directed* y = &d[2 * k + 1];
      p->across[3 * x->face + corner_of(p, x->face, x->from)] = y->face;
      p->across[3 * y->face + corner_of(p, y->face, y->from)] = x->face;
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
  if (status == MYB_OK) {
    corner_normals(p);
    status = build_tree(p);
  }
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

/* ---- The mesh as a surface (myb/geometry.h) --------------------------------------------------- */

size_t myb_polyhedron_face_count(const myb_polyhedron* polyhedron) {
  MYB_CHECK(polyhedron != NULL);
  return polyhedron->face_count;
}

void myb_polyhedron_face_normal(const myb_polyhedron* polyhedron, size_t face, double normal[3]) {
  MYB_CHECK(polyhedron != NULL && face < polyhedron->face_count && normal != NULL);
  memcpy(normal, polyhedron->normals + 3 * face, 3 * sizeof(double));
}

double myb_polyhedron_solid_angle(const myb_polyhedron* polyhedron, const double r[3]) {
  MYB_CHECK(polyhedron != NULL && r != NULL);
  const myb_polyhedron* p = polyhedron;
  double omega = 0.0;
  /* A. van Oosterom, J. Strackee, IEEE Trans. Biomed. Eng. 30, 125 (1983), as in the field */
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
    omega += 2.0 * myb_atan2(num, den);
  }
  return omega;
}

/* Where the ray origin + t d meets the face f (T. Moller, B. Trumbore, J. Graphics Tools 2,
 * 21, 1997); 0 if it does not. */
static int face_hit(const myb_polyhedron* p, size_t f, const double origin[3], const double d[3],
                    double* distance) {
  static const double margin = 1.0e-12;   /* of a face, in barycentric coordinates */
  static const double parallel = 1.0e-14; /* of the product of two edges */
  const double* v0 = p->vertices + 3 * p->faces[3 * f];
  double e1[3];
  double e2[3];
  double s[3];
  sub3(p->vertices + 3 * p->faces[3 * f + 1], v0, e1);
  sub3(p->vertices + 3 * p->faces[3 * f + 2], v0, e2);
  double h[3];
  myb_vec3_cross(d, e2, h);
  const double det = myb_vec3_dot(e1, h);
  if (!(fabs(det) > parallel * myb_vec3_norm(e1) * myb_vec3_norm(e2))) {
    return 0;
  }
  sub3(origin, v0, s);
  const double u = myb_vec3_dot(s, h) / det;
  if (u < -margin || u > 1.0 + margin) {
    return 0;
  }
  double q[3];
  myb_vec3_cross(s, e1, q);
  const double v = myb_vec3_dot(d, q) / det;
  if (v < -margin || u + v > 1.0 + margin) {
    return 0;
  }
  const double t = myb_vec3_dot(e2, q) / det;
  if (!(t >= 0.0)) {
    return 0;
  }
  *distance = t;
  return 1;
}

/* Does the ray meet the box of a node before the distance `limit`? The slabs of the three
 * axes; a direction along an axis gives infinities, and 0 times infinity is not a number:
 * the comparisons then leave the interval as it was. */
static int box_hit(const tree_node* node, const double origin[3], const double inverse[3],
                   double limit) {
  double t0 = 0.0;
  double t1 = limit;
  for (int i = 0; i < 3; ++i) {
    double near = (node->lo[i] - origin[i]) * inverse[i];
    double far = (node->hi[i] - origin[i]) * inverse[i];
    if (near > far) {
      const double swap = near;
      near = far;
      far = swap;
    }
    if (near > t0) {
      t0 = near;
    }
    if (far < t1) {
      t1 = far;
    }
    if (t0 > t1) {
      return 0;
    }
  }
  return 1;
}

int myb_polyhedron_ray(const myb_polyhedron* polyhedron, const double origin[3],
                       const double direction[3], double* distance, double point[3],
                       size_t* face) {
  MYB_CHECK(polyhedron != NULL && origin != NULL && direction != NULL);
  const myb_polyhedron* p = polyhedron;
  if (!myb_vec3_finite(origin) || !myb_vec3_finite(direction)) {
    return 0;
  }
  const double length = myb_vec3_norm(direction);
  if (!(length > 0.0)) {
    return 0;
  }
  const double d[3] = {direction[0] / length, direction[1] / length, direction[2] / length};
  const double inverse[3] = {1.0 / d[0], 1.0 / d[1], 1.0 / d[2]};
  int hit = 0;
  double best = (double)INFINITY;
  size_t best_face = 0;
  size_t stack[TREE_STACK];
  size_t top = 0;
  stack[top++] = 0;
  /* bound: the nodes of the tree */
  while (top > 0) {
    const tree_node* node = &p->nodes[stack[--top]];
    if (!box_hit(node, origin, inverse, best)) {
      continue;
    }
    if (node->count == 0) {
      MYB_CHECK(top + 2 <= TREE_STACK);
      stack[top++] = node->left;
      stack[top++] = node->right;
      continue;
    }
    for (size_t k = 0; k < node->count; ++k) {
      const size_t f = p->order[node->first + k];
      double t = 0.0;
      /* of faces met at the same distance — an edge they share — the first in the mesh */
      if (face_hit(p, f, origin, d, &t) && (!hit || t < best || (t == best && f < best_face))) {
        hit = 1;
        best = t;
        best_face = f;
      }
    }
  }
  if (!hit) {
    return 0;
  }
  if (distance != NULL) {
    *distance = best;
  }
  if (point != NULL) {
    for (int i = 0; i < 3; ++i) {
      point[i] = origin[i] + best * d[i];
    }
  }
  if (face != NULL) {
    *face = best_face;
  }
  return 1;
}

/* What part of a face holds the point nearest to a position: the inside of the face, an edge
 * (from corner j to corner j + 1) or a corner. */
enum { NEAR_FACE = 0, NEAR_EDGE = 1, NEAR_CORNER = 4 };

/* The point q of face f nearest to r, the part of the face it lies in (NEAR_EDGE + j,
 * NEAR_CORNER + j) and the squared distance: the regions of the corners, the edges and the
 * inside told apart by the projections on the edges (C. Ericson, Real-Time Collision
 * Detection, 2005, §5.1.5). */
static double face_nearest(const myb_polyhedron* p, size_t f, const double r[3], double q[3],
                           int* part) {
  const double* a = p->vertices + 3 * p->faces[3 * f];
  const double* b = p->vertices + 3 * p->faces[3 * f + 1];
  const double* c = p->vertices + 3 * p->faces[3 * f + 2];
  double ab[3];
  double ac[3];
  double ar[3];
  double br[3];
  double cr[3];
  sub3(b, a, ab);
  sub3(c, a, ac);
  sub3(r, a, ar);
  sub3(r, b, br);
  sub3(r, c, cr);
  const double d1 = myb_vec3_dot(ab, ar);
  const double d2 = myb_vec3_dot(ac, ar);
  const double d3 = myb_vec3_dot(ab, br);
  const double d4 = myb_vec3_dot(ac, br);
  const double d5 = myb_vec3_dot(ab, cr);
  const double d6 = myb_vec3_dot(ac, cr);
  const double vc = d1 * d4 - d3 * d2;
  const double vb = d5 * d2 - d1 * d6;
  const double va = d3 * d6 - d5 * d4;
  if (d1 <= 0.0 && d2 <= 0.0) {
    memcpy(q, a, 3 * sizeof *q);
    *part = NEAR_CORNER;
  } else if (d3 >= 0.0 && d4 <= d3) {
    memcpy(q, b, 3 * sizeof *q);
    *part = NEAR_CORNER + 1;
  } else if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double s = d1 / (d1 - d3);
    for (int i = 0; i < 3; ++i) {
      q[i] = a[i] + s * ab[i];
    }
    *part = NEAR_EDGE;
  } else if (d6 >= 0.0 && d5 <= d6) {
    memcpy(q, c, 3 * sizeof *q);
    *part = NEAR_CORNER + 2;
  } else if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double s = d2 / (d2 - d6);
    for (int i = 0; i < 3; ++i) {
      q[i] = a[i] + s * ac[i];
    }
    *part = NEAR_EDGE + 2;
  } else if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0) {
    const double s = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    for (int i = 0; i < 3; ++i) {
      q[i] = b[i] + s * (c[i] - b[i]);
    }
    *part = NEAR_EDGE + 1;
  } else {
    const double sum = va + vb + vc;
    const double v = vb / sum;
    const double w = vc / sum;
    for (int i = 0; i < 3; ++i) {
      q[i] = a[i] + v * ab[i] + w * ac[i];
    }
    *part = NEAR_FACE;
  }
  const double d[3] = {r[0] - q[0], r[1] - q[1], r[2] - q[2]};
  return myb_vec3_dot(d, d);
}

/* The squared distance from r to a box of the tree (0 inside it). */
static double box_distance2(const tree_node* node, const double r[3]) {
  double sum = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double below = node->lo[i] - r[i];
    const double above = r[i] - node->hi[i];
    const double d = below > 0.0 ? below : (above > 0.0 ? above : 0.0);
    sum += d * d;
  }
  return sum;
}

myb_status myb_polyhedron_nearest(const myb_polyhedron* polyhedron, const double position[3],
                                  double point[3], double normal[3], double* distance,
                                  size_t* face) {
  MYB_CHECK(polyhedron != NULL && position != NULL);
  const myb_polyhedron* p = polyhedron;
  if (!myb_vec3_finite(position)) {
    return MYB_E_DOMAIN;
  }
  double best = (double)INFINITY;
  size_t best_face = 0;
  int best_part = NEAR_FACE;
  double best_point[3] = {0.0, 0.0, 0.0};
  size_t stack[TREE_STACK];
  size_t top = 0;
  stack[top++] = 0;
  /* bound: the nodes of the tree; a box farther than the best face so far is not opened, and
   * of two children the nearer is opened first */
  while (top > 0) {
    const tree_node* node = &p->nodes[stack[--top]];
    if (!(box_distance2(node, position) <= best)) {
      continue;
    }
    if (node->count == 0) {
      MYB_CHECK(top + 2 <= TREE_STACK);
      const int left_first = box_distance2(&p->nodes[node->left], position) <=
                             box_distance2(&p->nodes[node->right], position);
      stack[top++] = left_first ? node->right : node->left;
      stack[top++] = left_first ? node->left : node->right;
      continue;
    }
    for (size_t k = 0; k < node->count; ++k) {
      const size_t f = p->order[node->first + k];
      double q[3];
      int part = NEAR_FACE;
      const double d2 = face_nearest(p, f, position, q, &part);
      /* of faces equally near — an edge or a corner they share — the first in the mesh */
      if (d2 < best || (d2 == best && f < best_face)) {
        best = d2;
        best_face = f;
        best_part = part;
        memcpy(best_point, q, sizeof q);
      }
    }
  }
  if (!(best < (double)INFINITY)) {
    return MYB_E_DOMAIN;
  }
  /* inside or outside: the side of the angle-weighted normal of the part nearest
   * (J. A. Baerentzen, H. Aanaes, IEEE Trans. Vis. Comput. Graph. 11, 243, 2005) */
  const double* nf = p->normals + 3 * best_face;
  double pseudo[3];
  if (best_part == NEAR_FACE) {
    memcpy(pseudo, nf, sizeof pseudo);
  } else if (best_part < NEAR_CORNER) {
    const double* ng = p->normals + 3 * p->across[3 * best_face + (size_t)(best_part - NEAR_EDGE)];
    for (int i = 0; i < 3; ++i) {
      pseudo[i] = nf[i] + ng[i];
    }
  } else {
    const size_t vertex = p->faces[3 * best_face + (size_t)(best_part - NEAR_CORNER)];
    memcpy(pseudo, p->corner_normals + 3 * vertex, sizeof pseudo);
  }
  const double d[3] = {position[0] - best_point[0], position[1] - best_point[1],
                       position[2] - best_point[2]};
  const double length = myb_vec3_norm(d);
  const int inside = myb_vec3_dot(d, pseudo) < 0.0;
  if (point != NULL) {
    memcpy(point, best_point, sizeof best_point);
  }
  if (distance != NULL) {
    *distance = inside ? -length : length;
  }
  if (face != NULL) {
    *face = best_face;
  }
  if (normal != NULL) {
    /* the face itself: its normal; an edge or a corner: the line to the position, turned
     * outward — or the angle-weighted normal for a position on it */
    const double pn = myb_vec3_norm(pseudo);
    if (best_part == NEAR_FACE || (!(length > 0.0) && !(pn > 0.0))) {
      memcpy(normal, nf, 3 * sizeof *normal);
    } else if (length > 0.0) {
      const double k = (inside ? -1.0 : 1.0) / length;
      for (int i = 0; i < 3; ++i) {
        normal[i] = k * d[i];
      }
    } else {
      for (int i = 0; i < 3; ++i) {
        normal[i] = pseudo[i] / pn;
      }
    }
  }
  return MYB_OK;
}

/* ---- Light on the mesh (myb/surface.h) --------------------------------------------------------- */

/* The force of the light on a part of a face of area `area` at the cosine `cosine` of the
 * angle of incidence, as a plate (myb_force_radiation_plates). */
static void plate_force(const myb_plate* surface, double pressure, double area, double cosine,
                        const double u[3], const double n[3], double f[3]) {
  const double absorbed = 1.0 - surface->specular - surface->diffuse;
  const double along = (1.0 - surface->specular) * pressure * area * cosine;
  const double normal = pressure * area * cosine *
                        (2.0 * (surface->specular * cosine + surface->diffuse / 3.0) +
                         (surface->thermal ? 2.0 / 3.0 * absorbed : 0.0));
  for (int i = 0; i < 3; ++i) {
    f[i] = along * u[i] - normal * n[i];
  }
}

/* d f / d u of plate_force, added to g (rows): with c = -n.u,
 *   P A [(1 - rho_s) (c 1 - u n^T) + (N(c) + 2 rho_s c) n n^T],
 * N the factor of the normal part. */
static void plate_gradient(const myb_plate* surface, double pressure, double area, double cosine,
                           const double u[3], const double n[3], double g[9]) {
  const double absorbed = 1.0 - surface->specular - surface->diffuse;
  const double along = (1.0 - surface->specular) * pressure * area;
  const double normal = pressure * area *
                        (2.0 * (surface->specular * cosine + surface->diffuse / 3.0) +
                         (surface->thermal ? 2.0 / 3.0 * absorbed : 0.0) +
                         2.0 * surface->specular * cosine);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      g[3 * i + j] += along * ((i == j ? cosine : 0.0) - u[i] * n[j]) + normal * n[i] * n[j];
    }
  }
}

static myb_status radiation_mesh(const myb_polyhedron* shape, const myb_plate* surface,
                                 double pressure, const double u[3], int shadows,
                                 double force[3], double torque[3], double* gradient);

enum { MAX_LAG_PHASES = 256 }; /* steps of a rotation of the late re-emission */

myb_status myb_force_radiation_mesh(const myb_polyhedron* shape, const myb_plate* surface,
                                    double pressure, const double u[3], int shadows,
                                    double force[3], double torque[3]) {
  return radiation_mesh(shape, surface, pressure, u, shadows, force, torque, NULL);
}

/* The light of the moment without the re-emission, plus each element's re-emission answered
 * late (myb_force_radiation_mesh_lagged). */
myb_status myb_force_radiation_mesh_lagged(const myb_polyhedron* shape, const myb_plate* surface,
                                           double pressure, const double u[3],
                                           const double axis[3], double theta, int phases,
                                           int shadows, double force[3], double torque[3]) {
  MYB_CHECK(shape != NULL && surface != NULL && u != NULL && axis != NULL && force != NULL &&
            torque != NULL && shadows >= 0);
  const myb_polyhedron* p = shape;
  if (!(fabs(myb_vec3_dot(axis, axis) - 1.0) <= 1e-12) || !(theta >= 0.0) || !isfinite(theta) ||
      phases < 8 || phases > MAX_LAG_PHASES) {
    return MYB_E_DOMAIN;
  }
  myb_plate cold = *surface;
  cold.thermal = 0;
  myb_status status = radiation_mesh(p, &cold, pressure, u, shadows, force, torque, NULL);
  if (status != MYB_OK || !surface->thermal) {
    return status;
  }
  const double absorbed = 1.0 - surface->specular - surface->diffuse;
  /* the kernel of the late answer: emitted(0) = sum_j kernel_j absorbed(phi_j), with
   * kernel_j = (1/N) [1 + 2 sum_k Re(H_k e^(-i k phi_j)) (+ the Nyquist term)],
   * H_k = 1 / (1 + theta sqrt(k) (1 + i) / 2); for theta = 0 it picks the moment alone */
  double kernel[MAX_LAG_PHASES];
  const int top = (phases - 1) / 2; /* the harmonics below the Nyquist one */
  for (int step = 0; step < phases; ++step) {
    double sum = 1.0;
    for (int h = 1; h <= top; ++h) {
      const double d_re = 1.0 + 0.5 * theta * sqrt((double)h);
      const double d_im = 0.5 * theta * sqrt((double)h);
      const double dd = d_re * d_re + d_im * d_im;
      const double angle = 6.283185307179586 * (double)((h * step) % phases) / (double)phases;
      /* Re((d_re - i d_im) / dd * (cos - i sin)): the harmonic k of the past answered now */
      sum += 2.0 * (d_re * myb_cos(angle) - d_im * myb_sin(angle)) / dd;
    }
    if (phases % 2 == 0) {
      const int h = phases / 2;
      const double d_re = 1.0 + 0.5 * theta * sqrt((double)h);
      const double d_im = 0.5 * theta * sqrt((double)h);
      const double dd = d_re * d_re + d_im * d_im;
      sum += (step % 2 == 0 ? 1.0 : -1.0) * d_re / dd;
    }
    kernel[step] = sum / (double)phases;
  }
  double emitted_total[3] = {0.0, 0.0, 0.0};
  double emitted_moment[3] = {0.0, 0.0, 0.0};
  const int parts = shadows > 0 ? shadows : 1;
  const double k = (double)parts;
  /* the direction of the source at each phase: the body has turned by phi, so the source is
   * seen turned by -phi about the axis (Rodrigues) */
  const double s0[3] = {-u[0], -u[1], -u[2]};
  double along = myb_vec3_dot(axis, s0);
  double across[3];
  double third_dir[3];
  for (int q = 0; q < 3; ++q) {
    across[q] = s0[q] - along * axis[q];
  }
  myb_vec3_cross(axis, across, third_dir);
  for (size_t f = 0; f < p->face_count; ++f) {
    const double* n = p->normals + 3 * f;
    const double* v0 = p->vertices + 3 * p->faces[3 * f];
    double e1[3];
    double e2[3];
    double normal[3];
    sub3(p->vertices + 3 * p->faces[3 * f + 1], v0, e1);
    sub3(p->vertices + 3 * p->faces[3 * f + 2], v0, e2);
    myb_vec3_cross(e1, e2, normal);
    const double area = 0.5 * myb_vec3_norm(normal);
    const double part_area = area / (k * k);
    const double lift = 1.0e-9 * fmax(myb_vec3_norm(e1), myb_vec3_norm(e2));
    for (int i = 0; i < parts; ++i) {
      for (int j = 0; j < parts - i; ++j) {
        for (int turned = 0; turned < 2; ++turned) {
          if (turned && i + j > parts - 2) {
            continue;
          }
          const double third = turned ? 2.0 / 3.0 : 1.0 / 3.0;
          const double a = shadows > 0 ? ((double)i + third) / k : 1.0 / 3.0;
          const double b = shadows > 0 ? ((double)j + third) / k : 1.0 / 3.0;
          double centre[3];
          double start[3];
          for (int q = 0; q < 3; ++q) {
            centre[q] = v0[q] + a * e1[q] + b * e2[q];
            start[q] = centre[q] + lift * n[q];
          }
          /* what the element emits now: its absorbed power over the rotation weighted by the
           * kernel of the late answer */
          double now = 0.0;
          for (int step = 0; step < phases; ++step) {
            const double phi = 6.283185307179586 * (double)step / (double)phases;
            const double cp = myb_cos(phi);
            const double sp = myb_sin(phi);
            double to_source[3];
            for (int q = 0; q < 3; ++q) {
              to_source[q] = along * axis[q] + cp * across[q] - sp * third_dir[q];
            }
            const double cosine = myb_vec3_dot(n, to_source);
            if (!(cosine > 0.0)) {
              continue;
            }
            if (shadows > 0 && myb_polyhedron_ray(p, start, to_source, NULL, NULL, NULL)) {
              continue;
            }
            now += kernel[step] * absorbed * pressure * part_area * cosine;
          }
          /* the recoil of a Lambertian radiator: (2/3) of the power along -n */
          double push[3];
          double lever[3];
          for (int q = 0; q < 3; ++q) {
            push[q] = -2.0 / 3.0 * now * n[q];
          }
          myb_vec3_cross(centre, push, lever);
          for (int q = 0; q < 3; ++q) {
            emitted_total[q] += push[q];
            emitted_moment[q] += lever[q];
          }
        }
      }
    }
  }
  for (int q = 0; q < 3; ++q) {
    force[q] += emitted_total[q];
    torque[q] += emitted_moment[q];
  }
  return MYB_OK;
}

myb_status myb_internal_radiation_mesh_gradient(const myb_polyhedron* shape,
                                                const myb_plate* surface, double pressure,
                                                const double u[3], int shadows, double force[3],
                                                double gradient[9]) {
  MYB_CHECK(gradient != NULL);
  memset(gradient, 0, 9 * sizeof(double));
  return radiation_mesh(shape, surface, pressure, u, shadows, force, NULL, gradient);
}

/* The force and the torque of the light on the mesh; with `gradient`, also d force / d u for
 * the parts that are lit now (which parts are lit changes in steps: the derivative is that of
 * the force between two of them). */
static myb_status radiation_mesh(const myb_polyhedron* shape, const myb_plate* surface,
                                 double pressure, const double u[3], int shadows,
                                 double force[3], double torque[3], double* gradient) {
  MYB_CHECK(shape != NULL && surface != NULL && u != NULL && shadows >= 0);
  const myb_polyhedron* p = shape;
  const double u2 = myb_vec3_dot(u, u);
  if (!isfinite(pressure) || !(fabs(u2 - 1.0) <= 1e-12) || !(surface->specular >= 0.0) ||
      !(surface->diffuse >= 0.0) || !(surface->specular + surface->diffuse <= 1.0)) {
    return MYB_E_DOMAIN;
  }
  const double to_source[3] = {-u[0], -u[1], -u[2]};
  double total[3] = {0.0, 0.0, 0.0};
  double moment[3] = {0.0, 0.0, 0.0};
  for (size_t f = 0; f < p->face_count; ++f) {
    const double* n = p->normals + 3 * f;
    const double cosine = -myb_vec3_dot(n, u);
    if (!(cosine > 0.0)) {
      continue;
    }
    const double* v0 = p->vertices + 3 * p->faces[3 * f];
    double e1[3];
    double e2[3];
    double normal[3];
    sub3(p->vertices + 3 * p->faces[3 * f + 1], v0, e1);
    sub3(p->vertices + 3 * p->faces[3 * f + 2], v0, e2);
    myb_vec3_cross(e1, e2, normal);
    const double area = 0.5 * myb_vec3_norm(normal);
    if (shadows == 0) {
      double part[3];
      double centre[3];
      double lever[3];
      plate_force(surface, pressure, area, cosine, u, n, part);
      if (gradient != NULL) {
        plate_gradient(surface, pressure, area, cosine, u, n, gradient);
      }
      for (int i = 0; i < 3; ++i) {
        centre[i] = v0[i] + (e1[i] + e2[i]) / 3.0;
        total[i] += part[i];
      }
      myb_vec3_cross(centre, part, lever);
      for (int i = 0; i < 3; ++i) {
        moment[i] += lever[i];
      }
      continue;
    }
    /* k^2 equal triangles: those standing as the face and, between them, those turned over */
    const double k = (double)shadows;
    const double lift = 1.0e-9 * fmax(myb_vec3_norm(e1), myb_vec3_norm(e2));
    const double part_area = area / (k * k);
    for (int i = 0; i < shadows; ++i) {
      for (int j = 0; j < shadows - i; ++j) {
        for (int turned = 0; turned < 2; ++turned) {
          if (turned && i + j > shadows - 2) {
            continue;
          }
          const double third = turned ? 2.0 / 3.0 : 1.0 / 3.0;
          const double a = ((double)i + third) / k;
          const double b = ((double)j + third) / k;
          double centre[3];
          double start[3];
          for (int q = 0; q < 3; ++q) {
            centre[q] = v0[q] + a * e1[q] + b * e2[q];
            start[q] = centre[q] + lift * n[q];
          }
          if (myb_polyhedron_ray(p, start, to_source, NULL, NULL, NULL)) {
            continue; /* in the shadow of another part */
          }
          double part[3];
          double lever[3];
          plate_force(surface, pressure, part_area, cosine, u, n, part);
          if (gradient != NULL) {
            plate_gradient(surface, pressure, part_area, cosine, u, n, gradient);
          }
          myb_vec3_cross(centre, part, lever);
          for (int q = 0; q < 3; ++q) {
            total[q] += part[q];
            moment[q] += lever[q];
          }
        }
      }
    }
  }
  if (force != NULL) {
    memcpy(force, total, sizeof total);
  }
  if (torque != NULL) {
    memcpy(torque, moment, sizeof moment);
  }
  return MYB_OK;
}

myb_status myb_force_yorp(const myb_polyhedron* shape, double pressure, double sun_colatitude,
                          double theta, int shadows, int phases, double torque[3]) {
  MYB_CHECK(shape != NULL && torque != NULL);
  const myb_polyhedron* p = shape;
  if (!isfinite(pressure) || !isfinite(sun_colatitude) || !(theta >= 0.0) || !isfinite(theta) ||
      shadows < 0 || phases < 8 || phases > 100000) {
    return MYB_E_DOMAIN;
  }
  const double sine = myb_sin(sun_colatitude);
  const double cosine_axis = myb_cos(sun_colatitude);
  double along = 0.0;    /* the mean torque along the axis */
  double across_re = 0.0; /* the first harmonic of the torque across it, in the inertial axes */
  double across_im = 0.0;
  for (int step = 0; step < phases; ++step) {
    const double phi = 6.283185307179586 * ((double)step + 0.5) / (double)phases;
    const double cp = myb_cos(phi);
    const double sp = myb_sin(phi);
    /* the body has turned by phi: the Sun is seen at -phi */
    const double to_source[3] = {sine * cp, -sine * sp, cosine_axis};
    double moment[3] = {0.0, 0.0, 0.0};
    for (size_t f = 0; f < p->face_count; ++f) {
      const double* n = p->normals + 3 * f;
      const double cosine = myb_vec3_dot(n, to_source);
      if (!(cosine > 0.0)) {
        continue;
      }
      const double* v0 = p->vertices + 3 * p->faces[3 * f];
      double e1[3];
      double e2[3];
      double normal[3];
      sub3(p->vertices + 3 * p->faces[3 * f + 1], v0, e1);
      sub3(p->vertices + 3 * p->faces[3 * f + 2], v0, e2);
      myb_vec3_cross(e1, e2, normal);
      const double area = 0.5 * myb_vec3_norm(normal);
      const int parts = shadows > 0 ? shadows : 1;
      const double k = (double)parts;
      const double lift = 1.0e-9 * fmax(myb_vec3_norm(e1), myb_vec3_norm(e2));
      /* the recoil of what a part re-emits: (2/3) p A cos along -n */
      const double recoil = 2.0 / 3.0 * pressure * (area / (k * k)) * cosine;
      for (int i = 0; i < parts; ++i) {
        for (int j = 0; j < parts - i; ++j) {
          for (int turned = 0; turned < 2; ++turned) {
            if (turned && i + j > parts - 2) {
              continue;
            }
            const double third = turned ? 2.0 / 3.0 : 1.0 / 3.0;
            const double a = ((double)i + third) / k;
            const double b = ((double)j + third) / k;
            double centre[3];
            double start[3];
            for (int q = 0; q < 3; ++q) {
              centre[q] = v0[q] + a * e1[q] + b * e2[q];
              start[q] = centre[q] + lift * n[q];
            }
            if (shadows > 0 && myb_polyhedron_ray(p, start, to_source, NULL, NULL, NULL)) {
              continue;
            }
            double lever[3];
            myb_vec3_cross(centre, n, lever);
            for (int q = 0; q < 3; ++q) {
              moment[q] -= recoil * lever[q];
            }
          }
        }
      }
    }
    along += moment[2];
    /* body components to inertial ones: times e^(i phi) */
    across_re += moment[0] * cp - moment[1] * sp;
    across_im += moment[0] * sp + moment[1] * cp;
  }
  along /= (double)phases;
  across_re /= (double)phases;
  across_im /= (double)phases;
  /* the first harmonic answers late: times the conjugate of 1 / (1 + theta (1 + i) / 2) */
  const double re = 1.0 + 0.5 * theta;
  const double im = 0.5 * theta;
  const double norm = re * re + im * im;
  const double answer_re = re / norm;
  const double answer_im = im / norm;
  torque[0] = across_re * answer_re - across_im * answer_im;
  torque[1] = across_re * answer_im + across_im * answer_re;
  torque[2] = along;
  return MYB_OK;
}
