/* SPDX-License-Identifier: BSD-3-Clause */
/* Magnetic fields from Gauss coefficients, the Lorentz force, co-rotation, plasma streams. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/force.h"
#include "myb/grav.h"
#include "myb/io.h"
#include "num/vec3.h"

struct myb_magnetic_field {
  myb_allocator allocator;
  myb_gravity_field* potential; /* V / a with GM = a^2, C = g / sqrt(2n + 1) */
  int degree;
};

static size_t coefficient_count(int degree) {
  return (size_t)(degree + 1) * (size_t)(degree + 2) / 2;
}

myb_status myb_magnetic_field_create(double radius, int degree, const double* g, const double* h,
                                     const myb_allocator* allocator,
                                     myb_magnetic_field** field) {
  MYB_CHECK(g != NULL && h != NULL && field != NULL);
  if (degree < 1 || degree > 2700 || !(radius > 0.0) || !isfinite(radius)) {
    return MYB_E_DOMAIN;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  const size_t count = coefficient_count(degree);
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, 2 * count, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  double* c = block;
  double* s = c + count;
  c[0] = 0.0;
  s[0] = 0.0;
  for (int n = 1; n <= degree && st == MYB_OK; ++n) {
    const double k = 1.0 / sqrt(2.0 * n + 1.0);
    for (int m = 0; m <= n; ++m) {
      const size_t i = (size_t)n * (size_t)(n + 1) / 2 + (size_t)m;
      if (!isfinite(g[i]) || !isfinite(h[i])) {
        st = MYB_E_DOMAIN;
        break;
      }
      c[i] = g[i] * k;
      s[i] = m == 0 ? 0.0 : h[i] * k;
    }
  }
  void* fblock = NULL;
  if (st == MYB_OK) {
    st = myb_internal_alloc_array(&a, 1, sizeof(myb_magnetic_field), &fblock);
  }
  myb_gravity_field* potential = NULL;
  if (st == MYB_OK) {
    st = myb_gravity_field_create(radius * radius, radius, degree, c, s, allocator, &potential);
  }
  myb_internal_release(&a, block);
  if (st != MYB_OK) {
    myb_internal_release(&a, fblock);
    return st;
  }
  myb_magnetic_field* f = fblock;
  f->allocator = a;
  f->potential = potential;
  f->degree = degree;
  *field = f;
  return MYB_OK;
}

void myb_magnetic_field_destroy(myb_magnetic_field* field) {
  if (field == NULL) {
    return;
  }
  const myb_allocator a = field->allocator;
  myb_gravity_field_destroy(field->potential);
  myb_internal_release(&a, field);
}

int myb_magnetic_field_degree(const myb_magnetic_field* field) {
  MYB_CHECK(field != NULL);
  return field->degree;
}

myb_status myb_magnetic_field_eval(const myb_magnetic_field* field, myb_gravity_work* work,
                                   const double r[3], double b[3]) {
  MYB_CHECK(field != NULL && work != NULL && r != NULL && b != NULL);
  double grad[3];
  const myb_status st = myb_gravity_field_eval(field->potential, work, r, 1, -1, NULL, grad);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    b[i] = -grad[i];
  }
  return MYB_OK;
}

/* ---- Coefficient tables in the IGRF format -------------------------------------------------- */

enum { MAX_EPOCHS = 64, MAX_TABLE_DEGREE = 30 };

/* Next white-space separated token of the line at *p (the end of the line excluded) if it is a
 * number as a whole; 0 at the end of the line or for a token that is not a number. */
static int next_number(const char** p, const char* end, double* value) {
  while (*p < end && (**p == ' ' || **p == '\t' || **p == '\r')) { /* bound: the line */
    ++*p;
  }
  const char* start = *p;
  const char* stop = start;
  while (stop < end && *stop != ' ' && *stop != '\t' && *stop != '\r') { /* bound: the line */
    ++stop;
  }
  if (stop == start) {
    return 0;
  }
  size_t used = 0;
  if (myb_parse_double(start, (size_t)(stop - start), value, &used) != MYB_OK ||
      used != (size_t)(stop - start)) {
    return 0;
  }
  *p = stop;
  return 1;
}

myb_status myb_magnetic_field_load_table(const char* text, size_t size, double radius,
                                         double year, const myb_allocator* allocator,
                                         myb_magnetic_field** field, myb_diag* diag) {
  MYB_CHECK(text != NULL && field != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  const size_t count = coefficient_count(MAX_TABLE_DEGREE);
  void* block = NULL;
  myb_status st = myb_internal_alloc_array(&a, 2 * count, sizeof(double), &block);
  if (st != MYB_OK) {
    return st;
  }
  double* g = block;
  double* h = g + count;
  memset(g, 0, 2 * count * sizeof(double));
  double epochs[MAX_EPOCHS];
  int nepochs = 0;
  int degree = 0;
  long line_number = 0;
  const char* p = text;
  const char* end_text = text + size;
  while (p < end_text) { /* bound: the lines of the text */
    const char* eol = memchr(p, '\n', (size_t)(end_text - p));
    const char* end = eol != NULL ? eol : end_text;
    ++line_number;
    const char* q = p;
    while (q < end && (*q == ' ' || *q == '\t')) { /* bound: the line */
      ++q;
    }
    if (end - q >= 3 && strncmp(q, "g/h", 3) == 0) {
      /* "g/h n m e1 e2 ... eN SV": the numbers are the epochs */
      const char* r = q + 3;
      while (r < end && (*r == ' ' || *r == '\t' || *r == 'n' || *r == 'm')) { /* bound */
        ++r;
      }
      double e;
      nepochs = 0;
      while (nepochs < MAX_EPOCHS && next_number(&r, end, &e)) { /* bound: MAX_EPOCHS */
        epochs[nepochs++] = e;
      }
    } else if (end - q >= 2 && (*q == 'g' || *q == 'h') && (q[1] == ' ' || q[1] == '\t')) {
      if (nepochs < 1) {
        myb_internal_diag(diag, line_number, "coefficients before the header line g/h n m ...");
        st = MYB_E_FORMAT;
        goto fail;
      }
      const char* r = q + 1;
      double nn = 0.0;
      double mm = 0.0;
      if (!next_number(&r, end, &nn) || !next_number(&r, end, &mm) || nn != floor(nn) ||
          mm != floor(mm) || nn < 1.0 || nn > (double)MAX_TABLE_DEGREE || mm < 0.0 || mm > nn) {
        myb_internal_diag(diag, line_number, "invalid degree or order");
        st = MYB_E_FORMAT;
        goto fail;
      }
      double values[MAX_EPOCHS + 1];
      int nv = 0;
      while (nv < MAX_EPOCHS + 1 && next_number(&r, end, &values[nv])) { /* bound */
        ++nv;
      }
      if (nv != nepochs + 1) {
        myb_internal_diag(diag, line_number,
                          "expected %d values (the epochs and the secular variation)",
                          nepochs + 1);
        st = MYB_E_FORMAT;
        goto fail;
      }
      if (!(year >= epochs[0])) {
        myb_internal_diag(diag, 0, "year %.3f before the first epoch %.1f", year, epochs[0]);
        st = MYB_E_RANGE;
        goto fail;
      }
      double v;
      if (year >= epochs[nepochs - 1]) {
        v = values[nepochs - 1] + values[nepochs] * (year - epochs[nepochs - 1]);
      } else {
        int k = 0;
        while (k + 2 < nepochs && year >= epochs[k + 1]) { /* bound: nepochs */
          ++k;
        }
        const double w = (year - epochs[k]) / (epochs[k + 1] - epochs[k]);
        v = values[k] + w * (values[k + 1] - values[k]);
      }
      const int n = (int)nn;
      const int m = (int)mm;
      const size_t idx = (size_t)n * (size_t)(n + 1) / 2 + (size_t)m;
      if (*q == 'g') {
        g[idx] = v * 1e-9;
      } else {
        h[idx] = v * 1e-9;
      }
      if (n > degree) {
        degree = n;
      }
    }
    p = end < end_text ? end + 1 : end_text;
  }
  if (degree < 1) {
    myb_internal_diag(diag, 0, "no coefficients");
    st = MYB_E_FORMAT;
    goto fail;
  }
  st = myb_magnetic_field_create(radius, degree, g, h, allocator, field);
fail:
  myb_internal_release(&a, block);
  return st;
}

/* ---- Forces ------------------------------------------------------------------------------- */

void myb_lorentz_acceleration(double charge_to_mass, const double e[3], const double v[3],
                              const double b[3], double a[3]) {
  MYB_CHECK(e != NULL && v != NULL && b != NULL && a != NULL);
  double vxb[3];
  myb_vec3_cross(v, b, vxb);
  for (int i = 0; i < 3; ++i) {
    a[i] = charge_to_mass * (e[i] + vxb[i]);
  }
}

void myb_corotation_field(const double omega[3], const double r[3], const double b[3],
                          double e[3]) {
  MYB_CHECK(omega != NULL && r != NULL && b != NULL && e != NULL);
  double u[3];
  myb_vec3_cross(omega, r, u);
  double uxb[3];
  myb_vec3_cross(u, b, uxb);
  for (int i = 0; i < 3; ++i) {
    e[i] = -uxb[i];
  }
}

void myb_force_plasma_stream(double rho0, double r0, double speed, double cd_area_over_mass,
                             const double r[3], const double v[3], double a[3]) {
  MYB_CHECK(r != NULL && v != NULL && a != NULL);
  const double rr = myb_vec3_norm(r);
  double w[3];
  for (int i = 0; i < 3; ++i) {
    w[i] = speed * r[i] / rr - v[i];
  }
  const double rho = rho0 * (r0 / rr) * (r0 / rr);
  const double f = 0.5 * cd_area_over_mass * rho * myb_vec3_norm(w);
  for (int i = 0; i < 3; ++i) {
    a[i] = f * w[i];
  }
}
