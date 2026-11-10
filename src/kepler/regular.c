/* SPDX-License-Identifier: BSD-3-Clause */
/* Kustaanheimo-Stiefel and Sundman regularisations. */
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "core/check.h"
#include "myb/regular.h"
#include "num/vec3.h"

/* L(u) w */
static void ks_l(const double u[4], const double w[4], double out[4]) {
  out[0] = u[0] * w[0] - u[1] * w[1] - u[2] * w[2] + u[3] * w[3];
  out[1] = u[1] * w[0] + u[0] * w[1] - u[3] * w[2] - u[2] * w[3];
  out[2] = u[2] * w[0] + u[3] * w[1] + u[0] * w[2] + u[1] * w[3];
  out[3] = u[3] * w[0] - u[2] * w[1] + u[1] * w[2] - u[0] * w[3];
}

/* L(u)^T w */
static void ks_lt(const double u[4], const double w[4], double out[4]) {
  out[0] = u[0] * w[0] + u[1] * w[1] + u[2] * w[2] + u[3] * w[3];
  out[1] = -u[1] * w[0] + u[0] * w[1] + u[3] * w[2] - u[2] * w[3];
  out[2] = -u[2] * w[0] - u[3] * w[1] + u[0] * w[2] + u[1] * w[3];
  out[3] = u[3] * w[0] - u[2] * w[1] + u[1] * w[2] - u[0] * w[3];
}

void myb_ks_from_cartesian(double mu, double t, const double x[3], const double v[3],
                           double y[10]) {
  MYB_CHECK(x != NULL && v != NULL && y != NULL);
  const double r = myb_vec3_norm(x);
  double u[4] = {0.0, 0.0, 0.0, 0.0};
  if (r > 0.0) {
    if (x[0] >= 0.0) {
      u[0] = sqrt(0.5 * (r + x[0]));
      u[1] = x[1] / (2.0 * u[0]);
      u[2] = x[2] / (2.0 * u[0]);
    } else {
      u[1] = sqrt(0.5 * (r - x[0]));
      u[0] = x[1] / (2.0 * u[1]);
      u[3] = x[2] / (2.0 * u[1]);
    }
  }
  const double w[4] = {v[0], v[1], v[2], 0.0};
  double up[4];
  ks_lt(u, w, up);
  for (int i = 0; i < 4; ++i) {
    y[i] = u[i];
    y[4 + i] = 0.5 * up[i];
  }
  y[8] = r > 0.0 ? mu / r - 0.5 * myb_vec3_dot(v, v) : 0.0;
  y[9] = t;
}

void myb_ks_to_cartesian(const double y[10], double* t, double x[3], double v[3]) {
  MYB_CHECK(y != NULL && x != NULL && v != NULL);
  const double* u = y;
  double xx[4];
  double vv[4];
  ks_l(u, u, xx);
  ks_l(u, y + 4, vv);
  const double r = u[0] * u[0] + u[1] * u[1] + u[2] * u[2] + u[3] * u[3];
  for (int i = 0; i < 3; ++i) {
    x[i] = xx[i];
    v[i] = r > 0.0 ? 2.0 * vv[i] / r : 0.0;
  }
  if (t != NULL) {
    *t = y[9];
  }
}

myb_status myb_ks_rhs(void* problem, double s, const double* y, double* dy) {
  (void)s;
  MYB_CHECK(problem != NULL && y != NULL && dy != NULL);
  const myb_regular_problem* p = problem;
  const double* u = y;
  const double* up = y + 4;
  const double h = y[8];
  const double r = u[0] * u[0] + u[1] * u[1] + u[2] * u[2] + u[3] * u[3];
  double lp[4] = {0.0, 0.0, 0.0, 0.0}; /* L(u)^T P */
  if (p->perturbation != NULL && r > 0.0) {
    double x[3];
    double v[3];
    myb_ks_to_cartesian(y, NULL, x, v);
    double a[3];
    const myb_status st = p->perturbation(p->context, y[9], x, v, a);
    if (st != MYB_OK) {
      return st;
    }
    const double w[4] = {a[0], a[1], a[2], 0.0};
    ks_lt(u, w, lp);
  }
  double hp = 0.0;
  for (int i = 0; i < 4; ++i) {
    dy[i] = up[i];
    dy[4 + i] = -0.5 * h * u[i] + 0.5 * r * lp[i];
    hp += up[i] * lp[i];
  }
  dy[8] = -2.0 * hp;
  dy[9] = r;
  return MYB_OK;
}

myb_status myb_sundman_rhs(void* problem, double s, const double* y, double* dy) {
  (void)s;
  MYB_CHECK(problem != NULL && y != NULL && dy != NULL);
  const myb_regular_problem* p = problem;
  const double r = myb_vec3_norm(y);
  if (!(r > 0.0)) {
    return MYB_E_DOMAIN;
  }
  double g = r; /* r^alpha for the supported exponents 1, 3/2, 2 */
  if (p->alpha == 1.5) {
    g = r * sqrt(r);
  } else if (p->alpha == 2.0) {
    g = r * r;
  } else if (p->alpha != 1.0) {
    return MYB_E_DOMAIN;
  }
  double a[3] = {0.0, 0.0, 0.0};
  if (p->perturbation != NULL) {
    const myb_status st = p->perturbation(p->context, y[6], y, y + 3, a);
    if (st != MYB_OK) {
      return st;
    }
  }
  const double c = p->mu / (r * r * r);
  for (int i = 0; i < 3; ++i) {
    dy[i] = g * y[3 + i];
    dy[3 + i] = g * (a[i] - c * y[i]);
  }
  dy[6] = g;
  return MYB_OK;
}

typedef struct time_event {
  int index; /* of t in y */
  double t1;
} time_event;

static double event_time(void* context, double s, const double* y) {
  (void)s;
  const time_event* e = context;
  return y[e->index] - e->t1;
}

myb_status myb_regular_propagate(myb_regular_kind kind, const myb_regular_problem* problem,
                                 const myb_dop853_options* options, double t0, const double x0[3],
                                 const double v0[3], double t1, double x1[3], double v1[3],
                                 int64_t* evaluations) {
  MYB_CHECK(problem != NULL && options != NULL && x0 != NULL && v0 != NULL && x1 != NULL &&
            v1 != NULL && (kind == MYB_REGULAR_KS || kind == MYB_REGULAR_SUNDMAN));
  if (!(problem->mu > 0.0) || !isfinite(t0) || !isfinite(t1)) {
    return MYB_E_DOMAIN;
  }
  if (evaluations != NULL) {
    *evaluations = 0;
  }
  if (t1 == t0) {
    memcpy(x1, x0, 3 * sizeof *x0);
    memcpy(v1, v0, 3 * sizeof *v0);
    return MYB_OK;
  }
  const int n = kind == MYB_REGULAR_KS ? 10 : 7;
  double y0[10];
  myb_ode_function f = NULL;
  if (kind == MYB_REGULAR_KS) {
    myb_ks_from_cartesian(problem->mu, t0, x0, v0, y0);
    f = myb_ks_rhs;
  } else {
    memcpy(y0, x0, 3 * sizeof *x0);
    memcpy(y0 + 3, v0, 3 * sizeof *v0);
    y0[6] = t0;
    f = myb_sundman_rhs;
  }
  myb_dop853* d = NULL;
  myb_status st = myb_dop853_create((size_t)n, options, NULL, &d);
  if (st != MYB_OK) {
    return st;
  }
  st = myb_dop853_start(d, f, (void*)(uintptr_t)problem, 0.0, y0);
  const double direction = t1 > t0 ? 1.0 : -1.0;
  time_event ev = {.index = n - 1, .t1 = t1};
  const int64_t budget = options->max_steps > 0 ? options->max_steps : INT64_C(1000000000);
  int found = 0;
  double y[10];
  for (int64_t k = 0; k < budget && st == MYB_OK && !found; ++k) {
    st = myb_dop853_step(d, direction * 1e300);
    if (st != MYB_OK) {
      break;
    }
    const double* ys = myb_dop853_state(d);
    if ((ys[n - 1] - t1) * direction >= 0.0) {
      double s_event;
      st = myb_dop853_locate(d, event_time, &ev, &found, &s_event, y);
      if (st == MYB_OK && !found) {
        /* t reached t1 exactly at the end of the step */
        memcpy(y, ys, (size_t)n * sizeof *ys);
        found = 1;
      }
    }
  }
  if (st == MYB_OK && !found) {
    st = MYB_E_NOCONV;
  }
  if (st == MYB_OK) {
    if (kind == MYB_REGULAR_KS) {
      myb_ks_to_cartesian(y, NULL, x1, v1);
    } else {
      memcpy(x1, y, 3 * sizeof *x1);
      memcpy(v1, y + 3, 3 * sizeof *v1);
    }
  }
  if (evaluations != NULL) {
    *evaluations = myb_dop853_statistics(d).evaluations;
  }
  myb_dop853_destroy(d);
  return st;
}
