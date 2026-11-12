/* SPDX-License-Identifier: BSD-3-Clause */
/* Quasi-periodic orbits of the circular restricted three-body problem (myb/cr3bp.h):
 * invariant circles of the flow over a time T, their families, their stability and their
 * manifolds. */
#include <math.h> /* sqrt, fabs, fmax */
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/cr3bp.h"
#include "myb/model.h"
#include "myb/num.h"
#include "myb/uncert.h"
#include "num/linalg_internal.h"

enum {
  MAX_POINTS = MYB_CR3BP_MAX_CIRCLE,
  EXTRA = 4,      /* unknowns after the states: T, rho and the two multipliers */
  NEWTON = 20,
  HALVINGS = 10,
  INVERSE_STEPS = 4
};

static const double two_pi = 6.283185307179586476925286766559;

typedef struct solver {
  double mu;
  myb_force_model* model;
  myb_variational* var;
  size_t points;
  size_t unknowns; /* 6 points + EXTRA */
  myb_cr3bp_torus_family family;
  double kept;           /* the value the family keeps */
  double* turn;          /* points x points: R(-rho) */
  double* turn_rate;     /* its derivative by rho */
  double* by_angle;      /* points x points: d / d theta of the interpolant */
  double* ends;          /* 6 points: the states after T */
  double* phis;          /* 36 points: their transition matrices */
  double* reference;     /* 6 points: the circle the phases are counted from */
  double* along_circle;  /* 6 points: its derivative by the angle */
  double* along_flow;    /* 6 points: the flow on it */
  double* unfold_jacobi; /* 6 points: the gradient of the Jacobi integral on it, unit */
  double* unfold_action; /* 6 points: the gradient of the action of the curve, unit */
  double* predicted;     /* unknowns */
  double* tangent;       /* unknowns */
  double* jacobian;      /* unknowns x unknowns */
  double* f;
  double* f_trial;
  double* z;
  double* z_trial;
  double* step;
} solver;

static void flow_at(double mu, const double y[6], double f[6]) {
  f[0] = y[3];
  f[1] = y[4];
  f[2] = y[5];
  myb_cr3bp_acceleration(mu, y, y + 3, f + 3);
}

/* The gradient of the Jacobi integral C = 2 Omega - |v|^2. */
static void jacobi_gradient(double mu, const double y[6], double g[6]) {
  double a[3];
  myb_cr3bp_acceleration(mu, y, y + 3, a);
  g[0] = 2.0 * (a[0] - 2.0 * y[4]);
  g[1] = 2.0 * (a[1] + 2.0 * y[3]);
  g[2] = 2.0 * a[2];
  g[3] = -2.0 * y[3];
  g[4] = -2.0 * y[4];
  g[5] = -2.0 * y[5];
}

/* W d for the matrix W of the symplectic form in (x, v): the canonical momenta are
 * p = v + (-y, x, 0), so that the flow is P grad E with P = [[0, 1], [-1, 2 B]],
 * B = [[0, 1, 0], [-1, 0, 0], [0, 0, 0]], and W = P^-1 = [[2 B, -1], [1, 0]]. The action of
 * a closed curve u(theta) changes by the integral of du . W u' under a change du. */
static void symplectic(const double d[6], double w[6]) {
  w[0] = 2.0 * d[1] - d[3];
  w[1] = -2.0 * d[0] - d[4];
  w[2] = -d[5];
  w[3] = d[0];
  w[4] = d[1];
  w[5] = d[2];
}

/* The kernel of the trigonometric interpolant on n equal angles, D(x) = (1 / n) sum over
 * |k| <= (n - 1) / 2 of cos(k x), and its derivative. */
static void kernel(size_t n, double x, double* value, double* slope) {
  double sum = 1.0;
  double rate = 0.0;
  for (size_t k = 1; k <= (n - 1) / 2; ++k) {
    sum += 2.0 * myb_cos((double)k * x);
    rate -= 2.0 * (double)k * myb_sin((double)k * x);
  }
  *value = sum / (double)n;
  *slope = rate / (double)n;
}

/* R(-rho): the values at the angles of the interpolant turned back by rho; its derivative
 * by rho when wanted. */
static void turn_matrix(size_t n, double rho, double* turn, double* rate) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double value = 0.0;
      double slope = 0.0;
      kernel(n, two_pi * ((double)i - (double)j) / (double)n - rho, &value, &slope);
      turn[i * n + j] = value;
      if (rate != NULL) {
        rate[i * n + j] = -slope;
      }
    }
  }
}

/* out_i = sum_j m_ij in_j for points of six components. */
static void mix(size_t n, const double* m, const double* in, double* out) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t c = 0; c < 6; ++c) {
      double sum = 0.0;
      for (size_t j = 0; j < n; ++j) {
        sum += m[i * n + j] * in[6 * j + c];
      }
      out[6 * i + c] = sum;
    }
  }
}

static double mean_jacobi(double mu, size_t n, const double* states) {
  double sum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sum += myb_cr3bp_jacobi(mu, states + 6 * i);
  }
  return sum / (double)n;
}

/* The weight of an unknown in the length of a step: the mean over the points for the
 * states, one for T and rho, none for the multipliers. */
static double weight_of(const solver* s, size_t k) {
  const size_t states = 6 * s->points;
  return k < states ? 1.0 / (double)s->points : (k < states + 2 ? 1.0 : 0.0);
}

/* The conditions f for the unknowns z = (states, T, rho, two multipliers), and their
 * Jacobian when wanted. */
static myb_status evaluate(solver* s, const double* z, double* f, double* jacobian) {
  const size_t n = s->points;
  const size_t states = 6 * n;
  const size_t m = s->unknowns;
  const double period = z[states];
  const double rho = z[states + 1];
  if (!(period > 0.0) || !isfinite(period) || !isfinite(rho)) {
    return MYB_E_DOMAIN;
  }
  turn_matrix(n, rho, s->turn, jacobian != NULL ? s->turn_rate : NULL);
  for (size_t j = 0; j < n; ++j) {
    myb_status st = myb_variational_start(s->var, 0.0, z + 6 * j);
    if (st == MYB_OK) {
      st = myb_variational_integrate(s->var, period);
    }
    if (st != MYB_OK) {
      return st;
    }
    myb_variational_get(s->var, s->ends + 6 * j, s->phis + 36 * j, NULL);
  }
  /* the circle after T, turned back, against itself; the multipliers */
  mix(n, s->turn, s->ends, f);
  for (size_t k = 0; k < states; ++k) {
    f[k] += -z[k] + z[states + 2] * s->unfold_jacobi[k] + z[states + 3] * s->unfold_action[k];
  }
  double phase_circle = 0.0;
  double phase_flow = 0.0;
  for (size_t k = 0; k < states; ++k) {
    phase_circle += (z[k] - s->reference[k]) * s->along_circle[k];
    phase_flow += (z[k] - s->reference[k]) * s->along_flow[k];
  }
  f[states] = phase_circle / (double)n;
  f[states + 1] = phase_flow / (double)n;
  f[states + 2] = s->family == MYB_CR3BP_TORUS_JACOBI
                      ? mean_jacobi(s->mu, n, z) - s->kept
                      : (s->family == MYB_CR3BP_TORUS_PERIOD ? period - s->kept : rho - s->kept);
  double arc = 0.0;
  for (size_t k = 0; k < m; ++k) {
    arc += weight_of(s, k) * (z[k] - s->predicted[k]) * s->tangent[k];
  }
  f[states + 3] = arc;
  if (jacobian == NULL) {
    return MYB_OK;
  }
  memset(jacobian, 0, m * m * sizeof *jacobian);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      const double r = s->turn[i * n + j];
      const double dr = s->turn_rate[i * n + j];
      const double* phi = s->phis + 36 * j;
      double rate[6];
      flow_at(s->mu, s->ends + 6 * j, rate);
      for (size_t a = 0; a < 6; ++a) {
        double* row = jacobian + (6 * i + a) * m;
        for (size_t b = 0; b < 6; ++b) {
          row[6 * j + b] = r * phi[6 * a + b];
        }
        row[states] += r * rate[a];
        row[states + 1] += dr * s->ends[6 * j + a];
      }
    }
  }
  for (size_t k = 0; k < states; ++k) {
    double* row = jacobian + k * m;
    row[k] -= 1.0;
    row[states + 2] = s->unfold_jacobi[k];
    row[states + 3] = s->unfold_action[k];
    jacobian[states * m + k] = s->along_circle[k] / (double)n;
    jacobian[(states + 1) * m + k] = s->along_flow[k] / (double)n;
  }
  double* family_row = jacobian + (states + 2) * m;
  if (s->family == MYB_CR3BP_TORUS_JACOBI) {
    for (size_t i = 0; i < n; ++i) {
      double g[6];
      jacobi_gradient(s->mu, z + 6 * i, g);
      for (size_t c = 0; c < 6; ++c) {
        family_row[6 * i + c] = g[c] / (double)n;
      }
    }
  } else {
    family_row[s->family == MYB_CR3BP_TORUS_PERIOD ? states : states + 1] = 1.0;
  }
  for (size_t k = 0; k < m; ++k) {
    jacobian[(states + 3) * m + k] = weight_of(s, k) * s->tangent[k];
  }
  return MYB_OK;
}

static double largest(const double* f, size_t m) {
  double worst = 0.0;
  for (size_t k = 0; k < m; ++k) {
    worst = fmax(worst, fabs(f[k]));
  }
  return worst;
}

/* A vector of six components per point made unit. */
static void make_unit(size_t count, double* v) {
  double sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    sum += v[k] * v[k];
  }
  const double norm = sqrt(sum);
  for (size_t k = 0; k < count; ++k) {
    v[k] = norm > 0.0 ? v[k] / norm : 0.0;
  }
}

/* The force model of the problem (its context the double *mu, which must outlive it) and
 * the variational equations on it. */
static myb_status flow_open(double* mu, myb_force_model** model, myb_variational** var) {
  myb_status st = myb_force_model_create(NULL, model);
  if (st == MYB_OK) {
    st = myb_force_model_add(*model, "cr3bp", myb_cr3bp_term, myb_cr3bp_term_jacobian, mu, NULL);
  }
  const myb_variational_options o = {.epsilon = 1e-12};
  if (st == MYB_OK) {
    st = myb_variational_create(*model, &o, NULL, var);
  }
  return st;
}

static myb_status solver_open(solver* s, double mu, size_t points, void** block) {
  memset(s, 0, sizeof *s);
  *block = NULL;
  s->mu = mu;
  s->points = points;
  s->unknowns = 6 * points + EXTRA;
  myb_status st = flow_open(&s->mu, &s->model, &s->var);
  const size_t n = points;
  const size_t m = s->unknowns;
  const size_t doubles = 3 * n * n + 6 * n + 36 * n + 5 * 6 * n + 7 * m + m * m;
  const myb_allocator a = myb_internal_allocator(NULL);
  if (st == MYB_OK) {
    st = myb_internal_alloc_array(&a, doubles, sizeof(double), block);
  }
  if (st != MYB_OK) {
    return st;
  }
  double* w = *block;
  s->turn = w;
  s->turn_rate = s->turn + n * n;
  s->by_angle = s->turn_rate + n * n;
  s->ends = s->by_angle + n * n;
  s->phis = s->ends + 6 * n;
  s->reference = s->phis + 36 * n;
  s->along_circle = s->reference + 6 * n;
  s->along_flow = s->along_circle + 6 * n;
  s->unfold_jacobi = s->along_flow + 6 * n;
  s->unfold_action = s->unfold_jacobi + 6 * n;
  s->predicted = s->unfold_action + 6 * n;
  s->tangent = s->predicted + m;
  s->f = s->tangent + m;
  s->f_trial = s->f + m;
  s->z = s->f_trial + m;
  s->z_trial = s->z + m;
  s->step = s->z_trial + m;
  s->jacobian = s->step + m;
  return MYB_OK;
}

static void solver_close(solver* s, void* block) {
  const myb_allocator a = myb_internal_allocator(NULL);
  myb_internal_release(&a, block);
  myb_variational_destroy(s->var);
  myb_force_model_destroy(s->model);
}

static int torus_valid(double mu, const myb_cr3bp_torus* torus) {
  if (!(mu > 0.0 && mu <= 0.5) || torus->points < 5 || torus->points > MAX_POINTS ||
      torus->points % 2 == 0 || !(torus->period > 0.0) || !isfinite(torus->period) ||
      !isfinite(torus->rotation)) {
    return 0;
  }
  for (size_t k = 0; k < 6 * torus->points; ++k) {
    if (!isfinite(torus->states[k])) {
      return 0;
    }
  }
  return 1;
}

myb_status myb_cr3bp_torus_continue(double mu, myb_cr3bp_torus* torus,
                                    myb_cr3bp_torus_family family, double step,
                                    double tolerance, int max_iterations) {
  MYB_CHECK(torus != NULL);
  if (!torus_valid(mu, torus) || !torus->has_tangent || !isfinite(step) ||
      !(tolerance >= 0.0) || max_iterations < 0 ||
      (family != MYB_CR3BP_TORUS_JACOBI && family != MYB_CR3BP_TORUS_PERIOD &&
       family != MYB_CR3BP_TORUS_ROTATION)) {
    return MYB_E_DOMAIN;
  }
  const double limit = tolerance > 0.0 ? tolerance : 1e-10;
  const int steps = max_iterations > 0 ? max_iterations : NEWTON;
  solver s;
  void* block = NULL;
  myb_status st = solver_open(&s, mu, torus->points, &block);
  if (st != MYB_OK) {
    solver_close(&s, block);
    return st;
  }
  const size_t n = s.points;
  const size_t states = 6 * n;
  const size_t m = s.unknowns;
  s.family = family;
  s.kept = family == MYB_CR3BP_TORUS_JACOBI
               ? mean_jacobi(mu, n, torus->states)
               : (family == MYB_CR3BP_TORUS_PERIOD ? torus->period : torus->rotation);
  /* the prediction along the tangent */
  memcpy(s.tangent, torus->tangent, (states + 2) * sizeof *s.tangent);
  s.tangent[states + 2] = s.tangent[states + 3] = 0.0;
  for (size_t k = 0; k < states; ++k) {
    s.predicted[k] = torus->states[k] + step * s.tangent[k];
  }
  s.predicted[states] = torus->period + step * s.tangent[states];
  s.predicted[states + 1] = torus->rotation + step * s.tangent[states + 1];
  s.predicted[states + 2] = s.predicted[states + 3] = 0.0;
  memcpy(s.z, s.predicted, m * sizeof *s.z);
  /* what the phases and the multipliers are taken from: the predicted circle */
  memcpy(s.reference, s.predicted, states * sizeof *s.reference);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double value = 0.0;
      double slope = 0.0;
      kernel(n, two_pi * ((double)i - (double)j) / (double)n, &value, &slope);
      s.by_angle[i * n + j] = slope;
    }
  }
  mix(n, s.by_angle, s.reference, s.along_circle);
  for (size_t i = 0; i < n; ++i) {
    flow_at(mu, s.reference + 6 * i, s.along_flow + 6 * i);
    jacobi_gradient(mu, s.reference + 6 * i, s.unfold_jacobi + 6 * i);
    symplectic(s.along_circle + 6 * i, s.unfold_action + 6 * i);
  }
  make_unit(states, s.unfold_jacobi);
  make_unit(states, s.unfold_action);
  /* Newton's method with the step halved while the conditions do not come closer */
  st = evaluate(&s, s.z, s.f, s.jacobian);
  double norm = st == MYB_OK ? largest(s.f, m) : 0.0;
  /* bound: `steps` iterations */
  for (int it = 0; it < steps && st == MYB_OK && norm > limit; ++it) {
    for (size_t k = 0; k < m; ++k) {
      s.step[k] = -s.f[k];
    }
    if (myb_internal_lu_solve(m, s.jacobian, s.step, 1) != MYB_OK) {
      st = MYB_E_NOCONV;
      break;
    }
    double scale = 1.0;
    int accepted = 0;
    /* bound: HALVINGS */
    for (int half = 0; half <= HALVINGS && !accepted; ++half) {
      for (size_t k = 0; k < m; ++k) {
        s.z_trial[k] = s.z[k] + scale * s.step[k];
      }
      /* with the Jacobian, which the solution of the step has used up: the transition
       * matrices come with the states, and a step that stands needs it next */
      accepted = evaluate(&s, s.z_trial, s.f_trial, s.jacobian) == MYB_OK &&
                 largest(s.f_trial, m) < norm;
      scale *= 0.5;
    }
    if (!accepted) {
      st = MYB_E_NOCONV;
      break;
    }
    memcpy(s.z, s.z_trial, m * sizeof *s.z);
    memcpy(s.f, s.f_trial, m * sizeof *s.f);
    norm = largest(s.f, m);
  }
  if (st == MYB_OK && !(norm <= limit)) {
    st = MYB_E_NOCONV;
  }
  if (st == MYB_OK) {
    /* the tangent of the family: across every condition but the last, whose row is the
     * tangent before */
    memset(s.step, 0, m * sizeof *s.step);
    s.step[m - 1] = 1.0;
    double length = 0.0;
    if (myb_internal_lu_solve(m, s.jacobian, s.step, 1) == MYB_OK) {
      for (size_t k = 0; k < m; ++k) {
        length += weight_of(&s, k) * s.step[k] * s.step[k];
      }
      length = sqrt(length);
    }
    if (!(length > 0.0) || !isfinite(length)) {
      st = MYB_E_NOCONV;
    } else {
      memcpy(torus->states, s.z, states * sizeof *s.z);
      torus->period = s.z[states];
      torus->rotation = s.z[states + 1];
      torus->defect = fmax(fabs(s.z[states + 2]), fabs(s.z[states + 3]));
      /* the last condition makes the new tangent point the way of the old one */
      for (size_t k = 0; k < states + 2; ++k) {
        torus->tangent[k] = s.step[k] / length;
      }
    }
  }
  solver_close(&s, block);
  return st;
}

myb_status myb_cr3bp_torus_from_orbit(double mu, const double y0[6], double period,
                                      size_t points, int pair, double size,
                                      myb_cr3bp_torus* torus) {
  MYB_CHECK(y0 != NULL && torus != NULL);
  if (!(mu > 0.0 && mu <= 0.5) || !(period > 0.0) || !isfinite(period) || points < 5 ||
      points > MAX_POINTS || points % 2 == 0 || pair < 0 || !(size > 0.0) || !isfinite(size)) {
    return MYB_E_DOMAIN;
  }
  double mono[36];
  myb_status st = myb_cr3bp_monodromy(mu, y0, period, mono);
  if (st != MYB_OK) {
    return st;
  }
  double copy[36];
  double re[6];
  double im[6];
  memcpy(copy, mono, sizeof copy);
  st = myb_eigenvalues(6, copy, re, im);
  if (st != MYB_OK) {
    return st;
  }
  /* the pairs on the unit circle away from 1, in the order of their angle */
  double angles[3];
  int count = 0;
  for (int k = 0; k < 6; ++k) {
    const double modulus = sqrt(re[k] * re[k] + im[k] * im[k]);
    const double angle = myb_atan2(im[k], re[k]);
    if (im[k] > 0.0 && fabs(modulus - 1.0) < 1e-5 && angle > 1e-3 && count < 3) {
      int at = count++;
      /* bound: at most two shifts */
      for (; at > 0 && angles[at - 1] > angle; --at) {
        angles[at] = angles[at - 1];
      }
      angles[at] = angle;
    }
  }
  if (pair >= count) {
    return MYB_E_DOMAIN;
  }
  const double rho = angles[pair];
  /* a + i b by inverse iteration on the real form of M - exp(i rho'), rho' beside rho */
  const double beside = rho + 1e-7;
  const double c = myb_cos(beside);
  const double sn = myb_sin(beside);
  double x[12];
  for (int k = 0; k < 12; ++k) {
    x[k] = 1.0 / (1.0 + (double)((7 * k) % 12)); /* any vector with a part along the pair */
  }
  /* bound: INVERSE_STEPS */
  for (int it = 0; it < INVERSE_STEPS; ++it) {
    double big[144] = {0};
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        const double entry = mono[6 * i + j] - (i == j ? c : 0.0);
        big[12 * i + j] = entry;
        big[12 * (6 + i) + 6 + j] = entry;
      }
      big[12 * i + 6 + i] = sn;
      big[12 * (6 + i) + i] = -sn;
    }
    if (myb_lu_solve(12, big, x, 1) != MYB_OK) {
      return MYB_E_NOCONV;
    }
    double norm = 0.0;
    for (int k = 0; k < 12; ++k) {
      norm += x[k] * x[k];
    }
    norm = sqrt(0.5 * norm); /* (|a|^2 + |b|^2) / 2 = 1: a circle of unit mean square */
    if (!(norm > 0.0) || !isfinite(norm)) {
      return MYB_E_NOCONV;
    }
    for (int k = 0; k < 12; ++k) {
      x[k] /= norm;
    }
  }
  const double* a = x;
  const double* b = x + 6;
  /* M a = cos rho a - sin rho b must hold */
  double miss = 0.0;
  for (int i = 0; i < 6; ++i) {
    double ma = 0.0;
    for (int j = 0; j < 6; ++j) {
      ma += mono[6 * i + j] * a[j];
    }
    miss = fmax(miss, fabs(ma - (myb_cos(rho) * a[i] - myb_sin(rho) * b[i])));
  }
  if (!(miss < 1e-6)) {
    return MYB_E_NOCONV;
  }
  memset(torus, 0, sizeof *torus);
  torus->points = points;
  torus->period = period;
  torus->rotation = rho;
  for (size_t i = 0; i < points; ++i) {
    const double theta = two_pi * (double)i / (double)points;
    const double ct = myb_cos(theta);
    const double st_ = myb_sin(theta);
    for (size_t k = 0; k < 6; ++k) {
      const double out = ct * a[k] - st_ * b[k];
      torus->states[6 * i + k] = y0[k] + size * out;
      torus->tangent[6 * i + k] = out;
    }
  }
  torus->has_tangent = 1;
  return MYB_OK;
}

/* The trigonometric interpolant of values of six components at the n points, at an angle. */
static void interpolate(size_t n, const double* values, double angle, double y[6]) {
  memset(y, 0, 6 * sizeof *y);
  for (size_t j = 0; j < n; ++j) {
    double value = 0.0;
    double slope = 0.0;
    kernel(n, angle - two_pi * (double)j / (double)n, &value, &slope);
    for (size_t k = 0; k < 6; ++k) {
      y[k] += value * values[6 * j + k];
    }
  }
}

void myb_cr3bp_torus_state(const myb_cr3bp_torus* torus, double angle, double y[6]) {
  MYB_CHECK(torus != NULL && y != NULL && torus->points >= 1 && torus->points <= MAX_POINTS &&
            torus->points % 2 == 1);
  interpolate(torus->points, torus->states, angle, y);
}

/* ---- Stability and manifolds -------------------------------------------------------------------
 * The linearised map A = R(-rho) diag(Phi_j) of myb/cr3bp.h is assembled from the transition
 * matrices of the points over T (the variational equations of the corrector, IAS15 at 1e-12)
 * and its 6 n eigenvalues found by the QR iteration of myb/num.h in its form without the
 * limit on the order (n = 35: a matrix of 210, 0.35 MB, two copies). The eigenvalues are
 * sorted by modulus, from the largest.
 *
 * The classes come from the circles themselves (classify below): an eigenvalue whose turn by
 * +-rho is not in the spectrum, within the tolerance, is one the points do not resolve.
 * Measured: about the L1 Lyapunov orbit (25 points, sizes 1e-3 .. 4e-3) the members of a
 * circle away from its ends are turns of each other to 1e-9 .. 1e-6, those at the ends and
 * the spurious ones to 1e-5 .. 3e-3; the Jordan blocks of 1 split radially by 2.6e-6.
 *
 * The hyperbolic pair: of the resolved real eigenvalues off the unit circle, on each side the
 * farthest from it — the k = 0 member of its circle, the only real one since n is odd (but
 * for a turn k rho of pi, which makes -lambda real too: then the one whose eigenvector, by
 * inverse iteration on A and taken to the Fourier coefficients of its six components, has
 * the lowest mean frequency |k|, weighted by the energy of the harmonics). A resolved circle
 * off the unit circle with no real member on its side is no spectrum of a torus of the
 * problem (MYB_E_NOCONV). */

static const double default_tolerance = 1e-4; /* on the unit circle */

typedef struct spectrum_work {
  double* turn;        /* n x n: R(-rho) */
  double* phis;        /* 36 n: the transition matrices of the points over T */
  double* map;         /* 6 n x 6 n: A */
  double* work;        /* 6 n x 6 n: the copy the eigenvalues destroy, then inverse iteration */
  double* trial;       /* 6 n: a candidate eigenvector (and the reflectors of the QR iteration) */
  double* image;       /* 6 n: A v */
  unsigned char* kind; /* 6 n: the classes of the eigenvalues */
  size_t* parent;      /* 6 n: their sets */
  size_t* size;        /* 6 n */
} spectrum_work;

/* The Fourier content of a real function of the angle given at the n points (six components
 * each): the mean frequency |k| weighted by the energy of the harmonics, and the share of the
 * energy above half the highest frequency (n - 1) / 2. */
static void frequency_content(size_t n, const double* v, double* mean, double* tail) {
  double c[MAX_POINTS];
  double s[MAX_POINTS];
  for (size_t j = 0; j < n; ++j) {
    c[j] = myb_cos(two_pi * (double)j / (double)n);
    s[j] = myb_sin(two_pi * (double)j / (double)n);
  }
  const size_t top = (n - 1) / 2;
  double total = 0.0;
  double weighted = 0.0;
  double upper = 0.0;
  for (size_t k = 0; k <= top; ++k) {
    double energy = 0.0;
    for (size_t a = 0; a < 6; ++a) {
      double re = 0.0;
      double im = 0.0;
      for (size_t j = 0; j < n; ++j) {
        const size_t at = (k * j) % n; /* exp(-i k theta_j) */
        re += v[6 * j + a] * c[at];
        im -= v[6 * j + a] * s[at];
      }
      energy += re * re + im * im;
    }
    if (k > 0) {
      energy *= 2.0; /* with the harmonic -k of a real function */
    }
    total += energy;
    weighted += (double)k * energy;
    if (2 * k > top) {
      upper += energy;
    }
  }
  *mean = total > 0.0 ? weighted / total : 0.0;
  *tail = total > 0.0 ? upper / total : 0.0;
}

/* The direction turned to a positive mean x component and scaled to unit root mean square
 * over the points; |A v - lambda v| / (|lambda| |v|). */
static double finish_direction(size_t n, const double* map, double lambda, double* v,
                               double* image) {
  const size_t m = 6 * n;
  double along_x = 0.0;
  double sum = 0.0;
  for (size_t j = 0; j < n; ++j) {
    along_x += v[6 * j];
  }
  for (size_t k = 0; k < m; ++k) {
    sum += v[k] * v[k];
  }
  const double scale = (along_x < 0.0 ? -1.0 : 1.0) * sqrt((double)n / sum);
  for (size_t k = 0; k < m; ++k) {
    v[k] *= scale;
  }
  double miss = 0.0;
  double size = 0.0;
  for (size_t r = 0; r < m; ++r) {
    double s = 0.0;
    for (size_t c = 0; c < m; ++c) {
      s += map[r * m + c] * v[c];
    }
    image[r] = s - lambda * v[r];
    miss += image[r] * image[r];
    size += v[r] * v[r];
  }
  return sqrt(miss / size) / fabs(lambda);
}

/* Whether the eigenvalue f is the eigenvalue e turned by +rho or -rho (cos rho = c, sin
 * rho = sn), within near |lambda_e|. */
static int is_turn(const myb_cr3bp_torus_spectrum* s, size_t e, size_t f, double c, double sn,
                   double near) {
  const double re = s->re[e];
  const double im = s->im[e];
  const double limit = near * myb_hypot(re, im);
  return myb_hypot(s->re[f] - (re * c - im * sn), s->im[f] - (re * sn + im * c)) <= limit ||
         myb_hypot(s->re[f] - (re * c + im * sn), s->im[f] - (im * c - re * sn)) <= limit;
}

enum { KIND_UNRESOLVED = 0, KIND_OFF = 1, KIND_ON = 2, KIND_TRIVIAL = 3 };

/* The root of the set of e (union-find with path halving). */
static size_t root_of(size_t* parent, size_t e) {
  /* bound: the path halves at every step; at most as many steps as eigenvalues */
  while (parent[e] != e) {
    parent[e] = parent[parent[e]];
    e = parent[e];
  }
  return e;
}

/* The classes of the eigenvalues, from their circles: two eigenvalues are joined when one is
 * the other turned by +-rho within the tolerance, and the sets so joined are the circles —
 * near their ends the members drift along the circle together, each still a turn of its
 * neighbour. A set of more than half the points is resolved (the four trivial circles make
 * one set, the others one each); the smaller ones are the spurious eigenvalues of high
 * frequency and the loose ends. A resolved set holding 1 is trivial; the others are on the
 * unit circle (elliptic) or off it (hyperbolic). parent, size: count values each. */
static void classify(myb_cr3bp_torus_spectrum* s, size_t points, double rho, double near,
                     unsigned char* kind, size_t* parent, size_t* size) {
  const size_t m = s->count;
  const double c = myb_cos(rho);
  const double sn = myb_sin(rho);
  for (size_t e = 0; e < m; ++e) {
    parent[e] = e;
    size[e] = 0;
    kind[e] = KIND_ON;
  }
  for (size_t e = 0; e < m; ++e) {
    for (size_t f = e + 1; f < m; ++f) {
      if (is_turn(s, e, f, c, sn, near)) {
        const size_t a = root_of(parent, e);
        const size_t b = root_of(parent, f);
        if (a != b) {
          parent[a > b ? a : b] = a > b ? b : a;
        }
      }
    }
  }
  for (size_t e = 0; e < m; ++e) {
    const size_t r = root_of(parent, e);
    ++size[r];
    if (myb_hypot(s->re[e] - 1.0, s->im[e]) <= near) {
      kind[r] = KIND_TRIVIAL; /* the mark of the set, at its root */
    }
  }
  /* at the roots, the code of the set: 0 unresolved, 1 resolved, 2 resolved and trivial */
  for (size_t r = 0; r < m; ++r) {
    if (parent[r] == r) {
      size[r] = 2 * size[r] > points ? (kind[r] == KIND_TRIVIAL ? 2u : 1u) : 0u;
    }
  }
  s->trivial = s->elliptic = s->hyperbolic = s->unresolved = 0;
  for (size_t e = 0; e < m; ++e) {
    const size_t code = size[root_of(parent, e)];
    const double modulus = myb_hypot(s->re[e], s->im[e]);
    unsigned char k = KIND_UNRESOLVED;
    if (code == 2) {
      k = KIND_TRIVIAL;
    } else if (code == 1) {
      k = fabs(modulus - 1.0) <= near ? KIND_ON : KIND_OFF;
    }
    kind[e] = k;
    s->unresolved += k == KIND_UNRESOLVED;
    s->hyperbolic += k == KIND_OFF;
    s->elliptic += k == KIND_ON;
    s->trivial += k == KIND_TRIVIAL;
  }
}

static myb_status spectrum_compute(myb_variational* var, const myb_cr3bp_torus* torus, double near,
                                   const spectrum_work* w, myb_cr3bp_torus_spectrum* out) {
  const size_t n = torus->points;
  const size_t m = 6 * n;
  for (size_t j = 0; j < n; ++j) {
    myb_status st = myb_variational_start(var, 0.0, torus->states + 6 * j);
    if (st == MYB_OK) {
      st = myb_variational_integrate(var, torus->period);
    }
    if (st != MYB_OK) {
      return st;
    }
    myb_variational_get(var, NULL, w->phis + 36 * j, NULL);
  }
  turn_matrix(n, torus->rotation, w->turn, NULL);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      const double r = w->turn[i * n + j];
      const double* phi = w->phis + 36 * j;
      for (size_t a = 0; a < 6; ++a) {
        double* row = w->map + (6 * i + a) * m + 6 * j;
        for (size_t b = 0; b < 6; ++b) {
          row[b] = r * phi[6 * a + b];
        }
      }
    }
  }
  memcpy(w->work, w->map, m * m * sizeof *w->work);
  myb_status st = myb_internal_eigenvalues(m, w->work, out->re, out->im, w->trial);
  if (st != MYB_OK) {
    return st;
  }
  out->count = m;
  /* by modulus, from the largest (insertion: m^2 / 2 comparisons at most) */
  for (size_t k = 1; k < m; ++k) {
    const double re = out->re[k];
    const double im = out->im[k];
    const double modulus = myb_hypot(re, im);
    size_t at = k;
    /* bound: at most k shifts */
    for (; at > 0 && myb_hypot(out->re[at - 1], out->im[at - 1]) < modulus; --at) {
      out->re[at] = out->re[at - 1];
      out->im[at] = out->im[at - 1];
    }
    out->re[at] = re;
    out->im[at] = im;
  }
  classify(out, n, torus->rotation, near, w->kind, w->parent, w->size);
  int off_above = 0;
  int off_below = 0;
  for (size_t e = 0; e < m; ++e) {
    if (w->kind[e] == KIND_OFF) {
      off_above = off_above || myb_hypot(out->re[e], out->im[e]) > 1.0;
      off_below = off_below || myb_hypot(out->re[e], out->im[e]) < 1.0;
    }
  }
  /* the representatives of the hyperbolic pair: real and resolved, the farthest from the
   * unit circle on each side; of two of one modulus (lambda and -lambda, when a turn k rho
   * is pi) the one of the lowest frequency */
  double lowest[2] = {HUGE_VAL, HUGE_VAL}; /* [0]: stable, [1]: unstable */
  double tails[2] = {0.0, 0.0};
  double lambdas[2] = {0.0, 0.0};
  double* directions[2] = {out->stable_direction, out->unstable_direction};
  memset(out->stable_direction, 0, sizeof out->stable_direction);
  memset(out->unstable_direction, 0, sizeof out->unstable_direction);
  for (size_t e = 0; e < m; ++e) {
    if (out->im[e] != 0.0 || w->kind[e] != KIND_OFF) {
      continue;
    }
    const double modulus = fabs(out->re[e]);
    const int side = modulus > 1.0 ? 1 : 0;
    const double best = fabs(lambdas[side]);
    const double margin = near * fmax(modulus, best);
    const int farther =
        lambdas[side] == 0.0 || (side == 1 ? modulus > best + margin : modulus < best - margin);
    const int level = !farther && fabs(modulus - best) <= margin;
    if (!farther && !level) {
      continue;
    }
    if (myb_internal_eigenvector_real(m, w->map, out->re[e], w->trial, w->work) != MYB_OK) {
      continue; /* not an eigenvalue inverse iteration can resolve */
    }
    double mean = 0.0;
    double tail = 0.0;
    frequency_content(n, w->trial, &mean, &tail);
    if (farther || mean < lowest[side]) {
      lowest[side] = mean;
      tails[side] = tail;
      lambdas[side] = out->re[e];
      memcpy(directions[side], w->trial, m * sizeof *w->trial);
    }
  }
  if ((off_above && lambdas[1] == 0.0) || (off_below && lambdas[0] == 0.0)) {
    return MYB_E_NOCONV;
  }
  out->stable = lambdas[0];
  out->unstable = lambdas[1];
  out->residual = 0.0;
  out->tail = fmax(tails[0], tails[1]);
  for (int side = 0; side < 2; ++side) {
    if (lambdas[side] != 0.0) {
      out->residual = fmax(out->residual,
                           finish_direction(n, w->map, lambdas[side], directions[side], w->image));
    }
  }
  const double big = fabs(out->unstable);
  out->index = big > 0.0 ? 0.5 * (big + 1.0 / big) : 1.0;
  return MYB_OK;
}

myb_status myb_cr3bp_torus_stability(double mu, const myb_cr3bp_torus* torus, double tolerance,
                                     myb_cr3bp_torus_spectrum* spectrum) {
  MYB_CHECK(torus != NULL && spectrum != NULL);
  if (!torus_valid(mu, torus) || !(tolerance >= 0.0) || !isfinite(tolerance)) {
    return MYB_E_DOMAIN;
  }
  const double near = tolerance > 0.0 ? tolerance : default_tolerance;
  /* the neighbours of a member of a circle must stand apart from it */
  if (!(fabs(2.0 * myb_sin(0.5 * torus->rotation)) > 10.0 * near)) {
    return MYB_E_DOMAIN;
  }
  const size_t n = torus->points;
  const size_t m = 6 * n;
  /* doubles, then the indices of the sets, then the classes */
  const size_t indices = (2 * m * sizeof(size_t) + sizeof(double) - 1) / sizeof(double);
  double context = mu;
  myb_force_model* model = NULL;
  myb_variational* var = NULL;
  void* block = NULL;
  const myb_allocator a = myb_internal_allocator(NULL);
  myb_status st = flow_open(&context, &model, &var);
  if (st == MYB_OK) {
    st = myb_internal_alloc_array(&a, n * n + 36 * n + 2 * m * m + 2 * m + indices + (m + 7) / 8,
                                  sizeof(double), &block);
  }
  if (st == MYB_OK) {
    spectrum_work w;
    w.turn = block;
    w.phis = w.turn + n * n;
    w.map = w.phis + 36 * n;
    w.work = w.map + m * m;
    w.trial = w.work + m * m;
    w.image = w.trial + m;
    w.parent = (size_t*)(void*)(w.image + m);
    w.size = w.parent + m;
    w.kind = (unsigned char*)(w.image + m + indices);
    st = spectrum_compute(var, torus, near, &w, spectrum);
  }
  myb_internal_release(&a, block);
  myb_variational_destroy(var);
  myb_force_model_destroy(model);
  return st;
}

myb_status myb_cr3bp_torus_manifold(const myb_cr3bp_torus* torus,
                                    const myb_cr3bp_torus_spectrum* spectrum, double angle,
                                    int unstable, int sign, double epsilon, double y[6]) {
  MYB_CHECK(torus != NULL && spectrum != NULL && y != NULL);
  const size_t n = torus->points;
  if (n < 5 || n > MAX_POINTS || n % 2 == 0 || spectrum->count != 6 * n || !isfinite(angle) ||
      !(epsilon > 0.0) || !isfinite(epsilon) || (sign != 1 && sign != -1)) {
    return MYB_E_DOMAIN;
  }
  const double lambda = unstable ? spectrum->unstable : spectrum->stable;
  if (!(fabs(lambda) > 0.0) || !isfinite(lambda)) {
    return MYB_E_DOMAIN;
  }
  double base[6];
  double d[6];
  interpolate(n, torus->states, angle, base);
  interpolate(n, unstable ? spectrum->unstable_direction : spectrum->stable_direction, angle, d);
  const double norm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  if (!(norm > 0.0) || !isfinite(norm)) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < 6; ++k) {
    y[k] = base[k] + (double)sign * epsilon * d[k] / norm;
  }
  return MYB_OK;
}
