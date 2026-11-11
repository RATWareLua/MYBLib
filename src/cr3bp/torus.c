/* SPDX-License-Identifier: BSD-3-Clause */
/* Quasi-periodic orbits of the circular restricted three-body problem (myb/cr3bp.h):
 * invariant circles of the flow over a time T, and their families. */
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

static myb_status solver_open(solver* s, double mu, size_t points, void** block) {
  memset(s, 0, sizeof *s);
  *block = NULL;
  s->mu = mu;
  s->points = points;
  s->unknowns = 6 * points + EXTRA;
  myb_status st = myb_force_model_create(NULL, &s->model);
  if (st == MYB_OK) {
    st = myb_force_model_add(s->model, "cr3bp", myb_cr3bp_term, myb_cr3bp_term_jacobian, &s->mu,
                             NULL);
  }
  const myb_variational_options o = {.epsilon = 1e-12};
  if (st == MYB_OK) {
    st = myb_variational_create(s->model, &o, NULL, &s->var);
  }
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

void myb_cr3bp_torus_state(const myb_cr3bp_torus* torus, double angle, double y[6]) {
  MYB_CHECK(torus != NULL && y != NULL && torus->points >= 1 && torus->points <= MAX_POINTS &&
            torus->points % 2 == 1);
  const size_t n = torus->points;
  memset(y, 0, 6 * sizeof *y);
  for (size_t j = 0; j < n; ++j) {
    double value = 0.0;
    double slope = 0.0;
    kernel(n, angle - two_pi * (double)j / (double)n, &value, &slope);
    for (size_t k = 0; k < 6; ++k) {
      y[k] += value * torus->states[6 * j + k];
    }
  }
}
