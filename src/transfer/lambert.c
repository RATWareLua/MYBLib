/* SPDX-License-Identifier: BSD-3-Clause */
/* Lambert's problem after D. Izzo (2015). Notation of the paper:
 *   c = |r2 - r1|, s = (r1 + r2 + c) / 2, lambda^2 = 1 - c/s (sign by geometry and direction),
 *   T = sqrt(2 mu / s^3) tof,  y(x) = sqrt(1 - lambda^2 (1 - x^2)),
 *   T(x) = [ (psi + M pi) / sqrt|1 - x^2| - x + lambda y ] / (1 - x^2)          (eq. 18)
 *     cos psi = x y + lambda (1 - x^2) (ellipse), sinh psi = (y - x lambda) sqrt(x^2 - 1)
 * (hyperbola); near x = 1 (M = 0): eta = y - lambda x, S1 = (1 - lambda - x eta)/2, T = (eta^3 Q +
 * 4 lambda eta) / 2, Q = 4/3 2F1(3, 1; 5/2; S1)                (eq. 20) T'  = (3 T x - 2 + 2
 * lambda^3 x / y) / (1 - x^2)                              (eq. 22) T'' = (3 T + 5 x T' + 2 (1 -
 * lambda^2) lambda^3 / y^3) / (1 - x^2) velocities: gamma = sqrt(mu s / 2), rho = (r1 - r2)/c,
 * sigma = sqrt(1 - rho^2), Vr1 = gamma [(lambda y - x) - rho (lambda y + x)] / r1, Vr2 = -gamma
 * [(lambda y - x) + rho (lambda y + x)] / r2, Vt1 = gamma sigma (y + lambda x) / r1,   Vt2 = gamma
 * sigma (y + lambda x) / r2. */
#include "myb/lambert.h"

#include <math.h> /* sqrt, fabs */

#include "core/base.h"
#include "myb/cmath.h"
#include "num/root.h"
#include "num/vec3.h"

enum { MAX_ITERATIONS = 200, MAX_SERIES_TERMS = 400, MAX_DOUBLINGS = 200 };

static const double pi = 0x1.921fb54442d18p+1;

static int finite(double x) {
  return x - x == 0.0;
}

typedef struct time_equation {
  double lambda;
  int revolutions;
} time_equation;

/* Gauss hypergeometric 2F1(3, 1; 5/2; z) by its series (|z| < 1 where it is used). */
static double hyp2f1b(double z) {
  double sum = 1.0;
  double term = 1.0;
  /* bound: MAX_SERIES_TERMS; terms decrease geometrically like z^k for |z| < 1 */
  for (int k = 0; k < MAX_SERIES_TERMS; ++k) {
    term *= (3.0 + k) / (2.5 + k) * z;
    const double next = sum + term;
    if (next == sum) {
      break;
    }
    sum = next;
  }
  return sum;
}

static double y_of(const time_equation* q, double x) {
  return sqrt(1.0 - q->lambda * q->lambda * (1.0 - x * x));
}

/* Non-dimensional time of flight T(x). */
static double time_of_flight(const time_equation* q, double x) {
  const double l = q->lambda;
  const double y = y_of(q, x);
  if (q->revolutions == 0 && x > sqrt(0.6) && x < sqrt(1.4)) {
    const double eta = y - l * x;
    const double s1 = 0.5 * (1.0 - l - x * eta);
    const double big_q = 4.0 / 3.0 * hyp2f1b(s1);
    return 0.5 * (eta * eta * eta * big_q + 4.0 * l * eta);
  }
  const double one_minus_x2 = 1.0 - x * x;
  double psi = 0.0;
  if (x < 1.0) {
    psi = myb_acos(x * y + l * one_minus_x2);
  } else {
    psi = myb_asinh((y - x * l) * sqrt(x * x - 1.0));
  }
  return ((psi + q->revolutions * pi) / sqrt(fabs(one_minus_x2)) - x + l * y) / one_minus_x2;
}

static double first_derivative(const time_equation* q, double x, double t) {
  const double l = q->lambda;
  const double y = y_of(q, x);
  return (3.0 * t * x - 2.0 + 2.0 * l * l * l * x / y) / (1.0 - x * x);
}

static double second_derivative(const time_equation* q, double x, double t, double dt) {
  const double l = q->lambda;
  const double y = y_of(q, x);
  return (3.0 * t + 5.0 * x * dt + 2.0 * (1.0 - l * l) * l * l * l / (y * y * y)) / (1.0 - x * x);
}

/* ---- Root finding on each branch ----------------------------------------------------------- */

typedef struct branch_problem {
  const time_equation* q;
  double target; /* T to reach */
  double sign;   /* +1: f = T - target (T increasing), -1: f = target - T (T decreasing) */
} branch_problem;

static void branch_eval(const void* context, double x, double* f, double* df) {
  const branch_problem* p = context;
  const double t = time_of_flight(p->q, x);
  *f = p->sign * (t - p->target);
  *df = p->sign * first_derivative(p->q, x, t);
}

/* Minimum of T for M >= 1: root of T'(x) on (-1, 1), where T is convex. */
typedef struct minimum_problem {
  const time_equation* q;
} minimum_problem;

static void minimum_eval(const void* context, double x, double* f, double* df) {
  const minimum_problem* p = context;
  const double t = time_of_flight(p->q, x);
  const double d1 = first_derivative(p->q, x, t);
  *f = d1;
  *df = second_derivative(p->q, x, t, d1);
}

static const double x_edge = 1.0 - 0x1p-40; /* |x| < 1 for multi-revolution ellipses */

static myb_status solve_branch(const time_equation* q, double target, double sign, double lo,
                               double hi, double start, double* x) {
  const branch_problem p = {.q = q, .target = target, .sign = sign};
  const myb_internal_root_function function = {.eval = branch_eval, .context = &p};
  return myb_internal_newton_bracketed(&function, lo, hi, start, MAX_ITERATIONS, x);
}

/* Zero revolutions: T decreases from +inf (x -> -1) to 0 (x -> inf). */
static myb_status solve_zero_revolutions(const time_equation* q, double target, double* x) {
  double hi = 1.0;
  /* bound: MAX_DOUBLINGS; T(x) -> 0 as x grows */
  for (int i = 0; i < MAX_DOUBLINGS && time_of_flight(q, hi) > target; ++i) {
    hi = 2.0 * hi;
  }
  if (!(time_of_flight(q, hi) <= target)) {
    return MYB_E_NOCONV;
  }
  /* Izzo's starter (eq. 30). */
  const double l = q->lambda;
  const double t00 = myb_acos(l) + l * sqrt(1.0 - l * l);
  const double t1 = 2.0 * (1.0 - l * l * l) / 3.0;
  double start = 0.0;
  if (target >= t00) {
    start = myb_exp(myb_log(t00 / target) * (2.0 / 3.0)) - 1.0;
  } else if (target < t1) {
    start = 2.5 * t1 / target * (t1 - target) / (1.0 - l * l * l * l * l) + 1.0;
  } else {
    start = myb_exp(myb_log(t00 / target) * (myb_log(t1 / t00) / myb_log(2.0))) - 1.0;
  }
  return solve_branch(q, target, -1.0, -x_edge, hi, start, x);
}

/* Geometry of the transfer: unit radial and transverse directions and the constants of the
 * velocity reconstruction. */
typedef struct geometry {
  double lambda;
  double gamma;
  double rho;
  double sigma;
  double r1n;
  double r2n;
  double ir1[3];
  double ir2[3];
  double it1[3];
  double it2[3];
} geometry;

static void velocities(const geometry* g, double x, double y, myb_lambert_solution* out) {
  const double l = g->lambda;
  const double vr1 = g->gamma * ((l * y - x) - g->rho * (l * y + x)) / g->r1n;
  const double vr2 = -g->gamma * ((l * y - x) + g->rho * (l * y + x)) / g->r2n;
  const double vt1 = g->gamma * g->sigma * (y + l * x) / g->r1n;
  const double vt2 = g->gamma * g->sigma * (y + l * x) / g->r2n;
  myb_vec3_combine(g->ir1, vr1, g->it1, vt1, out->v1);
  myb_vec3_combine(g->ir2, vr2, g->it2, vt2, out->v2);
}

myb_status myb_lambert(double mu, const double r1[3], const double r2[3], double tof,
                       int retrograde, int max_revolutions, myb_lambert_solution* solutions,
                       size_t capacity, size_t* count) {
  MYB_CHECK(r1 != NULL && r2 != NULL && solutions != NULL && count != NULL);
  MYB_CHECK(max_revolutions >= 0 && max_revolutions <= 100000);
  MYB_CHECK(capacity >= (size_t)(1 + 2 * max_revolutions));
  *count = 0;
  if (!(mu > 0.0) || !finite(mu) || !(tof > 0.0) || !finite(tof) || !myb_vec3_finite(r1) ||
      !myb_vec3_finite(r2)) {
    return MYB_E_DOMAIN;
  }
  const double r1n = myb_vec3_norm(r1);
  const double r2n = myb_vec3_norm(r2);
  const double chord[3] = {r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};
  const double c = myb_vec3_norm(chord);
  if (!(r1n > 0.0) || !(r2n > 0.0) || !(c > 0.0)) {
    return MYB_E_DOMAIN;
  }
  const double s = 0.5 * (r1n + r2n + c);
  geometry g = {.r1n = r1n, .r2n = r2n};
  double ih[3];
  myb_vec3_combine(r1, 1.0 / r1n, r1, 0.0, g.ir1);
  myb_vec3_combine(r2, 1.0 / r2n, r2, 0.0, g.ir2);
  myb_vec3_cross(g.ir1, g.ir2, ih);
  const double ih_norm = myb_vec3_norm(ih);
  if (!(ih_norm > 0x1p-40)) {
    return MYB_E_DOMAIN; /* collinear: the plane of the transfer is undefined */
  }
  myb_vec3_combine(ih, 1.0 / ih_norm, ih, 0.0, ih);
  double lambda = sqrt(1.0 - c / s);
  if (ih[2] < 0.0) {
    lambda = -lambda;
    myb_vec3_cross(g.ir1, ih, g.it1);
    myb_vec3_cross(g.ir2, ih, g.it2);
  } else {
    myb_vec3_cross(ih, g.ir1, g.it1);
    myb_vec3_cross(ih, g.ir2, g.it2);
  }
  if (retrograde) {
    lambda = -lambda;
    myb_vec3_combine(g.it1, -1.0, g.it1, 0.0, g.it1);
    myb_vec3_combine(g.it2, -1.0, g.it2, 0.0, g.it2);
  }
  const double target = sqrt(2.0 * mu / (s * s * s)) * tof;
  g.lambda = lambda;
  g.gamma = sqrt(0.5 * mu * s);
  g.rho = (r1n - r2n) / c;
  g.sigma = sqrt(1.0 - g.rho * g.rho);

  /* Zero revolutions: always one solution. */
  time_equation q = {.lambda = lambda, .revolutions = 0};
  double x = 0.0;
  myb_status status = solve_zero_revolutions(&q, target, &x);
  if (status != MYB_OK) {
    return status;
  }
  myb_lambert_solution* out = &solutions[0];
  *out = (myb_lambert_solution){.revolutions = 0, .branch = 0};
  velocities(&g, x, y_of(&q, x), out);
  size_t n = 1;

  /* M >= 1: two solutions while T >= T_min(M). T(x) on (-1, 1) is convex with one minimum. */
  const int feasible = (int)(target / pi);
  const int last = feasible < max_revolutions ? feasible : max_revolutions;
  for (int m = 1; m <= last; ++m) {
    q.revolutions = m;
    const minimum_problem mp = {.q = &q};
    const myb_internal_root_function mf = {.eval = minimum_eval, .context = &mp};
    double x_min = 0.0;
    status = myb_internal_newton_bracketed(&mf, -x_edge, x_edge, 0.0, MAX_ITERATIONS, &x_min);
    if (status != MYB_OK) {
      return status;
    }
    if (time_of_flight(&q, x_min) > target) {
      break; /* no solution with m (and more) revolutions */
    }
    /* Izzo's starters for the two branches. */
    const double left_base = myb_exp(myb_log((m * pi + pi) / (8.0 * target)) * (2.0 / 3.0));
    const double right_base = myb_exp(myb_log(8.0 * target / (m * pi)) * (2.0 / 3.0));
    const double left_start = (left_base - 1.0) / (left_base + 1.0);
    const double right_start = (right_base - 1.0) / (right_base + 1.0);
    double x_left = 0.0;
    double x_right = 0.0;
    status = solve_branch(&q, target, -1.0, -x_edge, x_min, left_start, &x_left);
    if (status == MYB_OK) {
      status = solve_branch(&q, target, 1.0, x_min, x_edge, right_start, &x_right);
    }
    if (status != MYB_OK) {
      return status;
    }
    out = &solutions[n++];
    *out = (myb_lambert_solution){.revolutions = m, .branch = -1};
    velocities(&g, x_left, y_of(&q, x_left), out);
    out = &solutions[n++];
    *out = (myb_lambert_solution){.revolutions = m, .branch = 1};
    velocities(&g, x_right, y_of(&q, x_right), out);
  }
  *count = n;
  return MYB_OK;
}
