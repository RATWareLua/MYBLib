/* SPDX-License-Identifier: BSD-3-Clause */
/* Optimal transfers with low thrust: the indirect method against the conditions it must meet
 * (the Hamiltonian constant, the arrival where asked, the engine full on or off), against
 * Hohmann's transfer where the thrust is high, the least time against the least propellant;
 * the solver of nonlinear programs on problems with known answers; the method of Sims and
 * Flanagan against the indirect one on the same transfer. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/env.h"
#include "myb/kepler.h"
#include "myb/model.h"
#include "myb/optimal.h"

static const double MU = 3.986004418e14;
static const double PI = 3.141592653589793;

/* The state on a circle of radius r at the angle given from the x-axis, in a plane turned
 * about that axis by `tilt`. With the exact functions of the library: the solvers follow
 * their input to the last bit, and the sine of the C library is not the same everywhere. */
static void on_tilted_circle(double r, double angle, double tilt, double y[6]) {
  const double v = sqrt(MU / r);
  y[0] = r * myb_cos(angle);
  y[1] = r * myb_sin(angle) * myb_cos(tilt);
  y[2] = r * myb_sin(angle) * myb_sin(tilt);
  y[3] = -v * myb_sin(angle);
  y[4] = v * myb_cos(angle) * myb_cos(tilt);
  y[5] = v * myb_cos(angle) * myb_sin(tilt);
}

static void on_circle(double r, double angle, double y[6]) {
  on_tilted_circle(r, angle, 0.0, y);
}

enum { SAMPLES = 401 };

/* What a solution looks like along its way. */
typedef struct profile {
  double miss_r;      /* from the target at the arrival (m) */
  double miss_v;      /* m/s */
  double hamiltonian; /* the largest change along the way */
  double between;     /* the part of the samples with the engine neither off nor full on */
  double burning;     /* the mean throttle */
  double mass1;
  double worst_switch; /* the largest switching function */
} profile;

static profile look(const myb_indirect_problem* problem, const myb_indirect_solution* solution,
                    const double target[6]) {
  static double times[SAMPLES];
  static myb_indirect_point points[SAMPLES];
  for (int k = 0; k < SAMPLES; ++k) {
    times[k] = problem->t0 + (solution->t1 - problem->t0) * k / (SAMPLES - 1);
  }
  times[SAMPLES - 1] = solution->t1;
  MT_CHECK(myb_indirect_sample(problem, solution, times, SAMPLES, NULL, points) == MYB_OK);
  profile p = {0};
  const double* end = points[SAMPLES - 1].y;
  p.miss_r = sqrt((end[0] - target[0]) * (end[0] - target[0]) +
                  (end[1] - target[1]) * (end[1] - target[1]) +
                  (end[2] - target[2]) * (end[2] - target[2]));
  p.miss_v = sqrt((end[3] - target[3]) * (end[3] - target[3]) +
                  (end[4] - target[4]) * (end[4] - target[4]) +
                  (end[5] - target[5]) * (end[5] - target[5]));
  p.mass1 = end[6];
  p.worst_switch = -1e300;
  int between = 0;
  for (int k = 0; k < SAMPLES; ++k) {
    p.hamiltonian = fmax(p.hamiltonian, fabs(points[k].hamiltonian - points[0].hamiltonian));
    between += points[k].throttle > 0.01 && points[k].throttle < 0.99;
    p.burning += points[k].throttle / (double)SAMPLES;
    p.worst_switch = fmax(p.worst_switch, points[k].switching);
    const double* d = points[k].direction;
    MT_CHECK(fabs(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - 1.0) < 1e-12);
    MT_CHECK(points[k].throttle >= 0.0 && points[k].throttle <= 1.0);
    if (k > 0) {
      MT_CHECK(points[k].y[6] <= points[k - 1].y[6]);
    }
  }
  p.between = (double)between / (double)SAMPLES;
  return p;
}

static void test_indirect_fuel(void) {
  /* from a circle of 7000 km to a point of a circle of 8000 km in a plane 0.06 rad away, in
   * 4500 s with 0.3 m/s^2: where the engine burns matters, so it is either full on or off */
  myb_indirect_problem problem = {.mu = MU, .thrust = 300.0, .exhaust = 3000.0, .t0 = 100.0,
                                  .mass0 = 1000.0, .t1 = 4600.0, .goal = MYB_INDIRECT_FUEL};
  on_circle(7.0e6, 0.0, problem.start);
  on_tilted_circle(8.0e6, 4.4, 0.06, problem.target);
  /* the least energy first: smooth */
  problem.smoothing = 1.0;
  myb_indirect_solution smooth;
  MT_CHECK(myb_indirect_solve(&problem, NULL, NULL, &smooth) == MYB_OK);
  const profile ps = look(&problem, &smooth, problem.target);
  printf("indirect, 7000 -> 8000 km and 0.06 rad of plane in 4500 s at 0.3 m/s^2\n"
         "  least energy: %d first guesses, %.3f kg left, misses %.1e m and %.1e m/s, Hamiltonian constant to %.1e\n",
         smooth.tries, smooth.mass1, ps.miss_r, ps.miss_v, ps.hamiltonian);
  MT_CHECK(ps.miss_r < 0.05 && ps.miss_v < 5e-5);
  MT_CHECK(ps.hamiltonian < 1e-8);
  MT_CHECK(smooth.stages == 0 && smooth.smoothing == 1.0);
  /* then the least propellant */
  problem.smoothing = 0.0;
  myb_indirect_solution sharp;
  MT_CHECK(myb_indirect_solve(&problem, NULL, NULL, &sharp) == MYB_OK);
  const profile pf = look(&problem, &sharp, problem.target);
  const double dv = problem.exhaust * log(problem.mass0 / sharp.mass1);
  const double v1 = sqrt(MU / 7.0e6);
  const double v2 = sqrt(MU / 8.0e6);
  const double hohmann = v1 * (sqrt(2.0 * 8.0 / 15.0) - 1.0) + v2 * (1.0 - sqrt(2.0 * 7.0 / 15.0));
  printf("  least propellant (eps = %.0e, %d steps from 1): %.3f kg left, %.1f m/s (%.1f "
         "for Hohmann's without the turn of the plane); the engine on %.1f %% of the time, between "
         "off and full in %.1f %% of the samples; misses %.1e m, %.1e m/s; Hamiltonian to "
         "%.1e\n",
         sharp.smoothing, sharp.stages, sharp.mass1, dv, hohmann, 100.0 * pf.burning,
         100.0 * pf.between, pf.miss_r, pf.miss_v, pf.hamiltonian);
  MT_CHECK(sharp.smoothing == 1e-4 && sharp.stages > 3);
  MT_CHECK(pf.miss_r < 0.05 && pf.miss_v < 5e-5);
  MT_CHECK(pf.hamiltonian < 1e-7);
  MT_CHECK(pf.between < 0.02);                /* the engine full on or off */
  MT_CHECK(sharp.mass1 > smooth.mass1);       /* and less is spent than for the least energy */
  MT_CHECK(dv > hohmann && dv < 2.0 * hohmann);
  /* from its own costates the solution is found at once; from those of a neighbour too */
  myb_indirect_solution again;
  MT_CHECK(myb_indirect_solve(&problem, sharp.costate, NULL, &again) == MYB_OK);
  MT_CHECK(fabs(again.mass1 - sharp.mass1) < 1e-3 && again.tries == 0);
  myb_indirect_problem neighbour = problem;
  on_tilted_circle(8.0e6, 4.41, 0.06, neighbour.target);
  MT_CHECK(myb_indirect_solve(&neighbour, sharp.costate, NULL, &again) == MYB_OK);
  printf("  the target 0.01 rad further on, from the costates of the first: %.3f kg left\n",
         again.mass1);
  MT_CHECK(fabs(again.mass1 - sharp.mass1) < 2.0);
  /* the field as a model gives the same */
  myb_point_mass_term centre = {MU};
  myb_force_model* field = NULL;
  MT_CHECK(myb_force_model_create(NULL, &field) == MYB_OK);
  MT_CHECK(myb_force_model_add(field, "centre", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &centre, NULL) == MYB_OK);
  myb_indirect_problem modelled = problem;
  modelled.field = field;
  MT_CHECK(myb_indirect_solve(&modelled, sharp.costate, NULL, &again) == MYB_OK);
  double apart = 0.0;
  for (int k = 0; k < 8; ++k) {
    apart = fmax(apart, fabs(again.costate[k] - sharp.costate[k]));
  }
  printf("  with the field as a model: costates within %.1e, mass within %.1e kg\n", apart,
         fabs(again.mass1 - sharp.mass1));
  MT_CHECK(apart < 1e-6 && fabs(again.mass1 - sharp.mass1) < 1e-4);
  myb_force_model_destroy(field);
}

/* High thrust: Hohmann's transfer with the two burns centred on its ends. */
static void test_indirect_hohmann(void) {
  const double r1 = 7.0e6;
  const double r2 = 8.0e6;
  const double v1 = sqrt(MU / r1);
  const double v2 = sqrt(MU / r2);
  const double at = 0.5 * (r1 + r2);
  const double dv1 = v1 * (sqrt(r2 / at) - 1.0);
  const double dv2 = v2 * (1.0 - sqrt(r1 / at));
  const double coast = PI * sqrt(at * at * at / MU);
  const double thrust = 1000.0; /* 1 m/s^2 on 1000 kg */
  const double b1 = dv1 * 1000.0 / thrust;
  const double b2 = dv2 * 1000.0 / thrust;
  myb_indirect_problem problem = {.mu = MU, .thrust = thrust, .exhaust = 3000.0,
                                  .mass0 = 1000.0, .t1 = coast + 0.5 * (b1 + b2),
                                  .goal = MYB_INDIRECT_FUEL, .smoothing = 1e-3};
  on_circle(r1, -0.5 * b1 * v1 / r1, problem.start);
  on_circle(r2, PI + 0.5 * b2 * v2 / r2, problem.target);
  myb_indirect_solution solution;
  MT_CHECK(myb_indirect_solve(&problem, NULL, NULL, &solution) == MYB_OK);
  const profile p = look(&problem, &solution, problem.target);
  const double dv = problem.exhaust * log(problem.mass0 / solution.mass1);
  printf("indirect at 1 m/s^2 over Hohmann's half ellipse: %.2f m/s against %.2f of the two "
         "impulses (%.2f %% more), the engine on %.1f %% of the time; %d first guesses, %d "
         "steps\n",
         dv, dv1 + dv2, 100.0 * (dv / (dv1 + dv2) - 1.0), 100.0 * p.burning, solution.tries,
         solution.stages);
  MT_CHECK(p.miss_r < 0.05 && p.miss_v < 5e-5);
  /* never below the impulses, and above them by the loss of burns that last: about 1 % */
  MT_CHECK(dv > dv1 + dv2 && dv < 1.02 * (dv1 + dv2));
}

static void test_indirect_time(void) {
  /* the least time to a target that moves on the circle of 8000 km */
  myb_indirect_problem problem = {.mu = MU, .thrust = 300.0, .exhaust = 3000.0, .t0 = 0.0,
                                  .mass0 = 1000.0, .t1 = 2500.0, .goal = MYB_INDIRECT_TIME};
  on_circle(7.0e6, 0.0, problem.start);
  on_circle(8.0e6, 0.5, problem.target);
  myb_indirect_solution fastest;
  MT_CHECK(myb_indirect_solve(&problem, NULL, NULL, &fastest) == MYB_OK);
  myb_pv from;
  myb_pv to;
  memcpy(from.r, problem.target, sizeof from.r);
  memcpy(from.v, problem.target + 3, sizeof from.v);
  MT_CHECK(myb_kepler_propagate(MU, &from, fastest.t1 - problem.t0, &to) == MYB_OK);
  double there[6];
  memcpy(there, to.r, sizeof to.r);
  memcpy(there + 3, to.v, sizeof to.v);
  const profile p = look(&problem, &fastest, there);
  const double spent = problem.mass0 - fastest.mass1;
  printf("indirect, the least time to a target on the circle of 8000 km: %.3f s (%d first "
         "guesses), %.3f kg spent (%.3f at full thrust all the way); misses %.1e m, %.1e m/s; "
         "Hamiltonian to %.1e; the switching function below %.2e\n",
         fastest.t1 - problem.t0, fastest.tries, spent,
         problem.thrust / problem.exhaust * (fastest.t1 - problem.t0), p.miss_r, p.miss_v,
         p.hamiltonian, p.worst_switch);
  MT_CHECK(p.miss_r < 0.05 && p.miss_v < 5e-5);
  MT_CHECK(p.hamiltonian < 1e-8);
  MT_CHECK(p.worst_switch < 0.0); /* full thrust is what the principle asks everywhere */
  MT_CHECK(fabs(spent - problem.thrust / problem.exhaust * (fastest.t1 - problem.t0)) < 1e-6);
  /* the least propellant in that very time needs the engine nearly all the way; in 3 % less
   * the target cannot be reached at all */
  myb_indirect_problem same = problem;
  same.goal = MYB_INDIRECT_FUEL;
  same.t1 = fastest.t1;
  same.smoothing = 1e-3;
  memcpy(same.target, there, sizeof there);
  myb_indirect_solution frugal;
  MT_CHECK(myb_indirect_solve(&same, NULL, NULL, &frugal) == MYB_OK);
  printf("  the least propellant in that time: %.3f kg spent\n", problem.mass0 - frugal.mass1);
  MT_CHECK(fabs((problem.mass0 - frugal.mass1) / spent - 1.0) < 0.01);
  myb_indirect_problem shorter = same;
  shorter.t1 = problem.t0 + 0.97 * (fastest.t1 - problem.t0);
  shorter.tries = 8;
  MT_CHECK(myb_kepler_propagate(MU, &from, shorter.t1 - problem.t0, &to) == MYB_OK);
  memcpy(shorter.target, to.r, sizeof to.r);
  memcpy(shorter.target + 3, to.v, sizeof to.v);
  MT_CHECK(myb_indirect_solve(&shorter, NULL, NULL, &frugal) == MYB_E_NOCONV);
  /* refused */
  myb_indirect_problem bad = problem;
  bad.thrust = 0.0;
  MT_CHECK(myb_indirect_solve(&bad, NULL, NULL, &frugal) == MYB_E_DOMAIN);
  bad = problem;
  bad.t1 = problem.t0;
  MT_CHECK(myb_indirect_solve(&bad, NULL, NULL, &frugal) == MYB_E_DOMAIN);
  const double none[8] = {0};
  MT_CHECK(myb_indirect_solve(&problem, none, NULL, &frugal) == MYB_E_DOMAIN);
  const double late = fastest.t1 + 10.0;
  myb_indirect_point point;
  MT_CHECK(myb_indirect_sample(&problem, &fastest, &late, 1, NULL, &point) == MYB_E_DOMAIN);
}

/* ---- The solver of nonlinear programs -------------------------------------------------------- */

/* Problem 71 of W. Hock and K. Schittkowski (Test Examples for Nonlinear Programming Codes,
 * 1981): min x1 x4 (x1 + x2 + x3) + x3 with x1 x2 x3 x4 >= 25, the sum of the squares 40,
 * 1 <= x <= 5. */
static myb_status hs71(void* context, const double* x, double* f, double* c, double* gradient,
                       double* jacobian) {
  (void)context;
  *f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
  c[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 40.0;
  c[1] = 25.0 - x[0] * x[1] * x[2] * x[3];
  if (gradient != NULL) {
    gradient[0] = x[3] * (2.0 * x[0] + x[1] + x[2]);
    gradient[1] = x[0] * x[3];
    gradient[2] = x[0] * x[3] + 1.0;
    gradient[3] = x[0] * (x[0] + x[1] + x[2]);
    for (int j = 0; j < 4; ++j) {
      jacobian[j] = 2.0 * x[j];
    }
    jacobian[4] = -x[1] * x[2] * x[3];
    jacobian[5] = -x[0] * x[2] * x[3];
    jacobian[6] = -x[0] * x[1] * x[3];
    jacobian[7] = -x[0] * x[1] * x[2];
  }
  return MYB_OK;
}

/* Rosenbrock's valley, to be left at a bound. */
static myb_status valley(void* context, const double* x, double* f, double* c, double* gradient,
                         double* jacobian) {
  (void)context;
  (void)c;
  (void)jacobian;
  const double a = 1.0 - x[0];
  const double b = x[1] - x[0] * x[0];
  *f = a * a + 100.0 * b * b;
  if (gradient != NULL) {
    gradient[0] = -2.0 * a - 400.0 * x[0] * b;
    gradient[1] = 200.0 * b;
  }
  return MYB_OK;
}

/* min x + y on the circle of radius sqrt 2. */
static myb_status on_a_circle(void* context, const double* x, double* f, double* c,
                              double* gradient, double* jacobian) {
  (void)context;
  *f = x[0] + x[1];
  c[0] = x[0] * x[0] + x[1] * x[1] - 2.0;
  if (gradient != NULL) {
    gradient[0] = gradient[1] = 1.0;
    jacobian[0] = 2.0 * x[0];
    jacobian[1] = 2.0 * x[1];
  }
  return MYB_OK;
}

/* A linear program: max x + y with x + 2 y <= 4 in a box. */
static myb_status linear(void* context, const double* x, double* f, double* c, double* gradient,
                         double* jacobian) {
  (void)context;
  *f = -x[0] - x[1];
  c[0] = x[0] + 2.0 * x[1] - 4.0;
  if (gradient != NULL) {
    gradient[0] = gradient[1] = -1.0;
    jacobian[0] = 1.0;
    jacobian[1] = 2.0;
  }
  return MYB_OK;
}

/* Fails away from the start: the solver must come back. */
static myb_status fenced(void* context, const double* x, double* f, double* c, double* gradient,
                         double* jacobian) {
  (void)context;
  (void)c;
  (void)jacobian;
  if (x[0] > 3.0) {
    return MYB_E_DOMAIN;
  }
  *f = (x[0] - 2.5) * (x[0] - 2.5);
  if (gradient != NULL) {
    gradient[0] = 2.0 * (x[0] - 2.5);
  }
  return MYB_OK;
}

static void test_nlp(void) {
  const double lower4[4] = {1.0, 1.0, 1.0, 1.0};
  const double upper4[4] = {5.0, 5.0, 5.0, 5.0};
  const myb_nlp_problem p71 = {.n = 4, .equalities = 1, .inequalities = 1, .lower = lower4,
                               .upper = upper4, .function = hs71};
  double x[4] = {1.0, 5.0, 5.0, 1.0};
  double multipliers[2];
  myb_nlp_result result;
  MT_CHECK(myb_nlp_solve(&p71, NULL, x, multipliers, NULL, &result) == MYB_OK);
  const double best[4] = {1.0, 4.74299963, 3.82114998, 1.37940829};
  double apart = 0.0;
  for (int i = 0; i < 4; ++i) {
    apart = fmax(apart, fabs(x[i] - best[i]));
  }
  printf("nonlinear programs\n  Hock and Schittkowski 71: f = %.7f (17.0140173), x within "
         "%.1e, constraints to %.1e, multipliers %.5f and %.5f; %d updates, %d steps, %d "
         "evaluations\n",
         result.objective, apart, result.violation, multipliers[0], multipliers[1], result.outer,
         result.inner, result.evaluations);
  MT_CHECK(fabs(result.objective - 17.0140173) < 1e-6);
  MT_CHECK(apart < 1e-5);
  MT_CHECK(result.violation <= 1e-8 && result.stationarity <= 1e-6);
  /* the multipliers of the solution: 0.1614686 for the sum of the squares and 0.5522937 for
   * the product (in the signs of the constraints as written here) */
  MT_CHECK(fabs(multipliers[0] - 0.1614686) < 1e-5 && fabs(multipliers[1] - 0.5522937) < 1e-5);
  /* the valley with its floor outside the box: the least value on the wall x = 0.5 */
  const double lo2[2] = {-2.0, -2.0};
  const double hi2[2] = {0.5, 2.0};
  const myb_nlp_problem walled = {.n = 2, .lower = lo2, .upper = hi2, .function = valley};
  double y[2] = {-1.2, 1.0};
  MT_CHECK(myb_nlp_solve(&walled, NULL, y, NULL, NULL, &result) == MYB_OK);
  printf("  Rosenbrock's valley behind a wall: (%.8f, %.8f), f = %.8f; %d steps\n", y[0], y[1],
         result.objective, result.inner);
  MT_CHECK(y[0] == 0.5 && fabs(y[1] - 0.25) < 1e-6 && fabs(result.objective - 0.25) < 1e-9);
  /* without the wall: the floor of the valley */
  const myb_nlp_problem open = {.n = 2, .function = valley};
  y[0] = -1.2;
  y[1] = 1.0;
  MT_CHECK(myb_nlp_solve(&open, NULL, y, NULL, NULL, &result) == MYB_OK);
  printf("  and without it: (%.8f, %.8f); %d steps\n", y[0], y[1], result.inner);
  MT_CHECK(fabs(y[0] - 1.0) < 1e-5 && fabs(y[1] - 1.0) < 1e-5);
  /* an equality alone */
  const myb_nlp_problem circle = {.n = 2, .equalities = 1, .function = on_a_circle};
  y[0] = 0.3;
  y[1] = -1.7;
  MT_CHECK(myb_nlp_solve(&circle, NULL, y, multipliers, NULL, &result) == MYB_OK);
  printf("  the lowest point of x + y on a circle: (%.9f, %.9f), multiplier %.9f (0.5)\n", y[0],
         y[1], multipliers[0]);
  MT_CHECK(fabs(y[0] + 1.0) < 1e-7 && fabs(y[1] + 1.0) < 1e-7);
  MT_CHECK(fabs(multipliers[0] - 0.5) < 1e-6);
  /* a linear program: the corner (3, 0.5) */
  const double lo0[2] = {0.0, 0.0};
  const double hi3[2] = {3.0, 3.0};
  const myb_nlp_problem lp = {.n = 2, .inequalities = 1, .lower = lo0, .upper = hi3,
                              .function = linear};
  y[0] = 1.0;
  y[1] = 1.0;
  MT_CHECK(myb_nlp_solve(&lp, NULL, y, multipliers, NULL, &result) == MYB_OK);
  printf("  a linear program: (%.9f, %.9f), multiplier %.9f (0.5)\n", y[0], y[1],
         multipliers[0]);
  MT_CHECK(y[0] == 3.0 && fabs(y[1] - 0.5) < 1e-7 && fabs(multipliers[0] - 0.5) < 1e-6);
  /* a function that fails beyond x = 3: the steps are shortened */
  const myb_nlp_problem fence = {.n = 1, .function = fenced};
  double z = -40.0;
  MT_CHECK(myb_nlp_solve(&fence, NULL, &z, NULL, NULL, &result) == MYB_OK);
  MT_CHECK(fabs(z - 2.5) < 1e-6);
  z = 4.0;
  MT_CHECK(myb_nlp_solve(&fence, NULL, &z, NULL, NULL, &result) == MYB_E_DOMAIN);
  /* refused */
  const myb_nlp_problem crossed = {.n = 2, .lower = hi2, .upper = lo2, .function = valley};
  MT_CHECK(myb_nlp_solve(&crossed, NULL, y, NULL, NULL, &result) == MYB_E_DOMAIN);
  const myb_nlp_problem none = {.n = 0, .function = valley};
  MT_CHECK(myb_nlp_solve(&none, NULL, y, NULL, NULL, &result) == MYB_E_DOMAIN);
}

/* ---- Sims and Flanagan ----------------------------------------------------------------------- */

static void test_sims_flanagan(void) {
  /* the transfer of test_indirect_fuel */
  myb_sims_flanagan leg = {.mu = MU, .thrust = 300.0, .exhaust = 3000.0, .t0 = 100.0,
                           .mass0 = 1000.0, .t1 = 4600.0, .segments = 30};
  on_circle(7.0e6, 0.0, leg.start);
  on_tilted_circle(8.0e6, 4.4, 0.06, leg.target);
  /* the Jacobian of the mismatch against differences, at throttles of all sizes */
  enum { N = 30 };
  static double throttles[3 * N];
  static double jacobian[6 * 3 * N];
  mt_rng rng = {19};
  for (int i = 0; i < N; ++i) {
    const double size = i % 5 == 0 ? 0.0 : ((double)(mt_next(&rng) >> 11) * 0x1p-53);
    const double a = 6.28 * ((double)(mt_next(&rng) >> 11) * 0x1p-53);
    const double b = 1.4 * (((double)(mt_next(&rng) >> 11) * 0x1p-53) - 0.5);
    throttles[3 * i] = size * myb_cos(b) * myb_cos(a);
    throttles[3 * i + 1] = size * myb_cos(b) * myb_sin(a);
    throttles[3 * i + 2] = size * myb_sin(b);
  }
  double mismatch[6];
  double mass1 = 0.0;
  MT_CHECK(myb_sims_flanagan_mismatch(&leg, throttles, mismatch, &mass1, jacobian) == MYB_OK);
  double worst = 0.0;
  for (int j = 0; j < 3 * N; ++j) {
    if (throttles[3 * (j / 3)] == 0.0 && throttles[3 * (j / 3) + 1] == 0.0) {
      continue; /* a throttle at zero: no derivative of its size */
    }
    const double keep = throttles[j];
    const double h = 1e-6;
    double up[6];
    double down[6];
    throttles[j] = keep + h;
    MT_CHECK(myb_sims_flanagan_mismatch(&leg, throttles, up, NULL, NULL) == MYB_OK);
    throttles[j] = keep - h;
    MT_CHECK(myb_sims_flanagan_mismatch(&leg, throttles, down, NULL, NULL) == MYB_OK);
    throttles[j] = keep;
    for (int r = 0; r < 6; ++r) {
      const double scale = r < 3 ? 1.0e5 : 1.0e2; /* what a throttle moves: m, m/s */
      worst = fmax(worst, fabs(jacobian[r * 3 * N + j] - (up[r] - down[r]) / (2.0 * h)) / scale);
    }
  }
  printf("Sims and Flanagan, 30 segments\n  the Jacobian of the mismatch against differences: "
         "%.1e of its size; the mismatch of random throttles %.0f km, %.1f kg left\n",
         worst, sqrt(mismatch[0] * mismatch[0] + mismatch[1] * mismatch[1] +
                     mismatch[2] * mismatch[2]) / 1e3, mass1);
  MT_CHECK(worst < 1e-6);
  /* the least propellant, from the first guess built in */
  memset(throttles, 0, sizeof throttles);
  myb_sims_flanagan_result result;
  MT_CHECK(myb_sims_flanagan_solve(&leg, NULL, throttles, NULL, &result) == MYB_OK);
  int full = 0;
  int off = 0;
  for (int i = 0; i < N; ++i) {
    const double size = sqrt(throttles[3 * i] * throttles[3 * i] +
                             throttles[3 * i + 1] * throttles[3 * i + 1] +
                             throttles[3 * i + 2] * throttles[3 * i + 2]);
    MT_CHECK(size <= 1.0 + 1e-12);
    full += size > 0.999;
    off += size < 0.001;
  }
  /* the same transfer by the indirect method */
  myb_indirect_problem problem = {.mu = MU, .thrust = 300.0, .exhaust = 3000.0, .t0 = 100.0,
                                  .mass0 = 1000.0, .t1 = 4600.0, .goal = MYB_INDIRECT_FUEL};
  memcpy(problem.start, leg.start, sizeof leg.start);
  memcpy(problem.target, leg.target, sizeof leg.target);
  myb_indirect_solution continuous;
  MT_CHECK(myb_indirect_solve(&problem, NULL, NULL, &continuous) == MYB_OK);
  const double spent = leg.mass0 - result.mass1;
  const double spent_continuous = leg.mass0 - continuous.mass1;
  printf("  the least propellant: %.3f kg spent against %.3f of the indirect method (%+.2f "
         "%%); %d segments full on, %d off; the halves meet to %.1e m and %.1e m/s; %d "
         "updates, %d steps, %d evaluations\n",
         spent, spent_continuous, 100.0 * (spent / spent_continuous - 1.0), full, off,
         sqrt(result.mismatch[0] * result.mismatch[0] + result.mismatch[1] * result.mismatch[1] +
              result.mismatch[2] * result.mismatch[2]),
         sqrt(result.mismatch[3] * result.mismatch[3] + result.mismatch[4] * result.mismatch[4] +
              result.mismatch[5] * result.mismatch[5]),
         result.nlp.outer, result.nlp.inner, result.nlp.evaluations);
  MT_CHECK(fabs(spent / spent_continuous - 1.0) < 0.02);
  MT_CHECK(full + off >= N - 6); /* the engine is full on or off but for a few segments */
  MT_CHECK(fabs(result.mismatch[0]) + fabs(result.mismatch[1]) + fabs(result.mismatch[2]) < 0.1);
  /* what the solver returns meets by the public function too */
  MT_CHECK(myb_sims_flanagan_mismatch(&leg, throttles, mismatch, &mass1, NULL) == MYB_OK);
  MT_CHECK(fabs(mismatch[0] - result.mismatch[0]) < 1e-6 && mass1 == result.mass1);
  /* from its own answer the solver has nothing to do */
  myb_sims_flanagan_result again;
  MT_CHECK(myb_sims_flanagan_solve(&leg, NULL, throttles, NULL, &again) == MYB_OK);
  MT_CHECK(fabs(again.mass1 - result.mass1) < 1e-3);
  /* twice the segments: closer to the continuous optimum */
  leg.segments = 60;
  static double finer[3 * 60];
  memset(finer, 0, sizeof finer);
  MT_CHECK(myb_sims_flanagan_solve(&leg, NULL, finer, NULL, &again) == MYB_OK);
  printf("  with 60 segments: %.3f kg spent (%+.2f %%); %d steps\n", leg.mass0 - again.mass1,
         100.0 * ((leg.mass0 - again.mass1) / spent_continuous - 1.0), again.nlp.inner);
  MT_CHECK(fabs((leg.mass0 - again.mass1) / spent_continuous - 1.0) < 0.01);
  /* refused */
  leg.segments = 1;
  MT_CHECK(myb_sims_flanagan_solve(&leg, NULL, finer, NULL, &again) == MYB_E_DOMAIN);
  leg.segments = 30;
  throttles[0] = 2.0;
  MT_CHECK(myb_sims_flanagan_mismatch(&leg, throttles, mismatch, NULL, NULL) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_indirect_fuel();
  test_indirect_hohmann();
  test_indirect_time();
  test_nlp();
  test_sims_flanagan();
  return mt_finish("test_optimal");
}
