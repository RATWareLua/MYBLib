/* SPDX-License-Identifier: BSD-3-Clause */
/* The indirect method for transfers with low thrust (myb/optimal.h): the state and the
 * costates integrated together, the conditions at the arrival solved for the costates at the
 * start. */
#include <math.h> /* sqrt, fabs, fmax, floor */
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/integ.h"
#include "myb/kepler.h"
#include "myb/optimal.h"
#include "myb/transfer.h"
#include "myb/uncert.h"
#include "num/vec3.h"

enum {
  STATE = 14,        /* x, v, m and their costates */
  COSTATES = 8,      /* the multiplier of the cost and the seven costates */
  MAX_UNKNOWNS = 10, /* and the time of the arrival for the least time, the phase on an orbit */
  NEWTON_STEPS = 25,
  FIRST_STEPS = 50, /* for a first guess: one that leads somewhere may start slowly */
  MAX_STAGES = 120,
  STEP_BUDGET = 200000,
  PHASE_SAMPLES = 72 /* of the target orbit, for the first guess of the phase */
};

/* A problem in the units of the solver, with the integrator. */
typedef struct shooter {
  myb_indirect_goal goal;
  myb_indirect_arrival arrival;
  double period; /* ORBIT: of the target orbit */
  myb_force_model* field;
  double length; /* m */
  double time;   /* s */
  double mass;   /* kg */
  double t0;     /* s */
  double thrust;
  double exhaust;
  double start[7];
  double target[6];
  double duration; /* FUEL: of the transfer */
  double epsilon;
  double lambda0; /* of the propagation in hand */
  myb_dop853* integrator;
} shooter;

/* The field at (t, r, v) and, when wanted, its derivatives by r and v (row-major). */
static myb_status field_at(const shooter* s, double t, const double r[3], const double v[3],
                           double g[3], double dgdr[9], double dgdv[9]) {
  if (s->field == NULL) {
    const double d = myb_vec3_norm(r);
    if (!(d > 0.0) || !isfinite(d)) {
      return MYB_E_DOMAIN;
    }
    const double k = 1.0 / (d * d * d);
    for (int i = 0; i < 3; ++i) {
      g[i] = -k * r[i];
    }
    if (dgdr != NULL) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          dgdr[3 * i + j] = k * (3.0 * r[i] * r[j] / (d * d) - (i == j ? 1.0 : 0.0));
        }
      }
      memset(dgdv, 0, 9 * sizeof *dgdv);
    }
    return MYB_OK;
  }
  const double speed = s->length / s->time;
  const double x[3] = {r[0] * s->length, r[1] * s->length, r[2] * s->length};
  const double w[3] = {v[0] * speed, v[1] * speed, v[2] * speed};
  const double at = s->t0 + t * s->time;
  double a[3];
  myb_status st = myb_force_model_acceleration(s->field, at, x, w, a);
  if (st != MYB_OK) {
    return st;
  }
  const double scale = s->time * s->time / s->length;
  for (int i = 0; i < 3; ++i) {
    g[i] = a[i] * scale;
  }
  if (dgdr != NULL) {
    st = myb_force_model_jacobian(s->field, at, x, w, dgdr, dgdv, NULL);
    if (st != MYB_OK) {
      return st;
    }
    for (int i = 0; i < 9; ++i) {
      dgdr[i] *= s->time * s->time;
      dgdv[i] *= s->time;
    }
  }
  return MYB_OK;
}

typedef struct control {
  double direction[3];
  double throttle;
  double switching;
} control;

/* The thrust that makes the Hamiltonian least at the state and costates z. */
static void control_at(const shooter* s, const double z[STATE], control* c) {
  const double* lv = z + 10;
  const double size = myb_vec3_norm(lv);
  for (int i = 0; i < 3; ++i) {
    c->direction[i] = size > 0.0 ? -lv[i] / size : 0.0;
  }
  const double gain = (s->exhaust * size / z[6] + z[13]) / s->lambda0;
  if (s->goal == MYB_INDIRECT_TIME) {
    c->switching = -gain;
    c->throttle = 1.0;
    return;
  }
  c->switching = 1.0 - gain;
  const double eps = s->epsilon;
  if (c->switching > eps) {
    c->throttle = 0.0;
  } else if (c->switching < -eps) {
    c->throttle = 1.0;
  } else {
    c->throttle = (eps - c->switching) / (2.0 * eps);
  }
}

static myb_status rhs(void* context, double t, const double* z, double* dz) {
  const shooter* s = context;
  const double m = z[6];
  if (!(m > 0.0)) {
    return MYB_E_DOMAIN; /* the propellant is gone */
  }
  double g[3];
  double dgdr[9];
  double dgdv[9];
  const myb_status st = field_at(s, t, z, z + 3, g, dgdr, dgdv);
  if (st != MYB_OK) {
    return st;
  }
  control c;
  control_at(s, z, &c);
  const double push = s->thrust * c.throttle / m;
  const double* lv = z + 10;
  for (int i = 0; i < 3; ++i) {
    dz[i] = z[3 + i];
    dz[3 + i] = g[i] + push * c.direction[i];
  }
  dz[6] = -s->thrust * c.throttle / s->exhaust;
  for (int j = 0; j < 3; ++j) {
    double by_r = 0.0;
    double by_v = 0.0;
    for (int i = 0; i < 3; ++i) {
      by_r += dgdr[3 * i + j] * lv[i];
      by_v += dgdv[3 * i + j] * lv[i];
    }
    dz[7 + j] = -by_r;
    dz[10 + j] = -z[7 + j] - by_v;
  }
  dz[13] = -myb_vec3_norm(lv) * push / m;
  return MYB_OK;
}

static myb_status hamiltonian_at(const shooter* s, double t, const double z[STATE], double* h) {
  double g[3];
  const myb_status st = field_at(s, t, z, z + 3, g, NULL, NULL);
  if (st != MYB_OK) {
    return st;
  }
  control c;
  control_at(s, z, &c);
  const double u = c.throttle;
  double value = myb_vec3_dot(z + 7, z + 3) + myb_vec3_dot(z + 10, g) -
                 myb_vec3_norm(z + 10) * s->thrust * u / z[6] - z[13] * s->thrust * u / s->exhaust;
  if (s->goal == MYB_INDIRECT_TIME) {
    value += s->lambda0;
  } else {
    value += s->lambda0 * s->thrust / s->exhaust * (u - s->epsilon * u * (1.0 - u));
  }
  *h = value;
  return MYB_OK;
}

/* The state and the costates after `duration` from the costates u at the start. */
static myb_status fly(shooter* s, const double* u, double duration, double z[STATE]) {
  if (!(u[0] > 0.0) || !(duration > 0.0) || !isfinite(duration)) {
    return MYB_E_DOMAIN;
  }
  s->lambda0 = u[0];
  double z0[STATE];
  memcpy(z0, s->start, 7 * sizeof *z0);
  memcpy(z0 + 7, u + 1, 7 * sizeof *z0);
  myb_status st = myb_dop853_start(s->integrator, rhs, s, 0.0, z0);
  if (st == MYB_OK) {
    st = myb_dop853_integrate(s->integrator, duration);
  }
  if (st == MYB_OK) {
    memcpy(z, myb_dop853_state(s->integrator), STATE * sizeof *z);
  }
  return st;
}

/* The unknowns: the costates, then the time of the arrival (TIME), then the phase (ORBIT). */
static size_t unknowns(const shooter* s) {
  return (s->goal == MYB_INDIRECT_TIME ? 9u : 8u) + (s->arrival == MYB_INDIRECT_ORBIT ? 1u : 0u);
}

/* Where the phase is among the unknowns, and its condition among the conditions. */
static size_t phase_index(const shooter* s) {
  return s->goal == MYB_INDIRECT_TIME ? 9u : 8u;
}

/* The state of the target conic a time `along` after `target`, and its acceleration. */
static myb_status target_at(const shooter* s, double along, double y[6], double fall[3]) {
  myb_pv from;
  myb_pv to;
  memcpy(from.r, s->target, sizeof from.r);
  memcpy(from.v, s->target + 3, sizeof from.v);
  const myb_status st = myb_kepler_propagate(1.0, &from, along, &to);
  if (st != MYB_OK) {
    return st;
  }
  memcpy(y, to.r, sizeof to.r);
  memcpy(y + 3, to.v, sizeof to.v);
  const double d = myb_vec3_norm(y);
  const double k = -1.0 / (d * d * d);
  for (int i = 0; i < 3; ++i) {
    fall[i] = k * y[i];
  }
  return MYB_OK;
}

/* The conditions for the unknowns u: the state at the arrival, the costate of the mass
 * there, the norm of the costates, for the least time the Hamiltonian, and on an orbit the
 * transversality condition. */
static myb_status conditions(shooter* s, const double* u, double* f) {
  const double duration = s->goal == MYB_INDIRECT_TIME ? u[8] : s->duration;
  const int orbit = s->arrival == MYB_INDIRECT_ORBIT;
  double z[STATE];
  myb_status st = fly(s, u, duration, z);
  if (st != MYB_OK) {
    return st;
  }
  double target[6];
  double fall[3] = {0.0, 0.0, 0.0};
  memcpy(target, s->target, sizeof target);
  if (s->goal == MYB_INDIRECT_TIME || orbit) {
    st = target_at(s, orbit ? u[phase_index(s)] : duration, target, fall);
    if (st != MYB_OK) {
      return st;
    }
  }
  for (int i = 0; i < 6; ++i) {
    f[i] = z[i] - target[i];
  }
  f[6] = z[13];
  double norm = 0.0;
  for (int k = 0; k < COSTATES; ++k) {
    norm += u[k] * u[k];
  }
  f[7] = norm - 1.0;
  if (s->goal == MYB_INDIRECT_TIME) {
    double h = 0.0;
    st = hamiltonian_at(s, duration, z, &h);
    if (st != MYB_OK) {
      return st;
    }
    /* a target that moves: H = l_x . v_T + l_v . a_T; an orbit that does not: H = 0 */
    f[8] = orbit ? h : h - myb_vec3_dot(z + 7, target + 3) - myb_vec3_dot(z + 10, fall);
  }
  if (orbit) {
    /* transversality: the costates orthogonal to the tangent of the orbit */
    f[phase_index(s)] = myb_vec3_dot(z + 7, target + 3) + myb_vec3_dot(z + 10, fall);
  }
  return MYB_OK;
}

/* The conditions and their Jacobian by central differences (myb_shoot_function). */
static myb_status shoot_function(void* context, const double* u, double* f, double* jacobian) {
  shooter* s = context;
  const size_t n = unknowns(s);
  myb_status st = conditions(s, u, f);
  if (st != MYB_OK || jacobian == NULL) {
    return st;
  }
  double trial[MAX_UNKNOWNS];
  double up[MAX_UNKNOWNS];
  double down[MAX_UNKNOWNS];
  for (size_t j = 0; j < n; ++j) {
    memcpy(trial, u, n * sizeof *trial);
    const double h = 1e-6 * fmax(1.0, fabs(u[j]));
    trial[j] = u[j] + h;
    st = conditions(s, trial, up);
    if (st != MYB_OK) {
      return st;
    }
    trial[j] = u[j] - h;
    st = conditions(s, trial, down);
    if (st != MYB_OK) {
      return st;
    }
    for (size_t i = 0; i < n; ++i) {
      jacobian[i * n + j] = (up[i] - down[i]) / (2.0 * h);
    }
  }
  return MYB_OK;
}

static myb_status newton(shooter* s, double* u, double tolerance, int steps, double* work,
                         double* residual) {
  const size_t n = unknowns(s);
  const myb_shoot_options options = {.tolerance = tolerance, .max_iterations = steps};
  return myb_shoot(shoot_function, s, n, n, u, &options, work, NULL, residual);
}

static int problem_valid(const myb_indirect_problem* p) {
  double sum = p->mu + p->thrust + p->exhaust + p->t0 + p->mass0 + p->t1 + p->smoothing +
               p->tolerance + p->phase;
  for (int i = 0; i < 6; ++i) {
    sum += p->start[i] + p->target[i];
  }
  return isfinite(sum) && p->mu > 0.0 && p->thrust > 0.0 && p->exhaust > 0.0 &&
         p->mass0 > 0.0 && p->t1 > p->t0 && p->smoothing >= 0.0 && p->smoothing <= 1.0 &&
         p->tolerance >= 0.0 && p->tries >= 0 &&
         (p->goal == MYB_INDIRECT_FUEL || p->goal == MYB_INDIRECT_TIME) &&
         (p->arrival == MYB_INDIRECT_RENDEZVOUS || p->arrival == MYB_INDIRECT_ORBIT) &&
         myb_vec3_norm(p->start) > 0.0 && myb_vec3_norm(p->target) > 0.0;
}

static myb_status shooter_open(shooter* s, const myb_indirect_problem* p,
                               const myb_allocator* allocator) {
  memset(s, 0, sizeof *s);
  if (!problem_valid(p)) {
    return MYB_E_DOMAIN;
  }
  s->goal = p->goal;
  s->field = p->field;
  s->length = myb_vec3_norm(p->start);
  s->time = sqrt(s->length * s->length * s->length / p->mu);
  s->mass = p->mass0;
  s->t0 = p->t0;
  s->thrust = p->thrust * s->time * s->time / (s->length * s->mass);
  s->exhaust = p->exhaust * s->time / s->length;
  for (int i = 0; i < 3; ++i) {
    s->start[i] = p->start[i] / s->length;
    s->start[3 + i] = p->start[3 + i] * s->time / s->length;
    s->target[i] = p->target[i] / s->length;
    s->target[3 + i] = p->target[3 + i] * s->time / s->length;
  }
  s->start[6] = 1.0;
  s->duration = (p->t1 - p->t0) / s->time;
  s->epsilon = 1.0;
  s->lambda0 = 1.0;
  s->arrival = p->arrival;
  if (s->arrival == MYB_INDIRECT_ORBIT) {
    /* the target orbit must close: its period bounds the phase */
    const double energy =
        0.5 * myb_vec3_dot(s->target + 3, s->target + 3) - 1.0 / myb_vec3_norm(s->target);
    if (!(energy < 0.0)) {
      return MYB_E_DOMAIN;
    }
    const double a = -0.5 / energy;
    s->period = 2.0 * 3.141592653589793 * a * sqrt(a);
  }
  const myb_dop853_options options = {.rtol = 1e-12, .atol = 1e-12, .max_steps = STEP_BUDGET};
  return myb_dop853_create(STATE, &options, allocator, &s->integrator);
}

/* A first guess: a point of the unit sphere of the costates with a positive multiplier of
 * the cost. */
static void random_guess(myb_rng* rng, double u[COSTATES]) {
  double norm = 0.0;
  for (int k = 0; k < COSTATES; ++k) {
    u[k] = myb_rng_normal(rng);
    norm += u[k] * u[k];
  }
  norm = sqrt(norm);
  u[0] = fabs(u[0]);
  for (int k = 0; k < COSTATES; ++k) {
    u[k] = norm > 0.0 ? u[k] / norm : (k == 0 ? 1.0 : 0.0);
  }
}

/* The phase (units of the solver) of the point of the target orbit whose direction from the
 * centre is nearest to e: the best of PHASE_SAMPLES points spread over a period. */
static myb_status nearest_phase(const shooter* s, const double e[3], double* phase) {
  double best = -HUGE_VAL;
  *phase = 0.0;
  for (int k = 0; k < PHASE_SAMPLES; ++k) {
    const double along = s->period * (double)k / (double)PHASE_SAMPLES;
    double y[6];
    double fall[3];
    const myb_status st = target_at(s, along, y, fall);
    if (st != MYB_OK) {
      return st;
    }
    const double cosine = myb_vec3_dot(y, e) / myb_vec3_norm(y);
    if (cosine > best) {
      best = cosine;
      *phase = along;
    }
  }
  return MYB_OK;
}

/* The first guess of the phase when there is no guess of the costates: the point of the orbit
 * in the direction halfway between where the start would coast in `duration` and where the
 * orbit carries its own point nearest to the start in that time — between the slower and the
 * faster of the two motions, as a spiral from one to the other is. */
static myb_status guess_phase(const shooter* s, double duration, double* phase) {
  myb_pv from;
  myb_pv to;
  memcpy(from.r, s->start, sizeof from.r);
  memcpy(from.v, s->start + 3, sizeof from.v);
  myb_status st = myb_kepler_propagate(1.0, &from, duration, &to);
  double nearest = 0.0;
  if (st == MYB_OK) {
    st = nearest_phase(s, s->start, &nearest);
  }
  double y[6];
  double fall[3];
  if (st == MYB_OK) {
    st = target_at(s, nearest + duration, y, fall);
  }
  if (st != MYB_OK) {
    return st;
  }
  const double a = myb_vec3_norm(to.r);
  const double b = myb_vec3_norm(y);
  double e[3];
  for (int i = 0; i < 3; ++i) {
    e[i] = to.r[i] / a + y[i] / b;
  }
  if (!(myb_vec3_norm(e) > 1e-3)) {
    memcpy(e, to.r, sizeof e); /* the two opposite: no halfway, the coast decides */
  }
  return nearest_phase(s, e, phase);
}

myb_status myb_indirect_solve(const myb_indirect_problem* problem, const double* guess,
                              const myb_allocator* allocator, myb_indirect_solution* solution) {
  MYB_CHECK(problem != NULL && solution != NULL);
  memset(solution, 0, sizeof *solution);
  shooter s;
  myb_status st = shooter_open(&s, problem, allocator);
  if (st != MYB_OK) {
    return st;
  }
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  st = myb_internal_alloc_array(&a, myb_shoot_work_size(MAX_UNKNOWNS, MAX_UNKNOWNS),
                                sizeof(double), &block);
  if (st != MYB_OK) {
    myb_dop853_destroy(s.integrator);
    return st;
  }
  double* work = block;
  const int fuel = problem->goal == MYB_INDIRECT_FUEL;
  const int orbit = problem->arrival == MYB_INDIRECT_ORBIT;
  const double tolerance = problem->tolerance > 0.0 ? problem->tolerance : 1e-9;
  const double asked = fuel ? (problem->smoothing > 0.0 ? problem->smoothing : 1e-4) : 1.0;
  double u[MAX_UNKNOWNS] = {0};
  double residual = 0.0;
  double epsilon = asked;
  int tries = 0;
  int stages = 0;
  int found = 0;
  double phase = 0.0; /* ORBIT: the first guess of the phase */
  if (orbit) {
    if (guess != NULL) {
      phase = problem->phase / s.time;
    } else {
      st = guess_phase(&s, s.duration, &phase);
    }
  }
  if (st == MYB_OK && guess != NULL) {
    double norm = 0.0;
    for (int k = 0; k < COSTATES; ++k) {
      norm += guess[k] * guess[k];
    }
    norm = sqrt(norm);
    if (!(norm > 0.0) || !isfinite(norm) || !(guess[0] > 0.0)) {
      st = MYB_E_DOMAIN;
    } else {
      for (int k = 0; k < COSTATES; ++k) {
        u[k] = guess[k] / norm;
      }
      u[8] = s.duration;
      if (orbit) {
        u[phase_index(&s)] = phase;
      }
      s.epsilon = asked;
      found = newton(&s, u, tolerance, NEWTON_STEPS, work, &residual) == MYB_OK;
      st = found ? MYB_OK : MYB_E_NOCONV;
    }
  } else if (st == MYB_OK) {
    myb_rng rng;
    myb_rng_seed(&rng, UINT64_C(0x6d79622d696e6469));
    const int limit = problem->tries > 0 ? problem->tries : 100;
    /* bound: the first guesses allowed */
    for (int k = 0; k < limit && !found; ++k) {
      random_guess(&rng, u);
      u[8] = s.duration;
      if (orbit) {
        u[phase_index(&s)] = phase;
      }
      s.epsilon = 1.0;
      found = newton(&s, u, tolerance, FIRST_STEPS, work, &residual) == MYB_OK;
      tries = k + 1;
    }
    st = found ? MYB_OK : MYB_E_NOCONV;
    epsilon = 1.0;
    if (found && fuel) {
      double factor = 0.5;
      double before[MAX_UNKNOWNS]; /* the solution one step back, at the eps `behind` times this */
      double behind = 0.0;
      memcpy(before, u, sizeof before);
      /* bound: MAX_STAGES attempts, fewer when the steps stop shrinking */
      for (int attempt = 0; attempt < MAX_STAGES && epsilon > asked && factor < 0.999;
           ++attempt) {
        const double next = fmax(asked, epsilon * factor);
        double trial[MAX_UNKNOWNS];
        double left = 0.0;
        s.epsilon = next;
        int solved = 0;
        if (behind > 0.0) {
          /* along the secant of the last two solutions in the logarithm of eps */
          const double along = myb_log(next / epsilon) / myb_log(epsilon / behind);
          for (int k = 0; k < MAX_UNKNOWNS; ++k) {
            trial[k] = u[k] + along * (u[k] - before[k]);
          }
          solved = newton(&s, trial, tolerance, NEWTON_STEPS, work, &left) == MYB_OK;
        }
        if (!solved) {
          memcpy(trial, u, sizeof trial);
          solved = newton(&s, trial, tolerance, NEWTON_STEPS, work, &left) == MYB_OK;
        }
        if (solved) {
          memcpy(before, u, sizeof before);
          behind = epsilon;
          memcpy(u, trial, sizeof u);
          epsilon = next;
          residual = left;
          ++stages;
        } else {
          factor = sqrt(factor); /* a shorter step in the logarithm */
        }
      }
      if (epsilon > asked) {
        st = MYB_E_NOCONV;
      }
    }
  }
  if (found) {
    const double duration = fuel ? s.duration : u[8];
    double z[STATE];
    s.epsilon = epsilon;
    const myb_status flown = fly(&s, u, duration, z);
    if (flown == MYB_OK) {
      memcpy(solution->costate, u, sizeof solution->costate);
      solution->t1 = problem->t0 + duration * s.time;
      solution->mass1 = z[6] * s.mass;
      solution->smoothing = fuel ? epsilon : 0.0;
      solution->residual = residual;
      solution->tries = tries;
      solution->stages = stages;
      if (orbit) {
        const double along = u[phase_index(&s)];
        solution->phase = (along - s.period * floor(along / s.period)) * s.time;
      }
    } else {
      st = flown;
    }
  }
  myb_internal_release(&a, block);
  myb_dop853_destroy(s.integrator);
  return st;
}

myb_status myb_indirect_sample(const myb_indirect_problem* problem,
                               const myb_indirect_solution* solution, const double* times,
                               size_t count, const myb_allocator* allocator,
                               myb_indirect_point* points) {
  MYB_CHECK(problem != NULL && solution != NULL &&
            (count == 0 || (times != NULL && points != NULL)));
  shooter s;
  myb_status st = shooter_open(&s, problem, allocator);
  if (st != MYB_OK) {
    return st;
  }
  const int fuel = problem->goal == MYB_INDIRECT_FUEL;
  const double duration = (solution->t1 - problem->t0) / s.time;
  if (!(solution->costate[0] > 0.0) || !(duration > 0.0) || !isfinite(duration) ||
      (fuel && !(solution->smoothing > 0.0 && solution->smoothing <= 1.0))) {
    myb_dop853_destroy(s.integrator);
    return MYB_E_DOMAIN;
  }
  s.epsilon = fuel ? solution->smoothing : 1.0;
  s.lambda0 = solution->costate[0];
  double z[STATE];
  memcpy(z, s.start, 7 * sizeof *z);
  memcpy(z + 7, solution->costate + 1, 7 * sizeof *z);
  st = myb_dop853_start(s.integrator, rhs, &s, 0.0, z);
  double reached = 0.0;
  for (size_t k = 0; k < count && st == MYB_OK; ++k) {
    double tau = (times[k] - problem->t0) / s.time;
    if (!(tau >= reached) || !(tau <= duration * (1.0 + 1e-12))) {
      st = MYB_E_DOMAIN;
      break;
    }
    tau = tau > duration ? duration : tau;
    if (tau > reached) {
      st = myb_dop853_integrate(s.integrator, tau);
      if (st != MYB_OK) {
        break;
      }
      memcpy(z, myb_dop853_state(s.integrator), sizeof z);
      reached = tau;
    }
    myb_indirect_point* p = &points[k];
    for (int i = 0; i < 3; ++i) {
      p->y[i] = z[i] * s.length;
      p->y[3 + i] = z[3 + i] * s.length / s.time;
    }
    p->y[6] = z[6] * s.mass;
    control c;
    control_at(&s, z, &c);
    memcpy(p->direction, c.direction, sizeof p->direction);
    p->throttle = c.throttle;
    p->switching = c.switching;
    memcpy(p->costate, z + 7, sizeof p->costate);
    st = hamiltonian_at(&s, tau, z, &p->hamiltonian);
  }
  myb_dop853_destroy(s.integrator);
  return st;
}
