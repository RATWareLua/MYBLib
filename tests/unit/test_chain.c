/* SPDX-License-Identifier: BSD-3-Clause */
/* Chains of legs of Sims and Flanagan joined by flybys: the derivatives of the program against
 * differences, in the field of the Sun and in a field that changes with time; an Earth - Venus
 * - Mars chain against the same chain kept far from Venus and against the direct transfer.
 * The planets move on circles. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mtest.h"
#include "myb/cmath.h"
#include "myb/env.h"
#include "myb/model.h"
#include "myb/optimal.h"
#include "myb/time.h"

static const double SUN = 1.32712440018e20;
static const double AU = 1.495978707e11;
static const double DAY = 86400.0;

/* A planet on a circle about the Sun, in a plane turned about the x-axis by `tilt`. */
typedef struct planet {
  double radius;
  double phase; /* at J2000 */
  double tilt;
} planet;

static void planet_at(const planet* p, double t, double y[6]) {
  const double n = sqrt(SUN / (p->radius * p->radius * p->radius));
  const double v = n * p->radius;
  const double a = p->phase + n * t;
  const double c = myb_cos(p->tilt);
  const double s = myb_sin(p->tilt);
  y[0] = p->radius * myb_cos(a);
  y[1] = p->radius * myb_sin(a) * c;
  y[2] = p->radius * myb_sin(a) * s;
  y[3] = -v * myb_sin(a);
  y[4] = v * myb_cos(a) * c;
  y[5] = v * myb_cos(a) * s;
}

static myb_status planet_state(void* context, myb_time t, double y[6]) {
  planet_at(context, myb_time_seconds(t), y);
  return MYB_OK;
}

static double uniform(mt_rng* rng) {
  return (double)(mt_next(rng) >> 11) * 0x1p-53;
}

static double norm3(const double* a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/* The largest difference of the analytic derivatives of the chain from central differences,
 * each per the size of its row and of its unknown (the solver's units: an astronomical unit,
 * the time sqrt(au^3 / mu); r_p of a powered flyby per the least distance). The flybys of a
 * two-leg chain: powered or not, as the chain has it. */
static double derivatives_apart(const myb_chain* chain, double* x, size_t n, size_t m,
                                double* worst_objective) {
  static double jacobian[64 * 256];
  static double gradient[256];
  static double c[64];
  static double up[64];
  static double down[64];
  double f = 0.0;
  MT_CHECK(myb_chain_evaluate(chain, x, &f, c, gradient, jacobian, NULL) == MYB_OK);
  const size_t legs = chain->legs;
  const double time = sqrt(AU * AU * AU / SUN);
  const double speed = AU / time;
  const int powered = legs == 2 && chain->impulse_exhaust[1] > 0.0;
  const size_t throttles = 7 * legs + 1 + 3 * (chain->segments[0] + chain->segments[1]);
  double worst = 0.0;
  *worst_objective = 0.0;
  for (size_t j = 0; j < n; ++j) {
    const int is_time = j < 1 + legs;
    const int is_excess = !is_time && j < 7 * legs + 1;
    const int is_flyby = powered && j >= throttles;
    double scale = is_time ? time : (is_excess ? speed : 1.0);
    if (is_flyby) {
      scale = j == throttles ? chain->bodies[1].periapsis : speed;
    } else if (!is_time && !is_excess) {
      const double* tau = x + 7 * legs + 1 + 3 * ((j - 7 * legs - 1) / 3);
      if (tau[0] == 0.0 && tau[1] == 0.0 && tau[2] == 0.0) {
        continue; /* a throttle at zero: no derivative of its size */
      }
    }
    const double h = (is_time ? 1e-5 : 1e-6) * scale;
    const double keep = x[j];
    double f_up = 0.0;
    double f_down = 0.0;
    x[j] = keep + h;
    MT_CHECK(myb_chain_evaluate(chain, x, &f_up, up, NULL, NULL, NULL) == MYB_OK);
    x[j] = keep - h;
    MT_CHECK(myb_chain_evaluate(chain, x, &f_down, down, NULL, NULL, NULL) == MYB_OK);
    x[j] = keep;
    for (size_t r = 0; r < m; ++r) {
      double unit = speed; /* the excess speeds, the speeds of a powered flyby */
      if (r < 6 * legs) {
        unit = r % 6 < 3 ? AU : speed;
      } else if (!powered && r >= 7 * legs - 1 && r < 8 * legs - 2) {
        unit = 1.0; /* the turns */
      } else if (powered && r == 6 * legs) {
        unit = 1.0; /* the turn of the powered flyby */
      }
      const double difference = (up[r] - down[r]) / (2.0 * h);
      worst = fmax(worst, fabs(jacobian[r * n + j] - difference) * scale / unit);
    }
    const double difference = (f_up - f_down) / (2.0 * h);
    *worst_objective =
        fmax(*worst_objective, fabs(gradient[j] - difference) * scale / chain->mass0);
  }
  return worst;
}

/* The Earth, Venus and Mars on circles, Mars in its own plane. */
static planet earth = {AU, 0.0, 0.0};
static planet venus = {0.723332 * AU, 250.0 * 3.141592653589793 / 180.0, 0.0};
static planet mars = {1.523679 * AU, 270.0 * 3.141592653589793 / 180.0, 0.0323};

/* Two legs from the Earth by Venus to a rendezvous with Mars, the times within `window` days
 * of 120, 170 and 220 days after J2000; 1 N on 1000 kg, 30 km/s, 2.5 km/s from the launcher. */
static myb_chain evm_chain(size_t segments, double window) {
  myb_chain chain = {.mu = SUN,
                     .thrust = 1.0,
                     .exhaust = 30000.0,
                     .mass0 = 1000.0,
                     .legs = 2,
                     .launch = 2500.0,
                     .arrival = 0.0};
  chain.segments[0] = chain.segments[1] = segments;
  chain.bodies[0] = (myb_chain_body){planet_state, &earth, 3.986004418e14, 6678.0e3};
  chain.bodies[1] = (myb_chain_body){planet_state, &venus, 3.24858592e14, 6351.8e3};
  chain.bodies[2] = (myb_chain_body){planet_state, &mars, 4.282837e13, 3696.0e3};
  chain.departure[0] = (120.0 - window) * DAY;
  chain.departure[1] = (120.0 + window) * DAY;
  chain.flight[0][0] = (170.0 - window) * DAY;
  chain.flight[0][1] = (170.0 + window) * DAY;
  chain.flight[1][0] = (220.0 - window) * DAY;
  chain.flight[1][1] = (220.0 + window) * DAY;
  return chain;
}

/* Random throttles (one in seven at zero) and the excess velocities moved off the conics. */
static void stir(double* x, size_t legs, size_t total, uint64_t seed) {
  mt_rng rng = {seed};
  for (size_t i = 1 + legs; i < 7 * legs + 1; ++i) {
    x[i] += 300.0 * (uniform(&rng) - 0.5);
  }
  for (size_t i = 0; i < total; ++i) {
    const double size = i % 7 == 3 ? 0.0 : uniform(&rng);
    const double a = 6.28 * uniform(&rng);
    const double b = 1.4 * (uniform(&rng) - 0.5);
    double* tau = x + 7 * legs + 1 + 3 * i;
    tau[0] = size * myb_cos(b) * myb_cos(a);
    tau[1] = size * myb_cos(b) * myb_sin(a);
    tau[2] = size * myb_sin(b);
  }
}

static void test_derivatives(void) {
  myb_chain chain = evm_chain(10, 20.0);
  size_t n = 0;
  size_t equalities = 0;
  size_t inequalities = 0;
  MT_CHECK(myb_chain_size(&chain, &n, &equalities, &inequalities) == MYB_OK);
  MT_CHECK(n == 15 + 60 && equalities == 13 && inequalities == 3);
  static double x[256];
  /* the guess: the conics of Lambert, which meet without the engine when nothing cuts their
   * excess speeds at the ends */
  myb_chain open = chain;
  open.launch = open.arrival = 1e5;
  MT_CHECK(myb_chain_guess(&open, x) == MYB_OK);
  double c[16];
  double f = 0.0;
  MT_CHECK(myb_chain_evaluate(&open, x, &f, c, NULL, NULL, NULL) == MYB_OK);
  double meet = 0.0;
  for (int r = 0; r < 12; ++r) {
    meet = fmax(meet, fabs(c[r]) / (r % 6 < 3 ? 1.0 : 1e-3));
  }
  printf(
      "chains, Earth - Venus - Mars on circles\n  the first guess: the halves of the conics "
      "of Lambert meet to %.1e m (mm/s) and spend nothing; the excess speeds %.0f m/s at the "
      "departure, %.0f and %.0f at Venus, %.0f at Mars\n",
      meet, norm3(x + 3), norm3(x + 6), norm3(x + 9), norm3(x + 12));
  MT_CHECK(meet < 1e-2 && f == 0.0);
  MT_CHECK(myb_chain_guess(&chain, x) == MYB_OK);
  MT_CHECK(fabs(norm3(x + 3) - chain.launch) < 1e-9);
  MT_CHECK(x[12] == 0.0 && x[13] == 0.0 && x[14] == 0.0); /* a rendezvous */
  stir(x, 2, 20, 41);
  double worst_objective = 0.0;
  const double worst = derivatives_apart(&chain, x, n, equalities + inequalities, &worst_objective);
  printf(
      "  the derivatives against central differences, per the size of the row and of the "
      "unknown: constraints %.1e, propellant %.1e\n",
      worst, worst_objective);
  MT_CHECK(worst < 1e-6 && worst_objective < 1e-8);
}

/* A field that changes with time: the Sun and a body of a tenth of its mass on a circle of
 * 3 au (the Jacobian of the latter left to the differences of the model). */
typedef struct companion {
  planet orbit;
  double gm;
} companion;

static myb_status companion_term(void* context, double t, const double x[3], const double v[3],
                                 double a[3]) {
  (void)v;
  const companion* b = context;
  double y[6];
  planet_at(&b->orbit, t, y);
  const double d[3] = {y[0] - x[0], y[1] - x[1], y[2] - x[2]};
  const double dd = norm3(d);
  const double r = norm3(y);
  for (int i = 0; i < 3; ++i) {
    a[i] = b->gm * (d[i] / (dd * dd * dd) - y[i] / (r * r * r));
  }
  return MYB_OK;
}

static void test_field(void) {
  myb_chain chain = evm_chain(6, 20.0);
  myb_point_mass_term sun = {SUN};
  companion other = {{3.0 * AU, 1.0, 0.1}, 0.1 * SUN};
  myb_force_model* field = NULL;
  MT_CHECK(myb_force_model_create(NULL, &field) == MYB_OK);
  MT_CHECK(myb_force_model_add(field, "Sun", myb_term_point_mass, myb_term_point_mass_jacobian,
                               &sun, NULL) == MYB_OK);
  size_t index = 0;
  MT_CHECK(myb_force_model_add(field, "companion", companion_term, NULL, &other, &index) == MYB_OK);
  myb_force_model_set_enabled(field, index, 0);
  size_t n = 0;
  size_t equalities = 0;
  size_t inequalities = 0;
  MT_CHECK(myb_chain_size(&chain, &n, &equalities, &inequalities) == MYB_OK);
  static double x[256];
  MT_CHECK(myb_chain_guess(&chain, x) == MYB_OK);
  stir(x, 2, 12, 43);
  /* the Sun as a model: the conics again, to the tolerance of the integration */
  double by_conics[16];
  double by_model[16];
  double f_conics = 0.0;
  double f_model = 0.0;
  MT_CHECK(myb_chain_evaluate(&chain, x, &f_conics, by_conics, NULL, NULL, NULL) == MYB_OK);
  myb_chain modelled = chain;
  modelled.field = field;
  MT_CHECK(myb_chain_evaluate(&modelled, x, &f_model, by_model, NULL, NULL, NULL) == MYB_OK);
  double apart = 0.0;
  for (int r = 0; r < 12; ++r) {
    apart = fmax(apart, fabs(by_model[r] - by_conics[r]) / (r % 6 < 3 ? AU : 3e4));
  }
  MT_CHECK(apart < 1e-12 && f_model == f_conics);
  /* the companion: the field moves with time, and so must the derivatives by the times */
  myb_force_model_set_enabled(field, index, 1);
  MT_CHECK(myb_chain_evaluate(&modelled, x, NULL, by_model, NULL, NULL, NULL) == MYB_OK);
  double moved = 0.0;
  for (int r = 0; r < 3; ++r) {
    moved = fmax(moved, fabs(by_model[r] - by_conics[r]));
  }
  double worst_objective = 0.0;
  const double worst =
      derivatives_apart(&modelled, x, n, equalities + inequalities, &worst_objective);
  printf(
      "  in a force model: the Sun alone meets as the conics to %.1e of the orbit; with a "
      "companion of a tenth of its mass at 3 au the halves meet %.0f km elsewhere, and the "
      "derivatives agree with differences to %.1e (propellant %.1e)\n",
      apart, moved / 1e3, worst, worst_objective);
  MT_CHECK(moved > 1e6);
  MT_CHECK(worst < 1e-5 && worst_objective < 1e-8);
  myb_force_model_destroy(field);
}

static void print_result(const char* what, const myb_chain* chain, const double* x,
                         const myb_chain_result* r) {
  const size_t legs = chain->legs;
  printf("  %s: %.3f kg left (%.1f m/s); departure %.2f d, flights", what, r->mass[legs],
         chain->exhaust * log(chain->mass0 / r->mass[legs]), x[0] / DAY);
  for (size_t k = 0; k < legs; ++k) {
    printf(" %.2f d", x[1 + k] / DAY);
  }
  printf("; excess %.1f m/s at the departure", r->excess[0][1]);
  for (size_t b = 1; b < legs; ++b) {
    printf(", %.3f and %.3f at the flyby, turned by %.4f rad at %.1f km", r->excess[b][0],
           r->excess[b][1], r->turn[b], r->periapsis[b] / 1e3);
  }
  printf(
      "; %d + %d updates, %d + %d steps, %d + %d evaluations; violation %.1e, stationarity "
      "%.1e\n",
      r->first.outer, r->nlp.outer, r->first.inner, r->nlp.inner, r->first.evaluations,
      r->nlp.evaluations, r->nlp.violation, r->nlp.stationarity);
}

/* The mass that the chain by Venus leaves (test_venus), for test_powered. */
static double g_by_venus;

static void test_venus(void) {
  myb_chain chain = evm_chain(10, 20.0);
  static double x[256];
  memset(x, 0, sizeof x);
  myb_chain_result with;
  MT_CHECK(myb_chain_solve(&chain, NULL, x, NULL, &with) == MYB_OK);
  print_result("by Venus", &chain, x, &with);
  g_by_venus = with.mass[2];
  /* the flyby holds: the excess speed kept, the turn within the least distance; the launcher
   * gives no more than it can and Mars is met */
  double c[16];
  double f = 0.0;
  MT_CHECK(myb_chain_evaluate(&chain, x, &f, c, NULL, NULL, NULL) == MYB_OK);
  double miss_r = 0.0;
  double miss_v = 0.0;
  for (int r = 0; r < 12; ++r) {
    if (r % 6 < 3) {
      miss_r = fmax(miss_r, fabs(c[r]));
    } else {
      miss_v = fmax(miss_v, fabs(c[r]));
    }
  }
  const double allowed =
      2.0 *
      myb_asin(1.0 / (1.0 + 6351.8e3 * with.excess[1][0] * with.excess[1][1] / 3.24858592e14));
  printf(
      "    the legs meet to %.1e m and %.1e m/s, the excess speeds at Venus %.1e m/s apart, "
      "the turn %.4f rad of the %.4f allowed, the excess speed at the departure %+.1e m/s "
      "from what the launcher gives\n",
      miss_r, miss_v, fabs(c[12]), with.turn[1], allowed, c[14]);
  MT_CHECK(with.turn[1] <= allowed * (1.0 + 1e-9));
  MT_CHECK(miss_r < 150.0 && miss_v < 3e-5); /* 1e-9 of the units of the solver */
  MT_CHECK(fabs(c[12]) < 1e-4);              /* m/s */
  MT_CHECK(c[13] < 1e-9 && c[14] < 1e-4);    /* the turn, the launcher */
  MT_CHECK(with.periapsis[1] > 6351.8e3 * (1.0 - 1e-6));
  MT_CHECK(fabs(with.excess[1][0] - with.excess[1][1]) < 1e-4 && with.excess[2][0] == 0.0);
  MT_CHECK(fabs(1000.0 - f - with.mass[2]) < 1e-9);
  MT_CHECK(with.mass[0] == 1000.0 && with.mass[1] < 1000.0 && with.mass[2] < with.mass[1]);
  /* from its own answer the solver has little to do */
  myb_chain_result again;
  static double y[256];
  memcpy(y, x, sizeof y);
  MT_CHECK(myb_chain_solve(&chain, NULL, y, NULL, &again) == MYB_OK);
  MT_CHECK(fabs(again.mass[2] - with.mass[2]) < 1e-3);
  /* the same chain kept 10^6 km from Venus, where it turns the excess velocity by a degree at
   * most: from the guess built in and from the chain by Venus, the same answer, converged.
   * (Before the solver let the slope decide where its values could not tell, the
   * stationarity of this solve stopped at 2e-6 .. 1e-4 and it ended with MYB_E_NOCONV.) */
  myb_chain far = chain;
  far.bodies[1].periapsis = 1e9;
  static double z[256];
  memset(z, 0, sizeof z);
  myb_chain_result distant;
  MT_CHECK(myb_chain_solve(&far, NULL, z, NULL, &distant) == MYB_OK);
  print_result("10^6 km from Venus", &far, z, &distant);
  memcpy(z, x, sizeof z);
  myb_chain_result distant2;
  MT_CHECK(myb_chain_solve(&far, NULL, z, NULL, &distant2) == MYB_OK);
  printf("    and from the chain by Venus: %.3f kg left (%.1e kg apart, stationarity %.1e); the "
         "flyby saves %.1f kg\n",
         distant2.mass[2], fabs(distant.mass[2] - distant2.mass[2]), distant2.nlp.stationarity,
         with.mass[2] - distant.mass[2]);
  MT_CHECK(distant.nlp.violation < 1e-9 && distant2.nlp.violation < 1e-9);
  MT_CHECK(distant.nlp.stationarity < 1e-6 && distant2.nlp.stationarity < 1e-6);
  MT_CHECK(fabs(distant.mass[2] - distant2.mass[2]) < 1e-4);
  MT_CHECK(distant.turn[1] < 0.0196);
  MT_CHECK(distant.mass[2] < with.mass[2] - 50.0); /* the flyby saves */
}

/* ---- Powered flybys -------------------------------------------------------------------------- */

/* The impulse and the turn of a powered flyby against the closed form evaluated at 50 digits
 * (mpmath) from the same doubles: Venus deep and far, the Earth with excess speeds 3e-6 m/s
 * apart (where the difference of the two speeds at the periapsis loses all but 8 digits, and
 * the form of the library keeps them), Mars with a large change; the derivatives against
 * central differences. */
static void test_flyby(void) {
  static const double cases[4][6] = {
      /* gm, v_in, v_out, r_p, dv, turn */
      {324858592000000.0, 5842.6999999999998, 6500.0, 7622160.0, 365.16137188688765,
       1.1147943487857491},
      {324858592000000.0, 6034.5, 4100.25, 1000000000.0, 1909.3676046844877,
       0.027799997005043948},
      {398600441800000.0, 3000.0, 3000.0000030000001, 6678000.0, 7.9432565602204766e-07,
       2.1062499407508304},
      {42828370000000.0, 2650.5, 9000.75, 3696000.0, 4711.788462660541, 0.79748213417128187}};
  /* the derivatives of the closed form by mpmath (of |dv|: the second has v_out < v_in) */
  static const double exact[4][6] = {
      {-0.53475145953627321, 0.57567098794194022, 1.655093156684853e-5,
       -0.00010163991593474199, -8.8951165839458405e-5, -7.6883331941890267e-8},
      {0.99119666401146573, -0.98121965655748747, 2.4381383746432116e-8,
       -2.904712569777457e-6, -9.072934466370697e-6, -2.7364893774029257e-11},
      {-0.26477520918558382, 0.26477520943179673, 5.5303888282945215e-14,
       -0.00015338849435800137, -0.00015338849430705213, -6.8907679429119947e-8},
      {-0.4823024662201102, 0.881795160992514, 0.00026335046616970065,
       -0.0002265723519349056, -2.4518939441832499e-5, -1.1109562540375286e-7}};
  double worst_dv = 0.0;
  double worst_turn = 0.0;
  double worst_partial = 0.0;
  double worst_exact = 0.0;
  for (int k = 0; k < 4; ++k) {
    const double* q = cases[k];
    double dv = 0.0;
    double turn = 0.0;
    double partials[6];
    MT_CHECK(myb_chain_flyby(q[0], q[1], q[2], q[3], &dv, &turn, partials) == MYB_OK);
    worst_dv = fmax(worst_dv, fabs(dv / q[4] - 1.0));
    worst_turn = fmax(worst_turn, fabs(turn / q[5] - 1.0));
    for (int i = 0; i < 6; ++i) {
      worst_exact = fmax(worst_exact, fabs(partials[i] / exact[k][i] - 1.0));
    }
    /* d / d v_in, v_out, r_p against central differences of 1e-6 of each; for the impulse of
     * the speeds as near as 3e-6 m/s, 1e-6 of their difference (not to step across the
     * corner of |.| at equal speeds) */
    for (int i = 0; i < 3; ++i) {
      for (int which = 0; which < 2; ++which) {
        double args[3] = {q[1], q[2], q[3]};
        const double span = which == 0 && i < 2 ? fmin(args[i], fabs(q[2] - q[1])) : args[i];
        const double h = 1e-6 * span;
        const double keep = args[i];
        double value_up[2] = {0.0, 0.0};
        double value_down[2] = {0.0, 0.0};
        args[i] = keep + h;
        const double up = args[i];
        MT_CHECK(myb_chain_flyby(q[0], args[0], args[1], args[2], &value_up[0], &value_up[1],
                                 NULL) == MYB_OK);
        args[i] = keep - h;
        const double down = args[i];
        MT_CHECK(myb_chain_flyby(q[0], args[0], args[1], args[2], &value_down[0],
                                 &value_down[1], NULL) == MYB_OK);
        /* by the steps the doubles take (3e-12 m/s is a few of their units at 3 km/s) */
        const double difference = (value_up[which] - value_down[which]) / (up - down);
        /* by the speeds the impulse per m/s (of the order of one), by r_p per its relative
         * change and the size of the impulse; the turn per the relative change */
        const double argument = i == 2 ? q[3] : q[1 + i];
        const double unit =
            which == 0 ? (i == 2 ? argument / fmax(1.0, q[4]) : 1.0) : argument;
        worst_partial =
            fmax(worst_partial, fabs(partials[3 * which + i] - difference) * unit);
      }
    }
  }
  /* the unpowered flyby is the case of equal speeds: no impulse, the turn 2 asin(1 / e) */
  double dv = 1.0;
  double turn = 0.0;
  MT_CHECK(myb_chain_flyby(3.24858592e14, 5842.7, 5842.7, 7622.16e3, &dv, &turn, NULL) == MYB_OK);
  const double e = 1.0 + 7622.16e3 * 5842.7 * 5842.7 / 3.24858592e14;
  printf(
      "powered flybys\n  against the closed form at 50 digits: the impulse to %.1e, the turn to "
      "%.1e, their derivatives to %.1e (relative); the derivatives against differences %.1e; "
      "equal speeds: dv = %g, the turn of the unpowered flyby to %.1e\n",
      worst_dv, worst_turn, worst_exact, worst_partial, dv, fabs(turn - 2.0 * myb_asin(1.0 / e)));
  MT_CHECK(worst_dv < 4e-15 && worst_turn < 4e-15 && worst_exact < 1e-14);
  MT_CHECK(worst_partial < 1e-8);
  MT_CHECK(dv == 0.0 && fabs(turn - 2.0 * myb_asin(1.0 / e)) < 1e-15);
  /* refused */
  MT_CHECK(myb_chain_flyby(0.0, 1.0, 1.0, 1.0, &dv, &turn, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_chain_flyby(1.0, 0.0, 1.0, 1.0, &dv, &turn, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_chain_flyby(1.0, 1.0, 1.0, (double)NAN, &dv, &turn, NULL) == MYB_E_DOMAIN);
}

/* The chain of test_derivatives with the flyby of Venus powered: its sizes, the guess (the
 * turn of the conics met by r_p, the impulse that of the two hyperbolas), the derivatives of the
 * program by r_p and dv too, through the masses of the second leg. */
static void test_powered_derivatives(void) {
  myb_chain chain = evm_chain(10, 20.0);
  chain.impulse_exhaust[1] = 3000.0;
  size_t n = 0;
  size_t equalities = 0;
  size_t inequalities = 0;
  MT_CHECK(myb_chain_size(&chain, &n, &equalities, &inequalities) == MYB_OK);
  MT_CHECK(n == 15 + 60 + 2 && equalities == 13 && inequalities == 4);
  static double x[256];
  MT_CHECK(myb_chain_guess(&chain, x) == MYB_OK);
  double c[17];
  double f = 0.0;
  MT_CHECK(myb_chain_evaluate(&chain, x, &f, c, NULL, NULL, NULL) == MYB_OK);
  double dv = 0.0;
  double turn = 0.0;
  MT_CHECK(myb_chain_flyby(3.24858592e14, norm3(x + 6), norm3(x + 9), x[75], &dv, &turn, NULL) ==
           MYB_OK);
  printf(
      "  the chain by a powered flyby of Venus: the conics of the guess (%.0f and %.0f m/s) "
      "turn %.3f rad more than two hyperbolas can at the least distance, r_p = %.1f km, where "
      "dv = %.3f m/s (the closed form %.3f); the two sides of the impulse %.1e, %.1e m/s; the "
      "propellant %.3f kg (the rocket equation %.3f)\n",
      norm3(x + 6), norm3(x + 9), c[12], x[75] / 1e3, x[76], dv, c[13], c[16], f,
      1000.0 * (1.0 - exp(-dv / 3000.0)));
  MT_CHECK(x[75] == 6351.8e3 && c[12] > 0.0 && x[76] == dv);
  MT_CHECK(fabs(c[13]) < 1e-9 && fabs(c[16] + 2.0 * dv) < 1e-9);
  MT_CHECK(fabs(f - 1000.0 * (1.0 - exp(-dv / 3000.0))) < 1e-9);
  stir(x, 2, 20, 41);
  x[75] *= 1.3;
  x[76] += 50.0;
  double worst_objective = 0.0;
  const double worst = derivatives_apart(&chain, x, n, equalities + inequalities, &worst_objective);
  printf(
      "    the derivatives against central differences, r_p and dv among them: constraints "
      "%.1e, propellant %.1e\n",
      worst, worst_objective);
  MT_CHECK(worst < 1e-6 && worst_objective < 1e-8);
  /* refused: an impulse at the departure, a negative one, a negative exhaust velocity */
  myb_chain bad = chain;
  bad.impulse_exhaust[0] = 3000.0;
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  bad = chain;
  bad.impulse_exhaust[1] = -1.0;
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  x[76] = -1.0;
  MT_CHECK(myb_chain_evaluate(&chain, x, &f, NULL, NULL, NULL, NULL) == MYB_E_DOMAIN);
}

/* The flyby of a powered chain: the turn of its excess velocities that of the two hyperbolas
 * of r_p (*turn_miss, rad), the impulse that of the closed form (*dv_miss, m/s), r_p no less
 * than allowed. The solver holds them to 1e-9 in its units: 1e-9 rad, 3e-5 m/s. */
static void flyby_miss(const myb_chain* chain, const double* x, const myb_chain_result* r,
                       double* turn_miss, double* dv_miss) {
  double dv = 0.0;
  double turn = 0.0;
  MT_CHECK(myb_chain_flyby(chain->bodies[1].gm, r->excess[1][0], r->excess[1][1],
                           r->periapsis[1], &dv, &turn, NULL) == MYB_OK);
  MT_CHECK(r->periapsis[1] >= chain->bodies[1].periapsis);
  MT_CHECK(r->periapsis[1] == x[75] && r->impulse[1] == x[76]);
  *turn_miss = fabs(turn - r->turn[1]);
  *dv_miss = fabs(dv - r->impulse[1]);
}

/* Where an impulse at the periapsis pays and where it does not. The chain by Venus with an
 * impulse of 3000 m/s (a chemical engine, ten times the propellant of the electric one for a
 * change of speed): the electric engine does better everywhere, and the solver drives dv to 0
 * — the chain of test_venus. The same chain with 0.25 N in place of 1 N: the electric engine
 * cannot meet Mars within the dates (the flights at their longest, 190 and 240 days, and the
 * constraints still far off); an impulse at Venus can, and from the guess built in and from
 * where the unpowered chain stopped the solver finds the same chain. */
static void test_powered(void) {
  myb_chain chain = evm_chain(10, 20.0);
  chain.impulse_exhaust[1] = 3000.0;
  static double x[256];
  memset(x, 0, sizeof x);
  myb_chain_result r;
  MT_CHECK(myb_chain_solve(&chain, NULL, x, NULL, &r) == MYB_OK);
  print_result("powered by Venus, 3000 m/s", &chain, x, &r);
  double turn_miss = 0.0;
  double dv_miss = 0.0;
  flyby_miss(&chain, x, &r, &turn_miss, &dv_miss);
  printf(
      "    dv = %.1e m/s, %.1e kg from the unpowered chain; the turn of the hyperbolas to "
      "%.1e rad, the impulse of the closed form to %.1e m/s\n",
      r.impulse[1], fabs(r.mass[2] - g_by_venus), turn_miss, dv_miss);
  MT_CHECK(r.impulse[1] < 1e-3 && fabs(r.mass[2] - g_by_venus) < 1e-3);
  MT_CHECK(turn_miss < 2e-9 && dv_miss < 6e-5);
  /* a weak electric engine */
  myb_chain weak = evm_chain(10, 20.0);
  weak.thrust = 0.25;
  static double y[256];
  memset(y, 0, sizeof y);
  myb_chain_result unpowered;
  const myb_nlp_options few = {.feasibility = 1e-9, .optimality = 1e-6, .max_outer = 5};
  MT_CHECK(myb_chain_solve(&weak, &few, y, NULL, &unpowered) == MYB_E_NOCONV);
  printf("  unpowered with 0.25 N: flights %.2f and %.2f d, the constraints %.1e off\n",
         y[1] / DAY, y[2] / DAY, unpowered.nlp.violation);
  MT_CHECK(unpowered.nlp.violation > 1e-3);
  MT_CHECK(fabs(y[1] - 190.0 * DAY) < 1.0 && fabs(y[2] - 240.0 * DAY) < 1.0);
  weak.impulse_exhaust[1] = 3000.0;
  static double z[256];
  memset(z, 0, sizeof z);
  myb_chain_result powered;
  MT_CHECK(myb_chain_solve(&weak, NULL, z, NULL, &powered) == MYB_OK);
  print_result("powered with 0.25 N", &weak, z, &powered);
  flyby_miss(&weak, z, &powered, &turn_miss, &dv_miss);
  /* from where the unpowered chain stopped, the impulse at 0 */
  y[75] = weak.bodies[1].periapsis;
  y[76] = 0.0;
  myb_chain_result again;
  MT_CHECK(myb_chain_solve(&weak, NULL, y, NULL, &again) == MYB_OK);
  double c[17];
  double f = 0.0;
  MT_CHECK(myb_chain_evaluate(&weak, z, &f, c, NULL, NULL, NULL) == MYB_OK);
  printf(
      "    dv = %.3f m/s at r_p = %.1f km; the turn of the hyperbolas to %.1e rad, the impulse "
      "of the closed form to %.1e m/s; from the unpowered chain %.6f kg (%.1e kg apart); the "
      "propellant %.1e kg from what the program counts\n",
      powered.impulse[1], powered.periapsis[1] / 1e3, turn_miss, dv_miss, again.mass[2],
      fabs(again.mass[2] - powered.mass[2]), fabs(1000.0 - f - powered.mass[2]));
  MT_CHECK(powered.impulse[1] > 100.0 && turn_miss < 2e-9 && dv_miss < 6e-5);
  MT_CHECK(fabs(again.mass[2] - powered.mass[2]) < 1e-3);
  MT_CHECK(fabs(1000.0 - f - powered.mass[2]) < 1e-9);
  /* dv is the impulse: one side of it holds with equality, the other is 2 dv inside */
  MT_CHECK(fmin(fabs(c[13]), fabs(c[16])) < 6e-5);
  MT_CHECK(fmax(fabs(c[13]), fabs(c[16])) > powered.impulse[1]);
}

static void test_refused(void) {
  myb_chain chain = evm_chain(10, 20.0);
  size_t n = 0;
  myb_chain bad = chain;
  bad.legs = 0;
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  bad = chain;
  bad.bodies[1].gm = 0.0; /* a flyby without a field */
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  bad = chain;
  bad.flight[1][0] = bad.flight[1][1] + 1.0;
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  bad = chain;
  bad.segments[0] = bad.segments[1] = 40; /* 15 + 4 x 80 unknowns for the solver */
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  bad = chain;
  bad.bodies[2].state = NULL;
  MT_CHECK(myb_chain_size(&bad, &n, NULL, NULL) == MYB_E_DOMAIN);
  static double x[256];
  MT_CHECK(myb_chain_guess(&chain, x) == MYB_OK);
  x[15] = 1.5; /* a throttle longer than 1 */
  double f = 0.0;
  MT_CHECK(myb_chain_evaluate(&chain, x, &f, NULL, NULL, NULL, NULL) == MYB_E_DOMAIN);
  MT_CHECK(myb_chain_guess(&chain, x) == MYB_OK);
  x[1] = -DAY; /* a flight backward */
  MT_CHECK(myb_chain_evaluate(&chain, x, &f, NULL, NULL, NULL, NULL) == MYB_E_DOMAIN);
  const myb_nlp_options heavy = {.weight = 1e13};
  myb_chain_result result;
  MT_CHECK(myb_chain_guess(&chain, x) == MYB_OK);
  MT_CHECK(myb_chain_solve(&chain, &heavy, x, NULL, &result) == MYB_E_DOMAIN);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_derivatives();
  test_field();
  test_flyby();
  test_powered_derivatives();
  test_venus();
  test_powered();
  test_refused();
  return mt_finish("test_chain");
}
