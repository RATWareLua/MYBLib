/* SPDX-License-Identifier: BSD-3-Clause */
#include "integ/event.h"

#include <math.h> /* fabs */
#include <string.h>

#include "core/base.h"

/* Every third iteration bisects, and 2 * 1100 halvings shrink any bracket of finite doubles
 * to adjacent doubles (1075 halvings at most on each side of zero). */
enum { ILLINOIS_ITERATIONS = 3 * 2200 };

static int sign_of(double v) {
  return (v > 0.0) - (v < 0.0);
}

static double larger_of(double a, double b) {
  return a > b ? a : b;
}

static double smaller_of(double a, double b) {
  return a < b ? a : b;
}

myb_status myb_internal_locate_event(const myb_internal_step_view* step, myb_event_function g,
                                     void* context, int* found, double* t_event,
                                     double* y_event) {
  MYB_CHECK(step != NULL && g != NULL && found != NULL && t_event != NULL);
  *found = 0;
  double a = step->t0;
  double b = step->t1;
  double ga = g(context, a, step->y0);
  double gb = g(context, b, step->y1);
  if (gb == 0.0) {
    *found = 1;
    *t_event = b;
    if (y_event != NULL) {
      memcpy(y_event, step->y1, step->m * sizeof(double));
    }
    return MYB_OK;
  }
  if (sign_of(ga) * sign_of(gb) >= 0) {
    return MYB_OK;
  }
  /* |g| at the ends, before the Illinois halvings of the retained end's value. */
  double abs_a = fabs(ga);
  double abs_b = fabs(gb);
  double best = 0.0;
  int exact = 0;
  int side = 0;
  /* bound: ILLINOIS_ITERATIONS (see above); the loop ends when the bracket is one ulp wide */
  for (int iteration = 0; iteration < ILLINOIS_ITERATIONS; ++iteration) {
    double c = (iteration % 3 == 2) ? 0.5 * (a + b) : (ga * b - gb * a) / (ga - gb);
    if (!(c > smaller_of(a, b) && c < larger_of(a, b))) {
      c = 0.5 * (a + b);
    }
    if (!(c > smaller_of(a, b) && c < larger_of(a, b))) {
      break; /* a and b are adjacent doubles */
    }
    const myb_status status = step->dense(step->integrator, c, step->scratch);
    if (status != MYB_OK) {
      return status;
    }
    const double gc = g(context, c, step->scratch);
    if (gc == 0.0) {
      best = c;
      exact = 1;
      break;
    }
    if (sign_of(gc) == sign_of(gb)) {
      b = c;
      gb = gc;
      abs_b = fabs(gc);
      if (side == -1) {
        ga *= 0.5;
      }
      side = -1;
    } else {
      a = c;
      ga = gc;
      abs_a = fabs(gc);
      if (side == 1) {
        gb *= 0.5;
      }
      side = 1;
    }
  }
  if (!exact) {
    best = abs_a <= abs_b ? a : b;
  }
  *found = 1;
  *t_event = best;
  if (y_event != NULL) {
    return step->dense(step->integrator, best, y_event);
  }
  return MYB_OK;
}
