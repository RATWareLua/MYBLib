/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/budget.h — the budget of effects of a scenario: how much each term of a force model moves
 * the object over the interval, so that what is switched off can be shown to be negligible and
 * what is on can be shown to matter (docs/PLAN.en.md §4: an effect is modelled when its
 * contribution exceeds the noise; the budget warns when a disabled effect exceeds the
 * tolerance).
 *
 * For each term, the scenario is propagated with that term switched (on if it was off, off if
 * it was on) and compared with the reference propagation of the model as configured: the
 * effect is the largest position difference over the interval (checked at `samples` equally
 * spaced times) — the full non-linear effect, not a linear estimate. The largest acceleration
 * of the term along the reference trajectory is reported as well. */
#ifndef MYB_BUDGET_H_
#define MYB_BUDGET_H_

#include <stddef.h>

#include "myb/core.h"
#include "myb/model.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct myb_budget_entry {
  char name[MYB_MODEL_NAME];
  int enabled;              /* as configured */
  double max_acceleration;  /* along the reference trajectory, m s^-2 */
  double effect;            /* max |r_switched - r_reference| over the samples, m */
  int flag;                 /* MYB_BUDGET_* */
} myb_budget_entry;

enum {
  MYB_BUDGET_OK = 0,
  MYB_BUDGET_MISSING = 1,   /* disabled, but its effect exceeds the tolerance */
  MYB_BUDGET_NEGLIGIBLE = 2 /* enabled, effect below the tolerance (could be dropped) */
};

typedef struct myb_budget_options {
  double tolerance;  /* m; > 0 */
  double epsilon;    /* IAS15 precision; 0 = 1e-11 */
  size_t samples;    /* comparison times over the interval (>= 1; 0 = 100) */
} myb_budget_options;

/* entries: myb_force_model_count(model) values, in the order of the terms. The model is left as
 * configured. Errors: MYB_E_DOMAIN (t1 == t0, invalid options), the statuses of the propagation. */
MYB_API MYB_NODISCARD myb_status myb_budget(myb_force_model* model, double t0, const double y0[6],
                                            double t1, const myb_budget_options* options,
                                            myb_budget_entry* entries);

/* The same budget to first order, from one propagation instead of one per term: along the
 * trajectory of the model as configured, every term k drives its own variation,
 *   d_k'' = (da/dx) d_k + (da/dv) d_k' + a_k(t, x, v),   d_k(t0) = d_k'(t0) = 0,
 * that is d_k(t) = int Phi(t, s) a_k(s) ds with the transition matrix of the model — what
 * switching the term changes, while the change stays small against the orbit; `effect` is the
 * largest |d_k| over the samples. The error is of the second order of the effect: right for
 * the small terms, which are the ones a budget has to decide, and the order of magnitude only
 * for those that change the orbit itself. Unlike the difference of two propagations it has no
 * floor of round-off, so that terms far below the accuracy of the integration are still
 * measured. The work is one propagation with the Jacobian of the model
 * (myb_force_model_jacobian) and every term once more: it grows as the number of terms, not as
 * its square. Less than myb_budget from about five terms on when they bring their Jacobians
 * (all the ready-made ones do); a term differentiated by differences costs 13 values per value
 * of the model, so that a few such terms are cheaper to switch one by one. Errors: as
 * myb_budget. */
MYB_API MYB_NODISCARD myb_status myb_budget_linear(myb_force_model* model, double t0,
                                                   const double y0[6], double t1,
                                                   const myb_budget_options* options,
                                                   myb_budget_entry* entries);

#ifdef __cplusplus
}
#endif

#endif /* MYB_BUDGET_H_ */
