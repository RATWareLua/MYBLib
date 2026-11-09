/* SPDX-License-Identifier: BSD-3-Clause */
/* myb/budget.h — the budget of effects of a scenario: how much each term of a force model moves
 * the object over the interval, so that what is switched off can be shown to be negligible and
 * what is on can be shown to matter (docs/PLAN.md §4: an effect is modelled when its
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

#ifdef __cplusplus
}
#endif

#endif /* MYB_BUDGET_H_ */
