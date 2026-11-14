/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_common.c — what every part of the C half of the Go binding uses: the trampolines that
 * call Go functions back, the providers of states and of orientations described by data. */
#include <string.h>

#include "_cgo_export.h" /* the Go functions called back */
#include "glue_internal.h"

/* ---- Go functions ----------------------------------------------------------------------------
 * The exported Go function finds the function by its handle, calls it and returns its result.
 * The arguments are copied, so that Go receives pointers it may write through without casting
 * away const. */

myb_status mybgo_go_term(void* context, double t, const double x[3], const double v[3],
                         double a[3]) {
  const mybgo_go_context* c = context;
  double xc[3] = {x[0], x[1], x[2]};
  double vc[3] = {v[0], v[1], v[2]};
  return (myb_status)mybgoTerm(c->handle, t, xc, vc, a);
}

double mybgo_go_event(void* context, double t, const double* y) {
  const mybgo_go_context* c = context;
  double yc[6];
  memcpy(yc, y, sizeof yc);
  return mybgoEvent(c->handle, t, yc);
}

static myb_status go_state(const mybgo_go_context* c, myb_time t, double state[6]) {
  return (myb_status)mybgoState(c->handle, t.s, t.f, state);
}

myb_status mybgo_go_action(void* context, double t, const myb_system_view* view,
                           double acceleration[3], double torque[3]) {
  const mybgo_go_context* c = context;
  myb_system_view copy = *view;
  return (myb_status)mybgoAction(c->handle, t, &copy, acceleration, torque);
}

myb_status mybgo_go_shoot(void* context, const double* u, double* f, double* jacobian) {
  const mybgo_go_shoot_context* c = context;

  enum { MAX_UNKNOWNS = 64 }; /* the most myb_shoot takes */

  double uc[MAX_UNKNOWNS];
  if (c->n > MAX_UNKNOWNS) {
    return MYB_E_DOMAIN;
  }
  memcpy(uc, u, c->n * sizeof(double));
  return (myb_status)mybgoShoot(c->handle, uc, f, jacobian);
}

myb_status mybgo_go_gaussian(void* context, double t, double mean[6], double covariance[36]) {
  const mybgo_go_context* c = context;
  return (myb_status)mybgoGaussian(c->handle, t, mean, covariance);
}

/* ---- Providers of states ----------------------------------------------------------------------
 */

myb_status mybgo_provider_state(void* context, myb_time t, double state[6]) {
  const mybgo_provider* p = context;
  switch (p->kind) {
    case MYBGO_PROVIDER_SPK: {
      if (p->spk == NULL) {
        return MYB_E_DOMAIN;
      }
      myb_spk_body body = {
          .spk = p->spk, .target = p->target, .center = p->center, .frame = p->frame};
      return myb_spk_body_state(&body, t, state);
    }
    case MYBGO_PROVIDER_KEPLER: {
      myb_pv pv0;
      memcpy(pv0.r, p->state, sizeof pv0.r);
      memcpy(pv0.v, p->state + 3, sizeof pv0.v);
      myb_pv pv;
      const myb_status status =
          myb_kepler_propagate(p->mu, &pv0, myb_dd_to_double(myb_time_diff(t, p->epoch)), &pv);
      if (status == MYB_OK) {
        memcpy(state, pv.r, sizeof pv.r);
        memcpy(state + 3, pv.v, sizeof pv.v);
      }
      return status;
    }
    case MYBGO_PROVIDER_GO: {
      const mybgo_go_context c = {.handle = p->handle};
      return go_state(&c, t, state);
    }
    default:
      return MYB_E_DOMAIN;
  }
}

myb_status mybgo_provider_states(const mybgo_provider* provider, const myb_time* t, size_t count,
                                 double* states, size_t* failed) {
  mybgo_provider p = *provider;
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = mybgo_provider_state(&p, t[k], states + 6 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

/* ---- Providers of orientations -----------------------------------------------------------------
 */

myb_status mybgo_orientation_of(const mybgo_orientation* d, mybgo_orientation_store* s,
                                myb_orientation_function* function, void** context) {
  *function = NULL;
  *context = NULL;
  s->base = d->base;
  const myb_rotation* base = d->has_base ? &s->base : NULL;
  switch (d->kind) {
    case MYBGO_ORIENTATION_NONE:
      return MYB_OK;
    case MYBGO_ORIENTATION_FIXED:
      s->fixed = d->fixed;
      *function = myb_orientation_fixed;
      *context = &s->fixed;
      return MYB_OK;
    case MYBGO_ORIENTATION_CIP: {
      if (d->cip == NULL || d->eop == NULL || d->timesys == NULL || d->tt_scale == NULL) {
        return MYB_E_DOMAIN;
      }
      myb_scale tt = {0};
      const myb_status status = myb_timesys_scale(d->timesys, d->tt_scale, &tt);
      if (status != MYB_OK) {
        return status;
      }
      s->cip = (myb_cip_orientation){.cip = d->cip, .eop = d->eop, .timesys = d->timesys, .tt = tt};
      *function = myb_orientation_cip;
      *context = &s->cip;
      return MYB_OK;
    }
    case MYBGO_ORIENTATION_PCK:
      if (d->pck == NULL) {
        return MYB_E_DOMAIN;
      }
      s->pck = (myb_pck_orientation){.pck = d->pck, .id = d->pck_id, .base = base};
      *function = myb_orientation_pck;
      *context = &s->pck;
      return MYB_OK;
    case MYBGO_ORIENTATION_WGCCRE:
      if (d->wgccre == NULL) {
        return MYB_E_DOMAIN;
      }
      s->wgccre = (myb_wgccre_orientation){.model = d->wgccre, .base = base};
      *function = myb_orientation_wgccre;
      *context = &s->wgccre;
      return MYB_OK;
    case MYBGO_ORIENTATION_CK:
      if (d->ck == NULL || d->sclk == NULL) {
        return MYB_E_DOMAIN;
      }
      s->ck = (myb_ck_orientation){.ck = d->ck,
                                   .sclk = d->sclk,
                                   .instrument = d->instrument,
                                   .base = d->ck_base,
                                   .tolerance = d->tolerance};
      *function = myb_orientation_ck;
      *context = &s->ck;
      return MYB_OK;
    default:
      return MYB_E_DOMAIN;
  }
}

myb_status mybgo_orientation_sample(const mybgo_orientation* orientation, const myb_time* t,
                                    size_t count, myb_orientation* out, size_t* failed) {
  mybgo_orientation_store store;
  myb_orientation_function function = NULL;
  void* context = NULL;
  myb_status status = mybgo_orientation_of(orientation, &store, &function, &context);
  if (status != MYB_OK || function == NULL) {
    return mybgo_end(status == MYB_OK ? MYB_E_DOMAIN : status, 0, failed);
  }
  return myb_orientation_many(function, context, t, count, out, failed);
}
