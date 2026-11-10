/* SPDX-License-Identifier: BSD-3-Clause */
/* Orientation providers: the rotation models of myb/frames.h in the form of
 * myb_orientation_function. */
#include "core/base.h"
#include "myb/frames.h"

static void finish(const myb_rotation* model, const myb_rotation* base, myb_time t,
                   myb_orientation* out) {
  if (base != NULL) {
    myb_rotation_compose(model, base, &out->rotation);
  } else {
    out->rotation = *model;
  }
  out->ut = t;
  out->xp = 0.0;
  out->yp = 0.0;
}

myb_status myb_orientation_cip(void* context, myb_time t, myb_orientation* orientation) {
  const myb_cip_orientation* c = context;
  MYB_CHECK(c != NULL && c->cip != NULL && c->eop != NULL && c->timesys != NULL &&
            orientation != NULL);
  myb_eop_values eop;
  myb_time ut1;
  const myb_status status = myb_eop_at(c->eop, c->timesys, t, c->tt, &eop, &ut1);
  if (status != MYB_OK) {
    return status;
  }
  myb_cip_terrestrial(c->cip, t, ut1, &eop, &orientation->rotation);
  orientation->ut = ut1;
  orientation->xp = eop.xp;
  orientation->yp = eop.yp;
  return MYB_OK;
}

myb_status myb_orientation_pck(void* context, myb_time t, myb_orientation* orientation) {
  const myb_pck_orientation* c = context;
  MYB_CHECK(c != NULL && c->pck != NULL && orientation != NULL);
  myb_rotation model;
  const myb_status status = myb_pck_rotation(c->pck, c->id, t, &model, NULL);
  if (status != MYB_OK) {
    return status;
  }
  finish(&model, c->base, t, orientation);
  return MYB_OK;
}

myb_status myb_orientation_wgccre(void* context, myb_time t, myb_orientation* orientation) {
  const myb_wgccre_orientation* c = context;
  MYB_CHECK(c != NULL && c->model != NULL && orientation != NULL);
  myb_rotation model;
  myb_wgccre_rotation(c->model, t, &model);
  finish(&model, c->base, t, orientation);
  return MYB_OK;
}

myb_status myb_orientation_fixed(void* rotation, myb_time t, myb_orientation* orientation) {
  const myb_rotation* r = rotation;
  MYB_CHECK(r != NULL && orientation != NULL);
  finish(r, NULL, t, orientation);
  return MYB_OK;
}
