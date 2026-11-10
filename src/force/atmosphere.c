/* SPDX-License-Identifier: BSD-3-Clause */
/* Atmospheres in the form of myb_density_function. */
#include <math.h>

#include "core/base.h"
#include "myb/cmath.h"
#include "myb/force.h"

myb_status myb_density_exponential(void* atmosphere, myb_time t, const double r[3],
                                   double* density) {
  const myb_exponential_atmosphere* atm = atmosphere;
  MYB_CHECK(atm != NULL && r != NULL && density != NULL);
  (void)t;
  if (atm->count < 1 || atm->count > MYB_ATMOSPHERE_LAYERS) {
    return MYB_E_DOMAIN;
  }
  for (size_t k = 0; k < atm->count; ++k) {
    if (!isfinite(atm->base[k]) || !isfinite(atm->density[k]) || !isfinite(atm->scale[k]) ||
        !(atm->density[k] >= 0.0) || !(atm->scale[k] > 0.0) ||
        (k > 0 && !(atm->base[k] > atm->base[k - 1]))) {
      return MYB_E_DOMAIN;
    }
  }
  double longitude = 0.0;
  double latitude = 0.0;
  double height = 0.0;
  const myb_status status =
      myb_geodetic_from_cartesian(&atm->ellipsoid, r, &longitude, &latitude, &height);
  if (status != MYB_OK) {
    return status;
  }
  size_t layer = 0;
  /* bound: the layers */
  while (layer + 1 < atm->count && height >= atm->base[layer + 1]) {
    ++layer;
  }
  *density = atm->density[layer] * myb_exp(-(height - atm->base[layer]) / atm->scale[layer]);
  return MYB_OK;
}
