/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_observe.c — attitude laws and observables over slices for the Go binding
 * (glue_tasks.h). */
#include <stdlib.h>
#include <string.h>

#include "glue_internal.h"

/* ---- Attitude ----------------------------------------------------------------------------------
 */

myb_status mybgo_attitudes(mybgo_model* model, const mybgo_attitude* law, const double* t,
                           const double* y, size_t count, double* m, size_t* failed) {
  myb_attitude_function function = NULL;
  void* context = NULL;
  double fixed[9];
  myb_orbital_attitude orbital = {0};
  myb_sun_attitude sun = {0};
  myb_spin_attitude spin = {0};
  const size_t bodies = model != NULL ? mybgo_model_bodies(model) : 0;
  switch (law->kind) {
    case MYBGO_ATTITUDE_FIXED:
      memcpy(fixed, law->m, sizeof fixed);
      function = myb_attitude_fixed;
      context = fixed;
      break;
    case MYBGO_ATTITUDE_ORBITAL:
      if (law->body >= bodies) {
        return mybgo_end(MYB_E_DOMAIN, 0, failed);
      }
      orbital = (myb_orbital_attitude){.env = mybgo_model_env(model), .body = law->body};
      function = myb_attitude_orbital;
      context = &orbital;
      break;
    case MYBGO_ATTITUDE_SUN_NADIR:
      if (law->body >= bodies || law->source >= bodies) {
        return mybgo_end(MYB_E_DOMAIN, 0, failed);
      }
      sun = (myb_sun_attitude){
          .env = mybgo_model_env(model), .body = law->body, .source = law->source};
      function = myb_attitude_sun_nadir;
      context = &sun;
      break;
    case MYBGO_ATTITUDE_SPIN:
      memcpy(spin.m0, law->m, sizeof spin.m0);
      spin.rate = law->rate;
      spin.phase = law->phase;
      function = myb_attitude_spin;
      context = &spin;
      break;
    default:
      return mybgo_end(MYB_E_DOMAIN, 0, failed);
  }
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = function(context, t[k], y + 6 * k, y + 6 * k + 3, m + 9 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_free_rotation(const double moments[3], const double q0[4], const double l0[3],
                               const double* t, size_t count, double* q, double* l) {
  myb_free_rotation* rotation = malloc(sizeof *rotation);
  if (rotation == NULL) {
    return MYB_E_NOMEM;
  }
  const myb_status status = myb_free_rotation_init(moments, q0, l0, rotation);
  for (size_t k = 0; k < count && status == MYB_OK; ++k) {
    myb_free_rotation_state(rotation, t[k], q + 4 * k, l + 3 * k);
  }
  free(rotation);
  return status;
}

/* ---- Observables --------------------------------------------------------------------------------
 */

/* The light model of the records (the bodies' providers copied to `copies`). */
static myb_status light_of(const mybgo_light_model* model, const double* gm,
                           const mybgo_provider* bodies, size_t body_count, myb_light_body** list,
                           mybgo_provider** copies, myb_light_model* out) {
  *list = NULL;
  *copies = NULL;
  if (body_count > 0) {
    *list = calloc(body_count, sizeof **list);
    *copies = calloc(body_count, sizeof **copies);
    if (*list == NULL || *copies == NULL) {
      return MYB_E_NOMEM;
    }
  }
  for (size_t k = 0; k < body_count; ++k) {
    (*copies)[k] = bodies[k];
    (*list)[k] =
        (myb_light_body){.gm = gm[k], .state = mybgo_provider_state, .context = &(*copies)[k]};
  }
  *out =
      (myb_light_model){.c = model->c, .gamma = model->gamma, .bodies = *list, .count = body_count};
  return MYB_OK;
}

myb_status mybgo_light_times(const mybgo_light_model* model, const double* gm,
                             const mybgo_provider* bodies, size_t body_count,
                             const mybgo_provider* target, const myb_time* t,
                             const double* observers, size_t count, int direction,
                             myb_light_path* paths, size_t* failed) {
  myb_light_body* list = NULL;
  mybgo_provider* copies = NULL;
  myb_light_model light;
  myb_status status = light_of(model, gm, bodies, body_count, &list, &copies, &light);
  mybgo_provider aim = *target;
  size_t k = 0;
  for (; k < count && status == MYB_OK; ++k) {
    status = myb_light_time(&light, mybgo_provider_state, &aim, t[k], observers + 6 * k,
                            (myb_light_direction)direction, &paths[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  free(copies);
  free(list);
  return mybgo_end(status, k, failed);
}

myb_status mybgo_round_trips(const mybgo_light_model* model, const double* gm,
                             const mybgo_provider* bodies, size_t body_count,
                             const mybgo_provider* station, const mybgo_provider* target,
                             const myb_time* t, size_t count, double* total, myb_time* bounce,
                             myb_time* transmit, size_t* failed) {
  myb_light_body* list = NULL;
  mybgo_provider* copies = NULL;
  myb_light_model light;
  myb_status status = light_of(model, gm, bodies, body_count, &list, &copies, &light);
  mybgo_provider from = *station;
  mybgo_provider aim = *target;
  size_t k = 0;
  for (; k < count && status == MYB_OK; ++k) {
    status = myb_round_trip(&light, mybgo_provider_state, &from, mybgo_provider_state, &aim, t[k],
                            &total[k], &bounce[k], &transmit[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  free(copies);
  free(list);
  return mybgo_end(status, k, failed);
}

void mybgo_aberrations(const double* p, const double* v, size_t count, double c, double* out) {
  for (size_t k = 0; k < count; ++k) {
    myb_aberration(p + 3 * k, v + 3 * k, c, out + 3 * k);
  }
}

void mybgo_deflections(double gm, double c, double gamma, const double* p, const double* q,
                       const double* e, const double* em, size_t count, double dlim, double* out) {
  for (size_t k = 0; k < count; ++k) {
    myb_light_deflection(gm, c, gamma, p + 3 * k, q + 3 * k, e + 3 * k, em[k], dlim, out + 3 * k);
  }
}

void mybgo_shapiro_delays(double gm, double c, double gamma, const double* x1, const double* x2,
                          const double* xb, size_t count, double* delay) {
  for (size_t k = 0; k < count; ++k) {
    delay[k] = myb_shapiro_delay(gm, c, gamma, x1 + 3 * k, x2 + 3 * k, xb + 3 * k);
  }
}

void mybgo_radecs(const double* p, size_t count, double* ra, double* dec) {
  for (size_t k = 0; k < count; ++k) {
    myb_direction_radec(p + 3 * k, &ra[k], &dec[k]);
  }
}

myb_status mybgo_plasma_columns(const myb_plasma_profile* profile, const double* p1,
                                const double* p2, size_t count, double* column, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status = myb_plasma_column(profile, p1 + 3 * k, p2 + 3 * k, &column[k]);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}
