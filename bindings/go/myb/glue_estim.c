/* SPDX-License-Identifier: BSD-3-Clause */
/* glue_estim.c — estimation and uncertainty for the Go binding (glue_tasks.h): the filter and
 * the smoother over a batch of observations, measurements over slices, the covariance carried
 * along an orbit, conjunctions of two orbits. */
#include <stdlib.h>
#include <string.h>

#include "glue_internal.h"

myb_status mybgo_filter_run(mybgo_model* model, const mybgo_filter_options* options, double t0,
                            const double* y0, const double* covariance, double t_end,
                            const double* bias, const mybgo_station* stations, size_t station_count,
                            const mybgo_light* light, const mybgo_observation* observations,
                            size_t count, double* innovations, size_t* steps, double* step_t,
                            double* step_unknowns, double* step_covariance, double* final_t,
                            double* final_y, double* final_covariance, double* final_bias,
                            size_t* failed) {
  *steps = 0;
  myb_force_model* force = mybgo_model_force(model);
  const size_t unknowns = 6 + myb_force_model_parameter_count(force) + options->bias_count;
  myb_status status = MYB_E_NOMEM;
  size_t k = 0;
  myb_filter* filter = NULL;
  myb_station* st = calloc(station_count > 0 ? station_count : 1, sizeof *st);
  myb_light_measure* lm = calloc(station_count > 0 ? station_count : 1, sizeof *lm);
  if (st == NULL || lm == NULL) {
    goto cleanup;
  }
  status = mybgo_stations_of(model, stations, station_count, light, st, lm);
  if (status != MYB_OK) {
    goto cleanup;
  }
  const myb_filter_options fo = {.epsilon = options->epsilon,
                                 .process_noise = options->process_noise,
                                 .bias_count = options->bias_count,
                                 .history = options->smooth ? count + 1 : 0};
  status = myb_filter_create(force, &fo, NULL, &filter);
  if (status == MYB_OK) {
    status = myb_filter_start(filter, t0, y0, covariance, t_end);
  }
  if (status == MYB_OK && options->bias_count > 0) {
    status = myb_filter_set_bias(filter, bias);
  }
  size_t value = 0;
  for (; k < count && status == MYB_OK; ++k) {
    myb_observation o = {0};
    status = mybgo_measure_of(&observations[k], st, light != NULL ? lm : NULL, station_count, &o);
    if (status == MYB_OK) {
      status = myb_filter_update(filter, &o, innovations + value);
    }
    if (status != MYB_OK) {
      break;
    }
    value += o.count;
  }
  if (status == MYB_OK && options->smooth) {
    status = myb_filter_smooth(filter);
  }
  if (filter != NULL) {
    myb_filter_get(filter, final_t, final_y, final_covariance);
    if (options->bias_count > 0) {
      myb_filter_bias(filter, final_bias);
    }
    const size_t kept = myb_filter_steps(filter);
    for (size_t j = 0; j < kept && j <= count; ++j) {
      myb_filter_step(filter, j, &step_t[j], step_unknowns + unknowns * j,
                      step_covariance + unknowns * unknowns * j);
    }
    *steps = kept <= count + 1 ? kept : count + 1;
  }
cleanup:
  myb_filter_destroy(filter);
  free(lm);
  free(st);
  return mybgo_end(status, k, failed);
}

myb_status mybgo_measure(mybgo_model* model, int kind, const mybgo_station* station,
                         const mybgo_light* light, const double* t, const double* y, size_t count,
                         double* values, double* partials, size_t* failed) {
  myb_station st = {0};
  myb_light_measure lm = {0};
  const size_t stations = station != NULL ? 1 : 0;
  myb_status status = mybgo_stations_of(model, station, stations, light, &st, &lm);
  if (status != MYB_OK) {
    return mybgo_end(status, 0, failed);
  }
  const mybgo_observation record = {.kind = kind, .station = 0};
  myb_observation o = {0};
  status = mybgo_measure_of(&record, &st, light != NULL ? &lm : NULL, stations, &o);
  if (status != MYB_OK) {
    return mybgo_end(status, 0, failed);
  }
  size_t k = 0;
  for (; k < count; ++k) {
    status = o.measure(o.context, t[k], y + 6 * k, values + o.count * k,
                       partials != NULL ? partials + 6 * o.count * k : NULL);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

/* ---- Orbits with their variational equations ------------------------------------------------- */

static myb_status start_orbit(mybgo_model* model, const mybgo_orbit_start* start,
                              myb_orbit** orbit) {
  myb_orbit_options options = start->options;
  options.variational = 1;
  myb_status status = myb_orbit_create(mybgo_model_force(model), &options, NULL, orbit);
  if (status == MYB_OK) {
    status = myb_orbit_start(*orbit, start->t0, start->y0, start->t_end);
  }
  return status;
}

myb_status mybgo_covariance_propagate(mybgo_model* model, const mybgo_orbit_start* start,
                                      const double* times, size_t count, double* states,
                                      double* covariances, size_t* failed) {
  myb_orbit* orbit = NULL;
  myb_status status = start_orbit(model, start, &orbit);
  size_t k = 0;
  for (; k < count && status == MYB_OK; ++k) {
    double phi[36];
    status = myb_orbit_partials(orbit, times[k], states + 6 * k, phi, NULL);
    if (status != MYB_OK) {
      break;
    }
    myb_covariance_transform(6, phi, start->covariance, covariances + 36 * k);
  }
  myb_orbit_destroy(orbit);
  return mybgo_end(status, k, failed);
}

myb_status mybgo_covariance_rtn(const double* states, const double* covariances, size_t count,
                                double* out, size_t* failed) {
  myb_status status = MYB_OK;
  size_t k = 0;
  for (; k < count; ++k) {
    status =
        myb_covariance_rtn(states + 6 * k, states + 6 * k + 3, covariances + 36 * k, out + 36 * k);
    if (status != MYB_OK) {
      break;
    }
  }
  return mybgo_end(status, k, failed);
}

myb_status mybgo_conjunctions(mybgo_model* first, const mybgo_orbit_start* first_start,
                              mybgo_model* second, const mybgo_orbit_start* second_start,
                              double radius, double t1, double step, myb_conjunction** found,
                              size_t* count, double* probability, double* error) {
  *count = 0;
  if (found != NULL) {
    *found = NULL;
  }
  myb_orbit* a = NULL;
  myb_orbit* b = NULL;
  myb_status status = start_orbit(first, first_start, &a);
  if (status == MYB_OK) {
    status = start_orbit(second, second_start, &b);
  }
  if (status == MYB_OK && probability != NULL) {
    status = myb_conjunction_long(a, first_start->covariance, b, second_start->covariance, radius,
                                  t1, step, probability, error);
  } else if (status == MYB_OK) {
    /* the count may exceed the room: then once more from the start with room for all */
    size_t capacity = 16;
    myb_conjunction* list = NULL;
    for (int pass = 0; pass < 2 && status == MYB_OK; ++pass) {
      free(list);
      list = calloc(capacity, sizeof *list);
      if (list == NULL) {
        status = MYB_E_NOMEM;
        break;
      }
      size_t n = 0;
      status = myb_conjunction_find(a, first_start->covariance, b, second_start->covariance, radius,
                                    t1, step, list, capacity, &n);
      if (status != MYB_OK || n <= capacity) {
        *count = n;
        break;
      }
      capacity = n;
      myb_orbit_destroy(a);
      myb_orbit_destroy(b);
      a = NULL;
      b = NULL;
      status = start_orbit(first, first_start, &a);
      if (status == MYB_OK) {
        status = start_orbit(second, second_start, &b);
      }
    }
    if (status == MYB_OK) {
      *found = list;
    } else {
      free(list);
      *count = 0;
    }
  }
  myb_orbit_destroy(b);
  myb_orbit_destroy(a);
  return status;
}

myb_status mybgo_collision_long(uintptr_t handle, double t0, double t1, double step, double radius,
                                double* probability, double* error) {
  mybgo_go_context context = {.handle = handle};
  return myb_collision_probability_long(mybgo_go_gaussian, &context, t0, t1, step, radius,
                                        probability, error);
}
