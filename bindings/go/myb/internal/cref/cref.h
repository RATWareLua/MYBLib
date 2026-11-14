/* SPDX-License-Identifier: BSD-3-Clause */
/* cref.h — the C reference of the tests of the Go binding (cref.c). Every function returns the
 * status of the library. */
#ifndef MYBGO_CREF_H_
#define MYBGO_CREF_H_

#include <stddef.h>

#include "myb/myb.h"

/* A system of point masses with radii (0: none) through the times, mergers into the heavier
 * body: states 2 n each, events 7 each (kind, a, b, t, distance, survivor, merger energy). */
int cref_system_run(const double* gm, const double* radius, size_t bodies, const double* y0,
                    const double* times, size_t count, double* states, double* events,
                    size_t* event_count, double* y_end);
int cref_invariants(const double* gm, size_t bodies, const double* y, size_t count,
                    double* energy_hi, double* momentum);
int cref_system_partials(const double* gm, size_t bodies, const double* y0, double t1, double* y1,
                         double* phi, double* sensitivity);
int cref_wh(const double* gm, size_t massive, const double* x, const double* v, double h,
            const double* times, size_t count, double* xs, double* vs);
int cref_filter(double gm, const double* y0, const double* covariance, double t_end,
                const double* times, const double* positions, size_t count, double sigma,
                double* innovations, double* step_t, double* step_unknowns, size_t* steps);
int cref_range(double gm, double longitude, double latitude, double radius, const double* t,
               const double* y, size_t count, double* values, double* partials);
int cref_cr3bp_family(double mu, const double* y0, double period, size_t patches, double step,
                      size_t count, double* periods, double* states);
int cref_torus(double mu, const double* y0, double period, size_t points, double size, double step,
               double* rotations, double* normal);
int cref_manifold(double mu, const double* y0, double period, double phase, double* y);
int cref_lambert_cell(double mu, const double* from, const double* to, double departure,
                      double arrival, double* v1, double* v2);
int cref_conjunctions(double gm, const double* y1, const double* y2, const double* covariance,
                      double radius, double t1, double step, double* first_t,
                      double* first_probability, size_t* count);
int cref_qlaw_fly(double mu, const double* target, const double* weight, double thrust,
                  double exhaust, const double* y0, double t_limit, double step, double* y,
                  double* t);
int cref_chain(double mu, const double* states, const double* gm, const double* periapsis,
               const double* impulse_exhaust, const double* bounds, double* x, size_t unknowns,
               double* propellant, double* constraints, size_t constraint_count);
int cref_sims_flanagan(double mu, const double* start, const double* target, double t1,
                       size_t segments, const double* throttles, double* mismatch, double* mass1);
int cref_sun_nadir(double gm, const double* sun, double sun_gm, const double* t, const double* y,
                   size_t count, double* m);
int cref_light_time(double c, double mu, const double* target, const double* observer, double t,
                    double* light_time, double* range);
int cref_ck_pointing(const char* ck_path, int instrument, double ticks, double* q);
int cref_sclk_time(const char* kernel, long long id, double ticks, long long* s, double* f);

#endif /* MYBGO_CREF_H_ */
