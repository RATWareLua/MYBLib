/* SPDX-License-Identifier: BSD-3-Clause */
/* myb_bench — timings of the hot paths of myb, one row of a table each (bench.h):
 *
 *   kepler      Kepler's equation of the ellipse; two-body propagation over up to ten periods
 *   spk         states of DE440 (data/download/de440.bsp): the Moon from the Earth, Mars from
 *               the barycentre; without the file, a segment of tests/data/spk
 *   frames      the celestial-to-terrestrial rotation of IAU 2006/2000A, from the series and
 *               from their table (myb_cip_tabulate)
 *   field       the potential and the acceleration of a spherical-harmonic field of degree 8,
 *               70 and 360 (coefficients made up: the time does not depend on them)
 *   nbody       all accelerations of 10^3 and 10^4 bodies: the direct sum on 1 .. 4 threads
 *               against the fast multipole method (order 8, theta 0.5) given the same pool —
 *               the method itself decides whether it uses the threads
 *   integrator  a propagation of six hours of a LAGEOS-like orbit by IAS15, DOP853 and RKN12 on
 *               the force model below; the note gives the time of one step
 *   orbit       two days of that orbit by IAS15 with the edges of the shadows of the Earth and
 *               the Moon as breaks of the model, and without them
 *   estim       the least-squares fit of the initial state to six hours of positions
 *
 * The force model: the Earth's field to degree 20 (GOCO06s from data/download if there, else
 * made up) in the frame of IAU 2006/2000A, the Sun and the Moon from DE440, and the pressure
 * of sunlight on a sphere with the shadows of the Earth and the Moon. Without DE440 the
 * cases that need it say so and are left out.
 *
 *   myb_bench [--smoke] [--repeats=N] [--min-time=SECONDS] [--threads=N] [--filter=TEXT]
 *             [--data=DIR]
 *
 * --smoke is the short run of the test "bench_smoke": small sizes, short spans, few repeats,
 * no downloaded data. The orbit cases take their options from one table (orbit_variants) so
 * that a new way of handling the edges of a shadow is timed by adding a row there. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bench.h"
#include "myb/myb.h"

#ifndef MYB_DATA_DIR
#define MYB_DATA_DIR "data"
#endif

enum { SAMPLES = 1024, MAX_THREADS = 64 };

static const double pi = 3.141592653589793;
static const double gm_earth = 3.986004418e14;
static const double gm_sun = 1.32712440041279419e20;
static const double gm_moon = 4.9028000661637961e12;

/* ---- Deterministic numbers ----------------------------------------------------------------- */

static uint64_t next_random(uint64_t* state) {
  uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

/* Uniform in [0, 1). */
static double uniform(uint64_t* state) {
  return (double)(next_random(state) >> 11) * 0x1p-53;
}

static void path_of(char* out, size_t size, const char* dir, const char* file) {
  snprintf(out, size, "%s/%s", dir, file);
}

static int file_exists(const char* path) {
  FILE* f = fopen(path, "rb");
  if (f == NULL) {
    return 0;
  }
  fclose(f);
  return 1;
}

/* ---- kepler ----------------------------------------------------------------------------------- */

typedef struct kepler_case {
  double m[SAMPLES];
  double e[SAMPLES];
  myb_pv start;
  double dt[SAMPLES];
} kepler_case;

static int run_elliptic(void* context, long iterations) {
  const kepler_case* k = context;
  double sum = 0.0;
  for (long i = 0; i < iterations; ++i) {
    const size_t j = (size_t)i % SAMPLES;
    double anomaly = 0.0;
    if (myb_kepler_elliptic(k->m[j], k->e[j], &anomaly) != MYB_OK) {
      return 1;
    }
    sum += anomaly;
  }
  return isfinite(sum) ? 0 : 1;
}

static int run_propagate(void* context, long iterations) {
  const kepler_case* k = context;
  for (long i = 0; i < iterations; ++i) {
    myb_pv pv;
    if (myb_kepler_propagate(gm_earth, &k->start, k->dt[(size_t)i % SAMPLES], &pv) != MYB_OK) {
      return 1;
    }
  }
  return 0;
}

static void bench_kepler(const bench_options* o) {
  kepler_case* k = malloc(sizeof *k);
  if (k == NULL) {
    return;
  }
  uint64_t random = 1;
  for (size_t j = 0; j < SAMPLES; ++j) {
    k->m[j] = (uniform(&random) * 2.0 - 1.0) * 20.0;
    k->e[j] = 0.99 * uniform(&random);
    k->dt[j] = uniform(&random) * 10.0 * 5828.5; /* up to ten periods of the orbit below */
  }
  k->start = (myb_pv){.r = {7000e3, 0.0, 0.0}, .v = {0.0, 7546.05 * 1.1 * 0.8, 7546.05 * 0.3}};
  bench_measure(o, "kepler", "elliptic E - e sin E = M", run_elliptic, k, "e in [0, 0.99)");
  bench_measure(o, "kepler", "propagate (universal)", run_propagate, k, "dt <= 10 periods");
  free(k);
}

/* ---- spk ------------------------------------------------------------------------------------- */

typedef struct spk_case {
  myb_spk* spk;
  int32_t target;
  int32_t center;
  double t[SAMPLES];
} spk_case;

static int run_spk(void* context, long iterations) {
  const spk_case* c = context;
  for (long i = 0; i < iterations; ++i) {
    double state[6];
    if (myb_spk_state(c->spk, c->target, c->center,
                      myb_time_from_seconds(c->t[(size_t)i % SAMPLES]), state, NULL) != MYB_OK) {
      return 1;
    }
  }
  return 0;
}

static void bench_spk(const bench_options* o) {
  if (!bench_selected(o, "spk", "")) {
    return;
  }
  spk_case* c = malloc(sizeof *c);
  if (c == NULL || myb_spk_create(NULL, &c->spk) != MYB_OK) {
    free(c);
    return;
  }
  char path[512];
  path_of(path, sizeof path, o->data, "download/de440.bsp");
  uint64_t random = 2;
  if (!o->smoke && file_exists(path) && myb_spk_load_file(c->spk, path, NULL) == MYB_OK) {
    for (size_t j = 0; j < SAMPLES; ++j) {
      c->t[j] = (uniform(&random) * 2.0 - 1.0) * 1.5e9; /* 1952 .. 2047 */
    }
    c->target = 301;
    c->center = 399;
    bench_measure(o, "spk", "DE440 Moon from Earth", run_spk, c, "two chains of segments");
    c->target = 4;
    c->center = 0;
    bench_measure(o, "spk", "DE440 Mars barycentre", run_spk, c, "one segment");
  } else {
    path_of(path, sizeof path, o->data, "../tests/data/spk/t02_03.bsp");
    myb_spk_segment_info info;
    if (myb_spk_load_file(c->spk, path, NULL) == MYB_OK && myb_spk_segment_count(c->spk) > 0) {
      myb_spk_segment(c->spk, 0, &info);
      for (size_t j = 0; j < SAMPLES; ++j) {
        c->t[j] = info.start + uniform(&random) * (info.stop - info.start);
      }
      c->target = info.target;
      c->center = info.center;
      bench_measure(o, "spk", "test segment (no DE440)", run_spk, c, "tests/data/spk/t02_03");
    }
  }
  myb_spk_destroy(c->spk);
  free(c);
}

/* ---- frames ------------------------------------------------------------------------------------ */

typedef struct frames_case {
  myb_cip* cip;
  myb_time t[SAMPLES];
} frames_case;

static int run_c2t(void* context, long iterations) {
  const frames_case* c = context;
  for (long i = 0; i < iterations; ++i) {
    const myb_time tt = c->t[(size_t)i % SAMPLES];
    myb_rotation rotation;
    myb_cip_terrestrial(c->cip, tt, myb_time_add(tt, -69.184), NULL, &rotation);
  }
  return 0;
}

static myb_cip* load_cip(const char* data) {
  char path[512];
  path_of(path, sizeof path, data, "frames/iau2006-cip.tk");
  myb_kpool* pool = NULL;
  myb_cip* cip = NULL;
  if (myb_kpool_create(NULL, &pool) == MYB_OK && myb_kpool_load_file(pool, path, NULL) == MYB_OK &&
      myb_cip_create(pool, "IAU2006A", NULL, &cip, NULL) != MYB_OK) {
    cip = NULL;
  }
  myb_kpool_destroy(pool);
  return cip;
}

static void bench_frames(const bench_options* o) {
  if (!bench_selected(o, "frames", "")) {
    return;
  }
  frames_case* c = malloc(sizeof *c);
  if (c == NULL) {
    return;
  }
  c->cip = load_cip(o->data);
  if (c->cip == NULL) {
    printf("| frames | (no data/frames/iau2006-cip.tk) | | | | | | |\n");
    free(c);
    return;
  }
  uint64_t random = 3;
  const double span = 30.0 * 86400.0;
  for (size_t j = 0; j < SAMPLES; ++j) {
    c->t[j] = myb_time_from_seconds(8.0e8 + uniform(&random) * span);
  }
  bench_measure(o, "frames", "celestial to terrestrial, series", run_c2t, c,
                "IAU 2006/2000A, no EOP");
  if (myb_cip_tabulate(c->cip, myb_time_from_seconds(8.0e8 - 3600.0),
                       myb_time_from_seconds(8.0e8 + span + 3600.0), 3600.0) == MYB_OK) {
    bench_measure(o, "frames", "celestial to terrestrial, table", run_c2t, c,
                  "myb_cip_tabulate, 1 h");
  }
  myb_cip_destroy(c->cip);
  free(c);
}

/* ---- field ----------------------------------------------------------------------------------- */

typedef struct field_case {
  myb_gravity_field* field;
  myb_gravity_work* work;
  double r[SAMPLES][3];
} field_case;

/* A field of the given degree with coefficients of the size of the Earth's (made up). */
static myb_gravity_field* made_up_field(int degree, uint64_t seed) {
  const size_t count = (size_t)(degree + 1) * (size_t)(degree + 2) / 2;
  double* c = calloc(count, sizeof(double));
  double* s = calloc(count, sizeof(double));
  myb_gravity_field* field = NULL;
  if (c != NULL && s != NULL) {
    uint64_t random = seed;
    c[0] = 1.0;
    for (int n = 2; n <= degree; ++n) {
      const double size = 1e-5 / ((double)n * (double)n);
      for (int m = 0; m <= n; ++m) {
        const size_t i = (size_t)n * (size_t)(n + 1) / 2 + (size_t)m;
        c[i] = size * (uniform(&random) * 2.0 - 1.0);
        s[i] = m == 0 ? 0.0 : size * (uniform(&random) * 2.0 - 1.0);
      }
    }
    if (degree >= 2) {
      c[3] = -4.841652e-4; /* C20 */
    }
    if (myb_gravity_field_create(gm_earth, 6378136.3, degree, c, s, NULL, &field) != MYB_OK) {
      field = NULL;
    }
  }
  free(s);
  free(c);
  return field;
}

static int run_field(void* context, long iterations) {
  field_case* c = context;
  for (long i = 0; i < iterations; ++i) {
    double potential = 0.0;
    double a[3];
    if (myb_gravity_field_eval(c->field, c->work, c->r[(size_t)i % SAMPLES], 0, -1, &potential,
                               a) != MYB_OK) {
      return 1;
    }
  }
  return 0;
}

static void bench_field(const bench_options* o) {
  if (!bench_selected(o, "field", "")) {
    return;
  }
  static const int degrees[3] = {8, 70, 360};
  const int count = o->smoke ? 2 : 3;
  field_case* c = malloc(sizeof *c);
  if (c == NULL) {
    return;
  }
  uint64_t random = 4;
  for (size_t j = 0; j < SAMPLES; ++j) {
    const double z = uniform(&random) * 2.0 - 1.0;
    const double phi = uniform(&random) * 2.0 * pi;
    const double rho = sqrt(1.0 - z * z);
    const double radius = 6.778e6 + uniform(&random) * 3.0e5;
    c->r[j][0] = radius * rho * myb_cos(phi);
    c->r[j][1] = radius * rho * myb_sin(phi);
    c->r[j][2] = radius * z;
  }
  for (int k = 0; k < count; ++k) {
    c->field = made_up_field(degrees[k], 5 + (uint64_t)k);
    c->work = NULL;
    if (c->field != NULL && myb_gravity_work_create(degrees[k], NULL, &c->work) == MYB_OK) {
      char name[64];
      snprintf(name, sizeof name, "potential and acceleration, degree %d", degrees[k]);
      bench_measure(o, "field", name, run_field, c, "LEO, any latitude");
    }
    myb_gravity_work_destroy(c->work);
    myb_gravity_field_destroy(c->field);
  }
  free(c);
}

/* ---- nbody ----------------------------------------------------------------------------------- */

typedef struct nbody_case {
  myb_nbody* nbody;
  double* x;
  double* a;
} nbody_case;

static int run_nbody(void* context, long iterations) {
  nbody_case* c = context;
  for (long i = 0; i < iterations; ++i) {
    if (myb_nbody_acceleration(c->nbody, 0.0, c->x, c->a) != MYB_OK) {
      return 1;
    }
  }
  return 0;
}

static void bench_nbody(const bench_options* o) {
  if (!bench_selected(o, "nbody", "")) {
    return;
  }
  const size_t sizes[2] = {o->smoke ? 256 : 1000, o->smoke ? 1024 : 10000};
  int most = o->threads < 1 ? 1 : o->threads > MAX_THREADS ? MAX_THREADS : o->threads;
  if (o->smoke && most > 2) {
    most = 2;
  }
  for (int k = 0; k < 2; ++k) {
    const size_t n = sizes[k];
    nbody_case c = {.nbody = NULL, .x = malloc(3 * n * sizeof(double)),
                    .a = malloc(3 * n * sizeof(double))};
    double* gm = malloc(n * sizeof(double));
    if (c.x == NULL || c.a == NULL || gm == NULL) {
      free(gm);
      free(c.a);
      free(c.x);
      return;
    }
    /* a uniform ball of 1 au: no close pairs to speak of, the cost the sums themselves */
    uint64_t random = 6;
    for (size_t i = 0; i < n; ++i) {
      double p[3];
      double r2 = 2.0;
      /* bound: rejection sampling, accepts with probability pi / 6 */
      for (int tries = 0; tries < 64 && r2 > 1.0; ++tries) {
        for (int j = 0; j < 3; ++j) {
          p[j] = uniform(&random) * 2.0 - 1.0;
        }
        r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
      }
      for (int j = 0; j < 3; ++j) {
        c.x[3 * i + (size_t)j] = 1.496e11 * p[j];
      }
      gm[i] = 1e12 * (1.0 + uniform(&random));
    }
    if (myb_nbody_create(n, 0, gm, NULL, &c.nbody) == MYB_OK) {
      for (int method = 0; method < 2; ++method) {
        if (method == 1 && myb_nbody_set_approximation(c.nbody, 8, 0.5) != MYB_OK) {
          break;
        }
        for (int threads = 1; threads <= most; ++threads) {
          myb_pool* pool = NULL;
          if (threads > 1 && myb_pool_create(threads, NULL, &pool) != MYB_OK) {
            break;
          }
          myb_nbody_set_pool(c.nbody, pool);
          char name[64];
          char note[64];
          snprintf(name, sizeof name, "%zu bodies, %s, %d thread%s", n,
                   method == 0 ? "direct" : "FMM p 8", threads, threads > 1 ? "s" : "");
          snprintf(note, sizeof note, "one evaluation of all %zu", n);
          bench_measure(o, "nbody", name, run_nbody, &c, note);
          myb_nbody_set_pool(c.nbody, NULL);
          myb_pool_destroy(pool);
        }
      }
      myb_nbody_destroy(c.nbody);
    }
    free(gm);
    free(c.a);
    free(c.x);
  }
}

/* ---- The force model ---------------------------------------------------------------------------- */

typedef struct force_setup {
  myb_timesys* timesys;
  myb_cip* cip;
  myb_spk* spk;
  myb_env* env;
  myb_force_model* model;
  myb_gravity_field* field;
  myb_gravity_work* work;
  myb_spk_body sun_state;
  myb_spk_body moon_state;
  myb_object object;
  myb_field_term gravity;
  myb_third_body_term third[2];
  myb_radiation_term light;
  myb_eclipse_event edges[4];
  int bodies; /* the Sun and the Moon are in the model */
  int breaks; /* the edges of the shadows are breaks of the model */
} force_setup;

/* The Earth's orientation without EOP: UT1 = TT - 69.184 s, no polar motion. */
static myb_status earth_orientation(void* context, myb_time t, myb_orientation* orientation) {
  const myb_cip* cip = context;
  const myb_time ut1 = myb_time_add(t, -69.184);
  myb_cip_terrestrial(cip, t, ut1, NULL, &orientation->rotation);
  orientation->ut = ut1;
  orientation->xp = 0.0;
  orientation->yp = 0.0;
  return MYB_OK;
}

static void teardown(force_setup* s) {
  myb_force_model_destroy(s->model);
  myb_env_destroy(s->env);
  myb_gravity_work_destroy(s->work);
  myb_gravity_field_destroy(s->field);
  myb_spk_destroy(s->spk);
  myb_cip_destroy(s->cip);
  myb_timesys_destroy(s->timesys);
  memset(s, 0, sizeof *s);
}

static myb_timesys* load_timesys(const char* data) {
  static const char* const files[3] = {"time/time-scales.tk", "time/tdb-tt-fb1990.tk",
                                       "time/utc-tai.tk"};
  myb_kpool* pool = NULL;
  myb_timesys* system = NULL;
  int ok = myb_kpool_create(NULL, &pool) == MYB_OK;
  for (int k = 0; k < 3 && ok; ++k) {
    char path[512];
    path_of(path, sizeof path, data, files[k]);
    ok = myb_kpool_load_file(pool, path, NULL) == MYB_OK;
  }
  if (!ok || myb_timesys_create(pool, NULL, &system, NULL) != MYB_OK) {
    system = NULL;
  }
  myb_kpool_destroy(pool);
  return system;
}

/* Builds the model of the integrator, orbit and estim cases at 2025-06-01 TT over `span`
 * seconds; with ephemerides = 0 the field alone. Returns 0 on failure. */
static int setup_forces(force_setup* s, const bench_options* o, int ephemerides, int breaks,
                        double span) {
  memset(s, 0, sizeof *s);
  const myb_calendar epoch_label = {2025, 6, 1, 0, 0, 0.0};
  myb_time epoch;
  if (myb_calendar_to_time(&epoch_label, &epoch) != MYB_OK) {
    return 0;
  }
  s->timesys = load_timesys(o->data);
  s->cip = load_cip(o->data);
  if (s->timesys == NULL || s->cip == NULL ||
      myb_cip_tabulate(s->cip, myb_time_add(epoch, -7200.0), myb_time_add(epoch, span + 7200.0),
                       3600.0) != MYB_OK) {
    teardown(s);
    return 0;
  }
  char path[512];
  path_of(path, sizeof path, o->data, "download/GOCO06s.gfc");
  const int degree = o->smoke ? 8 : 20;
  if (o->smoke || !file_exists(path) ||
      myb_gravity_field_load_icgem_file(path, degree, NULL, &s->field, NULL) != MYB_OK) {
    s->field = made_up_field(degree, 7);
  }
  const myb_env_options options = {.epoch = epoch, .timesys = s->timesys, .scale = "TT"};
  if (s->field == NULL || myb_gravity_work_create(degree, NULL, &s->work) != MYB_OK ||
      myb_env_create(&options, NULL, &s->env) != MYB_OK ||
      myb_force_model_create(NULL, &s->model) != MYB_OK) {
    teardown(s);
    return 0;
  }
  const myb_env_body earth = {.id = 399,
                              .gm = gm_earth,
                              .orientation = earth_orientation,
                              .orientation_context = s->cip,
                              .radius = 6378137.0};
  size_t index[3] = {0, 0, 0};
  int ok = myb_env_add_body(s->env, &earth, &index[0]) == MYB_OK;
  s->gravity = (myb_field_term){.env = s->env, .body = index[0], .field = s->field,
                                .work = s->work, .n_min = 0, .n_max = -1};
  size_t term = 0;
  ok = ok && myb_force_model_add(s->model, "earth", myb_term_field, myb_term_field_jacobian,
                                 &s->gravity, &term) == MYB_OK;
  myb_force_model_set_position_only(s->model, term, 1);
  path_of(path, sizeof path, o->data, "download/de440.bsp");
  if (ok && ephemerides && file_exists(path) && myb_spk_create(NULL, &s->spk) == MYB_OK &&
      myb_spk_load_file(s->spk, path, NULL) == MYB_OK) {
    s->sun_state = (myb_spk_body){.spk = s->spk, .target = 10, .center = 399, .frame = 1};
    s->moon_state = (myb_spk_body){.spk = s->spk, .target = 301, .center = 399, .frame = 1};
    const myb_env_body sun = {.id = 10, .gm = gm_sun, .state = myb_spk_body_state,
                              .state_context = &s->sun_state, .state_scale = "TDB",
                              .radius = 695700000.0, .luminosity = 3.828e26,
                              .limb = {0.6, 0.0}};
    const myb_env_body moon = {.id = 301, .gm = gm_moon, .state = myb_spk_body_state,
                               .state_context = &s->moon_state, .state_scale = "TDB",
                               .radius = 1737400.0};
    ok = myb_env_add_body(s->env, &sun, &index[1]) == MYB_OK &&
         myb_env_add_body(s->env, &moon, &index[2]) == MYB_OK;
    static const char* const names[2] = {"sun", "moon"};
    for (int b = 0; b < 2 && ok; ++b) {
      s->third[b] = (myb_third_body_term){s->env, index[1 + b]};
      ok = myb_force_model_add(s->model, names[b], myb_term_third_body,
                               myb_term_third_body_jacobian, &s->third[b], &term) == MYB_OK;
      myb_force_model_set_position_only(s->model, term, 1);
    }
    s->object = (myb_object){.mass = 406.965, .area = pi * 0.3 * 0.3, .cr = 1.13};
    s->light = (myb_radiation_term){.env = s->env, .source = index[1], .object = &s->object,
                                    .c = 299792458.0, .occulters = {index[0], index[2]},
                                    .occulter_count = 2};
    ok = ok && myb_force_model_add(s->model, "sunlight", myb_term_radiation,
                                   myb_term_radiation_jacobian, &s->light, &term) == MYB_OK;
    myb_force_model_set_position_only(s->model, term, 1);
    if (ok && breaks) {
      ok = myb_term_radiation_breaks(&s->light, s->model, s->edges) == MYB_OK;
      s->breaks = ok;
    }
    s->bodies = ok;
  }
  if (!ok) {
    teardown(s);
  }
  return ok;
}

/* A LAGEOS-like orbit: a = 12 270 km, circular, inclination 109.84 degrees, the ascending node
 * at 78.4 degrees, where on 2025-06-01 the Sun lies near the plane of the orbit: the Earth
 * eclipses it once a revolution (3.8 h), four edges of shadow each time. */
static void lageos_state(double y[6]) {
  const double a = 12270e3;
  const double v = sqrt(gm_earth / a);
  const double i = 109.84 * pi / 180.0;
  const double node = 78.4 * pi / 180.0;
  const double cn = myb_cos(node);
  const double sn = myb_sin(node);
  const double state[6] = {a * cn, a * sn, 0.0, -v * sn * myb_cos(i), v * cn * myb_cos(i),
                           v * myb_sin(i)};
  memcpy(y, state, sizeof state);
}

/* ---- integrator and orbit ---------------------------------------------------------------------- */

typedef struct orbit_case {
  force_setup* forces;
  myb_orbit_options options;
  double span;
  int64_t steps;
  int64_t breaks;
} orbit_case;

static int run_orbit(void* context, long iterations) {
  orbit_case* c = context;
  for (long i = 0; i < iterations; ++i) {
    myb_orbit* orbit = NULL;
    double y[6];
    lageos_state(y);
    if (myb_orbit_create(c->forces->model, &c->options, NULL, &orbit) != MYB_OK) {
      return 1;
    }
    const int ok = myb_orbit_start(orbit, 0.0, y, c->span) == MYB_OK &&
                   myb_orbit_state(orbit, c->span, y) == MYB_OK;
    c->steps = myb_orbit_steps(orbit);
    c->breaks = myb_orbit_breaks(orbit);
    myb_orbit_destroy(orbit);
    if (!ok) {
      return 1;
    }
  }
  return 0;
}

/* Times one propagation and notes its steps and the time of one. */
static void measure_orbit(const bench_options* o, const char* group, const char* name,
                          orbit_case* c) {
  if (!bench_selected(o, group, name)) {
    return;
  }
  /* one run first for the counts of the note */
  if (run_orbit(c, 1) != 0) {
    printf("| %s | %s | failed | | | | | |\n", group, name);
    return;
  }
  char note[96];
  snprintf(note, sizeof note, "%lld steps, %lld breaks", (long long)c->steps,
           (long long)c->breaks);
  const double median = bench_measure(o, group, name, run_orbit, c, note);
  if (median > 0.0 && c->steps > 0) {
    char per_step[16];
    bench_format_time(median / (double)c->steps, per_step);
    printf("| %s | (the same, per step) | %s | | | | | median / %lld steps |\n", group, per_step,
           (long long)c->steps);
  }
}

static void bench_integrators(const bench_options* o) {
  if (!bench_selected(o, "integrator", "")) {
    return;
  }
  force_setup* s = malloc(sizeof *s);
  const double span = o->smoke ? 7200.0 : 6.0 * 3600.0;
  if (s == NULL || !setup_forces(s, o, !o->smoke, 0, span)) {
    printf("| integrator | (no model: data/ missing) | | | | | | |\n");
    free(s);
    return;
  }
  static const struct {
    const char* name;
    myb_orbit_integrator integrator;
  } integrators[3] = {{"IAS15, epsilon 1e-10", MYB_ORBIT_IAS15},
                      {"DOP853, rtol 1e-12", MYB_ORBIT_DOP853},
                      {"RKN12, rtol 1e-12", MYB_ORBIT_RKN12}};
  for (int k = 0; k < 3; ++k) {
    orbit_case c = {.forces = s,
                    .options = {.epsilon = 1e-10, .integrator = integrators[k].integrator},
                    .span = span};
    char name[96];
    snprintf(name, sizeof name, "%s, %s", integrators[k].name,
             s->bodies ? "field 20 + Sun, Moon, sunlight" : "field alone");
    measure_orbit(o, "integrator", name, &c);
  }
  teardown(s);
  free(s);
}

/* The ways of propagating across the edges of a shadow that are timed: add a row for a new
 * one (an option of myb_orbit_options, or a model built differently). */
typedef struct orbit_variant {
  const char* name;
  int breaks;                 /* the edges are breaks of the model */
  myb_orbit_options options;
} orbit_variant;

static const orbit_variant orbit_variants[] = {
    {"IAS15, edges as breaks (the step is taken again up to the edge)", 1,
     {.epsilon = 1e-10}},
    {"IAS15, no breaks (the steps cross the edges)", 0, {.epsilon = 1e-10}},
};

static void bench_orbit(const bench_options* o) {
  if (!bench_selected(o, "orbit", "") || o->smoke) {
    return;
  }
  const double span = 2.0 * 86400.0;
  for (size_t k = 0; k < sizeof orbit_variants / sizeof orbit_variants[0]; ++k) {
    force_setup* s = malloc(sizeof *s);
    if (s == NULL) {
      return;
    }
    if (!setup_forces(s, o, 1, orbit_variants[k].breaks, span) || !s->bodies) {
      printf("| orbit | (needs data/download/de440.bsp) | | | | | | |\n");
      if (s->model != NULL) {
        teardown(s);
      }
      free(s);
      return;
    }
    orbit_case c = {.forces = s, .options = orbit_variants[k].options, .span = span};
    measure_orbit(o, "orbit", orbit_variants[k].name, &c);
    teardown(s);
    free(s);
  }
}

/* ---- estim ----------------------------------------------------------------------------------- */

enum { OBSERVATIONS = 24 };

typedef struct fit_case {
  force_setup* forces;
  myb_observation observations[OBSERVATIONS];
  double truth[6];
  int iterations;
} fit_case;

static int run_fit(void* context, long iterations) {
  fit_case* c = context;
  for (long i = 0; i < iterations; ++i) {
    double y0[6];
    memcpy(y0, c->truth, sizeof y0);
    const double offset[6] = {100.0, -50.0, 30.0, 0.05, -0.02, 0.01};
    for (int j = 0; j < 6; ++j) {
      y0[j] += offset[j];
    }
    myb_fit_result result;
    const myb_fit_options options = {.max_iterations = 10, .epsilon = 1e-10};
    if (myb_fit_orbit(c->forces->model, 0.0, y0, c->observations, OBSERVATIONS, &options, NULL,
                      &result, NULL, NULL) != MYB_OK) {
      return 1;
    }
    c->iterations = result.iterations;
  }
  return 0;
}

static void bench_estim(const bench_options* o) {
  if (!bench_selected(o, "estim", "")) {
    return;
  }
  fit_case* c = malloc(sizeof *c);
  force_setup* s = malloc(sizeof *s);
  const double span = o->smoke ? 3600.0 : 6.0 * 3600.0;
  if (c == NULL || s == NULL || !setup_forces(s, o, !o->smoke, 0, span)) {
    printf("| estim | (no model: data/ missing) | | | | | | |\n");
    free(s);
    free(c);
    return;
  }
  c->forces = s;
  lageos_state(c->truth);
  myb_orbit* orbit = NULL;
  double times[OBSERVATIONS];
  double states[6 * OBSERVATIONS];
  for (int k = 0; k < OBSERVATIONS; ++k) {
    times[k] = span * (double)(k + 1) / (double)OBSERVATIONS;
  }
  const myb_orbit_options options = {.epsilon = 1e-10};
  int ok = myb_orbit_create(s->model, &options, NULL, &orbit) == MYB_OK &&
           myb_orbit_start(orbit, 0.0, c->truth, span * 1.01) == MYB_OK &&
           myb_orbit_states(orbit, times, OBSERVATIONS, states) == MYB_OK;
  myb_orbit_destroy(orbit);
  /* positions with errors of about a metre (the sum of 12 uniform numbers less 6: nearly
   * normal, sigma 1), so that the fit stops at the noise as a real one does */
  uint64_t random = 9;
  for (int k = 0; k < OBSERVATIONS && ok; ++k) {
    c->observations[k] = (myb_observation){.t = times[k], .measure = myb_measure_position,
                                           .count = 3, .sigma = {1.0, 1.0, 1.0}};
    for (int j = 0; j < 3; ++j) {
      double noise = -6.0;
      for (int q = 0; q < 12; ++q) {
        noise += uniform(&random);
      }
      c->observations[k].value[j] = states[6 * k + j] + noise;
    }
  }
  if (ok && run_fit(c, 1) == 0) {
    char note[96];
    snprintf(note, sizeof note, "%d positions over %.0f h, %d iterations, %s", OBSERVATIONS,
             span / 3600.0, c->iterations, s->bodies ? "full model" : "field alone");
    bench_measure(o, "estim", "least squares, initial state", run_fit, c, note);
  } else {
    printf("| estim | least squares, initial state | failed | | | | | |\n");
  }
  teardown(s);
  free(s);
  free(c);
}

/* ---- Main -------------------------------------------------------------------------------------- */

static int option_value(const char* arg, const char* name, const char** value) {
  const size_t n = strlen(name);
  if (strncmp(arg, name, n) != 0) {
    return 0;
  }
  *value = arg + n;
  return 1;
}

int main(int argc, char** argv) {
  bench_options o = {.repeats = 7, .min_seconds = 0.2, .smoke = 0, .filter = NULL,
                     .threads = 4, .data = MYB_DATA_DIR};
  /* bound: the arguments */
  for (int i = 1; i < argc; ++i) {
    const char* value = NULL;
    if (strcmp(argv[i], "--smoke") == 0) {
      o.smoke = 1;
    } else if (option_value(argv[i], "--repeats=", &value)) {
      o.repeats = atoi(value);
    } else if (option_value(argv[i], "--min-time=", &value)) {
      o.min_seconds = atof(value);
    } else if (option_value(argv[i], "--threads=", &value)) {
      o.threads = atoi(value);
    } else if (option_value(argv[i], "--filter=", &value)) {
      o.filter = value;
    } else if (option_value(argv[i], "--data=", &value)) {
      o.data = value;
    } else {
      fprintf(stderr, "usage: %s [--smoke] [--repeats=N] [--min-time=SECONDS] [--threads=N] "
                      "[--filter=TEXT] [--data=DIR]\n", argv[0]);
      return 2;
    }
  }
  if (o.smoke) {
    o.repeats = 3;
    o.min_seconds = 0.005;
  }
  printf("myb %ld bench: median of %d repeats of at least %.3f s; %d threads at most; "
         "processors %d\n\n", myb_version_number(), o.repeats, o.min_seconds, o.threads,
         myb_platform_processor_count());
  bench_header();
  bench_kepler(&o);
  bench_spk(&o);
  bench_frames(&o);
  bench_field(&o);
  bench_nbody(&o);
  bench_integrators(&o);
  bench_orbit(&o);
  bench_estim(&o);
  return 0;
}
