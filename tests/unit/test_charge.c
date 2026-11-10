/* SPDX-License-Identifier: BSD-3-Clause */
/* Charge and plasma: IGRF-14 from the NOAA table against ppigrf (independent code, 160 points,
 * 2000-2025); the field of g10 against the dipole formula; a charge in uniform fields against the
 * exact helix and the E x B drift; in a dipole, the kinetic energy and the canonical angular
 * momentum p_phi = rho v_phi + (q/m) rho A_phi conserved; co-rotation field; plasma stream. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtest.h"
#include "myb/force.h"
#include "myb/grav.h"
#include "myb/integ.h"

#define DOWNLOAD MYB_DATA_DIR "/download/"

static const double PI = 3.141592653589793;

static double norm(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static char* read_all(const char* path, size_t* size) {
  FILE* f = fopen(path, "rb");
  if (f == NULL) {
    return NULL;
  }
  fseek(f, 0, SEEK_END);
  const long n = ftell(f);
  fseek(f, 0, SEEK_SET);
  char* text = malloc((size_t)n + 1);
  if (text != NULL && fread(text, 1, (size_t)n, f) != (size_t)n) {
    free(text);
    text = NULL;
  }
  fclose(f);
  if (text != NULL) {
    text[n] = '\0';
    *size = (size_t)n;
  }
  return text;
}

static void test_igrf(void) {
  size_t size = 0;
  char* table = read_all(DOWNLOAD "igrf14coeffs.txt", &size);
  if (table == NULL) {
    printf("igrf14coeffs.txt not downloaded: IGRF skipped\n");
    return;
  }
  FILE* f = fopen(MYB_TEST_DATA_DIR "/magnetic/igrf_reference.txt", "r");
  MT_CHECK(f != NULL);
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_work_create(13, NULL, &work) == MYB_OK);
  char line[512];
  double worst = 0.0;
  int rows = 0;
  double current_year = -1.0;
  myb_magnetic_field* field = NULL;
  while (f != NULL && fgets(line, sizeof line, f) != NULL) {
    if (line[0] == '#') {
      continue;
    }
    double year, r, th, ph, br, bt, bp;
    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &year, &r, &th, &ph, &br, &bt, &bp) != 7) {
      continue;
    }
    if (year != current_year) {
      myb_magnetic_field_destroy(field);
      field = NULL;
      myb_diag diag = {0};
      MT_CHECK(myb_magnetic_field_load_table(table, size, 6371.2e3, year, NULL, &field, &diag) ==
               MYB_OK);
      current_year = year;
      MT_CHECK(field != NULL && myb_magnetic_field_degree(field) == 13);
    }
    if (field == NULL) {
      continue;
    }
    const double t = th * PI / 180.0;
    const double p = ph * PI / 180.0;
    const double x[3] = {r * 1e3 * sin(t) * cos(p), r * 1e3 * sin(t) * sin(p), r * 1e3 * cos(t)};
    double b[3];
    MT_CHECK(myb_magnetic_field_eval(field, work, x, b) == MYB_OK);
    /* spherical components (nT): r, theta (south), phi (east) */
    const double er[3] = {sin(t) * cos(p), sin(t) * sin(p), cos(t)};
    const double et[3] = {cos(t) * cos(p), cos(t) * sin(p), -sin(t)};
    const double ep[3] = {-sin(p), cos(p), 0.0};
    const double ours[3] = {1e9 * (b[0] * er[0] + b[1] * er[1] + b[2] * er[2]),
                            1e9 * (b[0] * et[0] + b[1] * et[1] + b[2] * et[2]),
                            1e9 * (b[0] * ep[0] + b[1] * ep[1] + b[2] * ep[2])};
    const double ref[3] = {br, bt, bp};
    const double scale = sqrt(br * br + bt * bt + bp * bp);
    for (int i = 0; i < 3; ++i) {
      worst = fmax(worst, fabs(ours[i] - ref[i]) / scale);
    }
    ++rows;
  }
  if (f != NULL) {
    fclose(f);
  }
  printf("IGRF-14 vs ppigrf (%d points, 2000-2025, 1-4 Earth radii): %.2e of |B|\n", rows, worst);
  MT_CHECK(rows == 160);
  MT_CHECK(worst < 1e-12);
  /* errors */
  myb_magnetic_field* none = NULL;
  myb_diag diag = {0};
  MT_CHECK(myb_magnetic_field_load_table(table, size, 6371.2e3, 1850.0, NULL, &none, &diag) ==
           MYB_E_RANGE);
  MT_CHECK(myb_magnetic_field_load_table("g 1 0 1.0\n", 10, 6371.2e3, 2000.0, NULL, &none,
                                         &diag) == MYB_E_FORMAT);
  myb_magnetic_field_destroy(field);
  myb_gravity_work_destroy(work);
  free(table);
}

static void test_dipole(void) {
  /* g10 only: B = (a^3 / r^3) g10 (3 (z.rhat) rhat - z) ... with the sign of V = a (a/r)^2 g10 cos */
  const double a = 6.0e7;
  double g[3] = {0.0, 4.2e-4, 0.0};
  double h[3] = {0.0, 0.0, 0.0};
  myb_magnetic_field* field = NULL;
  MT_CHECK(myb_magnetic_field_create(a, 1, g, h, NULL, &field) == MYB_OK);
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_work_create(1, NULL, &work) == MYB_OK);
  double worst = 0.0;
  mt_rng rng = {5};
  for (int k = 0; k < 200; ++k) {
    double x[3];
    for (int i = 0; i < 3; ++i) {
      x[i] = ((double)(mt_next(&rng) >> 11) * 0x1p-53 - 0.5) * 8.0 * a;
    }
    if (norm(x) < a) {
      continue;
    }
    double b[3];
    MT_CHECK(myb_magnetic_field_eval(field, work, x, b) == MYB_OK);
    /* V = g10 a^3 z / r^3, B = -grad V = g10 a^3 (3 z x / r^5 - zhat / r^3) */
    const double r = norm(x);
    const double r3 = r * r * r;
    const double expected[3] = {g[1] * a * a * a * 3.0 * x[2] * x[0] / (r3 * r * r),
                                g[1] * a * a * a * 3.0 * x[2] * x[1] / (r3 * r * r),
                                g[1] * a * a * a * (3.0 * x[2] * x[2] / (r * r) - 1.0) / r3};
    for (int i = 0; i < 3; ++i) {
      worst = fmax(worst, fabs(b[i] - expected[i]) / norm(expected));
    }
  }
  printf("dipole from g10 vs formula: %.2e\n", worst);
  MT_CHECK(worst < 1e-15);
  myb_gravity_work_destroy(work);
  myb_magnetic_field_destroy(field);
}

/* ---- Motion of a charge ---------------------------------------------------------------------- */

typedef struct charge {
  double qm;
  double e[3];
  double b[3];
  const myb_magnetic_field* field; /* NULL: uniform b */
  myb_gravity_work* work;
} charge;

static myb_status charge_accel(void* context, double t, const double* x, const double* v,
                               double* a) {
  (void)t;
  const charge* c = context;
  double b[3] = {c->b[0], c->b[1], c->b[2]};
  if (c->field != NULL) {
    const myb_status st = myb_magnetic_field_eval(c->field, c->work, x, b);
    if (st != MYB_OK) {
      return st;
    }
  }
  myb_lorentz_acceleration(c->qm, c->e, v, b, a);
  return MYB_OK;
}

static void integrate(charge* c, const double y0[6], double t1, double y[6]) {
  const myb_ias15_options o = {.epsilon = 1e-9};
  myb_ias15* d = NULL;
  MT_CHECK(myb_ias15_create(3, &o, NULL, &d) == MYB_OK);
  MT_CHECK(myb_ias15_start(d, charge_accel, c, 0.0, y0) == MYB_OK);
  MT_CHECK(myb_ias15_integrate(d, t1) == MYB_OK);
  memcpy(y, myb_ias15_state(d), 6 * sizeof *y);
  myb_ias15_destroy(d);
}

static void test_uniform(void) {
  /* helix: B along z, q/m, gyrofrequency w = (q/m) B, x = rL sin(w t) ..., 50 turns */
  charge c = {.qm = 2.0, .b = {0.0, 0.0, 0.5}};
  const double w = c.qm * c.b[2];
  const double vperp = 3.0;
  const double vpar = 0.7;
  const double y0[6] = {0.0, 0.0, 0.0, vperp, 0.0, vpar};
  const double t1 = 50.0 * 2.0 * PI / w;
  double y[6];
  integrate(&c, y0, t1, y);
  /* dv/dt = w v x z: vx = vperp cos(w t), vy = -vperp sin(w t) */
  const double ex[3] = {vperp / w * sin(w * t1), vperp / w * (cos(w * t1) - 1.0), vpar * t1};
  double err = 0.0;
  for (int i = 0; i < 3; ++i) {
    err = fmax(err, fabs(y[i] - ex[i]));
  }
  printf("uniform B, 50 gyrations: position vs the helix %.2e (radius %.3f)\n", err, vperp / w);
  MT_CHECK(err < 1e-11);
  /* E x B drift: E along x, B along z: guiding centre at v_d = E x B / B^2 = (0, -E/B, 0) */
  charge d = {.qm = 2.0, .e = {0.2, 0.0, 0.0}, .b = {0.0, 0.0, 0.5}};
  const double vd = -d.e[0] / d.b[2];
  const double y1[6] = {0.0, 0.0, 0.0, 0.0, vd, 0.0}; /* started on the drift: no gyration */
  integrate(&d, y1, t1, y);
  const double drift_err = fmax(fabs(y[1] - vd * t1), fabs(y[0]));
  printf("E x B drift: %.2e m off the straight drift line over %.1f s\n", drift_err, t1);
  MT_CHECK(drift_err < 1e-10);
}

static void test_dipole_invariants(void) {
  /* a proton-like q/m in a dipole of Jupiter scale: the magnetic force does no work, and the
   * field is axisymmetric: |v| and p_phi = (x vy - y vx) + (q/m) rho A_phi are constant */
  const double a = 7.0e7;
  double g[3] = {0.0, 4.2e-4, 0.0};
  double h[3] = {0.0, 0.0, 0.0};
  myb_magnetic_field* field = NULL;
  MT_CHECK(myb_magnetic_field_create(a, 1, g, h, NULL, &field) == MYB_OK);
  myb_gravity_work* work = NULL;
  MT_CHECK(myb_gravity_work_create(1, NULL, &work) == MYB_OK);
  charge c = {.qm = 9.58e7 * 1e-6, .field = field, .work = work};
  const double y0[6] = {3.0 * a, 0.0, 0.2 * a, 1.0e4, 2.0e5, 3.0e4};
  double y[6];
  integrate(&c, y0, 2000.0, y);
  const double v0 = norm(y0 + 3);
  const double v1 = norm(y + 3);
  /* B = (3 (mu.rhat) rhat - mu) / r^3 with mu = g10 a^3 z: A = mu x r / r^3, A_phi = g10 a^3 rho / r^3 */
  double pphi[2];
  const double* s[2] = {y0, y};
  for (int k = 0; k < 2; ++k) {
    const double* x = s[k];
    const double rho2 = x[0] * x[0] + x[1] * x[1];
    const double r = norm(x);
    const double aphi_rho = g[1] * a * a * a * rho2 / (r * r * r);
    pphi[k] = (x[0] * x[4] - x[1] * x[3]) + c.qm * aphi_rho;
  }
  printf("dipole, 2000 s: |v| change %.2e, p_phi change %.2e (relative)\n",
         fabs(v1 - v0) / v0, fabs(pphi[1] - pphi[0]) / fabs(pphi[0]));
  MT_CHECK(fabs(v1 - v0) / v0 < 1e-12);
  MT_CHECK(fabs(pphi[1] - pphi[0]) / fabs(pphi[0]) < 1e-11);
  myb_gravity_work_destroy(work);
  myb_magnetic_field_destroy(field);
}

static void test_fields(void) {
  /* co-rotation: E = -(w x r) x B is perpendicular to B, and the E x B drift is w x r */
  const double omega[3] = {0.0, 0.0, 1.76e-4};
  const double r[3] = {2.1e8, -0.7e8, 0.3e8};
  const double b[3] = {1e-7, 3e-8, -4e-7};
  double e[3];
  myb_corotation_field(omega, r, b, e);
  const double edb = e[0] * b[0] + e[1] * b[1] + e[2] * b[2];
  const double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
  const double drift[3] = {(e[1] * b[2] - e[2] * b[1]) / b2, (e[2] * b[0] - e[0] * b[2]) / b2,
                           (e[0] * b[1] - e[1] * b[0]) / b2};
  const double wxr[3] = {omega[1] * r[2] - omega[2] * r[1], omega[2] * r[0] - omega[0] * r[2],
                         omega[0] * r[1] - omega[1] * r[0]};
  /* the drift is the component of w x r perpendicular to B */
  const double wb = (wxr[0] * b[0] + wxr[1] * b[1] + wxr[2] * b[2]) / b2;
  double worst = fabs(edb) / (norm(e) * norm(b));
  for (int i = 0; i < 3; ++i) {
    worst = fmax(worst, fabs(drift[i] - (wxr[i] - wb * b[i])) / norm(wxr));
  }
  MT_CHECK(worst < 1e-15);
  /* plasma stream: 5 cm^-3 protons at 400 km/s at 1 au, at rest object: rho v^2 cd A/m / 2 */
  const double rho = 5e6 * 1.67262192e-27;
  const double au = 149597870700.0;
  const double x[3] = {au, 0.0, 0.0};
  const double v0[3] = {0.0, 0.0, 0.0};
  double acc[3];
  myb_force_plasma_stream(rho, au, 4e5, 2.0 * 0.01, x, v0, acc);
  const double expected = 0.5 * 2.0 * 0.01 * rho * 4e5 * 4e5;
  printf("solar wind at 1 au, A/m = 0.01 m^2/kg: %.3e m/s^2 (radiation pressure ~4.6e-8)\n",
         acc[0]);
  MT_CHECK(fabs(acc[0] - expected) < 1e-15 * expected && acc[1] == 0.0 && acc[2] == 0.0);
  const double x2[3] = {2.0 * au, 0.0, 0.0};
  double acc2[3];
  myb_force_plasma_stream(rho, au, 4e5, 2.0 * 0.01, x2, v0, acc2);
  MT_CHECK(fabs(acc2[0] - 0.25 * expected) < 1e-15 * expected);
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0);
  test_igrf();
  test_dipole();
  test_uniform();
  test_dipole_invariants();
  test_fields();
  return mt_finish("test_charge");
}
