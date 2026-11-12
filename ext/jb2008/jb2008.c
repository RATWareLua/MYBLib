/* SPDX-License-Identifier: BSD-3-Clause */
/* The Jacchia-Bowman 2008 model of the thermosphere (myb_jb2008.h): written from the papers of
 * Bowman et al. (2008), Bowman et al. (JB2006, 2008) and Jacchia (1970, 1971); the choices they
 * leave open were settled against the outputs of the program of SET run as a black box. */
#include <math.h> /* fabs, floor, sqrt */
#include <stdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "jb2008_internal.h"
#include "myb_jb2008.h"
#include "platform/platform.h"

enum {
  GAUSS = 8,       /* nodes of the Gauss-Legendre rule on every panel */
  MAX_PANELS = 512 /* of one integral: more than the heights the model is used at need */
};

static const double pi = 3.14159265358979323846;
static const double degree = 0.017453292519943295;

/* ---- Constants of the static model (Jacchia 1970) ------------------------------------------- */

static const double gas_constant = 8.31432;  /* J K^-1 mol^-1 */
static const double avogadro = 6.02257e23;   /* mol^-1 */
static const double g0 = 9.80665;            /* m s^-2 */
static const double earth_radius = 6356.766; /* km, of the gravity law */
static const double t0 = 183.0;              /* K at 90 km */
static const double z0 = 90.0;
static const double zx = 125.0;              /* the inflection of the profile */
static const double z_homopause = 105.0;
static const double rho0 = 3.46e-9;          /* g cm^-3 at 90 km */
static const double mean_mass0 = 28.960;     /* at sea level */

/* N2, O2, O, Ar, He, H */
enum { N2, O2, O1, AR, HE, H1, SPECIES };
static const double mass[SPECIES] = {28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797};
static const double thermal_diffusion[SPECIES] = {0.0, 0.0, 0.0, 0.0, -0.38, 0.0};
static const double sea_level[SPECIES] = {0.78110, 0.20955, 0.0, 0.00934, 0.00001289, 0.0};

/* the mean molecular mass between 90 and 105 km: a polynomial in z - 100 */
static const double mean_mass_poly[7] = {28.15204, -0.085586, 1.2840e-4, -1.0056e-5,
                                         -1.0210e-5, 1.5044e-6, 9.9826e-8};

/* Gauss-Legendre, 8 nodes on [-1, 1] (the halves) */
static const double gl_node[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                  0.9602898564975363};
static const double gl_weight[4] = {0.3626837833783620, 0.3137066458778873,
                                    0.2223810344533745, 0.1012285362903763};

/* ---- The corrections of JB2006 (Table 4) ------------------------------------------------------ */

static const double b_coefficient[19] = {
    -0.457512297e+01, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03, 0.703316291e+03,
    -0.194349234e+04, 0.110651308e+04,  -0.174378996e+03, 0.188594601e+04, -0.709371517e+04,
    0.922454523e+04,  -0.384508073e+04, -0.645841789e+01, 0.409703319e+02, -0.482006560e+03,
    0.181870931e+04,  -0.237389204e+04, 0.996703815e+03,  0.361416936e+02};
static const double c_coefficient[23] = {
    -0.155986211e+02, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03,  0.703316291e+03,
    -0.194349234e+04, 0.110651308e+04,  -0.220835117e+03, 0.143256989e+04,  -0.318481844e+04,
    0.328981513e+04,  -0.135332119e+04, 0.199956489e+02,  -0.127093998e+02, 0.212825156e+02,
    -0.275555432e+01, 0.110234982e+02,  0.148881951e+03,  -0.751640284e+03, 0.637876542e+03,
    0.127093998e+02,  -0.212825156e+02, 0.275555432e+01};

static int finite(double x) {
  return x - x == 0.0;
}

/* x^p for x >= 0 */
static double power(double x, double p) {
  return x > 0.0 ? myb_exp(p * myb_log(x)) : 0.0;
}

/* JB2006 eq. 7 (300-600 km): H = z / 100 */
static double band_upper(double f, double t, double c, double z) {
  const double* b = b_coefficient;
  const double h = z / 100.0;
  const double t2 = t * t;
  const double t3 = t2 * t;
  const double t4 = t3 * t;
  const double t5 = t4 * t;
  return b[0] + f * (b[1] + b[2] * t + b[3] * t2 + b[4] * t3 + b[5] * t4 + b[6] * t5) +
         c * (b[7] * t + b[8] * t2 + b[9] * t3 + b[10] * t4 + b[11] * t5) +
         c * h * (b[12] + b[13] * t + b[14] * t2 + b[15] * t3 + b[16] * t4 + b[17] * t5) +
         b[18] * c;
}

/* JB2006 eq. 8 (200-240 km): H = (z - 200) / 50 */
static double band_lower(double f, double t, double c, double z) {
  const double* k = c_coefficient;
  const double h = (z - 200.0) / 50.0;
  const double t2 = t * t;
  const double t3 = t2 * t;
  const double t4 = t3 * t;
  const double t5 = t4 * t;
  return h * k[0] + h * f * (k[1] + k[2] * t + k[3] * t2 + k[4] * t3 + k[5] * t4 + k[6] * t5) +
         h * c *
             (k[7] * t + k[8] * t2 + k[9] * t3 + k[10] * t4 + k[11] * t5 + k[12] + k[13] * f +
              k[14] * f * t + k[15] * f * t2) +
         k[16] + c * (k[17] * t + k[18] * t2 + k[19] * t3 + k[20] * f + k[21] * f * t +
                      k[22] * f * t2);
}

/* The correction of the exospheric temperature at the height z (km) for f = (F10 - 100) / 100,
 * t the local solar time in days and c the cosine of the latitude: the equations in their
 * bands, joined by cubics of Hermite (the slopes per unit of their own H, as the program has
 * them). */
static double temperature_correction(double f, double t, double c, double z) {
  if (z < 120.0 || z >= 800.0) {
    return 0.0;
  }
  const double lower_slope = band_lower(f, t, c, 250.0) - band_lower(f, t, c, 200.0);
  const double upper_slope = band_upper(f, t, c, 400.0) - band_upper(f, t, c, 300.0);
  if (z < 200.0) {
    const double s = (z - 120.0) / 80.0;
    return (3.0 - 2.0 * s) * s * s * band_lower(f, t, c, 200.0) + (s - 1.0) * s * s * lower_slope;
  }
  if (z <= 240.0) {
    return band_lower(f, t, c, z);
  }
  if (z < 300.0) {
    const double s = (z - 240.0) / 60.0;
    const double s2 = s * s;
    const double s3 = s2 * s;
    return (2.0 * s3 - 3.0 * s2 + 1.0) * band_lower(f, t, c, 240.0) +
           (s3 - 2.0 * s2 + s) * lower_slope + (3.0 * s2 - 2.0 * s3) * band_upper(f, t, c, 300.0) +
           (s3 - s2) * upper_slope;
  }
  if (z <= 600.0) {
    return band_upper(f, t, c, z);
  }
  const double s = (z - 600.0) / 200.0;
  const double s2 = s * s;
  const double s3 = s2 * s;
  return (2.0 * s3 - 3.0 * s2 + 1.0) * band_upper(f, t, c, 600.0) +
         (s3 - 2.0 * s2 + s) * 2.0 * upper_slope;
}

/* ---- The temperature profile (Jacchia 1970, eq. 9-13) ---------------------------------------- */

typedef struct profile {
  double tinf; /* the exospheric temperature */
  double tx;   /* at 125 km */
  double gx;   /* the gradient there, K/km */
  double a;    /* 2 (tinf - tx) / pi */
} profile;

static profile make_profile(double tinf) {
  profile p;
  p.tinf = tinf;
  p.tx = 444.3807 + 0.02385 * tinf - 392.8292 * myb_exp(-0.0021357 * tinf);
  p.gx = 1.9 * (p.tx - t0) / (zx - z0);
  p.a = 2.0 / pi * (tinf - p.tx);
  return p;
}

static double temperature_at(const profile* p, double z) {
  if (z <= zx) {
    const double u = (z - zx) / (zx - z0);
    return p->tx + (p->tx - t0) * u * (1.9 - u * u * (1.7 + 0.8 * u));
  }
  const double h = z - zx;
  return p->tx + p->a * myb_atan(p->gx / p->a * h * (1.0 + 4.5e-6 * power(h, 2.5)));
}

/* The acceleration of gravity, m s^-2 */
static double gravity(double z) {
  const double q = 1.0 + z / earth_radius;
  return g0 / (q * q);
}

static double mean_mass(double z) {
  const double x = z - 100.0;
  double sum = mean_mass_poly[6];
  for (int k = 5; k >= 0; --k) {
    sum = sum * x + mean_mass_poly[k];
  }
  return sum;
}

/* The integral over [a, b] (km, a <= b) of g / (R T), times the mean molecular mass (in kg/mol)
 * where `with_mass`: per km of height, in kg/mol^-1 units such that a mass in g/mol times it is
 * the exponent of the barometric or diffusion law. The panels: narrow near the bottom where
 * the temperature turns, wider above; their ends at 105 and 125 km. */
static double integral(const profile* p, double a, double b, int with_mass) {
  if (!(b > a)) {
    return 0.0;
  }
  double sum = 0.0;
  double lo = a;
  for (int panel = 0; panel < MAX_PANELS && lo < b; ++panel) { /* bound: MAX_PANELS */
    double width = lo < zx ? 5.0 : (lo < 200.0 ? 15.0 : (lo < 600.0 ? 40.0 : 100.0));
    double hi = lo + width;
    if (lo < z_homopause && hi > z_homopause) {
      hi = z_homopause;
    } else if (lo < zx && hi > zx) {
      hi = zx;
    }
    if (hi > b || panel == MAX_PANELS - 1) {
      hi = b;
    }
    const double centre = 0.5 * (lo + hi);
    const double half = 0.5 * (hi - lo);
    double part = 0.0;
    for (int k = 0; k < GAUSS; ++k) {
      const double z = centre + (k < 4 ? -1.0 : 1.0) * half * gl_node[k % 4];
      double f = gravity(z) / (gas_constant * temperature_at(p, z));
      if (with_mass) {
        f *= mean_mass(z) * 1e-3;
      }
      part += gl_weight[k % 4] * f;
    }
    sum += half * part;
    lo = hi;
  }
  return sum;
}

/* The static density (g cm^-3) at the height z >= 90 km for the exospheric temperature tinf;
 * number (may be NULL) receives the number densities (cm^-3) of N2, O2, O, Ar, He, H above the
 * homopause (zeros below it). */
static double static_density(double tinf, double z, double* number) {
  if (number != NULL) {
    for (int s = 0; s < SPECIES; ++s) {
      number[s] = 0.0;
    }
  }
  const profile p = make_profile(tinf);
  const double zb = z < z_homopause ? z : z_homopause;
  /* barometric from 90 km: the integral is of M g / (R T) per km times 1000 m */
  const double rho_b = rho0 * (mean_mass(zb) / mean_mass(z0)) *
                       (t0 / temperature_at(&p, zb)) *
                       myb_exp(-1000.0 * integral(&p, z0, zb, 1));
  if (z <= z_homopause) {
    return rho_b;
  }
  /* the number densities at the homopause (cm^-3) */
  const double mh = mean_mass(z_homopause);
  const double total = avogadro * rho_b / mh;
  double n[SPECIES] = {0.0};
  n[N2] = sea_level[N2] * mh / mean_mass0 * total;
  n[AR] = sea_level[AR] * mh / mean_mass0 * total;
  n[HE] = sea_level[HE] * mh / mean_mass0 * total;
  n[O1] = 2.0 * total * (1.0 - mh / mean_mass0);
  n[O2] = total * (mh / mean_mass0 * (1.0 + sea_level[O2]) - 1.0);
  /* diffusion above it: per g/mol of mass the exponent is 1000 m/km * 1e-3 kg/g = 1 */
  const double base = integral(&p, z_homopause, z, 0);
  const double ratio = temperature_at(&p, z_homopause) / temperature_at(&p, z);
  double rho = 0.0;
  for (int s = 0; s < H1; ++s) {
    const double here = n[s] * power(ratio, 1.0 + thermal_diffusion[s]) * myb_exp(-mass[s] * base);
    rho += here * mass[s];
    if (number != NULL) {
      number[s] = here;
    }
  }
  /* hydrogen from 500 km, also below it (Jacchia 1971) */
  const double t500 = temperature_at(&p, 500.0);
  const double lt = myb_log(t500) / myb_log(10.0);
  const double n500 = myb_exp((73.13 - 39.40 * lt + 5.5 * lt * lt) * myb_log(10.0));
  const double from500 = z >= 500.0 ? integral(&p, 500.0, z, 0) : -integral(&p, z, 500.0, 0);
  const double nh = n500 * (t500 / temperature_at(&p, z)) * myb_exp(-mass[H1] * from500);
  rho += nh * mass[H1];
  if (number != NULL) {
    number[H1] = nh;
  }
  return rho / avogadro;
}

double myb_internal_jb2008_static(double tinf, double z, double* temperature, double number[6]) {
  const profile p = make_profile(tinf);
  *temperature = temperature_at(&p, z);
  return static_density(tinf, z, number);
}

/* ---- The time ------------------------------------------------------------------------------- */

/* The MJD of day `doy` (1 on January 1) of the Gregorian `year`. */
static int64_t mjd_of(int64_t year, int64_t doy) {
  const int64_t yp = year - 1; /* January counts in the March-based year before */
  const int64_t era = (yp >= 0 ? yp : yp - 399) / 400;
  const int64_t yoe = yp - era * 400;
  const int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + 306;
  return era * 146097 + doe - 678881 + doy - 1;
}

/* The year and the day of the year (1 on January 1) of the day with the modified Julian date
 * mjd. */
static void year_and_day(int64_t mjd, int64_t* year, int64_t* day) {
  const int64_t z = mjd + 678881; /* days from 0000-03-01 */
  const int64_t era = (z >= 0 ? z : z - 146096) / 146097;
  const int64_t doe = z - era * 146097;
  const int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const int64_t doy_march = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const int64_t mp = (5 * doy_march + 2) / 153;
  const int64_t y = yoe + era * 400 + (mp >= 10 ? 1 : 0); /* January and February */
  *year = y;
  *day = mjd - mjd_of(y, 1) + 1;
}

/* ---- The model ------------------------------------------------------------------------------ */

myb_status myb_jb2008(double mjd, double sun_right_ascension, double sun_declination,
                      double right_ascension, double latitude, double height,
                      const myb_jb2008_indices* indices, myb_jb2008_output* output) {
  MYB_CHECK(indices != NULL && output != NULL);
  const myb_jb2008_indices* x = indices;
  const double values[9] = {x->f10, x->f10_average, x->s10, x->s10_average, x->m10,
                            x->m10_average, x->y10, x->y10_average, x->dtc};
  for (int k = 0; k < 9; ++k) {
    if (!finite(values[k]) || (k < 8 && !(values[k] >= 0.0))) {
      return MYB_E_DOMAIN;
    }
  }
  if (!finite(mjd) || !finite(sun_right_ascension) || !finite(sun_declination) ||
      !finite(right_ascension) || !finite(latitude) || !finite(height) || !(height >= 90.0) ||
      fabs(mjd) > 1e8) {
    return MYB_E_DOMAIN;
  }
  /* the nighttime minimum of the global exospheric temperature (JB2008 eq. 1-2) */
  const double weight = power(x->f10_average / 240.0, 0.25);
  const double mixed = x->f10_average * weight + x->s10_average * (1.0 - weight);
  const double tc = 392.4 + 3.227 * mixed + 0.298 * (x->f10 - x->f10_average) +
                    2.259 * (x->s10 - x->s10_average) + 0.312 * (x->m10 - x->m10_average) +
                    0.178 * (x->y10 - x->y10_average);
  /* the diurnal variation (Jacchia 1970 eq. 15-17: m = 2.5, n = 3, R = 0.31) */
  const double hour_angle = right_ascension - sun_right_ascension;
  const double eta = 0.5 * fabs(latitude - sun_declination);
  const double theta = 0.5 * fabs(latitude + sun_declination);
  double tau = hour_angle + (-37.0 + 6.0 * myb_sin(hour_angle + 43.0 * degree)) * degree;
  tau -= 2.0 * pi * floor((tau + pi) / (2.0 * pi));
  const double c = power(myb_cos(eta), 2.5);
  const double s = power(myb_sin(theta), 2.5);
  const double night = tc * (1.0 + 0.31 * s);
  const double amplitude = 0.31 * (c - s) / (1.0 + 0.31 * s);
  const double day = power(fabs(myb_cos(0.5 * tau)), 3.0);
  const double tl = night * (1.0 + amplitude * day);
  const double tinf = tl + x->dtc;
  /* the correction of the temperature at the height (JB2006) */
  double solar_time = 12.0 + hour_angle / degree / 15.0;
  solar_time -= 24.0 * floor(solar_time / 24.0);
  const double tinf_here = tinf + temperature_correction((x->f10 - 100.0) / 100.0,
                                                         solar_time / 24.0, myb_cos(latitude),
                                                         height);
  const profile p = make_profile(tinf_here);
  double log_rho = myb_log(static_density(tinf_here, height, NULL)) / myb_log(10.0);
  /* the seasonal-latitudinal variation of the lower thermosphere: the amplitude of Jacchia
   * 1970 eq. 24, the time of Jacchia 1971 eq. 23 */
  if (latitude != 0.0) {
    const double dz = height - 90.0;
    const double phase = 2.0 * pi * (mjd - 36204.0) / 365.2422 + 1.72;
    const double sl = myb_sin(latitude);
    log_rho += 0.02 * dz * myb_exp(-0.045 * dz) * (latitude > 0.0 ? 1.0 : -1.0) * sl * sl *
               myb_sin(phase);
  }
  /* the semiannual variation (JB2008 eq. 3-7), F(z) not below zero */
  {
    int64_t year = 0;
    int64_t day_number = 0;
    const double whole = floor(mjd);
    year_and_day((int64_t)whole, &year, &day_number);
    (void)year;
    const double t = (double)day_number + (mjd - whole);
    const double w = 2.0 * pi * (t - 1.0) / 365.0;
    const double fsmj = x->f10_average - 0.70 * x->s10_average - 0.04 * x->m10_average;
    const double zz = height / 1000.0;
    double fz = 2.69e-01 - 1.18e-02 * fsmj + 2.78e-02 * zz * fsmj - 2.78e-02 * zz * zz * fsmj +
                3.47e-04 * zz * fsmj * fsmj;
    fz = fz > 0.0 ? fz : 0.0;
    const double fsm = x->f10_average - 0.75 * x->s10_average - 0.37 * x->m10_average;
    const double sw = myb_sin(w);
    const double cw = myb_cos(w);
    const double s2w = myb_sin(2.0 * w);
    const double c2w = myb_cos(2.0 * w);
    const double gt = -3.63e-01 + 8.51e-02 * sw + 2.40e-01 * cw - 1.90e-01 * s2w - 2.55e-01 * c2w +
                      fsm * (-1.79e-02 + 5.65e-04 * sw - 6.41e-04 * cw - 3.42e-03 * s2w -
                             1.25e-03 * c2w);
    log_rho += fz * gt;
  }
  double rho = myb_exp(log_rho * myb_log(10.0));
  /* above 1000 km (JB2006 eq. 10-11) */
  if (height > 1000.0) {
    const double c1 = 2.200e-01;
    const double c2 = -2.000e-03;
    const double c3 = 1.150e-03;
    const double c4 = -2.110e-06;
    const double fbar = x->f10_average;
    double factor = 0.0;
    if (height >= 1500.0) {
      factor = c1 + c2 * fbar + c3 * height + c4 * height * fbar;
    } else {
      const double h = (height - 1000.0) / 500.0;
      const double f1500 = c1 + c2 * fbar + c3 * 1500.0 + c4 * 1500.0 * fbar;
      const double slope = 500.0 * (c3 + c4 * fbar);
      factor = 1.0 + (3.0 * f1500 - slope - 3.0) * h * h + (slope - 2.0 * f1500 + 2.0) * h * h * h;
    }
    rho *= factor;
  }
  output->density = rho * 1000.0; /* g cm^-3 -> kg m^-3 */
  output->exospheric_temperature = tinf;
  output->temperature = temperature_at(&p, height);
  return finite(output->density) ? MYB_OK : MYB_E_DOMAIN;
}

/* ---- The files of indices ------------------------------------------------------------------- */

struct myb_jb2008_solar {
  myb_allocator allocator;
  int64_t first; /* MJD of the first day */
  size_t days;
  double* values; /* [day][8]: F10 F81c S10 S81c M10 M81c Y10 Y81c */
};

struct myb_jb2008_storm {
  myb_allocator allocator;
  int64_t first;
  size_t days;
  double* values; /* [day][24] */
};

/* The next line of the text; 0 at the end. */
static int next_line(const char* text, size_t size, size_t* at, const char** line,
                     size_t* length) {
  if (*at >= size) {
    return 0;
  }
  size_t end = *at;
  while (end < size && text[end] != '\n') { /* bound: the characters of the text */
    ++end;
  }
  *line = text + *at;
  *length = end - *at;
  *at = end + 1;
  while (*length > 0 && ((*line)[*length - 1] == '\r' || (*line)[*length - 1] == ' ' ||
                         (*line)[*length - 1] == '\t')) {
    --*length;
  }
  return 1;
}

/* Up to `capacity` numbers separated by blanks at the start of a line (the first that is not a
 * number ends them): their count. */
static int numbers(const char* line, size_t length, double* out, int capacity) {
  int n = 0;
  size_t at = 0;
  while (at < length && n < capacity) { /* bound: the characters of the line */
    if (line[at] == ' ' || line[at] == '\t') {
      ++at;
      continue;
    }
    size_t end = at;
    while (end < length && line[end] != ' ' && line[end] != '\t') {
      ++end;
    }
    double value = 0.0;
    size_t used = 0;
    if (myb_parse_double(line + at, end - at, &value, &used) != MYB_OK || used != end - at ||
        !finite(value)) {
      break;
    }
    out[n++] = value;
    at = end;
  }
  return n;
}

/* Reads records of `width` values after a date (year, day of year) — `skip` leading fields
 * between the date and the values — into a table of whole days. */
static myb_status read_days(const char* text, size_t size, const myb_allocator* allocator,
                            int prefix, int skip, int width, int64_t* first, size_t* days,
                            double** values, myb_diag* diag) {
  const myb_allocator a = myb_internal_allocator(allocator);
  /* count the records */
  size_t at = 0;
  const char* line = NULL;
  size_t length = 0;
  size_t records = 0;
  while (next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    size_t k = 0;
    while (k < length && (line[k] == ' ' || line[k] == '\t')) {
      ++k;
    }
    if (k < length && line[k] != '#') {
      ++records;
    }
  }
  if (records == 0) {
    myb_internal_diag(diag, 0, "no records");
    return MYB_E_FORMAT;
  }
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, records * (size_t)width, sizeof(double), &block);
  if (status != MYB_OK) {
    return status;
  }
  double* out = block;
  at = 0;
  long number = 0;
  size_t day = 0;
  int64_t start = 0;
  while (next_line(text, size, &at, &line, &length)) { /* bound: the lines */
    ++number;
    size_t k = 0;
    while (k < length && (line[k] == ' ' || line[k] == '\t')) {
      ++k;
    }
    if (k == length || line[k] == '#') {
      continue;
    }
    if (prefix) { /* "DTC" */
      if (length - k < 3 || line[k] != 'D' || line[k + 1] != 'T' || line[k + 2] != 'C') {
        myb_internal_diag(diag, number, "a record that does not begin with DTC");
        myb_internal_release(&a, out);
        return MYB_E_FORMAT;
      }
      k += 3;
    }
    double field[2 + 2 + 24];
    const int want = 2 + skip + width;
    if (numbers(line + k, length - k, field, want) != want || field[0] != floor(field[0]) ||
        field[1] != floor(field[1]) || !(field[1] >= 1.0 && field[1] <= 366.0) ||
        fabs(field[0]) > 1e6) {
      myb_internal_diag(diag, number, "expected the year, the day and %d values", width);
      myb_internal_release(&a, out);
      return MYB_E_FORMAT;
    }
    const int64_t mjd = mjd_of((int64_t)field[0], (int64_t)field[1]);
    if (day == 0) {
      start = mjd;
    } else if (mjd != start + (int64_t)day) {
      myb_internal_diag(diag, number, "a day out of sequence");
      myb_internal_release(&a, out);
      return MYB_E_FORMAT;
    }
    memcpy(out + day * (size_t)width, field + 2 + skip, (size_t)width * sizeof(double));
    ++day;
  }
  *first = start;
  *days = day;
  *values = out;
  return MYB_OK;
}

myb_status myb_jb2008_solar_load(const char* text, size_t size, const myb_allocator* allocator,
                                 myb_jb2008_solar** table, myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && table != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_jb2008_solar), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_jb2008_solar* t = block;
  memset(t, 0, sizeof *t);
  t->allocator = a;
  status = read_days(text, size, allocator, 0, 1, 8, &t->first, &t->days, &t->values, diag);
  if (status != MYB_OK) {
    myb_internal_release(&a, t);
    return status;
  }
  *table = t;
  return MYB_OK;
}

void myb_jb2008_solar_destroy(myb_jb2008_solar* table) {
  if (table == NULL) {
    return;
  }
  const myb_allocator a = table->allocator;
  myb_internal_release(&a, table->values);
  myb_internal_release(&a, table);
}

static const double* solar_day(const myb_jb2008_solar* table, int64_t mjd) {
  if (mjd < table->first || mjd >= table->first + (int64_t)table->days) {
    return NULL;
  }
  return table->values + (size_t)(mjd - table->first) * 8;
}

myb_status myb_jb2008_solar_at(const myb_jb2008_solar* table, double mjd,
                               myb_jb2008_indices* indices) {
  MYB_CHECK(table != NULL && indices != NULL);
  if (!finite(mjd) || fabs(mjd) > 1e8) {
    return MYB_E_DOMAIN;
  }
  const int64_t today = (int64_t)floor(mjd);
  const double* one = solar_day(table, today - 1);
  const double* two = solar_day(table, today - 2);
  const double* five = solar_day(table, today - 5);
  if (one == NULL || two == NULL || five == NULL) {
    return MYB_E_RANGE;
  }
  indices->f10 = one[0];
  indices->f10_average = one[1];
  indices->s10 = one[2];
  indices->s10_average = one[3];
  indices->m10 = two[4];
  indices->m10_average = two[5];
  indices->y10 = five[6];
  indices->y10_average = five[7];
  return MYB_OK;
}

myb_status myb_jb2008_storm_load(const char* text, size_t size, const myb_allocator* allocator,
                                 myb_jb2008_storm** table, myb_diag* diag) {
  MYB_CHECK((text != NULL || size == 0) && table != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_jb2008_storm), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_jb2008_storm* t = block;
  memset(t, 0, sizeof *t);
  t->allocator = a;
  status = read_days(text, size, allocator, 1, 0, 24, &t->first, &t->days, &t->values, diag);
  if (status != MYB_OK) {
    myb_internal_release(&a, t);
    return status;
  }
  *table = t;
  return MYB_OK;
}

void myb_jb2008_storm_destroy(myb_jb2008_storm* table) {
  if (table == NULL) {
    return;
  }
  const myb_allocator a = table->allocator;
  myb_internal_release(&a, table->values);
  myb_internal_release(&a, table);
}

myb_status myb_jb2008_storm_at(const myb_jb2008_storm* table, double mjd, double* dtc) {
  MYB_CHECK(table != NULL && dtc != NULL);
  if (!finite(mjd) || fabs(mjd) > 1e8) {
    return MYB_E_DOMAIN;
  }
  const double whole = floor(mjd);
  const int64_t today = (int64_t)whole;
  if (today < table->first || today >= table->first + (int64_t)table->days) {
    return MYB_E_RANGE;
  }
  int hour = (int)floor((mjd - whole) * 24.0);
  hour = hour < 0 ? 0 : (hour > 23 ? 23 : hour);
  *dtc = table->values[(size_t)(today - table->first) * 24 + (size_t)hour];
  return MYB_OK;
}

myb_status myb_jb2008_solar_load_file(const char* path, const myb_allocator* allocator,
                                      myb_jb2008_solar** table, myb_diag* diag) {
  MYB_CHECK(path != NULL && table != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_jb2008_solar_load((const char*)map.data, map.size, allocator, table, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

myb_status myb_jb2008_storm_load_file(const char* path, const myb_allocator* allocator,
                                      myb_jb2008_storm** table, myb_diag* diag) {
  MYB_CHECK(path != NULL && table != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  myb_platform_mapping map = {0};
  const myb_status status = myb_platform_map_file(path, &a, &map);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot map '%s'", path);
    return status;
  }
  const myb_status loaded =
      myb_jb2008_storm_load((const char*)map.data, map.size, allocator, table, diag);
  myb_platform_unmap_file(&map);
  return loaded;
}

/* ---- The atmosphere of a force model ---------------------------------------------------------- */

myb_status myb_jb2008_density(void* atmosphere, myb_time t, const double r[3], double* density) {
  const myb_jb2008_atmosphere* a = atmosphere;
  MYB_CHECK(a != NULL && a->timesys != NULL && r != NULL && density != NULL);
  myb_calendar utc;
  myb_status status = myb_timesys_to_utc(a->timesys, t, a->scale, &utc);
  if (status != MYB_OK) {
    return status;
  }
  if (utc.year < -100000 || utc.year > 100000) {
    return MYB_E_DOMAIN;
  }
  /* the MJD of the UTC day and its fraction */
  int64_t day_number = 0;
  {
    static const int before[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    const int64_t y = utc.year;
    const int leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    day_number = before[utc.month - 1] + utc.day + (leap && utc.month > 2 ? 1 : 0);
  }
  const double seconds = 3600.0 * utc.hour + 60.0 * utc.minute + utc.second;
  const double day_fraction = (seconds < 86400.0 ? seconds : 86399.999) / 86400.0;
  const double mjd = (double)mjd_of(utc.year, day_number) + day_fraction;
  myb_jb2008_indices indices = a->fixed;
  if (a->solar != NULL) {
    status = myb_jb2008_solar_at(a->solar, mjd, &indices);
    if (status != MYB_OK) {
      return status;
    }
  }
  if (a->storm != NULL) {
    status = myb_jb2008_storm_at(a->storm, mjd, &indices.dtc);
    if (status != MYB_OK) {
      return status;
    }
  }
  double longitude = 0.0;
  double latitude = 0.0;
  double height = 0.0;
  status = myb_geodetic_from_cartesian(&a->ellipsoid, r, &longitude, &latitude, &height);
  if (status != MYB_OK) {
    return status;
  }
  const double geocentric = myb_atan2(r[2], sqrt(r[0] * r[0] + r[1] * r[1]));
  /* the low-precision Sun of the Astronomical Almanac and the mean sidereal time of 1982 */
  const double n = mjd - 51544.5;
  const double mean_longitude = (280.460 + 0.9856474 * n) * degree;
  const double anomaly = (357.528 + 0.9856003 * n) * degree;
  const double ecliptic = mean_longitude + (1.915 * myb_sin(anomaly) +
                                            0.020 * myb_sin(2.0 * anomaly)) * degree;
  const double obliquity = (23.439 - 0.0000004 * n) * degree;
  const double sun_ra = myb_atan2(myb_cos(obliquity) * myb_sin(ecliptic), myb_cos(ecliptic));
  const double sun_dec = myb_asin(myb_sin(obliquity) * myb_sin(ecliptic));
  const double centuries = n / 36525.0;
  double gmst = 280.46061837 + 360.98564736629 * n + 0.000387933 * centuries * centuries -
                centuries * centuries * centuries / 38710000.0;
  gmst -= 360.0 * floor(gmst / 360.0);
  myb_jb2008_output out;
  status = myb_jb2008(mjd, sun_ra, sun_dec, gmst * degree + longitude, geocentric,
                      height / 1000.0, &indices, &out);
  if (status == MYB_OK) {
    *density = out.density;
  }
  return status;
}
