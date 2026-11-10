/* SPDX-License-Identifier: MIT */
/* SGP4/SDP4 — a C port of the Vallado revision of Spacetrack Report #3 (see NOTICE): the
 * routines dpper, dscom, dsinit, dspace, initl, sgp4init, sgp4 and gstime keep the order of the
 * operations of the original, so that results agree to rounding. Changes: correctly rounded
 * elementary functions (myb_sin, ...); fmod computed exactly with one fused multiply-add;
 * pow(x, 2/3) as cbrt(x^2) and integer and half-integer powers by products and sqrt; the deep-space
 * resonance integrator restarts from the epoch on every call (its steps of 720 minutes are fixed,
 * so the cached restart point of the original gave the same values); the per-call fields of the
 * original (con41, x1mth2, x7thm1, xlcof, aycof in deep space) are local. */
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "core/alloc.h"
#include "core/check.h"
#include "myb/cmath.h"
#include "myb/io.h"
#include "myb_sgp4.h"

#define SGP4_PI 3.14159265358979323846

static const double PI = SGP4_PI;
static const double TWOPI = 2.0 * SGP4_PI;
static const double DEG2RAD = SGP4_PI / 180.0;
static const double X2O3 = 2.0 / 3.0;

/* fmod(x, y) of C for y > 0 and |x / y| < 2^52: exact, as every fmod. */
static double fmod_exact(double x, double y) {
  double q = trunc(x / y);
  double r = myb_fma(-q, y, x);
  if (x >= 0.0) {
    if (r < 0.0) {
      q -= 1.0;
      r = myb_fma(-q, y, x);
    } else if (r >= y) {
      q += 1.0;
      r = myb_fma(-q, y, x);
    }
  } else {
    if (r > 0.0) {
      q += 1.0;
      r = myb_fma(-q, y, x);
    } else if (r <= -y) {
      q -= 1.0;
      r = myb_fma(-q, y, x);
    }
  }
  return r;
}

static double pow23(double x) { /* x^(2/3), x > 0 */
  return myb_cbrt(x * x);
}

/* ---- Lunar and solar periodics (dpper) -------------------------------------------------------- */

static void dpper(const myb_sgp4* s, double t, double* ep, double* inclp, double* nodep,
                  double* argpp, double* mp) {
  const double zns = 1.19459e-5;
  const double zes = 0.01675;
  const double znl = 1.5835218e-4;
  const double zel = 0.05490;
  double zm = s->zmos + zns * t;
  double zf = zm + 2.0 * zes * myb_sin(zm);
  double sinzf = myb_sin(zf);
  double f2 = 0.5 * sinzf * sinzf - 0.25;
  double f3 = -0.5 * sinzf * myb_cos(zf);
  const double ses = s->se2 * f2 + s->se3 * f3;
  const double sis = s->si2 * f2 + s->si3 * f3;
  const double sls = s->sl2 * f2 + s->sl3 * f3 + s->sl4 * sinzf;
  const double sghs = s->sgh2 * f2 + s->sgh3 * f3 + s->sgh4 * sinzf;
  const double shs = s->sh2 * f2 + s->sh3 * f3;
  zm = s->zmol + znl * t;
  zf = zm + 2.0 * zel * myb_sin(zm);
  sinzf = myb_sin(zf);
  f2 = 0.5 * sinzf * sinzf - 0.25;
  f3 = -0.5 * sinzf * myb_cos(zf);
  const double sel = s->ee2 * f2 + s->e3 * f3;
  const double sil = s->xi2 * f2 + s->xi3 * f3;
  const double sll = s->xl2 * f2 + s->xl3 * f3 + s->xl4 * sinzf;
  const double sghl = s->xgh2 * f2 + s->xgh3 * f3 + s->xgh4 * sinzf;
  const double shll = s->xh2 * f2 + s->xh3 * f3;
  double pe = ses + sel;
  double pinc = sis + sil;
  double pl = sls + sll;
  double pgh = sghs + sghl;
  double ph = shs + shll;
  /* init = 'n' (the call with init = 'y' changes nothing) */
  pe = pe - s->peo;
  pinc = pinc - s->pinco;
  pl = pl - s->plo;
  pgh = pgh - s->pgho;
  ph = ph - s->pho;
  *inclp = *inclp + pinc;
  *ep = *ep + pe;
  const double sinip = myb_sin(*inclp);
  const double cosip = myb_cos(*inclp);
  if (*inclp >= 0.2) {
    ph = ph / sinip;
    pgh = pgh - cosip * ph;
    *argpp = *argpp + pgh;
    *nodep = *nodep + ph;
    *mp = *mp + pl;
  } else {
    /* Lyddane modification */
    const double sinop = myb_sin(*nodep);
    const double cosop = myb_cos(*nodep);
    double alfdp = sinip * sinop;
    double betdp = sinip * cosop;
    const double dalf = ph * cosop + pinc * cosip * sinop;
    const double dbet = -ph * sinop + pinc * cosip * cosop;
    alfdp = alfdp + dalf;
    betdp = betdp + dbet;
    *nodep = fmod_exact(*nodep, TWOPI);
    if (*nodep < 0.0 && s->mode == MYB_SGP4_AFSPC) {
      *nodep = *nodep + TWOPI;
    }
    const double xls = *mp + *argpp + pl + pgh + (cosip - pinc * sinip) * *nodep;
    const double xnoh = *nodep;
    *nodep = myb_atan2(alfdp, betdp);
    if (*nodep < 0.0 && s->mode == MYB_SGP4_AFSPC) {
      *nodep = *nodep + TWOPI;
    }
    if (fabs(xnoh - *nodep) > PI) {
      if (*nodep < xnoh) {
        *nodep = *nodep + TWOPI;
      } else {
        *nodep = *nodep - TWOPI;
      }
    }
    *mp = *mp + pl;
    *argpp = xls - *mp - cosip * *nodep;
  }
}

/* ---- Deep-space common terms (dscom) ---------------------------------------------------------- */

typedef struct dscom_out {
  double snodm, cnodm, sinim, cosim, sinomm, cosomm, day, em, emsq, gam, rtemsq;
  double s1, s2, s3, s4, s5, s6, s7;
  double ss1, ss2, ss3, ss4, ss5, ss6, ss7;
  double sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33;
  double z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33;
  double nm;
} dscom_out;

static void dscom(double epoch, double ep, double argpp, double tc, double inclp, double nodep,
                  double np, myb_sgp4* s, dscom_out* o) {
  const double zes = 0.01675;
  const double zel = 0.05490;
  const double c1ss = 2.9864797e-6;
  const double c1l = 4.7968065e-7;
  const double zsinis = 0.39785416;
  const double zcosis = 0.91744867;
  const double zcosgs = 0.1945905;
  const double zsings = -0.98088458;

  o->nm = np;
  o->em = ep;
  o->snodm = myb_sin(nodep);
  o->cnodm = myb_cos(nodep);
  o->sinomm = myb_sin(argpp);
  o->cosomm = myb_cos(argpp);
  o->sinim = myb_sin(inclp);
  o->cosim = myb_cos(inclp);
  o->emsq = o->em * o->em;
  const double betasq = 1.0 - o->emsq;
  o->rtemsq = sqrt(betasq);

  s->peo = 0.0;
  s->pinco = 0.0;
  s->plo = 0.0;
  s->pgho = 0.0;
  s->pho = 0.0;
  o->day = epoch + 18261.5 + tc / 1440.0;
  const double xnodce = fmod_exact(4.5236020 - 9.2422029e-4 * o->day, TWOPI);
  const double stem = myb_sin(xnodce);
  const double ctem = myb_cos(xnodce);
  const double zcosil = 0.91375164 - 0.03568096 * ctem;
  const double zsinil = sqrt(1.0 - zcosil * zcosil);
  const double zsinhl = 0.089683511 * stem / zsinil;
  const double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
  o->gam = 5.8351514 + 0.0019443680 * o->day;
  double zx = 0.39785416 * stem / zsinil;
  const double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
  zx = myb_atan2(zx, zy);
  zx = o->gam + zx - xnodce;
  const double zcosgl = myb_cos(zx);
  const double zsingl = myb_sin(zx);

  /* solar terms, then lunar */
  double zcosg = zcosgs;
  double zsing = zsings;
  double zcosi = zcosis;
  double zsini = zsinis;
  double zcosh = o->cnodm;
  double zsinh = o->snodm;
  double cc = c1ss;
  const double xnoi = 1.0 / o->nm;

  for (int lsflg = 1; lsflg <= 2; ++lsflg) {
    const double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
    const double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
    const double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
    const double a8 = zsing * zsini;
    const double a9 = zsing * zsinh + zcosg * zcosi * zcosh;
    const double a10 = zcosg * zsini;
    const double a2 = o->cosim * a7 + o->sinim * a8;
    const double a4 = o->cosim * a9 + o->sinim * a10;
    const double a5 = -o->sinim * a7 + o->cosim * a8;
    const double a6 = -o->sinim * a9 + o->cosim * a10;

    const double x1 = a1 * o->cosomm + a2 * o->sinomm;
    const double x2 = a3 * o->cosomm + a4 * o->sinomm;
    const double x3 = -a1 * o->sinomm + a2 * o->cosomm;
    const double x4 = -a3 * o->sinomm + a4 * o->cosomm;
    const double x5 = a5 * o->sinomm;
    const double x6 = a6 * o->sinomm;
    const double x7 = a5 * o->cosomm;
    const double x8 = a6 * o->cosomm;

    o->z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
    o->z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
    o->z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
    o->z1 = 3.0 * (a1 * a1 + a2 * a2) + o->z31 * o->emsq;
    o->z2 = 6.0 * (a1 * a3 + a2 * a4) + o->z32 * o->emsq;
    o->z3 = 3.0 * (a3 * a3 + a4 * a4) + o->z33 * o->emsq;
    o->z11 = -6.0 * a1 * a5 + o->emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
    o->z12 = -6.0 * (a1 * a6 + a3 * a5) +
             o->emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
    o->z13 = -6.0 * a3 * a6 + o->emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
    o->z21 = 6.0 * a2 * a5 + o->emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
    o->z22 = 6.0 * (a4 * a5 + a2 * a6) +
             o->emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
    o->z23 = 6.0 * a4 * a6 + o->emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
    o->z1 = o->z1 + o->z1 + betasq * o->z31;
    o->z2 = o->z2 + o->z2 + betasq * o->z32;
    o->z3 = o->z3 + o->z3 + betasq * o->z33;
    o->s3 = cc * xnoi;
    o->s2 = -0.5 * o->s3 / o->rtemsq;
    o->s4 = o->s3 * o->rtemsq;
    o->s1 = -15.0 * o->em * o->s4;
    o->s5 = x1 * x3 + x2 * x4;
    o->s6 = x2 * x3 + x1 * x4;
    o->s7 = x2 * x4 - x1 * x3;

    if (lsflg == 1) {
      o->ss1 = o->s1;
      o->ss2 = o->s2;
      o->ss3 = o->s3;
      o->ss4 = o->s4;
      o->ss5 = o->s5;
      o->ss6 = o->s6;
      o->ss7 = o->s7;
      o->sz1 = o->z1;
      o->sz2 = o->z2;
      o->sz3 = o->z3;
      o->sz11 = o->z11;
      o->sz12 = o->z12;
      o->sz13 = o->z13;
      o->sz21 = o->z21;
      o->sz22 = o->z22;
      o->sz23 = o->z23;
      o->sz31 = o->z31;
      o->sz32 = o->z32;
      o->sz33 = o->z33;
      zcosg = zcosgl;
      zsing = zsingl;
      zcosi = zcosil;
      zsini = zsinil;
      zcosh = zcoshl * o->cnodm + zsinhl * o->snodm;
      zsinh = o->snodm * zcoshl - o->cnodm * zsinhl;
      cc = c1l;
    }
  }

  s->zmol = fmod_exact(4.7199672 + 0.22997150 * o->day - o->gam, TWOPI);
  s->zmos = fmod_exact(6.2565837 + 0.017201977 * o->day, TWOPI);

  /* solar terms */
  s->se2 = 2.0 * o->ss1 * o->ss6;
  s->se3 = 2.0 * o->ss1 * o->ss7;
  s->si2 = 2.0 * o->ss2 * o->sz12;
  s->si3 = 2.0 * o->ss2 * (o->sz13 - o->sz11);
  s->sl2 = -2.0 * o->ss3 * o->sz2;
  s->sl3 = -2.0 * o->ss3 * (o->sz3 - o->sz1);
  s->sl4 = -2.0 * o->ss3 * (-21.0 - 9.0 * o->emsq) * zes;
  s->sgh2 = 2.0 * o->ss4 * o->sz32;
  s->sgh3 = 2.0 * o->ss4 * (o->sz33 - o->sz31);
  s->sgh4 = -18.0 * o->ss4 * zes;
  s->sh2 = -2.0 * o->ss2 * o->sz22;
  s->sh3 = -2.0 * o->ss2 * (o->sz23 - o->sz21);

  /* lunar terms */
  s->ee2 = 2.0 * o->s1 * o->s6;
  s->e3 = 2.0 * o->s1 * o->s7;
  s->xi2 = 2.0 * o->s2 * o->z12;
  s->xi3 = 2.0 * o->s2 * (o->z13 - o->z11);
  s->xl2 = -2.0 * o->s3 * o->z2;
  s->xl3 = -2.0 * o->s3 * (o->z3 - o->z1);
  s->xl4 = -2.0 * o->s3 * (-21.0 - 9.0 * o->emsq) * zel;
  s->xgh2 = 2.0 * o->s4 * o->z32;
  s->xgh3 = 2.0 * o->s4 * (o->z33 - o->z31);
  s->xgh4 = -18.0 * o->s4 * zel;
  s->xh2 = -2.0 * o->s2 * o->z22;
  s->xh3 = -2.0 * o->s2 * (o->z23 - o->z21);
}

/* ---- Deep-space initialisation (dsinit), at t = tc = 0 ---------------------------------------- */

static void dsinit(myb_sgp4* s, const dscom_out* o, double xpidot, double eccsq) {
  const double q22 = 1.7891679e-6;
  const double q31 = 2.1460748e-6;
  const double q33 = 2.2123015e-7;
  const double root22 = 1.7891679e-6;
  const double root44 = 7.3636953e-9;
  const double root54 = 2.1765803e-9;
  const double rptim = 4.37526908801129966e-3; /* 7.29211514668855e-5 rad/s in rad/min */
  const double root32 = 3.7393792e-7;
  const double root52 = 1.1428639e-7;
  const double znl = 1.5835218e-4;
  const double zns = 1.19459e-5;

  const double cosim = o->cosim;
  const double sinim = o->sinim;
  double em = o->em;
  double emsq = o->emsq;
  const double nm = o->nm;
  const double inclm = s->inclo;

  /* resonance */
  s->irez = 0;
  if ((nm < 0.0052359877) && (nm > 0.0034906585)) {
    s->irez = 1;
  }
  if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5)) {
    s->irez = 2;
  }

  /* solar terms */
  const double ses = o->ss1 * zns * o->ss5;
  const double sis = o->ss2 * zns * (o->sz11 + o->sz13);
  const double sls = -zns * o->ss3 * (o->sz1 + o->sz3 - 14.0 - 6.0 * emsq);
  const double sghs = o->ss4 * zns * (o->sz31 + o->sz33 - 6.0);
  double shs = -zns * o->ss2 * (o->sz21 + o->sz23);
  if ((inclm < 5.2359877e-2) || (inclm > PI - 5.2359877e-2)) {
    shs = 0.0;
  }
  if (sinim != 0.0) {
    shs = shs / sinim;
  }
  const double sgs = sghs - cosim * shs;

  /* lunar terms */
  s->dedt = ses + o->s1 * znl * o->s5;
  s->didt = sis + o->s2 * znl * (o->z11 + o->z13);
  s->dmdt = sls - znl * o->s3 * (o->z1 + o->z3 - 14.0 - 6.0 * emsq);
  const double sghl = o->s4 * znl * (o->z31 + o->z33 - 6.0);
  double shll = -znl * o->s2 * (o->z21 + o->z23);
  if ((inclm < 5.2359877e-2) || (inclm > PI - 5.2359877e-2)) {
    shll = 0.0;
  }
  s->domdt = sgs + sghl;
  s->dnodt = shs;
  if (sinim != 0.0) {
    s->domdt = s->domdt - cosim / sinim * shll;
    s->dnodt = s->dnodt + shll / sinim;
  }

  /* resonance effects */
  const double theta = fmod_exact(s->gsto, TWOPI); /* gsto + tc * rptim, tc = 0 */
  if (s->irez == 0) {
    return;
  }
  const double aonv = pow23(nm / s->xke);
  if (s->irez == 2) {
    /* geopotential resonance for 12-hour orbits */
    const double cosisq = cosim * cosim;
    const double emo = em;
    em = s->ecco;
    const double emsqo = emsq;
    emsq = eccsq;
    const double eoc = em * emsq;
    const double g201 = -0.306 - (em - 0.64) * 0.440;
    double g211, g310, g322, g410, g422, g520, g521, g532, g533;
    if (em <= 0.65) {
      g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
      g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
      g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
      g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
      g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
      g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
    } else {
      g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
      g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
      g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
      g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
      g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
      if (em > 0.715) {
        g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
      } else {
        g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq;
      }
    }
    if (em < 0.7) {
      g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21 * eoc;
      g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
      g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4 * eoc;
    } else {
      g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
      g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
      g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
    }
    const double sini2 = sinim * sinim;
    const double f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
    const double f221 = 1.5 * sini2;
    const double f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
    const double f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
    const double f441 = 35.0 * sini2 * f220;
    const double f442 = 39.3750 * sini2 * sini2;
    const double f522 = 9.84375 * sinim *
                        (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
                         0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
    const double f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq) +
                                 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
    const double f542 = 29.53125 * sinim *
                        (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
    const double f543 = 29.53125 * sinim *
                        (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
    const double xno2 = nm * nm;
    const double ainv2 = aonv * aonv;
    double temp1 = 3.0 * xno2 * ainv2;
    double temp = temp1 * root22;
    s->d2201 = temp * f220 * g201;
    s->d2211 = temp * f221 * g211;
    temp1 = temp1 * aonv;
    temp = temp1 * root32;
    s->d3210 = temp * f321 * g310;
    s->d3222 = temp * f322 * g322;
    temp1 = temp1 * aonv;
    temp = 2.0 * temp1 * root44;
    s->d4410 = temp * f441 * g410;
    s->d4422 = temp * f442 * g422;
    temp1 = temp1 * aonv;
    temp = temp1 * root52;
    s->d5220 = temp * f522 * g520;
    s->d5232 = temp * f523 * g532;
    temp = 2.0 * temp1 * root54;
    s->d5421 = temp * f542 * g521;
    s->d5433 = temp * f543 * g533;
    s->xlamo = fmod_exact(s->mo + s->nodeo + s->nodeo - theta - theta, TWOPI);
    s->xfact = s->mdot + s->dmdt + 2.0 * (s->nodedot + s->dnodt - rptim) - s->no_unkozai;
    em = emo;
    emsq = emsqo;
  }
  if (s->irez == 1) {
    /* synchronous resonance */
    const double g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
    const double g310 = 1.0 + 2.0 * emsq;
    const double g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
    const double f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
    const double f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
    double f330 = 1.0 + cosim;
    f330 = 1.875 * f330 * f330 * f330;
    s->del1 = 3.0 * nm * nm * aonv * aonv;
    s->del2 = 2.0 * s->del1 * f220 * g200 * q22;
    s->del3 = 3.0 * s->del1 * f330 * g300 * q33 * aonv;
    s->del1 = s->del1 * f311 * g310 * q31 * aonv;
    s->xlamo = fmod_exact(s->mo + s->nodeo + s->argpo - theta, TWOPI);
    s->xfact = s->mdot + xpidot - rptim + s->dmdt + s->domdt + s->dnodt - s->no_unkozai;
  }
}

/* ---- Deep-space secular effects and resonance integration (dspace) ---------------------------- */

static void dspace(const myb_sgp4* s, double t, double* em, double* argpm, double* inclm,
                   double* mm, double* nodem, double* nm) {
  const double fasx2 = 0.13130908;
  const double fasx4 = 2.8843198;
  const double fasx6 = 0.37448087;
  const double g22 = 5.7686396;
  const double g32 = 0.95240898;
  const double g44 = 1.8014998;
  const double g52 = 1.0508330;
  const double g54 = 4.4108898;
  const double rptim = 4.37526908801129966e-3;
  const double stepp = 720.0;
  const double stepn = -720.0;
  const double step2 = 259200.0;

  const double tc = t;
  const double theta = fmod_exact(s->gsto + tc * rptim, TWOPI);
  *em = *em + s->dedt * t;
  *inclm = *inclm + s->didt * t;
  *argpm = *argpm + s->domdt * t;
  *nodem = *nodem + s->dnodt * t;
  *mm = *mm + s->dmdt * t;

  if (s->irez == 0) {
    return;
  }
  /* Euler-Maclaurin integration of the resonance terms from the epoch in steps of 720 min */
  double atime = 0.0;
  double xni = s->no_unkozai;
  double xli = s->xlamo;
  const double delt = t > 0.0 ? stepp : stepn;
  double ft = 0.0;
  double xndt = 0.0;
  double xldot = 0.0;
  double xnddt = 0.0;
  int done = 0;
  const double steps = fabs(t) / stepp + 2.0; /* bound: |t| / 720 steps, then the last point */
  for (double k = 0.0; k < steps && !done; k += 1.0) {
    if (s->irez != 2) {
      /* near-synchronous resonance terms */
      xndt = s->del1 * myb_sin(xli - fasx2) + s->del2 * myb_sin(2.0 * (xli - fasx4)) +
             s->del3 * myb_sin(3.0 * (xli - fasx6));
      xldot = xni + s->xfact;
      xnddt = s->del1 * myb_cos(xli - fasx2) + 2.0 * s->del2 * myb_cos(2.0 * (xli - fasx4)) +
              3.0 * s->del3 * myb_cos(3.0 * (xli - fasx6));
      xnddt = xnddt * xldot;
    } else {
      /* near half-day resonance terms */
      const double xomi = s->argpo + s->argpdot * atime;
      const double x2omi = xomi + xomi;
      const double x2li = xli + xli;
      xndt = s->d2201 * myb_sin(x2omi + xli - g22) + s->d2211 * myb_sin(xli - g22) +
             s->d3210 * myb_sin(xomi + xli - g32) + s->d3222 * myb_sin(-xomi + xli - g32) +
             s->d4410 * myb_sin(x2omi + x2li - g44) + s->d4422 * myb_sin(x2li - g44) +
             s->d5220 * myb_sin(xomi + xli - g52) + s->d5232 * myb_sin(-xomi + xli - g52) +
             s->d5421 * myb_sin(xomi + x2li - g54) + s->d5433 * myb_sin(-xomi + x2li - g54);
      xldot = xni + s->xfact;
      xnddt = s->d2201 * myb_cos(x2omi + xli - g22) + s->d2211 * myb_cos(xli - g22) +
              s->d3210 * myb_cos(xomi + xli - g32) + s->d3222 * myb_cos(-xomi + xli - g32) +
              s->d5220 * myb_cos(xomi + xli - g52) + s->d5232 * myb_cos(-xomi + xli - g52) +
              2.0 * (s->d4410 * myb_cos(x2omi + x2li - g44) + s->d4422 * myb_cos(x2li - g44) +
                     s->d5421 * myb_cos(xomi + x2li - g54) +
                     s->d5433 * myb_cos(-xomi + x2li - g54));
      xnddt = xnddt * xldot;
    }
    if (fabs(t - atime) >= stepp) {
      xli = xli + xldot * delt + xndt * step2;
      xni = xni + xndt * delt + xnddt * step2;
      atime = atime + delt;
    } else {
      ft = t - atime;
      done = 1;
    }
  }
  MYB_CHECK(done);
  *nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
  const double xl = xli + xldot * ft + xndt * ft * ft * 0.5;
  if (s->irez != 1) {
    *mm = xl - 2.0 * *nodem + 2.0 * theta;
  } else {
    *mm = xl - *nodem - *argpm + theta;
  }
  /* dndt = nm - no; nm = no + dndt: the same value up to rounding, kept as in the code */
  const double dndt = *nm - s->no_unkozai;
  *nm = s->no_unkozai + dndt;
}

/* ---- Sidereal time (gstime) ------------------------------------------------------------------- */

double myb_sgp4_gmst(double jd1, double jd2) {
  const double tut1 = ((jd1 - 2451545.0) + jd2) / 36525.0;
  double temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
                (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841; /* s */
  temp = fmod_exact(temp * DEG2RAD / 240.0, TWOPI); /* 360 / 86400 = 1 / 240 */
  if (temp < 0.0) {
    temp += TWOPI;
  }
  return temp;
}

/* ---- Initialisation (initl, sgp4init) --------------------------------------------------------- */

void myb_sgp4_constants_of(myb_sgp4_gravity gravity, myb_sgp4_constants* c) {
  MYB_CHECK(c != NULL && (gravity == MYB_SGP4_WGS72_OLD || gravity == MYB_SGP4_WGS72 ||
                          gravity == MYB_SGP4_WGS84));
  double mu = 0.0; /* km^3 / s^2 */
  switch (gravity) {
  case MYB_SGP4_WGS72_OLD:
    c->radius = 6378.135;
    c->xke = 0.0743669161;
    c->j2 = 0.001082616;
    c->j3 = -0.00000253881;
    c->j4 = -0.00000165597;
    return;
  case MYB_SGP4_WGS72:
    mu = 398600.8;
    c->radius = 6378.135;
    c->j2 = 0.001082616;
    c->j3 = -0.00000253881;
    c->j4 = -0.00000165597;
    break;
  case MYB_SGP4_WGS84:
    mu = 398600.5;
    c->radius = 6378.137;
    c->j2 = 0.00108262998905;
    c->j3 = -0.00000253215306;
    c->j4 = -0.00000161098761;
    break;
  }
  c->xke = 60.0 / sqrt(c->radius * c->radius * c->radius / mu);
}

myb_status myb_sgp4_init_elements(const myb_sgp4_constants* constants, myb_sgp4_mode mode,
                                  double epoch, double bstar, double eccentricity,
                                  double argument_of_perigee, double inclination,
                                  double mean_anomaly, double mean_motion, double node,
                                  myb_sgp4* sgp4) {
  MYB_CHECK(constants != NULL && sgp4 != NULL &&
            (mode == MYB_SGP4_IMPROVED || mode == MYB_SGP4_AFSPC));
  const myb_sgp4_constants* c = constants;
  if (!(c->radius > 0.0) || !(c->xke > 0.0) || !isfinite(c->radius) || !isfinite(c->xke) ||
      !isfinite(c->j2) || !(c->j2 != 0.0) || !isfinite(c->j3) || !isfinite(c->j4)) {
    return MYB_E_DOMAIN;
  }
  if (!isfinite(epoch) || !isfinite(bstar) || !isfinite(argument_of_perigee) ||
      !isfinite(inclination) || !isfinite(mean_anomaly) || !isfinite(node) ||
      !(eccentricity >= 0.0 && eccentricity < 1.0) || !(mean_motion > 0.0) ||
      !isfinite(mean_motion)) {
    return MYB_E_DOMAIN;
  }
  myb_sgp4* s = sgp4;
  memset(s, 0, sizeof *s);
  s->radius = c->radius;
  s->xke = c->xke;
  s->j2 = c->j2;
  s->j3 = c->j3;
  s->j4 = c->j4;
  s->tumin = 1.0 / s->xke;
  s->j3oj2 = s->j3 / s->j2;
  s->mode = mode == MYB_SGP4_AFSPC ? MYB_SGP4_AFSPC : MYB_SGP4_IMPROVED;
  s->epoch = epoch;
  s->bstar = bstar;
  s->ecco = eccentricity;
  s->argpo = argument_of_perigee;
  s->inclo = inclination;
  s->mo = mean_anomaly;
  s->no_kozai = mean_motion;
  s->nodeo = node;

  const double temp4 = 1.5e-12;
  const double ss = 78.0 / s->radius + 1.0;
  const double qzms2ttemp = (120.0 - 78.0) / s->radius;
  const double qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;

  /* initl: un-Kozai the mean motion */
  const double eccsq = s->ecco * s->ecco;
  const double omeosq = 1.0 - eccsq;
  const double rteosq = sqrt(omeosq);
  const double cosio = myb_cos(s->inclo);
  const double cosio2 = cosio * cosio;
  const double ak = pow23(s->xke / s->no_kozai);
  const double d1 = 0.75 * s->j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
  double del = d1 / (ak * ak);
  const double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
  del = d1 / (adel * adel);
  s->no_unkozai = s->no_kozai / (1.0 + del);
  const double ao = pow23(s->xke / s->no_unkozai);
  const double sinio = myb_sin(s->inclo);
  const double po = ao * omeosq;
  const double con42 = 1.0 - 5.0 * cosio2;
  s->con41 = -con42 - cosio2 - cosio2;
  const double posq = po * po;
  const double rp = ao * (1.0 - s->ecco);
  if (s->mode == MYB_SGP4_AFSPC) {
    const double ts70 = s->epoch - 7305.0;
    const double ds70 = floor(ts70 + 1.0e-8);
    const double tfrac = ts70 - ds70;
    const double c1 = 1.72027916940703639e-2;
    const double thgr70 = 1.7321343856509374;
    const double fk5r = 5.07551419432269442e-15;
    const double c1p2p = c1 + TWOPI;
    s->gsto = fmod_exact(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, TWOPI);
    if (s->gsto < 0.0) {
      s->gsto = s->gsto + TWOPI;
    }
  } else {
    s->gsto = myb_sgp4_gmst(s->epoch + 2433281.5, 0.0);
  }
  if (!(rp > 0.0)) {
    return MYB_E_DOMAIN;
  }

  /* sgp4init */
  s->isimp = 0;
  if (rp < (220.0 / s->radius + 1.0)) {
    s->isimp = 1;
  }
  double sfour = ss;
  double qzms24 = qzms2t;
  const double perige = (rp - 1.0) * s->radius;
  /* for perigees below 156 km, s and qoms2t are altered */
  if (perige < 156.0) {
    sfour = perige - 78.0;
    if (perige < 98.0) {
      sfour = 20.0;
    }
    const double qzms24temp = (120.0 - sfour) / s->radius;
    qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
    sfour = sfour / s->radius + 1.0;
  }
  const double pinvsq = 1.0 / posq;
  const double tsi = 1.0 / (ao - sfour);
  s->eta = ao * s->ecco * tsi;
  const double etasq = s->eta * s->eta;
  const double eeta = s->ecco * s->eta;
  const double psisq = fabs(1.0 - etasq);
  const double tsi2 = tsi * tsi;
  const double coef = qzms24 * (tsi2 * tsi2);
  const double coef1 = coef / (psisq * psisq * psisq * sqrt(psisq));
  const double cc2 =
      coef1 * s->no_unkozai *
      (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
       0.375 * s->j2 * tsi / psisq * s->con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
  s->cc1 = s->bstar * cc2;
  double cc3 = 0.0;
  if (s->ecco > 1.0e-4) {
    cc3 = -2.0 * coef * tsi * s->j3oj2 * s->no_unkozai * sinio / s->ecco;
  }
  s->x1mth2 = 1.0 - cosio2;
  s->cc4 = 2.0 * s->no_unkozai * coef1 * ao * omeosq *
           (s->eta * (2.0 + 0.5 * etasq) + s->ecco * (0.5 + 2.0 * etasq) -
            s->j2 * tsi / (ao * psisq) *
                (-3.0 * s->con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                 0.75 * s->x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
                     myb_cos(2.0 * s->argpo)));
  s->cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
  const double cosio4 = cosio2 * cosio2;
  const double temp1 = 1.5 * s->j2 * pinvsq * s->no_unkozai;
  const double temp2 = 0.5 * temp1 * s->j2 * pinvsq;
  const double temp3 = -0.46875 * s->j4 * pinvsq * pinvsq * s->no_unkozai;
  s->mdot = s->no_unkozai + 0.5 * temp1 * rteosq * s->con41 +
            0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
  s->argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
               temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
  const double xhdot1 = -temp1 * cosio;
  s->nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) *
                            cosio;
  const double xpidot = s->argpdot + s->nodedot;
  s->omgcof = s->bstar * cc3 * myb_cos(s->argpo);
  s->xmcof = 0.0;
  if (s->ecco > 1.0e-4) {
    s->xmcof = -X2O3 * coef * s->bstar / eeta;
  }
  s->nodecf = 3.5 * omeosq * xhdot1 * s->cc1;
  s->t2cof = 1.5 * s->cc1;
  /* division by zero for an inclination of 180 degrees */
  if (fabs(cosio + 1.0) > 1.5e-12) {
    s->xlcof = -0.25 * s->j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
  } else {
    s->xlcof = -0.25 * s->j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
  }
  s->aycof = -0.5 * s->j3oj2 * sinio;
  const double delmotemp = 1.0 + s->eta * myb_cos(s->mo);
  s->delmo = delmotemp * delmotemp * delmotemp;
  s->sinmao = myb_sin(s->mo);
  s->x7thm1 = 7.0 * cosio2 - 1.0;

  /* deep-space initialisation */
  if ((2.0 * PI / s->no_unkozai) >= 225.0) {
    s->deep = 1;
    s->isimp = 1;
    dscom_out o;
    dscom(s->epoch, s->ecco, s->argpo, 0.0, s->inclo, s->nodeo, s->no_unkozai, s, &o);
    dsinit(s, &o, xpidot, eccsq);
  }

  /* drag terms of higher order when the perigee is above 220 km */
  if (s->isimp != 1) {
    const double cc1sq = s->cc1 * s->cc1;
    s->d2 = 4.0 * ao * tsi * cc1sq;
    const double temp = s->d2 * tsi * s->cc1 / 3.0;
    s->d3 = (17.0 * ao + sfour) * temp;
    s->d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s->cc1;
    s->t3cof = s->d2 + 2.0 * cc1sq;
    s->t4cof = 0.25 * (3.0 * s->d3 + s->cc1 * (12.0 * s->d2 + 10.0 * cc1sq));
    s->t5cof = 0.2 * (3.0 * s->d4 + 12.0 * s->cc1 * s->d3 + 6.0 * s->d2 * s->d2 +
                      15.0 * cc1sq * (2.0 * s->d2 + cc1sq));
  }
  return MYB_OK;
}

/* ---- Propagation (sgp4) ----------------------------------------------------------------------- */

myb_status myb_sgp4_propagate_km(const myb_sgp4* sgp4, double minutes, double r[3], double v[3],
                                 myb_sgp4_error* error) {
  MYB_CHECK(sgp4 != NULL && r != NULL && v != NULL);
  const myb_sgp4* s = sgp4;
  if (error != NULL) {
    *error = MYB_SGP4_OK;
  }
  if (!isfinite(minutes)) {
    return MYB_E_DOMAIN;
  }
  const double temp4 = 1.5e-12;
  const double vkmpersec = s->radius * s->xke / 60.0;
  const double t = minutes;

  /* secular gravity and atmospheric drag */
  const double xmdf = s->mo + s->mdot * t;
  const double argpdf = s->argpo + s->argpdot * t;
  const double nodedf = s->nodeo + s->nodedot * t;
  double argpm = argpdf;
  double mm = xmdf;
  const double t2 = t * t;
  double nodem = nodedf + s->nodecf * t2;
  double tempa = 1.0 - s->cc1 * t;
  double tempe = s->bstar * s->cc4 * t;
  double templ = s->t2cof * t2;

  if (s->isimp != 1) {
    const double delomg = s->omgcof * t;
    const double delmtemp = 1.0 + s->eta * myb_cos(xmdf);
    const double delm = s->xmcof * (delmtemp * delmtemp * delmtemp - s->delmo);
    const double temp = delomg + delm;
    mm = xmdf + temp;
    argpm = argpdf - temp;
    const double t3 = t2 * t;
    const double t4 = t3 * t;
    tempa = tempa - s->d2 * t2 - s->d3 * t3 - s->d4 * t4;
    tempe = tempe + s->bstar * s->cc5 * (myb_sin(mm) - s->sinmao);
    templ = templ + s->t3cof * t3 + t4 * (s->t4cof + t * s->t5cof);
  }

  double nm = s->no_unkozai;
  double em = s->ecco;
  double inclm = s->inclo;
  if (s->deep) {
    dspace(s, t, &em, &argpm, &inclm, &mm, &nodem, &nm);
  }

  if (nm <= 0.0) {
    if (error != NULL) {
      *error = MYB_SGP4_MEAN_MOTION;
    }
    return MYB_E_DOMAIN;
  }
  const double am = pow23(s->xke / nm) * tempa * tempa;
  nm = s->xke / (am * sqrt(am));
  em = em - tempe;

  if ((em >= 1.0) || (em < -0.001)) {
    if (error != NULL) {
      *error = MYB_SGP4_MEAN_ECCENTRICITY;
    }
    return MYB_E_DOMAIN;
  }
  if (em < 1.0e-6) {
    em = 1.0e-6;
  }
  mm = mm + s->no_unkozai * templ;
  double xlm = mm + argpm + nodem;
  nodem = fmod_exact(nodem, TWOPI);
  argpm = fmod_exact(argpm, TWOPI);
  xlm = fmod_exact(xlm, TWOPI);
  mm = fmod_exact(xlm - argpm - nodem, TWOPI);

  /* lunar and solar periodics */
  double ep = em;
  double xincp = inclm;
  double argpp = argpm;
  double nodep = nodem;
  double mp = mm;
  double sinip = myb_sin(inclm);
  double cosip = myb_cos(inclm);
  double aycof = s->aycof;
  double xlcof = s->xlcof;
  double con41 = s->con41;
  double x1mth2 = s->x1mth2;
  double x7thm1 = s->x7thm1;
  if (s->deep) {
    dpper(s, t, &ep, &xincp, &nodep, &argpp, &mp);
    if (xincp < 0.0) {
      xincp = -xincp;
      nodep = nodep + PI;
      argpp = argpp - PI;
    }
    if ((ep < 0.0) || (ep > 1.0)) {
      if (error != NULL) {
        *error = MYB_SGP4_PERTURBED_ECCENTRICITY;
      }
      return MYB_E_DOMAIN;
    }
    /* long-period periodics */
    sinip = myb_sin(xincp);
    cosip = myb_cos(xincp);
    aycof = -0.5 * s->j3oj2 * sinip;
    if (fabs(cosip + 1.0) > 1.5e-12) {
      xlcof = -0.25 * s->j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
    } else {
      xlcof = -0.25 * s->j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
    }
  }
  const double axnl = ep * myb_cos(argpp);
  double temp = 1.0 / (am * (1.0 - ep * ep));
  const double aynl = ep * myb_sin(argpp) + temp * aycof;
  const double xl = mp + argpp + nodep + temp * xlcof * axnl;

  /* Kepler's equation */
  const double u = fmod_exact(xl - nodep, TWOPI);
  double eo1 = u;
  double tem5 = 9999.9;
  double sineo1 = 0.0;
  double coseo1 = 0.0;
  for (int ktr = 1; fabs(tem5) >= 1.0e-12 && ktr <= 10; ++ktr) { /* bound: 10 iterations */
    sineo1 = myb_sin(eo1);
    coseo1 = myb_cos(eo1);
    tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
    tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
    if (fabs(tem5) >= 0.95) {
      tem5 = tem5 > 0.0 ? 0.95 : -0.95;
    }
    eo1 = eo1 + tem5;
  }

  /* short-period periodics */
  const double ecose = axnl * coseo1 + aynl * sineo1;
  const double esine = axnl * sineo1 - aynl * coseo1;
  const double el2 = axnl * axnl + aynl * aynl;
  const double pl = am * (1.0 - el2);
  if (pl < 0.0) {
    if (error != NULL) {
      *error = MYB_SGP4_SEMILATUS_RECTUM;
    }
    return MYB_E_DOMAIN;
  }
  const double rl = am * (1.0 - ecose);
  const double rdotl = sqrt(am) * esine / rl;
  const double rvdotl = sqrt(pl) / rl;
  const double betal = sqrt(1.0 - el2);
  temp = esine / (1.0 + betal);
  const double sinu = am / rl * (sineo1 - aynl - axnl * temp);
  const double cosu = am / rl * (coseo1 - axnl + aynl * temp);
  double su = myb_atan2(sinu, cosu);
  const double sin2u = (cosu + cosu) * sinu;
  const double cos2u = 1.0 - 2.0 * sinu * sinu;
  temp = 1.0 / pl;
  const double temp1 = 0.5 * s->j2 * temp;
  const double temp2 = temp1 * temp;

  if (s->deep) {
    const double cosisq = cosip * cosip;
    con41 = 3.0 * cosisq - 1.0;
    x1mth2 = 1.0 - cosisq;
    x7thm1 = 7.0 * cosisq - 1.0;
  }
  const double mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) + 0.5 * temp1 * x1mth2 * cos2u;
  su = su - 0.25 * temp2 * x7thm1 * sin2u;
  const double xnode = nodep + 1.5 * temp2 * cosip * sin2u;
  const double xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
  const double mvt = rdotl - nm * temp1 * x1mth2 * sin2u / s->xke;
  const double rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / s->xke;

  /* orientation vectors */
  const double sinsu = myb_sin(su);
  const double cossu = myb_cos(su);
  const double snod = myb_sin(xnode);
  const double cnod = myb_cos(xnode);
  const double sini = myb_sin(xinc);
  const double cosi = myb_cos(xinc);
  const double xmx = -snod * cosi;
  const double xmy = cnod * cosi;
  const double ux = xmx * sinsu + cnod * cossu;
  const double uy = xmy * sinsu + snod * cossu;
  const double uz = sini * sinsu;
  const double vx = xmx * cossu - cnod * sinsu;
  const double vy = xmy * cossu - snod * sinsu;
  const double vz = sini * cossu;

  /* position and velocity (km, km/s) */
  r[0] = (mrt * ux) * s->radius;
  r[1] = (mrt * uy) * s->radius;
  r[2] = (mrt * uz) * s->radius;
  v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
  v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
  v[2] = (mvt * uz + rvdot * vz) * vkmpersec;

  if (mrt < 1.0 && error != NULL) {
    *error = MYB_SGP4_DECAYED;
  }
  return MYB_OK;
}

myb_status myb_sgp4_propagate(const myb_sgp4* sgp4, double minutes, double r[3], double v[3],
                              myb_sgp4_error* error) {
  double rk[3];
  double vk[3];
  const myb_status st = myb_sgp4_propagate_km(sgp4, minutes, rk, vk, error);
  if (st != MYB_OK) {
    return st;
  }
  for (int i = 0; i < 3; ++i) {
    r[i] = 1000.0 * rk[i];
    v[i] = 1000.0 * vk[i];
  }
  return MYB_OK;
}

/* ---- Two-line elements ------------------------------------------------------------------------ */

static size_t line_length(const char* line) {
  size_t n = 0;
  while (n < 80 && line[n] != '\0' && line[n] != '\n' && line[n] != '\r') { /* bound: 80 */
    ++n;
  }
  return n;
}

/* A number from the columns [from, to) with blanks removed; MYB_E_FORMAT if none. */
static myb_status field_double(const char* line, size_t from, size_t to, double* value) {
  char text[24];
  size_t n = 0;
  for (size_t i = from; i < to && n + 1 < sizeof text; ++i) {
    if (line[i] != ' ') {
      text[n++] = line[i];
    }
  }
  if (n == 0) {
    return MYB_E_FORMAT;
  }
  size_t used = 0;
  const myb_status st = myb_parse_double(text, n, value, &used);
  return st == MYB_OK && used == n ? MYB_OK : MYB_E_FORMAT;
}

static myb_status field_long(const char* line, size_t from, size_t to, long* value) {
  long v = 0;
  long sign = 1;
  int digits = 0;
  for (size_t i = from; i < to; ++i) {
    const char c = line[i];
    if (c == ' ') {
      continue;
    }
    if ((c == '-' || c == '+') && digits == 0) {
      sign = c == '-' ? -1 : 1;
      continue;
    }
    if (c < '0' || c > '9') {
      return MYB_E_FORMAT;
    }
    v = 10 * v + (c - '0');
    ++digits;
  }
  *value = sign * v;
  return MYB_OK;
}

static int checksum_ok(const char* line, size_t length) {
  if (length < 69 || line[68] < '0' || line[68] > '9') {
    return 1; /* no checksum */
  }
  int sum = 0;
  for (size_t i = 0; i < 68; ++i) {
    if (line[i] >= '0' && line[i] <= '9') {
      sum += line[i] - '0';
    } else if (line[i] == '-') {
      sum += 1;
    }
  }
  return sum % 10 == line[68] - '0';
}

/* 10^e for |e| <= 9, exact for e >= 0 and correctly rounded for e < 0, like pow(10.0, e). */
static double ten_to(long e) {
  static const double POS[10] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9};
  static const double NEG[10] = {1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9};
  return e >= 0 ? POS[e] : NEG[-e];
}

static myb_status bad(myb_diag* diag, int line, size_t column, const char* what) {
  myb_internal_diag(diag, line, "column %zu: %s", column + 1, what);
  return MYB_E_FORMAT;
}

myb_status myb_tle_parse(const char* line1, const char* line2, unsigned flags, myb_tle* tle,
                         myb_diag* diag) {
  MYB_CHECK(line1 != NULL && line2 != NULL && tle != NULL);
  myb_internal_diag_clear(diag);
  const size_t n1 = line_length(line1);
  const size_t n2 = line_length(line2);
  if (n1 < 64 || line1[0] != '1' || line1[1] != ' ' || line1[8] != ' ' || line1[23] != '.' ||
      line1[32] != ' ' || line1[34] != '.' || line1[43] != ' ' || line1[52] != ' ' ||
      line1[61] != ' ' || line1[63] != ' ') {
    return bad(diag, 1, 0, "not the first line of a TLE");
  }
  if (n2 < 68 || line2[0] != '2' || line2[1] != ' ' || line2[7] != ' ' || line2[11] != '.' ||
      line2[16] != ' ' || line2[20] != '.' || line2[25] != ' ' || line2[33] != ' ' ||
      line2[37] != '.' || line2[42] != ' ' || line2[46] != '.' || line2[51] != ' ') {
    return bad(diag, 2, 0, "not the second line of a TLE");
  }
  const int verify = (flags & MYB_TLE_NO_CHECKSUM) == 0;
  if (verify && !checksum_ok(line1, n1)) {
    return bad(diag, 1, 68, "checksum");
  }
  if (verify && !checksum_ok(line2, n2)) {
    return bad(diag, 2, 68, "checksum");
  }
  if (memcmp(line1 + 2, line2 + 2, 5) != 0) {
    return bad(diag, 2, 2, "catalogue numbers of the two lines differ");
  }
  myb_tle t;
  memset(&t, 0, sizeof t);
  memcpy(t.satnum, line1 + 2, 5);
  t.satnum[5] = '\0';
  /* Alpha-5: a leading letter A..Z without I and O stands for 10..33 */
  long first = 0;
  const char c0 = line1[2];
  if (c0 >= 'A' && c0 <= 'Z' && c0 != 'I' && c0 != 'O') {
    first = 10 + (c0 - 'A') - (c0 > 'I') - (c0 > 'O');
  } else if (c0 >= '0' && c0 <= '9') {
    first = c0 - '0';
  } else if (c0 != ' ') {
    return bad(diag, 1, 2, "catalogue number");
  }
  long rest = 0;
  if (field_long(line1, 3, 7, &rest) != MYB_OK || rest < 0) {
    return bad(diag, 1, 3, "catalogue number");
  }
  t.catalog_number = first * 10000 + rest;
  t.classification = line1[7] == ' ' ? 'U' : line1[7];
  memcpy(t.designator, line1 + 9, 8);
  t.designator[8] = '\0';
  for (int i = 7; i >= 0 && t.designator[i] == ' '; --i) {
    t.designator[i] = '\0';
  }
  long year2 = 0;
  if (field_long(line1, 18, 20, &year2) != MYB_OK || year2 < 0) {
    return bad(diag, 1, 18, "epoch year");
  }
  t.epoch_year = (int)(year2 < 57 ? year2 + 2000 : year2 + 1900);
  if (field_double(line1, 20, 32, &t.epoch_day) != MYB_OK || !(t.epoch_day >= 1.0) ||
      !(t.epoch_day < 367.0)) {
    return bad(diag, 1, 20, "epoch day");
  }
  if (field_double(line1, 33, 43, &t.ndot) != MYB_OK) {
    return bad(diag, 1, 33, "first derivative of the mean motion");
  }
  /* "[sign]ddddd[sign]d": 0.ddddd x 10^e */
  char text[8];
  long exponent = 0;
  double mantissa = 0.0;
  text[0] = line1[44] == ' ' ? '+' : line1[44];
  text[1] = '.';
  memcpy(text + 2, line1 + 45, 5);
  if (field_double(text, 0, 7, &mantissa) != MYB_OK ||
      field_long(line1, 50, 52, &exponent) != MYB_OK || exponent < -9 || exponent > 9) {
    return bad(diag, 1, 44, "second derivative of the mean motion");
  }
  t.nddot = mantissa * ten_to(exponent);
  text[0] = line1[53] == ' ' ? '+' : line1[53];
  memcpy(text + 2, line1 + 54, 5);
  if (field_double(text, 0, 7, &mantissa) != MYB_OK ||
      field_long(line1, 59, 61, &exponent) != MYB_OK || exponent < -9 || exponent > 9) {
    return bad(diag, 1, 53, "drag term");
  }
  t.bstar = mantissa * ten_to(exponent);
  t.ephemeris_type = line1[62] >= '0' && line1[62] <= '9' ? line1[62] - '0' : 0;
  long number = 0;
  if (field_long(line1, 64, 68 < n1 ? 68 : n1, &number) != MYB_OK) {
    return bad(diag, 1, 64, "element set number");
  }
  t.element_number = (int)number;

  char ecc[10] = "0.";
  for (int i = 0; i < 7; ++i) {
    ecc[2 + i] = line2[26 + i] == ' ' ? '0' : line2[26 + i];
  }
  if (field_double(line2, 8, 16, &t.inclination) != MYB_OK) {
    return bad(diag, 2, 8, "inclination");
  }
  if (field_double(line2, 17, 25, &t.node) != MYB_OK) {
    return bad(diag, 2, 17, "right ascension of the node");
  }
  if (field_double(ecc, 0, 9, &t.eccentricity) != MYB_OK) {
    return bad(diag, 2, 26, "eccentricity");
  }
  if (field_double(line2, 34, 42, &t.argument_of_perigee) != MYB_OK) {
    return bad(diag, 2, 34, "argument of perigee");
  }
  if (field_double(line2, 43, 51, &t.mean_anomaly) != MYB_OK) {
    return bad(diag, 2, 43, "mean anomaly");
  }
  if (field_double(line2, 52, 63, &t.mean_motion) != MYB_OK) {
    return bad(diag, 2, 52, "mean motion");
  }
  long revolution = 0;
  if (field_long(line2, 63, 68, &revolution) != MYB_OK) {
    return bad(diag, 2, 63, "revolution number");
  }
  t.revolution = revolution;
  *tle = t;
  return MYB_OK;
}

static int leap_year(int y) {
  return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
}

/* The epoch in days since 1949 December 31 0h as the reference code computes it (days2mdhms,
 * jday, then (jd + jdFrac) - 2433281.5): its rounding at the scale of a Julian date (2^-31 day,
 * 40 microseconds) is part of the verification output — the lunar and solar terms of highly
 * eccentric orbits move by kilometres per second of epoch. */
static double reference_epoch(int year, double days) {
  static const int LMONTH[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const int dayofyr = (int)floor(days);
  int i = 1;
  int inttemp = 0;
  for (; dayofyr > inttemp + (LMONTH[i] + (i == 2 && year % 4 == 0)) && i < 12; ++i) {
    inttemp = inttemp + LMONTH[i] + (i == 2 && year % 4 == 0);
  }
  const int mon = i;
  const int day = dayofyr - inttemp;
  double temp = (days - dayofyr) * 24.0;
  const int hr = (int)floor(temp);
  temp = (temp - hr) * 60.0;
  const int minute = (int)floor(temp);
  const double sec = (temp - minute) * 60.0;
  double jd = 367.0 * year - floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) +
              floor(275 * mon / 9.0) + day + 1721013.5;
  double jd_frac = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;
  if (fabs(jd_frac) > 1.0) {
    const double dtt = floor(jd_frac);
    jd = jd + dtt;
    jd_frac = jd_frac - dtt;
  }
  return (jd + jd_frac) - 2433281.5;
}

void myb_tle_epoch(const myb_tle* tle, myb_calendar* utc) {
  MYB_CHECK(tle != NULL && utc != NULL);
  static const int MONTH[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const double whole = floor(tle->epoch_day);
  int day = (int)whole;
  const double seconds = (tle->epoch_day - whole) * 86400.0;
  int month = 0;
  for (; month < 11; ++month) {
    const int length = MONTH[month] + (month == 1 && leap_year(tle->epoch_year));
    if (day <= length) {
      break;
    }
    day -= length;
  }
  const double hours = floor(seconds / 3600.0);
  const double minutes = floor((seconds - 3600.0 * hours) / 60.0);
  *utc = (myb_calendar){.year = tle->epoch_year,
                        .month = month + 1,
                        .day = day,
                        .hour = (int)hours,
                        .minute = (int)minutes,
                        .second = seconds - 3600.0 * hours - 60.0 * minutes};
}

myb_status myb_sgp4_init(const myb_tle* tle, myb_sgp4_gravity gravity, myb_sgp4_mode mode,
                         myb_sgp4* sgp4) {
  MYB_CHECK(tle != NULL && sgp4 != NULL);
  const double xpdotp = 1440.0 / (2.0 * PI); /* 229.1831180523293 */
  const double epoch = reference_epoch(tle->epoch_year, tle->epoch_day);
  myb_sgp4_constants c;
  myb_sgp4_constants_of(gravity, &c);
  return myb_sgp4_init_elements(&c, mode, epoch, tle->bstar, tle->eccentricity,
                                tle->argument_of_perigee * DEG2RAD, tle->inclination * DEG2RAD,
                                tle->mean_anomaly * DEG2RAD, tle->mean_motion / xpdotp,
                                tle->node * DEG2RAD, sgp4);
}

/* ---- SPK type 10 ------------------------------------------------------------------------------ */

/* Generic segment meta data (NAIF SGPARAM): the last word of the segment holds the number of
 * items (15, 16 or 17), which precede it; base addresses are relative to the segment start. */
enum {
  META_CONBAS = 0,
  META_NCON = 1,
  META_REFBAS = 5,
  META_NREF = 6,
  META_PKTBAS = 10,
  META_NPKT = 11,
  META_PKTSZ = 14,
  META_PKTOFF = 15,
  META_ITEMS = 17,
  GEOPHS = 8,
  ELEMENTS = 10
};

typedef struct segment10 {
  const myb_spk* spk;
  size_t index;
  size_t meta[META_ITEMS];
  double geophs[GEOPHS];
} segment10;

static int as_count(double v, size_t limit, size_t* out) {
  if (!(v >= 0.0 && v <= (double)limit) || v != floor(v)) {
    return 0;
  }
  *out = (size_t)v;
  return 1;
}

static myb_status open_segment10(const myb_spk* spk, size_t index, segment10* g) {
  myb_spk_segment_info info;
  myb_spk_segment(spk, index, &info);
  if (info.type != 10) {
    return MYB_E_FORMAT;
  }
  const size_t size = myb_spk_segment_size(spk, index);
  double last = 0.0;
  size_t count = 0;
  if (size < 16 || myb_spk_segment_data(spk, index, size - 1, 1, &last) != MYB_OK ||
      !as_count(last, size, &count) || count < 15) {
    return MYB_E_FORMAT;
  }
  /* 15 items are stored with the count after them (16 words); more than 17: the first 17 */
  const size_t stored = count == 15 ? 16 : count;
  const size_t used = count > META_ITEMS ? META_ITEMS : stored;
  double meta[META_ITEMS] = {0};
  if (stored > size || myb_spk_segment_data(spk, index, size - stored, used, meta) != MYB_OK) {
    return MYB_E_FORMAT;
  }
  g->spk = spk;
  g->index = index;
  for (size_t i = 0; i < META_ITEMS; ++i) {
    g->meta[i] = 0;
    if (i < used && i < 16 && !as_count(meta[i], size, &g->meta[i])) {
      return MYB_E_FORMAT;
    }
  }
  if (count <= 16) {
    g->meta[META_PKTOFF] = 0;
  }
  const size_t* m = g->meta;
  if (m[META_NCON] != GEOPHS || m[META_NREF] == 0 || m[META_NPKT] != m[META_NREF] ||
      (m[META_PKTSZ] != 14 && m[META_PKTSZ] != 10) || m[META_PKTOFF] > 16 ||
      m[META_CONBAS] + GEOPHS > size || m[META_REFBAS] + m[META_NREF] > size ||
      m[META_PKTBAS] + m[META_NPKT] * (m[META_PKTSZ] + m[META_PKTOFF]) > size) {
    return MYB_E_FORMAT;
  }
  return myb_spk_segment_data(spk, index, m[META_CONBAS], GEOPHS, g->geophs) == MYB_OK
             ? MYB_OK
             : MYB_E_FORMAT;
}

static double epoch10(const segment10* g, size_t k) { /* 0-based */
  double t = 0.0;
  const myb_status st = myb_spk_segment_data(g->spk, g->index, g->meta[META_REFBAS] + k, 1, &t);
  MYB_CHECK(st == MYB_OK); /* within the bounds checked by open_segment10 */
  return t;
}

/* TEME state (km, km/s) of set k at et. */
static myb_status set_state(const segment10* g, size_t k, const myb_timesys* system,
                            myb_scale tdb, myb_time et, double state[6]) {
  double e[ELEMENTS];
  /* fixed packets of PKTSZ words, each preceded by PKTOFF words (SGFPKT) */
  const size_t at = g->meta[META_PKTBAS] + k * (g->meta[META_PKTSZ] + g->meta[META_PKTOFF]) +
                    g->meta[META_PKTOFF];
  if (myb_spk_segment_data(g->spk, g->index, at, ELEMENTS, e) != MYB_OK) {
    return MYB_E_FORMAT;
  }
  /* the epoch as a UTC Julian date in one double, like SPICE (ttrans TDB -> JDUTC) */
  myb_calendar utc;
  myb_status st = myb_timesys_to_utc(system, myb_time_from_seconds(e[9]), tdb, &utc);
  if (st != MYB_OK) {
    return st;
  }
  /* days from 2000 January 1.5 of the date, then the seconds of the day */
  const int64_t y = utc.year;
  const int64_t mo = utc.month;
  const int64_t a = (14 - mo) / 12;
  const int64_t yy = y + 4800 - a;
  const int64_t mm = mo + 12 * a - 3;
  const int64_t jdn = utc.day + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
  const double utc_seconds = (double)(jdn - 2451545) * 86400.0 - 43200.0 +
                             (utc.hour * 3600.0 + utc.minute * 60.0 + utc.second);
  const double jd_utc = 2451545.0 + utc_seconds / 86400.0;
  const myb_sgp4_constants c = {.radius = g->geophs[6],
                                .xke = g->geophs[3],
                                .j2 = g->geophs[0],
                                .j3 = g->geophs[1],
                                .j4 = g->geophs[2]};
  myb_sgp4 s;
  st = myb_sgp4_init_elements(&c, MYB_SGP4_AFSPC, jd_utc - 2433281.5, e[2], e[5], e[6], e[3],
                              e[7], e[8], e[4], &s);
  if (st != MYB_OK) {
    return st;
  }
  const double minutes = myb_dd_to_double(myb_time_diff(et, myb_time_from_seconds(e[9]))) / 60.0;
  return myb_sgp4_propagate_km(&s, minutes, state, state + 3, NULL);
}

myb_status myb_sgp4_spk10_state(const myb_spk* spk, size_t index, const myb_timesys* system,
                                const myb_equinox* equinox, myb_time et, double state[6]) {
  MYB_CHECK(spk != NULL && system != NULL && equinox != NULL && state != NULL &&
            index < myb_spk_segment_count(spk));
  segment10 g;
  myb_status st = open_segment10(spk, index, &g);
  if (st != MYB_OK) {
    return st;
  }
  myb_scale tdb;
  st = myb_timesys_scale(system, "TDB", &tdb);
  if (st != MYB_OK) {
    return st;
  }
  /* the bracketing sets */
  const size_t n = g.meta[META_NREF];
  const double x = myb_time_seconds(et);
  size_t first = 0;
  size_t second = 0;
  if (x <= epoch10(&g, 0)) {
    first = second = 0;
  } else if (x >= epoch10(&g, n - 1)) {
    first = second = n - 1;
  } else {
    /* epoch[lo] < x < epoch[hi] with hi = lo + 1 at the end */
    size_t lo = 0;
    size_t hi = n - 1;
    while (hi - lo > 1) { /* bound: log2(n) halvings */
      const size_t mid = lo + (hi - lo) / 2;
      if (epoch10(&g, mid) <= x) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    first = lo;
    second = hi;
  }
  double s1[6];
  st = set_state(&g, first, system, tdb, et, s1);
  if (st != MYB_OK) {
    return st;
  }
  double teme[6];
  if (first == second) {
    memcpy(teme, s1, sizeof teme);
  } else {
    double s2[6];
    st = set_state(&g, second, system, tdb, et, s2);
    if (st != MYB_OK) {
      return st;
    }
    /* W(t) = 0.5 + 0.5 cos(pi (t - t1) / (t2 - t1)), as SPKE10 */
    const double t1 = epoch10(&g, first);
    const double t2 = epoch10(&g, second);
    const double numer = myb_dd_to_double(myb_time_diff(et, myb_time_from_seconds(t1)));
    const double denom = t2 - t1;
    const double arg = numer * PI / denom;
    const double dargdt = PI / denom;
    const double w = myb_cos(arg) * 0.5 + 0.5;
    const double dwdt = myb_sin(arg) * -0.5 * dargdt;
    for (int i = 0; i < 6; ++i) {
      teme[i] = w * s1[i] + (1.0 - w) * s2[i];
    }
    for (int i = 0; i < 3; ++i) {
      teme[3 + i] = teme[3 + i] + (dwdt * s1[i] + -dwdt * s2[i]);
    }
  }
  /* TEME of et -> J2000 */
  myb_rotation to_teme;
  myb_equinox_rotation(equinox, MYB_EQUINOX_TEME, et, &to_teme);
  myb_rotation from_teme;
  myb_rotation_invert(&to_teme, &from_teme);
  double j2000[6];
  myb_rotation_apply_state(&from_teme, teme, j2000);
  for (int i = 0; i < 6; ++i) {
    state[i] = 1000.0 * j2000[i];
  }
  return MYB_OK;
}

/* ---- Writing element sets ---------------------------------------------------------------------- */

/* The catalogue number in 5 characters: digits, or a letter for the tens of thousands above
 * 99 999 (Alpha-5: A = 10 .. Z = 33 without I and O). */
static myb_status satnum_of(long number, char out[6]) {
  static const char letters[] = "ABCDEFGHJKLMNPQRSTUVWXYZ";
  if (number < 0 || number > 339999) {
    return MYB_E_DOMAIN;
  }
  long rest = number % 10000;
  const long high = number / 10000;
  out[0] = high < 10 ? (char)('0' + high) : letters[high - 10];
  for (int i = 4; i >= 1; --i) {
    out[i] = (char)('0' + rest % 10);
    rest /= 10;
  }
  out[5] = '\0';
  return MYB_OK;
}

/* value with `decimals` digits after the point, right-aligned in `width` characters padded
 * with `pad`; MYB_E_DOMAIN if it does not fit. */
static myb_status fixed_field(double value, int decimals, size_t width, char pad, char* out) {
  char text[64];
  size_t length = 0;
  if (!isfinite(value) || myb_format_fixed(value, decimals, text, sizeof text, &length) != MYB_OK ||
      length > width) {
    return MYB_E_DOMAIN;
  }
  memset(out, pad, width - length);
  memcpy(out + width - length, text, length);
  return MYB_OK;
}

/* +-.ddddd x 10^e in 8 characters: sign, five digits, the sign and the digit of e. A zero is
 * written with the exponent sign the catalogue uses in that field ("00000-0" for the second
 * derivative of the mean motion, "00000+0" for the drag term). */
static myb_status exponent_field(double value, char zero_sign, char out[8]) {
  if (!isfinite(value)) {
    return MYB_E_DOMAIN;
  }
  out[0] = value < 0.0 ? '-' : ' ';
  if (value == 0.0) {
    memcpy(out + 1, "00000", 5);
    out[6] = zero_sign;
    out[7] = '0';
    return MYB_OK;
  }
  char text[32];
  size_t length = 0;
  if (myb_format_scientific(fabs(value), 5, text, sizeof text, &length) != MYB_OK || length < 10) {
    return MYB_E_DOMAIN;
  }
  /* d.dddde+XX = .ddddd x 10^(XX + 1) */
  int exponent = 0;
  for (size_t i = 8; i < length; ++i) {
    exponent = 10 * exponent + (text[i] - '0');
  }
  exponent = (text[7] == '-' ? -exponent : exponent) + 1;
  if (exponent < -9 || exponent > 9) {
    return MYB_E_DOMAIN;
  }
  out[1] = text[0];
  memcpy(out + 2, text + 2, 4);
  out[6] = exponent > 0 ? '+' : '-';
  out[7] = (char)('0' + (exponent < 0 ? -exponent : exponent));
  return MYB_OK;
}

static myb_status integer_field(long value, size_t width, char* out) {
  if (value < 0) {
    return MYB_E_DOMAIN;
  }
  memset(out, ' ', width);
  size_t at = width;
  long rest = value;
  /* bound: the width of the field */
  do {
    if (at == 0) {
      return MYB_E_DOMAIN;
    }
    out[--at] = (char)('0' + rest % 10);
    rest /= 10;
  } while (rest != 0);
  return MYB_OK;
}

static char checksum_of(const char* line) {
  int sum = 0;
  for (int i = 0; i < 68; ++i) {
    if (line[i] >= '0' && line[i] <= '9') {
      sum += line[i] - '0';
    } else if (line[i] == '-') {
      sum += 1;
    }
  }
  return (char)('0' + sum % 10);
}

myb_status myb_tle_format(const myb_tle* tle, char line1[70], char line2[70]) {
  MYB_CHECK(tle != NULL && line1 != NULL && line2 != NULL);
  char satnum[6];
  if (strlen(tle->satnum) == 5) {
    memcpy(satnum, tle->satnum, 6);
  } else if (satnum_of(tle->catalog_number, satnum) != MYB_OK) {
    return MYB_E_DOMAIN;
  }
  const size_t designator = strlen(tle->designator);
  if (designator > 8 || tle->classification < 32 || tle->classification > 126 ||
      tle->epoch_year < 1957 || tle->epoch_year > 2056 || !(tle->epoch_day >= 1.0) ||
      !(tle->epoch_day < 367.0) || !(fabs(tle->ndot) < 1.0) || tle->ephemeris_type < 0 ||
      tle->ephemeris_type > 9 || !(tle->eccentricity >= 0.0) || !(tle->eccentricity < 1.0)) {
    return MYB_E_DOMAIN;
  }
  memset(line1, ' ', 69);
  memset(line2, ' ', 69);
  line1[69] = '\0';
  line2[69] = '\0';
  line1[0] = '1';
  memcpy(line1 + 2, satnum, 5);
  line1[7] = tle->classification;
  memcpy(line1 + 9, tle->designator, designator);
  line1[18] = (char)('0' + (tle->epoch_year % 100) / 10);
  line1[19] = (char)('0' + tle->epoch_year % 10);
  char ndot[12];
  myb_status status = fixed_field(tle->epoch_day, 8, 12, '0', line1 + 20);
  if (status == MYB_OK) {
    status = fixed_field(fabs(tle->ndot), 8, 10, ' ', ndot);
  }
  if (status == MYB_OK && ndot[0] != '0') {
    status = MYB_E_DOMAIN; /* rounds to 1 */
  }
  if (status == MYB_OK) {
    /* "0.dddddddd" without its zero, the sign in its place */
    line1[33] = tle->ndot < 0.0 ? '-' : ' ';
    memcpy(line1 + 34, ndot + 1, 9);
    status = exponent_field(tle->nddot, '-', line1 + 44);
  }
  if (status == MYB_OK) {
    status = exponent_field(tle->bstar, '+', line1 + 53);
  }
  if (status == MYB_OK) {
    line1[62] = (char)('0' + tle->ephemeris_type);
    status = integer_field(tle->element_number, 4, line1 + 64);
  }
  line2[0] = '2';
  memcpy(line2 + 2, satnum, 5);
  char eccentricity[12];
  if (status == MYB_OK) {
    status = fixed_field(tle->inclination, 4, 8, ' ', line2 + 8);
  }
  if (status == MYB_OK) {
    status = fixed_field(tle->node, 4, 8, ' ', line2 + 17);
  }
  if (status == MYB_OK) {
    status = fixed_field(tle->eccentricity, 7, 9, ' ', eccentricity);
  }
  if (status == MYB_OK && eccentricity[0] != '0') {
    status = MYB_E_DOMAIN; /* rounds to 1 */
  }
  if (status == MYB_OK) {
    memcpy(line2 + 26, eccentricity + 2, 7);
    status = fixed_field(tle->argument_of_perigee, 4, 8, ' ', line2 + 34);
  }
  if (status == MYB_OK) {
    status = fixed_field(tle->mean_anomaly, 4, 8, ' ', line2 + 43);
  }
  if (status == MYB_OK) {
    status = fixed_field(tle->mean_motion, 8, 11, ' ', line2 + 52);
  }
  if (status == MYB_OK) {
    status = integer_field(tle->revolution, 5, line2 + 63);
  }
  if (status != MYB_OK) {
    return status;
  }
  if (tle->inclination < 0.0 || tle->node < 0.0 || tle->argument_of_perigee < 0.0 ||
      tle->mean_anomaly < 0.0 || !(tle->mean_motion > 0.0)) {
    return MYB_E_DOMAIN;
  }
  line1[68] = checksum_of(line1);
  line2[68] = checksum_of(line2);
  return MYB_OK;
}

static int day_of_year(const myb_calendar* c) {
  static const int before[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  const int leap = (c->year % 4 == 0 && c->year % 100 != 0) || c->year % 400 == 0;
  return before[c->month - 1] + c->day + (leap && c->month > 2 ? 1 : 0);
}

myb_status myb_tle_from_omm(const myb_omm* omm, myb_tle* tle) {
  MYB_CHECK(omm != NULL && tle != NULL);
  if (!omm->has_mean_motion || !omm->has_tle || omm->epoch.month < 1 || omm->epoch.month > 12 ||
      omm->epoch.year < 1957 || omm->epoch.year > 2056) {
    return MYB_E_DOMAIN;
  }
  myb_tle t;
  memset(&t, 0, sizeof t);
  if (satnum_of(omm->norad_cat_id, t.satnum) != MYB_OK) {
    return MYB_E_DOMAIN;
  }
  t.catalog_number = omm->norad_cat_id;
  t.classification = omm->classification_type != '\0' ? omm->classification_type : 'U';
  /* the international designator YYYY-NNNP[PP] as YYNNNP[PP] */
  const char* id = omm->object_id;
  const size_t n = strlen(id);
  if (n >= 9 && n <= 11 && id[4] == '-') {
    int ok = 1;
    for (size_t i = 0; i < 8 && ok; ++i) {
      ok = i == 4 || (id[i] >= '0' && id[i] <= '9');
    }
    if (ok) {
      memcpy(t.designator, id + 2, 2);
      memcpy(t.designator + 2, id + 5, n - 5);
      t.designator[n - 3] = '\0';
    }
  }
  t.epoch_year = (int)omm->epoch.year;
  t.epoch_day = (double)day_of_year(&omm->epoch) +
                ((double)(omm->epoch.hour * 3600 + omm->epoch.minute * 60) + omm->epoch.second) /
                    86400.0;
  t.ndot = omm->mean_motion_dot;
  t.nddot = omm->mean_motion_ddot;
  t.bstar = omm->bstar;
  t.ephemeris_type = omm->ephemeris_type;
  t.element_number = (int)omm->element_set_no;
  t.inclination = omm->inclination;
  t.node = omm->ra_of_asc_node;
  t.eccentricity = omm->eccentricity;
  t.argument_of_perigee = omm->arg_of_pericenter;
  t.mean_anomaly = omm->mean_anomaly;
  t.mean_motion = omm->mean_motion;
  t.revolution = omm->rev_at_epoch;
  *tle = t;
  return MYB_OK;
}

myb_status myb_omm_from_tle(const myb_tle* tle, myb_omm* omm) {
  MYB_CHECK(tle != NULL && omm != NULL);
  memset(omm, 0, sizeof *omm);
  myb_tle_epoch(tle, &omm->epoch);
  omm->creation_date = omm->epoch;
  memcpy(omm->originator, "MYB", 4);
  memcpy(omm->center_name, "EARTH", 6);
  memcpy(omm->ref_frame, "TEME", 5);
  memcpy(omm->time_system, "UTC", 4);
  memcpy(omm->mean_element_theory, "SGP4", 5);
  const size_t satnum = strlen(tle->satnum);
  memcpy(omm->object_name, tle->satnum, satnum + 1);
  /* YYNNNP[PP] as YYYY-NNNP[PP]; the catalogue number where there is no designator */
  const size_t n = strlen(tle->designator);
  int digits = n >= 6;
  for (size_t i = 0; i < 5 && digits; ++i) {
    digits = tle->designator[i] >= '0' && tle->designator[i] <= '9';
  }
  if (digits) {
    const int two = 10 * (tle->designator[0] - '0') + (tle->designator[1] - '0');
    const int year = two >= 57 ? 1900 + two : 2000 + two;
    omm->object_id[0] = (char)('0' + year / 1000);
    omm->object_id[1] = (char)('0' + (year / 100) % 10);
    omm->object_id[2] = tle->designator[0];
    omm->object_id[3] = tle->designator[1];
    omm->object_id[4] = '-';
    memcpy(omm->object_id + 5, tle->designator + 2, n - 2);
    omm->object_id[n + 3] = '\0';
  } else {
    memcpy(omm->object_id, tle->satnum, satnum + 1);
  }
  omm->has_mean_motion = 1;
  omm->mean_motion = tle->mean_motion;
  omm->eccentricity = tle->eccentricity;
  omm->inclination = tle->inclination;
  omm->ra_of_asc_node = tle->node;
  omm->arg_of_pericenter = tle->argument_of_perigee;
  omm->mean_anomaly = tle->mean_anomaly;
  omm->has_tle = 1;
  omm->ephemeris_type = tle->ephemeris_type;
  omm->classification_type = tle->classification;
  omm->norad_cat_id = tle->catalog_number;
  omm->element_set_no = tle->element_number;
  omm->rev_at_epoch = tle->revolution;
  omm->bstar = tle->bstar;
  omm->mean_motion_dot = tle->ndot;
  omm->mean_motion_ddot = tle->nddot;
  return MYB_OK;
}

/* ---- Writing SPK type 10 ------------------------------------------------------------------------
 * The generic segment that SPKW10 builds (NAIF SGPARAM): the constants; the packets, each
 * after one word that holds its epoch; the epochs once more as the reference values; every
 * hundredth of them as their directory; and the seventeen words of meta data that say where
 * all that is. */

enum { PACKET_WORDS = 14, DIRECTORY_STEP = 100, EXPLICIT_CLOSEST = 4 };

void myb_sgp4_spk10_constants(const myb_sgp4_constants* model, double constants[8]) {
  MYB_CHECK(model != NULL && constants != NULL);
  constants[0] = model->j2;
  constants[1] = model->j3;
  constants[2] = model->j4;
  constants[3] = model->xke;
  constants[4] = 120.0;
  constants[5] = 78.0;
  constants[6] = model->radius;
  constants[7] = 1.0;
}

myb_status myb_sgp4_spk10_packet(const myb_tle* tle, const myb_timesys* system,
                                 const myb_equinox* equinox, double packet[14]) {
  MYB_CHECK(tle != NULL && system != NULL && packet != NULL);
  myb_scale tdb;
  myb_status st = myb_timesys_scale(system, "TDB", &tdb);
  if (st != MYB_OK) {
    return st;
  }
  myb_calendar utc;
  myb_tle_epoch(tle, &utc);
  myb_time epoch;
  st = myb_timesys_from_utc(system, &utc, tdb, &epoch);
  if (st != MYB_OK) {
    return st;
  }
  const double minutes_per_day = 1440.0;
  const double rad = PI / 180.0;
  packet[0] = tle->ndot * (2.0 * PI) / (minutes_per_day * minutes_per_day);
  packet[1] = tle->nddot * (2.0 * PI) / (minutes_per_day * minutes_per_day * minutes_per_day);
  packet[2] = tle->bstar;
  packet[3] = tle->inclination * rad;
  packet[4] = tle->node * rad;
  packet[5] = tle->eccentricity;
  packet[6] = tle->argument_of_perigee * rad;
  packet[7] = tle->mean_anomaly * rad;
  packet[8] = tle->mean_motion * (2.0 * PI) / minutes_per_day;
  packet[9] = myb_time_seconds(epoch);
  packet[10] = packet[11] = packet[12] = packet[13] = 0.0;
  if (equinox != NULL) {
    double nutation[2];
    double rates[2];
    myb_equinox_nutation(equinox, epoch, nutation, rates);
    packet[10] = nutation[1]; /* in obliquity */
    packet[11] = nutation[0]; /* in longitude */
    packet[12] = rates[1];
    packet[13] = rates[0];
  }
  return MYB_OK;
}

size_t myb_sgp4_spk10_size(size_t count) {
  if (count == 0 || count > SIZE_MAX / 256) {
    return 0;
  }
  return GEOPHS + count * (PACKET_WORDS + 1) + count + (count - 1) / DIRECTORY_STEP + META_ITEMS;
}

myb_status myb_sgp4_spk10_words(const double constants[8], const double* packets, size_t count,
                                double* words, size_t size) {
  MYB_CHECK(constants != NULL && packets != NULL && words != NULL);
  const size_t total = myb_sgp4_spk10_size(count);
  if (total == 0) {
    return MYB_E_DOMAIN;
  }
  if (size < total) {
    return MYB_E_OVERFLOW;
  }
  for (size_t k = 0; k < GEOPHS; ++k) {
    if (!isfinite(constants[k])) {
      return MYB_E_DOMAIN;
    }
  }
  for (size_t k = 0; k < count; ++k) {
    for (size_t i = 0; i < PACKET_WORDS; ++i) {
      if (!isfinite(packets[PACKET_WORDS * k + i])) {
        return MYB_E_DOMAIN;
      }
    }
    if (k > 0 && !(packets[PACKET_WORDS * k + 9] > packets[PACKET_WORDS * (k - 1) + 9])) {
      return MYB_E_DOMAIN;
    }
  }
  double* w = words;
  memcpy(w, constants, GEOPHS * sizeof(double));
  w += GEOPHS;
  for (size_t k = 0; k < count; ++k) {
    *w++ = packets[PACKET_WORDS * k + 9];
    memcpy(w, packets + PACKET_WORDS * k, PACKET_WORDS * sizeof(double));
    w += PACKET_WORDS;
  }
  const size_t references = (size_t)(w - words);
  for (size_t k = 0; k < count; ++k) {
    *w++ = packets[PACKET_WORDS * k + 9];
  }
  const size_t directory = (size_t)(w - words);
  const size_t entries = (count - 1) / DIRECTORY_STEP;
  for (size_t k = 1; k <= entries; ++k) {
    *w++ = packets[PACKET_WORDS * (k * DIRECTORY_STEP - 1) + 9];
  }
  const double meta[META_ITEMS] = {
      [META_CONBAS] = 0.0,
      [META_NCON] = (double)GEOPHS,
      [2] = (double)directory, /* the directory of the reference values */
      [3] = (double)entries,
      [4] = (double)EXPLICIT_CLOSEST,
      [META_REFBAS] = (double)references,
      [META_NREF] = (double)count,
      [META_PKTBAS] = (double)GEOPHS,
      [META_NPKT] = (double)count,
      [META_PKTSZ] = (double)PACKET_WORDS,
      [META_PKTOFF] = 1.0,
      [16] = (double)META_ITEMS};
  memcpy(w, meta, sizeof meta);
  w += META_ITEMS;
  MYB_CHECK((size_t)(w - words) == total);
  return MYB_OK;
}
