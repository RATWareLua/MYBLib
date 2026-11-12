/* SPDX-License-Identifier: BSD-3-Clause */
/* jb2008_internal.h — the static model of the module alone, for the tests (not exported). */
#ifndef MYB_EXT_JB2008_INTERNAL_H_
#define MYB_EXT_JB2008_INTERNAL_H_

/* The static model of Jacchia (1970) at the exospheric temperature tinf (K) and the height z
 * (km, >= 90): the density (g cm^-3), the temperature there and the number densities (cm^-3)
 * of N2, O2, O, Ar, He, H above the homopause (zeros below it). */
double myb_internal_jb2008_static(double tinf, double z, double* temperature, double number[6]);

#endif /* MYB_EXT_JB2008_INTERNAL_H_ */
