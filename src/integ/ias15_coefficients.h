/* SPDX-License-Identifier: BSD-3-Clause */
/* Gauss-Radau constants of IAS15 (H. Rein, D. S. Spiegel, MNRAS 446, 1424, 2015; E. Everhart,
 * in Dynamics of Comets, 1985). Computed with MPFR at 512 bits and
 * rounded to nearest double once (tests/mpfr/test_mpfr_radau.c recomputes every entry and
 * compares bit for bit); not copied from any implementation.
 *
 * Nodes h_0 = 0 < h_1 < ... < h_7 < 1: h_0 and the roots of P_7(2h - 1) + P_8(2h - 1) mapped
 * to [0, 1] (Radau quadrature with the left end fixed, exact to degree 14).
 * Within a step the acceleration is the polynomial
 *   a(h) = a_0 + b_0 h + b_1 h^2 + ... + b_6 h^7
 *        = a_0 + g_0 h + g_1 h (h - h_1) + ... + g_6 h (h - h_1) ... (h - h_6),
 * the g being divided differences of the accelerations at the nodes. Then
 *   b_k = sum_j C[k][j] g_j   (C[k][j]: coefficient of h^(k+1) in h (h - h_1) ... (h - h_j)),
 *   g_j = sum_k D[j][k] b_k   (D = C^-1),
 *   R[n][k] = 1 / (h_n - h_k), k < n. */
#ifndef MYB_INTEG_IAS15_COEFFICIENTS_H_
#define MYB_INTEG_IAS15_COEFFICIENTS_H_

// clang-format off
static const double myb_ias15_h[8] = {0x0p+0, 0x1.cce7242fd9812p-5, 0x1.7122082358d27p-3, 0x1.6916742f0c9c5p-2, 0x1.18248525f4802p-1, 0x1.77ea65770fa49p-1, 0x1.c548c982467cdp-1, 0x1.f47d94f4d2bb2p-1};
/* C[k][j]: b_k = sum_j C[k][j] g_j */
static const double myb_ias15_c[7][7] = {
    {0x1p+0, -0x1.cce7242fd9812p-5, 0x1.4c4b384151c1ep-7, -0x1.d4b3379a1992dp-9, 0x1.007370e57b89dp-9, -0x1.7893e98a2a997p-10, 0x1.4d64672378b3dp-10},
    {0x0p+0, 0x1p+0, -0x1.e45bd12f4f32cp-3, 0x1.7f2161e152586p-4, -0x1.c08e5e7efcf23p-5, 0x1.595cd739ee74fp-5, -0x1.3d865cb429919p-5},
    {0x0p+0, 0x0p+0, 0x1p+0, -0x1.2da22e635a1aep-1, 0x1.a9dcc2e57fbbfp-2, -0x1.70bdf2d91dd44p-2, 0x1.71a0161ae20fep-2},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, -0x1.22e359c4a74d8p+0, 0x1.4009e0866d488p+0, -0x1.7785b93619434p+0},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, -0x1.ded88c802f1fcp+0, 0x1.73fc45dfb9889p+1},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, -0x1.60be78a0a92f1p+1},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0},
};
/* D = C^-1: g_j = sum_k D[j][k] b_k */
static const double myb_ias15_d[7][7] = {
    {0x1p+0, 0x1.cce7242fd9812p-5, 0x1.9ee7bd5f249p-9, 0x1.757f83047d1ecp-13, 0x1.50392d238405bp-17, 0x1.2eab281d37645p-21, 0x1.107669a72b0ebp-25},
    {0x0p+0, 0x1p+0, 0x1.e45bd12f4f32cp-3, 0x1.7722d9e44e93cp-5, 0x1.144b8b5c061cdp-7, 0x1.9106110af09c6p-10, 0x1.21b6dca35a7eap-12},
    {0x0p+0, 0x0p+0, 0x1p+0, 0x1.2da22e635a1aep-1, 0x1.039e6f5c0b72fp-2, 0x1.90baa1e56027dp-4, 0x1.27255aa2a25aap-5},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x1.22e359c4a74d8p+0, 0x1.c0216120faaf2p-1, 0x1.27498c3d6d93p-1},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x1.ded88c802f1fcp+0, 0x1.1fd1c1aad40f1p+1},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x1.60be78a0a92f1p+1},
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0},
};
/* R[n][k] = 1 / (h_n - h_k), k < n */
static const double myb_ias15_r[8][8] = {
    {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
    {0x1.1c618574b517bp+4, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
    {0x1.6314ac07f6a9ap+2, 0x1.021c2b5e4af6fp+3, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
    {0x1.6afdfcc2a4aefp+1, 0x1.afe76c59caa8ap+1, 0x1.73439c0b4562bp+2, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
    {0x1.d3e03b89d73cfp+0, 0x1.04c014a3946fbp+1, 0x1.5cdb4a576d57ap+1, 0x1.48fffc3fdf2e6p+2, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
    {0x1.5cac8b52e27b1p+0, 0x1.799c3c56b4e9dp+0, 0x1.ce1e8b86f5ae8p+0, 0x1.4f714afea15ccp+1, 0x1.56247c4b878b4p+2, 0x0p+0, 0x0p+0, 0x0p+0},
    {0x1.212921ce4baffp+0, 0x1.34c8b7028f418p+0, 0x1.6b1448c75e4d5p+0, 0x1.e092f6d7c8a42p+0, 0x1.7a82c7145ec17p+1, 0x1.a787c6413556bp+2, 0x0p+0, 0x0p+0},
    {0x1.05e3166485188p+0, 0x1.15e18175af101p+0, 0x1.41177c7d3aff5p+0, 0x1.99ab118ef2216p+0, 0x1.296bab4a3e0d2p+1, 0x1.0709d7d8d5e3cp+2, 0x1.5b12a584540e2p+3, 0x0p+0},
};

// clang-format on

#endif /* MYB_INTEG_IAS15_COEFFICIENTS_H_ */
