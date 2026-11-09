/* SPDX-License-Identifier: BSD-3-Clause */
/* dyn/wh_coefficients.h — coefficients of the Wisdom-Holman splitting schemes and of the
 * symplectic correctors, as the doubles nearest to the exact values (checked against MPFR by
 * tests/mpfr/test_mpfr_wh.c).
 *
 * Schemes (fractions of the step h, in the order of application):
 *   SABAn (Laskar & Robutel 2001): drifts c1, c2, ..., c1 around n kicks at the Gauss-Legendre
 *     nodes of [0, 1], the kicks weighted by the Gauss weights:
 *       SABA2: c1 = 1/2 - sqrt3/6, c2 = sqrt3/3; d = 1/2
 *       SABA3: c1 = 1/2 - sqrt15/10, c2 = sqrt15/10; d1 = 5/18, d2 = 4/9
 *       SABA4: c1 = 1/2 - u/70, c2 = (u - v)/70, c3 = v/35, u, v = sqrt(525 +- 70 sqrt30);
 *              d1 = 1/4 - sqrt30/72, d2 = 1/4 + sqrt30/72
 *   SBABn: kicks at the n+1 Gauss-Lobatto nodes, drifts between them:
 *       SBAB1: d = 1/2; c = 1                       SBAB2: d1 = 1/6, d2 = 2/3; c = 1/2
 *       SBAB3: c2 = 1/2 - sqrt5/10, c3 = sqrt5/5; d1 = 1/12, d2 = 5/12
 *       SBAB4: c2 = 1/2 - sqrt21/14, c3 = sqrt21/14; d1 = 1/20, d2 = 49/180, d3 = 16/45
 *
 * Symplectic correctors of the DKD map (Wisdom, Holman & Touma 1996), derived here. To first
 * order in the perturbation B, the DKD map is exp(h A + phi(h ad_A) h B) with
 * phi(t) = (t/2) / sinh(t/2), so the corrector exp(chi) that conjugates it to the exact flow
 * up to O(e h^(2m+2)) has chi = -h ((1 - phi(t)) / t)(t = h ad_A) B, an odd series
 * sum_n gamma_(2n-1) t^(2n-1), gamma_(2n-1) = (2^(2n) - 2) B_(2n) / (2^(2n) (2n)!) (Bernoulli
 * numbers). It is built from Z(a, b) = A(a) B(b) A(-2a) B(-b) A(a), whose first-order
 * generator is 2 b h sinh(a h ad_A) B, as
 *   C = Z(a_m, b_m) ... Z(a_1, b_1) Z(-a_1, -b_1) ... Z(-a_m, -b_m),   C^-1 = C with -b,
 * with a_i = i alpha, alpha = sqrt(7/40) (the choice of Wisdom et al.), and b_i solving
 *   sum_i b_i a_i^k = k! gamma_k / 4 = (2^(2n) - 2) B_(2n) / (4 2^(2n) 2n),  k = 2n - 1,
 * for n = 1 .. m: the corrector of order 2m + 1. For m = 1, b_1 = 1/(96 alpha) =
 * 0.024900596027799867..., the value published by Wisdom et al. */
#ifndef MYB_SRC_DYN_WH_COEFFICIENTS_H_
#define MYB_SRC_DYN_WH_COEFFICIENTS_H_

enum { MYB_WH_SCHEMES = 8, MYB_WH_MAX_KICKS = 5, MYB_WH_CORRECTOR_STAGES = 8 };

typedef struct myb_wh_table {
  int kick_first; /* 0: A B A ... B A (DKD, SABA); 1: B A B ... A B (SBAB) */
  int kicks;      /* kicks per step; drifts = kicks + 1 (SABA) or kicks - 1 (SBAB) */
  double drift[MYB_WH_MAX_KICKS];
  double kick[MYB_WH_MAX_KICKS];
} myb_wh_table;

/* Indexed by myb_wh_scheme. */
static const myb_wh_table myb_wh_tables[MYB_WH_SCHEMES] = {
    /* MYB_WH_DKD */
    {.kick_first = 0, .kicks = 1, .drift = {0.5, 0.5}, .kick = {1.0}},
    /* MYB_WH_SABA2 */
    {.kick_first = 0,
     .kicks = 2,
     .drift = {0x1.b0cb174df99c7p-3, 0x1.279a74590331cp-1, 0x1.b0cb174df99c7p-3},
     .kick = {0.5, 0.5}},
    /* MYB_WH_SABA3 */
    {.kick_first = 0,
     .kicks = 3,
     .drift = {0x1.cda042f0236e1p-4, 0x1.8c97ef43f7248p-2, 0x1.8c97ef43f7248p-2,
               0x1.cda042f0236e1p-4},
     .kick = {0x1.1c71c71c71c72p-2, 0x1.c71c71c71c71cp-2, 0x1.1c71c71c71c72p-2}},
    /* MYB_WH_SABA4 */
    {.kick_first = 0,
     .kicks = 4,
     .drift = {0x1.1c6490c2719ecp-4, 0x1.0ad4dd0079a87p-2, 0x1.5c23fd9dd3dfcp-2,
               0x1.0ad4dd0079a87p-2, 0x1.1c6490c2719ecp-4},
     .kick = {0x1.64340f7e7b66bp-3, 0x1.4de5f840c24cap-2, 0x1.4de5f840c24cap-2,
              0x1.64340f7e7b66bp-3}},
    /* MYB_WH_SBAB1 */
    {.kick_first = 1, .kicks = 2, .drift = {1.0}, .kick = {0.5, 0.5}},
    /* MYB_WH_SBAB2 */
    {.kick_first = 1,
     .kicks = 3,
     .drift = {0.5, 0.5},
     .kick = {0x1.5555555555555p-3, 0x1.5555555555555p-1, 0x1.5555555555555p-3}},
    /* MYB_WH_SBAB3 */
    {.kick_first = 1,
     .kicks = 4,
     .drift = {0x1.1b06d1d200913p-2, 0x1.c9f25c5bfedd9p-2, 0x1.1b06d1d200913p-2},
     .kick = {0x1.5555555555555p-4, 0x1.aaaaaaaaaaaabp-2, 0x1.aaaaaaaaaaaabp-2,
              0x1.5555555555555p-4}},
    /* MYB_WH_SBAB4 */
    {.kick_first = 1,
     .kicks = 5,
     .drift = {0x1.61a277d8695abp-3, 0x1.4f2ec413cb52bp-2, 0x1.4f2ec413cb52bp-2,
               0x1.61a277d8695abp-3},
     .kick = {0x1.999999999999ap-5, 0x1.16c16c16c16c1p-2, 0x1.6c16c16c16c17p-2,
              0x1.16c16c16c16c1p-2, 0x1.999999999999ap-5}},
};

/* a_i = i sqrt(7/40), i = 1 .. 8. */
static const double myb_wh_corrector_a[MYB_WH_CORRECTOR_STAGES] = {
    0x1.ac5eb3f7ab2f8p-2, 0x1.ac5eb3f7ab2f8p-1, 0x1.414706f9c063ap+0, 0x1.ac5eb3f7ab2f8p+0,
    0x1.0bbb307acafdbp+1, 0x1.414706f9c063ap+1, 0x1.76d2dd78b5c99p+1, 0x1.ac5eb3f7ab2f8p+1};

/* Row m - 1: b_1 .. b_m of the corrector of order 2m + 1. */
static const double myb_wh_corrector_b[MYB_WH_CORRECTOR_STAGES][MYB_WH_CORRECTOR_STAGES] = {
    {0x1.97f8ab65c79afp-6},
    {0x1.53f9e42a26567p-5, -0x1.0ffb1cee8511fp-7},
    {0x1.ba138c281ad38p-5, -0x1.2b59ce73fceaap-6, 0x1.46b87ff974c34p-9},
    {0x1.065069a8159b4p-4, -0x1.d0745cc41db0cp-6, 0x1.be6533f3cdfe7p-8, -0x1.796145496f7bcp-11},
    {0x1.2957b105a15f8p-4, -0x1.384ad181052b2p-5, 0x1.935808ff91397p-7, -0x1.33da044f97665p-9,
     0x1.ab0365fa770ecp-13},
    {0x1.4817a7d2ba70bp-4, -0x1.852aba81c3d62p-5, 0x1.302be5d61c2b5p-6, -0x1.3decd11851642p-8,
     0x1.94ef0bee1c79ep-11, -0x1.dd16b718cabd1p-15},
    {0x1.63c49c0bd971cp-4, -0x1.cef7f06f6bd91p-5, 0x1.9edfb6ba982fcp-6, -0x1.0a4f712ba8b66p-7,
     0x1.d6d61b85ce47bp-10, -0x1.01d10e30f1fd2p-12, 0x1.083d9f127631fp-16},
    {0x1.7d2865a643682p-4, -0x1.0b07ac0fe3df1p-4, 0x1.099a47793a306p-5, -0x1.8b9144f7f2b3fp-7,
     0x1.b2467afd31964p-9, -0x1.4d7273c9751c8p-11, 0x1.4098e956e85c7p-14,
     -0x1.22e64af0557ffp-18},
};

#endif /* MYB_SRC_DYN_WH_COEFFICIENTS_H_ */
