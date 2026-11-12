/* SPDX-License-Identifier: BSD-3-Clause */
/* Weights of the symmetric compositions of the leapfrog (H. Yoshida, Phys. Lett. A 150, 262,
 * 1990, the "triple jump"): S_(p+2)(h) = S_p(z1 h) S_p(z0 h) S_p(z1 h) with
 * z1 = 1 / (2 - 2^(1/(p+1))), z0 = 1 - 2 z1, for p = 2, 4, 6. Each weight is a product of such
 * factors, computed with MPFR at 512 bits and rounded to nearest once
 * (tests/mpfr/test_mpfr_yoshida.c recomputes every entry and compares bit for bit).
 *
 * The compositions of 7 stages, order 6 (H. Yoshida 1990, solution A) and of 15 stages, order 8
 * (E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration, 2nd ed., 2006, V.3, the
 * method p8 s15, after W. Kahan, R.-C. Li, Math. Comp. 66, 1089, 1997, and R. I. McLachlan,
 * SIAM J. Sci. Comput. 16, 151, 1995): the published values refined by Newton's method on the
 * order conditions in MPFR and rounded to nearest once (tests/mpfr/test_mpfr_compositions.c
 * derives the conditions from the formula of Baker, Campbell and Hausdorff, refines the
 * weights at 512 bits from these doubles and compares bit for bit). */
#ifndef MYB_INTEG_SYMPLECTIC_COEFFICIENTS_H_
#define MYB_INTEG_SYMPLECTIC_COEFFICIENTS_H_

// clang-format off
static const double myb_symplectic_order2[1] = {1.0};
static const double myb_symplectic_order4[3] = {
    0x1.59e8b6eb96338p+0, -0x1.b3d16dd72c671p+0, 0x1.59e8b6eb96338p+0,
};
static const double myb_symplectic_order6[9] = {
    0x1.96545f73f687p+0, -0x1.fff1751765b91p+0, 0x1.96545f73f687p+0,
    -0x1.d2c007fc56da7p+0, 0x1.2608be2bcf859p+1, -0x1.d2c007fc56da7p+0,
    0x1.96545f73f687p+0, -0x1.fff1751765b91p+0, 0x1.96545f73f687p+0,
};
static const double myb_symplectic_order8[27] = {
    0x1.c589c3f9eaa63p+0, -0x1.1db60c8b87e9cp+1, 0x1.c589c3f9eaa63p+0,
    -0x1.047d403d5b813p+1, 0x1.483226c05243dp+1, -0x1.047d403d5b813p+1,
    0x1.c589c3f9eaa63p+0, -0x1.1db60c8b87e9cp+1, 0x1.c589c3f9eaa63p+0,
    -0x1.f4bf287fdec56p+0, 0x1.3b735e8b5cf7p+1, -0x1.f4bf287fdec56p+0,
    0x1.1f9a7c7c8b953p+1, -0x1.6a5b8f54d5022p+1, 0x1.1f9a7c7c8b953p+1,
    -0x1.f4bf287fdec56p+0, 0x1.3b735e8b5cf7p+1, -0x1.f4bf287fdec56p+0,
    0x1.c589c3f9eaa63p+0, -0x1.1db60c8b87e9cp+1, 0x1.c589c3f9eaa63p+0,
    -0x1.047d403d5b813p+1, 0x1.483226c05243dp+1, -0x1.047d403d5b813p+1,
    0x1.c589c3f9eaa63p+0, -0x1.1db60c8b87e9cp+1, 0x1.c589c3f9eaa63p+0,
};

static const double myb_symplectic_order6_s7[7] = {
    0x1.91abc49889363p-1, 0x1.e2743579895ddp-3, -0x1.2d7c6f7933b97p+0, 0x1.50b00cfb7be54p+0,
    -0x1.2d7c6f7933b97p+0, 0x1.e2743579895ddp-3, 0x1.91abc49889363p-1,
};
static const double myb_symplectic_order8_s15[15] = {
    0x1.7bbc37ce98b91p-1, -0x1.a2eb53ae4f982p-2, 0x1.86aa67d23e0d1p-3, -0x1.25d14d43b0cd8p-1,
    0x1.323de1786ab7dp-2, 0x1.56a7ea212931ep-2, 0x1.42dc3141feca6p-2, -0x1.9801b22c31f99p-1,
    0x1.42dc3141feca6p-2, 0x1.56a7ea212931ep-2, 0x1.323de1786ab7dp-2, -0x1.25d14d43b0cd8p-1,
    0x1.86aa67d23e0d1p-3, -0x1.a2eb53ae4f982p-2, 0x1.7bbc37ce98b91p-1,
};

// clang-format on

#endif /* MYB_INTEG_SYMPLECTIC_COEFFICIENTS_H_ */
