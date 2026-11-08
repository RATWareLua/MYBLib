/* SPDX-License-Identifier: BSD-3-Clause */
/* Weights of the symmetric compositions of the leapfrog (H. Yoshida, Phys. Lett. A 150, 262,
 * 1990, the "triple jump"): S_(p+2)(h) = S_p(z1 h) S_p(z0 h) S_p(z1 h) with
 * z1 = 1 / (2 - 2^(1/(p+1))), z0 = 1 - 2 z1, for p = 2, 4, 6. Each weight is a product of such
 * factors, computed with MPFR at 512 bits and rounded to nearest once
 * (tests/mpfr/test_mpfr_yoshida.c recomputes every entry and compares bit for bit). */
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

// clang-format on

#endif /* MYB_INTEG_SYMPLECTIC_COEFFICIENTS_H_ */
