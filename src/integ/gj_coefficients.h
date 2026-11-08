/* SPDX-License-Identifier: BSD-3-Clause */
/* Coefficients of the Gauss-Jackson method of order 8 in summed ordinate form (cf. M. M. Berry,
 * L. M. Healy, J. Astronaut. Sci. 52, 331, 2004), derived here from the operator identities
 *   h v - s = P1(D) h^2 a,   P1(D) = 1/D - 1/(e^D - 1) = 1/2 - D/12 + D^3/720 - ...,
 *   x - S   = P2(D) h^2 a,   P2(D) = 1/D^2 - e^D/(e^D - 1)^2 = 1/12 - D^2/240 + ...,
 * with the sums s(n+1) = s(n) + h^2 a(n), S(n+1) = S(n) + s(n+1) and D = h d/dt (Bernoulli
 * numbers), applied to the degree-8 interpolant of the accelerations at the window nodes
 * u = -4..4. Row j + 4 gives the correction at window position j (j = 5: one step beyond the
 * window, the predictor; j = 4: the newest node, the corrector). The predictor row reproduces
 * the published one (3250433/53222400, ..., 103798439/159667200). Exact rationals rounded to
 * nearest double once; tests/mpfr/test_mpfr_gauss_jackson.c recomputes every entry. */
#ifndef MYB_INTEG_GJ_COEFFICIENTS_H_
#define MYB_INTEG_GJ_COEFFICIENTS_H_

// clang-format off
/* L[k][m]: coefficient of u^m in the Lagrange basis polynomial of node k - 4 */
static const double myb_gj_lagrange[9][9] = {
    {0x0.0p+0, 0x1.d41d41d41d41dp-9, -0x1.d41d41d41d41dp-11, -0x1.3e93e93e93e94p-8, 0x1.3e93e93e93e94p-10, 0x1.6c16c16c16c17p-10, -0x1.6c16c16c16c17p-12, -0x1.a01a01a01a01ap-14, 0x1.a01a01a01a01ap-16},
    {0x0.0p+0, -0x1.3813813813814p-5, 0x1.a01a01a01a01ap-7, 0x1.999999999999ap-5, -0x1.1111111111111p-6, -0x1.999999999999ap-7, 0x1.1111111111111p-8, 0x1.3813813813814p-11, -0x1.a01a01a01a01ap-13},
    {0x0.0p+0, 0x1.999999999999ap-3, -0x1.999999999999ap-4, -0x1.e0b60b60b60b6p-3, 0x1.e0b60b60b60b6p-4, 0x1.27d27d27d27d2p-5, -0x1.27d27d27d27d2p-6, -0x1.6c16c16c16c17p-10, 0x1.6c16c16c16c17p-11},
    {0x0.0p+0, -0x1.999999999999ap-1, 0x1.999999999999ap-1, 0x1.5b05b05b05b06p-2, -0x1.5b05b05b05b06p-2, -0x1.49f49f49f49f5p-5, 0x1.49f49f49f49f5p-5, 0x1.6c16c16c16c17p-10, -0x1.6c16c16c16c17p-10},
    {0x1.0000000000000p+0, 0x0.0p+0, -0x1.6c71c71c71c72p+0, 0x0.0p+0, 0x1.e555555555555p-2, 0x0.0p+0, -0x1.aaaaaaaaaaaabp-5, 0x0.0p+0, 0x1.c71c71c71c71cp-10},
    {0x0.0p+0, 0x1.999999999999ap-1, 0x1.999999999999ap-1, -0x1.5b05b05b05b06p-2, -0x1.5b05b05b05b06p-2, 0x1.49f49f49f49f5p-5, 0x1.49f49f49f49f5p-5, -0x1.6c16c16c16c17p-10, -0x1.6c16c16c16c17p-10},
    {0x0.0p+0, -0x1.999999999999ap-3, -0x1.999999999999ap-4, 0x1.e0b60b60b60b6p-3, 0x1.e0b60b60b60b6p-4, -0x1.27d27d27d27d2p-5, -0x1.27d27d27d27d2p-6, 0x1.6c16c16c16c17p-10, 0x1.6c16c16c16c17p-11},
    {0x0.0p+0, 0x1.3813813813814p-5, 0x1.a01a01a01a01ap-7, -0x1.999999999999ap-5, -0x1.1111111111111p-6, 0x1.999999999999ap-7, 0x1.1111111111111p-8, -0x1.3813813813814p-11, -0x1.a01a01a01a01ap-13},
    {0x0.0p+0, -0x1.d41d41d41d41dp-9, -0x1.d41d41d41d41dp-11, 0x1.3e93e93e93e94p-8, 0x1.3e93e93e93e94p-10, -0x1.6c16c16c16c17p-10, -0x1.6c16c16c16c17p-12, 0x1.a01a01a01a01ap-14, 0x1.a01a01a01a01ap-16},
};
/* V[j + 4][k]: velocity correction at window position j = -4..5 for node k - 4 */
static const double myb_gj_velocity[10][9] = {
    {0x1.6d118de5ab278p-1, -0x1.2d94005c778e9p-1, 0x1.ed9359ebaa1b5p-1, -0x1.3daafde311645p+0, 0x1.23fc06dcdf951p+0, -0x1.711f8cfe80431p-1, 0x1.3100d00d00d01p-2, -0x1.28f294d9a8a04p-4, 0x1.029f8aaeb64c2p-7},
    {0x1.029f8aaeb64c2p-7, 0x1.48b31e651984dp-1, -0x1.383484b462079p-2, 0x1.3443e7ccb5aefp-2, -0x1.f8af8af8af8b0p-3, 0x1.2b37d2c721116p-3, -0x1.dae26f930ff27p-5, 0x1.c1aa810e774dep-7, -0x1.7fc63546e2a71p-10},
    {-0x1.7fc63546e2a71p-10, 0x1.592f434f3aa40p-6, 0x1.2db72ea61d951p-1, -0x1.748ef6723f517p-3, 0x1.dd81664f4882cp-4, -0x1.fba1ba1ba1ba2p-5, 0x1.7aedfe84e29ddp-6, -0x1.591b9d1a87b3ep-8, 0x1.1eb851eb851ecp-11},
    {0x1.1eb851eb851ecp-11, -0x1.a280e97aae6c5p-8, 0x1.4ddf4fbc18335p-5, 0x1.16320fedcba98p-1, -0x1.cee07c40a798bp-4, 0x1.8687eb56e2f10p-5, -0x1.069f9d2d06039p-6, 0x1.c2f512df66da0p-9, -0x1.68c40f191f14dp-12},
    {-0x1.68c40f191f14dp-12, 0x1.dd8aa577243f1p-9, -0x1.338e82dccd16cp-6, 0x1.1d4ffcd24a4c7p-4, 0x1.0000000000000p-1, -0x1.1d4ffcd24a4c7p-4, 0x1.338e82dccd16cp-6, -0x1.dd8aa577243f1p-9, 0x1.68c40f191f14dp-12},
    {0x1.68c40f191f14dp-12, -0x1.c2f512df66da0p-9, 0x1.069f9d2d06039p-6, -0x1.8687eb56e2f10p-5, 0x1.cee07c40a798bp-4, 0x1.d39be02468acfp-2, -0x1.4ddf4fbc18335p-5, 0x1.a280e97aae6c5p-8, -0x1.1eb851eb851ecp-11},
    {-0x1.1eb851eb851ecp-11, 0x1.591b9d1a87b3ep-8, -0x1.7aedfe84e29ddp-6, 0x1.fba1ba1ba1ba2p-5, -0x1.dd81664f4882cp-4, 0x1.748ef6723f517p-3, 0x1.a491a2b3c4d5ep-2, -0x1.592f434f3aa40p-6, 0x1.7fc63546e2a71p-10},
    {0x1.7fc63546e2a71p-10, -0x1.c1aa810e774dep-7, 0x1.dae26f930ff27p-5, -0x1.2b37d2c721116p-3, 0x1.f8af8af8af8b0p-3, -0x1.3443e7ccb5aefp-2, 0x1.383484b462079p-2, 0x1.6e99c335ccf67p-2, -0x1.029f8aaeb64c2p-7},
    {-0x1.029f8aaeb64c2p-7, 0x1.28f294d9a8a04p-4, -0x1.3100d00d00d01p-2, 0x1.711f8cfe80431p-1, -0x1.23fc06dcdf951p+0, 0x1.3daafde311645p+0, -0x1.ed9359ebaa1b5p-1, 0x1.2d94005c778e9p-1, 0x1.25dce434a9b10p-2},
    {0x1.25dce434a9b10p-2, -0x1.4b9b2045ed9d7p+1, 0x1.4cea65e4f2386p+3, -0x1.8675eec552bb9p+4, 0x1.2709eed7d10b5p+5, -0x1.2a6550dabe06ep+5, 0x1.958c9b634fce9p+4, -0x1.6971b659f988dp+3, 0x1.95fd80d25ccacp+1},
};
/* X[j + 4][k]: position correction at window position j = -4..5 for node k - 4 */
static const double myb_gj_position[10][9] = {
    {0x1.f44ea13ac1712p-5, 0x1.9b657df8b6b02p-4, -0x1.be7466c206c2dp-3, 0x1.35dd7e09869b8p-2, -0x1.26106b0bdcb37p-2, 0x1.7a2a32ab710cfp-3, -0x1.3bc6b0017f516p-4, 0x1.359e34225a744p-6, -0x1.0f07422ab2714p-9},
    {-0x1.0f07422ab2714p-9, 0x1.46615b3962e87p-4, 0x1.a9f54e22b7c30p-6, -0x1.6aeb7fa7f238dp-5, 0x1.5892c73bf7a07p-5, -0x1.b4545e9d50bfdp-6, 0x1.670abd366d84dp-7, -0x1.5bd0b22ed23f6p-9, 0x1.2d826c946d386p-12},
    {0x1.2d826c946d386p-12, -0x1.311cfe28d6a85p-8, 0x1.70c7b27e42445p-4, 0x1.e3a1f9fe86910p-10, -0x1.087c73175ad56p-7, 0x1.7e3322cee0e7ep-8, -0x1.44c980d40b2f2p-9, 0x1.3d8030f72a574p-11, -0x1.13bf00faebff7p-14},
    {-0x1.13bf00faebff7p-14, 0x1.ccf81764801b8p-11, -0x1.cc386eb5fb680p-8, 0x1.87665d92d7a05p-4, -0x1.a5f78b6e5efaap-8, 0x1.bce477e95e856p-13, 0x1.44871858b289ap-12, -0x1.d253f73835f84p-14, 0x1.d253f73835f84p-17},
    {0x1.d253f73835f84p-17, -0x1.902ebb8d145b6p-13, 0x1.69a3a93a0f3b9p-10, -0x1.0c5d1aa29a216p-7, 0x1.8e92883044f4dp-4, -0x1.0c5d1aa29a216p-7, 0x1.69a3a93a0f3b9p-10, -0x1.902ebb8d145b6p-13, 0x1.d253f73835f84p-17},
    {0x1.d253f73835f84p-17, -0x1.d253f73835f84p-14, 0x1.44871858b289ap-12, 0x1.bce477e95e856p-13, -0x1.a5f78b6e5efaap-8, 0x1.87665d92d7a05p-4, -0x1.cc386eb5fb680p-8, 0x1.ccf81764801b8p-11, -0x1.13bf00faebff7p-14},
    {-0x1.13bf00faebff7p-14, 0x1.3d8030f72a574p-11, -0x1.44c980d40b2f2p-9, 0x1.7e3322cee0e7ep-8, -0x1.087c73175ad56p-7, 0x1.e3a1f9fe86910p-10, 0x1.70c7b27e42445p-4, -0x1.311cfe28d6a85p-8, 0x1.2d826c946d386p-12},
    {0x1.2d826c946d386p-12, -0x1.5bd0b22ed23f6p-9, 0x1.670abd366d84dp-7, -0x1.b4545e9d50bfdp-6, 0x1.5892c73bf7a07p-5, -0x1.6aeb7fa7f238dp-5, 0x1.a9f54e22b7c30p-6, 0x1.46615b3962e87p-4, -0x1.0f07422ab2714p-9},
    {-0x1.0f07422ab2714p-9, 0x1.359e34225a744p-6, -0x1.3bc6b0017f516p-4, 0x1.7a2a32ab710cfp-3, -0x1.26106b0bdcb37p-2, 0x1.35dd7e09869b8p-2, -0x1.be7466c206c2dp-3, 0x1.9b657df8b6b02p-4, 0x1.f44ea13ac1712p-5},
    {0x1.f44ea13ac1712p-5, -0x1.1a7b41f337821p-1, 0x1.1bd7771951849p+1, -0x1.4d42b48e94ef8p+2, 0x1.f84eb84b31f3cp+2, -0x1.fede6d6694369p+2, 0x1.5bb171af275bfp+2, -0x1.3553811d2d3bdp+1, 0x1.4cd8ea7023a5ap-1},
};

// clang-format on

#endif /* MYB_INTEG_GJ_COEFFICIENTS_H_ */
