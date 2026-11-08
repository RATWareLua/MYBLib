/* SPDX-License-Identifier: BSD-3-Clause */
/* Coefficients of the Runge-Kutta-Nystrom pair RKN12(10) of J. R. Dormand, M. E. A. El-Mikkawy
 * and P. J. Prince (IMA J. Numer. Anal. 7, 423, 1987), 17 stages, for x'' = a(t, x).
 * Source of the numbers: the exact fractions of DPRKN12Tableau in OrdinaryDiffEq.jl
 * (lib/OrdinaryDiffEqRKN/src/rkn_tableaus.jl, MIT licence; see THIRD_PARTY_NOTICES.md),
 * rounded to nearest double once; the fractions are kept in tests/data/rkn12/dprkn12.txt and
 * tests/mpfr/test_mpfr_rkn12.c compares every entry bit for bit and re-checks the order
 * conditions in exact arithmetic. Rows 14-17 of A are 50-digit values (row sums equal c^2/2 to
 * 1e-53); everything else is exact.
 *
 * Stage i (0-based, c[0] = 0): X_i = x + h (c_i v + h sum_j a_ij k_j), k_i = a(t + c_i h, X_i);
 *   x_new = x + h v + h^2 sum b_i k_i,   v_new = v + h sum bp_i k_i      (order 12),
 *   error: h^2 sum bt_i k_i and h sum bpt_i k_i (difference to the order-10 solution). */
#ifndef MYB_INTEG_RKN12_COEFFICIENTS_H_
#define MYB_INTEG_RKN12_COEFFICIENTS_H_

enum { MYB_RKN12_STAGES = 17 };

// clang-format off
static const double myb_rkn12_c[MYB_RKN12_STAGES] = {0x0.0p+0, 0x1.47ae147ae147bp-6, 0x1.47ae147ae147bp-5, 0x1.999999999999ap-4, 0x1.1111111111111p-3, 0x1.47ae147ae147bp-3, 0x1.999999999999ap-5, 0x1.999999999999ap-3, 0x1.0000000000000p-2, 0x1.5555555555555p-2, 0x1.0000000000000p-1, 0x1.1c71c71c71c72p-1, 0x1.8000000000000p-1, 0x1.b6db6db6db6dbp-1, 0x1.e3f361749d2afp-1, 0x1.0000000000000p+0, 0x1.0000000000000p+0};
static const double myb_rkn12_a[17][17] = {
    {0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.a36e2eb1c432dp-13, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.179ec9cbd821ep-12, 0x1.179ec9cbd821ep-11, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.7e4b17e4b17e5p-9, -0x1.1111111111111p-8, 0x1.999999999999ap-8, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.af83440e53dbcp-10, 0x0.0p+0, 0x1.67980e0bf08c7p-8, 0x1.cc47c00f4861ep-10, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.fe07017c01026p-10, 0x0.0p+0, 0x1.d4b26350ce3dap-8, 0x1.7dc65dee09bcfp-9, 0x1.9e2855e5a8473p-11, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.28f5c28f5c28fp-11, 0x0.0p+0, 0x1.cde233788cde2p-11, -0x1.ca252adb363bfp-12, 0x1.63796ac9dfd13p-12, -0x1.a12f684bda12fp-14, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.94237fa89e60fp-9, 0x0.0p+0, 0x0.0p+0, 0x1.d208a5a912e32p-10, 0x1.61e4f765fd8aep-9, 0x1.9dbcc48676f31p-10, 0x1.63e15dec2a59ap-7, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.dea7655eba68cp-9, 0x0.0p+0, 0x1.03dc8ad5c043cp-8, 0x1.a3122f5a30b3cp-9, 0x1.0d66aa9aa246dp-7, -0x1.583828be35fbfp-10, 0x1.4032033a9b1abp-7, 0x1.ec46b2a5b9d40p-9, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.e5ea5e4c0ca7bp-9, 0x0.0p+0, 0x1.4d0e94bdf96f6p-8, 0x1.33f151ce8894dp-10, -0x1.5a7b96b7e57c5p-6, 0x1.ec608f457e72fp-5, 0x1.49699052bcfa0p-6, -0x1.d07faa261596ep-6, 0x1.e792def12339ap-7, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.1fbf40f2d12efp-5, 0x0.0p+0, -0x1.1a522483073a7p-7, -0x1.7b8c5ff9be6d7p-8, 0x1.f20fe184b8ac4p+0, -0x1.b7b2238e37a66p+1, -0x1.bfb8b55ca214cp-4, 0x1.2cc0f2bd84becp+1, -0x1.8313aa4b8ee13p-1, 0x1.c0a1739a7f84cp-4, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.5053cd66e2441p-6, 0x0.0p+0, -0x1.dd864a8eb77e5p-8, -0x1.15438e4faba11p-9, 0x1.daebdec5edc22p-1, -0x1.a6fbfc225cb71p+0, -0x1.595e186e4c1bbp-6, 0x1.34df925789ee1p+0, -0x1.a7a4c49267c14p-2, 0x1.73e9612146663p-4, 0x1.5efb434573ff2p-8, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {-0x1.255b6d39436dep-3, 0x0.0p+0, 0x1.9a8a627bbcdd1p-7, 0x1.bf59991234e61p-8, -0x1.332bea64d5450p+2, 0x1.6cb645d0fc2b7p+2, 0x1.82bc5291f185fp-1, -0x1.053b7dffaad9ep-3, -0x1.f5e9659fbbfefp+0, 0x1.d64d9d726b004p-1, -0x1.e9106c35aba0cp-3, 0x1.45dbe3e5dd845p-3, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {0x1.9be7acfcae340p-1, 0x0.0p+0, -0x1.10eeed3cf334bp-6, -0x1.5ee084b5191eap-6, 0x1.0d3c5bbdcf850p+4, -0x1.6587de03d53e0p+3, -0x1.b046c17f7c157p+1, -0x1.e7c954fb60a30p+3, 0x1.12e09b7064988p+4, -0x1.5c03978ef7f04p+2, 0x1.634b432747ab6p+0, -0x1.2f67025eb1d9fp-1, 0x1.e507a4c6c308cp-6, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {-0x1.d39ba24103f09p-1, 0x0.0p+0, 0x1.3c0ce1b59dbacp-9, 0x1.214f88dea628fp-6, -0x1.db40b709f16ddp+3, 0x1.145928496694ep+1, 0x1.ffbbb2a4d47d0p+1, 0x1.c6f2492877f39p+4, -0x1.94265b3ce38e9p+4, 0x1.eef7ddb4aba54p+2, -0x1.e42c6d3033a3cp+0, 0x1.006149e51cb91p+0, 0x1.302a687e0cb9fp-8, 0x1.6f9dbcf602f76p-7, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0},
    {-0x1.19cd0ee1e9cd5p-2, 0x0.0p+0, 0x1.2becb21e70b85p-5, 0x1.40c87173c6be6p-7, -0x1.8960c448e8a67p+3, 0x1.c6a1995f50f06p+3, 0x1.962e8b05bd72ap+0, 0x1.3a9852b816a2dp+1, -0x1.1ded1b539c663p+3, 0x1.17ea41086bf0cp+2, -0x1.d5b01725a326dp+0, 0x1.28c1e3e0e56c2p+0, -0x1.1b4893f8a1f25p-6, 0x1.3ca306c4062fdp-6, 0x1.55141035f03adp-8, 0x0.0p+0, 0x0.0p+0},
    {0x1.4ec1710d6d9afp+0, 0x0.0p+0, 0x1.1c7e5a46063ffp-6, -0x1.2fd5167f900f5p-6, 0x1.d9f7fd0de4cf7p+3, 0x1.2c42faf631185p+3, -0x1.4e9e88ec93ed7p+2, -0x1.879c38f6deb2ep+5, 0x1.326073e2ba058p+5, -0x1.52cbd3f212f68p+3, 0x1.377418509c483p+1, -0x1.0b9b711c778f6p+0, 0x1.25fbaa39196a7p-4, 0x1.1b67c2d03791fp-9, 0x1.cb61808bb13b1p-8, 0x0.0p+0, 0x0.0p+0},
};
static const double myb_rkn12_b[17] = {0x1.8d67ef536454dp-7, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.617972b3b9033p-4, 0x1.029bab774bf0fp-2, -0x1.94503eaf7161ap-3, 0x1.a02076a35a7e1p-3, -0x1.54640e0f78f76p-6, 0x1.c13dc5072eebdp-4, 0x1.37d45e6a35e63p-5, 0x1.7d39a48459aa8p-7, 0x1.3144c53e19f25p-8, 0x0.0p+0, 0x0.0p+0};
static const double myb_rkn12_bp[17] = {0x1.8d67ef536454dp-7, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.74140cf313979p-4, 0x1.434296551eed3p-2, -0x1.0d8ad474f6411p-2, 0x1.381858fa83de9p-2, -0x1.54640e0f78f76p-5, 0x1.f9657da814c94p-3, 0x1.37d45e6a35e63p-3, 0x1.4d926ff3ce753p-4, 0x1.5c43ebc18c31bp-4, -0x1.2bff4c0661eacp-7, 0x1.999999999999ap-6};
static const double myb_rkn12_bt[17] = {-0x1.3fdecbc9589eep-8, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.cc00574a1efaep-7, -0x1.e9630ab0cd5f3p-4, 0x1.a29dc275ce951p-3, -0x1.0ee29a79b355dp-3, 0x1.c4bbd4fa104e2p-4, -0x1.46abfb34930b8p-4, 0x1.6fcc822deb584p-7, -0x1.3228637ba2290p-8, 0x1.c1dc5c340011cp-11, 0x0.0p+0, 0x0.0p+0};
static const double myb_rkn12_bpt[17] = {-0x1.3fdecbc9589eep-8, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.e43640efb4d23p-7, -0x1.31dde6ae805b8p-3, 0x1.1713d6f934636p-2, -0x1.9653e7b68d00bp-3, 0x1.c4bbd4fa104e2p-3, -0x1.6f817a9b256cfp-3, 0x1.6fcc822deb584p-5, -0x1.0be3570c2de3ep-5, 0x1.009c8db9471b6p-6, -0x1.ddadba7e123d1p-6, 0x1.999999999999ap-6};

// clang-format on

#endif /* MYB_INTEG_RKN12_COEFFICIENTS_H_ */
