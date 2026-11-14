/* SPDX-License-Identifier: BSD-3-Clause */
/* pn_flux_test_terms.h: the gravitational-wave energy flux of a compact binary for general
 * orbits, F = (8/15) nu^2 / r^4 {f_N + f_1PN + f_S1.5 + f_S2.5} (G m = c = 1), of L. Blanchet,
 * A. Buonanno, G. Faye, PRD 74, 104034 (2006), gr-qc/0605140, eqs. flux1PN, flux15PN, flux25PN,
 * printed by machine from its LaTeX source (sympy, a script kept outside the repository), for
 * tests/unit/test_pn_sorr.c. The spin terms are in the spin variables of that paper; with the
 * conserved-norm ones the 2.5PN term gets -f_S1.5(dS), dS_1 = -(1/2)(v_1.S_1) v_1 + (m_2 / r) S_1
 * (eq. S1K). nsv = n.(S x v), nzv = n.(Sigma x v). */
#ifndef MYB_TESTS_PN_FLUX_TEST_TERMS_H_
#define MYB_TESTS_PN_FLUX_TEST_TERMS_H_

static double pn_flux_fN(double rd, double v2) {
  const double rd2 = rd * rd, rd4 = rd2 * rd2, v22 = v2 * v2;
  (void)rd2; (void)rd4; (void)v22;
  return v2 * (12.0) +
    rd2 * (-11.0);
}

static double pn_flux_f1PN(double nu, double ri, double rd, double v2) {
  const double rd2 = rd * rd, rd4 = rd2 * rd2, v22 = v2 * v2;
  (void)rd2; (void)rd4; (void)v22;
  const double ri2 = ri * ri, ri3 = ri2 * ri;
  (void)ri2; (void)ri3;
  return v22 * (785.0 / 28.0 - 213.0 / 7.0 * nu) +
    rd2 * v2 * (-1487.0 / 14.0 + 696.0 / 7.0 * nu) +
    rd4 * (2061.0 / 28.0 - 465.0 / 7.0 * nu) +
    ri * v2 * (-680.0 / 7.0 + 40.0 / 7.0 * nu) +
    ri * rd2 * (734.0 / 7.0 - 30.0 / 7.0 * nu) +
    ri2 * (4.0 / 7.0 - 16.0 / 7.0 * nu);
}

static double pn_flux_fS15(double dm, double ri, double rd, double v2, double nsv, double nzv) {
  const double rd2 = rd * rd, rd4 = rd2 * rd2, v22 = v2 * v2;
  (void)rd2; (void)rd4; (void)v22;
  const double ri2 = ri * ri, ri3 = ri2 * ri;
  (void)ri2; (void)ri3;
  return nzv * ri * v2 * (43.0 * dm) +
    nzv * ri * rd2 * (-51.0 * dm) +
    nzv * ri2 * (-4.0 * dm) +
    nsv * ri * v2 * (80.0) +
    nsv * ri * rd2 * (-78.0) +
    nsv * ri2 * (8.0);
}

static double pn_flux_fS25(double nu, double dm, double ri, double rd, double v2, double nsv, double nzv) {
  const double rd2 = rd * rd, rd4 = rd2 * rd2, v22 = v2 * v2;
  (void)rd2; (void)rd4; (void)v22;
  const double ri2 = ri * ri, ri3 = ri2 * ri;
  (void)ri2; (void)ri3;
  return nzv * ri * v22 * (2603.0 / 28.0 * dm - 1040.0 / 7.0 * nu * dm) +
    nzv * ri * rd2 * v2 * (-2364.0 / 7.0 * dm + 3621.0 / 7.0 * nu * dm) +
    nzv * ri * rd4 * (7941.0 / 28.0 * dm - 2676.0 / 7.0 * nu * dm) +
    nzv * ri2 * v2 * (-4785.0 / 14.0 * dm + 65.0 * nu * dm) +
    nzv * ri2 * rd2 * (6613.0 / 14.0 * dm - 1031.0 / 7.0 * nu * dm) +
    nzv * ri3 * (109.0 / 7.0 * dm - 18.0 * nu * dm) +
    nsv * ri * v22 * (1207.0 / 7.0 - 1810.0 / 7.0 * nu) +
    nsv * ri * rd2 * v2 * (-3519.0 / 7.0 + 5004.0 / 7.0 * nu) +
    nsv * ri * rd4 * (2244.0 / 7.0 - 3144.0 / 7.0 * nu) +
    nsv * ri2 * v2 * (-3504.0 / 7.0 + 20.0 * nu) +
    nsv * ri2 * rd2 * (2866.0 / 7.0 - 170.0 / 7.0 * nu) +
    nsv * ri3 * (-972.0 / 7.0 - 166.0 / 7.0 * nu);
}

#endif /* MYB_TESTS_PN_FLUX_TEST_TERMS_H_ */
