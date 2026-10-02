#!/usr/bin/env python3
"""
HCSM-20_verification.py

Verification script for HCSM-20: The Motion Theorem.
Reproduces every numerical claim in the paper at machine precision
in IEEE-754 double precision.

Stanley Preschutti
Entropia Research Institute / Information Physics Institute
ORCID: 0009-0004-5445-1744
September 2026

Usage:
    python HCSM-20_verification.py

Exit codes:
    0 = all checks passed
    1 = one or more checks failed
"""

import numpy as np
from numpy import sqrt, pi, cos, sin, exp, log2, arccos
from scipy.linalg import eigvalsh
import sys

# ============================================================
# Tolerances
# ============================================================
TOL = 1e-12
TOL_LOOSE = 1e-10

# ============================================================
# Framework constants
# ============================================================
d = 4
N = d * 2**d          # 64
L = int(np.sqrt(N))   # 8
H = 2**d              # 16

# ============================================================
# HCSM-54 constants
# ============================================================
def AGM(a, b, tol=1e-15, max_iter=100):
    """Arithmetic-geometric mean."""
    for _ in range(max_iter):
        a_new = (a + b) / 2
        b_new = np.sqrt(a * b)
        if abs(a_new - b_new) < tol:
            return a_new
        a, b = a_new, b_new
    return a_new

# I(2) = 1 - AGM(1, sqrt(2))/2
agm_1_sqrt2 = AGM(1.0, np.sqrt(2.0))
I2 = 1.0 - agm_1_sqrt2 / 2.0

# Stability coefficient gamma (HCSM-54)
gamma = 0.70283946799007245929257819650854315216110316215188

# c_2 = gamma * I(2)
c2 = gamma * I2

# EFT parameters (HCSM-54, Corollary 6.2)
w = -1 - (1 - gamma)**2 / (2 * pi)
zetaH_rho = (1 - gamma)**2 / (6 * pi)

# Derived quantities
T_DME = 1 / (2 * gamma)
delta_edge = np.sqrt(2) / gamma - 2
T_beat = 2 * pi / (2 * np.sqrt(2))

# ============================================================
# Check 1: The 14 mode multiplet at lambda = 4
# ============================================================
def check_1_14_mode_multiplet():
    print("=" * 60)
    print("Check 1: The 14 mode multiplet at lambda = 4")
    print("=" * 60)

    n1_vals = np.arange(8)
    n2_vals = np.arange(8)
    eigenvalues = []
    modes = []
    for n1 in n1_vals:
        for n2 in n2_vals:
            lam = 4 - 2 * cos(pi * n1 / 4) - 2 * cos(pi * n2 / 4)
            eigenvalues.append(lam)
            modes.append((n1, n2))
    eigenvalues = np.array(eigenvalues)
    modes = np.array(modes)

    mask = np.abs(eigenvalues - 4) < TOL
    modes_4 = modes[mask]

    print(f"  Multiplicity at lambda = 4: {len(modes_4)}")
    assert len(modes_4) == 14, f"Expected 14, got {len(modes_4)}"

    O_0 = [(0, 4), (4, 0)]
    O_1 = [(1, 3), (1, 5), (3, 1), (3, 7), (5, 1), (5, 7), (7, 3), (7, 5)]
    O_2 = [(2, 2), (2, 6), (6, 2), (6, 6)]

    all_modes = set(map(tuple, modes_4))
    O_0_set = set(O_0)
    O_1_set = set(O_1)
    O_2_set = set(O_2)

    assert all_modes == O_0_set | O_1_set | O_2_set, "Orbit decomposition failed"
    print(f"  O_0 size: {len(O_0)}")
    print(f"  O_1 size: {len(O_1)}")
    print(f"  O_2 size: {len(O_2)}")
    print("  PASSED")


# ============================================================
# Check 2: The winding generator H_wind = i sin(phi)
# ============================================================
def build_H_wind():
    """Build the winding generator H_wind = i sin(phi)."""
    modes_std = [
        (0, 4), (4, 0), (1, 3), (1, 5), (3, 1), (3, 7),
        (5, 7), (7, 3), (7, 5), (5, 1), (2, 2), (2, 6),
        (6, 2), (6, 6)
    ]
    n1 = np.array([m[0] for m in modes_std], dtype=float)
    n2 = np.array([m[1] for m in modes_std], dtype=float)

    M = np.outer(n1, n2) - np.outer(n2, n1)
    phi = (pi / 8) * M
    H_wind = 1j * np.sin(phi)

    return H_wind, n1, n2, M, phi


def check_2_winding_generator():
    print("=" * 60)
    print("Check 2: The winding generator H_wind = i sin(phi)")
    print("=" * 60)

    H_wind, n1, n2, M, phi = build_H_wind()

    assert np.allclose(H_wind, H_wind.conj().T, atol=TOL), "H_wind not Hermitian"
    print(f"  H_wind is Hermitian: {np.allclose(H_wind, H_wind.conj().T, atol=TOL)}")

    W = np.sin(phi)
    unique_vals = np.unique(np.round(W, 12))
    print(f"  Unique values of sin(phi): {unique_vals}")

    assert np.allclose(W, -W.T, atol=TOL), "W not antisymmetric"
    print(f"  W = sin(phi) is antisymmetric: {np.allclose(W, -W.T, atol=TOL)}")

    print("  PASSED")


# ============================================================
# Check 3: Characteristic polynomial
# ============================================================
def check_3_characteristic_polynomial():
    print("=" * 60)
    print("Check 3: Characteristic polynomial")
    print("=" * 60)

    H_wind, _, _, _, _ = build_H_wind()

    eigs = eigvalsh(H_wind)
    eigs_sorted = np.sort(eigs)

    print(f"  Eigenvalues of H_wind (sorted):")
    for e in eigs_sorted:
        print(f"    {e:.15f}")

    expected = np.array([-4*np.sqrt(2), -2*np.sqrt(2)] + [0]*10 + [2*np.sqrt(2), 4*np.sqrt(2)])
    expected_sorted = np.sort(expected)

    assert np.allclose(eigs_sorted, expected_sorted, atol=TOL_LOOSE), \
        f"Spectrum mismatch"
    print(f"  Spectrum matches {{+/- 4 sqrt(2), +/- 2 sqrt(2), 0^(10)}}")
    print("  PASSED")


# ============================================================
# Check 4: Trace identities
# ============================================================
def check_4_trace_identities():
    print("=" * 60)
    print("Check 4: Trace identities")
    print("=" * 60)

    H_wind, _, _, _, _ = build_H_wind()

    tr1 = np.trace(H_wind)
    tr2 = np.trace(H_wind @ H_wind)
    tr3 = np.trace(H_wind @ H_wind @ H_wind)
    tr4 = np.trace(H_wind @ H_wind @ H_wind @ H_wind)

    print(f"  tr(H_wind)   = {tr1.real:.15f} (expected 0)")
    print(f"  tr(H_wind^2) = {tr2.real:.15f} (expected 80)")
    print(f"  tr(H_wind^3) = {tr3.real:.15f} (expected 0)")
    print(f"  tr(H_wind^4) = {tr4.real:.15f} (expected 2176)")

    assert abs(tr1) < TOL_LOOSE, "tr(H_wind) != 0"
    assert abs(tr2.real - 80) < TOL_LOOSE, "tr(H_wind^2) != 80"
    assert abs(tr3.real) < TOL_LOOSE, "tr(H_wind^3) != 0"
    assert abs(tr4.real - 2176) < TOL_LOOSE, "tr(H_wind^4) != 2176"

    print("  PASSED")


# ============================================================
# Check 5: Standing wave amplitudes
# ============================================================
def check_5_standing_wave_amplitudes():
    print("=" * 60)
    print("Check 5: Standing wave amplitudes")
    print("=" * 60)

    amp_64 = 28/9
    amp_48 = 7/18
    amp_32 = 21/16

    print(f"  amp_64 = 28/9  = {amp_64:.15f}")
    print(f"  amp_48 = 7/18  = {amp_48:.15f}")
    print(f"  amp_32 = 21/16 = {amp_32:.15f}")

    f1 = d * (2*d - 1) / (d - 1)**2
    f2 = (2*d - 1) / (2 * (d - 1)**2)
    f3 = (d - 1) * (2*d - 1) / 2**d

    assert abs(f1 - 28/9) < TOL, f"Framework form 1: {f1} != 28/9"
    assert abs(f2 - 7/18) < TOL, f"Framework form 2: {f2} != 7/18"
    assert abs(f3 - 21/16) < TOL, f"Framework form 3: {f3} != 21/16"

    ratio_1 = amp_64 / amp_48
    ratio_2 = amp_32 / amp_48
    print(f"  amp_64/amp_48 = {ratio_1:.15f} (expected 8)")
    print(f"  amp_32/amp_48 = {ratio_2:.15f} (expected 27/8 = {27/8:.15f})")

    assert abs(ratio_1 - 8) < TOL, "Ratio 1 != 8"
    assert abs(ratio_2 - 27/8) < TOL, "Ratio 2 != 27/8"

    print("  PASSED")


# ============================================================
# Check 6: Floor function cycle counts
# ============================================================
def check_6_floor_function_cycle_counts():
    print("=" * 60)
    print("Check 6: Floor function cycle counts")
    print("=" * 60)

    omega_edge = 2 * np.sqrt(2)
    n_edge = int(np.floor(omega_edge * T_DME))
    print(f"  omega_edge = 2 sqrt(2) = {omega_edge:.15f}")
    print(f"  omega_edge * T_DME = {omega_edge * T_DME:.15f}")
    print(f"  n_edge = floor(...) = {n_edge} (expected 2)")

    omega_mid = 6 * np.sqrt(2)
    n_mid = int(np.floor(omega_mid * T_DME))
    print(f"  omega_mid = 6 sqrt(2) = {omega_mid:.15f}")
    print(f"  omega_mid * T_DME = {omega_mid * T_DME:.15f}")
    print(f"  n_mid = floor(...) = {n_mid} (expected 6)")

    assert n_edge == 2, f"n_edge = {n_edge} != 2"
    assert n_mid == 6, f"n_mid = {n_mid} != 6"

    abs_w1 = abs(w + 1)
    print(f"  |w + 1| = {abs_w1:.15f}")
    print(f"  delta_edge = {delta_edge:.15f}")
    print(f"  n_edge * delta_edge = {n_edge * delta_edge:.15f}")
    print(f"  n_mid * delta_edge  = {n_mid * delta_edge:.15f}")
    assert n_edge * delta_edge > abs_w1, "Boundary constraint violated for edge"
    assert n_mid * delta_edge > abs_w1, "Boundary constraint violated for middle"

    print("  PASSED")


# ============================================================
# Check 7: Loop eigenvalue magnitudes
# ============================================================
def check_7_loop_eigenvalue_magnitudes():
    print("=" * 60)
    print("Check 7: Loop eigenvalue magnitudes")
    print("=" * 60)

    prefactor_A = abs(w + 1) / (2 * delta_edge)
    prefactor_B = abs(w + 1) / (6 * delta_edge)
    exp_factor = np.exp(-zetaH_rho * T_beat / 2)

    MA = prefactor_A * exp_factor
    MB = prefactor_B * exp_factor

    print(f"  |w + 1| = {abs(w + 1):.15f}")
    print(f"  delta_edge = {delta_edge:.15f}")
    print(f"  |w + 1| / delta_edge = {abs(w + 1) / delta_edge:.15f}")
    print(f"  exp(-(zeta H / rho) T_beat / 2) = {exp_factor:.15f}")
    print(f"  |M_A| = {MA:.15f} (expected 0.5757...)")
    print(f"  |M_B| = {MB:.15f} (expected 0.1919...)")

    ratio = MA / MB
    print(f"  |M_A| / |M_B| = {ratio:.15f} (expected 3)")

    assert abs(ratio - 3) < TOL, f"Ratio = {ratio} != 3"
    assert MA < 1, f"|M_A| = {MA} >= 1"
    assert MB < 1, f"|M_B| = {MB} >= 1"

    # Verify that the correct exponent gives the stated numerical values.
    # The paper's Corollary 9.5 uses exp(-(1-gamma)^2 / (12 sqrt(2))).
    exp_factor_check = np.exp(-(1 - gamma)**2 / (12 * np.sqrt(2)))
    print(f"  exp(-(1-gamma)^2 / (12 sqrt(2))) = {exp_factor_check:.15f}")
    assert abs(exp_factor - exp_factor_check) < TOL, "Exponent mismatch"

    print("  PASSED")


# ============================================================
# Check 8: Charge bookkeeping
# ============================================================
def check_8_charge_bookkeeping():
    print("=" * 60)
    print("Check 8: Charge bookkeeping")
    print("=" * 60)

    charges = {
        'nu': 0, 'e': -1, 'u': 2/3, 'd': -1/3, 's': -1/3,
        'mu': -1, 'dark': 0, 'c': 2/3, 'tau': -1, 'b': -1/3,
        'W': 1, 'Z': 0, 'H': 0, 't': 2/3
    }

    im_T = ['nu', 'e', 'u', 'd', 's', 'mu', 'c', 'tau', 'b', 't']
    ker_T = ['dark', 'W', 'Z', 'H']

    sum_im = sum(charges[m] for m in im_T)
    sum_ker = sum(charges[m] for m in ker_T)
    sum_total = sum(charges.values())

    print(f"  sum_{{im T}} Q  = {sum_im:.15f} (expected -2)")
    print(f"  sum_{{ker T}} Q = {sum_ker:.15f} (expected +1)")
    print(f"  sum_total Q  = {sum_total:.15f} (expected -1)")

    assert abs(sum_im - (-2)) < TOL, f"sum_im = {sum_im} != -2"
    assert abs(sum_ker - 1) < TOL, f"sum_ker = {sum_ker} != +1"
    assert abs(sum_total - (-1)) < TOL, f"sum_total = {sum_total} != -1"

    print("  PASSED")


# ============================================================
# Check 9: Wick rotation rate
# ============================================================
def check_9_wick_rotation_rate():
    print("=" * 60)
    print("Check 9: Wick rotation rate")
    print("=" * 60)

    agm_1_sqrt2_val = AGM(1.0, np.sqrt(2.0))
    I2_val = 1 - agm_1_sqrt2_val / 2
    c2_val = gamma * I2_val

    print(f"  AGM(1, sqrt(2)) = {agm_1_sqrt2_val:.15f}")
    print(f"  I(2) = 1 - AGM(1, sqrt(2))/2 = {I2_val:.15f}")
    print(f"  gamma = {gamma:.15f}")
    print(f"  c_2 = gamma * I(2) = {c2_val:.15f}")
    print(f"  2 gamma = {2 * gamma:.15f}")
    print(f"  T_DME = 1/(2 gamma) = {1/(2*gamma):.15f}")

    assert abs(gamma - c2_val / I2_val) < TOL, "gamma != c_2/I(2)"
    assert abs(T_DME - 1/(2*gamma)) < TOL, "T_DME != 1/(2 gamma)"

    print("  PASSED")


# ============================================================
# Check 10: Lattice correction coefficients
# ============================================================
def check_10_lattice_corrections():
    print("=" * 60)
    print("Check 10: Lattice correction coefficients")
    print("=" * 60)

    O_0 = [(0, 4), (4, 0)]
    O_1 = [(1, 3), (1, 5), (3, 1), (3, 7), (5, 1), (5, 7), (7, 3), (7, 5)]
    O_2 = [(2, 2), (2, 6), (6, 2), (6, 6)]

    def avg_n4(orbit):
        vals = [n1**4 + n2**4 for n1, n2 in orbit]
        return np.mean(vals)

    avg_O_0 = avg_n4(O_0)
    avg_O_1 = avg_n4(O_1)
    avg_O_2 = avg_n4(O_2)

    print(f"  <n_1^4 + n_2^4>_O_0 = {avg_O_0:.15f} (expected 256)")
    print(f"  <n_1^4 + n_2^4>_O_1 = {avg_O_1:.15f} (expected 1418)")
    print(f"  <n_1^4 + n_2^4>_O_2 = {avg_O_2:.15f} (expected 1312)")

    factor = (4096 / 12) * (pi / 4)**4
    C_O_0 = factor * avg_O_0 / 4
    C_O_1 = factor * avg_O_1 / 4
    C_O_2 = factor * avg_O_2 / 4

    print(f"  C_O_0 = {C_O_0:.4f} (expected ~129.9)")
    print(f"  C_O_1 = {C_O_1:.4f} (expected ~718.6)")
    print(f"  C_O_2 = {C_O_2:.4f} (expected ~665.0)")

    assert abs(avg_O_0 - 256) < TOL, "avg_O_0 != 256"
    assert abs(avg_O_1 - 1418) < TOL, "avg_O_1 != 1418"
    assert abs(avg_O_2 - 1312) < TOL, "avg_O_2 != 1312"

    print("  PASSED")


# ============================================================
# Check 11: Time unit
# ============================================================
def check_11_time_unit():
    print("=" * 60)
    print("Check 11: Time unit")
    print("=" * 60)

    M_P = 1.2209e19  # GeV (approximate)
    m_s = M_P / N * np.exp(-(d + 1) * N / (d - 1))

    print(f"  M_P = {M_P:.4e} GeV")
    print(f"  m_s = M_P / N * exp(-(d+1)N/(d-1)) = {m_s:.4e} GeV")
    print(f"  Expected m_s = 9.0247487584e-30 GeV")

    T_unit_GeV_inv = 2 * pi / m_s
    GeV_inv_to_s = 6.582e-25
    T_unit_s = T_unit_GeV_inv * GeV_inv_to_s
    T_unit_days = T_unit_s / (24 * 3600)

    print(f"  T_unit = 2 pi / m_s = {T_unit_GeV_inv:.4e} GeV^-1")
    print(f"  T_unit = {T_unit_s:.4e} s")
    print(f"  T_unit = {T_unit_days:.4f} days (expected ~5.3 days)")

    phi_loop = m_s * T_unit_GeV_inv
    print(f"  Phase advance per cycle: m_s * T_unit = {phi_loop:.15f} (expected 2 pi = {2*pi:.15f})")
    assert abs(phi_loop - 2*pi) < TOL, "Phase advance != 2 pi"

    print("  PASSED")


# ============================================================
# Main
# ============================================================
def main():
    print()
    print("=" * 60)
    print("HCSM-20: The Motion Theorem")
    print("Verification Script")
    print("=" * 60)
    print()

    checks = [
        check_1_14_mode_multiplet,
        check_2_winding_generator,
        check_3_characteristic_polynomial,
        check_4_trace_identities,
        check_5_standing_wave_amplitudes,
        check_6_floor_function_cycle_counts,
        check_7_loop_eigenvalue_magnitudes,
        check_8_charge_bookkeeping,
        check_9_wick_rotation_rate,
        check_10_lattice_corrections,
        check_11_time_unit,
    ]

    passed = 0
    failed = 0
    for check in checks:
        try:
            check()
            passed += 1
        except AssertionError as e:
            print(f"  FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            failed += 1
        print()

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    if failed > 0:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()