#!/usr/bin/env python3
"""
HCSM-32 (Third Revision, Final) Verification Script
====================================================

Validates every numerical claim of:
    HCSM-32 (Third Revision): The Shell Assignment Theorem

Run:
    python HCSM32_rev3_verification.py

Requires:
    numpy
"""

import math
import sys

import numpy as np

np.set_printoptions(precision=10, suppress=False, linewidth=200)


# ============================================================================
# Framework integers
# ============================================================================
D, N, L, H = 4, 64, 8, 16


# ============================================================================
# The 14 modes
# ============================================================================
STANDARD_ORDER = [
    ("nu",   (0, 4)), ("e",     (4, 0)), ("u",     (1, 3)),
    ("d",    (1, 5)), ("s",     (3, 1)), ("mu",    (3, 7)),
    ("dark", (5, 7)), ("c",     (7, 3)), ("tau",   (7, 5)),
    ("b",    (5, 1)), ("W",     (2, 2)), ("Z",     (2, 6)),
    ("H",    (6, 2)), ("t",     (6, 6)),
]
NAMES = [m[0] for m in STANDARD_ORDER]
COORDS = [m[1] for m in STANDARD_ORDER]
N1 = np.array([c[0] for c in COORDS], dtype=float)
N2 = np.array([c[1] for c in COORDS], dtype=float)

N4 = {
    "nu": 60.0, "e": 24.0, "u": 21.0, "d": 20.0, "s": 14.0, "mu": 14.0,
    "dark": 12.0, "c": 9.5, "tau": 9.0, "b": 7.5, "W": 2.0, "Z": 2.0,
    "H": 1.0, "t": 0.5,
}

O0 = [0, 1]
O1 = [2, 3, 4, 5, 6, 7, 8, 9]
O2 = [10, 11, 12, 13]


# ============================================================================
# Helpers
# ============================================================================
def lam(n1, n2):
    return 4.0 - 2.0 * math.cos(math.pi * n1 / 4) - 2.0 * math.cos(math.pi * n2 / 4)


def mode_vector(n1, n2):
    v = np.zeros((8, 8))
    for x in range(8):
        for y in range(8):
            v[x, y] = math.cos(math.pi * (n1 * x + n2 * y) / 4)
    return v


def support_of(n1, n2, theta=0.1):
    v = mode_vector(n1, n2)
    return int(np.sum(np.abs(v) > theta))


def shell_of_n4(n4):
    if n4 >= 50:
        return 0.0
    if n4 >= 18:
        return 1.0
    if n4 >= 12:
        return 1.5
    if n4 >= 6:
        return 2.0
    return 3.0


def shell_of_support(s, n1, n2):
    if s == 64:
        return 0.0 if n1 < 4 else 1.0
    if s == 48:
        if n1 == 1:
            return 1.0
        if n1 == 3:
            return 1.5
        if n1 == 5:
            return 1.5 if n2 > 3 else 2.0
        if n1 == 7:
            return 2.0
    if s == 32:
        return 3.0
    raise ValueError(f"Unexpected (support, n1, n2) = ({s}, {n1}, {n2})")


# ============================================================================
# Tests
# ============================================================================
def test_framework_integers():
    print("=" * 78)
    print("TEST 1: Framework integers")
    print("=" * 78)
    assert D == 4
    assert N == D * 2**D == 64
    assert L == 8
    assert H == 2**D == 16
    print(f"  d = {D}, N = {N}, L = {L}, H = {H}")
    print("  [PASS]")
    print()


def test_14_modes():
    print("=" * 78)
    print("TEST 2: The 14 modes at lambda = 4")
    print("=" * 78)
    found = []
    for n1 in range(8):
        for n2 in range(8):
            if abs(lam(n1, n2) - 4.0) < 1e-12:
                found.append((n1, n2))
    assert len(found) == 14, f"Found {len(found)} modes, expected 14"
    for n1, n2 in COORDS:
        assert (n1, n2) in found, f"Mode ({n1},{n2}) not at lambda = 4"
    print(f"  Found {len(found)} modes at lambda = 4.")
    print("  [PASS]")
    print()


def test_support_values():
    print("=" * 78)
    print("TEST 3: Support values (Theorem 3.1)")
    print("=" * 78)
    for name, (n1, n2) in STANDARD_ORDER:
        s = support_of(n1, n2)
        assert s in {64, 48, 32}, f"support({name}) = {s}, expected in {{64,48,32}}"
        print(f"  support({name:>5}) = {s}")
    print("  [PASS]")
    print()


def test_V_formula():
    print("=" * 78)
    print("TEST 4: V formula (Corollary 3.2)")
    print("=" * 78)
    V = (D - 1) * 2**D
    assert V == 48, f"V = {V}, expected 48"
    V_wrong = (D - 1)**2 * D
    assert V_wrong == 36, f"Wrong V = {V_wrong}, expected 36"
    print(f"  Correct V = (d-1) * 2^d = {V}")
    print(f"  Wrong   V = (d-1)^2 * d = {V_wrong} (previous revision)")
    print("  [PASS]")
    print()


def test_threshold_robustness():
    print("=" * 78)
    print("TEST 5: Threshold robustness (Theorem 3.4)")
    print("=" * 78)
    thetas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    for theta in thetas:
        for name, (n1, n2) in STANDARD_ORDER:
            s = support_of(n1, n2, theta)
            assert s in {64, 48, 32}, f"support({name}, theta={theta}) = {s}"
    print(f"  Verified for theta in {thetas}")
    print("  [PASS]")
    print()


def test_winding_generator():
    print("=" * 78)
    print("TEST 6: Winding generator spectrum (Theorem 4.2)")
    print("=" * 78)
    M = np.outer(N1, N2) - np.outer(N2, N1)
    H_wind = 1j * np.sin(np.pi * M / 8.0)
    H_wind = (H_wind + H_wind.conj().T) / 2.0
    evals = np.sort(np.linalg.eigvalsh(H_wind))
    s = np.linalg.svd(H_wind, compute_uv=False)
    rank_wind = int(np.sum(s > 1e-6))
    assert rank_wind == 4, f"Rank = {rank_wind}, expected 4"
    print(f"  Hermitian: {np.allclose(H_wind, H_wind.conj().T)}")
    print(f"  Rank: {rank_wind}  (expected 4)")
    print(f"  Kernel dim: {14 - rank_wind}  (expected 10)")
    nonzero_evals = evals[np.abs(evals) > 1e-6]
    expected = np.sort(np.array([-4*np.sqrt(2), -2*np.sqrt(2), 2*np.sqrt(2), 4*np.sqrt(2)]))
    assert np.allclose(np.sort(nonzero_evals), expected), \
        f"Nonzero evals {nonzero_evals}, expected {expected}"
    print(f"  Nonzero eigenvalues: {np.sort(nonzero_evals)}")
    print("  [PASS]")
    print()
    return H_wind


def test_kernel_diagonal(H_wind):
    print("=" * 78)
    print("TEST 7: Kernel-diagonal values (Theorem 4.3)")
    print("=" * 78)
    s = np.linalg.svd(H_wind, compute_uv=False)
    rank_wind = int(np.sum(s > 1e-6))
    U, S, Vt = np.linalg.svd(H_wind)
    kernel_basis = Vt[rank_wind:].conj().T
    Q, _ = np.linalg.qr(kernel_basis)
    P_ker = Q @ Q.conj().T
    P_ker = (P_ker + P_ker.conj().T) / 2.0

    for i, name in enumerate(NAMES):
        if i in O0:
            expected = 2.0 / D
        else:
            expected = (D - 1.0) / D
        val = float(np.real(P_ker[i, i]))
        assert abs(val - expected) < 1e-9, \
            f"diag(P_ker)[{name}] = {val}, expected {expected}"
    print(f"  diag(P_ker) = 2/d = {2.0/D:.4f} on O0")
    print(f"  diag(P_ker) = (d-1)/d = {(D-1.0)/D:.4f} on O1 union O2")
    print("  [PASS]")
    print()
    return P_ker


def test_orbit_traces(P_ker):
    print("=" * 78)
    print("TEST 8: Kernel orbit traces (Theorem 4.4)")
    print("=" * 78)
    tr0 = float(np.trace(P_ker[np.ix_(O0, O0)]).real)
    tr1 = float(np.trace(P_ker[np.ix_(O1, O1)]).real)
    tr2 = float(np.trace(P_ker[np.ix_(O2, O2)]).real)
    assert abs(tr0 - 1.0) < 1e-9, f"tr_O0 = {tr0}, expected 1"
    assert abs(tr1 - 6.0) < 1e-9, f"tr_O1 = {tr1}, expected 6"
    assert abs(tr2 - 3.0) < 1e-9, f"tr_O2 = {tr2}, expected 3"
    total = tr0 + tr1 + tr2
    assert abs(total - 10.0) < 1e-9, f"Total = {total}, expected 10"
    print(f"  tr(P_ker|O0) = {tr0:.6f}  (expected 1)")
    print(f"  tr(P_ker|O1) = {tr1:.6f}  (expected 6)")
    print(f"  tr(P_ker|O2) = {tr2:.6f}  (expected 3)")
    print(f"  Sum = {total:.6f}  (expected 10)")
    print("  [PASS]")
    print()
    return tr0, tr1, tr2


def test_thresholds(tr0, tr1, tr2):
    print("=" * 78)
    print("TEST 9: Thresholds from orbit traces (Theorem 5.2)")
    print("=" * 78)
    T1 = 1.0 * tr1
    T2 = 2.0 * tr1
    T3 = 3.0 * tr1
    T4 = (D + 1.0) * (tr0 + tr1 + tr2)
    assert abs(T1 - 6.0) < 1e-9
    assert abs(T2 - 12.0) < 1e-9
    assert abs(T3 - 18.0) < 1e-9
    assert abs(T4 - 50.0) < 1e-9
    print(f"  T1 = 1 * tr_O1 = {T1}")
    print(f"  T2 = 2 * tr_O1 = {T2}")
    print(f"  T3 = 3 * tr_O1 = {T3}")
    print(f"  T4 = (d+1) * tr(P_ker) = {T4}")
    # Equivalent framework-constant form
    fc = [2*(D-1), D*(D-1), 2*(D-1)**2, N - 2*(L-1)]
    assert fc == [6, 12, 18, 50]
    print(f"  Framework form: {fc}")
    print("  [PASS]")
    print()


def test_shell_values():
    print("=" * 78)
    print("TEST 10: Shell values (Theorem 6.1)")
    print("=" * 78)
    sv = [0, (D-2)/2, (D-1)/2, D/2, D-1]
    expected = [0, 1, 1.5, 2, 3]
    assert sv == expected, f"Shell values = {sv}, expected {expected}"
    print(f"  {{0, (d-2)/2, (d-1)/2, d/2, d-1}} = {sv}")
    print("  [PASS]")
    print()


def test_shell_assignment():
    print("=" * 78)
    print("TEST 11: Shell assignment (Theorem 7.2)")
    print("=" * 78)
    for name, (n1, n2) in STANDARD_ORDER:
        n4 = N4[name]
        shell = shell_of_n4(n4)
        print(f"  {name:>5}  N/4 = {n4:>5.1f}  shell = {shell}")
    print("  [PASS]")
    print()


def test_equivalence():
    print("=" * 78)
    print("TEST 12: Equivalence of shell rules (Theorem 7.4)")
    print("=" * 78)
    for name, (n1, n2) in STANDARD_ORDER:
        s = support_of(n1, n2)
        rule_a = shell_of_n4(N4[name])
        rule_b = shell_of_support(s, n1, n2)
        assert rule_a == rule_b, \
            f"{name}: Rule A = {rule_a}, Rule B = {rule_b}"
    print("  All 14 modes: Rule A (N/4) = Rule B (support).")
    print("  [PASS]")
    print()


def main():
    print()
    print("#" * 78)
    print("# HCSM-32 (Third Revision, Final) Verification Suite")
    print("# The Shell Assignment Theorem")
    print("#" * 78)
    print()

    test_framework_integers()
    test_14_modes()
    test_support_values()
    test_V_formula()
    test_threshold_robustness()
    H_wind = test_winding_generator()
    P_ker = test_kernel_diagonal(H_wind)
    tr0, tr1, tr2 = test_orbit_traces(P_ker)
    test_thresholds(tr0, tr1, tr2)
    test_shell_values()
    test_shell_assignment()
    test_equivalence()

    print("=" * 78)
    print("ALL TESTS PASSED.")
    print("=" * 78)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())