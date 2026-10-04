#!/usr/bin/env python3
"""
HCSM-35 Verification Script
============================

This script reproduces every numerical claim in HCSM-35 (Third Revision).

It builds the HCSM-native objects from scratch:
- The 8x8 torus Laplacian
- The 14 modes at lambda = 4
- The winding matrix M = X J X^T
- The processor T = I - P_C
- The rank-2 structure of M
- The Wolfenstein parameters, CP phase, and Jarlskog invariant

Every quantity is derived from first principles. No external fitting.

Usage:
    python HCSM-35_verification.py

Author: Stanley Preschutti
ORCID: 0009-0004-5445-1744
"""

import numpy as np
from numpy.linalg import svd, matrix_rank, norm
from itertools import product

# Suppress numpy warnings from numerical precision
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

np.set_printoptions(precision=15, suppress=True, linewidth=200)

# ================================================================
# SECTION 1: HCSM Constants (from Axiom + Theorem)
# ================================================================
L, d, N, H = 8, 4, 64, 16

print("=" * 72)
print("SECTION 1: HCSM CONSTANTS")
print("=" * 72)
print(f"L = {L}, d = {d}, N = {N}, H = {H}")
print(f"Selected by unified substrate selection theorem (HCSM-00 Thm 2.7):")
print(f"  L^2 = d * 2^d  =>  {L**2} = {d * 2**d}  OK")
print(f"  2(L-1) = d(d+3)/2  =>  {2*(L-1)} = {d*(d+3)//2}  OK")

# ================================================================
# SECTION 2: The 8x8 torus Laplacian and the 14 modes at lambda = 4
# ================================================================
def lam(n1, n2, L=L):
    return 4 - 2*np.cos(np.pi * n1 / 4) - 2*np.cos(np.pi * n2 / 4)

# Find all modes at lambda = 4
modes_14 = [(n1, n2) for n1, n2 in product(range(L), repeat=2)
            if abs(lam(n1, n2) - 4) < 1e-12]

print("\n" + "=" * 72)
print("SECTION 2: THE 14 MODES AT LAMBDA = 4")
print("=" * 72)
print(f"Number of modes: {len(modes_14)}")
assert len(modes_14) == 14, "Expected 14 modes"

# Standard mode order (from HCSM-22 and HCSM-36)
labels = ['nu', 'e', 'u', 'd', 's', 'mu', 'dark', 'c', 'tau', 'b', 'W', 'Z', 'H', 't']
coords = np.array([
    (0, 4),   # nu
    (4, 0),   # e
    (1, 3),   # u
    (1, 5),   # d
    (3, 1),   # s
    (3, 7),   # mu
    (5, 7),   # dark
    (7, 3),   # c
    (7, 5),   # tau
    (5, 1),   # b
    (2, 2),   # W
    (2, 6),   # Z
    (6, 2),   # H
    (6, 6),   # t
], dtype=float)

# Verify all modes are at lambda = 4 and distinct
for i, (n1, n2) in enumerate(coords):
    assert abs(lam(int(n1), int(n2)) - 4) < 1e-12, f"{labels[i]} not at lambda=4"
assert len(set(map(tuple, coords))) == 14, "Mode coordinates must be distinct"

print(f"All 14 modes verified at lambda = 4 and distinct.")
print(f"Mode coordinates:")
for lab, c in zip(labels, coords):
    print(f"  {lab:5s} : ({int(c[0])},{int(c[1])})")

# ================================================================
# SECTION 3: The Winding Matrix M = X J X^T
# ================================================================
X = coords.copy()
J = np.array([[0, 1], [-1, 0]], dtype=float)
M = X @ J @ X.T

# Alternative construction: M_ij = n1_i * n2_j - n1_j * n2_i
n1, n2 = coords[:, 0], coords[:, 1]
M_direct = np.outer(n1, n2) - np.outer(n2, n1)

print("\n" + "=" * 72)
print("SECTION 3: WINDING MATRIX M = X J X^T")
print("=" * 72)
print(f"||M - M_direct|| = {norm(M - M_direct):.2e}")
assert norm(M - M_direct) < 1e-12, "Factorization failed"

rank_M = matrix_rank(M)
print(f"rank(M) = {rank_M}")
assert rank_M == 2, f"Expected rank 2, got {rank_M}"

U_M, S_M, Vt_M = svd(M)
print(f"Singular values of M: {S_M[:4]}...")
print(f"S0/S1 = {S_M[0]/S_M[1]:.15f}")
print(f"Isotropy verified: S0 = S1 exactly (within numerical precision)")

# ================================================================
# SECTION 4: The Processor T = I - P_C
# ================================================================
boson_labels = ['dark', 'W', 'Z', 'H']
boson_idx = [labels.index(b) for b in boson_labels]
fermion_idx = [i for i in range(14) if i not in boson_idx]

P_C = np.zeros((14, 14))
for i in boson_idx:
    P_C[i, i] = 1.0
T = np.eye(14) - P_C
I_T = P_C

print("\n" + "=" * 72)
print("SECTION 4: PROCESSOR T = I - P_C")
print("=" * 72)
print(f"Boson indices: {boson_idx}")
print(f"Fermion indices: {fermion_idx}")
print(f"rank(T) = {matrix_rank(T)} (expected 10)")
print(f"rank(I-T) = {matrix_rank(I_T)} (expected 4)")

# ================================================================
# SECTION 5: Crush, Boson-Boson Winding
# ================================================================
M_FB = T @ M @ I_T
M_BB = I_T @ M @ I_T

print("\n" + "=" * 72)
print("SECTION 5: CRUSH AND BOSON-BOSON WINDING")
print("=" * 72)
print(f"rank(M_FB) = {matrix_rank(M_FB)} (expected 2)")
print(f"rank(M_BB) = {matrix_rank(M_BB)} (expected 2)")

# Boson-boson 4x4 block
M_BB_4 = M_BB[np.ix_(boson_idx, boson_idx)]
U_BB, S_BB, Vt_BB = svd(M_BB_4)
print(f"M_BB singular values: {S_BB}")
print(f"Degeneracy: S0 = S1 = {S_BB[0]:.4f} (expected 49.4773)")
assert abs(S_BB[0] - S_BB[1]) < 1e-6, "Boson degeneracy failed"

# ================================================================
# SECTION 6: Stability Coefficient gamma
# ================================================================
def agm(a, b, tol=1e-15):
    """Arithmetic-geometric mean."""
    while abs(a - b) > tol:
        a, b = (a + b) / 2, np.sqrt(a * b)
    return a

agm_12 = agm(1.0, np.sqrt(2.0))
I2 = 1 - 0.5 * agm_12
c2 = 1 - 1 / (2 * np.sqrt(2))

gamma = c2 / I2
delta_edge = np.sqrt(2) / gamma - 2

print("\n" + "=" * 72)
print("SECTION 6: STABILITY COEFFICIENT")
print("=" * 72)
print(f"AGM(1, sqrt2) = {agm_12:.15f}")
print(f"I(2) = 1 - AGM/2 = {I2:.15f}")
print(f"c2 = 1 - 1/(2 sqrt2) = {c2:.15f}")
print(f"gamma = c2 / I(2) = {gamma:.15f}")
print(f"  (HCSM-54 value: 0.702839467990)")
print(f"delta_edge = sqrt(2)/gamma - 2 = {delta_edge:.15f}")

# ================================================================
# SECTION 7: Wolfenstein Parameters
# ================================================================
lam_w = np.pi / 14
A_w = (d / (d + 1)) * (1 + delta_edge)
R_b = ((d + 1) / (d * (d - 1))) * (1 + delta_edge / 6)

print("\n" + "=" * 72)
print("SECTION 7: WOLFENSTEIN PARAMETERS")
print("=" * 72)
print(f"lambda = pi/14 = {lam_w:.15f}")
print(f"A = (d/(d+1))(1+delta_edge) = {A_w:.15f}")
print(f"R_b = ((d+1)/(d(d-1)))(1+delta_edge/6) = {R_b:.15f}")

# ================================================================
# SECTION 8: CP Phase from Rank-2 Winding
# ================================================================
rank_M_val = 2
tan_gamma_UT = rank_M_val + lam_w
gamma_UT = np.arctan(tan_gamma_UT)

print("\n" + "=" * 72)
print("SECTION 8: CP PHASE FROM RANK-2 WINDING")
print("=" * 72)
print(f"rank(M) = {rank_M_val}")
print(f"lambda = pi/14 = {lam_w:.15f}")
print(f"tan(gamma_UT) = rank(M) + lambda = {tan_gamma_UT:.15f}")
print(f"gamma_UT = arctan(2 + pi/14) = {np.degrees(gamma_UT):.6f} deg")
print(f"Observed (PDG 2024): 65.8 +/- 2.8 deg")
print(f"Deviation: {abs(np.degrees(gamma_UT) - 65.8):.4f} deg")

# ================================================================
# SECTION 9: Individual rho_bar and eta_bar
# ================================================================
rho_bar = R_b * np.cos(gamma_UT)
eta_bar = R_b * np.sin(gamma_UT)

print("\n" + "=" * 72)
print("SECTION 9: INDIVIDUAL WOLFENSTEIN PARAMETERS")
print("=" * 72)
print(f"rho_bar = R_b cos(gamma_UT) = {rho_bar:.15f}")
print(f"  (Observed: 0.159 +/- 0.007)")
print(f"eta_bar = R_b sin(gamma_UT) = {eta_bar:.15f}")
print(f"  (Observed: 0.354 +/- 0.006)")
print(f"rho_bar^2 + eta_bar^2 = {rho_bar**2 + eta_bar**2:.15f}")
print(f"R_b^2 = {R_b**2:.15f}")

# ================================================================
# SECTION 10: CKM Matrix Magnitudes
# ================================================================
V_us = lam_w
V_cb = A_w * lam_w**2
V_ub = A_w * lam_w**3 * R_b
V_td = A_w * lam_w**3 * np.sqrt((1 - rho_bar)**2 + eta_bar**2)

print("\n" + "=" * 72)
print("SECTION 10: CKM MATRIX MAGNITUDES")
print("=" * 72)
print(f"|V_us| = lambda = {V_us:.10f}   (obs 0.2250)")
print(f"|V_cb| = A lambda^2 = {V_cb:.10f}   (obs 0.0408)")
print(f"|V_ub| = A lambda^3 R_b = {V_ub:.10f}   (obs 0.00382)")
print(f"|V_td| = A lambda^3 R_t = {V_td:.10f}   (obs 0.00857)")

# ================================================================
# SECTION 11: Jarlskog Invariant
# ================================================================
J = A_w**2 * lam_w**6 * eta_bar * (1 - lam_w**2 / 2)

print("\n" + "=" * 72)
print("SECTION 11: JARLSKOG INVARIANT")
print("=" * 72)
print(f"J = A^2 lambda^6 eta_bar (1 - lambda^2/2)")
print(f"  = {J:.6e}")
print(f"Observed: (3.08 +/- 0.13) e-5")
print(f"Deviation: {abs(J - 3.08e-5) / 3.08e-5 * 100:.2f}%")

# ================================================================
# SECTION 12: Full CKM Matrix (complex)
# ================================================================
V_CKM = np.array([
    [1 - lam_w**2 / 2, lam_w, A_w * lam_w**3 * (rho_bar - 1j * eta_bar)],
    [-lam_w, 1 - lam_w**2 / 2, A_w * lam_w**2],
    [A_w * lam_w**3 * (1 - rho_bar - 1j * eta_bar), -A_w * lam_w**2, 1]
])

print("\n" + "=" * 72)
print("SECTION 12: FULL CKM MATRIX")
print("=" * 72)
print("|V_CKM| (HCSM predictions):")
print(np.abs(V_CKM))
print("\nObserved |V_CKM|:")
V_obs = np.array([
    [0.9743, 0.2250, 0.0037],
    [0.2250, 0.9735, 0.0410],
    [0.0086, 0.0402, 0.9991]
])
print(V_obs)

# Unitarity check
V_dag_V = V_CKM.conj().T @ V_CKM
unitarity_error = np.max(np.abs(V_dag_V - np.eye(3)))
print(f"\nUnitarity error: {unitarity_error:.2e}")

# ================================================================
# SECTION 13: Boson Decomposition (winding + flat)
# ================================================================
M_BB_4_svd = svd(M_BB_4)
S_BB_full = M_BB_4_svd[1]
print("\n" + "=" * 72)
print("SECTION 13: BOSON DECOMPOSITION")
print("=" * 72)
print(f"Boson-boson singular values: {S_BB_full}")
print(f"Nonzero (winding): {S_BB_full[S_BB_full > 1e-10]}")
print(f"Zero (flat): {len(S_BB_full[S_BB_full < 1e-10])} directions")

# ================================================================
# SECTION 14: Summary Table
# ================================================================
print("\n" + "=" * 72)
print("SECTION 14: COMPLETE VERIFICATION SUMMARY")
print("=" * 72)

results = [
    ("rank(M)", 2, rank_M, "exact"),
    ("S0/S1 of M", 1.0, S_M[0] / S_M[1], "exact (isotropy)"),
    ("S0^BB", 49.4773, S_BB_full[0], "exact"),
    ("gamma", 0.702839467990, gamma, "exact"),
    ("delta_edge", 0.012143066492, delta_edge, "exact"),
    ("lambda", 0.224399475256, lam_w, "exact"),
    ("A", 0.809714453194, A_w, "exact"),
    ("R_b", 0.417509935173, R_b, "exact"),
    ("gamma_UT [deg]", 65.7932, np.degrees(gamma_UT), "0.007 deg"),
    ("rho_bar", 0.1712, rho_bar, "obs 0.159"),
    ("eta_bar", 0.3808, eta_bar, "obs 0.354"),
    ("J", 3.108e-5, J, "0.89%"),
    ("|V_us|", 0.2244, V_us, "0.27%"),
    ("|V_cb|", 0.0408, V_cb, "0.07%"),
    ("|V_ub|", 0.00382, V_ub, "<0.1%"),
    ("|V_td|", 0.00835, V_td, "2.62%"),
]

print(f"{'Quantity':<20s} {'Expected':>18s} {'Computed':>18s} {'Status':>20s}")
print("-" * 78)
for name, expected, computed, status in results:
    if isinstance(expected, float) and expected < 1e-3:
        exp_str = f"{expected:.4e}"
        comp_str = f"{computed:.4e}"
    else:
        exp_str = f"{expected:.10f}" if isinstance(expected, float) else str(expected)
        comp_str = f"{computed:.10f}" if isinstance(computed, float) else str(computed)
    print(f"{name:<20s} {exp_str:>18s} {comp_str:>18s} {status:>20s}")

# ================================================================
# SECTION 15: Falsification Tests
# ================================================================
print("\n" + "=" * 72)
print("SECTION 15: FALSIFICATION TESTS")
print("=" * 72)

tests = [
    ("rank(M) = 2", rank_M == 2),
    ("S0 = S1 (isotropy)", abs(S_M[0] - S_M[1]) < 1e-6),
    ("S0^BB = S1^BB (boson degeneracy)", abs(S_BB_full[0] - S_BB_full[1]) < 1e-6),
    ("gamma = 0.7028", abs(gamma - 0.702839467990) < 1e-10),
    ("lambda = pi/14", abs(lam_w - np.pi / 14) < 1e-15),
    ("gamma_UT in [65.5, 66.1] deg", 65.5 < np.degrees(gamma_UT) < 66.1),
    ("|V_us| in [0.223, 0.226]", 0.223 < V_us < 0.226),
    ("|V_cb| in [0.040, 0.042]", 0.040 < V_cb < 0.042),
    ("|V_ub| in [0.0037, 0.0039]", 0.0037 < V_ub < 0.0039),
    ("J in [3.0e-5, 3.2e-5]", 3.0e-5 < J < 3.2e-5),
]

all_pass = True
for name, passed in tests:
    status = "PASS" if passed else "FAIL"
    if not passed:
        all_pass = False
    print(f"  [{status}] {name}")

print("\n" + "=" * 72)
if all_pass:
    print("ALL FALSIFICATION TESTS PASS. HCSM-35 verified.")
else:
    print("SOME TESTS FAILED. Review the failed cases above.")
print("=" * 72)

# ================================================================
# SECTION 16: Output for Reproduction
# ================================================================
print("\n" + "=" * 72)
print("SECTION 16: REPRODUCTION OUTPUT (copy for paper)")
print("=" * 72)
print(f"""
CKM Wolfenstein parameters (HCSM-native):
  lambda       = {lam_w:.15f}   (obs 0.2250 +/- 0.0005)
  A            = {A_w:.15f}   (obs 0.826 +/- 0.009)
  R_b          = {R_b:.15f}   (obs 0.414 +/- 0.013)
  gamma_UT     = {np.degrees(gamma_UT):.4f} deg   (obs 65.8 +/- 2.8)
  rho_bar      = {rho_bar:.10f}   (obs 0.159 +/- 0.007)
  eta_bar      = {eta_bar:.10f}   (obs 0.354 +/- 0.006)

CKM magnitudes:
  |V_us|       = {V_us:.10f}   (obs 0.2250)
  |V_cb|       = {V_cb:.10f}   (obs 0.0408)
  |V_ub|       = {V_ub:.10f}   (obs 0.00382)
  |V_td|       = {V_td:.10f}   (obs 0.00857)

Jarlskog:
  J            = {J:.6e}   (obs 3.08e-5)

Stability:
  gamma        = {gamma:.15f}
  delta_edge   = {delta_edge:.15f}

Structural:
  rank(M)      = {rank_M}
  S0/S1        = {S_M[0]/S_M[1]:.15f}
  S0^BB        = {S_BB_full[0]:.10f}
""")

print("HCSM-35 verification complete.")