#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HCSM-21 Validation Script
=========================

Numerical verification of every claim in HCSM-21: The Time Theorem.

All computations are performed at 50-digit precision using mpmath.

Author: Stanley Preschutti
Date:   October 1, 2026
"""

import mpmath as mp

# ============================================================================
# SETUP
# ============================================================================

mp.mp.dps = 50  # 50 decimal places

def header(title):
    print("=" * 78)
    print(title)
    print("=" * 78)

def check(name, computed, expected, tol=mp.mpf('1e-40')):
    """Compare computed to expected; report match."""
    diff = abs(computed - expected)
    status = "OK " if diff < tol else "FAIL"
    print(f"  [{status}] {name}")
    print(f"        computed = {mp.nstr(computed, 40)}")
    print(f"        expected = {mp.nstr(expected, 40)}")
    if diff >= tol:
        print(f"        |diff|   = {mp.nstr(diff, 10)}")
    return diff < tol

# ============================================================================
# SECTION 1: STABILITY COEFFICIENT GAMMA
# ============================================================================

header("SECTION 1: Stability coefficient gamma")

# AGM(1, sqrt(2))
a = mp.mpf(1)
b = mp.sqrt(2)
for _ in range(200):
    a, b = (a + b) / 2, mp.sqrt(a * b)
AGM = a

# I(2) = 1 - (1/2) AGM(1, sqrt(2))
I2 = 1 - AGM / 2

# I(alpha) = (1/(2 alpha)) [B(a, 1/2) - 1/a], a = -1/2 + 1/(2 alpha)
def I(alpha):
    aa = -mp.mpf(1) / 2 + 1 / (2 * alpha)
    return (mp.beta(aa, mp.mpf(1) / 2) - 1 / aa) / (2 * alpha)

# c_2 = 2 I''(2) via central finite difference
h = mp.mpf('1e-25')
I_plus  = I(2 + h)
I_minus = I(2 - h)
I_double_prime = (I_plus - 2 * I2 + I_minus) / h**2
c2 = 2 * I_double_prime

gamma = c2 / I2

print(f"  AGM(1, sqrt(2)) = {mp.nstr(AGM, 40)}")
print(f"  I(2)            = {mp.nstr(I2, 40)}")
print(f"  c_2             = {mp.nstr(c2, 40)}")
print(f"  gamma           = {mp.nstr(gamma, 40)}")

check("gamma", gamma, mp.mpf('0.70283946799007245929257819650854315216110316215188'))

# ============================================================================
# SECTION 2: DERIVED CONSTANTS
# ============================================================================

header("SECTION 2: Derived constants")

T_DME = 1 / (2 * gamma)
delta_edge = mp.sqrt(2) / gamma - 2
w = -1 - (1 - gamma)**2 / (2 * mp.pi)
zeta_H_rho = (1 - gamma)**2 / (6 * mp.pi)
T_beat = mp.pi / mp.sqrt(2)
L_EFT = mp.exp(-zeta_H_rho * T_beat / 2)

M_P = mp.mpf('1.2209e19')
N = 64
d = 4
m_s = M_P / N * mp.exp(-(d + 1) * N / (d - 1))

print(f"  T_DME        = {mp.nstr(T_DME, 40)}")
print(f"  delta_edge   = {mp.nstr(delta_edge, 40)}")
print(f"  w            = {mp.nstr(w, 40)}")
print(f"  zeta H / rho = {mp.nstr(zeta_H_rho, 40)}")
print(f"  T_beat       = {mp.nstr(T_beat, 40)}")
print(f"  L_EFT        = {mp.nstr(L_EFT, 40)}")
print(f"  m_s          = {mp.nstr(m_s, 40)} GeV")

check("T_DME",      T_DME,      mp.mpf('0.71140000351696586214085618566222702281882679896640'))
check("delta_edge", delta_edge, mp.mpf('0.012143066491921425517407564301260126573205825055250'))
check("w",          w,          mp.mpf('-1.0140540788576648688014217438329986838970812421309'))
check("zetaH/rho",  zeta_H_rho, mp.mpf('0.0046846929525549563089897441470455395854206540938412'))
check("L_EFT",      L_EFT,      mp.mpf('0.99481012837355545749133011302003631191369417632946'))

# ============================================================================
# SECTION 3: LEVEL 1 - SUBSTRATE
# ============================================================================

header("SECTION 3: Level 1 - Substrate")

omega_A = 2 * mp.sqrt(2)
omega_B = 6 * mp.sqrt(2)

Phi_A = omega_A * T_DME
Phi_B = omega_B * T_DME

print(f"  omega_A = {mp.nstr(omega_A, 40)}")
print(f"  omega_B = {mp.nstr(omega_B, 40)}")
print(f"  Phi_A   = {mp.nstr(Phi_A, 40)}")
print(f"  Phi_B   = {mp.nstr(Phi_B, 40)}")

check("Phi_A", Phi_A, mp.mpf('2.0121430664919214255174075643012601265732058250552'))
check("Phi_B", Phi_B, mp.mpf('6.0364291994757642765522226929037803797196174751657'))
check("Phi_A - 2 = delta_edge", Phi_A - 2, delta_edge)
check("Phi_B - 6 = 3 delta_edge", Phi_B - 6, 3 * delta_edge)
check("Phi