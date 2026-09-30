<div align="center">

# HCSM — The Hodge Complex Standard Model

**A structural derivation of the Standard Model from a single discrete substrate**

[![ORCID](https://img.shields.io/badge/ORCID-0009--0004--5445--1744-green?logo=orcid)](https://orcid.org/0009-0004-5445-1744)
[![License](https://img.shields.io/badge/License-All%20Rights%20Reserved-red.svg)](#license)
[![Papers](https://img.shields.io/badge/Papers-58%20of%2060-orange.svg)](#repository-map)
[![Rank%20Blockers](https://img.shields.io/badge/Rank%20Blockers-7%2F7%20Closed-brightgreen.svg)](#rank-closure-status)
[![Anchors](https://img.shields.io/badge/Particle--Sector%20Anchors-Zero-brightgreen.svg)](#what-hcsm-derives)
[![Confirmed](https://img.shields.io/badge/Empirical%20Matches-7-brightgreen.svg)](#evidence--confirmed)
[![Predictions](https://img.shields.io/badge/Active%20Predictions-8-blue.svg)](#evidence--falsifiable)

**Physics is what the machine cannot cancel. The Hodge complex is the machine.**

Stanley Preschutti · Entropia Research Institute / Information Physics Institute

</div>

---

## TL;DR

- The **Standard Model** — gauge structure, matter content, particle spectrum, mixing matrices, spacetime dimension, graviton, Newton constant, and vacuum stability — is derived from a single discrete substrate: the periodic `8 × 8` torus `Λ = ℤ₈ × ℤ₈`.
- The substrate itself is selected by a **single theorem** (Unified Substrate Selection), from two framework identities: `N = d · 2ᵈ` and `mult_{λ=4}(L) = |F| + |B|`. The unique even-integer solution is `(L, d) = (8, 4)`, `N = 64`.
- The framework has **one dimensionful anchor** (`M_P`) and **zero empirical anchors in the particle sector**.
- **All seven rank blockers** of the derivation program are closed; the **need-to-derive list is empty**.
- **Seven derived quantities** already match observation at sub-percent precision. **Eight active predictions** are testable at upcoming facilities.

---

## What HCSM Derives

Every quantity below is derived from the substrate. No quantity is fitted.

| Domain | Result |
|---|---|
| **Substrate** | `(L, d) = (8, 4)`, `N = 64` — Unified Substrate Selection Theorem |
| **Spacetime** | `d = 4`; Lorentzian signature `(−, +, +, +)` from Hodge star on `ker D` |
| **Hodge complex** | `D = d + d*` on `Ω⁰ ⊕ Ω¹ ⊕ Ω²`; `D² = Δ` exactly; `dim = 256` |
| **Particle content** | 14 modes at `λ = 4`; three `D₄` orbits `{2, 8, 4}` |
| **Gauge algebra** | `su(3) ⊕ su(2) ⊕ u(1)⁵` from the `D₄` commutant |
| **Charge** | `q₃ = F(n₁, σ)`; bijection to `{−3, −1, 0, +2, +3}`; total charge `−1e` |
| **Mass** | `m = A · exp(−B · N/4)`, RMS `0.0449 dex` across 13 orders of magnitude |
| **Flavor** | `sin θ₁₂ = π/14`; Wolfenstein `√(ρ̄²+η̄²) = 5/12`; PMNS `(60, 60, 58)` |
| **Higgs** | `λ_tree = ln(kL)/[(d−1)π²]`; `v_EW = M_P exp(−kL_phys)`; `m_H = 125.14 GeV` |
| **Gravity** | `G_N = (d+1)² / (8π M_P²) = 25 / (8π M_P²)`; graviton exactly massless |
| **Cosmology** | `kL_bare = 192/5 = 38.4`; `ρ_DE = ((d+1)/2)(1 + cδ) m_ν⁴`, `c = 5/9`, `w = −1` |

---

## How It Works

Each level uses only the level above it. Every arrow is a theorem.

| Level | Object | Function | Result |
|:---:|---|---|---|
| **0** | Unified substrate selection | Single theorem | `(L, d) = (8, 4)`, `N = 64` |
| **1** | 8×8 torus `Λ = ℤ₈ × ℤ₈` | Substrate | `N = L² = 64`; 5-point Laplacian forced |
| **2** | Hodge complex `Ω⁰ ⊕ Ω¹ ⊕ Ω²` | Spectral arena | `D = d + d*`; `D² = Δ` exactly; `dim = 256` |
| **3** | 14-mode multiplet at `λ = 4` | Particle content | `D₄` orbits `{2, 8, 4}`; supports `{64, 48, 32}` |
| **4** | `D₄` representation theory | Gauge algebra | `su(3) ⊕ su(2) ⊕ u(1)⁵` |
| **5** | `ker D`, Betti `(1, 2, 1)` | Emergent spacetime | 4D Minkowski; signature `(−, +, +, +)` |
| **6** | `Sym²(Ω¹)|(λ=4)` | Gravity sector | Graviton `54B₁ ⊕ 54B₂`; exact masslessness |
| **7** | Hidden `B₁` sector | Normalization | `(d+1)² = 25` |
| **8** | Charge theorem | Electric charge | Bijection to five SM charge eigenvalues |
| **9** | `N/4` labels | Mass ordering | Mass map; RMS `0.0449 dex` |
| **10** | `D₄` breaking pattern | Flavor | CKM hierarchical in `π/14`; PMNS `(60, 60, 58)` |
| **11** | `kL` self-consistency | Cosmology | `ρ_DE`, `w = −1`, `S₀` |

---

## The Substrate Selection

The two identities that select the substrate are not new — they are the framework's own site-count relation and its own mode-content relation.

| Identity | Statement |
|---|---|
| **Site-count relation** | `N = d · 2ᵈ` |
| **Mode-content relation** | `mult_{λ=4}(L) = |F| + |B|` |

Combined:
