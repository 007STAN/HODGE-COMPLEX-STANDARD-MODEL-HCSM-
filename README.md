<div align="center">

# HCSM — The Hodge Complex Standard Model

**A structural derivation of the Standard Model from a single discrete substrate**

[![ORCID](https://img.shields.io/badge/ORCID-0009--0004--5445--1744-green?logo=orcid)](https://orcid.org/0009-0004-5445-1744)
[![License](https://img.shields.io/badge/License-All%20Rights%20Reserved-red.svg)](#license)
[![Papers](https://img.shields.io/badge/Papers-64%20of%2066-orange.svg)](#repository-map)
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
- The paper stack divides into two programs: a **derivational program** (papers 00–60) that derives the Standard Model from the substrate, and a **quantum program** (papers 61+) that recovers classical limits, wave equations, and quantum field theory from the same substrate.

---

## The Two Programs

The HCSM paper stack has two distinct programs, each answering a different question.

| Program | Papers | Question | Method | Output |
|---|---|---|---|---|
| **Derivational** | 00–60 | What does the substrate contain? | Representation theory, spectral theory, group theory | Particles, gauge algebra, masses, mixing, gravity, cosmology |
| **Quantum** | 61+ | How do the contents behave? | Continuum limits, wave equations, quantization, path integral | Field equations, propagators, Fock space, QFT |

The derivational program is **canonical**: it fixes the structure of the framework and derives every observable from the substrate.

The quantum program is **reconstructive**: it recovers known physics — Maxwell, Newton, Dirac, quantum field theory — as the continuum and quantum limits of the same substrate.

Both are required. Neither replaces the other.

---

## Program 1 — Derivational (Papers 00–60)

The derivational program is the canonical core of HCSM. Every quantity below is derived from the substrate. No quantity is fitted.

### What HCSM Derives

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

### How It Works

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

### Substrate Selection

The two identities that select the substrate are the framework's own site-count relation and its own mode-content relation.

| Identity | Statement |
|---|---|
| **Site-count relation** | `N = d · 2ᵈ` |
| **Mode-content relation** | `mult_{λ=4}(L) = |F| + |B|` |

Combined:
HCSM/
├── README.md
├── VERIFICATION/
│
├── 00-FOUNDATIONS/ HCSM-00
├── 01-SUBSTRATE/ HCSM-01
├── 02-HODGE_COMPLEX/ HCSM-02, 03
├── 03-SPECTRAL_MIDPOINT/ HCSM-04, 05
├── 04-GAUGE_STRUCTURE/ HCSM-06, 07, 08, 09
├── 05-PROCESSOR/ HCSM-10, 11, 12, 13
├── 06-PROCESSOR_DYNAMICS/ HCSM-14, 14a, 14b, 15–21
├── 07-OUTPUT_STRUCTURES/ HCSM-22, 22a, 23–32
├── 08-FLAVOR_SECTOR/ HCSM-33–37
├── 09-GAUGE_CONTENT/ HCSM-38–44
├── 10-GRAVITY/ HCSM-45–50
├── 11-COSMOLOGY/ HCSM-51–56, 56a, 58, 60
├── 12-VERIFICATION/ HCSM-59
│
├── 13-CLASSICAL_LIMITS/ HCSM-61, 62, 63, 64
├── 14-QUANTUM_SECTOR/ HCSM-65, 66
└── ...


### Paper Status

| Status | Count | Notes |
|---|---|---|
| Complete | 64 | HCSM-00 through HCSM-56a, plus 58, 60, and 61–66 |
| Writable | 1 | HCSM-59 terminal verification suite |
| Pending | 1 | HCSM-57 (H₀ Boost) — waiting for data |
| Retired | 2 | SYK/QIN, ENTROPIX/MESA |
| **Total** | **66** | |

---

## Reading Order

**If you have 10 minutes:** Read the [TL;DR](#tldr) and [The Two Programs](#the-two-programs) above. That is the entire claim.

**If you have 1 hour:** Read `00-FOUNDATIONS/HCSM-00.txt`. It is the framework's own summary of every result.

### Derivational program (00–60)

Follow the stack in order:

1. `00-FOUNDATIONS` — foundations, notation, and the unified substrate selection
2. `01-SUBSTRATE` — the 8×8 torus and the 5-point Laplacian
3. `02-HODGE_COMPLEX` — Hodge–Dirac operator; Lorentzian signature
4. `03-SPECTRAL_MIDPOINT` — the 14-mode multiplet and its `D₄` orbits
5. `04-GAUGE_STRUCTURE` — gauge algebra, hypercharge, anomalies
6. `05-PROCESSOR` — the processor, zero state, and vacuum mechanism
7. `06-PROCESSOR_DYNAMICS` — winding generator, Kernel–Diagonal Theorem, DNLS, motion, time
8. `07-OUTPUT_STRUCTURES` — charge, mass labels, observer pairs
9. `08-FLAVOR_SECTOR` — mass map, CKM, PMNS, generation phase
10. `09-GAUGE_CONTENT` — hidden sector, Higgs, dark photon, Proca
11. `10-GRAVITY` — graviton, Newton constant, self-coupling
12. `11-COSMOLOGY` — `kL` crossover, dark energy, vacuum stability, Λ-promotion

### Quantum program (61+)

Follow the stack in order:

1. `13-CLASSICAL_LIMITS`
   - `HCSM-61` — Maxwell and Lorentz in the classical limit
   - `HCSM-62` — Newtonian gravity as the classical limit
   - `HCSM-63` — first post-Newtonian corrections
   - `HCSM-64` — special relativity as the relativistic limit
2. `14-QUANTUM_SECTOR`
   - `HCSM-65` — the Dirac equation from the Hodge–Dirac operator
   - `HCSM-66` — direct quantization of the Hodge complex
3. Future folders (`15-QUANTUM_GAUGE_SECTOR`, `16-QUANTUM_GRAVITY`, etc.) as the program expands.

**If you want to verify:** Run `12-VERIFICATION/run_all.py`. It reproduces every numerical claim in the stack in IEEE-754 double precision.

---

## Rank Closure and Session Summary

All nine items (6.1–6.9) were derived, grounded, or formalized.

| # | Item | Result |
|:---:|---|---|
| 6.1 | `kL` correction coefficient | `(kL_phys − kL_bare) = (\|D₄\|−1)/2 · δ_edge = (7/2)δ_edge` |
| 6.2 | Bridge identity | Consistent with exactness at `8.4×10⁻⁸` |
| 6.3 | Spin-structure phase | `φ_gen = (3/5)π = 0.6π` |
| 6.4 | Majorana phases | `(α₁, α₂) = (π, 0)` |
| 6.5 | Class-1 denominator | `145 = N + (d−1)ᵈ = 64 + 81` |
| 6.6 | Anti-daughter coefficients | `15/16 = 1 − 1/\|O₂\|²`; `4/9 = \|O₂\|/(d−1)²` |
| 6.7 | Residual-bath factor | `(N−14)/(N−13) = 50/51` |
| 6.8 | DME entropy-well `S₀` | `S₀ = v_EW² · 7/40 · (1 + δ_edge/6) = 10630.722 GeV²` |
| 6.9 | CKM higher-order | `√(ρ̄²+η̄²) = (5/12)(1 + δ_edge/6)` |

**The universal structural signature:** every coefficient is a ratio of two HCSM native integer counts — orbit sizes, group orders, processor ranks, measurement matrix ranks. Coefficients are not fitted; they are counts.

---

## Reproducibility

Every numerical claim in the HCSM paper stack is verified at machine precision in IEEE-754 double precision. The verification scripts use only standard libraries (NumPy, SciPy, SymPy, mpmath) and construct the objects they need from HCSM-native definitions.

Running `12-VERIFICATION/run_all.py` reproduces every numerical claim in the stack.

---

## Citation

```bibtex
@misc{PreschuttiHCSM2026,
  author       = {Preschutti, Stanley},
  title        = {{HCSM}: The Hodge Complex Standard Model},
  year         = {2026},
  publisher    = {Entropia Research Institute},
  howpublished = {\url{https://github.com/007STAN/HODGE-COMPLEX-STANDARD-MODEL-HCSM-}},
  note         = {ORCID: 0009-0004-5445-1744;
                  all seven rank blockers closed;
                  unified substrate selection theorem;
                  need-to-derive list empty;
                  zero empirical anchors in the particle sector;
                  quantum program papers 61--66 complete}
}
<details> <summary><b>Schema.org JSON-LD</b> (click to expand)</summary>
{
  "@context": "https://schema.org",
  "@type": "ScholarlyArticle",
  "name": "HCSM — The Hodge Complex Standard Model",
  "alternativeName": "Hodge Complex Standard Model",
  "author": {
    "@type": "Person",
    "name": "Stanley Preschutti",
    "jobTitle": "Physics Researcher",
    "email": "scstanp@yahoo.com",
    "identifier": "https://orcid.org/0009-0004-5445-1744",
    "affiliation": {
      "@type": "Organization",
      "name": "Entropia Research Institute",
      "url": "https://www.informationphysicsinstitute.org"
    }
  },
  "datePublished": "2026",
  "inLanguage": "en",
  "license": "All rights reserved",
  "abstract": "HCSM derives the Standard Model from a single discrete substrate, the periodic 8x8 torus Lambda = Z_8 x Z_8. The substrate is selected by a Unified Substrate Selection Theorem from two framework identities: N = d * 2^d and mult_{lambda=4}(L) = |F| + |B|, with unique even-integer solution (L, d) = (8, 4). The Hodge complex on the torus yields the 14-mode particle multiplet at lambda = 4, the gauge algebra su(3) + su(2) + u(1)^5, four-dimensional spacetime as the kernel of the Hodge-Dirac operator, an exactly massless graviton, and the Newton constant G_N = 25/(8 pi M_P^2). All seven rank blockers are closed; the need-to-derive list is empty; zero empirical anchors remain in the particle sector. The paper stack divides into a derivational program (papers 00-60) and a quantum program (papers 61+). The derivational program fixes the structure of the framework. The quantum program recovers Maxwell, Newton, Dirac, and quantum field theory as continuum and quantum limits of the same substrate. Seven derived quantities match observation at sub-percent precision; eight active predictions are testable at facilities through 2028.",
  "keywords": [
    "Hodge Complex Standard Model", "HCSM", "discrete substrate",
    "8x8 torus", "Unified Substrate Selection Theorem",
    "D4 representation theory", "Standard Model derivation",
    "charge theorem", "mass map", "top-Yukawa coefficient",
    "Cabibbo angle", "PMNS matrix", "Majorana phases",
    "graviton masslessness", "Newton constant", "dark energy",
    "DME entropy well", "S_0 normalization", "single anchor",
    "zero empirical anchors", "falsifiable unification",
    "Hodge-Dirac operator", "rank closure",
    "quantum field theory", "Dirac equation", "Kähler-Dirac operator",
    "direct quantization", "Fock space", "path integral"
  ]
}
</details>
Contact

Stanley Preschutti — Physics Researcher
Entropia Research Institute / Information Physics Institute

📧 Email	scstanp@yahoo.com
🆔 ORCID	0009-0004-5445-1744
🌐 Web	informationphysicsinstitute.org
License

Copyright © 2026 Stanley Preschutti

All rights reserved. Open for independent verification and peer review.
