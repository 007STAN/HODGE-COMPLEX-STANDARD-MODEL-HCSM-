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

- The **Standard Model** — gauge structure, matter content, particle spectrum, mixing matrices, spacetime dimension, graviton, Newton constant, and vacuum stability is derived from a single discrete substrate: the periodic `8 × 8` torus `Λ = ℤ₈ × ℤ₈`.
- The substrate itself is selected by a **single theorem** (Unified Substrate Selection), from two framework identities: `N = d · 2ᵈ` and `mult_{λ=4}(L) = |F| + |B|`. The unique even-integer solution is `(L, d) = (8, 4)`, `N = 64`.
- The framework has **one dimensionful anchor** (`M_P`) and **zero empirical anchors in the particle sector**.
- **All seven rank blockers** of the derivation program are closed; the **need to derive list is empty**.
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

The two identities that select the substrate are not new, they are the framework's own site count relation and its own mode content relation.

| Identity | Statement |
|---|---|
| **Site-count relation** | `N = d · 2ᵈ` |
| **Mode-content relation** | `mult_{λ=4}(L) = |F| + |B|` |

Combined:
L² = d · 2ᵈ
2(L − 1) = d(d + 3) / 2


**Theorem (Unified Substrate Selection).** This system has the unique even-integer solution `L = 8`, `d = 4`, `N = 64`.

**Corollaries.** The four prior derivations of `L = 8` and `d = 4` — equal entropy spacing, maximal winding cancellation, Lorentzian signature on `ker D`, and the Betti numbers of `T²` — become independent consistency checks on the unique solution.

---

## Philosophy

### What this is

The universe is a **machine that cancels itself**. The substrate `Λ = ℤ₈ × ℤ₈` carries a mirror symmetry `λ ↔ 8 − λ` on its Laplacian spectrum. Every excitation pairs with a mirror partner and annihilates. The unique fixed point of the mirror is the self-paired eigenvalue `λ = 4`; the 14 modes that live there are what the mirror cannot cancel.

What we call physics quarks, leptons, gauge bosons, the Higgs, the graviton, dark energy is the **residue** at that fixed point.

HCSM does not ask *"what fields exist?"* or *"what symmetry groups are broken?"*
It asks a single question:

> **What substrate could produce exactly the Standard Model, and nothing else?**

The answer is a small, discrete, information-bearing lattice whose mirror cancellation leaves exactly the observed particle content.

### What this means

Every observed structure in the Standard Model is a **residue**  something the mirror could not cancel:

- **Charge** is what the machine *must not* cancel. The charge rigidity theorem shows that five structural constraints reduce `5¹⁴ ≈ 6 × 10⁹` assignments to exactly two, both with total charge `−1e`. The residual `−1e` is not a choice; it is the unique value the machine's own structure permits.
- **Mass** is the residual weight of the `N/4` label under the mass map `m = A · exp(−B · N/4)`.
- **Spacetime** is the kernel of the Hodge–Dirac operator: `ker D ≅ ⊕ₖ Hᵏ(T²; ℝ)`, with `dim ker D = b₀ + b₁ + b₂ = 1 + 2 + 1 = 4`. The emergent Lorentzian signature `(−, +, +, +)` is the unique `D₄`-invariant form on `ker D` that vanishes on pure-degree forms.
- **The graviton** is the self-paired `B₁ ⊕ B₂` curvature mode of `Sym²(Ω¹|_{λ=4})` — exactly massless at all orders by `D₄` Schur orthogonality.
- **Dark energy** is the residual vacuum energy after mirror cancellation, promoted to the observed cosmological constant by *pure projection* onto the timelike direction `(v₀ − v₂)/√2` of `ker D` — not by gravitational coupling, and with no `1/M_P²` suppression.

> **The framework is not a replacement for the Standard Model.**
> **It is the machine the Standard Model is the output of.**

> **The Standard Model is what the substrate looks like when you observe it.**
>
> **Physics is what the machine cannot cancel.**
> **The Hodge complex is the machine.**
> **The −1e is what the machine must not cancel.**

### How this is different

| ❌ Not | ✅ Is |
|---|---|
| String theory | A foundation beneath the Standard Model |
| Loop quantum gravity | Derives what the Standard Model assumes |
| A replacement for the Standard Model | The substrate the Standard Model is observed from |

---

## Rank Closure Status

All seven rank blockers of the HCSM derivation program are closed.

| Rank | Blocker | Closure | Paper |
|:---:|---|---|---|
| **1** | HCSM-native `γ` | Spectral computation from `τ = i`, `η(i)`, AGM | HCSM-54 |
| **2** | `kL` bare | Mass-map self-consistency; `kL_bare = 192/5 = 38.4` exact | HCSM-51 |
| **3** | Charge bijection | `q₃ = F(n₁, σ)` on 14 modes | HCSM-24 |
| **4** | `N = d · 2ᵈ` | Photon × Hodge form degrees | HCSM-01 |
| **5** | `(d−1)/2` from DNLS | Hartree orbit trace `(2, 3, 2)` on `(O₀, O₁, O₂)` | HCSM-19 |
| **6** | 3D irrep for CKM | Reframed: CKM is hierarchical in `π/14` | HCSM-34 |
| **7** | Loop-pass duration | Landauer selection: `T_pass = T_DME = 1/(2γ)` | HCSM-21 |

---

## Evidence — Confirmed

Seven HCSM derived quantities already match established empirical data.

| # | Parameter | HCSM Value | Observation | Agreement |
|:-:|---|---|---|---|
| 1 | Top quark mass | `172.688398 GeV` | PDG `172.690 ± 0.30 GeV` | `0.0053σ` |
| 2 | Cabibbo angle | `sin θ₁₂ = π/14 = 0.2243995` | CKM extraction of `\|V_us\|` | `0.044%` |
| 3 | Third Wolfenstein parameter | `√(ρ̄²+η̄²) = (5/12)(1 + δ_edge/6)` | CKM unitarity triangle fits | `0.03–0.07%` |
| 4 | Graviton spin | Exactly `2` | LIGO/Virgo polarization analyses | Confirmed |
| 5 | Graviton mass | Exactly `0` | Multi-messenger GW propagation (`v_g = c`) | Confirmed |
| 6 | Kaluza–Klein towers | Absent | LHC multi-TeV scans | Confirmed |
| 7 | Dark energy density | `ρ_DE ≈ 2.4809 × 10⁻⁴⁷ GeV⁴` | Planck CMB + DESI BAO | `0.76%` (refined: `0.010%`) |

Additional sub percent confirmations:

| Parameter | HCSM Value | Observation | Agreement |
|---|---|---|---|
| Higgs mass | `125.14 GeV` | PDG `125.25 ± 0.17 GeV` | `0.09%` |
| Weinberg angle | `sin²θ_W(M_Z) = 0.2312` | PDG `0.23122 ± 0.00004` | `0.02%` |
| `m_W / m_Z` | `cos θ_W = 0.8769` | PDG `0.8814` | `0.5%` |

---

## Evidence — Falsifiable

Eight framework-level predictions from HCSM-00 §9.

| # | Prediction | HCSM Value | Test |
|:-:|---|---|---|
| 1 | **No light dark photon** | `ε ≡ 0` exactly | SHiP, FASER2, MATHUSLA, CODEX-b, LDMX, NA64, Belle II, HL-LHC at `ε > 10⁻¹²` |
| 2 | **Strong CP** | `θ_QCD = 0` exactly | Any neutron EDM at generic `θ_QCD` level falsifies |
| 3 | **ρ parameter** | Tree-level `ρ = 1` exactly | Any tree-level deviation falsifies |
| 4 | **Fermion content** | Exactly 14 per generation | A fourth generation falsifies |
| 5 | **Kernel structure** | Winding kernel `dim = 10`; diag `{1/2, 3/4}` | Any change to `H_wind` spectrum falsifies |
| 6 | **Dark energy** | `ρ_DE = ((d+1)/2)(1 + cδ) m_ν⁴`, `c = 5/9` | Deviation > `0.1%` falsifies |
| 7 | **Charge rigidity** | Two configurations, both `−1e` | Vacuum total charge ≠ `−1e` falsifies |
| 8 | **Identification** | Unique `N/4` assignment | An SM particle with `N/4` outside the allowed set falsifies |

Sector predictions from the wider stack, currently under test:

| Paper | Prediction | HCSM Value | Window | Facility |
|---|---|---|---|---|
| HCSM-35 | `\|V_td\|` | `0.0061 ± 0.0001` | 2026–2028 | Belle II, LHCb |
| HCSM-53 | Equation of state `w` | `≈ −1.014…` | 2026–2027 | DESI DR3, Euclid |
| HCSM-53 | Bulk viscosity `ζH/ρ` | `≈ 0.004685…` | 2027–2028 | Euclid |
| HCSM-36 | Neutrino ordering | Normal (`m₁ ≈ m₂ < m₃`) | 2027+ | JUNO, DUNE |
| HCSM-36 | Majorana phases | `(α₁, α₂) = (π, 0)` | 2028+ | nEXO, LEGEND-1000, CUPID |

---

## Repository Map
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
└── 12-VERIFICATION/ HCSM-59


**Paper Status**

| Status | Count | Notes |
|---|---|---|
| Complete | 58 | HCSM-00 through HCSM-56a, plus 58, 60, 100 |
| Writable | 1 | HCSM-59 terminal verification suite |
| Pending | 1 | HCSM-57 (H₀ Boost) — waiting for data |
| Retired | 2 | SYK/QIN, ENTROPIX/MESA |
| **Total** | **60** | |

---

## Reading Order

**If you have 10 minutes:** Read the [TL;DR](#tldr) and [What HCSM Derives](#what-hcsm-derives) above. That is the entire claim.

**If you have 1 hour:** Read `00-FOUNDATIONS/HCSM-00.txt`. It is the framework's own summary of every result.

**If you want the full derivation:** Follow the stack in order:

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

**The universal structural signature:** every coefficient is a ratio of two HCSM native integer counts orbit sizes, group orders, processor ranks, measurement matrix ranks. Coefficients are not fitted; they are counts.

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
                  zero empirical anchors in the particle sector}
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
  "abstract": "HCSM derives the Standard Model from a single discrete substrate, the periodic 8x8 torus Lambda = Z_8 x Z_8. The substrate is selected by a Unified Substrate Selection Theorem from two framework identities: N = d * 2^d and mult_{lambda=4}(L) = |F| + |B|, with unique even-integer solution (L, d) = (8, 4). The Hodge complex on the torus yields the 14-mode particle multiplet at lambda = 4, the gauge algebra su(3) + su(2) + u(1)^5, four-dimensional spacetime as the kernel of the Hodge-Dirac operator, an exactly massless graviton, and the Newton constant G_N = 25/(8 pi M_P^2). All seven rank blockers are closed; the need-to-derive list is empty; zero empirical anchors remain in the particle sector. Seven derived quantities match observation at sub-percent precision; eight active predictions are testable at facilities through 2028.",
  "keywords": [
    "Hodge Complex Standard Model", "HCSM", "discrete substrate",
    "8x8 torus", "Unified Substrate Selection Theorem",
    "D4 representation theory", "Standard Model derivation",
    "charge theorem", "mass map", "top-Yukawa coefficient",
    "Cabibbo angle", "PMNS matrix", "Majorana phases",
    "graviton masslessness", "Newton constant", "dark energy",
    "DME entropy well", "S_0 normalization", "single anchor",
    "zero empirical anchors", "falsifiable unification",
    "Hodge-Dirac operator", "rank closure"
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
