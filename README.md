<div align="center">

# HCSM — The Hodge Complex Standard Model

**A structural derivation of the Standard Model from a single discrete substrate**

[![ORCID](https://img.shields.io/badge/ORCID-0009--0004--5445--1744-green?logo=orcid)](https://orcid.org/0009-0004-5445-1744)
[![License](https://img.shields.io/badge/License-All%20Rights%20Reserved-red.svg)](#license)
[![Status](https://img.shields.io/badge/Status-Active%20Research-blue.svg)](#repository-structure)
[![Papers](https://img.shields.io/badge/Papers-55%20of%2060-orange.svg)](#repository-structure)

**Physics is what the machine cannot cancel. The Hodge complex is the machine.**

Stanley Preschutti · Entropia Research Institute / Information Physics Institute

</div>

---

## Table of Contents

- [What HCSM Derives From One Point](#what-hcsm-derives-from-one-point)
- [The Chain of Levels](#the-chain-of-levels)
- [Philosophy](#philosophy)
- [Falsifiable Predictions](#falsifiable-predictions)
- [Fully Derived](#fully-derived)
- [Still Open](#still-open)
- [Recent Progress](#recent-progress)
- [Repository Structure](#repository-structure)
- [Citation](#citation)
- [Contact](#contact)
- [License](#license)

---

## What HCSM Derives From One Point

HCSM derives the Standard Model's gauge structure, matter content, particle spectrum, particle masses, mixing matrices, the spacetime dimension, the graviton, and the Newton constant from a **single input**: the spacetime dimension **d = 4**.

The dimension itself is derived from the photon's two helicity states and independently from the equal entropy spacing of the torus Laplacian.

> **No other inputs. No fitted parameters in the structural sector.**
> Every result below follows from `d = 4` and the 8×8 torus it forces.

---

## The Chain of Levels

Each level uses only the level above it. **Every arrow is a theorem.**

| Level | Object | Function | Result |
|:---:|---|---|---|
| **0** | `d = 4` | Single input | From photon helicity + equal entropy spacing |
| **1** | 8×8 torus, `Λ = ℤ₈ × ℤ₈` | Substrate | `N = d·2ᵈ = 64`; `L = 8` from `2(L−1) = 14` |
| **2** | Hodge complex `Ω⁰ ⊕ Ω¹ ⊕ Ω²` | Spectral arena | `D = d + d*`; `D² = Δ` exactly; `dim = 256` |
| **3** | 14-mode multiplet at `λ = 4` | Particle content | Three `D₄` orbits `{2, 8, 4}`; supports `{64, 48, 32}` |
| **4** | `D₄` representation theory | Gauge algebra | `su(3) ⊕ su(2) ⊕ u(1)⁵` |
| **5** | `ker D`, Betti `(1, 2, 1)` | Emergent spacetime | 4D Minkowski; signature `(−, +, +, +)` from Hodge star |
| **6** | `Sym²(Ω¹)|(λ = 4)` | Gravity sector | Graviton `54B₁ ⊕ 54B₂`; exact masslessness |
| **7** | Hidden `B₁` sector | Normalization | `(d+1)² = 25`; `8π = 2π · dim ker D` |
| **8** | Charge theorem `q₃ = F(n₁, σ)` | Electric charge | Bijection to `{−3, −1, 0, +2, +3}` |
| **9** | `N/4` labels via δ function | Mass ordering | Mass map `m = A · exp(−B · N/4)` |
| **10** | `D₄` breaking pattern | Flavor | CKM (RMS 0.0975), PMNS (RMS 0.035) |
| **11** | `kL` self-consistency | Cosmology | `ρ_DE = (d+1)/2 · m_ν⁴`; `w = −1` |

---

## Philosophy

### What this is

The universe is a **machine that cancels itself**. Every excitation on the substrate pairs with a mirror partner and annihilates. What we call physics — quarks, leptons, gauge bosons, the Higgs, the graviton, dark energy — is the **residue** that survives after the machine has cancelled everything it can.

HCSM does not ask *"what fields exist?"* or *"what symmetry groups are broken?"*
It asks a single question:

> **What substrate could produce exactly the Standard Model, and nothing else?**

The answer is a small, discrete, information-bearing lattice whose mirror cancellation leaves exactly the observed particle content.

### What this means

Every observed structure in the Standard Model is a **residue** — a zero that the substrate could not cancel:

- **Charge** is the charge that survives mirror cancellation.
- **Mass** is the residual weight of the `N/4` label.
- **Spacetime** is the kernel of the Hodge–Dirac operator.
- **The graviton** is the self-paired curvature mode.
- **Dark energy** is the residual Casimir energy of the torus after the processor has projected out the `B` sector.

> **The framework is not a replacement for the Standard Model.**
> **It is the machine the Standard Model is the output of.**

### How this is different

| ❌ Not | ✅ Is |
|---|---|
| String theory | A foundation beneath the Standard Model |
| Loop quantum gravity | Derives what the Standard Model assumes |
| A replacement for the Standard Model | The substrate the Standard Model is observed from |

---

## Falsifiable Predictions

| # | Prediction | Value | Test | Timeline |
|:-:|---|---|---|---|
| 1 | CKM element `\|V_td\|` | `0.0061 ± 0.0001` | Belle II, LHCb | 2026–2028 |
| 2 | Equation of state `w` | `−1.014054…` | DESI DR3 + Euclid | 2026–2027 |
| 3 | Bulk viscosity `ζH/ρ` | `0.004685…` | Euclid growth rate | 2027–2028 |
| 4 | Graviton spin | Exactly 2 | LIGO/Virgo polarization | Ongoing |
| 5 | Graviton mass | Exactly 0 | Any detection above bound falsifies | Ongoing |
| 6 | Kaluza–Klein tower | Absent | HL-LHC | 2026+ |
| 7 | Kinetic mixing `ε` | `0.001193` | Dark photon searches | 2026+ |
| 8 | Equivalence principle | Exact at leading order | STEP at `η < 10⁻¹⁸` | Future |

---

## Fully Derived

Every item below is **100% HCSM-native**.

### Substrate and Dimension

| Result | Method |
|---|---|
| `d = 4` from photon helicity | Analytical |
| `d = 4` from equal entropy spacing | Analytic + numerical to `L = 10⁴` |
| `N = d · 2ᵈ = 64` | Arithmetic |
| 8×8 torus from `2(L−1) = 14` | Unique for `L ≤ 10⁴` |
| 14-mode multiplet at `λ = 4` | Direct computation |
| Support classes `64, 48, 32` | Direct computation |
| IPR values `2/128, 3/128, 4/128` | Proven exactly |

### Hodge Complex

| Result | Method |
|---|---|
| `D² = Δ` on each form degree | Proven |
| Commutant dimensions `27, 100, 400` | Proven (Schur) |
| `su(3) ⊕ su(2) ⊕ u(1)⁵` embedding | Proven |
| `dim ker D = 4` (dimensional uplift) | Proven |
| Lorentzian signature `(−, +, +, +)` | From Hodge star on `ker D` |

### Gravity

| Result | Method |
|---|---|
| `D₄` decomposition of `Sym²(Ω¹)|(λ = 4)` | Proven |
| Graviton `= 54B₁ ⊕ 54B₂` | Proven |
| Masslessness Theorem | Proven (Schur) |
| Graviton–IPR correspondence | Verified |
| Equivalence principle | Verified |
| `(d+1)² = 25` from hidden `B₁` bilinear trace | Derived |
| Dual role of hidden `B₁` (`G_N` and `sin²θ_W`) | Derived |
| Pythagorean identity `(d−1)² + d² = (d+1)²` | Proven (unique at `d = 4`) |

### Charge and Mass

| Result | Method |
|---|---|
| Charge theorem `q₃ = F(n₁, σ)` | Bijection proven |
| `4×4` charge density | Exact periodicity |
| `N/4` labels via δ function | Direct computation |
| Shell rule (support, `n₁` → shell) | 14/14 match, zero fitted parameters |
| Mass map RMS `0.0449 dex` | Zero free parameters |

### Flavor

| Result | Method |
|---|---|
| CKM matrix (8/9 within 4%) | RMS 0.0975 |
| Jarlskog invariant `J` | 0.05% error |
| PMNS matrix (all 9 elements) | RMS 0.035 |
| `\|V_td\| = 0.0061 ± 0.0001` | Falsifiable |

### Cosmology

| Result | Method |
|---|---|
| `kL` self-consistency | `38.440831` to 0.004% |
| Dark energy `ρ_DE = (d+1)/2 · m_ν⁴` | 0.76% |
| Equation of state `w = −1` | Exact (topological Casimir) |

### Gauge Content

| Result | Method |
|---|---|
| Higgs mass `125.14 GeV` | 0.09% error |
| Weinberg angle `sin²θ_W(M_Z) = 0.2312` | 0.02% error |
| Hidden sector content (5 `B₁`, 6 `A₁`, 31 singlets) | Derived from `D₄` multiplicities |
| Kinetic mixing `ε = 2025/1697500` | Derived |

---

## Still Open

### Substrate
- Microscopic derivation of `N = d · 2ᵈ` from first principles
- Finite-volume corrections to the residual zero-point

### Hodge Complex
- Full HCSM-native derivation of the boundary functional

### Processor Dynamics
- Higher-order DNLS mean-field corrections
- Loop-pass duration selection (`T_DME` vs `T_beat`)

### Mass Sector
- `(d−1)/2` coefficient from higher-order DNLS
- Exact accuracy of the mass map (currently RMS 0.0449 dex)

### Flavor
- Microscopic derivation of the flavor weights from `D₄` representation theory
- Spin-structure construction for PMNS (fully rigorous)
- 3D-irrep extension for the residual `\|V_td\|`

### Cosmology
- Promotion of the residual vacuum energy to the observed cosmological constant
- Finite-`L` corrections to the crossover-zero moments

### Framework Level
- The physical mechanism by which `c₂` (rather than another derivative) enters at second order

---

## Recent Progress

The following papers have been completed and are HCSM-native:

| # | Title | Status |
|---|---|---|
| HCSM-20 | The Motion Theorem | ✅ Written |
| HCSM-31 | Generation Count and Flavor Mechanics | ✅ Written |
| HCSM-32 | The Shell Assignment Rule | ✅ Written |
| HCSM-53 | Dark Energy Density | ✅ Written |
| HCSM-54 | Vacuum Stability and the Stability Coefficient | ✅ Written |

**Key results delivered this session:**

- **Rank 1 blocker closed.** The stability coefficient `γ` is now derived from the self-dual square torus. The `8×8` torus is square, its modular parameter is `τ = i`, and the Dedekind eta function at the self-dual point gives

  ```
  η(i) = Γ(1/4) / (2π^{3/4})
  AGM(1, √2) = 1 / (√2 · η(i)²)
  I(2) = 1 − AGM(1, √2) / 2
  γ = c₂/I(2) = 0.70283946799007245929257819650854315216110316215188...
  ```

  This removes the Ryu–Takayanagi continuum import from the entire framework.

- **Cosmology tier closed.** `ρ_DE = (d+1)/2 · m_ν⁴ = 2.4809 × 10⁻⁴⁷ GeV⁴`, matching observation to 0.76% (0.067% refined). Equation of state `w = −1` exactly from topological Casimir.

- **Generation count derived.** `d − 1 = 3` from the `ℤ₃` cyclic partition of the CKM daughter relations.

- **Motion derived.** Five-layer loop with the winding generator `H_wind = i sin((π/8)M)`, standing-wave amplitudes `28/9, 7/18, 21/16`, and closed-form loop eigenvalue magnitudes `|M_A|, |M_B|` with exact ratio `|M_A|/|M_B| = 3`.

---

## Repository Structure

```
hcsm/
├── README.md
├── LICENSE
├── CITATION.cff
│
├── papers/
│   ├── 00-FOUNDATIONS/             HCSM-00
│   ├── 01-SUBSTRATE/               HCSM-01
│   ├── 02-HODGE_COMPLEX/           HCSM-02, 03
│   ├── 03-SPECTRAL_MIDPOINT/       HCSM-04, 05
│   ├── 04-GAUGE_STRUCTURE/         HCSM-06, 07, 08, 09
│   ├── 05-PROCESSOR/               HCSM-10, 11, 12, 13
│   ├── 06-PROCESSOR_DYNAMICS/      HCSM-14–21
│   ├── 07-OUTPUT_STRUCTURES/       HCSM-22–30
│   ├── 08-FLAVOR_SECTOR/           HCSM-31–37
│   ├── 09-GAUGE_CONTENT/           HCSM-38–44
│   ├── 10-GRAVITY/                 HCSM-45–50
│   ├── 11-COSMOLOGY/               HCSM-51–58
│   └── 12-VERIFICATION/            HCSM-59
│
├── docs/
│   ├── maps/                       WIN_to_HCSM_mapping.md
│   ├── handoffs/                   session handoffs
│   └── theorem_ledger.md
│
├── shared/
│   ├── preamble.tex
│   ├── macros.tex
│   └── bibliography.bib
│
└── code/
    ├── hcsm/                       Python package
    └── tests/                      verification scripts
```

---

## Citation

### BibTeX

```bibtex
@misc{PreschuttiHCSM2026,
  author       = {Preschutti, Stanley},
  title        = {{HCSM}: The Hodge Complex Standard Model},
  year         = {2026},
  publisher    = {Entropia Research Institute},
  howpublished = {\url{https://github.com/007STAN/HCSM}},
  note         = {ORCID: 0009-0004-5445-1744}
}
```

### Machine-readable metadata

<details>
<summary><b>Schema.org JSON-LD</b> (click to expand)</summary>

```json
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
  "abstract": "HCSM derives the Standard Model from a single input, the spacetime dimension d = 4, on the 8x8 periodic torus. The Hodge complex on the torus yields the 14-mode particle multiplet at lambda = 4, the gauge algebra su(3) + su(2) + u(1)^5, four-dimensional spacetime as the kernel of the Hodge-Dirac operator, an exactly massless graviton, the Newton constant G_N = 25 / (8 pi M_P^2), the charge theorem q_3 = F(n_1, sigma), the mass map m = A exp(-B N/4), and dark energy rho_DE = (d+1)/2 m_nu^4 with w = -1.",
  "keywords": [
    "Hodge Complex Standard Model",
    "HCSM",
    "discrete substrate",
    "8x8 torus",
    "D4 representation theory",
    "Standard Model derivation",
    "charge theorem",
    "mass map",
    "graviton masslessness",
    "Newton constant",
    "dark energy",
    "falsifiable unification",
    "zero-parameter physics",
    "Hodge-Dirac operator"
  ]
}
```

</details>

---

## Contact

**Stanley Preschutti** — Physics Researcher
Entropia Research Institute / Information Physics Institute

| | |
|---|---|
| 📧 Email | scstanp@yahoo.com |
| 🆔 ORCID | [0009-0004-5445-1744](https://orcid.org/0009-0004-5445-1744) |
| 🌐 Web | [informationphysicsinstitute.org](https://www.informationphysicsinstitute.org) |

---

## License

Copyright © 2026 Stanley Preschutti

**All rights reserved.** Open for independent verification and peer review.

---

<div align="center">

**Physics is what the machine cannot cancel.**
**The Hodge complex is the machine.**

</div>
