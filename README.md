<div align="center">

# HCSM — The Hodge Complex Standard Model

**A structural derivation of the Standard Model from a single discrete substrate**

[![ORCID](https://img.shields.io/badge/ORCID-0009--0004--5445--1744-green?logo=orcid)](https://orcid.org/0009-0004-5445-1744)
[![License](https://img.shields.io/badge/License-All%20Rights%20Reserved-red.svg)](#license)
[![Papers](https://img.shields.io/badge/Papers-64%20of%2066-orange.svg)](#repository-map)
[![Anchors](https://img.shields.io/badge/Particle--Sector%20Anchors-Zero-brightgreen.svg)](#what-hcsm-derives)
[![Confirmed](https://img.shields.io/badge/Empirical%20Matches-7-brightgreen.svg)](#evidence--confirmed)
[![Predictions](https://img.shields.io/badge/Active%20Predictions-8-blue.svg)](#evidence--falsifiable)

**Physics is what the mirror cannot cancel. The Hodge complex is the mirror.**

'If you cannot derive from 1st principles the input and the output it is NOT physics that is referred to as creative writing."

Stanley Preschutti · Entropia Research Institute / Information Physics Institute

</div>

---

## TL;DR

- The **Standard Model** — gauge structure, matter content, particle spectrum, mixing matrices, spacetime dimension, graviton, Newton constant, and vacuum stability — is derived from a single discrete substrate: the periodic `8 × 8` torus `Λ = ℤ₈ × ℤ₈`.
- The substrate itself is selected by a **single theorem** (Unified Substrate Selection), from two framework identities: `N = d · 2ᵈ` and `mult_{λ=4}(L) = |F| + |B|`. The unique even-integer solution is `(L, d) = (8, 4)`, `N = 64`.
- The framework has **one dimensionful anchor** (`M_P`) and **zero empirical anchors in the particle sector**.
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
| **Flavor** | `sin θ₁₂ = π/14`; Wolfenstein `√(ρ̃²+η̃²) = 5/12`; PMNS `(60, 60, 58)` |
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
L² = d · 2ᵈ
2(L − 1) = d(d + 3) / 2

**Theorem (Unified Substrate Selection).** This system has the unique even-integer solution `L = 8`, `d = 4`, `N = 64`.

**Corollaries.** The four prior derivations of `L = 8` and `d = 4` — equal entropy spacing, maximal winding cancellation, Lorentzian signature on `ker D`, and the Betti numbers of `T²` — become independent consistency checks on the unique solution.

### Evidence — Confirmed

Seven derived quantities already match established empirical data.

| # | Parameter | HCSM Value | Observation | Agreement |
|:-:|---|---|---|---|
| 1 | Top quark mass | `172.688398 GeV` | PDG `172.690 ± 0.30 GeV` | `0.0053σ` |
| 2 | Cabibbo angle | `sin θ₁₂ = π/14 = 0.2243995` | CKM extraction of `\|V_us\|` | `0.044%` |
| 3 | Third Wolfenstein parameter | `√(ρ̃²+η̃²) = (5/12)(1 + δ_edge/6)` | CKM unitarity triangle fits | `0.03–0.07%` |
| 4 | Graviton spin | Exactly `2` | LIGO/Virgo polarization analyses | Confirmed |
| 5 | Graviton mass | Exactly `0` | Multi-messenger GW propagation | Confirmed |
| 6 | Kaluza–Klein towers | Absent | LHC multi-TeV scans | Confirmed |
| 7 | Dark energy density | `ρ_DE ≈ 2.4809 × 10⁻⁴⁷ GeV⁴` | Planck CMB + DESI BAO | `0.76%` (refined `0.010%`) |

Additional sub-percent confirmations:

| Parameter | HCSM Value | Observation | Agreement |
|---|---|---|---|
| Higgs mass | `125.14 GeV` | PDG `125.25 ± 0.17 GeV` | `0.09%` |
| Weinberg angle | `sin²θ_W(M_Z) = 0.2312` | PDG `0.23122 ± 0.00004` | `0.02%` |
| `m_W / m_Z` | `cos θ_W = 0.8769` | PDG `0.8814` | `0.5%` |

### Evidence — Falsifiable

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

Sector predictions currently under test:

| Paper | Prediction | HCSM Value | Window | Facility |
|---|---|---|---|---|
| HCSM-35 | `\|V_td\|` | `0.0061 ± 0.0001` | 2026–2028 | Belle II, LHCb |
| HCSM-53 | Equation of state `w` | `≈ −1.014…` | 2026–2027 | DESI DR3, Euclid |
| HCSM-53 | Bulk viscosity `ζH/ρ` | `≈ 0.004685…` | 2027–2028 | Euclid |
| HCSM-36 | Neutrino ordering | Normal (`m₁ ≈ m₂ < m₃`) | 2027+ | JUNO, DUNE |
| HCSM-36 | Majorana phases | `(α₁, α₂) = (π, 0)` | 2028+ | nEXO, LEGEND-1000, CUPID |

---

## Program 2 — Quantum (Papers 61+)

The quantum program recovers the classical and quantum limits of the substrate. Each paper is a continuum-limit or quantization theorem applied to the derivational core.

### Papers 61–66 — Classical Limits and Quantization

| Paper | Title | Subject | Status |
|---|---|---|---|
| **HCSM-61** | Maxwell and Lorentz in the Classical Limit | Electromagnetic field equations; Lorentz covariance | Complete |
| **HCSM-62** | Newtonian Gravity as the Classical Limit | Poisson equation; geodesic equation; equivalence principle; `λ_L` | Complete |
| **HCSM-63** | First Post-Newtonian Corrections | EIH Lagrangian; PPN parameters; `γ = β = 1` | Complete |
| **HCSM-64** | Special Relativity as the Relativistic Limit | Four-momentum; mass-shell condition; `E² = p²c² + m²c⁴` | Complete |
| **HCSM-65** | The Dirac Equation from the Hodge–Dirac Operator | Kähler–Dirac correspondence; Clifford algebra; chirality | Complete |
| **HCSM-66** | Direct Quantization of the Hodge Complex | Hilbert space; Fock space; propagator; path integral | Complete |

### What the Quantum Program Derives

| Paper | Result |
|---|---|
| **HCSM-61** | `D = d + d*` → Maxwell equations; Lorentzian signature from Hodge star; `O(ε²)` lattice corrections; Tensor Product Theorem `V₅₆ ≅ V₁₄ ⊗ ker D`; Bundle Construction Theorem; explicit `S_YM^(2) = −(ε²/48)∫Tr(∂_λ F_μν ∂^λ F^μν)` |
| **HCSM-62** | `φ(r) = −G_N M / r`; `∇²φ = 4πG_N ρ_m`; `G_N = 25/(8πM_P²)`; `λ_L = 5M_P/1536`; geodesic equation; equivalence principle as theorem; speed of gravity `c_g = c_γ = c` |
| **HCSM-63** | 1PN metric; EIH Lagrangian; `γ = β = 1`; `α_i = ζ_i = ξ = 0`; `η_N = 0`; lattice corrections `O(ε²)` |
| **HCSM-64** | Four-momentum `p^μ = (E/c, p)`; mass-shell condition; `E² = p²c² + m²c⁴`; DNLS coupling `g = q = 33/64` |
| **HCSM-65** | `D/a → iγ^μ ∂_μ`; Clifford algebra `Cl(3,1)` from substrate operators `A₁, A₂`; chirality `γ₅ = −iA₁A₂`; mode-dependent Dirac equation `(D₅₆ + m_i)ψ_i = 0`; five Dirac bilinears; Majorana reduction |
| **HCSM-66** | `H₂₅₆` Hilbert space (dim 256); `D` as free Hamiltonian; full spectrum; `ℤ₂` grading → statistics; Fock space `F = Sym(H₊) ⊗ Λ(H₋)` with `dim Λ(H₋) = 2¹²⁸`; propagator `G(p) = (M(p) − iε)⁻¹`; unique path-integral measure; Ward identities |

### Quantum Program — Research Directions

| Folder | Program | Example Papers |
|---|---|---|
| **15-QUANTUM_GAUGE_SECTOR** | Gauge quantization, BRST, anomalies | HCSM-67–72 |
| **16-QUANTUM_SYMMETRY_BREAKING** | Quantum EWSB, chiral condensate | HCSM-73–75 |
| **17-QUANTUM_GRAVITY** | Graviton quantization, black holes | HCSM-76–80 |
| **18-QUANTUM_FOUNDATIONS** | Entanglement, Bell, measurement | HCSM-81–86 |
| **19-QUANTUM_STRUCTURE** | SUSY, topology, RG flow | HCSM-87–92 |
| **20-QUANTUM_PHENOMENOLOGY** | `g−2`, precision EW, neutrinos | HCSM-93–97 |

---

## Philosophy

### What this is

Most frameworks in theoretical physics take the Standard Model as given. They assume its gauge group, its particle content, and its parameters, then ask what lies beyond it — a larger symmetry, an extra dimension, a deeper structure. The Standard Model is the starting point.

HCSM takes the opposite direction. It does not postulate the Standard Model. It asks what substrate produces exactly this structure, and nothing else. The Standard Model is the output, not the input.

The answer is a small, discrete substrate: the periodic torus `Λ = ℤ₈ × ℤ₈`. On this substrate, the Laplacian spectrum carries a mirror symmetry `λ ↔ 8 − λ`. Every excitation pairs with a mirror partner and cancels. The unique fixed point of the mirror is the self-paired eigenvalue `λ = 4`. The 14 modes that live there are what the mirror cannot cancel.

What we call physics — quarks, leptons, gauge bosons, the Higgs, the graviton, dark energy — is the **residue** at that fixed point.

HCSM does not ask *"what fields exist?"* or *"what symmetry groups are broken?"* It asks:

> **What substrate produces exactly the Standard Model, and nothing else?**

### What this means

Every observed structure in the Standard Model is a **residue** — something the mirror could not cancel:

- **Charge** is what the mirror cannot cancel. The charge rigidity theorem shows that five structural constraints reduce `5¹⁴ ≈ 6 × 10⁹` assignments to exactly two, both with total charge `−1e`. The residual `−1e` is not a choice; it is the unique value the substrate's own structure permits.
- **Mass** is the residual weight of the `N/4` label under the mass map `m = A · exp(−B · N/4)`.
- **Spacetime** is the kernel of the Hodge–Dirac operator: `ker D ≅ ⊕ₖ Hᵏ(T²; ℝ)`, with `dim ker D = b₀ + b₁ + b₂ = 1 + 2 + 1 = 4`. The emergent Lorentzian signature `(−, +, +, +)` is the unique `D₄`-invariant form on `ker D` that vanishes on pure-degree forms.
- **The graviton** is the self-paired `B₁ ⊕ B₂` curvature mode of `Sym²(Ω¹|_{λ=4})` — exactly massless at all orders by `D₄` Schur orthogonality.
- **Dark energy** is the residual vacuum energy after mirror cancellation, promoted to the observed cosmological constant by *pure projection* onto the timelike direction `(v₀ − v₂)/√2` of `ker D` — not by gravitational coupling, and with no `1/M_P²` suppression.

> **The framework is not a replacement for the Standard Model.**
> **It is the substrate from which the Standard Model emerges.**

> **The Standard Model is what the substrate looks like when you observe it.**
>
> **Physics is what the mirror cannot cancel.**
> **The Hodge complex is the mirror.**
> **The −1e is what the mirror must not cancel.**

### How this is different

| ❌ Not | ✅ Is |
|---|---|
| A replacement for the Standard Model | The substrate from which the Standard Model emerges |
| A BSM extension (adds particles, symmetries, dimensions) | A pre-geometric framework (derives SM content) |
| A model that fits SM parameters | A substrate that derives SM parameters |

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
3. Future folders as the program expands.

**If you want to verify:** Run `12-VERIFICATION/run_all.py`. It reproduces every numerical claim in the stack in IEEE-754 double precision.

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
                  unified substrate selection theorem;
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
    "email": "scstannp@yahoo.com",
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
  "abstract": "HCSM derives the Standard Model from a single discrete substrate, the periodic 8x8 torus Lambda = Z_8 x Z_8. The substrate is selected by a Unified Substrate Selection Theorem from two framework identities: N = d * 2^d and mult_{lambda=4}(L) = |F| + |B|, with unique even-integer solution (L, d) = (8, 4). The Hodge complex on the torus yields the 14-mode particle multiplet at lambda = 4, the gauge algebra su(3) + su(2) + u(1)^5, four-dimensional spacetime as the kernel of the Hodge-Dirac operator, an exactly massless graviton, and the Newton constant G_N = 25/(8 pi M_P^2). Zero empirical anchors remain in the particle sector. The paper stack divides into a derivational program (papers 00-60) and a quantum program (papers 61+). The derivational program fixes the structure of the framework. The quantum program recovers Maxwell, Newton, Dirac, and quantum field theory as continuum and quantum limits of the same substrate. Seven derived quantities match observation at sub-percent precision; eight active predictions are testable at facilities through 2028.",
  "keywords": [
    "Hodge Complex Standard Model", "HCSM", "discrete substrate",
    "8x8 torus", "Unified Substrate Selection Theorem",
    "D4 representation theory", "Standard Model derivation",
    "charge theorem", "mass map", "top-Yukawa coefficient",
    "Cabibbo angle", "PMNS matrix", "Majorana phases",
    "graviton masslessness", "Newton constant", "dark energy",
    "DME entropy well", "S_0 normalization", "single anchor",
    "zero empirical anchors", "falsifiable unification",
    "Hodge-Dirac operator",
    "quantum field theory", "Dirac equation", "Kähler-Dirac operator",
    "direct quantization", "Fock space", "path integral",
    "Maxwell equations", "Lorentz covariance", "tensor product theorem",
    "bundle construction", "O(epsilon^2) corrections"
  ]
}
</details>
Contact

Stanley Preschutti — Physics Researcher
Entropia Research Institute / Information Physics Institute

📧 Email	scstannp@yahoo.com
🆔 ORCID	0009-0004-5445-1744
🌐 Web	informationphysicsinstitute.org
License

Copyright © 2026 Stanley Preschutti

All rights reserved. Open for independent verification and peer review.
