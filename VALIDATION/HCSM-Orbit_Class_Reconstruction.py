"""
HCSM-28 Validation Module
=========================

Validates the two central theorems of HCSM-28:

    Theorem 1 (Orbit-Class Reconstruction Theorem)
    Theorem 2 (Distinguished-Pair Graph Decomposition Theorem)

All claims are checked against HCSM-native definitions at machine
precision in IEEE-754 double precision.

Author:  Stanley Preschutti
ORCID:   0009-0004-5445-1744
Series:  HCSM White Paper Series, Paper 28
Date:    September 2026
"""

from itertools import combinations
from typing import FrozenSet, List, Tuple


# ---------------------------------------------------------------------------
# Section 2: Setup
# ---------------------------------------------------------------------------

MODE_ORDER: Tuple[str, ...] = (
    "nu", "e", "u", "d", "s", "mu", "dark", "c", "tau", "b",
    "W", "Z", "H", "t",
)

MODE_COORDS: dict = {
    "nu":   (0, 4),
    "e":    (4, 0),
    "u":    (1, 3),
    "d":    (1, 5),
    "s":    (3, 1),
    "mu":   (3, 7),
    "dark": (5, 7),
    "c":    (7, 3),
    "tau":  (7, 5),
    "b":    (5, 1),
    "W":    (2, 2),
    "Z":    (2, 6),
    "H":    (6, 2),
    "t":    (6, 6),
}

DIHEDRAL_ORBITS: Tuple[FrozenSet[str], ...] = (
    frozenset({"nu", "e"}),
    frozenset({"u", "d", "s", "mu", "dark", "c", "tau", "b"}),
    frozenset({"W", "Z", "H", "t"}),
)

WINDING_CLASSES: Tuple[FrozenSet[str], ...] = (
    frozenset({"nu", "e", "t"}),                       # A
    frozenset({"u", "d", "s", "mu", "c", "tau", "b"}), # B
    frozenset({"dark", "W", "Z", "H"}),                # C
)

FERMION_SECTOR: FrozenSet[str] = frozenset({
    "nu", "e", "u", "d", "s", "mu", "c", "tau", "b", "t",
})

BOSON_SECTOR: FrozenSet[str] = frozenset({"dark", "W", "Z", "H"})

CHARGE_LINES: dict = {
    -1: frozenset({"nu", "d", "s", "dark", "c", "Z"}),
     0: frozenset({"W", "t"}),
    +1: frozenset({"e", "mu", "tau", "b", "H", "u"}),
}

# Ten edges of the distinguished-pair graph (Section 2.7).
DISTINGUISHED_EDGES: Tuple[FrozenSet[str], ...] = (
    frozenset({"nu", "e"}),
    frozenset({"u", "s"}),
    frozenset({"u", "dark"}),
    frozenset({"d", "b"}),
    frozenset({"mu", "c"}),
    frozenset({"s", "tau"}),
    frozenset({"dark", "tau"}),
    frozenset({"dark", "t"}),
    frozenset({"Z", "H"}),
    frozenset({"W", "t"}),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mode_to_orbit() -> dict:
    out = {}
    for i, orbit in enumerate(DIHEDRAL_ORBITS):
        for m in orbit:
            out[m] = i
    return out


def _mode_to_class() -> dict:
    out = {}
    for i, cls in enumerate(WINDING_CLASSES):
        for m in cls:
            out[m] = i
    return out


def _mode_to_fb() -> dict:
    out = {}
    for m in FERMION_SECTOR:
        out[m] = "F"
    for m in BOSON_SECTOR:
        out[m] = "B"
    return out


def _mode_to_sigma() -> dict:
    out = {}
    for sigma, line in CHARGE_LINES.items():
        for m in line:
            out[m] = sigma
    return out


# ---------------------------------------------------------------------------
# Theorem 1: Orbit-Class Reconstruction
# ---------------------------------------------------------------------------

def check_lemma_3_1() -> bool:
    """Lemma 3.1: partitions agree on exactly 12 modes, differ on {t, dark}."""
    mo = _mode_to_orbit()
    mc = _mode_to_class()
    agreed, differ = [], []
    for m in MODE_ORDER:
        if mo[m] == mc[m]:
            agreed.append(m)
        else:
            differ.append(m)
    ok = len(agreed) == 12 and set(differ) == {"t", "dark"}
    print(f"Lemma 3.1: agreed={len(agreed)}, differ={sorted(differ)}  "
          f"[{'PASS' if ok else 'FAIL'}]")
    return ok


def check_theorem_3_3() -> bool:
    """Theorem 3.3: enumerate all (3, 7, 4) partitions that are
    F/B-homogeneous and agree with the dihedral orbit partition on
    at least 12 modes. Verify uniqueness."""
    mo = _mode_to_orbit()
    mfb = _mode_to_fb()

    # Step 1: C is forced to be the boson sector.
    C = BOSON_SECTOR
    fermionic_pool = sorted(FERMION_SECTOR)

    # Step 2: A must be a 3-subset of F, B = F \ A.
    # Step 3 + 4: count agreements with dihedral orbits.
    candidates: List[Tuple[FrozenSet[str], FrozenSet[str], FrozenSet[str]]] = []
    for A in combinations(fermionic_pool, 3):
        A_set = frozenset(A)
        B_set = FERMION_SECTOR - A_set
        # Compute agreements with dihedral correspondence:
        # A <-> O_0, B <-> O_1, C <-> O_2.
        agreements = 0
        for m in MODE_ORDER:
            orbit = mo[m]
            if orbit == 0 and m in A_set:
                agreements += 1
            elif orbit == 1 and m in B_set:
                agreements += 1
            elif orbit == 2 and m in C:
                agreements += 1
        if agreements >= 12:
            candidates.append((A_set, B_set, C))

    unique = (
        len(candidates) == 1
        and candidates[0][0] == frozenset({"nu", "e", "t"})
        and candidates[0][1] == frozenset({"u", "d", "s", "mu", "c", "tau", "b"})
        and candidates[0][2] == C
    )
    print(f"Theorem 3.3: {len(candidates)} candidate(s) satisfy both "
          f"hypotheses  [{'PASS' if unique else 'FAIL'}]")
    if candidates:
        A, B, _ = candidates[0]
        print(f"             A = {sorted(A)}")
        print(f"             B = {sorted(B)}")
    return unique


# ---------------------------------------------------------------------------
# Theorem 2: Distinguished-Pair Graph Decomposition
# ---------------------------------------------------------------------------

def _adjacency() -> dict:
    adj = {m: set() for m in MODE_ORDER}
    for edge in DISTINGUISHED_EDGES:
        a, b = tuple(edge)
        adj[a].add(b)
        adj[b].add(a)
    return adj


def _components(adj: dict) -> List[FrozenSet[str]]:
    visited = set()
    comps = []
    for m in MODE_ORDER:
        if m in visited:
            continue
        stack = [m]
        comp = set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            visited.add(x)
            stack.extend(adj[x] - comp)
        comps.append(frozenset(comp))
    return comps


def check_theorem_4_1_components() -> bool:
    """Theorem 4.1, Part 1: component decomposition."""
    adj = _adjacency()
    comps = _components(adj)
    iso = {c for c in comps if len(c) == 2}
    core = {c for c in comps if len(c) == 6}
    expected_iso = {
        frozenset({"nu", "e"}),
        frozenset({"d", "b"}),
        frozenset({"mu", "c"}),
        frozenset({"Z", "H"}),
    }
    expected_core = frozenset({"u", "s", "dark", "tau", "W", "t"})
    ok = (
        len(comps) == 5
        and iso == expected_iso
        and core == {expected_core}
    )
    print(f"Theorem 4.1 (components): {len(iso)} isolated edges, "
          f"{len(core)} six-vertex core  [{'PASS' if ok else 'FAIL'}]")
    return ok


def check_theorem_4_1_homogeneity() -> bool:
    """Theorem 4.1, Part 2: isolated-edge homogeneity."""
    mo = _mode_to_orbit()
    mc = _mode_to_class()
    mfb = _mode_to_fb()
    ms = _mode_to_sigma()

    iso_edges = [
        frozenset({"nu", "e"}),
        frozenset({"d", "b"}),
        frozenset({"mu", "c"}),
        frozenset({"Z", "H"}),
    ]
    ok = True
    for edge in iso_edges:
        a, b = tuple(edge)
        hom = (mo[a] == mo[b]) and (mc[a] == mc[b]) and (mfb[a] == mfb[b])
        het = (ms[a] != ms[b])
        if not (hom and het):
            ok = False
            print(f"             edge {sorted(edge)}: hom={hom}, het={het}")
    print(f"Theorem 4.1 (homogeneity): all isolated edges homogeneous in "
          f"orbit/class/F-B and heterogeneous in sigma  "
          f"[{'PASS' if ok else 'FAIL'}]")
    return ok


def check_theorem_4_1_core() -> bool:
    """Theorem 4.1, Part 3: core heterogeneity, degree, cycle rank."""
    mo = _mode_to_orbit()
    mc = _mode_to_class()
    mfb = _mode_to_fb()
    ms = _mode_to_sigma()

    core = {"u", "s", "dark", "tau", "W", "t"}
    core_edges = [
        frozenset({"u", "s"}),
        frozenset({"u", "dark"}),
        frozenset({"s", "tau"}),
        frozenset({"dark", "tau"}),
        frozenset({"dark", "t"}),
        frozenset({"W", "t"}),
    ]

    # Heterogeneity: orbit, class, F/B, sigma must all be non-constant.
    orbits = {mo[m] for m in core}
    classes = {mc[m] for m in core}
    fbs = {mfb[m] for m in core}
    sigmas = {ms[m] for m in core}
    het = (
        len(orbits) >= 2
        and len(classes) == 3
        and len(fbs) == 2
        and sigmas == {-1, 0, 1}
    )
    print(f"Theorem 4.1 (core heterogeneity): orbits={sorted(orbits)}, "
          f"classes={sorted(classes)}, F/B={sorted(fbs)}, "
          f"sigma={sorted(sigmas)}  [{'PASS' if het else 'FAIL'}]")

    # Degree distribution within the core.
    deg = {m: 0 for m in core}
    for edge in core_edges:
        a, b = tuple(edge)
        deg[a] += 1
        deg[b] += 1
    max_deg_vertices = [m for m, d in deg.items() if d == max(deg.values())]
    unique_hub = (max_deg_vertices == ["dark"] and deg["dark"] == 3)
    print(f"Theorem 4.1 (degree): degrees={deg}, unique hub = "
          f"{max_deg_vertices}  [{'PASS' if unique_hub else 'FAIL'}]")

    # Cycle rank of the core: E - V + 1.
    V = len(core)
    E = len(core_edges)
    rank = E - V + 1
    rank_ok = rank == 1
    print(f"Theorem 4.1 (cycle rank): V={V}, E={E}, rank={rank}  "
          f"[{'PASS' if rank_ok else 'FAIL'}]")

    # Explicit cycle.
    explicit_cycle = ["dark", "u", "s", "tau"]
    cycle_ok = True
    for i in range(len(explicit_cycle)):
        a = explicit_cycle[i]
        b = explicit_cycle[(i + 1) % len(explicit_cycle)]
        if frozenset({a, b}) not in {frozenset(e) for e in core_edges}:
            cycle_ok = False
    print(f"Theorem 4.1 (explicit cycle): dark -> u -> s -> tau -> dark  "
          f"[{'PASS' if cycle_ok else 'FAIL'}]")

    return het and unique_hub and rank_ok and cycle_ok


# ---------------------------------------------------------------------------
# Cross-checks
# ---------------------------------------------------------------------------

def check_graph_edge_count() -> bool:
    """Definition 2.2: the distinguished-pair graph has ten edges."""
    ok = len(DISTINGUISHED_EDGES) == 10
    print(f"Definition 2.2: {len(DISTINGUISHED_EDGES)} distinguished edges  "
          f"[{'PASS' if ok else 'FAIL'}]")
    return ok


def check_partition_consistency() -> bool:
    """Verify the dihedral, winding, and F/B partitions are well-formed
    partitions of the fourteen modes."""
    universe = set(MODE_ORDER)
    ok = True

    # Dihedral orbits partition the multiplet.
    if set().union(*DIHEDRAL_ORBITS) != universe:
        ok = False
    for a, b in combinations(DIHEDRAL_ORBITS, 2):
        if a & b:
            ok = False

    # Winding classes partition the multiplet.
    if set().union(*WINDING_CLASSES) != universe:
        ok = False
    for a, b in combinations(WINDING_CLASSES, 2):
        if a & b:
            ok = False

    # F/B split partitions the multiplet.
    if FERMION_SECTOR | BOSON_SECTOR != universe:
        ok = False
    if FERMION_SECTOR & BOSON_SECTOR:
        ok = False

    # Charge lines partition the multiplet.
    if set().union(*CHARGE_LINES.values()) != universe:
        ok = False
    for a, b in combinations(CHARGE_LINES.values(), 2):
        if a & b:
            ok = False

    print(f"Partition consistency: dihedral, winding, F/B, sigma lines all "
          f"partition the 14 modes  [{'PASS' if ok else 'FAIL'}]")
    return ok


def check_FB_correspondence() -> bool:
    """HCSM-13 correspondence: A union B = F, C = B_bos."""
    A, B, C = WINDING_CLASSES
    ok = (A | B == FERMION_SECTOR) and (C == BOSON_SECTOR)
    print(f"HCSM-13 correspondence: A union B = F, C = B_bos  "
          f"[{'PASS' if ok else 'FAIL'}]")
    return ok


# ---------------------------------------------------------------------------
# Top-level runner
# ---------------------------------------------------------------------------

def run_all() -> bool:
    print("=" * 70)
    print("HCSM-28 Validation Module")
    print("Orbit-Class Reconstruction")
    print("=" * 70)
    print()

    results = [
        ("Partition consistency",        check_partition_consistency()),
        ("HCSM-13 correspondence",       check_FB_correspondence()),
        ("Distinguished edge count",     check_graph_edge_count()),
        ("Lemma 3.1",                    check_lemma_3_1()),
        ("Theorem 3.3",                  check_theorem_3_3()),
        ("Theorem 4.1 (components)",     check_theorem_4_1_components()),
        ("Theorem 4.1 (homogeneity)",    check_theorem_4_1_homogeneity()),
        ("Theorem 4.1 (core)",           check_theorem_4_1_core()),
    ]

    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    n_pass = sum(1 for _, ok in results if ok)
    n_total = len(results)
    for name, ok in results:
        flag = "PASS" if ok else "FAIL"
        print(f"  [{flag}] {name}")
    print()
    print(f"  {n_pass} of {n_total} checks passed.")
    print("=" * 70)
    return n_pass == n_total


if __name__ == "__main__":
    import sys
    sys.exit(0 if run_all() else 1)