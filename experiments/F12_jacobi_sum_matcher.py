#!/usr/bin/env python3
"""
F12: Jacobi Sum Matcher

Computes Gauss sums g(chi^k) and Jacobi sums J(chi^a, chi^b) for characters
on F_p*, then matches them against eigenvalues from the Witt carry matrices
(untwisted and character-twisted from F9).

Key theoretical point: the Fermat curve X^p + Y^p = 1 degenerates in
characteristic p (becoming (X+Y)^p = 1, i.e. a line). The relevant variety
is NOT the Fermat curve itself but the Artin-Schreier variety z^p - z = c_2
(explored in F11). Nevertheless, the Jacobi/Gauss sums of F_p* are the
natural candidates for eigenvalues at the Weil bound.

Gauss sum: g(chi^k) = sum_{x=1}^{p-1} chi^k(x) * omega^x
           where chi is a character of order m=p-1, omega = exp(2*pi*i/p)

Jacobi sum: J(chi^a, chi^b) = sum_{x in F_p} chi^a(x) * chi^b(1-x)
            Relation: J(chi^a, chi^b) = g(chi^a) * g(chi^b) / g(chi^{a+b})
            when chi^a, chi^b, chi^{a+b} are all non-trivial.

This script:
  1. Computes all Gauss sums and non-trivial Jacobi sums for p = 3, 5, 7
  2. Computes the Witt carry eigenvalues (untwisted + all twists)
  3. Matches carry eigenvalues against the full Gauss/Jacobi catalog
  4. Reports which sums the carry operator "sees"
"""

import math
import cmath
from itertools import product as cartesian_product
import numpy as np

PI = math.pi


# =====================================================================
# UTILITY FUNCTIONS (from F8c/F9)
# =====================================================================

def witt_mul_exact(x, y, p, n):
    x_int, y_int = list(x), list(y)
    def ghost(a, k):
        return sum(p**i * a[i]**(p**(k-i)) for i in range(k+1))
    ghost_prod = [ghost(x_int, k) * ghost(y_int, k) for k in range(n)]
    s = []
    for k in range(n):
        rhs = ghost_prod[k] - sum(p**i * s[i]**(p**(k-i)) for i in range(k))
        s.append((rhs // p**k) % p)
    return tuple(s)


def build_level2_exponent_matrix(p):
    n = 3
    states = list(cartesian_product(range(p), repeat=2))
    state_idx = {s: i for i, s in enumerate(states)}
    D = p * p
    E = np.zeros((D, D), dtype=int)
    for x0 in range(p):
        for y0 in range(p):
            for x1 in range(p):
                for y1 in range(p):
                    x = (x0, x1, 0)
                    y = (y0, y1, 0)
                    pr = witt_mul_exact(x, y, p, n)
                    g2x = x0**(p**2) + p * x1**p
                    g2y = y0**(p**2) + p * y1**p
                    g2_prod = g2x * g2y
                    partial = g2_prod - pr[0]**(p**2) - p * pr[1]**p
                    c2 = (partial // (p**2)) % p
                    i_from = state_idx[(x0, y0)]
                    j_to = state_idx[(x1, y1)]
                    E[j_to, i_from] = c2
    return E, states, state_idx


def primitive_root(p):
    for g in range(2, p):
        order = 1
        val = g % p
        while val != 1:
            val = (val * g) % p
            order += 1
        if order == p - 1:
            return g
    return None


def build_character_table(p):
    g = primitive_root(p)
    zeta = np.exp(2j * PI / (p - 1))
    m = p - 1
    dlog = {}
    val = 1
    for j in range(m):
        dlog[val] = j
        val = (val * g) % p
    chars = {}
    for k in range(m):
        chi = {}
        chi[0] = 0
        for a in range(1, p):
            chi[a] = zeta ** (k * dlog[a])
        chars[k] = chi
    return chars, g, dlog


# =====================================================================
# GAUSS AND JACOBI SUM COMPUTATION
# =====================================================================

def compute_gauss_sums(p):
    """Compute Gauss sums g(chi^k) = sum_{x=1}^{p-1} chi^k(x) omega^x.

    Returns dict: k -> g(chi^k), for k = 1, ..., p-2 (non-trivial characters).
    """
    g_root = primitive_root(p)
    m = p - 1
    omega = np.exp(2j * PI / p)
    zeta_m = np.exp(2j * PI / m)

    dlog = {}
    val = 1
    for j in range(m):
        dlog[val] = j
        val = (val * g_root) % p

    gauss = {}
    for k in range(1, m):
        G = sum(zeta_m ** (k * dlog[x]) * omega ** x for x in range(1, p))
        gauss[k] = G

    return gauss


def compute_jacobi_sums(p):
    """Compute all non-trivial Jacobi sums J(chi^a, chi^b).

    Valid pairs: 1 <= a, b <= p-2, and a+b not divisible by p-1.
    Returns dict: (a, b) -> J(chi^a, chi^b).
    """
    g_root = primitive_root(p)
    m = p - 1
    zeta_m = np.exp(2j * PI / m)

    dlog = {}
    val = 1
    for j in range(m):
        dlog[val] = j
        val = (val * g_root) % p

    def chi_power(a, x):
        if x == 0:
            return 0.0
        return zeta_m ** (a * dlog[x])

    jacobi = {}
    for a in range(1, m):
        for b in range(1, m):
            if (a + b) % m == 0:
                continue
            J = sum(chi_power(a, x) * chi_power(b, (1 - x) % p) for x in range(p))
            jacobi[(a, b)] = J

    return jacobi


def compute_witt_eigenvalues_all_twists(p):
    """Compute Witt carry eigenvalues for untwisted + all multiplicative/additive twists.

    Returns list of (eigenvalue_array, label_string).
    """
    E, states, state_idx = build_level2_exponent_matrix(p)
    omega = np.exp(2j * PI / p)
    D = p * p
    chars, g, dlog = build_character_table(p)
    results = []

    # Untwisted
    M_base = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            M_base[i, j] = omega ** E[i, j]
    results.append((np.linalg.eigvals(M_base), "untwisted"))

    # Multiplicative twists
    for k in range(1, p - 1):
        chi = chars[k]
        order = (p - 1) // math.gcd(k, p - 1)
        M = np.zeros((D, D), dtype=complex)
        for x0 in range(p):
            for y0 in range(p):
                chi_val = chi[(x0 * y0) % p]
                i_from = state_idx[(x0, y0)]
                for x1 in range(p):
                    for y1 in range(p):
                        j_to = state_idx[(x1, y1)]
                        M[j_to, i_from] = chi_val * omega ** E[j_to, i_from]
        results.append((np.linalg.eigvals(M), f"mul chi_{k} (ord {order})"))

    # Additive twists
    for k in range(1, p - 1):
        chi = chars[k]
        order = (p - 1) // math.gcd(k, p - 1)
        M = np.zeros((D, D), dtype=complex)
        for x0 in range(p):
            for y0 in range(p):
                chi_val = chi[(x0 + y0) % p]
                i_from = state_idx[(x0, y0)]
                for x1 in range(p):
                    for y1 in range(p):
                        j_to = state_idx[(x1, y1)]
                        M[j_to, i_from] = chi_val * omega ** E[j_to, i_from]
        results.append((np.linalg.eigvals(M), f"add chi_{k} (ord {order})"))

    return results


def extract_weil(eigs, p, tol_frac=0.005):
    sqrt_p = math.sqrt(p)
    return [e for e in eigs if abs(abs(e) - sqrt_p) < tol_frac * sqrt_p and abs(e) > 1e-6]


# =====================================================================
# MATCHING ENGINE
# =====================================================================

def build_candidate_catalog(gauss, jacobi, p):
    """Build a catalog of all candidate complex numbers at |z| = sqrt(p).

    Includes: Gauss sums, Jacobi sums, negated versions, and conjugates.
    """
    sqrt_p = math.sqrt(p)
    catalog = {}

    for k, G in gauss.items():
        if abs(abs(G) - sqrt_p) < 0.01 * sqrt_p:
            catalog[f"g(chi^{k})"] = G
            catalog[f"-g(chi^{k})"] = -G

    for (a, b), J in jacobi.items():
        if abs(abs(J) - sqrt_p) < 0.01 * sqrt_p:
            catalog[f"J(chi^{a},chi^{b})"] = J

    return catalog


# =====================================================================
# MAIN
# =====================================================================

print("=" * 78)
print("  F12: JACOBI SUM MATCHER")
print("  Matching Witt carry eigenvalues against Gauss/Jacobi sums")
print("=" * 78)

for p in [3, 5, 7]:
    print(f"\n{'='*78}")
    print(f"  PRIME p = {p}")
    m = p - 1
    g_genus = (p - 1) * (p - 2) // 2
    sqrt_p = math.sqrt(p)
    print(f"  |F_p*| = {m}, primitive root = {primitive_root(p)}")
    print(f"  Fermat genus (smooth, char != p): g = {g_genus}")
    print(f"  NOTE: X^p + Y^p = 1 is degenerate in char p (becomes X+Y=1)")
    print(f"{'='*78}")

    # --- Step 1: Gauss sums ---
    gauss = compute_gauss_sums(p)
    print(f"\n  GAUSS SUMS g(chi^k) = sum chi^k(x) omega^x:")
    print(f"  {'k':>4s}  {'|g|':>10s}  {'arg/pi':>10s}  {'Re':>12s}  {'Im':>12s}  {'|g|=sqrt(p)?':>12s}")
    print(f"  {'-'*4:>4s}  {'-'*10:>10s}  {'-'*10:>10s}  {'-'*12:>12s}  {'-'*12:>12s}  {'-'*12:>12s}")
    for k in sorted(gauss):
        G = gauss[k]
        at_weil = "YES" if abs(abs(G) - sqrt_p) < 0.001 * sqrt_p else "no"
        print(f"  {k:>4d}  {abs(G):10.6f}  {cmath.phase(G)/PI:+10.6f}  "
              f"{G.real:+12.6f}  {G.imag:+12.6f}  {at_weil:>12s}")

    # --- Step 2: Jacobi sums ---
    jacobi = compute_jacobi_sums(p)
    n_weil_jacobi = sum(1 for J in jacobi.values() if abs(abs(J) - sqrt_p) < 0.01 * sqrt_p)
    print(f"\n  JACOBI SUMS J(chi^a, chi^b) = sum chi^a(x) chi^b(1-x):")
    print(f"  Valid pairs: {len(jacobi)}, at Weil bound: {n_weil_jacobi}")

    if jacobi:
        print(f"  {'(a,b)':>8s}  {'|J|':>10s}  {'arg/pi':>10s}  {'Re':>12s}  {'Im':>12s}  {'Weil?':>6s}")
        print(f"  {'-'*8:>8s}  {'-'*10:>10s}  {'-'*10:>10s}  {'-'*12:>12s}  {'-'*12:>12s}  {'-'*6:>6s}")
        for (a, b) in sorted(jacobi):
            J = jacobi[(a, b)]
            at_weil = "YES" if abs(abs(J) - sqrt_p) < 0.01 * sqrt_p else ""
            print(f"  ({a},{b})    {abs(J):10.6f}  {cmath.phase(J)/PI:+10.6f}  "
                  f"{J.real:+12.6f}  {J.imag:+12.6f}  {at_weil:>6s}")
    else:
        print(f"  (No valid Jacobi sum pairs for p={p} — "
              f"all (a,b) with 1<=a,b<={m-1} have a+b ≡ 0 mod {m})")

    # --- Step 3: Build candidate catalog ---
    catalog = build_candidate_catalog(gauss, jacobi, p)
    print(f"\n  Candidate catalog (|z| = sqrt({p})): {len(catalog)} entries")

    # --- Step 4: Compute Witt eigenvalues and match ---
    print(f"\n  Computing Witt carry eigenvalues (all twists)...")
    all_twist_results = compute_witt_eigenvalues_all_twists(p)

    all_matches = []

    for eigs, label in all_twist_results:
        weil = extract_weil(eigs, p)
        if not weil:
            continue
        print(f"\n  Twist: {label} — {len(weil)} Weil eigenvalues")

        for e in weil:
            best_name = None
            best_err = float('inf')
            best_val = None
            for name, z in catalog.items():
                err = abs(e - z) / sqrt_p
                if err < best_err:
                    best_err = err
                    best_name = name
                    best_val = z
            match_quality = ("EXACT" if best_err < 1e-6 else
                             "CLOSE" if best_err < 0.01 else
                             "ROUGH" if best_err < 0.05 else "NONE")
            print(f"    mu = {e.real:+10.6f}{e.imag:+10.6f}i  "
                  f"|mu|={abs(e):.6f}  arg/pi={cmath.phase(e)/PI:+.6f}  "
                  f"-> {best_name or 'N/A':>20s}  err={best_err:.2e}  {match_quality}")
            all_matches.append((label, e, best_name, best_val, best_err, match_quality))

    # --- Step 5: Summary ---
    print(f"\n  {'='*70}")
    print(f"  SUMMARY for p = {p}")
    print(f"  {'='*70}")

    exact_matches = [m for m in all_matches if m[5] in ("EXACT", "CLOSE")]
    rough_matches = [m for m in all_matches if m[5] == "ROUGH"]
    no_matches = [m for m in all_matches if m[5] == "NONE"]

    print(f"  Total Weil eigenvalues found: {len(all_matches)}")
    print(f"  Exact/close matches to Gauss/Jacobi: {len(exact_matches)}")
    print(f"  Rough matches (1-5%): {len(rough_matches)}")
    print(f"  No match: {len(no_matches)}")

    if exact_matches:
        matched_names = set(m[2] for m in exact_matches)
        print(f"  Matched sums: {sorted(matched_names)}")

    # Check Galois structure: are matched Gauss sums related by conjugation?
    matched_gauss = [m for m in exact_matches if m[2] and m[2].startswith("g(") or
                     (m[2] and m[2].startswith("-g("))]
    matched_jacobi = [m for m in exact_matches if m[2] and m[2].startswith("J(")]
    if matched_gauss:
        print(f"  Matched Gauss sums: {len(matched_gauss)}")
    if matched_jacobi:
        print(f"  Matched Jacobi sums: {len(matched_jacobi)}")


# =====================================================================
# FINAL SUMMARY
# =====================================================================
print(f"\n{'='*78}")
print(f"  F12 FINAL SUMMARY")
print(f"{'='*78}")
print("""
  This experiment computes all Gauss sums g(chi^k) and Jacobi sums
  J(chi^a, chi^b) for multiplicative characters on F_p*, then matches
  them against the Weil-bound eigenvalues of Witt carry matrices.

  THEORETICAL NOTE:
  The Fermat curve X^p + Y^p = 1 degenerates in characteristic p
  (it becomes X + Y = 1, a line). The standard Jacobi sum formula
  for Fermat curve Frobenius does NOT apply directly.

  The actual variety whose cohomology the carry operator computes
  is the Artin-Schreier variety z^p - z = c_2(x0,x1,y0,y1),
  which is explored in F11.

  Nevertheless, Gauss sums g(chi^k) with |g| = sqrt(p) provide
  the complete set of "available" complex numbers at the Weil bound,
  and any Frobenius eigenvalue of a variety over F_p must be
  expressible in terms of these.
""")
