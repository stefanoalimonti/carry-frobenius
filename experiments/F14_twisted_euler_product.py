#!/usr/bin/env python3
"""
F14: Twisted Euler Product Coherence

Tests whether the Legendre-twisted Witt carry Frobenius traces assemble
into a coherent L-function across primes.

F9 showed that character twists recover Weil-bound eigenvalues at p=5,13.
F10 showed that the UNTWISTED carry traces are cross-prime INCONSISTENT.
Question: are the TWISTED traces consistent?

Method:
  For each prime p = 3, 5, 7, 11, 13, 17, 19, 23:
    - Build the p^2 x p^2 Witt carry matrix (level 2, mul)
    - Twist by additive Legendre character chi_2(x0 + y0 mod p)
    - Extract Weil-bound eigenvalues |mu| = sqrt(p)
    - Compute Frobenius trace a_p = -sum(Weil eigenvalues)
  Assemble (a_3, a_5, ..., a_23) and search for pattern / LMFDB match.
  Compute partial Euler product.
"""

import math
import numpy as np
from itertools import product as cartesian_product

PI = math.pi


def witt_mul_exact(x, y, p, n):
    """Witt multiplication via ghost components."""
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
    """Build the p^2 x p^2 exponent matrix."""
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
                    pr = witt_mul_exact(x, y, p, 3)
                    g2x = x0**(p**2) + p * x1**p
                    g2y = y0**(p**2) + p * y1**p
                    partial = g2x * g2y - pr[0]**(p**2) - p * pr[1]**p
                    c2 = (partial // (p**2)) % p
                    i_from = state_idx[(x0, y0)]
                    j_to = state_idx[(x1, y1)]
                    E[j_to, i_from] = c2
    return E, states, state_idx


def primitive_root(p):
    for g in range(2, p):
        order, val = 1, g % p
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
        chi = {0: 0}
        for a in range(1, p):
            chi[a] = zeta ** (k * dlog[a])
        chars[k] = chi
    return chars, g, dlog


def build_additive_twisted_matrix(E, p, chi, states, state_idx):
    """M(chi)[j,i] = chi(x0 + y0 mod p) * omega^{c_2}"""
    omega = np.exp(2j * PI / p)
    D = p * p
    M = np.zeros((D, D), dtype=complex)
    for x0 in range(p):
        for y0 in range(p):
            sum_mod_p = (x0 + y0) % p
            chi_val = chi[sum_mod_p]
            i_from = state_idx[(x0, y0)]
            for x1 in range(p):
                for y1 in range(p):
                    j_to = state_idx[(x1, y1)]
                    M[j_to, i_from] = chi_val * omega ** E[j_to, i_from]
    return M


def build_multiplicative_twisted_matrix(E, p, chi, states, state_idx):
    """M(chi)[j,i] = chi(x0 * y0 mod p) * omega^{c_2}"""
    omega = np.exp(2j * PI / p)
    D = p * p
    M = np.zeros((D, D), dtype=complex)
    for x0 in range(p):
        for y0 in range(p):
            prod_mod_p = (x0 * y0) % p
            chi_val = chi[prod_mod_p]
            i_from = state_idx[(x0, y0)]
            for x1 in range(p):
                for y1 in range(p):
                    j_to = state_idx[(x1, y1)]
                    M[j_to, i_from] = chi_val * omega ** E[j_to, i_from]
    return M


def extract_weil_eigenvalues(M, p, tol_frac=0.05):
    """Extract eigenvalues at the Weil bound |mu| = sqrt(p)."""
    eigs = np.linalg.eigvals(M)
    target = math.sqrt(p)
    weil = [e for e in eigs if abs(abs(e) - target) < tol_frac * target and abs(e) > 0.01]
    return weil


# =====================================================================
print("=" * 78)
print("  F14: TWISTED EULER PRODUCT COHERENCE")
print("=" * 78)

primes = [3, 5, 7, 11, 13, 17, 19, 23]
results = {}

for p in primes:
    print(f"\n{'='*78}")
    print(f"  Prime p = {p} (matrix size {p**2}x{p**2})")
    print(f"{'='*78}")

    E, states, state_idx = build_level2_exponent_matrix(p)
    chars, g, dlog = build_character_table(p)

    # Find the Legendre character (quadratic, order 2)
    # It's chi_{(p-1)/2}
    legendre_idx = (p - 1) // 2
    chi_leg = chars[legendre_idx]

    print(f"  Primitive root: {g}")
    print(f"  Legendre character index: {legendre_idx}")

    # Test untwisted first
    omega = np.exp(2j * PI / p)
    D = p * p
    M_plain = np.zeros((D, D), dtype=complex)
    for i_col in range(D):
        for j_row in range(D):
            M_plain[j_row, i_col] = omega ** E[j_row, i_col]

    weil_plain = extract_weil_eigenvalues(M_plain, p)
    print(f"\n  Untwisted: {len(weil_plain)} Weil eigenvalues at |mu|=sqrt({p})={math.sqrt(p):.4f}")
    if weil_plain:
        trace_plain = -sum(e for e in weil_plain).real
        print(f"    a_p(untwisted) = {trace_plain:.4f}")
        for e in sorted(weil_plain, key=lambda x: -abs(x)):
            print(f"      mu = {e.real:+.6f} {e.imag:+.6f}i  |mu| = {abs(e):.6f}")

    # Additive Legendre twist
    M_add_leg = build_additive_twisted_matrix(E, p, chi_leg, states, state_idx)
    weil_add = extract_weil_eigenvalues(M_add_leg, p)
    print(f"\n  Additive Legendre twist: {len(weil_add)} Weil eigenvalues")
    trace_add = 0.0
    if weil_add:
        trace_add = -sum(e for e in weil_add).real
        print(f"    a_p(add_Leg) = {trace_add:.4f}")
        for e in sorted(weil_add, key=lambda x: -abs(x)):
            print(f"      mu = {e.real:+.6f} {e.imag:+.6f}i  |mu| = {abs(e):.6f}")

    # Multiplicative Legendre twist
    M_mul_leg = build_multiplicative_twisted_matrix(E, p, chi_leg, states, state_idx)
    weil_mul = extract_weil_eigenvalues(M_mul_leg, p)
    print(f"\n  Multiplicative Legendre twist: {len(weil_mul)} Weil eigenvalues")
    trace_mul = 0.0
    if weil_mul:
        trace_mul = -sum(e for e in weil_mul).real
        print(f"    a_p(mul_Leg) = {trace_mul:.4f}")
        for e in sorted(weil_mul, key=lambda x: -abs(x)):
            print(f"      mu = {e.real:+.6f} {e.imag:+.6f}i  |mu| = {abs(e):.6f}")

    # Scan ALL character twists (additive)
    best_add = []
    for k in range(p - 1):
        chi_k = chars[k]
        M_tw = build_additive_twisted_matrix(E, p, chi_k, states, state_idx)
        weil_tw = extract_weil_eigenvalues(M_tw, p)
        if weil_tw:
            tr = -sum(e for e in weil_tw).real
            best_add.append((k, len(weil_tw), tr))

    if best_add:
        print(f"\n  Productive additive twists:")
        for k, n_weil, tr in best_add:
            order = (p - 1) // math.gcd(k, p - 1) if k > 0 else 1
            print(f"    chi_{k} (order {order}): {n_weil} Weil eigs, a_p = {tr:.4f}")

    results[p] = {
        'n_weil_plain': len(weil_plain),
        'n_weil_add_leg': len(weil_add),
        'n_weil_mul_leg': len(weil_mul),
        'trace_plain': -sum(e for e in weil_plain).real if weil_plain else None,
        'trace_add_leg': trace_add if weil_add else None,
        'trace_mul_leg': trace_mul if weil_mul else None,
        'best_add': best_add,
    }

# =====================================================================
# SUMMARY TABLE
# =====================================================================
print("\n" + "=" * 78)
print("  SUMMARY: Cross-Prime Trace Table")
print("=" * 78)

print(f"\n  {'p':>3s}  {'#Weil(plain)':>12s}  {'a_p(plain)':>12s}  "
      f"{'#Weil(add_Leg)':>14s}  {'a_p(add_Leg)':>14s}")
for p in primes:
    r = results[p]
    ap = f"{r['trace_plain']:.2f}" if r['trace_plain'] is not None else "—"
    al = f"{r['trace_add_leg']:.2f}" if r['trace_add_leg'] is not None else "—"
    print(f"  {p:3d}  {r['n_weil_plain']:12d}  {ap:>12s}  "
          f"{r['n_weil_add_leg']:14d}  {al:>14s}")

# =====================================================================
# EULER PRODUCT
# =====================================================================
print("\n" + "=" * 78)
print("  Partial Euler Product")
print("=" * 78)

# Assemble L(s) = prod_p (1 - a_p * p^{-s} + p^{1-2s})^{-1}
# for s values on Re(s) = 1

print("\n  Using additive Legendre traces (where available):")
s_values = [1.0, 1.5, 2.0]
for s in s_values:
    L = 1.0
    for p in primes:
        r = results[p]
        ap = r['trace_add_leg']
        if ap is None:
            ap = r['trace_plain']
        if ap is None:
            ap = 0.0  # no Weil eigenvalues = trivial factor
        factor = 1.0 / (1.0 - ap * p**(-s) + p**(1 - 2*s))
        L *= factor
    print(f"  s = {s:.1f}: L(s) = {L:.6f}")

print("\n" + "=" * 78)
print("  F14 COMPLETE")
print("=" * 78)
