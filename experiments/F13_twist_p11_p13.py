#!/usr/bin/env python3
"""
F13: Extended Twist Scan at p=11 and p=13

Extends the F9 character twist methodology to p=11 (121x121 matrix)
and p=13 (169x169 matrix). For each prime, tests ALL multiplicative
and additive character twists to find Weil-bound eigenvalues.

Goal: collect the trace sequence (a_3, a_5, a_7, a_11, a_13) for
potential LMFDB matching.

Known results from F8c/F9:
  p=3:  a_3 = -3 (untwisted, Theorem A)
  p=5:  a_5 from additive chi_2 twist (2 eigenvalues at sqrt(5))
  p=7:  a_7 = 0 (untwisted, 6 eigenvalues at sqrt(7))
  p=11: ? (untwisted: NO Weil eigenvalues — can twisting rescue?)
  p=13: ? (untwisted: NO Weil eigenvalues — can twisting rescue?)
"""

import math
import cmath
import time
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


def extract_weil(eigs, p, tol_frac=0.005):
    sqrt_p = math.sqrt(p)
    return [e for e in eigs if abs(abs(e) - sqrt_p) < tol_frac * sqrt_p and abs(e) > 1e-6]


# =====================================================================
# MAIN
# =====================================================================

print("=" * 78)
print("  F13: EXTENDED TWIST SCAN — p=11 AND p=13")
print("=" * 78)

results_by_prime = {}

for p in [3, 5, 7, 11, 13]:
    print(f"\n{'='*78}")
    print(f"  PRIME p = {p}  (matrix: {p*p}x{p*p}, {p-1} characters)")
    print(f"{'='*78}")

    t0 = time.time()
    E, states, state_idx = build_level2_exponent_matrix(p)
    t_build = time.time() - t0
    print(f"  Exponent matrix built in {t_build:.1f}s")

    chars, g, dlog = build_character_table(p)
    omega = np.exp(2j * PI / p)
    D = p * p
    sqrt_p = math.sqrt(p)

    prime_results = []

    # --- Untwisted ---
    M_base = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            M_base[i, j] = omega ** E[i, j]
    eigs_base = np.linalg.eigvals(M_base)
    weil_base = extract_weil(eigs_base, p)
    n_weil = len(weil_base)

    if n_weil > 0:
        traces = [e for e in weil_base]
        a_p = -sum(e for e in weil_base).real
        print(f"  Untwisted: {n_weil} Weil eigenvalues, a_p = {a_p:+.4f}")
        for e in weil_base:
            print(f"    mu = {e.real:+10.6f}{e.imag:+10.6f}i  |mu|={abs(e):.6f}  arg/pi={cmath.phase(e)/PI:+.6f}")
        prime_results.append(("untwisted", n_weil, weil_base))
    else:
        print(f"  Untwisted: 0 Weil eigenvalues")
        prime_results.append(("untwisted", 0, []))

    # --- Multiplicative twists ---
    print(f"\n  --- Multiplicative twists chi_k(x0*y0) ---")
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
        eigs = np.linalg.eigvals(M)
        weil = extract_weil(eigs, p)
        label = f"mul chi_{k} (ord {order})"
        if weil:
            print(f"  {label}: {len(weil)} Weil eigenvalues")
            for e in weil:
                print(f"    mu = {e.real:+10.6f}{e.imag:+10.6f}i  |mu|={abs(e):.6f}  arg/pi={cmath.phase(e)/PI:+.6f}")
        prime_results.append((label, len(weil), weil))

    # --- Additive twists ---
    print(f"\n  --- Additive twists chi_k(x0+y0) ---")
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
        eigs = np.linalg.eigvals(M)
        weil = extract_weil(eigs, p)
        label = f"add chi_{k} (ord {order})"
        if weil:
            print(f"  {label}: {len(weil)} Weil eigenvalues")
            for e in weil:
                print(f"    mu = {e.real:+10.6f}{e.imag:+10.6f}i  |mu|={abs(e):.6f}  arg/pi={cmath.phase(e)/PI:+.6f}")
        prime_results.append((label, len(weil), weil))

    # --- Summary table for this prime ---
    print(f"\n  Summary for p={p}:")
    total_weil = 0
    productive_twists = []
    for (label, n_weil, weil_eigs) in prime_results:
        if n_weil > 0:
            total_weil += n_weil
            productive_twists.append((label, n_weil, weil_eigs))

    if productive_twists:
        print(f"  Total Weil eigenvalues found: {total_weil} across {len(productive_twists)} twists")
        for (label, n, w) in productive_twists:
            a_trace = -sum(e for e in w)
            print(f"    {label:>30s}: {n:2d} eigenvalues, trace = {a_trace.real:+.4f}{a_trace.imag:+.4f}i")
    else:
        print(f"  NO Weil eigenvalues found from ANY twist!")

    t_total = time.time() - t0
    print(f"  Total time for p={p}: {t_total:.1f}s")

    results_by_prime[p] = prime_results


# =====================================================================
# GLOBAL SUMMARY AND LMFDB MATCHING
# =====================================================================

print(f"\n{'='*78}")
print(f"  F13: GLOBAL TRACE SEQUENCE")
print(f"{'='*78}")
print(f"\n  Assembling (a_3, a_5, a_7, a_11, a_13) from best-matching twists:")

trace_sequence = {}
for p in [3, 5, 7, 11, 13]:
    results = results_by_prime[p]
    productive = [(l, n, w) for (l, n, w) in results if n > 0]

    if not productive:
        print(f"\n  p={p}: NO productive twists — a_{p} = UNKNOWN")
        trace_sequence[p] = None
        continue

    # For the "best" trace, use the untwisted result if available,
    # otherwise the twist with the most Weil eigenvalues
    best = None
    for (l, n, w) in productive:
        if l == "untwisted":
            best = (l, n, w)
            break
    if best is None:
        best = max(productive, key=lambda x: x[1])

    label, n_weil, weil_eigs = best
    a_p = -sum(e for e in weil_eigs)

    # Check if a_p is a real integer
    if abs(a_p.imag) < 0.01 and abs(a_p.real - round(a_p.real)) < 0.01:
        a_p_int = int(round(a_p.real))
        print(f"\n  p={p} ({label}): a_{p} = {a_p_int}")
        trace_sequence[p] = a_p_int
    else:
        print(f"\n  p={p} ({label}): a_{p} = {a_p.real:+.4f}{a_p.imag:+.4f}i (non-integer!)")
        trace_sequence[p] = a_p

print(f"\n  Trace sequence:")
for p in [3, 5, 7, 11, 13]:
    a = trace_sequence.get(p)
    if a is None:
        print(f"    a_{p:2d} = UNKNOWN")
    elif isinstance(a, int):
        print(f"    a_{p:2d} = {a:+d}")
    else:
        print(f"    a_{p:2d} = {a.real:+.4f}{a.imag:+.4f}i")

# --- LMFDB: brute-force search for matching elliptic curves ---
print(f"\n{'='*78}")
print(f"  LMFDB-STYLE SEARCH: Elliptic curves matching trace sequence")
print(f"{'='*78}")

def count_curve_points(A, B, p):
    """Count #E(F_p) for E: y^2 = x^3 + Ax + B, return a_p = p + 1 - #E."""
    if (4 * A**3 + 27 * B**2) % p == 0:
        return None  # singular
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x**3 + A * x + B) % p
        if rhs == 0:
            count += 1
        else:
            # Legendre symbol
            leg = pow(rhs, (p - 1) // 2, p)
            if leg == 1:
                count += 2
    return p + 1 - count

known_traces = {p: t for p, t in trace_sequence.items() if isinstance(t, int)}
target_primes = sorted(known_traces.keys())

if len(target_primes) >= 2:
    print(f"\n  Searching y^2 = x^3 + Ax + B with |A|,|B| <= 20...")
    print(f"  Must match: {', '.join(f'a_{p}={known_traces[p]:+d}' for p in target_primes)}")

    candidates = []
    for A in range(-20, 21):
        for B in range(-20, 21):
            match = True
            trace_values = {}
            for p in target_primes:
                a = count_curve_points(A, B, p)
                if a is None:
                    trace_values[p] = None
                    if p in known_traces:
                        match = False
                        break
                else:
                    trace_values[p] = a
                    if p in known_traces and a != known_traces[p]:
                        match = False
                        break
            if match:
                # Compute a_p for all test primes
                full_traces = {}
                for p_test in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
                    full_traces[p_test] = count_curve_points(A, B, p_test)
                candidates.append((A, B, full_traces))

    if candidates:
        print(f"\n  Found {len(candidates)} candidate curves:")
        header = "  y^2 = x^3 + Ax + B"
        for p_test in [2, 3, 5, 7, 11, 13, 17, 19]:
            header += f"  a_{p_test:2d}"
        print(header)
        print("  " + "-" * len(header))
        for (A, B, ft) in candidates[:20]:
            line = f"  A={A:+3d}, B={B:+3d}         "
            for p_test in [2, 3, 5, 7, 11, 13, 17, 19]:
                v = ft.get(p_test)
                if v is None:
                    line += f"  bad"
                else:
                    line += f"  {v:+3d}"
            print(line)
    else:
        print(f"\n  No matching curves found in range |A|,|B| <= 20")
else:
    print(f"\n  Not enough integer traces to search. Need at least 2.")


# =====================================================================
# FINAL SUMMARY
# =====================================================================
print(f"\n{'='*78}")
print(f"  F13 FINAL SUMMARY")
print(f"{'='*78}")
print(f"""
  Extended twist scan for the Witt carry operator.

  For each prime p, tested:
  - {'{'}p-1{'}'} multiplicative twists chi_k(x0*y0 mod p)
  - {'{'}p-1{'}'} additive twists chi_k(x0+y0 mod p)

  Twist scan results:
""")
for p in [3, 5, 7, 11, 13]:
    results = results_by_prime[p]
    n_productive = sum(1 for (_, n, _) in results if n > 0)
    n_total_weil = sum(n for (_, n, _) in results)
    productive_labels = [l for (l, n, _) in results if n > 0]
    a = trace_sequence.get(p)
    a_str = f"{a:+d}" if isinstance(a, int) else ("UNKNOWN" if a is None else f"{a}")
    print(f"  p={p:2d}: {n_productive:2d} productive twists, "
          f"{n_total_weil:3d} total Weil eigenvalues, a_{p} = {a_str}")
    if productive_labels:
        for l in productive_labels:
            print(f"         - {l}")
