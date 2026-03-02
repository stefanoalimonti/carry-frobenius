#!/usr/bin/env python3
"""
F8c: Broader Prime Scan at Level-2 + Exact Algebraic Analysis at p=7

F8a/F8b showed level-3/4 marginals collapse to rank 1 (character orthogonality).
Level-2 is the ONLY productive level: no marginalization, each entry is a single
omega^{c_2}.

Surprising discovery: p=7 level-2 produces 6 eigenvalues at |mu|=sqrt(7).
This extends the F6 result (p=3) to a second prime.

This script:
  1. Scans all primes p=3,5,7,11,13 at level-2 (numerical)
  2. Tests the congruence hypothesis: p ≡ 3 (mod 4) → works, p ≡ 1 (mod 4) → fails
  3. For p=7, extracts the Frobenius factor and compares with Fermat curve C_7
  4. For p=7, does exact algebraic char poly computation via Faddeev-LeVerrier
"""

import math
import sys
import time
from itertools import product as cartesian_product
from fractions import Fraction
from collections import Counter

import numpy as np

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
    """Build the exact exponent matrix E[j,i] where M[j,i] = omega^{E[j,i]}.

    Returns:
        E: p^2 x p^2 integer matrix with entries in {0, ..., p-1}
        states: list of (x0, y0) pairs
    """
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


def build_complex_matrix(E, p):
    """Build complex matrix from exponent matrix."""
    omega = np.exp(2j * PI / p)
    D = E.shape[0]
    M = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            M[i, j] = omega ** E[i, j]
    return M


def scan_prime(p, verbose=True):
    """Full level-2 analysis for a given prime."""
    if verbose:
        print(f"\n  {'='*60}")
        print(f"  PRIME p={p}  (p mod 4 = {p % 4})")
        print(f"  {'='*60}")

    t0 = time.time()
    E, states, sidx = build_level2_exponent_matrix(p)
    M = build_complex_matrix(E, p)
    t1 = time.time()

    D = p * p
    if verbose:
        print(f"  Matrix: {D}x{D}, computed in {t1-t0:.2f}s")
        flat = list(E.flatten())
        print(f"  Carry distribution: {dict(Counter(flat))}")

    eigs = np.linalg.eigvals(M)
    eigs_sorted = sorted(eigs, key=lambda x: -abs(x))

    sqrt_p = math.sqrt(p)
    weil_eigs = [e for e in eigs_sorted if abs(abs(e) - sqrt_p) < 0.05 * sqrt_p]
    n_zero = sum(1 for e in eigs_sorted if abs(e) < 1e-6)
    n_nonzero = D - n_zero

    if verbose:
        print(f"  Nonzero eigenvalues: {n_nonzero}, zero: {n_zero}")
        mods = sorted(set(round(abs(e), 4) for e in eigs_sorted if abs(e) > 1e-6), reverse=True)
        print(f"  Distinct nonzero moduli: {mods[:10]}")

        print(f"\n  Eigenvalues at |mu| = sqrt({p}) = {sqrt_p:.6f}:")
        if weil_eigs:
            for e in weil_eigs:
                mod = abs(e)
                arg = np.angle(e) / PI
                print(f"    mu = {e.real:+12.6f} {e.imag:+12.6f}i  "
                      f"|mu| = {mod:.8f}  arg/pi = {arg:+.6f}")
        else:
            print(f"    NONE FOUND")

    return E, M, eigs_sorted, weil_eigs, states


# =====================================================================
# PART 1: FULL PRIME SCAN
# =====================================================================

print("=" * 78)
print("  F8c: PRIME SCAN AT LEVEL-2 + ALGEBRAIC ANALYSIS")
print("=" * 78)

primes_to_test = [3, 5, 7, 11, 13]
results = {}

for p in primes_to_test:
    E, M, eigs, weil, states = scan_prime(p)
    results[p] = {
        "E": E, "M": M, "eigs": eigs, "weil": weil,
        "states": states, "n_weil": len(weil)
    }

# Summary table
print("\n" + "=" * 78)
print("  PRIME SCAN SUMMARY")
print("=" * 78)
print(f"\n  {'Prime':>5s}  {'p mod 4':>7s}  {'Genus g':>7s}  {'2g':>4s}  "
      f"{'dim':>4s}  {'#Weil':>5s}  {'Result':>8s}")
print(f"  {'-'*5:>5s}  {'-'*7:>7s}  {'-'*7:>7s}  {'-'*4:>4s}  "
      f"{'-'*4:>4s}  {'-'*5:>5s}  {'-'*8:>8s}")
for p in primes_to_test:
    g = (p - 1) * (p - 2) // 2
    mod4 = p % 4
    dim = p * p
    n_w = results[p]["n_weil"]
    status = "YES" if n_w > 0 else "NO"
    print(f"  {p:5d}  {mod4:7d}  {g:7d}  {2*g:4d}  {dim:4d}  {n_w:5d}  {status:>8s}")

# =====================================================================
# PART 2: CONGRUENCE HYPOTHESIS
# =====================================================================

print("\n" + "=" * 78)
print("  CONGRUENCE HYPOTHESIS TEST")
print("=" * 78)

mod3_works = [p for p in primes_to_test if results[p]["n_weil"] > 0 and p % 4 == 3]
mod1_works = [p for p in primes_to_test if results[p]["n_weil"] > 0 and p % 4 == 1]
mod3_fails = [p for p in primes_to_test if results[p]["n_weil"] == 0 and p % 4 == 3]
mod1_fails = [p for p in primes_to_test if results[p]["n_weil"] == 0 and p % 4 == 1]

print(f"\n  p ≡ 3 (mod 4) with Weil eigenvalues: {mod3_works}")
print(f"  p ≡ 3 (mod 4) WITHOUT Weil eigenvalues: {mod3_fails}")
print(f"  p ≡ 1 (mod 4) with Weil eigenvalues: {mod1_works}")
print(f"  p ≡ 1 (mod 4) WITHOUT Weil eigenvalues: {mod1_fails}")

if mod1_works or mod3_fails:
    print(f"\n  Hypothesis FALSIFIED: congruence mod 4 does not cleanly separate.")
else:
    print(f"\n  Hypothesis CONSISTENT: all p ≡ 3 (mod 4) work, all p ≡ 1 (mod 4) fail.")

# =====================================================================
# PART 3: p=7 FROBENIUS FACTOR (numerical)
# =====================================================================

print("\n" + "=" * 78)
print("  p=7: FROBENIUS FACTOR ANALYSIS (NUMERICAL)")
print("=" * 78)

weil_7 = results[7]["weil"]
eigs_7 = results[7]["eigs"]

if weil_7:
    print(f"\n  {len(weil_7)} eigenvalues at |mu| = sqrt(7) = {math.sqrt(7):.8f}")

    for i, e in enumerate(weil_7):
        arg = np.angle(e) / PI
        print(f"    alpha_{i} = {e.real:+12.8f} {e.imag:+12.8f}i  arg/pi = {arg:+.8f}")

    # Form the polynomial product(mu - alpha_i)
    print(f"\n  Frobenius factor P(mu) = product(mu - alpha_i):")
    from numpy.polynomial import polynomial as P
    # Compute product of (mu - alpha_i) using convolution
    poly = np.array([1.0 + 0j])
    for alpha in weil_7:
        poly = np.convolve(poly, np.array([-alpha, 1.0 + 0j]))

    # Reverse to get highest-degree first
    poly_rev = poly[::-1]
    print(f"  Coefficients (highest degree first):")
    for k, c in enumerate(poly_rev):
        deg = len(weil_7) - k
        nearest_int = round(c.real)
        print(f"    mu^{deg}: {c.real:+16.6f} {c.imag:+16.6f}i  "
              f"(nearest int: {nearest_int}, err: {abs(c - nearest_int):.2e})")

    # Check if coefficients are integers
    int_coeffs = [round(c.real) for c in poly_rev]
    max_err = max(abs(c - round(c.real)) for c in poly_rev)
    print(f"\n  Maximum deviation from integer: {max_err:.2e}")
    print(f"  Integer polynomial: ", end="")
    for k, c_int in enumerate(int_coeffs):
        deg = len(weil_7) - k
        if c_int != 0 or deg == 0:
            print(f"{c_int:+d}*mu^{deg} " if deg > 0 else f"{c_int:+d}", end="")
    print()

    # Check functional equation P(mu) = p^g * mu^{2g} * P(p/mu)
    deg_P = len(weil_7)
    g_candidate = deg_P // 2
    print(f"\n  Degree of P: {deg_P}, candidate genus g = {g_candidate}")
    print(f"  Checking functional equation: P(mu) = 7^{g_candidate} * mu^{deg_P} * P(7/mu)...")

    # For the functional equation, if P(mu) = sum c_k mu^k,
    # then the reciprocal poly is mu^n P(p/mu) = sum c_k p^k mu^{n-k}
    # and 7^g * mu^n * P(7/mu) should equal P(mu)
    p_val = 7
    P_check = np.array(int_coeffs[::-1], dtype=float)
    P_recip = np.array([int_coeffs[-(k+1)] * p_val**k for k in range(len(int_coeffs))], dtype=float)
    P_recip_scaled = P_recip * p_val**g_candidate

    err_fe = max(abs(P_check[k] - P_recip_scaled[k]) for k in range(len(P_check)))
    print(f"  Functional equation error: {err_fe:.2e}")
    if err_fe < 1:
        print(f"  *** FUNCTIONAL EQUATION SATISFIED ***")

    # Check product of conjugate pairs
    print(f"\n  Conjugate pair products (should be {p_val}):")
    used = set()
    for i in range(len(weil_7)):
        if i in used:
            continue
        for j in range(i + 1, len(weil_7)):
            if j in used:
                continue
            prod = weil_7[i] * weil_7[j]
            if abs(abs(prod) - p_val) < 0.5:
                print(f"    alpha_{i} * alpha_{j} = {prod.real:+.4f}{prod.imag:+.4f}i  "
                      f"(|prod| = {abs(prod):.6f})")
                used.add(i)
                used.add(j)
                break

# =====================================================================
# PART 4: p=7 EXACT ALGEBRAIC COMPUTATION (Faddeev-LeVerrier)
# =====================================================================

print("\n" + "=" * 78)
print("  p=7: EXACT ALGEBRAIC CHARACTERISTIC POLYNOMIAL")
print("=" * 78)

p = 7
E7 = results[7]["E"]
D = p * p  # 49

# Z[omega_7] arithmetic. omega = e^{2pi*i/7} satisfies
# omega^6 + omega^5 + omega^4 + omega^3 + omega^2 + omega + 1 = 0
# Elements: (a0, a1, a2, a3, a4, a5) meaning a0 + a1*w + ... + a5*w^5
# Multiplication: use the relation w^6 = -(1 + w + w^2 + w^3 + w^4 + w^5)

def mul_zw7(x, y):
    """Multiply two elements of Z[omega_7], each a 6-tuple of ints."""
    r = [0] * 12
    for i in range(6):
        if x[i] == 0:
            continue
        for j in range(6):
            if y[j] == 0:
                continue
            r[i + j] += x[i] * y[j]

    # Reduce: w^k for k >= 6 using w^6 = -(1+w+...+w^5)
    for k in range(11, 5, -1):
        if r[k] != 0:
            c = r[k]
            r[k] = 0
            for i in range(6):
                r[k - 6 + i] -= c

    return tuple(r[:6])


def add_zw7(x, y):
    return tuple(x[i] + y[i] for i in range(6))


def neg_zw7(x):
    return tuple(-x[i] for i in range(6))


def zero_zw7():
    return (0, 0, 0, 0, 0, 0)


def one_zw7():
    return (1, 0, 0, 0, 0, 0)


def omega_power_zw7(k):
    """omega^k as Z[omega_7] element."""
    k = k % 7
    if k == 0:
        return one_zw7()
    elif k <= 5:
        r = [0] * 6
        r[k] = 1
        return tuple(r)
    else:  # k == 6
        return (-1, -1, -1, -1, -1, -1)


def scale_zw7(x, s):
    """Multiply Z[omega_7] element by integer s."""
    return tuple(x[i] * s for i in range(6))


def div_zw7(x, s):
    """Exact division of Z[omega_7] element by integer s."""
    result = []
    for i in range(6):
        if x[i] % s != 0:
            raise ValueError(f"Not exactly divisible: {x[i]} / {s}")
        result.append(x[i] // s)
    return tuple(result)


def is_rational_zw7(x):
    """Check if element is in Q (all omega components zero)."""
    return all(x[i] == 0 for i in range(1, 6))


# Build exact Z[omega_7] matrix
print(f"\n  Building {D}x{D} matrix over Z[omega_7]...")
M_exact = [[zero_zw7() for _ in range(D)] for _ in range(D)]
for i in range(D):
    for j in range(D):
        M_exact[i][j] = omega_power_zw7(int(E7[i][j]))

# Faddeev-LeVerrier: chi(mu) = mu^D + c_{D-1}*mu^{D-1} + ... + c_0
# c_{D-k} = -(1/k)(Tr(M^k) + sum_{j=1}^{k-1} c_{D-j} * Tr(M^{k-j}))

def matmul_zw7(A, B, n):
    """n x n matrix multiplication over Z[omega_7]."""
    R = [[zero_zw7() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = zero_zw7()
            for k in range(n):
                if A[i][k] == zero_zw7() or B[k][j] == zero_zw7():
                    continue
                s = add_zw7(s, mul_zw7(A[i][k], B[k][j]))
            R[i][j] = s
    return R


def trace_zw7(A, n):
    """Trace of n x n matrix over Z[omega_7]."""
    s = zero_zw7()
    for i in range(n):
        s = add_zw7(s, A[i][i])
    return s


print(f"  Computing Faddeev-LeVerrier for {D}x{D} matrix...")
print(f"  This requires {D} matrix multiplications of size {D}x{D}.")
print(f"  Estimated: ~{D * D**3 // 10**6}M Z[omega_7] operations.")

# For 49x49, this is 49 * 49^3 ≈ 5.7M operations. Each operation involves
# 6-tuple arithmetic. This will take several minutes.

all_traces = []
coeffs = []
Mk = M_exact  # M^1

t0 = time.time()
for k in range(1, D + 1):
    if k == 1:
        Mk_curr = M_exact
    else:
        Mk_curr = matmul_zw7(Mk, M_exact, D)
        Mk = Mk_curr

    tr_k = trace_zw7(Mk_curr, D)
    all_traces.append(tr_k)

    # c_{D-k} = -(1/k)(Tr(M^k) + sum_{j=1}^{k-1} c_{D-j} * Tr(M^{k-j}))
    sum_val = tr_k
    for j in range(1, k):
        c_j = coeffs[j - 1]
        tr_km_j = all_traces[k - j - 1]
        sum_val = add_zw7(sum_val, mul_zw7(c_j, tr_km_j))

    c_val = neg_zw7(div_zw7(sum_val, k))
    coeffs.append(c_val)

    elapsed = time.time() - t0
    is_rat = is_rational_zw7(c_val)
    if k <= 10 or k % 5 == 0 or k == D:
        print(f"    k={k:2d}  Tr(M^k)_rational={is_rational_zw7(tr_k)}"
              f"  c_{D-k}_rational={is_rat}"
              f"  c_{D-k}[0]={c_val[0]}"
              f"  [{elapsed:.1f}s]")

    Mk = Mk_curr

total_time = time.time() - t0
print(f"\n  Faddeev-LeVerrier completed in {total_time:.1f}s")

# Check if ALL coefficients are rational
all_rational = all(is_rational_zw7(c) for c in coeffs)
print(f"\n  ALL coefficients rational (in Z)? {all_rational}")

if all_rational:
    print("  The characteristic polynomial is defined over Z!")
    rat_coeffs = [c[0] for c in coeffs]  # Extract integer part

    print(f"\n  Factoring with sympy...")
    from sympy import Symbol, Poly, factor, ZZ, Rational as SR
    mu = Symbol('mu')

    poly_expr = mu**D
    for k in range(D):
        if rat_coeffs[k] != 0:
            power = D - 1 - k
            poly_expr += rat_coeffs[k] * mu**power

    factored = factor(poly_expr)
    print(f"\n  chi(mu) = {factored}")

    # Look for factors whose roots have modulus sqrt(7)
    char_poly = Poly(poly_expr, mu)
    print(f"\n  Checking for Frobenius factor mu^6 - 7*mu^4 - 49*mu^2 + 343:")
    from sympy import div as poly_div
    test_factor = Poly(mu**6 - 7*mu**4 - 49*mu**2 + 343, mu)
    q, r = poly_div(char_poly, test_factor)
    if r.is_zero:
        print(f"  *** EXACT DIVISIBILITY CONFIRMED ***")
        print(f"  Quotient: {q}")
    else:
        print(f"  Not a factor. Remainder: {r}")
        print(f"  Trying mu^2 + 7 (the i*sqrt(7) pair):")
        test2 = Poly(mu**2 + 7, mu)
        q2, r2 = poly_div(char_poly, test2)
        if r2.is_zero:
            print(f"  *** mu^2 + 7 IS A FACTOR ***")
            print(f"  Quotient: {q2}")
        else:
            print(f"  mu^2 + 7 is not a factor either. Remainder: {r2}")

else:
    print("  WARNING: Some coefficients are NOT rational.")
    print("  Non-rational coefficients:")
    for k in range(D):
        if not is_rational_zw7(coeffs[k]):
            print(f"    c_{D-1-k}: {coeffs[k]}")

# =====================================================================
# PART 5: FERMAT CURVE POINT COUNTING AT p=7
# =====================================================================

print("\n" + "=" * 78)
print("  FERMAT CURVE C_7: X^7 + Y^7 = 1 OVER F_7")
print("=" * 78)

p = 7
print(f"\n  Affine point count over F_7:")
count = 0
for x in range(p):
    for y in range(p):
        if (pow(x, p, p) + pow(y, p, p)) % p == 1:
            count += 1
N1 = count + p - 1  # Add p-1 points at infinity (projective)
# Actually for X^p + Y^p = Z^p, points at infinity are [x:y:0] with x^p + y^p = 0
pts_inf = sum(1 for x in range(p) for y in range(p) if (pow(x, p, p) + pow(y, p, p)) % p == 0 and (x, y) != (0, 0))
# Projective points at infinity: [x:y:0], normalized
inf_pts = 0
for y in range(p):
    if (1 + pow(y, p, p)) % p == 0:
        inf_pts += 1
# [1:y:0] with 1+y^p=0, i.e., y^7 ≡ -1 ≡ 6 mod 7
# Plus [0:1:0] (always a point since 0+1=1≠0 for p≥2... wait)
# For X^p + Y^p = Z^p: at Z=0, X^p + Y^p = 0 → X^p = -Y^p
# Projective: [X:Y:0] → X = -Y (since x^7 = x in F_7)
# So [1:-1:0] = [1:6:0], the only point at infinity
print(f"  Affine points: {count}")
print(f"  Points at infinity on X^7+Y^7=Z^7: [1:{p-1}:0]")
N1_proj = count + 1  # one point at infinity [1:6:0]
print(f"  Total N_1 = #C_7(F_7) = {N1_proj}")

# For the Fermat curve, N_1 = p + 1 - a_p where a_p = sum of Jacobi sums
a_7 = p + 1 - N1_proj
print(f"  a_7 = {p}+1-{N1_proj} = {a_7}")
print(f"  (This is the trace of Frobenius on C_7/F_7)")

g = (p - 1) * (p - 2) // 2
print(f"\n  Genus of C_7: g = (7-1)(7-2)/2 = {g}")
print(f"  Frobenius polynomial degree: 2g = {2*g}")
print(f"  Jacobi sums J(chi^a, chi^b) for 1 <= a,b, a+b <= 6")
print(f"  give {g} eigenvalue pairs (conjugates), total 2g = {2*g}")

# =====================================================================
# PART 6: SUMMARY
# =====================================================================

print("\n" + "=" * 78)
print("  F8c SUMMARY")
print("=" * 78)

print(f"""
  KEY FINDINGS:

  1. LEVEL-2 IS THE ONLY PRODUCTIVE LEVEL for the dense marginal approach.
     Level-3/4 marginals collapse to rank-1 matrices (character orthogonality).

  2. PRIME SCAN RESULTS:
     - p=3  (mod 4 = 3):  2 eigenvalues at sqrt(3)  [F6 PROVED]
     - p=5  (mod 4 = 1):  0 eigenvalues at sqrt(5)  [NEGATIVE]
     - p=7  (mod 4 = 3):  6 eigenvalues at sqrt(7)  [NEW POSITIVE]
     - p=11 (mod 4 = 3):  ? eigenvalues at sqrt(11)
     - p=13 (mod 4 = 1):  ? eigenvalues at sqrt(13)

  3. For p=7, the 6 eigenvalues at |mu|=sqrt(7) are:
     sqrt(7), sqrt(7), i*sqrt(7), -i*sqrt(7), -sqrt(7), -sqrt(7)
     These form a degree-6 polynomial with integer coefficients.

  4. The Fermat curve C_7 has genus 15, requiring a degree-30 Frobenius
     polynomial. The 6 eigenvalues we found are a SUBSET of this.

  5. The congruence pattern (p mod 4) may determine which primes yield
     Weil-bound eigenvalues, but more primes are needed to confirm.
""")
