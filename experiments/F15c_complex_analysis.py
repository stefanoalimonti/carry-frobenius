"""
F15c: D-odd Binary Carry Eigenvalues — Full Complex Analysis
=============================================================

The D-odd conditioned 4×4 binary carry matrix has X↔Y symmetry: M[1,j] ↔ M[2,j].
This decomposes into:
  - Symmetric sector: 2×2 reduced matrix S (complex conjugate eigenvalues)
  - Anti-symmetric sector: scalar eigenvalue f - g (real)
  - Zero eigenvalue from the null row (state a1=c1=1 never D-odd)

CLOSED-FORM for carry_{D-2}:
    carry_{D-2} = (Q >> (2K-3)) - (X >> (K-2) & 1) - (Y >> (K-2) & 1) - 2
"""

import numpy as np
import math
import time
import sys

PI = math.pi
PHI = (1 + math.sqrt(5)) / 2

MAX_K = 16
if len(sys.argv) > 1:
    MAX_K = int(sys.argv[1])

# =====================================================================
# PART 1: Compute corrected matrices and FULL complex eigenvalues
# =====================================================================
print("=" * 78)
print(f"  Part 1: D-odd carry matrices K=4..{MAX_K}")
print("=" * 78)

results = []

for K in range(4, MAX_K + 1):
    t0 = time.time()
    lo = 1 << (K - 1)
    hi = 1 << K
    D = 2 * K - 1
    shift_dm2 = 2 * K - 3
    shift_a0 = K - 3
    shift_a1 = K - 2

    threshold = np.uint64(1) << np.uint64(D)
    half_threshold = np.uint64(1) << np.uint64(D - 1)

    M_flat = np.zeros(16, dtype=np.float64)
    N_flat = np.zeros(16, dtype=np.float64)

    Y_all = np.arange(lo, hi, dtype=np.uint64)

    for X_int in range(lo, hi):
        X = np.uint64(X_int)
        Q = X * Y_all

        dodd = (Q >= half_threshold) & (Q < threshold)
        Q_d = Q[dodd]
        Y_d = Y_all[dodd]
        if len(Q_d) == 0:
            continue

        a0 = int((X_int >> shift_a0) & 1)
        a1 = int((X_int >> shift_a1) & 1)

        c0 = ((Y_d >> np.uint64(shift_a0)) & np.uint64(1)).astype(np.int32)
        c1 = ((Y_d >> np.uint64(shift_a1)) & np.uint64(1)).astype(np.int32)

        carry = (Q_d >> np.uint64(shift_dm2)).astype(np.int32) - a1 - c1 - 2
        val = (2 * carry - 1).astype(np.float64)

        cell = ((a1 * 2 + c1) * 4 + a0 * 2 + c0).astype(np.intp)

        M_flat += np.bincount(cell, weights=val, minlength=16)
        N_flat += np.bincount(cell, minlength=16).astype(np.float64)

    M = M_flat.reshape(4, 4)
    N = N_flat.reshape(4, 4)

    for col in range(4):
        tc = N[:, col].sum()
        if tc > 0:
            M[:, col] /= tc

    eigs_full = np.linalg.eigvals(M)
    eigs_full_sorted = sorted(eigs_full, key=lambda z: -abs(z))

    # Decompose using X↔Y symmetry
    a_val = M[0, 0]
    b_val = M[0, 1]  # = M[0,2] by symmetry
    d_val = M[0, 3]
    e_val = M[1, 0]  # = M[2,0]
    f_val = M[1, 1]  # diagonal
    g_val = M[1, 2]  # off-diagonal

    # Symmetric sector: 2x2 with eigenvector [x, y, y, 0]
    S = np.array([[a_val, 2 * b_val],
                  [e_val, f_val + g_val]])
    eigs_sym = np.linalg.eigvals(S)
    tr_S = S[0, 0] + S[1, 1]
    det_S = S[0, 0] * S[1, 1] - S[0, 1] * S[1, 0]
    disc = tr_S**2 - 4 * det_S

    # Anti-symmetric eigenvalue
    lam_anti = f_val - g_val

    elapsed = time.time() - t0
    n_dodd = int(N.sum())

    results.append({
        'K': K,
        'M': M.copy(),
        'N': N.copy(),
        'eigs_full': eigs_full_sorted,
        'eigs_sym': eigs_sym,
        'S': S.copy(),
        'tr_S': tr_S,
        'det_S': det_S,
        'disc': disc,
        'lam_anti': lam_anti,
        'n_dodd': n_dodd,
        'elapsed': elapsed,
    })

# =====================================================================
# PART 2: Report eigenvalues with complex structure
# =====================================================================
print("\n" + "=" * 78)
print("  Part 2: Full eigenvalue spectrum (complex!)")
print("=" * 78)

print(f"\n  {'K':>3s}  {'Re(lam_dom)':>12s}  {'Im(lam_dom)':>12s}  "
      f"{'|lam_dom|':>12s}  {'arg(deg)':>10s}  {'lam_anti':>12s}  "
      f"{'disc':>12s}  {'n_dodd':>12s}")

for r in results:
    K = r['K']
    e = r['eigs_sym']
    # Dominant symmetric eigenvalue
    idx_dom = np.argmax(np.abs(e))
    lam_dom = e[idx_dom]
    re_dom = lam_dom.real
    im_dom = lam_dom.imag
    mod_dom = abs(lam_dom)
    arg_dom = math.degrees(math.atan2(im_dom, re_dom)) if mod_dom > 1e-10 else 0

    lam_a = r['lam_anti']
    disc = r['disc']

    print(f"  {K:3d}  {re_dom:+12.8f}  {im_dom:+12.8f}  "
          f"{mod_dom:12.8f}  {arg_dom:+10.4f}  {lam_a:+12.8f}  "
          f"{disc:+12.6f}  {r['n_dodd']:12d}")

# =====================================================================
# PART 3: Convergence of modulus, real part, imaginary part, phase
# =====================================================================
print("\n" + "=" * 78)
print("  Part 3: Convergence analysis")
print("=" * 78)

print("\n  --- Modulus |lambda_dom| ---")
mods = []
for r in results:
    e = r['eigs_sym']
    mod = max(abs(e[0]), abs(e[1]))
    mods.append((r['K'], mod))

print(f"  {'K':>3s}  {'|lam|':>14s}  {'Delta':>14s}  {'Ratio':>10s}")
for i, (K, m) in enumerate(mods):
    if i == 0:
        print(f"  {K:3d}  {m:14.10f}")
    elif i == 1:
        d = m - mods[i-1][1]
        print(f"  {K:3d}  {m:14.10f}  {d:+14.10f}")
    else:
        d = m - mods[i-1][1]
        pd = mods[i-1][1] - mods[i-2][1]
        r_val = d / pd if abs(pd) > 1e-15 else float('nan')
        print(f"  {K:3d}  {m:14.10f}  {d:+14.10f}  {r_val:10.6f}")

# Aitken on modulus
if len(mods) >= 3:
    vals = [m for (_, m) in mods]
    print(f"\n  Aitken on |lam| (last 3 points):")
    s0, s1, s2 = vals[-3], vals[-2], vals[-1]
    den = s2 - 2*s1 + s0
    if abs(den) > 1e-15:
        aitken = s0 - (s1 - s0)**2 / den
        print(f"    Aitken limit = {aitken:.10f}")

    if len(vals) >= 5:
        s0, s1, s2 = vals[-5], vals[-4], vals[-3]
        den = s2 - 2*s1 + s0
        if abs(den) > 1e-15:
            aitken = s0 - (s1 - s0)**2 / den
            print(f"    Aitken limit (5 back) = {aitken:.10f}")

print("\n  --- Real part Re(lambda_dom) ---")
reals = []
for r in results:
    e = r['eigs_sym']
    idx_dom = np.argmax(np.abs(e))
    reals.append((r['K'], e[idx_dom].real))

print(f"  {'K':>3s}  {'Re(lam)':>14s}  {'Delta':>14s}  {'Ratio':>10s}")
for i, (K, re_v) in enumerate(reals):
    if i == 0:
        print(f"  {K:3d}  {re_v:+14.10f}")
    elif i == 1:
        d = re_v - reals[i-1][1]
        print(f"  {K:3d}  {re_v:+14.10f}  {d:+14.10f}")
    else:
        d = re_v - reals[i-1][1]
        pd = reals[i-1][1] - reals[i-2][1]
        r_val = d / pd if abs(pd) > 1e-15 else float('nan')
        print(f"  {K:3d}  {re_v:+14.10f}  {d:+14.10f}  {r_val:10.6f}")

# Aitken on Re
if len(reals) >= 3:
    vals = [r for (_, r) in reals]
    s0, s1, s2 = vals[-3], vals[-2], vals[-1]
    den = s2 - 2*s1 + s0
    if abs(den) > 1e-15:
        aitken_re = s0 - (s1 - s0)**2 / den
        print(f"\n  Aitken on Re (last 3): {aitken_re:.10f}")

print("\n  --- Imaginary part Im(lambda_dom) ---")
imags = []
for r in results:
    e = r['eigs_sym']
    idx_dom = np.argmax(np.abs(e))
    imags.append((r['K'], abs(e[idx_dom].imag)))

print(f"  {'K':>3s}  {'|Im(lam)|':>14s}  {'Delta':>14s}  {'Ratio':>10s}")
for i, (K, im_v) in enumerate(imags):
    if i == 0:
        print(f"  {K:3d}  {im_v:14.10f}")
    elif i == 1:
        d = im_v - imags[i-1][1]
        print(f"  {K:3d}  {im_v:14.10f}  {d:+14.10f}")
    else:
        d = im_v - imags[i-1][1]
        pd = imags[i-1][1] - imags[i-2][1]
        r_val = d / pd if abs(pd) > 1e-15 else float('nan')
        print(f"  {K:3d}  {im_v:14.10f}  {d:+14.10f}  {r_val:10.6f}")

if len(imags) >= 3:
    vals = [i for (_, i) in imags]
    s0, s1, s2 = vals[-3], vals[-2], vals[-1]
    den = s2 - 2*s1 + s0
    if abs(den) > 1e-15:
        aitken_im = s0 - (s1 - s0)**2 / den
        print(f"\n  Aitken on |Im| (last 3): {aitken_im:.10f}")

print("\n  --- Anti-symmetric eigenvalue (f - g) ---")
antis = [(r['K'], r['lam_anti']) for r in results]
print(f"  {'K':>3s}  {'lam_anti':>14s}  {'Delta':>14s}")
for i, (K, la) in enumerate(antis):
    if i == 0:
        print(f"  {K:3d}  {la:+14.10f}")
    else:
        d = la - antis[i-1][1]
        print(f"  {K:3d}  {la:+14.10f}  {d:+14.10f}")

if len(antis) >= 3:
    vals = [a for (_, a) in antis]
    s0, s1, s2 = vals[-3], vals[-2], vals[-1]
    den = s2 - 2*s1 + s0
    if abs(den) > 1e-15:
        aitken_anti = s0 - (s1 - s0)**2 / den
        print(f"\n  Aitken on lam_anti (last 3): {aitken_anti:.10f}")

# =====================================================================
# PART 4: Candidate matching for |lam|, Re, |Im|, and lam_anti
# =====================================================================
print("\n" + "=" * 78)
print("  Part 4: Candidate matching")
print("=" * 78)

last_mod = mods[-1][1]
last_re = reals[-1][1]
last_im = imags[-1][1]
last_anti = antis[-1][1]

# Aitken-extrapolated values
aitken_values = {}
for name, series in [('mod', mods), ('re', reals), ('im', imags), ('anti', antis)]:
    vals = [v for (_, v) in series]
    if len(vals) >= 3:
        s0, s1, s2 = vals[-3], vals[-2], vals[-1]
        den = s2 - 2*s1 + s0
        if abs(den) > 1e-15:
            aitken_values[name] = s0 - (s1 - s0)**2 / den

candidates_mod = [
    ("1/sqrt(5) = phi-1", 1/math.sqrt(5)),
    ("sqrt(1/5)", 1/math.sqrt(5)),
    ("(sqrt(5)-1)/2 / phi", (math.sqrt(5)-1)/2 / PHI),
    ("1/phi^2 = 2-phi", 2 - PHI),
    ("1/(1+phi) = (3-sqrt(5))/2", (3 - math.sqrt(5))/2),
    ("phi/2 - 1/phi", PHI/2 - 1/PHI),
    ("sqrt(phi/pi)", math.sqrt(PHI/PI)),
    ("3/(2*pi)", 3/(2*PI)),
    ("ln(phi)", math.log(PHI)),
    ("1/phi = phi-1", PHI - 1),
    ("pi/7", PI/7),
    ("sqrt(2)/pi", math.sqrt(2)/PI),
    ("(1+phi)/5", (1+PHI)/5),
]

candidates_re = [
    ("-3/10", -0.3),
    ("-1/pi", -1/PI),
    ("-pi/10", -PI/10),
    ("-3/pi^2", -3/PI**2),
    ("-ln(2)", -math.log(2)),
    ("-1/e", -1/math.e),
    ("-ln(phi)", -math.log(PHI)),
    ("-(3-sqrt(5))/2", -(3-math.sqrt(5))/2),
    ("-phi/2+1/phi", -PHI/2 + 1/PHI),
    ("-(2-phi)", -(2-PHI)),
]

candidates_anti = [
    ("-1/2*ln(2)", -math.log(2)/2),
    ("-pi/20", -PI/20),
    ("-1/(2*pi)", -1/(2*PI)),
    ("-(3-sqrt(5))/4", -(3-math.sqrt(5))/4),
    ("-3/20", -0.15),
    ("-(phi-1)/4", -(PHI-1)/4),
]

for label, val_best, cands in [
    ("|lambda_dom|", aitken_values.get('mod', last_mod), candidates_mod),
    ("Re(lambda_dom)", aitken_values.get('re', last_re), candidates_re),
    ("lam_anti", aitken_values.get('anti', last_anti), candidates_anti),
]:
    print(f"\n  --- {label} (best estimate: {val_best:.10f}) ---")
    ranked = sorted(cands, key=lambda x: abs(x[1] - abs(val_best)) if 'mod' in label.lower() else abs(x[1] - val_best))
    for name, val in ranked[:5]:
        delta = val - val_best
        rel = abs(delta / val_best) * 100 if abs(val_best) > 1e-15 else float('inf')
        print(f"    {name:>30s}  {val:+14.10f}  delta={delta:+14.10f}  ({rel:.4f}%)")

# =====================================================================
# PART 5: Phase angle analysis
# =====================================================================
print("\n" + "=" * 78)
print("  Part 5: Phase angle of complex eigenvalue")
print("=" * 78)

print(f"\n  {'K':>3s}  {'phase (deg)':>12s}  {'phase/pi':>12s}")
for r in results:
    K = r['K']
    e = r['eigs_sym']
    if max(abs(e[0].imag), abs(e[1].imag)) < 1e-10:
        print(f"  {K:3d}  real eigenvalues (disc > 0)")
        continue
    idx_dom = np.argmax(np.abs(e))
    lam = e[idx_dom]
    phase_rad = math.atan2(abs(lam.imag), lam.real)
    phase_deg = math.degrees(phase_rad)
    phase_over_pi = phase_rad / PI
    print(f"  {K:3d}  {phase_deg:12.6f}  {phase_over_pi:12.8f}")

# =====================================================================
# PART 6: det(S) = |lambda|^2 analysis
# =====================================================================
print("\n" + "=" * 78)
print("  Part 6: det(S) = |lambda_dom|^2 convergence")
print("=" * 78)

dets = [(r['K'], r['det_S']) for r in results]
print(f"\n  {'K':>3s}  {'det(S)':>14s}  {'|lam|^2':>14s}  {'sqrt(det)':>14s}")
for K, d in dets:
    print(f"  {K:3d}  {d:14.10f}  {d:14.10f}  {math.sqrt(abs(d)):14.10f}")

det_candidates = [
    ("1/5", 0.2),
    ("1/phi^2 = 3-phi", 3 - PHI),
    ("(3-sqrt(5))/2", (3-math.sqrt(5))/2),
    ("1/(2*phi)", 1/(2*PHI)),
    ("1/e^2", 1/math.e**2),
    ("phi/8", PHI/8),
    ("pi/15", PI/15),
    ("2/pi^2", 2/PI**2),
    ("ln(2)/pi", math.log(2)/PI),
]

det_best = aitken_values.get('mod', mods[-1][1])**2
print(f"\n  Best estimate det(S) -> {det_best:.10f}")
print(f"  {'Candidate':>25s}  {'Value':>14s}  {'delta':>14s}")
for name, val in sorted(det_candidates, key=lambda x: abs(x[1] - det_best)):
    delta = val - det_best
    print(f"  {name:>25s}  {val:14.10f}  {delta:+14.10f}")

print("\n" + "=" * 78)
print("  Part 7: 2x2 Symmetric sector matrix at largest K")
print("=" * 78)

r = results[-1]
print(f"\n  K={r['K']}:")
S = r['S']
print(f"    S = [[{S[0,0]:+.8f}  {S[0,1]:+.8f}]")
print(f"         [{S[1,0]:+.8f}  {S[1,1]:+.8f}]]")
print(f"    tr(S) = {r['tr_S']:.10f}")
print(f"    det(S) = {r['det_S']:.10f}")
print(f"    disc = {r['disc']:.10f}")
print(f"    Eigenvalues: {r['eigs_sym'][0]:.10f}, {r['eigs_sym'][1]:.10f}")
print(f"    |lam| = {abs(r['eigs_sym'][0]):.10f}")
print(f"    Anti-sym lam = {r['lam_anti']:.10f}")

print("\n" + "=" * 78)
print("  DONE")
print("=" * 78)
