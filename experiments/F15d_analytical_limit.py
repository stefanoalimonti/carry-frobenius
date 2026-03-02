"""
F15d: Analytical Derivation of the Limit Matrix S∞
===================================================

As K → ∞, X = 2^{K-1} U, Y = 2^{K-1} V with U, V uniform on [1,2).
D-odd condition: UV < 2.

Bits partition [1,2) into subintervals:
  I(0,0) = [1, 5/4),   I(0,1) = [5/4, 3/2),
  I(1,0) = [3/2, 7/4), I(1,1) = [7/4, 2)
where I(a1, a0) is the interval for (bit K-2 = a1, bit K-3 = a0).

The val observable:
  val = +1 iff a1 = c1 = 0 AND UV ≥ 3/2
  val = -1 otherwise (for D-odd pairs)
  (a1=c1=1 is never D-odd)

Each matrix entry M[j,i] (to-state j, from-state i) is computed from
definite integrals over the product region.

Goal: compute S∞ in EXACT closed form and identify the eigenvalue constants.
"""

from sympy import (Rational, log, sqrt, pi, Symbol, Matrix, simplify,
                   nsimplify, N, Abs, arg, atan2, I, re, im, conjugate,
                   GoldenRatio)
import numpy as np

phi = GoldenRatio  # (1+sqrt(5))/2

# Interval endpoints: I(a1, a0) = [lo, hi)
def interval(a1, a0):
    """Return (lo, hi) for the bit pattern (a1, a0)."""
    # U = 1 + a1/2 + a0/4 + ... ∈ [1 + a1/2 + a0/4, 1 + a1/2 + a0/4 + 1/4)
    lo = 1 + Rational(a1, 2) + Rational(a0, 4)
    hi = lo + Rational(1, 4)
    return (lo, hi)

# Verify intervals
print("=" * 70)
print("  Interval partition of [1, 2)")
print("=" * 70)
for a1 in [0, 1]:
    for a0 in [0, 1]:
        lo, hi = interval(a1, a0)
        print(f"  I({a1},{a0}) = [{lo}, {hi})")

# =====================================================================
# Compute the integrals
# =====================================================================
# For region U ∈ [u1, u2), V ∈ [v1, v2), condition UV < T:
# Area = ∫_{u1}^{u2} ∫_{v1}^{min(v2, T/u)} dv du
#
# For region UV ∈ [S, T) with same U, V bounds:
# Area = ∫_{u1}^{u2} ∫_{max(v1, S/u)}^{min(v2, T/u)} dv du

def area_uv_less_than(u1, u2, v1, v2, T):
    """Exact area of {(u,v) : u ∈ [u1,u2), v ∈ [v1,v2), uv < T}."""
    # Cutoff: T/u = v2 when u = T/v2
    u_cut = T / v2

    if u_cut >= u2:
        # T/u ≥ v2 for all u in [u1, u2): full rectangle is included
        return (u2 - u1) * (v2 - v1)
    elif u_cut <= u1:
        # T/u ≤ v2 for all u: need T/u > v1, i.e., u < T/v1
        u_top = T / v1
        if u_top <= u1:
            return Rational(0)
        u_eff = min(u2, u_top)
        # ∫_{u1}^{u_eff} (T/u - v1) du = T*ln(u_eff/u1) - v1*(u_eff - u1)
        return T * log(u_eff / u1) - v1 * (u_eff - u1)
    else:
        # Split at u_cut
        part1 = (u_cut - u1) * (v2 - v1)
        # ∫_{u_cut}^{u2} (T/u - v1) du, but need T/u > v1 i.e. u < T/v1
        u_top = T / v1
        u_eff = min(u2, u_top)
        if u_eff <= u_cut:
            part2 = Rational(0)
        else:
            part2 = T * log(u_eff / u_cut) - v1 * (u_eff - u_cut)
        return part1 + part2


def area_uv_in_range(u1, u2, v1, v2, S, T):
    """Exact area of {(u,v) : u ∈ [u1,u2), v ∈ [v1,v2), S ≤ uv < T}."""
    return area_uv_less_than(u1, u2, v1, v2, T) - area_uv_less_than(u1, u2, v1, v2, S)


# =====================================================================
# Build the exact 4×4 matrix
# =====================================================================
print("\n" + "=" * 70)
print("  Computing exact limit matrix entries")
print("=" * 70)

M_num = [[None]*4 for _ in range(4)]   # val-weighted area
N_num = [[None]*4 for _ in range(4)]   # D-odd count area

for i_from in range(4):
    a0 = (i_from >> 1) & 1
    c0 = i_from & 1
    for j_to in range(4):
        a1 = (j_to >> 1) & 1
        c1 = j_to & 1

        u1, u2 = interval(a1, a0)
        v1, v2 = interval(c1, c0)

        # D-odd area: UV < 2
        dodd_area = area_uv_less_than(u1, u2, v1, v2, Rational(2))
        dodd_area = simplify(dodd_area)
        N_num[j_to][i_from] = dodd_area

        if a1 == 0 and c1 == 0:
            # val = +1 when UV ≥ 3/2, val = -1 when UV < 3/2
            area_plus = area_uv_in_range(u1, u2, v1, v2, Rational(3, 2), Rational(2))
            area_minus = area_uv_in_range(u1, u2, v1, v2, Rational(1), Rational(3, 2))
            area_plus = simplify(area_plus)
            area_minus = simplify(area_minus)
            val_weighted = area_plus - area_minus
        elif a1 + c1 >= 2:
            # a1=c1=1: never D-odd
            val_weighted = Rational(0)
        else:
            # a1+c1=1: val = -1 always
            val_weighted = -dodd_area

        M_num[j_to][i_from] = simplify(val_weighted)

# Print raw (unnormalized) areas
print("\n  D-odd areas N[j][i]:")
for j in range(4):
    row = "    ["
    for i in range(4):
        row += f"  {N_num[j][i]}"
        if i < 3: row += ","
    print(row + " ]")

print("\n  Val-weighted areas M_num[j][i]:")
for j in range(4):
    row = "    ["
    for i in range(4):
        row += f"  {M_num[j][i]}"
        if i < 3: row += ","
    print(row + " ]")

# Column totals (for normalization)
col_totals = []
for i in range(4):
    total = sum(N_num[j][i] for j in range(4))
    col_totals.append(simplify(total))
    print(f"\n  Column {i} total: {total} = {N(total, 15)}")

# Total D-odd fraction
total_dodd = sum(col_totals)
print(f"\n  Total D-odd fraction: {simplify(total_dodd)} = {N(total_dodd, 15)}")
print(f"  Expected 2*ln(2) - 1: {N(2*log(2) - 1, 15)}")

# =====================================================================
# Normalized matrix
# =====================================================================
print("\n" + "=" * 70)
print("  Exact normalized matrix M∞")
print("=" * 70)

M_exact = [[None]*4 for _ in range(4)]
for j in range(4):
    row_str = "    ["
    for i in range(4):
        if col_totals[i] == 0:
            M_exact[j][i] = Rational(0)
        else:
            M_exact[j][i] = simplify(M_num[j][i] / col_totals[i])
        row_str += f"  {M_exact[j][i]}"
        if i < 3: row_str += ","
    print(row_str + " ]")

# Numerical verification
print("\n  Numerical values:")
for j in range(4):
    row_str = f"    ["
    for i in range(4):
        row_str += f"  {float(N(M_exact[j][i], 15)):+.10f}"
        if i < 3: row_str += ","
    print(row_str + " ]")

# =====================================================================
# Extract the 2×2 symmetric sector and antisymmetric eigenvalue
# =====================================================================
print("\n" + "=" * 70)
print("  2×2 Symmetric sector S∞")
print("=" * 70)

a_val = M_exact[0][0]
b_val = M_exact[0][1]  # = M_exact[0][2] by symmetry
e_val = M_exact[1][0]  # = M_exact[2][0]
f_val = M_exact[1][1]  # diagonal of rows 1,2
g_val = M_exact[1][2]  # off-diagonal of rows 1,2

print(f"\n  a = M[0,0] = {a_val}")
print(f"  b = M[0,1] = M[0,2] = {b_val}")
print(f"  e = M[1,0] = M[2,0] = {e_val}")
print(f"  f = M[1,1] = M[2,2] = {f_val}")
print(f"  g = M[1,2] = M[2,1] = {g_val}")

# Verify symmetry
print(f"\n  Symmetry check: M[0,1] - M[0,2] = {simplify(M_exact[0][1] - M_exact[0][2])}")
print(f"  Symmetry check: M[1,0] - M[2,0] = {simplify(M_exact[1][0] - M_exact[2][0])}")
print(f"  Symmetry check: M[1,1] - M[2,2] = {simplify(M_exact[1][1] - M_exact[2][2])}")

# 2×2 symmetric matrix
S_a = a_val
S_b = 2 * b_val
S_c = e_val
S_d = f_val + g_val

S = Matrix([[S_a, S_b], [S_c, S_d]])
print(f"\n  S = [[{S_a},  {S_b}],")
print(f"       [{S_c},  {S_d}]]")

tr_S = simplify(S_a + S_d)
det_S = simplify(S_a * S_d - S_b * S_c)
disc_S = simplify(tr_S**2 - 4 * det_S)

print(f"\n  tr(S) = {tr_S}")
print(f"       = {N(tr_S, 20)}")
print(f"\n  det(S) = {det_S}")
print(f"        = {N(det_S, 20)}")
print(f"\n  disc = tr² - 4·det = {disc_S}")
print(f"       = {N(disc_S, 20)}")

# Anti-symmetric eigenvalue
lam_anti = simplify(f_val - g_val)
print(f"\n  λ_anti = f - g = {lam_anti}")
print(f"        = {N(lam_anti, 20)}")

# =====================================================================
# Eigenvalues
# =====================================================================
print("\n" + "=" * 70)
print("  Eigenvalues of S∞")
print("=" * 70)

if N(disc_S) < 0:
    print("\n  Discriminant < 0 → COMPLEX eigenvalues!")
    re_lam = tr_S / 2
    im_lam = sqrt(-disc_S) / 2
    mod_lam = sqrt(det_S)
    phase_lam = atan2(im_lam, re_lam)

    print(f"\n  Re(λ) = tr/2 = {simplify(re_lam)}")
    print(f"        = {N(re_lam, 20)}")
    print(f"\n  |Im(λ)| = √(-disc)/2 = {simplify(im_lam)}")
    print(f"          = {N(im_lam, 20)}")
    print(f"\n  |λ| = √(det) = {simplify(mod_lam)}")
    print(f"      = {N(mod_lam, 20)}")
    print(f"\n  |λ|² = det = {simplify(det_S)}")
    print(f"       = {N(det_S, 20)}")

    print(f"\n  Phase = {N(phase_lam, 20)} rad")
    print(f"        = {N(phase_lam * 180 / pi, 20)} degrees")
    print(f"        = {N(phase_lam / pi, 20)} × π")
else:
    print("\n  Discriminant ≥ 0 → REAL eigenvalues")
    lam1 = (tr_S + sqrt(disc_S)) / 2
    lam2 = (tr_S - sqrt(disc_S)) / 2
    print(f"  λ₁ = {simplify(lam1)} = {N(lam1, 20)}")
    print(f"  λ₂ = {simplify(lam2)} = {N(lam2, 20)}")

# =====================================================================
# Try to identify constants
# =====================================================================
print("\n" + "=" * 70)
print("  Constant identification")
print("=" * 70)

det_val = N(det_S, 30)
tr_val = N(tr_S, 30)
re_val = N(tr_S/2, 30)
anti_val = N(lam_anti, 30)

# Check det against known constants
candidates_det = {
    "1/5": Rational(1, 5),
    "2*ln(2)-1": 2*log(2)-1,
    "ln(2)^2": log(2)**2,
    "pi/15": pi/15,
    "2/pi^2": 2/pi**2,
    "(2*ln2-1)/2": (2*log(2)-1)/2,
    "ln(3/2)": log(Rational(3,2)),
    "2*ln(3/2)": 2*log(Rational(3,2)),
    "ln(2)*ln(3/2)": log(2)*log(Rational(3,2)),
    "ln(9/8)": log(Rational(9,8)),
    "(ln2)^2/2": log(2)**2/2,
    "1-ln(2)": 1-log(2),
    "(1-ln2)^2": (1-log(2))**2,
    "ln(4/3)*ln(3/2)": log(Rational(4,3))*log(Rational(3,2)),
}

print(f"\n  det(S∞) = {det_val}")
for name, val in sorted(candidates_det.items(), key=lambda x: abs(float(N(x[1])) - float(det_val))):
    delta = float(det_val) - float(N(val))
    if abs(delta) < 0.05:
        print(f"    {name:25s} = {N(val, 15):>20s}  delta = {delta:+.10f}")

print(f"\n  tr(S∞)/2 = Re(λ) = {re_val}")
candidates_re = {
    "-3/10": Rational(-3, 10),
    "-ln(2)+1/4": -log(2) + Rational(1,4),
    "-3/pi^2": -3/pi**2,
    "-(2*ln2-1)/2": -(2*log(2)-1)/2,
    "-ln(3/2)": -log(Rational(3,2)),
    "1/2-ln(2)": Rational(1,2) - log(2),
    "-ln(2)/2": -log(2)/2,
}
for name, val in sorted(candidates_re.items(), key=lambda x: abs(float(N(x[1])) - float(re_val))):
    delta = float(re_val) - float(N(val))
    if abs(delta) < 0.05:
        print(f"    {name:25s} = {N(val, 15):>20s}  delta = {delta:+.10f}")

print(f"\n  λ_anti = {anti_val}")
candidates_anti = {
    "-(phi-1)/4": -(phi-1)/4,
    "-pi/20": -pi/20,
    "-1/(2*pi)": -1/(2*pi),
    "-(sqrt(5)-1)/8": -(sqrt(5)-1)/8,
    "-ln(2)/4": -log(2)/4,
    "-3/20": Rational(-3,20),
    "-(2*ln2-1)/4": -(2*log(2)-1)/4,
}
for name, val in sorted(candidates_anti.items(), key=lambda x: abs(float(N(x[1])) - float(anti_val))):
    delta = float(anti_val) - float(N(val))
    if abs(delta) < 0.05:
        print(f"    {name:25s} = {N(val, 15):>20s}  delta = {delta:+.10f}")

# =====================================================================
# Numerical comparison with F15c K=18 data
# =====================================================================
print("\n" + "=" * 70)
print("  Comparison with numerical F15c data (K=18)")
print("=" * 70)

# K=18 values from F15c
k18_re = -0.3017760194
k18_im = 0.3397095065
k18_mod = 0.4543911472
k18_anti = -0.1552633930
k18_det = 0.2064713147
k18_tr = 2 * k18_re

print(f"\n  {'Quantity':>20s}  {'Analytical':>20s}  {'K=18 numerical':>20s}  {'Delta':>15s}")
print(f"  {'':>20s}  {'':>20s}  {'':>20s}  {'':>15s}")
for name, exact, k18 in [
    ("Re(λ)", tr_S/2, k18_re),
    ("|λ|²=det(S)", det_S, k18_det),
    ("tr(S)", tr_S, k18_tr),
    ("λ_anti", lam_anti, k18_anti),
]:
    exact_f = float(N(exact, 20))
    delta = k18 - exact_f
    print(f"  {name:>20s}  {exact_f:+20.15f}  {k18:+20.15f}  {delta:+15.2e}")

print("\n" + "=" * 70)
print("  DONE — Exact analytical matrix S∞ derived")
print("=" * 70)
