"""
F15g: Formal Proof of the Closed-Form Carry Formula
====================================================

THEOREM. For K-bit integers X, Y with X, Y ∈ [2^{K-1}, 2^K),
let Q = X·Y and D = 2K-1.

If Q has exactly D bits (i.e., 2^{D-1} ≤ Q < 2^D, the "D-odd" condition),
then the carry entering position D-2 = 2K-3 of the schoolbook multiplication is:

    carry_{D-2} = ⌊Q / 2^{D-2}⌋ - x_{K-2} - y_{K-2} - 2

where x_{K-2} = ⌊X/2^{K-2}⌋ mod 2 and y_{K-2} = ⌊Y/2^{K-2}⌋ mod 2.

PROOF. Consider the schoolbook multiplication of two K-bit integers.

Let D = 2K-1. The product Q = X·Y can be written as:
    Q = Σ_{j=0}^{D-1} q_j · 2^j

where q_j are the bits of Q.

The carry entering position j of the addition is defined recursively:
    c_0 = 0
    For j ≥ 1: c_j = ⌊(s_j + c_{j-1}) / 2⌋

where s_j = Σ_{i+k=j} x_i · y_k is the partial product sum at position j.

But by the definition of binary representation:
    Q = Σ_{j=0}^{D-1} q_j · 2^j = Σ_{j=0}^{D-1} s_j · 2^j

The key observation: the carry c_j and bit q_j at position j satisfy:
    s_j + c_{j-1} = q_j + 2·c_j

Summing from j=0 to j=J-1:
    Σ_{j=0}^{J-1} s_j + Σ_{j=0}^{J-1} c_{j-1} = Σ_{j=0}^{J-1} q_j + 2·Σ_{j=0}^{J-1} c_j

Note: Σ c_{j-1} from j=0..J-1 = c_{-1} + c_0 + ... + c_{J-2} = Σ c_j from j=-1..J-2.
With c_{-1} = 0: = Σ c_j from j=0..J-2.

So: Σ s_j + Σ_{j=0}^{J-2} c_j = Σ q_j + 2·Σ_{j=0}^{J-1} c_j

Rearranging: Σ s_j - Σ q_j = 2·Σ_{j=0}^{J-1} c_j - Σ_{j=0}^{J-2} c_j
            = Σ_{j=0}^{J-2} c_j + 2·c_{J-1}
            = c_{J-1} + (c_{J-1} + Σ_{j=0}^{J-2} c_j)

Actually, let me use a cleaner approach.

DIRECT APPROACH:
    ⌊Q / 2^J⌋ = Σ_{j=J}^{D-1} q_j · 2^{j-J}

For the top bits of Q:
    ⌊Q / 2^{D-2}⌋ = q_{D-1} · 2 + q_{D-2} = 2 + q_{D-2}

since q_{D-1} = 1 (MSB of D-bit number).

Now, at position D-2 = 2K-3:
    s_{D-2} + c_{D-3} = q_{D-2} + 2·c_{D-2}

where s_{D-2} = Σ_{i+k=D-2} x_i · y_k.

For K-bit integers, the non-zero terms have i ∈ [0, K-1] and k = D-2-i ∈ [0, K-1].
Since D-2 = 2K-3, we need both i ≤ K-1 and 2K-3-i ≤ K-1, i.e., i ≥ K-2.
So i ∈ {K-2, K-1}, giving:

    s_{D-2} = x_{K-2}·y_{K-1} + x_{K-1}·y_{K-2}

Since x_{K-1} = y_{K-1} = 1 (MSBs of K-bit numbers):
    s_{D-2} = y_{K-2} + x_{K-2}

At position D-1 = 2K-2:
    s_{D-1} = x_{K-1}·y_{K-1} = 1
    s_{D-1} + c_{D-2} = q_{D-1} + 2·c_{D-1}
    1 + c_{D-2} = 1 + 2·c_{D-1}    (since q_{D-1} = 1)
    c_{D-1} = c_{D-2} / 2

For D-odd (exactly D bits): c_{D-1} = 0, so c_{D-2} ∈ {0, 1}.
(And c_{D-2} = 0 iff Q < 2^D, which is guaranteed by D-odd.)

Wait, c_{D-1} is the carry OUT of the last position. For a D-bit result, c_{D-1} = 0.
So c_{D-2} must be 0 or 1, and from 1 + c_{D-2} = 1 + 0, c_{D-2} = 0... 

No wait. c_{D-1} is the carry INTO position D-1, and c_{D} is the carry OUT.
Let me be more precise.

Position D-2: s_{D-2} + c_{D-3} = q_{D-2} + 2·c_{D-2}
Position D-1: s_{D-1} + c_{D-2} = q_{D-1} + 2·c_{D-1}

c_{D-1} is the final carry out. For a D-bit result: c_{D-1} = 0.

From position D-1: 1 + c_{D-2} = 1 + 0, so c_{D-2} = 0?

That can't be right — c_{D-2} is the carry entering position D-1.

Hmm, let me recheck. If s_{D-1} = 1 and c_{D-2} is the carry into position D-1:
    1 + c_{D-2} = q_{D-1} + 2·c_{D-1}

For D-odd: the result has exactly D bits, so q_{D-1} = 1 and c_{D-1} = 0.
    1 + c_{D-2} = 1 + 0 = 1
    c_{D-2} = 0

But our data shows carry_{D-2} is NOT always 0! There's a confusion about indexing.

Let me clarify. The "carry at position D-2" in our code is carry ENTERING position D-2,
i.e., c_{D-3} (carry propagated from position D-3 to position D-2).

In the code: for each position j, we compute carry_out = (s_j + carry_in) >> 1.
The carry_in for position D-2 is the carry_out from position D-3.

So carry_{D-2} in our observable is actually c_{D-3} (the carry entering position D-2),
NOT c_{D-2} (the carry entering position D-1).

Let me redo with correct indexing.

CORRECTED NOTATION:
Let c_j = carry entering position j (carry out of position j-1). c_0 = 0.
At position j: s_j + c_j = q_j + 2·c_{j+1}

At position D-2 = 2K-3:
    s_{D-2} + c_{D-2} = q_{D-2} + 2·c_{D-1}

s_{D-2} = x_{K-2} + y_{K-2} (as computed above).

At position D-1 = 2K-2:
    s_{D-1} + c_{D-1} = q_{D-1} + 2·c_D
    1 + c_{D-1} = q_{D-1} + 2·c_D

D-odd: q_{D-1} = 1 and c_D = 0 (no bits beyond D-1).
So: 1 + c_{D-1} = 1, hence c_{D-1} = 0.

From position D-2:
    (x_{K-2} + y_{K-2}) + c_{D-2} = q_{D-2} + 2·0 = q_{D-2}

So: c_{D-2} = q_{D-2} - x_{K-2} - y_{K-2}

And q_{D-2} is the (D-2)-th bit of Q = the (2K-3)-th bit:
    q_{D-2} = ⌊Q / 2^{D-2}⌋ mod 2

But also: ⌊Q / 2^{D-2}⌋ = 2·q_{D-1} + q_{D-2} = 2 + q_{D-2}

So: q_{D-2} = ⌊Q / 2^{D-2}⌋ - 2

Therefore: c_{D-2} = ⌊Q / 2^{D-2}⌋ - 2 - x_{K-2} - y_{K-2}    ∎

This is EXACTLY the closed-form formula! QED.

Now verify numerically:
"""

import numpy as np

print("=" * 70)
print("  Formal Verification of Carry Formula")
print("=" * 70)
print()
print("  THEOREM: For K-bit X, Y with D-odd product Q = X·Y:")
print("    c_{D-2} = ⌊Q / 2^{D-2}⌋ - x_{K-2} - y_{K-2} - 2")
print()
print("  PROOF SKETCH:")
print("    1. s_{D-2} = x_{K-2} + y_{K-2} (only 2 terms in convolution)")
print("    2. s_{D-1} = x_{K-1}·y_{K-1} = 1 (both MSBs are 1)")
print("    3. D-odd ⟹ c_{D-1} = 0 (carry into top position is zero)")
print("    4. Position D-2: s_{D-2} + c_{D-2} = q_{D-2} + 2·c_{D-1}")
print("       ⟹ c_{D-2} = q_{D-2} - x_{K-2} - y_{K-2}")
print("    5. q_{D-2} = ⌊Q/2^{D-2}⌋ - 2 (since q_{D-1} = 1)")
print("    6. Combining: c_{D-2} = ⌊Q/2^{D-2}⌋ - 2 - x_{K-2} - y_{K-2}  ∎")

# Numerical verification
print(f"\n  Numerical verification (brute-force vs formula):")

errors = 0
total = 0
for K in range(4, 13):
    lo = 1 << (K - 1)
    hi = 1 << K
    D = 2 * K - 1

    for X in range(lo, hi):
        x_bits = [(X >> b) & 1 for b in range(K)]
        for Y in range(lo, hi):
            Q = X * Y
            if Q.bit_length() != D:
                continue

            carry = 0
            target_carry = -1
            for j in range(D):
                if j == D - 2:
                    target_carry = carry
                s = 0
                for i in range(max(0, j - K + 1), min(j + 1, K)):
                    s += x_bits[i] * ((Y >> (j - i)) & 1)
                carry = (s + carry) >> 1

            shift = D - 2
            formula_carry = (Q >> shift) - ((X >> (K-2)) & 1) - ((Y >> (K-2)) & 1) - 2

            if target_carry != formula_carry:
                errors += 1
                if errors <= 3:
                    print(f"    MISMATCH: K={K}, X={X}, Y={Y}, brute={target_carry}, formula={formula_carry}")
            total += 1
    print(f"    K={K:2d}: {total:8d} pairs verified, {errors} errors")

print(f"\n  Total verified: {total} D-odd pairs, {errors} mismatches")
print(f"  VERDICT: {'PROOF VERIFIED ✓' if errors == 0 else 'ERRORS FOUND ✗'}")

# =====================================================================
# Prove the range of c_{D-2}
# =====================================================================
print(f"\n{'='*70}")
print(f"  Range analysis of c_{{D-2}}")
print(f"{'='*70}")
print()
print("  Since x_{K-2}, y_{K-2} ∈ {0, 1} and ⌊Q/2^{D-2}⌋ ∈ {2, 3}:")
print("  (⌊Q/2^{D-2}⌋ = 2 when q_{D-2}=0, = 3 when q_{D-2}=1)")
print()
print("  c_{D-2} = ⌊Q/2^{D-2}⌋ - 2 - x_{K-2} - y_{K-2}")
print("          ∈ {0-2, 0-1, 0-0, 1-2, 1-1, 1-0}")
print("          = {-2, -1, 0, -1, 0, 1}")
print()
print("  Wait — can c_{D-2} be negative or > 1?")
print("  By definition, carry is always 0 or 1 in binary.")
print("  Let's verify the range...")

carries_seen = set()
K = 10
lo = 1 << (K - 1)
hi = 1 << K
D = 2 * K - 1
for X in range(lo, hi):
    for Y in range(lo, hi):
        Q = X * Y
        if Q.bit_length() != D:
            continue
        c = (Q >> (D-2)) - ((X >> (K-2)) & 1) - ((Y >> (K-2)) & 1) - 2
        carries_seen.add(c)

print(f"  Observed carry values: {sorted(carries_seen)}")
print(f"  Carry is always in {{0, 1}} ✓" if carries_seen <= {0, 1} else
      f"  UNEXPECTED: carry values outside {{0, 1}}!")

print()
print("  Since c_{D-2} ∈ {0, 1}, the formula constraints are:")
print("    c_{D-2} = 0: ⌊Q/2^{D-2}⌋ = 2 + x_{K-2} + y_{K-2}")
print("    c_{D-2} = 1: ⌊Q/2^{D-2}⌋ = 3 + x_{K-2} + y_{K-2} − 2 = 1 + x_{K-2} + y_{K-2}")
print()
print("  This means:")
print("    ⌊Q/2^{D-2}⌋ ∈ {2, 3} for D-odd Q")
print("    When x_{K-2} + y_{K-2} = 0: c = ⌊Q/2^{D-2}⌋ - 2 ∈ {0, 1}")
print("    When x_{K-2} + y_{K-2} = 1: c = ⌊Q/2^{D-2}⌋ - 3 ∈ {-1, 0}")
print("    When x_{K-2} + y_{K-2} = 2: c = ⌊Q/2^{D-2}⌋ - 4 ∈ {-2, -1}")
print()
print("  The last two give c < 0, which is impossible.")
print("  This means: when x_{K-2} + y_{K-2} ≥ 1 AND D-odd,")
print("  ⌊Q/2^{D-2}⌋ must be large enough to keep c ∈ {0,1}.")
print("  Specifically: ⌊Q/2^{D-2}⌋ ≥ 2 + x_{K-2} + y_{K-2}")
print()
print("  This is automatically satisfied by the structure of multiplication.")

print(f"\n{'='*70}")
print(f"  FORMAL PROOF COMPLETE")
print(f"{'='*70}")
print(f"""
  THEOREM (Closed-Form Carry at D-2).
  Let X, Y be K-bit positive integers (x_{{K-1}} = y_{{K-1}} = 1).
  Let Q = X·Y with exactly D = 2K-1 bits (D-odd condition).
  Then the carry entering position D-2 of schoolbook multiplication is:

      c_{{D-2}} = ⌊Q / 2^{{D-2}}⌋ - x_{{K-2}} - y_{{K-2}} - 2

  where x_{{K-2}} and y_{{K-2}} are the second-most-significant bits of X and Y.

  PROOF. At position D-2 = 2K-3 of the schoolbook algorithm:
    (a) The partial product sum s_{{D-2}} = x_{{K-2}} + y_{{K-2}}
        (only terms i=K-2, k=K-1 and i=K-1, k=K-2 contribute).
    (b) The carry into position D-1 is c_{{D-1}} = 0
        (from position D-1: 1 + c_{{D-1}} = 1 + 2·c_D, and c_D = 0 for D-bit result).
    (c) From the recurrence at D-2: s_{{D-2}} + c_{{D-2}} = q_{{D-2}} + 2·c_{{D-1}} = q_{{D-2}}.
    (d) The bit q_{{D-2}} = ⌊Q/2^{{D-2}}⌋ - 2 (since q_{{D-1}} = 1 for D-bit Q).
    (e) Substituting: c_{{D-2}} = ⌊Q/2^{{D-2}}⌋ - 2 - x_{{K-2}} - y_{{K-2}}.  ∎

  COROLLARY. The observable val = 2c_{{D-2}} - 1 ∈ {{-1, +1}} is:
      val = 2⌊Q/2^{{D-2}}⌋ - 2x_{{K-2}} - 2y_{{K-2}} - 5

  COMPUTATIONAL IMPACT.
  This replaces an O(K²) carry propagation with an O(1) formula:
  one right-shift, two bit-extractions, and three subtractions.
  Enables computation to K=18+ in seconds (vs hours for brute-force).
""")
