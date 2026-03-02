# carry-frobenius

**Frobenius Eigenvalues and Gauss Sums from Witt Vector Carry Arithmetic**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

This repository contains the paper and computational verification for the discovery that the Frobenius eigenvalues and Gauss sums of algebraic varieties over finite fields are *algebraically encoded* in the carry structure of Witt vector multiplication.

## Main Results

**Theorem A** (*p* = 3). The 9 × 9 multiplicative Witt carry matrix contains the exact Frobenius polynomial μ² + 3μ + 3 of the supersingular elliptic curve *E*: y² = x³ + 2x + 1 over F₃.

**Theorem B** (*p* = 7). The 49 × 49 carry matrix yields a genus-3 Frobenius factor (μ² − 7)²(μ² + 7) with 6 eigenvalues at |μ| = √7. The supersingular component μ² + 7 has roots ±i√7, which are exactly the quadratic Gauss sum g(χ³) of the Legendre symbol.

The carry exponent c₂ is equidistributed on F*_p, and Weil-bound eigenvectors are orthogonal to the all-ones vector: the Frobenius is invisible in the total exponential sum but emerges in Tr(M^k). We conjecture the carry operator computes the Frobenius on the Artin-Schreier variety z^p − z = c₂(x₀, x₁, y₀, y₁).

## Status

- **Proved:** Theorem A (`p=3`) and Theorem B (`p=7`) factor extractions.
- **Open target:** Artin-Schreier carry conjecture for full cross-prime factor coverage.

The operator language is compatible with the P1 Step-3 character-resolvent framework (real `p=2` channel vs phase-bearing `p>=3` Frobenius pairs). The base-2 Witt–carry bridge is developed further in [L].

## Cross-Prime Evidence

| *p* | genus | 2*g* | Weil eigs | Method | Status |
|-----|-------|------|-----------|--------|--------|
| 2   | 0     | 0    | 0         | —      | Ihara zeta (golden ratio) |
| 3   | 1     | 2    | 2         | untwisted | **Proved** (Theorem A) |
| 5   | 6     | 12   | 2         | χ₂ twist | Partial (§5.4) |
| 7   | 15    | 30   | 6         | untwisted | **Proved** (Theorem B) |
| 11  | 45    | 90   | 2         | χ₂ ord 5  | Partial (§5.8) |
| 13  | 66    | 132  | 2         | χ₄ ord 3  | Partial (§5.8) |
| 17  | 120   | 240  | 14        | multiple  | Rich (§5.8) |
| 19  | 153   | 306  | 23+       | untwisted + twists | Rich (§5.8) |
| 23  | 231   | 462  | 14+       | ord 11    | Partial (§5.8) |

Every prime tested produces Weil eigenvalues with appropriate character twists. No single character gives a coherent cross-prime *L*-function (§5.8, §8.2).

## Repository Structure

```
paper/carry_frobenius.md                    The paper
experiments/F01_verify_theorem_a.py         Computational proof of Theorem A
experiments/F11_artin_schreier_points.py    Artin-Schreier analysis and carry balance
experiments/F12_jacobi_sum_matcher.py       Gauss/Jacobi sum matching
experiments/F13_twist_p11_p13.py            Extended twist scan (p=3,5,7,11,13)
experiments/F14_twisted_euler_product.py    Cross-prime coherence (p=3..23)
experiments/F15c_complex_analysis.py        D-odd binary carry: complex eigenvalue analysis
experiments/F15d_analytical_limit.py        Analytical derivation of limit matrix S∞
experiments/F15g_carry_formula_proof.py     Formal proof of O(1) carry formula
CITATION.cff                                Citation metadata
LICENSE                                     CC BY 4.0 (paper) + MIT (code)
```

## Reproduction

```bash
pip install numpy sympy
python experiments/F01_verify_theorem_a.py    # Theorem A proof
python experiments/F12_jacobi_sum_matcher.py  # Gauss sum identification
python experiments/F14_twisted_euler_product.py  # Cross-prime scan (p=3..23)
```

F01 verifies the characteristic polynomial factorization, Weil bound, and point counts. F12 matches carry eigenvalues against Gauss/Jacobi sums. F15d derives the exact limit matrix analytically. F15g proves the closed-form carry formula (2.16M pairs, 0 mismatches).

## Dependencies

- Python ≥ 3.9
- [NumPy](https://numpy.org/) (numerical eigenvalue verification)
- [SymPy](https://www.sympy.org/) (polynomial factorization)

All core computations use exact arithmetic (`fractions.Fraction` and manual Q(ω) representation). NumPy is used only for independent numerical cross-checks.

### Citation

```bibtex
@article{alimonti2026witt_frobenius,
  author  = {Alimonti, Stefano},
  title   = {Frobenius Eigenvalues and Gauss Sums from Witt Vector Carry Arithmetic},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-frobenius}
}
```

## License

Paper: CC BY 4.0. Code: MIT — see [LICENSE](LICENSE).
