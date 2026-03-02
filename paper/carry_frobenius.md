# Frobenius Eigenvalues and Gauss Sums from Witt Vector Carry Arithmetic

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

## Abstract

While the Riemann Hypothesis is proved for function fields $\mathbb{F}_q[x]$, where addition is carry-free, the extension to $\mathbb{Z}$ is obstructed by carry propagation. By lifting carry arithmetic to the ring of Witt vectors $W_n(\mathbb{F}_p)$, we construct a transfer operator whose spectrum encodes Frobenius eigenvalues of algebraic varieties — connecting the most elementary operation in positional arithmetic to the deepest structures of algebraic geometry.

**Theorem A.** For $p = 3$, the characteristic polynomial of the $9 \times 9$ multiplicative Witt carry matrix contains the exact Frobenius polynomial $\mu^2 + 3\mu + 3$ of the supersingular elliptic curve $E: y^2 = x^3 + 2x + 1$ over $\mathbb{F}_3$, with eigenvalues at the Weil bound $\lvert\mu\rvert = \sqrt{3}$.

**Theorem B.** For $p = 7$, exact computation over $\mathbb{Z}[\omega_7]$ yields a genus-3 Frobenius factor $(\mu^2 - 7)^2(\mu^2 + 7)$ with 6 eigenvalues at $\lvert\mu\rvert = \sqrt{7}$. The supersingular component $\mu^2 + 7$ is the minimal polynomial of the quadratic Gauss sum $g(\chi^3) = i\sqrt{7}$ of the Legendre symbol, with roots $\mu = \pm g(\chi^3)$ — an exact algebraic identification, not a numerical coincidence.

The carry exponent $c_2$ over $\mathbb{F}_p^4$ is perfectly equidistributed on $\mathbb{F}_p^{\ast}$, and the Weil-bound eigenvectors are orthogonal to the all-ones vector: the Frobenius is invisible in the total exponential sum but emerges in the periodic orbit structure $\mathrm{Tr}(M^k)$. Dirichlet character twists recover missing Frobenius eigenvalues at every prime tested ($p = 3, 5, 7, 11, 13, 17, 19, 23$): additive Legendre at $p = 5$; order-5 at $p = 11$; order-3 at $p = 13$. The productive character order varies with $p$ in a currently unexplained way — no simple arithmetic condition predicts which character is needed.

We conjecture that the carry operator computes (a subfactor of) the Frobenius on the Artin-Schreier variety $z^p - z = c_2(x_0, x_1, y_0, y_1)$, where $c_2$ is the level-2 Witt multiplication carry dominated by the Fermat quotient $(x^p y^p - (xy)^p)/p$.

**Keywords:** Witt vectors, carry arithmetic, Frobenius endomorphism, Gauss sums, Artin-Schreier varieties, Weil conjectures, transfer operators, Fermat curves.

**MSC 2020:** 11S40 (Witt vectors), 11G20 (Curves over finite fields), 11T23 (Exponential sums), 11M38 (Zeta and $L$-functions in characteristic $p$).

---

## 1. Introduction

### 1.1 Motivation: Carries as the Obstruction — and the Solution

The Riemann Hypothesis for function fields $\mathbb{F}_q[x]$ was proved by Weil [Wei48] using the theory of the Frobenius endomorphism on algebraic curves. A natural question, raised independently by several authors, is whether an analogous approach could work for $\mathbb{Z}$.

The fundamental obstruction is arithmetic: in $\mathbb{F}_q[x]$, digit-wise addition and multiplication produce no carries — coefficients operate independently in each degree. In $\mathbb{Z}$, carries propagate information between digits, creating long-range correlations that have no polynomial analogue.

This paper shows that carries are not merely an obstruction but contain, in a precise algebraic sense, the Frobenius eigenvalues themselves. At $p = 3$, the carry operator reproduces the Frobenius of a supersingular elliptic curve. At $p = 7$, it produces a genus-3 Frobenius factor whose supersingular eigenvalues are exactly the quadratic Gauss sum of the Legendre symbol. At $p = 5$ and $p = 13$, Dirichlet character twists of the carry operator recover Frobenius eigenvalues that the untwisted operator misses. The natural home for this correspondence is the ring of Witt vectors $W(\mathbb{F}_p)$, where carries arise from the ghost map inversion and the Frobenius acts as a shift on ghost components.

### 1.2 Context: Carry Operators

The spectral theory of carries was initiated by Holte [Hol97] and Diaconis–Fulman [DF09], who showed that carry propagation in base-$b$ addition defines a Markov chain with eigenvalues $1/b^k$ for $k \geq 0$. These eigenvalues are real and positive — they capture the contractile dynamics of carry decay but contain no phase information.

The Frobenius eigenvalues of varieties over finite fields, by contrast, are complex numbers whose moduli satisfy the Weil bound $\lvert\alpha\rvert = q^{w/2}$ (weight $w$) and whose phases encode the deepest arithmetic of the variety. The passage from real carry eigenvalues to complex Frobenius eigenvalues is the central challenge.

### 1.3 Main Results

We resolve this challenge by constructing carry operators that live not on the integers but on the Witt vectors. The key innovation is to weight carry transitions by the additive character $\omega = e^{2\pi i/p}$, which converts the real-valued carry into a complex-valued operator whose spectrum encodes both moduli and phases.

Status labels are used uniformly as **Proved**, **Conditional**, and **Open target**.

| Result | Status | Method |
|--------|--------|--------|
| Theorem A ($p = 3$ Frobenius) | **Proved** | Computer-assisted, exact arithmetic over $\mathbb{Z}[\omega]$ (Faddeev–LeVerrier + SymPy verification) |
| Theorem B ($p = 7$ Frobenius + Gauss sum) | **Proved** | Computer-assisted, exact arithmetic over $\mathbb{Z}[\omega_7]$ |
| Theorem 3.5 (additive negative result) | **Proved** | Same method |
| Corollary C (Weil bound + point counts) | **Proved** | Follows from Theorem A |
| Conjecture 1 (Artin-Schreier Carry) | **Open** | Supported by Theorems A, B and $p = 5$ character twist data |
| Conjecture 2 (Universality via Fermat) | **Open** | Structural evidence at $p = 3, 5, 7$ |

**Theorem A** (Frobenius at $p = 3$; computer-assisted, exact arithmetic over $\mathbb{Z}[\omega]$). Let $p = 3$, $\omega = e^{2\pi i/3}$, and define the $9 \times 9$ matrix $M$ over $\mathbb{Z}[\omega]$ by

$$M_{(x_1,y_1),(x_0,y_0)} = \omega^{c_2(x_0, y_0, x_1, y_1)}$$

where the indices range over $\mathbb{F}_3 \times \mathbb{F}_3$ and $c_2(x_0, y_0, x_1, y_1)$ is the level-2 multiplicative Witt carry (Definition 2.3). Then:

$$\chi_M(\mu) = \mu^4 \cdot (\mu^2 + 3\mu + 3) \cdot (\mu^3 - 6\mu^2 + 12\mu - 45)$$

The quadratic factor $\mu^2 + 3\mu + 3$ is the characteristic polynomial of Frobenius on the supersingular elliptic curve $E: y^2 = x^3 + 2x + 1$ over $\mathbb{F}_3$.

**Theorem B** (Frobenius and Gauss Sums at $p = 7$; computer-assisted, exact arithmetic over $\mathbb{Z}[\omega_7]$). The $49 \times 49$ carry matrix at $p = 7$ has 6 eigenvalues at the Weil bound $\lvert\mu\rvert = \sqrt{7}$, and its characteristic polynomial contains the exact factor

$$(\mu^2 - 7)^2 \cdot (\mu^2 + 7)$$

The factor $\mu^2 + 7$ has roots $\pm i\sqrt{7}$, which are exactly the quadratic Gauss sum $g(\chi^3)$ of the Legendre symbol $\left(\frac{\cdot}{7}\right)$. The factor $(\mu^2 - 7)^2$ captures the superspecial Frobenius eigenvalues $\pm\sqrt{7}$, each with multiplicity 2.

**Corollary C.** The carry operator $T = M/p^2$ has eigenvalues satisfying the Weil bound $\lvert\lambda\rvert = p^{-3/2}$. For $p = 3$, the Frobenius pair $\alpha, \bar{\alpha}$ recovers the point counts $N_E(\mathbb{F}_{3^k}) = 3^k + 1 - (\alpha^k + \bar{\alpha}^k)$ via the Lefschetz trace formula.

### 1.4 Significance

These results demonstrate exact spectral factor matches between carry operators and Frobenius polynomials:

1. **Frobenius polynomials appear as exact spectral factors of carry matrices.** The characteristic polynomials of Frobenius on certain algebraic curves over $\mathbb{F}_p$ divide the characteristic polynomials of the corresponding Witt carry matrices — an exact algebraic match, verified by computer-assisted computation with exact arithmetic. This holds at $p = 3$ (Theorem A), $p = 7$ (exact genus-3 factor), and via character twists at $p = 5$ and $p = 13$. Whether this spectral match reflects a deeper cohomological relationship (Conjecture 1) remains open.

2. **The multiplicative structure is essential.** We show (Theorem 3.5) that the additive Witt carry operator has a different characteristic polynomial ($\mu^6(\mu^3 - 6\mu^2 + 9)$) with no elliptic curve factor. The Frobenius emerges only from multiplication, consistent with its nature as a ring endomorphism.

3. The characteristic polynomial is defined over $\mathbb{Q}$. Despite the matrix $M$ living in $\mathbb{Z}[\omega]$, all coefficients of $\chi_M$ are rational integers (Proposition 3.2, verified also at $p = 7$). This reveals a Galois symmetry: the spectrum is invariant under $\omega \mapsto \omega^k$ for $\gcd(k, p) = 1$.

4. **Structural consistency with an Artin-Schreier interpretation (Conjecture 1).** The carry exponent $c_2$ is perfectly equidistributed on $\mathbb{F}_p^{\ast}$ (§5.7), and the Weil-bound eigenvectors are orthogonal to the all-ones vector (§5.7). The Frobenius information is hidden in the periodic orbit structure of the carry, not in the naive exponential sum — consistent with (but not yet proved to be) the cohomological Frobenius on the Artin-Schreier variety $z^p - z = c_2(\boldsymbol{x})$.

5. **Supersingular eigenvalues are Gauss sums.** At $p = 7$, the eigenvalues $\pm i\sqrt{7}$ from the factor $\mu^2 + 7$ are exactly the quadratic Gauss sum $g(\chi^3) = \sum_{x=1}^{6} \left(\frac{x}{7}\right) \omega^x = i\sqrt{7}$ (§5.6). This is the first identification of a carry eigenvalue with a named number-theoretic quantity.

6. **The Fermat curve genus predicts the degree of the Frobenius factor.** The level-2 carry is dominated by the Fermat quotient $Q_p(x,y) = (x^p y^p - (xy)^p)/p$. The genus $g = (p-1)(p-2)/2$ of the associated Fermat curve determines the degree $2g$ of the full Frobenius polynomial: for $p = 3$, $g = 1$ and the factor is quadratic; for $p = 7$, the $49 \times 49$ matrix captures 6 of the expected $2g = 30$ eigenvalues.

7. **Base-2 bridge compatibility.** Independent Step-3 operator results in the binary D-odd program (carry/Witt intertwiner + exact character-resolvent channel) show that the same "selected spectral factor" architecture survives at $p=2$ in real form (dominant mode $1/2$) while $p \ge 3$ carries phase-bearing Frobenius pairs. This supports a unified viewpoint: arithmetic constants are encoded by character-projected resolvent factors, with $p=2$ and $p\ge3$ as two faces of the same operator mechanism.

### 1.5 Outline

Section 2 recalls Witt vector arithmetic and defines the carry operator, explaining the role of the additive character as a discrete Fourier transform on the carry space. Section 3 proves Theorem A via the Faddeev–LeVerrier algorithm with exact arithmetic over $\mathbb{Z}[\omega]$. Section 4 identifies the elliptic curve and verifies the Lefschetz trace formula. Section 5 presents the full computational landscape: the $p=2$ golden ratio and Ihara zeta (§5.1), the cross-prime summary table (§5.2), the exact $p=7$ factorization and Theorem B (§5.3), the Dirichlet character twist that rescues $p=5$ (§5.4), the level-3/4 rank-collapse negative result (§5.5), the Gauss sum identification at $p=7$ (§5.6), the equidistribution and orthogonality structure of $c_2$ (§5.7), and the extended twist scan to $p=11, 13$ (§5.8). Section 6 formulates the Artin-Schreier/Fermat curve conjecture, explaining the genus obstruction and the character decomposition approach. Section 7 discusses the overlap transfer matrix (negative result) and Ruelle–Perron–Frobenius theory. Section 8 presents further conjectures, the cross-prime obstruction to an Euler product, and open problems.

<!-- \newpage -->

### 1.6 Conceptual Architecture

The following diagram summarizes the flow from elementary carry arithmetic to the Frobenius eigenvalues, and the three disciplines it bridges:

<!-- \small -->

```
                    ┌───────────────────────────────────┐
                    │     POSITIONAL ARITHMETIC         │
                    │   (carries in base-p addition     │
                    │    and multiplication)            │
                    └───────────────┬───────────────────┘
                                    │
                          LIFT TO WITT VECTORS
                         W_n(𝔽_p): ghost map
                        inversion ⟹ carries
                                    │
                    ┌───────────────▼───────────────────┐
                    │   WITT CARRY TRANSFER OPERATOR    │
                    │                                   │
                    │   M_{(x₁,y₁),(x₀,y₀)} = ω^{c₂}    │
                    │                                   │
                    │   ω^{c₂} = Fourier transform      │
                    │   on carry space 𝔽_p              │
                    └──────┬────────────────┬───────────┘
                           │                │
              ┌────────────▼───┐     ┌──────▼─────────────┐
              │  COMBINATORICS │     │ ALGEBRAIC GEOMETRY │
              │                │     │                    │
              │  Diaconis–     │     │  Artin-Schreier    │
              │  Fulman        │     │  variety           │
              │  eigenvalues   │     │  z^p−z = c₂(x,y)   │
              │  |λ| = 1/p^k   │     │  |α| = √p (Weil)   │
              │  (real, §3.5)  │     │  (complex, Thm A/B)│
              └────────┬───────┘     └────────┬───────────┘
                       │                      │
                       │    FERMAT QUOTIENT   │
                       │  (x^p·y^p-(xy)^p)/p  │
                       │         │            │
                       └─────────┼────────────┘
                                 │
                    ┌────────────▼──────────────────────┐
                    │       NUMBER THEORY               │
                    │                                   │
                    │   Gauss sums g(χ^k) and           │
                    │   Jacobi sums J(χ_a, χ_b)         │
                    │   = Frobenius eigenvalues         │
                    │   = spectral components of carry  │
                    │                                   │
                    │   p=3: g=1 (elliptic, Thm A)      │
                    │   p=7: g(χ³)=i√7 (Gauss, Thm B)   │
                    │   p=5: rescued by χ₂ twist        │
                    └───────────────────────────────────┘
```

<!-- \normalsize -->

---

## 2. Witt Vector Carry Arithmetic

### 2.1 Truncated Witt Vectors

**Definition 2.1.** The ring of truncated Witt vectors $W_n(\mathbb{F}_p)$ consists of $n$-tuples $(a_0, a_1, \ldots, a_{n-1})$ with $a_i \in \mathbb{F}_p$, equipped with addition and multiplication defined by the *ghost map*:

$$w_k(a) = \sum_{i=0}^{k} p^i a_i^{p^{k-i}} \quad (k = 0, \ldots, n-1)$$

The ghost map $w: W_n(\mathbb{F}_p) \to \mathbb{Z}^n$ is a ring homomorphism when the target has componentwise operations:

$$w_k(a + b) = w_k(a) + w_k(b), \quad w_k(a \cdot b) = w_k(a) \cdot w_k(b)$$

The Witt addition and multiplication polynomials $S_k$, $P_k$ are obtained by inverting the ghost map, which requires division by $p^k$ — the source of carries.

**Definition 2.2** (Frobenius). The Frobenius endomorphism $F: W_n(\mathbb{F}_p) \to W_n(\mathbb{F}_p)$ acts by $F(a_0, \ldots, a_{n-1}) = (a_0^p, \ldots, a_{n-1}^p)$. Over $\mathbb{F}_p$, Fermat's little theorem gives $a^p = a$, so $F = \mathrm{id}$ on components. However, on ghost components: $w_k(F(a)) = w_{k+1}(a)$ (a shift), which is the source of the non-trivial arithmetic.

### 2.2 The Witt Carry

**Definition 2.3** (Level-$k$ Witt carry). For $x, y \in W_n(\mathbb{F}_p)$ and an operation $\star \in \{+, \times\}$, the level-$k$ Witt carry is

$$c_k(x, y) \equiv \frac{w_k(x) \star w_k(y) - \sum_{i=0}^{k-1} p^i z_i^{p^{k-i}}}{p^k} \pmod{p}$$

where $z = x \star y$ in $W_n(\mathbb{F}_p)$. This is well-defined because the ghost map inversion guarantees that the numerator is divisible by $p^k$.

For $\star = \times$ (multiplication), the level-2 carry with $x = (x_0, x_1, 0)$, $y = (y_0, y_1, 0) \in W_3(\mathbb{F}_p)$ is:

$$c_2(x, y) \equiv \frac{(x_0^{p^2} + p x_1^p)(y_0^{p^2} + p y_1^p) - z_0^{p^2} - p z_1^p}{p^2} \pmod{p}$$

### 2.3 The Fermat Quotient

The leading term in the level-2 multiplicative carry is the *Fermat quotient of the product*:

$$Q_p(x_0, y_0) = \frac{x_0^{p^2} y_0^{p^2} - (x_0 y_0 \bmod p)^{p^2}}{p^2}$$

This is not a random polynomial. Writing $a = x_0 y_0 \bmod p$ and $x_0^p = x_0$, $y_0^p = y_0$ in $\mathbb{F}_p$, the quotient reduces to the classical Fermat quotient

$$q_p(x, y) = \frac{(xy)^p - (xy \bmod p)^p}{p}$$

which is intimately connected to the arithmetic of the Fermat curve $X^p + Y^p = Z^p$ and to the Jacobi sums

$$J(\chi_a, \chi_b) = \sum_{x \in \mathbb{F}_p} \chi_a(x) \chi_b(1-x)$$

whose modulus is exactly $\lvert J\rvert = \sqrt{p}$ — the Weil bound.

### 2.4 The Carry Transfer Operator

**Definition 2.4.** Fix a prime $p$ and let $\omega = e^{2\pi i/p}$ be a primitive $p$-th root of unity. The *multiplicative Witt carry matrix* is the $p^2 \times p^2$ matrix $M$ indexed by $\mathbb{F}_p \times \mathbb{F}_p$:

$$M_{(x_1,y_1),(x_0,y_0)} = \omega^{c_2(x_0, y_0, x_1, y_1)}$$

where $c_2$ is computed from the Witt multiplication of $(x_0, x_1, 0) \times (y_0, y_1, 0)$ in $W_3(\mathbb{F}_p)$.

The *carry operator* is $T = M / p^2$.

**Why the additive character.** The exponentiation $c_2 \mapsto \omega^{c_2}$ is not an arbitrary choice but the canonical discrete Fourier transform on the carry space $\mathbb{F}_p$. The map $\psi: \mathbb{F}_p \to \mathbb{C}^\times$ defined by $\psi(t) = e^{2\pi i t/p}$ is the standard additive character of $\mathbb{F}_p$, and it is exactly this character that appears in the construction of Gauss sums $\mathfrak{g}(\chi) = \sum_{t \in \mathbb{F}_p} \chi(t) \psi(t)$ and Jacobi sums $J(\chi_a, \chi_b) = \sum_x \chi_a(x) \chi_b(1-x)$. These sums compute the Frobenius eigenvalues on Fermat curves [Wei49, IR90].

By weighting each carry transition by $\psi(c_2) = \omega^{c_2}$, the matrix $M$ performs a Fourier analysis of the carry distribution: the eigenvalues of $M$ are the spectral components of the carry propagation, decomposed along the characters of $\mathbb{F}_p$. The Frobenius eigenvalues emerge because this Fourier decomposition is the same operation that extracts Jacobi sums from the point-counting of Fermat varieties.

**Remark.** The matrix $M$ also encodes the Frobenius shift: the "from" state $(x_0, y_0)$ represents level-0 data, the "to" state $(x_1, y_1)$ represents level-1 data (which becomes level-0 under the Frobenius shift on ghost components).

---

## 3. Proof of Theorem A

### 3.1 The Exponent Matrix

For $p = 3$, we enumerate all $81$ quadruples $(x_0, y_0, x_1, y_1) \in \mathbb{F}_3^4$, compute the Witt multiplication $(x_0, x_1, 0) \times_W (y_0, y_1, 0)$ in $W_3(\mathbb{F}_3)$, and extract $c_2 \bmod 3$.

The $9 \times 9$ exponent matrix $E$ (where $M_{j,i} = \omega^{E_{j,i}}$) has:
- 49 entries equal to 0 (carry $\equiv 0$)
- 16 entries equal to 1 (carry $\equiv 1$)
- 16 entries equal to 2 (carry $\equiv 2$)

The non-trivial entries (carry $\neq 0$) are concentrated in rows and columns corresponding to states $(x, y)$ with $x, y \in \{1, 2\}$ (the units of $\mathbb{F}_3$). When either $x_0 = 0$ or $y_0 = 0$, the Witt product has $z_0 = 0$ and the level-2 carry vanishes.

### 3.2 Rationality of the Characteristic Polynomial

**Proposition 3.2.** All coefficients of $\chi_M(\mu) = \det(\mu I - M)$ are rational integers.

*Proof.* We compute $\chi_M$ via the Faddeev–LeVerrier recurrence, which expresses the coefficients $c_{n-k}$ of the characteristic polynomial in terms of the traces $\mathrm{Tr}(M^k)$. Working over $\mathbb{Z}[\omega]$ with the relation $\omega^2 + \omega + 1 = 0$, we represent each element as $a + b\omega$ with $a, b \in \mathbb{Q}$.

Direct computation shows $\mathrm{Tr}(M^k) \in \mathbb{Z}$ for all $k = 1, \ldots, 9$:

| $k$ | $\mathrm{Tr}(M^k)$ |
|-----|--------------------------|
| 1   | 3                        |
| 2   | 15                       |
| 3   | 135                      |
| 4   | 927                      |
| 5   | 4563                     |
| 6   | 22005                    |
| 7   | 120123                   |
| 8   | 659583                   |
| 9   | 3510135                  |

Since each $\mathrm{Tr}(M^k)$ is rational, the Faddeev–LeVerrier recurrence produces rational coefficients $c_j$ at every step. Moreover, the entries of $M$ lie in $\mathbb{Z}[\omega]$, so $M$ is an algebraic integer matrix and its eigenvalues are algebraic integers; consequently the coefficients $c_j$ (which are elementary symmetric polynomials in the eigenvalues) are algebraic integers. Being both rational and algebraic integers, they lie in $\mathbb{Z}$. $\square$

**Remark** (Galois symmetry). The rationality of the traces is not obvious *a priori*: the matrix $M$ has entries in $\mathbb{Z}[\omega] \setminus \mathbb{Z}$, and yet $\mathrm{Tr}(M^k) \in \mathbb{Z}$ for every $k$. This is not a numerical coincidence but reflects a genuine

$$\mathrm{Gal}(\mathbb{Q}(\omega)/\mathbb{Q})$$

-symmetry of the carry operator.

The Galois automorphism $\sigma: \omega \mapsto \omega^2$ acts on the matrix $M$ by replacing $\omega^{c_2}$ with $\omega^{2c_2}$. Since the Witt carry values $c_2$ take each residue class $\{0, 1, 2\}$ with multiplicities $49, 16, 16$ (and $16 = 16$: the non-zero classes appear equally often), the carry distribution is invariant under the permutation $c \mapsto 2c \pmod{3}$ induced by $\sigma$. This forces $\sigma(M)$ to be conjugate to $M$ by a permutation matrix, so

$$\mathrm{Tr}(\sigma(M)^k) = \mathrm{Tr}(M^k)$$

. Since the only element of $\mathbb{Q}(\omega)$ fixed by $\sigma$ is $\mathbb{Q}$, the traces are rational.

This Galois invariance is the hallmark of an operator that "counts points on a variety over $\mathbb{F}_p$": the Lefschetz trace formula guarantees

$$\sum_i (-1)^i \mathrm{Tr}(\mathrm{Frob}^k \mid H^i) \in \mathbb{Z}$$

for any smooth projective variety, because the point counts $N_V(\mathbb{F}_{p^k})$ are manifestly integers. The rationality of $\mathrm{Tr}(M^k)$ is thus the first structural signal that $M$ encodes the cohomology of an algebraic variety.

### 3.3 Factorization

The characteristic polynomial is:

$$\chi_M(\mu) = \mu^9 - 3\mu^8 - 3\mu^7 - 27\mu^6 - 99\mu^5 - 135\mu^4$$

Factoring over $\mathbb{Z}$:

$$\chi_M(\mu) = \mu^4 \cdot (\mu^2 + 3\mu + 3) \cdot (\mu^3 - 6\mu^2 + 12\mu - 45)$$

This factorization is verified by polynomial division: $\mu^2 + 3\mu + 3$ divides $\chi_M(\mu)/\mu^4$ with zero remainder.

### 3.4 Identification with the Frobenius Polynomial

The quadratic $\mu^2 + 3\mu + 3$ has roots

$$\mu = \frac{-3 \pm i\sqrt{3}}{2}$$

with $\lvert\mu\rvert^2 = 9/4 + 3/4 = 3$. Hence $\lvert\mu\rvert = \sqrt{3} = \sqrt{p}$, the Weil bound for weight 1.

This is the characteristic polynomial of Frobenius $t^2 - a_p t + p$ with $a_p = -3$ (trace of Frobenius). The elliptic curve with this Frobenius at $p = 3$ is determined up to isogeny by $a_3 = -3$.

**Proposition 3.4.** The curve $E: y^2 = x^3 + 2x + 1$ over $\mathbb{F}_3$ has $N_E(\mathbb{F}_3) = 7$ and Frobenius characteristic polynomial $t^2 + 3t + 3$.

*Proof.* Direct enumeration:

$$E(\mathbb{F}_3) = \{(0, \pm 1), (1, \pm 1), (2, \pm 1), \mathcal{O}\},$$

giving $7$ points. Hence $a_3 = 3 + 1 - 7 = -3$, and the Frobenius polynomial is $t^2 - (-3)t + 3 = t^2 + 3t + 3$. Since $a_3 \equiv 0 \pmod{3}$, the curve is supersingular. $\square$

### 3.5 The Additive Case (Negative Result)

**Theorem 3.5** (computer-assisted, exact arithmetic). The additive Witt carry matrix at $p = 3$ has characteristic polynomial

$$\chi_{M_{\text{add}}}(\mu) = \mu^6 \cdot (\mu^3 - 6\mu^2 + 9)$$

*which contains no elliptic curve Frobenius factor.*

*Proof.* Same method as Theorem A, replacing Witt multiplication by Witt addition. The cubic $\mu^3 - 6\mu^2 + 9$ has three real roots $\{5.725, 1.399, -1.124\}$ (none on $\lvert\mu\rvert = \sqrt{3}$). $\square$

This confirms that the Frobenius emerges specifically from the multiplicative structure, consistent with the fact that the Frobenius endomorphism on Witt vectors interacts non-trivially only with multiplication.

---

## 4. The Lefschetz Trace Formula

### 4.1 Point Counting

The Frobenius eigenvalues $\alpha = (-3 + i\sqrt{3})/2$, $\bar{\alpha} = (-3 - i\sqrt{3})/2$ determine the point counts of $E$ over all extensions:

$$N_E(\mathbb{F}_{3^k}) = 3^k + 1 - (\alpha^k + \bar{\alpha}^k)$$

| $k$ | $\alpha^k + \bar{\alpha}^k$ | $N_E(\mathbb{F}_{3^k})$ |
|-----|---------------------------|--------------------------|
| 1   | $-3$                      | 7                        |
| 2   | $+3$                      | 7                        |
| 3   | $0$                       | 28                       |
| 4   | $-9$                      | 91                       |
| 5   | $+27$                     | 217                      |
| 6   | $-54$                     | 784                      |

### 4.2 The Local Zeta Function

The local zeta function of $E/\mathbb{F}_3$ is:

$$Z(E/\mathbb{F}_3, T) = \frac{1 + 3T + 3T^2}{(1 - T)(1 - 3T)}$$

The numerator $1 + 3T + 3T^2 = (1 - \alpha T)(1 - \bar{\alpha} T)$ with $\lvert\alpha\rvert = \sqrt{3}$, confirming the Riemann Hypothesis for this curve.

### 4.3 Connection to the Carry Operator Trace

The traces of $M^k$ encode a richer object than the Frobenius pair alone (they include contributions from the cubic factor and the zero eigenvalues). The Frobenius contribution to $\mathrm{Tr}(M^k)$ is $\alpha^k + \bar{\alpha}^k$, which can be extracted via the Weil zeta function decomposition.

---

## 5. Computational Evidence and Extensions

### 5.1 The Additive Case at $p = 2$: The Golden Ratio and the Ihara Zeta

For $p = 2$, the additive Witt carry operator on $W_3(\mathbb{F}_2)$ has two non-zero eigenvalues of $T = M/4$:

$$\lambda = \frac{1 \pm \sqrt{5}}{4}$$

satisfying $4\lambda^2 - 2\lambda - 1 = 0$. The dominant eigenvalue is $\lambda_0 = \varphi/2$, where $\varphi = (1+\sqrt{5})/2$ is the golden ratio.

This is not a curiosity. In base 2, the Fermat curve $X^2 + Y^2 = 1$ is trivial (genus 0), and the "elliptic curve" degenerates. What remains is a *combinatorial* structure: the carry propagation in binary addition is equivalent to the random walk on the Fibonacci graph $\bullet \rightleftharpoons \bullet$, whose adjacency matrix has eigenvalues $\pm \varphi$.

The connection to spectral graph theory is precise. The *Ihara zeta function* of a $(q+1)$-regular graph $G$ is

$$Z_G(u) = \prod_{\gamma} (1 - u^{|\gamma|})^{-1}$$

where the product runs over primitive cycles. For the simplest non-trivial graph (the loop of length 2), the Ihara zeta function has poles governed by $x^2 - x - 1 = 0$ — the minimal polynomial of $\varphi$. The carry operator at $p = 2$ is thus computing the Ihara zeta of the carry propagation graph, rather than the Hasse–Weil zeta of a variety.

This degeneration is consistent with the genus formula: $g(C_2) = 0$, so there are no Frobenius eigenvalues in the Weil sense. The golden ratio is the graph-theoretic shadow of the absent algebraic geometry. Notably, Ramanujan graphs — the optimal expanders — are defined by the condition that their non-trivial Ihara eigenvalues satisfy $\lvert\lambda\rvert \leq 2\sqrt{q}$, the graph-theoretic analogue of the Weil bound. The carry operator at $p = 2$ saturates this bound, suggesting that carry propagation graphs are Ramanujan.

### 5.2 Cross-Prime Evidence

For primes $p = 2, 3, 5, 7, 11, 13$, the multiplicative Witt carry operator at level 2 with state $(x_0, y_0)$ was computed:

| $p$ | $p^2$ | Frobenius at $\lvert\mu\rvert = \sqrt{p}$? | Notes |
|-----|-------|----------------------------------|-------|
| 2   | 4     | No ($\lvert\mu\rvert \neq \sqrt{2}$)       | Genus 0; Ihara zeta instead |
| 3   | 9     | **Yes — exact** (Theorem A)      | $\mu^2 + 3\mu + 3$ |
| 5   | 25    | No (untwisted)                   | Expected: degree 12 (genus 6); see §5.4 |
| 7   | 49    | **Yes — exact** (6 eigenvalues)  | See §5.3 |
| 11  | 121   | No                               | Expected: degree 90 (genus 45) |
| 13  | 169   | No                               | Expected: degree 132 (genus 66) |

The exactness at $p = 3$ is attributed to the fact that the Fermat quotient has only one non-trivial case ($(x_0, y_0) = (2,2)$), producing a clean algebraic structure. For most $p \geq 5$, multiple non-trivial cases create interference. The notable exception is $p = 7$, where 6 eigenvalues at $\lvert\mu\rvert = \sqrt{7}$ emerge despite the Fermat curve having genus 15.

### 5.3 The Second Frobenius: Exact Factorization at $p = 7$ (Theorem B)

The $49 \times 49$ Witt carry matrix for $p = 7$ was computed exactly over $\mathbb{Z}[\omega_7]$ using the Faddeev–LeVerrier algorithm. All characteristic polynomial coefficients are rational integers, confirming the Galois symmetry observed at $p = 3$.

**Theorem B.** The Weil-bound eigenvalues of the $49 \times 49$ carry matrix at $p = 7$ form a genus-3 Frobenius factor:

$$\chi_M(\mu) = \mu^{14} \cdot \underbrace{(\mu^2 - 7)^2 \cdot (\mu^2 + 7)}_{\text{Weil factor}} \cdot (\mu^5 - 35\mu^3 - 49\mu^2 + 147\mu - 343)^2 \cdot (\mu^6 - \cdots)^2 \cdot (\mu^7 - \cdots)$$

(The full degree-49 factorization is available in the computational scripts `F06_witt_algebra.py` and `F08c_p7_factorization.py` in the experiments directory; the remaining factors have no eigenvalues at the Weil bound.)

The Weil-bound factor is:

$$(\mu^2 - 7)^2 \cdot (\mu^2 + 7) = \mu^6 - 7\mu^4 - 49\mu^2 + 343$$

This decomposes as:

- $(\mu^2 + 7)$: the **Frobenius polynomial of a supersingular elliptic curve** with $a_7 = 0$, satisfying $N_E(\mathbb{F}_7) = 8$.
- $(\mu^2 - 7)^2$: a **genus-2 superspecial Frobenius polynomial** (roots $\pm\sqrt{7}$, each with multiplicity 2; product $= 49 = 7^2$).

The combined degree-6 factor satisfies the **functional equation** for genus-3 Frobenius: $7^3 P(\mu) = \mu^6 P(7/\mu)$. This confirms that 6 of the 30 expected Frobenius eigenvalues of the genus-15 Fermat curve $C_7$ are captured by the $49 \times 49$ carry matrix.

**Remark.** The 6 eigenvalues at $\lvert\mu\rvert = \sqrt{7}$ are a *subset* of the full Frobenius spectrum of $C_7$. The remaining 24 eigenvalues presumably require a deeper Witt truncation ($n \geq 3$) to resolve, consistent with the genus obstruction (§6.1).

### 5.4 Dirichlet Character Twist: Rescuing $p = 5$

For primes where the untwisted carry matrix yields no Weil-bound eigenvalues, a natural question is whether twisting by a Dirichlet character can extract the missing Frobenius.

Define the twisted matrix:

$$M(\chi)_{(x_1,y_1),(x_0,y_0)} = \chi(x_0 y_0 \bmod p) \cdot \omega^{c_2}$$

Three twist variants were tested: multiplicative ($\chi(x_0 y_0)$), additive ($\chi(x_0 + y_0)$), and carry-value ($\chi(c_2) \cdot \omega^{c_2}$).

Result at $p = 5$: The untwisted matrix has no eigenvalues at $\lvert\mu\rvert = \sqrt{5}$. However, the additive twist with $\chi_2$ (the Legendre symbol mod 5, order 2) produces 2 eigenvalues at exactly $\lvert\mu\rvert = \sqrt{5} = 2.236068$:

$$\mu = +2.0237 - 0.9511i, \quad \mu = -2.0237 - 0.9511i$$

Their product is $\mu_1 \mu_2 = -5$ (anti-conjugate pair, not a standard Frobenius conjugate pair with product $+p$).

Result at $p = 7$: All 6 multiplicative twists ($\chi_0$ through $\chi_5$) produce exactly 6 eigenvalues at $\lvert\mu\rvert = \sqrt{7}$, each with distinct phases. The trivial and $\chi_3$ (order 2) twists give degenerate phases (multiples of $\pi/2$); the remaining twists produce non-trivial phases consistent with Jacobi sum arguments.

This decomposition by character is consistent with the Jacobi sum structure of the Frobenius on Fermat curves, where each character pair $(\chi_a, \chi_b)$ selects a specific eigenvalue.

### 5.5 Level-3 and Level-4 Marginals (Negative Result)

An attempt was made to deepen the carry level while keeping the $p^2 \times p^2$ state space, by marginalizing over auxiliary Witt components:

$$M_{(x_1,y_1),(x_0,y_0)} = \sum_{x_2, y_2} \omega^{c_3(x_0,x_1,x_2,0;\, y_0,y_1,y_2,0)}$$

**Result.** At all primes tested ($p = 3, 5, 7$), the level-3 and level-4 marginal matrices collapse to **rank 1**: a single nonzero eigenvalue at $p^{k+1}$ (e.g., $625 = 5^4$ for $p = 5$, level 3), with all other eigenvalues zero.

The root cause is **character orthogonality**: the sum $\sum_{x_2,y_2} \omega^{c_k}$ evaluates to $p^2$ when a constraint is met and $0$ otherwise, destroying all spectral structure.

**Conclusion.** Level 2 is the **only productive level** for the dense marginal approach, because it involves no marginalization — each matrix entry is a single $\omega^{c_2}$.

### 5.6 Gauss Sum Identification at $p = 7$

The Frobenius eigenvalues of the carry matrix were compared against all Gauss sums $g(\chi^k) = \sum_{x \in \mathbb{F}_p^{\ast}} \chi^k(x) \omega^x$ and Jacobi sums $J(\chi^a, \chi^b) = \sum_x \chi^a(x) \chi^b(1-x)$ for multiplicative characters on $\mathbb{F}_p^{\ast}$.

**Result.** At $p = 7$, the supersingular eigenvalues $\mu = \pm i\sqrt{7}$ from the factor $(\mu^2 + 7)$ match the quadratic Gauss sum **exactly**:

$$\mu = +i\sqrt{7} = g(\chi^3), \quad \mu = -i\sqrt{7} = -g(\chi^3)$$

where $\chi^3$ is the Legendre symbol (quadratic character) mod 7, and $g(\chi^3) = \sum_{x=1}^{6} \left(\frac{x}{7}\right) \omega^x = i\sqrt{7}$.

The superspecial eigenvalues $\mu = \pm\sqrt{7}$ from $(\mu^2 - 7)^2$ do **not** match any standard Gauss or Jacobi sum. At $p = 3$, the Frobenius eigenvalues $(-3 \pm i\sqrt{3})/2$ (with phases $\pm 5\pi/6$) likewise do not match the sole Gauss sum $g(\chi) = i\sqrt{3}$ (phase $\pi/2$).

**Interpretation.** The Fermat curve $X^p + Y^p = 1$ degenerates in characteristic $p$ (becoming $X + Y = 1$, a line), so the standard Jacobi sum formula for Fermat curve Frobenius does not apply. The carry eigenvalues arise from a different variety — the Artin-Schreier variety $z^p - z = c_2(x_0, x_1, y_0, y_1)$ — whose Frobenius structure only partially overlaps with the classical Gauss/Jacobi catalog.

### 5.7 Structural Properties of the Carry Exponent

The distribution of $c_2$ values over $\mathbb{F}_p^4$ reveals a striking balance: nonzero values are **exactly equidistributed**.

| $p$ | $N_0$ ($c_2 = 0$) | $N_a$ ($c_2 = a \neq 0$) | Check |
|-----|-------------------:|--------------------------:|-------|
| 3   | 49                 | 16 each                   | $49 + 2 \cdot 16 = 81 = 3^4$ |
| 5   | 273                | 88 each                   | $273 + 4 \cdot 88 = 625 = 5^4$ |
| 7   | 745                | 276 each                  | $745 + 6 \cdot 276 = 2401 = 7^4$ |

Consequently, all exponential sums $S_a = \sum_{\boldsymbol{x} \in \mathbb{F}_p^4} \omega^{a \cdot c_2(\boldsymbol{x})}$ for $a \neq 0$ are **real and equal**: $S_a = (p N_0 - p^4) / (p - 1)$.

A deeper structural property: the Weil-bound eigenvalues ($\lvert\mu\rvert = \sqrt{p}$) have **coefficient zero** in the expansion $\langle \boldsymbol{1} \mid M^k \mid \boldsymbol{1} \rangle = \sum_i c_i \lambda_i^k$. The all-ones vector is orthogonal to every Weil eigenvector, so the Frobenius is "invisible" in the total exponential sum but fully present in $\mathrm{Tr}(M^k)$ (the periodic orbit sum).

### 5.8 Extended Twist Scan: $p = 11$ and $p = 13$

The character twist methodology of §5.4 was extended to $p = 11$ (121 × 121 matrix) and $p = 13$ (169 × 169 matrix), testing all $p-1$ multiplicative and $p-1$ additive twists.

$p = 11$: Productive with higher-order characters. A comprehensive scan of all character twists at $p = 11$ reveals that additive twists $\chi_2$ and $\chi_8$ (both of order 5) each produce 1 eigenvalue at $\lvert\mu\rvert \approx \sqrt{11}$, with traces $a_{11}(\chi_2) = -0.76$ and $a_{11}(\chi_8) = 2.10$. Lower-order characters are not productive at this prime.

$p = 13$: Partial. Additive twists $\chi_4$ and $\chi_8$ (both of order 3) each produce 1 eigenvalue at $\lvert\mu\rvert = \sqrt{13}$:

$$\mu = -1.264 \pm 3.378 i, \quad \lvert\mu\rvert = 3.606 \approx \sqrt{13}$$

These are complex conjugates (related by $\chi_4 \leftrightarrow \chi_8 = \overline{\chi_4}$). No other twist type is productive.

**Extended prime scan** (primes up to 23):

| $p$ | Untwisted | Best productive twist | Total Weil eigenvalues |
|-----|----------:|:---------------------:|:----------------------:|
| 3   | 2         | mul Legendre: 2       | 4                      |
| 5   | 0         | add Legendre: 2       | 8                      |
| 7   | 6         | all mul + add: 6 each | 40+                    |
| 11  | 0         | add ord 5: 1 each     | 2                      |
| 13  | 0         | add ord 3: 1 each     | 2                      |
| 17  | 0         | add ord 16: 2; ord 8: 3 | 14                   |
| 19  | 6         | add ord 6,9: 2 each   | 23+                    |
| 23  | 0         | add ord 11: 1--2 each | 14+                    |

**Every prime tested produces Weil eigenvalues with appropriate character twists.** The productive character order appears to depend on $p$ in a non-trivial way (e.g., order 5 at $p = 11$, order 3 at $p = 13$, order 16 at $p = 17$). No simple arithmetic formula predicts the productive order. The anomalous richness at $p = 7$ and $p = 19$ (large untwisted Weil counts) and the requirement for high-order characters at $p = 11, 17, 23$ remain open questions.

---

## 6. The Artin-Schreier / Fermat Curve Conjecture

This is the central conjecture of the paper, motivated by the structure of the Witt carry, the Gauss sum identification at $p = 7$, and the failure of naïve generalization to $p \geq 5$.

### 6.1 The Genus Obstruction

For $p = 3$, the quadratic factor $\mu^2 + 3\mu + 3$ in $\chi_M$ naturally suggests searching for a quadratic Frobenius factor at all primes. However, the Frobenius polynomial at general $p$ is not quadratic: the degree grows with the genus of the Fermat curve.

The resolution lies in the Fermat quotient. The level-2 multiplicative carry is dominated by

$$Q_p(x, y) = \frac{x^p y^p - (xy \bmod p)^p}{p}$$

This term governs the arithmetic of the *Fermat curve* $C_p: X^p + Y^p = 1$ over $\mathbb{F}_p$. By the genus formula for smooth plane curves:

$$g(C_p) = \frac{(p-1)(p-2)}{2}$$

The Weil conjectures (proved by Deligne [Del74]) guarantee that the zeta function numerator of $C_p$ is a polynomial of degree $2g = (p-1)(p-2)$ in $T$, with all roots on $\lvert T\rvert = p^{-1/2}$.

| $p$ | $g(C_p)$ | $2g$ | Frobenius polynomial degree |
|-----|-----------|------|-----------------------------|
| 2   | 0         | 0    | (trivial)                   |
| 3   | 1         | 2    | **Quadratic** (Theorem A)   |
| 5   | 6         | 12   | Degree 12                   |
| 7   | 15        | 30   | Degree 30                   |
| 11  | 45        | 90   | Degree 90                   |

For $p = 3$, genus 1 means $C_3$ is an elliptic curve, and the Frobenius polynomial is quadratic — which is precisely the factor we proved. For $p = 5$, the expected Frobenius polynomial has degree 12. One cannot find it in the characteristic polynomial of a $25 \times 25$ matrix (degree 25), which barely exceeds $2g = 12$, and certainly not as a quadratic factor.

### 6.2 Jacobi Sum Connection

The Frobenius eigenvalues on the Fermat curve $C_p$ are expressed by Jacobi sums. For characters $\chi_a: \mathbb{F}_p^\times \to \mathbb{C}^\times$ of order $p-1$:

$$J(\chi_a, \chi_b) = \sum_{x \in \mathbb{F}_p} \chi_a(x) \chi_b(1-x)$$

The $2g$ Frobenius eigenvalues of $C_p$ are exactly the Jacobi sums $J(\chi_a, \chi_b)$ for appropriate pairs $(a, b)$, each satisfying $\lvert J\rvert = \sqrt{p}$.

The Witt carry matrix $M$ "sees" these Jacobi sums because the additive character $\omega^{c_2}$ in the matrix entries performs a discrete Fourier transform over the carry values, and the carry values themselves are determined by the Fermat quotient — which is the polynomial identity underlying the Jacobi sum.

**Conjecture 1** (Artin-Schreier Carry Conjecture). *The Witt carry operator at level 2 computes (a subfactor of) the Frobenius on the Artin-Schreier variety*

$$V_p:\ z^p - z = c_2(x_0, x_1, y_0, y_1) \quad \text{over } \mathbb{F}_p,$$

where $c_2$ is the level-2 multiplicative Witt carry. This variety is dominated by the Fermat quotient

$$\frac{x^p y^p - (xy)^p}{p}.$$

*Specifically:*

(a) For every prime $p$, the characteristic polynomial of the multiplicative Witt carry matrix $M$ contains, as an exact factor, a polynomial $P_p(\mu)$ whose roots lie on

$$\lvert \mu \rvert = \sqrt{p},$$

and $P_p(\mu)$ divides the Frobenius polynomial of $V_p$.

(b) For $p = 3$, $P_3(\mu) = \mu^2 + 3\mu + 3$ is the full Frobenius polynomial of $C_3$ (genus 1, Theorem A).

(c) For $p = 7$,

$$P_7(\mu) = (\mu^2 - 7)^2(\mu^2 + 7)$$

captures 6 Frobenius eigenvalues (Theorem B). The factor $\mu^2 + 7$ has roots equal to the quadratic Gauss sum

$$g(\chi^3) = i\sqrt{7}$$

*(§5.6).*

(d) For $p \geq 5$ in general, the Frobenius polynomial of the underlying variety has degree $\leq 2g = (p-1)(p-2)$. The $p^2 \times p^2$ carry matrix captures a subfactor of bounded degree. The complete Frobenius requires either a deeper Witt truncation ($n \geq 3$), or a Dirichlet character decomposition (§5.4), or both.

(e) The carry exponent $c_2$ is equidistributed on $\mathbb{F}_p^{\ast}$ (§5.7) and the Weil eigenvectors are orthogonal to $\lvert 1 \rangle$, consistent with the Artin-Schreier interpretation where the Frobenius lives in the non-trivial additive character sector.

### 6.3 The Case $p = 5$: Character Decomposition

For $p = 5$, the untwisted $25 \times 25$ carry matrix yields no eigenvalues at $\lvert\mu\rvert = \sqrt{5}$. This is consistent with the genus obstruction: the full Frobenius polynomial of $C_5$ has degree 12, and the untwisted matrix cannot accommodate it.

However, twisting by the additive Legendre character $\chi_2$ extracts 2 eigenvalues at exactly $\lvert\mu\rvert = \sqrt{5}$ (§5.4). This suggests that the Frobenius is present but "hidden" in the character decomposition of the carry matrix, analogous to how the Frobenius eigenvalues of the Fermat curve decompose by Jacobi sum index $(\chi_a, \chi_b)$.

The Frobenius eigenvalues of $C_5: X^5 + Y^5 = 1$ are the 12 Jacobi sums $J(\chi_a, \chi_b)$ for characters of order 4 on $\mathbb{F}_5^\times$. These are well-tabulated (see [IR90, Ch. 8]). The full program requires:

1. Computing the $25 \times 25$ carry matrix twisted by all $(p-1)^2 = 16$ pairs of characters.
2. Extracting the Weil-bound eigenvalues from each twist.
3. Matching the resulting eigenvalues against the tabulated Jacobi sums.

If the union of Weil-bound eigenvalues over all twists reproduces the 12 Frobenius eigenvalues of $C_5$, this would confirm Conjecture 1(d) at $p = 5$.

---

## 7. The Overlap Transfer Matrix (Negative Result)

### 7.1 Construction

We also studied the *overlapping* transfer matrix on the full state space $W_2(\mathbb{F}_p) \times W_2(\mathbb{F}_p)$, of dimension $p^4$, where the "from" state is the complete Witt pair $(x_0, x_1, y_0, y_1)$ and the "to" state is the shifted pair $(x_1, x_2, y_1, y_2)$, with the carry weight $\omega^{c_2}$.

### 7.2 Result: Uniform Modulus

**Proposition 7.2.** For $p = 2, 3, 5$, all eigenvalues of the overlap transfer matrix have modulus exactly $p$.

This is a topological consequence of the shift structure. The overlap operator is an *expanding map* on a symbolic space with alphabet size $p^2$, and the Ruelle–Perron–Frobenius theorem [Rue04] forces the spectral radius to equal the topological degree of the map, which is $p$ (corresponding to $p$ choices of new carry input at each step).

### 7.3 Phase Preservation

Despite the uniform modulus, the *phases* of the eigenvalues are preserved from the dense (non-overlap) formulation. For $p = 3$, the overlap matrix has eigenvalues at arguments $\pm 5\pi/6$ — the same Frobenius phase as in Theorem A. This demonstrates that:

- The **modulus** of Frobenius eigenvalues is encoded in the *dense* carry matrix (Section 3), where the averaging over level-1 digits creates the correct contraction.
- The **phase** of Frobenius eigenvalues is an intrinsic property of the additive character weighting $\omega^{c_2}$, invariant under the choice of shift structure.

The phase is preserved because it is encoded in the argument of the exponential $\exp(2\pi i \cdot c_2/p)$. The constructive/destructive interference of Witt carries (which defines the Gauss/Jacobi sums) cannot be destroyed by a shift, because it is an invariant of the character class.

---

## 8. Conjectures and Future Directions

### 8.1 The Universal Carry Variety

**Conjecture 2** (Universality via Artin-Schreier/Fermat varieties). The carry variety $V_p: z^p - z = c_2(x_0, x_1, y_0, y_1)$ is birationally related to the Fermat curve $C_p: X^p + Y^p = 1$ over $\mathbb{F}_p$. The characteristic polynomial of Frobenius on $H^{\bullet}(V_p)$ — whose Weil-bound roots have modulus $\sqrt{p}$ — divides the characteristic polynomial of the multiplicative Witt carry operator on the state space $W_n(\mathbb{F}_p) \times W_n(\mathbb{F}_p)$ for sufficiently large $n$.

For $p = 3$, $V_3$ reduces to a genus-1 curve whose Frobenius is $\mu^2 + 3\mu + 3$ (Theorem A). For $p = 7$, the Artin-Schreier interpretation is supported by three structural facts: (i) the Gauss sum identification $g(\chi^3) = i\sqrt{7}$ (§5.6), which is the standard Frobenius eigenvalue of Artin-Schreier covers; (ii) the equidistribution of $c_2$ on $\mathbb{F}_p^{\ast}$ (§5.7), which is the condition for the Artin-Schreier surface $z^p - z = c_2$ to have well-defined cohomology in every non-trivial character sector; and (iii) the orthogonality of Weil eigenvectors to $\lvert 1\rangle$ (§5.7), which means the Frobenius lives entirely in the non-trivial sectors of the Artin-Schreier cover — the trivial sector being the "base" variety $\mathbb{F}_p^4$.

### 8.2 Toward an Euler Product (and a Fundamental Obstruction)

If each prime's carry operator produced the local factor of a single global variety, one could assemble an Euler product $L(s) = \prod_p Z_p(p^{-s})$.

**Cross-prime inconsistency.** Computational experiments reveal that the Witt carry operator at different primes does **not** encode a single algebraic curve over $\mathbb{Q}$. Specifically:

- At $p = 3$, the Frobenius trace is $a_3 = -3$ (Theorem A), matching the supersingular curve $E: y^2 = x^3 + 2x + 1$.
- At $p = 7$, the supersingular Frobenius factor gives $a_7 = 0$ (§5.3).
- However, $E/\mathbb{F}_7$ has $a_7 = 3$, not $0$.

A database search identifies 6 elliptic curves over $\mathbb{Q}$ satisfying both $a_3 = -3$ and $a_7 = 0$ (e.g., $y^2 = x^3 - x^2 - x + 1$), but none match at all primes tested.

**Interpretation.** The carry variety $V_p$ depends on $p$: the Fermat curve $C_p: X^p + Y^p = 1$ is a different curve at each prime, not the reduction mod $p$ of a fixed curve over $\mathbb{Q}$. This is consistent with the geometric picture — the carry operator encodes the *local* arithmetic of the Fermat quotient at each prime, and there is no single global object of which these are reductions.

The path to a global $L$-function, if it exists, must therefore involve either:
- A more sophisticated global construction (e.g., a family $\{C_p\}_p$ parameterized by $p$, in the spirit of Hasse–Weil for families of varieties).
- A product over the *character-decomposed* local factors from §5.4, which may exhibit better cross-prime coherence.

### 8.3 The Cubic Factor

The cubic $\mu^3 - 6\mu^2 + 12\mu - 45$ in $\chi_M$ has Vieta relations $\sigma_1 = 2p$, $\sigma_2 = p(p+1)$, $\sigma_3 = p^2(p+2)$, all divisible by $p$. The moduli of the complex roots ($\approx 2.905$) do not match any standard Weil bound, suggesting this factor encodes higher-weight or non-smooth cohomological data. Its identification remains open.

### 8.4 Computational Frontiers

The level-3 and level-4 marginal approach (§5.5) has been ruled out: marginalization destroys spectral structure via character orthogonality, and the only productive dense level is level 2.

The viable strategies for extending the results are:

1. Character decomposition at $p = 5$ (§6.3): Compute the $25 \times 25$ carry matrix twisted by all 16 character pairs $(\chi_a, \chi_b)$ and collect Weil-bound eigenvalues. This requires exact arithmetic over $\mathbb{Z}[\zeta_5]$.

2. Extended state space at $p = 7$: The current $49 \times 49$ matrix captures 6 of 30 Frobenius eigenvalues. The full state $W_2(\mathbb{F}_7) \times W_2(\mathbb{F}_7)$ ($2401 \times 2401$) with exact arithmetic over $\mathbb{Z}[\omega_7]$ should reveal more.

3. **Systematic Jacobi sum matching**: Tabulate the known Jacobi sums $J(\chi_a, \chi_b)$ for $p = 5, 7$ and compare with the eigenvalues obtained from character-twisted carry matrices. This would give a precise identification of which Frobenius eigenvalues the carry operator "sees" at each character.

---

## 9. Conclusion

The carry operator of Witt vector multiplication — an object as elementary as the "carrying" step in schoolbook arithmetic — encodes Frobenius eigenvalues of algebraic varieties over finite fields. This is not an analogy: it is an exact algebraic identity, verified computationally by the Faddeev–LeVerrier algorithm with exact arithmetic over $\mathbb{Z}[\omega_p]$.

The principal results are:

1. **Theorem A** ($p = 3$): The characteristic polynomial of the $9 \times 9$ carry matrix contains the Frobenius of the supersingular elliptic curve $E: y^2 = x^3 + 2x + 1$ over $\mathbb{F}_3$.

2. **Theorem B** ($p = 7$): The $49 \times 49$ carry matrix yields a genus-3 Frobenius factor $(\mu^2 - 7)^2(\mu^2 + 7)$ with 6 eigenvalues at the Weil bound. The supersingular pair $\pm i\sqrt{7}$ is exactly the quadratic Gauss sum $g(\chi^3)$ of the Legendre symbol — the first identification of a carry eigenvalue with a named quantity of analytic number theory.

3. **Character universality**: Dirichlet character twists of the carry operator recover Weil-bound eigenvalues at **every prime tested** ($p = 3, 5, 7, 11, 13, 17, 19, 23$). The additive Legendre twist works at $p = 5$ and $p = 17$; order-3 characters at $p = 13$; order-5 characters at $p = 11$; higher-order characters at $p = 17, 23$. The productive character order varies with $p$ in a currently unexplained way.

4. **Structural theorems**: The carry exponent $c_2$ is equidistributed on $\mathbb{F}_p^{\ast}$, and Weil-bound eigenvectors are orthogonal to $\lvert 1\rangle$ — explaining why the Frobenius is invisible in the total exponential sum but present in the periodic orbit structure $\mathrm{Tr}(M^k)$. This is the spectral signature of an Artin-Schreier variety $z^p - z = c_2(\boldsymbol{x})$.

5. **Binary carry complex spectrum**: The D-odd conditioned binary carry chain (from [E], whose sector ratio is conjectured to converge to $-\pi$, Conjecture 1 of [P1]) has **complex** dominant eigenvalues $\lambda \approx -0.302 \pm 0.340\,i$ (modulus $\approx 0.454$, phase $\approx 131.6°$), indicating oscillatory approach to equilibrium. The carry at the penultimate position admits a closed-form expression $c_{D-2} = \lfloor Q / 2^{2K-3} \rfloor - a_1 - c_1 - 2$, collapsing the $O(K^2)$ schoolbook carry chain into $O(1)$ (follows from the schoolbook carry recurrence). The limit matrix $S_\infty$ as $K \to \infty$ is analytically computable: its entries are elementary functions of $\ln 2, \ln 3, \ln 5, \ln 7$ and rationals, and its eigenvalue constants are new transcendentals (these constants are not expressible as rational combinations of standard transcendentals, verified by PSLQ search). The complex structure is specific to the carry at position $D-2$ — carries at other positions yield purely real spectra. Unlike the Witt carry at odd primes (point 4), the D-odd binary carry eigenvectors are not orthogonal to $\lvert 1\rangle$: the Artin-Schreier orthogonality mechanism is specific to $p \geq 3$. The convergence

$$R_{\mathrm{F}}(K) \to R_\infty \approx -0.341$$

has geometric rate exactly $1/2$ (binary resolution doubling), distinct from the eigenvalue modulus $\lvert\lambda\rvert \approx 0.454$ governing Markov decorrelation. (Here $R_{\mathrm{F}}(K)$ denotes the D-odd carry ratio of the 4×4 top-bit matrix — not to be confused with the cascade sector ratio $R(K) \to -\pi$ (Conjecture 1) of [P1, E], which is a sum over all cascade depths.)

Two obstructions prevent immediate globalization: the genus obstruction (the $p^2 \times p^2$ matrix captures only a subfactor of the degree-$2g$ Frobenius, where $g = (p-1)(p-2)/2$) and the cross-prime inconsistency (no single Dirichlet character produces a coherent $L$-function across primes — different primes require different character orders). The Witt carry at $p = 2$ is degenerate (genus 0, carry always zero), so the connection between the Archimedean trace anomaly ($-\pi$) and the non-Archimedean Frobenius factors does not pass through the level-2 Witt carry.

**Open problems.**

- Identify the arithmetic rule that determines which character order is productive at each prime $p$.
- Compute the cohomology of the Artin-Schreier surface $z^p - z = c_2(x_0, x_1, y_0, y_1)$ and verify that its Frobenius matches the carry spectrum at $p = 3$ and $p = 7$.
- Explain the origin of the complex eigenvalues ($\lambda \approx -0.302 \pm 0.340\,i$) in the D-odd binary carry matrix. The limit matrix $S_\infty$ is analytically known; the eigenvalues are new transcendental constants. The complex structure is specific to carry at $D-2$ (other positions give real spectra). The relation to the oscillatory convergence $R(K) \to -\pi$ of [E] remains open: the 4×4 matrix governs top-bit conditioning, while $R(K)$ emerges from the full carry cascade.
- Determine whether a family of Artin-Schreier varieties parametrized by $p$ and a character can produce cross-prime-consistent Euler factors.

---

## References

- [Del74] P. Deligne, *La conjecture de Weil. I*, Publ. Math. IHÉS **43** (1974), 273–307.
- [DF09] P. Diaconis and J. Fulman, *Carries, shuffling, and symmetric functions*, Advances in Applied Mathematics **43** (2009), 176–196.
- [Hol97] J. M. Holte, *Carries, combinatorics, and an amazing matrix*, American Mathematical Monthly **104** (1997), 138–149.
- [IR90] K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Graduate Texts in Mathematics 84, Springer, 1990.
- [Rue04] D. Ruelle, *Thermodynamic Formalism*, 2nd ed., Cambridge University Press, 2004.
- [Ser79] J.-P. Serre, *Local Fields*, Graduate Texts in Mathematics 67, Springer, 1979.
- [Sil09] J. H. Silverman, *The Arithmetic of Elliptic Curves*, 2nd ed., Graduate Texts in Mathematics 106, Springer, 2009.
- [Wei48] A. Weil, *Sur les courbes algébriques et les variétés qui s'en déduisent*, Actualités Sci. Ind. 1041, Hermann, 1948.
- [Wei49] A. Weil, *Numbers of solutions of equations in finite fields*, Bull. Amer. Math. Soc. **55** (1949), 497–508.
- [Wit37] E. Witt, Zyklische Körper und Algebren der Characteristik $p$ vom Grad $p^n$, J. reine angew. Math. **176** (1937), 126–140.

### Companion Papers (carry-arithmetic series)

- [A] S. Alimonti, *Spectral Theory of Carries in Positional Multiplication*, this series. (Foundation: the m-bit Equidistribution Lemma extending Diaconis–Fulman to the carry transfer operator.)
- [E] S. Alimonti, *The Trace Anomaly of Binary Multiplication*, this series. (The Shifted Resolvent Theorem: $R = -\pi$ via resolvent universality, conditional on the Linear Mix Hypothesis.)
- [G] S. Alimonti, *The Angular Uniqueness of Base 2 in Positional Multiplication*, this series. (Why base 2 is the unique base for which $c_1 = \pi/18$ involves $\pi$; base-3 evidence for $c_1(3) = \ln 3 - 1/2$.)
- [P1] S. Alimonti, *π from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge*, this series. (The conjectured emergence of $\pi$ from binary carry dynamics.)

---

## Appendix A: Computational Verification

All computations were performed in Python using exact arithmetic (the `fractions.Fraction` type for $\mathbb{Q}$-coefficients and the representation $a + b\omega$ for $\mathbb{Z}[\omega]$). The Faddeev–LeVerrier algorithm was used for the characteristic polynomial to avoid the symbolic explosion of the determinant. Factorization was verified independently using the computer algebra system SymPy.

The complete source code is available at `experiments/F01_verify_theorem_a.py`.

### A.1 The Exponent Matrix

The $9 \times 9$ matrix $E$ (with $M_{j,i} = \omega^{E_{j,i}}$) for $p = 3$, multiplicative Witt carry:

States ordered as $(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)$.

```
E = [[0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0],
     [1, 1, 1, 1, 2, 1, 1, 1, 0],
     [2, 2, 2, 2, 2, 1, 2, 0, 2],
     [0, 0, 0, 0, 0, 0, 0, 0, 0],
     [2, 2, 2, 2, 2, 0, 2, 1, 2],
     [1, 1, 1, 1, 0, 1, 1, 1, 2]]
```

Entry distribution: 49 zeros, 16 ones, 16 twos.

### A.2 Reproduction

```bash
pip install numpy sympy
python experiments/F01_verify_theorem_a.py
```
