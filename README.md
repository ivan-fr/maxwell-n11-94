# Maxwell's Electrostatic Conjecture: $M^+(11) \ge 94$ and $M^+(2k+1) \ge 20k - 6$

[![Certified Rust Proof CI](https://github.com/ivan-fr/maxwell-n11-94/actions/workflows/ci.yml/badge.svg)](https://github.com/ivan-fr/maxwell-n11-94/actions/workflows/ci.yml)

**Certified Pure-Rust Proof via `pandsysroots 0.34.0` and Pandrosion V31 Order-4 Projective Radical Enclosures**

Ivan Besevic · ESIEE-IT, France · ORCID `0009-0004-4124-1465`

## Main Theorem

There exists an explicit, exact-rational configuration of **eleven strictly positive point charges** $(a_j, q_j)_{j=1}^{11} \subset \mathbb{R}^3 \times \mathbb{Q}_{>0}$ (with coordinates and weights in $\mathbb{Z}/10^{90}$) whose three-dimensional Coulomb potential
$$V(x) = \sum_{j=1}^{11} \frac{q_j}{\|x - a_j\|}$$
admits **at least $94$ distinct, finite, non-degenerate critical points** ($c_1 = 42$ index-$1$ saddles and $c_2 = 52$ index-$2$ saddles, matching the exact Poincaré–Hopf identity $c_2 - c_1 = 11 - 1 = 10$).

Consequently:
- **$M^+(11) \ge 94$**, improving the previous best known all-positive lower bound $M^+(11) \ge 88$ (obtained by applying the Arathoon–Ball–Kvalheim $+20$ axial-pair lift to Kintali's $M^+(9) \ge 68$) by **$+6$**.
- **Infinite Odd Family Corollary:** By the Arathoon–Ball–Kvalheim $+20$ axial-pair lift,
  $$M^+(2k+1) \ge 20k - 6 \quad \text{for every integer } k \ge 5$$
  (giving $M^+(13) \ge 114$, $M^+(15) \ge 134$, $M^+(17) \ge 154$, $\dots$, improving every odd bound for $n \ge 11$ by $+6$).

## Four-Scale Harmonic Architecture

| Scale | Mechanism | Equilibria | Morse $(c_1, c_2)$ | $c_2 - c_1$ | Radial Range $\|x\|$ |
|---|---|---:|---:|---:|---|
| 1 | $100\%$-Bézout $C_{3v}$ quartic core $P_{27}$ ($27/27$ real) | $27$ | $(16, 11)$ | $-5$ | $[7.20\times 10^{-46},\; 1.85\times 10^{-13}]$ |
| 2 | Quintic harmonic nesting ($P_{27} + \varepsilon_5 H_5$, $52/64$ real) | $+25$ | $(10, 15)$ | $+5$ | $[1.92\times 10^{-12},\; 1.56\times 10^{-11}]$ |
| 3 | Sextic nesting against 11-charge 5-design $H_6$ ($71/125$ real) | $+19$ | $(9, 10)$ | $+1$ | $[8.33\times 10^{-9},\; 5.65\times 10^{-8}]$ |
| 4 | Multi-shell $n=11$ positive Kelvin 5-design outer equilibria | $+23$ | $(7, 16)$ | $+9$ | $[3.60\times 10^{-2},\; 8.01\times 10^{-1}]$ |
| **Total** | **Exact 11-charge positive Coulomb potential $V(x)$** | **$94$** | **$(42, 52)$** | **$+10$** | **All $94$ disjoint, min pole dist $0.033859$** |

## Verification Methodology (Zero `f64`/`libm` Trust)

1. **Limiting Polynomial Certification (`pandsysroots 0.34.0`)**:
   - Certifies the degree-$4$ core $P_{27}$ (`27/27` certified, `complete = true`, all `27` real) and the degree-$6$ nested harmonic limit $P_{71}$ (`125/125` certified, `complete = true`, `71` real, `c1 = 35, c2 = 36`).
2. **Direct 3D Coulomb Certification (512-bit Dyadic Intervals + Pandrosion V31 + Krawczyk)**:
   - Every rational charge $(a_j, q_j) \in (\mathbb{Z}/10^{90})^4$ is enclosed in 512-bit outward-rounded dyadic intervals with certified $q_j^{\mathrm{lo}} > 0$.
   - Every Euclidean distance $r_j = \sqrt{\|x - a_j\|^2}$ is enclosed via the **Pandrosion V31 projective order-4 two-sided radical bounds**
     $$\frac{4 s S (s^2 + S)}{S^2 + 6 s^2 S + s^4} \le \sqrt{S} \le \frac{S^2 + 6 s^2 S + s^4}{4 s (s^2 + S)},$$
     double-checked by exact integer squaring ($r_{\mathrm{lo}}^2 \le S_{\mathrm{lo}}$ and $r_{\mathrm{hi}}^2 \ge S_{\mathrm{hi}}$).
   - Every one of the $94$ critical points passes strict 3D interval Krawczyk inclusion $K(B_k) \subset \operatorname{int}(B_k)$ with contraction norm $\|I - Y_k \operatorname{Hess} V(B_k)\|_\infty \le 1.386 \times 10^{-21} < 1$, strict non-zero determinant $\det \operatorname{Hess} V(B_k)$ certifying exact Morse index, pairwise disjointness of all $94$ boxes, and strict separation from all $11$ charge poles ($\ge 0.033859$).

## Reproduce in < 1 Second

```bash
cargo test --release -- --nocapture
cargo run --release
```

## Files

- [`src/dyadic.rs`](src/dyadic.rs): 512-bit exact outward-rounded dyadic intervals and Pandrosion V31 order-4 projective radical enclosure.
- [`src/coulomb.rs`](src/coulomb.rs): 3D Coulomb gradient, Hessian, and interval Krawczyk certifier.
- [`src/limiting_poly.rs`](src/limiting_poly.rs): `pandsysroots 0.34.0` complete certification of $P_{27}$ and $P_{71}$.
- [`data/witness_n11_94.json`](data/witness_n11_94.json): Exact 90-decimal rational coordinates of the 11 positive charges and all 94 critical-point centers.
- [`paper/Maxwell_Eleven_Charges_94.tex`](paper/Maxwell_Eleven_Charges_94.tex) & [`paper/Maxwell_Eleven_Charges_94.pdf`](paper/Maxwell_Eleven_Charges_94.pdf): Full manuscript.
