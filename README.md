# Dual-LP Certified Lower Bounds and Spectral Matrix Invariants for Giuga's Primality Conjecture

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24539270.svg)](https://doi.org/10.5281/zenodo.24539270)

**Authors:** Shubham Barmase, Gemini 4 Neon  
**Date:** November 25, 2026

## Main Result

We establish an unconditional, self-certifying lower bound showing that any composite counterexample $n = \prod_{i=1}^k p_i$ to **Giuga's Primality Conjecture** ($\sum_{j=1}^{n-1} j^{n-1} \equiv -1 \pmod n$) must satisfy:

- **$k \ge 8{,}138$ distinct prime factors**
- **Largest prime factor $p_k > 265{,}000$**
- **$n > 10^{40{,}730}$** (at least **40,731 decimal digits**)

This improves the standing rigorous lower bound of Borwein, Maitland, and Skerritt (2013) ($k \ge 4{,}771$ prime factors, $19{,}908$ decimal digits) by more than a factor of two, and surpasses the heuristic estimate of Sorini (2001) ($10^{36{,}067}$).

## Contents

- [`giuga_40730.pdf`](giuga_40730.pdf) — Complete 8-page `amsart` paper with full proofs, tables, figures, and references.
- [`giuga_40730.tex`](giuga_40730.tex) — LaTeX source.
- [`GiugaLPBound.lean`](GiugaLPBound.lean) — Pure, zero-dependency Lean 4 (`v4.36.0`) kernel formalization of the theoretical spine (`0 errors`, `0 warnings`, `0 sorry`).
- [`verify_giuga_40730_strict.py`](verify_giuga_40730_strict.py) — Standalone 50-digit directed-rounding interval-arithmetic + exact arbitrary-precision integer Weak-LP-Duality certificate verifier (zero IEEE-754 float trust).
- [`generate_figures.py`](generate_figures.py) & `figures/` — Reproducible vector figures.

## Independent Verification

### 1. Lean 4 Kernel Check (~5 seconds)
```bash
lean GiugaLPBound.lean
```
Verifies `giuga_odd`, `no_divisor_chain`, `odd_k_forces_M_ge_2`, `three_div_forces_k_odd`, and `weak_lp_duality_certificate` with zero external dependencies.

### 2. 50-Digit Directed-Rounding Dual-LP Certificate (~15 seconds)
```bash
pip install numpy scipy
python verify_giuga_40730_strict.py
```
Independently solves the LP on the 27,685-prime, 53,058-edge prime divisibility DAG ($X = 321{,}000$) and verifies every dual inequality in exact integer arithmetic scaled by $10^{12}$.
