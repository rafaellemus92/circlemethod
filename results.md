# Results — the circle method for binary Goldbach at N = 10^6, made visible

Exploration and visualisation only: nothing below is, or is meant as, a proof or a proof strategy. All numbers are produced by `python run_all.py` (one run, 342 s on 4 cores); every figure is in `figures/`.

## 0. Anchors and quadrature certificates

| quantity | value | check |
|---:|---:|---:|
| π(10^6) | 78,498 | = 78,498 ✓ |
| G(10^6) (direct count) | 5,402 | = 5,402 ✓ |
| ψ(N) = S(0) | 999,586.597 | ψ/N − 1 = -4.13e-04 (within 0.2% ✓) |
| prime-power weight Σ_{k≥2} θ(N^{1/k}) | 1,102.4 | O(√N), 0.11% of ψ |
| Σ Λ(n)² | 1.280299e7 | N(log N − 1) = 1.281551e7 |
| R2(N) direct | 1,756,567.536732 | R2/(SS2·N) = 0.997805 |
| SS2(N) = 2C2·(5−1)/(5−2), C2 = literature value | 1.7604315089 | Σ_{q≤10^6} = 1.7604315083 |
| C2 (our product over p ≤ 10^7) | 0.6601618197 | literature 0.6601618158 |
| Parseval on grid (1/M)Σ\|S\|² | 1.280298648413e7 | rel. err 0.0e+00 |
| (1/M)Σ S² e(−Nk/M) | 1,756,567.536732 | = R2, rel. err 1.3e-16 |
| same with the wrong FFT sign | -2.40e-10 | 0 to rounding (R2 ≈ 1.8×10⁶) |
| R3(999,999) via FFT convolution | 7.279938865e11 | grid ∫S³e(−N'α): rel. err 6.7e-16 |
| R3(1,000,001) via FFT convolution | 1.150440137e12 | grid ∫S³e(−N'α): rel. err 0.0e+00 |

**FFT sign convention.** `numpy.fft.fft` uses exp(−2πi nk/M). We need S(k/M) = Σ Λ(n) e(+nk/M), so S = M·ifft(Λ) = conj(fft(Λ)), evaluated with the half-length `rfft` (`expsum.S_grid`). S(0) = ψ, S(1/2) = 2⌊log₂N⌋log 2 − ψ = -999,560.2579 (predicted -999,560.2579) and S(−α) = conj S(α) hold for either sign, so they do not settle it. Two checks do: Im S(1/M) = 357,037.4019 on the grid vs 357,037.4019 by direct summation (same sign), and the binary grid integral returns R2(N) exactly with the chosen sign but -2.4e-10 with the other (the wrong sign evaluates the autoconvolution at M − N > 2N, where it vanishes). One quadrature rule is used throughout: the rectangle rule on k/M, k = 0..M−1, M = 2^23 = 8,388,608. On the full circle it is exact for R2 (needs M > N), for R3 (needs M > 2N') and for Parseval; on arcs it is a Riemann sum with membership decided exactly in integers, so I_M + I_m = R2 to rounding (max |closure| = 5.3e-16).

## Headline table

All (N, Q, c) evaluated. B3/main is the ternary ratio for N' = N−1 / N+1 (odd), using the same S. For N = 10^6 all 18 grid pairs (Q, c) satisfy 2cQ² < N (largest: 2·20·150² = 900,000), so nothing was skipped; the rows (158, 20) and (223, 10) for N = 10^6 are the Qmax pairs of §7. Rows marked * are the smooth-N probe (§7).

| N | Q | c | R2 | SS2·N | I_M | I_m | P_m | P_m/R2 | sup_m\|S\| | B3/main N−1 / N+1 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000,000 | 5 | 5 | 1.7566e6 | 1.7604e6 | 1.9577e6 | -2.012e+05 | 1.011e7 | 5.756 | 4.998e5 | 8.789 / 5.562 |
| 1,000,000 | 5 | 10 | 1.7566e6 | 1.7604e6 | 1.9772e6 | -2.207e+05 | 1.008e7 | 5.740 | 4.998e5 | 8.789 / 5.562 |
| 1,000,000 | 5 | 20 | 1.7566e6 | 1.7604e6 | 1.9873e6 | -2.307e+05 | 1.007e7 | 5.732 | 4.998e5 | 8.789 / 5.562 |
| 1,000,000 | 10 | 5 | 1.7566e6 | 1.7604e6 | 1.9295e6 | -1.730e+05 | 9.212e6 | 5.244 | 1.667e5 | 2.932 / 1.855 |
| 1,000,000 | 10 | 10 | 1.7566e6 | 1.7604e6 | 1.9487e6 | -1.921e+05 | 9.175e6 | 5.223 | 1.667e5 | 2.932 / 1.855 |
| 1,000,000 | 10 | 20 | 1.7566e6 | 1.7604e6 | 1.9585e6 | -2.020e+05 | 9.155e6 | 5.212 | 1.667e5 | 2.932 / 1.855 |
| 1,000,000 | 20 | 5 | 1.7566e6 | 1.7604e6 | 1.8165e6 | -5.993e+04 | 8.629e6 | 4.913 | 1.252e5 | 2.201 / 1.393 |
| 1,000,000 | 20 | 10 | 1.7566e6 | 1.7604e6 | 1.8346e6 | -7.805e+04 | 8.583e6 | 4.886 | 1.252e5 | 2.201 / 1.393 |
| 1,000,000 | 20 | 20 | 1.7566e6 | 1.7604e6 | 1.8437e6 | -8.717e+04 | 8.551e6 | 4.868 | 1.252e5 | 2.201 / 1.393 |
| 1,000,000 | 50 | 5 | 1.7566e6 | 1.7604e6 | 1.7316e6 | +2.493e+04 | 7.646e6 | 4.353 | 5.723e4 | 1.007 / 0.637 |
| 1,000,000 | 50 | 10 | 1.7566e6 | 1.7604e6 | 1.7491e6 | +7.426e+03 | 7.563e6 | 4.306 | 5.170e4 | 0.909 / 0.575 |
| 1,000,000 | 50 | 20 | 1.7566e6 | 1.7604e6 | 1.7568e6 | -1.928e+02 | 7.460e6 | 4.247 | 5.170e4 | 0.909 / 0.575 |
| 1,000,000 | 100 | 5 | 1.7566e6 | 1.7604e6 | 1.7165e6 | +4.010e+04 | 6.908e6 | 3.933 | 5.723e4 | 1.007 / 0.637 |
| 1,000,000 | 100 | 10 | 1.7566e6 | 1.7604e6 | 1.7334e6 | +2.322e+04 | 6.704e6 | 3.817 | 3.219e4 | 0.566 / 0.358 |
| 1,000,000 | 100 | 20 | 1.7566e6 | 1.7604e6 | 1.7409e6 | +1.567e+04 | 6.332e6 | 3.605 | 3.219e4 | 0.566 / 0.358 |
| 1,000,000 | 150 | 5 | 1.7566e6 | 1.7604e6 | 1.7148e6 | +4.172e+04 | 6.329e6 | 3.603 | 5.723e4 | 1.007 / 0.637 |
| 1,000,000 | 150 | 10 | 1.7566e6 | 1.7604e6 | 1.7310e6 | +2.553e+04 | 5.908e6 | 3.363 | 3.213e4 | 0.565 / 0.358 |
| 1,000,000 | 150 | 20 | 1.7566e6 | 1.7604e6 | 1.7393e6 | +1.725e+04 | 5.046e6 | 2.873 | 2.237e4 | 0.393 / 0.249 |
| 10,000 | 15 | 20 | 1.6855e4 | 1.7604e4 | 1.7926e4 | -1.071e+03 | 3.204e4 | 1.901 | 1.281e3 | 1.381 / 0.911 |
| 10,000 | 20 | 10 | 1.6855e4 | 1.7604e4 | 1.7856e4 | -1.001e+03 | 3.293e4 | 1.954 | 1.281e3 | 1.381 / 0.911 |
| 10,000 | 22 | 10 | 1.6855e4 | 1.7604e4 | 1.7765e4 | -9.103e+02 | 2.985e4 | 1.771 | 1.281e3 | 1.381 / 0.911 |
| 30,000 | 20 | 10 | 1.0368e5 | 1.0563e5 | 1.0494e5 | -1.255e+03 | 1.448e5 | 1.396 | 3.796e3 | 1.019 / 1.022 |
| 30,000 | 27 | 20 | 1.0368e5 | 1.0563e5 | 1.0297e5 | +7.110e+02 | 1.087e5 | 1.048 | 3.796e3 | 1.019 / 1.022 |
| 30,000 | 38 | 10 | 1.0368e5 | 1.0563e5 | 1.0639e5 | -2.704e+03 | 1.005e5 | 0.969 | 2.605e3 | 0.699 / 0.702 |
| 100,000 | 20 | 10 | 1.7657e5 | 1.7604e5 | 1.8312e5 | -6.545e+03 | 6.211e5 | 3.517 | 1.263e4 | 1.730 / 1.165 |
| 100,000 | 49 | 20 | 1.7657e5 | 1.7604e5 | 1.7665e5 | -8.110e+01 | 3.995e5 | 2.262 | 5.277e3 | 0.723 / 0.487 |
| 100,000 | 50 | 10 | 1.7657e5 | 1.7604e5 | 1.7561e5 | +9.592e+02 | 4.754e5 | 2.692 | 5.277e3 | 0.723 / 0.487 |
| 100,000 | 70 | 10 | 1.7657e5 | 1.7604e5 | 1.7684e5 | -2.706e+02 | 3.783e5 | 2.143 | 4.287e3 | 0.587 / 0.396 |
| 300,000 | 20 | 10 | 1.0452e6 | 1.0563e6 | 1.0506e6 | -5.363e+03 | 2.206e6 | 2.111 | 3.779e4 | 1.316 / 1.279 |
| 300,000 | 50 | 10 | 1.0452e6 | 1.0563e6 | 1.0649e6 | -1.965e+04 | 1.867e6 | 1.786 | 1.531e4 | 0.533 / 0.518 |
| 300,000 | 86 | 20 | 1.0452e6 | 1.0563e6 | 1.0551e6 | -9.906e+03 | 1.313e6 | 1.256 | 9.938e3 | 0.346 / 0.336 |
| 300,000 | 122 | 10 | 1.0452e6 | 1.0563e6 | 1.0475e6 | -2.317e+03 | 1.247e6 | 1.193 | 9.184e3 | 0.320 / 0.311 |
| 1,000,000 | 158 | 20 | 1.7566e6 | 1.7604e6 | 1.7381e6 | +1.850e+04 | 4.812e6 | 2.740 | 2.237e4 | 0.393 / 0.249 |
| 1,000,000 | 223 | 10 | 1.7566e6 | 1.7604e6 | 1.7322e6 | +2.439e+04 | 4.594e6 | 2.615 | 3.213e4 | 0.565 / 0.358 |
| 1,531,530* | 20 | 10 | 8.3614e6 | 8.3662e6 | 6.3733e6 | +1.988e+06 | 1.381e7 | 1.651 | 1.921e5 | 1.443 / 1.444 |
| 1,531,530* | 50 | 10 | 8.3614e6 | 8.3662e6 | 7.3738e6 | +9.876e+05 | 1.227e7 | 1.467 | 7.766e4 | 0.583 / 0.584 |
| 1,531,530* | 195 | 20 | 8.3614e6 | 8.3662e6 | 8.0775e6 | +2.838e+05 | 7.608e6 | 0.910 | 3.399e4 | 0.255 / 0.256 |
| 1,531,530* | 276 | 10 | 8.3614e6 | 8.3662e6 | 8.1175e6 | +2.439e+05 | 7.303e6 | 0.873 | 4.638e4 | 0.348 / 0.349 |
| 2,000,000 | 20 | 10 | 3.5243e6 | 3.5209e6 | 3.6782e6 | -1.538e+05 | 1.857e7 | 5.269 | 2.514e5 | 1.482 / 2.213 |
| 2,000,000 | 50 | 10 | 3.5243e6 | 3.5209e6 | 3.5143e6 | +9.994e+03 | 1.657e7 | 4.702 | 1.014e5 | 0.598 / 0.893 |
| 2,000,000 | 223 | 20 | 3.5243e6 | 3.5209e6 | 3.5071e6 | +1.719e+04 | 1.008e7 | 2.861 | 3.373e4 | 0.199 / 0.297 |
| 2,000,000 | 316 | 10 | 3.5243e6 | 3.5209e6 | 3.4907e6 | +3.360e+04 | 9.712e6 | 2.756 | 5.868e4 | 0.346 / 0.517 |

## 1–3. Spikes, zooms, phasors

* Squarefree q ≤ 30: max | |S(a/q)| / (ψ/φ(q)) − 1 | = 2.70e-02 (worst a/q = 27/29); sign of Re S(a/q) equals μ(q) for every one: True.
* Non-squarefree q ≤ 30: max |S(a/q)|/ψ = 1.73e-03.
* S(1/4) = 11.7835 -835.8223 i. Real part = ψ(N;4,0) − ψ(N;4,2) = 11.7835 (powers of 2); imaginary part = ψ(N;4,1) − ψ(N;4,3) = -835.8223. So |Im S(1/4)| is exactly the mod-4 race ψ(N;4,3) − ψ(N;4,1) = 835.8. For primes alone the race is θ(N;4,3) − θ(N;4,1) = 1816.6 ≈ 1.8√N. About θ(√N) ≈ 956 of that is Chebyshev's bias proper — the prime squares p² ≡ 1 (mod 4) that θ leaves out; Λ puts them back and so removes the bias. What remains in the ψ race, 836 ≈ 0.84√N, is the oscillation driven by the zeros of L(s, χ₄), which happens to be positive at N = 10^6 and has no preferred sign.
* Zoom around 1/3: on |Nβ| ≤ 10 the measured S differs from the model −v(β)/2 by at most 4.56e-03·(ψ/2).
* Golden-ratio window: max |S| = 5,506 = 1.54·√ΣΛ², at N·β = -8.13, on the crest of 144/233 (model ψ/φ(q) = 4,309, N·β = -8.24). 2 model crests in the window exceed √ΣΛ² = 3,578. Exact S(a/q) at the 20 squarefree fractions with crest ≥ 10³: real part has the sign of μ(q) in 19 cases — S(144/233) = -5,397 -206i (μ = -1); S(377/610) = -3,206 +910i (μ = -1); S(343/555) = -3,459 +1,130i (μ = -1); S(233/377) = +3,415 -637i (μ = +1); S(775/1254) = +3,000 +42i (μ = +1).
* Phasor resultants: S(1/3) = -499,774 +619i; S(1/4) = 12 -836i; S(1/5) = -249,800 -72i; S(1/6) = 499,760 +620i; S(1/30) = -124,627 -391i; S(frac(golden ratio)) = -2,837 +325i; S(frac(sqrt 2)) = -570 -565i.

## 4. Major arcs vs truth (binary, N = 10^6)

(i) I_M from the grid; (ii) the same with S replaced on each arc by (μ(q)/φ(q))·v(β); (iii) SS2_Q(N)·N; (iv) SS2(N)·N from the product; (v) R2(N) exact. Relative errors are against R2. J(c)/N is the fraction of ∫v² e(−Nβ) inside |β| ≤ c/N: c=5: 0.97978, c=10: 0.98987, c=20: 0.99493 (≈ 1 − 1/(π²c)).

| Q | c | \|M\| | (i) I_M | (ii) model | (iii) SS2_Q·N | (iv) SS2·N | (v) R2 | (i)/R2−1 | (ii)/R2−1 | (iii)/R2−1 | (iv)/R2−1 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 5 | 0.0001 | 1.95773e6 | 1.95956e6 | 2.00000e6 | 1.76043e6 | 1.75657e6 | +0.1145 | +0.1156 | +0.1386 | +0.0022 |
| 5 | 10 | 0.0002 | 1.97724e6 | 1.97975e6 | 2.00000e6 | 1.76043e6 | 1.75657e6 | +0.1256 | +0.1271 | +0.1386 | +0.0022 |
| 5 | 20 | 0.0004 | 1.98728e6 | 1.98987e6 | 2.00000e6 | 1.76043e6 | 1.75657e6 | +0.1313 | +0.1328 | +0.1386 | +0.0022 |
| 10 | 5 | 0.0003 | 1.92955e6 | 1.93234e6 | 1.97222e6 | 1.76043e6 | 1.75657e6 | +0.0985 | +0.1001 | +0.1228 | +0.0022 |
| 10 | 10 | 0.0006 | 1.94872e6 | 1.95225e6 | 1.97222e6 | 1.76043e6 | 1.75657e6 | +0.1094 | +0.1114 | +0.1228 | +0.0022 |
| 10 | 20 | 0.0013 | 1.95854e6 | 1.96223e6 | 1.97222e6 | 1.76043e6 | 1.75657e6 | +0.1150 | +0.1171 | +0.1228 | +0.0022 |
| 20 | 5 | 0.0013 | 1.81650e6 | 1.82043e6 | 1.85801e6 | 1.76043e6 | 1.75657e6 | +0.0341 | +0.0364 | +0.0577 | +0.0022 |
| 20 | 10 | 0.0026 | 1.83461e6 | 1.83919e6 | 1.85801e6 | 1.76043e6 | 1.75657e6 | +0.0444 | +0.0470 | +0.0577 | +0.0022 |
| 20 | 20 | 0.0051 | 1.84374e6 | 1.84860e6 | 1.85801e6 | 1.76043e6 | 1.75657e6 | +0.0496 | +0.0524 | +0.0577 | +0.0022 |
| 50 | 5 | 0.0077 | 1.73164e6 | 1.73793e6 | 1.77380e6 | 1.76043e6 | 1.75657e6 | -0.0142 | -0.0106 | +0.0098 | +0.0022 |
| 50 | 10 | 0.0155 | 1.74914e6 | 1.75583e6 | 1.77380e6 | 1.76043e6 | 1.75657e6 | -0.0042 | -0.0004 | +0.0098 | +0.0022 |
| 50 | 20 | 0.0310 | 1.75676e6 | 1.76481e6 | 1.77380e6 | 1.76043e6 | 1.75657e6 | +0.0001 | +0.0047 | +0.0098 | +0.0022 |
| 100 | 5 | 0.0304 | 1.71647e6 | 1.72587e6 | 1.76149e6 | 1.76043e6 | 1.75657e6 | -0.0228 | -0.0175 | +0.0028 | +0.0022 |
| 100 | 10 | 0.0609 | 1.73335e6 | 1.74365e6 | 1.76149e6 | 1.76043e6 | 1.75657e6 | -0.0132 | -0.0074 | +0.0028 | +0.0022 |
| 100 | 20 | 0.1218 | 1.74090e6 | 1.75257e6 | 1.76149e6 | 1.76043e6 | 1.75657e6 | -0.0089 | -0.0023 | +0.0028 | +0.0022 |
| 150 | 5 | 0.0686 | 1.71484e6 | 1.72343e6 | 1.75900e6 | 1.76043e6 | 1.75657e6 | -0.0238 | -0.0189 | +0.0014 | +0.0022 |
| 150 | 10 | 0.1372 | 1.73103e6 | 1.74118e6 | 1.75900e6 | 1.76043e6 | 1.75657e6 | -0.0145 | -0.0088 | +0.0014 | +0.0022 |
| 150 | 20 | 0.2743 | 1.73931e6 | 1.75009e6 | 1.75900e6 | 1.76043e6 | 1.75657e6 | -0.0098 | -0.0037 | +0.0014 | +0.0022 |

Agreement is at the few-percent level once Q ≥ 20 and c ≥ 10 (and ≤ 1.5% for Q ≥ 50), as expected. The difference (i) − (ii) is the error of the major-arc model on the major arcs themselves (-0.0066 to -0.0010 of R2; part of it is simply ψ ≠ N, (ψ/N)² − 1 = -8.3e-04). SS2_Q(N) for Q = 5, 10, 20, 50, 100, 150: 2.000000, 1.972222, 1.858007, 1.773796, 1.761490, 1.758996.

## 5. The minor arcs — the headline

P_m = ∫_m |S|² (mass) and I_m = ∫_m S² e(−Nα) (signed), both summed directly over the minor-arc grid points. ‘heur.’ is the no-fit estimate P_m/R2 ≈ (log N − log Q − 2.33)/SS2(N), with 2.33 = 1 + γ + Σ_p log p/(p(p−1)) from Σ_{q≤Q} μ²(q)/φ(q) ≈ log Q + 1.33; ‘heur. exact’ uses the exact partial sum and the kernel capture: (ΣΛ² − J(c)·Σ_{q≤Q} μ²(q)/φ(q))/R2. ‘leak’ = (SS2 − SS2_Q)(N−1) + SS2_Q·(N−1−J(c)) is the deterministic part of I_m that the major-arc model predicts: the q > Q terms of the singular series (spikes the arcs leave out) plus the Dirichlet-kernel tails that spill past c/N. ‘resid’ = I_m − leak. ‘model err’ = I_M − I_M(model) from §4. σ_rp = P_m/√(N|m|) is a random-phase yardstick: the size of a sum of N|m| independent coherence cells of width 1/N carrying the mass P_m with random signs. sup_m|S| is the grid maximum refined by direct summation within ±1/M (the grid alone can miss a crest by ≤ 0.6% on a main lobe, ≤ 1.8% on a sidelobe); ‘at’ names the rational whose model spike explains it.

| Q | c | P_m/R2 | heur. | heur. exact | I_m/R2 | \|I_m\|/√P_m | leak/R2 | resid/R2 | model err/R2 | σ_rp/R2 | sup_m\|S\| | at a/q | N·β |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 5 | 5.756 | 5.609 | 5.755 | -0.11452 | 63.3 | -0.11336 | -0.00116 | -0.00104 | 0.0058 | 4.998e5 | 1/6 | -0.00 |
| 5 | 10 | 5.740 | 5.609 | 5.739 | -0.12563 | 69.5 | -0.12486 | -0.00077 | -0.00143 | 0.0057 | 4.998e5 | 1/6 | -0.00 |
| 5 | 20 | 5.732 | 5.609 | 5.731 | -0.13134 | 72.7 | -0.13062 | -0.00072 | -0.00148 | 0.0057 | 4.998e5 | 1/6 | -0.00 |
| 10 | 5 | 5.244 | 5.215 | 5.243 | -0.09848 | 57.0 | -0.09787 | -0.00061 | -0.00159 | 0.0052 | 1.667e5 | 9/14 | -0.00 |
| 10 | 10 | 5.223 | 5.215 | 5.222 | -0.10939 | 63.4 | -0.10920 | -0.00019 | -0.00201 | 0.0052 | 1.667e5 | 9/14 | -0.00 |
| 10 | 20 | 5.212 | 5.215 | 5.212 | -0.11498 | 66.8 | -0.11488 | -0.00010 | -0.00210 | 0.0052 | 1.667e5 | 9/14 | -0.00 |
| 20 | 5 | 4.913 | 4.821 | 4.913 | -0.03412 | 20.4 | -0.03416 | +0.00004 | -0.00224 | 0.0049 | 1.252e5 | 11/30 | +0.00 |
| 20 | 10 | 4.886 | 4.821 | 4.888 | -0.04443 | 26.6 | -0.04484 | +0.00041 | -0.00261 | 0.0049 | 1.252e5 | 11/30 | +0.00 |
| 20 | 20 | 4.868 | 4.821 | 4.876 | -0.04963 | 29.8 | -0.05019 | +0.00057 | -0.00277 | 0.0049 | 1.252e5 | 11/30 | +0.00 |
| 50 | 5 | 4.353 | 4.301 | 4.363 | +0.01419 | 9.0 | +0.01281 | +0.00138 | -0.00358 | 0.0044 | 5.723e4 | 0/1 | +5.48 |
| 50 | 10 | 4.306 | 4.301 | 4.333 | +0.00423 | 2.7 | +0.00262 | +0.00161 | -0.00381 | 0.0043 | 5.170e4 | 49/66 | +0.00 |
| 50 | 20 | 4.247 | 4.301 | 4.318 | -0.00011 | 0.1 | -0.00249 | +0.00238 | -0.00458 | 0.0043 | 5.170e4 | 49/66 | +0.00 |
| 100 | 5 | 3.933 | 3.907 | 3.992 | +0.02283 | 15.3 | +0.01968 | +0.00315 | -0.00535 | 0.0040 | 5.723e4 | 0/1 | +5.48 |
| 100 | 10 | 3.817 | 3.907 | 3.958 | +0.01322 | 9.0 | +0.00955 | +0.00367 | -0.00586 | 0.0039 | 3.219e4 | 61/102 | +0.00 |
| 100 | 20 | 3.605 | 3.907 | 3.941 | +0.00892 | 6.2 | +0.00448 | +0.00444 | -0.00664 | 0.0038 | 3.219e4 | 61/102 | +0.00 |
| 150 | 5 | 3.603 | 3.677 | 3.754 | +0.02375 | 16.6 | +0.02107 | +0.00269 | -0.00489 | 0.0037 | 5.723e4 | 0/1 | +5.48 |
| 150 | 10 | 3.363 | 3.677 | 3.717 | +0.01454 | 10.5 | +0.01096 | +0.00358 | -0.00578 | 0.0036 | 3.213e4 | 0/1 | +10.49 |
| 150 | 20 | 2.873 | 3.677 | 3.699 | +0.00982 | 7.7 | +0.00589 | +0.00393 | -0.00613 | 0.0034 | 2.237e4 | 149/210 | +0.01 |

* **The L² wall.** P_m/R2 runs from 2.87 to 5.76 over the 18 grid pairs, and down to
  2.62 at the Qmax pairs (158, 20) and (223, 10) of §7: any estimate of the minor arcs by their magnitude
  alone, |I_m| ≤ P_m, is too weak by that factor. For Q ≤ 20 the exact-sum heuristic is within 0.008
  (the asymptotic form within 0.15, its o(1) error at small Q). As the arcs widen (|M| from 3% to 27%)
  P_m falls increasingly below both, because wide arcs swallow part of the unresolved background as well as the
  spikes (2.87 against 3.70 at (150, 20)); it is still nearly 3× the signal. Raising Q to (log N)^B would only subtract
  B·log log N = 2.63·B from log N − 2.33 = 11.49.
* **What the minor arcs contribute.** |I_m|/R2 is 3–13% at Q ≤ 20 and 1.1e-04–2.4e-02 for
  Q ≥ 50, against a mass of 2.9–5.8·R2. At small Q nearly all of it is the predictable leakage (the
  q > Q terms of the singular series and the kernel tails); |I_m|/√P_m is 20–73 there, well above log N, because
  that part is not random.
* **What is left, and what it is not.** Because I_m = R2 − I_M, the residual satisfies the identity
  resid = (R2 − SS2·N) − (I_M − I_M(model)) + O(10⁻⁶·R2). The first term is the global Hardy–Littlewood error,
  -0.0022·R2 in every row; the second is the major-arc model error (-0.0066 to
  -0.0010·R2). So the residual is the minor-arc share of the singular-series model's error, and at this
  N it cannot be separated from the model's error on the major arcs. Its size for Q ≥ 50
  (1.4e-03–4.4e-03·R2) is comparable to the random-phase yardstick σ_rp
  (3.4e-03–4.4e-03·R2): the data are consistent with square-root cancellation of the minor-arc
  mass, but they do not measure it, and nothing here certifies it for any other N.

## 6. Ternary control

| N' | factorisation | R3 exact | SS3 (product) | SS3_150 | main = SS3·N'²/2 | R3/main |
|---:|---:|---:|---:|---:|---:|---:|
| 999,999 | 3^3·7·11·13·37 | 7.279939e11 | 1.455984 | 1.455295 | 7.279904e11 | 1.000005 |
| 1,000,001 | 101·9901 | 1.150440e12 | 2.300729 | 2.300761 | 1.150367e12 | 1.000064 |

SS3_Q truncations (Q = 5, 10, 20, 50, 100, 150, 1000, 10^4): N' = 999,999: 1.76562, 1.50347, 1.45526, 1.45205, 1.45530, 1.45529, 1.45600, 1.45598; N' = 1,000,001: 2.14062, 2.28588, 2.29446, 2.30016, 2.30072, 2.30076, 2.30073, 2.30073.

B3 = sup_m|S| · ∫₀¹|S|² (the textbook minor-arc bound, ∫|S|² = ΣΛ² = 1.2803e7), against the binary P_m/R2 at the same (Q, c):

| Q | c | sup_m\|S\| | B3 | B3/main N'=999,999 | B3/main N'=1,000,001 | I3_M/main (999,999) | I3_m/main (999,999) | binary P_m/R2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 5 | 4.998e5 | 6.398e12 | 8.789 | 5.562 | 1.21258 | -2.13e-01 | 5.756 |
| 5 | 10 | 4.998e5 | 6.398e12 | 8.789 | 5.562 | 1.21262 | -2.13e-01 | 5.740 |
| 5 | 20 | 4.998e5 | 6.398e12 | 8.789 | 5.562 | 1.21262 | -2.13e-01 | 5.732 |
| 10 | 5 | 1.667e5 | 2.134e12 | 2.932 | 1.855 | 1.03257 | -3.26e-02 | 5.244 |
| 10 | 10 | 1.667e5 | 2.134e12 | 2.932 | 1.855 | 1.03260 | -3.26e-02 | 5.223 |
| 10 | 20 | 1.667e5 | 2.134e12 | 2.932 | 1.855 | 1.03261 | -3.26e-02 | 5.212 |
| 20 | 5 | 1.252e5 | 1.602e12 | 2.201 | 1.393 | 0.99949 | +5.11e-04 | 4.913 |
| 20 | 10 | 1.252e5 | 1.602e12 | 2.201 | 1.393 | 0.99953 | +4.79e-04 | 4.886 |
| 20 | 20 | 1.252e5 | 1.602e12 | 2.201 | 1.393 | 0.99953 | +4.78e-04 | 4.868 |
| 50 | 5 | 5.723e4 | 7.327e11 | 1.007 | 0.637 | 0.99734 | +2.66e-03 | 4.353 |
| 50 | 10 | 5.170e4 | 6.619e11 | 0.909 | 0.575 | 0.99737 | +2.63e-03 | 4.306 |
| 50 | 20 | 5.170e4 | 6.619e11 | 0.909 | 0.575 | 0.99738 | +2.63e-03 | 4.247 |
| 100 | 5 | 5.723e4 | 7.327e11 | 1.007 | 0.637 | 0.99960 | +4.02e-04 | 3.933 |
| 100 | 10 | 3.219e4 | 4.122e11 | 0.566 | 0.358 | 0.99963 | +3.77e-04 | 3.817 |
| 100 | 20 | 3.219e4 | 4.122e11 | 0.566 | 0.358 | 0.99963 | +3.74e-04 | 3.605 |
| 150 | 5 | 5.723e4 | 7.327e11 | 1.007 | 0.637 | 0.99957 | +4.32e-04 | 3.603 |
| 150 | 10 | 3.213e4 | 4.114e11 | 0.565 | 0.358 | 0.99959 | +4.13e-04 | 3.363 |
| 150 | 20 | 2.237e4 | 2.864e11 | 0.393 | 0.249 | 0.99958 | +4.29e-04 | 2.873 |

The bound closes when sup_m|S| < main/ΣΛ² = 56,861 = ψ/17.6 for N' = 999,999, and
89,851 = ψ/11.1 for N' = 1,000,001. Two kinds of stray crest on the minor arcs set sup_m|S|:

* **Excluded squarefree q with small φ(q)** — these need Q to grow. The table's sup sits at 1/6 (ψ/2) for Q = 5,
  9/14 (ψ/6) for Q = 10, 11/30 (ψ/8) for Q = 20, 49/66 (ψ/20) for Q = 50 with c ≥ 10, 61/102 (ψ/32) for Q = 100 with
  c ≥ 10, and 149/210 (ψ/48, measured 22,373) at (150, 20).
* **The first Dirichlet sidelobe just outside the α = 0 (and α = 1/2) arc**, |S| ≈ ψ/(π(c+½)) — this needs c to
  grow. Measured 57,232 at c = 5 and 32,131 at c = 10 (model 57,851 and
  30,303). No arc with q ≤ Q can cover it, so at c = 5 it keeps B3/main(999,999) at
  1.007 for every Q ≥ 50; at c = 10 it (32,131) is almost exactly as tall as the
  61/102 crest (32,192), so B3/main stays at 0.565–0.566
  for Q = 100, 150 (and 223 in §7).

At the largest grid pair (150, 20) — and at the Qmax pair (158, 20), which has the same sup — the bound closes with
room to spare: B3/main = 0.393 (N' = 999,999) and 0.249 (N' = 1,000,001), against
binary P_m/R2 = 2.873 in the same row. Either way the sup is set by near-rational α. In the
asymptotic theory the arcs have width ~Q/N, so Q and the distance to the arc edge grow together (Vinogradov: Q a
power of log N). The actual ternary minor-arc share |I3_m|/main is ≤ 3×10⁻³ for Q ≥ 20.

## 7. Scaling in N

M = smallest power of two ≥ 8N for each N (alias-free for R3, ≥ 8 samples per 1/N). Series: (Q, c) = (20, 10); (50, 10) where admissible (N ≥ 10^5); and Q = Qmax(N, c) = the largest Q with 2cQ² < N, for c = 10 and 20. The extra row N = 1,531,530 = 2·3²·5·7·11·13·17 is a probe at the largest SS2 for N ≤ 2·10^6 (5.4627, shared with 510,510 and 1,021,020). ‘resid’ and σ_rp as in §5.

| N | series | Q | c | SS2(N) | P_m/R2 | P_m/N | I_m/R2 | resid/R2 | σ_rp/R2 | sup_m\|S\| | at | B3/main N−1 | B3/main N+1 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | Q=20, c=10 | 20 | 10 | 1.760 | 1.954 | 3.293 | -0.0594 | -0.0126 | 0.0227 | 1.281e3 | 23/30 | 1.381 | 0.911 |
| 10,000 | Q=Qmax, c=10 | 22 | 10 | 1.760 | 1.771 | 2.985 | -0.0540 | -0.0090 | 0.0212 | 1.281e3 | 23/30 | 1.381 | 0.911 |
| 10,000 | Q=Qmax, c=20 | 15 | 20 | 1.760 | 1.901 | 3.204 | -0.0635 | -0.0070 | 0.0225 | 1.281e3 | 23/30 | 1.381 | 0.911 |
| 30,000 | Q=20, c=10 | 20 | 10 | 3.521 | 1.396 | 4.826 | -0.0121 | -0.0153 | 0.0084 | 3.796e3 | 19/30 | 1.019 | 1.022 |
| 30,000 | Q=Qmax, c=10 | 38 | 10 | 3.521 | 0.969 | 3.350 | -0.0261 | -0.0093 | 0.0067 | 2.605e3 | 13/42 | 0.699 | 0.702 |
| 30,000 | Q=Qmax, c=20 | 27 | 20 | 3.521 | 1.048 | 3.622 | +0.0069 | -0.0006 | 0.0073 | 3.796e3 | 19/30 | 1.019 | 1.022 |
| 100,000 | Q=20, c=10 | 20 | 10 | 1.760 | 3.517 | 6.211 | -0.0371 | +0.0075 | 0.0113 | 1.263e4 | 1/30 | 1.730 | 1.165 |
| 100,000 | Q=50, c=10 | 50 | 10 | 1.760 | 2.692 | 4.754 | +0.0054 | +0.0028 | 0.0093 | 5.277e3 | 5/66 | 0.723 | 0.487 |
| 100,000 | Q=Qmax, c=10 | 70 | 10 | 1.760 | 2.143 | 3.783 | -0.0015 | -0.0097 | 0.0081 | 4.287e3 | 53/78 | 0.587 | 0.396 |
| 100,000 | Q=Qmax, c=20 | 49 | 20 | 1.760 | 2.262 | 3.995 | -0.0005 | +0.0020 | 0.0086 | 5.277e3 | 5/66 | 0.723 | 0.487 |
| 300,000 | Q=20, c=10 | 20 | 10 | 3.521 | 2.111 | 7.355 | -0.0051 | -0.0084 | 0.0039 | 3.779e4 | 13/30 | 1.316 | 1.279 |
| 300,000 | Q=50, c=10 | 50 | 10 | 3.521 | 1.786 | 6.222 | -0.0188 | -0.0081 | 0.0033 | 1.531e4 | 5/66 | 0.533 | 0.518 |
| 300,000 | Q=Qmax, c=10 | 122 | 10 | 3.521 | 1.193 | 4.157 | -0.0022 | -0.0044 | 0.0026 | 9.184e3 | 0/1 | 0.320 | 0.311 |
| 300,000 | Q=Qmax, c=20 | 86 | 20 | 3.521 | 1.256 | 4.375 | -0.0095 | -0.0025 | 0.0027 | 9.938e3 | 23/102 | 0.346 | 0.336 |
| 1,000,000 | Q=20, c=10 | 20 | 10 | 1.760 | 4.886 | 8.583 | -0.0444 | +0.0004 | 0.0049 | 1.252e5 | 11/30 | 2.201 | 1.393 |
| 1,000,000 | Q=50, c=10 | 50 | 10 | 1.760 | 4.306 | 7.563 | +0.0042 | +0.0016 | 0.0043 | 5.170e4 | 49/66 | 0.909 | 0.575 |
| 1,000,000 | Q=Qmax, c=10 | 223 | 10 | 1.760 | 2.615 | 4.594 | +0.0139 | +0.0034 | 0.0031 | 3.213e4 | 0/1 | 0.565 | 0.358 |
| 1,000,000 | Q=Qmax, c=20 | 158 | 20 | 1.760 | 2.740 | 4.812 | +0.0105 | +0.0045 | 0.0033 | 2.237e4 | 149/210 | 0.393 | 0.249 |
| 1,531,530* | Q=20, c=10 | 20 | 10 | 5.463 | 1.651 | 9.015 | +0.2378 | -0.0011 | 0.0013 | 1.921e5 | 19/30 | 1.443 | 1.444 |
| 1,531,530* | Q=50, c=10 | 50 | 10 | 5.463 | 1.467 | 8.010 | +0.1181 | -0.0013 | 0.0012 | 7.766e4 | 47/66 | 0.583 | 0.584 |
| 1,531,530* | Q=Qmax, c=10 | 276 | 10 | 5.463 | 0.873 | 4.769 | +0.0292 | -0.0033 | 0.0008 | 4.638e4 | 0/1 | 0.348 | 0.349 |
| 1,531,530* | Q=Qmax, c=20 | 195 | 20 | 5.463 | 0.910 | 4.967 | +0.0339 | -0.0034 | 0.0009 | 3.399e4 | 43/210 | 0.255 | 0.256 |
| 2,000,000 | Q=20, c=10 | 20 | 10 | 1.760 | 5.269 | 9.284 | -0.0437 | +0.0010 | 0.0037 | 2.514e5 | 1/30 | 1.482 | 2.213 |
| 2,000,000 | Q=50, c=10 | 50 | 10 | 1.760 | 4.702 | 8.286 | +0.0028 | +0.0002 | 0.0033 | 1.014e5 | 47/66 | 0.598 | 0.893 |
| 2,000,000 | Q=Qmax, c=10 | 316 | 10 | 1.760 | 2.756 | 4.856 | +0.0095 | -0.0014 | 0.0023 | 5.868e4 | 0/1 | 0.346 | 0.517 |
| 2,000,000 | Q=Qmax, c=20 | 223 | 20 | 1.760 | 2.861 | 5.042 | +0.0049 | -0.0005 | 0.0024 | 3.373e4 | 0/1 | 0.199 | 0.297 |

Quadrature certificates per N:

| N | M | Parseval rel. err | R2 grid rel. err | R3 grid rel. err | R2/(SS2·N) | R3/main N−1 / N+1 |
|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | 2^17 | 1.8e-16 | 2.2e-16 | 0.0e+00 | 0.9574 | 1.00057 / 0.99261 |
| 30,000 | 2^18 | 0.0e+00 | 0.0e+00 | 3.5e-16 | 0.9816 | 0.99697 / 0.99762 |
| 100,000 | 2^20 | 0.0e+00 | 0.0e+00 | 3.7e-16 | 1.0030 | 0.99947 / 0.99940 |
| 300,000 | 2^22 | 0.0e+00 | 1.1e-16 | 0.0e+00 | 0.9895 | 0.99956 / 0.99903 |
| 1,000,000 | 2^23 | 0.0e+00 | 1.3e-16 | 6.7e-16 | 0.9978 | 1.00000 / 1.00006 |
| 2,000,000 | 2^24 | 0.0e+00 | 1.3e-16 | 3.2e-16 | 1.0010 | 1.00008 / 1.00000 |
| 1,531,530 | 2^24 | 1.8e-16 | 1.1e-16 | 1.8e-16 | 0.9994 | 0.99998 / 0.99991 |

* **Binary.** At fixed (20, 10), P_m/N rises from 3.29 (N = 10^4) to
  9.28 (N = 2·10^6), roughly linearly in log N; P_m/R2 zigzags because SS2(N) is 3.52
  at N = 3·10^4, 3·10^5 and 1.76 elsewhere. With Q = Qmax(N) (Q = 22 … 316, arcs covering
  ≈ 30% of the circle) P_m/N still climbs, 2.98 → 4.86: a rise
  of 1.9 against 6.0 at fixed (20, 10), while log N rises by 5.3. So Q ≈ √(N/2c) removes
  about 69% of the fixed-Q growth (slope 0.35 per unit log N against
  1.13; theory for Q ∝ √N: ½), more than ½ because arcs this wide also absorb background.
  P_m/R2 at Qmax is 1.77, 0.97, 2.14, 1.19, 2.62, 2.76 — it dips to
  0.97 < 1 at N = 3·10^4, where log N is small and SS2 is doubled. What has to be
  beaten is SS2(N), which never exceeds 5.46 for N ≤ 2·10^6 and grows at most like log log N, while P_m/N grows like
  log N, so the gap widens. The smooth probe shows the other side of the coin: at N = 1,531,530 (SS2 = 5.46)
  P_m/R2 = 1.65 (Q=20, c=10), 1.47 (Q=50, c=10), 0.87 (Q=Qmax, c=10), 0.91 (Q=Qmax, c=20). In
  3 tabulated rows the minor-arc mass is below the major-arc integral, so the crude inequality
  R2 ≥ I_M − P_m is already positive: N = 30,000 at (38, 10), I_M − P_m = +0.06·R2, N = 1,531,530 at (276, 10), I_M − P_m = +0.10·R2, N = 1,531,530 at (195, 20), I_M − P_m = +0.06·R2.
  All are N with an unusually large SS2. That says nothing about other N: typical N of the same size sit at
  P_m/R2 ≈ 2.6–2.9 at Qmax, and log N outruns log log N.
* **Ternary — this differs from the expected behaviour.** At fixed (20, 10), B3/main does **not** stay below 1:
  it is 1.38, 1.02, 1.73, 1.32, 2.20, 1.48 (N' = N−1) and creeps up with log N.
  The sup column shows why: with Q fixed, sup_m|S| is pinned at the q = 30 crest ψ/φ(30) ≈ N/8, a fixed fraction of
  N, so B3/main ≈ 2(log N − 1)/(8·SS3(N')) grows like log N, just like the binary ratio. What differs is the response
  to the arcs: grown as far as disjointness allows, the ternary ratio falls —
  1.38, 1.02, 0.72, 0.35, 0.39, 0.20 at (Qmax, 20) — because it scales like the tallest
  stray crest (1/φ(q*) for the first excluded q*, or 1/(π(c+½)) at the arc edge), while the binary ratio scales like
  log N − log Q. The asymptotic theory needs Q and c growing like powers of log N; the constraint 2cQ² < N leaves only
  Q ≤ √(N/2c), so at these N the decay cannot be displayed cleanly — only its mechanism.

## 8. Random control (Cramér model)

Main draw (seed 20261016): n ≥ 3 kept with probability 1/log n, weight log n; 78,482 members, ΣW = 998,067, ΣW² = 1.2788e7 (primes: 1.2803e7). Binary count R2_R(N) = 981,854.1 (exact; E R2_R = 1,000,007; grid integral agrees to 3.6e-16), i.e. singular series ≈ 1, against primes R2 = 1,756,568 ≈ SS2·N with SS2 = 1.7604.

| Q | c | I_M/R2_R | P_m/R2_R | I_m/R2_R | sup_m\|S_R\| | primes P_m/R2 |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 5 | 0.9943 | 12.029 | +0.0057 | 6.024e4 | 5.756 |
| 5 | 10 | 1.0053 | 12.017 | -0.0053 | 3.157e4 | 5.740 |
| 5 | 20 | 1.0101 | 12.010 | -0.0101 | 1.646e4 | 5.732 |
| 10 | 5 | 0.9945 | 12.027 | +0.0055 | 6.024e4 | 5.244 |
| 10 | 10 | 1.0052 | 12.012 | -0.0052 | 3.157e4 | 5.223 |
| 10 | 20 | 1.0096 | 12.000 | -0.0096 | 1.646e4 | 5.212 |
| 20 | 5 | 0.9948 | 12.015 | +0.0052 | 6.024e4 | 4.913 |
| 20 | 10 | 1.0049 | 11.989 | -0.0049 | 3.157e4 | 4.886 |
| 20 | 20 | 1.0095 | 11.954 | -0.0095 | 1.646e4 | 4.868 |
| 50 | 5 | 0.9957 | 11.937 | +0.0043 | 6.024e4 | 4.353 |
| 50 | 10 | 1.0068 | 11.834 | -0.0068 | 3.157e4 | 4.306 |
| 50 | 20 | 1.0125 | 11.641 | -0.0125 | 1.646e4 | 4.247 |
| 100 | 5 | 0.9930 | 11.668 | +0.0070 | 6.024e4 | 3.933 |
| 100 | 10 | 1.0044 | 11.291 | -0.0044 | 3.157e4 | 3.817 |
| 100 | 20 | 1.0117 | 10.555 | -0.0117 | 1.646e4 | 3.605 |
| 150 | 5 | 0.9916 | 11.209 | +0.0084 | 6.024e4 | 3.603 |
| 150 | 10 | 0.9990 | 10.376 | +0.0010 | 3.157e4 | 3.363 |
| 150 | 20 | 1.0063 | 8.728 | -0.0063 | 1.646e4 | 2.873 |

Five more independent draws (seeds 1–5) plus the main draw:

| seed | R2_R/N | P_m/R2_R (20, 10) | I_m/R2_R (20, 10) | P_m/R2_R (150, 20) | I_m/R2_R (150, 20) |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.9942 | 11.862 | +0.0047 | 8.623 | +0.0025 |
| 2 | 0.9770 | 12.020 | -0.0042 | 8.741 | -0.0122 |
| 3 | 1.0103 | 11.612 | +0.0309 | 8.453 | +0.0140 |
| 4 | 1.0311 | 11.434 | +0.0414 | 8.327 | +0.0181 |
| 5 | 0.9740 | 12.044 | -0.0056 | 8.768 | -0.0017 |
| 20261016 | 0.9819 | 11.989 | -0.0049 | 8.728 | -0.0063 |

Across the six draws R2_R/N = 0.974–1.031 and P_m/R2_R = 11.43–12.04
at (20, 10), 8.33–8.77 at (150, 20). The random set has the same L² wall, larger in
fact, because its energy ΣW² ≈ N log N is the same while nothing beyond the α = 0 arc is phase-locked and its R2 is
smaller (≈ N). So the wall belongs to the density 1/log n and to the method, not to arithmetic. The random model
escapes it — R2_R ≈ N ± O(√N log N) — by independence: its minor-arc phases really are random. The primes are
deterministic and offer no such independence; that is the gap. (sup_m|S_R| is always at the edge of the α = 0 arc:
there are no other crests.)

## 9. Goldbach's comet

G(N) for all even N ≤ 200,000 by one FFT autoconvolution of the prime indicator (spot checks against direct counting: True). G(N)/prediction does not scatter around 1 with no trend: its mean rises steadily toward 1 while the scatter shrinks — binned means 0.947 (1,000–5,000), 0.968 (5,000–10,000), 0.977 (10,000–20,000), 0.985 (20,000–50,000), 0.990 (50,000–100,000), 0.993 (100,000–150,000), 0.994 (150,000–200,000), and 5402 / 5417.1 = 0.9972 at N = 10^6. A slow approach to the asymptotic prediction, observed and not explained here.

## Where the oscillator picture is literal, and where it is only an analogy

* **Literal.** S(α) is exactly a sum of phasors Λ(n)e(nα). "Resultant", "channel", "constructive /
  destructive", "random walk" describe the arithmetic exactly, and the Ramanujan sum c_q(a) is literally
  the resultant of the φ(q) reduced-class phasors e(ar/q).
* **Not literal: there is no coupling.** The oscillators do not act on each other. Nothing in S(α)
  plays the role of a coupling strength K, nothing evolves in α, and there is no Adler/Kuramoto
  threshold |Δω| < K. "Phase-locking" at a/q means commensurability: the phase 2πna/q depends only on
  n mod q, so the oscillators fall into channels because of arithmetic (primes spread evenly over the
  reduced residue classes), not because they pull each other into step.
* **What plays the role of a locking window** is the coherence window |β| ≲ 1/N: detuned by β, the
  channel phases fan out over 2πNβ across n ≤ N and the resultant decays along the Dirichlet kernel. That
  is a finite-observation-length effect (the linewidth of a pulse train watched for N periods), not a
  balance between detuning and coupling.
* **Look-alikes.** The Farey forest resembles the Arnold-tongue diagram of a circle map, and F(α) in
  Fig F resembles a devil's staircase: structure at every rational, weaker as q grows. The mechanism is
  different — tongue widths come from coupling, spike heights here from μ(q)/φ(q).
* **"Beat"** in the captions means demodulation: S² is multiplied by the reference e(−Nα) and averaged
  over α, which keeps exactly the pairs n₁ + n₂ = N.

## Where the results differ from what the brief expected

1. **Ternary at fixed (Q, c) = (20, 10) does not stay below 1, nor flat** (§7): B3/main is
   1.38, 1.02, 1.73, 1.32, 2.20, 1.48 for N' = N − 1 and grows with log N, because
   sup_m|S| is pinned at the left-out 1/30 crest ≈ N/8. It falls below 1 once the arcs grow with N (Qmax series)
   — and only if c is large enough: at c = 5 the arc-edge sidelobe keeps B3/main(10^6 − 1) at 1.007 for every Q.
2. **P_m/R2 falls below the 3.5–7 band** at (150, 10) = 3.36, (150, 20) =
   2.87 and at the N = 10^6 Qmax pairs (158, 20), (223, 10) =
   2.74, 2.62: these arcs cover 14–30% of the
   circle and absorb background as well as spikes. It is still > 2.6.
3. **|I_m|/R2 ≈ 10⁻² only for Q ≥ 50.** At Q ≤ 20 it is 3–13%, and |I_m|/√P_m is 20–73 (also 15–17 at c = 5 with
   Q = 100, 150), above log N = 13.8: almost all of it is the deterministic leakage of left-out spikes and kernel
   tails (§5). What remains after removing it is the singular-series model's own error, not a clean measurement
   of random cancellation.
4. **Binary P_m/R2 is not > 1 for every tabulated N and (Q, c)**: it is 0.97 at N = 3·10^4
   with Q = Qmax = 38, and 0.87–0.91 at the smooth
   probe N = 1,531,530 with Qmax — both N with an unusually large singular series. For N = 10^6 every tabulated
   pair gives P_m/R2 ≥ 2.62.
5. **B3/main at (150, 20)** is 0.39 for N' = 10^6 − 1 (inside the expected 0.3–0.6) but
   0.25 for N' = 10^6 + 1, whose singular series is larger (2.30 vs 1.46).
6. **G(N)/prediction** sits below 1 with a slowly rising mean (0.947 at 10³–5·10³ up to 0.994
   at 1.5–2·10^5, 0.9972 at 10^6) rather than scattering around 1 with no trend.
7. **Memory**: peak RSS 1,432 MB (the N = 2·10^6 and 1,531,530 runs use M = 2^24).

## Figure captions

**Fig A (fig_A_farey_forest.png)** — Each prime power n ≤ 10^6 is an oscillator with amplitude Λ(n) and phase 2πnα; |S(α)| is the length of their resultant. At a rational a/q every phase depends only on n mod q, so the oscillators lock into the φ(q) reduced channels, each carrying ≈ ψ/φ(q), and the channels add up to the Ramanujan sum: a spike of height |μ(q)|ψ/φ(q) (orange ticks), which is where the measured crests (blue) stand. For non-squarefree q (red ×) the locked channel phasors sum to zero and only a residue of 10²–10³ survives, below the sea. Between the spikes the bin medians (≈ 2,073) are close to the Rayleigh median √(ln2·⟨|S|²⟩) = 2,040 of a sea with mean energy ⟨|S|²⟩ ≈ 6.0e+06 (Fig E window), about half of ΣΛ² = 1.28e+07: Parseval's energy is shared between forest and sea, and around 1/2, 1/4, 3/4 the forest thins because a small-q rational has no low-q neighbours nearby (Farey gaps).

**Fig B (fig_B_spike_heights.png)** — (a) S(a/q) for every reduced a/q with q ≤ 30, summed exactly over residue classes, S(a/q) = Σ_r e(ar/q)·ψ(N;q,r). For squarefree q the length matches |μ(q)|ψ/φ(q) to within 2.7% (worst case 27/29); non-squarefree q sit at ≤ 1.7e-03·ψ, below even the random-walk level √ΣΛ²/ψ = 3.6e-03, because their locked channel phasors cancel. (b) The resultant also points the predicted way: φ(q)·S(a/q)/ψ lies on the real axis at μ(q), i.e. the φ(q) channels add to exactly one channel, forwards or backwards according to the parity of the number of prime factors of q. The residue at 1/4 is literally the mod-4 race: S(1/4) = [ψ(N;4,0) − ψ(N;4,2)] + i[ψ(N;4,1) − ψ(N;4,3)] = 11.8 -835.8i.

**Fig C (fig_C_zooms.png)** — Zooms of half-width 20/N around 1/3, 1/6 and 1/4 (direct summation, 60 samples per 1/N). Detuning by β makes oscillator n precess at rate nβ, so a locked channel's phases fan out over 2πNβ; the collective phasor decoheres along the Dirichlet kernel v(β), whose zeros at Nβ = ±1, ±2, … are the detunings at which the fan closes into whole turns. The model (μ(q)/φ(q))·v(β) (orange) follows the measured modulus through ~20 sidelobes on each side (at ±20/N they are still 2–3× the sea level √ΣΛ²), and in the frame co-rotating with the carrier e((N+1)β/2) the phase flips between 0 and π at every zero, starting at π for 1/3 (μ = −1) and at 0 for 1/6 (μ = +1). At 1/4 the model is identically zero; what is left is sea plus the mod-4 residue, and its phase wanders.

**Fig D (fig_D_golden_zoom.png)** — Same width and scale, centred on frac(golden ratio) = 0.6180339…, the number hardest to approximate by rationals. Even here the convergents and their neighbours with q ≈ √N carry μ(q)-locked crests ψ/φ(q) as tall as the sea: 2 stand above √ΣΛ² = 3,578, the tallest 144/233 at 4,309. The window maximum |S| = 5,506 sits on that crest, where exactly S(144/233) = -5,397 -206i points backwards as μ(233) = -1 predicts, and the exact sums at the other sizeable crests point along μ(q) too (19 of 20). Between the crests the co-rotating phase drifts freely: at N = 10^6 even the deepest minor arc is weak phase-locking at the level of the sea, not pure noise.

**Fig P (fig_P_phasors.png)** — Λ-weighted histograms of the phases 2πnα (360 bins) with the resultant S(α) drawn to the same radial scale. At 1/3 the primes fill two channels (n ≡ ±1 mod 3) of ψ/2 each at ±120°, and at 1/6 two channels at ±60°: in both cases the resultant has length ψ/2 (-499,774 and +499,760), pointing backwards at 1/3 (μ(3) = −1) and forwards at 1/6 (μ(6) = +1). At 1/4 the channels meet head-on at 90° and 270° and leave 836 ≈ 8e-04·ψ, the imbalance of the mod-4 race. 1/5 and 1/30 have 4 and 8 channels whose sum is exactly one channel reversed (-249,800 ≈ −ψ/4, -124,627 ≈ −ψ/8). At the irrationals the bins form a flat ring and the resultant is random-walk sized (2,856 and 803 against √ΣΛ² = 3,578).

**Fig S (fig_S_singular_series.png)** — SS2_Q(10^6) against Q. Each squarefree q adds the coherent beat of two locked channel sets, μ(q)²c_q(N)/φ(q)², a product over its primes of c_p(N)/(p−1)²: +1 for p = 2 and +1/4 for p = 5 (they divide N), −1/(p−1)² for every other p — so the kick is positive when q has an even number of prime factors outside {2, 5} (q = 21 adds +1/144) and negative when it has an odd number. The partial sums swing — 1, 2, 1.75, 2, 1.75, 1.72, … — and settle: SS2_20 = 1.8580, SS2_150 = 1.75900, and the sum to Q = 10^6 equals 2C2·4/3 (C2 to 20 digits) = 1.7604315 to 6e-10. Panel (b): the signed error changes sign only at Q = 2, 3, 5, 6, 10, 130, then stays negative and shrinks from ~10⁻³ to ~10⁻⁸; beyond Q ≈ 19,914 the ~10⁻⁸ tail oscillates about 0.

**Fig E (fig_E_minor_stretch.png)** — A stretch of width 300/N at α₀ = 0.300216, inside the minor arcs of every (Q, c) tabulated at N = 10^6. (a) The mass density |S|² is mostly speckle — coherence cells of width ~1/N, window mean 6.00e+06 — and (b) the integrand Re[S²e(−Nα)] has the same envelope, but its phase (the beat of S² against the reference e(Nα)) turns once per ~1/N, so neighbouring cells enter with effectively random signs. (c) Mass piles up linearly to 1,801 over the window while the signed integral wanders and ends at +18: a random walk of ~300 steps. The fractions with q ≤ 300 inside the window (88/293, 85/283) have model crests only 1.4× the sea's rms √⟨|S|²⟩ = 2,450, and the two highest peaks (5.7× and 5.6× the mean) are speckle, ≥ 134/N from any of them; still, the exact S(a/q) point along μ(q) in 2 of 2 cases — weakly locked crests that the major arcs leave on the minor arcs.

**Fig F (fig_F_staircase.png)** — The running integral F(α) = ∫₀^α S²e(−Nt)dt. Each major-arc spike is a burst of phase-locked beat between S² and e(Nα) that adds the phasor (N μ(q)²/φ(q)²)·e(−aN/q): real jump N μ(q)²cos(2πaN/q)/φ(q)², so +N at 1/2, −N/8 at 1/3, 2/3, 1/6 and 5/6, +N/16 at each of 1/5…4/5, nothing at 1/4, 3/4, and smaller steps at q = 7, 10, … by the same rule. Measured jumps over ±20/N agree to 5.8e-03·N (the kernel tail beyond 20/N), and between jumps the plateaus wander at the 10⁻³N level (inset: a 1.5e-03·N range over α ∈ [0.777, 0.782]) — small spikes plus minor-arc sea. The imaginary jumps of a/q and (q−a)/q cancel in pairs, so Im F is symmetric about 1/2 and returns to -6.7e-11; the real part ends at R2(N)/N = 1.756568, exactly R2.

**Fig T (fig_T_binary_vs_ternary.png)** — (a) Magnitude-only minor-arc bounds divided by the main term, for every grid pair (Q, c) at N = 10^6. Binary: P_m/R2 = 2.9–5.8, always above 1 — the decoherent energy left on the minor arcs exceeds the whole signal. Ternary: sup_m|S|·∫|S|² ÷ SS3·N'²/2 drops below 1 only when both kinds of stray crest are short enough — the excluded squarefree q (needs Q) and the arc-edge sidelobe ψ/(π(c+½)) (needs c; at c = 5 it keeps N' = 10^6−1 at 1.007 for every Q) — reaching 0.39 (N' = 10^6−1, solid) and 0.25 (N' = 10^6+1, hollow) at (150, 20). (b) What the minor arcs really contribute: binary |I_m|/R2 of 3–13% at Q ≤ 20, almost all predictable leakage of left-out spikes; what is left after removing it is the model's own error, comparable in size to the random-phase yardstick P_m/√(N|m|) (grey) — the phases cancel to a fraction of a percent of a mass several times larger than the signal.

**Fig G (fig_G_scaling.png)** — N from 10^4 to 2·10^6, four (Q, c) series, no fits; stars mark the smooth-N probe N = 1,531,530 (SS2 = 5.46). (a, b) At fixed (Q, c) the binary minor-arc energy P_m/N climbs roughly linearly in log N, and even Q = Qmax(N) (the most arcs 2cQ² < N allows) only lowers the slope; the singular series it has to beat (black: its largest value at each size) barely moves. (c) The ternary ratio at fixed (Q, c) also creeps up with log N, because the sup is pinned by the left-out crest at 1/30 (≈ N/8) — contrary to the brief's expectation that it stays flat — while with Q = Qmax it falls, from 0.9–1.4 at 10^4 to 0.2–0.5 at 2·10^6. The ternary bound only needs the strongest stray resultant on the minor arcs to be small, which more and wider arcs can arrange; the binary bound needs their total energy to be small, which they cannot.

**Fig H (fig_H_random_model.png)** — (a) Max |S| per bin for the primes (blue) and, drawn on top, a Cramér random set with the same density (orange): randomising destroys every phase-locking except the trivial one at α = 0, so only the primes' forest sticks out above the random band. (b) The median per bin: the random sea sits uniformly higher, because it keeps all of its energy ΣW² instead of lending about half to spikes. (c) At the rationals themselves the random resultant is of random-walk size. Its binary count is R2_R = 981,854 ≈ N (singular series 1) against 1,756,568 ≈ 1.76N for the primes, yet its minor-arc energy is P_m/R2_R = 8.7–12.0 (table in §8): the same L² wall, higher — the density's wall, not arithmetic's, which the random set gets past by an independence the primes do not provide.

**Fig K (fig_K_comet.png)** — G(N) for every even N ≤ 2·10^5. (a) The comet splits into bands by N mod 3: if 3 | N, the two channels n ≡ 1, 2 (mod 3) pair up with each other in both orders, doubling the count (factor (3−1)/(3−2) = 2), while N ≡ 1 and N ≡ 2 share the lower bands. (b) Coloured by log SS2(N), each band is one value of the singular series — one pattern of phase-locked channels mod the small primes dividing N — and the Hardy–Littlewood prediction (orange) runs through each band. (c) Observed over predicted is not trend-free: its binned mean (orange) rises from 0.947 at N = 10³–5·10³ to 0.994 at 1.5–2·10^5 (0.9972 at 10^6) while the scatter shrinks — a slow approach to the prediction.

## What this does and does not show

**Does show (for N = 10^6 and the listed N, in floating point):**

1. The circle-method identities hold numerically to machine precision: Parseval, R2 = ∫S²e(−Nα) and R3 = ∫S³e(−N'α)
   on a uniform grid with M > 2N' (relative errors ≤ 1e-15), with the FFT sign pinned by a check that fails loudly.
2. The major-arc picture is quantitatively right at this N: spike heights |μ(q)|ψ/φ(q) to ≤ 2.7% for q ≤ 30,
   the Dirichlet-kernel shape and phase through ~20 sidelobes, and the major-arc integral within a few percent of
   R2 for Q ≥ 20, c ≥ 10 (≤ 1.5% for Q ≥ 50). The singular series predicts R2(10^6) to 0.22%
   and R3(N') to ≤ 7×10⁻⁵.
3. The binary obstruction, measured: at N = 10^6 the minor arcs carry 2.6–5.8× the signal in L²
   mass for every tabulated (Q, c), while their signed contribution is 1e-04–2e-02 of R2 for
   Q ≥ 50 and is accounted for, to within the major-arc model's own error (≤ 0.7% of R2),
   by the predictable leakage. A magnitude estimate cannot see that cancellation. P_m/N grows like log N, while
   SS2(N), which it must stay below for the L² bound to close, grows at most like log log N.
4. The ternary contrast: with S³ one factor can be bounded by its sup, and at (Q, c) = (150, 20) the bound
   sup_m|S|·∫|S|² is 0.25–0.39 of the main term. It closes because Q absorbs the
   tall excluded crests and c pushes the arc-edge sidelobes down; at fixed small Q, or at c = 5 for N' = 10^6 − 1,
   it does not.
5. The Cramér model (six draws) has the same L² wall but no Farey forest and singular series 1: the wall comes from
   density and method, the forest from arithmetic.

**Does not show:**

1. Anything asymptotic. N = 10^6 is pre-asymptotic: the textbook arc width Q/N with Q = (log N)^B (B ≥ 2) would make
   2Q³ > N, disjointness caps Q at √(N/2c), and all results are tabulated over (Q, c) rather than at one canonical
   choice. No trend has been fitted to minor-arc data; the heuristic curves contain no fitted parameters.
2. That the minor-arc cancellation persists, or how large it is beyond the model error. The residual after leakage
   is identically the Hardy–Littlewood error minus the major-arc model error; it is ≤ 4.4e-03·R2 at
   N = 10^6 (Q ≥ 50) and up to 1.5e-02·R2 at N ≤ 3·10^5. Nothing here bounds it for any other N, and that
   is precisely what would be needed.
3. Rigorous bounds. sup_m|S| is a grid maximum refined locally by direct summation, not a certified supremum;
   arithmetic is IEEE double precision, not interval arithmetic.
4. A proof strategy, binary or ternary. The ternary closure is shown for two specific N'; the actual theorem needs
   sup_m|S| ≪ N/(log N)^A uniformly (Vinogradov, with Siegel–Walfisz on the major arcs), which no finite computation
   displays.
5. An explanation of the 1–5% deficit of G(N) against the Hardy–Littlewood integral at N ≤ 2·10^5, shrinking with N.
