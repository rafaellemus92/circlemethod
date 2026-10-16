# goldbach_circle — the Hardy–Littlewood circle method for binary Goldbach at N = 10^6, made visible

> **This is an exploration and visualisation project, not a proof and not an attempt at one.**
> It computes, at one finite N (and a handful of others), every object the circle method talks about —
> the prime exponential sum on the whole circle, its spikes at rationals, the major-arc model, the
> singular series, the minor-arc integrand — so that the reason the method proves the *ternary*
> problem but not the *binary* one can be seen and measured. N = 10^6 is far from asymptotic; the
> results are tabulated over the arc parameters (Q, c) instead of hiding that sensitivity.

## The picture in one paragraph (phasor language)

Each n ≤ N with Λ(n) ≠ 0 is an oscillator of amplitude Λ(n) whose phase at "frequency" α is 2πnα;
S(α) = Σ Λ(n) e(nα) is the resultant. At a rational a/q the phase depends only on n mod q, so the
oscillators **phase-lock** into the φ(q) reduced residue channels, each carrying ≈ ψ(N)/φ(q); the
channels add up to the Ramanujan sum, giving a spike of height |μ(q)|ψ/φ(q) (Farey forest). Detuning
by β makes the channels precess and the resultant decoheres along the Dirichlet kernel (width 1/N).
Away from rationals with small q the phases are scattered and the resultant is a random walk of size
√ΣΛ² (minor arcs). The binary count R2(N) = ∫ S(α)² e(−Nα) dα is the beat of S² against the reference
e(Nα): the major arcs give the singular-series prediction, while the minor arcs carry **2.6–6 times
more energy than the whole signal** and yet contribute only 10⁻⁴–2·10⁻² of it (for Q ≥ 50), because the
phases of that energy cancel. A bound that sees only magnitudes cannot see that cancellation (the L² wall).
With S³ (ternary) one factor can be bounded by its supremum, and that bound closes once the arcs swallow
the tall near-rational spikes (Q large) and push the arc-edge sidelobes down (c large).

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

## Layout

| file | contents |
|---|---|
| `sieve.py` | module 0: sieve, Λ(n), μ, φ, ψ(N), ΣΛ², R2(N) by direct summation, G(N), FFT autoconvolution (R3, the comet), residue-class sums ψ(N; q, r), the Cramér random set |
| `expsum.py` | module 1: S on the grid by FFT (sign convention documented in the module docstring), cached as `.npz`; direct and exact-rational evaluation of S; Dirichlet kernel; grid integrals; FFT sign checks |
| `arcs.py` | major/minor arcs on the grid (exact integer membership, disjointness asserted), major-arc model, arc integrals, refined sup over the minor arcs, the per-(N, Q, c) table rows |
| `singular_series.py` | Ramanujan sums (definition, divisor formula, Hölder), SS2 and SS3 (truncated sums and Euler products), C2, C3, the Hardy–Littlewood integral |
| `figures.py` | all figures (A–F as named in the brief, plus P phasors, S singular series, T binary vs ternary, G scaling, H random model, K comet) |
| `run_all.py` | runs modules 0–9, writes `figures/*.png` and `results.md` |
| `tests/` | pytest: π(10^6); G(10^6); Parseval on the grid; grid integral = R2; FFT sign; R3 grid = convolution; c_q vs divisor formula and multiplicativity; SS2 sum vs product to 1e-4; SS3 likewise; disjointness assertion; major + minor = whole |
| `results.md` | headline table, per-module tables, a caption per figure, deviations from the brief, and "What this does and does not show" |

## Running

```bash
pip install numpy scipy matplotlib pytest
python -m pytest -q tests          # ~10 s
python run_all.py                  # ~5 min on 4 cores, peak RSS ~1.4 GB; writes figures/ and results.md
python run_all.py --figures        # redraw figures and results.md from cache/ctx.pkl
```

`cache/` holds S on each grid (`S_<set>_N<N>_M<M>.npz`, 134 MB for M = 2^23) and the pickled results.

## Definitions and conventions

* e(x) = exp(2πix); Λ = von Mangoldt, prime powers included (their total weight is
  Σ_{k≥2} θ(N^{1/k}) ≈ 1.1×10³ at N = 10^6, i.e. O(√N)).
* S(α) = Σ_{n≤N} Λ(n) e(nα); ψ(N) = S(0); R2, R3, G, c_q, v(β), SS2, SS3 exactly as in the brief.
* **Grid.** α_k = k/M, k = 0..M−1, M = 2^23 for N = 10^6 (M = smallest power of two ≥ 8N for the
  other N). One quadrature rule everywhere: the rectangle rule (1/M)Σ_k. On the full circle it is
  exact for Parseval and R2 (needs M > N) and for R3 (needs M > 2N'), so those integrals are
  certificates of the quadrature, not approximations.
* **FFT sign convention.** numpy's `fft` is X[k] = Σ x[n] exp(−2πi nk/M). Hence
  S(k/M) = M·ifft(Λ)[k] = conj(fft(Λ)[k]); we compute it with the half-length `rfft`. The checks
  S(0) = ψ, S(1/2) = 2⌊log₂N⌋log 2 − ψ and S(−α) = conj S(α) hold for either sign and so do not decide
  it; what decides it is (a) Im S(1/M) against a direct sum and (b) the binary grid integral, which
  returns R2(N) exactly with the right sign and ≈ 2×10⁻¹⁰ (zero to rounding, against R2 ≈ 1.8×10⁶) with
  the wrong one, because the wrong sign evaluates the autoconvolution at M − N > 2N.
* **Arcs.** M(Q, c) = ∪_{q≤Q, (a,q)=1} {‖α − a/q‖ ≤ c/N}; a grid point is in the arc iff
  |kqN − aMN| ≤ cMq (integers). `MajorArcs` asserts 2cQ² < N (sufficient for disjointness, as the
  brief asks; the sharp condition is 2cQ(Q−1) < N) and re-checks disjointness on the grid. Width c/N
  rather than the textbook Q/N, because at N = 10^6 the asymptotic choice Q = (log N)^B with B ≥ 2 would
  make 2Q³ > N; the kernel tail outside c/N carries ≈ 1/(π²c) of each spike.
* **sup over the minor arcs** = grid maximum, then the top candidates are re-evaluated by direct
  summation on a ±1/M mesh restricted to the minor arcs (the grid alone can miss a crest by ≤ 0.6% on a
  main lobe and ≤ 1.8% on a sidelobe).

## What to look at first

`results.md` → the headline table and §5 (the minor arcs), then `figures/fig_E_minor_stretch.png`
(mass vs integrand vs running integral) and `figures/fig_T_binary_vs_ternary.png`.
