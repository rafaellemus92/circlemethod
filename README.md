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
| `escapes.py`, `run_escapes.py` | the escape-route experiments E0–E5 (below); `run_escapes.py` writes `cache/escapes.pkl`, `figures/fig_X*.png`, `escapes.md` |
| `figures_escapes.py`, `escapes_report.py` | figures X0–X5 and the `escapes.md` writer (every number in it is read from the pickle) |
| `tests/test_escapes.py` | pytest on a small grid (X = 10^4): Bessel identity in N, 4th moment = additive energy, transmission certificate, G and R_K grid = direct, theorem-shaped bound ≥ joint bound, almost-prime and smooth counts, the z ≤ 7 closed form, η constants, window bound monotone |
| `escapes.md` | the escape-route results: certificates, one section per experiment, captions, "What this does and does not show", references |
| `fractals.py`, `run_fractals.py` | the fractal lens M1–M5, A1–A3 (below): estimators (D_s, f(α), Hill, structure functions), no-fit models (Beckmann, Rice-FG, Farey staircase, spike-sum tail, fourth moment, triangular envelope), the singular-series surrogate H_R, Mandelbrot-set bulbs; `run_fractals.py` writes `cache/fractals.pkl`, `figures/fig_M*.png`, `fractals.md` |
| `figures_fractals.py`, `fractals_report.py` | figures M1–M6 and the `fractals.md` writer (every number in it is read from the pickle) |
| `tests/test_fractals.py` | pytest on X = 10^4: estimators on uniform and cascade measures, Farey count and staircase, sinc integrals, all-a FFT, Beckmann vs Rice, FG and q*, collapse, spike CCDF asymptote, Hill, sup prediction, structure functions, triangular envelope, spectral identity, surrogate identities, uncertainty products, bulb centres, size law, mediant counts |
| `fractals.md` | the fractal-lens results: certificates, one section per block, the physics dictionary, literal versus analogy, synthesis, captions, references |

## Running

Needs Python 3.11 or newer (tested on 3.11 and 3.13), about 2 GB of free RAM and ~1 GB of disk for `cache/`.

```bash
pip install numpy scipy matplotlib pytest
python -m pytest -q tests          # ~30 s
python run_all.py                  # ~5 min on 4 cores, peak RSS ~1.4 GB; writes figures/ and results.md
python run_all.py --figures        # redraw figures and results.md from cache/ctx.pkl
python run_escapes.py              # escape routes E0–E5: ~4 min, peak RSS ~1.4 GB; needs run_all.py's cache
python run_escapes.py --report     # redraw figures/fig_X*.png and escapes.md from cache/escapes.pkl
python run_fractals.py             # fractal lens M1–M5, A1–A3: ~10 min, peak RSS ~1.6 GB; needs run_all.py's and run_escapes.py's caches
python run_fractals.py --report    # redraw figures/fig_M*.png and fractals.md from cache/fractals.pkl
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

## Escape routes: what gets through the wall, and what that does not mean (`escapes.md`)

The main run shows the obstacle: at N = 10^6 the minor arcs carry 2.6–5.8 times more energy than the whole count, so
any bound that sees only |S| there fails (as it does for a typical N; E5 below finds the rare exceptions). Getting past that for every large N means proving that the minor-arc phases
cancel, which is the open problem; nothing here attempts it. `escapes.md` instead measures the obstacle for one N and
tests, at the same X = 10^6, the classical relaxations under which the *magnitude* inequality does close:

* **E0, the obstacle.** The 6.1 M minor-arc phasors at (150,20) form a walk 2.9× the signal long that ends 1% of the way
  out; every magnitude-only refinement is longer still. Phase-summed windows ≈ 21/N wide would suffice at this N (W < 0.99).
* **E1, average over N.** One energy budget Σ_N I_m(N)² = ∫_m|S|⁴ is shared by all targets, so at most 274 of the
  250,001 even N in [X/2, X] can fail at (150,20) (15,058 with the textbook sup form sup²·ΣΛ² of the 1937–38 proofs);
  none does.
* **E5, single-N transmission.** With no relaxation, 21 smooth N close to X (SS2 ≈ 5) have P_m < I_M(N) at (223,10),
  16 of them below the 0.99 margin; the share of targets with P_m < I_M(N) falls at every tabulated X, from 2.8% at
  X = 10^4 to 0.0014% at 2·10^6.
* **E2, Linnik (add powers of 2).** One power of 2 already closes the joint magnitude bound at Q ≥ 20 — generic
  averaging, as a random-shift control shows; the substitutes proofs can use push K to 4–7 at (150,20).
* **E3, almost-primes.** Replacing one prime by m with Ω(m) ≤ 2 closes at (150,20), Ω ≤ 3–4 elsewhere; random sets of
  the same density never close. Chen's and Rényi's theorems are sieve theorems; r = 2 here is not Chen's 2.
* **E4, smooth weight (control).** Suppresses the arc-edge sidelobes (to 0.9% of the peak just past the c = 5 edge,
  against 5.7% for the sharp cutoff), which removes the 1.007 plateau of the ternary bound, and leaves the binary wall
  where it was.

Every K\*, r\* and z\* there is a finite-X statistic, with a decision margin (closes < 0.99 < undecided < 1.01 < fails),
not a constant comparable to a theorem's.

## The fractal lens: Mandelbrot's principles as measurements (`fractals.md`)

`fractals.md` applies three of Mandelbrot's principles to the same objects at X = 10^6: scaling only between cutoffs
(self-similarity is not self-affinity), roughness H versus intermittency D_s and f(α), and mild versus wild randomness
(Noah versus Joseph). Every model is no-fit, every exponent is quoted with its window, and every result is checked
against a surrogate made of the singular-series spikes q ≤ 1000 plus a Gaussian sea of matched power.

* **One fact read five ways.** The Farey count #{a/q : μ(q) ≠ 0, φ(q) ≤ y} ~ y²/2 (constant exactly ½), applied to
  spikes of height ψ/φ(q) and width 1/N, gives the minor-arc tail of index 1 in |S|², the coarse-grained straight-line
  spectrum f(α) = α (offset −ln(2L)/ln(1/ε), L = 10.82), the large-value count exponent 2(1 − σ), the kurtosis, and equal
  energy per octave of q (the L² wall). These are not independent confirmations.
* **M1.** The spikes are a data collapse (heights scaled by μ(q)/φ(q), widths fixed at 1/N), not self-similar. The
  collapse ends at q* = 521 (513 predicted, with no fit, from the variance of primes in progressions), which is not the
  arc range √(N/2c).
* **M3.** Minor arcs at Q ≤ 20 carry a tail that is wild in the untruncated idealisation; at Q ≳ 50 they are neither
  wild nor exponential (residual spikes on a sea).
* **M4.** The minor-arc walk is smooth below 1/N (generic), then diffusive with H₂ ≈ 0.52–0.53; heavy tails (Noah) show
  in H_½ ≈ 0.78 at (20,10); no long memory.
* **M5.** The Mandelbrot set's bulbs share the Farey index set, circle order and a q^(−2) size law with the spikes
  (through Diophantine approximation, not dynamics); the mediant rule holds for bulbs in 140/140 neighbour pairs of F_21
  and for spikes in 56/140.
* A physics dictionary sorts the vocabulary (phasors, conservation and leakage, uncertainty principle, diffraction,
  equipartition: literal; spin, time, fluids: not).
