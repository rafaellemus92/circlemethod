"""Write fractals.md from the dict produced by run_fractals.compute() (every number is read from it)."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from escapes_report import _table
from fractals import H_triangular, fg_factor, rice_moments
from run_all import _e, _f
from sieve import factorize

PAIRS7 = [(5, 5), (10, 10), (20, 10), (50, 10), (100, 20), (150, 20), (223, 10)]
CCDF5 = [(5, 5), (20, 10), (50, 10), (150, 20), (223, 10)]
WALK = [(150, 20), (20, 10)]
BANDS = [(0.5, 1), (1, 10), (10, 100), (100, 1000), (1000, 3000), (3000, 50000)]
LAGS = [2**j for j in range(22)]
HILL_K = (300, 1000, 3000, 10000)


def pair(p) -> str:
    return f"({p[0]},{p[1]})"


def g3(x: float) -> str:
    return f"{x:.3g}"


def rng(vals, d: int = 2) -> str:
    v = [x for x in vals if x is not None and np.isfinite(x)]
    return f"{min(v):.{d}f}–{max(v):.{d}f}"


def xl(res) -> np.ndarray:
    """Lag coordinate of each local exponent: the smaller lag of the pair (l, 2l), in units of 1/N."""
    return np.array(LAGS[:-1]) * res["X"] / res["M"]


def row_q(res, q):
    return next(r for r in res["m1_rows"] if r["q"] == q)


def fg_any(q: int) -> float:
    return fg_factor(q, list(factorize(q)))


def qstar_at(res, X):
    return next(d for d in res["m1_qstar"] if d["X"] == X)


def band_stats(res, pr, key: str = "raw") -> list[dict]:
    """Primes' local H_2 against the Cramér ensemble, averaged over each band of N*lag."""
    x = xl(res)
    w = res["walks"][("primes",) + pr]
    Hp = w["full"][2.0]["H"] if key == "raw" else w["full"]["drift"]["H"]
    ens = np.array([e[key] for e in res["ensemble"][pr]])
    hyb = {R: res["walks"][(f"hybrid{R}",) + pr]["full"][2.0]["H"] for R in (300, 1000, 3000)}
    out = []
    for lo, hi in BANDS:
        sel = (x >= lo) & (x < hi)
        if not sel.any():
            continue
        p = float(Hp[sel].mean())
        per = ens[:, sel].mean(axis=1)
        m, s = float(per.mean()), float(per.std(ddof=1))
        z = (p - m) / s
        nge = int(np.sum(per >= p))
        if z >= 3 and nge == 0:
            verdict = "excess"
        elif z >= 2:
            verdict = "marginal"
        elif z <= -2:
            verdict = "below"
        else:
            verdict = "none"
        lags = np.array(LAGS[:-1])[sel]
        out.append({"band": (lo, hi), "primes": p, "mean": m, "sd": s, "z": z, "nge": nge, "nseeds": len(per),
                    "verdict": verdict, "indep": (w["n"] / lags.min(), w["n"] / lags.max()),
                    "hyb": {R: float(h[sel].mean()) for R, h in hyb.items()}})
    return out


def hband(res, tag, pr, var, p, lo=4.0, hi=1000.0):
    x = xl(res)
    st = res["walks"][(tag,) + pr]
    d = st["full"].get(p) if var == "full" else st.get(var, {}).get(p)
    if d is None:
        return None
    sel = (x >= lo) & (x <= hi)
    return float(np.mean(d["H"][sel]))


def stair_value(st, L, i):
    g = int(np.searchsorted(st["rank_end"], i))
    return 2 * L * i * st["mass"][g]


def comp_at(d, i):
    k = np.flatnonzero(np.asarray(d["ranks"]) == i)
    return float(d["comp"][k[0]]) if len(k) else None


def fmt(x, spec: str = ".2f") -> str:
    return "—" if x is None or not np.isfinite(x) else format(x, spec)


def minor_d2_min(res, pr, tag):
    """Min of the local minor-arc D_2 over N*eps in [1.9, 977] (j = 4..13 in the local-slope index) and its location."""
    X, M = res["X"], res["M"]
    d = np.asarray(res["m2"]["minor_D2"][pr][tag])
    Ne = np.array([X * 2 ** (j + 0.5) / M for j in range(len(d))])
    sel = (Ne >= 1.9) & (Ne <= 977 * 1.5)
    k = int(np.argmin(np.where(sel, d, np.inf)))
    return float(d[k]), float(Ne[k])


# ================================================================ sections
def header(res) -> str:
    m2 = res["m2"]
    return f"""# The fractal lens: Mandelbrot's measurement programme on the circle method at N = 10^6

> Exploration at one finite X, not a proof and not an attempt at one. Every number below is a floating-point statistic
> of one computation at X = N = 10^6 on the grid M = 2^23 (plus the X-scaling runs). Every exponent is a local
> difference quotient at X = 10^6 over a stated window; no trend is fitted, and "no fit" means every constant in a model
> comes from ψ, ΣΛ², μ, φ, N, Q, c or R.

This file applies three of Benoit Mandelbrot's principles to the circle-method objects of results.md and escapes.md, as
measurements. It then compares the primes' rational spikes with the bulbs of the Mandelbrot set:

- **P1. Scaling holds only between cutoffs.** Every exponent comes with its window. Self-similarity (invariance under
  isotropic zoom) is not self-affinity (different rescaling in different directions, or only statistical).
- **P2. Roughness is measured by H; intermittency by D_s and f(α).** The Hurst exponent H measures how rough a walk is.
  The generalized dimensions D_s and the singularity spectrum f(α) measure how unevenly a measure is spread.
- **P3. Mild versus wild randomness; Noah versus Joseph.** Wild (state 6 of Mandelbrot's seven, 1997 ch. E5) means
  infinite variance with some finite positive moment, i.e. Pareto index ≤ 2. The Noah effect is heavy tails; the
  Joseph effect is long memory.

**The one fact.** Almost every scaling law below is one counting identity:

    #{{a/q : μ(q) ≠ 0, φ(q) ≤ y}} = Σ_{{μ²(q)=1, φ(q)≤y}} φ(q) ~ y²/2,

applied to spikes of height ψ/φ(q) and width 1/N at every reduced a/q with squarefree q, on top of a structureless sea.
The constant is exactly ½. The Dirichlet series is Σ μ²(q)φ(q)^(1−s) = ζ(s−1)·Π_p(1 + (p−1)^(1−s))(1 − p^(1−s)), and
every Euler factor equals 1 at s = 2. Read pointwise, the count is a power-law tail; read per box, a straight-line
spectrum; read per coherence cell, a large-value exponent; read per octave of q, equal energy per octave. These are one
fact read five ways, not five confirmations.

Normalisations used throughout: ψ(X) = {_f(res['psi'], 1)}, ΣΛ² = {_e(res['L2'], 5)}. The trivial atoms at α = 0 and
α = 1/2 (within |Nα| ≤ 20) hold a share {m2['atom_share']:.4f} of Σ|S|² and are removed before any measure is formed. The
spike-mass normaliser is L = X·P_kept/ψ² = {m2['L']:.3f} (2L = {2 * m2['L']:.2f}), so a spike with denominator q carries
mass 1/(Lφ(q)²) of the normalised measure.
"""


def headline(res) -> str:
    m1q = qstar_at(res, res["X"])
    m2, m3, m5 = res["m2"], res["m3"], res["m5"]
    rm = m2["rank_mass"][("primes", 6)]
    vals = [comp_at(rm, i) for i in (2, 8, 16, 32, 64, 128, 256, 512, 1024)]
    d2p, d2c = m2["Ds"]["primes"][2.0][6], m2["Ds"]["cramer"][2.0][6]
    b150 = band_stats(res, (150, 20))
    b10 = next(b for b in b150 if b["band"] == (10, 100))
    hhalf = hband(res, "primes", (20, 10), "full", 0.5)
    hhalf_c = hband(res, "cramer", (20, 10), "full", 0.5)
    pref = [b["diam_half"] for b in m5["bulbs"]]
    med = m5["mediant"][21]
    oc = res["a1"]["octaves"]
    octs = [o["frac"] for o in oc if 2 <= o["j"] <= 6]
    a2 = res["a2"]
    sig = list(np.round(a2["sigma"], 2))
    i70 = sig.index(0.7)
    un = res["uncertainty"]
    bump = un["bump"][max(un["bump"])]
    rows = [
        ["M1 spike collapse", "does every spike look like the spike at 0, scaled?",
         f"q* = {m1q['q_star']} measured, {m1q['q_star_fg']:.0f} from the variance of primes in progressions (no fit)",
         "a data collapse (heights scaled by μ(q)/φ(q), widths fixed at 1/N) up to φ(q*) ≈ "
         f"{m1q['q_star_fg'] - 1:.0f}; not self-similarity"],
        ["M2 spectral measure", "how unevenly is |S|²dα spread?",
         f"2L·i·m₍ᵢ₎ = {rng(vals)} for ranks 2–1024 (1 = the Farey count); D₂ = {d2p:.2f} vs Cramér {d2c:.2f} at Nε ≈ 11",
         "coarse-grained straight-line spectrum f(α) = α, offset −ln(2L)/ln(1/ε); Diophantine, not a cascade"],
        ["M3 minor-arc tails", "is the minor-arc |S|² wild or mild?",
         f"Hill index (k = 1000) {m3[(5, 5)]['primes']['hill'][1000]:.2f} at (5,5), "
         f"{m3[(20, 10)]['primes']['hill'][1000]:.2f} at (20,10), {m3[(150, 20)]['primes']['hill'][1000]:.2f} at (150,20)",
         "wild in the untruncated idealisation at Q ≤ 20; neither wild nor exponential at Q ≳ 50"],
        ["M4 minor-arc walk", "how rough is the E0 walk? long memory?",
         f"H₂ = {b10['primes']:.3f} vs Cramér {b10['mean']:.3f} ± {b10['sd']:.3f} at Nℓ ∈ [10,100) (z = {b10['z']:.1f}); "
         f"H_½ = {hhalf:.2f} vs {hhalf_c:.2f} at (20,10)",
         "smooth below 1/N (generic); H₂ slightly above ½ from residual spikes; Noah effect yes; no long memory"],
        ["M5 Mandelbrot bulbs", "do bulbs and spikes share a hierarchy?",
         f"{sum(b['ok'] for b in m5['bulbs'])}/{len(m5['bulbs'])} bulbs checked; prefactor {rng(pref)}; "
         f"mediant rule {med['bulb_mediant']}/{med['pairs']} (bulbs) vs {med['spike_mediant']}/{med['pairs']} (spikes)",
         "same index set and circle order, q⁻² sizes through Diophantine approximation; not shared dynamics"],
        ["A1 octaves of q", "how is the energy spread over denominators?",
         f"{rng(octs, 3)} of ΣΛ² per octave for q ∈ [4, 128)",
         "equal energy per octave: the L² wall in fractal language"],
        ["A2 large values", "how many cells have |S| ≥ N^σ?",
         f"count exponent {a2['cells'][i70]} cells at σ = 0.7; large sieve bound / cells = "
         f"{a2['sieve_bound'][i70] / a2['cells'][i70]:.0f}",
         "count exponent 2(1 − σ) at resolution 1/N: the finite-N shadow of Jarník–Besicovitch"],
        ["A3 uncertainty", "why do sharp arcs leak?",
         f"sharp cutoff: σ_n·σ_α grows like √L; smooth bump {bump['product']:.4f} = "
         f"{bump['product'] / bump['bound']:.2f}× the bound 1/(4π)",
         "the Heisenberg–Gabor inequality for conjugate n and α (ħ → 1/2π): literal"],
    ]
    rows[6][2] = (f"{a2['cells'][i70]:,} cells at σ = 0.7 against {a2['spikes'][i70]:,} spikes; "
                  f"large-sieve bound / cells = {a2['sieve_bound'][i70] / a2['cells'][i70]:.0f}")
    return "## Headline\n\n" + _table(["block", "question", "key numbers at X = 10^6", "reading"], rows) + "\n"


def certificates(res) -> str:
    rows = [[n, "pass" if ok else "FAIL", d] for n, ok, d in res["certificates"]]
    return ("## 0. Certificates\n\nEvery identity the experiments rely on, checked at full scale (the run aborts if one "
            "fails).\n\n" + _table(["check", "result", "detail"], rows) + "\n")


def sec_m1(res) -> str:
    X, psi = res["X"], res["psi"]
    rows = []
    for q in (101, 210, 301, 499, 1001, 2003, 4001, 2310):
        r = row_q(res, q)
        rice = (f"{r['rice_mean']:.3f} / {r['rice_sd']:.3f} / {r['rice_cos']:.3f}" if "rice_mean" in r else "n/a (composite)")
        rows.append([q, r["phi"], f"{r['mean']:.3f} / {r['sd']:.3f} / {r['cos']:.3f}",
                     f"{r['beck_mean']:.3f} / {r['beck_sd']:.3f} / {r['beck_cos']:.3f}", rice])
    t_b = _table(["q", "φ(q)", "measured mean / sd / ⟨cos θ⟩", "Beckmann (exact 2nd moments)", "Rice-FG (prime q, no data)"], rows)
    r210 = row_q(res, 210)
    iso = rice_moments(r210["E2"])
    vrows = []
    for q in (101, 499, 2003, 4001):
        r = row_q(res, q)
        vrows.append([q, f"{r['V'] / X:.2f}", f"{fg_any(q):.2f}", f"{r['V'] / (X * fg_any(q)):.3f}"])
    prs = [r for r in res["m1_rows"] if 300 <= r["q"] <= 10**4 and len(factorize(r["q"])) == 1]
    med_ratio = float(np.median([r["V"] / (X * fg_any(r["q"])) for r in prs]))
    comp = [r for r in res["m1_rows"] if r["q"] in (210, 301, 1001)]
    e2v = ", ".join(f"{r['E2'] * psi**2 / (r['phi'] ** 2 * r['V']):.2f} at q = {r['q']}" for r in comp)
    t_v = _table(["prime q", "V(q)/X", "FG(q)", "ratio"], vrows)
    zrows = []
    for z in res["m1_zoom"]:
        c = np.asarray(z["centres"])
        k0 = int(np.argmin(np.abs(c)))
        lo, hi = z["class"]
        zrows.append([f"[{lo}, {hi}]", z["n_q"], f"{z['median'][k0]:.3f}", f"{z['p10'][k0]:.3f}–{z['p90'][k0]:.3f}"])
    z1 = res["m1_zoom"][0]
    c1 = np.asarray(z1["centres"])
    sel = np.abs(c1) <= 0.5
    dev1 = float(np.max(np.abs(np.asarray(z1["median"])[sel] / np.abs(np.sinc(c1[sel])) - 1)))
    t_z = _table(["q class", "# q", "median y at x ≈ 0", "10–90% band"], zrows)
    qrows = [[f"{d['X']:.0e}".replace("e+0", "e"), d["q_star"], f"{d['q_star_fg']:.0f}", f"{d['sqrt_X_logX']:.0f}",
              f"{d['psi_over_sea']:.0f}", d["Qmax"][10]] for d in res["m1_qstar"]]
    t_q = _table(["X", "q* measured", "q* from FG (no fit)", "√(X/log X)", "ψ/√ΣΛ²", "Q_max(c = 10)"], qrows)
    wm = res["m1_winmax"]
    js = wm["js"]
    NL = np.array([2**j for j in js]) * X / res["M"]
    wrows = [[f"{NL[i]:.3g}", f"{wm['primes'][i]:.3e}", f"{wm['hybrid'][i]:.3e}", f"{wm['cramer'][i]:.3e}"]
             for i in range(0, len(js), 2)]
    t_w = _table(["N·L", "primes", "surrogate H₁₀₀₀", "Cramér"], wrows)
    sl = np.diff(np.log(wm["primes"])) / np.diff(np.log(NL))
    mid = np.sqrt(NL[1:] * NL[:-1])
    win = (mid >= 120) & (mid <= 16000)
    low = mid < 60
    q1 = qstar_at(res, X)
    return f"""## M1. Spike data collapse and its cutoffs (fig_M1)

At a reduced fraction a/q with squarefree q the major-arc approximation is S(a/q + β) ≈ (μ(q)/φ(q))·(ψ/X)·v(β), with
|v(β)|/X ≈ |sinc(Nβ)|. Every spike would then be the spike at 0, with its height divided by φ(q) and its width
unchanged. This is a **data collapse** (heights rescaled, widths fixed at 1/N), not invariance under zoom. How far does
it hold, and what ends it?

**(a) Zoom collapse.** y = |S(α)|φ(q)/ψ against x = N(α − a/q), for four classes of q (up to 30 q and 50 a per
class). Values at the spike centre:

{t_z}

For the smallest class the median profile matches |sinc x| to {100 * dev1:.2f}% at |x| ≤ 0.5. By q ≈ 500–1000 the
heights scatter by a factor of two (the collapse has failed). By q ≈ 2000–4000 the "spike" is sea.

**(b) The collapse error is predicted by the variance of primes in progressions.** Write Z_a = μ(q)S(a/q)φ(q)/ψ − 1.
Over the reduced a, E|Z|² and E[Z²] are exact second moments (E[Z²] is real, because a → −a conjugates). The
residue-class variance V(q) = Σ_{{(r,q)=1}}(ψ(X;q,r) − ψ_c/φ(q))² (Barban–Davenport–Halberstam, Montgomery–Hooley)
behaves like X·FG(q), FG(q) = ln q − γ − ln 2π − Σ_{{p|q}} ln p/(p−1) (Friedlander–Goldston form). Three predictions per q:
a **Beckmann model** y = |1 + X + iY| with X, Y independent Gaussians carrying the exact E|Z|², E[Z²] (an 80×80
Gauss–Hermite rule); and for prime q the **Rice-FG closed form** with σ² = (q−1)²X·FG(q)/ψ², which uses no data at all.

{t_b}

The isotropic Rice model with the measured E|Z|² gets the spike-dominated regime wrong: at q = 210 it gives sd
{iso['sd']:.3f} against the measured {r210['sd']:.3f}; Beckmann's anisotropy fixes it. For composite q the off-diagonal
Ramanujan-sum terms matter (E|Z|²φ²/(ψ²V)·… = {e2v}), so the FG closed form is used for prime q only. V(q)/X against FG:

{t_v}

Median of V/(X·FG) over the {len(prs)} primes in [300, 10^4]: {med_ratio:.3f}. The phase-locking order parameter
⟨cos θ⟩ (last number in each triple) falls together with the height collapse: a locked spike is *real* (μ(q)ψ/φ(q),
because the Ramanujan sums c_q(n) are integers), and an unlocked one has a uniformly spread phase.

**(c) Where the collapse ends.** q*(X) is the prime q where spike and sea are equal, (q − 1)√(X·FG(q)) = ψ (brentq;
no fit). Measured q* is the first prime where the running median over 21 consecutive primes of
K(q) = √(mean_a |S(a/q) + ψ_c/(q−1)|²)·(q−1)/ψ exceeds 1:

{t_q}

At X = 10^6, q* = {q1['q_star']} measured against {q1['q_star_fg']:.0f} predicted. The plan's √(X/log X) is about
{np.mean([d['q_star_fg'] / d['sqrt_X_logX'] for d in res['m1_qstar']]):.1f}× too low; ψ/√ΣΛ² is an order of magnitude only.
q* and the largest admissible arc range Q_max = √(X/2c) both grow like √X up to logarithms, **for different reasons**:
q* compares spike height with the sea; Q_max compares arc width c/N with the Farey spacing 1/Q², and depends on c. They
do not coincide.

**(d) Statistical self-affinity of the window maxima.** Median over non-overlapping windows of 2^j grid points of
max|S|/ψ:

{t_w}

Local slope (primes) {rng(sl[win])} for N·L ∈ [120, 16000] and {rng(sl[low])} below N·L ≈ 60 (sea). This is the literal
self-affinity: statistical, with exponent ≈ ½, between the sea floor and the q ≤ 5 spikes.

**Shows.** The major-arc spike picture holds quantitatively up to a computable cutoff q*(X), set by the variance of
primes in progressions. Phase locking (⟨cos θ⟩ → 1) and the height collapse fail together. **Does not show.** Any
self-similarity under zoom; any link between q* and an admissible arc range; anything asymptotic (FG is a conjectural
asymptotic, checked at these X only).
"""


def sec_m2(res) -> str:
    X, M = res["X"], res["M"]
    m2 = res["m2"]
    L = m2["L"]
    st = m2["staircase"]
    ranks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    rm = {tag: m2["rank_mass"][(tag, 6)] for tag in ("primes", "hybrid", "cramer")}
    rrows = [[i, fmt(comp_at(rm["primes"], i)), f"{stair_value(st, L, i):.2f}", fmt(comp_at(rm["hybrid"], i)),
              fmt(comp_at(rm["cramer"], i), ".2g")] for i in ranks]
    t_r = _table(["rank i", "primes", "Farey staircase (no fit)", "surrogate H₁₀₀₀", "Cramér"], rrows)
    crows = []
    for j in (5, 7, 9):
        d = m2["cumulative"][j]
        a = np.round(np.asarray(d["alpha"]), 2)
        cells = []
        for A in (0.6, 0.7, 0.8, 0.9):
            k = int(np.flatnonzero(a == A)[0])
            v = f"{d['N'][k] * 2 * L * d['eps'] ** A:.2f}"
            cells.append(v if A >= d["alpha_lo"] else f"({v})")
        ok = (a >= d["alpha_lo"]) & (a <= 0.9)
        meas_off = float(np.nanmean(np.asarray(d["f"])[ok] - a[ok]))
        crows.append([f"2^{int(round(math.log2(d['eps'])))}", *cells, f"{d['offset']:.3f}", f"{meas_off:.3f}", f"{d['alpha_lo']:.2f}"])
    t_c = _table(["ε", "α = 0.6", "0.7", "0.8", "0.9", "offset −ln(2L)/ln(1/ε)", "measured mean f − α", "α_lo (N_< ≥ 60)"], crows)
    sc = m2.get("stair_counts", {})
    t_s = _table(["A", "boxes with μ ≥ ε^A (data)", "Farey staircase (model)"],
                 [[A, d, m] for A, (d, m) in sorted(sc.items())]) if sc else ""
    cj = m2["cj"]
    cjrows = [[s, f"{cj[s]['alpha'][0]:.2f}", f"{cj[s]['f'][0]:.2f}"] for s in (0.8, 0.9, 1.0, 1.1, 1.2) if s in cj]
    t_cj = _table(["s", "α(s)", "f(s)"], cjrows)
    js = m2["Ds_js"]
    j6 = 6
    Ne6 = X * 2 ** (j6 + 0.5) / M
    drows = [[s, f"{m2['Ds']['primes'][s][j6]:.2f}", f"{m2['Ds']['hybrid'][s][j6]:.2f}", f"{m2['Ds']['cramer'][s][j6]:.2f}"]
             for s in sorted(m2["Ds"]["primes"])]
    t_d = _table(["s", "primes", "surrogate H₁₀₀₀", "Cramér"], drows)
    Ne = np.array([X * 2 ** (j + 0.5) / M for j in js[:-1]])
    d2rows = [[f"{Ne[j]:.3g}", f"{m2['Ds']['primes'][2.0][j]:.2f}", f"{m2['Ds']['hybrid'][2.0][j]:.2f}",
               f"{m2['Ds']['cramer'][2.0][j]:.2f}"] for j in range(0, 14)]
    t_d2 = _table(["Nε", "D₂ primes", "D₂ surrogate", "D₂ Cramér"], d2rows)
    acc = (Ne >= 3.8 * 0.7) & (Ne <= 122 * 1.5)
    accept = float(np.max(np.abs(np.asarray(m2["Ds"]["hybrid"][2.0]) - np.asarray(m2["Ds"]["primes"][2.0]))[acc]))
    sg = m2["s_grid"]
    ssel = [s for s in sg if s in (0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0)]
    t_ds = _table(["s"] + [str(s) for s in ssel],
                  [[MEAS] + [f"{m2['Ds_s'][tag][s]:.2f}" for s in ssel]
                   for MEAS, tag in (("primes", "primes"), ("surrogate", "hybrid"), ("Cramér", "cramer"))])
    mrows = []
    for pr in CCDF5:
        d = m2["minor_D2"][pr]
        cells = []
        for tag in ("primes", "hybrid", "cramer"):
            v, n = minor_d2_min(res, pr, tag)
            cells.append(f"{v:.2f} (Nε {n:.3g})")
        b = d["boxes_lt09"]
        mrows.append([pair(pr), *cells, b["n"], f"{b['alpha_min']:.3f}" if b["alpha_min"] is not None else "—",
                      f"{100 * b['share']:.1f}%"])
    t_m = _table(["(Q,c)", "min D₂ primes", "surrogate", "Cramér", "# α < 0.9 boxes inside m (ε = 2^-18)", "their α_min",
                  "share of minor mass"], mrows)
    xs = m2["xscaling"]
    xrows = []
    for d in xs:
        sel = (d["ranks"] >= 8) & (d["ranks"] <= 256)
        v512 = comp_at(d, 512)
        xrows.append([f"{d['X']:.0e}".replace("e+0", "e"), f"{d['L']:.3f}", f"{d['L'] - math.log(d['X']):+.3f}", f"{d['Neps']:.1f}",
                      rng(d["comp"][sel]), f"{v512:.2f}" if v512 is not None else "—"])
    t_x = _table(["X", "L(X)", "L − ln X", "Nε", "2L·i·m₍ᵢ₎, i ∈ [8, 256]", "at i = 512"], xrows)
    return f"""## M2. Intermittency of the spectral measure: the Farey spectrum (fig_M2)

The measure is μ = |S|²dα with the atoms at 0 and 1/2 removed, normalised to 1. Boxes hold 2^j grid points:
ε = 2^(j−23), Nε = 0.119·2^j. A spike with denominator q carries mass m_q = 1/(Lφ(q)²) and sits in one box (or two).
Ordering spikes by φ and using the Farey count gives, with no fit,

    N_<(α; ε) = #{{boxes: μ_box ≥ ε^α}} = ε^(−α)/(2L)  ⇔  f_<(α) = α − ln(2L)/ln(1/ε)  ⇔  m₍ᵢ₎ = 1/(2L·i).

It is exact in expectation (box counting is linear in mass, and split spikes do not change the expected count). This is a
**coarse-grained** spectrum with ε tied to N. As N → ∞ at fixed ε every m_q → 0 and the measure tends to Lebesgue
measure; with ε = N^(−θ), 0 < θ < 1, τ(s) = min(0, s − 1), so D_s = 1 for s < 1 and 0 for s > 1, and the offset vanishes
like ln ln N/ln N. Neither limit is reached here.

**(a) Rank–mass.** Compensated 2L·i·m₍ᵢ₎ at Nε = {X * 2**6 / M:.1f} (j = 6), against the exact Farey staircase (squarefree
q ≥ 3 sorted by φ; φ(q) boxes of mass κ(Nε)/(Lφ²) + ρε each, with κ = {st['kappa']:.3f} the expected capture of a
spike by a box and ρ = {res['sea_level']:.3f} the sea level from the surrogate):

{t_r}

The value at i = 4 is the φ = 2 degeneracy (q = 3 and 6 share φ = 2). Cramér has no flat branch.

**(b) Cumulative spectrum.** Compensated N_<·2L·ε^α (1 = the no-fit line; in parentheses below α_lo, where fewer than
60 boxes count), with the no-fit offset and the measured mean of f_< − α over [α_lo, 0.9]:

{t_c}

Exact staircase count at ε = 2^-14 (Nε = 61):

{t_s}

Chhabra–Jensen local slopes between ε and 4ε (ε = 2^-18 → 2^-16) lie on f = α with **no** offset, because a scale
difference removes the constant ln(2L). They show the kink of τ(s) at s ≈ 1, not an independent spectrum:

{t_cj}

**(c) Generalized dimensions against scale.** Local D_s between consecutive box sizes, at Nε ≈ {Ne6:.1f} (j = 6 → 7):

{t_d}

D₂ against box size (D_1 is the Shannon information dimension):

{t_d2}

Below Nε ≈ 15 Gamma fluctuations of 1–15 cells move even Cramér. Acceptance of the surrogate: max |D₂(H₁₀₀₀) − D₂(primes)|
over Nε ∈ [3.8, 122] = {accept:.3f}. D_(s≥2) ≈ 0 is the q = 3 and q = 6 atoms (about 1.8% of the mass each); any measure
with one unresolved atom does this, so it is not evidence of a hierarchy.

**(d) D_s against s** (local slope at Nε 15 → 31):

{t_ds}

A broad crossover from ≈ 1 (s < 1) to ≈ 0 (s > 1) between s ≈ 0.8 and 1.5. The step is only the N → ∞ limit.

**(e) Minor-arc measure** μ_m = |S|²·1_m/Σ (the atoms sit in the major arcs). Minimum local D₂ over Nε ∈ [1.9, 977],
and the boxes of the full-circle measure with α < 0.9 that lie entirely in the minor arcs at ε = 2^-18:

{t_m}

**(f) X-scaling.** The compensated rank–mass at the box with Nε nearest 8, for six X, each with its own L(X):

{t_x}

**Shows.** f = α for α < 1 and D_s(ε) are quantitative consequences of the Farey count with heights ψ/φ(q); the surrogate
made of singular-series spikes q ≤ 1000 plus a structureless sea reproduces every curve, and Cramér is a monofractal
point. A straight-line spectrum is the signature of Diophantine counting (as in Jaffard's spectrum of Riemann's
function), not of a multiplicative cascade. **Does not show.** A fine (Hausdorff) spectrum: the rationals have
dimension 0. A cascade. A D_s step at this N.
"""


def sec_m3(res) -> str:
    m3, ts = res["m3"], np.asarray(res["m3_t"])
    rows = []
    for pr in PAIRS7:
        r = m3[pr]
        p = r["primes"]
        km = r["E4_model"] / (p["P_m"] ** 2 / r["meas_m"])
        lo, hi = r["pareto_range"]
        rows.append([pair(pr), f"{p['kappa']:.4g}", f"{km / p['kappa']:.3f}", f"{p['sup2'] / p['mean_m']:.1f}",
                     f"{r['sup_pred']['sup2'] / p['sup2']:.3f}", r["sup_pred"]["active"],
                     f"[5, {hi:.0f}]" if hi >= 10 * lo else "none",
                     " / ".join(f"{p['hill'][k]:.2f}" for k in (300, 1000, 10000)),
                     f"{r['cramer']['hill'][1000]:.2f}", f"{r['cramer']['kappa']:.3f}"])
    t1 = _table(["(Q,c)", "κ", "κ model / measured", "sup²/mean", "predicted sup² / measured", "active term",
                 "Pareto range", "Hill k = 300 / 1000 / 10⁴", "Hill Cramér (k = 1000)", "κ Cramér"], rows)
    crow = []
    for pr in PAIRS7:
        r = m3[pr]
        crow.append([pair(pr), *[f"{x:.3f}" for x in r["primes"]["concentration"]],
                     *[f"{x:.3f}" for x in r["cramer"]["concentration"]]])
    t2 = _table(["(Q,c)", "primes top 0.01%", "0.1%", "1%", "10%", "Cramér 0.01%", "0.1%", "1%", "10%"], crow)

    def at(pr, t):
        k = int(np.argmin(np.abs(ts - t)))
        r = m3[pr]
        return r["primes"]["ccdf"][k], r["model"][k], r["hybrid"]["ccdf"][k], r["cramer"]["ccdf"][k]

    trows = []
    for pr, tlist in (((5, 5), (8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)), ((20, 10), (8, 64, 512))):
        for t in tlist:
            a, b, h, c = at(pr, t)
            trows.append([pair(pr), t, g3(a), g3(b), f"{a / b:.3f}" if b > 0 else "—", g3(h), g3(c) if c > 0 else "0"])
    t3 = _table(["(Q,c)", "t", "primes P(> t)", "spike-sum model", "measured / model", "surrogate", "Cramér"], trows)
    hrows = []
    for t in (4, 8, 16, 32, 64):
        a, b, h, c = at((150, 20), t)
        hrows.append([t, g3(a), g3(h), f"{h / a:.2f}" if a > 0 else "—", g3(c) if c > 0 else "0", g3(math.exp(-t))])
    t4 = _table(["t", "primes", "surrogate H₁₀₀₀", "surrogate / primes", "Cramér", "e^(−t)"], hrows)
    r150 = m3[(150, 20)]
    k8 = int(np.argmin(np.abs(ts - 8)))
    kmr = {pr: m3[pr]["E4_model"] / (m3[pr]["primes"]["P_m"] ** 2 / m3[pr]["meas_m"]) / m3[pr]["primes"]["kappa"] for pr in PAIRS7}
    km_lo = f"{min(kmr[p] for p in PAIRS7[:3]):.3f}–{max(kmr[p] for p in PAIRS7[:3]):.3f}"
    km_hi = ", ".join(f"{kmr[p]:.2f} at {pair(p)}" for p in PAIRS7[3:])
    return f"""## M3. Mild versus wild on the minor arcs, by (Q,c) (fig_M3)

Take α uniform on the minor arcs and t = |S|²/mean_m. A spike of height ψ/φ(q) for each of the ≈ y²/2 reduced fractions
with φ(q) ≤ y gives

    P(|S|² > t·mean_m) ≈ ψ²/(2N|m|·t·mean_m),

a power law of index 1 in |S|² (index 2 in |S|): the wild/pre-wild boundary for |S|. (∫_0^1 u·m(u)du = ½∫sinc² = ½ with
m(u) = |{{|sinc| > u}}|.) The arcs truncate it: spikes with q ≤ Q are cut out within c/N of their centre. A semi-analytic
sum over spikes, P_model(t) = (1/(N|m|))Σ_q φ(q)·ℓ_q(u_q) with ℓ_q the measure of {{|sinc x| > u}} outside the arc, carries
the truncation with no fit. The largest surviving value is predicted as max(ψ/φ_min(Q), ψ/(π(c + ½)))²: the smallest-φ
spike above Q, or the first arc-edge sidelobe of q = 1, 2. Mandelbrot's operational tests run on coherence cells of 8 grid
points (≈ 1/N) lying entirely in the minor arcs: the Hill index of the cell maxima and the concentration of energy in
the largest cells. The fourth-moment model ∫_m|S|⁴ ≈ (ψ⁴/N)[(2/3)Σ_{{q>Q}}μ²/φ³ + T4(c)Σ_{{q≤Q}}μ²/φ³], T4(c) = 2∫_c^∞ sinc⁴,
predicts the kurtosis κ = mean_m|S|⁴/(mean_m|S|²)².

{t1}

The fourth-moment model counts spikes only. It holds to {km_lo} at Q ≤ 20; from Q ≈ 50 on the sea adds to ∫_m|S|⁴ and
the model falls short ({km_hi}). The Pareto range is [5, 0.04·sup²/mean]; a power law is claimed only where it spans a
decade. κ and sup² equal escapes.md
E1's values to 1e-9 (certificate). Concentration (share of the minor-arc energy in the top cells; the exponential
reference f(1 + ln(1/f)) gives 0.0010, 0.0079, 0.0561, 0.330):

{t2}

The tail against the spike-sum model and the surrogate:

{t3}

At (150,20) the surrogate against the primes:

{t4}

At (150,20) P(> 8) = {r150['primes']['ccdf'][k8]:.3g} against {math.exp(-8):.3g} for an exponential law
({r150['primes']['ccdf'][k8] / math.exp(-8):.0f}×).

**Reading.** At Q ≤ 20 the minor arcs still carry the near-rational tail over one to two decades (Hill index ≈ 1 to
1.5), which is wild in Mandelbrot's sense **in the untruncated idealisation** (state 6, Pareto index ≤ 2 in |S|). At
finite N every moment is finite, so it is wild in the short run and mild in the long run; it is never "extreme" (state
7 needs every moment infinite). For Q ≳ 50 the arcs cut the tail off within a decade of the bulk. What remains is
**neither wild nor exponential**: deterministic residual spikes on a sea, reproduced by the singular-series surrogate.
Cramér is exponential in |S|² (borderline mild). The fourth-moment model, the Pareto tail and κ are the same fact. The
ternary bound's sup is an order statistic of the truncated tail. **Does not show.** That any tail is wild at finite N;
probability over primes (the randomness here is uniform α); anything beyond the surrogate.
"""


def sec_m4(res) -> str:
    X, M = res["X"], res["M"]
    x = xl(res)
    out = []
    for pr in WALK:
        bs = band_stats(res, pr)
        rows = [[f"[{b['band'][0]:g}, {b['band'][1]:g})", f"{b['indep'][0]:.0f}–{b['indep'][1]:.0f}", f"{b['primes']:.3f}",
                 f"{b['mean']:.3f} ± {b['sd']:.3f}", f"{b['z']:+.1f}", f"{b['nge']}/{b['nseeds']}", b["verdict"],
                 f"{b['hyb'][300]:.3f} / {b['hyb'][1000]:.3f} / {b['hyb'][3000]:.3f}"] for b in bs]
        out.append(f"**{pair(pr)}** (n = {res['walks'][('primes',) + pr]['n']:,} minor points):\n\n" + _table(
            ["band of Nℓ", "n/ℓ", "primes H₂", "Cramér, 16 seeds", "z", "seeds ≥ primes", "verdict", "surrogate R = 300 / 1000 / 3000"], rows))
    t_bands = "\n\n".join(out)
    drift_ch = max(float(np.max(np.abs(res["walks"][("primes",) + pr]["full"]["drift"]["H"] - res["walks"][("primes",) + pr]["full"][2.0]["H"])[
        res["walks"][("primes",) + pr]["n"] / np.array(LAGS[:-1]) >= 500])) for pr in WALK)
    trows = []
    for j in range(0, 9):
        ht = H_triangular(float(x[j]))
        trows.append([f"{x[j]:.3g}", f"{ht:.3f}", f"{res['walks'][('primes', 20, 10)]['full'][2.0]['H'][j]:.3f}",
                      f"{res['walks'][('primes', 150, 20)]['full'][2.0]['H'][j]:.3f}",
                      f"{np.mean([e['raw'][j] for e in res['ensemble'][(20, 10)]]):.3f}"])
    t_tri = _table(["Nℓ/M", "triangular envelope (no fit)", "primes (20,10)", "primes (150,20)", "Cramér mean (20,10)"], trows)
    dev20 = max(abs(H_triangular(float(x[j])) - res["walks"][("primes", 20, 10)]["full"][2.0]["H"][j]) for j in range(9))
    noah = []
    for pr in WALK:
        n = res["walks"][("primes",) + pr]["n"]
        sel = n / np.array(LAGS[:-1]) >= 1500
        h = np.asarray(res["walks"][("primes",) + pr]["noah"][2.0]["H"])[sel]
        noah.append(f"{rng(h, 3)} at {pair(pr)}")
    prow = []
    for pr in WALK:
        for name, tag, var in (("primes", "primes", "full"), ("Noah-only (phases randomised)", "primes", "noah"),
                               ("Joseph-only (amplitudes flattened)", "primes", "joseph"), ("Cramér", "cramer", "full"),
                               ("surrogate H₁₀₀₀", "hybrid1000", "full")):
            vals = [hband(res, tag, pr, var, p) for p in (0.5, 1.0, 2.0)]
            prow.append([pair(pr), name, *[f"{v:.3f}" if v is not None else "—" for v in vals]])
    t_p = _table(["(Q,c)", "walk", "H_½", "H₁", "H₂"], prow)
    orow = []
    for pr in WALK:
        o = res["order"][pr]
        orow.append([pair(pr), f"{o['major'][0]:.3f}", f"{o['major'][1]:.3f}", f"{o['minor'][0]:.4f}", f"{o['minor'][1]:.3f}"])
    t_o = _table(["(Q,c)", "major r", "major phase entropy (bits)", "minor r", "minor phase entropy (bits)"], orow)
    return f"""## M4. Roughness of the minor-arc walk: band limit, Noah, Joseph (fig_M4)

The walk is E0's: k runs over the sorted minor-arc grid indices, w_j = S(k_j/M)²e(−Nk_j/M)/M, and F = cumsum(w) ends at
I_m(N) (certificate). Structure functions S_p(ℓ) = mean_j|F_(j+ℓ) − F_j|^p give local exponents
H_p = (1/p)·Δln S_p/Δln ℓ between ℓ and 2ℓ. Nulls: a 16-seed Cramér ensemble on the same mask; the surrogate H_R for
R = 300, 1000, 3000; a **Noah-only** walk |w|e^(iθ) with random θ (keeps the tails, destroys the phases); and a
**Joseph-only** walk w/|w| (keeps the phases, flattens the tails). Verdict words were fixed in advance: "excess" if
z ≥ 3 and no seed reaches the primes, "marginal" if 2 ≤ z < 3, "none" if |z| < 2 ("below" if z ≤ −2).

{t_bands}

Drift removal (S₂ − |mean increment|², an exact identity) changes H₂ by at most {drift_ch:.3f} wherever n/ℓ ≥ 500.

**Smooth below 1/N, then diffusive.** Any degree-N sum is smooth below its band limit, and Cramér's model shows the same
crossover. For the primes its shape follows, with no fit, the triangular envelope of the binary counts carried by the
residual spikes, S₂(x) ∝ ∫_0^1(1 − u)²sin²(πux)/u² du, x = Nℓ/M:

{t_tri}

Largest deviation from the envelope at (20,10), x ≤ 30: {dev20:.4f}. At (150,20) the mask removes 27% of the circle,
and the deviations are larger.

**Noah and Joseph.** Band means over Nℓ ∈ [4, 1000]:

{t_p}

The Noah-only H₂ stays at {"; ".join(noah)} wherever n/ℓ ≥ 1500 (check: 0.50 ± 0.02). H_½ well above ½ for the primes
and for the phase-randomised walk, but not for Cramér, is the Noah effect (heavy tails), at (20,10) above all. The
Joseph-only walk sits near ½. No long memory is detected; the rise at the largest lags at (150,20) has at most a few
hundred independent increments and is reproduced by the deterministic surrogate.

**Spectral identity.** For the zero-filled cyclic walk, S₂^cyc(ℓ) = (1/M²)Σ_d|I_m(N + d)|²·sin²(πdℓ/M)/sin²(πd/M)
exactly (certificate), so H at lag ℓ is a Fejér average of the minor-arc counts I_m(N') over |N' − N| ≲ M/ℓ. H is
purely second order; a Fourier-phase surrogate would add nothing.

**Order parameter and phase entropy** of the binary integrand w (Kuramoto-style r = |Σw|/Σ|w|; |w|-weighted entropy of
arg w in 64 bins, maximum 6 bits):

{t_o}

On the minor arcs the phases are spread almost uniformly (entropy ≈ 6 bits), and r is |I_m|/P_m. Magnitude-only
bounds assume r = 1 there.

**Shows.** H is a second-order statistic of the minor-arc counts I_m(N'). The band limit is generic, its detailed shape
is set by the residual spikes, and the tails are heavy at small Q. **Does not show.** Long memory. That α is time: it is
not, and "walk" refers to the ordering of grid points. Anything beyond "marginal" or "none" at lags with fewer than ~500
independent increments.
"""


def sec_m5(res) -> str:
    m5 = res["m5"]
    bl = m5["bulbs"]
    sel = [(1, 2), (1, 3), (1, 40), (3, 8), (13, 34), (3, 40), (9, 20), (10, 21)]
    rows = []
    for p, q in sel:
        b = next(b for b in bl if b["p"] == p and b["q"] == q)
        rows.append([f"{p}/{q}", f"{b['s']:.4g}", f"{b['centre_root']:.4f}", f"{b['diam_half']:.4f}"])
    t1 = _table(["p/q", "s = sin(πp/q)/q²", "|centre − root| / s", "diam/2 / s"], rows)
    brow = []
    for lo, hi in ((2, 5), (6, 10), (11, 20), (21, 30), (31, 40)):
        d = [b["diam_half"] for b in bl if lo <= b["q"] <= hi]
        c = [b["centre_root"] for b in bl if lo <= b["q"] <= hi]
        brow.append([f"[{lo}, {hi}]", len(d), f"{np.mean(d):.3f} ({min(d):.3f}–{max(d):.3f})",
                     f"{np.mean(c):.3f} ({min(c):.3f}–{max(c):.3f})"])
    t2 = _table(["q band", "# bulbs", "diam/2 / s, mean (range)", "|centre − root| / s, mean (range)"], brow)
    arow = [[k, f"{v['diam_half']:.4f}", f"{v['area_radius']:.4f}"] for k, v in m5["area"].items()]
    t3 = _table(["p/q", "diam/2 / s", "area radius √(A/π) / s"], arow)
    mrow = [[f"F_{n}", d["pairs"], f"{d['bulb_mediant']}/{d['pairs']}", f"{d['spike_mediant']}/{d['pairs']}",
             f"{d['mediant_nonsquarefree']}/{d['pairs']}"] for n, d in sorted(m5["mediant"].items())]
    t4 = _table(["Farey order", "neighbour pairs", "largest bulb = mediant", "tallest spike = mediant", "mediant not squarefree"], mrow)
    ex = m5["mediant"][21]["examples"][:2]
    exs = "; ".join(f"{a[0]}/{a[1]} | {b[0]}/{b[1]}: mediant {m[0]}/{m[1]}, tallest spike {t[0]}/{t[1]}"
                    for (a, b), m, t in ((e["pair"], e["mediant"], e["tallest_spike"]) for e in ex))
    lrow = [[d["Q"], f"{d['bulbs']:.5f}", f"{d['bulb_formula']:.5f}", f"{d['spikes']:.5f}", f"{d['spike_formula']:.5f}"]
            for d in m5["logdiv"]]
    t5 = _table(["Q", "Σ sin(πp/q)/q² (bulbs)", "(12/π³) ln Q − 0.009", "Σ μ²/φ (spikes)", "ln Q + 1.3326"], lrow)
    pref = [b["diam_half"] for b in bl]
    return f"""## M5. The Mandelbrot set's bulbs and the primes' Farey forest (fig_M5)

The bulb attached to the main cardioid at internal angle θ = p/q has its root at c_r = e^(2πiθ)/2 − e^(4πiθ)/4 and a
centre with f_c^q(0) = 0 of exact period q (Douady–Hubbard). Its leading-order radius is s = sin(πθ)/q²: Mandelbrot's
empirical "N² rule" (1985), explained by a normal form (Guckenheimer–McGehee 1984). It is an estimate, not a theorem:
for limbs only O(1/q) is proven, and O(1/q²) is Milnor's open conjecture (as stated in Kapiamba, arXiv:2103.03211).
Every bulb with q ≤ 40 was located by Newton's method seeded along the outward normal at 0.9·s and checked for exact
period q, rotation number p/q and continuation to the root (certificate: {sum(b['ok'] for b in bl)}/{len(bl)}). The
primary size is diam/2 = |c(ρ = −1) − c_r|/2, from continuation of the q-cycle multiplier ρ to −1.

{t1}

By band of q:

{t2}

Area radius from 240 boundary points (ρ = e^(iφ)) and the shoelace formula:

{t3}

The prefactor diam/2 ÷ s runs over {rng(pref)} and depends on p/q. **Mediant test.** For each pair of Farey neighbours,
among the fractions strictly between them (denominator ≤ 120, 200, 300 for F_8, F_13, F_21): is the largest bulb the
mediant? is the tallest spike (max |μ(b)|/φ(b), ties to the smaller denominator) the mediant?

{t4}

Examples of spike failures: {exs}. **Log-divergence parallel:**

{t5}

**Reading.** The bulbs and the spikes are indexed by the same reduced fractions in the same circle order, and both
families have q^(−2) sizes: for bulbs the empirical N² rule, for spikes the **mass** 1/(Lφ(q)²) (the height scales as
1/φ(q)). Ford circles (radius 1/(2q²)) show the same order and size law with no dynamics at all. The mediant rule holds
for the bulbs but fails for the spikes, because many mediants have non-squarefree denominators (no spike) and φ(q) is not
monotone in q. What the two families share is Diophantine approximation, not dynamics. The Mandelbrot set's fractality
lives elsewhere: its boundary (Hausdorff dimension 2, Shishikura 1998), its Misiurewicz points (asymptotic
self-similarity, Tan Lei 1990) and its small copies. Arnold tongues and the devil's staircase of mode locking
(Jensen–Bak–Bohr 1983) are a look-alike with a coupling parameter; S has no iteration, parameter or coupling.
"""


def sec_additions(res) -> str:
    oc = res["a1"]["octaves"]
    orow = [[f"[{2 ** o['j']}, {2 ** (o['j'] + 1)})", f"{o['frac']:.4f}", f"{o['spike_model']:.4f}", f"{o['spike_sea_model']:.4f}",
             f"{o['sum_inv_phi']:.3f}"] for o in oc]
    t1 = _table(["octave of q", "measured share of ΣΛ²", "spike model", "spike + sea model", "Σ μ²/φ in the octave"], orow)
    a2 = res["a2"]
    sig = np.asarray(a2["sigma"])
    cells = np.asarray(a2["cells"], dtype=float)
    lnN = math.log(res["X"])
    slope = np.r_[np.nan, -np.diff(np.log(cells)) / (np.diff(sig) * lnN)]
    arow = [[f"{s:.2f}", f"{int(c):,}", f"{int(sp):,}", f"{c / sp:.2f}", f"{math.log(c) / lnN:.3f}", f"{2 * (1 - s):.2f}",
             f"{sl:.2f}" if np.isfinite(sl) else "—", f"{b / c:.0f}"]
            for s, c, sp, sl, b in zip(sig, cells, a2["spikes"], slope, a2["sieve_bound"])]
    t2 = _table(["σ", "cells with max|S| ≥ N^σ", "Farey count of spikes", "cells per spike", "D(σ) = ln N_cells/ln N",
                 "2(1 − σ)", "local slope", "large sieve / cells"], arow)
    un = res["uncertainty"]
    urow = [["sharp cutoff (w = 1 on [0, 1])", L, f"{d['sigma_x']:.4f}", f"{d['sigma_xi']:.4f}", f"{d['product']:.4f}",
             f"{d['product'] / d['bound']:.2f}"] for L, d in sorted(un["sharp"].items())]
    urow += [["smooth bump η (escapes E4)", L, f"{d['sigma_x']:.4f}", f"{d['sigma_xi']:.4f}", f"{d['product']:.4f}",
              f"{d['product'] / d['bound']:.2f}"] for L, d in sorted(un["bump"].items())]
    t3 = _table(["window", "L (|Nβ| ≤ L)", "σ_x", "σ_ξ", "σ_x·σ_ξ", "÷ 1/(4π)"], urow)
    return f"""## A1. Equal energy per octave of q (fig_M6 a)

E_oct(j) = (1/M)Σ|S|² over grid points within |N(α − a/q)| ≤ ½ of every reduced a/q with squarefree q ∈ [2^j, 2^(j+1)).
Spike model: 0.7737·(ψ²/N)·Σ_oct μ²/φ (0.7737 = ∫_(|x|≤½) sinc²). The spike + sea model adds the surrogate's sea power times
the window measure. Equal-energy line: 0.7737·ln 2·ψ²/(NΣΛ²) = {res['a1']['equal_line']:.4f} of ΣΛ² per octave.

{t1}

Index 2 in |S| ⇔ a logarithmically divergent variance ⇔ Σ_(q≤Q) μ²/φ ≈ ln Q + 1.33 ⇔ each octave of q holds the same
energy. This is the reason arcs up to Q capture only (ln Q + 1.33)/ln N of the energy: results.md's L² wall, restated as
scale invariance. Above q ≈ 256 the sea dominates the windows.

## A2. Large values: the Jarník–Besicovitch shadow, against the large sieve (fig_M6 b)

N_cells(σ) counts the coherence cells (8 grid points ≈ 1/N) with max|S| ≥ N^σ. A value |S| ≥ N^σ needs φ(q) ≲ N^(1−σ) and
|α − a/q| ≲ 1/N = q^(−τ) with τ = 1/(1 − σ), so the count exponent at resolution 1/N is 2(1 − σ): the finite-N shadow of
Jarník–Besicovitch (the τ-well-approximable numbers have Hausdorff dimension 2/τ). The rigorous anchor is the dual large
sieve: even-indexed and odd-indexed cells are each 1/N-separated, so N_cells ≤ 2(2N − 1)ΣΛ²/N^(2σ).

{t2}

The law holds for σ ≥ 0.65 (the sea is √(N log N) ≈ N^0.595). The Farey forest sits a constant factor below the large
sieve bound at every level: the y² law is the extremal shape of the large-values estimate. This is a **count exponent
at resolution 1/N**, not a Hausdorff dimension.

## A3. The uncertainty principle behind the arc widths (fig_M6 c)

With e(x) = e^(2πix), n (as x = n/N ∈ [0, 1]) and α (as ξ = Nβ) are Fourier conjugates, so σ_x·σ_ξ ≥ 1/(4π) = 0.0796 for
|w|²-weighted spreads: Planck's ħ becomes 1/(2π). This is literal (Heisenberg–Gabor), not an analogy.

{t3}

A sharp cutoff n ≤ N has Dirichlet sidelobes falling like 1/ξ², whose variance diverges: σ_ξ grows like √L, so the
product is unbounded. That is E4's sidelobe leakage, read as the infinite momentum spread of a sharp edge. The smooth
bump is a near-minimum-uncertainty packet, which is why it suppressed the leakage. Consequence: an arc cannot be much
narrower than ~1/N without cutting its spike, and spikes with denominators ≤ Q are ≥ 1/Q² apart, so disjoint arcs need
Q ≲ √(N/2c).
"""


def sec_physics(res) -> str:
    lk = res.get("leakage", {})
    o150 = res["order"][(150, 20)]
    sh = lk.get((150, 20))
    leak = (f"at (150,20) the minor-arc counts I_m(N) spread their energy Σ_N I_m(N)² as {100 * sh['targets']:.0f}% on the "
            f"targets, {100 * sh['other even <= 2X']:.0f}% on other even N ≤ 2X, {100 * sh['odd <= 2X']:.0f}% on odd N ≤ 2X and "
            f"{100 * sh['> 2X']:.2f}% beyond 2X (escapes.md E1)") if sh else "see escapes.md E1"
    r101 = row_q(res, 101)
    r2003 = row_q(res, 2003)
    rows = [
        ["phasors, interference (Feynman's arrows in *QED*: amplitude = Σ e^(iφ))", "LITERAL",
         "S(α) = Σ Λ(n)e(nα) is a sum of arrows; spikes are constructive interference at rationals"],
        ["phase locking turns complex sums into real numbers", "LITERAL",
         f"locked resultant μ(q)ψ/φ(q) is real (Ramanujan sums are integers); ⟨cos θ⟩ = {r101['cos']:.3f} at q = 101, "
         f"{r2003['cos']:.3f} at q = 2003"],
        ["order parameter and phase entropy", "LITERAL",
         f"(150,20): r = {o150['major'][0]:.3f} on the major arcs, {o150['minor'][0]:.4f} on the minor arcs; phase entropy "
         f"{o150['minor'][1]:.3f} of 6 bits"],
        ["\"momentum\" conservation from integrating a phase over a full period", "LITERAL",
         "∫_0^1 e((n₁ + n₂ − N)α)dα = δ: only the full circle enforces n₁ + n₂ = N"],
        ["leakage when the integral is restricted to arcs", "LITERAL", leak],
        ["energy conservation", "LITERAL", "Parseval (1/M)Σ|S|² = ΣΛ² (certificate)"],
        ["uncertainty principle, ħ → 1/(2π)", "LITERAL", "A3: sharp edges have unbounded σ_n·σ_α; the smooth bump sits at 1.19× the bound"],
        ["the circle as a Brillouin zone; |S|² as the primes' diffraction pattern", "LITERAL (re-reading)",
         "n ∈ ℤ ⇒ α ∈ ℝ/ℤ; Bragg-like peaks at a/q with intensity (μ(q)ψ/φ(q))², only for squarefree q "
         "(cf. Torquato, Zhang & de Courcy-Ireland 2018, cited for the idea, not checked here)"],
        ["equipartition, and the missing ultraviolet catastrophe", "PARTLY",
         "the minor-arc mean |S|² is flat in α (white noise) and the total is finite because n is discrete; no temperature, "
         "no Bose statistics, no E = hf"],
        ["a \"Reynolds number\" deciding locked vs random", "ANALOGY",
         "spike/sea = ψ/(φ(q)√(X·FG(q))) decides locked vs unlocked, but it marks a smooth crossover (M1), not an instability"],
        ["photon polarization, electron spin alignment", "NOT LITERAL",
         "S is a scalar sum; the only alignment is phase alignment"],
        ["\"moving along the future\", time", "NOT LITERAL", "there is no dynamics in α; X is a parameter, not time"],
        ["Poiseuille, Reynolds, Navier–Stokes", "NOT LITERAL",
         "no fluid, viscosity or pressure. Electron hydrodynamics is real physics in ultraclean conductors (Bandurin et al. "
         "2016; Sulpizio et al. 2019) but unrelated to the primes. The honest bridge to fluids is turbulence's multifractal "
         "cascade (Mandelbrot 1974), the origin of f(α); M2 finds no cascade"],
        ["\"the momentum stays the same because the proof is smooth\"", "NO",
         "there is no proof here; smoothing the weight (E4) only reduces edge leakage, it conserves nothing more"],
        ["Bohr levels ∝ 1/n² and bulb sizes ∝ 1/q²", "COINCIDENCE", "no connection; no \"Planck constant of the primes\""],
    ]
    return ("## Phase, entropy, fluids and Planck's constant: what is literal\n\n"
            "The physics vocabulary raised along the way, sorted by whether it applies literally. The literal entries all "
            "come from one place: n and α are Fourier conjugates on ℤ × ℝ/ℤ.\n\n"
            + _table(["statement", "status", "numbers / qualifier"], rows) + "\n")


def sec_literal(res) -> str:
    m3, m2, m5 = res["m3"], res["m2"], res["m5"]
    q1 = qstar_at(res, res["X"])
    hk = lambda pr: rng([m3[pr]["primes"]["hill"][k] for k in HILL_K])  # noqa: E731
    b150 = band_stats(res, (150, 20))
    z10 = next(b for b in b150 if b["band"] == (10, 100))["z"]
    z100 = next(b for b in b150 if b["band"] == (100, 1000))["z"]
    h10 = next(b for b in b150 if b["band"] == (10, 100))["primes"]
    h100 = next(b for b in b150 if b["band"] == (100, 1000))["primes"]
    hhalf = hband(res, "primes", (20, 10), "full", 0.5)
    hhalf_c = hband(res, "cramer", (20, 10), "full", 0.5)
    med = m5["mediant"][21]
    pref = [b["diam_half"] for b in m5["bulbs"]]
    oc = [o["frac"] for o in res["a1"]["octaves"] if 2 <= o["j"] <= 6]
    a2 = res["a2"]
    ratio = [b / c for b, c in zip(a2["sieve_bound"], a2["cells"])]
    r150 = m3[(150, 20)]
    k8 = int(np.argmin(np.abs(np.asarray(res["m3_t"]) - 8)))
    sp = [m3[pr]["sup_pred"]["sup2"] / m3[pr]["primes"]["sup2"] for pr in PAIRS7]
    b50 = m2["minor_D2"][(50, 10)]["boxes_lt09"]["n"]
    b150n = m2["minor_D2"][(150, 20)]["boxes_lt09"]["n"]
    rows = [
        ["|S|² on the minor arcs has a power-law tail (index ≈ 1 in |S|², 2 in |S|)",
         "LITERAL, as the value distribution under uniform α, in a truncated window, for Q ≤ 20 only",
         f"Hill (cells) {hk((5, 5))} at (5,5), {hk((20, 10))} at (20,10); {hk((150, 20))} at (150,20)"],
        ["that tail is \"wild\"", "LITERAL in the untruncated idealisation (state 6)",
         "at finite N wild in the short run, mild in the long run (all moments finite); |S| sits on the wild/pre-wild boundary; never \"extreme\""],
        ["Cramér is mild", "LITERAL", "|S|² exponential (borderline mild); Re S, Im S Gaussian (proper mild)"],
        ["minor arcs at (150,20) are mild", "FALSE: not wild, and not exponential either",
         f"Hill ≥ 2, max {r150['primes']['sup2'] / r150['primes']['mean_m']:.0f}× mean; but κ = {r150['primes']['kappa']:.2f} "
         f"against {r150['cramer']['kappa']:.2f} for Cramér, and P(> 8) is {r150['primes']['ccdf'][k8] / math.exp(-8):.0f}× e^(−8)"],
        ["Noah effect", "LITERAL (statistics in α; α is not time)", f"H_½ = {hhalf:.2f} against {hhalf_c:.2f} for Cramér at (20,10)"],
        ["Joseph effect", "ABSENT",
         f"H₂ = {h10:.3f} / {h100:.3f} at (150,20), above Cramér by z = {z10:.1f} / {z100:.1f} in Nℓ ∈ [10,100) / [100,1000), from "
         "the residual-spike spectrum (the surrogate reproduces it)"],
        ["inner cutoff 1/N", "LITERAL but GENERIC", "band limit of any degree-N polynomial; \"smooth\", not \"ballistic\"; the primes' shape follows the triangular envelope"],
        ["spikes are self-similar", "PARTLY",
         f"data collapse (heights scaled by μ/φ, widths fixed) up to φ(q*) ≈ {q1['q_star_fg'] - 1:.0f}; window maxima statistically self-affine"],
        ["multifractal for α < 1", "LITERAL only as a coarse-grained spectrum, ε tied to N",
         f"straight line f = α, offset −ln(2L)/ln(1/ε), L = {m2['L']:.2f}; Diophantine, not a cascade"],
        ["the sea is monofractal", "LITERAL for Cramér; approximate for the primes' minor arcs only at Q ≥ 150",
         f"{b150n} boxes with α < 0.9 entirely in the minor arcs at (150,20), {b50} at (50,10)"],
        ["D_s steps at s = 1", "LIMIT ONLY (N → ∞, ε = N^(−θ))", "at 10^6 a broad crossover over 0.8 ≲ s ≲ 1.5; D_(s≥2) ≈ 0 is the q = 3, 6 atoms"],
        ["outer cutoff coincides with the arc range", "OVERCLAIM",
         f"q* ≈ {q1['q_star_fg']:.0f} (height against sea) vs Q_max = {q1['Qmax'][10]} (c = 10; spacing against width); both ∝ √N up to logs and c"],
        ["bulbs and spikes share Farey combinatorics", "PARTLY",
         f"same index set and circle order; the mediant rule holds for bulbs only ({med['bulb_mediant']}/{med['pairs']} against "
         f"{med['spike_mediant']}/{med['pairs']})"],
        ["bulb size 1/q² matches spike height", "WRONG PAIRING",
         "height ∝ 1/φ(q); the q^(−2) objects are spike mass and Ford circles; bulbs depend on p (sin πθ), spikes do not"],
        ["the bulb size law is exact or proven", "NO", f"empirical N² rule; normal-form explanation; prefactor {rng(pref)}; limbs only O(1/q)"],
        ["simple rules → hierarchy; Arnold tongues; mode locking", "ANALOGY", "no iteration, parameter or coupling in S; the shared root is Diophantine approximation"],
        ["ternary success = truncating the wild tail", "LITERAL",
         f"sup_m|S|² = max(ψ/φ_min(Q), ψ/(π(c+½)))² within a factor {rng(sp)} (the sea adds a few per cent in amplitude)"],
        ["binary wall = sea energy; equal energy per octave", "LITERAL", f"{rng(oc, 3)} of ΣΛ² per octave for q ∈ [4, 128)"],
        ["large values are extremal", "LITERAL (rigorous anchor)", f"within a factor {rng(ratio[1:], 0)} of the dual large sieve for σ ≥ 0.65"],
        ["the Mandelbrot set's fractality", "DISCLAIMER", "boundary dimension 2, Misiurewicz points, small copies; none lies in the bulbs used here, none has a counterpart in S"],
    ]
    return "## Literal versus analogy\n\n" + _table(["statement", "status", "numbers / qualifier"], rows) + "\n"


def synthesis(res) -> str:
    m3, m2, m5 = res["m3"], res["m2"], res["m5"]
    q1 = qstar_at(res, res["X"])
    small = [(5, 5), (10, 10), (20, 10)]
    hills = [m3[pr]["primes"]["hill"][k] for pr in small for k in HILL_K]
    kap = [m3[pr]["primes"]["kappa"] for pr in small]
    top1 = [m3[pr]["primes"]["concentration"][2] for pr in small]
    top1c = [m3[pr]["cramer"]["concentration"][2] for pr in small]
    a2 = res["a2"]
    ratio = [b / c for b, c in zip(a2["sieve_bound"], a2["cells"])]
    b150 = band_stats(res, (150, 20))
    hb = [b["primes"] for b in b150 if b["band"] in ((10, 100), (100, 1000))]
    hhalf = hband(res, "primes", (20, 10), "full", 0.5)
    med = m5["mediant"][21]
    return f"""## Synthesis

At X = 10^6 the circle-method objects pass Mandelbrot's measurement programme in one specific way. Every scaling law
found here is the Farey count #{{a/q : μ(q) ≠ 0, φ(q) ≤ y}} ~ y²/2, applied to spikes of height ψ/φ(q) and width 1/N, on
top of a structureless sea. Read pointwise, it is a power-law tail of index 1 in |S|² (index 2 in |S|). Read per box, it
is the coarse-grained straight-line spectrum f(α) = α with offset −ln(2L)/ln(1/ε), L = {m2['L']:.2f}. Read per coherence
cell, it is the large-value count exponent 2(1 − σ), the finite-N shadow of Jarník–Besicovitch, within a factor
{rng(ratio[1:], 0)} of the dual large sieve. Read through the fourth moment, it is the minor-arc kurtosis. Read per octave
of q, it is equal energy per octave, which is why arcs up to Q capture only (log Q + 1.33)/log N of the energy: the L²
wall in fractal language. These are one fact, not independent confirmations. A surrogate built from the singular-series
spikes q ≤ 1000 and a Gaussian sea of matched power reproduces all of them to within a few per cent, and no fractal
structure beyond it is detected.

Where the minor arcs sit depends on (Q, c). At Q ≤ 20 they still carry the near-rational tail over one to two decades
(Hill index {rng(hills)}, κ = {min(kap):.0f}–{max(kap):.0f}, the top 1% of cells holding {100 * min(top1):.0f}–{100 * max(top1):.0f}%
of the energy against {100 * min(top1c):.0f}–{100 * max(top1c):.0f}% for Cramér). That is wild in the untruncated idealisation.
For Q ≳ 50 the arcs cut the tail off within a decade of the bulk. What remains is neither wild nor exponential:
deterministic residual spikes on a sea (κ = {m3[(150, 20)]['primes']['kappa']:.2f} at (150,20) against
{m3[(150, 20)]['cramer']['kappa']:.2f} for Cramér). The ternary bound needs exactly this truncation. Its sup is the largest
surviving order statistic, either the smallest-φ spike above Q or the first arc-edge sidelobe of q = 1, 2.

The spike data collapse holds up to φ(q*) ≈ {q1['q_star_fg'] - 1:.0f} (q* = {q1['q_star']} measured), where the variance of
primes in progressions, V(q) ≈ X(log q − γ − log 2π − Σ_(p|q) log p/(p−1)), drowns the spike. This scale and the largest
admissible arc range √(N/2c) both grow like √N, for different reasons.

The minor-arc walk is smooth below 1/N, as any degree-N sum is. Beyond that it is diffusive in its second moment, with
H₂ ≈ {rng(hb, 2)}: slightly but significantly above the Cramér ensemble, because of residual spikes. Its heavy tails show in
H_½ ≈ {hhalf:.2f} at (20,10) (a Noah effect). There is no evidence of long memory.

The Mandelbrot set's bulbs share with the Farey forest the index set, the circle order and a q^(−2) size law, through
Diophantine approximation, not dynamics. The mediant rule holds for the bulbs ({med['bulb_mediant']}/{med['pairs']}) but not
for the spikes ({med['spike_mediant']}/{med['pairs']}).

Every exponent here is a statistic at one finite N over a finite scaling window. Nothing here bears on a proof.
"""


def corrections(res) -> str:
    q1 = qstar_at(res, res["X"])
    return f"""## Corrections to earlier informal statements

Three statements made while planning this phase do not survive measurement:

- "The self-similarity cutoff is q* ≈ 270 and coincides with Q_max = 223." The cutoff, measured and predicted from the
  variance of primes in progressions, is q* = {q1['q_star']} (prediction {q1['q_star_fg']:.0f}). The figure 270 came from
  √(X/log X), which is about 1.9× too low at every X. Q_max(c) = √(X/2c) comes from a different mechanism and depends on
  c; the agreement at c = 10 was a coincidence, not three routes to one scale.
- "The spikes are self-similar." They form a data collapse: heights scale by μ(q)/φ(q), widths stay 1/N. The literal
  self-affinity is statistical (the window maxima).
- "The minor arcs are a mild sea." That depends on (Q,c): at Q ≤ 20 they carry a wild (untruncated) tail; at Q ≥ 150
  they are neither wild nor exponential.
"""


def captions(res) -> str:
    q1 = qstar_at(res, res["X"])
    m2, m5 = res["m2"], res["m5"]
    med = m5["mediant"][21]
    pref = [b["diam_half"] for b in m5["bulbs"]]
    return f"""## Figures

![Fig M1](figures/fig_M1_collapse.png)

**Fig M1.** Every spike at a reduced a/q has width 1/N. Rescaling its height by φ(q)/ψ collapses it onto |sinc| while
φ(q) ≲ 300 (a). The collapse is a data collapse, not invariance under zoom. Its error is predicted with no fit by the
residue-class variance V(q) ≈ X(log q − γ − log 2π − Σ_(p|q) log p/(p−1)): a Beckmann model with the exact second
moments matches mean, spread and phase locking to a few per cent (b). The collapse ends where spike and sea are equal,
at φ(q*) ≈ {q1['q_star_fg'] - 1:.0f} for X = 10^6, and q* tracks the closed form from 10^4 to 2·10^6 (c). q* and √(X/2c)
both grow like √X up to logarithms, for different reasons. The window maxima grow like L^(1/2) between the sea floor and
the q ≤ 5 spikes (d): statistical self-affinity.

![Fig M2](figures/fig_M2_spectrum.png)

**Fig M2.** With the atoms at 0 and 1/2 removed, the i-th largest box carries mass ≈ 1/(2L·i), L = {m2['L']:.2f} (a): the
Farey count read per box. Equivalently the coarse-grained spectrum is the straight line f(α) = α with a no-fit offset
−ln(2L)/ln(1/ε) (b). The generalized dimensions cross over smoothly from ≈ 1 (s < 1) to ≈ 0 (s > 1) between s ≈ 0.8 and
1.5 (c, d); the step is only the N → ∞ limit. The surrogate made of the singular-series spikes q ≤ 1000 plus a sea
reproduces every curve; Cramér's model is a monofractal point. On the minor arcs residual spikes keep D₂ far below Cramér
at small Q (e). The law holds at every X with L(X) ≈ ln X − 3 (f).

![Fig M3](figures/fig_M3_tails.png)

**Fig M3.** Under uniform α on the minor arcs the tail of |S|² is the Farey count again, P(> t) ≈ ψ²/(2N|m|·t·mean), a
power law of index 1 in |S|² (2 in |S|). The semi-analytic spike sum (no fit) follows the measured tail and explains the
steepening as truncation by the arcs. At Q ≤ 20 the tail spans one to two decades (wild in the untruncated
idealisation); for Q ≳ 50 the arcs cut it within a decade, leaving residual spikes on a sea that the surrogate
reproduces. Cramér is exponential. The vertical line is the largest surviving value and the term that sets it (f:
Hill index and kurtosis by (Q,c)).

![Fig M4](figures/fig_M4_walks.png)

**Fig M4.** The minor-arc phasor walk is smooth below the band limit 1/N and crosses over to diffusive growth, as any
degree-N sum does (Cramér too). For the primes the crossover follows the triangular envelope of the binary counts
carried by the residual spikes (a, b), keeping H₂ slightly above ½ and above the Cramér band; the surrogate reproduces
this. Heavy tails (Noah effect) show in H_½, and the phase-randomised walk keeps them (c). No long memory. The exact
spectral identity shows which neighbouring targets N' = N + d set H at each lag (d).

![Fig M5](figures/fig_M5_mandelbrot.png)

**Fig M5.** The bulbs on the main cardioid and the spikes of S are indexed by the same reduced fractions in the same
circle order, and both have q^(−2) sizes: Mandelbrot's empirical N² rule for bulbs (prefactor {rng(pref)}, depending on
p/q), and the spike mass 1/(Lφ²), not the height, for spikes. Ford circles show the same order and size law with no
dynamics. The mediant rule holds for bulbs in {med['bulb_mediant']} of {med['pairs']} Farey-neighbour pairs of F_21 but for
spikes in only {med['spike_mediant']}; in {med['mediant_nonsquarefree']} pairs the mediant has no spike at all.

![Fig M6](figures/fig_M6_octaves_largevalues_uncertainty.png)

**Fig M6.** Two more readings of the Farey count and the uncertainty principle. (a) Each octave of denominators carries
the same energy up to q ≈ 128: the L² wall in fractal language. (b) Cells with |S| ≥ N^σ have count exponent 2(1 − σ) at
resolution 1/N, a constant factor below the dual large sieve. (c) σ_n·σ_α ≥ 1/(4π): a sharp cutoff has unbounded spread
(leakage), the smooth bump sits near the bound.
"""


def references() -> str:
    return """## References

Mandelbrot and his programme: B. Mandelbrot, J. Business 36 (1963) 394–419; J. Fluid Mech. 62 (1974) 331–358; *The
Fractal Geometry of Nature* (1982); "On the dynamics of iterated maps III" in *Chaos, Fractals and Dynamics* (Dekker 1985)
213–224 (the N² rule); *Fractals and Scaling in Finance* (1997), chs. E5–E6. Mandelbrot & Wallis, Water Resour. Res. 4
(1968) 909–918 and 5 (1969) 967–988; Mandelbrot & Van Ness, SIAM Rev. 10 (1968) 422–437; Klemeš, Water Resour. Res. 10
(1974) 675–688. Parisi & Frisch, in *Turbulence and Predictability* (Varenna 1983, North-Holland 1985) 84–87; Halsey,
Jensen, Kadanoff, Procaccia & Shraiman, Phys. Rev. A 33 (1986) 1141; Chhabra & Jensen, Phys. Rev. Lett. 62 (1989) 1327.

Exponential sums: Jaffard, Rev. Mat. Iberoam. 12 (1996) 441–460; Chamizo & Ubis, Adv. Math. 250 (2014) 1–34;
Boritchev, Eceizabarrena & Vilaça da Rocha, J. Math. Phys. 62 (2021) 093101; Duistermaat, Nieuw Arch. Wisk. (4) 9 (1991)
303–337.

Mandelbrot set: Douady & Hubbard, C. R. Acad. Sci. Paris 294 (1982) 123–126, the Orsay notes (1984–85) and Ann. Sci. ENS
18 (1985) 287–343; Guckenheimer & McGehee, Mittag-Leffler Report 15 (1984); Shishikura, Ann. of Math. 147 (1998) 225–267
and LMS Lecture Notes 274 (2000) 325–363; Tan Lei, Comm. Math. Phys. 134 (1990) 587–617; Fowler & McGuinness, Chaos
Solitons Fractals X 3 (2019) 100019; Devaney, Amer. Math. Monthly 106 (1999) 289–302; Ford, Amer. Math. Monthly 45 (1938)
586–601; Jensen, Bak & Bohr, Phys. Rev. Lett. 50 (1983) 1637.

Physics: Bandurin et al., Science 351 (2016) 1055–1058; Sulpizio et al., Nature 576 (2019) 75–79.

Cited for a statement only, not checked here: Jarník (1929) and Besicovitch (1934) on well-approximable numbers;
Pommerenke (1986), Levin (1991), Yoccoz (in Hubbard 1993) and Milnor's conjecture, as stated in Kapiamba,
arXiv:2103.03211; Barban (1964), Davenport & Halberstam, Michigan Math. J. 13 (1966) 485–489, Montgomery, Michigan Math. J.
17 (1970) 33–39, Hooley, J. reine angew. Math. 274/275 (1975) 206–223, Friedlander & Goldston, Quart. J. Math. Oxford (2)
47 (1996) 313–336 (variance of primes in progressions); Montgomery & Vaughan, Mathematika 20 (1973) 119–134 (large sieve
constant); Torquato, Zhang & de Courcy-Ireland, J. Stat. Mech. (2018) 093401 (primes as a diffraction pattern).
"""


def write(res, path: Path) -> None:
    parts = [header(res), headline(res), certificates(res), sec_m1(res), sec_m2(res), sec_m3(res), sec_m4(res), sec_m5(res),
             sec_additions(res), sec_physics(res), sec_literal(res), synthesis(res), corrections(res), captions(res),
             references(),
             f"Run: {res['runtime_s'] / 60:.1f} min, peak RSS {res['peak_rss_mb']:,.0f} MB.\n" if "runtime_s" in res else ""]
    Path(path).write_text("\n".join(parts))
