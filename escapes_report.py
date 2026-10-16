"""Write escapes.md from the dict produced by run_escapes.compute() (all numbers are read from it)."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from run_all import _e, _f

NUMERIC = set("0123456789.,%+−-eE ()/–·<>≈")


def _table(head: list[str], rows: list[list[str]]) -> str:
    """Markdown table; | inside cells is escaped. Columns whose cells are all number-like are right-aligned."""
    esc = lambda cells: [str(c).replace("|", "\\|") for c in cells]  # noqa: E731
    num = [all(set(str(r[j])) <= NUMERIC for r in rows) for j in range(len(head))]
    out = ["| " + " | ".join(esc(head)) + " |", "|" + "|".join("---:" if nm else ":---" for nm in num) + "|"]
    out += ["| " + " | ".join(esc(r)) + " |" for r in rows]
    return "\n".join(out)
from sieve import factorize
from singular_series import SS2_product

SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
SHORT5 = [(150, 20), (100, 20), (50, 10), (20, 10), (5, 5)]


def lab(x: float) -> str:
    """Decision margin for a ratio computed from Riemann sums on the 2^23 grid."""
    if x is None or not np.isfinite(x):
        return "n/a"
    return "closes" if x < 0.99 else ("undecided" if x <= 1.01 else "fails")


def r3(x: float) -> str:
    """Ratio to 3 decimals, flagged when it falls inside the undecided band."""
    if x is None or not np.isfinite(x):
        return "—"
    s = f"{x:.3f}" if x >= 0.0095 else f"{x:.2e}"
    return s + " (undecided)" if lab(x) == "undecided" else s


def first_close(vals, start: int = 0):
    for K in range(start, len(vals)):
        if vals[K] is not None and np.isfinite(vals[K]) and vals[K] < 0.99:
            return K
    return None


def kstr(K, vals=None, kmax: int = 14) -> str:
    if K is None:
        return f"> {kmax}"
    return f"{K} ({vals[K]:.3f})" if vals is not None else str(K)


def fact(n: int) -> str:
    return "·".join(f"{p}" + (str(e).translate(SUP) if e > 1 else "") for p, e in sorted(factorize(n).items()))


def pc(x: float, n: int) -> str:
    return f"{x:,.0f} ({100 * x / n:.2g}%)" if x < 10 * n else f"{x:.3g} ({100 * x / n:,.0f}%)"


def pair(p) -> str:
    return f"({p[0]},{p[1]})"


# ====================================================================== sections

def certificates(res) -> str:
    R2, X = res["R2"], res["X"]
    e1, e2, e3, e4 = res["e1"], res["e2"], res["e3"], res["e4"]
    if "sum_r2_sq" in res:
        sr2 = res["sum_r2_sq"]
    else:  # older cache: recompute (2 s)
        from sieve import autoconvolve, von_mangoldt
        sr2 = math.fsum(autoconvolve(von_mangoldt(X)) ** 2)
    full4 = e1[(150, 20)]["full4"]
    bes = max(o["bessel_relerr"] for o in e1.values())
    rk = np.abs(e2["RK_grid"].real - e2["RK_direct"]) / e2["RK_direct"]
    rki = np.abs(e2["RK_grid"].imag) / e2["RK_direct"]
    e3err = max(abs(s["grid"] - s["direct"]) / s["direct"] for s in e3)
    e3im = max(abs(s["grid_imag"]) / s["direct"] for s in e3)
    eb = abs(e4["R_eta"] - e4["R_eta_grid"]) / e4["R_eta"]
    et = max(abs(t["R3"] - t["R3_grid"]) / t["R3"] for t in e4["ternary"].values())
    e0 = res["e0"]
    walk = max(abs(o["I_m"].real - (R2 - o["I_M"])) / R2 for o in e0.values())
    walki = max(abs(o["I_m"].imag) / R2 for o in e0.values())
    pm = max(abs(o["P_m"] - e1[k]["P_m"]) / o["P_m"] for k, o in e0.items())
    cert = all(o["certificates_ok"] for o in e1.values()) and all(s["certificates_ok"] for s in res["e5_scaling"])
    rows = [
        ["Bessel identity Σ_N I_m(N)² = ∫_m|S|⁴, all 20 pairs", f"max rel. error {_e(bes, 1)}"],
        ["full-circle 4th moment ∫|S|⁴ = Σ_N r2(N)² (additive energy)", f"{_e(full4, 6)} vs {_e(sr2, 6)}, rel. {_e(abs(full4 - sr2) / sr2, 1)}"],
        ["E0: walk endpoint = R2 − I_M; path length = E1's P_m", f"{_e(walk, 1)}·R2 (imag ≤ {_e(walki, 1)}·R2); {_e(pm, 1)}"],
        ["E1/E5: every transmitted N has r2(N) ≥ I_M(N) − P_m > 0 (20 pairs + 12 scaling runs)", "all hold" if cert else "**FAILED**"],
        ["E2: R_K grid = direct, K ≤ 14 (grid exact)", f"max rel. {_e(rk[:15].max(), 1)}, imag ≤ {_e(rki[:15].max(), 1)}"],
        ["E2: K = 15, 16 (past the exactness limit; the aliased tuples are negligible)", f"rel. {_e(rk[15], 1)}, {_e(rk[16], 1)}"],
        [f"E3: R^(A)(X) grid = direct, all {len(e3)} sets", f"max rel. {_e(e3err, 1)}, imag ≤ {_e(e3im, 1)}"],
        ["E4: R_η(X) and R3_η(N') grid = direct (4 targets)", f"{_e(eb, 1)}; max {_e(et, 1)}"],
        ["E4: R_η / main, R3_η / main", f"{e4['R_eta'] / e4['main_bin']:.5f}; "
         + ", ".join(f"{t['R3'] / t['main_eta']:.5f}" for t in e4["ternary"].values())],
    ]
    tr = [[f"{t['N']:,}", fact(t["N"]), f"2^{int(math.log2(t['M']))}", pair((t["Q"], t["c"])), f"{t['P_m'] / t['I_M']:.3f}",
           lab(t["ratio"])] for t in res["e5_trunc"]]
    return ("## 0. Certificates\n\nEvery identity the experiments rely on, checked at full scale (M = 2^23 unless stated). "
            "The FFT sign convention is the one documented and checked in results.md.\n\n"
            + _table(["check", "result"], rows)
            + "\n\nE5 single-N certificates recomputed with S truncated at N itself (P_m(S_N) < P_m(S_X), so the X = 10^6 "
              "arrays are conservative for N < X), arcs (Q_max(N, 10), 10):\n\n"
            + _table(["N", "factorisation", "M", "(Q,c)", "P_m / I_M(N)", "verdict"], tr)
            + f"\n\nRun: {res['runtime_s'] / 60:.1f} min, peak RSS {res['peak_rss_mb']:,.0f} MB.\n")


def sec_e0(res) -> tuple[str, dict]:
    R2 = res["R2"]
    out = {}
    rows = []
    for k, o in res["e0"].items():
        nd = [d["N_Delta"] for d in o["dial"]]
        W = [d["W"] for d in o["dial"]]
        dstar = None
        for j in range(1, len(W)):
            if W[j - 1] >= 1 > W[j]:
                t = math.log(W[j - 1]) / (math.log(W[j - 1]) - math.log(W[j]))
                dstar = nd[j - 1] * (nd[j] / nd[j - 1]) ** t
        best = min(o["model"], key=lambda m: m["P_E"])
        bnd = [m["bound"] / o["I_M"] for m in o["model"]]
        out[k] = dict(dstar=dstar, best=best, bnd=(min(bnd), max(bnd)))
        rows.append([pair(k), f"{o['P_m'] / R2:.3f}", f"{o['max_abs_F'] / R2:.4f}", f"{o['I_m'].real / R2:+.4f}",
                     f"{o['I_M'] / R2:.3f}", f"{o['random_phase_rms'] / R2:.4f}", f"{o['rp_scale'] / R2:.4f}",
                     f"{o['holder'] / R2:.2f}", f"{o['sup_L1'] / R2:.2f}",
                     f"{best['P_E'] / R2:.2f} (Q' = {best['Qp']})", f"{min(bnd):.2f}–{max(bnd):.2f}",
                     f"{o['dial_floor']:.3f} ({o['n_components']:,})", f"≈ {dstar:.0f}/N" if dstar else "—"])
    a, b = res["e0"][(150, 20)], res["e0"][(20, 10)]
    oa, ob = out[(150, 20)], out[(20, 10)]
    dial_rows = []
    for da, db in zip(a["dial"], b["dial"]):
        dial_rows.append([f"{da['B']:,}", f"{da['N_Delta']:.2f}", f"{da['W']:.3f}", f"{db['W']:.3f}"])
    dial_rows.append(["whole components", "∞", f"{a['dial_floor']:.3f}", f"{b['dial_floor']:.3f}"])
    lad = [[f"{m['Qp']}", f"{m['P_E'] / R2:.3f}", f"{m['P_star'] / R2:.3f}", f"{m['J_star'].real / R2:+.3f}",
            f"{m['bound'] / a['I_M']:.2f}"] for m in a["model"]]
    indep = (a["P_m"] / a["I_M"]) ** 2
    txt = f"""## 1. E0 — the obstacle, measured

For one target N = X, the minor arcs add M·|m| phasors w_k = S(k/M)² e(−Nk/M)/M to the count. Laid head to tail they
form a walk F_j = Σ_(i≤j) w_(k_i): its **path length** is P_m = ∫_m|S|², its **endpoint** is I_m(N) = R2(N) − I_M(N).
A bound that knows only magnitudes must allow every step to point the same way, so it can never beat P_m.
Every column below is in units of R2(X) = {_f(R2, 1)} except the last four.

{_table(["(Q,c)", "path P_m", "max |F|", "endpoint I_m", "I_M", "rms, 1/N cells", "σ_rp = P_m/√(N|m|)", "Hölder √(|m|∫_m|S|⁴)", "sup_m|S|·∫_m|S|", "sea after model, best", "model bound / I_M", "W(∞) (components)", "Δ*"], rows)}

**Bound ladder at (150,20).** Every magnitude-only refinement we tried is at least as long as the path: Hölder is
{a['holder'] / a['P_m']:.2f}·P_m, sup × L¹ is {a['sup_L1'] / a['P_m']:.1f}·P_m. Subtracting a major-arc model S\\* (spikes
(μ(q)/φ(q))·v(α − a/q) on disjoint windows of width 1/(qQ′), squarefree q ≤ Q′) and bounding the direction that
could spoil the count, −I_m ≤ P_E + 2√(P\\*·P_E) − Re J\\*, gives {oa['bnd'][0]:.2f}–{oa['bnd'][1]:.2f}·I_M. The decoherent remainder
P_E = ∫_m|S − S\\*|² bottoms out at {oa['best']['P_E'] / R2:.2f}·R2 near Q′ = {oa['best']['Qp']} and rises again beyond,
where the spikes are no longer resolved by windows of width 1/(qQ′). That remainder is not itself a bound; it is still
about twice the signal.

{_table(["Q′", "P_E / R2", "P\\* / R2", "Re J\\* / R2", "bound / I_M"], lad)}

**Partial-phase dial.** Cut each connected component of the minor arcs into windows of B grid points (width Δ = B/M,
N·Δ = B·X/M) and add |window integral| over all windows: W(B) = Σ|∫_window S²e(−Nα)| / I_M. B = 1 knows only magnitudes
(W = P_m/I_M); whole components know all the phase inside each gap between arcs.

{_table(["B", "N·Δ", "W at (150,20)", "W at (20,10)"], dial_rows)}

W drops below 1 at N·Δ ≈ {oa['dstar']:.0f} at (150,20) and ≈ {ob['dstar']:.0f} at (20,10) (log–log interpolation between
the tabulated B). Independent cells would need only (P_m/I_M)² = {indep:.1f}. The rms of the cell sums at width 1/N is
{a['random_phase_rms'] / a['rp_scale']:.1f}× the random-phase yardstick σ_rp at (150,20) and
{b['random_phase_rms'] / b['rp_scale']:.1f}× at (20,10): cells of width 1/N are not independent, because the minor arcs
still contain unresolved spikes at a/q with q > Q, each coherent over about 1/N. At (20,10) those spikes (q = 21…) carry
coherent local integrals of size about N/φ(q)² that cancel only between different a/q with the same q (a Ramanujan sum,
which is non-local), so the windows must be much wider before W < 1.

**Reading.** At one N, with numbers: the binary wall is a statement about phases. The walk's path is
{a['P_m'] / a['I_M']:.1f}× the radius it would have to leave to spoil the count, yet it never strays more than
{100 * a['max_abs_F'] / a['I_M']:.1f}% of the way out (at (150,20); {b['P_m'] / b['I_M']:.1f}× and
{100 * b['max_abs_F'] / b['I_M']:.1f}% at (20,10)). Knowing the phase-summed integral over windows about
{oa['dstar']:.0f}/N wide would be enough at (150,20). That is a descriptive statistic of this S, not a quantity any
theorem controls: proving it for every large N would need a minor-arc estimate uniform in N, which is the open problem.
"""
    return txt, out


def sec_e1(res) -> tuple[str, dict]:
    e1 = res["e1"]
    n = e1[(150, 20)]["n_targets"]
    rows = []
    for p in res["pairs"]:
        o = e1[p]
        w = o["worst"][0]
        rows.append([pair(p), _e(o["E4m"], 3), f"{o['kurtosis']:.3g}", pc(o["bound_4th"], n), pc(o["bound_natural_sup"], n),
                     pc(o["bound_textbook"], n), pc(o["bound_oracle"], n), f"{o['F']}",
                     f"{w[1]:.3f} ({w[0]:,}; {'harmless' if w[2] > 0 else 'harmful'} sign)",
                     f"{o['pm_ratio_min']:.3f}–{o['pm_ratio_median']:.2f}–{o['pm_ratio_max']:.2f}",
                     f"{len(o['transmitted'])}"])
    a = e1[(150, 20)]
    q = a["quantiles"]
    sh = a["shares"]
    lo = min(o["shares"]["targets"] for o in e1.values())
    hi = max(o["shares"]["targets"] for o in e1.values())
    grid_pairs = [p for p in res["pairs"] if p not in ((158, 20), (223, 10))]
    worst_grid = sorted({e1[p]["worst"][0][0] for p in grid_pairs})
    worst_qmax = [(p, e1[p]["worst"][0]) for p in ((158, 20), (223, 10))]
    s = {"lo": lo, "hi": hi}
    txt = f"""## 2. E1 — average over N (the failure count)

The set of targets is 𝒩 = even N in [X/2, X], |𝒩| = {n:,}. For every N at once, I_M(N) = ∫_M S²e(−Nα) is one
length-M forward FFT of S²·1_M (numpy's sign is exactly the demodulation by e(−Nα)), and I_m(N) = R2(N) − I_M(N) with R2
from the exact autoconvolution (M > 2X, so R2(N) = ∫S²e(−Nα) exactly for N ≤ X). Discrete Parseval in N gives the
**Bessel identity** Σ_(N mod M) I_m(N)² = ∫_m|S|⁴ (checked in §0). With the **failure count**
F = #{{N ∈ 𝒩 : |I_m(N)| ≥ I_M(N)}} (failures of the method's inequality, not numbers that fail to be sums of two primes),
Chebyshev gives F ≤ U / min_𝒩 I_M(N)², with four choices of U:

* **textbook**: U = sup_m|S|²·ΣΛ² (the form the 1937–38 proofs use, with Vinogradov's sup bound);
* **natural sup**: U = sup_m|S|²·P_m, valid because ∫_m|S|⁴ ≤ sup_m|S|²·∫_m|S|² (the ternary device);
* **true 4th moment**: U = ∫_m|S|⁴, measured — no theorem gives it;
* **oracle**: U = Σ_(N∈𝒩) I_m(N)², Chebyshev applied to the actual values — no theorem gives it either.

sup is the grid maximum of |S|² on the minor arcs. κ = ∫_m|S|⁴/(P_m²/|m|) is the minor-arc kurtosis (2 for a Rayleigh
sea). "harmful sign" means I_m(N) < 0.

{_table(["(Q,c)", "∫_m|S|⁴", "κ", "F ≤ (4th moment)", "F ≤ (sup²·P_m)", "F ≤ (sup²·ΣΛ²)", "F ≤ (oracle)", "F", "max |I_m|/I_M (N)", "P_m/I_M(N) min–median–max", "# P_m < I_M(N)"], rows)}

**Reading.** F = 0 at every pair: no target in 𝒩 has |I_m(N)| ≥ I_M(N). At (150,20) the median |I_m|/I_M is
{q[0.5]:.4f}, the 99th percentile {q[0.99]:.4f}, the 99.9th {q[0.999]:.4f}, the worst {q[1.0]:.4f}. Yet for a single
target the magnitude bound P_m/I_M(N) is {a['pm_ratio_min']:.2f}–{a['pm_ratio_max']:.2f} there: each receiver alone sees
the wall. Taken together, the receivers share one energy budget, so few can fail: the true 4th moment allows
{a['bound_4th']:,.0f} failures ({100 * a['bound_4th'] / n:.2f}%), the sup form {a['bound_natural_sup']:,.0f}
({100 * a['bound_natural_sup'] / n:.1f}%), the textbook form {a['bound_textbook']:,.0f}
({100 * a['bound_textbook'] / n:.1f}%). At (20,10) the sup forms allow more failures than there are targets while the
4th moment still allows {100 * e1[(20, 10)]['bound_4th'] / n:.1f}%; at Q = 5 even the 4th moment allows
{100 * e1[(5, 5)]['bound_4th'] / n:.0f}%, because the tall spikes at q = 6, 7, 10, … left on the minor arcs make the sea very
non-Gaussian (κ ≈ {e1[(5, 5)]['kurtosis']:.0f}).

Where the budget goes: only {100 * lo:.0f}–{100 * hi:.0f}% of ∫_m|S|⁴ lands on 𝒩 ({100 * sh['targets']:.1f}% at (150,20));
other even N ≤ 2X take {100 * sh['other even <= 2X']:.1f}%, odd N ≤ 2X {100 * sh['odd <= 2X']:.1f}%, N > 2X
{100 * sh['> 2X']:.2f}%. At the 18 grid pairs the worst target is always highly composite
({", ".join(f"{w:,} = {fact(w)}" for w in worst_grid)}), where the singular-series terms with q > Q that the arcs leave
out are largest. At the two Q_max pairs that deterministic leakage is small and the worst targets are not special:
{"; ".join(f"{w[0]:,} = {fact(w[0])} at {pair(p)} ({w[1]:.3f}, I_m {'< 0' if w[2] < 0 else '> 0'})" for p, w in worst_qmax)}.

**What this is.** It reproduces, at one X, the mean-square (Bessel + Chebyshev) step of the exceptional-set theorems of
van der Corput, Chudakov and Estermann (1937–38), who proved E(X) ≪ X(log X)^(−A) using Vinogradov's minor-arc bound.
Those proofs use the sup form, which at (150,20) gives {100 * a['bound_natural_sup'] / n:.1f}%, not
{100 * a['bound_4th'] / n:.2f}%, and they need a provable major-arc model; here I_M is computed from the true S.
Montgomery–Vaughan (1975) proved E(X) ≪ X^(1−δ); their new input is on the major arcs (moduli up to X^δ, zero-density
estimates, a possible exceptional zero), which this experiment does not test. One X cannot tell X^(1−δ) from
X(log X)^(−A): Q = 150 is both X^0.36 and (log X)^1.9. And since every even N ≤ 4·10^18 is known to be a sum of two
primes, F counts failures of the method, not of Goldbach.
"""
    return txt, s


def sec_e5(res) -> str:
    e1 = res["e1"]
    o = e1[(223, 10)]
    tr = sorted(o["transmitted"], key=lambda t: t[1])
    rows = [[f"{N:,}", fact(N), f"{ss:.3f}", f"{r:.3f}"] for N, r, ss in tr]
    o158 = e1[(158, 20)]
    pr = o["pm_ratio_targets"].astype(float)
    nontr = float(pr[pr >= 1.0].min())
    sc = res["e5_scaling"]
    srows = []
    for Xs in sorted({r["X"] for r in sc}):
        cells = [f"{Xs:,}"]
        for c in (10, 20):
            r = next(r for r in sc if r["X"] == Xs and r["c"] == c)
            cells += [f"{r['Q']}", f"{r['n_transmitted']} / {r['n_targets']:,} ({100 * r['n_transmitted'] / r['n_targets']:.2g}%)",
                      f"{r['need']:.2f}"]
        cells.append(f"{next(r['max_SS2'] for r in sc if r['X'] == Xs):.2f}")
        srows.append(cells)
    a = e1[(150, 20)]
    nX = SS2_product(res["X"])
    p1, p2 = 9_699_690, 223_092_870
    first = next(r for r in sc if r["c"] == 10 and r["X"] == min(x["X"] for x in sc))
    last = next(r for r in sc if r["c"] == 10 and r["X"] == max(x["X"] for x in sc))
    minpm_grid = min(e1[p]["pm_ratio_min"] for p in res["pairs"] if p not in ((158, 20), (223, 10)))
    return f"""## 3. E5 — single-N transmission with no relaxation

P_m does not depend on N; the signal I_M(N) ≈ SS2_Q(N)·N does. Whenever I_M(N) > P_m, the plain magnitude inequality
(the triangle inequality |I_m(N)| ≤ P_m) certifies R2(N) ≥ I_M(N) − P_m > 0 for that one N (in floating point; checked
against the exact r2 in §0).
need(X) = P_m·SS2(X)/I_M(X) is the singular series a target near X would need to get through.

At X = 10^6 no grid pair transmits anything (closest: P_m/I_M(N) = {minpm_grid:.3f}). At the two Q_max pairs the arcs
cover more of the circle and the wall drops below the signal of the most composite targets:
**{len(o158['transmitted'])}** targets at (158,20) ({", ".join(f"{t[0]:,}" for t in sorted(o158['transmitted']))}) and
**{len(tr)}** at (223,10). The closest non-transmitted target at (223,10) has P_m/I_M(N) = {nontr:.3f}. For N = 10^6
itself (SS2 = {nX:.2f}) the ratio is {float(a['pm_ratio_targets'][-1]):.2f} at (150,20) and {float(pr[-1]):.2f} at
(223,10).

Transmitted targets at (223,10), best margin first:

{_table(["N", "factorisation", "SS2(N)", "P_m / I_M(N)"], rows)}

**Scaling.** The same experiment for each cached S_X, arcs (Q_max(X, c), c), targets even N in [X/2, X]:

{_table(["X", "Q (c=10)", "transmitted (c=10)", "need (c=10)", "Q (c=20)", "transmitted (c=20)", "need (c=20)", "max SS2 on 𝒩_X"], srows)}

**Reading.** The transmitted fraction falls from {100 * first['n_transmitted'] / first['n_targets']:.1f}% at
X = {first['X']:,} to {100 * last['n_transmitted'] / last['n_targets']:.4f}% at X = {last['X']:,} (c = 10). need(X)
rises at every step of the table, as the no-fit heuristic of results.md §5 suggests (P_m/R2 ≈ (log N − log Q − 2.33)/SS2
with Q = Q_max ≈ √(N/2c), so the wall grows like ½·log N). The largest singular series available below X
grows only like log log X, in steps at primorial-type N (SS2 = {SS2_product(p1):.2f} at {p1:,} and
{SS2_product(p2):.2f} at {p2:,}). We do not extrapolate where the two curves part for good. This corrects the premise that
nothing gets through at X = 10^6, and shows why it is not a route to Goldbach: it works only for a vanishing set of very
smooth N, with arcs covering about 30% of the circle and I_M computed from the true S rather than a provable model.
"""


def sec_e2(res) -> tuple[str, dict]:
    e2 = res["e2"]
    L = e2["L"]
    R2 = res["R2"]
    rk = e2["RK_direct"]
    rel = np.abs(e2["RK_grid"].real - rk) / rk
    rho = rk / (float(L) ** np.arange(len(rk)) * R2)
    krows = [[f"{K}", _e(rk[K], 4), _e(rel[K], 1), f"{rho[K]:.3f}"] for K in range(0, 15)]
    jrows = []
    for p in [(150, 20), (100, 20), (50, 10), (20, 10), (10, 5), (5, 5), (223, 10)]:
        r = e2["rows"][p]
        jrows.append([pair(p)] + [r3(v) for v in r["joint_ratio"][:7]])
    trows = []
    out = {}
    for p in res["pairs"]:
        r = e2["rows"][p]
        th = r["theorem"]
        ks = {C: first_close([v[0] for v in th[C]["split"]], 2) for C in (1, 2, 4)}
        kl = {C: first_close([v[0] for v in th[C]["literature"]], 2) for C in (1, 2, 4)}
        kn = {C: first_close(th[C]["no_split"], 2) for C in (1, 2, 4)}
        ke = first_close([v[0] for v in th[1]["split_exact"]], 2)
        kj, ksup, kh = first_close(r["joint_ratio"]), first_close(r["sup_ratio"]), first_close(r["holder_ratio"])
        out[p] = dict(joint=kj, sup=ksup, holder=kh, split=ks, lit=kl, nosplit=kn, exact=ke)
        trows.append([pair(p), kstr(kj, r["joint_ratio"]), kstr(ksup, r["sup_ratio"]), kstr(kh, r["holder_ratio"]),
                      " / ".join(kstr(ks[C]) for C in (1, 2, 4)), " / ".join(kstr(kl[C]) for C in (1, 2, 4)),
                      " / ".join(kstr(kn[C]) for C in (1, 2, 4)), kstr(ke)])
    a = e2["rows"][(150, 20)]
    b = e2["rows"][(20, 10)]
    Ks = [1, 2, 3, 4, 8, 14]
    sp = e2["sparse_control"]
    th = a["theorem"]
    lam = {C: th[C]["split"][out[(150, 20)]["split"][C]][1] for C in (1, 2, 4) if out[(150, 20)]["split"][C] is not None}
    pp = [(a["joint_ratio"][K + 1] / a["joint_ratio"][K]) for K in range(0, 6)]
    I2 = a["I_MK"][2]
    sup_s = b["sup_ratio"][14]
    txt = f"""## 4. E2 — add K powers of 2 (Linnik's problem)

G(α) = Σ_(k=1..{L}) e(2^k α) (powers 2…2^{L}; k ≥ 1 keeps the parity), and
R_K(N) = Σ_(k_1..k_K) R2(N − 2^(k_1) − … − 2^(k_K)) = ∫S²G^K e(−Nα) at N = 10^6. The grid is exact for K ≤ 14 (the
frequencies reach 2X + K·2^{L} and alias onto N + M only from K = 15). ρ_K = R_K/(L^K·R2).

{_table(["K", "R_K(10^6) (direct)", "grid vs direct", "ρ_K"], krows)}

Minor-arc magnitude bounds, divided by the major part I_M,K = ∫_M S²G^K e(−Nα):

* **joint** B_K = ∫_m|S|²|G|^K — the sharpest magnitude-only bound;
* **sup**: (sup_m|G|)^K·P_m, with grid sup_m|G| = {a['sup_G']:.3f} = {a['sup_G'] / L:.3f}·L for Q ≥ 10;
* **Hölder**: (∫_m|S|⁴)^(1/2)·(∫_m|G|^(2K))^(1/2);
* **theorem-shaped** (the substitutes the proofs actually use):
  B^thm_K(λ, C) = (λL)^(K−2)·[L·ΣΛ² + C·OD] + L^K·sup_m|S|²·meas(E_λ ∩ m), minimised over λ ∈ [0.30, 0.84], where
  E_λ = {{|G| ≥ λL}}, L·ΣΛ² is the diagonal of ∫|S|²|G|² (known asymptotically) and OD = ∫|S|²|G|² − L·ΣΛ² is the
  off-diagonal (correlations ΣΛ(n)Λ(n+h), h = 2^k − 2^l), for which a proof pays a sieve constant C. Variants: the
  literature form puts C on the whole ∫|S|²|G|²; "no split" takes λL = sup_m|G|; "exact split" replaces
  sup_m|S|²·meas(E_λ ∩ m) by the measured ∫_(E_λ∩m)|S|² (C = 1).

Joint bound B_K / I_M,K:

{_table(["(Q,c)", "K=0", "K=1", "K=2", "K=3", "K=4", "K=5", "K=6"], jrows)}

Smallest closing K (ratio < 0.99) at each pair; K ≤ 14, the value in brackets is the ratio there:

{_table(["(Q,c)", "joint", "sup", "Hölder", "theorem-shaped, split, C = 1 / 2 / 4", "literature form, C = 1 / 2 / 4", "no split, C = 1 / 2 / 4", "exact split, C = 1"], trows)}

Diagnostics at (150,20):

* the correlation ratio B_K/(P_m·⟨|G|^K⟩_m) is {", ".join(f"{a['corr_ratio'][K]:.3f}" for K in Ks)} for
  K = {", ".join(str(K) for K in Ks)}: |S| and |G| are nearly independent on the minor arcs for small K;
* the Rayleigh prediction (P_m/R2)·Γ(1+K/2)·L^(−K/2)/ρ_K gives {", ".join(f"{a['rayleigh_ratio'][K]:.3f}" for K in (1, 2, 3, 4))}
  for K = 1–4, against measured {", ".join(f"{a['joint_ratio'][K]:.3f}" for K in (1, 2, 3, 4))};
* the cost per extra power, (B_(K+1)/B_K)/(I_M,K+1/I_M,K), is {", ".join(f"{v:.2f}" for v in pp)} for K = 0…5,
  rising toward sup_m|G|/L = {a['sup_G'] / L:.3f};
* inputs of the theorem-shaped bound: L·ΣΛ²/I_M,2 = {L * res['L2'] / I2:.3f}, OD/I_M,2 = {e2['OD'] / I2:.3f}; optimal
  λ = {", ".join(f"{lam[C]:.2f} (C = {C})" for C in lam)};
* meas(E_λ ∩ m) = {", ".join(f"{_e(v, 1)} (λ = {l})" for l, v in a['measE'].items())};
* share of B_K from minor points within R/N of α = 0 or 1/2: {"; ".join(f"R = {R}: " + ", ".join(f"{100 * a['shoulder_share'][R][K]:.1f}%" for K in (1, 2, 6, 10, 14)) for R in (60, 200, 1000))}
  for K = 1, 2, 6, 10, 14. For small K the weight is spread over the whole minor arc; only for K ≳ 10 do the shoulders
  of the arcs at 0 and 1/2 carry a large part of it;
* **sparse-target control**: 19 random even shifts t_k, scale-matched to 2^k (t_k from (2^(k−1), 2^(k+1)]), 8 seeds:
  B_1/I_M,1 = {np.mean(sp[(150, 20)]):.3f} (max {np.max(sp[(150, 20)]):.3f}) at (150,20) against
  {a['joint_ratio'][1]:.3f} for powers of 2; {np.mean(sp[(20, 10)]):.3f} (max {np.max(sp[(20, 10)]):.3f}) at (20,10)
  against {b['joint_ratio'][1]:.3f}.

**Reading.** At N = 10^6 the magnitude inequality closes with **one** power of 2 at every Q ≥ 20 (joint bound). The
sparse control shows this is generic averaging over 19 targets of the right sizes, not something special about powers of
2. In the Rayleigh model the K = 1 ratio scales like (P_m/R2)/√L, i.e. like √(log N), so K = 1 is a finite-N phenomenon.
Replacing exact quantities by what proofs can use moves K\\* to {kstr(out[(150, 20)]['split'][1])} /
{kstr(out[(150, 20)]['split'][2])} / {kstr(out[(150, 20)]['split'][4])} at (150,20) for C = 1 / 2 / 4 and to
{kstr(out[(20, 10)]['split'][1])} / {kstr(out[(20, 10)]['split'][2])} / {kstr(out[(20, 10)]['split'][4])} at (20,10). The
theorems' K also pays for major-arc and asymptotic losses (power-size moduli, zero-density estimates) that one X cannot
show, so no K in these tables is comparable to a theorem's K; that a 7 or an 8 appears is a coincidence. (The sup bound
at (20,10) reaches {sup_s:.3f} at K = 14.)

**Literature.** Linnik proved, on GRH (1951) and unconditionally (1953), that some absolute K makes every large even
N = p₁ + p₂ + 2^(v₁) + … + 2^(v_K); he gave no value. Gallagher (1975) gave a simpler proof. Heath-Brown and Puchta (2002)
obtained K = 13 unconditionally and 7 on GRH; Pintz and Ruzsa obtained 7 on GRH (2003) and **K = 8 unconditionally**
(2020), the best published value.
"""
    return txt, out


def sec_e3(res) -> tuple[str, dict]:
    e3 = {s["name"]: s for s in res["e3"]}
    names = [f"Omega<={r}" for r in range(1, 7)] + ["odd Omega<=2", "Omega<=2, (m,10)=1", "Chen-shaped: Omega<=2, (m,6)=1",
                                                     "odd Omega<=3", "Omega<=3, (m,10)=1", "101-rough", "1000-rough"]
    names += [n for n in e3 if n.startswith("Bernoulli")]
    nice = lambda n: (n.replace("Omega<=", "Ω ≤ ").replace("Chen-shaped: ", "Chen-shaped ")  # noqa: E731
                      .replace("[control for Ω ≤ 2]", "control").replace("[control for Ω ≤ 3]", "control"))
    rows = []
    for nme in names:
        s = e3[nme]
        rows.append([nice(nme), f"{s['density']:.4f}"]
                    + [f"{r3(s['rows'][p]['joint_ratio'])} / {r3(s['rows'][p]['CS_ratio'])}" for p in SHORT5]
                    + [f"{s['rows'][(150, 20)]['f']:.3f}"])
    rstar = []
    out = {}
    for p in SHORT5:
        j = [e3[f"Omega<={r}"]["rows"][p]["joint_ratio"] for r in range(1, 7)]
        c = [e3[f"Omega<={r}"]["rows"][p]["CS_ratio"] for r in range(1, 7)]
        rj, rc = first_close(j), first_close(c)
        zs = sorted(int(n.split("-")[0]) for n in e3 if n.endswith("-rough"))
        zz = {}
        for key in ("joint_ratio", "CS_ratio"):
            last = None
            for z in zs:
                v = e3[f"{z}-rough"]["rows"][p][key]
                if v < 0.99:
                    last = (z, v)
                else:
                    nxt = (z, v)
                    break
            zz[key] = (last, nxt)
        out[p] = dict(rj=None if rj is None else rj + 1, rc=None if rc is None else rc + 1, z=zz)
        fz = lambda t: f"{t[0][0]} ({t[0][1]:.3f}); next z = {t[1][0]}: {t[1][1]:.3f} {lab(t[1][1])}" if t[0] else "none"  # noqa: E731
        rstar.append([pair(p), f"{out[p]['rj']} ({j[rj]:.3f})" if rj is not None else "> 6",
                      f"{out[p]['rc']} ({c[rc]:.3f})" if rc is not None else "> 6",
                      fz(zz["joint_ratio"]), fz(zz["CS_ratio"])])
    zs = sorted(int(n.split("-")[0]) for n in e3 if n.endswith("-rough"))
    zrows = [[f"{z}", f"{e3[f'{z}-rough']['density']:.4f}"]
             + [f"{r3(e3[f'{z}-rough']['rows'][p]['joint_ratio'])} / {r3(e3[f'{z}-rough']['rows'][p]['CS_ratio'])}"
                for p in ((150, 20), (50, 10), (20, 10))] for z in zs]
    a2 = e3["Omega<=2"]["rows"][(150, 20)]
    bern = next(e3[n] for n in e3 if n.startswith("Bernoulli") and "Omega<=2" in n)
    chen = e3["Chen-shaped: Omega<=2, (m,6)=1"]["rows"]
    chen_close = [pair(p) for p in SHORT5 if chen[p]["joint_ratio"] < 0.99]
    txt = f"""## 5. E3 — relax one prime to an almost-prime (the Rényi / Chen direction)

R^(A)(N) = Σ_(n+m=N) Λ(n)·1_A(m) = ∫S·T_A·e(−Nα), T_A = Σ_(m∈A) e(mα) (unweighted, 2 ≤ m ≤ X), at N = X. Minor-arc
magnitude bounds: **joint** ∫_m|S||T_A| and **Cauchy–Schwarz** √(P_m·∫_m|T_A|²); f = ∫_m|T_A|²/|A| is the share of T_A's
energy left on the minor arcs. Sets: A_r = {{Ω(m) ≤ r}} (A_1 is exactly the primes, 78,498 of them, unweighted);
odd A_r (an m = N − p with p odd is odd, so even members are dead mass); A_r coprime to 10 (N = 2⁶·5⁶);
**Chen-shaped** = Ω ≤ 2 with (m,6) = 1, i.e. all prime factors ≥ N^(1/10) ≈ 3.98; z-rough B_z = no prime factor < z;
and Bernoulli random sets with the densities of A_2 and A_3 (controls). Cells: joint / CS, divided by the major part.

{_table(["set", "density", "(150,20)", "(100,20)", "(50,10)", "(20,10)", "(5,5)", "f at (150,20)"], rows)}

Smallest closing r and largest closing z (finite-X statistics; z runs over every prime 3…101, then 317, 1000):

{_table(["(Q,c)", "r\\* joint", "r\\* CS", "z\\* joint", "z\\* CS"], rstar)}

<details><summary>z-sweep (joint / CS)</summary>

{_table(["z", "density", "(150,20)", "(50,10)", "(20,10)"], zrows)}

</details>

**Reading.** With one prime replaced by a member of A_2, the joint bound closes at (150,20) ({a2['joint_ratio']:.3f});
with A_3 or A_4 it closes at every short-list pair except (5,5) at r = 3. A random set of the density of A_2 never closes
({bern['rows'][(150, 20)]['joint_ratio']:.2f} at (150,20)): its energy is sea from the start (f =
{bern['rows'][(150, 20)]['f']:.3f} against {a2['f']:.3f} for Ω ≤ 2, matching the random prediction (1 − δ)(1 − |M|)). So
what closes the bound is the arithmetic concentration of the second family on the major arcs, plus density. For z ≤ 7,
T_z is exactly a sum of Dirichlet kernels at j/d, d | P(z) ∈ {{2, 6, 30}}, which the arcs swallow once Q ≥ d (z = 5 at
Q = 5 shows the exception: the spike at 1/6 sits on the minor arcs and the CS bound fails). As z grows the lines spread to
larger q, and at z ≈ √X the set is the primes again (1000-rough reproduces Ω ≤ 1). The Chen-shaped set closes only at
{", ".join(chen_close) if chen_close else "no short-list pair"}.

**r\\* and z\\* reflect how dense these sets are at X = 10^6, not the arithmetic of almost-primes.** For fixed r the density
of A_r tends to 0 like (log log X)^(r−1)/((r−1)!·log X), so every fixed-r ratio eventually fails. Do not read
r\\* = 2 at (150,20) as Chen's 2. The theorems in this direction are sieve theorems with no circle-method inequality:
Rényi (1948) proved N = p + P_K for some absolute K with Brun's sieve and Linnik's large sieve; Chen (announced 1966,
published 1973) proved N = p + P_2 with a weighted lower-bound sieve, Bombieri–Vinogradov at level 1/2 and the switching
principle. The parity barrier is why the sieve stops at P_2.
"""
    return txt, out


def sec_e4(res) -> str:
    e4 = res["e4"]
    e1 = res["e1"]
    R2 = res["R2"]
    X = res["X"]
    Nt = res["nters"]
    rows = e4["rows"]
    brow = []
    for p in res["pairs"]:
        r = rows[p]
        a, q, nb = r["sup_at"]
        brow.append([pair(p), f"{e1[p]['P_m'] / R2:.3f}", f"{r['P_m_over_R']:.3f}", f"{r['I_m_over_R']:+.4f}",
                     f"{a}/{q}", _e(r["sup"], 3)])
    trow = []
    for Q in (5, 10, 20, 50, 100, 150):
        cells = [f"{Q}"]
        for Np in (Nt[0], Nt[2]):
            cells.append(r3(rows[(Q, 5)][f"B3_eta_{Np}"]))
            cells.append(" / ".join(r3(rows[(Q, c)][f"B3_sharp_{Np}"]) for c in (5, 10, 20)))
        trow.append(cells)
    cs = e4["constants"]
    ter = e4["ternary"]
    side = e4["sidelobes"]
    srow = [["max_(ξ>c) |η̂(ξ)|/η̂(0)"] + [_e(side[c]["max_sidelobe"], 1) for c in (5, 10, 20)],
            ["sharp sidelobe 1/(πc)"] + [_e(side[c]["sharp_sidelobe"], 1) for c in (5, 10, 20)],
            ["binary-kernel tail outside c, η"] + [_e(side[c]["kernel_tail"], 1) for c in (5, 10, 20)],
            ["binary-kernel tail outside c, sharp ≈ 1/(π²c)"] + [_e(side[c]["sharp_tail"], 1) for c in (5, 10, 20)]]
    pe = [rows[p]["P_m_over_R"] for p in res["pairs"]]
    return f"""## 6. E4 — control: a smooth weight cures the arc-edge leakage, not the wall

w(n) = Λ(n)η(n/X) with η(x) = exp(1 − 1/(1 − (2x−1)²)) on (0,1): a smooth compactly supported weight (standard practice;
Helfgott's ternary proof uses weights of a different form). Constants: ∫η = {cs['int_eta']:.7f},
∫η² = {cs['int_eta2']:.7f}, (η∗η∗η)(1) = {cs['eta3_at_1']:.7f}, (η∗η∗η)(1.5) = {cs['eta3_at_1.5']:.7f};
ψ_η = {_e(e4['psi_eta'], 4)}, Σw² = {_e(e4['sum_w2'], 4)}. Binary: R_η(X) = {_f(e4['R_eta'], 1)}, main term
SS2(X)·X·∫η² = {_f(e4['main_bin'], 1)}. Ternary: B3 = sup_m|S_w|·Σw² against SS3(N′)·X²·(η∗η∗η)(N′/X) (sharp:
SS3(N′)·X²·vol(N′/X)). At the same spike, η pays the geometric factor ∫η·∫η²/((η∗η∗η)(t)/vol(t)) =
{ter[Nt[0]]['geom_factor']:.3f} at t = 1 and {ter[Nt[2]]['geom_factor']:.3f} at t = 1.5 relative to the sharp cutoff;
that is bookkeeping and must be divided out before crediting smoothing.

**Binary wall** (P_m / R, sharp against η), and where sup_m|S_η| sits:

{_table(["(Q,c)", "sharp P_m/R2", "η: P_m/R_η", "η: I_m/R_η", "sup_m|S_η| at", "sup_m|S_η|"], brow)}

**Ternary B3 / main.** η's sup does not depend on c (it is always an unresolved spike a/q with q > Q), so one η column
per N′; sharp for c = 5 / 10 / 20:

{_table(["Q", f"η, N′ = {Nt[0]:,}", f"sharp c = 5 / 10 / 20, N′ = {Nt[0]:,}", f"η, N′ = {Nt[2]:,}", f"sharp c = 5 / 10 / 20, N′ = {Nt[2]:,}"], trow)}

**Sidelobes** of the window function (|η̂| against the Dirichlet kernel of the sharp cutoff):

{_table(["", "c = 5", "c = 10", "c = 20"], srow)}

**Reading.** Smoothing removes the c = 5 plateau of the sharp cutoff (1.007 at N′ = X−1 for every Q ≥ 50, which is
the arc-edge sidelobe at 0/1): with η, (150,5) drops to {rows[(150, 5)][f'B3_eta_{Nt[0]}']:.3f} at X−1 and to
{rows[(150, 5)][f'B3_eta_{Nt[2]}']:.3f} at 1.5X−1. At X ± 1 it pays the factor {ter[Nt[0]]['geom_factor']:.2f}
everywhere else, so (50,·) fails with η ({rows[(50, 10)][f'B3_eta_{Nt[0]}']:.3f}) where sharp (50,10) closes. Moving
N′ to 1.5X helps both weights through the geometry of the singular integral and is not credited to smoothing. The
sidelobes do not vanish: at the c = 5 edge they are about 1% of the peak. The binary P_m/R_η stays at
{min(pe):.2f}–{max(pe):.2f}: the main term and Σw² both scale with ∫η², so the decoherent power is unchanged. Leakage is
an edge effect; the wall is a power budget.
"""


def captions(res, e0x, e1x, e2x, e3x) -> str:
    a = res["e0"][(150, 20)]
    R2 = res["R2"]
    e1 = res["e1"]
    n = e1[(150, 20)]["n_targets"]
    b = e1[(150, 20)]
    tr = e1[(223, 10)]["transmitted"]
    e2 = res["e2"]
    r = e2["rows"][(150, 20)]
    e3 = {s["name"]: s for s in res["e3"]}
    bern = next(e3[nm] for nm in e3 if nm.startswith("Bernoulli") and "Omega<=2" in nm)
    e4 = res["e4"]
    sc = [s for s in res["e5_scaling"] if s["c"] == 10]
    return f"""## Figures (captions in phasor language)

![Fig X0](figures/fig_X0_obstacle.png)

**Fig X0 — the obstacle.** Every minor-arc grid point adds one phasor S(α)²e(−Nα)dα to the binary count. Laid head to
tail, the {a['n_minor'] / 1e6:.1f} million phasors at (150,20) form a path {a['P_m'] / a['I_M']:.1f} times longer than the
radius of the circle they would have to leave to spoil the count, yet the walk never strays more than
{100 * a['max_abs_F'] / a['I_M']:.1f}% of the way out (inset ×50). A bound that knows only the step lengths must allow
the straight path, and every magnitude refinement in panel (b) is at least that long; even after the resolvable locked
channels are stripped off (Q′ ≈ {e0x[(150, 20)]['best']['Qp']}) a decoherent sea of
{e0x[(150, 20)]['best']['P_E'] / R2:.2f}·R2 remains. Panel (c) measures how much phase knowledge would close the gap:
window integrals about {e0x[(150, 20)]['dstar']:.0f}/N wide at (150,20), but about {e0x[(20, 10)]['dstar']:.0f}/N at
(20,10), where the left-out spikes cancel only across different a/q.

![Fig X1](figures/fig_X1_exceptional.png)

**Fig X1 — averaging over N.** Each even target N is a receiver: it demodulates S² with its own reference e(Nα), and
I_m(N) is its minor-arc output. Alone, a receiver must allow every minor-arc phasor to line up (P_m/I_M(N) =
{b['pm_ratio_min']:.2f}–{b['pm_ratio_max']:.2f} at (150,20)). Together, the receivers share one energy budget
Σ_N I_m(N)² = ∫_m|S|⁴, so the decoherent power cannot be in phase with many references at once: Chebyshev lets at most
{b['bound_4th']:,.0f} of the {n:,} targets fail at (150,20), the sup form (the ternary device) {b['bound_natural_sup']:,.0f}.
None fails, and the median |I_m|/I_M is {b['quantiles'][0.5]:.3f}. Panel (c): only about {100 * e1x['lo']:.0f}–{100 * e1x['hi']:.0f}%
of the budget lands on the targets; the rest goes to odd N and to N outside [X/2, X].

![Fig X5](figures/fig_X5_transmission.png)

**Fig X5 — single-N transmission.** The minor-arc energy P_m is the same for every target: the decoherent power of S²
does not know which reference will demodulate it. The signal does: at a spike a/q the reference samples the phase
e(−aN/q), and summed over the φ(q) spikes with denominator q that is the Ramanujan sum c_q(N) — φ(q) when q | N, −1 when
q is a prime not dividing N. So SS2(N) counts how many locked spikes the reference catches in phase. Multiples of
2·3·5·7·11 collect about 3.6× the signal of N = 2p against the same wall, and {len(tr)} of them get through at (223,10).
Panel (b): the wall need(X) rises at every step ({", ".join(f"{s['need']:.2f}" for s in sc)}) while the best available
signal climbs only in rare steps, so the transmitted counts (labels) thin out.

![Fig X2](figures/fig_X2_linnik.png)

**Fig X2 — Linnik.** A second bank of L = {e2['L']} oscillators at frequencies 2, 4, …, 2^{e2['L']}. Near α = 0 and 1/2
they all point the same way; over most of the circle they behave like a random walk (⟨|G|⟩ = {e2['G_moments'][1]:.2f},
∫|G|² = L). Weighting the minor-arc energy by |G|^K costs about {r['joint_ratio'][1] / r['joint_ratio'][0]:.2f} against the
main term for the first power and {r['joint_ratio'][2] / r['joint_ratio'][1]:.2f}–{r['joint_ratio'][6] / r['joint_ratio'][5]:.2f}
for each further one, enough for K = 1 at every Q ≥ 20 (panel a). As K grows, |G|^K moves the weight toward the
shoulders of the arcs at 0 and 1/2, where G's slow oscillators are still in phase and S still has Dirichlet sidelobes:
the share of B_K within 60/N of those two points is {100 * r['shoulder_share'][60][1]:.1f}% at K = 1,
{100 * r['shoulder_share'][60][6]:.0f}% at K = 6 and {100 * r['shoulder_share'][60][14]:.0f}% at K = 14 at (150,20)
(panel b), and the joint bound slowly turns into the sup bound. G is large only on a tiny part of the minor arcs
(panel c), which is what the theorem-shaped λ-split exploits.

![Fig X3](figures/fig_X3_almost_primes.png)

**Fig X3 — almost-primes.** Replace one prime by a member of a set whose oscillators lock into fewer, lower-q channels.
z-rough numbers are the primes with the sieve stopped at z: for z ≤ 7 their spectrum is a few Dirichlet lines at j/d,
d | 30, which the arcs swallow; as z grows the lines spread to larger q and build their own sea, and at z ≈ √X they are
the primes again. A random set of the same density never closes ({bern['rows'][(150, 20)]['joint_ratio']:.2f} at
(150,20)): its energy is sea from the start. What closes is arithmetic concentration on the major arcs, plus density.

![Fig X4](figures/fig_X4_window.png)

**Fig X4 — smooth weight (control).** Apodization tapers each oscillator's amplitude to zero at both ends of the
observation window. It removes the ringing of the rectangular window (the Dirichlet sidelobes, panel a) that held the
ternary bound at 1.007 for c = 5 at every Q ≥ 50 (panel b). It does not change the power of the decoherent sea: the binary
P_m/R_η stays at {min(v['P_m_over_R'] for v in e4['rows'].values()):.2f}–{max(v['P_m_over_R'] for v in e4['rows'].values()):.2f}
(panel c).
"""


def write(res, path: Path) -> None:
    R2 = res["R2"]
    e1 = res["e1"]
    n = e1[(150, 20)]["n_targets"]
    pmr = [e1[p]["P_m"] / R2 for p in res["pairs"]]
    s0 = certificates(res)
    s1, e0x = sec_e0(res)
    s2, e1x = sec_e1(res)
    s3 = sec_e5(res)
    s4, e2x = sec_e2(res)
    s5, e3x = sec_e3(res)
    s6 = sec_e4(res)
    a0 = res["e0"][(150, 20)]
    b1 = e1[(150, 20)]
    tr = e1[(223, 10)]["transmitted"]
    best_tr = min(t[1] for t in tr)
    r2 = res["e2"]["rows"][(150, 20)]
    e3 = {s["name"]: s for s in res["e3"]}
    bern = next(e3[nm] for nm in e3 if nm.startswith("Bernoulli") and "Omega<=2" in nm)
    e4 = res["e4"]
    sc10 = [s for s in res["e5_scaling"] if s["c"] == 10]
    head = [
        ["E0 obstacle", "how far is the binary wall from closing, for one N?",
         f"P_m/I_M = {a0['P_m'] / a0['I_M']:.2f} at (150,20); best refinement {e0x[(150, 20)]['bnd'][0]:.1f}·I_M",
         f"fails at every (Q,c); window phases over ≈ {e0x[(150, 20)]['dstar']:.0f}/N would suffice", "the open problem"],
        ["E1 average over N", f"how many of the {n:,} even N in [X/2, X] can the inequality lose?",
         f"F ≤ {b1['bound_4th']:,.0f} (true 4th moment), ≤ {b1['bound_natural_sup']:,.0f} (sup form); actual F = 0",
         f"closes for all but ≤ {100 * b1['bound_4th'] / n:.2f}% (sup form {100 * b1['bound_natural_sup'] / n:.1f}%)",
         "van der Corput, Chudakov, Estermann 1937–38"],
        ["E5 transmission", "does any single N get through with no relaxation?",
         f"{len(tr)} N at (223,10), best P_m/I_M(N) = {best_tr:.3f}",
         f"closes for N with SS2 ≈ 5; fraction {100 * sc10[0]['n_transmitted'] / sc10[0]['n_targets']:.1f}% at 10^4 → "
         f"{100 * sc10[-1]['n_transmitted'] / sc10[-1]['n_targets']:.4f}% at 2·10^6", "finite checks to 4·10^18"],
        ["E2 Linnik", "N = p + p′ + 2^(k₁) + … + 2^(k_K)",
         f"joint B_1/I_M,1 = {r2['joint_ratio'][1]:.3f}; theorem-shaped K\\* = "
         + " / ".join(kstr(e2x[(150, 20)]['split'][C]) for C in (1, 2, 4)),
         "closes at K = 1 (joint bound), generic averaging", "Linnik 1951/53; best published K = 8"],
        ["E3 almost-primes", "N = p + m with Ω(m) ≤ r",
         f"Ω ≤ 2: joint {e3['Omega<=2']['rows'][(150, 20)]['joint_ratio']:.3f} at (150,20); random set of equal density "
         f"{bern['rows'][(150, 20)]['joint_ratio']:.2f}", f"closes at r = {e3x[(150, 20)]['rj']}–4 depending on (Q,c)",
         "Rényi 1948, Chen 1973 (sieve theorems)"],
        ["E4 smooth weight", "does smoothing move the wall? (control)",
         f"binary P_m/R_η = {min(v['P_m_over_R'] for v in e4['rows'].values()):.2f}–"
         f"{max(v['P_m_over_R'] for v in e4['rows'].values()):.2f}",
         "binary wall unchanged; ternary c = 5 plateau removed", "standard practice"],
    ]
    txt = f"""# Relaxations under which the magnitude inequality closes at N = 10^6

> Exploration at one finite X, not a proof and not an attempt at one. Every number below is a floating-point
> computation on the grid M = 2^23 at X = N = 10^6 (plus the scaling runs); none is a theorem, and none of the
> experiments gets past the binary wall for a single typical N. Each one changes the question.

At N = 10^6 the magnitude-only bound fails for every tabulated (Q,c) (P_m/R2 = {min(pmr):.2f}–{max(pmr):.2f}, results.md).
Proving R2(N) > 0 for every large N this way would require showing that the minor-arc phases cancel uniformly in N; that
is the open problem. For a single N, R2(N) > 0 is a finite check (done to 4·10^18). This file does two things: it
measures the obstacle for one N (E0), and it tests, at the same X, the relaxations under which the circle method's
*magnitude* inequality does close — averaging over N (E1), adding powers of 2 (E2), relaxing one prime to an almost-prime
(E3) — plus the targets that get through with no relaxation (E5) and a control that should not help (E4).

**What "closes" means.** For a target N and arcs M(Q,c) (half-width c/N, exact integer membership, as in results.md),
R(N) = I_M + I_m. A bound on |I_m| that uses magnitudes only is computed and divided by I_M. Because the arc integrals
are Riemann sums, each ratio carries a decision margin: **closes** < 0.99, **undecided** 0.99–1.01, **fails** > 1.01.
"Closes" means the method's inequality alone gives R(N) > 0 for that target, numerically, at this X. Every K\\*, r\\* and
z\\* below is a finite-X statistic of this computation, not comparable to a constant in a theorem.

## Headline

{_table(["experiment", "question", "key numbers", "verdict at X = 10^6", "literature direction"], head)}

{s0}
{s1}
{s2}
{s3}
{s4}
{s5}
{s6}
{captions(res, e0x, e1x, e2x, e3x)}
## What this does and does not show

**Shows.**
* At one N, with numbers, that the binary wall is about phases, not sizes: the minor-arc walk is
  {a0['P_m'] / a0['I_M']:.1f}× the signal long and ends within {100 * abs(a0['I_m'].real) / a0['I_M']:.1f}% of the origin,
  and every magnitude-only refinement tried (Hölder, sup × L¹, model subtraction) is longer than the path (E0).
* That each classical relaxation turns the binary problem into one where a magnitude inequality closes at this X, and
  where the slack comes from: a shared energy budget across targets (E1), a second bank of oscillators that is mostly
  decoherent on the minor arcs (E2), a second family whose energy sits on the major arcs (E3).
* That a few very smooth N get through with no relaxation at all, and that this dies out as X grows (E5).
* That a smooth weight fixes the arc-edge leakage of the ternary bound but leaves the binary wall where it was (E4).
* Controls: random sets of the same density (E3) and scale-matched random shifts (E2) separate arithmetic from density
  and from generic averaging.

**Does not show.**
* Anything about a single large N: no experiment here bounds the minor-arc phases, and that is the open problem.
* Asymptotics: X = 10^6 is far from asymptotic, Q = 150 is simultaneously X^0.36 and (log X)^1.9, and I_M is computed
  from the true S, not from a provable major-arc model.
* That any K\\*, r\\* or z\\* matches a theorem: the theorems pay for major-arc and asymptotic losses that one X cannot show,
  and Rényi's and Chen's theorems are sieve theorems that do not use this inequality at all. r\\* = 2 is not Chen's 2.
* Trends: no fits are made; where a quantity is said to grow, the table shows each value.

## References

* J. G. van der Corput, Sur l'hypothèse de Goldbach pour presque tous les nombres pairs, Acta Arith. 2 (1937).
* N. G. Chudakov, On the density of the set of even numbers which are not representable as a sum of two odd primes,
  Izv. Akad. Nauk SSSR (1938).
* T. Estermann, On Goldbach's problem: proof that almost all even positive integers are sums of two primes,
  Proc. London Math. Soc. 44 (1938).
* H. L. Montgomery, R. C. Vaughan, The exceptional set in Goldbach's problem, Acta Arith. 27 (1975).
* Yu. V. Linnik, Prime numbers and powers of two, Trudy Mat. Inst. Steklov 38 (1951); Addition of prime numbers with
  powers of one and the same number, Mat. Sb. 32 (1953).
* P. X. Gallagher, Primes and powers of 2, Invent. Math. 29 (1975).
* D. R. Heath-Brown, J.-C. Puchta, Integers represented as a sum of primes and powers of two, Asian J. Math. 6 (2002).
* J. Pintz, I. Z. Ruzsa, On Linnik's approximation to Goldbach's problem I, Acta Arith. 109 (2003); II, Acta Math.
  Hungar. 161 (2020).
* A. Rényi, On the representation of an even number as the sum of a single prime and a single almost-prime number,
  Izv. Akad. Nauk SSSR Ser. Mat. 12 (1948).
* J. R. Chen, On the representation of a larger even integer as the sum of a prime and the product of at most two primes,
  Sci. Sinica 16 (1973).
* T. Oliveira e Silva, S. Herzog, S. Pardi, Empirical verification of the even Goldbach conjecture and computation of
  prime gaps up to 4·10^18, Math. Comp. 83 (2014).

## Reproduce

    python run_all.py              # the main project (writes the S caches this file reuses)
    python run_escapes.py          # E0–E5: ~{res['runtime_s'] / 60:.0f} min, peak ~{res['peak_rss_mb'] / 1024:.1f} GB
    python run_escapes.py --report # redraw figures/fig_X*.png and this file from cache/escapes.pkl
    python -m pytest tests/        # includes tests/test_escapes.py (small grid X = 10^4)
"""
    Path(path).write_text(txt)
