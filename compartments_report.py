"""Figure fig_W1_compartments.png and compartments.md for the body-water lens (every number read from the pickle)."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import figures as base  # noqa: F401  (shared rcParams / palette)
from escapes_report import _table
from figures import AQUA, BLUE, GRID, INK, INK2, MAGENTA, MUTED, ORANGE, RED, VIOLET, YELLOW, _save, plt

# Physiology, for comparison only (percent of body weight).  Fetal and infant values: Abdominal Key, "Body Composition in
# the Fetus and Newborn: Effects of Intrauterine Growth Aberration" (ch. 2 of Nephrology and Fluid/Electrolyte Physiology:
# Neonatology Questions and Controversies).  Adult: StatPearls, "Physiology, Water Balance" (70 kg man: 42 L water, 28 L
# intracellular, 14 L extracellular = 10.5 L interstitial + 3.5 L plasma).
PHYSIOLOGY = [
    {"stage": "first trimester", "TBW": 94, "ECW": None, "ICW": None, "solids": 6},
    {"stage": "24 weeks", "TBW": 86, "ECW": 60, "ICW": 26, "solids": 14},
    {"stage": "term", "TBW": 78, "ECW": 44, "ICW": 34, "solids": 22},
    {"stage": "1 year", "TBW": 70, "ECW": 28, "ICW": 42, "solids": 30},
    {"stage": "adult man", "TBW": 60, "ECW": 20, "ICW": 40, "solids": 40},
]
SOURCE_FETAL = ("https://abdominalkey.com/body-composition-in-the-fetus-and-newborn-effects-of-intrauterine-growth-aberration/")
SOURCE_ADULT = "https://www.ncbi.nlm.nih.gov/books/NBK541059/"


def _sh(x, tot):
    return x / tot


def figure(res, out: Path) -> None:
    out = Path(out)
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    # (a) physiology
    ax = axes[0, 0]
    xs = np.arange(len(PHYSIOLOGY))
    for i, p in enumerate(PHYSIOLOGY):
        if p["ECW"] is None:
            ax.bar(i, p["TBW"], 0.6, color=YELLOW, alpha=0.45, edgecolor="white")
            ax.text(i, p["TBW"] / 2, "water\n(not split)", ha="center", va="center", fontsize=7)
            ax.bar(i, p["solids"], 0.6, bottom=p["TBW"], color=MUTED, edgecolor="white")
        else:
            ax.bar(i, p["ECW"], 0.6, color=YELLOW, edgecolor="white", label="extracellular water" if i == 1 else None)
            ax.bar(i, p["ICW"], 0.6, bottom=p["ECW"], color=VIOLET, edgecolor="white", label="intracellular water" if i == 1 else None)
            ax.bar(i, p["solids"], 0.6, bottom=p["TBW"], color=MUTED, edgecolor="white", label="solids" if i == 1 else None)
        ax.text(i, 101, f"water {p['TBW']}%", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(xs, [p["stage"] for p in PHYSIOLOGY], fontsize=8)
    ax.set_ylim(0, 110)
    ax.set_ylabel("% of body weight")
    ax.set_title("(a) physiology, for comparison: body water by age (Abdominal Key ch. 2; StatPearls); colours not paired with (b)")
    ax.legend(loc="lower right", fontsize=7)
    # (b) development of S's compartments with X
    ax = axes[0, 1]
    dev = res["development"]
    X = np.array([r["X"] for r in dev], dtype=float)
    sh = lambda key, sub=None: np.array([(r[sub][key] if sub else r[key]) / r["L2"] for r in dev])  # noqa: E731
    for key, col, lab in (("A", RED, "atoms (q = 1, 2)"), ("F", BLUE, "forest (3 ≤ q ≤ q*)"), ("sea", AQUA, "sea (the rest)")):
        ax.plot(X, sh(key, "local_centre"), "o", color=col, ms=6, label=f"{lab}, measured")
    pa = np.array([r["pred"]["atoms"] for r in dev])
    pf = np.array([r["pred"]["forest"] for r in dev])
    ps = np.array([r["pred"]["sea"] for r in dev])
    for y, col in ((pa, RED), (pf, BLUE), (ps, AQUA)):
        ax.plot(X, y, color=col, lw=1.0)
    ax.plot([], [], color=INK2, lw=1.0, label="no-fit model (Ramanujan sums, Farey count)")
    Xf = np.logspace(4, 40, 200)
    asy = [_asym(x) for x in Xf]
    for i, col in zip(range(3), (RED, BLUE, AQUA)):
        ax.plot(Xf, [a[i] for a in asy], color=col, lw=0.8, ls=(0, (4, 2)))
    ax.plot([], [], color=INK2, lw=0.8, ls=(0, (4, 2)), label="leading-order law, extended to X = 10⁴⁰")
    ax.axhline(0.5, color=INK, lw=0.8, ls=(0, (1, 1.5)))
    ax.text(2e4, 0.51, "limit ½", fontsize=7, color=INK2)
    ax.set_xscale("log")
    ax.set_ylim(0, 0.62)
    ax.set_xlabel("X (the 'growth' parameter; not time)")
    ax.set_ylabel("share of ΣΛ² (Parseval: the conserved total)")
    ax.set_title("(b) the energy of S by compartment, X = 10⁴ … 2·10⁶, and its X → ∞ limit (0, ½, ½)")
    ax.legend(loc="center right", fontsize=6.6)
    # (c) the sea is thinner near the forest
    ax = axes[1, 0]
    row = next(r for r in dev if r["X"] == 10**6)
    bands = [b for b in row["sea_bands"] if b["q"][0] >= 31]
    xb = np.arange(len(bands))
    wdt = 0.2
    ax.bar(xb - 1.5 * wdt, [b["eff"] / row["X"] for b in bands], wdt, color=BLUE, label="measured in the windows (reading − capture × spike)")
    ax.bar(xb - 0.5 * wdt, [b["centre"] / row["X"] for b in bands], wdt, color=AQUA, label="measured at the spike centres")
    ax.bar(xb + 0.5 * wdt, [b["fg"] / row["X"] for b in bands], wdt, color=ORANGE, label="X·FG(q): variance of primes in progressions (no fit)")
    ax.axhline(row["rho"] / row["X"], color=INK, lw=1.0, ls=(0, (4, 2)), label="uniform sea (what plain dilution assumes)")
    ax.set_xticks(xb, [f"q ∈ [{b['q'][0]}, {b['q'][1]}]" for b in bands], fontsize=8)
    ax.set_ylabel("sea level per unit α, in units of X")
    ax.set_title("(c) the sea is thinner near the forest (X = 10⁶): the uniform-mixing assumption fails")
    ax.legend(loc="upper left", fontsize=7)
    # (d) growth aberration: four sets
    ax = axes[1, 1]
    ab = res["aberration"]
    y = np.arange(len(ab))[::-1]
    for yi, s in zip(y, ab):
        T = s["L2"]
        A, F = s["A"] / T, max(s["F"], 0.0) / T
        ax.barh(yi, A, 0.55, color=RED, edgecolor="white")
        ax.barh(yi, F, 0.55, left=A, color=BLUE, edgecolor="white")
        ax.barh(yi, 1 - A - F, 0.55, left=A + F, color=AQUA, edgecolor="white")
        ax.plot([s["pred_A"] / T, (s["pred_A"] + s["pred_F"]) / T], [yi, yi], "|", color=INK, ms=16, mew=1.6)
        ax.text(1.01, yi, f"forest {F:.1%}", va="center", fontsize=7.5)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    hd = [Patch(color=RED, label="atoms"), Patch(color=BLUE, label="forest"), Patch(color=AQUA, label="sea"),
          Line2D([], [], marker="|", color=INK, ms=10, mew=1.6, ls="none", label="no-fit model (compartment boundaries)")]
    ax.legend(handles=hd, loc="upper center", ncol=4, fontsize=7)
    ax.set_yticks(y, [s["name"] for s in ab], fontsize=8)
    ax.set_ylim(-0.6, len(ab) - 0.1)
    ax.set_xlim(0, 1.18)
    ax.set_xlabel("share of the set's own ΣΛ² (uniform-sea reading)")
    ax.set_title("(d) 'growth aberration': less congruence structure, less forest (X = 10⁶)")
    fig.suptitle("Fig W1. The energy of S as body-water-style compartments: one conserved total, nested indicators, capture corrections",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_W1_compartments.png")


def _asym(X: float, lnX: float | None = None) -> tuple:
    """Leading-order shares at X (or at ln X = lnX, for X beyond floats) with q* from q = sqrt(X / FG(q)) by fixed point
    (FG with the prime-q form ln q - gamma - ln 2 pi) and L_N ~ ln X - 1."""
    lnX = math.log(X) if lnX is None else lnX
    lq = 0.5 * lnX
    for _ in range(60):
        fg = max(lq - 0.5772156649 - math.log(2 * math.pi), 0.5)
        lq = 0.5 * (lnX - math.log(fg))
    LN = lnX - 1.0
    f = (lq + 1.3326 - 2.0) / LN
    return 2.0 / LN, f, 1.0 - 2.0 / LN - f


def write(res, path: Path) -> None:
    dev = res["development"]
    r6 = next(r for r in dev if r["X"] == 10**6)
    T6 = r6["L2"]
    lc = r6["local_centre"]
    lf = r6["local_fg"]
    p6 = r6["pred"]
    pct = lambda x: f"{100 * x:.2f}%"  # noqa: E731

    cert = _table(["check", "result", "detail"], [[n, "pass" if ok else "FAIL", d] for n, ok, d in res["certificates"]])
    pq = _table(["q", "window energy / (capture × model)", "window capture, min–max"],
                [[r["q"], f"{r['ratio']:.4f}", f"{r['cap_min']:.4f}–{r['cap_max']:.4f}"] for r in r6["per_q"]])
    b1 = _table(["indicator", "what it reads", "share of ΣΛ²"], [
        ["whole circle (Parseval)", "the conserved total ΣΛ²", "100%"],
        ["windows at q = 1, 2", "raw reading", pct(r6["E_A"] / T6)],
        ["windows at 3 ≤ q ≤ q*", "raw reading", pct(r6["E_F"] / T6)],
    ])
    b1s = _table(["compartment", "one average capture 0.7737, uniform sea", "exact capture, uniform sea", "exact capture, sea at X·FG(q)",
                  "exact capture, sea measured at the centres", "no-fit model"], [
        ["atoms (q = 1, 2)", pct(r6["naive_A"] / T6), pct(r6["A"] / T6), pct(lf["A"] / T6), pct(lc["A"] / T6), pct(p6["atoms"])],
        ["forest (3 ≤ q ≤ q*)", pct(r6["naive_F"] / T6), pct(r6["F"] / T6), pct(lf["F"] / T6), pct(lc["F"] / T6), pct(p6["forest"])],
        ["sea (the rest)", pct(1 - (r6["naive_A"] + r6["naive_F"]) / T6), pct(r6["sea"] / T6), pct(lf["sea"] / T6),
         pct(lc["sea"] / T6), pct(p6["sea"])],
    ])
    bands = [b for b in r6["sea_bands"] if b["q"][0] >= 31]
    bt = _table(["windows of", "# q", "measured in the windows", "measured at the centres", "X·FG(q) (no fit)", "uniform sea"],
                [[f"q ∈ [{b['q'][0]}, {b['q'][1]}]", b["n_q"], f"{b['eff'] / r6['X']:.2f} X", f"{b['centre'] / r6['X']:.2f} X",
                  f"{b['fg'] / r6['X']:.2f} X", f"{r6['rho'] / r6['X']:.2f} X"] for b in bands])
    dt = _table(["X", "L_N", "q*", "atoms: measured / model", "forest: measured / model", "sea: measured / model",
                 "forest with a uniform sea"],
                [[f"{r['X']:.0e}".replace("e+0", "e"), f"{r['pred']['L_N']:.2f}", f"{r['q_star']:.0f}",
                  f"{r['local_centre']['A'] / r['L2']:.4f} / {r['pred']['atoms']:.4f}",
                  f"{r['local_centre']['F'] / r['L2']:.4f} / {r['pred']['forest']:.4f}",
                  f"{r['local_centre']['sea'] / r['L2']:.4f} / {r['pred']['sea']:.4f}", f"{r['F'] / r['L2']:.4f}"] for r in dev])
    asy = _table(["X", "atoms", "forest", "sea"],
                 [[f"10^{e}", *[f"{v:.3f}" for v in _asym(0.0, e * math.log(10))]] for e in (6, 12, 20, 40, 100, 1000)])
    derived = {("1 year", "ECW"), ("first trimester", "solids"), ("adult man", "solids")}
    cell = lambda p, k: ("—" if p[k] is None else f"{'≈ ' if (p['stage'], k) in derived else ''}{p[k]}%")  # noqa: E731
    phys = _table(["stage", "total water", "extracellular", "intracellular", "solids"],
                  [[p["stage"], cell(p, "TBW"), cell(p, "ECW"), cell(p, "ICW"), cell(p, "solids")] for p in PHYSIOLOGY])
    ab = res["aberration"]
    abt = _table(["set (X = 10^6, same total weight ≈ ψ)", "what it keeps", "atoms: measured / model", "forest: measured / model",
                  "sea", "# denominators with spikes ≤ q*"],
                 [[s["name"], s["desc"], f"{s['A'] / s['L2']:.4f} / {s['pred_A'] / s['L2']:.4f}",
                   f"{max(s['F'], 0) / s['L2']:.4f} / {s['pred_F'] / s['L2']:.4f}", f"{s['sea'] / s['L2']:.4f}", s["n_qF"]] for s in ab])
    a0, a1 = dev[0], dev[-1]
    errs = {r["X"]: max(abs(r["local_centre"]["A"] / r["L2"] / r["pred"]["atoms"] - 1),
                        abs(r["local_centre"]["F"] / r["L2"] / r["pred"]["forest"] - 1)) for r in dev}
    big = [x for x in errs if x >= 3 * 10**4]
    err_big, err_small = 100 * max(errs[x] for x in big), 100 * errs[dev[0]["X"]]
    fg_gap = [1 - b["fg"] / b["eff"] for b in bands]
    coh0 = (a0["local_centre"]["A"] + a0["local_centre"]["F"]) / a0["L2"]
    coh1 = (a1["local_centre"]["A"] + a1["local_centre"]["F"]) / a1["L2"]
    text = f"""# Compartments: the body-water lens on the energy of S at N = 10^6

> Exploration at finite X, not a proof. Every number below is read from `cache/compartments.pkl` (grids M = 2^17 … 2^24,
> X = 10^4 … 2·10^6). Body-water physiology is used as a teacher for bookkeeping, not as a mechanism: nothing in S moves,
> mixes or equilibrates.

**What physiology does.** Total body water (TBW) is about 60% of an adult man's weight (42 L of 70 kg). Two thirds of it
is intracellular (28 L) and one third extracellular (14 L: 10.5 L interstitial and 3.5 L plasma) (StatPearls, *Physiology,
Water Balance*). Compartments are measured by **indicator dilution**: give a known amount I of an indicator that spreads
only through the compartment of interest, let it equilibrate, read its plasma concentration C, and set V = I/C. Heavy
water (D₂O, H₂¹⁸O) or antipyrine measures TBW; bromide, sucrose or inulin measures extracellular water; Evans blue
measures plasma. The rest follows by subtraction: intracellular = TBW − extracellular, interstitial = extracellular −
plasma, solids = body weight − TBW (Abdominal Key, ch. 2). The composition changes with growth. TBW falls from 94% of body
weight in the first trimester to 86% at 24 weeks, 78% at term and about 70% at one year, as solids accrete. Extracellular
water falls (60% → 44% → ≈ 28%) and intracellular water rises (26% → 34% → 42%). It also depends on growth conditions:
infants small for gestational age have more water and less fat; infants large for gestational age have more fat and
less water.

**What the circle has.** One conserved total, by Parseval: (1/M)Σ_k |S(k/M)|² = ΣΛ² exactly. The "indicators" are windows
|N(α − a/q)| ≤ ½ around the reduced fractions a/q of a set of denominators. A window "reads" the energy it covers; it
catches a known share (its capture) of the Dirichlet spike centred in it, and it also covers some sea. Three compartments:

- **atoms**: the spikes at q = 1, 2 (α = 0 and ½);
- **forest**: the spikes with 3 ≤ q ≤ q*(X), where q* is the denominator at which spike and sea are equal (fractals.md, M1);
- **sea**: the rest, i.e. the structureless part plus the spikes with q > q* that it drowns.

The no-fit spike model is S(a/q + β) ≈ μ(q)(ψ/φ(q))·v(β)/X, with v the Dirichlet kernel of energy X. So the φ(q) spikes of
denominator q carry ψ²μ²(q)/(Xφ(q)) together, and the shares are atoms = 2/L_N and forest = (Σ_{{3≤q≤q*}} μ²/φ)/L_N, with
L_N = XΣΛ²/ψ². **L_N is a dilution volume**: the total energy divided by the energy of one reference spike (the one at 0,
ψ²/X). At X = 10^6, L_N = {p6['L_N']:.2f}: the energy of S is spread as if over {p6['L_N']:.1f} copies of its main spike.

## 0. Certificates

{cert}

Each window's energy against its capture times the model spike, at X = 10^6 (small q, where the sea is negligible):

{pq}

## B1. Compartments by indicator dilution (X = 10^6)

{b1}

The readings are turned into compartments by the three-indicator solution. Each window reads its capture of its own
spikes plus the sea it covers; the rest of the circle holds the spikes' leakage plus the remaining sea:

{b1s}

Two corrections matter, and both have close counterparts in tracer practice:

1. **The indicator's distribution space.** A window's capture depends on where the fraction falls between grid points: at
   X = 10^6 a window holds 8 or 9 grid points and catches 0.755 or 0.803 of its spike, not the average 0.7737. Using one
   average capture misreads the atoms by {100 * (r6['naive_A'] / r6['A'] - 1):+.1f}%. The same kind of error happens when a
   tracer reaches more than its compartment (in standard tracer practice, not from the two sources above: heavy water
   exchanges with non-aqueous hydrogen, and bromide enters cells a little). This was also the A1 octave bug that the fractal-lens review caught.
2. **The equilibration assumption.** Plain dilution assumes the sea is uniformly mixed. It is not (next table). With a
   uniform sea the forest reads {100 * (r6['F'] / T6 / p6['forest'] - 1):+.1f}% against the model. With the sea level
   measured at the spike centres it reads {100 * (lc['F'] / T6 / p6['forest'] - 1):+.1f}%; with the no-fit X·FG(q) it reads
   {100 * (lf['F'] / T6 / p6['forest'] - 1):+.1f}%.

**The sea is thinner near the forest.** Sea level per unit α in the windows of each band of denominators, at X = 10^6:

{bt}

Near a fraction with denominator q, the structureless part of S is set by the variance of the primes among the residue
classes mod q: about X·FG(q) = X(ln q − γ − ln 2π − Σ_{{p|q}} ln p/(p−1)) (Friedlander–Goldston form of the Barban–Davenport–Halberstam /
Montgomery–Hooley variance). That is well below the circle's average sea level. The level climbs back toward the
average as q grows, and X·FG(q) falls {100 * min(fg_gap):.0f}–{100 * max(fg_gap):.0f}% below the window values. (Bands with q ≤ 30 are left out: there the sea is too
small a fraction of the window to measure.)

## B2. Development: the composition of S as X grows

{dt}

For comparison, the body (percent of body weight; fetal and infant values from the Abdominal Key chapter, adult from
StatPearls; "≈" marks a value obtained by subtraction):

{phys}

The bookkeeping matches closely: one conserved total redistributing among compartments as a growth parameter
increases, toward a limit. The **atoms** fall like 2/L_N ≈ 2/(ln X − 1), from {dev[0]['local_centre']['A'] / dev[0]['L2']:.1%}
at X = 10^4 to {dev[-1]['local_centre']['A'] / dev[-1]['L2']:.1%} at 2·10^6. Their energy 2ψ²/X grows only like X, while
the total grows like X ln X. The **forest** stays near 0.44. The **sea** grows. As X → ∞ the shares tend to **0, ½, ½**:
q* grows like √X up to logarithms, so the forest holds (ln q* + 1.33 − 2)/(ln X − 1) → ½ of the energy. The visible forest
ends up holding half the logarithm. The approach is slow (leading-order law):

{asy}

**Where the analogy stops.** The coherent share of S (atoms plus forest) *falls* with X, from {coh0:.1%} to {coh1:.1%}
(limit ½), while a fetus *gains* solids. Which compartment of S "corresponds" to extracellular or intracellular water is
a choice with no physical basis. Only the bookkeeping and the derived limit are literal. X is a parameter, not time, and
there is no counterpart to birth or to the postnatal contraction of extracellular water.

## B3. "Growth aberration": the arithmetic of the set sets its composition (X = 10^6)

Four sets with the same total weight (≈ ψ) but different congruence structure. The forest windows run to the primes' q*
({res['aberration_Q']}), placed wherever each set's own no-fit model has spikes. That model is the set's weight spread
evenly over the residue classes it allows mod lcm(D, q), summed with e(ab/q).

{abt}

Less congruence structure means less forest. Cramér's model has one spike (q = 1). Restricting it to odd numbers adds the
spike at ½, and nothing more. The primes have the whole Farey forest. Keeping only the primes ≡ 1 (mod 4) adds spikes at
q = 4·(odd squarefree): {ab[-1]['n_qF']} denominators instead of {ab[-2]['n_qF']}. The forest grows to
{max(ab[-1]['F'], 0) / ab[-1]['L2']:.1%} of a total energy that has doubled, so the atoms' share halves. These forest
readings use the uniform sea, so they share the primes' {100 * (r6['F'] / T6 / p6['forest'] - 1):+.1f}% bias; the
comparison between sets is unaffected. In body-composition terms (SGA: more water, less fat; LGA: more fat, less water)
the sea plays "water" and the forest plays "solids". That pairing is a reading device: the measurement is literal, the
nutrition language is not.

## Literal versus analogy

{_table(["body water", "on the circle", "status"], [
        ["total body water, conserved and split into compartments", "Parseval: ΣΛ² split into atoms, forest and sea", "LITERAL bookkeeping"],
        ["dilution V = I/C", "compartment = (window reading − sea covered)/capture; L_N = total/(energy of one spike)", "LITERAL algebra"],
        ["indicators with different spaces (D₂O, bromide/inulin, Evans blue)", "windows around different sets of denominators; nested subtraction",
         "LITERAL nesting"],
        ["indicator corrections (tracer outside its space)", "per-window capture 0.755/0.803 instead of 0.7737; otherwise "
         f"{100 * (r6['naive_A'] / r6['A'] - 1):+.1f}% on the atoms", "LITERAL, same kind of correction"],
        ["equilibration (the tracer mixes evenly)", "uniform sea: fails near rationals, where the sea is X·FG(q)", "MEASURED failure"],
        ["development: TBW 94 → 86 → 78 → 70 → 60%, extracellular ↓, intracellular ↑",
         "atoms 2/L_N ↓, forest ≈ 0.44, sea ↑; limits 0, ½, ½", "ANALOGY (X is not time; compartment pairing is a choice)"],
        ["growth aberration (SGA more water, LGA more fat)", "less / more congruence structure, less / more forest", "LITERAL measurement, ANALOGY wording"],
        ["osmotic equilibrium across membranes", "equidistribution of primes over reduced classes (V(q) is the imbalance)",
         "ANALOGY (no flow, no membrane)"],
        ["postnatal contraction of extracellular water, diuresis", "none", "NO COUNTERPART"],
        ["60-40-20 rule", "no fixed fractions at finite X; only the limit (0, ½, ½) is constant", "NO"],
    ])}

## What this does and does not show

It shows that the energy budget of S can be measured the way a physiologist measures body water, and that the
corrections physiology learned the hard way are needed here too. The windows need their exact capture. The sea is not
evenly mixed: it is thinner near the rationals, at the level of the variance of primes in progressions. Once both are
handled, the atoms and the forest match the no-fit spike model to within {err_big:.1f}% at X ≥ 3·10^4 ({err_small:.1f}% at
10^4). As X grows the
composition drifts toward the limit (0, ½, ½) like 1/ln X.

It does not show any mechanism shared with physiology, and it gives no new information about Goldbach's problem: the
forest's energy is the major-arc energy of results.md, and the sea is the minor-arc obstacle of escapes.md.

## References

- A. Tobias, B. D. Ballard, S. S. Mohiuddin, "Physiology, Water Balance", StatPearls (2022), NCBI Bookshelf NBK541059,
  {SOURCE_ADULT} (the full text could not be fetched from this environment; the figures used are the standard ones it
  states for a 70 kg man).
- "Body Composition in the Fetus and Newborn: Effects of Intrauterine Growth Aberration", ch. 2 of *Nephrology and
  Fluid/Electrolyte Physiology: Neonatology Questions and Controversies*, via {SOURCE_FETAL} (open part only).
- Variance of primes in progressions: Barban (1964); Davenport & Halberstam (1966); Montgomery (1970); Hooley (1975);
  Friedlander & Goldston (1996), as in fractals.md.

Run: {res['runtime_s'] / 60:.1f} min, peak RSS {res['peak_rss_mb']:,.0f} MB.
"""
    Path(path).write_text(text)
