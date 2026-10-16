"""Figure fig_W1_compartments.png and compartments.md for the body-water lens (every number read from the pickle)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

import figures as base  # noqa: F401  (shared rcParams / palette)
from escapes_report import _table
from figures import AQUA, BLUE, INK, INK2, MUTED, ORANGE, RED, VIOLET, YELLOW, _save, plt

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
DERIVED = {("first trimester", "solids"), ("24 weeks", "ICW"), ("24 weeks", "solids"), ("1 year", "ECW"), ("adult man", "solids")}
SOURCE_FETAL = "https://abdominalkey.com/body-composition-in-the-fetus-and-newborn-effects-of-intrauterine-growth-aberration/"
SOURCE_ADULT = "https://www.ncbi.nlm.nih.gov/books/NBK541059/"


def figure(res, out: Path) -> None:
    out = Path(out)
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    # (a) physiology
    ax = axes[0, 0]
    xs = np.arange(len(PHYSIOLOGY))
    for i, p in enumerate(PHYSIOLOGY):
        if p["ECW"] is None:
            ax.bar(i, p["TBW"], 0.6, color="white", edgecolor=MUTED, hatch="//", label="water (not split)")
            ax.bar(i, p["solids"], 0.6, bottom=p["TBW"], color=MUTED, edgecolor="white")
        else:
            ax.bar(i, p["ECW"], 0.6, color=YELLOW, edgecolor="white", label="extracellular water" if i == 1 else None)
            ax.bar(i, p["ICW"], 0.6, bottom=p["ECW"], color=VIOLET, edgecolor="white", label="intracellular water" if i == 1 else None)
            ax.bar(i, p["solids"], 0.6, bottom=p["TBW"], color=MUTED, edgecolor="white", label="solids" if i == 1 else None)
        ax.text(i, 101, f"water {p['TBW']}%", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(xs, [p["stage"] for p in PHYSIOLOGY], fontsize=8)
    ax.set_ylim(0, 128)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("% of body weight")
    ax.set_title("(a) physiology, for comparison: body water by age (colours not paired with (b)–(d))")
    ax.legend(loc="upper center", ncol=4, fontsize=7)
    # (b) development of S's compartments with X
    ax = axes[0, 1]
    dev = res["development"]
    X = np.array([r["X"] for r in dev], dtype=float)
    meas = {k: np.array([r["local_centre"][k] / r["L2"] for r in dev]) for k in ("A", "F", "sea")}
    pred = {"A": np.array([r["pred"]["atoms"] for r in dev]), "F": np.array([r["pred"]["forest"] for r in dev]),
            "sea": np.array([r["pred"]["sea"] for r in dev])}
    lc = res["limits_curve"]
    Xf = 10.0 ** np.asarray(lc["log10X"])
    for i, (key, col, lab) in enumerate((("A", RED, "atoms"), ("F", BLUE, "forest"), ("sea", AQUA, "sea"))):
        ax.plot(Xf, np.asarray(lc["shares"])[:, i], color=col, lw=0.9, ls=(0, (4, 2)))
        ax.plot(X, meas[key], "o", mfc="none", mec=col, ms=5, mew=1.2, label=f"{lab}, measured")
    ax.plot([], [], color=INK2, lw=0.9, ls=(0, (4, 2)), label="leading-order law (a model), to X = 10⁴⁰")
    ax.axhline(0.5, color=INK, lw=0.8, ls=(0, (1, 1.5)))
    ax.text(2e4, 0.51, "½", fontsize=8, color=INK2)
    ax.set_xscale("log")
    ax.set_ylim(0, 0.62)
    ax.set_xlabel("X (the 'growth' parameter; not time)")
    ax.set_ylabel("share of ΣΛ² (Parseval: the conserved total)")
    ax.set_title("(b) the energy of S by compartment, and the leading-order law toward (0, ½, ½)")
    ax.legend(loc="center right", fontsize=6.6)
    ins = ax.inset_axes([0.27, 0.2, 0.33, 0.28])
    for key, col in (("A", RED), ("F", BLUE), ("sea", AQUA)):
        ins.plot(X, pred[key], color=col, lw=1.0)
        ins.plot(X, meas[key], "o", mfc="none", mec=col, ms=4, mew=1.0)
    ins.set_xscale("log")
    ins.set_title("X = 10⁴ … 2·10⁶: measured (○) vs no-fit model (—)", fontsize=6.5, loc="left")
    ins.tick_params(labelsize=6)
    # (c) the sea is thinner near the forest
    ax = axes[1, 0]
    row = next(r for r in dev if r["X"] == 10**6)
    bands = [b for b in row["sea_bands"] if b["q"][0] >= 31]
    xb = np.arange(len(bands))
    wdt = 0.2
    ax.bar(xb - 1.5 * wdt, [b["eff"] / row["X"] for b in bands], wdt, color=INK2, label="measured in the windows (reading − capture × spike)")
    ax.bar(xb - 0.5 * wdt, [b["centre"] / row["X"] for b in bands], wdt, color=MUTED, label="measured at the spike centres")
    ax.bar(xb + 0.5 * wdt, [b["fg"] / row["X"] for b in bands], wdt, color=ORANGE,
           label="divisor sum of X·FG(d), d | q: variances of primes in progressions (no fit)")
    ax.axhline(row["rho"] / row["X"], color=INK, lw=1.0, ls=(0, (4, 2)), label="uniform sea (what plain dilution assumes)")
    ax.set_xticks(xb, [f"q ∈ [{b['q'][0]}, {b['q'][1]}]" for b in bands], fontsize=8)
    ax.set_ylim(0, 7.6)
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
    hd = [Patch(color=RED, label="atoms (exact copies of the main spike)"), Patch(color=BLUE, label="forest"),
          Patch(color=AQUA, label="sea"), Line2D([], [], marker="|", color=INK, ms=10, mew=1.6, ls="none", label="no-fit model")]
    ax.legend(handles=hd, loc="upper center", ncol=4, fontsize=6.8)
    ax.set_yticks(y, [s["name"] for s in ab], fontsize=8)
    ax.set_ylim(-0.6, len(ab) - 0.1)
    ax.set_xlim(0, 1.18)
    ax.set_xlabel("share of the set's own ΣΛ² (uniform-sea reading)")
    ax.set_title("(d) 'growth aberration': removing congruence structure removes forest (X = 10⁶)")
    fig.suptitle("Fig W1. The energy of S as body-water-style compartments: one conserved total, nested windows, recovery and mixing corrections",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_W1_compartments.png")


def write(res, path: Path) -> None:
    dev = res["development"]
    r6 = next(r for r in dev if r["X"] == 10**6)
    T6 = r6["L2"]
    lc, lf, p6 = r6["local_centre"], r6["local_fg"], r6["pred"]
    pct = lambda x: f"{100 * x:.2f}%"  # noqa: E731
    rel = lambda a, b: f"{100 * (a / b - 1):+.1f}%"  # noqa: E731

    cert = _table(["check", "result", "detail"], [[n, "pass" if ok else "FAIL", d] for n, ok, d in res["certificates"]])
    pq = _table(["q", "window energy / (capture × model)", "window capture, min–max"],
                [[r["q"], f"{r['ratio']:.4f}", f"{r['cap_min']:.4f}–{r['cap_max']:.4f}"] for r in r6["per_q"]])
    b1 = _table(["indicator", "what it reads", "share of ΣΛ²"], [
        ["whole circle (Parseval)", "the conserved total ΣΛ²", "100%"],
        ["windows at q = 1, 2", "raw reading", pct(r6["E_A"] / T6)],
        [f"windows at 3 ≤ q ≤ Q = {r6['Q']}", "raw reading", pct(r6["E_F"] / T6)],
    ])
    b1s = _table(["compartment", "one average capture, uniform sea", "exact capture, uniform sea",
                  "exact capture, sea from Σ_{d|q} of X·FG(d) (no fit)", "exact capture, sea measured at the centres", "no-fit model"], [
        ["atoms (q = 1, 2)", pct(r6["naive_A"] / T6), pct(r6["A"] / T6), pct(lf["A"] / T6), pct(lc["A"] / T6), pct(p6["atoms"])],
        ["forest (3 ≤ q ≤ Q)", pct(r6["naive_F"] / T6), pct(r6["F"] / T6), pct(lf["F"] / T6), pct(lc["F"] / T6), pct(p6["forest"])],
        ["sea (the rest)", pct(1 - (r6["naive_A"] + r6["naive_F"]) / T6), pct(r6["sea"] / T6), pct(lf["sea"] / T6),
         pct(lc["sea"] / T6), pct(p6["sea"])],
    ])
    bands = [b for b in r6["sea_bands"] if b["q"][0] >= 31]
    bt = _table(["windows of", "# q", "measured in the windows", "measured at the centres", "divisor sum of X·FG(d) (no fit)",
                 "uniform sea"],
                [[f"q ∈ [{b['q'][0]}, {b['q'][1]}]", b["n_q"], f"{b['eff'] / r6['X']:.2f} X", f"{b['centre'] / r6['X']:.2f} X",
                  f"{b['fg'] / r6['X']:.2f} X", f"{r6['rho'] / r6['X']:.2f} X"] for b in bands])
    fg_gap = [1 - b["fg"] / b["eff"] for b in bands]
    dt = _table(["X", "L_N", "q* (Q = ⌊q*⌋)", "atoms: measured / model", "forest: measured / model", "sea: measured / model",
                 "forest with a uniform sea"],
                [[f"{r['X']:.0e}".replace("e+0", "e"), f"{r['pred']['L_N']:.2f}", f"{r['q_star']:.1f} ({r['Q']})",
                  f"{r['local_centre']['A'] / r['L2']:.4f} / {r['pred']['atoms']:.4f}",
                  f"{r['local_centre']['F'] / r['L2']:.4f} / {r['pred']['forest']:.4f}",
                  f"{r['local_centre']['sea'] / r['L2']:.4f} / {r['pred']['sea']:.4f}", f"{r['F'] / r['L2']:.4f}"] for r in dev])
    errs = {r["X"]: max(abs(r["local_centre"]["A"] / r["L2"] / r["pred"]["atoms"] - 1),
                        abs(r["local_centre"]["F"] / r["L2"] / r["pred"]["forest"] - 1)) for r in dev}
    err_big = 100 * max(v for x, v in errs.items() if x >= 3 * 10**4)
    err_small = 100 * errs[dev[0]["X"]]
    asy = _table(["X", "atoms", "forest", "sea"],
                 [[f"10^{d['log10X']}", f"{d['atoms']:.3f}", f"{d['forest']:.3f}", f"{d['sea']:.3f}"] for d in res["limits"]])
    cell = lambda p, k: ("—" if p[k] is None else f"{'≈ ' if (p['stage'], k) in DERIVED else ''}{p[k]}%")  # noqa: E731
    phys = _table(["stage", "total water", "extracellular", "intracellular", "solids"],
                  [[p["stage"], cell(p, "TBW"), cell(p, "ECW"), cell(p, "ICW"), cell(p, "solids")] for p in PHYSIOLOGY])
    ab = res["aberration"]
    abr = []
    for s in ab:
        T = s["L2"]
        if s["n_qF"] == 0:
            F_txt = f"no forest windows (model has none); control in the primes' windows: {s['F_control'] / T:+.4f}"
        else:
            F_txt = f"{s['F'] / T:.4f} / {s['pred_F'] / T:.4f}"
        abr.append([s["name"], s["desc"], ", ".join(map(str, s["atom_qs"])), f"{s['A'] / T:.4f} / {s['pred_A'] / T:.4f}", F_txt,
                    f"{s['sea'] / T:.4f}", s["n_qF"]])
    abt = _table(["set (X = 10^6, same total weight ≈ ψ)", "what it keeps", "atom denominators (full-height spikes)",
                  "atoms: measured / model", "forest: measured / model", "sea", "# forest denominators with model spikes"], abr)
    m4 = ab[-1]
    pr_ = next(s for s in ab if s["name"] == "primes")
    a0, a1 = dev[0], dev[-1]
    coh0 = (a0["local_centre"]["A"] + a0["local_centre"]["F"]) / a0["L2"]
    coh1 = (a1["local_centre"]["A"] + a1["local_centre"]["F"]) / a1["L2"]
    pp = max(r["prime_power_max_dev"] for r in dev)
    text = f"""# Compartments: the body-water lens on the energy of S at N = 10^6

> Exploration at finite X, not a proof. Every number below is read from `cache/compartments.pkl` (grids M = 2^17 … 2^24,
> X = 10^4 … 2·10^6); the physiology numbers come from the two sources listed at the end. Body-water physiology is used as
> a teacher for bookkeeping, not as a mechanism: nothing in S moves, mixes or equilibrates.

**What physiology does.** Total body water (TBW) is about 60% of an adult man's weight (42 L of 70 kg). Two thirds of it
is intracellular (28 L) and one third extracellular (14 L: 10.5 L interstitial and 3.5 L plasma) (StatPearls, *Physiology,
Water Balance*). Compartments are measured by **indicator dilution**: give a known amount I of an indicator that spreads
only through the compartment of interest, let it equilibrate, read its plasma concentration C, and set V = I/C. Heavy
water (D₂O, H₂¹⁸O) or antipyrine measures TBW; bromide, sucrose or inulin measures extracellular water; Evans blue
measures plasma. The rest follows by subtraction: intracellular = TBW − extracellular, interstitial = extracellular −
plasma, solids = body weight − TBW (Abdominal Key, ch. 2). The composition changes with growth. TBW falls from 94% of body
weight in the first trimester to about 86% at 24 weeks, 78% at term and about 70% at one year, as solids accrete.
Extracellular water falls (60% → 44% → ≈ 28%) and intracellular water rises (≈ 26% → 34% → 42%). It also depends on growth
conditions: infants small for gestational age have more water and less fat; infants large for gestational age have more
fat and less water.

**What the circle has.** One conserved total, by Parseval: (1/M)Σ_k |S(k/M)|² = ΣΛ² exactly. The windows
|N(α − a/q)| ≤ ½ around the reduced fractions a/q of a set of denominators play the indicators. A window reads the energy it
covers: a known share (its capture) of the Dirichlet spike centred in it, plus the sea it covers. Three compartments:

- **atoms**: the spikes at q = 1, 2 (α = 0 and ½);
- **forest**: the spikes with 3 ≤ q ≤ Q = ⌊q*(X)⌋, where q* is the denominator at which spike and sea are equal
  (fractals.md, M1);
- **sea**: the rest, i.e. the structureless part plus the spikes with q > Q that it drowns.

The no-fit spike model is S(a/q + β) ≈ μ(q)(ψ/φ(q))·v(β)/X, with v the Dirichlet kernel of energy X. So the φ(q) spikes of
denominator q carry ψ²μ²(q)/(Xφ(q)) together, and the shares are atoms = 2/L_N and forest = (Σ_{{3≤q≤Q}} μ²/φ)/L_N, with
L_N = XΣΛ²/ψ². L_N is the total energy divided by the energy of one reference spike (the one at 0, ψ²/X), a dilution
volume in the literal sense of a ratio. At X = 10^6, L_N = {p6['L_N']:.2f}: the energy of S is spread as if over
{p6['L_N']:.1f} copies of its main spike.

## 0. Certificates

{cert}

The divisor identity is exact once the variances are taken over n coprime to q. With the prime powers p^k (p | q) that
S(a/q) actually contains, the centre level moves by up to {100 * pp:.0f}% for a few q at small X; the local-sea reading
uses the actual centre level.

Each window's energy against its capture times the model spike, at X = 10^6 (small q, where the sea is negligible):

{pq}

## B1. Compartments at X = 10^6

{b1}

The readings become compartments by a three-window solution. Each window reads its capture of its own spikes plus the
sea it covers; the rest of the circle holds the spikes' leakage and the remaining sea:

{b1s}

Two corrections matter, and each resembles one that tracer practice needs:

1. **Incomplete recovery.** A window catches only part of its spike, and how much depends on where the fraction falls
   between grid points. At X = 10^6 a window holds 8 or 9 grid points and catches 0.755 or 0.803 of its spike, not the
   average 0.7737. Using one average capture misreads the atoms by {rel(r6['naive_A'], r6['A'])}. The error changes sign
   with the grid offset, so it is a recovery correction rather than a tracer leaking one way.
2. **The tracer outside its space, and uneven mixing.** A window also reads sea that is not its spike's: that is the
   counterpart of an indicator that reaches beyond its compartment (in standard tracer practice, not from the two
   sources: heavy water exchanges with non-aqueous hydrogen, and bromide enters cells a little). Plain dilution assumes
   that sea is evenly mixed. It is not (next table). With a uniform sea the forest reads {rel(r6['F'], T6 * p6['forest'])}
   against the model. With the sea level measured at the spike centres it reads {rel(lc['F'], T6 * p6['forest'])}; with the
   no-fit divisor form below it reads {rel(lf['F'], T6 * p6['forest'])}.

**The sea is thinner near the forest.** Sea level per unit α in the windows of each band of denominators, at X = 10^6:

{bt}

Near a reduced fraction a/q (q squarefree) the structureless part of S has mean square

    (1/φ(q)) Σ_{{d|q}} μ(q/d) · d · V(X; d),

where V(X; d) is the variance of the primes over the reduced residue classes mod d (V(X; 1) = V(X; 2) = 0). This follows
from Σ_{{a reduced}} e(a(r − s)/q) = c_q(r − s) = Σ_{{d|q, d|r−s}} d μ(q/d); the certificate above checks it to rounding. For prime
q it is V(X; q)·q/(q − 1), about X·FG(q) with FG(q) = ln q − γ − ln 2π − ln q/(q − 1) (Friedlander–Goldston form of the
Barban–Davenport–Halberstam / Montgomery–Hooley variance). For composite q the divisor sum is larger. Replacing each
V(X; d) by X·FG(d) gives the no-fit column, which falls {100 * min(fg_gap):.0f}–{100 * max(fg_gap):.0f}% below the window values.
Every band sits well below the circle's average sea level, and the level climbs back toward it as q grows. (Bands with
q ≤ 30 are left out: there the sea is too small a fraction of the window to measure.)

## B2. Development: the composition of S as X grows

{dt}

With the centre-measured sea, the atoms and the forest match the no-fit model to within {err_big:.1f}% at X ≥ 3·10^4
({err_small:.1f}% at 10^4).

For comparison, the body (percent of body weight; fetal and infant values from the Abdominal Key chapter, adult from
StatPearls; "≈" marks a value obtained by subtraction):

{phys}

Both show one conserved total redistributing among compartments as a growth parameter increases. The **atoms** fall like
2/L_N ≈ 2/(ln X − 1), from {a0['local_centre']['A'] / a0['L2']:.1%} at X = 10^4 to {a1['local_centre']['A'] / a1['L2']:.1%} at
2·10^6: their energy 2ψ²/X grows only like X, while the total grows like X ln X. The **forest** stays near 0.44 and the **sea**
grows.

The leading-order law below is a model, not a theorem: it uses the Friedlander–Goldston asymptotic, the spike model up to
q* and q* from q = √(X/FG(q)). In that model the shares tend to **0, ½, ½**, because q* grows like √X up to logarithms and the
forest holds (ln q* + 1.33 − 2)/(ln X − 1). The atoms approach 0 like 2/ln X, but the forest and the sea approach ½ only like
ln ln X/(2 ln X), and the sea passes ½ near X ≈ 10^38 and returns to it from above:

{asy}

**Where the analogy stops.** The coherent share of S (atoms plus forest) *falls* with X, from {coh0:.1%} to {coh1:.1%}
(model limit ½), while a fetus *gains* solids. Which compartment of S "corresponds" to extracellular or intracellular water
is a choice with no physical basis. Only the bookkeeping is literal; the limit belongs to the leading-order model. X is a
parameter, not time, and there is no counterpart to birth or to the postnatal contraction of extracellular water.

## B3. "Growth aberration": the arithmetic of the set sets its composition (X = 10^6)

Four sets with the same total weight (≈ ψ) but different congruence structure. Each set's no-fit model spreads its weight
evenly over the residue classes it allows mod lcm(D, q) and sums them with e(ab/q). **Atoms** are the set's exact copies of
its main spike (model height equal to the height at q = 1); the forest windows sit where the set's model has spikes, up to
the primes' Q = {res['aberration_Q']}.

{abt}

Removing congruence structure removes forest. Cramér's model has one atom (q = 1) and no forest. Restricting it to odd
numbers adds the atom at ½, and nothing more. The primes have the whole Farey forest. Read in the primes' forest windows,
the two Cramér sets give forests of {ab[0]['F_control'] / ab[0]['L2']:+.4f} and {ab[1]['F_control'] / ab[1]['L2']:+.4f}: that is a
measurement, not just the model.

Adding a congruence condition does **not** add forest. The primes ≡ 1 (mod 4) live on one class mod 4, so
S(α + ¼) = i·S(α): the spikes at ¼ and ¾ are exact copies of the atoms at 0 and ½ and count as atoms. The set then has the
same atoms as the primes ({m4['A'] / m4['L2']:.1%}) but a smaller forest ({m4['F'] / m4['L2']:.1%} against
{pr_['F'] / pr_['L2']:.1%}; model {m4['pred_F'] / m4['L2']:.1%} against {pr_['pred_F'] / pr_['L2']:.1%}) and a larger sea. Filing
the q = 4 copies under the forest instead would give {m4['q12']['F'] / m4['L2']:.1%}, which is the artefact a first draft of
this file reported. The forest readings here use the uniform sea, so they share the primes'
{rel(r6['F'], T6 * p6['forest'])} bias; the comparison between sets is unaffected.

In body-composition terms (small for gestational age: more water, less fat) only the first direction carries over: less
structure, less forest, more sea. The "large for gestational age" direction has no counterpart here. The measurement is
literal; the nutrition language is not.

## Literal versus analogy

{_table(["body water", "on the circle", "status"], [
        ["total body water, conserved and split into compartments", "Parseval: ΣΛ² split into atoms, forest and sea", "LITERAL bookkeeping"],
        ["dilution volume V = I/C", "L_N = total energy / energy of one spike", "LITERAL ratio"],
        ["V = I/C as the measuring rule", "compartment = (window reading − sea covered)/capture, a recovery correction",
         "ANALOGY (same arithmetic form; nothing is given or diluted)"],
        ["indicators with different spaces (D₂O, bromide/inulin, Evans blue)", "windows around different sets of denominators; nested subtraction",
         "LITERAL nesting"],
        ["incomplete recovery of an indicator", "a window catches 0.755–0.803 of its spike; one average capture errs by "
         f"{rel(r6['naive_A'], r6['A'])} on the atoms, with a sign that depends on the grid offset", "ANALOGY"],
        ["an indicator reaching beyond its space", "the sea a window covers", "ANALOGY"],
        ["equilibration (the indicator mixes evenly)", "uniform sea: fails near rationals, where the sea is the divisor sum of the "
         "residue-class variances", "MEASURED failure (with an exact identity)"],
        ["development: TBW 94 → 86 → 78 → 70 → 60%, extracellular ↓, intracellular ↑ to one year",
         "atoms 2/L_N ↓, forest ≈ 0.44, sea ↑; leading-order limit 0, ½, ½", "ANALOGY (X is not time; compartment pairing is a choice)"],
        ["growth aberration: small for gestational age (more water, less fat)", "less congruence structure: less forest, more sea",
         "LITERAL measurement, ANALOGY wording"],
        ["growth aberration: large for gestational age (more fat)", "an extra congruence condition adds atoms (copies), not forest",
         "NO COUNTERPART"],
        ["osmotic equilibrium across membranes", "equidistribution of primes over reduced classes (V(q) is the imbalance)",
         "ANALOGY (no flow, no membrane)"],
        ["postnatal contraction of extracellular water, diuresis", "none", "NO COUNTERPART"],
        ["60-40-20 rule", "no fixed fractions at finite X; only the model limit (0, ½, ½) is constant", "NO"],
    ])}

## What this does and does not show

It shows that the energy budget of S can be measured with body-water bookkeeping, and that two corrections physiology
knows are needed here too. Each window needs its exact capture. The sea is not evenly mixed: near a/q it sits at the
divisor sum of the variances of primes in progressions, well below its average. With both handled, the atoms and the
forest match the no-fit spike model to within {err_big:.1f}% at X ≥ 3·10^4. Removing congruence structure from the set
removes forest; adding a condition adds copies of the main spike instead.

It does not show any mechanism shared with physiology, and it gives no new information about Goldbach's problem. The
forest is the same kind of quantity as the major-arc mass of results.md (larger Q, narrower arcs, without q = 1, 2), and
the sea is the minor-arc obstacle of escapes.md.

## References

- A. Tobias, B. D. Ballard, S. S. Mohiuddin, "Physiology, Water Balance", StatPearls (2022), NCBI Bookshelf NBK541059,
  {SOURCE_ADULT} (the full text could not be fetched from this environment; the figures used are the standard ones it
  states for a 70 kg man).
- "Body Composition in the Fetus and Newborn: Effects of Intrauterine Growth Aberration", ch. 2 of *Nephrology and
  Fluid/Electrolyte Physiology: Neonatology Questions and Controversies*, via {SOURCE_FETAL} (open part only).
- Variance of primes in progressions: Barban (1964); Davenport & Halberstam (1966); Montgomery (1970); Hooley (1975);
  Friedlander & Goldston (1996), as in fractals.md.

Run: {res['runtime_s']:.0f} s, peak RSS {res['peak_rss_mb']:,.0f} MB.
"""
    Path(path).write_text(text)
