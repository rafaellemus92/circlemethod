"""All figures.  Every function takes the context dict produced by run_all.compute()."""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# ---------------------------------------------------------------- style (reference palette, light)
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
MAGENTA, GREEN, VIOLET, RED = "#e87ba4", "#008300", "#6250d6", "#e34948"
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seqblue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "text.color": INK, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False, "legend.fontsize": 8,
    "lines.linewidth": 1.4, "lines.solid_capstyle": "round", "savefig.dpi": 160,
})


def _save(fig, out: Path, name: str) -> None:
    fig.savefig(out / name, bbox_inches="tight")
    plt.close(fig)


def _sci(x: float, digits: int = 2) -> str:
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10**e
    return f"{m:.{digits}f}×10$^{{{e}}}$"


def farey_list(qmax: int) -> list[tuple[int, int]]:
    return sorted({(a, q) for q in range(1, qmax + 1) for a in range(q) if math.gcd(a, q) == 1},
                  key=lambda t: t[0] / t[1])


# ================================================================ Fig A
def fig_A(ctx, out):
    psi, L2 = ctx["psi"], ctx["L2"]
    env, med = ctx["env_max"], ctx["env_med"]
    x = (np.arange(len(env)) + 0.5) / len(env)
    exact = {(r["a"], r["q"]): r for r in ctx["farey"]}
    fig, ax = plt.subplots(figsize=(16, 6.2))
    ax.semilogy(x, env, color=BLUE, lw=0.7, label=f"max |S| in each bin of width 1/{len(env)} ($\\approx$61/N)")
    ax.semilogy(x, med, color=MUTED, lw=0.6, label="median |S| in each bin (the incoherent sea)")
    ax.axhline(math.sqrt(L2), color=INK2, lw=0.8,
               label=f"random-walk scale $\\sqrt{{\\Sigma\\Lambda^2}}$ = {math.sqrt(L2):.0f}")
    first_pred = first_zero = True
    for a, q in farey_list(12):
        r = exact[(a, q)]
        if r["mu"] != 0:
            h = r["pred"]
            ax.plot([a / q], [h], marker="_", ms=11, mew=2.2, color=ORANGE,
                    label="predicted |μ(q)|ψ/φ(q)" if first_pred else None)
            first_pred = False
            ax.text(a / q, h * 1.35, f"{a}/{q}", rotation=90, ha="center", va="bottom", fontsize=6.5, color=INK)
        else:
            h = r["absS"]
            ax.plot([a / q], [h], marker="x", ms=6, mew=1.4, color=RED,
                    label="non-squarefree q: prediction 0, exact |S(a/q)| shown" if first_zero else None)
            first_zero = False
            ax.text(a / q, h / 1.5, f"{a}/{q}", rotation=90, ha="center", va="top", fontsize=6.5, color=INK2)
    ax.set_ylim(80, 4e6)
    ax.set_xlim(-0.005, 1.005)
    ax.set_xlabel("α")
    ax.set_ylabel("|S(α)|  (log scale)")
    ax.set_title(f"Fig A. The Farey forest: |S(α)| for N = 10$^6$ on the M = 2$^{{23}}$ grid "
                 f"(ψ(N) = {psi:,.0f}); every a/q with q ≤ 12 marked")
    ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.09))
    _save(fig, out, "fig_A_farey_forest.png")


# ================================================================ Fig B
def fig_B(ctx, out):
    psi, L2 = ctx["psi"], ctx["L2"]
    rows = ctx["farey"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.2))
    for q in range(1, 31):
        rq = [r for r in rows if r["q"] == q]
        xs = q + np.linspace(-0.22, 0.22, len(rq)) if len(rq) > 1 else np.array([q])
        sf = rq[0]["mu"] != 0
        ys = np.array([r["absS"] for r in rq]) / psi
        if sf:
            a1.plot([q - 0.46, q + 0.46], [1 / rq[0]["phi"]] * 2, color=ORANGE, lw=1.6, solid_capstyle="butt",
                    zorder=4)
            a1.scatter(xs, ys, s=8, color=BLUE, zorder=3, lw=0)
        else:
            a1.scatter(xs, ys, s=22, color=RED, marker="x", zorder=3, linewidth=1.2)
        re = np.array([r["S"].real for r in rq]) * rq[0]["phi"] / psi
        im = np.array([r["S"].imag for r in rq]) * rq[0]["phi"] / psi
        a2.plot([q - 0.46, q + 0.46], [rq[0]["mu"]] * 2, color=ORANGE, lw=1.6, solid_capstyle="butt", zorder=4)
        a2.scatter(xs, re, s=8 if sf else 20, color=BLUE if sf else RED, marker="o" if sf else "x", zorder=3,
                   linewidth=0 if sf else 1.2)
        a2.scatter(xs, im, s=5, color=MUTED, zorder=2)
    a1.axhline(math.sqrt(L2) / psi, color=INK2, lw=0.8)
    a1.text(30.6, math.sqrt(L2) / psi * 1.12, "$\\sqrt{\\Sigma\\Lambda^2}/ψ$", ha="right", fontsize=7.5, color=INK2)
    a1.set_yscale("log")
    a1.set_xlabel("q")
    a1.set_ylabel("|S(a/q)| / ψ(N)")
    a1.set_title("(a) exact |S(a/q)| for every reduced a/q, q ≤ 30")
    a1.plot([], [], color=ORANGE, lw=1.6, label="prediction |μ(q)|/φ(q)")
    a1.scatter([], [], s=12, color=BLUE, label="measured, squarefree q")
    a1.scatter([], [], s=22, color=RED, marker="x", label="measured, non-squarefree q (prediction 0)")
    a1.legend(loc="upper right")
    m4 = ctx["mod4"]
    a1.annotate(f"S(1/4) = {m4['S_quarter'].real:.1f} − {abs(m4['S_quarter'].imag):.1f} i\n"
                f"Im S(1/4) = ψ(N;4,1) − ψ(N;4,3)  (the mod-4 prime race)", xy=(4, abs(m4["S_quarter"]) / psi),
                xytext=(5.0, 1.0e-2), fontsize=7.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    a1.set_xlim(0.3, 30.7)
    a2.set_xlim(0.3, 30.7)
    a2.set_ylim(-1.35, 1.35)
    a2.set_xlabel("q")
    a2.set_ylabel("φ(q) · S(a/q) / ψ(N)")
    a2.set_title("(b) the resultant is real with sign μ(q): φ(q)·Re S(a/q)/ψ ≈ μ(q)")
    a2.plot([], [], color=ORANGE, lw=1.6, label="μ(q)")
    a2.scatter([], [], s=12, color=BLUE, label="φ(q)·Re S(a/q)/ψ")
    a2.scatter([], [], s=5, color=MUTED, label="φ(q)·Im S(a/q)/ψ")
    a2.legend(loc="lower left", ncol=3)
    fig.suptitle("Fig B. Spike heights at rationals: measured by exact residue-class sums vs the "
                 "Ramanujan-sum prediction", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_B_spike_heights.png")


# ================================================================ Fig C / D
def _zoom_panel(axm, axp, beta, S, N, model=None, title="", ylim=(30, 1.2e6), marks=(), label_y=None):
    x = N * beta
    carrier = np.exp(-1j * np.pi * (N + 1) * beta)
    axm.semilogy(x, np.abs(S), color=BLUE, lw=1.0, label="measured |S|")
    if model is not None and np.any(model != 0):
        axm.semilogy(x, np.abs(model), color=ORANGE, lw=1.0, ls=(0, (4, 2)), label="model |μ(q)/φ(q)·v(β)|")
    axm.set_ylim(*ylim)
    axm.set_title(title)
    for xm, h, lab in marks:
        axm.plot([xm], [h], marker="_", ms=10, mew=2, color=ORANGE)
        if lab:
            ly = label_y if label_y else h * 1.3
            if label_y:
                axm.plot([xm, xm], [h * 1.15, ly / 1.1], color=GRID, lw=0.6, zorder=0)
            axm.text(xm, ly, lab, ha="center", va="bottom", fontsize=6.5, rotation=90, color=INK2)
    # wrap into [-pi/2, 3pi/2) so the model's two levels 0 and pi sit mid-axis instead of on the cut
    wrap = lambda p: (p + np.pi / 2) % (2 * np.pi) - np.pi / 2  # noqa: E731
    ph = wrap(np.angle(S * carrier))
    axp.scatter(x, ph, s=1.2, color=BLUE, lw=0, label="measured", zorder=3)
    if model is not None and np.any(model != 0):
        amp = (model * carrier).real                      # the model is real in the co-rotating frame
        mp = np.where(amp >= 0, 0.0, np.pi)
        mp[1:][np.diff(mp) != 0] = np.nan                 # no vertical connectors at sign flips
        axp.plot(x, mp, color=ORANGE, lw=1.6, label="model (0 or π)", zorder=2)
    axp.set_ylim(-np.pi / 2 - 0.2, 3 * np.pi / 2 + 0.2)
    axp.set_yticks([-np.pi / 2, 0, np.pi / 2, np.pi, 3 * np.pi / 2], ["−π/2", "0", "π/2", "π", "3π/2"])
    axp.set_xlabel("N·β   (β = α − centre)")


def fig_C(ctx, out):
    N, psi = ctx["N"], ctx["psi"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 6.4), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    for j, z in enumerate(ctx["zooms"]):
        a, q = z["a"], z["q"]
        mu = {3: -1, 6: 1, 4: 0}[q]
        phi = {3: 2, 6: 2, 4: 2}[q]
        pred = f"μ/φ = {mu:+d}/{phi}" if mu else "μ(4) = 0: model ≡ 0"
        _zoom_panel(axes[0, j], axes[1, j], z["beta"], z["S"], N, z["model"],
                    title=f"α = {a}/{q} + β    ({pred})")
        axes[0, j].axhline(math.sqrt(ctx["L2"]), color=INK2, lw=0.7)
    axes[0, 2].text(19.5, math.sqrt(ctx["L2"]) * 1.15, "sea level $\\sqrt{\\Sigma\\Lambda^2}$", ha="right",
                    va="bottom", fontsize=7.5, color=INK2)
    axes[0, 0].set_ylabel("|S(α)|")
    axes[1, 0].set_ylabel("arg S·e(−(N+1)β/2)")
    axes[0, 0].legend(loc="upper left", fontsize=7.5)
    axes[1, 0].legend(loc="upper left", fontsize=7, markerscale=4, ncol=2)
    fig.suptitle("Fig C. Zooms of width 40/N around 1/3, 1/6, 1/4: measured S (direct summation) vs the "
                 "major-arc model (μ(q)/φ(q))·v(β); phase shown in the frame co-rotating with the carrier "
                 "e((N+1)β/2)", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_C_zooms.png")


def fig_D(ctx, out):
    N = ctx["N"]
    z = ctx["zoom_golden"]
    fig, axes = plt.subplots(2, 1, figsize=(6.2, 6.4), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    fr = sorted([f for f in z["fractions"] if f["height"] > 0], key=lambda f: -f["height"])
    top = {(f["a"], f["q"]) for f in fr[:6]}
    marks = [(f["N_beta"], f["height"], f"{f['a']}/{f['q']}" if (f["a"], f["q"]) in top else "")
             for f in fr if f["height"] >= 1000]
    _zoom_panel(axes[0], axes[1], z["beta"], z["S"], N, None,
                title="α = frac(golden ratio) + β", marks=marks, label_y=1.4e4)
    axes[0].axhline(math.sqrt(ctx["L2"]), color=INK2, lw=0.7)
    axes[0].text(19.5, math.sqrt(ctx["L2"]) * 1.12, "$\\sqrt{\\Sigma\\Lambda^2}$", ha="right", fontsize=7.5,
                 color=INK2)
    axes[0].set_ylabel("|S(α)|")
    axes[1].set_ylabel("arg S·e(−(N+1)β/2)")
    axes[0].plot([], [], marker="_", ms=10, mew=2, color=ORANGE, ls="none",
                 label="model crest ψ|μ(q)|/φ(q) ≥ 10$^3$ of each a/q (q ≤ 2000)\nin the window; tallest six labelled")
    axes[0].legend(loc="upper left", fontsize=7)
    fig.suptitle("Fig D. Same width and scale, the deepest minor arc:\nweak μ(q)-locked crests at sea level",
                 x=0.01, ha="left", fontsize=10.5,
                 fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_D_golden_zoom.png")


# ================================================================ phasors
def fig_P(ctx, out):
    psi, L2 = ctx["psi"], ctx["L2"]
    fig = plt.figure(figsize=(16, 8.6))
    th = np.deg2rad(np.arange(360))
    mu = {3: -1, 4: 0, 5: -1, 6: 1, 30: -1}
    phi = {3: 2, 4: 2, 5: 4, 6: 2, 30: 8}
    blurb = {3: "2 channels at ±120°: resultant ψ/2, backwards",
             4: "2 channels at 90°, 270°: head-on cancellation",
             5: "4 channels: resultant = −1 channel",
             6: "2 channels at ±60°: resultant ψ/2, forwards",
             30: "8 channels: resultant = −1 channel"}
    for i, p in enumerate(ctx["phasors"]):
        ax = fig.add_subplot(2, 4, i + 1, projection="polar")
        h = p["hist"]
        S = p["S"]
        ax.bar(th, h, width=np.deg2rad(1.0), color=BLUE, edgecolor="none", alpha=0.9)
        rmax = 1.18 * max(h.max(), abs(S))
        ax.set_ylim(0, rmax)
        ax.annotate("", xy=(np.angle(S), abs(S)), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.2, mutation_scale=14))
        ax.set_yticklabels([])
        ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
        ax.tick_params(labelsize=7, colors=MUTED)
        ax.grid(color=GRID, lw=0.6)
        ax.spines["polar"].set_color(AXIS)
        if p["q"] is not None:
            q = p["q"]
            ax.set_title(f"α = {p['name']}", loc="center", pad=12)
            pred = mu[q] * psi / phi[q]
            txt = (f"{blurb[q]}\nS = {S.real:+,.0f} {S.imag:+,.0f}i\n"
                   f"predicted μ(q)ψ/φ(q) = {pred:+,.0f}" if mu[q] else
                   f"{blurb[q]}\nS = {S.real:+,.0f} {S.imag:+,.0f}i  (|S| = {abs(S) / psi:.1e}ψ)\n"
                   f"predicted 0; residue = mod-4 prime race")
        else:
            ax.set_title(f"α = {p['name']}", loc="center", pad=12)
            txt = (f"phases spread uniformly: a flat ring\n|S| = {abs(S):,.0f} = {abs(S) / psi:.4f}ψ\n"
                   f"random-walk scale $\\sqrt{{\\Sigma\\Lambda^2}}$ = {math.sqrt(L2):,.0f}")
        ax.text(0.5, -0.2, txt, transform=ax.transAxes, ha="center", va="top", fontsize=7.5, color=INK)
    ax = fig.add_subplot(2, 4, 8)
    ax.axis("off")
    ax.text(0.0, 0.95,
            "How to read a panel\n\n"
            "Each n ≤ N with Λ(n) ≠ 0 is an oscillator of\namplitude Λ(n) and phase 2π·n·α.\n\n"
            "Blue bars: total amplitude Λ landing in each 1°\nphase bin (360 bins, radial scale per panel).\n\n"
            "Orange arrow: the resultant S(α) — the vector\nsum of all bars, drawn on the same radial scale.\n\n"
            "At α = a/q the phase depends only on n mod q,\nso the oscillators lock into the φ(q) reduced\n"
            "channels with ≈ ψ/φ(q) each; their resultant is\nψ·c_q(a)/φ(q) = μ(q)ψ/φ(q) (Ramanujan sum).",
            va="top", fontsize=8.2, color=INK, linespacing=1.35)
    fig.suptitle("Fig P. Phasor pictures: Λ-weighted histograms of the phases 2π·n·α (mod 2π) and their resultant",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.subplots_adjust(hspace=0.75, wspace=0.35, top=0.88, bottom=0.12)
    _save(fig, out, "fig_P_phasors.png")


# ================================================================ singular series
def fig_S(ctx, out):
    part = ctx["SS2_part"]
    SS2 = ctx["SS2"]
    Q = np.arange(len(part))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6))
    a1.step(Q[1:10001], part[1:10001], where="post", color=BLUE, lw=1.0)
    a1.axhline(SS2, color=ORANGE, lw=1.2)
    a1.text(9000, SS2 + 0.006, f"SS2(N) = 2C$_2$·4/3 = {SS2:.6f}", ha="right", fontsize=8, color=INK)
    for q in (5, 10, 20, 50, 100, 150):
        a1.plot([q], [part[q]], "o", ms=5, color=BLUE, mec=SURF, mew=1.2)
        a1.annotate(f"Q={q}", (q, part[q]), xytext=(4, -12) if q == 150 else (4, 6), textcoords="offset points",
                    fontsize=7, color=INK2)
    a1.set_xscale("log")
    a1.set_ylim(0.95, 2.05)
    a1.set_xlabel("Q")
    a1.set_ylabel("SS2$_Q$(10$^6$)")
    a1.set_title("(a) truncated singular series Σ$_{q≤Q}$ μ(q)$^2$ c$_q$(N)/φ(q)$^2$")
    d = part[1:] - SS2
    a2.step(Q[1:], d, where="post", color=BLUE, lw=0.6)
    a2.axhline(0, color=AXIS, lw=0.8)
    a2.set_xscale("log")
    a2.set_yscale("symlog", linthresh=1e-5)
    a2.set_ylim(-1.5, 1.5)
    a2.set_xlabel("Q")
    a2.set_ylabel("SS2$_Q$ − SS2")
    a2.set_title("(b) signed error SS2$_Q$ − SS2 (symlog): each squarefree q adds μ(q)$^2$c$_q$(N)/φ(q)$^2$")
    fig.suptitle("Fig S. The binary singular series for N = 10$^6$ = 2$^6$5$^6$ converges; the sign of each kick "
                 "is set by how many prime factors of q do not divide N", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_S_singular_series.png")


# ================================================================ Fig E
def fig_E(ctx, out):
    N, psi = ctx["N"], ctx["psi"]
    W = ctx["window"]
    x = N * (W["alpha"] - W["x0"])
    fig, axes = plt.subplots(3, 1, figsize=(13, 8.6), sharex=True)
    mass, integ = W["mass"] / 1e7, W["integrand"].real / 1e7
    lim = 1.65 * mass.max()          # top 35% of panel (a) is a label strip
    axes[0].plot(x, mass, color=BLUE, lw=0.8)
    axes[0].axhline(np.mean(mass), color=INK2, lw=0.8)
    axes[0].text(x[0], lim * 0.97, f"horizontal line: window mean |S|$^2$ = {np.mean(W['mass']):.2e}\n"
                 f"(circle mean ΣΛ$^2$ = {ctx['L2']:.2e})", ha="left", va="top", fontsize=7.5, color=INK2)
    axes[0].set_ylim(0, lim)
    axes[0].set_ylabel("|S|$^2$  (×10$^7$)")
    axes[0].set_title("(a) mass density |S(α)|$^2$ — always positive")
    axes[1].plot(x, integ, color=BLUE, lw=0.6)
    axes[1].set_ylim(-mass.max() * 1.08, mass.max() * 1.08)
    axes[1].axhline(0, color=AXIS, lw=0.8)
    axes[1].set_ylabel("Re S$^2$e(−Nα)  (×10$^7$)")
    axes[1].set_title("(b) integrand Re[S(α)$^2$ e(−Nα)] — the same envelope, but the phase spins")
    axes[2].plot(x, W["cum_mass"], color=ORANGE, lw=1.4, label="cumulative mass ∫|S|$^2$ dα")
    axes[2].plot(x, W["cum"], color=BLUE, lw=1.4, label="cumulative integral ∫Re S$^2$e(−Nα) dα")
    axes[2].axhline(0, color=AXIS, lw=0.8)
    axes[2].set_ylabel("running integral")
    axes[2].set_xlabel(f"N·(α − α$_0$),   α$_0$ = {W['x0']:.6f}   (window of width 300/N inside the minor "
                       "arcs of every (Q, c) tabulated at N = 10$^6$)")
    axes[2].set_title("(c) mass piles up linearly; the signed integral wanders like a random walk")
    axes[2].legend(loc="upper left")
    for i, f in enumerate(sorted(W["fractions"], key=lambda f: f["alpha"])):
        xf = N * (f["alpha"] - W["x0"])
        axes[0].vlines(xf, 0, mass.max() * 1.05, color=MUTED, lw=0.8, zorder=0)
        for ax in axes[1:]:
            ax.axvline(xf, color=MUTED, lw=0.8, zorder=0)
        crest = f"model crest$^2$ = {f['height'] ** 2 / 1e7:.1f}×10$^7$" if f["height"] > 0 else "μ(q) = 0: no crest"
        axes[0].text(xf, lim * (0.97 - 0.11 * (i % 3)), f"{f['a']}/{f['q']}: {crest}", ha="center",
                     va="top", fontsize=7, color=INK2)
    fig.suptitle("Fig E. One stretch of minor arc: mass vs integrand vs running integral", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_E_minor_stretch.png")


# ================================================================ Fig F
def fig_F(ctx, out):
    N = ctx["N"]
    st = ctx["stair"]
    al, F = np.r_[st["alpha"], 1.0], np.r_[st["F"], st["final"]]   # F(1) = full circle = R2
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(14, 8.4), sharex=True, gridspec_kw={"height_ratios": [1.6, 1]})
    a1.plot(al, F.real / N, color=BLUE, lw=0.9)
    a1.axhline(ctx["R2"] / N, color=ORANGE, lw=1.0)
    a1.text(0.005, ctx["R2"] / N + 0.03, f"final height R2(N)/N = {ctx['R2'] / N:.6f}", fontsize=8, color=INK)
    a2.plot(al, F.imag / N, color=BLUE, lw=0.9)
    a2.axhline(0, color=AXIS, lw=0.8)
    for j in st["jumps"]:
        a, q = j["a"], j["q"]
        if a == 0:
            continue
        y = np.interp(a / q, al, F.real / N)
        lab = f"{a}/{q}\npred {j['pred'] / N:+.4f}\nmeas {j['meas'] / N:+.4f}"
        off = (10, -8) if (a, q) == (1, 2) else (0, 26 if j["pred"] > 0 else -42)
        a1.annotate(lab, (a / q, y), xytext=off, textcoords="offset points",
                    ha="left" if (a, q) == (1, 2) else "center", fontsize=6.5, color=INK2,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    a1.set_ylabel("Re F(α) / N")
    a2.set_ylabel("Im F(α) / N")
    a2.set_xlabel("α")
    a1.set_title("(a) real part: a staircase of jumps N μ(q)$^2$cos(2πaN/q)/φ(q)$^2$ at each a/q "
                 "(labels: predicted vs measured jump, q ≤ 6; smaller steps at q = 7, 10, … follow the same rule)")
    a2.set_title("(b) imaginary part: jumps −N μ(q)$^2$sin(2πaN/q)/φ(q)$^2$ cancel in pairs a, q−a; "
                 "the curve is symmetric about 1/2 and returns to 0")
    # inset: a plateau
    ins = a1.inset_axes([0.62, 0.08, 0.34, 0.36])
    lo, hi = 0.777, 0.782                    # a stretch with no jump larger than ~1e-3 N
    sel = (al >= lo) & (al <= hi)
    ins.plot(al[sel], F.real[sel] / N, color=BLUE, lw=0.8)
    ins.set_title(f"plateau zoom, α in [{lo}, {hi}]: wander ~10$^{{-3}}$N", fontsize=7.5, loc="left")
    ins.tick_params(labelsize=6.5)
    a1.indicate_inset_zoom(ins, edgecolor=MUTED, lw=0.6)
    fig.suptitle("Fig F. F(α) = ∫$_0^α$ S(t)$^2$e(−Nt) dt over the whole circle (N = 10$^6$; the α = 0 "
                 "spike straddles 0 = 1 mod 1, so half of it sits at each end)", x=0.01, ha="left", fontsize=10.5,
                 fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_F_staircase.png")


# ================================================================ binary vs ternary (module 6)
def fig_T(ctx, out):
    rows = ctx["rows"]
    lab = [f"({r['Q']},{r['c']})" for r in rows]
    x = np.arange(len(rows))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 4.9))
    a1.plot(x, [r["P_m_over_R2"] for r in rows], "o-", color=BLUE, ms=6, mec=SURF, mew=1.2, lw=1.2,
            label="binary: P$_m$/R2 = ∫$_m$|S|$^2$ / R2(N)")
    a1.plot(x, [r["B3_main_999999"] for r in rows], "s-", color=ORANGE, ms=6, mec=SURF, mew=1.2, lw=1.2,
            label="ternary: B3/main, N' = 10$^6$−1")
    a1.plot(x, [r["B3_main_1000001"] for r in rows], "s", color=ORANGE, ms=6, mfc=SURF, mew=1.2, lw=0.8,
            ls=(0, (3, 2)), label="ternary: B3/main, N' = 10$^6$+1")
    a1.axhline(1.0, color=INK, lw=1.0)
    a1.text(len(rows) - 0.6, 1.06, "bound = signal", ha="right", fontsize=7.5)
    a1.set_yscale("log")
    a1.set_xticks(x, lab, rotation=60, fontsize=7.5)
    a1.set_xlabel("(Q, c)")
    a1.set_title("(a) the bounds: magnitude-only estimates of the minor arcs ÷ main term")
    a1.set_ylabel("bound ÷ main term (log)")
    a2.plot(x, [abs(r["I_m"]) / r["R2"] for r in rows], "o-", color=BLUE, ms=6, mec=SURF, mew=1.2, lw=1.2,
            label="binary: |I$_m$|/R2 (actual signed minor-arc contribution)")
    a2.plot(x, [abs(r["I_m_resid"]) / r["R2"] for r in rows], "o-", color=VIOLET, ms=5, mec=SURF, mew=1.0,
            lw=1.0, label="binary: |I$_m$ − leakage|/R2 = |(R2 − SS2·N) − (I$_M$ − model)|/R2 (model error, not a "
            "clean cancellation measurement)")
    a2.plot(x, [r["rp_scale"] / r["R2"] for r in rows], "-", color=MUTED, lw=1.0,
            label="random-phase yardstick P$_m$/√(N|m|) ÷ R2")
    a2.plot(x, [abs(r["I3_m_999999"]) / r["main3_999999"] for r in rows], "s-", color=ORANGE, ms=6, mec=SURF,
            mew=1.2, lw=1.2, label="ternary: |I3$_m$|/main, N' = 10$^6$−1")
    a2.set_yscale("log")
    a2.set_xticks(x, lab, rotation=60, fontsize=7.5)
    a2.set_xlabel("(Q, c)")
    a2.set_title("(b) the truth: what the minor arcs actually contribute")
    a2.set_ylabel("|signed contribution| ÷ main term (log)")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=7.5, ncol=1)
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=7.5, ncol=1)
    fig.suptitle("Fig T. Binary vs ternary at N = 10$^6$: the binary L$^2$ bound never closes; the ternary sup×L$^2$ "
                 "bound closes once Q swallows the tall excluded spikes and c pushes the arc-edge sidelobes down",
                 x=0.01, ha="left", fontsize=10.5,
                 fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_T_binary_vs_ternary.png")


# ================================================================ scaling (module 7)
def fig_G(ctx, out):
    rows = [r for r in ctx["scaling_rows"] if not r.get("probe")]
    probe = [r for r in ctx["scaling_rows"] if r.get("probe")]
    series = [("Q=20, c=10", BLUE, "o"), ("Q=50, c=10", ORANGE, "s"), ("Q=Qmax, c=10", AQUA, "D"),
              ("Q=Qmax, c=20", VIOLET, "^")]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.4))
    for name, col, mk in series:
        rs = sorted([r for r in rows if r["series"] == name], key=lambda r: r["N"])
        if not rs:
            continue
        Ns = np.array([r["N"] for r in rs])
        kw = dict(marker=mk, color=col, ms=6, mec=SURF, mew=1.2, lw=1.2)
        axes[0].plot(Ns, [r["P_m_over_R2"] for r in rs], label=name, **kw)
        axes[1].plot(Ns, [r["P_m"] / r["N"] for r in rs], label=name, **kw)
        axes[2].plot(Ns, [r[f"B3_main_{r['N'] - 1}"] for r in rs], label=f"{name}, N' = N−1", **kw)
        axes[2].plot(Ns, [r[f"B3_main_{r['N'] + 1}"] for r in rs], marker=mk, color=col, ms=6, mfc=SURF,
                     mew=1.2, lw=0.8, ls=(0, (3, 2)), label=f"{name}, N' = N+1")
        pr = [r for r in probe if r["series"] == name]
        for r in pr:
            axes[0].plot([r["N"]], [r["P_m_over_R2"]], marker="*", ms=13, color=col, mec=INK, mew=0.6, ls="none")
            axes[1].plot([r["N"]], [r["P_m"] / r["N"]], marker="*", ms=13, color=col, mec=INK, mew=0.6, ls="none")
    if probe:
        Np = probe[0]["N"]
        axes[0].plot([], [], marker="*", ms=11, color=MUTED, mec=INK, mew=0.6, ls="none",
                     label=f"smooth-N probe N = {Np:,} (SS2 = {probe[0]['SS2']:.2f})")
        axes[0].annotate(f"N = {Np:,} = 2·3²·5·7·11·13·17", (Np, min(r["P_m_over_R2"] for r in probe)),
                         xytext=(0.30, 0.04), textcoords="axes fraction", fontsize=7, color=INK2,
                         arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    from singular_series import SS2_even_array
    Ns_all = np.array(sorted({r["N"] for r in rows}))
    arr = SS2_even_array(int(Ns_all.max()))
    axes[1].plot(Ns_all, [arr[: Nn + 1].max() for Nn in Ns_all], color=INK, lw=1.0,
                 label="largest SS2(n) over even n ≤ N")
    axes[1].plot(Ns_all, [arr[Nn] for Nn in Ns_all], color=INK2, lw=0.8, ls=(0, (1, 2)),
                 label="SS2(N) at the tabulated N")
    for ax in (axes[0], axes[2]):
        ax.axhline(1.0, color=INK, lw=1.0)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("N (log scale)")
    axes[0].set_ylabel("P$_m$ / R2(N)")
    axes[0].set_title("(a) binary L$^2$ wall: ∫$_m$|S|$^2$ ÷ R2(N)")
    axes[1].set_ylabel("P$_m$ / N")
    axes[1].set_title("(b) arithmetic removed: P$_m$/N vs SS2, which it must beat")
    axes[2].set_ylabel("B3 / main")
    axes[2].set_title("(c) ternary: sup$_m$|S|·∫|S|$^2$ ÷ SS3(N')N'$^2$/2")
    for ax, nc in zip(axes, (1, 1, 2)):
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=7.2, ncol=nc)
    fig.suptitle("Fig G. Scaling in N (no fits). Fixed (Q, c): both ratios creep up with log N. Q grown as far as "
                 "2cQ$^2$ < N allows: the ternary ratio falls (once c is large enough), the binary one does not",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_G_scaling.png")


# ================================================================ random model (module 8)
def fig_H(ctx, out):
    R = ctx["rand"]
    env_p, env_r = ctx["env_max"], R["env_max"]
    x = (np.arange(len(env_p)) + 0.5) / len(env_p)
    fig, (a1, a3, a2) = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw={"width_ratios": [2.0, 1.2, 1]})
    # the random envelope is drawn on top: the prime forest shows as the blue spikes sticking out above it
    a1.semilogy(x, env_p, color=BLUE, lw=0.6, label="primes: max |S| per bin")
    a1.semilogy(x, env_r, color=ORANGE, lw=0.6, label="Cramér random set: max |S$_R$| per bin")
    a1.set_ylim(80, 4e6)
    a1.set_xlabel("α")
    a1.set_ylabel("max |S| per bin of width 1/16384 (log)")
    a1.set_title("(a) overlay: only the α = 0 spike survives randomisation")
    a1.legend(loc="upper center", ncol=2)
    a3.plot(x, ctx["env_med"], color=BLUE, lw=0.6, label="primes")
    a3.plot(x, R["env_med"], color=ORANGE, lw=0.6, label="random")
    a3.set_xlabel("α")
    a3.set_ylabel("median |S| per bin")
    a3.set_title("(b) the sea: the random set sits uniformly higher")
    a3.legend(loc="upper center", ncol=2)
    psi = ctx["psi"]
    fr = [r for r in ctx["farey"] if r["q"] <= 12]
    xs = np.array([r["a"] / r["q"] for r in fr])
    a2.scatter(xs, [r["absS"] / psi for r in fr], s=16, color=BLUE, zorder=3, label="primes |S(a/q)|/ψ",
               edgecolor=SURF, linewidth=0.6)
    a2.scatter(xs, [r["absS"] / R["psi"] for r in R["farey"]], s=16, color=ORANGE, marker="s", zorder=3,
               label="random |S$_R$(a/q)|/ΣW", edgecolor=SURF, linewidth=0.6)
    a2.axhline(math.sqrt(R["L2"]) / R["psi"], color=INK2, lw=0.8)
    a2.text(0.99, math.sqrt(R["L2"]) / R["psi"] * 1.25, "random-walk level", ha="right", fontsize=7, color=INK2,
            transform=a2.get_yaxis_transform())
    a2.set_yscale("log")
    a2.set_xlabel("a/q  (q ≤ 12)")
    a2.set_ylabel("|S(a/q)| ÷ total weight")
    a2.set_title("(c) at the rationals themselves")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=7.5, ncol=1)
    fig.suptitle(f"Fig H. Cramér model (n in R with probability 1/log n, weight log n): R2$_R$(N) = {R['R2']:,.0f} ≈ N, "
                 f"vs primes R2(N) = {ctx['R2']:,.0f} ≈ SS2(N)·N", x=0.01, ha="left", fontsize=10.5,
                 fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_H_random_model.png")


# ================================================================ comet (module 9)
def fig_K(ctx, out):
    c = ctx["comet"]
    N, G, SS2, pred = c["N"], c["G"], c["SS2"], c["pred"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    # draw all classes in one interleaved pass, so N = 1 and N = 2 (mod 3), which share the
    # same bands, both stay visible instead of one painting over the other
    cls_col = np.array([BLUE, ORANGE, AQUA])
    perm = np.random.default_rng(3).permutation(len(N))
    axes[0].scatter(N[perm], G[perm], s=0.6, c=cls_col[N[perm] % 3], lw=0, rasterized=True)
    for r, col in enumerate(cls_col):
        axes[0].scatter([], [], s=6, color=col, label=f"N ≡ {r} (mod 3)")
    axes[0].set_xlabel("even N")
    axes[0].set_ylabel("G(N) = # {p ≤ p′ : p + p′ = N}")
    axes[0].set_title("(a) the comet, coloured by N mod 3")
    axes[0].legend(loc="upper left", markerscale=1.5, fontsize=7.5)
    order = np.argsort(SS2)
    sc = axes[1].scatter(N[order], G[order], s=0.6, c=np.log(SS2[order]), cmap=SEQ_BLUE, lw=0, rasterized=True)
    axes[1].scatter(N, pred, s=0.06, color=ORANGE, lw=0, rasterized=True, alpha=0.5)
    axes[1].scatter([], [], s=4, color=ORANGE, label="prediction SS2(N)·½∫$_2^{N-2}$dt/(log t log(N−t))")
    axes[1].legend(loc="upper left", markerscale=3, fontsize=7.5)
    cb = fig.colorbar(sc, ax=axes[1], pad=0.01)
    cb.set_label("log SS2(N)", color=INK2)
    cb.outline.set_edgecolor(AXIS)
    axes[1].set_xlabel("even N")
    axes[1].set_title("(b) coloured by log SS2(N), prediction overlaid (orange)")
    ratio = G / pred
    s = N >= 1000
    axes[2].scatter(N[s], ratio[s], s=0.6, color=BLUE, lw=0, rasterized=True)
    axes[2].axhline(1.0, color=INK, lw=1.0)
    edges = np.array([1e3, 5e3, 1e4, 2e4, 5e4, 1e5, 1.5e5, 2.0e5 + 1])
    centres, means = [], []
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        b = (N >= lo_) & (N < hi_)
        centres.append(np.median(N[b]))
        means.append(ratio[b].mean())
    axes[2].plot(centres, means, "o-", color=ORANGE, ms=6, mec=SURF, mew=1.2, lw=1.4, label="binned mean (no fit)")
    axes[2].legend(loc="upper right", fontsize=7.5)
    ylo, yhi = 0.80, 1.20
    out_mask = (ratio[s] < ylo) | (ratio[s] > yhi)
    clipped = int(np.sum(out_mask))
    axes[2].set_ylim(ylo, yhi)
    if clipped:
        axes[2].text(0.99, 0.02, f"{clipped} points outside [{ylo}, {yhi}], all at N ≤ {N[s][out_mask].max():,}",
                     ha="right", transform=axes[2].transAxes, fontsize=7, color=INK2)
    axes[2].set_xlabel("even N ≥ 1000")
    axes[2].set_ylabel("G(N) / prediction")
    axes[2].set_title("(c) observed ÷ predicted")
    kfmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1000:.0f}k" if v else "0")
    for ax in axes:
        ax.set_xticks([0, 50000, 100000, 150000, 200000])
        ax.xaxis.set_major_formatter(kfmt)
    fig.suptitle("Fig K. Goldbach's comet for even N ≤ 2·10$^5$ (FFT autoconvolution of the prime indicator)",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_K_comet.png")


# ================================================================
ALL = [fig_A, fig_B, fig_C, fig_D, fig_P, fig_S, fig_E, fig_F, fig_T, fig_G, fig_H, fig_K]


def make_all(ctx, out: Path) -> None:
    out = Path(out)
    out.mkdir(exist_ok=True)
    for f in ALL:
        f(ctx, out)
