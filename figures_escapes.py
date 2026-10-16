"""Figures for the escape-route experiments (E0-E5). Takes the dict produced by run_escapes.compute()."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import figures as base  # applies the shared rcParams / palette
from figures import AQUA, AXIS, BLUE, GRID, INK, INK2, MUTED, ORANGE, SURF, VIOLET, YELLOW, _save, plt

PAIR_COL = {(150, 20): BLUE, (50, 10): ORANGE, (20, 10): AQUA, (223, 10): VIOLET}


def _plain_log(ax, axis="y"):
    """Log axis labelled with plain numbers (1, 2, 3, 5 per decade) instead of 2×10^0."""
    from matplotlib.ticker import FuncFormatter

    def minor(v, _):
        if v <= 0:
            return ""
        lead = round(v / 10 ** math.floor(math.log10(v)))
        return f"{v:g}" if lead in (2, 3, 5) else ""
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_minor_formatter(FuncFormatter(minor))


def _undecided(ax):
    ax.axhspan(0.99, 1.01, color=GRID, zorder=0)
    ax.axhline(1.0, color=INK, lw=1.0)


# ================================================================ E0
def fig_X0(res, out):
    R2 = res["R2"]
    e0 = res["e0"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw={"width_ratios": [1.15, 1.1, 1.1]})
    ax = axes[0]
    o = e0[(150, 20)]
    wk = o["walk"] / R2
    IM = o["I_M"] / R2
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(IM * np.cos(th), IM * np.sin(th), color=MUTED, lw=1.0, label="|z| = I_M (to spoil the count)")
    ax.plot([0, -o["P_m"] / R2], [0, 0], color=ORANGE, lw=2.0, label=f"adversarial path: same step lengths, aligned (P_m = {o['P_m'] / R2:.2f}·R2)")
    ax.plot(wk.real, wk.imag, color=BLUE, lw=1.0, label=f"actual walk (6.1 M phasors), ends at {o['I_m'].real / R2:+.4f}·R2")
    ax.set_aspect("equal")
    ax.set_xlim(-3.1, 1.2)
    ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("Re (units of R2)")
    ax.set_ylabel("Im (units of R2)")
    ax.set_title("(a) minor-arc phasors head to tail, (Q, c) = (150, 20)")
    ins = ax.inset_axes([0.13, 0.06, 0.36, 0.36])
    ins.plot(wk.real * 100, wk.imag * 100, color=BLUE, lw=0.6)
    ins.plot([0], [0], "o", color=INK, ms=3.5)
    ins.plot([wk[-1].real * 100], [wk[-1].imag * 100], "s", color=ORANGE, ms=4)
    ins.set_title("walk ×100 (● start, ■ end)", fontsize=6.5, loc="left", pad=2)
    ins.set_aspect("equal")
    ins.tick_params(labelsize=6)
    ax.legend(loc="upper left", fontsize=6.8)
    ax = axes[1]
    best = min(o["model"], key=lambda m: m["P_E"])
    bestb = min(o["model"], key=lambda m: m["bound"])
    bars = [("actual |I_m|", abs(o["I_m"].real)), ("cells of width 1/N, random phases: rms", o["random_phase_rms"]),
            (f"sea after subtracting spikes, Q' = {best['Qp']} (not a bound)", best["P_E"]),
            ("P_m = ∫_m|S|$^2$ (magnitude bound)", o["P_m"]), ("Hölder √(|m|∫_m|S|$^4$)", o["holder"]),
            (f"best model-subtraction bound (Q' = {bestb['Qp']})", bestb["bound"]), ("sup_m|S| · ∫_m|S|", o["sup_L1"])]
    y = np.arange(len(bars))
    ax.barh(y, [b[1] / R2 for b in bars], color=[AQUA, AQUA, MUTED, BLUE, BLUE, BLUE, BLUE], height=0.55)
    ax.axvline(IM, color=ORANGE, lw=1.6)
    ax.text(IM * 1.1, -0.45, "I_M (signal)", color=INK, fontsize=8, va="center")
    ax.set_yticks(y, [b[0] for b in bars], fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel("units of R2 (log)")
    ax.set_title("(b) the bound ladder at (150, 20): every magnitude-only bound (blue) exceeds I_M")
    ax = axes[2]
    for (Q, c), col in ((((150, 20)), BLUE), (((20, 10)), AQUA)):
        o = e0[(Q, c)]
        nd = [d["N_Delta"] for d in o["dial"]]
        W = [d["W"] for d in o["dial"]]
        ax.plot(nd, W, "o-", color=col, ms=5, mec=SURF, mew=1.0, label=f"({Q}, {c}): Σ|window integrals| ÷ I_M")
        pm = o["P_m"] / o["I_M"]
        xs = np.geomspace(1.0, max(nd), 50)          # the equal-cell yardstick only makes sense for N·Δ >= 1
        ax.plot(xs, pm / np.sqrt(xs), color=col, lw=0.8, ls=(0, (3, 2)))
        ax.axhline(o["dial_floor"], color=col, lw=0.6)
    _undecided(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("window width N·Δ (Δ in units of 1/N)")
    ax.set_ylabel("W(Δ)")
    ax.set_title("(c) how much phase knowledge would suffice")
    ax.plot([], [], color=INK2, lw=0.8, ls=(0, (3, 2)), label="equal-cell yardstick (P_m/I_M)/√(NΔ), pair colours")
    ax.plot([], [], color=INK2, lw=0.6, label="thin line: W(∞), whole gaps between arcs")
    ax.legend(loc="lower left", fontsize=7)
    a, b = e0[(150, 20)], e0[(20, 10)]
    lo, hi = sorted([a["P_m"] / a["I_M"], b["P_m"] / b["I_M"]])
    elo, ehi = sorted([abs(a["I_m"].real) / a["I_M"], abs(b["I_m"].real) / b["I_M"]])
    fig.suptitle(f"Fig X0. The obstacle, measured at N = 10$^6$: the minor arcs' path is {lo:.1f}–{hi:.1f}× the signal, "
                 f"its endpoint {100 * elo:.0f}–{100 * ehi:.0f}% of it", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_X0_obstacle.png")


# ================================================================ E1
def fig_X1(res, out):
    e1 = res["e1"]
    pairs = res["pairs"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [1.1, 1.3, 1.0]})
    ax = axes[0]
    edges = np.linspace(-6, 0.5, 131)
    ctr = 10 ** (0.5 * (edges[1:] + edges[:-1]))
    for pr in ((20, 10), (50, 10), (150, 20)):
        h = e1[pr]["ratio_hist"]
        ax.step(ctr, h, where="mid", color=PAIR_COL[pr], lw=1.2, label=f"{pr}: worst {e1[pr]['quantiles'][1.0]:.3f}")
    ax.axvline(1.0, color=INK, lw=1.0)
    ax.set_xscale("log")
    ax.set_xlim(1e-5, 3)
    ax.set_xlabel("|I_m(N)| / I_M(N)  over the 250,001 even N in [X/2, X]")
    ax.set_ylabel("number of targets per bin")
    ax.set_title("(a) every target's minor-arc output is far below its signal")
    ax.legend(loc="upper left", fontsize=7.5)
    ax = axes[1]
    x = np.arange(len(pairs))
    n = e1[pairs[0]]["n_targets"]
    for key, col, mk, lab in (("bound_textbook", MUTED, "v", "sup$^2$·ΣΛ$^2$ (textbook: the 1937–38 proofs)"),
                              ("bound_natural_sup", ORANGE, "s", "sup$^2$·P_m (sup pulled out, minor-arc L$^2$)"),
                              ("bound_4th", BLUE, "o", "∫_m|S|$^4$ (true 4th moment)"),
                              ("bound_oracle", AQUA, "D", "Σ_targets I_m(N)$^2$ (oracle, no theorem)")):
        ax.plot(x, [max(e1[p][key], 0.5) for p in pairs], mk + "-", color=col, ms=5, mec=SURF, mew=1.0, lw=1.0, label=lab)
    ax.axhline(n, color=INK, lw=1.0)
    ax.text(len(pairs) - 0.5, n * 1.2, "number of targets", ha="right", fontsize=7.5)
    ax.set_yscale("log")
    ax.set_xticks(x, [f"({p[0]},{p[1]})" for p in pairs], rotation=60, fontsize=7)
    ax.set_ylabel("allowed failures F (log); actual F = 0 everywhere")
    ax.set_title("(b) Bessel + Chebyshev: how many targets the method may lose")
    ax.legend(loc="lower left", fontsize=7)
    ax = axes[2]
    keys = ["targets", "other even <= 2X", "odd <= 2X", "> 2X"]
    names = {"targets": "targets (even N in [X/2, X])", "other even <= 2X": "other even N ≤ 2X", "odd <= 2X": "odd N ≤ 2X",
             "> 2X": "N > 2X"}
    cols = [BLUE, AQUA, ORANGE, MUTED]
    sp = [(150, 20), (100, 20), (50, 10), (20, 10), (5, 5), (223, 10)]
    bottom = np.zeros(len(sp))
    for kk, col in zip(keys, cols):
        vals = np.array([e1[p]["shares"][kk] for p in sp])
        ax.bar(np.arange(len(sp)), vals, bottom=bottom, color=col, width=0.6, edgecolor=SURF, linewidth=1.0, label=names[kk])
        bottom += vals
    ax.set_xticks(np.arange(len(sp)), [f"({p[0]},{p[1]})" for p in sp], fontsize=7.5)
    ax.set_ylabel("share of Σ_N I_m(N)$^2$ = ∫_m|S|$^4$")
    ax.set_title("(c) where the shared energy budget goes")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=7.2)
    fig.suptitle("Fig X1. Averaging over N: one energy budget Σ_N I_m(N)$^2$ = ∫_m|S|$^4$ is shared by all targets, "
                 "so few can fail", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_X1_exceptional.png")


# ================================================================ E5
def fig_X5(res, out):
    o = res["e1"][(223, 10)]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4))
    ax = axes[0]
    ss2, pr = o["SS2_targets"], o["pm_ratio_targets"]
    ax.scatter(ss2, pr, s=0.5, color=BLUE, lw=0, rasterized=True)
    tr = o["transmitted"]
    cl = [t for t in tr if t[1] < 0.99]
    ud = [t for t in tr if t[1] >= 0.99]
    ax.scatter([t[2] for t in cl], [t[1] for t in cl], s=18, color=ORANGE, zorder=3, edgecolor=SURF, linewidth=0.6,
               label=f"{len(cl)} targets close (P_m/I_M(N) < 0.99): R2(N) > 0 by magnitudes alone")
    ax.scatter([t[2] for t in ud], [t[1] for t in ud], s=18, facecolor="none", edgecolor=ORANGE, linewidth=1.0, zorder=3,
               label=f"{len(ud)} undecided with ratio in [0.99, 1) (still grid-certified)")
    hi = (pr >= 1) & (pr <= 1.01)
    ax.scatter(ss2[hi], pr[hi], s=16, marker="s", facecolor="none", edgecolor=ORANGE, linewidth=0.8, zorder=3,
               label=f"{int(hi.sum())} undecided with ratio in [1, 1.01] (no grid certificate)")
    for t in res["e5_trunc"]:
        from singular_series import SS2_product
        ax.scatter([SS2_product(t["N"])], [t["ratio"]], s=70, marker="D", facecolor="none", edgecolor=INK, linewidth=1.0,
                   zorder=4)
        if t["ratio"] < 1:
            ax.annotate(f"{t['N']:,}: S cut at N, Q = {t['Q']}", (SS2_product(t["N"]), t["ratio"]), xytext=(-8, -16),
                        ha="right", textcoords="offset points", fontsize=7, color=INK2,
                        bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.85))
        else:                                   # far from X: lead the label out of the cloud into the empty corner
            ax.annotate(f"{t['N']:,}: S cut at N, Q = {t['Q']}", (SS2_product(t["N"]), t["ratio"]), xytext=(1.45, 1.35),
                        textcoords="data", ha="left", fontsize=7, color=INK2,
                        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    _undecided(ax)
    ax.set_yscale("log")
    ax.set_xlabel("SS2(N)")
    ax.set_ylabel("P_m / I_M(N)  at (Q, c) = (223, 10)")
    _plain_log(ax)
    ax.plot([], [], "D", mfc="none", mec=INK, ls="none", label="open diamonds: S truncated at N, arcs (Qmax(N, 10), 10)")
    ax.set_title("(a) the wall is the same for every target; the signal grows with SS2(N)")
    ax.legend(loc="upper right", fontsize=7.5)
    ax = axes[1]
    sc = res["e5_scaling"]
    for c, col, mk in ((10, AQUA, "D"), (20, VIOLET, "^")):
        rs = [r for r in sc if r["c"] == c]
        xs = [r["X"] for r in rs]
        ax.plot(xs, [r["need"] for r in rs], mk + "-", color=col, ms=6, mec=SURF, mew=1.0, lw=1.2,
                label=f"SS2 a target near X would need (Q = Qmax, c = {c})")
        for r in rs:
            ax.annotate(f"{r['n_transmitted']}", (r["X"], r["need"]), xytext=(0, -13 if c == 10 else 7),
                        textcoords="offset points", fontsize=7, color=col, ha="center")
    xs = sorted({r["X"] for r in sc})
    ax.plot(xs, [max(r["max_SS2"] for r in sc if r["X"] == x_) for x_ in xs], "s", color=INK, ms=6, ls="none",
            label="largest SS2(N) among the targets at each X")
    ax.set_xscale("log")
    ax.set_xlabel("X")
    ax.set_ylabel("singular series")
    ax.set_title("(b) the SS2 a target needs rises at every X; the largest available rises slowly\n(labels: targets with P_m < I_M)")
    ax.legend(loc="upper left", fontsize=7.5)
    fig.suptitle("Fig X5. Single-N transmission with no relaxation: a few very smooth N get through, a falling fraction as X grows",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_X5_transmission.png")


# ================================================================ E2
def fig_X2(res, out):
    e2 = res["e2"]
    L = e2["L"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    ax = axes[0]
    r = e2["rows"][(150, 20)]
    K = np.arange(len(r["joint_ratio"]))
    ax.plot(K, r["joint_ratio"], "o-", color=BLUE, ms=6, mec=SURF, mew=1.0, lw=1.2, label="joint ∫_m|S|$^2$|G|$^K$ (sharpest magnitude bound)")
    ax.plot(K, r["rayleigh_ratio"], color=MUTED, lw=1.0, label="Rayleigh prediction (|S|, |G| independent)")
    ax.plot(K, r["holder_ratio"], "s-", color=YELLOW, ms=5, mec=SURF, mew=1.0, lw=1.0, label="Hölder")
    ax.plot(K, r["sup_ratio"], "v-", color=ORANGE, ms=5, mec=SURF, mew=1.0, lw=1.0, label="sup_m|G|$^K$·P_m")
    for C, ls in ((1, "-"), (2, (0, (4, 2))), (4, (0, (1, 2)))):
        vals = [v[0] for v in r["theorem"][C]["split"]]
        ax.plot(K[2:], vals[2:], color=VIOLET, lw=1.2, ls=ls, label=f"theorem-shaped (λ-split, sieve constant C = {C})")
    r2 = e2["rows"][(20, 10)]
    ax.plot(K, r2["joint_ratio"], "o", color=AQUA, mfc=SURF, ms=5, mew=1.0, ls=(0, (3, 2)), lw=0.8, label="joint at (20, 10)")
    _undecided(ax)
    ax.set_yscale("log")
    ax.set_xlabel("K (number of powers of 2 added)")
    ax.set_ylabel("bound ÷ major part I_M,K")
    ax.set_title("(a) (Q, c) = (150, 20): one power already closes the joint bound")
    ax.legend(loc="lower left", fontsize=6.6)
    # (b) how much of B_K comes from the shoulders of the arcs at 0 and 1/2
    ax = axes[1]
    for pr, ls in (((150, 20), "-"), ((20, 10), (0, (3, 2)))):
        sh = e2["rows"][pr]["shoulder_share"]
        for R, col in ((60, VIOLET), (200, ORANGE), (1000, YELLOW)):
            ax.plot(K, 100 * sh[R], "o", color=col, ms=4.5, mec=SURF, mew=0.8, ls=ls, lw=1.1,
                    label=f"within {R}/N of 0 or 1/2, {pr}")
    ax.set_xlabel("K (number of powers of 2 added)")
    ax.set_ylabel("share of B_K = ∫_m|S|$^2$|G|$^K$ (%)")
    ax.set_ylim(0, 100)
    ax.set_title("(b) as K grows, B_K's weight moves to the shoulders of the arcs at 0 and 1/2")
    ax.legend(loc="upper left", fontsize=6.8, ncol=2)
    # (c) large values of G
    ax = axes[2]
    lams = sorted(e2["rows"][(150, 20)]["measE"])
    for pr, col in (((150, 20), BLUE), ((20, 10), AQUA)):
        me = e2["rows"][pr]["measE"]
        ax.plot(lams, [me[l_] for l_ in lams], "o-", color=col, ms=6, mec=SURF, mew=1.0, label=f"meas{{α in minor arcs : |G| ≥ λL}}, {pr}")
    ax.set_yscale("log")
    ax.set_xlabel("λ")
    ax.set_ylabel("measure (log)")
    ax.set_title("(c) G is large only on a tiny part of the minor arcs")
    ax.legend(loc="upper right", fontsize=7.2)
    fig.suptitle("Fig X2. Linnik's relaxation: N = p + p′ + 2$^{k_1}$ + … + 2$^{k_K}$ at N = 10$^6$", x=0.01, ha="left",
                 fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_X2_linnik.png")


# ================================================================ E3
def fig_X3(res, out):
    e3 = {s["name"]: s for s in res["e3"]}
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.4))
    ax = axes[0]
    rs = np.arange(1, 7)
    for pr in ((150, 20), (50, 10), (20, 10)):
        col = PAIR_COL[pr]
        j = [e3[f"Omega<={r}"]["rows"][pr]["joint_ratio"] for r in rs]
        cs = [e3[f"Omega<={r}"]["rows"][pr]["CS_ratio"] for r in rs]
        ax.plot(rs, j, "o-", color=col, ms=6, mec=SURF, mew=1.0, lw=1.2, label=f"Ω ≤ r, joint, {pr}")
        ax.plot(rs, cs, "o", color=col, mfc=SURF, ms=6, mew=1.0, ls=(0, (3, 2)), lw=0.8, label=f"Ω ≤ r, Cauchy–Schwarz, {pr}")
    pr = (150, 20)
    for name, mk, dx in (("Chen-shaped: Omega<=2, (m,6)=1", "*", 0.0), ("odd Omega<=2", "P", 0.08), ("Omega<=2, (m,10)=1", "X", 0.16),
                         ("odd Omega<=3", "P", 1.08), ("Omega<=3, (m,10)=1", "X", 1.16)):
        ax.plot([2 + dx], [e3[name]["rows"][pr]["joint_ratio"]], mk, color=BLUE, ms=9, mec=INK, mew=0.5, ls="none")
    for name in [n for n in e3 if n.startswith("Bernoulli")]:
        rr = 2 if "Omega<=2" in name else 3
        ax.plot([rr - 0.12], [e3[name]["rows"][pr]["joint_ratio"]], "x", color=MUTED, ms=9, mew=1.6, ls="none")
    ax.plot([], [], "*", color=BLUE, mec=INK, mew=0.5, ms=9, ls="none", label="Chen-shaped Ω ≤ 2, (m,6) = 1 (150, 20)")
    ax.plot([], [], "P", color=BLUE, mec=INK, mew=0.5, ms=8, ls="none", label="odd members only (150, 20)")
    ax.plot([], [], "X", color=BLUE, mec=INK, mew=0.5, ms=8, ls="none", label="coprime to 10 = N's primes (150, 20)")
    ax.plot([], [], "x", color=MUTED, ms=9, mew=1.6, ls="none", label="random set of the same density (150, 20)")
    _undecided(ax)
    ax.set_yscale("log")
    ax.set_xlabel("r  (second summand has at most r prime factors)")
    ax.set_ylabel("minor-arc bound ÷ major part")
    _plain_log(ax)
    ax.set_title("(a) N = p + m, Ω(m) ≤ r: joint closes at r = 2–3, Cauchy–Schwarz at 3–4; random sets never do")
    ax.legend(loc="upper right", fontsize=6.6, ncol=2)
    ax = axes[1]
    zs = sorted(int(n.split("-")[0]) for n in e3 if n.endswith("-rough"))
    for pr in ((150, 20), (20, 10)):
        col = PAIR_COL[pr]
        ax.plot(zs, [e3[f"{z}-rough"]["rows"][pr]["joint_ratio"] for z in zs], "o-", color=col, ms=4.5, mec=SURF, mew=0.8,
                lw=1.1, label=f"joint, {pr}")
        ax.plot(zs, [e3[f"{z}-rough"]["rows"][pr]["CS_ratio"] for z in zs], "o", color=col, mfc=SURF, ms=4.5, mew=0.9,
                ls=(0, (3, 2)), lw=0.8, label=f"Cauchy–Schwarz, {pr}")
    _undecided(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("z  (second summand has no prime factor below z; z ≈ √X gives the primes again)")
    ax.set_ylabel("minor-arc bound ÷ major part")
    from escapes_report import largest_closing

    def zstar(pr):
        last, _ = largest_closing([(z, e3[f"{z}-rough"]["rows"][pr]["joint_ratio"]) for z in zs])
        return last[0] if last else "none"
    ax.set_title(f"(b) stopping the sieve at z: the joint bound closes up to z = {zstar((150, 20))} at (150, 20), "
                 f"{zstar((20, 10))} at (20, 10)")
    ax.legend(loc="lower right", fontsize=7.2)
    fig.suptitle("Fig X3. Relaxing one prime to an almost-prime (the direction of Rényi and Chen — whose theorems need the sieve)",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_X3_almost_primes.png")


# ================================================================ E4
def fig_X4(res, out):
    e4 = res["e4"]
    X = res["X"]
    Nt = res["nters"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    pf = e4["profile"]
    ax.semilogy(pf["Nbeta"], pf["S"], color=BLUE, lw=1.0, label="sharp cutoff: |S(β)|/ψ")
    ax.semilogy(pf["Nbeta"], pf["S_eta"], color=INK, lw=1.0, label="smooth bump: |S_η(β)|/ψ_η")
    far = np.asarray(pf["Nbeta"]) >= 10
    rip = float(np.sqrt(np.mean(np.asarray(pf["S_eta"])[far] ** 2)))
    ax.plot([10, 30], [rip, rip], color=INK2, lw=0.9, ls=(0, (1, 2)),
            label="rms |S_η|, Nβ ≥ 10 (primes' fluctuation)")
    for c in (5, 10, 20):
        ax.axvline(c, color=MUTED, lw=0.8)
        ax.text(c + 0.2, 1.4e-4, f"c = {c}", fontsize=7, color=INK2)
    ax.set_ylim(1e-4, 2)
    ax.set_xlabel("N·β next to α = 0")
    ax.set_ylabel("normalised |S| (log)")
    ax.set_title("(a) the sharp window rings; the bump's sidelobes are far lower")
    ax.legend(loc="upper right", fontsize=7.5, frameon=True, framealpha=0.92, edgecolor="none", facecolor=SURF)
    ax = axes[1]
    rows = e4["rows"]
    Qs = [5, 10, 20, 50, 100, 150]
    styles = ((Nt[0], "-", "N' = X−1"), (Nt[2], (0, (1, 2)), "N' = 1.5X−1"))
    for c, col in ((5, BLUE), (10, ORANGE), (20, AQUA)):
        for Np, ls, lab in styles:
            ax.plot(Qs, [rows[(Q, c)][f"B3_sharp_{Np}"] for Q in Qs], "o", color=col, ms=5, mec=SURF, mew=1.0, ls=ls,
                    lw=1.1, label=f"sharp, c = {c}" if Np == Nt[0] else None)
    for Np, ls, lab in styles:   # the smooth sup does not depend on c: one curve per N'
        ax.plot(Qs, [rows[(Q, 5)][f"B3_eta_{Np}"] for Q in Qs], "s", color=INK, mfc=SURF, ms=6, mew=1.2, ls=ls, lw=1.0,
                label="smooth η (every c)" if Np == Nt[0] else None)
    ax.plot([], [], color=INK2, ls="-", label="solid: N' = X−1")
    ax.plot([], [], color=INK2, ls=(0, (1, 2)), label="dotted: N' = 1.5X−1")
    _undecided(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Q")
    ax.set_ylabel("ternary B3 / main")
    ax.set_title("(b) ternary: smoothing removes the c = 5 plateau")
    ax.text(0.02, 0.02, "sharp: all three c coincide for Q ≤ 20;\nc = 10 and 20 coincide up to Q = 100 (orange hidden)",
            transform=ax.transAxes, fontsize=7, color=INK2, va="bottom")
    ax.legend(loc="upper right", fontsize=6.8, ncol=2)
    ax = axes[2]
    for c, col in ((10, ORANGE), (20, AQUA)):
        Qs = [5, 10, 20, 50, 100, 150]
        ax.plot(Qs, [res["e1"][(Q, c)]["P_m"] / res["R2"] for Q in Qs], "o-", color=col, ms=5, mec=SURF, mew=1.0, lw=1.2,
                label=f"sharp, c = {c}")
        ax.plot(Qs, [rows[(Q, c)]["P_m_over_R"] for Q in Qs], "s", color=col, mfc=SURF, ms=5, mew=1.0, ls=(0, (3, 2)), lw=0.9,
                label=f"smooth η, c = {c}")
    ax.axhline(1.0, color=INK, lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("Q")
    ax.set_ylabel("binary P_m / R")
    ax.set_title("(c) binary: the wall does not move")
    ax.legend(loc="upper right", fontsize=7.5)
    fig.suptitle("Fig X4. Control: a smooth weight (apodization) cures the arc-edge leakage, not the L$^2$ wall",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    _save(fig, out, "fig_X4_window.png")


ALL = [fig_X0, fig_X1, fig_X5, fig_X2, fig_X3, fig_X4]


def make_all(res, out: Path) -> None:
    out = Path(out)
    out.mkdir(exist_ok=True)
    for f in ALL:
        f(res, out)
