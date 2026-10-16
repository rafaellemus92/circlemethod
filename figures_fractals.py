"""Figures for the fractal lens (M1-M5, A1-A3). Takes the dict produced by run_fractals.compute()."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import figures as base  # noqa: F401  (applies the shared rcParams / palette)
from figures import AQUA, AXIS, BLUE, GRID, INK, INK2, MAGENTA, MUTED, ORANGE, RED, SURF, VIOLET, YELLOW, _save, plt

LAGS = [2**j for j in range(22)]
PAIR_COL = {(5, 5): VIOLET, (10, 10): MAGENTA, (20, 10): AQUA, (50, 10): ORANGE, (100, 20): YELLOW, (150, 20): BLUE,
            (223, 10): RED}
MEAS_STYLE = {"primes": dict(color=BLUE, ls="-"), "hybrid": dict(color=ORANGE, ls=(0, (4, 2))),
              "cramer": dict(color=MUTED, ls=(0, (1, 1.5)))}
MEAS_LABEL = {"primes": "primes", "hybrid": "surrogate H₁₀₀₀ (spikes q ≤ 1000 + sea)", "cramer": "Cramér model"}


def _suptitle(fig, text):
    fig.suptitle(text, x=0.01, ha="left", fontsize=10.5, fontweight="bold")


def _xl(res):
    return np.array(LAGS[:-1]) * res["X"] / res["M"]


# ================================================================ M1
def fig_M1(res, out):
    X, psi = res["X"], res["psi"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    ax = axes[0, 0]
    xs = np.linspace(-10, 10, 2001)
    ax.plot(xs, np.abs(np.sinc(xs)), color=INK, lw=1.0, label="|sinc x| (the spike at 0)")
    cols = [BLUE, AQUA, ORANGE, RED]
    for z, col in zip(res["m1_zoom"], cols):
        lo, hi = z["class"]
        ax.plot(z["centres"], z["median"], color=col, lw=1.1, label=f"q in [{lo}, {hi}] (median of {z['n_q']} q)")
        ax.fill_between(z["centres"], z["p10"], z["p90"], color=col, alpha=0.15, lw=0)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 20)
    ax.set_xlabel("x = N(α − a/q)   (the width is never rescaled)")
    ax.set_ylabel("|S(α)| φ(q)/ψ   (height rescaled by φ(q)/ψ)")
    ax.set_title("(a) data collapse onto the spike at 0, and its failure (10–90% bands)")
    ax.legend(loc="upper right", fontsize=7.5)
    ax = axes[0, 1]
    rows = [r for r in res["m1_rows"] if r["q"] >= 3]
    ph = np.array([r["phi"] for r in rows], dtype=float)
    ax.scatter(ph, [r["mean"] for r in rows], s=3, color=BLUE, lw=0, label="mean of |S(a/q)| φ/ψ over a (measured)", rasterized=True)
    ax.scatter(ph, [r["sd"] for r in rows], s=3, color=ORANGE, lw=0, label="sd over a (measured)", rasterized=True)
    sub = rows[:: max(1, len(rows) // 400)]
    ax.scatter([r["phi"] for r in sub], [r["beck_mean"] for r in sub], s=12, facecolor="none", edgecolor=INK, lw=0.5,
               label="Beckmann model (exact second moments, no fit)")
    ax.scatter([r["phi"] for r in sub], [r["beck_sd"] for r in sub], s=12, facecolor="none", edgecolor=INK, lw=0.5)
    pr = sorted([r for r in rows if "rice_mean" in r], key=lambda r: r["phi"])
    ax.plot([r["phi"] for r in pr], [r["rice_mean"] for r in pr], color=INK2, lw=1.0, label="Rice-FG closed form (prime q, no data)")
    ax.plot([r["phi"] for r in pr], [r["rice_sd"] for r in pr], color=INK2, lw=1.0)
    qs = res["m1_qstar"][[q["X"] for q in res["m1_qstar"]].index(X)]["q_star_fg"]
    ax.axvline(qs - 1, color=MUTED, lw=0.8)
    ax.text(qs * 1.05, 0.012, f"φ(q*) ≈ {qs - 1:.0f}: spike = sea", fontsize=7.5, color=INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 30)
    ax.set_xlabel("φ(q)   (squarefree q ≤ 10⁴)")
    ax.set_ylabel("rescaled spike height")
    ax.set_title("(b) collapse error predicted by the variance of primes in progressions")
    ax.legend(loc="upper left", fontsize=7)
    ins = ax.inset_axes([0.62, 0.08, 0.36, 0.33])
    ins.scatter(ph, [r["cos"] for r in rows], s=2, color=VIOLET, lw=0, rasterized=True)
    ins.plot([r["phi"] for r in pr], [r["rice_cos"] for r in pr], color=INK2, lw=0.8)
    ins.set_xscale("log")
    ins.set_ylim(-0.05, 1.05)
    ins.set_title("⟨cos θ⟩: locked phase = real", fontsize=7, loc="left")
    ins.tick_params(labelsize=6)
    ax = axes[1, 0]
    qsx = res["m1_qstar"]
    Xs = np.array([q["X"] for q in qsx], dtype=float)
    ax.plot(Xs, [q["q_star"] for q in qsx], "o", color=BLUE, ms=7, mec=SURF, mew=1.0, ls="none", label="q* measured (running median of 21 primes)")
    ax.plot(Xs, [q["q_star_fg"] for q in qsx], "-", color=BLUE, lw=1.2, label="q* from (q−1)√(X·FG(q)) = ψ (no fit)")
    ax.plot(Xs, [q["sqrt_X_logX"] for q in qsx], color=MUTED, lw=0.9, ls=(0, (4, 2)), label="√(X/log X): order of magnitude only")
    ax.plot(Xs, [q["psi_over_sea"] for q in qsx], color=INK2, lw=0.9, ls=(0, (4, 2)), label="ψ/√ΣΛ²: order of magnitude only")
    for c, col in ((1, AQUA), (10, ORANGE), (20, VIOLET)):
        ax.plot(Xs, [q["Qmax"][c] for q in qsx], color=col, lw=0.9, ls=(0, (1, 1.5)), label=f"Q_max = √(X/2c), c = {c}: arc spacing, different mechanism")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("X")
    ax.set_ylabel("denominator q")
    ax.set_title("(c) where the collapse ends, from X = 10⁴ to 2·10⁶")
    ax.legend(loc="upper left", fontsize=7)
    ax = axes[1, 1]
    wm = res["m1_winmax"]
    NL = np.array([2**j for j in wm["js"]]) * X / res["M"]
    for tag in ("primes", "hybrid", "cramer"):
        ax.plot(NL, wm[tag], marker="o", ms=3.5, lw=1.1, label=MEAS_LABEL[tag], **MEAS_STYLE[tag])
    ax.axhline(math.sqrt(res["L2"]) / psi, color=MUTED, lw=0.8)
    ax.text(NL[0], math.sqrt(res["L2"]) / psi * 1.15, "sea rms √ΣΛ²/ψ", fontsize=7, color=INK2)
    ax.axvspan(120, 16000, color=GRID, alpha=0.6, zorder=0)
    x0, y0 = 300.0, 0.01
    ax.plot([x0, x0 * 10, x0 * 10, x0], [y0, y0, y0 * math.sqrt(10), y0], color=INK2, lw=0.8)
    ax.text(x0 * 11, y0 * 1.5, "slope ½", fontsize=7, color=INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("window width N·L")
    ax.set_ylabel("median over windows of max |S|/ψ")
    ax.set_title("(d) statistical self-affinity of the window maxima (shaded: NL ∈ [120, 16000])")
    ax.legend(loc="upper left", fontsize=7.5)
    _suptitle(fig, "Fig M1. The spikes are a data collapse (heights scale by μ(q)/φ(q), widths stay 1/N) up to a computable cutoff q*")
    fig.tight_layout()
    _save(fig, out, "fig_M1_collapse.png")


# ================================================================ M2
def fig_M2(res, out):
    m2 = res["m2"]
    L, M, X = m2["L"], res["M"], res["X"]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    ax = axes[0, 0]
    for j, col in ((3, VIOLET), (5, BLUE), (6, AQUA), (9, ORANGE)):
        d = m2["rank_mass"][("primes", j)]
        ax.plot(d["ranks"], d["comp"], color=col, lw=1.1, label=f"primes, Nε = {X * 2**j / M:.1f}")
    st = m2["staircase"]
    r_end = st["rank_end"]
    r_start = np.r_[1, r_end[:-1] + 1]
    xs_st = np.ravel(np.column_stack([r_start, r_end]))
    ys_st = np.ravel(np.column_stack([2 * L * r_start * st["mass"], 2 * L * r_end * st["mass"]]))
    sel = xs_st <= 1e5
    ax.plot(xs_st[sel], ys_st[sel], color=INK, lw=0.7, label=f"exact Farey staircase, Nε = {X * 2**6 / M:.1f} (no fit)")
    for tag in ("hybrid", "cramer"):
        d = m2["rank_mass"][(tag, 6)]
        ax.plot(d["ranks"], d["comp"], lw=1.0, label=f"{MEAS_LABEL[tag]}, Nε = 7.6", **{k: v for k, v in MEAS_STYLE[tag].items()})
    ax.axhline(1.0, color=INK2, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.3, 30)
    ax.set_xlabel("rank i of the box (largest first)")
    ax.set_ylabel("2L · i · m₍ᵢ₎   (1 = the Farey count)")
    ax.set_title(f"(a) rank–mass: m₍ᵢ₎ ≈ 1/(2L·i), L = {L:.2f}")
    ax.legend(loc="upper right", fontsize=7)
    ax = axes[0, 1]
    for (j, d), col in zip(sorted(m2["cumulative"].items()), (BLUE, AQUA, ORANGE)):
        a = d["alpha"]
        ok = a >= d["alpha_lo"]
        ax.plot(a[ok], d["f"][ok], color=col, lw=1.4, label=f"f₍<₎(α), ε = 2^{int(round(math.log2(d['eps'])))}")
        ax.plot(a, a + d["offset"], color=col, lw=0.8, ls=(0, (4, 2)))
    cj = m2["cj"]
    ax.plot([cj[s]["alpha"][0] for s in cj], [cj[s]["f"][0] for s in cj], "o", ms=3.5, color=INK, ls="none",
            label="Chhabra–Jensen local slopes, s ∈ [0.8, 1.26]")
    ax.plot([0.3, 1.0], [0.3, 1.0], color=INK2, lw=0.8)
    ax.plot([], [], color=INK2, lw=0.8, ls=(0, (4, 2)), label="no-fit lines α − ln(2L)/ln(1/ε)")
    ax.set_xlim(0.3, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("coarse Hölder exponent α")
    ax.set_ylabel("f(α)")
    ax.set_title("(b) coarse-grained spectrum: the straight line f = α, shifted")
    ax.legend(loc="upper left", fontsize=7)
    ax = axes[0, 2]
    js = m2["Ds_js"]
    Ne = np.array([X * 2 ** (j + 0.5) / M for j in js[:-1]])
    for s, col in ((0.5, AQUA), (1.0, BLUE), (1.5, ORANGE), (2.0, RED)):
        for tag in ("primes", "hybrid", "cramer"):
            st_ = MEAS_STYLE[tag]
            ax.plot(Ne, m2["Ds"][tag][s], color=col, ls=st_["ls"], lw=1.1 if tag == "primes" else 0.9)
        ax.plot([], [], color=col, lw=1.1, label=f"D_{s:g}")
    for tag in ("primes", "hybrid", "cramer"):
        ax.plot([], [], color=INK2, ls=MEAS_STYLE[tag]["ls"], label=MEAS_LABEL[tag])
    ax.axvspan(Ne[0] * 0.8, 15, color=GRID, alpha=0.6, zorder=0)
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("box size Nε")
    ax.set_ylabel("local D_s")
    ax.set_title("(c) generalized dimensions against scale (grey: Nε < 15)")
    ax.legend(loc="lower right", fontsize=6.6, ncol=2)
    ax = axes[1, 0]
    sg = m2["s_grid"]
    for tag in ("primes", "hybrid", "cramer"):
        ax.plot(sg, [m2["Ds_s"][tag][s] for s in sg], marker="o", ms=3.5, lw=1.1, label=MEAS_LABEL[tag], **MEAS_STYLE[tag])
    ax.plot([0.5, 1, 1, 3], [1, 1, 0, 0], color=INK2, lw=0.8, ls=(0, (4, 2)), label="N → ∞ limit (ε = N^(−θ)): not reached")
    ax.set_xlabel("moment order s")
    ax.set_ylabel("D_s at Nε 15 → 31")
    ax.set_title("(d) a broad crossover in s, not a step")
    ax.legend(loc="upper right", fontsize=7)
    ax = axes[1, 1]
    Nm = np.array([X * 2 ** (j + 0.5) / M for j in range(14)])
    for pr, d in m2["minor_D2"].items():
        col = PAIR_COL[pr]
        ax.plot(Nm, d["primes"], color=col, lw=1.2, label=f"{pr}")
        ax.plot(Nm, d["cramer"], color=col, lw=0.8, ls=(0, (1, 1.5)))
        ax.plot(Nm, d["hybrid"], color=col, lw=0.8, ls=(0, (4, 2)))
    ax.plot([], [], color=INK2, lw=1.2, label="primes (solid)")
    ax.plot([], [], color=INK2, lw=0.8, ls=(0, (4, 2)), label="surrogate (dashed)")
    ax.plot([], [], color=INK2, lw=0.8, ls=(0, (1, 1.5)), label="Cramér (dotted)")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("box size Nε")
    ax.set_ylabel("local D₂ of the minor-arc measure")
    ax.set_title("(e) on the minor arcs: residual spikes keep D₂ low at small Q")
    ax.legend(loc="lower right", fontsize=6.8, ncol=2)
    ax = axes[1, 2]
    cols = [VIOLET, MAGENTA, AQUA, ORANGE, BLUE, RED]
    for d, col in zip(m2["xscaling"], cols):
        ax.plot(d["ranks"], d["comp"], color=col, lw=1.1, label=f"X = {d['X']:.0e}, L = {d['L']:.2f}, Nε = {d['Neps']:.1f}")
    ax.axhline(1.0, color=INK2, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.3, 30)
    ax.set_xlabel("rank i")
    ax.set_ylabel("2L(X) · i · m₍ᵢ₎")
    ax.set_title("(f) the same law at every X, each with its own L(X)")
    ax.legend(loc="upper left", fontsize=7)
    _suptitle(fig, "Fig M2. The spectral measure |S|²dα: the Farey count read per box (coarse-grained, ε tied to N; Diophantine, not a cascade)")
    fig.tight_layout()
    _save(fig, out, "fig_M2_spectrum.png")


# ================================================================ M3
def fig_M3(res, out):
    m3, ts = res["m3"], res["m3_t"]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    pairs = [(5, 5), (20, 10), (50, 10), (150, 20), (223, 10)]
    for ax, pr in zip(axes.ravel()[:5], pairs):
        r = m3[pr]
        tmax = r["primes"]["sup2"] / r["primes"]["mean_m"]
        for tag in ("primes", "hybrid", "cramer"):
            y = np.where(r[tag]["ccdf"] > 0, r[tag]["ccdf"], np.nan)
            ax.plot(ts, y, lw=1.3 if tag == "primes" else 1.0, label=MEAS_LABEL[tag], **MEAS_STYLE[tag])
        sel = ts >= 4
        ax.plot(ts[sel], np.where(r["model"][sel] > 0, r["model"][sel], np.nan), color=INK, lw=0.9,
                label="spike sum, semi-analytic (no fit)")
        sa = ts <= tmax
        ax.plot(ts[sa], r["asymptote"][sa], color=INK2, lw=0.8, ls=(0, (4, 2)), label="asymptote ψ²/(2N|m|·t·mean)")
        ax.plot(ts, np.exp(-ts), color=MUTED, lw=0.6, ls=(0, (2, 2)), label="e^(−t) (exponential)")
        ax.axvline(tmax, color=PAIR_COL[pr], lw=0.9)
        ax.text(tmax * 0.95, 0.3, f"max: {r['sup_pred']['active']}", fontsize=6.6, ha="right", color=INK2, rotation=90, va="top")
        lo, hi = r["pareto_range"]
        if hi >= 10 * lo:
            ax.axvspan(lo, hi, color=GRID, alpha=0.6, zorder=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-8, 1.5)
        ax.set_xlim(1, ts[-1])
        ax.set_xlabel("t = |S|²/mean over the minor arcs")
        ax.set_ylabel("P(> t)")
        ax.set_title(f"({'abcde'[pairs.index(pr)]}) (Q, c) = {pr}" + ("  (shaded: Pareto range)" if hi >= 10 * lo else ""))
        if pr == (5, 5):
            ax.legend(loc="lower left", fontsize=6.8)
    ax = axes[1, 2]
    tp = res.get("m3_pairs_order") or [(5, 5), (10, 10), (20, 10), (50, 10), (100, 20), (150, 20), (223, 10)]
    x = np.arange(len(tp))
    for tag, mk in (("primes", "o"), ("hybrid", "o"), ("cramer", "x")):
        y = [m3[p][tag]["hill"][1000] for p in tp]
        if tag == "hybrid":
            ax.plot(x, y, mk, mfc="none", mec=ORANGE, ms=7, ls="none", label="Hill index (k = 1000), surrogate")
        elif tag == "primes":
            ax.plot(x, y, mk, color=BLUE, ms=7, ls="none", label="Hill index (k = 1000), primes")
        else:
            ax.plot(x, y, mk, color=MUTED, ms=7, mew=1.5, ls="none", label="Hill index (k = 1000), Cramér")
    ax.axhline(1.0, color=INK, lw=0.9)
    ax.text(len(tp) - 0.6, 1.05, "index 1 in |S|² = 2 in |S|: wild/pre-wild boundary for |S|", fontsize=6.6, ha="right")
    ax.axhline(2.0, color=INK2, lw=0.8, ls=(0, (4, 2)))
    ax.text(len(tp) - 0.6, 2.08, "above 2: variance of |S|² finite", fontsize=6.6, ha="right")
    ax.set_xticks(x, [f"({p[0]},{p[1]})" for p in tp], fontsize=7.5)
    ax.set_ylabel("Hill tail index of |S|² (cell maxima)")
    ax.set_ylim(0, 8)
    ax2 = ax.twinx()
    ax2.plot(x, [m3[p]["primes"]["kappa"] for p in tp], "s", color=VIOLET, ms=6, ls="none", label="κ measured")
    km = [m3[p]["E4_model"] / (m3[p]["primes"]["P_m"] ** 2 / m3[p]["meas_m"]) for p in tp]
    ax2.plot(x, km, "s", mfc="none", mec=VIOLET, ms=9, ls="none", label="κ from the fourth-moment model")
    ax2.set_yscale("log")
    ax2.set_ylabel("minor-arc kurtosis κ (log)")
    ax2.grid(False)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=6.8)
    ax.set_title("(f) wild at small Q, neither wild nor exponential at large Q")
    _suptitle(fig, "Fig M3. Mild versus wild on the minor arcs, by (Q, c): the tail of |S|² is the Farey count again, truncated by the arcs")
    fig.tight_layout()
    _save(fig, out, "fig_M3_tails.png")


# ================================================================ M4
def _band(ens, key):
    a = np.array([e[key] for e in ens])
    return a.mean(axis=0), a.std(axis=0, ddof=1)


def fig_M4(res, out):
    from fractals import H_bridge, H_triangular

    X, M = res["X"], res["M"]
    xl = _xl(res)
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    for ax, pr in zip(axes[0], ((150, 20), (20, 10))):
        w = res["walks"][("primes",) + pr]
        n = w["n"]
        H = w["full"][2.0]["H"]
        good = n / np.array(LAGS[:-1]) >= 500
        mean, sd = _band(res["ensemble"][pr], "raw")
        ax.fill_between(xl, mean - 2 * sd, mean + 2 * sd, color=GRID, alpha=0.9, lw=0, label="Cramér ensemble, 16 seeds (±2 sd)")
        ax.plot(xl, H, color=BLUE, lw=1.3)
        ax.plot(xl[good], H[good], "o", color=BLUE, ms=4, ls="none", label="primes H₂ (solid: ≥ 500 independent increments)")
        ax.plot(xl[~good], H[~good], "o", mfc="none", mec=BLUE, ms=4, ls="none")
        ax.plot(xl, w["full"]["drift"]["H"], color=BLUE, lw=0.6, ls=(0, (1, 1)), label="primes, drift removed")
        for R, col in ((300, AQUA), (1000, ORANGE), (3000, VIOLET)):
            ax.plot(xl, res["walks"][(f"hybrid{R}",) + pr]["full"][2.0]["H"], color=col, lw=0.9, ls=(0, (4, 2)),
                    label=f"surrogate, spikes q ≤ {R}")
        ht = [H_triangular(x) for x in xl[xl < 200]]
        ax.plot(xl[xl < 200], ht, color=INK, lw=0.9, label="triangular envelope (no fit)")
        hb = [H_bridge(L, n) for L in LAGS[:-1] if 2 * L < n]
        ax.plot(xl[: len(hb)], hb, color=MUTED, lw=0.7, ls=(0, (2, 2)), label="pinned-walk (bridge) reference")
        sp = res["spectral"][("primes",) + pr]
        Hc = 0.5 * np.diff(np.log(sp["S2_spectral"])) / math.log(2)
        ax.plot(xl, Hc, color=INK2, lw=0.6, ls=(0, (1, 2)), label="zero-filled cyclic walk (spectral identity)")
        ax.axhline(0.5, color=INK, lw=0.8)
        ax.set_xscale("log")
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel("lag N·ℓ/M (units of 1/N)")
        ax.set_ylabel("local H₂")
        ax.set_title(f"({'ab'[pr == (20, 10)]}) minor-arc walk at {pr}: smooth below 1/N, then ≈ ½")
        ax.legend(loc="lower left", fontsize=6.4)
    ax = axes[1, 0]
    groups = [("primes", "full"), ("primes", "noah"), ("primes", "joseph"), ("cramer", "full"), ("hybrid1000", "full")]
    names = ["primes", "Noah-only\n(phases randomised)", "Joseph-only\n(amplitudes flattened)", "Cramér", "surrogate\nH₁₀₀₀"]
    sel = (xl >= 4) & (xl <= 1000)
    width = 0.13
    k = 0
    for pr, hatch in (((150, 20), None), ((20, 10), "//")):
        for p, col in ((0.5, VIOLET), (1.0, AQUA), (2.0, BLUE)):
            vals = []
            for tag, var in groups:
                st = res["walks"][(tag,) + pr]
                d = st.get(var, {}).get(p) if var != "full" else st["full"].get(p)
                vals.append(float(np.mean(d["H"][sel])) if d is not None else np.nan)
            ax.bar(np.arange(len(groups)) + (k - 2.5) * width, vals, width, color=col, hatch=hatch, edgecolor=SURF, lw=0.6,
                   label=f"H_{p:g}, {pr}")
            k += 1
    ax.axhline(0.5, color=INK, lw=0.8)
    ax.set_xticks(np.arange(len(groups)), names, fontsize=7.5)
    ax.set_ylabel("H_p, mean over Nℓ ∈ [4, 1000]")
    ax.set_ylim(0, 1.0)
    ax.set_title("(c) heavy tails (Noah) show in H_½; no long memory (Joseph)")
    ax.legend(loc="upper right", fontsize=6.6, ncol=2)
    ax = axes[1, 1]
    sp = res["spectral"][("primes", 150, 20)]
    e = sp["d_edges"]
    ctr = np.sqrt(e[:-1] * e[1:])
    ax.plot(ctr, sp["P_binned"], color=BLUE, lw=1.2, label="primes")
    for tag in ("hybrid1000", "cramer"):
        if (tag, 150, 20) in res["spectral"]:
            st = MEAS_STYLE["hybrid" if tag.startswith("hybrid") else "cramer"]
            ax.plot(ctr, res["spectral"][(tag, 150, 20)]["P_binned"], lw=1.0, label=MEAS_LABEL["hybrid" if tag.startswith("hybrid") else "cramer"], **st)
    for xx in (1, 100, 3906):
        L = xx * M / X
        ax.axvline(M / L, color=MUTED, lw=0.7)
        ax.text(M / L * 1.1, ax.get_ylim()[1] if False else max(sp["P_binned"]) * 0.5, f"M/ℓ at Nℓ/M = {xx}", fontsize=6.6, rotation=90, color=INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("|d| = |N′ − N| (neighbouring targets)")
    ax.set_ylabel("binned |I_m(N + d)|²")
    ax.set_title("(d) H is a Fejér average of the minor-arc counts at nearby targets, (150, 20)")
    ax.legend(loc="lower left", fontsize=7)
    ins = ax.inset_axes([0.62, 0.62, 0.35, 0.3])
    ins.plot(xl, np.abs(sp["S2_spectral"][:-1] / sp["S2_direct"][:-1] - 1) + 1e-18, color=INK, lw=0.8)
    ins.set_xscale("log")
    ins.set_yscale("log")
    ins.set_title("|spectral/direct − 1|", fontsize=6.5, loc="left")
    ins.tick_params(labelsize=6)
    _suptitle(fig, "Fig M4. Roughness of the minor-arc walk: band limit (generic), residual spikes (H₂ ≈ 0.52), heavy tails (Noah), no long memory")
    fig.tight_layout()
    _save(fig, out, "fig_M4_walks.png")


# ================================================================ M5
def fig_M5(res, out):
    m5 = res["m5"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    ax = axes[0, 0]
    img = m5["render"]
    x0, x1, y0, y1 = m5["render_extent"]
    maxit = img.max()
    shade = np.where(img >= maxit, 0.0, np.log1p(img) / math.log1p(maxit))
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("m", [INK, "#2a3f5f", BLUE, "#cde2fb", SURF])
    ax.imshow(shade, extent=(x0, x1, y0, y1), cmap=cm, origin="upper", interpolation="nearest")
    for b in m5["bulbs"]:
        if b["q"] <= 8:
            from fractals import cardioid_root
            for sgn in (1, -1):
                r = cardioid_root(b["p"], b["q"])
                r = complex(r.real, sgn * r.imag)
                ax.plot(r.real, r.imag, "o", ms=3, color=ORANGE)
                ax.annotate(f"{b['p'] if sgn > 0 else b['q'] - b['p']}/{b['q']}", (r.real, r.imag), xytext=(4, 4),
                            textcoords="offset points", fontsize=6.5, color=ORANGE)
    ax.set_xlabel("Re c")
    ax.set_ylabel("Im c")
    ax.set_title("(a) the Mandelbrot set; bulb roots on the main cardioid labelled by p/q (q ≤ 8)")
    ax.grid(False)
    ax = axes[0, 1]
    bl = m5["bulbs"]
    sc = ax.scatter([b["theta"] for b in bl], [b["diam_half"] for b in bl], c=[b["q"] for b in bl], cmap="viridis", s=14,
                    label="diam/2 ÷ sin(πp/q)/q²")
    ax.scatter([b["theta"] for b in bl], [b["centre_root"] for b in bl], c=[b["q"] for b in bl], cmap="viridis", s=14,
               marker="o", facecolors="none", label="|centre − root| ÷ sin(πp/q)/q²")
    ax.axhline(1.0, color=INK, lw=0.8)
    fig.colorbar(sc, ax=ax, label="q")
    ax.set_xlabel("θ = p/q")
    ax.set_ylabel("measured size ÷ leading-order estimate")
    ax.set_title("(b) Mandelbrot's empirical N² rule: prefactor 0.98–1.16, depending on p/q")
    ax.legend(loc="upper left", fontsize=7.5)
    ax = axes[1, 0]
    sz = m5["sizes"]
    q = np.array([s["q"] for s in sz])
    ax.plot(q, [s["bulb"] for s in sz], "o", color=BLUE, ms=5, ls="none", label="largest bulb with denominator q (diam/2)")
    sq = np.array([s["squarefree"] for s in sz])
    ax.plot(q[sq], [s["spike_mass"] for s in sz if s["squarefree"]], "s", color=ORANGE, ms=5, ls="none",
            label="spike mass 1/(L φ(q)²), squarefree q")
    ax.plot(q, [s["ford"] for s in sz], color=INK2, lw=0.9, label="Ford circle radius 1/(2q²) (no dynamics)")
    for qq in q[~sq]:
        ax.axvline(qq, color=GRID, lw=2.0, zorder=0)
    ax.plot([], [], color=GRID, lw=4, label="non-squarefree q: bulbs, but no spike")
    ax.plot([3, 30], [0.2, 0.2 / 100], color=MUTED, lw=0.8, ls=(0, (4, 2)), label="slope −2 guide")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("denominator q")
    ax.set_ylabel("size")
    ax.set_title("(c) three q⁻² families: bulbs, spike masses (not heights), Ford circles")
    ax.legend(loc="lower left", fontsize=7)
    ax = axes[1, 1]
    from fractals import farey
    pts = [(a, b) for a, b in farey(21) if 0 < a / b <= 0.5 and b >= 2]
    sizes = {(b["p"], b["q"]): b["diam_half_abs"] for b in bl}
    L = res["m2"]["L"]
    phi = None
    for a, b in pts:
        th = a / b
        bs = sizes.get((a, b), None)
        if bs:
            ax.vlines(th, 0, bs / sizes[(1, 2)], color=BLUE, lw=1.2)
    mus = {s["q"]: s for s in m5["sizes"]}
    for a, b in pts:
        th = a / b
        s = mus.get(b)
        if s and s["squarefree"]:
            ax.vlines(th, 0, -s["spike_mass"] / (1.0 / (L * 1.0)), color=ORANGE, lw=1.2)
        else:
            ax.plot(th, 0, "x", color=RED, ms=4)
    ax.axhline(0, color=INK, lw=0.8)
    med = m5["mediant"][21]
    ax.text(0.02, 0.92, f"mediant rule (F₂₁): bulbs {med['bulb_mediant']}/{med['pairs']}, spikes {med['spike_mediant']}/{med['pairs']}; "
            f"{med['mediant_nonsquarefree']} mediants non-squarefree", transform=ax.transAxes, fontsize=7.5)
    ax.plot([], [], color=BLUE, lw=1.2, label="bulb size (÷ the 1/2 bulb)")
    ax.plot([], [], color=ORANGE, lw=1.2, label="spike mass (÷ the q = 1, 2 mass), downwards")
    ax.plot([], [], "x", color=RED, ls="none", label="non-squarefree q: no spike")
    ax.set_xlim(0, 0.52)
    ax.set_xlabel("θ = a/q (denominators ≤ 21)")
    ax.set_ylabel("relative size")
    ax.set_title("(d) same index set and circle order; the mediant rule holds for bulbs only")
    ax.legend(loc="lower right", fontsize=7)
    _suptitle(fig, "Fig M5. The Mandelbrot set's bulbs and the primes' Farey forest: shared Diophantine bookkeeping, not shared dynamics")
    fig.tight_layout()
    _save(fig, out, "fig_M5_mandelbrot.png")


# ================================================================ M6 (A1, A2, A3)
def fig_M6(res, out):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))
    ax = axes[0]
    oc = res["a1"]["octaves"]
    x = np.arange(len(oc))
    ax.bar(x, [o["frac"] for o in oc], 0.6, color=BLUE, label="measured energy in |N(α − a/q)| ≤ ½ (share of ΣΛ²)")
    ax.plot(x, [o["spike_model"] for o in oc], "o", color=ORANGE, ms=6, ls="none", label="spike model 0.7737·(ψ²/N)Σμ²/φ")
    ax.plot(x, [o["spike_sea_model"] for o in oc], "s", mfc="none", mec=INK, ms=7, ls="none", label="spike + sea model")
    ax.axhline(res["a1"]["equal_line"], color=INK2, lw=0.9, ls=(0, (4, 2)), label="equal energy per octave")
    ax.set_xticks(x, [f"[{2**o['j']},{2**(o['j'] + 1)})" for o in oc], rotation=45, fontsize=7)
    ax.set_xlabel("octave of denominators q")
    ax.set_ylabel("share of ΣΛ²")
    ax.set_title("(a) equal energy per octave of q: the L² wall in fractal language")
    ax.legend(loc="upper left", fontsize=7)
    ax = axes[1]
    a2 = res["a2"]
    V = np.array([res["X"] ** s for s in a2["sigma"]])
    ax.plot(V, a2["cells"], "o-", color=BLUE, ms=5, label="cells (width ≈ 1/N) with max |S| ≥ N^σ")
    ax.plot(V, a2["spikes"], "s-", color=ORANGE, ms=5, label="Farey count of spikes taller than N^σ")
    ax.plot(V, a2["sieve_bound"], color=INK2, lw=1.0, ls=(0, (4, 2)), label="dual large sieve bound (rigorous)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("level V = N^σ")
    ax.set_ylabel("count")
    ax.set_title("(b) large values: count exponent 2(1 − σ) (Jarník–Besicovitch shadow)")
    ax.legend(loc="lower left", fontsize=7)
    ins = ax.inset_axes([0.6, 0.58, 0.37, 0.36])
    sig = np.array(a2["sigma"])
    cells = np.array(a2["cells"], dtype=float)
    slope = -np.diff(np.log(cells)) / (np.diff(sig) * math.log(res["X"]))
    ins.plot(0.5 * (sig[1:] + sig[:-1]), slope, "o-", color=BLUE, ms=3)
    ins.axhline(2.0, color=INK2, lw=0.8)
    ins.set_title("local slope vs 2", fontsize=6.5, loc="left")
    ins.tick_params(labelsize=6)
    ax = axes[2]
    un = res["uncertainty"]
    Ls = sorted(un["sharp"])
    ax.plot(Ls, [un["sharp"][L]["product"] for L in Ls], "o-", color=RED, ms=6, label="sharp cutoff n ≤ N: σ_x·σ_ξ over |Nβ| ≤ L")
    Lb = sorted(un["bump"])
    ax.plot(Lb, [un["bump"][L]["product"] for L in Lb], "s-", color=AQUA, ms=6, label="smooth bump η: σ_x·σ_ξ")
    ax.axhline(1 / (4 * math.pi), color=INK, lw=1.0, ls=(0, (4, 2)), label="Heisenberg–Gabor bound 1/(4π) (ħ → 1/2π)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("frequency window L (units of 1/N)")
    ax.set_ylabel("σ_n · σ_α")
    ax.set_title("(c) the uncertainty principle: sharp sums have infinite spread (leakage)")
    ax.legend(loc="upper left", fontsize=7)
    _suptitle(fig, "Fig M6. Two more readings of the Farey count, and the uncertainty principle behind the arc widths")
    fig.tight_layout()
    _save(fig, out, "fig_M6_octaves_largevalues_uncertainty.png")


ALL = [fig_M1, fig_M2, fig_M3, fig_M4, fig_M5, fig_M6]


def make_all(res, out: Path) -> None:
    out = Path(out)
    out.mkdir(exist_ok=True)
    for f in ALL:
        f(res, out)
