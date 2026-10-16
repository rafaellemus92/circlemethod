"""The fractal lens (M1-M5, A1-A3): Mandelbrot's measurement programme on the circle-method objects at X = 10^6.

    python run_fractals.py            # full run (~12 min on 4 cores, peak RSS ~1.6 GB); needs run_all.py's and run_escapes.py's caches
    python run_fractals.py --report   # redraw figures/fig_M*.png and fractals.md from cache/fractals.pkl
"""
from __future__ import annotations

import argparse
import math
import pickle
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fractals as fr  # noqa: E402
from arcs import MajorArcs  # noqa: E402
from expsum import S_grid, S_rational, e_minus  # noqa: E402
from sieve import cramer_weights, mobius_phi, primes_upto, von_mangoldt  # noqa: E402

CACHE = HERE / "cache"
FIGS = HERE / "figures"
X = 10**6
M = 2**23
WALK_PAIRS = [(150, 20), (20, 10)]
CCDF_PAIRS = [(5, 5), (20, 10), (50, 10), (150, 20), (223, 10)]
TAIL_PAIRS = [(5, 5), (10, 10), (20, 10), (50, 10), (100, 20), (150, 20), (223, 10)]
XS = [(10**4, 2**17), (3 * 10**4, 2**18), (10**5, 2**20), (3 * 10**5, 2**22), (10**6, 2**23), (2 * 10**6, 2**24)]
RS = [300, 1000, 3000]
LAGS = [2**j for j in range(22)]
N_SEEDS = 16
T0 = time.time()


def log(msg: str) -> None:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[{time.time() - T0:7.1f}s, peak {peak:6.0f} MB] {msg}", flush=True)


def load_S(tag: str, Xs: int = X, Ms: int = M) -> np.ndarray:
    return np.load(CACHE / f"S_{tag}_N{Xs}_M{Ms}.npz")["S"]


def certify(name: str, ok: bool, detail: str, res: dict) -> None:
    res["certificates"].append((name, bool(ok), detail))
    if not ok:
        raise AssertionError(f"certificate failed: {name}: {detail}")


# ---------------------------------------------------------------- stage 1: spike collapse (no grid needed)
def stage_collapse(res, lam, psi, mu, phi):
    n = np.flatnonzero(lam)
    w = lam[n]
    qs = [q for q in range(1, 10001) if mu[q] != 0]
    rows = fr.spike_collapse(n, w, psi, qs, mu, phi)
    pset = set(int(p) for p in primes_upto(10001))
    for r in rows:
        E2, EZ2 = r["E2"], r["EZ2"]
        b = fr.beckmann_moments(0.5 * (E2 + EZ2), 0.5 * (E2 - EZ2))
        r.update(beck_mean=b["mean"], beck_sd=b["sd"], beck_cos=b["cos"])
        q = r["q"]
        if q in pset and fr.fg_factor(q, [q]) > 0:            # the FG asymptotic is positive only for q >= 13
            s2 = (q - 1) ** 2 * X * fr.fg_factor(q, [q]) / psi**2
            rc = fr.rice_moments(s2)
            r.update(fg=fr.fg_factor(q, [q]), rice_mean=rc["mean"], rice_sd=rc["sd"], rice_cos=rc["cos"])
    res["m1_rows"] = rows
    # certificates 3 and 4
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(20):
        q = int(rng.integers(2, 10001))
        res_q = np.bincount(n % q, weights=w, minlength=q)
        a = int(rng.integers(0, q))
        worst = max(worst, abs(q * np.fft.ifft(res_q)[a] - S_rational(a, q, res_q)))
    certify("all-a FFT = S_rational (20 random a/q, q <= 10^4)", worst <= 1e-12 * psi, f"max |diff| = {worst / psi:.1e} psi", res)
    for q in (101, 210, 1001):
        res_q = np.bincount(n % q, weights=w, minlength=q)
        r = np.arange(q)
        red = np.gcd(r, q) == 1
        delta = np.where(red, res_q - res_q[red].sum() / phi[q], 0.0)
        lhs = np.mean([abs(np.sum(np.exp(2j * np.pi * a * r / q) * delta)) ** 2 for a in r[red]])
        from singular_series import ramanujan_divisor
        cq = np.array([ramanujan_divisor(q, int(d)) for d in range(q)], dtype=float)
        dif = (r[:, None] - r[None, :]) % q
        rhs = float(delta @ (cq[dif] @ delta)) / phi[q]
        certify(f"second-moment identity q = {q}", abs(lhs - rhs) <= 1e-10 * abs(rhs), f"rel {abs(lhs - rhs) / abs(rhs):.1e}", res)
    log(f"M1(b): {len(rows)} squarefree q <= 10^4")


def stage_qstar(res):
    out = []
    for Xs, _ in XS:
        lam = von_mangoldt(Xs)
        n = np.flatnonzero(lam)
        w = lam[n]
        psi = float(lam.sum())
        hi = int(min(10**4, 4 * math.sqrt(Xs)))
        primes = [int(p) for p in primes_upto(hi) if p >= 7]
        qm, K = fr.q_star_measured(n, w, psi, Xs, primes)
        out.append({"X": Xs, "psi": psi, "q_star": qm, "q_star_fg": fr.q_star_fg(Xs, psi), "K": K,
                    "sqrt_X_logX": math.sqrt(Xs / math.log(Xs)), "psi_over_sea": psi / math.sqrt(float(np.sum(lam**2))),
                    "Qmax": {c: math.isqrt((Xs - 1) // (2 * c)) for c in (1, 10, 20)}})
    res["m1_qstar"] = out
    log("M1(c): q*(X)")


# ---------------------------------------------------------------- stage 2: primes grid
def zoom_collapse(res, absS, psi, mu, phi):
    classes = [(3, 30), (101, 211), (499, 1009), (2003, 4001)]
    bins = np.arange(-10.0, 10.0001, 0.1)
    out = []
    for lo, hi in classes:
        qs = [q for q in range(lo, hi + 1) if mu[q] != 0]                     # every squarefree q in the class
        xs, ys = [], []
        for q in qs:
            a_s = [a for a in range(1, q) if math.gcd(a, q) == 1]
            a_s = [a_s[int(i)] for i in np.linspace(0, len(a_s) - 1, min(50, len(a_s)))]
            for a in a_s:
                k0 = a * M / q
                k = np.arange(int(math.floor(k0 - 10 * M / X)), int(math.ceil(k0 + 10 * M / X)) + 1)
                x = X * (k / M - a / q)
                sel = np.abs(x) <= 10
                xs.append(x[sel])
                ys.append(absS[k[sel] % M] * phi[q] / psi)
        x, y = np.concatenate(xs), np.concatenate(ys)
        idx = np.digitize(x, bins) - 1
        med = np.full(len(bins) - 1, np.nan)
        p10 = med.copy()
        p90 = med.copy()
        for b in range(len(bins) - 1):
            v = y[idx == b]
            if len(v):
                med[b], p10[b], p90[b] = np.median(v), np.quantile(v, 0.1), np.quantile(v, 0.9)
        out.append({"class": (lo, hi), "n_q": len(qs), "centres": 0.5 * (bins[1:] + bins[:-1]), "median": med,
                    "p10": p10, "p90": p90})
    res["m1_zoom"] = out


def walk_stats(F: np.ndarray, ps=(0.5, 1.0, 2.0), drift: bool = True) -> dict:
    o = {}
    for p in ps:
        sf = fr.structure_function(F, LAGS, p)
        o[p] = {"S": sf, "H": fr.local_hurst(LAGS, sf, p)}
    if drift:
        sd = fr.structure_function(F, LAGS, 2.0, drift=True)
        o["drift"] = {"S": sd, "H": fr.local_hurst(LAGS, sd, 2.0)}
    return o


def stage_primes_grid(res, psi, mu, phi, keep):
    S = load_S("primes")
    a2 = S.real**2 + S.imag**2
    absS = np.sqrt(a2)
    zoom_collapse(res, absS, psi, mu, phi)
    res["m1_winmax"] = {"primes": fr.window_max_median(absS, range(3, 21), psi)}
    del absS
    esc = pickle.load(open(CACHE / "escapes.pkl", "rb"))
    res["walks"] = {}
    res["order"] = {}
    res["spectral"] = {}
    for (Q, c) in WALK_PAIRS:
        A = MajorArcs(X, M, Q, c, mu, phi)
        minor = ~A.mask()
        k = np.flatnonzero(minor).astype(np.int32)
        ph = e_minus(X, k, M)
        keep[(Q, c)] = (k, ph)
        w = S[k] ** 2 * ph / M
        F = np.cumsum(w)
        Im = esc["e0"][(Q, c)]["I_m"]
        certify(f"walk endpoint = escapes I_m at {(Q, c)}", abs(F[-1] - Im) <= 1e-9 * abs(Im) and len(k) == esc["e0"][(Q, c)]["n_minor"],
                f"{F[-1].real:.2f} vs {Im.real:.2f}, n = {len(k)}", res)
        st = {"full": walk_stats(F), "n": len(k)}
        st["noah"] = walk_stats(np.cumsum(fr.noah_only(w, 11)), drift=False)
        st["noah_seeds"] = [fr.local_hurst(LAGS, fr.structure_function(np.cumsum(fr.noah_only(w, sd)), LAGS, 2.0))
                            for sd in range(11, 27)]              # the random phases' seed matters when few |w| dominate
        st["joseph"] = walk_stats(np.cumsum(fr.joseph_only(w)), drift=False)
        res["walks"][("primes", Q, c)] = st
        wM = S[A.idx] ** 2 * e_minus(X, A.idx, M) / M
        res["order"][(Q, c)] = {"major": fr.order_parameter_entropy(wM), "minor": fr.order_parameter_entropy(w)}
        del F, wM
        # spectral identity of the zero-filled cyclic walk
        x = np.zeros(M, dtype=complex)
        x[k] = w
        del w
        Xh = np.fft.fft(x)
        S2s = fr.spectral_S2(Xh, LAGS)
        Ff = np.cumsum(x)
        del x
        tot = Ff[-1]
        S2d = np.empty(len(LAGS))
        for i, L in enumerate(LAGS):
            d = np.empty(M, dtype=complex)
            d[: M - L] = Ff[L:] - Ff[: M - L]
            d[M - L :] = Ff[: L] + tot - Ff[M - L :]
            S2d[i] = float(np.mean(d.real**2 + d.imag**2))
            del d
        del Ff
        rel = float(np.max(np.abs(S2s / S2d - 1)))
        certify(f"spectral identity, zero-filled cyclic walk {(Q, c)}", rel <= 1e-9, f"max rel {rel:.1e}", res)
        P = Xh.real**2 + Xh.imag**2
        del Xh
        dd = np.minimum(np.arange(M), M - np.arange(M))
        edges = np.unique(np.geomspace(1, M // 2, 41).astype(np.int64))
        binned = [float(P[(dd >= lo) & (dd < hi)].mean()) for lo, hi in zip(edges[:-1], edges[1:])]
        res["spectral"][("primes", Q, c)] = {"S2_spectral": S2s, "S2_direct": S2d, "d_edges": edges, "P_binned": binned}
        del P, dd
        log(f"M4 primes walk {(Q, c)}")
    del S
    return a2


def stage_other_walks(res, tag, S, keep, pairs=WALK_PAIRS, variants=True, ps=(0.5, 1.0, 2.0), spectral=False):
    for (Q, c) in pairs:
        k, ph = keep[(Q, c)]
        w = S[k] ** 2 * ph / M
        st = {"full": walk_stats(np.cumsum(w), ps=ps), "n": len(k)}
        if variants:
            st["noah"] = walk_stats(np.cumsum(fr.noah_only(w, 11)), drift=False)
            st["joseph"] = walk_stats(np.cumsum(fr.joseph_only(w)), drift=False)
        if spectral:
            x = np.zeros(M, dtype=complex)
            x[k] = w
            P = np.abs(np.fft.fft(x)) ** 2
            del x
            dd = np.minimum(np.arange(M), M - np.arange(M))
            edges = res["spectral"][("primes", Q, c)]["d_edges"]
            res["spectral"][(tag, Q, c)] = {"P_binned": [float(P[(dd >= lo) & (dd < hi)].mean())
                                                       for lo, hi in zip(edges[:-1], edges[1:])]}
            del P, dd
        res["walks"][(tag, Q, c)] = st
        del w


# ---------------------------------------------------------------- stage 6: measures on the grids
def spectrum_stats(res, a2p, a2c, a2h, psi, L2, mu, phi, sea_level):
    out = {}
    measures = {}
    for tag, a2 in (("primes", a2p), ("hybrid", a2h), ("cramer", a2c)):
        kept = fr.remove_atoms(a2, X)
        if tag == "primes":
            out["atom_share"] = 1.0 - kept.sum() / a2.sum()
            out["L"] = fr.kept_L(kept, X, psi)
        measures[tag] = kept / kept.sum()
        del kept
    L = out["L"]
    sea_level = (L2 - res["hybrid_sumsq"][1000]) / (L2 * (1 - out["atom_share"]))   # rho: sea mass density of mu
    res["sea_level"] = sea_level
    certify("atom share and L", abs(out["atom_share"] - 0.1553) <= 5e-4 and abs(L - 10.823) <= 5e-3,
            f"atoms {out['atom_share']:.4f}, L {L:.3f}", res)
    # (a) compensated rank-mass
    out["rank_mass"] = {}
    for tag in ("primes", "hybrid", "cramer"):
        for j in ((3, 5, 6, 9) if tag == "primes" else (6,)):
            m = fr.rank_mass(measures[tag], j, 10**5)        # every rank up to 10^5 (16384 boxes at j = 9)
            rj = np.arange(1, len(m) + 1)
            out["rank_mass"][(tag, j)] = {"ranks": rj, "comp": 2 * L * rj * m}
    eps6 = 2**6 / M
    kap = fr.box_capture(X * eps6)
    out["staircase"] = fr.farey_staircase(L, mu, phi, kappa=kap, sea=sea_level * eps6, qmax=20000)
    out["staircase"]["kappa"] = kap
    eps9 = 2**9 / M                                     # exact staircase count at eps = 2^-14 (Neps = 61)
    st9 = fr.farey_staircase(L, mu, phi, kappa=fr.box_capture(X * eps9), sea=sea_level * eps9, qmax=20000)
    b9 = fr.box_masses(measures["primes"], 9)
    out["stair_counts"] = {A: (int(np.sum(b9 >= eps9**A)), int(st9["phi"][st9["mass"] >= eps9**A].sum()))
                           for A in (0.5, 0.6, 0.7, 0.8, 0.9)}
    # (b) cumulative spectrum, CJ
    alphas = np.round(np.arange(0.30, 1.0001, 0.01), 2)
    out["cumulative"] = {}
    for j in (5, 7, 9):
        eps = 2**j / M
        nn, f = fr.cumulative_spectrum(fr.box_masses(measures["primes"], j), eps, alphas)
        out["cumulative"][j] = {"eps": eps, "alpha": alphas, "N": nn, "f": f, "offset": -math.log(2 * L) / math.log(1 / eps),
                                "alpha_lo": math.log(120 * L) / math.log(1 / eps)}   # model N_< = 60
    ss = np.round(np.arange(0.80, 1.2601, 0.02), 2)
    out["cj"] = fr.chhabra_jensen(measures["primes"], [5, 7], list(ss))
    out["offsets"] = {j: -math.log(2 * L) / math.log(1 / (2**j / M)) for j in (3, 5, 7, 9, 11, 13)}
    # (c, d) generalized dimensions
    s_list = [0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0]
    s_grid = sorted(set([round(x, 2) for x in np.arange(0.5, 1.5001, 0.05)] + [2.0, 3.0]))
    js = list(range(0, 23))
    out["Ds"] = {tag: fr.partition_dims(measures[tag], js, s_list) for tag in measures}
    out["Ds_s"] = {tag: {s: float(fr.partition_dims(measures[tag], [7, 8], [s])[s][0]) for s in s_grid} for tag in measures}
    out["Ds_js"] = js
    out["s_grid"] = s_grid
    res["m2"] = out
    return measures


def minor_arc_stats(res, a2p, a2c, a2h, measures, psi, L2, mu, phi, sea_level):
    m2e = {}
    m3 = {}
    esc = pickle.load(open(CACHE / "escapes.pkl", "rb"))
    ts = 2.0 ** (np.arange(61) / 4)
    levels = fr.SincLevels()
    full_boxes = fr.box_masses(measures["primes"], 5)
    alpha_full = np.full(len(full_boxes), np.inf)
    pos = full_boxes > 0
    alpha_full[pos] = np.log(full_boxes[pos]) / math.log(2**5 / M)
    for (Q, c) in TAIL_PAIRS:
        A = MajorArcs(X, M, Q, c, mu, phi)
        minor = ~A.mask()
        meas_m = minor.sum() / M
        row = {"meas_m": meas_m}
        for tag, a2 in (("primes", a2p), ("hybrid", a2h), ("cramer", a2c)):
            v = a2[minor]
            mean_m = float(v.mean())
            Pm = float(v.sum()) / M
            E4m = float(np.sum(v * v)) / M
            d = {"mean_m": mean_m, "P_m": Pm, "E4m": E4m, "kappa": E4m / (Pm**2 / meas_m), "sup2": float(v.max()),
                 "ccdf": fr.ccdf(v / mean_m, ts)}
            cm = fr.cell_maxima(a2 / mean_m, minor)
            d["hill"] = {kk: fr.hill(cm, kk) for kk in (30, 100, 300, 1000, 3000, 10000)}
            ce = (a2 * minor).reshape(-1, 8)[minor.reshape(-1, 8).all(axis=1)].sum(axis=1)
            d["concentration"] = fr.concentration(ce, [1e-4, 1e-3, 1e-2, 1e-1])
            row[tag] = d
            del v, cm, ce
        p = row["primes"]
        e1 = esc["e1"][(Q, c)]
        certify(f"kappa and sup^2 = escapes at {(Q, c)}", abs(p["kappa"] - e1["kurtosis"]) <= 1e-9 * e1["kurtosis"]
                and abs(p["sup2"] - e1["sup2"]) <= 1e-9 * e1["sup2"], f"kappa {p['kappa']:.4f}", res)
        sp = fr.sup_prediction(Q, c, psi, mu, phi)
        row["sup_pred"] = sp
        row["model"] = fr.spike_ccdf_model(ts, p["mean_m"], meas_m, X, psi, mu, phi, Q, c, levels=levels)
        row["asymptote"] = psi**2 / (2 * X * meas_m * ts * p["mean_m"])
        row["pareto_range"] = (5.0, 0.04 * p["sup2"] / p["mean_m"])
        row["E4_model"] = fr.fourth_moment_model(Q, c, psi, X, mu, phi)
        m3[(Q, c)] = row
        if (Q, c) in CCDF_PAIRS:
            dd = {}
            for tag, a2 in (("primes", a2p), ("hybrid", a2h), ("cramer", a2c)):
                mm = np.where(minor, a2, 0.0)
                mm /= mm.sum()
                dd[tag] = fr.partition_dims(mm, list(range(0, 15)), [2.0])[2.0]
                del mm
            # alpha < 0.9 boxes (eps = 2^-18) lying entirely in the minor arcs
            box_minor = minor.reshape(-1, 32).all(axis=1)
            sel = box_minor & (alpha_full < 0.9)
            mass_minor = a2p[minor].sum()
            boxes_a2 = a2p.reshape(-1, 32).sum(axis=1)
            dd["boxes_lt09"] = {"n": int(sel.sum()), "alpha_min": float(alpha_full[sel].min()) if sel.any() else None,
                                "share": float(boxes_a2[sel].sum() / mass_minor)}
            m2e[(Q, c)] = dd
        log(f"M3/M2e {(Q, c)}")
    res["m3"] = m3
    res["m3_t"] = ts
    res["m2"]["minor_D2"] = m2e


def additions(res, a2p, psi, L2, mu, phi, sea_level):
    oct_ = fr.octave_energy(a2p, X, mu, phi, jmax=10)
    sea_power = L2 - res["hybrid_sumsq"][1000]          # sea energy per unit of alpha in the surrogate
    for o in oct_:
        o["window_measure"] = o["measure"]
        o["spike_model"] = fr.SINC2_HALF * psi**2 / X * o["sum_inv_phi"] / L2
        o["spike_sea_model"] = o["spike_model"] + sea_power * o["window_measure"] / L2
        o["frac"] = o["energy"] / L2
    res["a1"] = {"octaves": oct_, "equal_line": fr.SINC2_HALF * math.log(2) * psi**2 / (X * L2)}
    absS = np.sqrt(a2p)
    sig = np.round(np.arange(0.60, 0.951, 0.05), 2)
    cells = fr.large_value_counts(absS, X, sig)
    spikes = [fr.spike_count(mu, phi, psi / X**s) for s in sig]
    res["a2"] = {"sigma": sig, "cells": cells, "spikes": spikes,
                 "sieve_bound": [2 * (2 * X - 1) * L2 / X ** (2 * s) for s in sig]}
    res["m1_winmax"]["js"] = list(range(3, 21))
    del absS


def stage_xscaling(res):
    out = []
    for Xs, Ms in XS:
        lam = von_mangoldt(Xs)
        psi = float(lam.sum())
        S = load_S("primes", Xs, Ms)
        a2 = S.real**2 + S.imag**2
        del S
        kept = fr.remove_atoms(a2, Xs)
        del a2
        L = fr.kept_L(kept, Xs, psi)
        mu_ = kept / kept.sum()
        del kept
        j = int(round(math.log2(8 * Ms / Xs)))
        m = fr.rank_mass(mu_, j, 20000)                     # every rank (only 1024 boxes at X = 10^4)
        ranks = np.arange(1, len(m) + 1)
        out.append({"X": Xs, "L": L, "j": j, "Neps": Xs * 2**j / Ms, "ranks": ranks, "comp": 2 * L * ranks * m})
        del mu_
        log(f"M2(f) X = {Xs}: L = {L:.3f}")
    res["m2"]["xscaling"] = out


def stage_mandelbrot(res, L, mu, phi):
    t = time.time()
    rows, allok = [], True
    for q in range(2, 41):
        for p in range(1, q // 2 + 1):
            if math.gcd(p, q) != 1:
                continue
            c, per, rot = fr.bulb_centre(p, q)
            croot = fr.cardioid_root(p, q)
            s = fr.leading_size(p, q)
            c_near_root = fr.bulb_continue(p, q, 0.999999 + 0j)
            ok = per == q and rot and abs(c_near_root - croot) <= 1e-2 * s
            allok &= ok
            dh = fr.bulb_diam_half(p, q)
            rows.append({"p": p, "q": q, "theta": p / q, "s": s, "centre": c, "centre_root": abs(c - croot) / s,
                         "diam_half": dh / s, "diam_half_abs": dh, "ok": ok})
    certify("Mandelbrot bulbs q <= 40, p/q <= 1/2 (the rest by conjugation): period, rotation number, root continuation", allok and len(rows) == 245,
            f"{sum(r['ok'] for r in rows)}/{len(rows)} in {time.time() - t:.1f} s", res)
    d12 = fr.bulb_diam_half(1, 2)
    certify("bulb 1/2: diam/2 = 1/4", abs(d12 - 0.25) <= 1e-10, f"{d12:.12f}", res)
    area = {f"{p}/{q}": {"area_radius": fr.bulb_area_radius(p, q) / fr.leading_size(p, q),
                         "diam_half": fr.bulb_diam_half(p, q) / fr.leading_size(p, q)}
            for p, q in [(1, 3), (1, 10), (7, 20), (1, 50), (21, 50), (13, 34)]}
    med = {n: fr.mediant_test(n, D, mu, phi) for n, D in ((8, 120), (13, 200), (21, 300))}
    sizes = []
    for q in range(2, 41):
        bulbs = [r["diam_half_abs"] for r in rows if r["q"] == q]
        sizes.append({"q": q, "bulb": max(bulbs), "spike_mass": (1.0 / (L * phi[q] ** 2)) if mu[q] != 0 else None,
                      "ford": 1.0 / (2 * q * q), "squarefree": bool(mu[q] != 0)})
    logdiv = []
    for Qd in (40, 400, 4000):
        sb = sum(math.sin(math.pi * p / q) / q**2 for q in range(2, Qd + 1) for p in range(1, q) if math.gcd(p, q) == 1)
        logdiv.append({"Q": Qd, "bulbs": sb, "bulb_formula": 12 / math.pi**3 * math.log(Qd) - 0.009,
                       "spikes": float(sum(1.0 / phi[q] for q in range(1, Qd + 1) if mu[q] != 0)),
                       "spike_formula": math.log(Qd) + 1.3326})
    img = fr.mandelbrot_escape((-2.0, 0.55), (-1.12, 1.12), 1600, 1400, 2000)
    res["m5"] = {"bulbs": rows, "area": area, "mediant": med, "sizes": sizes, "logdiv": logdiv, "render": img,
                 "render_extent": (-2.0, 0.55, -1.12, 1.12)}
    log("M5 Mandelbrot set")


def compute() -> dict:
    res: dict = {"X": X, "M": M, "certificates": []}
    lam = von_mangoldt(X)
    psi = float(lam.sum())
    L2 = float(np.sum(lam**2))
    mu, phi = mobius_phi(2 * 10**5)
    res.update(psi=psi, L2=L2)
    certify("spike-count constant: sum_{phi <= 10^4} mu^2 phi / (10^8/2)", 0.995 <= fr.spike_count(mu, phi, 1e4) / 5e7 <= 1.010,
            f"{fr.spike_count(mu, phi, 1e4) / 5e7:.4f}", res)
    log("stage 0")
    stage_collapse(res, lam, psi, mu, phi)
    stage_qstar(res)
    keep = {}
    a2p = stage_primes_grid(res, psi, mu, phi, keep)
    certify("Parseval, primes", abs(a2p.sum() / M - L2) <= 1e-9 * L2, f"{a2p.sum() / M:.6e}", res)
    Sc = load_S("cramer")
    a2c = Sc.real**2 + Sc.imag**2
    wc = cramer_weights(X)
    certify("Parseval, Cramer", abs(a2c.sum() / M - float(np.sum(wc**2))) <= 1e-9 * float(np.sum(wc**2)), "", res)
    del wc
    stage_other_walks(res, "cramer", Sc, keep, spectral=True)
    del Sc
    log("M4 cached Cramer walks")
    W = fr.ramanujan_weights(X, psi, mu, phi, RS)
    res["hybrid_sumsq"] = {R: float(np.sum(W[R] ** 2)) for R in RS}
    sea_level = None                                     # set in spectrum_stats once the atom share is measured
    a2h = None
    for R in RS:
        h = fr.hybrid_weights(W[R], L2, seed=7 + R)
        Sh = S_grid(h, M)
        certify(f"Parseval, hybrid R = {R}", abs(np.sum(Sh.real**2 + Sh.imag**2) / M - float(np.sum(h**2))) <= 1e-9 * float(np.sum(h**2)), "", res)
        stage_other_walks(res, f"hybrid{R}", Sh, keep, variants=False, ps=(0.5, 2.0) if R == 1000 else (2.0,),
                          spectral=(R == 1000))
        if R == 1000:
            a2h = Sh.real**2 + Sh.imag**2
        del Sh, h
    del W
    log("M4 hybrid walks")
    ens = {}
    for sd in range(N_SEEDS):
        Se = S_grid(cramer_weights(X, seed=1000 + sd), M)
        for (Q, c) in WALK_PAIRS:
            k, ph = keep[(Q, c)]
            F = np.cumsum(Se[k] ** 2 * ph / M)
            ens.setdefault((Q, c), []).append({"raw": fr.local_hurst(LAGS, fr.structure_function(F, LAGS, 2.0)),
                                               "drift": fr.local_hurst(LAGS, fr.structure_function(F, LAGS, 2.0, drift=True))})
            del F
        del Se
    res["ensemble"] = ens
    del keep
    log("M4 Cramer ensemble")
    measures = spectrum_stats(res, a2p, a2c, a2h, psi, L2, mu, phi, sea_level)
    sea_level = res["sea_level"]
    log("M2 spectrum")
    minor_arc_stats(res, a2p, a2c, a2h, measures, psi, L2, mu, phi, sea_level)
    del measures
    additions(res, a2p, psi, L2, mu, phi, sea_level)
    res["m1_winmax"]["hybrid"] = fr.window_max_median(np.sqrt(a2h), range(3, 21), psi)
    res["m1_winmax"]["cramer"] = fr.window_max_median(np.sqrt(a2c), range(3, 21), psi)
    del a2p, a2c, a2h
    log("A1, A2, window maxima")
    stage_xscaling(res)
    stage_mandelbrot(res, res["m2"]["L"], mu, phi)
    esc = pickle.load(open(CACHE / "escapes.pkl", "rb"))  # E1's spread of sum_N I_m(N)^2 over target classes
    res["leakage"] = {pr: dict(esc["e1"][pr]["shares"]) for pr in WALK_PAIRS}
    del esc
    res["uncertainty"] = {"sharp": {L: fr.uncertainty_product("sharp", L) for L in (5, 20, 80, 320)},
                          "bump": {L: fr.uncertainty_product("bump", L) for L in (10, 20, 40)}}
    res["runtime_s"] = time.time() - T0
    res["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    with open(CACHE / "fractals.pkl", "wb") as f:
        pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    log("saved cache/fractals.pkl")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="reuse cache/fractals.pkl")
    args = ap.parse_args()
    FIGS.mkdir(exist_ok=True)
    if args.report and (CACHE / "fractals.pkl").exists():
        res = pickle.load(open(CACHE / "fractals.pkl", "rb"))
    else:
        res = compute()
    import figures_fractals
    figures_fractals.make_all(res, FIGS)
    log("figures written")
    import fractals_report
    fractals_report.write(res, HERE / "fractals.md")
    log("fractals.md written")


if __name__ == "__main__":
    main()
