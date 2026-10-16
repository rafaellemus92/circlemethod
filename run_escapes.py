"""Escape routes from the binary L^2 wall: run E0-E5, write figures/fig_X*.png and escapes.md.

    python run_escapes.py            # full run (~4 min on 4 cores, ~1.4 GB peak); needs run_all.py's cache
    python run_escapes.py --report   # redraw figures and escapes.md from cache/escapes.pkl
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

import escapes as es  # noqa: E402
from expsum import e_minus, load_or_compute_S  # noqa: E402
from singular_series import SS2_even_array  # noqa: E402
from sieve import mobius_phi, von_mangoldt  # noqa: E402

CACHE = HERE / "cache"
FIGS = HERE / "figures"
X = 10**6
M = 2**23
QS = (5, 10, 20, 50, 100, 150)
CS = (5, 10, 20)
PAIRS = [(Q, c) for Q in QS for c in CS] + [(158, 20), (223, 10)]
SHORT = [(150, 20), (100, 20), (50, 10), (20, 10), (5, 5)]
NTERS = (X - 1, X + 1, 3 * X // 2 - 1, 3 * X // 2 + 1)
TRUNC_TARGETS = (510_510, 990_990, 999_998)
SCALING_X = (10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6, 2 * 10**6)
KEEP_TARGET_ARRAYS = ((150, 20), (223, 10))

T0 = time.time()


def log(msg: str) -> None:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[{time.time() - T0:7.1f}s, peak {peak:6.0f} MB] {msg}", flush=True)


def grid_size(N: int) -> int:
    return 1 << int(math.ceil(math.log2(8 * N)))


def compute() -> dict:
    res: dict = {"X": X, "M": M, "pairs": PAIRS, "short": SHORT, "nters": NTERS}
    ctx = es.Context(X, M, CACHE, PAIRS)
    res.update(R2=ctx.R2, L2=ctx.L2, psi=ctx.psi, sum_r2_sq=math.fsum(ctx.r2 ** 2))
    log("context ready")

    # ---------------- phase A: needs S itself
    absS = np.abs(ctx.S)
    W = ctx.S.copy()
    for s0 in range(0, M, 1 << 20):
        k = np.arange(s0, min(M, s0 + (1 << 20)), dtype=np.int64)
        W[k] *= e_minus(X, k, M)
    res["e3"] = es.e3_almost_primes(ctx, absS, W, es.almost_prime_sets(X), PAIRS, SHORT)
    del W, absS
    log(f"E3 done ({len(res['e3'])} sets)")
    res["e0"] = {}
    for Q, c in ((150, 20), (20, 10)):
        A = ctx.arc(Q, c)
        IMX = float(np.real(np.sum(ctx.S[A.idx] ** 2 * e_minus(X, A.idx, M)))) / M
        res["e0"][(Q, c)] = es.e0_obstacle(ctx, Q, c, IMX)
    log("E0 done")

    # ---------------- phase B: needs S^2 and |S|^2
    a2 = ctx.S.real ** 2 + ctx.S.imag ** 2
    S2 = ctx.S
    ctx.S = None
    np.multiply(S2, S2, out=S2)
    SS2e = SS2_even_array(X)
    res["e1"] = {}
    for Q, c in PAIRS:
        o = es.e1_exceptional(ctx, S2, a2, Q, c, SS2e)
        if (Q, c) not in KEEP_TARGET_ARRAYS:
            o.pop("SS2_targets")
            o.pop("pm_ratio_targets")
        res["e1"][(Q, c)] = o
    log("E1 + E5(a) done")
    res["e2"] = es.e2_linnik(ctx, S2, a2, PAIRS)
    del S2, a2
    log("E2 done")

    # ---------------- phase C: smooth weight
    main = pickle.load(open(CACHE / "ctx.pkl", "rb"))
    sharp = {(r["Q"], r["c"]): r["sup"] for r in main["rows"]}
    for r in main["scaling_rows"]:
        if r["N"] == X and not r.get("probe"):
            sharp.setdefault((r["Q"], r["c"]), r["sup"])
    res["e4"] = es.e4_smooth(ctx, PAIRS, sharp, NTERS)
    log("E4 done")

    # ---------------- phase D: single-N transmission beyond X = 10^6
    ctx.arcs.clear()
    res["e5_trunc"] = [es.e5_truncated(X, N, 10, ctx.mu, ctx.phi, CACHE) for N in TRUNC_TARGETS]
    SS2big = SS2_even_array(max(SCALING_X))
    scal = []
    for Xs in SCALING_X:
        lam = von_mangoldt(Xs)
        for c in (10, 20):
            S = load_or_compute_S(lam, grid_size(Xs), "primes", CACHE)
            scal.append(es.e5_scaling(Xs, S, lam, ctx.mu, ctx.phi, SS2big, c))
            del S
        log(f"E5 scaling X={Xs} done")
    res["e5_scaling"] = scal
    res["runtime_s"] = time.time() - T0
    res["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    with open(CACHE / "escapes.pkl", "wb") as f:
        pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="reuse cache/escapes.pkl")
    args = ap.parse_args()
    FIGS.mkdir(exist_ok=True)
    if args.report and (CACHE / "escapes.pkl").exists():
        res = pickle.load(open(CACHE / "escapes.pkl", "rb"))
    else:
        res = compute()
    import figures_escapes
    figures_escapes.make_all(res, FIGS)
    log("figures written")
    import escapes_report
    escapes_report.write(res, HERE / "escapes.md")
    log("escapes.md written")


if __name__ == "__main__":
    main()
