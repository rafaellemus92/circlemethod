"""Run every module, write figures/*.png and results.md.

    python run_all.py            # full run (~5 min on 4 cores, ~1.4 GB peak)
    python run_all.py --figures  # re-draw figures and results.md from cache/ctx.pkl
"""
from __future__ import annotations

import argparse
import math
import pickle
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import arcs as ar  # noqa: E402
import expsum as ex  # noqa: E402
import singular_series as ss  # noqa: E402
import sieve as sv  # noqa: E402

CACHE = HERE / "cache"
FIGS = HERE / "figures"

N0 = 10**6
M0 = 2**23
QS = (5, 10, 20, 50, 100, 150)
CS = (5, 10, 20)
NTER = (N0 - 1, N0 + 1)
SCALING_N = (10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6, 2 * 10**6)
SMOOTH_PROBE = 1531530        # = 2·3²·5·7·11·13·17: SS2 = 5.46, near the largest value for N <= 2e6
COMET_MAX = 2 * 10**5
GOLDEN = (math.sqrt(5.0) - 1.0) / 2.0          # frac of the golden ratio
SQRT2 = math.sqrt(2.0) - 1.0                   # frac of sqrt 2

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def grid_size(N: int) -> int:
    """Power of two >= 8N: M > 2N' (alias-free ternary) with ~8 samples per 1/N."""
    return 1 << int(math.ceil(math.log2(8 * N)))


def q_max(N: int, c: int) -> int:
    return math.isqrt((N - 1) // (2 * c))  # largest Q with 2 c Q^2 < N (up to rounding below)


# =============================================================== module 0 + 1

def base_numbers(ctx: dict) -> None:
    N = N0
    is_p = sv.prime_sieve(N)
    lam = sv.von_mangoldt(N)
    psi = sv.chebyshev_psi(lam)
    L2 = sv.sum_lambda_sq(lam)
    R2 = sv.R2_direct(lam, N)
    ctx.update(
        N=N, M=M0, pi_N=int(is_p.sum()), psi=psi, L2=L2, L2_pred=N * (math.log(N) - 1),
        pp_weight=sv.prime_power_weight(lam, is_p), G_N=sv.goldbach_G(N, is_p), R2=R2,
        SS2=ss.SS2_product(N), C2_check=ss.C2_product(), C3=ss.C3_product(),
    )
    log(f"pi={ctx['pi_N']} G={ctx['G_N']} psi={psi:.3f} sumL2={L2:.6e} R2={R2:.6f}")
    assert ctx["pi_N"] == 78498, "sieve broken: pi(10^6) != 78498"
    assert ctx["G_N"] == 5402, "sieve broken: G(10^6) != 5402"
    return is_p, lam


def grid_checks(ctx: dict, lam: np.ndarray, S: np.ndarray) -> None:
    N = ctx["N"]
    chk = ex.sign_convention_checks(S, lam, N, ctx["R2"])
    chk["parseval_grid"] = ex.grid_mass(S)
    chk["parseval_relerr"] = abs(chk["parseval_grid"] - ctx["L2"]) / ctx["L2"]
    chk["R2_grid_relerr"] = abs(chk["grid_R2_right_sign"] - ctx["R2"]) / ctx["R2"]
    ctx["checks"] = chk
    log(f"sign/Parseval/R2 checks: parseval {chk['parseval_relerr']:.1e}, R2 {chk['R2_grid_relerr']:.1e}, "
        f"wrong-sign R2 {chk['grid_R2_wrong_sign']:.2e}")


# =============================================================== module 2 + 3 data

def farey_data(ctx: dict, lam: np.ndarray, is_p: np.ndarray, mu, phi) -> None:
    psi = ctx["psi"]
    rows = []
    for q in range(1, 31):
        pr = sv.psi_residues(lam, q)
        for a in range(q):
            if math.gcd(a, q) != 1:
                continue
            Sv = ex.S_rational(a, q, pr)
            rows.append({"a": a, "q": q, "S": Sv, "absS": abs(Sv), "mu": int(mu[q]), "phi": int(phi[q]),
                         "pred": psi * abs(mu[q]) / phi[q], "pred_signed": psi * mu[q] / phi[q]})
    ctx["farey"] = rows
    p4 = sv.psi_residues(lam, 4)
    t4 = sv.theta_residues(lam, is_p, 4)
    ctx["mod4"] = {"psi41": p4[1], "psi43": p4[3], "psi40": p4[0], "psi42": p4[2],
                   "theta41": t4[1], "theta43": t4[3],
                   "S_quarter": ex.S_rational(1, 4, p4)}
    log(f"S(1/4) = {ctx['mod4']['S_quarter']:.3f};  psi(N;4,3)-psi(N;4,1) = {p4[3] - p4[1]:.3f}; "
        f"theta race = {t4[3] - t4[1]:.3f}")


def phasor_data(ctx: dict, lam: np.ndarray, n, w) -> None:
    out = []
    specs = [("1/3", 1, 3), ("1/4", 1, 4), ("1/5", 1, 5), ("1/6", 1, 6), ("1/30", 1, 30),
             ("frac(golden ratio)", None, GOLDEN), ("frac(sqrt 2)", None, SQRT2)]
    for name, a, q in specs:
        if a is not None:
            frac = ((n * a) % q) / q
            alpha = a / q
        else:
            alpha = q
            frac = (n.astype(float) * alpha) % 1.0
        b = np.rint(frac * 360).astype(np.int64) % 360
        hist = np.bincount(b, weights=w, minlength=360)
        Sv = complex(np.sum(w * np.exp(2j * np.pi * frac)))
        out.append({"name": name, "a": a, "q": q if a is not None else None, "alpha": alpha,
                    "hist": hist, "S": Sv})
    ctx["phasors"] = out


def zoom_data(ctx: dict, n, w, mu, phi) -> None:
    N, psi = ctx["N"], ctx["psi"]
    betas = np.linspace(-20 / N, 20 / N, 2401)
    zooms = []
    for a, q in ((1, 3), (1, 6), (1, 4)):
        Sv = ex.S_near_rational(a, q, betas, n, w)
        model = (mu[q] / phi[q]) * ex.dirichlet(betas, N)
        zooms.append({"a": a, "q": q, "beta": betas, "S": Sv, "model": model})
    Sg = ex.S_direct(GOLDEN + betas, n, w)
    conv = []
    for q in range(1, 2001):
        a = round(GOLDEN * q)
        d = GOLDEN - a / q
        if math.gcd(a, q) == 1 and abs(d) <= 20 / N:
            conv.append({"a": a, "q": q, "N_beta": -N * d, "height": psi * abs(mu[q]) / phi[q]})
    for f in conv:
        if f["height"] >= 1000:
            f["S_exact"] = complex(ex.S_near_rational(f["a"], f["q"], [0.0], n, w)[0])
            f["mu"] = int(mu[f["q"]])
    jmax = int(np.argmax(np.abs(Sg)))
    ctx["zooms"] = zooms
    ctx["zoom_golden"] = {"beta": betas, "S": Sg, "fractions": conv, "max_Nbeta": float(N * betas[jmax])}
    log("zooms done")


def envelope(absS: np.ndarray, bins: int) -> tuple[np.ndarray, np.ndarray]:
    r = absS.reshape(bins, -1)
    return r.max(axis=1), np.median(r, axis=1)


# =============================================================== module 5 figure data

def minor_window(ctx: dict, S: np.ndarray, mu, phi, n=None, w=None) -> None:
    """A 300/N stretch lying in the minor arcs of every (Q, c) tabulated at N = 10^6.

    That means clear of q <= 158 by 20/N (the grid pairs and (158, 20)) and of q <= 223 by
    10/N (the Qmax pair (223, 10)).
    """
    N, M, psi = ctx["N"], ctx["M"], ctx["psi"]
    half = 150 / N
    xs = np.arange(0.30, 0.45, 2e-6)
    ok = np.ones(len(xs), dtype=bool)
    for q in range(1, 224):
        x = xs * q
        ok &= np.abs(x - np.rint(x)) / q > half + (20 if q <= 158 else 10) / N
    cand = xs[ok]
    # Such windows always sit in a wide Farey gap next to some low-q rational, and the q <= 300
    # fractions inside are its "satellites".  Pick the window with the fewest of them (at least two
    # squarefree, so there are crests to look at), then the largest distance from the window edges.
    best, best_key = float(cand[len(cand) // 2]), None
    for x in cand:
        fr = [(a / q, q) for q in range(151, 301) for a in [round(x * q)]
              if math.gcd(a, q) == 1 and abs(x - a / q) <= half]
        n_sf = sum(mu[q] != 0 for _, q in fr)
        if n_sf < 2:
            continue
        margin = min(half - abs(f - x) for f, _ in fr)
        key = (len(fr), -margin)
        if best_key is None or key < best_key:
            best, best_key = float(x), key
    x0 = float(best)
    k0 = int(round(x0 * M))
    hw = int(round(half * M))
    k = np.arange(k0 - hw, k0 + hw + 1, dtype=np.int64)
    Sk = S[k]
    integrand = Sk * Sk * ex.e_minus(N, k, M)
    fracs = []
    for q in range(1, 301):
        a = round(x0 * q)
        if math.gcd(a, q) == 1 and abs(x0 - a / q) <= half:
            fr = {"a": a, "q": q, "alpha": a / q, "height": psi * abs(mu[q]) / phi[q], "mu": int(mu[q])}
            if n is not None:
                fr["S_exact"] = complex(ex.S_near_rational(a, q, [0.0], n, w)[0])
            fracs.append(fr)
    mass = np.abs(Sk) ** 2
    order = np.argsort(-mass)
    peaks = []                                   # the two tallest local peaks and their nearest fraction
    for j in order:
        if all(abs(int(j) - p) > 8 for p in peaks):
            peaks.append(int(j))
        if len(peaks) == 2:
            break
    peak_info = []
    for j in peaks:
        al = k[j] / M
        near = min(fracs, key=lambda f: abs(f["alpha"] - al)) if fracs else None
        peak_info.append({"alpha": al, "ratio": float(mass[j] / mass.mean()),
                          "near": (near["a"], near["q"]) if near else None,
                          "N_dist": float(N * abs(near["alpha"] - al)) if near else None})
    ctx["window"] = {"x0": x0, "alpha": k / M, "mass": mass, "integrand": integrand, "peaks": peak_info,
                     "cum": np.cumsum(integrand.real) / M, "cum_imag": np.cumsum(integrand.imag) / M,
                     "cum_mass": np.cumsum(np.abs(Sk) ** 2) / M, "fractions": fracs}
    log(f"minor window at {x0:.6f}, fractions q<=300 inside: {[(f['a'], f['q']) for f in fracs]}")


def staircase(ctx: dict, S: np.ndarray, mu, phi, stride: int = 16) -> None:
    N, M = ctx["N"], ctx["M"]
    F = np.empty(M, dtype=complex)
    for s in range(0, M, 1 << 20):
        k = np.arange(s, min(M, s + (1 << 20)), dtype=np.int64)
        F[k] = S[k] * S[k] * ex.e_minus(N, k, M)
    np.cumsum(F, out=F)
    F /= M
    jumps = []
    for q in range(1, 7):
        for a in range(q):
            if math.gcd(a, q) != 1:
                continue
            pred = N * mu[q] ** 2 * math.cos(2 * math.pi * ((a * N) % q) / q) / phi[q] ** 2
            pred_im = -N * mu[q] ** 2 * math.sin(2 * math.pi * ((a * N) % q) / q) / phi[q] ** 2
            if a == 0:
                meas = F[int(20 * M / N)] + (F[-1] - F[M - int(20 * M / N) - 1])
            else:
                kc = a * M / q
                lo, hi = int(math.floor(kc - 20 * M / N)), int(math.ceil(kc + 20 * M / N))
                meas = F[hi] - F[lo]
            jumps.append({"a": a, "q": q, "pred": pred, "pred_im": pred_im, "meas": meas.real,
                          "meas_im": meas.imag})
    ctx["stair"] = {"alpha": np.arange(0, M, stride) / M, "F": F[::stride].copy(), "final": F[-1],
                    "jumps": jumps}
    del F


# =============================================================== modules 4-6 driver

def module_main(ctx: dict) -> None:
    is_p, lam = base_numbers(ctx)
    N, M = ctx["N"], ctx["M"]
    n, w = ex.nonzero_terms(lam)
    S = ex.load_or_compute_S(lam, M, "primes", CACHE)
    log(f"S on grid M=2^{int(math.log2(M))} ready")
    grid_checks(ctx, lam, S)

    mu, phi = sv.mobius_phi(10**6)
    SS2_part = ss.SS2_partial(N, 10**6, mu, phi)
    ctx["SS2_part"] = SS2_part[: 10**5 + 1].copy()
    ctx["SS2_sum_1e6"] = float(SS2_part[-1])
    ctx["SS3"] = {Np: ss.SS3_product(Np) for Np in NTER}
    ctx["SS3_part"] = {Np: ss.SS3_partial(Np, 10**4, mu, phi) for Np in NTER}
    ctx["factor"] = {Np: sv.factorize(Np) for Np in (N, *NTER)}

    r2 = sv.autoconvolve(lam)
    ctx["R3"] = {Np: sv.R3_exact(lam, Np, r2) for Np in NTER}
    ctx["R3_grid"] = {Np: ex.grid_power(S, 3, Np) for Np in NTER}
    ctx["R2_from_conv"] = float(r2[N])
    del r2
    log("R3: " + ", ".join(f"{Np}: {ctx['R3'][Np]:.6e} (grid rel {abs(ctx['R3_grid'][Np] - ctx['R3'][Np]) / ctx['R3'][Np]:.1e})" for Np in NTER))

    pairs = [(Q, c) for Q in QS for c in CS]
    ctx["skipped_pairs"] = [(Q, c) for Q, c in pairs if not ar.admissible(N, Q, c)]
    J_cache: dict = {}
    ctx["rows"] = ar.arc_rows(S=S, N=N, pairs=pairs, mu=mu, phi=phi, SS2_part=SS2_part,
                              SS2_full=ctx["SS2"], R2=ctx["R2"], L2=ctx["L2"], psi=ctx["psi"], n=n, w=w,
                              Nter=NTER, R3=ctx["R3"], SS3=ctx["SS3"], J_cache=J_cache)
    ctx["J"] = {c: J_cache[(N, c)] for c in CS}
    log("arc table done")

    farey_data(ctx, lam, is_p, mu, phi)
    phasor_data(ctx, lam, n, w)
    zoom_data(ctx, n, w, mu, phi)
    absS = np.abs(S)
    ctx["env_max"], ctx["env_med"] = envelope(absS, 1 << 14)
    del absS
    minor_window(ctx, S, mu, phi, n, w)
    staircase(ctx, S, mu, phi)
    log("staircase done")
    del S

    # ---------------- module 8: Cramer random model on the same grid
    wR = sv.cramer_weights(N)
    SR = ex.load_or_compute_S(wR, M, "cramer", CACHE)
    nR, wRn = ex.nonzero_terms(wR)
    R2R = sv.R2_direct(wR, N)
    L2R = float(np.sum(wR ** 2))
    psiR = float(wR.sum())
    ctx["rand"] = {"R2": R2R, "L2": L2R, "psi": psiR, "count": int(np.count_nonzero(wR)),
                   "E_R2": (N - 6) + math.log(N / 2), "parseval": ex.grid_mass(SR),
                   "R2_grid": ex.grid_power(SR, 2, N).real}
    absSR = np.abs(SR)
    ctx["rand"]["env_max"], ctx["rand"]["env_med"] = envelope(absSR, 1 << 14)
    del absSR
    ctx["rand"]["farey"] = [{"a": r["a"], "q": r["q"], "absS": abs(SR[0]) if r["q"] == 1 else
                             abs(ex.S_direct([r["a"] / r["q"]], nR, wRn)[0])} for r in ctx["farey"] if r["q"] <= 12]
    ones = np.ones(10**5 + 1)
    ones[0] = 0.0
    ctx["rand"]["rows"] = ar.arc_rows(S=SR, N=N, pairs=pairs, mu=mu, phi=phi, SS2_part=ones, SS2_full=1.0,
                                      R2=R2R, L2=L2R, psi=psiR, n=nR, w=wRn, J_cache=J_cache, label="cramer")
    del SR
    # robustness of the random control: a few more independent draws (no sup search, two pairs)
    seeds = []
    for seed in (1, 2, 3, 4, 5):
        wS = sv.cramer_weights(N, seed=seed)
        SRs = ex.S_grid(wS, M)
        R2s = sv.R2_direct(wS, N)
        L2s = float(np.sum(wS ** 2))
        row = {"seed": seed, "R2": R2s, "L2": L2s}
        for Q, c in ((20, 10), (150, 20)):
            Ar = ar.MajorArcs(N, M, Q, c, mu, phi)
            Pm, Im = ar.minor_integrals(SRs, Ar.mask(), N)
            row[(Q, c)] = {"P_m_over_R2": Pm / R2s, "I_m_over_R2": Im.real / R2s}
        seeds.append(row)
        del SRs
    ctx["rand"]["seeds"] = seeds
    log("random model done")


# =============================================================== module 7

def module_scaling(ctx: dict) -> None:
    mu, phi = sv.mobius_phi(10**5)
    rows = []
    meta = []
    for N in (*SCALING_N, SMOOTH_PROBE):
        M = grid_size(N)
        lam = sv.von_mangoldt(N)
        n, w = ex.nonzero_terms(lam)
        S = ex.load_or_compute_S(lam, M, "primes", CACHE)
        psi, L2 = sv.chebyshev_psi(lam), sv.sum_lambda_sq(lam)
        R2 = sv.R2_direct(lam, N)
        nter = (N - 1, N + 1)
        r2 = sv.autoconvolve(lam)
        R3 = {Np: sv.R3_exact(lam, Np, r2) for Np in nter}
        del r2
        SS3 = {Np: ss.SS3_product(Np) for Np in nter}
        SS2_full = ss.SS2_product(N)
        Qmax = max(QS + (q_max(N, 10), q_max(N, 20)))
        SS2_part = ss.SS2_partial(N, Qmax, mu, phi)
        pairs = [(20, 10), (50, 10), (q_max(N, 10), 10), (q_max(N, 20), 20)]
        rr = ar.arc_rows(S=S, N=N, pairs=pairs, mu=mu, phi=phi, SS2_part=SS2_part, SS2_full=SS2_full,
                         R2=R2, L2=L2, psi=psi, n=n, w=w, Nter=nter, R3=R3, SS3=SS3)
        for r in rr:
            r["probe"] = N == SMOOTH_PROBE
            r["series"] = {(20, 10): "Q=20, c=10", (50, 10): "Q=50, c=10"}.get(
                (r["Q"], r["c"]), f"Q=Qmax, c={r['c']}")
            if (r["Q"], r["c"]) in ((20, 10), (50, 10)) and r["Q"] == q_max(N, r["c"]):
                r["series"] = f"Q=Qmax, c={r['c']}"
        rows += rr
        meta.append({"N": N, "M": M, "psi": psi, "L2": L2, "R2": R2, "SS2": SS2_full, "R3": R3, "SS3": SS3,
                     "parseval_relerr": abs(ex.grid_mass(S) - L2) / L2,
                     "R2_grid_relerr": abs(ex.grid_power(S, 2, N).real - R2) / R2,
                     "R3_grid_relerr": {Np: abs(ex.grid_power(S, 3, Np).real - R3[Np]) / R3[Np] for Np in nter},
                     "factor": {Np: sv.factorize(Np) for Np in (N, *nter)}})
        del S
        log(f"scaling N={N} done ({len(rr)} rows)")
    ctx["scaling_rows"] = rows
    ctx["scaling_meta"] = meta


# =============================================================== module 9

def module_comet(ctx: dict) -> None:
    Ns, G = sv.goldbach_comet(COMET_MAX)
    is_p = sv.prime_sieve(COMET_MAX)
    spot = [10, 100, 1000, 9998, 100000, 199998]
    spot_ok = all(G[(m - 4) // 2] == sv.goldbach_G(m, is_p) for m in spot)
    SS2a = ss.SS2_even_array(COMET_MAX)[Ns]
    keep = Ns >= 8
    Ns, G, SS2a = Ns[keep], G[keep], SS2a[keep]
    pred = SS2a * 0.5 * ss.hl_integral(Ns)
    ctx["comet"] = {"N": Ns, "G": G, "SS2": SS2a, "pred": pred, "spot_ok": spot_ok}
    log(f"comet done (spot checks vs direct count: {spot_ok})")


# =============================================================== results.md

ANALOGY = '## Where the oscillator picture is literal, and where it is only an analogy\n\n* **Literal.** S(α) is exactly a sum of phasors Λ(n)e(nα). "Resultant", "channel", "constructive /\n  destructive", "random walk" describe the arithmetic exactly, and the Ramanujan sum c_q(a) is literally\n  the resultant of the φ(q) reduced-class phasors e(ar/q).\n* **Not literal: there is no coupling.** The oscillators do not act on each other. Nothing in S(α)\n  plays the role of a coupling strength K, nothing evolves in α, and there is no Adler/Kuramoto\n  threshold |Δω| < K. "Phase-locking" at a/q means commensurability: the phase 2πna/q depends only on\n  n mod q, so the oscillators fall into channels because of arithmetic (primes spread evenly over the\n  reduced residue classes), not because they pull each other into step.\n* **What plays the role of a locking window** is the coherence window |β| ≲ 1/N: detuned by β, the\n  channel phases fan out over 2πNβ across n ≤ N and the resultant decays along the Dirichlet kernel. That\n  is a finite-observation-length effect (the linewidth of a pulse train watched for N periods), not a\n  balance between detuning and coupling.\n* **Look-alikes.** The Farey forest resembles the Arnold-tongue diagram of a circle map, and F(α) in\n  Fig F resembles a devil\'s staircase: structure at every rational, weaker as q grows. The mechanism is\n  different — tongue widths come from coupling, spike heights here from μ(q)/φ(q).\n* **"Beat"** in the captions means demodulation: S² is multiplied by the reference e(−Nα) and averaged\n  over α, which keeps exactly the pairs n₁ + n₂ = N.\n'

def _e(x: float, d: int = 3) -> str:
    return f"{x:.{d}e}".replace("e+0", "e").replace("e-0", "e-").replace("e+", "e")


def _f(x: float, d: int = 3) -> str:
    return f"{x:,.{d}f}"


def _b3(r: dict, d: int) -> float:
    """B3/main for N' = N + d."""
    return r[f"B3_main_{r['N'] + d}"]


def _table(head: list[str], rows: list[list[str]]) -> str:
    esc = lambda cells: [str(c).replace("|", "\\|") for c in cells]  # noqa: E731  (| inside cells)
    out = ["| " + " | ".join(esc(head)) + " |", "|" + "|".join("---:" for _ in head) + "|"]
    out += ["| " + " | ".join(esc(r)) + " |" for r in rows]
    return "\n".join(out)


def write_results(ctx: dict, path: Path) -> None:
    N, psi, L2, R2, SS2 = ctx["N"], ctx["psi"], ctx["L2"], ctx["R2"], ctx["SS2"]
    chk = ctx["checks"]
    rows = ctx["rows"]
    rmap = {(r["Q"], r["c"]): r for r in rows}
    big = rmap[(150, 20)]
    Nm, Np = NTER
    far = ctx["farey"]
    sf_dev = max(abs(r["absS"] / r["pred"] - 1) for r in far if r["mu"] != 0)
    sf_dev_q = max(far, key=lambda r: abs(r["absS"] / r["pred"] - 1) if r["mu"] else -1)
    sign_ok = all(np.sign(r["S"].real) == np.sign(r["mu"]) for r in far if r["mu"] != 0)
    nsf_max = max(r["absS"] / psi for r in far if r["mu"] == 0)
    m4 = ctx["mod4"]
    W = ctx["window"]
    win_mean = float(np.mean(W["mass"]))
    sea_med = float(np.median(ctx["env_med"]))
    zg = ctx["zoom_golden"]
    gold_max = float(np.max(np.abs(zg["S"])))
    gold_top = max(zg["fractions"], key=lambda f: f["height"])
    gold_above = [f for f in zg["fractions"] if f["height"] > math.sqrt(L2)]
    gold_locked = [f for f in zg["fractions"] if "S_exact" in f and f["mu"] != 0]
    gold_sign_ok = sum(np.sign(f["S_exact"].real) == np.sign(f["mu"]) for f in gold_locked)
    z13 = ctx["zooms"][0]
    lobe = np.abs(z13["beta"]) * N <= 10
    z_dev = float(np.max(np.abs(z13["S"][lobe] - z13["model"][lobe])) / (psi / 2))
    ph = {p["name"]: p["S"] for p in ctx["phasors"]}
    st = ctx["stair"]
    jmax = max(abs(j["meas"] - j["pred"]) / N for j in st["jumps"])
    sel = (st["alpha"] >= 0.777) & (st["alpha"] <= 0.782)
    plateau_range = float(np.ptp(st["F"].real[sel]) / N)
    R = ctx["rand"]
    rr = {(r["Q"], r["c"]): r for r in R["rows"]}
    sc = [r for r in ctx["scaling_rows"]]
    probe = [r for r in sc if r.get("probe")]
    com = ctx["comet"]
    ratio = com["G"] / com["pred"]
    edges = [1e3, 5e3, 1e4, 2e4, 5e4, 1e5, 1.5e5, 2.0e5 + 1]
    bins = [(lo, hi, float(ratio[(com["N"] >= lo) & (com["N"] < hi)].mean())) for lo, hi in zip(edges[:-1], edges[1:])]
    g6 = 5402 / (ss.SS2_product(10**6) * 0.5 * ss.hl_integral([10**6])[0])
    Pm_lo = min(r["P_m_over_R2"] for r in rows)
    Pm_hi = max(r["P_m_over_R2"] for r in rows)
    n6 = [r for r in sc if r["N"] == N0 and r["series"].startswith("Q=Qmax")]
    Pm_lo_all = min([Pm_lo] + [r["P_m_over_R2"] for r in n6])
    im_q50 = [abs(r["I_m"]) / R2 for r in rows if r["Q"] >= 50]
    res_q50 = [abs(r["I_m_resid"]) / R2 for r in rows if r["Q"] >= 50]
    rp_q50 = [r["rp_scale"] / R2 for r in rows if r["Q"] >= 50]
    me_all = [r["model_err"] / R2 for r in rows]
    thin = [r for r in rows if r["measure"] <= 0.07]
    heur_small = max(abs(r["P_m_over_R2"] - r["pred_Pm_over_R2"]) for r in rows if r["Q"] <= 20)
    heurx_small = max(abs(r["P_m_over_R2"] - r["pred_exact_Pm_over_R2"]) for r in rows if r["Q"] <= 20)
    qmx = {r["N"]: r for r in sc if r["series"] == "Q=Qmax, c=10" and not r.get("probe")}
    qmx20 = {r["N"]: r for r in sc if r["series"] == "Q=Qmax, c=20" and not r.get("probe")}
    fix20 = {r["N"]: r for r in sc if r["series"] == "Q=20, c=10" and not r.get("probe")}
    ns = sorted(qmx)
    small_resid = max(abs(r["I_m_resid"]) / r["R2"] for r in sc if r["N"] <= 3 * 10**5)
    pos_rows = [r for r in sc if r["I_M"] - r["P_m"] > 0]
    thr = {n: ss.SS3_product(n) * n * n / 2 / L2 for n in NTER}
    sup5, sup10 = rmap[(150, 5)]["sup"], rmap[(150, 10)]["sup"]
    d_ss = ctx["SS2_part"][1:] - SS2
    sg = np.sign(d_ss)
    flips = [int(i + 1) + 1 for i in np.flatnonzero(sg[1:] * sg[:-1] < 0)]   # Q at which the sign changes
    early_flips = [q for q in flips if q <= 10**4]
    late_flip = next((q for q in flips if q > 10**4), None)
    seeds = R.get("seeds", [])
    seed_R2 = [s["R2"] / N for s in seeds] + [R["R2"] / N]
    seed_pm20 = [s[(20, 10)]["P_m_over_R2"] for s in seeds] + [rr[(20, 10)]["P_m_over_R2"]]
    seed_pm150 = [s[(150, 20)]["P_m_over_R2"] for s in seeds] + [rr[(150, 20)]["P_m_over_R2"]]
    t = []
    A = t.append

    A("# Results — the circle method for binary Goldbach at N = 10^6, made visible\n")
    A("Exploration and visualisation only: nothing below is, or is meant as, a proof or a proof strategy. "
      "All numbers are produced by `python run_all.py` (one run, %.0f s on 4 cores); every figure is in `figures/`.\n"
      % ctx["runtime_s"])

    # ------------------------------------------------------------ anchors
    A("## 0. Anchors and quadrature certificates\n")
    A(_table(["quantity", "value", "check"], [
        ["π(10^6)", f"{ctx['pi_N']:,}", "= 78,498 ✓"],
        ["G(10^6) (direct count)", f"{ctx['G_N']:,}", "= 5,402 ✓"],
        ["ψ(N) = S(0)", _f(psi, 3), f"ψ/N − 1 = {psi / N - 1:+.2e} (within 0.2% ✓)"],
        ["prime-power weight Σ_{k≥2} θ(N^{1/k})", _f(ctx["pp_weight"], 1), "O(√N), 0.11% of ψ"],
        ["Σ Λ(n)²", _e(L2, 6), f"N(log N − 1) = {_e(ctx['L2_pred'], 6)}"],
        ["R2(N) direct", _f(R2, 6), f"R2/(SS2·N) = {R2 / (SS2 * N):.6f}"],
        ["SS2(N) = 2C2·(5−1)/(5−2), C2 = literature value", f"{SS2:.10f}",
         f"Σ_{{q≤10^6}} = {ctx['SS2_sum_1e6']:.10f}"],
        ["C2 (our product over p ≤ 10^7)", f"{ctx['C2_check']:.10f}", f"literature {ss.C2_KNOWN:.10f}"],
        ["Parseval on grid (1/M)Σ|S|²", _e(chk["parseval_grid"], 12), f"rel. err {chk['parseval_relerr']:.1e}"],
        ["(1/M)Σ S² e(−Nk/M)", _f(chk["grid_R2_right_sign"], 6), f"= R2, rel. err {chk['R2_grid_relerr']:.1e}"],
        ["same with the wrong FFT sign", _e(chk["grid_R2_wrong_sign"], 2), "0 to rounding (R2 ≈ 1.8×10⁶)"],
        [f"R3({Nm:,}) via FFT convolution", _e(ctx["R3"][Nm], 9),
         f"grid ∫S³e(−N'α): rel. err {abs(ctx['R3_grid'][Nm].real - ctx['R3'][Nm]) / ctx['R3'][Nm]:.1e}"],
        [f"R3({Np:,}) via FFT convolution", _e(ctx["R3"][Np], 9),
         f"grid ∫S³e(−N'α): rel. err {abs(ctx['R3_grid'][Np].real - ctx['R3'][Np]) / ctx['R3'][Np]:.1e}"],
    ]))
    A("\n**FFT sign convention.** `numpy.fft.fft` uses exp(−2πi nk/M). We need S(k/M) = Σ Λ(n) e(+nk/M), so "
      "S = M·ifft(Λ) = conj(fft(Λ)), evaluated with the half-length `rfft` (`expsum.S_grid`). S(0) = ψ, "
      f"S(1/2) = 2⌊log₂N⌋log 2 − ψ = {chk['S_half']:,.4f} (predicted {chk['S_half_pred']:,.4f}) and S(−α) = conj S(α) "
      "hold for either sign, so they do not settle it. Two checks do: "
      f"Im S(1/M) = {chk['Im_S_1_over_M_grid']:,.4f} on the grid vs {chk['Im_S_1_over_M_direct']:,.4f} by direct "
      "summation (same sign), and the binary grid integral returns R2(N) exactly with the chosen sign but "
      f"{chk['grid_R2_wrong_sign']:.1e} with the other (the wrong sign evaluates the autoconvolution at M − N > 2N, "
      "where it vanishes). One quadrature rule is used throughout: the rectangle rule on k/M, k = 0..M−1, "
      f"M = 2^23 = {M0:,}. On the full circle it is exact for R2 (needs M > N), for R3 (needs M > 2N') and "
      "for Parseval; on arcs it is a Riemann sum with membership decided exactly in integers, so "
      "I_M + I_m = R2 to rounding (max |closure| = %.1e).\n" % max(abs(r["closure"]) for r in rows))

    # ------------------------------------------------------------ headline
    A("## Headline table\n")
    A("All (N, Q, c) evaluated. B3/main is the ternary ratio for N' = N−1 / N+1 (odd), using the same S. "
      "For N = 10^6 all 18 grid pairs (Q, c) satisfy 2cQ² < N (largest: 2·20·150² = 900,000), so nothing was skipped; "
      "the rows (158, 20) and (223, 10) for N = 10^6 are the Qmax pairs of §7. Rows marked * are the smooth-N "
      "probe (§7).\n")
    hrows = []
    seen = set()
    for r in rows + sorted([r for r in sc], key=lambda r: (r["N"], r["Q"], r["c"])):
        key = (r["N"], r["Q"], r["c"])
        if key in seen:
            continue
        seen.add(key)
        hrows.append([f"{r['N']:,}{'*' if r.get('probe') else ''}", str(r["Q"]), str(r["c"]), _e(r["R2"], 4),
                      _e(r["SS2N"], 4), _e(r["I_M"], 4), f"{r['I_m']:+.3e}", _e(r["P_m"], 3),
                      f"{r['P_m_over_R2']:.3f}", _e(r["sup"], 3),
                      f"{_b3(r, -1):.3f} / {_b3(r, +1):.3f}"])
    A(_table(["N", "Q", "c", "R2", "SS2·N", "I_M", "I_m", "P_m", "P_m/R2", "sup_m|S|", "B3/main N−1 / N+1"],
             hrows))

    # ------------------------------------------------------------ module 2/3
    A("\n## 1–3. Spikes, zooms, phasors\n")
    A(f"* Squarefree q ≤ 30: max | |S(a/q)| / (ψ/φ(q)) − 1 | = {sf_dev:.2e} (worst a/q = {sf_dev_q['a']}/{sf_dev_q['q']}); "
      f"sign of Re S(a/q) equals μ(q) for every one: {sign_ok}.")
    A(f"* Non-squarefree q ≤ 30: max |S(a/q)|/ψ = {nsf_max:.2e}.")
    A(f"* S(1/4) = {m4['S_quarter'].real:.4f} {m4['S_quarter'].imag:+.4f} i. Real part = ψ(N;4,0) − ψ(N;4,2) = "
      f"{m4['psi40'] - m4['psi42']:.4f} (powers of 2); imaginary part = ψ(N;4,1) − ψ(N;4,3) = "
      f"{m4['psi41'] - m4['psi43']:.4f}. So |Im S(1/4)| is exactly the mod-4 race ψ(N;4,3) − ψ(N;4,1) = "
      f"{m4['psi43'] - m4['psi41']:.1f}. For primes alone the race is θ(N;4,3) − θ(N;4,1) = "
      f"{m4['theta43'] - m4['theta41']:.1f} ≈ 1.8√N. About θ(√N) ≈ 956 of that is Chebyshev's bias proper — the "
      "prime squares p² ≡ 1 (mod 4) that θ leaves out; Λ puts them back and so removes the bias. What remains in "
      f"the ψ race, {m4['psi43'] - m4['psi41']:.0f} ≈ 0.84√N, is the oscillation driven by the zeros of L(s, χ₄), "
      "which happens to be positive at N = 10^6 and has no preferred sign.")
    A(f"* Zoom around 1/3: on |Nβ| ≤ 10 the measured S differs from the model −v(β)/2 by at most "
      f"{z_dev:.2e}·(ψ/2).")
    A(f"* Golden-ratio window: max |S| = {gold_max:,.0f} = {gold_max / math.sqrt(L2):.2f}·√ΣΛ², at N·β = "
      f"{zg['max_Nbeta']:+.2f}, on the crest of {gold_top['a']}/{gold_top['q']} (model ψ/φ(q) = {gold_top['height']:,.0f}, "
      f"N·β = {gold_top['N_beta']:+.2f}). {len(gold_above)} model crests in the window exceed √ΣΛ² = {math.sqrt(L2):,.0f}. "
      f"Exact S(a/q) at the {len(gold_locked)} squarefree fractions with crest ≥ 10³: real part has the sign of μ(q) "
      f"in {gold_sign_ok} cases — "
      + "; ".join(f"S({f['a']}/{f['q']}) = {f['S_exact'].real:+,.0f} {f['S_exact'].imag:+,.0f}i (μ = {f['mu']:+d})"
                  for f in sorted(gold_locked, key=lambda f: -f['height'])[:5]) + ".")
    A("* Phasor resultants: " + "; ".join(f"S({k}) = {v.real:,.0f} {v.imag:+,.0f}i" for k, v in ph.items()) + ".\n")

    # ------------------------------------------------------------ module 4
    A("## 4. Major arcs vs truth (binary, N = 10^6)\n")
    A("(i) I_M from the grid; (ii) the same with S replaced on each arc by (μ(q)/φ(q))·v(β); (iii) SS2_Q(N)·N; "
      "(iv) SS2(N)·N from the product; (v) R2(N) exact. Relative errors are against R2. J(c)/N is the fraction of "
      "∫v² e(−Nβ) inside |β| ≤ c/N: " + ", ".join(f"c={c}: {ctx['J'][c] / N:.5f}" for c in CS)
      + " (≈ 1 − 1/(π²c)).\n")
    A(_table(["Q", "c", "|M|", "(i) I_M", "(ii) model", "(iii) SS2_Q·N", "(iv) SS2·N", "(v) R2", "(i)/R2−1",
              "(ii)/R2−1", "(iii)/R2−1", "(iv)/R2−1"],
             [[str(r["Q"]), str(r["c"]), f"{r['measure']:.4f}", _e(r["I_M"], 5), _e(r["I_M_model"], 5),
               _e(r["SS2QN"], 5), _e(r["SS2N"], 5), _e(r["R2"], 5), f"{r['I_M'] / R2 - 1:+.4f}",
               f"{r['I_M_model'] / R2 - 1:+.4f}", f"{r['SS2QN'] / R2 - 1:+.4f}", f"{r['SS2N'] / R2 - 1:+.4f}"]
              for r in rows]))
    A("\nAgreement is at the few-percent level once Q ≥ 20 and c ≥ 10 (and ≤ 1.5% for Q ≥ 50), as expected. "
      "The difference (i) − (ii) is the error of the major-arc model on the major arcs themselves "
      f"({min(me_all):+.4f} to {max(me_all):+.4f} of R2; part of it is simply ψ ≠ N, (ψ/N)² − 1 = {(psi / N) ** 2 - 1:+.1e}). "
      "SS2_Q(N) for Q = 5, 10, 20, 50, 100, 150: " + ", ".join(f"{ctx['SS2_part'][q]:.6f}" for q in QS) + ".\n")

    # ------------------------------------------------------------ module 5
    A("## 5. The minor arcs — the headline\n")
    A("P_m = ∫_m |S|² (mass) and I_m = ∫_m S² e(−Nα) (signed), both summed directly over the minor-arc grid points. "
      "‘heur.’ is the no-fit estimate P_m/R2 ≈ (log N − log Q − 2.33)/SS2(N), with 2.33 = 1 + γ + Σ_p log p/(p(p−1)) "
      "from Σ_{q≤Q} μ²(q)/φ(q) ≈ log Q + 1.33; ‘heur. exact’ uses the exact partial sum and the kernel capture: "
      "(ΣΛ² − J(c)·Σ_{q≤Q} μ²(q)/φ(q))/R2. ‘leak’ = (SS2 − SS2_Q)(N−1) + SS2_Q·(N−1−J(c)) is the deterministic part "
      "of I_m that the major-arc model predicts: the q > Q terms of the singular series (spikes the arcs leave out) "
      "plus the Dirichlet-kernel tails that spill past c/N. ‘resid’ = I_m − leak. ‘model err’ = I_M − I_M(model) "
      "from §4. σ_rp = P_m/√(N|m|) is a random-phase yardstick: the size of a sum of N|m| independent coherence "
      "cells of width 1/N carrying the mass P_m with random signs. sup_m|S| is the grid maximum refined by direct "
      "summation within ±1/M (the grid alone can miss a crest by ≤ 0.6% on a main lobe, ≤ 1.8% on a sidelobe); "
      "‘at’ names the rational whose model spike explains it.\n")
    A(_table(["Q", "c", "P_m/R2", "heur.", "heur. exact", "I_m/R2", "|I_m|/√P_m", "leak/R2", "resid/R2",
              "model err/R2", "σ_rp/R2", "sup_m|S|", "at a/q", "N·β"],
             [[str(r["Q"]), str(r["c"]), f"{r['P_m_over_R2']:.3f}", f"{r['pred_Pm_over_R2']:.3f}",
               f"{r['pred_exact_Pm_over_R2']:.3f}", f"{r['I_m_over_R2']:+.5f}", f"{abs(r['I_m']) / r['sqrtP_m']:.1f}",
               f"{r['leak'] / R2:+.5f}", f"{r['I_m_resid'] / R2:+.5f}", f"{r['model_err'] / R2:+.5f}",
               f"{r['rp_scale'] / R2:.4f}", _e(r["sup"], 3), f"{r['sup_a']}/{r['sup_q']}", f"{r['sup_Nbeta']:+.2f}"]
              for r in rows]))
    A(f"""
* **The L² wall.** P_m/R2 runs from {Pm_lo:.2f} to {Pm_hi:.2f} over the 18 grid pairs, and down to
  {Pm_lo_all:.2f} at the Qmax pairs (158, 20) and (223, 10) of §7: any estimate of the minor arcs by their magnitude
  alone, |I_m| ≤ P_m, is too weak by that factor. For Q ≤ 20 the exact-sum heuristic is within {heurx_small:.3f}
  (the asymptotic form within {heur_small:.2f}, its o(1) error at small Q). As the arcs widen (|M| from 3% to 27%)
  P_m falls increasingly below both, because wide arcs swallow part of the unresolved background as well as the
  spikes ({big['P_m_over_R2']:.2f} against {big['pred_exact_Pm_over_R2']:.2f} at (150, 20)); it is still nearly 3× the signal. Raising Q to (log N)^B would only subtract
  B·log log N = {math.log(math.log(N)):.2f}·B from log N − 2.33 = {math.log(N) - 2.33:.2f}.
* **What the minor arcs contribute.** |I_m|/R2 is 3–13% at Q ≤ 20 and {min(im_q50):.1e}–{max(im_q50):.1e} for
  Q ≥ 50, against a mass of {Pm_lo:.1f}–{Pm_hi:.1f}·R2. At small Q nearly all of it is the predictable leakage (the
  q > Q terms of the singular series and the kernel tails); |I_m|/√P_m is 20–73 there, well above log N, because
  that part is not random.
* **What is left, and what it is not.** Because I_m = R2 − I_M, the residual satisfies the identity
  resid = (R2 − SS2·N) − (I_M − I_M(model)) + O(10⁻⁶·R2). The first term is the global Hardy–Littlewood error,
  {(R2 - SS2 * N) / R2:+.4f}·R2 in every row; the second is the major-arc model error ({min(me_all):+.4f} to
  {max(me_all):+.4f}·R2). So the residual is the minor-arc share of the singular-series model's error, and at this
  N it cannot be separated from the model's error on the major arcs. Its size for Q ≥ 50
  ({min(res_q50):.1e}–{max(res_q50):.1e}·R2) is comparable to the random-phase yardstick σ_rp
  ({min(rp_q50):.1e}–{max(rp_q50):.1e}·R2): the data are consistent with square-root cancellation of the minor-arc
  mass, but they do not measure it, and nothing here certifies it for any other N.
""")

    # ------------------------------------------------------------ module 6
    A("## 6. Ternary control\n")
    A(_table(["N'", "factorisation", "R3 exact", "SS3 (product)", "SS3_150", "main = SS3·N'²/2", "R3/main"],
             [[f"{n:,}", "·".join(f"{p}^{e}" if e > 1 else str(p) for p, e in ctx["factor"][n].items()),
               _e(ctx["R3"][n], 6), f"{ctx['SS3'][n]:.6f}", f"{ctx['SS3_part'][n][150]:.6f}",
               _e(ctx["SS3"][n] * n * n / 2, 6), f"{ctx['R3'][n] / (ctx['SS3'][n] * n * n / 2):.6f}"]
              for n in NTER]))
    A("\nSS3_Q truncations (Q = 5, 10, 20, 50, 100, 150, 1000, 10^4): " + "; ".join(
        f"N' = {n:,}: " + ", ".join(f"{ctx['SS3_part'][n][q]:.5f}" for q in (5, 10, 20, 50, 100, 150, 1000, 10000))
        for n in NTER) + ".\n")
    A(f"B3 = sup_m|S| · ∫₀¹|S|² (the textbook minor-arc bound, ∫|S|² = ΣΛ² = {_e(L2, 4)}), against the binary "
      "P_m/R2 at the same (Q, c):\n")
    A(_table(["Q", "c", "sup_m|S|", "B3", f"B3/main N'={Nm:,}", f"B3/main N'={Np:,}", f"I3_M/main ({Nm:,})",
              f"I3_m/main ({Nm:,})", "binary P_m/R2"],
             [[str(r["Q"]), str(r["c"]), _e(r["sup"], 3), _e(r["B3"], 3), f"{r[f'B3_main_{Nm}']:.3f}",
               f"{r[f'B3_main_{Np}']:.3f}", f"{r[f'I3_M_{Nm}'] / r[f'main3_{Nm}']:.5f}",
               f"{r[f'I3_m_{Nm}'] / r[f'main3_{Nm}']:+.2e}", f"{r['P_m_over_R2']:.3f}"] for r in rows]))
    A(f"""
The bound closes when sup_m|S| < main/ΣΛ² = {thr[Nm]:,.0f} = ψ/{psi / thr[Nm]:.1f} for N' = {Nm:,}, and
{thr[Np]:,.0f} = ψ/{psi / thr[Np]:.1f} for N' = {Np:,}. Two kinds of stray crest on the minor arcs set sup_m|S|:

* **Excluded squarefree q with small φ(q)** — these need Q to grow. The table's sup sits at 1/6 (ψ/2) for Q = 5,
  9/14 (ψ/6) for Q = 10, 11/30 (ψ/8) for Q = 20, 49/66 (ψ/20) for Q = 50 with c ≥ 10, 61/102 (ψ/32) for Q = 100 with
  c ≥ 10, and {big['sup_a']}/{big['sup_q']} (ψ/48, measured {big['sup']:,.0f}) at (150, 20).
* **The first Dirichlet sidelobe just outside the α = 0 (and α = 1/2) arc**, |S| ≈ ψ/(π(c+½)) — this needs c to
  grow. Measured {sup5:,.0f} at c = 5 and {sup10:,.0f} at c = 10 (model {psi / (math.pi * 5.5):,.0f} and
  {psi / (math.pi * 10.5):,.0f}). No arc with q ≤ Q can cover it, so at c = 5 it keeps B3/main({Nm:,}) at
  {rmap[(150, 5)][f'B3_main_{Nm}']:.3f} for every Q ≥ 50; at c = 10 it ({sup10:,.0f}) is almost exactly as tall as the
  61/102 crest ({rmap[(100, 10)]['sup']:,.0f}), so B3/main stays at {rmap[(150, 10)][f'B3_main_{Nm}']:.3f}–{rmap[(100, 10)][f'B3_main_{Nm}']:.3f}
  for Q = 100, 150 (and 223 in §7).

At the largest grid pair (150, 20) — and at the Qmax pair (158, 20), which has the same sup — the bound closes with
room to spare: B3/main = {big[f'B3_main_{Nm}']:.3f} (N' = {Nm:,}) and {big[f'B3_main_{Np}']:.3f} (N' = {Np:,}), against
binary P_m/R2 = {big['P_m_over_R2']:.3f} in the same row. Either way the sup is set by near-rational α. In the
asymptotic theory the arcs have width ~Q/N, so Q and the distance to the arc edge grow together (Vinogradov: Q a
power of log N). The actual ternary minor-arc share |I3_m|/main is ≤ 3×10⁻³ for Q ≥ 20.
""")

    # ------------------------------------------------------------ module 7
    A("## 7. Scaling in N\n")
    A("M = smallest power of two ≥ 8N for each N (alias-free for R3, ≥ 8 samples per 1/N). Series: (Q, c) = (20, 10); "
      "(50, 10) where admissible (N ≥ 10^5); and Q = Qmax(N, c) = the largest Q with 2cQ² < N, for c = 10 and 20. "
      "The extra row N = 1,531,530 = 2·3²·5·7·11·13·17 is a probe at the largest SS2 for N ≤ 2·10^6 (5.4627, shared "
      "with 510,510 and 1,021,020). ‘resid’ and σ_rp as in §5.\n")
    A(_table(["N", "series", "Q", "c", "SS2(N)", "P_m/R2", "P_m/N", "I_m/R2", "resid/R2", "σ_rp/R2", "sup_m|S|", "at",
              "B3/main N−1", "B3/main N+1"],
             [[f"{r['N']:,}{'*' if r.get('probe') else ''}", r["series"], str(r["Q"]), str(r["c"]),
               f"{r['SS2']:.3f}", f"{r['P_m_over_R2']:.3f}", f"{r['P_m'] / r['N']:.3f}", f"{r['I_m_over_R2']:+.4f}",
               f"{r['I_m_resid'] / r['R2']:+.4f}", f"{r['rp_scale'] / r['R2']:.4f}",
               _e(r["sup"], 3), f"{r['sup_a']}/{r['sup_q']}", f"{_b3(r, -1):.3f}",
               f"{_b3(r, +1):.3f}"] for r in sorted(sc, key=lambda r: (r["N"], r["c"], r["Q"]))]))
    m_rows = []
    for m in ctx["scaling_meta"]:
        m_rows.append([f"{m['N']:,}", f"2^{int(math.log2(m['M']))}", f"{m['parseval_relerr']:.1e}",
                       f"{m['R2_grid_relerr']:.1e}", f"{max(m['R3_grid_relerr'].values()):.1e}",
                       f"{m['R2'] / (m['SS2'] * m['N']):.4f}",
                       " / ".join(f"{m['R3'][k] / (m['SS3'][k] * k * k / 2):.5f}" for k in sorted(m['R3']))])
    A("\nQuadrature certificates per N:\n")
    A(_table(["N", "M", "Parseval rel. err", "R2 grid rel. err", "R3 grid rel. err", "R2/(SS2·N)",
              "R3/main N−1 / N+1"], m_rows))
    rise_q = qmx[ns[-1]]["P_m"] / ns[-1] - qmx[ns[0]]["P_m"] / ns[0]
    rise_f = fix20[ns[-1]]["P_m"] / ns[-1] - fix20[ns[0]]["P_m"] / ns[0]
    dlog = math.log(ns[-1] / ns[0])
    A(f"""
* **Binary.** At fixed (20, 10), P_m/N rises from {fix20[ns[0]]['P_m'] / ns[0]:.2f} (N = 10^4) to
  {fix20[ns[-1]]['P_m'] / ns[-1]:.2f} (N = 2·10^6), roughly linearly in log N; P_m/R2 zigzags because SS2(N) is 3.52
  at N = 3·10^4, 3·10^5 and 1.76 elsewhere. With Q = Qmax(N) (Q = {qmx[ns[0]]['Q']} … {qmx[ns[-1]]['Q']}, arcs covering
  ≈ 30% of the circle) P_m/N still climbs, {qmx[ns[0]]['P_m'] / ns[0]:.2f} → {qmx[ns[-1]]['P_m'] / ns[-1]:.2f}: a rise
  of {rise_q:.1f} against {rise_f:.1f} at fixed (20, 10), while log N rises by {dlog:.1f}. So Q ≈ √(N/2c) removes
  about {1 - rise_q / rise_f:.0%} of the fixed-Q growth (slope {rise_q / dlog:.2f} per unit log N against
  {rise_f / dlog:.2f}; theory for Q ∝ √N: ½), more than ½ because arcs this wide also absorb background.
  P_m/R2 at Qmax is {', '.join(f"{qmx[n]['P_m_over_R2']:.2f}" for n in ns)} — it dips to
  {min(qmx[n]['P_m_over_R2'] for n in ns):.2f} < 1 at N = 3·10^4, where log N is small and SS2 is doubled. What has to be
  beaten is SS2(N), which never exceeds 5.46 for N ≤ 2·10^6 and grows at most like log log N, while P_m/N grows like
  log N, so the gap widens. The smooth probe shows the other side of the coin: at N = 1,531,530 (SS2 = 5.46)
  P_m/R2 = {', '.join(f"{r['P_m_over_R2']:.2f} ({r['series']})" for r in probe)}. In
  {len(pos_rows)} tabulated rows the minor-arc mass is below the major-arc integral, so the crude inequality
  R2 ≥ I_M − P_m is already positive: {', '.join(f"N = {r['N']:,} at ({r['Q']}, {r['c']}), I_M − P_m = {(r['I_M'] - r['P_m']) / r['R2']:+.2f}·R2" for r in pos_rows)}.
  All are N with an unusually large SS2. That says nothing about other N: typical N of the same size sit at
  P_m/R2 ≈ 2.6–2.9 at Qmax, and log N outruns log log N.
* **Ternary — this differs from the expected behaviour.** At fixed (20, 10), B3/main does **not** stay below 1:
  it is {', '.join(f"{fix20[n][f'B3_main_{n - 1}']:.2f}" for n in ns)} (N' = N−1) and creeps up with log N.
  The sup column shows why: with Q fixed, sup_m|S| is pinned at the q = 30 crest ψ/φ(30) ≈ N/8, a fixed fraction of
  N, so B3/main ≈ 2(log N − 1)/(8·SS3(N')) grows like log N, just like the binary ratio. What differs is the response
  to the arcs: grown as far as disjointness allows, the ternary ratio falls —
  {', '.join(f"{qmx20[n][f'B3_main_{n - 1}']:.2f}" for n in ns)} at (Qmax, 20) — because it scales like the tallest
  stray crest (1/φ(q*) for the first excluded q*, or 1/(π(c+½)) at the arc edge), while the binary ratio scales like
  log N − log Q. The asymptotic theory needs Q and c growing like powers of log N; the constraint 2cQ² < N leaves only
  Q ≤ √(N/2c), so at these N the decay cannot be displayed cleanly — only its mechanism.
""")

    # ------------------------------------------------------------ module 8
    A("## 8. Random control (Cramér model)\n")
    A(f"Main draw (seed 20261016): n ≥ 3 kept with probability 1/log n, weight log n; {R['count']:,} members, "
      f"ΣW = {R['psi']:,.0f}, ΣW² = {_e(R['L2'], 4)} (primes: {_e(L2, 4)}). Binary count R2_R(N) = {R['R2']:,.1f} "
      f"(exact; E R2_R = {R['E_R2']:,.0f}; grid integral agrees to {abs(R['R2_grid'] - R['R2']) / R['R2']:.1e}), "
      f"i.e. singular series ≈ 1, against primes R2 = {R2:,.0f} ≈ SS2·N with SS2 = {SS2:.4f}.\n")
    A(_table(["Q", "c", "I_M/R2_R", "P_m/R2_R", "I_m/R2_R", "sup_m|S_R|", "primes P_m/R2"],
             [[str(r["Q"]), str(r["c"]), f"{r['I_M'] / r['R2']:.4f}", f"{r['P_m_over_R2']:.3f}",
               f"{r['I_m_over_R2']:+.4f}", _e(r["sup"], 3), f"{rmap[(r['Q'], r['c'])]['P_m_over_R2']:.3f}"]
              for r in R["rows"]]))
    if seeds:
        A("\nFive more independent draws (seeds 1–5) plus the main draw:\n")
        A(_table(["seed", "R2_R/N", "P_m/R2_R (20, 10)", "I_m/R2_R (20, 10)", "P_m/R2_R (150, 20)", "I_m/R2_R (150, 20)"],
                 [[str(s_["seed"]), f"{s_['R2'] / N:.4f}", f"{s_[(20, 10)]['P_m_over_R2']:.3f}",
                   f"{s_[(20, 10)]['I_m_over_R2']:+.4f}", f"{s_[(150, 20)]['P_m_over_R2']:.3f}",
                   f"{s_[(150, 20)]['I_m_over_R2']:+.4f}"] for s_ in seeds]
                 + [["20261016", f"{R['R2'] / N:.4f}", f"{rr[(20, 10)]['P_m_over_R2']:.3f}",
                     f"{rr[(20, 10)]['I_m_over_R2']:+.4f}", f"{rr[(150, 20)]['P_m_over_R2']:.3f}",
                     f"{rr[(150, 20)]['I_m_over_R2']:+.4f}"]]))
    A(f"""
Across the six draws R2_R/N = {min(seed_R2):.3f}–{max(seed_R2):.3f} and P_m/R2_R = {min(seed_pm20):.2f}–{max(seed_pm20):.2f}
at (20, 10), {min(seed_pm150):.2f}–{max(seed_pm150):.2f} at (150, 20). The random set has the same L² wall, larger in
fact, because its energy ΣW² ≈ N log N is the same while nothing beyond the α = 0 arc is phase-locked and its R2 is
smaller (≈ N). So the wall belongs to the density 1/log n and to the method, not to arithmetic. The random model
escapes it — R2_R ≈ N ± O(√N log N) — by independence: its minor-arc phases really are random. The primes are
deterministic and offer no such independence; that is the gap. (sup_m|S_R| is always at the edge of the α = 0 arc:
there are no other crests.)
""")

    # ------------------------------------------------------------ module 9
    A("## 9. Goldbach's comet\n")
    A(f"G(N) for all even N ≤ {COMET_MAX:,} by one FFT autoconvolution of the prime indicator (spot checks against "
      f"direct counting: {com['spot_ok']}). G(N)/prediction does not scatter around 1 with no trend: its mean rises "
      "steadily toward 1 while the scatter shrinks — binned means "
      + ", ".join(f"{m:.3f} ({lo:,.0f}–{min(hi, 2e5):,.0f})" for lo, hi, m in bins)
      + f", and 5402 / {ss.SS2_product(10**6) * 0.5 * ss.hl_integral([10**6])[0]:.1f} = {g6:.4f} at N = 10^6. "
      "A slow approach to the asymptotic prediction, observed and not explained here.\n")

    A(ANALOGY)

    # ------------------------------------------------------------ deviations
    A("## Where the results differ from what the brief expected\n")
    A(f"""1. **Ternary at fixed (Q, c) = (20, 10) does not stay below 1, nor flat** (§7): B3/main is
   {', '.join(f"{fix20[n][f'B3_main_{n - 1}']:.2f}" for n in ns)} for N' = N − 1 and grows with log N, because
   sup_m|S| is pinned at the left-out 1/30 crest ≈ N/8. It falls below 1 once the arcs grow with N (Qmax series)
   — and only if c is large enough: at c = 5 the arc-edge sidelobe keeps B3/main(10^6 − 1) at 1.007 for every Q.
2. **P_m/R2 falls below the 3.5–7 band** at (150, 10) = {rmap[(150, 10)]['P_m_over_R2']:.2f}, (150, 20) =
   {big['P_m_over_R2']:.2f} and at the N = 10^6 Qmax pairs (158, 20), (223, 10) =
   {', '.join(f"{r['P_m_over_R2']:.2f}" for r in sorted(n6, key=lambda r: -r['c']))}: these arcs cover 14–30% of the
   circle and absorb background as well as spikes. It is still > 2.6.
3. **|I_m|/R2 ≈ 10⁻² only for Q ≥ 50.** At Q ≤ 20 it is 3–13%, and |I_m|/√P_m is 20–73 (also 15–17 at c = 5 with
   Q = 100, 150), above log N = 13.8: almost all of it is the deterministic leakage of left-out spikes and kernel
   tails (§5). What remains after removing it is the singular-series model's own error, not a clean measurement
   of random cancellation.
4. **Binary P_m/R2 is not > 1 for every tabulated N and (Q, c)**: it is {qmx[30000]['P_m_over_R2']:.2f} at N = 3·10^4
   with Q = Qmax = {qmx[30000]['Q']}, and {min(r['P_m_over_R2'] for r in probe):.2f}–{max(r['P_m_over_R2'] for r in probe if r['Q'] > 50):.2f} at the smooth
   probe N = 1,531,530 with Qmax — both N with an unusually large singular series. For N = 10^6 every tabulated
   pair gives P_m/R2 ≥ {Pm_lo_all:.2f}.
5. **B3/main at (150, 20)** is {big[f'B3_main_{Nm}']:.2f} for N' = 10^6 − 1 (inside the expected 0.3–0.6) but
   {big[f'B3_main_{Np}']:.2f} for N' = 10^6 + 1, whose singular series is larger (2.30 vs 1.46).
6. **G(N)/prediction** sits below 1 with a slowly rising mean ({bins[0][2]:.3f} at 10³–5·10³ up to {bins[-1][2]:.3f}
   at 1.5–2·10^5, {g6:.4f} at 10^6) rather than scattering around 1 with no trend.
7. **Memory**: peak RSS {ctx.get('peak_rss_mb', float('nan')):,.0f} MB (the N = 2·10^6 and 1,531,530 runs use M = 2^24).
""")

    # ------------------------------------------------------------ captions
    A("## Figure captions\n")
    caps = figure_captions(ctx, locals())
    for name, cap in caps:
        A(f"**{name}** — {cap}\n")

    # ------------------------------------------------------------ does / does not
    A("## What this does and does not show\n")
    A(f"""**Does show (for N = 10^6 and the listed N, in floating point):**

1. The circle-method identities hold numerically to machine precision: Parseval, R2 = ∫S²e(−Nα) and R3 = ∫S³e(−N'α)
   on a uniform grid with M > 2N' (relative errors ≤ 1e-15), with the FFT sign pinned by a check that fails loudly.
2. The major-arc picture is quantitatively right at this N: spike heights |μ(q)|ψ/φ(q) to ≤ {sf_dev:.1%} for q ≤ 30,
   the Dirichlet-kernel shape and phase through ~20 sidelobes, and the major-arc integral within a few percent of
   R2 for Q ≥ 20, c ≥ 10 (≤ 1.5% for Q ≥ 50). The singular series predicts R2(10^6) to {abs(R2 / (SS2 * N) - 1):.2%}
   and R3(N') to ≤ 7×10⁻⁵.
3. The binary obstruction, measured: at N = 10^6 the minor arcs carry {Pm_lo_all:.1f}–{Pm_hi:.1f}× the signal in L²
   mass for every tabulated (Q, c), while their signed contribution is {min(im_q50):.0e}–{max(im_q50):.0e} of R2 for
   Q ≥ 50 and is accounted for, to within the major-arc model's own error (≤ {max(abs(x) for x in me_all):.1%} of R2),
   by the predictable leakage. A magnitude estimate cannot see that cancellation. P_m/N grows like log N, while
   SS2(N), which it must stay below for the L² bound to close, grows at most like log log N.
4. The ternary contrast: with S³ one factor can be bounded by its sup, and at (Q, c) = (150, 20) the bound
   sup_m|S|·∫|S|² is {big[f'B3_main_{Np}']:.2f}–{big[f'B3_main_{Nm}']:.2f} of the main term. It closes because Q absorbs the
   tall excluded crests and c pushes the arc-edge sidelobes down; at fixed small Q, or at c = 5 for N' = 10^6 − 1,
   it does not.
5. The Cramér model (six draws) has the same L² wall but no Farey forest and singular series 1: the wall comes from
   density and method, the forest from arithmetic.

**Does not show:**

1. Anything asymptotic. N = 10^6 is pre-asymptotic: the textbook arc width Q/N with Q = (log N)^B (B ≥ 2) would make
   2Q³ > N, disjointness caps Q at √(N/2c), and all results are tabulated over (Q, c) rather than at one canonical
   choice. No trend has been fitted to minor-arc data; the heuristic curves contain no fitted parameters.
2. That the minor-arc cancellation persists, or how large it is beyond the model error. The residual after leakage
   is identically the Hardy–Littlewood error minus the major-arc model error; it is ≤ {max(res_q50):.1e}·R2 at
   N = 10^6 (Q ≥ 50) and up to {small_resid:.1e}·R2 at N ≤ 3·10^5. Nothing here bounds it for any other N, and that
   is precisely what would be needed.
3. Rigorous bounds. sup_m|S| is a grid maximum refined locally by direct summation, not a certified supremum;
   arithmetic is IEEE double precision, not interval arithmetic.
4. A proof strategy, binary or ternary. The ternary closure is shown for two specific N'; the actual theorem needs
   sup_m|S| ≪ N/(log N)^A uniformly (Vinogradov, with Siegel–Walfisz on the major arcs), which no finite computation
   displays.
5. An explanation of the 1–5% deficit of G(N) against the Hardy–Littlewood integral at N ≤ 2·10^5, shrinking with N.
""")
    path.write_text("\n".join(t))


def figure_captions(ctx: dict, v: dict) -> list[tuple[str, str]]:
    N, psi, L2, R2 = ctx["N"], ctx["psi"], ctx["L2"], ctx["R2"]
    m4, W, big = ctx["mod4"], ctx["window"], v["big"]
    ph = v["ph"]
    R = ctx["rand"]
    rr = v["rr"]
    fr = sorted(W["fractions"], key=lambda f: f["alpha"])
    sea_rms = math.sqrt(v["win_mean"])
    sea_rayleigh = math.sqrt(math.log(2) * v["win_mean"])
    crests = [f["height"] / sea_rms for f in fr if f["height"] > 0] or [0.0]
    locked = [f for f in fr if f.get("mu") and "S_exact" in f]
    lock_ok = sum(np.sign(f["S_exact"].real) == np.sign(f["mu"]) for f in locked)
    pk = W.get("peaks", [])
    on = [p for p in pk if p.get("near") and p["N_dist"] < 1.0]
    if on:
        pk_txt = ("yet the highest peaks sit on them (" + "; ".join(
            f"{p['ratio']:.1f}× the mean at {p['near'][0]}/{p['near'][1]}" for p in on) + ")")
    else:
        pk_txt = ("the two highest peaks (" + " and ".join(f"{p['ratio']:.1f}×" for p in pk)
                  + " the mean) are speckle, ≥ " + f"{min(p['N_dist'] for p in pk if p.get('N_dist') is not None):.0f}"
                  + "/N from any of them")
    st = ctx["stair"]
    com = ctx["comet"]
    zg = ctx["zoom_golden"]
    gt = v["gold_top"]
    flips = v["early_flips"]
    return [
        ("Fig A (fig_A_farey_forest.png)",
         "Each prime power n ≤ 10^6 is an oscillator with amplitude Λ(n) and phase 2πnα; |S(α)| is the length of "
         "their resultant. At a rational a/q every phase depends only on n mod q, so the oscillators lock into the "
         "φ(q) reduced channels, each carrying ≈ ψ/φ(q), and the channels add up to the Ramanujan sum: a spike of "
         "height |μ(q)|ψ/φ(q) (orange ticks), which is where the measured crests (blue) stand. For non-squarefree q "
         "(red ×) the locked channel phasors sum to zero and only a residue of 10²–10³ survives, below the sea. "
         f"Between the spikes the bin medians (≈ {v['sea_med']:,.0f}) are close to the Rayleigh median "
         f"√(ln2·⟨|S|²⟩) = {sea_rayleigh:,.0f} of a sea with mean energy ⟨|S|²⟩ ≈ {v['win_mean']:.1e} (Fig E window), "
         f"about half of ΣΛ² = {L2:.2e}: Parseval's energy is shared between forest and sea, and around 1/2, 1/4, "
         "3/4 the forest thins because a small-q rational has no low-q neighbours nearby (Farey gaps)."),
        ("Fig B (fig_B_spike_heights.png)",
         "(a) S(a/q) for every reduced a/q with q ≤ 30, summed exactly over residue classes, "
         "S(a/q) = Σ_r e(ar/q)·ψ(N;q,r). For squarefree q the length matches |μ(q)|ψ/φ(q) to within "
         f"{v['sf_dev']:.1%} (worst case {v['sf_dev_q']['a']}/{v['sf_dev_q']['q']}); non-squarefree q sit at ≤ "
         f"{v['nsf_max']:.1e}·ψ, below even the random-walk level √ΣΛ²/ψ = {math.sqrt(L2) / psi:.1e}, because their "
         "locked channel phasors cancel. (b) The resultant also points the predicted way: φ(q)·S(a/q)/ψ lies on the "
         "real axis at μ(q), i.e. the φ(q) channels add to exactly one channel, forwards or backwards according to "
         "the parity of the number of prime factors of q. The residue at 1/4 is literally the mod-4 race: "
         f"S(1/4) = [ψ(N;4,0) − ψ(N;4,2)] + i[ψ(N;4,1) − ψ(N;4,3)] = {m4['S_quarter'].real:.1f} "
         f"{m4['S_quarter'].imag:+.1f}i."),
        ("Fig C (fig_C_zooms.png)",
         "Zooms of half-width 20/N around 1/3, 1/6 and 1/4 (direct summation, 60 samples per 1/N). Detuning by β "
         "makes oscillator n precess at rate nβ, so a locked channel's phases fan out over 2πNβ; the collective "
         "phasor decoheres along the Dirichlet kernel v(β), whose zeros at Nβ = ±1, ±2, … are the detunings at which "
         "the fan closes into whole turns. The model (μ(q)/φ(q))·v(β) (orange) follows the measured modulus through "
         "~20 sidelobes on each side (at ±20/N they are still 2–3× the sea level √ΣΛ²), and in the frame "
         "co-rotating with the carrier e((N+1)β/2) the phase flips between 0 and π at every zero, starting at π for "
         "1/3 (μ = −1) and at 0 for 1/6 (μ = +1). At 1/4 the model is identically zero; what is left is sea plus the "
         "mod-4 residue, and its phase wanders."),
        ("Fig D (fig_D_golden_zoom.png)",
         "Same width and scale, centred on frac(golden ratio) = 0.6180339…, the number hardest to approximate by "
         f"rationals. Even here the convergents and their neighbours with q ≈ √N carry μ(q)-locked crests ψ/φ(q) "
         f"as tall as the sea: {len(v['gold_above'])} stand above √ΣΛ² = {math.sqrt(L2):,.0f}, the tallest "
         f"{gt['a']}/{gt['q']} at {gt['height']:,.0f}. The window maximum |S| = {v['gold_max']:,.0f} sits on that "
         f"crest, where exactly S({gt['a']}/{gt['q']}) = {gt['S_exact'].real:+,.0f} {gt['S_exact'].imag:+,.0f}i points "
         f"{'backwards' if gt['mu'] < 0 else 'forwards'} as μ({gt['q']}) = {gt['mu']:+d} predicts, and the exact sums "
         f"at the other sizeable crests point along μ(q) too ({v['gold_sign_ok']} of {len(v['gold_locked'])}). Between "
         "the crests the co-rotating phase drifts freely: at N = 10^6 even the deepest minor arc is weak phase-locking "
         "at the level of the sea, not pure noise."),
        ("Fig P (fig_P_phasors.png)",
         "Λ-weighted histograms of the phases 2πnα (360 bins) with the resultant S(α) drawn to the same radial "
         "scale. At 1/3 the primes fill two channels (n ≡ ±1 mod 3) of ψ/2 each at ±120°, and at 1/6 two channels "
         f"at ±60°: in both cases the resultant has length ψ/2 ({ph['1/3'].real:,.0f} and {ph['1/6'].real:+,.0f}), "
         "pointing backwards at 1/3 (μ(3) = −1) and forwards at 1/6 (μ(6) = +1). At 1/4 the channels meet head-on at "
         f"90° and 270° and leave {abs(ph['1/4']):,.0f} ≈ {abs(ph['1/4']) / psi:.0e}·ψ, the imbalance of the mod-4 race. "
         "1/5 and 1/30 have 4 and 8 channels whose sum is exactly one channel reversed "
         f"({ph['1/5'].real:,.0f} ≈ −ψ/4, {ph['1/30'].real:,.0f} ≈ −ψ/8). At the irrationals the bins form a flat ring "
         f"and the resultant is random-walk sized ({abs(ph['frac(golden ratio)']):,.0f} and "
         f"{abs(ph['frac(sqrt 2)']):,.0f} against √ΣΛ² = {math.sqrt(L2):,.0f})."),
        ("Fig S (fig_S_singular_series.png)",
         "SS2_Q(10^6) against Q. Each squarefree q adds the coherent beat of two locked channel sets, "
         "μ(q)²c_q(N)/φ(q)², a product over its primes of c_p(N)/(p−1)²: +1 for p = 2 and +1/4 for p = 5 (they divide "
         "N), −1/(p−1)² for every other p — so the kick is positive when q has an even number of prime factors "
         "outside {2, 5} (q = 21 adds +1/144) and negative when it has an odd number. The partial sums swing — 1, 2, "
         f"1.75, 2, 1.75, 1.72, … — and settle: SS2_20 = {ctx['SS2_part'][20]:.4f}, SS2_150 = {ctx['SS2_part'][150]:.5f}, "
         f"and the sum to Q = 10^6 equals 2C2·4/3 (C2 to 20 digits) = {ctx['SS2']:.7f} to "
         f"{abs(ctx['SS2_sum_1e6'] - ctx['SS2']):.0e}. Panel (b): the signed error changes sign only at Q = "
         f"{', '.join(str(q) for q in flips)}, then stays negative and shrinks from ~10⁻³ to ~10⁻⁸; beyond "
         f"Q ≈ {v['late_flip'] or 10**4:,} the ~10⁻⁸ tail oscillates about 0."),
        ("Fig E (fig_E_minor_stretch.png)",
         f"A stretch of width 300/N at α₀ = {W['x0']:.6f}, inside the minor arcs of every (Q, c) tabulated at "
         "N = 10^6. (a) The mass density |S|² is mostly speckle — coherence cells of width ~1/N, window mean "
         f"{v['win_mean']:.2e} — and (b) the integrand Re[S²e(−Nα)] has the same envelope, but its phase (the beat of "
         "S² against the reference e(Nα)) turns once per ~1/N, so neighbouring cells enter with effectively random "
         f"signs. (c) Mass piles up linearly to {W['cum_mass'][-1]:,.0f} over the window while the signed integral "
         f"wanders and ends at {W['cum'][-1]:+,.0f}: a random walk of ~300 steps. The fractions with q ≤ 300 inside "
         f"the window ({', '.join(f'{f_['a']}/{f_['q']}' for f_ in fr)}) have model crests only "
         f"{(f'{min(crests):.1f}' if max(crests) - min(crests) < 0.05 else f'{min(crests):.1f}–{max(crests):.1f}')}× the sea's rms √⟨|S|²⟩ = {sea_rms:,.0f}, and {pk_txt}; still, the exact "
         f"S(a/q) point along μ(q) in {lock_ok} of {len(locked)} cases — weakly locked crests that the major arcs "
         "leave on the minor arcs."),
        ("Fig F (fig_F_staircase.png)",
         "The running integral F(α) = ∫₀^α S²e(−Nt)dt. Each major-arc spike is a burst of phase-locked beat "
         "between S² and e(Nα) that adds the phasor (N μ(q)²/φ(q)²)·e(−aN/q): real jump N μ(q)²cos(2πaN/q)/φ(q)², "
         "so +N at 1/2, −N/8 at 1/3, 2/3, 1/6 and 5/6, +N/16 at each of 1/5…4/5, nothing at 1/4, 3/4, and smaller "
         f"steps at q = 7, 10, … by the same rule. Measured jumps over ±20/N agree to {v['jmax']:.1e}·N (the kernel "
         f"tail beyond 20/N), and between jumps the plateaus wander at the 10⁻³N level (inset: a {v['plateau_range']:.1e}·N "
         "range over α ∈ [0.777, 0.782]) — small spikes plus minor-arc sea. The imaginary jumps of a/q and (q−a)/q "
         "cancel in pairs, so Im F is symmetric about 1/2 and returns to "
         f"{st['final'].imag:.1e}; the real part ends at R2(N)/N = {st['final'].real / N:.6f}, exactly R2."),
        ("Fig T (fig_T_binary_vs_ternary.png)",
         "(a) Magnitude-only minor-arc bounds divided by the main term, for every grid pair (Q, c) at N = 10^6. "
         f"Binary: P_m/R2 = {v['Pm_lo']:.1f}–{v['Pm_hi']:.1f}, always above 1 — the decoherent energy left on the minor "
         "arcs exceeds the whole signal. Ternary: sup_m|S|·∫|S|² ÷ SS3·N'²/2 drops below 1 only when both kinds of "
         f"stray crest are short enough — the excluded squarefree q (needs Q) and the arc-edge sidelobe ψ/(π(c+½)) "
         f"(needs c; at c = 5 it keeps N' = 10^6−1 at 1.007 for every Q) — reaching {big['B3_main_999999']:.2f} "
         f"(N' = 10^6−1, solid) and {big['B3_main_1000001']:.2f} (N' = 10^6+1, hollow) at (150, 20). (b) What the minor "
         "arcs really contribute: binary |I_m|/R2 of 3–13% at Q ≤ 20, almost all predictable leakage of left-out "
         "spikes; what is left after removing it is the model's own error, comparable in size to the random-phase "
         "yardstick P_m/√(N|m|) (grey) — the phases cancel to a fraction of a percent of a mass several times larger "
         "than the signal."),
        ("Fig G (fig_G_scaling.png)",
         "N from 10^4 to 2·10^6, four (Q, c) series, no fits; stars mark the smooth-N probe N = 1,531,530 "
         "(SS2 = 5.46). (a, b) At fixed (Q, c) the binary minor-arc energy P_m/N climbs roughly linearly in log N, "
         "and even Q = Qmax(N) (the most arcs 2cQ² < N allows) only lowers the slope; the singular series it has to "
         "beat (black: its largest value at each size) barely moves. (c) The ternary ratio at fixed (Q, c) also "
         "creeps up with log N, because the sup is pinned by the left-out crest at 1/30 (≈ N/8) — contrary to the "
         "brief's expectation that it stays flat — while with Q = Qmax it falls, from 0.9–1.4 at 10^4 to 0.2–0.5 at "
         "2·10^6. The ternary bound only needs the strongest stray resultant on the minor arcs to be small, which "
         "more and wider arcs can arrange; the binary bound needs their total energy to be small, which they cannot."),
        ("Fig H (fig_H_random_model.png)",
         "(a) Max |S| per bin for the primes (blue) and, drawn on top, a Cramér random set with the same density "
         "(orange): randomising destroys every phase-locking except the trivial one at α = 0, so only the primes' "
         "forest sticks out above the random band. (b) The median per bin: the random sea sits uniformly higher, "
         "because it keeps all of its energy ΣW² instead of lending about half to spikes. (c) At the rationals "
         f"themselves the random resultant is of random-walk size. Its binary count is R2_R = {R['R2']:,.0f} ≈ N "
         f"(singular series 1) against {R2:,.0f} ≈ 1.76N for the primes, yet its minor-arc energy is "
         f"P_m/R2_R = {rr[(150, 20)]['P_m_over_R2']:.1f}–{rr[(5, 5)]['P_m_over_R2']:.1f} (table in §8): the same L² wall, "
         "higher — the density's wall, not arithmetic's, which the random set gets past by an independence the "
         "primes do not provide."),
        ("Fig K (fig_K_comet.png)",
         f"G(N) for every even N ≤ 2·10^5. (a) The comet splits into bands by N mod 3: if 3 | N, the two channels "
         "n ≡ 1, 2 (mod 3) pair up with each other in both orders, doubling the count (factor (3−1)/(3−2) = 2), while "
         "N ≡ 1 and N ≡ 2 share the lower bands. (b) Coloured by log SS2(N), each band is one value of the singular "
         "series — one pattern of phase-locked channels mod the small primes dividing N — and the Hardy–Littlewood "
         "prediction (orange) runs through each band. (c) Observed over predicted is not trend-free: its binned mean "
         f"(orange) rises from {v['bins'][0][2]:.3f} at N = 10³–5·10³ to {v['bins'][-1][2]:.3f} at 1.5–2·10^5 "
         f"({v['g6']:.4f} at 10^6) while the scatter shrinks — a slow approach to the prediction."),
    ]


# =============================================================== entry

def compute() -> dict:
    ctx: dict = {}
    module_main(ctx)
    module_scaling(ctx)
    module_comet(ctx)
    ctx["runtime_s"] = time.time() - T0
    import resource
    ctx["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    with open(CACHE / "ctx.pkl", "wb") as f:
        pickle.dump(ctx, f, protocol=pickle.HIGHEST_PROTOCOL)
    return ctx


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", action="store_true", help="reuse cache/ctx.pkl")
    args = ap.parse_args()
    CACHE.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)
    if args.figures and (CACHE / "ctx.pkl").exists():
        with open(CACHE / "ctx.pkl", "rb") as f:
            ctx = pickle.load(f)
    else:
        ctx = compute()
    import figures
    figures.make_all(ctx, FIGS)
    log("figures written")
    write_results(ctx, HERE / "results.md")
    log("results.md written")


if __name__ == "__main__":
    main()
