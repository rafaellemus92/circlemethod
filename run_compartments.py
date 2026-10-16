"""The body-water lens (B1-B3): the energy of S as compartments measured by indicator dilution.

    python run_compartments.py            # ~10 s, peak RSS ~0.8 GB; needs run_all.py's S caches
    python run_compartments.py --report   # redraw figures/fig_W1_compartments.png and compartments.md from the pickle
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

import compartments as cp  # noqa: E402
from expsum import S_grid  # noqa: E402
from fractals import fg_factor, q_star_fg  # noqa: E402
from sieve import factorize, mobius_phi, von_mangoldt  # noqa: E402

CACHE = HERE / "cache"
FIGS = HERE / "figures"
XS = [(10**4, 2**17), (3 * 10**4, 2**18), (10**5, 2**20), (3 * 10**5, 2**22), (10**6, 2**23), (2 * 10**6, 2**24)]
X0, M0 = 10**6, 2**23
T0 = time.time()


def log(msg: str) -> None:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[{time.time() - T0:7.1f}s, peak {peak:6.0f} MB] {msg}", flush=True)


def certify(name: str, ok: bool, detail: str, res: dict) -> None:
    res["certificates"].append((name, bool(ok), detail))
    if not ok:
        raise AssertionError(f"certificate failed: {name}: {detail}")


def local_sea(X: int, M: int, qs, lam: np.ndarray, psi: float) -> dict:
    """Sea level around the rationals a/q (per unit alpha), for squarefree q:
    'centre'   mean over reduced a of |S(a/q) - mu(q) psi_c/phi(q)|^2 at the spike centres (residue sums; measured);
    'centre_red' the same with the prime powers p^k (p | q) left out, i.e. reduced classes only;
    'divisor'  'centre_red' from the residue-class variances of the divisors d | q, taken over n coprime to q
               (cp.divisor_sea; an exact identity);
    'fg'       the divisor form with V(X; d) replaced by X * FG(d) (Friedlander-Goldston; 0 where negative; no fit).
    For prime q the divisor form is V(X; q) q/(q - 1), i.e. about X FG(q); for composite q it is larger.
    Also the window measure of each q, so that the sea a window covers is sum_q measure_q * level_q."""
    n = np.flatnonzero(lam)
    w = lam[n]
    out = {"fg": {}, "centre": {}, "centre_red": {}, "divisor": {}, "measure": {}}
    Vm, Vf = {1: 0.0, 2: 0.0}, {1: 0.0, 2: 0.0}

    def variances(d):
        if d not in Vm:
            res = np.bincount(n % d, weights=w, minlength=d)
            red = np.gcd(np.arange(d), d) == 1
            v = res[red]
            Vm[d] = float(np.sum((v - v.sum() / red.sum()) ** 2))
            Vf[d] = max(fg_factor(d, list(factorize(d))), 0.0) * X
        return Vm[d], Vf[d]

    for q in qs:
        ps = list(factorize(q))
        divs = [math.prod(p for i, p in enumerate(ps) if mask >> i & 1) for mask in range(1, 1 << len(ps))]
        for d in divs:
            variances(d)
        cop = np.gcd(n, q) == 1                     # the same variances over n coprime to q: the identity is then exact
        nq, wq = n[cop], w[cop]
        Vq = {1: 0.0}
        for d in divs:
            r_d = np.bincount(nq % d, weights=wq, minlength=d)
            red_d = np.gcd(np.arange(d), d) == 1
            Vq[d] = float(np.sum((r_d[red_d] - wq.sum() / red_d.sum()) ** 2)) if d > 2 else 0.0
        res = np.bincount(n % q, weights=w, minlength=q)
        Sa = q * np.fft.ifft(res)
        red = np.gcd(np.arange(q), q) == 1
        mu_q = (-1) ** len(ps)
        coh = mu_q * float(res[red].sum()) / int(red.sum())
        out["centre"][q] = float(np.mean(np.abs(Sa[red] - coh) ** 2))
        res_red = np.where(red, res, 0.0)                                   # reduced classes only: the identity is exact
        out["centre_red"][q] = float(np.mean(np.abs((q * np.fft.ifft(res_red))[red] - coh) ** 2))
        out["divisor"][q] = cp.divisor_sea(q, Vq, ps)
        out["fg"][q] = cp.divisor_sea(q, Vf, ps)
        a, k, ok, cap = cp.fraction_windows(M, X, q)
        out["measure"][q] = float(ok.sum()) / M
    return out


def budget(a2: np.ndarray, X: int, M: int, model: dict, Q: int, sea: dict | None = None, atom_qs=(1, 2)) -> dict:
    """Atoms (the denominators atom_qs) and forest (the other q <= Q with a nonzero model spike) read by their windows, then
    solved.  With `sea` (from local_sea), the forest is also re-read with the local sea level instead of a uniform one."""
    T = float(a2.sum()) / M
    tiny = 1e-9 * max(float(np.sum(model[q])) for q in model)          # rounding residue of an exact zero
    qa = [q for q in atom_qs if float(np.sum(model[q])) > tiny]
    qf = [q for q in range(3, Q + 1) if q not in atom_qs and float(np.sum(model[q])) > tiny]
    mA, kA = cp.compartment_indicator(M, X, qa, {q: model[q] for q in qa})
    mF, kF = cp.compartment_indicator(M, X, qf, {q: model[q] for q in qf})
    mF &= ~mA
    if not qa:
        kA = cp.SINC2_HALF
    if not qf:
        kF = cp.SINC2_HALF
    E_A, E_F = float(a2[mA].sum()) / M, float(a2[mF].sum()) / M
    m_A, m_F = float(mA.sum()) / M, float(mF.sum()) / M
    sol = cp.solve_compartments(E_A, E_F, T, m_A, m_F, kA, kF)
    pred_A = sum(float(np.sum(model[q])) for q in qa)
    pred_F = sum(float(np.sum(model[q])) for q in qf)
    naive = cp.solve_compartments(E_A, E_F, T, m_A, m_F, cp.SINC2_HALF, cp.SINC2_HALF)
    out = {"T": T, "E_A": E_A, "E_F": E_F, "m_A": m_A, "m_F": m_F, "k_A": kA, "k_F": kF, "n_qF": len(qf), **sol,
           "pred_A": pred_A, "pred_F": pred_F, "naive_A": naive["A"], "naive_F": naive["F"]}
    out["qa"] = qa
    if sea is not None:
        for kind in ("fg", "centre"):
            S_F = sum(sea["measure"][q] * sea[kind][q] for q in qf)
            F = (E_F - S_F) / kF
            A, rho = np.linalg.solve(np.array([[kA, m_A], [1.0 - kA, 1.0 - m_A - m_F]]),
                                     [E_A, T - E_A - E_F - (1.0 - kF) * F])
            out[f"local_{kind}"] = {"A": float(A), "F": float(F), "rho": float(rho), "sea": float(T - A - F), "S_F": S_F}
        out["sea_uniform_in_forest"] = sol["rho"] * m_F
    return out


def per_q_check(a2: np.ndarray, X: int, M: int, model: dict, qs) -> list[dict]:
    """Window energy of each denominator against capture x model (+ nothing for the sea: small q only)."""
    out = []
    for q in qs:
        a, k, ok, cap = cp.fraction_windows(M, X, q)
        E = sum(float(a2[k[i][ok[i]] % M].sum()) for i in range(len(a))) / M
        pred = float(np.sum(cap * model[q]))
        out.append({"q": q, "E": E, "pred": pred, "ratio": E / pred, "cap_min": float(cap.min()), "cap_max": float(cap.max())})
    return out


def stage_development(res, mu, phi):
    rows = []
    for Xs, Ms in XS:
        lam = von_mangoldt(Xs)
        psi, L2 = float(lam.sum()), float(np.sum(lam**2))
        S = np.load(CACHE / f"S_primes_N{Xs}_M{Ms}.npz")["S"]
        a2 = S.real**2 + S.imag**2
        del S
        certify(f"Parseval, primes, X = {Xs}", abs(a2.sum() / Ms - L2) <= 1e-9 * L2, "", res)
        qs = q_star_fg(Xs, psi)
        Q = int(qs)
        model = cp.model_energies(range(1, Q + 1), Xs, psi, 1, {0})
        qf = [q for q in range(3, Q + 1) if mu[q] != 0]
        sea = local_sea(Xs, Ms, qf, lam, psi)
        dev_id = max(abs(sea["divisor"][q] / sea["centre_red"][q] - 1) for q in qf)
        dev_pp = max(abs(sea["centre_red"][q] / sea["centre"][q] - 1) for q in qf)
        b = budget(a2, Xs, Ms, model, Q, sea)
        bands = []
        for lo, hi in ((3, 31), (31, 101), (101, 256), (256, Q + 1)):
            qq = [q for q in qf if lo <= q < hi]
            if not qq:
                continue
            meas = sum(sea["measure"][q] for q in qq)
            E = cap_model = 0.0
            for q in qq:
                a, k, ok, cap = cp.fraction_windows(Ms, Xs, q)
                E += sum(float(a2[k[i][ok[i]] % Ms].sum()) for i in range(len(a))) / Ms
                cap_model += float(np.sum(cap * model[q]))
            bands.append({"q": (min(qq), max(qq)), "n_q": len(qq), "E": E, "cap_model": cap_model, "measure": meas,
                          "eff": (E - cap_model) / meas, "fg": sum(sea["measure"][q] * sea["fg"][q] for q in qq) / meas,
                          "centre": sum(sea["measure"][q] * sea["centre"][q] for q in qq) / meas})
        pr = cp.prime_predictions(Xs, psi, L2, mu, phi, Q)
        row = {"X": Xs, "M": Ms, "psi": psi, "L2": L2, "q_star": qs, "Q": Q, **b, "pred": pr, "sea_bands": bands,
               "asym": cp.asymptotic_shares(Xs, qs, pr["L_N"])}
        row["divisor_identity_max_dev"], row["prime_power_max_dev"] = dev_id, dev_pp
        certify(f"centre sea (reduced classes) = divisor sum of residue variances, X = {Xs}", dev_id <= 1e-9,
                f"max |ratio - 1| = {dev_id:.1e}; prime powers p^k, p | q, move it by up to {dev_pp:.3f}", res)
        if Xs == X0:
            row["per_q"] = per_q_check(a2, Xs, Ms, model, [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 30])
            worst = max(abs(r["ratio"] - 1) for r in row["per_q"])
            certify("window energy = capture x model, q <= 30 at X = 10^6", worst <= 0.01, f"max |ratio - 1| = {worst:.4f}", res)
        rows.append(row)
        del a2
        log(f"B2 X = {Xs}: atoms {b['A'] / L2:.4f}, forest {b['F'] / L2:.4f}, sea {b['sea'] / L2:.4f}")
    res["development"] = rows
    res["limits"] = [{"log10X": e, **dict(zip(("atoms", "forest", "sea"), cp.leading_order_shares(e * math.log(10))))}
                     for e in (6, 12, 20, 40, 100, 1000)]
    lx = np.linspace(4, 40, 181)
    res["limits_curve"] = {"log10X": lx, "shares": np.array([cp.leading_order_shares(e * math.log(10)) for e in lx])}


def odd_cramer_weights(N: int, seed: int = 20261017) -> np.ndarray:
    """Cramer model restricted to odd n >= 3: member with probability 2/log n (capped at 1), weight log n, so E[w(n)] = 2
    on odd n and 0 on even n (the same total weight as the primes and Cramer's model)."""
    rng = np.random.default_rng(seed)
    n = np.arange(3, N + 1, 2)
    logn = np.log(n)
    member = rng.random(len(n)) < np.minimum(2.0 / logn, 1.0)
    w = np.zeros(N + 1)
    w[n[member]] = logn[member]
    return w


def stage_aberration(res, Q: int):
    lam = von_mangoldt(X0)
    psi = float(lam.sum())
    n = np.arange(X0 + 1)
    sets = []
    # (name, description, weights or cached S tag, model spec)
    S_c = np.load(CACHE / f"S_cramer_N{X0}_M{M0}.npz")["S"]
    from sieve import cramer_weights
    sets.append(("Cramér", "random n ≥ 3 with probability 1/ln n, weight ln n: no congruence structure",
                 cramer_weights(X0), S_c, dict(D=1, allowed={0}, scale=1.0, base="integers")))
    del S_c
    wo = odd_cramer_weights(X0)
    sets.append(("odd Cramér", "random odd n with probability 2/ln n, weight ln n: parity only", wo, None,
                 dict(D=2, allowed={1}, scale=2.0, base="integers")))
    sets.append(("primes", "Λ(n): every congruence condition of the primes", lam, None, dict(D=1, allowed={0}, scale=1.0)))
    w4 = np.where(n % 4 == 1, 2.0 * lam, 0.0)
    sets.append(("primes ≡ 1 (mod 4), weight 2Λ", "one extra congruence condition", w4, None,
                 dict(D=4, allowed={1}, scale=2.0)))
    out = []
    for name, desc, w, S, spec in sets:
        if S is None:
            S = S_grid(w, M0)
        a2 = S.real**2 + S.imag**2
        del S
        L2 = float(np.sum(w**2))
        certify(f"Parseval, {name}", abs(a2.sum() / M0 - L2) <= 1e-9 * L2, "", res)
        psi_set = float(np.sum(w))
        model = cp.model_energies(range(1, Q + 1), X0, psi_set, **spec)
        h1 = float(np.max(np.abs(cp.class_model_heights(1, psi_set, **spec))))
        full = tuple(q for q in range(1, Q + 1) if float(np.max(np.abs(cp.class_model_heights(q, psi_set, **spec)), initial=0)) >= 0.999 * h1)
        b = budget(a2, X0, M0, model, Q, atom_qs=full)                  # atoms = exact copies of the main spike
        b_q12 = budget(a2, X0, M0, model, Q, atom_qs=(1, 2))           # atoms = q in {1, 2}, for comparison
        entry = {"name": name, "desc": desc, "psi": psi_set, "L2": L2, **b, "atom_qs": full,
                 "q12": {k: b_q12[k] for k in ("A", "F", "sea", "pred_A", "pred_F")}}
        if b["n_qF"] == 0:                                               # control: read the primes' forest windows
            pm = cp.model_energies(range(1, Q + 1), X0, psi_set, 1, {0})
            ctl = budget(a2, X0, M0, pm, Q, atom_qs=full)
            entry["F_control"] = ctl["F"]
        out.append(entry)
        del a2
        log(f"B3 {name}: atoms {b['A'] / L2:.4f}, forest {b['F'] / L2:.4f}, sea {b['sea'] / L2:.4f}")
    res["aberration"] = out
    res["aberration_Q"] = Q


def compute() -> dict:
    res: dict = {"certificates": []}
    mu, phi = mobius_phi(2 * 10**5)
    log("start")
    stage_development(res, mu, phi)
    Q0 = next(r["Q"] for r in res["development"] if r["X"] == X0)
    stage_aberration(res, Q0)
    res["runtime_s"] = time.time() - T0
    res["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    with open(CACHE / "compartments.pkl", "wb") as f:
        pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    log("saved cache/compartments.pkl")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="reuse cache/compartments.pkl")
    args = ap.parse_args()
    FIGS.mkdir(exist_ok=True)
    if args.report and (CACHE / "compartments.pkl").exists():
        res = pickle.load(open(CACHE / "compartments.pkl", "rb"))
    else:
        res = compute()
    import compartments_report
    compartments_report.figure(res, FIGS)
    compartments_report.write(res, HERE / "compartments.md")
    log("figure and compartments.md written")


if __name__ == "__main__":
    main()
