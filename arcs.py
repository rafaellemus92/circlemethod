"""Major and minor arcs on the uniform grid alpha_k = k/M, and integrals over them.

Major arcs M(Q, c) = union over q <= Q, gcd(a, q) = 1 of {alpha : ||alpha - a/q|| <= c/N}.
They are pairwise disjoint if 2 c Q^2 < N (distinct a/q, a'/q' with q, q' <= Q are at
least 1/(q q') >= 1/Q^2 apart), which we assert as the brief asks.  (The sharp condition is
2 c Q (Q-1) < N, since the closest pair is 1/Q, 1/(Q-1).)  Minor arcs m = the complement.

One quadrature rule is used everywhere: the rectangle rule on the grid, with a grid
point belonging to an arc iff its exact distance to a/q is <= c/N (tested in integer
arithmetic).  Hence I_M + I_m = (full grid integral) = R2(N) identically.
"""
from __future__ import annotations

import math

import numpy as np

from expsum import S_direct, dirichlet, e_minus, kernel_capture


def admissible(N: int, Q: int, c: float) -> bool:
    return 2 * c * Q * Q < N


class MajorArcs:
    def __init__(self, N: int, M: int, Q: int, c: int, mu: np.ndarray, phi: np.ndarray):
        if not admissible(N, Q, c):
            raise AssertionError(f"(Q, c) not admissible: 2cQ^2 = {2 * c * Q * Q} >= N = {N}")
        if int(c) != c:
            raise ValueError("c must be an integer (exact integer membership test)")
        c = int(c)
        self.N, self.M, self.Q, self.c = N, M, Q, c
        idx, qq, aa, bnum = [], [], [], []
        for q in range(1, Q + 1):
            a = np.array([a for a in range(q) if math.gcd(a, q) == 1], dtype=np.int64)
            # |k/M - a/q| <= c/N  <=>  |k q N - a M N| <= c M q   (all integers)
            lo = -((-(a * M * N - c * M * q)) // (q * N))   # ceil
            hi = (a * M * N + c * M * q) // (q * N)         # floor
            for ai, l, h in zip(a, lo, hi):
                k = np.arange(l, h + 1, dtype=np.int64)
                idx.append(k % M)
                qq.append(np.full(len(k), q, dtype=np.int64))
                aa.append(np.full(len(k), ai, dtype=np.int64))
                bnum.append(k * q - ai * M)                 # beta = bnum / (q M), exact
        self.idx = np.concatenate(idx)
        self.q = np.concatenate(qq)
        self.a = np.concatenate(aa)
        self.beta = np.concatenate(bnum) / (self.q.astype(float) * M)
        self.n_arcs = int(sum(phi[1 : Q + 1]))
        # disjointness on the grid: no grid point claimed twice
        if len(np.unique(self.idx)) != len(self.idx):
            raise AssertionError("arcs overlap on the grid")
        self.mu_q = mu[self.q]
        self.phi_q = phi[self.q]

    @property
    def measure(self) -> float:
        return len(self.idx) / self.M

    def mask(self) -> np.ndarray:
        m = np.zeros(self.M, dtype=bool)
        m[self.idx] = True
        return m

    def model(self) -> np.ndarray:
        """Major-arc model (mu(q)/phi(q)) v(beta) at every major-arc grid point."""
        return (self.mu_q / self.phi_q) * dirichlet(self.beta, self.N)

    def contains(self, alpha) -> np.ndarray:
        """Exact-in-floats membership test for arbitrary alpha (vectorised)."""
        alpha = np.atleast_1d(np.asarray(alpha, dtype=float)) % 1.0
        inside = np.zeros(len(alpha), dtype=bool)
        for q in range(1, self.Q + 1):
            x = alpha * q
            inside |= np.abs(x - np.rint(x)) <= self.c * q / self.N
        return inside


def arc_integrals(S: np.ndarray, arcs: MajorArcs, N: int, Nter=()) -> dict:
    """Integrals over the major arcs (grid rule) of |S|^2, S^2 e(-N a), model^2 e(-N a), S^3 e(-N' a)."""
    k = arcs.idx
    Sk = S[k]
    M = len(S)
    eN = e_minus(N, k, M)
    model = arcs.model()
    out = {
        "P_M": float(np.sum(Sk.real ** 2 + Sk.imag ** 2) / M),
        "I_M": complex(np.sum(Sk * Sk * eN) / M),
        "I_M_model": complex(np.sum(model * model * eN) / M),
    }
    for Np in Nter:
        out[f"I3_M_{Np}"] = complex(np.sum(Sk ** 3 * e_minus(Np, k, M)) / M)
    return out


def minor_sup(S: np.ndarray, arcs: MajorArcs, n: np.ndarray, w: np.ndarray, n_cand: int = 12) -> dict:
    """sup over the minor arcs of |S|, from the grid then refined off-grid.

    The grid can miss a crest by up to 1/(2M) in alpha (a 0.6% dip on a main lobe, up to
    ~1.8% on a Dirichlet sidelobe, for M/N = 8.4), so the top grid candidates are re-evaluated by direct summation on a
    fine local mesh, keeping only points that are themselves in the minor arcs.
    """
    M = len(S)
    absS = np.abs(S)
    absS[arcs.idx] = -1.0
    grid_sup = float(absS.max())
    # top candidates that are distinct local peaks
    cand = np.argpartition(absS, -40 * n_cand)[-40 * n_cand :]
    cand = cand[np.argsort(-absS[cand])]
    picked = []
    for k in cand:
        if all(min(abs(int(k) - p), M - abs(int(k) - p)) > 4 for p in picked):
            picked.append(int(k))
        if len(picked) >= n_cand:
            break
    best, best_alpha = grid_sup, int(np.argmax(absS)) / M
    for k in picked:
        loc = (k + np.linspace(-1.0, 1.0, 41)) / M
        loc = loc[~arcs.contains(loc)]
        if len(loc) == 0:
            continue
        v = np.abs(S_direct(loc, n, w))
        j = int(np.argmax(v))
        if v[j] > best:
            best, best_alpha = float(v[j]), float(loc[j] % 1.0)
    return {"grid_sup": grid_sup, "sup": best, "alpha": best_alpha}


def explain_alpha(alpha: float, N: int, psi: float, mu: np.ndarray, phi: np.ndarray, qmax: int = 1000) -> dict:
    """Which rational's spike is responsible for |S(alpha)|?

    For every q <= qmax take the nearest a/q and the model height
    psi |mu(q)| / phi(q) * |v(beta)| / N; return the fraction with the largest one.
    """
    q = np.arange(1, qmax + 1)
    a = np.rint(alpha * q).astype(np.int64)
    beta = alpha - a / q
    ok = np.gcd(a, q) == 1
    h = np.where(ok, psi * np.abs(mu[q]) / phi[q] * np.abs(dirichlet(beta, N)) / N, 0.0)
    j = int(np.argmax(h))
    return {"a": int(a[j] % q[j]), "q": int(q[j]), "N_beta": float(N * beta[j]), "model_height": float(h[j])}


# ---------------------------------------------------------------- per-(N, Q, c) table rows
def minor_integrals(S: np.ndarray, mask: np.ndarray, N: int, chunk: int = 1 << 20) -> tuple[float, complex]:
    """(P_m, I_m) summed directly over the minor-arc grid points (not by subtraction)."""
    M = len(S)
    P, I = 0.0, 0.0j
    for s in range(0, M, chunk):
        sl = slice(s, min(M, s + chunk))
        keep = ~mask[sl]
        k = np.arange(sl.start, sl.stop, dtype=np.int64)[keep]
        Sk = S[sl][keep]
        P += float(np.sum(Sk.real ** 2 + Sk.imag ** 2))
        I += complex(np.sum(Sk * Sk * e_minus(N, k, M)))
    return P / M, I / M


def arc_rows(*, S, N, pairs, mu, phi, SS2_part, SS2_full, R2, L2, psi, n, w,
             Nter=(), R3=None, SS3=None, J_cache=None, label="primes"):
    """One dict per admissible (Q, c).  R3 / SS3: dicts keyed by N' (ternary columns)."""
    rows = []
    J_cache = {} if J_cache is None else J_cache
    for Q, c in pairs:
        if not admissible(N, Q, c):
            continue
        A = MajorArcs(N, len(S), Q, c, mu, phi)
        d = arc_integrals(S, A, N, Nter)
        P_m, I_m = minor_integrals(S, A.mask(), N)
        sup = minor_sup(S, A, n, w)
        ex = explain_alpha(sup["alpha"], N, psi, mu, phi)
        key = (N, c)
        if key not in J_cache:
            J_cache[key] = kernel_capture(N, c).real
        J = J_cache[key]
        SS2Q = float(SS2_part[Q])
        row = {
            "set": label, "N": N, "Q": Q, "c": c,
            "measure": A.measure, "n_arcs": A.n_arcs,
            "R2": R2, "SS2": SS2_full, "SS2N": SS2_full * N, "SS2Q": SS2Q, "SS2QN": SS2Q * N,
            "I_M": d["I_M"].real, "I_M_imag": d["I_M"].imag, "I_M_model": d["I_M_model"].real,
            "I_m": I_m.real, "I_m_imag": I_m.imag,
            "closure": (d["I_M"].real + I_m.real - R2) / R2,
            "P_M": d["P_M"], "P_m": P_m, "P_m_sub": L2 - d["P_M"],
            "P_m_over_R2": P_m / R2, "I_m_over_R2": I_m.real / R2,
            "sqrtP_m": math.sqrt(P_m), "P_m_over_sqrtN": P_m / math.sqrt(N),
            "J_over_N": J / N,
            # deterministic leakage into the minor arcs, normalised by the full-circle kernel
            # integral int v^2 e(-N b) = N - 1: unresolved spikes q > Q, and kernel tails past c/N
            "leak_tail_q": (SS2_full - SS2Q) * (N - 1),
            "leak_kernel": SS2Q * ((N - 1) - J),
            # error of the major-arc model on the major arcs themselves
            "model_err": d["I_M"].real - d["I_M_model"].real,
            # random-phase yardstick: N|m| independent coherence cells of width 1/N
            "rp_scale": P_m / math.sqrt(N * (1.0 - A.measure)),
            "sup": sup["sup"], "sup_grid": sup["grid_sup"], "sup_alpha": sup["alpha"],
            "sup_a": ex["a"], "sup_q": ex["q"], "sup_Nbeta": ex["N_beta"], "sup_model": ex["model_height"],
            "pred_Pm_over_R2": (math.log(N) - math.log(Q) - 2.3326) / SS2_full if SS2_full else float("nan"),
            # same heuristic with the exact sum_{q<=Q} mu^2/phi and the kernel capture J(c)
            "pred_exact_Pm_over_R2": (L2 - J * float(np.sum((mu[1:Q + 1] ** 2) / phi[1:Q + 1]))) / R2,
            "B3": sup["sup"] * L2,
        }
        row["leak"] = row["leak_tail_q"] + row["leak_kernel"]
        row["I_m_resid"] = row["I_m"] - row["leak"]
        for Np in Nter:
            main = SS3[Np] * Np ** 2 / 2
            I3M = d[f"I3_M_{Np}"].real
            row[f"main3_{Np}"] = main
            row[f"B3_main_{Np}"] = row["B3"] / main
            row[f"I3_M_{Np}"] = I3M
            row[f"I3_m_{Np}"] = R3[Np] - I3M
        rows.append(row)
    return rows
