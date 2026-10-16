"""Escape routes from the binary L^2 wall: relaxations under which the circle method's
magnitude-only minor-arc inequality *does* close at X = 10^6, plus controls that do not.

For one even N the binary wall is the statement that P_m = int_m |S|^2 exceeds R2(N): any bound that
uses only |S| on the minor arcs is consistent with |I_m(N)| = P_m (phases aligned against us), so it
cannot show R2(N) > 0.  Beating that for a single N means proving the phase cancellation -- the open
problem; nothing here does that.  What we test instead, numerically and for these N only, is where the
same inequality closes once the question is relaxed in one of the ways the literature relaxes it:

  E0  the obstacle, measured: the minor-arc phasor walk, the bound ladder, and how much phase knowledge would suffice;
  E1  average over N -> failure count (the Bessel + Chebyshev step of van der Corput, Chudakov, Estermann 1937-38);
  E5  single-N transmission: the few targets for which the plain inequality already closes, and how that dies out;
  E2  add K powers of 2 -> Linnik's problem (Linnik 1951-53; Gallagher 1975; Heath-Brown-Puchta 2002; Pintz-Ruzsa 2020);
  E3  one prime -> almost-prime            (direction of Renyi 1948 and Chen 1973 -- which are sieve theorems);
  E4  smooth weights (control): cures the arc-edge sidelobes, not the wall.

"Closes" always means: (magnitude-only bound on the minor-arc part) < (major-arc part), computed on the
same grid with the same rectangle rule as the rest of the project.  It is a floating-point statement about
these N, not a theorem.
"""
from __future__ import annotations

import math

import numpy as np

from arcs import MajorArcs
from expsum import S_grid, chunked_sum, e_minus

# ---------------------------------------------------------------- helpers on the grid


def minor_mask(arcs: MajorArcs) -> np.ndarray:
    return ~arcs.mask()


def grid_mean(values: np.ndarray, mask: np.ndarray | None = None) -> float | complex:
    """(1/M) sum_k values[k] over the mask (rectangle rule)."""
    M = len(values)
    if mask is None:
        return values.sum() / M
    return values[mask].sum() / M


def omega_array(X: int) -> np.ndarray:
    """Omega(n) = number of prime factors with multiplicity, n = 0..X (Omega(0) = Omega(1) = 0)."""
    from sieve import primes_upto

    om = np.zeros(X + 1, dtype=np.int8)
    for p in primes_upto(X):
        pk = int(p)
        while pk <= X:
            om[pk::pk] += 1
            pk *= int(p)
    return om


def rough_indicator(X: int, z: int) -> np.ndarray:
    """1 for 2 <= m <= X with no prime factor < z, else 0 (float64)."""
    from sieve import primes_upto

    ind = np.ones(X + 1)
    ind[:2] = 0.0
    for p in primes_upto(max(z - 1, 1)):
        ind[p::p] = 0.0
    return ind


def smooth_bump(x: np.ndarray) -> np.ndarray:
    """eta(x) = exp(1 - 1/(1 - (2x-1)^2)) on (0, 1), 0 outside; C-infinity, eta(1/2) = 1."""
    x = np.asarray(x, dtype=float)
    u = 2.0 * x - 1.0
    out = np.zeros_like(x)
    inside = np.abs(u) < 1.0
    out[inside] = np.exp(1.0 - 1.0 / (1.0 - u[inside] ** 2))
    return out


def powers_of_two_grid(M: int, L: int) -> np.ndarray:
    """G(k/M) = sum_{j=1}^{L} e(2^j k / M), phases reduced exactly in integers."""
    G = np.zeros(M, dtype=complex)
    k = np.arange(M, dtype=np.int64)
    for j in range(1, L + 1):
        G += np.exp(2j * np.pi * (((k << j) % M) / M))
    return G


def kfold_power_counts(L: int, K: int, length: int, dtype=np.int64) -> np.ndarray:
    """c_K[s] = #{(k_1..k_K) in {1..L}^K : 2^{k_1} + ... + 2^{k_K} = s} (ordered), s < length.

    Use dtype=float for K >= 15 (19^15 > 2^63)."""
    c = np.zeros(length, dtype=dtype)
    c[0] = 1
    for _ in range(K):
        nxt = np.zeros(length, dtype=dtype)
        for j in range(1, L + 1):
            s = 1 << j
            if s < length:
                nxt[s:] += c[: length - s]
        c = nxt
    return c


# ================================================================ shared context

class Context:
    """Everything the experiments share at X = N = 10^6 (loaded once)."""

    def __init__(self, X: int, M: int, cache_dir, pairs):
        from expsum import load_or_compute_S
        from sieve import autoconvolve, mobius_phi, sum_lambda_sq, von_mangoldt

        self.X, self.M = X, M
        self.lam = von_mangoldt(X)
        self.L2 = sum_lambda_sq(self.lam)
        self.psi = float(self.lam.sum())
        self.S = load_or_compute_S(self.lam, M, "primes", cache_dir)
        self.r2 = autoconvolve(self.lam)                    # R2(N) for N <= 2X, exact
        self.R2 = float(self.r2[X])
        self.mu, self.phi = mobius_phi(10**5)
        self.pairs = list(pairs)
        self.arcs = {}

    def arc(self, Q: int, c: int) -> MajorArcs:
        """MajorArcs for (Q, c), cached in slim form (only idx and the membership test are kept)."""
        if (Q, c) not in self.arcs:
            A = MajorArcs(self.X, self.M, Q, c, self.mu, self.phi)
            for attr in ("q", "a", "beta", "mu_q", "phi_q"):
                setattr(A, attr, None)
            A.idx = A.idx.astype(np.int32)                  # M = 2^23 < 2^31; e_minus casts back to int64
            self.arcs[(Q, c)] = A
        return self.arcs[(Q, c)]


def chunked_demod_sum(F, N: int, M: int, chunk: int = 1 << 20) -> complex:
    """(1/M) sum_k F(k) e(-N k/M), with F a function of an index slice."""
    tot = 0.0j
    for s in range(0, M, chunk):
        k = np.arange(s, min(M, s + chunk), dtype=np.int64)
        tot += complex(np.sum(F(k) * e_minus(N, k, M)))
    return tot / M


def assert_real(z: np.ndarray, rel: float = 1e-9) -> np.ndarray:
    scale = max(float(np.max(np.abs(z.real))), 1e-300)
    if float(np.max(np.abs(z.imag))) > rel * scale:
        raise AssertionError(f"imaginary part too large: {np.max(np.abs(z.imag)) / scale:.2e}")
    return z.real.copy()


# ================================================================ E3: one prime -> almost-prime

def almost_prime_sets(X: int, seed: int = 7):
    """Yield (name, indicator array) for the sets used in E3 (generated lazily: 8 MB each)."""
    from sieve import primes_upto

    om = omega_array(X)
    n = np.arange(X + 1)
    base = n >= 2
    for r in range(1, 7):
        yield f"Omega<={r}", (base & (om <= r)).astype(float)
    for r in (2, 3):
        yield f"odd Omega<={r}", (base & (om <= r) & (n % 2 == 1)).astype(float)
        yield f"Omega<={r}, (m,10)=1", (base & (om <= r) & (n % 2 == 1) & (n % 5 != 0)).astype(float)
    yield "Chen-shaped: Omega<=2, (m,6)=1", (base & (om <= 2) & (n % 2 == 1) & (n % 3 != 0)).astype(float)
    for z in [int(p) for p in primes_upto(101) if p >= 3] + [317, 1000]:
        yield f"{z}-rough", rough_indicator(X, z)
    rng = np.random.default_rng(seed)
    for r in (2, 3):
        dens = float((base & (om <= r))[2:].mean())
        ind = np.zeros(X + 1)
        ind[2:] = (rng.random(X - 1) < dens).astype(float)
        yield f"Bernoulli({dens:.4f}) [control for Omega<={r}]", ind


def e3_almost_primes(ctx: Context, absS: np.ndarray, W: np.ndarray, sets, pairs_full, pairs_short) -> list:
    """For each set A: R^(A)(X) = sum_n Lambda(n) 1_A(X-n) = int S T_A e(-X a), major part, joint and CS bounds.

    W = S * e(-X k / M) on the grid (so that int_M S T e(-X a) = (1/M) sum_{idx} W T).
    """
    X, M, lam = ctx.X, ctx.M, ctx.lam
    tot_a2 = float(np.sum(absS ** 2))
    out = []
    for name, ind in sets:
        T = S_grid(ind, M)
        size = float(ind.sum())
        direct = math.fsum(lam[1:X] * ind[X - 1 : 0 : -1])
        grid = complex(np.sum(W * T) / M)
        absT = np.abs(T)
        tot_joint = float(np.dot(absS, absT))
        tot_T2 = float(np.dot(absT, absT))
        is_r = name.startswith("Omega<=") and "," not in name and int(name[-1]) <= 4
        rows = {}
        for (Q, c) in (pairs_full if is_r else pairs_short):
            A = ctx.arc(Q, c)
            idx = A.idx
            IM = float(np.real(np.sum(W[idx] * T[idx])) / M)
            joint = (tot_joint - float(np.dot(absS[idx], absT[idx]))) / M
            T2m = (tot_T2 - float(np.dot(absT[idx], absT[idx]))) / M
            Pm = (tot_a2 - float(np.dot(absS[idx], absS[idx]))) / M
            rows[(Q, c)] = {"I_M": IM, "joint": joint, "CS": math.sqrt(Pm * T2m), "f": T2m / size,
                            "joint_ratio": joint / IM, "CS_ratio": math.sqrt(Pm * T2m) / IM}
        out.append({"name": name, "size": size, "density": size / (X - 1), "direct": direct,
                    "grid": grid.real, "grid_imag": grid.imag, "rows": rows})
        del T, absT
    return out


# ================================================================ E0: the obstacle, measured

def components_and_positions(k: np.ndarray):
    """For sorted minor-arc indices k: component id of each (runs of consecutive indices) and position in it."""
    brk = np.r_[False, np.diff(k) > 1]
    comp = np.cumsum(brk)
    starts = np.flatnonzero(np.r_[True, brk[1:]])
    pos = np.arange(len(k)) - starts[comp]
    return comp, pos


def model_spikes(X: int, M: int, Qp: int, mu: np.ndarray, phi: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """S*(k/M) = (mu(q)/phi(q)) v(k/M - a/q) on the disjoint windows |alpha - a/q| <= 1/(2 q Qp), squarefree q <= Qp.

    Returned as a dense complex array (zero outside the windows); `keep` (bool) restricts to minor-arc points.
    """
    from expsum import dirichlet

    Sst = np.zeros(M, dtype=complex)
    for q in range(1, Qp + 1):
        if mu[q] == 0:
            continue
        a = np.arange(q, dtype=np.int64)
        a = a[np.gcd(a, q) == 1]
        hw = M / (2.0 * q * Qp)
        w = int(math.floor(hw))
        off = np.arange(-w - 1, w + 2, dtype=np.int64)
        kk = np.rint(a[:, None] * M / q).astype(np.int64) + off[None, :]          # candidate grid points
        num = kk * q - a[:, None] * M                                             # beta = num / (q M)
        inside = np.abs(num) * 2 * Qp <= M                                        # |beta| <= 1/(2 q Qp)
        kk, num = kk[inside], num[inside]
        kk %= M
        sel = keep[kk]
        kk, num = kk[sel], num[sel]
        Sst[kk] = (mu[q] / phi[q]) * dirichlet(num / (q * float(M)), X)
    return Sst


def e0_obstacle(ctx: Context, Q: int, c: int, IM_X: float, Bs=(1, 8, 64, 256, 1024, 4096, 16384, 65536, 262144),
                Qps=(150, 300, 500, 707, 1000, 1414, 2000, 4000)) -> dict:
    """Minor-arc phasor walk, bound ladder and partial-phase dial at N = X (needs S itself)."""
    X, M, S = ctx.X, ctx.M, ctx.S
    A = ctx.arc(Q, c)
    minor = ~A.mask()
    k = np.flatnonzero(minor)
    ph = e_minus(X, k, M)                                   # e(-N k / M) on the minor indices, reused below
    w = S[k]
    np.multiply(w, w, out=w)
    w *= ph
    w /= M
    F = np.cumsum(w)
    Pm = float(np.abs(w).sum())
    out = {"Q": Q, "c": c, "P_m": Pm, "I_m": complex(F[-1]), "max_abs_F": float(np.max(np.abs(F))), "I_M": IM_X,
           "n_minor": int(len(k))}
    stride = max(1, len(F) // 100000)
    out["walk"] = np.r_[F[::stride], F[-1]]
    out["walk_adversarial_length"] = Pm
    del F
    meas_m = len(k) / M
    # ---- partial-phase dial: |sum over windows of B consecutive minor points|, windows nested inside components
    comp, pos = components_and_positions(k)
    dial = []
    for B in Bs:
        inv = np.cumsum(pos % B == 0) - 1          # windows are contiguous runs: a new one starts where pos = 0 mod B
        re = np.bincount(inv, weights=w.real)
        im = np.bincount(inv, weights=w.imag)
        dial.append({"B": B, "N_Delta": B * X / M, "W": float(np.sum(np.hypot(re, im)) / IM_X)})
        if B == 8:   # random-phase rms with cells of width ~1/N
            out["random_phase_rms"] = float(math.sqrt(np.sum(re ** 2 + im ** 2)))
        del inv
    re = np.bincount(comp, weights=w.real)
    im = np.bincount(comp, weights=w.imag)
    out["dial"] = dial
    out["dial_floor"] = float(np.sum(np.hypot(re, im)) / IM_X)
    out["n_components"] = int(comp.max() + 1)
    out["rp_scale"] = Pm / math.sqrt(X * meas_m)
    del w, comp, pos
    # ---- bound ladder (units of R2 are applied by the caller)
    absS_m = np.abs(S[k])
    E4m = float(np.sum(absS_m ** 4) / M)
    out["holder"] = math.sqrt(meas_m * E4m)
    out["sup_L1"] = float(absS_m.max() * absS_m.sum() / M)
    out["E4m"] = E4m
    del absS_m
    model = []
    for Qp in Qps:
        Sst = model_spikes(X, M, Qp, ctx.mu, ctx.phi, minor)
        Sm = Sst[k]
        del Sst
        E = S[k]
        E -= Sm
        PE = float(np.sum(E.real ** 2 + E.imag ** 2) / M)
        del E
        Pst = float(np.sum(Sm.real ** 2 + Sm.imag ** 2) / M)
        np.multiply(Sm, Sm, out=Sm)
        Sm *= ph
        Jst = complex(np.sum(Sm) / M)
        model.append({"Qp": Qp, "P_E": PE, "P_star": Pst, "J_star": Jst,
                      "bound": PE + 2 * math.sqrt(Pst * PE) - Jst.real})
        del Sm
    out["model"] = model
    return out


# ================================================================ E1 (+ E5a): average over N

def even_targets(X: int) -> np.ndarray:
    lo = X // 2 + (X // 2) % 2
    return np.arange(lo, X + 1, 2, dtype=np.int64)


def major_part_all_N(S2: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """I_M(N) = (1/M) sum_{k in M} S_k^2 e(-N k / M) for every N mod M: one forward FFT (numpy sign)."""
    import scipy.fft

    M = len(S2)
    x = np.zeros(M, dtype=complex)
    x[idx] = S2[idx]
    FT = scipy.fft.fft(x, overwrite_x=True, workers=2)
    del x
    IM = assert_real(FT, rel=1e-9) / M
    del FT
    return IM


def e1_exceptional(ctx: Context, S2: np.ndarray, a2: np.ndarray, Q: int, c: int, SS2_even: np.ndarray) -> dict:
    """Bessel + Chebyshev over the even targets N in [X/2, X], and the single-N transmitted set."""
    X, M = ctx.X, ctx.M
    A = ctx.arc(Q, c)
    idx = A.idx
    IM = major_part_all_N(S2, idx)
    r2full = np.zeros(M)
    r2full[: len(ctx.r2)] = ctx.r2
    Im = r2full - IM
    del r2full
    tot4 = float(np.sum(a2 ** 2) / M)
    E4m = tot4 - float(np.sum(a2[idx] ** 2) / M)
    Pm = (float(np.sum(a2)) - float(np.sum(a2[idx]))) / M
    a2m = a2.copy()
    a2m[idx] = -1.0
    sup2 = float(a2m.max())
    del a2m
    bessel = float(np.sum(Im ** 2))
    Ns = even_targets(X)
    IMN, ImN = IM[Ns], Im[Ns]
    allN = np.arange(M)
    in_N = np.zeros(M, dtype=bool)
    in_N[Ns] = True
    le2X = allN <= 2 * X
    shares = {
        "targets": float(np.sum(Im[in_N] ** 2)) / bessel,
        "other even <= 2X": float(np.sum(Im[le2X & (allN % 2 == 0) & ~in_N] ** 2)) / bessel,
        "odd <= 2X": float(np.sum(Im[le2X & (allN % 2 == 1)] ** 2)) / bessel,
        "> 2X": float(np.sum(Im[~le2X] ** 2)) / bessel,
    }
    del Im, IM, in_N, le2X, allN
    minIM = float(IMN.min())
    nT = len(Ns)
    ratio = np.abs(ImN) / IMN
    fail = np.abs(ImN) >= IMN
    order = np.argsort(-ratio)[:10]
    meas = len(idx) / M
    pm_ratio = Pm / IMN
    trans = np.flatnonzero(pm_ratio < 1.0)
    return {
        "Q": Q, "c": c, "n_targets": nT, "P_m": Pm, "E4m": E4m, "full4": tot4, "sup2": sup2, "meas": meas,
        "kurtosis": E4m / (Pm ** 2 / (1.0 - meas)),
        "bessel_relerr": abs(bessel - E4m) / E4m,
        "minIM": minIM, "minIM_N": int(Ns[int(np.argmin(IMN))]),
        "bound_textbook": sup2 * ctx.L2 / minIM ** 2,
        "bound_natural_sup": sup2 * Pm / minIM ** 2,
        "bound_4th": E4m / minIM ** 2,
        "bound_oracle": float(np.sum(ImN ** 2)) / minIM ** 2,
        "F": int(fail.sum()), "F_harmful": int(np.sum(ImN <= -IMN)), "F_harmless": int(np.sum(ImN >= IMN)),
        "shares": shares,
        "quantiles": {qq: float(np.quantile(ratio, qq)) for qq in (0.5, 0.9, 0.99, 0.999, 1.0)},
        "worst": [(int(Ns[j]), float(ratio[j]), float(ImN[j] / IMN[j]), float(SS2_even[Ns[j]])) for j in order],
        "pm_ratio_min": float(pm_ratio.min()), "pm_ratio_median": float(np.median(pm_ratio)),
        "pm_ratio_max": float(pm_ratio.max()), "pm_ratio_min_N": int(Ns[int(np.argmin(pm_ratio))]),
        "transmitted": [(int(Ns[j]), float(pm_ratio[j]), float(SS2_even[Ns[j]])) for j in trans],
        "certificates_ok": bool(np.all(ctx.r2[Ns[trans]] >= IMN[trans] - Pm)) if len(trans) else True,
        "ratio_hist": np.histogram(np.log10(np.maximum(ratio, 1e-8)), bins=np.linspace(-6, 0.5, 131))[0],
        "SS2_targets": SS2_even[Ns], "pm_ratio_targets": pm_ratio.astype(np.float32),
    }


def e5_truncated(lam_full_X: int, N: int, c: int, mu, phi, cache_dir) -> dict:
    """Single-N certificate with S truncated at N itself (its own grid M >= 4N), arcs (Qmax(N, c), c) of half-width c/N."""
    from expsum import S_grid as _Sg
    from sieve import von_mangoldt

    lam = von_mangoldt(N)
    M = 1 << int(math.ceil(math.log2(4 * N)))
    S = _Sg(lam, M)
    Q = math.isqrt((N - 1) // (2 * c))
    A = MajorArcs(N, M, Q, c, mu, phi)
    a2 = S.real ** 2 + S.imag ** 2
    Pm = (float(a2.sum()) - float(a2[A.idx].sum())) / M
    IMN = float(np.real(np.sum(S[A.idx] ** 2 * e_minus(N, A.idx, M))) / M)
    return {"N": N, "M": M, "Q": Q, "c": c, "P_m": Pm, "I_M": IMN, "ratio": Pm / IMN}


def e5_scaling(X: int, S: np.ndarray, lam: np.ndarray, mu, phi, SS2_even: np.ndarray, c: int) -> dict:
    """Transmitted fraction at X for the arcs (Qmax(X, c), c), using a cached S_X (squared in place)."""
    from sieve import autoconvolve

    M = len(S)
    Q = math.isqrt((X - 1) // (2 * c))
    A = MajorArcs(X, M, Q, c, mu, phi)
    a2 = S.real ** 2 + S.imag ** 2
    Pm = (float(a2.sum()) - float(a2[A.idx].sum())) / M
    del a2
    np.multiply(S, S, out=S)
    IM = major_part_all_N(S, A.idx)
    Ns = even_targets(X)
    IMN = IM[Ns]
    del IM
    sel = IMN > Pm
    trans = Ns[sel]
    ratio = Pm / IMN
    near = np.flatnonzero(ratio < 1.01)
    r2 = autoconvolve(lam)
    return {"X": X, "Q": Q, "c": c, "P_m": Pm, "n_targets": len(Ns), "n_transmitted": int(len(trans)),
            "n_closes": int(np.sum(ratio < 0.99)), "near": [(int(Ns[j]), float(ratio[j])) for j in near],
            "meas": len(A.idx) / M,
            "need": Pm * SS2_even[X] / float(IMN[-1]), "max_SS2": float(SS2_even[Ns].max()),
            "transmitted": [int(t) for t in trans],
            "certificates_ok": bool(np.all(r2[trans] >= IMN[sel] - Pm - 1e-9 * Pm) and np.all(r2[trans] > 0))}


# ================================================================ E2: add K powers of 2 (Linnik)

def pow2_indicator(X: int) -> np.ndarray:
    ind = np.zeros(X + 1)
    j = 1
    while (1 << j) <= X:
        ind[1 << j] = 1.0
        j += 1
    return ind


def e2_linnik(ctx: Context, S2: np.ndarray, a2: np.ndarray, pairs, Kmax: int = 14, Kfull: int = 16,
              lambdas=None, Cs=(1, 2, 4), sparse_seeds: int = 8) -> dict:
    """R_K(X) = int S^2 G^K e(-X a), its major part, and magnitude-only bounds on the minor part."""
    X, M, L2 = ctx.X, ctx.M, ctx.L2
    L = int(math.floor(math.log2(X)))
    G = S_grid(pow2_indicator(X), M)
    aG = np.abs(G)
    lambdas = np.round(np.arange(0.30, 0.845, 0.01), 2) if lambdas is None else np.asarray(lambdas)
    # full circle: R_K (grid), moments T_K = int |S|^2 |G|^K, and int |G|^j
    acc = np.zeros(Kfull + 1, dtype=complex)
    TK = np.zeros(Kmax + 1)
    Gmom = np.zeros(2 * Kmax + 1)
    for s0 in range(0, M, 1 << 20):
        ch = np.arange(s0, min(M, s0 + (1 << 20)), dtype=np.int64)
        p = S2[ch] * e_minus(X, ch, M)
        g = G[ch]
        for K in range(Kfull + 1):
            acc[K] += p.sum()
            p *= g
        a = a2[ch].copy()
        ag = aG[ch]
        for K in range(Kmax + 1):
            TK[K] += a.sum()
            a *= ag
        gp = np.ones_like(ag)
        for j in range(2 * Kmax + 1):
            Gmom[j] += gp.sum()
            gp *= ag
    RK_grid = acc / M
    TK /= M
    Gmom /= M
    rev = ctx.r2[X::-1][: X + 1]                       # rev[s] = r2[X - s]
    RK_dir = np.array([math.fsum(kfold_power_counts(L, K, X + 1, dtype=float) * rev) for K in range(Kfull + 1)])
    OD = TK[2] - L * L2
    out = {"L": L, "RK_grid": RK_grid, "RK_direct": RK_dir, "TK": TK, "G_moments": Gmom, "OD": OD,
           "rho": RK_dir / (float(L) ** np.arange(Kfull + 1) * ctx.R2), "lambdas": lambdas, "rows": {}}
    tot4 = float(np.sum(a2 ** 2))
    for (Q, c) in pairs:
        A = ctx.arc(Q, c)
        idx = A.idx
        p = S2[idx] * e_minus(X, idx, M)
        g = G[idx]
        IMK = np.zeros(Kmax + 1)
        for K in range(Kmax + 1):
            IMK[K] = float(np.real(p.sum())) / M
            p *= g
        a = a2[idx].copy()
        ag = aG[idx]
        TKmaj = np.zeros(Kmax + 1)
        for K in range(Kmax + 1):
            TKmaj[K] = float(a.sum()) / M
            a *= ag
        gp = np.ones_like(ag)
        Gmaj = np.zeros(2 * Kmax + 1)
        for j in range(2 * Kmax + 1):
            Gmaj[j] = float(gp.sum()) / M
            gp *= ag
        BK = TK - TKmaj
        Gm = Gmom - Gmaj
        minor = ~A.mask()
        # share of B_K from minor points within R/N of alpha = 0 or 1/2 (the shoulders of those two arcs)
        shoulder = {}
        for R in (60, 200, 1000):
            r = int(R * M // X)
            near = np.r_[np.arange(0, r + 1), np.arange(M - r, M), np.arange(M // 2 - r, M // 2 + r + 1)]
            near = near[minor[near]]
            a, ag = a2[near].copy(), aG[near]
            sh = np.zeros(Kmax + 1)
            for K in range(Kmax + 1):
                sh[K] = float(a.sum()) / M
                a *= ag
            shoulder[R] = sh / BK

        gm, a2m = aG[minor], a2[minor]
        n_minor = len(gm)
        supG, sup2 = float(gm.max()), float(a2m.max())
        Pm, E4m = float(BK[0]), (tot4 - float(np.sum(a2[idx] ** 2))) / M
        order = np.argsort(gm)
        gs = gm[order]
        tail_a2 = np.r_[np.cumsum(a2m[order][::-1])[::-1], 0.0] / M
        del gm, a2m, order, minor
        Ks = np.arange(Kmax + 1)
        meanGK_m = Gm[: Kmax + 1] / (n_minor / M)
        thm = {}
        for C in Cs:
            best_split, best_lit, best_exact = [], [], []
            for K in Ks:
                if K < 2:
                    best_split.append((np.inf, None)); best_lit.append((np.inf, None)); best_exact.append((np.inf, None))
                    continue
                vals_s, vals_l, vals_e = [], [], []
                for lam_ in lambdas:
                    j = int(np.searchsorted(gs, lam_ * L))
                    measE = (len(gs) - j) / M
                    base = (lam_ * L) ** (K - 2)
                    vals_s.append(base * (L * L2 + C * OD) + float(L) ** K * sup2 * measE)
                    vals_l.append(base * C * TK[2] + float(L) ** K * sup2 * measE)
                    vals_e.append(base * (L * L2 + C * OD) + float(L) ** K * tail_a2[j])
                for vals, best in ((vals_s, best_split), (vals_l, best_lit), (vals_e, best_exact)):
                    jj = int(np.argmin(vals))
                    best.append((vals[jj] / IMK[K], float(lambdas[jj])))
            thm[C] = {"split": best_split, "literature": best_lit, "split_exact": best_exact,
                      "no_split": [np.inf if K < 2 else supG ** (K - 2) * (L * L2 + C * OD) / IMK[K] for K in Ks]}
        rayleigh = [Pm * math.gamma(1 + K / 2) * L ** (K / 2) / IMK[K] for K in Ks]   # B_K ~ P_m E|G|^K, over I_M,K
        out["rows"][(Q, c)] = {
            "I_MK": IMK, "B_K": BK, "joint_ratio": BK / IMK, "sup_G": supG,
            "sup_ratio": supG ** Ks * Pm / IMK,
            "holder_ratio": np.sqrt(E4m) * np.sqrt(Gm[2 * Ks]) / IMK,
            "corr_ratio": BK / (Pm * meanGK_m), "rayleigh_ratio": np.array(rayleigh),
            "theorem": thm,
            "measE": {float(l): (len(gs) - int(np.searchsorted(gs, l * L))) / M for l in (0.5, 0.6, 0.7, 0.8)},
            "shoulder_share": shoulder,
        }
    # sparse-target control: 19 random even targets of matching scale instead of the powers of 2
    sparse = {}
    for (Q, c) in [p for p in pairs if p in ((150, 20), (20, 10))]:
        A = ctx.arc(Q, c)
        idx = A.idx
        ratios = []
        for sd in range(sparse_seeds):
            rng = np.random.default_rng([11, sd])
            chosen = set()
            for k in range(1, L + 1):
                lo, hi = max(2, (1 << (k - 1)) + 1), min(1 << L, 1 << (k + 1))
                cand = [t for t in range(lo + (lo % 2), hi + 1, 2) if t not in chosen]
                chosen.add(int(rng.choice(cand)))
            ind = np.zeros(X + 1)
            ind[sorted(chosen)] = 1.0
            H = S_grid(ind, M)
            aH = np.abs(H)
            B1 = (float(np.dot(a2, aH)) - float(np.dot(a2[idx], aH[idx]))) / M
            I1 = float(np.real(np.sum(S2[idx] * H[idx] * e_minus(X, idx, M)))) / M
            ratios.append(B1 / I1)
            del H, aH
        sparse[(Q, c)] = ratios
    out["sparse_control"] = sparse
    return out


# ================================================================ E4: smooth weight (control)

def eta_constants(h: float = 2e-4) -> dict:
    """int eta, int eta^2 = (eta*eta)(1), (eta*eta*eta)(1) and (eta*eta*eta)(1.5) by sampled convolution."""
    x = np.arange(0.0, 1.0 + h / 2, h)
    e = smooth_bump(x)
    ee = np.convolve(e, e) * h                 # (eta*eta)(t) at t = j h
    eee = np.convolve(ee, e) * h
    at = lambda arr, t: float(arr[int(round(t / h))])  # noqa: E731
    return {"int_eta": float(np.trapezoid(e, x)), "int_eta2": float(np.trapezoid(e * e, x)),
            "eta3_at_1": at(eee, 1.0), "eta3_at_1.5": at(eee, 1.5), "eta2_at_1": at(ee, 1.0)}


def eta_hat_table(xi_max: float = 300.0, step: float = 0.005, n_x: int = 4001) -> tuple[np.ndarray, np.ndarray]:
    """|eta^(xi)| / eta^(0) for xi in [0, xi_max], eta^(xi) = int_0^1 eta(x) e(-xi x) dx (trapezoid)."""
    x = np.linspace(0.0, 1.0, n_x)
    e = smooth_bump(x)
    wts = np.full(n_x, x[1] - x[0])
    wts[[0, -1]] *= 0.5
    xi = np.arange(0.0, xi_max + step / 2, step)
    vals = np.empty(len(xi))
    for s0 in range(0, len(xi), 500):
        ph = np.exp(-2j * np.pi * np.outer(xi[s0 : s0 + 500], x))
        vals[s0 : s0 + 500] = np.abs(ph @ (e * wts))
    return xi, vals / vals[0]


def sharp_volume(t: float) -> float:
    """Area of {t1 + t2 + t3 = t, 0 <= t_i <= 1} (the singular integral of the sharp ternary count)."""
    return t * t / 2 if t <= 1 else (-2 * t * t + 6 * t - 3) / 2


def e4_smooth(ctx: Context, pairs, sharp_sup: dict, Nters) -> dict:
    """Binary wall and ternary sup bound with w(n) = Lambda(n) eta(n/X), against the sharp cutoff."""
    from arcs import explain_alpha, minor_sup
    from expsum import nonzero_terms
    from singular_series import SS2_product, SS3_product
    from sieve import autoconvolve

    X, M = ctx.X, ctx.M
    w = ctx.lam * smooth_bump(np.arange(X + 1) / X)
    Se = S_grid(w, M)
    r2e = autoconvolve(w)
    consts = eta_constants()
    Rbin = float(r2e[X])
    Rbin_grid = chunked_demod_sum(lambda k: Se[k] ** 2, X, M).real
    sum_w2 = float(np.sum(w ** 2))
    psi_e = float(w.sum())
    ter = {}
    for Np in Nters:
        m = np.arange(max(0, Np - X), min(Np, len(r2e) - 1) + 1)
        R3 = math.fsum(r2e[m] * w[Np - m])
        R3_grid = chunked_demod_sum(lambda k: Se[k] ** 3, Np, M).real
        t = Np / X
        e3 = consts["eta3_at_1"] if abs(t - 1) < 1e-3 else consts["eta3_at_1.5"]
        ter[Np] = {"R3": R3, "R3_grid": R3_grid, "main_eta": SS3_product(Np) * X * X * e3,
                   "main_sharp": SS3_product(Np) * X * X * sharp_volume(t),
                   "geom_factor": consts["int_eta"] * consts["int_eta2"] / (e3 / sharp_volume(t))}
    nz, wn = nonzero_terms(w)
    a2e = Se.real ** 2 + Se.imag ** 2
    rows = {}
    for (Q, c) in pairs:
        A = ctx.arc(Q, c)
        idx = A.idx
        Pm = (float(a2e.sum()) - float(a2e[idx].sum())) / M
        IM = float(np.real(np.sum(Se[idx] ** 2 * e_minus(X, idx, M)))) / M
        sup = minor_sup(Se, A, nz, wn)
        ex = explain_alpha(sup["alpha"], X, psi_e, ctx.mu, ctx.phi)
        row = {"P_m_over_R": Pm / Rbin, "I_m_over_R": (Rbin - IM) / Rbin, "sup": sup["sup"],
               "sup_at": (ex["a"], ex["q"], ex["N_beta"])}
        for Np in Nters:
            row[f"B3_eta_{Np}"] = sup["sup"] * sum_w2 / ter[Np]["main_eta"]
            if (Q, c) in sharp_sup:
                row[f"B3_sharp_{Np}"] = sharp_sup[(Q, c)] * ctx.L2 / ter[Np]["main_sharp"]
        rows[(Q, c)] = row
    xi, eh = eta_hat_table()
    side = {}
    for c in (5, 10, 20):
        tail = float(np.sum(eh[xi > c] ** 2) / np.sum(eh ** 2))
        side[c] = {"max_sidelobe": float(eh[xi > c].max()), "kernel_tail": tail,
                   "sharp_sidelobe": 1 / (math.pi * c), "sharp_tail": 1 / (math.pi ** 2 * c)}
    # zoom profile near alpha = 0 for the figure: |S| and |S_eta| against N beta in [0, 30]
    from expsum import S_direct
    nS, wS = nonzero_terms(ctx.lam)
    betas = np.linspace(0, 30 / X, 1201)
    prof = {"Nbeta": betas * X, "S": np.abs(S_direct(betas, nS, wS)) / ctx.psi,
            "S_eta": np.abs(S_direct(betas, nz, wn)) / psi_e}
    return {"R_eta": Rbin, "R_eta_grid": Rbin_grid, "main_bin": SS2_product(X) * X * consts["int_eta2"],
            "sum_w2": sum_w2, "psi_eta": psi_e, "constants": consts, "ternary": ter, "rows": rows,
            "sidelobes": side, "profile": prof}
