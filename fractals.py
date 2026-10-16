"""The fractal lens: Mandelbrot's measurement programme applied to the circle-method objects at X = 10^6.

Three of Mandelbrot's principles, used as measurements rather than metaphors:

  scaling only between cutoffs   -- data collapse of the spikes, window maxima, and where each scaling window ends;
  roughness vs intermittency     -- Hurst exponents of the minor-arc walk; D_s and the coarse-grained f(alpha) of |S|^2 dalpha;
  mild vs wild randomness        -- the tail of |S|^2 on the minor arcs, by (Q, c).

What the measurements find is one fact read several ways: the Farey count #{a/q : mu(q) != 0, phi(q) <= y} ~ y^2/2 applied to
spikes of height psi/phi(q) and width 1/N, on top of a structureless sea.  A surrogate made of the singular-series spikes
and a Gaussian sea of matched power reproduces it.  The Mandelbrot set's p/q bulbs share the index set, the circle order and a
q^-2 size law with those spikes (through Diophantine approximation, not dynamics).

Every exponent is a statistic at one finite N over a finite window; nothing here bears on a proof.
"""
from __future__ import annotations

import math

import numpy as np

EULER_GAMMA = 0.5772156649015329
SINC2_HALF = 0.7737                       # int_{|x|<=1/2} sinc^2, for reference (computed exactly by box_capture(1))

# ================================================================ measures on the grid


def remove_atoms(a2: np.ndarray, X: int, radius: float = 20.0) -> np.ndarray:
    """Copy of |S|^2 on the grid with the trivial atoms at alpha = 0 and 1/2 (|N alpha| <= radius) set to zero."""
    M = len(a2)
    out = a2.copy()
    r = int(math.floor(radius * M / X))
    out[: r + 1] = 0.0
    out[M - r :] = 0.0
    out[M // 2 - r : M // 2 + r + 1] = 0.0
    return out


def kept_L(a2_kept: np.ndarray, X: int, psi: float) -> float:
    """Spike-mass normaliser L = X P_kept / psi^2 with P_kept = (1/M) sum of the kept |S|^2 (spike q has mass 1/(L phi^2))."""
    return float(X * (a2_kept.sum() / len(a2_kept)) / psi**2)


def box_masses(mu: np.ndarray, j: int) -> np.ndarray:
    """Masses of the boxes of 2^j consecutive grid points (len(mu) must be a power of two >= 2^j)."""
    return mu.reshape(-1, 1 << j).sum(axis=1)


def falpha_histogram(boxes: np.ndarray, eps: float, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Histogram estimate of the singularity spectrum at one box size eps (kept for tests; not used in the figures).

    alpha_box = log(mass) / log(eps);  N(alpha) = #boxes per bin;  f(alpha) = log N(alpha) / log(1/eps).
    """
    b = boxes[boxes > 0]
    al = np.log(b) / math.log(eps)
    h, _ = np.histogram(al, bins)
    f = np.full(len(h), np.nan)
    f[h > 0] = np.log(h[h > 0]) / math.log(1.0 / eps)
    return 0.5 * (bins[1:] + bins[:-1]), f, h


def rank_mass(mu: np.ndarray, j: int, imax: int) -> np.ndarray:
    """The imax largest box masses at box size 2^j, in decreasing order (m_(1) >= m_(2) >= ...)."""
    b = box_masses(mu, j)
    k = min(imax, len(b))
    top = np.partition(b, len(b) - k)[len(b) - k :]
    return np.sort(top)[::-1]


def cumulative_spectrum(boxes: np.ndarray, eps: float, alphas) -> tuple[np.ndarray, np.ndarray]:
    """N_<(alpha) = #{boxes with mass >= eps^alpha} and f_<(alpha) = ln N_< / ln(1/eps)."""
    b = np.sort(boxes[boxes > 0])
    thr = eps ** np.asarray(alphas, dtype=float)
    n = len(b) - np.searchsorted(b, thr, side="left")
    with np.errstate(divide="ignore"):
        f = np.where(n > 0, np.log(np.maximum(n, 1)) / math.log(1.0 / eps), np.nan)
    return n, f


def chhabra_jensen(mu: np.ndarray, js, s_list) -> dict:
    """Direct (Chhabra-Jensen 1989) alpha(s), f(s): local slopes, between consecutive box sizes, of the s-weighted averages
    A_j(s) = sum m_i log mu_i and F_j(s) = sum m_i log m_i (m_i = mu_i^s / sum mu^s) against log eps."""
    M = len(mu)
    out = {s: {"A": [], "F": []} for s in s_list}
    for j in js:
        b = box_masses(mu, j)
        b = b[b > 0]
        lb = np.log(b)
        for s in s_list:
            w = np.exp(s * lb - np.max(s * lb))
            m = w / w.sum()
            out[s]["A"].append(float(np.sum(m * lb)))
            out[s]["F"].append(float(np.sum(m[m > 0] * np.log(m[m > 0]))))
    le = np.log(np.array([2.0**j / M for j in js]))
    res = {}
    for s in s_list:
        A, F = np.array(out[s]["A"]), np.array(out[s]["F"])
        res[s] = {"alpha": np.diff(A) / np.diff(le), "f": np.diff(F) / np.diff(le), "log_eps_mid": 0.5 * (le[1:] + le[:-1])}
    return res


def partition_dims(mu: np.ndarray, js, s_list) -> dict:
    """Local generalized dimensions between consecutive box sizes: D_s = dlog Z_s / ((s - 1) dlog eps) with
    Z_s = sum mu_box^s; for s = 1 the information dimension D_1 = -dH/dlog eps with H = -sum mu log mu (Shannon)."""
    M = len(mu)
    Z = {s: [] for s in s_list}
    for j in js:
        b = box_masses(mu, j)
        b = b[b > 0]
        for s in s_list:
            Z[s].append(float(-np.sum(b * np.log(b))) if s == 1 else float(np.sum(b**s)))
    le = np.log(np.array([2.0**j / M for j in js]))
    out = {}
    for s in s_list:
        z = np.array(Z[s])
        out[s] = -np.diff(z) / np.diff(le) if s == 1 else np.diff(np.log(z)) / np.diff(le) / (s - 1.0)
    return out


# ================================================================ the Farey count and the sinc kernel


def sinc_F(t):
    """int_0^t sinc^2(x) dx = Si(2 pi t)/pi - sin^2(pi t)/(pi^2 t), sinc x = sin(pi x)/(pi x); odd in t, F(inf) = 1/2."""
    from scipy.special import sici

    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    nz = t != 0
    tn = t[nz]
    out[nz] = sici(2 * np.pi * tn)[0] / np.pi - np.sin(np.pi * tn) ** 2 / (np.pi**2 * tn)
    return out


def box_capture(w: float, n: int = 20001) -> float:
    """kappa(w) = (1/w) int_0^w du int_{-u}^{w-u} sinc^2 = (2/w) int_0^w F(u) du: the expected share of one spike's mass
    caught by a box of width w (in units of 1/N) whose position relative to the spike is uniform."""
    u = np.linspace(0.0, w, n)
    return float(2.0 / w * np.trapezoid(sinc_F(u), u))


def squarefree_by_phi(mu: np.ndarray, phi: np.ndarray, qmin: int = 3, qmax: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Squarefree q in [qmin, qmax] sorted by phi(q) (ties by q): arrays (q, phi(q))."""
    qmax = len(mu) - 1 if qmax is None else qmax
    q = np.arange(qmin, qmax + 1)
    q = q[mu[q] != 0]
    o = np.lexsort((q, phi[q]))
    return q[o], phi[q][o]


def farey_staircase(L: float, mu: np.ndarray, phi: np.ndarray, kappa: float = 1.0, sea: float = 0.0,
                    qmax: int = 20000) -> dict:
    """Exact Farey staircase for the rank-mass law: each squarefree q >= 3 contributes phi(q) boxes of mass
    kappa/(L phi^2) + sea, ordered by phi.  Returns per-q arrays: q, phi, rank_end (cumulative sum of phi), mass."""
    q, ph = squarefree_by_phi(mu, phi, 3, qmax)
    return {"q": q, "phi": ph, "rank_end": np.cumsum(ph), "mass": kappa / (L * ph.astype(float) ** 2) + sea}


def spike_count(mu: np.ndarray, phi: np.ndarray, y: float) -> int:
    """#{a/q : mu(q) != 0, phi(q) <= y} = sum over squarefree q with phi(q) <= y of phi(q) (q runs over the whole table)."""
    q = np.arange(1, len(mu))
    sel = (mu[q] != 0) & (phi[q] <= y)
    return int(phi[q][sel].sum())


# ================================================================ M1: spike data collapse


def fg_factor(q: int, primes_of_q) -> float:
    """FG(q) = ln q - gamma - ln 2 pi - sum_{p | q} ln p / (p - 1)  (Friedlander-Goldston form of V(q)/X)."""
    return math.log(q) - EULER_GAMMA - math.log(2 * math.pi) - sum(math.log(p) / (p - 1) for p in primes_of_q)


def spike_collapse(n: np.ndarray, w: np.ndarray, psi: float, qs, mu: np.ndarray, phi: np.ndarray) -> list[dict]:
    """For each squarefree q: S(a/q) for every reduced a by one inverse FFT of the residue sums psi(X; q, r).

    S(a/q) = sum_r psi(X; q, r) e(a r / q) = q ifft(psi_res)[a].  With Z_a = mu(q) S(a/q) phi(q)/psi - 1 the model predicts
    Z = 0; the sea makes Z a complex random variable whose exact second moments E|Z|^2, E[Z^2] (real, because a -> -a
    conjugates) feed the Beckmann model.  V(q) is the residue-class variance over reduced classes.
    """
    out = []
    for q in qs:
        q = int(q)
        res = np.bincount(n % q, weights=w, minlength=q)
        Sa = q * np.fft.ifft(res)
        a = np.arange(q)
        red = np.gcd(a, q) == 1 if q > 1 else np.array([True])
        Sr = Sa[red]
        y = np.abs(Sr) * phi[q] / psi
        Z = mu[q] * Sr * phi[q] / psi - 1.0
        E2 = float(np.mean(np.abs(Z) ** 2))
        EZ2 = float(np.mean(Z**2).real)
        pc = res[red].sum()
        V = float(np.sum((res[red] - pc / phi[q]) ** 2))
        out.append({"q": q, "mu": int(mu[q]), "phi": int(phi[q]), "n_a": int(red.sum()), "mean": float(y.mean()),
                    "sd": float(y.std()), "cos": float(np.mean((mu[q] * Sr).real / np.abs(Sr))),
                    "E2": E2, "EZ2": EZ2, "V": V, "psi_c": float(pc)})
    return out


def beckmann_moments(sx2: float, sy2: float, n: int = 80) -> dict:
    """y = |1 + X + iY|, X ~ N(0, sx2), Y ~ N(0, sy2) independent: mean, sd and <cos theta> = E[(1+X)/y] by an n x n
    probabilists' Gauss-Hermite rule; sd^2 = 1 + sx2 + sy2 - mean^2 uses only the exact second moment."""
    from numpy.polynomial.hermite_e import hermegauss

    x, wt = hermegauss(n)
    wt = wt / math.sqrt(2 * math.pi)
    XX = math.sqrt(max(sx2, 0.0)) * x[:, None]
    YY = math.sqrt(max(sy2, 0.0)) * x[None, :]
    W = wt[:, None] * wt[None, :]
    y = np.hypot(1 + XX, YY)
    mean = float(np.sum(W * y))
    return {"mean": mean, "sd": math.sqrt(max(1 + sx2 + sy2 - mean**2, 0.0)), "cos": float(np.sum(W * (1 + XX) / y))}


def rice_moments(sigma2: float) -> dict:
    """Isotropic Rice model, total variance sigma2 = E|Z|^2: y = |1 + Z|; mean, sd and <cos theta> in closed form, using
    exponentially scaled Bessel functions (stable for tiny sigma2, where scipy's rice moments overflow):
        mean = sc sqrt(pi/2) [(1 + x) I0e(x/2) + x I1e(x/2)],  x = 1/(2 sc^2) = 1/sigma2,  sd^2 = 2 sc^2 + 1 - mean^2."""
    from scipy.special import ive

    sc2 = sigma2 / 2
    x = 1.0 / sigma2
    mean = math.sqrt(sc2) * math.sqrt(math.pi / 2) * ((1 + x) * ive(0, x / 2) + x * ive(1, x / 2))
    cos = math.sqrt(math.pi * x) / 2 * (ive(0, x / 2) + ive(1, x / 2))
    return {"mean": float(mean), "sd": math.sqrt(max(2 * sc2 + 1 - mean**2, 0.0)), "cos": float(cos)}


def q_star_fg(X: float, psi: float) -> float:
    """Continuous root q of (q - 1) sqrt(X FG(q)) = psi with FG for prime q (sum_{p|q} = ln q/(q-1)): the spike-sea crossover."""
    from scipy.optimize import brentq

    def g(q):
        fg = math.log(q) - EULER_GAMMA - math.log(2 * math.pi) - math.log(q) / (q - 1)
        return (q - 1) * math.sqrt(X * fg) - psi

    return float(brentq(g, 20.0, 1e6))


def q_star_measured(n: np.ndarray, w: np.ndarray, psi: float, X: int, primes: np.ndarray, window: int = 21) -> tuple[int | None, dict]:
    """First prime q (in the given increasing list) where the running median over `window` primes of
    K(q) = sqrt(mean_a |S(a/q) + psi_c/(q - 1)|^2) (q - 1)/psi exceeds 1.  Returns (q*, {q: K})."""
    K = {}
    for q in primes:
        q = int(q)
        res = np.bincount(n % q, weights=w, minlength=q)
        Sa = q * np.fft.ifft(res)[1:]
        pc = res[1:].sum()
        K[q] = math.sqrt(float(np.mean(np.abs(Sa + pc / (q - 1)) ** 2))) * (q - 1) / psi
    qs = list(K)
    vals = np.array([K[q] for q in qs])
    h = window // 2
    for i in range(h, len(qs) - h):
        if np.median(vals[i - h : i + h + 1]) > 1.0:
            return qs[i], K
    return None, K


def window_max_median(absS: np.ndarray, js, scale: float) -> np.ndarray:
    """Median over non-overlapping windows of 2^j grid points of max |S| / scale, for each j."""
    return np.array([float(np.median(absS.reshape(-1, 1 << j).max(axis=1))) / scale for j in js])


# ================================================================ the singular-series surrogate H_R


def ramanujan_weights(X: int, psi: float, mu: np.ndarray, phi: np.ndarray, Rs) -> dict:
    """Truncated Ramanujan expansion w*_R(n) = (psi/X) sum_{q <= R, mu(q) != 0} mu(q) c_q(n)/phi(q), n = 0..X (w*[0] = 0),
    with c_q(n) = mu(q/g) phi(q)/phi(q/g), g = gcd(n, q).  Built incrementally in q; returns {R: w*_R}."""
    n = np.arange(X + 1, dtype=np.int64)
    acc = np.zeros(X + 1)
    out = {}
    Rs = sorted(Rs)
    for q in range(1, Rs[-1] + 1):
        if mu[q] != 0:
            g = np.gcd(n, q)
            qg = q // g
            acc += mu[q] * (mu[qg] * phi[q] / phi[qg]) / phi[q]
        if q in Rs:
            w = acc * (psi / X)
            w[0] = 0.0
            out[q] = w.copy()
    return out


def hybrid_weights(wstar: np.ndarray, L2: float, seed: int) -> np.ndarray:
    """h = w* + xi, xi iid Gaussian with variance (sum Lambda^2 - sum w*^2)/X per n (so sum h^2 ~ sum Lambda^2)."""
    X = len(wstar) - 1
    var = max((L2 - float(np.sum(wstar**2))) / X, 0.0)
    rng = np.random.default_rng(seed)
    h = wstar.copy()
    h[1:] += math.sqrt(var) * rng.standard_normal(X)
    return h


# ================================================================ M3: tails


def ccdf(values: np.ndarray, ts: np.ndarray) -> np.ndarray:
    """P(values > t) for each t."""
    v = np.sort(values)
    return 1.0 - np.searchsorted(v, ts, side="right") / len(v)


class SincLevels:
    """Measures of {x : |sinc x| > u} (x in units of 1/N), optionally restricted to |x| > c, from one sorted table."""

    def __init__(self, xmax: float = 400.0, dx: float = 1.0 / 400):
        self.dx = dx
        self.x = np.arange(dx / 2, xmax, dx)
        self.v = np.abs(np.sinc(self.x))
        self._cache = {}

    def _sorted(self, c: float):
        if c not in self._cache:
            self._cache[c] = np.sort(self.v[self.x > c])
        return self._cache[c]

    def measure(self, u, c: float = 0.0) -> np.ndarray:
        """|{x in R : |sinc x| > u, |x| > c}| (two-sided)."""
        s = self._sorted(c)
        return 2.0 * self.dx * (len(s) - np.searchsorted(s, np.asarray(u, dtype=float), side="right"))


def spike_ccdf_model(ts, mean_m: float, meas_m: float, N: int, psi: float, mu: np.ndarray, phi: np.ndarray, Q: int,
                     c: float, qcap: int = 20000, levels: SincLevels | None = None) -> np.ndarray:
    """Semi-analytic minor-arc CCDF of t = |S|^2/mean_m: every squarefree q <= qcap contributes phi(q) Dirichlet spikes of
    height psi/phi(q); a spike with q > Q is entirely on the minor arcs, one with q <= Q only beyond |x| = c.
        P(t) = (1/(N |m|)) sum_q phi(q) l_q(u_q),  u_q = sqrt(t mean_m) phi(q)/psi."""
    levels = levels or SincLevels()
    q = np.arange(1, qcap + 1)
    q = q[mu[q] != 0]
    ph = phi[q].astype(float)
    out = np.empty(len(ts))
    for i, t in enumerate(ts):
        u = math.sqrt(t * mean_m) * ph / psi
        lo = q <= Q
        ell = np.where(lo, levels.measure(u, c), levels.measure(u, 0.0))
        out[i] = float(np.sum(ph * ell)) / (N * meas_m)
    return out


def phi_min_above(Q: int, mu: np.ndarray, phi: np.ndarray, qlim: int = 100000) -> tuple[int, int]:
    """(phi_min, q) = smallest phi(q) over squarefree q > Q (q <= qlim), and the smallest such q."""
    q = np.arange(Q + 1, min(qlim, len(mu) - 1) + 1)
    q = q[mu[q] != 0]
    i = int(np.argmin(phi[q]))
    return int(phi[q][i]), int(q[i])


def sup_prediction(Q: int, c: float, psi: float, mu: np.ndarray, phi: np.ndarray) -> dict:
    """sup_m |S| ~ max(psi/phi_min(Q), psi/(pi (c + 1/2))): the smallest-phi spike above Q, or the first arc-edge sidelobe
    of q = 1, 2 (height ~ psi/(pi (c + 1/2)) at |x| ~ c + 1/2)."""
    pm, q = phi_min_above(Q, mu, phi)
    spike, side = psi / pm, psi / (math.pi * (c + 0.5))
    return {"sup2": max(spike, side) ** 2, "spike": spike, "sidelobe": side, "q": q, "phi": pm,
            "active": f"q = {q} (phi = {pm})" if spike >= side else "arc-edge sidelobe of q = 1, 2"}


def sinc4_tail(c: float, xmax: float = 4000.0, n: int = 4_000_001) -> float:
    """T4(c) = 2 int_c^inf sinc^4 (the tail beyond xmax is ~ 1/(3 pi^4 xmax^3), added analytically)."""
    x = np.linspace(c, xmax, n)
    return float(2 * (np.trapezoid(np.sinc(x) ** 4, x) + 1.0 / (8 * math.pi**4 * xmax**3)))   # tail: (3/8)/(3 pi^4 x^3)


def fourth_moment_model(Q: int, c: float, psi: float, N: int, mu: np.ndarray, phi: np.ndarray, qcap: int = 200000) -> float:
    """int_m |S|^4 ~ (psi^4/N) [ (2/3) sum_{q>Q} mu^2/phi^3 + T4(c) sum_{q<=Q} mu^2/phi^3 ]  (int sinc^4 = 2/3)."""
    q = np.arange(1, min(qcap, len(mu) - 1) + 1)
    sq = mu[q] != 0
    inv3 = np.where(sq, 1.0 / phi[q].astype(float) ** 3, 0.0)
    hi = float(inv3[q > Q].sum())
    lo = float(inv3[q <= Q].sum())
    return psi**4 / N * (2.0 / 3.0 * hi + sinc4_tail(c) * lo)


def cell_maxima(values: np.ndarray, minor: np.ndarray, cell: int = 8) -> np.ndarray:
    """Maxima over cells of `cell` consecutive grid points that lie entirely in the minor arcs."""
    v = values.reshape(-1, cell)
    full = minor.reshape(-1, cell).all(axis=1)
    return v[full].max(axis=1)


def hill(x: np.ndarray, k: int) -> float:
    """Hill tail index from the k + 1 largest values: k / sum_{i<k} ln(x_(i)/x_(k)) (0-based order statistics)."""
    top = np.sort(np.partition(x, len(x) - k - 1)[len(x) - k - 1 :])[::-1]
    return float(k / np.sum(np.log(top[:k] / top[k])))


def concentration(x: np.ndarray, fracs) -> list[float]:
    """Share of sum(x) held by the top fraction f of the values, for each f."""
    s = np.sort(x)[::-1]
    c = np.cumsum(s) / s.sum()
    return [float(c[max(int(round(f * len(s))), 1) - 1]) for f in fracs]


# ================================================================ M4: the minor-arc walk


def minor_walk(S: np.ndarray, minor: np.ndarray, X: int) -> tuple[np.ndarray, np.ndarray]:
    """Minor indices k (sorted) and phasors w_k = S(k/M)^2 e(-N k/M)/M (the E0 walk's steps)."""
    from expsum import e_minus

    M = len(S)
    k = np.flatnonzero(minor)
    w = S[k]
    w = w * w
    w *= e_minus(X, k, M)
    w /= M
    return k, w


def noah_only(w: np.ndarray, seed: int = 11) -> np.ndarray:
    """Keep |w|, randomise the phases (destroys order, keeps the heavy tail)."""
    rng = np.random.default_rng(seed)
    return np.abs(w) * np.exp(2j * np.pi * rng.random(len(w)))


def joseph_only(w: np.ndarray) -> np.ndarray:
    """Keep the phases, flatten the amplitudes (w/|w|)."""
    a = np.abs(w)
    return np.where(a > 0, w / np.where(a > 0, a, 1.0), 0.0)


def structure_function(F: np.ndarray, lags, p: float = 2.0, drift: bool = False) -> np.ndarray:
    """S_p(l) = mean_j |F[j + l] - F[j]|^p; with drift=True (p = 2 only) subtract |mean_j (F[j + l] - F[j])|^2."""
    out = np.empty(len(lags))
    for i, L in enumerate(lags):
        d = F[L:] - F[:-L]
        if p == 2.0:
            v = float(np.mean(d.real**2 + d.imag**2))
            if drift:
                v -= abs(complex(np.mean(d))) ** 2
        else:
            v = float(np.mean(np.abs(d) ** p))
        out[i] = v
    return out


def local_hurst(lags, sf: np.ndarray, p: float = 2.0) -> np.ndarray:
    """H_p between consecutive lags: (1/p) dlog S_p / dlog l."""
    return np.diff(np.log(sf)) / np.diff(np.log(np.asarray(lags, dtype=float))) / p


def walk_increments_fft(x: np.ndarray) -> np.ndarray:
    """DFT X_d = sum_k x_k e(-d k/M) of the zero-filled increments (length M)."""
    return np.fft.fft(x)


def spectral_S2(Xhat: np.ndarray, lags, chunk: int = 1 << 20) -> np.ndarray:
    """Cyclic structure function of the walk whose increments have DFT Xhat (length M):
    S2cyc(l) = (1/M) sum_j |F_{j+l} - F_j|^2 = (1/M^2) sum_d |Xhat_d|^2 sin^2(pi d l/M)/sin^2(pi d/M)  (l^2 at d = 0).
    Accumulated over chunks of d to keep memory flat."""
    M = len(Xhat)
    out = np.zeros(len(lags))
    for s0 in range(0, M, chunk):
        d = np.arange(s0, min(M, s0 + chunk), dtype=float)
        P = Xhat[s0 : s0 + len(d)]
        P = P.real**2 + P.imag**2
        s2 = np.sin(np.pi * d / M) ** 2
        for i, L in enumerate(lags):
            ker = np.empty(len(d))
            nz = d != 0
            ker[nz] = np.sin(np.pi * d[nz] * L / M) ** 2 / s2[nz]
            ker[~nz] = float(L) ** 2
            out[i] += float(np.dot(P, ker))
    return out / M**2


def H_triangular(x: float, n: int = 200001) -> float:
    """No-fit crossover shape: S_2(x) ~ int_0^1 (1 - u)^2 sin^2(pi u x)/u^2 du, H = (1/2) log2(S_2(2x)/S_2(x))."""
    u = np.linspace(1e-9, 1.0, n)

    def s2(xx):
        return np.trapezoid((1 - u) ** 2 * np.sin(np.pi * u * xx) ** 2 / u**2, u)

    return 0.5 * math.log2(s2(2 * x) / s2(x))


def H_bridge(l: float, n: float) -> float:
    """Local H between lags l and 2l of a random walk pinned at both ends (n steps): (1/2)[1 + log2((1-2l/n)/(1-l/n))]."""
    return 0.5 * (1 + math.log2((1 - 2 * l / n) / (1 - l / n)))


def order_parameter_entropy(w: np.ndarray, bins: int = 64) -> tuple[float, float]:
    """Order parameter r = |sum w| / sum |w| and the |w|-weighted entropy of arg w (bins, in bits; uniform = log2 bins)."""
    a = np.abs(w)
    r = abs(complex(w.sum())) / float(a.sum())
    h, _ = np.histogram(np.angle(w), bins=bins, range=(-np.pi, np.pi), weights=a)
    p = h / h.sum()
    return float(r), float(-np.sum(p[p > 0] * np.log2(p[p > 0])))


# ================================================================ additions A1, A2, A3


def octave_energy(a2: np.ndarray, X: int, mu: np.ndarray, phi: np.ndarray, jmax: int = 10, half: float = 0.5) -> list[dict]:
    """(1/M) sum |S|^2 over grid points within |N(k/M - a/q)| <= half of every reduced a/q, squarefree q in [2^j, 2^{j+1})."""
    M = len(a2)
    r = int(math.floor(half * M / X))
    off = np.arange(-r, r + 1, dtype=np.int64)
    out = []
    for j in range(jmax + 1):
        tot, inv_phi = 0.0, 0.0
        for q in range(max(1 << j, 1), 1 << (j + 1)):
            if mu[q] == 0:
                continue
            a = np.arange(q, dtype=np.int64)
            a = a[np.gcd(a, q) == 1] if q > 1 else np.array([0])
            centre = np.rint(a * M / q).astype(np.int64)
            k = (centre[:, None] + off[None, :]).ravel() % M
            tot += float(a2[k].sum())
            inv_phi += 1.0 / phi[q]
        out.append({"j": j, "energy": tot / M, "sum_inv_phi": inv_phi})
    return out


def large_value_counts(absS: np.ndarray, N: int, sigmas, cell: int = 8) -> list[int]:
    """#{cells of `cell` grid points with max |S| >= N^sigma}, for each sigma."""
    cm = absS.reshape(-1, cell).max(axis=1)
    return [int(np.sum(cm >= N**s)) for s in sigmas]


def uncertainty_product(kind: str, L: float, nx: int = 20001, nxi: int | None = None) -> dict:
    """Heisenberg-Gabor product for a window w on x = n/N in [0, 1] and its transform in xi = N beta (e(x) = e^{2 pi i x},
    so hbar -> 1/(2 pi) and sigma_x sigma_xi >= 1/(4 pi)): spreads of |w|^2 and of |w^|^2 over |xi| <= L.
    kind = 'sharp' (w = 1: closed forms; sigma_xi grows like sqrt(L)) or 'bump' (the smooth eta of escapes.smooth_bump)."""
    if kind == "sharp":
        sx = 1.0 / math.sqrt(12.0)
        F = float(sinc_F(np.array([L]))[0])
        sxi2 = (L - math.sin(2 * math.pi * L) / (2 * math.pi)) / (math.pi**2) / (2 * F)
        return {"sigma_x": sx, "sigma_xi": math.sqrt(sxi2), "product": sx * math.sqrt(sxi2), "bound": 1 / (4 * math.pi)}
    from escapes import smooth_bump

    x = np.linspace(0.0, 1.0, nx)
    dx = x[1] - x[0]
    wv = smooth_bump(x)
    p = wv**2 / np.trapezoid(wv**2, x)
    m = np.trapezoid(x * p, x)
    sx = math.sqrt(np.trapezoid((x - m) ** 2 * p, x))
    nxi = nxi or int(200 * L) + 1
    xi = np.linspace(-L, L, nxi)
    W = np.empty(nxi)
    for s in range(0, nxi, 1000):
        ph = np.exp(-2j * np.pi * np.outer(xi[s : s + 1000], x))
        W[s : s + 1000] = np.abs(ph @ (wv * dx)) ** 2
    q = W / np.trapezoid(W, xi)
    sxi = math.sqrt(np.trapezoid(xi**2 * q, xi))
    return {"sigma_x": sx, "sigma_xi": sxi, "product": sx * sxi, "bound": 1 / (4 * math.pi)}


# ================================================================ M5: the Mandelbrot set


def cardioid_root(p: int, q: int) -> complex:
    """Root of the p/q bulb on the main cardioid: c = e^{2 pi i t}/2 - e^{4 pi i t}/4, t = p/q."""
    t = p / q
    return complex(np.exp(2j * np.pi * t) / 2 - np.exp(4j * np.pi * t) / 4)


def leading_size(p: int, q: int) -> float:
    """Leading-order bulb radius s = sin(pi p/q)/q^2 (Mandelbrot's empirical N^2 rule, normal-form estimate)."""
    return math.sin(math.pi * p / q) / q**2


def _orbit_derivs(z0: complex, c: complex, q: int):
    """Iterate f_c^q from z0: returns z_q, d z_q/dz, d z_q/dc, the multiplier D = prod 2 z_k and dD/dz, dD/dc."""
    z, dz, dc = z0, 1 + 0j, 0j
    D, dDz, dDc = 1 + 0j, 0j, 0j
    for _ in range(q):
        dDz = 2 * dz * D + 2 * z * dDz
        dDc = 2 * dc * D + 2 * z * dDc
        D = 2 * z * D
        dc = 2 * z * dc + 1
        dz = 2 * z * dz
        z = z * z + c
    return z, dz, dc, D, dDz, dDc


def bulb_centre(p: int, q: int, iters: int = 200) -> tuple[complex, int, bool]:
    """Centre of the period-q component attached at p/q: Newton on f_c^q(0) = 0 seeded at root + 0.9 s along the outward
    normal lambda(1 - lambda)/|...|, lambda = e^{2 pi i p/q}.  Returns (centre, exact period, rotation number ok)."""
    lam = complex(np.exp(2j * np.pi * p / q))
    nrm = lam * (1 - lam)
    c = cardioid_root(p, q) + 0.9 * leading_size(p, q) * nrm / abs(nrm)
    for _ in range(iters):
        z, dz = 0j, 0j
        for _ in range(q):
            dz = 2 * z * dz + 1
            z = z * z + c
        step = z / dz
        c -= step
        if abs(step) < 1e-16:
            break
    orbit = [0j]
    z, per = 0j, 0
    for i in range(1, q + 1):
        z = z * z + c
        if abs(z) < 1e-9:
            per = i
            break
        orbit.append(z)
    rot_ok = False
    if per == q:
        za = (1 - np.sqrt(1 - 4 * c + 0j)) / 2
        ang = np.angle(np.array(orbit) - za)
        pos = np.empty(q, dtype=int)
        pos[np.argsort(ang)] = np.arange(q)
        shift = {(pos[(k + 1) % q] - pos[k]) % q for k in range(q)}
        rot_ok = shift in ({p % q}, {(q - p) % q})
    return c, per, rot_ok


def bulb_continue(p: int, q: int, rho_target: complex, steps: int = 40, iters: int = 60) -> complex:
    """Parameter c where the period-q cycle of the p/q bulb has multiplier rho_target: continuation from the centre
    (rho = 0, cycle through z = 0) by 2x2 complex Newton on f_c^q(z) = z, (f_c^q)'(z) = rho."""
    c, _, _ = bulb_centre(p, q)
    z = 0j
    for k in range(1, steps + 1):
        rho = rho_target * k / steps
        for _ in range(iters):
            zq, dz, dc, D, dDz, dDc = _orbit_derivs(z, c, q)
            F1, F2 = zq - z, D - rho
            a11, a12, a21, a22 = dz - 1, dc, dDz, dDc
            det = a11 * a22 - a12 * a21
            dzs = (F1 * a22 - a12 * F2) / det
            dcs = (a11 * F2 - a21 * F1) / det
            z -= dzs
            c -= dcs
            if abs(dzs) + abs(dcs) < 1e-15:
                break
    return c


def bulb_diam_half(p: int, q: int) -> float:
    """Half the distance from the root to the multiplier -1 point (primary size measure)."""
    return abs(bulb_continue(p, q, -1.0 + 0j) - cardioid_root(p, q)) / 2


def bulb_area_radius(p: int, q: int, n: int = 240) -> float:
    """sqrt(area/pi) of the bulb, boundary from the multiplier circle rho = e^{i phi} (shoelace formula)."""
    pts = [bulb_continue(p, q, complex(np.exp(1j * 2 * np.pi * (k + 0.5) / n)) * (1 - 1e-9)) for k in range(n)]
    xs = np.array([z.real for z in pts])
    ys = np.array([z.imag for z in pts])
    area = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))
    return math.sqrt(area / math.pi)


def farey(n: int) -> list[tuple[int, int]]:
    """Farey sequence F_n on [0, 1] as (a, q) pairs in increasing order."""
    a, b, c, d = 0, 1, 1, n
    out = [(0, 1)]
    while c <= n:
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
        out.append((a, b))
    return out


def mediant_test(n: int, D: int, mu: np.ndarray, phi: np.ndarray) -> dict:
    """For adjacent x < y in F_n, among reduced fractions strictly between them with denominator <= D: is the largest
    bulb (max sin(pi t')/b^2, t' = min(t, 1 - t)) the mediant?  Is the tallest spike (max |mu(b)|/phi(b), ties to the
    smaller b) the mediant?  Also count mediants with non-squarefree denominator."""
    F = farey(n)
    res = {"pairs": 0, "bulb_mediant": 0, "spike_mediant": 0, "mediant_nonsquarefree": 0, "examples": []}
    for (a1, q1), (a2, q2) in zip(F, F[1:]):
        lo, hi = a1 / q1, a2 / q2
        med = (a1 + a2, q1 + q2)
        cands = []
        for b in range(2, D + 1):
            for a in range(math.ceil(lo * b), math.floor(hi * b) + 1):
                if math.gcd(a, b) == 1 and lo < a / b < hi:
                    cands.append((a, b))
        if not cands:
            continue
        res["pairs"] += 1
        tb = max(cands, key=lambda t: (math.sin(math.pi * min(t[0] / t[1], 1 - t[0] / t[1])) / t[1] ** 2, -t[1]))
        ts = max(cands, key=lambda t: (abs(mu[t[1]]) / phi[t[1]], -t[1]))
        res["bulb_mediant"] += tb == med
        res["spike_mediant"] += ts == med
        res["mediant_nonsquarefree"] += mu[med[1]] == 0
        if ts != med and len(res["examples"]) < 4:
            res["examples"].append({"pair": ((a1, q1), (a2, q2)), "mediant": med, "tallest_spike": ts})
    return res


def mandelbrot_escape(re_range, im_range, width: int, height: int, max_iter: int = 400) -> np.ndarray:
    """Escape-time counts (max_iter = interior) on a width x height grid, rows from im_range[1] down; points inside the
    main cardioid and the period-2 disc are marked interior without iterating."""
    xs = np.linspace(re_range[0], re_range[1], width)
    ys = np.linspace(im_range[1], im_range[0], height)
    c = xs[None, :] + 1j * ys[:, None]
    x, y = c.real, c.imag
    qt = (x - 0.25) ** 2 + y**2
    inside = (qt * (qt + x - 0.25) <= y**2 / 4) | ((x + 1) ** 2 + y**2 <= 1 / 16)
    n = np.full(c.shape, max_iter, dtype=np.int32)
    alive = ~inside
    z = np.zeros_like(c)
    for i in range(max_iter):
        if not alive.any():
            break
        z[alive] = z[alive] ** 2 + c[alive]
        esc = alive & (z.real**2 + z.imag**2 > 4.0)
        n[esc] = i
        alive &= ~esc
    return n
