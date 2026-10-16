"""The fractal lens: Mandelbrot's principles used as measurements on the circle-method objects at X = 10^6.

Mandelbrot's programme, applied literally (as measurements, not metaphors):

  self-similarity with cutoffs   -- does the picture repeat under zoom, and between which scales?
  roughness has a dimension      -- singularity spectrum f(alpha), generalized dimensions D_s, Hurst exponent H;
  mild vs wild randomness        -- exponential tails vs Pareto power laws;
  rational hierarchies           -- the Mandelbrot set's p/q bulbs and the primes' a/q spikes share the Farey order.

Everything here is a finite-N statistic of one computation; the scaling ranges are bounded by 1/N (the coherence
length of S) below and by 1 (the circle) or sqrt(N/log N) (where the spikes drown in the sea) above.  Nothing here
bears on a proof.
"""
from __future__ import annotations

import math

import numpy as np

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


def box_masses(mu: np.ndarray, j: int) -> np.ndarray:
    """Masses of the boxes of 2^j consecutive grid points (len(mu) must be a power of two >= 2^j)."""
    return mu.reshape(-1, 1 << j).sum(axis=1)


def falpha_histogram(boxes: np.ndarray, eps: float, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Histogram estimate of the singularity spectrum at one box size eps.

    alpha_box = log(mass) / log(eps);  N(alpha) = #boxes per bin;  f(alpha) = log N(alpha) / log(1/eps).
    Returns (bin centres, f (nan where the bin is empty), counts).
    """
    b = boxes[boxes > 0]
    al = np.log(b) / math.log(eps)
    h, _ = np.histogram(al, bins)
    f = np.full(len(h), np.nan)
    f[h > 0] = np.log(h[h > 0]) / math.log(1.0 / eps)
    return 0.5 * (bins[1:] + bins[:-1]), f, h


def chhabra_jensen(mu: np.ndarray, js, s_list) -> dict:
    """Direct (Chhabra-Jensen 1989) estimate of alpha(s), f(s) from the s-weighted box measures, scale by scale.

    For each box size eps_j = 2^j / M:  m_i(s) = mu_i^s / sum mu^s,  A_j(s) = sum m_i log mu_i,  F_j(s) = sum m_i log m_i.
    alpha(s), f(s) are the slopes of A_j and F_j against log eps; returned as local slopes between consecutive j
    (difference quotients, not fits).
    """
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
    """Local generalized dimensions D_s = (1/(s-1)) dlog Z_s / dlog eps between consecutive box sizes, Z_s = sum mu_box^s."""
    M = len(mu)
    Z = {s: [] for s in s_list}
    for j in js:
        b = box_masses(mu, j)
        b = b[b > 0]
        for s in s_list:
            Z[s].append(float(np.sum(b**s)))
    le = np.log(np.array([2.0**j / M for j in js]))
    return {s: np.diff(np.log(Z[s])) / np.diff(le) / (s - 1.0) for s in s_list}


# ================================================================ tails


def ccdf(values: np.ndarray, ts: np.ndarray) -> np.ndarray:
    """P(values > t) for each t (values need not be sorted)."""
    v = np.sort(values)
    return 1.0 - np.searchsorted(v, ts, side="right") / len(v)


# ================================================================ roughness


def structure_function(F: np.ndarray, lags) -> np.ndarray:
    """<|F[j + l] - F[j]|^2> over all j, for each lag l (F complex or real)."""
    out = np.empty(len(lags))
    for i, L in enumerate(lags):
        d = F[L:] - F[:-L]
        out[i] = float(np.mean(d.real**2 + d.imag**2))
    return out


def local_hurst(lags, sf: np.ndarray) -> np.ndarray:
    """H between consecutive lags: (1/2) dlog S_2 / dlog l."""
    return 0.5 * np.diff(np.log(sf)) / np.diff(np.log(np.asarray(lags, dtype=float)))


# ================================================================ self-similarity of the spikes


def spike_collapse(lam: np.ndarray, qs, mu: np.ndarray, phi: np.ndarray, max_a: int = 400) -> list[dict]:
    """|S(a/q)| phi(q) / psi for reduced a mod q (exact, from the residue sums psi(N; q, r)); <= max_a values of a per q.

    For squarefree q the major-arc model predicts exactly 1 (the spike at a/q is the spike at 0 scaled by
    mu(q)/phi(q)); the sea adds a fluctuation of rms ~ sqrt(sum Lambda^2) phi(q) / psi.
    """
    from expsum import S_rational
    from sieve import psi_residues

    psi = float(lam.sum())
    sea = math.sqrt(float(np.sum(lam**2)))
    out = []
    for q in qs:
        pr = psi_residues(lam, q)
        a_s = [a for a in range(1, max(q, 2)) if math.gcd(a, q) == 1] if q > 1 else [0]
        if len(a_s) > max_a:
            a_s = a_s[:: int(math.ceil(len(a_s) / max_a))]
        r = np.array([abs(S_rational(a, q, pr)) * phi[q] / psi for a in a_s])
        out.append({"q": int(q), "mu": int(mu[q]), "phi": int(phi[q]), "n_a": len(a_s), "mean": float(r.mean()),
                    "sd": float(r.std()), "rms_dev": float(np.sqrt(np.mean((r - 1.0) ** 2))),
                    "min": float(r.min()), "max": float(r.max()), "sea_over_spike": sea * phi[q] / psi})
    return out


# ================================================================ the Mandelbrot set


def cardioid_root(p: int, q: int) -> complex:
    """Root of the p/q bulb on the main cardioid: c = e^{2 pi i t}/2 - e^{4 pi i t}/4, t = p/q."""
    t = p / q
    return complex(np.exp(2j * np.pi * t) / 2 - np.exp(4j * np.pi * t) / 4)


def bulb_centre(p: int, q: int, iters: int = 200) -> tuple[complex, int]:
    """Centre of the period-q hyperbolic component attached at p/q: Newton on f_c^q(0) = 0 seeded just outside the root.

    Returns (centre, exact period found by iterating 0 until it returns to 0; 0 if it does not within q steps).
    """
    root = cardioid_root(p, q)
    c = root + 0.5 / q**2 * (root / abs(root))
    for _ in range(iters):
        z, dz = 0j, 0j
        for _ in range(q):
            dz = 2 * z * dz + 1
            z = z * z + c
        step = z / dz
        c -= step
        if abs(step) < 1e-15:
            break
    z, per = 0j, 0
    for i in range(1, q + 1):
        z = z * z + c
        if abs(z) < 1e-9:
            per = i
            break
    return c, per


def mandelbrot_escape(re_range, im_range, width: int, height: int, max_iter: int = 400) -> np.ndarray:
    """Escape-time counts (max_iter = did not escape) on a width x height grid; rows run from im_range[1] down."""
    xs = np.linspace(re_range[0], re_range[1], width)
    ys = np.linspace(im_range[1], im_range[0], height)
    c = xs[None, :] + 1j * ys[:, None]
    z = np.zeros_like(c)
    n = np.full(c.shape, max_iter, dtype=np.int32)
    alive = np.ones(c.shape, dtype=bool)
    for i in range(max_iter):
        z[alive] = z[alive] ** 2 + c[alive]
        esc = alive & (z.real**2 + z.imag**2 > 4.0)
        n[esc] = i
        alive &= ~esc
        if not alive.any():
            break
    return n
