"""The body-water lens: the energy of S as a set of compartments measured by indicator dilution.

Physiology measures total body water and its compartments by dilution: give a known amount I of an indicator that
spreads only through the compartment of interest, wait for equilibrium, read its concentration C, and set V = I/C
(heavy water for total body water, bromide or inulin for extracellular water, Evans blue for plasma). The compartments
are then nested differences (intracellular = total - extracellular, interstitial = extracellular - plasma), and every
indicator needs a correction for the part of it that strays outside its space.

The circle has the same bookkeeping, with Parseval as the conserved total: (1/M) sum_k |S(k/M)|^2 = sum Lambda^2.
Here the "indicators" are windows |N(alpha - a/q)| <= 1/2 around the reduced fractions a/q of a set of denominators, a
window catches a known share (its capture) of the energy of a Dirichlet spike centred there, and the compartments are

    atoms   the spikes at q = 1, 2 (alpha = 0 and 1/2);
    forest  the spikes with 3 <= q <= Q (Q = q*(X), where spike and sea are equal; see fractals.md, M1);
    sea     everything else: the structureless part, plus the spikes with q > Q that it drowns.

The no-fit spike model is S(a/q + beta) ~ s(a/q) v(beta)/X with v the Dirichlet kernel, whose energy over the circle is
|s|^2/X.  For the primes s(a/q) = mu(q) psi/phi(q) (Ramanujan sums), so the phi(q) spikes of denominator q carry
psi^2 mu^2(q)/(X phi(q)) together; for any set given by allowed residue classes, s(a/q) follows from its class sums.
"""
from __future__ import annotations

import math

import numpy as np

SINC2_HALF = 0.7736950099                  # int_{|x| <= 1/2} sinc^2 = (2/pi) Si(pi) - 4/pi^2 (mean capture over offsets)


# ================================================================ windows (the indicators)
def fraction_windows(M: int, X: int, q: int):
    """Grid points of the windows |X(k/M - a/q)| <= 1/2 around every reduced a/q (exact integer test
    |2X(kq - aM)| <= Mq), with each window's capture under the rectangle rule on its actual grid offset:
        capture(a/q) = (X/M) sum_{k in window} sinc^2(X(k/M - a/q)).
    Returns (a, k, ok, capture): a[i] the numerators, k[i, j] candidate grid indices (unreduced), ok[i, j] membership."""
    a = np.arange(q, dtype=np.int64)
    a = a[np.gcd(a, q) == 1] if q > 1 else np.array([0], dtype=np.int64)
    r = int(math.ceil(0.5 * M / X)) + 1
    off = np.arange(-r, r + 1, dtype=np.int64)
    k = (a * M // q)[:, None] + off[None, :]
    ok = np.abs(2 * X * (k * q - a[:, None] * M)) <= M * q
    x = X * (k / M - a[:, None] / q)
    cap = (X / M) * np.sum(np.where(ok, np.sinc(x) ** 2, 0.0), axis=1)
    return a, k, ok, cap


def compartment_indicator(M: int, X: int, qs, weights=None):
    """Union mask of the windows of the denominators in qs, and the mean capture of those windows weighted by the model
    energy of their spikes: weights[q] is an array over the reduced a (or a scalar for all of them).  Without weights the
    returned capture is nan."""
    mask = np.zeros(M, dtype=bool)
    num, den = 0.0, 0.0
    for q in qs:
        a, k, ok, cap = fraction_windows(M, X, q)
        mask[k[ok] % M] = True
        if weights is not None:
            w = np.asarray(weights[q], dtype=float)
            w = np.broadcast_to(w, cap.shape) if w.ndim == 0 else w
            num += float(np.sum(w * cap))
            den += float(np.sum(w))
    return mask, (num / den if den > 0 else float("nan"))


def solve_compartments(E_A: float, E_F: float, T: float, m_A: float, m_F: float, k_A: float, k_F: float) -> dict:
    """Three-indicator solution for the spike energies A (atoms), F (forest) and a uniform sea density rho:
        E_A = k_A A + rho m_A,   E_F = k_F F + rho m_F,   T - E_A - E_F = (1 - k_A) A + (1 - k_F) F + rho (1 - m_A - m_F),
    i.e. what each window reads is its capture of its own spikes plus the sea it covers, and the rest of the circle holds
    the spikes' leakage and the remaining sea.  The sea compartment is then T - A - F (it includes the drowned spikes)."""
    m_O = 1.0 - m_A - m_F
    Mx = np.array([[k_A, 0.0, m_A], [0.0, k_F, m_F], [1.0 - k_A, 1.0 - k_F, m_O]])
    A, F, rho = np.linalg.solve(Mx, [E_A, E_F, T - E_A - E_F])
    return {"A": float(A), "F": float(F), "rho": float(rho), "sea": float(T - A - F), "T": float(T)}


# ================================================================ the no-fit congruence model of spike heights
def class_model_heights(q: int, psi: float, D: int, allowed, scale: float = 1.0, base: str = "primes") -> np.ndarray:
    """Spike heights s(a/q) for every reduced a, for a set built from `base` restricted to the classes `allowed` mod D and
    weighted by `scale` (total weight psi before the restriction).

    base = "primes": the weight is equidistributed over the phi(L) reduced classes mod L = lcm(D, q) (psi/phi(L) each);
    base = "integers": over all L classes (psi/L each; the Cramer model, which has no congruence structure).  The set keeps
    the classes b with b mod D in `allowed`, each with weight scale * (its share), and
        s(a/q) = sum_{b mod L} w_b e(ab/q).
    For primes with D = 1 this is mu(q) psi/phi(q) (the Ramanujan sum); for integers with D = 1 it is psi at q = 1 and 0
    otherwise.  Returns an array over the reduced a in increasing order."""
    L = D * q // math.gcd(D, q)
    b = np.arange(L, dtype=np.int64)
    if base == "primes":
        keep = (np.gcd(b, L) == 1) & np.isin(b % D, list(allowed))
        share = psi / _phi(L)
    elif base == "integers":
        keep = np.isin(b % D, list(allowed))
        share = psi / L
    else:
        raise ValueError(base)
    w = np.zeros(L)
    w[keep] = scale * share
    res = np.bincount(b % q, weights=w, minlength=q)          # fold the classes mod L onto residues mod q
    Sa = q * np.fft.ifft(res)                                  # s(a/q) = sum_c res_c e(ac/q) for every a
    a = np.arange(q)
    red = np.gcd(a, q) == 1 if q > 1 else np.array([True])
    return Sa[red]


def _phi(n: int) -> int:
    out, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            out -= out // p
        p += 1
    if m > 1:
        out -= out // m
    return out


def model_energies(qs, X: int, psi: float, D: int, allowed, scale: float = 1.0, base: str = "primes") -> dict:
    """Per-q model: the kernel energy |s(a/q)|^2 / X of every reduced fraction a/q (array over a), for q in qs."""
    return {q: np.abs(class_model_heights(q, psi, D, allowed, scale, base)) ** 2 / X for q in qs}


# ================================================================ predictions and limits
def prime_predictions(X: int, psi: float, L2: float, mu: np.ndarray, phi: np.ndarray, Q: float) -> dict:
    """No-fit shares of sum Lambda^2 for the primes: atoms 2/L_N, forest (sum_{3<=q<=Q} mu^2/phi)/L_N, with
    L_N = X sum Lambda^2 / psi^2 (the 'volume of distribution' of the energy in units of the spike at 0, psi^2/X)."""
    LN = X * L2 / psi**2
    sF = sum(1.0 / phi[q] for q in range(3, int(Q) + 1) if mu[q] != 0)
    return {"L_N": LN, "atoms": 2.0 / LN, "forest": sF / LN, "sea": 1.0 - (2.0 + sF) / LN}


def leading_order_shares(lnX: float) -> tuple[float, float, float]:
    """Leading-order (atoms, forest, sea) at ln X = lnX: q* from the fixed point q = sqrt(X / FG(q)) with the prime-q form
    FG(q) = ln q - gamma - ln 2 pi, sum_{q<=Q} mu^2/phi ~ ln Q + 1.3326, and L_N ~ ln X - 1.  This is a model (the FG
    asymptotic and the spike model up to q*), not a theorem.  The atoms fall like 2/ln X; the forest and the sea approach
    1/2 only like ln ln X / (2 ln X), and the sea passes 1/2 near X ~ 10^38 and returns to it from above."""
    lq = 0.5 * lnX
    for _ in range(80):
        fg = max(lq - 0.5772156649 - math.log(2 * math.pi), 0.5)
        lq = 0.5 * (lnX - math.log(fg))
    LN = lnX - 1.0
    f = (lq + 1.3326 - 2.0) / LN
    return 2.0 / LN, f, 1.0 - 2.0 / LN - f


def divisor_sea(q: int, V: dict, factors) -> float:
    """Mean over reduced a of |S(a/q) - coherent|^2 from the residue-class variances V[d] of the divisors d | q:
        (1/phi(q)) sum_{d|q} mu(q/d) d V[d]        (q squarefree; V[1] = V[2] = 0, a single reduced class)
    because sum_{a reduced} e(a(r - s)/q) = c_q(r - s) = sum_{d|q, d|r-s} d mu(q/d)."""
    ps = list(factors)
    tot = 0.0
    for mask in range(1 << len(ps)):
        d = 1
        for i, p in enumerate(ps):
            if mask >> i & 1:
                d *= p
        sign = (-1) ** (len(ps) - bin(mask).count("1"))
        tot += sign * d * V.get(d, 0.0)
    phq = 1
    for p in ps:
        phq *= p - 1
    return tot / phq


def asymptotic_shares(X: float, Q: float, LN: float | None = None) -> dict:
    """Leading-order shares with sum_{q<=Q} mu^2/phi ~ ln Q + 1.3326 and L_N ~ ln X - 1: atoms 2/L_N -> 0,
    forest (ln Q + 1.3326 - 2)/L_N -> 1/2 when Q ~ q*(X) ~ sqrt(X / FG), sea -> 1/2."""
    LN = LN if LN is not None else math.log(X) - 1.0
    f = (math.log(Q) + 1.3326 - 2.0) / LN
    return {"atoms": 2.0 / LN, "forest": f, "sea": 1.0 - 2.0 / LN - f}
