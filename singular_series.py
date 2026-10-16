"""Ramanujan sums, the binary and ternary singular series, and the Hardy-Littlewood integral.

Phasor picture: c_q(n) is the resultant of the phi(q) unit phasors e(a n / q),
gcd(a, q) = 1 - the "channels" that primes populate equally at the rational a/q.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from sieve import factorize, mobius_phi, primes_upto

# Twin-prime constant C2 = prod_{p>2} (1 - 1/(p-1)^2) (literature value, 20 digits).
C2_KNOWN = 0.66016181584686957393


@lru_cache(maxsize=None)
def _primes(P: int) -> np.ndarray:
    return primes_upto(P)


def C2_product(P: int = 10**7) -> float:
    p = _primes(P)[1:].astype(float)
    return float(np.exp(np.sum(np.log1p(-1.0 / (p - 1.0) ** 2))))


@lru_cache(maxsize=None)
def C3_product(P: int = 10**7) -> float:
    """C3 = prod_p (1 + 1/(p-1)^3) ~ 2.30096 (tail beyond P is < 1/(2 P^2))."""
    p = _primes(P).astype(float)
    return float(np.exp(np.sum(np.log1p(1.0 / (p - 1.0) ** 3))))


# ---------------------------------------------------------------- Ramanujan sums

def ramanujan_direct(q: int, n: int) -> float:
    """c_q(n) = sum_{1<=a<=q, gcd(a,q)=1} e(a n / q), summed as phasors (real part)."""
    a = np.array([a for a in range(1, q + 1) if math.gcd(a, q) == 1], dtype=np.int64)
    return float(np.sum(np.cos(2 * np.pi * ((a * n) % q) / q)))


def ramanujan_divisor(q: int, n: int) -> int:
    """c_q(n) = sum_{d | gcd(q, n)} mu(q/d) d."""
    g = math.gcd(q, n)
    mu, _ = mobius_phi(q)
    return int(sum(mu[q // d] * d for d in range(1, g + 1) if g % d == 0))


def ramanujan_holder(qs: np.ndarray, n: int, mu: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Vectorised c_q(n) via Holder: c_q(n) = mu(q/d) phi(q) / phi(q/d), d = gcd(q, n)."""
    qs = np.asarray(qs, dtype=np.int64)
    d = np.gcd(qs, n)
    return mu[qs // d] * (phi[qs] // phi[qs // d])


# ---------------------------------------------------------------- binary series

def SS2_terms(N: int, Qmax: int, mu=None, phi=None) -> np.ndarray:
    """t[q] = mu(q)^2 c_q(N) / phi(q)^2 for q = 0..Qmax (t[0] = 0)."""
    if mu is None:
        mu, phi = mobius_phi(Qmax)
    q = np.arange(1, Qmax + 1)
    t = np.zeros(Qmax + 1)
    cq = ramanujan_holder(q, N, mu, phi)
    t[1:] = (mu[q] ** 2) * cq / phi[q].astype(float) ** 2
    return t


def SS2_partial(N: int, Qmax: int, mu=None, phi=None) -> np.ndarray:
    """SS2_Q(N) for Q = 0..Qmax (cumulative sums of the terms)."""
    return np.cumsum(SS2_terms(N, Qmax, mu, phi))


def SS2_product(N: int) -> float:
    """SS2(N) = 2 C2 prod_{p | N, p > 2} (p-1)/(p-2) for even N; 0 for odd N."""
    if N % 2:
        return 0.0
    f = 2.0 * C2_KNOWN
    for p in factorize(N):
        if p > 2:
            f *= (p - 1) / (p - 2)
    return f


def SS2_even_array(Nmax: int) -> np.ndarray:
    """SS2(N) for all N = 0..Nmax (zero at odd N), by sieving the odd prime factors."""
    out = np.zeros(Nmax + 1)
    out[0::2] = 2.0 * C2_KNOWN
    for p in _primes(Nmax):
        if p > 2:
            out[p::p] *= (p - 1) / (p - 2)
    out[0] = 0.0
    return out


# ---------------------------------------------------------------- ternary series

def SS3_terms(Np: int, Qmax: int, mu=None, phi=None) -> np.ndarray:
    """t[q] = mu(q) c_q(N') / phi(q)^3."""
    if mu is None:
        mu, phi = mobius_phi(Qmax)
    q = np.arange(1, Qmax + 1)
    t = np.zeros(Qmax + 1)
    t[1:] = mu[q] * ramanujan_holder(q, Np, mu, phi) / phi[q].astype(float) ** 3
    return t


def SS3_partial(Np: int, Qmax: int, mu=None, phi=None) -> np.ndarray:
    return np.cumsum(SS3_terms(Np, Qmax, mu, phi))


def SS3_product(Np: int) -> float:
    """prod_{p | N'} (1 - 1/(p-1)^2) * prod_{p not | N'} (1 + 1/(p-1)^3)."""
    f = C3_product()
    for p in factorize(Np):
        f *= (1 - 1 / (p - 1) ** 2) / (1 + 1 / (p - 1) ** 3)
    return f


# ---------------------------------------------------------------- HL integral

_GL_X, _GL_W = np.polynomial.legendre.leggauss(160)


def hl_integral(Ns) -> np.ndarray:
    """I(N) = int_2^{N-2} dt / (log t log(N - t)), vectorised over N (N >= 8).

    Symmetric about N/2, so I = 2 int_2^{N/2}; substitute t = e^u and use
    160-point Gauss-Legendre in u (the integrand is smooth in u).
    """
    Ns = np.atleast_1d(np.asarray(Ns, dtype=float))
    lo = math.log(2.0)
    hi = np.log(Ns / 2.0)
    half = 0.5 * (hi - lo)
    u = lo + half[:, None] * (_GL_X[None, :] + 1.0)
    t = np.exp(u)
    f = t / (u * np.log(Ns[:, None] - t))
    return 2.0 * half * (f @ _GL_W)
