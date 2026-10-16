"""Module 0 - sieves, the von Mangoldt weights, and the exact counts.

Every number the circle method is later compared against is computed here by
plain counting, with no Fourier analysis: pi(N), psi(N), sum Lambda^2, the
weighted binary count R2(N), the unweighted pair count G(N), and residue-class
sums psi(N; q, r).
"""
from __future__ import annotations

import math

import numpy as np


# ---------------------------------------------------------------- sieves

def prime_sieve(N: int) -> np.ndarray:
    """Boolean array is_p[0..N] (Eratosthenes)."""
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if is_p[p]:
            is_p[p * p :: p] = False
    return is_p


def primes_upto(N: int) -> np.ndarray:
    return np.flatnonzero(prime_sieve(N))


def von_mangoldt(N: int) -> np.ndarray:
    """Lambda(n) for n = 0..N as float64; Lambda(p^k) = log p, else 0.

    Prime powers p^k (k >= 2) are included.  Their total weight is
    sum_{k>=2} theta(N^(1/k)) = O(sqrt N): about 1.0e3 out of psi(10^6) ~ 1e6.
    """
    lam = np.zeros(N + 1)
    p = primes_upto(N)
    logp = np.log(p.astype(float))
    lam[p] = logp
    pk, keep = p.copy(), np.ones(len(p), dtype=bool)
    while True:
        keep &= pk <= N // p          # p^(k+1) <= N, without overflow
        if not keep.any():
            break
        pk = np.where(keep, pk * p, pk)
        lam[pk[keep]] = logp[keep]
    return lam


def prime_power_weight(lam: np.ndarray, is_p: np.ndarray) -> float:
    """Total Lambda-weight carried by proper prime powers p^k, k >= 2."""
    n = np.flatnonzero(lam)
    return float(lam[n[~is_p[n]]].sum())


def mobius_phi(Q: int) -> tuple[np.ndarray, np.ndarray]:
    """(mu, phi) as int64 arrays indexed 0..Q (index 0 unused, set to 0)."""
    mu = np.ones(Q + 1, dtype=np.int64)
    phi = np.arange(Q + 1, dtype=np.int64)
    for p in primes_upto(Q):
        mu[p::p] *= -1
        mu[p * p :: p * p] = 0
        phi[p::p] -= phi[p::p] // p
    mu[0] = 0
    phi[0] = 0
    return mu, phi


def factorize(n: int) -> dict[int, int]:
    """Trial-division factorisation (fine for n <= 10^12)."""
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


# ---------------------------------------------------------------- exact counts

def chebyshev_psi(lam: np.ndarray) -> float:
    return math.fsum(lam)


def sum_lambda_sq(lam: np.ndarray) -> float:
    return math.fsum(lam[lam > 0] ** 2)


def R2_direct(lam: np.ndarray, N: int) -> float:
    """R2(N) = sum_{n1+n2=N} Lambda(n1) Lambda(n2) by direct summation (fsum)."""
    if N >= len(lam) + 1:
        raise ValueError("Lambda array too short")
    n1 = np.arange(1, N)
    n1 = n1[(lam[n1] > 0)]
    n1 = n1[N - n1 < len(lam)]
    return math.fsum(lam[n1] * lam[N - n1])


def goldbach_G(N: int, is_p: np.ndarray) -> int:
    """G(N) = #{p <= p' prime : p + p' = N} by direct counting."""
    p = np.flatnonzero(is_p[: N // 2 + 1])
    return int(np.count_nonzero(is_p[N - p]))


def autoconvolve(x: np.ndarray, length: int | None = None) -> np.ndarray:
    """r[m] = sum_j x[j] x[m-j] for m = 0..2(len(x)-1), via zero-padded rfft.

    Zero-padding to L >= 2 len(x) - 1 makes the circular convolution linear:
    no aliasing.  Float64 error is ~1e-16 * sum(x^2) per entry.
    """
    n = len(x)
    L = length or 1 << int(math.ceil(math.log2(2 * n - 1)))
    assert L >= 2 * n - 1
    X = np.fft.rfft(x, L)
    return np.fft.irfft(X * X, L)[: 2 * n - 1]


def R3_exact(lam: np.ndarray, Np: int, r2: np.ndarray | None = None) -> float:
    """R3(N') = sum_m R2(m) Lambda(N'-m), with R2(m) for all m by FFT autoconvolution."""
    if r2 is None:
        r2 = autoconvolve(lam)
    N = len(lam) - 1
    m = np.arange(max(0, Np - N), min(Np, len(r2) - 1) + 1)
    return math.fsum(r2[m] * lam[Np - m])


def goldbach_comet(Nmax: int) -> tuple[np.ndarray, np.ndarray]:
    """G(N) for every even N in [4, Nmax], via FFT autoconvolution of the prime indicator.

    Returns (N_values, G_values).  Ordered representations r(N) = sum 1_{p}(n) 1_{p}(N-n)
    count (p, p') and (p', p) separately and p = p' once, so G = (r + 1_{N/2 prime}) / 2.
    """
    is_p = prime_sieve(Nmax)
    r = np.rint(autoconvolve(is_p.astype(float))).astype(np.int64)
    Ns = np.arange(4, Nmax + 1, 2)
    G = (r[Ns] + is_p[Ns // 2]) // 2
    return Ns, G


# ---------------------------------------------------------------- residue classes

def psi_residues(lam: np.ndarray, q: int) -> np.ndarray:
    """psi(N; q, r) = sum_{n <= N, n = r mod q} Lambda(n), for r = 0..q-1."""
    n = np.flatnonzero(lam)
    return np.bincount(n % q, weights=lam[n], minlength=q)


def theta_residues(lam: np.ndarray, is_p: np.ndarray, q: int) -> np.ndarray:
    """theta(N; q, r): the same with primes only (no prime powers)."""
    n = np.flatnonzero(is_p[: len(lam)])
    return np.bincount(n % q, weights=lam[n], minlength=q)


# ---------------------------------------------------------------- random model

def cramer_weights(N: int, seed: int = 20261016) -> np.ndarray:
    """Cramer model: n >= 3 is a member with probability 1/log n; members weigh log n.

    Returns w[0..N] with w[n] = log n for members, 0 otherwise.  E[w(n)] = 1 exactly,
    so E[sum w] = N - 2 and the model's 'singular series' is 1: no congruence
    obstruction, no phase-locking at any a/q.
    """
    rng = np.random.default_rng(seed)
    n = np.arange(3, N + 1)
    logn = np.log(n)
    member = rng.random(len(n)) < 1.0 / logn
    w = np.zeros(N + 1)
    w[n[member]] = logn[member]
    return w
