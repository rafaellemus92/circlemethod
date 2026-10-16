"""Escape routes from the binary L^2 wall: relaxations under which the circle method's
magnitude-only minor-arc inequality *does* close at X = 10^6, plus controls that do not.

For one even N the binary wall is the statement that P_m = int_m |S|^2 exceeds R2(N): any bound that
uses only |S| on the minor arcs is consistent with |I_m(N)| = P_m (phases aligned against us), so it
cannot show R2(N) > 0.  Beating that for a single N means proving the phase cancellation -- the open
problem; nothing here does that.  What we test instead, numerically and for these N only, is where the
same inequality closes once the question is relaxed in one of the ways the literature relaxes it:

  E0  phase scrambling (no-go control): same |S| on the minor arcs, different phases;
  E1  average over N  -> exceptional set   (Chudakov, van der Corput, Estermann 1937-38; Montgomery-Vaughan 1975);
  E2  add K powers of 2 -> Linnik's problem (Linnik 1951-53; Gallagher 1975; Heath-Brown-Puchta 2002);
  E3  one prime -> almost-prime            (direction of Renyi 1948 and Chen 1973 -- which need the sieve);
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


def kfold_power_counts(L: int, K: int, length: int) -> np.ndarray:
    """c_K[s] = #{(k_1..k_K) in {1..L}^K : 2^{k_1} + ... + 2^{k_K} = s} (ordered), s < length."""
    u = np.zeros(length, dtype=np.int64)
    for j in range(1, L + 1):
        if (1 << j) < length:
            u[1 << j] = 1
    c = np.zeros(length, dtype=np.int64)
    c[0] = 1
    for _ in range(K):
        nxt = np.zeros(length, dtype=np.int64)
        for j in range(1, L + 1):
            s = 1 << j
            if s < length:
                nxt[s:] += c[: length - s]
        c = nxt
    return c
