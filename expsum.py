"""Module 1 - the prime exponential sum S(alpha) = sum_{n<=N} Lambda(n) e(n alpha).

Phasor picture: every n with Lambda(n) != 0 is an oscillator of amplitude
Lambda(n) whose phase at "frequency" alpha is 2 pi n alpha.  S(alpha) is their
resultant (the order parameter).

FFT SIGN CONVENTION (resolved and checked in `sign_convention_checks`)
---------------------------------------------------------------------
numpy.fft.fft computes  X[k] = sum_n x[n] exp(-2 pi i n k / M)   (minus sign),
numpy.fft.ifft computes x[n] = (1/M) sum_k X[k] exp(+2 pi i n k / M).
We want S(k/M) = sum_n Lambda(n) exp(+2 pi i n k / M), so

    S(k/M) = M * ifft(Lambda_padded)[k]  =  conj(fft(Lambda_padded)[k]).

Because the weights are real we evaluate it with the half-length rfft:
S(k/M) = conj(rfft[k]) for 0 <= k <= M/2 and S(k/M) = rfft[M-k] above (memory: one
complex array of length M plus one of length M/2).

S(0) = psi(N), S(1/2) = 2 floor(log2 N) log 2 - psi(N) and S(-alpha) = conj S(alpha)
hold under *either* sign, so they are sanity checks, not discriminators.  Two checks
do discriminate: (a) Im S(1/M) against the direct sum  sum Lambda(n) sin(2 pi n / M) > 0,
and (b) the binary grid integral (1/M) sum_k S(k/M)^2 e(-Nk/M), which returns R2(N)
under the right sign and exactly 0 under the wrong one (the wrong sign evaluates the
autoconvolution at M - N, which exceeds 2N, where it vanishes).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------- the grid

def S_grid(lam: np.ndarray, M: int) -> np.ndarray:
    """S(k/M) for k = 0..M-1 by one length-M inverse FFT of the zero-padded weights."""
    if M <= len(lam) - 1:
        raise ValueError("need M > N")
    if M % 2:
        raise ValueError("M must be even")
    x = np.zeros(M)
    x[: len(lam)] = lam
    R = np.fft.rfft(x)                 # R[k] = sum_n x[n] exp(-2 pi i n k / M)
    del x
    h = M // 2
    S = np.empty(M, dtype=complex)
    np.conjugate(R, out=S[: h + 1])    # S(k/M) = conj R[k]        (0 <= k <= M/2)
    S[h + 1 :] = R[1:h][::-1]          # S(k/M) = R[M - k]         (M/2 < k < M)
    return S


def load_or_compute_S(lam: np.ndarray, M: int, tag: str, cache_dir: Path) -> np.ndarray:
    """Cache S on the grid as .npz (uncompressed; a 2^23 complex128 array is 134 MB)."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    f = cache_dir / f"S_{tag}_N{len(lam) - 1}_M{M}.npz"
    if f.exists():
        with np.load(f) as z:
            if int(z["M"]) == M and int(z["N"]) == len(lam) - 1 and math.isclose(
                float(z["wsum"]), float(lam.sum()), rel_tol=1e-12
            ):
                return z["S"]
    S = S_grid(lam, M)
    np.savez(f, S=S, M=M, N=len(lam) - 1, wsum=lam.sum())
    return S


def e_minus(N: int, k: np.ndarray, M: int) -> np.ndarray:
    """e(-N k / M) with the phase reduced exactly in integers: (N k) mod M."""
    k = np.asarray(k, dtype=np.int64)
    return np.exp(-1j * TWO_PI * ((N * k) % M) / M)


def chunked_sum(f, M: int, chunk: int = 1 << 20):
    """sum_{k=0}^{M-1} f(k_slice) accumulated chunk by chunk (keeps memory flat)."""
    tot = 0.0
    for s in range(0, M, chunk):
        k = np.arange(s, min(M, s + chunk), dtype=np.int64)
        tot = tot + f(k)
    return tot


def grid_power(S: np.ndarray, power: int, Ntarget: int, k: np.ndarray | None = None) -> complex:
    """(1/M) sum_{k in set} S(k/M)^power e(-Ntarget k/M); full circle if k is None.

    On the full circle this is *exact* (no quadrature error) as soon as M exceeds
    the largest possible n1+...+n_power minus Ntarget, because the rectangle rule
    integrates trigonometric polynomials of degree < M exactly.
    """
    M = len(S)
    if k is None:
        return chunked_sum(lambda kk: np.sum(S[kk] ** power * e_minus(Ntarget, kk, M)), M) / M
    k = np.asarray(k, dtype=np.int64)
    return np.sum(S[k] ** power * e_minus(Ntarget, k, M)) / M


def grid_mass(S: np.ndarray, k: np.ndarray | None = None) -> float:
    """(1/M) sum_{k in set} |S(k/M)|^2; full circle (Parseval) if k is None."""
    M = len(S)
    if k is None:
        return float(chunked_sum(lambda kk: np.sum(S[kk].real ** 2 + S[kk].imag ** 2), M) / M)
    return float(np.sum(np.abs(S[k]) ** 2) / M)


# ---------------------------------------------------------------- direct evaluation

def nonzero_terms(lam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = np.flatnonzero(lam)
    return n.astype(np.int64), lam[n]


def S_direct(alphas, n: np.ndarray, w: np.ndarray, chunk: int = 48) -> np.ndarray:
    """S(alpha) = sum_j w_j e(n_j alpha) for arbitrary real alphas (phases reduced mod 1)."""
    alphas = np.atleast_1d(np.asarray(alphas, dtype=float))
    nf = n.astype(float)
    out = np.empty(len(alphas), dtype=complex)
    for s in range(0, len(alphas), chunk):
        ph = np.outer(alphas[s : s + chunk], nf)
        ph -= np.rint(ph)
        out[s : s + chunk] = np.exp(1j * TWO_PI * ph) @ w
    return out


def S_near_rational(a: int, q: int, betas, n: np.ndarray, w: np.ndarray, chunk: int = 48) -> np.ndarray:
    """S(a/q + beta): the rational part of every phase is reduced exactly as (n a mod q)/q."""
    betas = np.atleast_1d(np.asarray(betas, dtype=float))
    base = ((n * a) % q) / q
    nf = n.astype(float)
    out = np.empty(len(betas), dtype=complex)
    for s in range(0, len(betas), chunk):
        ph = base[None, :] + np.outer(betas[s : s + chunk], nf)
        ph -= np.rint(ph)
        out[s : s + chunk] = np.exp(1j * TWO_PI * ph) @ w
    return out


def S_rational(a: int, q: int, psi_res: np.ndarray) -> complex:
    """S(a/q) = sum_{r mod q} e(a r / q) psi(N; q, r), exact given the residue sums."""
    r = np.arange(q, dtype=np.int64)
    return complex(np.sum(np.exp(1j * TWO_PI * ((a * r) % q) / q) * psi_res))


# ---------------------------------------------------------------- the model

def dirichlet(beta, N: int) -> np.ndarray:
    """v(beta) = sum_{n=1}^N e(n beta) = e((N+1) beta / 2) sin(pi N beta) / sin(pi beta)."""
    beta = np.asarray(beta, dtype=float)
    s = np.sin(np.pi * beta)
    small = np.abs(s) < 1e-300
    ratio = np.where(small, float(N), np.sin(np.pi * N * beta) / np.where(small, 1.0, s))
    return np.exp(1j * np.pi * (N + 1) * beta) * ratio


def kernel_capture(N: int, c: float, n_pts: int = 200001) -> complex:
    """J(c) = int_{|beta| <= c/N} v(beta)^2 e(-N beta) d beta  (= N-1 when c -> N/2).

    v^2 e(-N beta) = e(beta) sin^2(pi N beta)/sin^2(pi beta); Simpson on a fine grid.
    The tail outside |beta| <= c/N carries about N/(pi^2 c).
    """
    b = np.linspace(-c / N, c / N, n_pts)
    f = dirichlet(b, N) ** 2 * np.exp(-1j * TWO_PI * N * b)
    h = b[1] - b[0]
    wts = np.ones(n_pts)
    wts[1:-1:2] = 4
    wts[2:-1:2] = 2
    return complex(np.sum(wts * f) * h / 3)


# ---------------------------------------------------------------- sign convention

def sign_convention_checks(S: np.ndarray, lam: np.ndarray, N: int, R2: float) -> dict:
    """Pin numpy's FFT sign: returns a dict of check values (see module docstring)."""
    M = len(S)
    psi = math.fsum(lam)
    n, w = nonzero_terms(lam)
    k_test = np.array([1, 7, 12345])
    direct = S_direct(k_test / M, n, w)
    s_half_pred = 2 * int(math.floor(math.log2(N))) * math.log(2) - psi
    kk = np.arange(1, 2049)
    sym = np.max(np.abs(S[M - kk] - np.conj(S[kk]))) / psi
    wrong = grid_power(np.conj(S), 2, N)
    right = grid_power(S, 2, N)
    return {
        "S0": S[0].real,
        "psi": psi,
        "S0_relerr": abs(S[0] - psi) / psi,
        "S_half": S[M // 2].real,
        "S_half_pred": s_half_pred,
        "S_half_abserr": abs(S[M // 2] - s_half_pred),
        "symmetry_maxrel": float(sym),
        "direct_vs_grid_maxrel": float(np.max(np.abs(direct - S[k_test])) / psi),
        "Im_S_1_over_M_grid": S[1].imag,
        "Im_S_1_over_M_direct": direct[0].imag,
        "grid_R2_right_sign": right.real,
        "grid_R2_wrong_sign": wrong.real,
        "R2": R2,
    }
