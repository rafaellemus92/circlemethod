"""pytest suite for compartments.py (the body-water lens) on a small grid (X = 10^4, M = 2^17)."""
from __future__ import annotations

import math

import numpy as np
import pytest

import compartments as cp
from expsum import S_grid
from sieve import mobius_phi, von_mangoldt

X = 10**4
M = 2**17


@pytest.fixture(scope="module")
def tables():
    return mobius_phi(2 * 10**5)


@pytest.fixture(scope="module")
def small():
    lam = von_mangoldt(X)
    S = S_grid(lam, M)
    return {"lam": lam, "a2": S.real**2 + S.imag**2, "psi": float(lam.sum()), "L2": float(np.sum(lam**2))}


def test_class_model_reduces_to_ramanujan(tables):
    mu, phi = tables
    psi = 1000.0
    for q in (1, 2, 3, 4, 6, 9, 10, 12, 30):
        s = cp.class_model_heights(q, psi, 1, {0})
        assert np.allclose(s, mu[q] * psi / phi[q], atol=1e-9 * psi)


def test_class_model_restricted_and_integer_bases():
    psi = 1000.0
    for q in (3, 5, 15):                                    # odd q: primes = 1 (mod 4) with weight 2 look like all primes
        assert np.allclose(cp.class_model_heights(q, psi, 4, {1}, 2.0), cp.class_model_heights(q, psi, 1, {0}), atol=1e-9)
    s4 = cp.class_model_heights(4, psi, 4, {1}, 2.0)        # e(n/4) = i for every n = 1 (mod 4): a spike of height psi
    assert np.allclose(np.abs(s4), psi)
    assert np.allclose(cp.class_model_heights(8, psi, 4, {1}, 2.0), 0.0, atol=1e-9)   # classes 1, 5 mod 8 cancel at odd a/8
    assert np.allclose(cp.class_model_heights(1, psi, 1, {0}, base="integers"), psi)
    for q in (2, 3, 7):                                     # Cramer: no congruence structure, no spike for q > 1
        assert np.allclose(cp.class_model_heights(q, psi, 1, {0}, base="integers"), 0.0, atol=1e-9)
    assert np.allclose(cp.class_model_heights(2, psi, 2, {1}, 2.0, base="integers"), -psi)   # odd numbers: spike at 1/2
    assert np.allclose(cp.class_model_heights(3, psi, 2, {1}, 2.0, base="integers"), 0.0, atol=1e-9)


def test_fraction_windows_membership_and_capture():
    k_all = np.arange(M)
    for q in (1, 2, 7, 30):
        a, k, ok, cap = cp.fraction_windows(M, X, q)
        for i, ai in enumerate(a):
            dist = np.abs(((k_all / M - ai / q + 0.5) % 1.0) - 0.5) * X
            brute = set(np.flatnonzero(dist <= 0.5 + 1e-12))
            assert set((k[i][ok[i]] % M).tolist()) == brute
            xs = X * (k[i][ok[i]] / M - ai / q)
            assert cap[i] == pytest.approx((X / M) * float(np.sum(np.sinc(xs) ** 2)), rel=1e-12)
            assert 0.70 < cap[i] < 0.85


def test_window_energy_matches_capture_times_model(small, tables):
    """For small q the sea is negligible: what a window reads is its capture of the model spike (to 1.5% at X = 10^4)."""
    mu, phi = tables
    a2, psi = small["a2"], small["psi"]
    for q in (1, 2, 3, 5, 6, 7):
        a, k, ok, cap = cp.fraction_windows(M, X, q)
        E = sum(float(a2[k[i][ok[i]] % M].sum()) for i in range(len(a))) / M
        model = np.abs(cp.class_model_heights(q, psi, 1, {0})) ** 2 / X
        assert E == pytest.approx(float(np.sum(cap * model)), rel=0.015)


def test_solve_compartments_recovers_synthetic_values():
    A, F, rho, m_A, m_F, k_A, k_F = 3.0, 7.0, 2.5, 0.01, 0.2, 0.8, 0.76
    E_A = k_A * A + rho * m_A
    E_F = k_F * F + rho * m_F
    T = A + F + rho * (1 - m_A - m_F) + rho * (m_A + m_F)       # spikes plus a uniform sea over the whole circle
    out = cp.solve_compartments(E_A, E_F, T, m_A, m_F, k_A, k_F)
    assert out["A"] == pytest.approx(A) and out["F"] == pytest.approx(F) and out["rho"] == pytest.approx(rho)
    assert out["sea"] == pytest.approx(rho)


def test_compartment_indicator_weighted_capture():
    qs = (3, 5)
    w = {q: np.arange(1, len(cp.fraction_windows(M, X, q)[0]) + 1, dtype=float) for q in qs}
    mask, k = cp.compartment_indicator(M, X, qs, w)
    num = sum(float(np.sum(w[q] * cp.fraction_windows(M, X, q)[3])) for q in qs)
    den = sum(float(np.sum(w[q])) for q in qs)
    assert k == pytest.approx(num / den, rel=1e-12)
    n = sum(int(cp.fraction_windows(M, X, q)[2].sum()) for q in qs)
    assert int(mask.sum()) == n                                    # disjoint windows at these q


def test_prime_predictions_and_limits(small, tables):
    mu, phi = tables
    pr = cp.prime_predictions(X, small["psi"], small["L2"], mu, phi, 75)
    assert pr["L_N"] == pytest.approx(X * small["L2"] / small["psi"] ** 2)
    assert pr["atoms"] == pytest.approx(2 / pr["L_N"])
    assert pr["atoms"] + pr["forest"] + pr["sea"] == pytest.approx(1.0)
    far = cp.asymptotic_shares(1e300, math.sqrt(1e300 / math.log(1e300)))
    assert far["atoms"] < 0.01 and abs(far["forest"] - 0.5) < 0.01 and abs(far["sea"] - 0.5) < 0.02


def test_constants_predictions_and_model_energies(small, tables):
    from scipy.special import sici

    mu, phi = tables
    assert cp.SINC2_HALF == pytest.approx(2 / math.pi * sici(math.pi)[0] - 4 / math.pi**2, rel=1e-9)
    pr = cp.prime_predictions(X, small["psi"], small["L2"], mu, phi, 75)
    assert pr["forest"] == pytest.approx(sum(1 / phi[q] for q in range(3, 76) if mu[q] != 0) / pr["L_N"], rel=1e-12)
    me = cp.model_energies([1, 4, 6], X, 1000.0, 4, {1}, 2.0)
    for q in (1, 4, 6):
        assert np.allclose(me[q], np.abs(cp.class_model_heights(q, 1000.0, 4, {1}, 2.0)) ** 2 / X)


def test_leading_order_law_values():
    a, f, s = cp.leading_order_shares(6 * math.log(10))
    assert (a, f, s) == pytest.approx((0.1561, 0.4346, 0.4093), abs=5e-4)
    a, f, s = cp.leading_order_shares(12 * math.log(10))
    assert (a, f, s) == pytest.approx((0.0751, 0.4500, 0.4749), abs=5e-4)
    a, f, s = cp.leading_order_shares(1000 * math.log(10))
    assert a < 0.001 and abs(f - 0.5) < 0.002 and s > 0.5            # the sea approaches 1/2 from above


def test_divisor_identity_exact(small):
    """Mean over reduced a of |S_red(a/q) - coherent|^2 = (1/phi(q)) sum_{d|q} mu(q/d) d V_q(d), with the variances over n
    coprime to q: an exact identity (Ramanujan sums)."""
    from sieve import factorize

    lam = small["lam"]
    n = np.flatnonzero(lam)
    w = lam[n]
    for q in (6, 30, 35, 42, 105):
        cop = np.gcd(n, q) == 1
        nq, wq = n[cop], w[cop]
        res = np.bincount(nq % q, weights=wq, minlength=q)
        red = np.gcd(np.arange(q), q) == 1
        ps = list(factorize(q))
        coh = (-1) ** len(ps) * res[red].sum() / red.sum()
        centre = float(np.mean(np.abs((q * np.fft.ifft(res))[red] - coh) ** 2))
        V = {1: 0.0}
        for mask in range(1, 1 << len(ps)):
            d = math.prod(p for i, p in enumerate(ps) if mask >> i & 1)
            r_d = np.bincount(nq % d, weights=wq, minlength=d)
            red_d = np.gcd(np.arange(d), d) == 1
            V[d] = float(np.sum((r_d[red_d] - wq.sum() / red_d.sum()) ** 2)) if d > 2 else 0.0
        assert cp.divisor_sea(q, V, ps) == pytest.approx(centre, rel=1e-9)


def test_budget_recovers_synthetic_compartments(tables):
    """Dirichlet spikes of known heights at q = 1, 2, 3 on a uniform sea: the driver's budget() recovers atoms, forest and
    the sea density."""
    import run_compartments as rc

    k = np.arange(M)
    h = {1: 900.0, 2: -900.0, 3: -450.0}
    rho = 2.5e3
    a2 = np.full(M, rho)
    for q, s in h.items():
        for a in ([0] if q == 1 else [b for b in range(q) if math.gcd(b, q) == 1]):
            beta = k / M - a / q
            num = np.sin(np.pi * X * beta) ** 2
            den = np.sin(np.pi * beta) ** 2
            D2 = np.where(den > 1e-30, num / np.where(den > 1e-30, den, 1.0), float(X) ** 2)
            a2 += s**2 * D2 / X**2
    model = {q: np.array([h[q] ** 2 / X] * (1 if q == 1 else sum(math.gcd(b, q) == 1 for b in range(q)))) for q in h}
    b = rc.budget(a2, X, M, model, 3)
    assert b["A"] == pytest.approx((h[1] ** 2 + h[2] ** 2) / X, rel=5e-3)
    assert b["F"] == pytest.approx(2 * h[3] ** 2 / X, rel=5e-3)
    assert b["rho"] == pytest.approx(rho, rel=5e-3)
