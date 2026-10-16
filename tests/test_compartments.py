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
