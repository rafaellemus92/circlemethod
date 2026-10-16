"""pytest suite for fractals.py: the estimators on measures with known spectra, then the circle-method objects on a small grid."""
from __future__ import annotations

import math

import numpy as np
import pytest

import fractals as fr


def _cascade(p: float, levels: int) -> np.ndarray:
    """Binomial (p-model) cascade on 2^levels cells: the textbook multifractal with a closed-form spectrum."""
    mu = np.array([1.0])
    for _ in range(levels):
        mu = np.stack([mu * p, mu * (1 - p)], axis=1).ravel()
    return mu


def test_estimators_on_uniform_measure():
    u = np.full(2**14, 2.0**-14)
    for s, d in fr.partition_dims(u, range(0, 8), [0.5, 2.0, 3.0]).items():
        assert np.allclose(d, 1.0, atol=1e-12)
    cj = fr.chhabra_jensen(u, range(0, 8), [-3.0, 0.5, 2.0])
    for s in cj:
        assert np.allclose(cj[s]["alpha"], 1.0, atol=1e-12) and np.allclose(cj[s]["f"], 1.0, atol=1e-12)


def test_estimators_on_binomial_cascade():
    p = 0.7
    mu = _cascade(p, 16)
    D2 = fr.partition_dims(mu, range(0, 10), [2.0])[2.0]
    assert np.allclose(D2, -math.log2(p * p + (1 - p) ** 2), atol=1e-12)
    cj = fr.chhabra_jensen(mu, range(0, 10), [40.0, -40.0, 1.0])
    assert cj[40.0]["alpha"][0] == pytest.approx(-math.log2(p), abs=1e-6)          # alpha_min
    assert cj[-40.0]["alpha"][0] == pytest.approx(-math.log2(1 - p), abs=1e-6)     # alpha_max
    H = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))                          # information dimension
    assert cj[1.0]["f"][0] == pytest.approx(H, abs=1e-12) and cj[1.0]["alpha"][0] == pytest.approx(H, abs=1e-12)


def test_falpha_histogram_on_cascade():
    """At level j the cascade has C(n, k) boxes of mass p^k (1-p)^(n-k): the histogram must count exactly those."""
    p, n = 0.7, 12
    mu = _cascade(p, n)
    eps = 2.0**-n
    k = np.arange(n + 1)
    al = (k * math.log(p) + (n - k) * math.log(1 - p)) / math.log(eps)
    bins = np.sort(np.r_[al - 1e-9, al + 1e-9])
    _, f, h = fr.falpha_histogram(mu, eps, bins)
    counts = h[::2][np.argsort(np.argsort(al))]
    assert sorted(counts.tolist()) == sorted(math.comb(n, int(x)) for x in k)


def test_ccdf_and_structure_function():
    v = np.arange(1, 101, dtype=float)
    assert fr.ccdf(v, np.array([0.0, 50.0, 100.0])).tolist() == [1.0, 0.5, 0.0]
    rng = np.random.default_rng(3)
    walk = np.cumsum(rng.standard_normal(2**18))                                  # Brownian: H = 1/2
    lags = [2**j for j in range(2, 9)]                                            # >= 1000 independent increments
    H = fr.local_hurst(lags, fr.structure_function(walk, lags))
    assert np.all(np.abs(H - 0.5) < 0.03)
    line = np.arange(2**12, dtype=float)                                          # ballistic: H = 1
    assert np.allclose(fr.local_hurst(lags[:5], fr.structure_function(line, lags[:5])), 1.0)


def test_remove_atoms():
    X, M = 1000, 2**14
    a2 = np.ones(M)
    out = fr.remove_atoms(a2, X, radius=20.0)
    r = int(20 * M / X)
    assert out[: r + 1].sum() == 0 and out[M - r :].sum() == 0 and out[M // 2 - r : M // 2 + r + 1].sum() == 0
    assert out.sum() == M - (r + 1) - r - (2 * r + 1)


@pytest.mark.parametrize("p, q, centre", [(1, 2, -1.0 + 0j), (1, 3, -0.1225611668766536 + 0.7448617666197442j)])
def test_bulb_centres_known_values(p, q, centre):
    c, per = fr.bulb_centre(p, q)
    assert per == q and abs(c - centre) < 1e-12


def test_bulb_size_law():
    """|centre - root| q^2 / sin(pi p/q) stays near 1 for the main-cardioid bulbs (the 1/q^2 law)."""
    for p, q in [(1, 4), (2, 5), (1, 7), (3, 8), (1, 13)]:
        c, per = fr.bulb_centre(p, q)
        assert per == q
        assert 0.95 < abs(c - fr.cardioid_root(p, q)) * q * q / math.sin(math.pi * p / q) < 1.1
