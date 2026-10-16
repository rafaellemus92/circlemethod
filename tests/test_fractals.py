"""pytest suite for fractals.py: estimators on measures with known spectra, the Farey count, the collapse models, the
tail model, the walk statistics, the surrogate and the Mandelbrot-set routines; circle-method objects on a small grid
(X = 10^4, M = 2^17, (Q, c) = (20, 10))."""
from __future__ import annotations

import math

import numpy as np
import pytest

import fractals as fr
from arcs import MajorArcs
from expsum import S_grid, S_rational
from sieve import mobius_phi, primes_upto, von_mangoldt

X = 10**4
M = 2**17


@pytest.fixture(scope="module")
def tables():
    mu, phi = mobius_phi(2 * 10**5)
    return mu, phi


@pytest.fixture(scope="module")
def small():
    lam = von_mangoldt(X)
    S = S_grid(lam, M)
    return {"lam": lam, "S": S, "psi": float(lam.sum()), "L2": float(np.sum(lam**2))}


def _cascade(p: float, levels: int) -> np.ndarray:
    """Binomial (p-model) cascade on 2^levels cells: the textbook multifractal with a closed-form spectrum."""
    mu = np.array([1.0])
    for _ in range(levels):
        mu = np.stack([mu * p, mu * (1 - p)], axis=1).ravel()
    return mu


# ------------------------------------------------------------------ estimators on known measures
def test_estimators_on_uniform_measure():
    u = np.full(2**14, 2.0**-14)
    for s, d in fr.partition_dims(u, range(0, 8), [0.5, 1.0, 2.0, 3.0]).items():
        assert np.allclose(d, 1.0, atol=1e-12)
    cj = fr.chhabra_jensen(u, range(0, 8), [-3.0, 0.5, 2.0])
    for s in cj:
        assert np.allclose(cj[s]["alpha"], 1.0, atol=1e-12) and np.allclose(cj[s]["f"], 1.0, atol=1e-12)


def test_estimators_on_binomial_cascade():
    p = 0.7
    mu = _cascade(p, 16)
    D = fr.partition_dims(mu, range(0, 10), [1.0, 2.0])
    H = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))                          # information dimension
    assert np.allclose(D[2.0], -math.log2(p * p + (1 - p) ** 2), atol=1e-12)
    assert np.allclose(D[1.0], H, atol=1e-12)
    cj = fr.chhabra_jensen(mu, range(0, 10), [40.0, -40.0, 1.0])
    assert cj[40.0]["alpha"][0] == pytest.approx(-math.log2(p), abs=1e-6)          # alpha_min
    assert cj[-40.0]["alpha"][0] == pytest.approx(-math.log2(1 - p), abs=1e-6)     # alpha_max
    assert cj[1.0]["f"][0] == pytest.approx(H, abs=1e-12) and cj[1.0]["alpha"][0] == pytest.approx(H, abs=1e-12)


def test_falpha_histogram_on_cascade():
    """At level n the cascade has C(n, k) boxes of mass p^k (1-p)^(n-k): the histogram must count exactly those."""
    p, n = 0.7, 12
    mu = _cascade(p, n)
    eps = 2.0**-n
    k = np.arange(n + 1)
    al = (k * math.log(p) + (n - k) * math.log(1 - p)) / math.log(eps)
    bins = np.sort(np.r_[al - 1e-9, al + 1e-9])
    _, f, h = fr.falpha_histogram(mu, eps, bins)
    assert sorted(h[::2].tolist()) == sorted(math.comb(n, int(x)) for x in k)


def test_rank_mass_and_cumulative_consistent():
    i = np.arange(1, 4097)
    m = (1.0 / i) / np.sum(1.0 / i)                                                # Zipf masses
    rng = np.random.default_rng(0)
    boxes = rng.permutation(m)
    top = fr.rank_mass(boxes, 0, 4096)
    assert np.allclose(i * top, i[0] * top[0], rtol=1e-12)
    eps = 1e-3
    al = np.array([0.5, 0.8, 1.0])
    n, _ = fr.cumulative_spectrum(boxes, eps, al)
    assert n.tolist() == [int(np.sum(boxes >= eps**a)) for a in al]


def test_atom_plus_lebesgue_Ds():
    """One atom a in one cell plus uniform (1 - a): Z_s(eps) = (a + (1-a)eps)^s + (1/eps - 1)((1-a)eps)^s exactly."""
    a, n = 0.2, 12
    mu = np.full(2**n, (1 - a) / 2**n)
    mu[0] += a
    js = list(range(0, 8))
    for s in (0.5, 2.0, 3.0):
        Z = [(a + (1 - a) * 2.0 ** (j - n)) ** s + (2.0 ** (n - j) - 1) * ((1 - a) * 2.0 ** (j - n)) ** s for j in js]
        le = np.log([2.0 ** (j - n) for j in js])
        ref = np.diff(np.log(Z)) / np.diff(le) / (s - 1)
        assert np.allclose(fr.partition_dims(mu, js, [s])[s], ref, rtol=1e-10)


def test_remove_atoms():
    Xs, Ms = 1000, 2**14
    out = fr.remove_atoms(np.ones(Ms), Xs, radius=20.0)
    r = int(20 * Ms / Xs)
    assert out[: r + 1].sum() == 0 and out[Ms - r :].sum() == 0 and out[Ms // 2 - r : Ms // 2 + r + 1].sum() == 0
    assert out.sum() == Ms - (r + 1) - r - (2 * r + 1)


# ------------------------------------------------------------------ the Farey count and the kernel
def test_spike_count_constant(tables):
    mu, phi = tables
    for y in (1e3, 1e4):
        assert 0.99 <= fr.spike_count(mu, phi, y) / (y * y / 2) <= 1.01
    q = np.arange(1, 10**5 + 1)
    assert np.sum(np.where(mu[q] != 0, phi[q], 0)) / (1e10 / 2) == pytest.approx(0.4283, abs=0.002)


def test_farey_staircase(tables):
    mu, phi = tables
    st = fr.farey_staircase(10.0, mu, phi, qmax=100)
    assert st["q"][:4].tolist() == [3, 6, 5, 10] and st["phi"][:4].tolist() == [2, 2, 4, 4]
    assert st["rank_end"][1] == 4 and st["rank_end"][3] == 12
    assert st["mass"][1] == pytest.approx(1 / (4 * 10.0)) and st["mass"][3] == pytest.approx(1 / (16 * 10.0))


def test_sinc_integrals_and_level_measure():
    assert 2 * float(fr.sinc_F(np.array([0.5]))[0]) == pytest.approx(0.7737, abs=1e-4)
    assert float(fr.sinc_F(np.array([1e6]))[0]) == pytest.approx(0.5, abs=1e-6)
    assert fr.sinc4_tail(0.0) == pytest.approx(2.0 / 3.0, abs=1e-4)
    from scipy.optimize import brentq
    lv = fr.SincLevels()
    for u in (0.3, 0.5, 0.9):
        xu = brentq(lambda x: np.sinc(x) - u, 1e-9, 1.0)
        assert lv.measure(u)[()] == pytest.approx(2 * xu, abs=2 * lv.dx)          # resolution of the table
    assert fr.box_capture(1e-3) == pytest.approx(1e-3, rel=1e-3)                   # a box of width w << 1/N catches ~ w
    assert 0 < fr.box_capture(1.0) < 0.7737 < fr.box_capture(50.0) < 1.0


# ------------------------------------------------------------------ M1: collapse
def test_all_a_fft_matches_S_rational(small, tables):
    mu, phi = tables
    lam, psi = small["lam"], small["psi"]
    n = np.flatnonzero(lam)
    for q in (7, 30, 101, 210):
        res = np.bincount(n % q, weights=lam[n], minlength=q)
        fft = q * np.fft.ifft(res)
        assert max(abs(fft[a] - S_rational(a, q, res)) for a in range(q)) <= 1e-9 * psi


def test_second_moment_quadratic_forms(small, tables):
    mu, phi = tables
    lam = small["lam"]
    n = np.flatnonzero(lam)
    from singular_series import ramanujan_divisor
    for q in (31, 30, 101):
        res = np.bincount(n % q, weights=lam[n], minlength=q)
        r = np.arange(q)
        red = np.gcd(r, q) == 1
        delta = np.where(red, res - res[red].sum() / phi[q], 0.0)
        z = np.array([np.sum(np.exp(2j * np.pi * a * r / q) * delta) for a in r[red]])
        cq = np.array([ramanujan_divisor(q, int(d)) for d in range(q)], dtype=float)
        assert np.mean(np.abs(z) ** 2) == pytest.approx(delta @ cq[(r[:, None] - r[None, :]) % q] @ delta / phi[q], rel=1e-10)
        assert np.mean(z**2).real == pytest.approx(delta @ cq[(r[:, None] + r[None, :]) % q] @ delta / phi[q], rel=1e-10)
        if q in (31, 101):
            V = float(np.sum(delta**2))
            assert np.mean(np.abs(z) ** 2) == pytest.approx(V * q / (q - 1), rel=1e-10)


def test_beckmann_rice_cos():
    for s2 in (0.1, 1.0, 5.0):
        b = fr.beckmann_moments(s2 / 2, s2 / 2)
        r = fr.rice_moments(s2)
        # 80-point Gauss-Hermite vs the exact Rice law: |1 + Z| is not smooth at Z = -1, so ~1e-4 (1e-3 at sigma^2 = 5)
        assert b["mean"] == pytest.approx(r["mean"], abs=1e-3) and b["sd"] == pytest.approx(r["sd"], abs=1e-3)
        assert b["cos"] == pytest.approx(r["cos"], abs=3e-3)
        rng = np.random.default_rng(1)
        Z = math.sqrt(s2 / 2) * (rng.standard_normal(200000) + 1j * rng.standard_normal(200000))
        y = 1 + Z
        assert r["cos"] == pytest.approx(float(np.mean(y.real / np.abs(y))), abs=0.005)
    b = fr.beckmann_moments(1e-8, 1e-8)
    assert b["mean"] == pytest.approx(1.0, abs=1e-6) and b["sd"] == pytest.approx(1e-4, rel=0.05)
    r = fr.rice_moments(6e-5)                                                         # tiny spread: no overflow
    assert r["mean"] == pytest.approx(1 + 6e-5 / 4, abs=1e-6) and r["cos"] == pytest.approx(1 - 6e-5 / 4, abs=1e-6)


def test_fg_and_qstar():
    assert fr.fg_factor(101, [101]) == pytest.approx(2.154, abs=1e-3)
    assert fr.q_star_fg(10**6, 999586.6) == pytest.approx(513, abs=2)


def test_collapse_small_grid(small, tables):
    mu, phi = tables
    lam, psi = small["lam"], small["psi"]
    n = np.flatnonzero(lam)
    rows = fr.spike_collapse(n, lam[n], psi, [q for q in range(3, 31) if mu[q] != 0], mu, phi)
    assert all(abs(r["mean"] - 1) < 0.02 for r in rows)
    assert all(r["cos"] > 0.9 for r in rows)                                         # q ~ 30 is near the X = 10^4 crossover


# ------------------------------------------------------------------ M3: tails
def test_spike_ccdf_asymptote(tables):
    """With Q = 0 and many spikes, P_model 2 N |m| t mean / psi^2 -> 1 where the relevant phi are 10^3 - 10^4."""
    mu, phi = tables
    N, psi, mean_m, meas = 10**6, 1e6, 1.0, 1.0
    ts = np.array([(psi / 5e3) ** 2, (psi / 2e3) ** 2])
    P = fr.spike_ccdf_model(ts, mean_m, meas, N, psi, mu, phi, 0, 0.0, qcap=200000)
    assert np.allclose(P * 2 * N * meas * ts * mean_m / psi**2, 1.0, atol=0.03)


def test_hill_and_concentration():
    rng = np.random.default_rng(2)
    x = 1.0 / rng.random(10**6)                                                     # Pareto, index 1
    assert 0.93 <= fr.hill(x, 1000) <= 1.07
    e = rng.exponential(size=10**6)
    assert fr.concentration(e, [0.01])[0] == pytest.approx(0.0561, abs=0.004)


def test_sup_prediction_small_grid(small, tables):
    mu, phi = tables
    S = small["S"]
    A = MajorArcs(X, M, 20, 10, mu, phi)
    minor = ~A.mask()
    sup2 = float(np.max(np.abs(S[minor]) ** 2))
    sp = fr.sup_prediction(20, 10, small["psi"], mu, phi)
    assert abs(sp["sup2"] / sup2 - 1) < 0.25 and sp["active"]


# ------------------------------------------------------------------ M4: walks
def test_structure_function_orders():
    rng = np.random.default_rng(3)
    walk = np.cumsum(rng.standard_normal(2**18))                                   # Brownian
    lags = [2**j for j in range(2, 9)]                                              # >= 1000 independent increments
    for p in (0.5, 1.0, 2.0):
        H = fr.local_hurst(lags, fr.structure_function(walk, lags, p), p)
        assert np.all(np.abs(H - 0.5) < 0.03)
    cauchy = np.cumsum(rng.standard_cauchy(2**18))                                  # Levy flight, index 1: H_1/2 = 1
    lg = [4, 8, 16, 32, 64]
    assert np.all(np.abs(fr.local_hurst(lg, fr.structure_function(cauchy, lg, 0.5), 0.5) - 1.0) < 0.05)
    line = np.arange(2**12, dtype=float)                                           # smooth: H = 1
    assert np.allclose(fr.local_hurst(lags[:5], fr.structure_function(line, lags[:5])), 1.0)


def test_bridge_H():
    rng = np.random.default_rng(4)
    n = 2**12
    lags = [n // 16, n // 8]
    S = np.zeros(2)
    for _ in range(200):
        w = rng.standard_normal(n)
        F = np.cumsum(w - w.mean())                                                # pinned at both ends
        S += fr.structure_function(F, lags)
    H = fr.local_hurst(lags, S)[0]
    assert abs(H - fr.H_bridge(lags[0], n)) < 0.03


def test_H_triangular_values():
    ref = [0.990, 0.962, 0.873, 0.722, 0.624, 0.570, 0.540, 0.523, 0.513]
    assert np.allclose([fr.H_triangular(0.119 * 2**j) for j in range(9)], ref, atol=2e-3)


def test_minor_walk_spectral_identity(small, tables):
    mu, phi = tables
    S = small["S"]
    A = MajorArcs(X, M, 20, 10, mu, phi)
    k, w = fr.minor_walk(S, ~A.mask(), X)
    x = np.zeros(M, dtype=complex)
    x[k] = w
    F = np.cumsum(x)
    lags = [1, 8, 64, 512]
    direct = []
    for L in lags:
        d = np.empty(M, dtype=complex)
        d[: M - L] = F[L:] - F[: M - L]
        d[M - L :] = F[:L] + F[-1] - F[M - L :]
        direct.append(np.mean(np.abs(d) ** 2))
    assert np.allclose(fr.spectral_S2(np.fft.fft(x), lags, chunk=1 << 14), direct, rtol=1e-9)


def test_order_parameter_entropy():
    w = np.exp(2j * np.pi * np.arange(6400) / 6400)                                 # uniform phases: r = 0, H = log2(64)
    r, H = fr.order_parameter_entropy(w)
    assert r < 1e-12 and H == pytest.approx(6.0, abs=1e-3)
    r, H = fr.order_parameter_entropy(np.full(100, 2.0 + 0j))
    assert r == pytest.approx(1.0) and H == pytest.approx(0.0)


# ------------------------------------------------------------------ the surrogate
def test_hybrid_small(small, tables):
    mu, phi = tables
    lam, psi, L2 = small["lam"], small["psi"], small["L2"]
    W = fr.ramanujan_weights(X, psi, mu, phi, [30])[30]
    # exact identity: sum_n w*(n) = (psi/X) sum_{q<=30} mu(q)/phi(q) sum_{d|q} d mu(q/d) floor(X/d)
    ref = sum(mu[q] / phi[q] * sum(d * mu[q // d] * (X // d) for d in range(1, q + 1) if q % d == 0)
              for q in range(1, 31) if mu[q] != 0) * psi / X
    assert W.sum() == pytest.approx(ref, rel=1e-9) and abs(W.sum() / psi - 1) < 1e-3
    h = fr.hybrid_weights(W, L2, seed=37)
    Sh = S_grid(h, M)
    assert np.sum(np.abs(Sh) ** 2) / M == pytest.approx(float(np.sum(h**2)), rel=1e-10)
    assert abs(np.sum(h**2) / L2 - 1) < 0.02


# ------------------------------------------------------------------ uncertainty (Planck, Fourier form)
def test_uncertainty_products():
    b = fr.uncertainty_product("bump", 20)
    assert b["product"] >= b["bound"] and b["product"] == pytest.approx(0.0947, abs=5e-4)
    s = [fr.uncertainty_product("sharp", L)["sigma_xi"] for L in (20, 80, 320)]
    assert s[1] / s[0] == pytest.approx(2.0, rel=0.02) and s[2] / s[1] == pytest.approx(2.0, rel=0.02)   # ~ sqrt(L)


# ------------------------------------------------------------------ the Mandelbrot set
@pytest.mark.parametrize("p, q, centre", [(1, 2, -1.0 + 0j), (1, 3, -0.1225611668766536 + 0.7448617666197442j)])
def test_bulb_centres_known_values(p, q, centre):
    c, per, rot = fr.bulb_centre(p, q)
    assert per == q and rot and abs(c - centre) < 1e-12


@pytest.mark.parametrize("p, q", [(1, 40), (1, 100)])
def test_bulb_centres_large_q(p, q):
    c, per, rot = fr.bulb_centre(p, q)
    assert per == q and rot


def test_bulb_size_law():
    for p, q in [(1, 4), (2, 5), (1, 7), (3, 8), (1, 13)]:
        s = fr.leading_size(p, q)
        assert 0.97 <= fr.bulb_diam_half(p, q) / s <= 1.16
        assert 0.98 <= abs(fr.bulb_centre(p, q)[0] - fr.cardioid_root(p, q)) / s <= 1.07


def test_bulb_all_q40():
    n = 0
    for q in range(2, 41):
        for p in range(1, q // 2 + 1):
            if math.gcd(p, q) != 1:
                continue
            c, per, rot = fr.bulb_centre(p, q)
            assert per == q and rot
            assert abs(fr.bulb_continue(p, q, 0.999999 + 0j) - fr.cardioid_root(p, q)) <= 1e-2 * fr.leading_size(p, q)
            n += 1
    assert n == 245
    assert fr.bulb_diam_half(1, 2) == pytest.approx(0.25, abs=1e-10)
    assert fr.bulb_diam_half(1, 3) / fr.leading_size(1, 3) == pytest.approx(0.9806, abs=1e-3)
    assert fr.bulb_diam_half(1, 40) / fr.leading_size(1, 40) == pytest.approx(1.0264, abs=1e-3)


def test_mediant_counts_F8(tables):
    mu, phi = tables
    r = fr.mediant_test(8, 120, mu, phi)
    assert (r["pairs"], r["bulb_mediant"], r["spike_mediant"], r["mediant_nonsquarefree"]) == (22, 22, 10, 8)
    assert fr.farey(5) == [(0, 1), (1, 5), (1, 4), (1, 3), (2, 5), (1, 2), (3, 5), (2, 3), (3, 4), (4, 5), (1, 1)]


# ------------------------------------------------------------------ direct tests of the functions behind the report's numbers
def test_spike_collapse_rows_direct(small, tables):
    """Every moment in a collapse row equals a brute-force computation from S(a/q) evaluated by S_rational."""
    mu, phi = tables
    lam, psi = small["lam"], small["psi"]
    n = np.flatnonzero(lam)
    for row in fr.spike_collapse(n, lam[n], psi, [30, 31, 101], mu, phi):
        q = row["q"]
        res = np.bincount(n % q, weights=lam[n], minlength=q)
        Sa = np.array([S_rational(a, q, res) for a in range(q) if math.gcd(a, q) == 1])
        y = np.abs(Sa) * phi[q] / psi
        Z = mu[q] * Sa * phi[q] / psi - 1.0
        red = np.gcd(np.arange(q), q) == 1
        assert row["mean"] == pytest.approx(y.mean(), rel=1e-12)
        assert row["sd"] == pytest.approx(y.std(), rel=1e-9)
        assert row["cos"] == pytest.approx(np.mean((mu[q] * Sa).real / np.abs(Sa)), rel=1e-12)
        assert row["E2"] == pytest.approx(np.mean(np.abs(Z) ** 2), rel=1e-9)
        assert row["EZ2"] == pytest.approx(np.mean(Z**2).real, rel=1e-9, abs=1e-15)
        assert row["V"] == pytest.approx(np.sum((res[red] - res[red].sum() / phi[q]) ** 2), rel=1e-12)


def test_spike_ccdf_truncation_branch(tables):
    """With Q > 0 the spikes q <= Q contribute only beyond |x| = c; the rest contribute in full."""
    mu, phi = tables
    lv = fr.SincLevels()
    N, psi, mean_m, meas, Q, c = 10**6, 1e6, 1.0, 0.5, 5, 3.0
    ts = np.array([1e6, 1e8, 1e10])
    P = fr.spike_ccdf_model(ts, mean_m, meas, N, psi, mu, phi, Q, c, qcap=2000, levels=lv)
    P0 = fr.spike_ccdf_model(ts, mean_m, meas, N, psi, mu, phi, 0, 0.0, qcap=2000, levels=lv)
    for t, p in zip(ts, P):
        tot = 0.0
        for q in range(1, 2001):
            if mu[q] != 0:
                u = math.sqrt(t * mean_m) * phi[q] / psi
                tot += phi[q] * float(lv.measure(u, c if q <= Q else 0.0))
        assert p == pytest.approx(tot / (N * meas), rel=1e-12)
    assert np.all(P <= P0) and P[-1] < 0.9 * P0[-1]


def test_structure_function_drift_identity():
    rng = np.random.default_rng(3)
    F = np.cumsum(rng.normal(size=5000) + 1j * rng.normal(size=5000) + 0.3)
    lags = [1, 4, 16]
    for L, v in zip(lags, fr.structure_function(F, lags, 2.0, drift=True)):
        d = F[L:] - F[:-L]
        assert v == pytest.approx(np.mean(np.abs(d - d.mean()) ** 2), rel=1e-10)


def test_octave_energy_exact_windows(tables):
    """The union of the windows |X(k/M - a/q)| <= 1/2, built independently from integer interval ends, equals the
    library's, point for point; from j = 7 on the windows overlap and are counted once."""
    mu, phi = tables
    a2 = np.random.default_rng(5).random(M)
    oc = fr.octave_energy(a2, X, mu, phi, jmax=8)
    for o in oc:
        mask = np.zeros(M, dtype=bool)
        qs = [q for q in range(max(1 << o["j"], 1), 1 << (o["j"] + 1)) if mu[q] != 0]
        for q in qs:
            for a in range(q) if q > 1 else [0]:
                if math.gcd(a, q) != 1:
                    continue
                lo, hi, den = 2 * X * a * M - M * q, 2 * X * a * M + M * q, 2 * X * q    # |2X(kq - aM)| <= Mq
                mask[np.arange(-((-lo) // den), hi // den + 1) % M] = True
        assert round(o["measure"] * M) == int(mask.sum())
        assert o["energy"] == pytest.approx(a2[mask].sum() / M, rel=1e-12)
        nwin = sum(int(phi[q]) for q in qs)
        if o["j"] <= 5:
            assert abs(o["measure"] - nwin / X) <= nwin / M        # disjoint windows of M/X +- 1 grid points
        if o["j"] >= 7:
            assert o["measure"] < nwin / X - 0.1                     # overlapping windows are counted once


def test_fourth_moment_model_direct(tables):
    from scipy.integrate import quad

    mu, phi = tables
    Q, c, psi, N = 20, 10.0, 1e6, 10**6
    got = fr.fourth_moment_model(Q, c, psi, N, mu, phi, qcap=500)
    hi = sum(1.0 / phi[q] ** 3 for q in range(Q + 1, 501) if mu[q] != 0)
    lo = sum(1.0 / phi[q] ** 3 for q in range(1, Q + 1) if mu[q] != 0)
    T4 = 2 * sum(quad(lambda x: np.sinc(x) ** 4, x0, x0 + 1)[0] for x0 in range(10, 4000)) + 2 / (8 * math.pi**4 * 4000**3)
    assert got == pytest.approx(psi**4 / N * (2 / 3 * hi + T4 * lo), rel=1e-6)


def test_cells_large_values_and_window_maxima():
    v = np.arange(32.0)
    minor = np.ones(32, dtype=bool)
    minor[9] = False
    assert list(fr.cell_maxima(v, minor)) == [7.0, 23.0, 31.0]
    assert fr.large_value_counts(np.arange(16.0), 4, [0.5, 1.0, 1.5], cell=4) == [4, 3, 2]
    assert list(fr.window_max_median(np.arange(16.0), [1, 2], 2.0)) == [4.0, 4.5]


def test_q_star_measured_small_grid(small):
    """K(q) equals the direct S_rational computation, and q* at X = 10^4 lands near the FG prediction (75)."""
    lam, psi = small["lam"], small["psi"]
    n = np.flatnonzero(lam)
    primes = [int(p) for p in primes_upto(400) if p >= 7]
    qs, K = fr.q_star_measured(n, lam[n], psi, X, primes)
    q = 31
    res = np.bincount(n % q, weights=lam[n], minlength=q)
    Sa = np.array([S_rational(a, q, res) for a in range(1, q)])
    assert K[q] == pytest.approx(math.sqrt(np.mean(np.abs(Sa + res[1:].sum() / (q - 1)) ** 2)) * (q - 1) / psi, rel=1e-10)
    assert qs is not None and 50 <= qs <= 110


def test_mediant_detail_consistent(tables):
    mu, phi = tables
    d = fr.mediant_test(8, 120, mu, phi)
    assert len(d["detail"]) == d["pairs"]
    assert sum(x["bulb_ok"] for x in d["detail"]) == d["bulb_mediant"]
    assert sum(x["spike_ok"] for x in d["detail"]) == d["spike_mediant"]
    assert sum(x["nonsquarefree"] for x in d["detail"]) == d["mediant_nonsquarefree"]
