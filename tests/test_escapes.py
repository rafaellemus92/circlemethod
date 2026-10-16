"""pytest suite for escapes.py on a small grid: X = 10^4, M = 2^17, (Q, c) = (20, 10), L = 13.

2cQ^2 = 8000 < X, so the arcs are admissible; M > 3X, so binary and ternary counts are exact on the grid;
2X + K 2^13 < X + M for K <= 14, so R_K is exact on the grid for K <= 14.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

import escapes as es
from expsum import S_grid, e_minus
from sieve import autoconvolve
from singular_series import SS2_even_array

X = 10**4
M = 2**17
QC = (20, 10)
L = 13


@pytest.fixture(scope="module")
def small(tmp_path_factory):
    ctx = es.Context(X, M, tmp_path_factory.mktemp("cache"), [QC])
    S2 = ctx.S * ctx.S
    a2 = ctx.S.real ** 2 + ctx.S.imag ** 2
    A = ctx.arc(*QC)
    return {"ctx": ctx, "S2": S2, "a2": a2, "A": A, "IM": es.major_part_all_N(S2, A.idx)}


@pytest.fixture(scope="module")
def linnik(small):
    return es.e2_linnik(small["ctx"], small["S2"], small["a2"], [QC], Kmax=14, Kfull=14, sparse_seeds=2)


def test_bessel_identity_in_N(small):
    ctx, S2, a2, A, IM = small["ctx"], small["S2"], small["a2"], small["A"], small["IM"]
    r2full = np.zeros(M)
    r2full[: len(ctx.r2)] = ctx.r2
    Im = r2full - IM
    minor = ~A.mask()
    E4m = float(np.sum(a2[minor] ** 2)) / M
    assert abs(float(np.sum(Im ** 2)) - E4m) <= 1e-12 * E4m
    # I_M(N) from the FFT is the demodulated arc integral, and I_M + I_m(direct) = R2(N) for N <= X
    k = np.flatnonzero(minor)
    for N in (2, 1000, 5002, X):
        direct_M = np.sum(S2[A.idx] * e_minus(N, A.idx, M)) / M
        direct_m = np.sum(S2[k] * e_minus(N, k, M)) / M
        assert abs(direct_M.real - IM[N]) <= 1e-9 * ctx.R2
        assert abs(direct_M.imag) <= 1e-9 * ctx.R2
        assert abs(IM[N] + direct_m.real - ctx.r2[N]) <= 1e-9 * ctx.R2


def test_fourth_moment_is_additive_energy(small):
    ctx, a2 = small["ctx"], small["a2"]
    full4 = float(np.sum(a2 ** 2)) / M
    energy = math.fsum(ctx.r2 ** 2)
    assert abs(full4 - energy) <= 1e-12 * energy


def test_transmission_certificate(small):
    ctx, S2, a2, A, IM = small["ctx"], small["S2"], small["a2"], small["A"], small["IM"]
    Pm = (float(a2.sum()) - float(a2[A.idx].sum())) / M
    r2full = np.zeros(M)
    r2full[: len(ctx.r2)] = ctx.r2
    assert np.all(np.abs(r2full - IM) <= Pm * (1 + 1e-12))            # triangle inequality, every N mod M
    Ns = es.even_targets(X)
    trans = Ns[IM[Ns] > Pm]
    assert np.all(ctx.r2[trans] > 0)
    o = es.e1_exceptional(ctx, S2, a2, *QC, SS2_even_array(X))
    assert o["certificates_ok"] and o["bessel_relerr"] < 1e-12
    assert o["bound_oracle"] <= o["bound_4th"] <= o["bound_natural_sup"] <= o["bound_textbook"]
    assert o["F"] <= o["bound_oracle"]


def test_G_constructions():
    G = S_grid(es.pow2_indicator(X), M)
    G_ref = es.powers_of_two_grid(M, L)
    assert np.max(np.abs(G - G_ref)) <= 1e-12 * L
    assert abs(float(np.mean(np.abs(G) ** 2)) - L) <= 1e-12 * L
    assert int(es.pow2_indicator(X).sum()) == L


def test_kfold_power_counts():
    for K in range(4):
        c = es.kfold_power_counts(L, K, 2**15)
        assert c.sum() == L ** K                                          # every tuple lands below 2^15 for K <= 3
    c2 = es.kfold_power_counts(L, 2, 2**15)
    assert c2[6] == 2 and c2[8] == 1                                      # 6 = 2 + 4 = 4 + 2; 8 = 4 + 4


def test_RK_grid_equals_direct(linnik):
    rk = linnik["RK_direct"]
    rel = np.abs(linnik["RK_grid"].real - rk) / rk
    assert np.all(rel[:7] <= 1e-12)
    assert np.all(rel <= 1e-11)                                           # K <= 14: still exact on the grid
    assert np.all(np.abs(linnik["RK_grid"].imag) <= 1e-12 * rk)
    assert linnik["rho"][0] == pytest.approx(1.0)


def test_theorem_shaped_bound_dominates(linnik):
    r = linnik["rows"][QC]
    for K in range(2, 7):
        assert r["theorem"][1]["split"][K][0] >= r["joint_ratio"][K] * (1 - 1e-12)
        assert r["theorem"][1]["split_exact"][K][0] >= r["joint_ratio"][K] * (1 - 1e-12)
        assert r["sup_ratio"][K] >= r["joint_ratio"][K] * (1 - 1e-12)
        assert r["holder_ratio"][K] >= r["joint_ratio"][K] * (1 - 1e-12)


def test_almost_prime_counts(small):
    ctx, A = small["ctx"], small["A"]
    W = ctx.S * e_minus(X, np.arange(M, dtype=np.int64), M)
    om = es.omega_array(X)
    n = np.arange(X + 1)
    sets = [("Omega<=2", ((n >= 2) & (om <= 2)).astype(float)),
            ("odd Omega<=2", ((n >= 2) & (om <= 2) & (n % 2 == 1)).astype(float)),
            ("7-rough", es.rough_indicator(X, 7))]
    out = es.e3_almost_primes(ctx, np.abs(ctx.S), W, sets, [QC], [QC])
    for s in out:
        assert abs(s["grid"] - s["direct"]) <= 1e-12 * s["direct"]
        r = s["rows"][QC]
        assert r["CS"] >= r["joint"] * (1 - 1e-12)                        # Cauchy-Schwarz dominates the joint bound
    assert int(om[2 * 3 * 5]) == 3 and int(om[64]) == 6 and int(om[97]) == 1


@pytest.mark.parametrize("z, P", [(3, 2), (5, 6), (7, 30)])
def test_rough_closed_form(z, P):
    """T_z(a) = sum_{d | P(z)} mu(d) D_{floor(X/d)}(d a) - e(a), D_K(b) = sum_{j=1}^K e(j b)."""
    ind = es.rough_indicator(X, z)
    T = S_grid(ind, M)
    k = np.arange(M, dtype=np.int64)
    ref = -np.exp(2j * np.pi * k / M)
    for d in [d for d in range(1, P + 1) if P % d == 0]:
        mu = (-1) ** len([p for p in (2, 3, 5) if d % p == 0])
        K = X // d
        r = (d * k) % M
        ph = np.exp(2j * np.pi * r / M)
        num = 1 - np.exp(2j * np.pi * ((K * r) % M) / M)
        with np.errstate(invalid="ignore", divide="ignore"):
            D = np.where(r == 0, K, ph * num / (1 - ph))
        ref += mu * D
    assert np.max(np.abs(T - ref)) <= 1e-9 * ind.sum()


def test_smooth_counts(small):
    ctx = small["ctx"]
    w = ctx.lam * es.smooth_bump(np.arange(X + 1) / X)
    Se = S_grid(w, M)
    r2e = autoconvolve(w)
    grid2 = es.chunked_demod_sum(lambda k: Se[k] ** 2, X, M)
    assert abs(grid2.real - r2e[X]) <= 1e-12 * r2e[X]
    Np = 3 * X // 2 - 1
    m = np.arange(Np - X, Np + 1)                                         # w[Np - m] needs 0 <= Np - m <= X
    R3 = math.fsum(r2e[m] * w[Np - m])
    grid3 = es.chunked_demod_sum(lambda k: Se[k] ** 3, Np, M)
    assert abs(grid3.real - R3) <= 1e-12 * R3
    assert abs(grid3.imag) <= 1e-12 * R3


def test_eta_constants():
    c = es.eta_constants()
    assert c["int_eta"] == pytest.approx(0.6034502, abs=1e-6)
    assert c["int_eta2"] == pytest.approx(0.4916904, abs=1e-6)
    assert c["eta2_at_1"] == pytest.approx(c["int_eta2"], abs=1e-9)
    assert c["eta3_at_1"] == pytest.approx(0.0949085, abs=1e-6)
    assert c["eta3_at_1.5"] == pytest.approx(0.2446696, abs=1e-6)
    assert es.sharp_volume(1.0) == 0.5 and es.sharp_volume(1.5) == 0.75


def test_window_bound_monotone(small):
    ctx, A, IM = small["ctx"], small["A"], small["IM"]
    o = es.e0_obstacle(ctx, *QC, float(IM[X]), Qps=(20, 50))
    W = [d["W"] for d in o["dial"]]
    assert all(W[j + 1] <= W[j] * (1 + 1e-12) for j in range(len(W) - 1))
    assert W[0] == pytest.approx(o["P_m"] / o["I_M"], rel=1e-12)
    assert o["dial_floor"] <= W[-1] * (1 + 1e-12)
    assert o["dial_floor"] >= abs(o["I_m"]) / o["I_M"] * (1 - 1e-12)
    assert o["holder"] >= o["P_m"] * (1 - 1e-12)
    assert abs(o["I_m"].real - (ctx.R2 - o["I_M"])) <= 1e-9 * ctx.R2
    for m in o["model"]:                                                   # valid bound on the harmful direction
        assert m["bound"] >= -o["I_m"].real


def test_model_spikes_values(small):
    """S* = (mu(q)/phi(q)) v(alpha - a/q) inside the windows |alpha - a/q| <= 1/(2 q Q'), squarefree q <= Q'; 0 elsewhere."""
    from expsum import dirichlet

    ctx = small["ctx"]
    Qp = 30
    Sst = es.model_spikes(X, M, Qp, ctx.mu, ctx.phi, np.ones(M, dtype=bool))
    rng = np.random.default_rng(5)
    for k in rng.integers(0, M, 3000):
        hit = None
        for q in range(1, Qp + 1):
            if ctx.mu[q] == 0:
                continue
            a = int(round(k * q / M)) % q
            num = k * q - round(k * q / M) * M                     # beta = num / (q M), nearest a/q
            if math.gcd(a, q) == 1 and abs(num) * 2 * Qp <= M:
                hit = (q, num)
                break
        if hit is None:
            assert Sst[k] == 0
        else:
            q, num = hit
            ref = (ctx.mu[q] / ctx.phi[q]) * dirichlet(np.array([num / (q * float(M))]), X)[0]
            assert abs(Sst[k] - ref) <= 1e-9 * abs(ref) + 1e-12


def test_e5_scaling_counts(small):
    ctx = small["ctx"]
    c = 10
    SS2e = SS2_even_array(X)
    o = es.e5_scaling(X, ctx.S.copy(), ctx.lam, ctx.mu, ctx.phi, SS2e, c)
    Q = math.isqrt((X - 1) // (2 * c))
    A = es.MajorArcs(X, M, Q, c, ctx.mu, ctx.phi)
    a2 = np.abs(ctx.S) ** 2
    Pm = (a2.sum() - a2[A.idx].sum()) / M
    S2M = ctx.S[A.idx] ** 2
    Ns = es.even_targets(X)
    IM = np.array([np.sum(S2M * e_minus(int(N), A.idx, M)).real / M for N in Ns])   # direct, no FFT
    assert o["Q"] == Q and abs(o["P_m"] - Pm) <= 1e-12 * Pm
    assert o["n_transmitted"] == int(np.sum(IM > Pm)) and o["n_closes"] == int(np.sum(Pm / IM < 0.99))
    assert o["need"] == pytest.approx(Pm * SS2e[X] / IM[-1], rel=1e-9)
    assert o["certificates_ok"]


def test_report_decision_margin():
    import escapes_report as rep

    assert [rep.lab(x) for x in (0.989, 0.99, 1.0, 1.01, 1.0101)] == ["closes", "undecided", "undecided", "undecided", "fails"]
    assert rep.first_close([1.2, 0.995, 0.98, 0.5]) == 2
    assert rep.first_close([0.5, 0.98], start=1) == 1 and rep.first_close([1.0, 1.1]) is None
    assert rep.r3(0.995).endswith("(undecided)") and rep.r3(0.0072) == "0.0072"
