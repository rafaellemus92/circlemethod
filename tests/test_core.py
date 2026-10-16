"""pytest suite: anchors, quadrature certificates, Ramanujan sums, singular series, arc geometry."""
from __future__ import annotations

import math

import numpy as np
import pytest

import arcs as ar
import expsum as ex
import singular_series as ss
import sieve as sv

N = 10**6
M = 2**23


@pytest.fixture(scope="module")
def base():
    lam = sv.von_mangoldt(N)
    S = ex.S_grid(lam, M)
    return {"lam": lam, "S": S, "R2": sv.R2_direct(lam, N), "L2": sv.sum_lambda_sq(lam)}


# ------------------------------------------------------------------ module 0
def test_pi_1e6():
    assert int(sv.prime_sieve(N).sum()) == 78498


def test_goldbach_G_1e6():
    assert sv.goldbach_G(N, sv.prime_sieve(N)) == 5402


def test_von_mangoldt_small():
    lam = sv.von_mangoldt(30)
    expect = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3, 11: 11, 13: 13, 16: 2, 17: 17, 19: 19, 23: 23,
              25: 5, 27: 3, 29: 29}
    for n in range(31):
        assert lam[n] == pytest.approx(math.log(expect[n]) if n in expect else 0.0)


def test_comet_matches_direct_count():
    Ns, G = sv.goldbach_comet(20000)
    is_p = sv.prime_sieve(20000)
    for n in (4, 6, 100, 1000, 9998, 20000):
        assert G[(n - 4) // 2] == sv.goldbach_G(n, is_p)


# ------------------------------------------------------------------ module 1
def test_parseval_on_grid(base):
    assert ex.grid_mass(base["S"]) == pytest.approx(base["L2"], rel=1e-10)


def test_grid_integral_equals_R2(base):
    assert ex.grid_power(base["S"], 2, N).real == pytest.approx(base["R2"], rel=1e-10)


def test_fft_sign_convention(base):
    S, lam = base["S"], base["lam"]
    chk = ex.sign_convention_checks(S, lam, N, base["R2"])
    assert chk["S0_relerr"] < 1e-12
    assert chk["S_half_abserr"] < 1e-6 * chk["psi"]
    assert chk["symmetry_maxrel"] < 1e-12
    assert chk["direct_vs_grid_maxrel"] < 1e-10
    assert np.sign(chk["Im_S_1_over_M_grid"]) == np.sign(chk["Im_S_1_over_M_direct"])
    # the discriminating check: the wrong sign kills the binary integral
    assert abs(chk["grid_R2_wrong_sign"]) < 1e-6 * base["R2"]


def test_ternary_grid_integral_equals_convolution():
    n0 = 10**4
    lam = sv.von_mangoldt(n0)
    S = ex.S_grid(lam, 2**17)            # M > 2N' needed for an alias-free cube
    for Np in (n0 - 1, n0 + 1):
        assert ex.grid_power(S, 3, Np).real == pytest.approx(sv.R3_exact(lam, Np), rel=1e-10)


def test_S_rational_matches_direct(base):
    lam = base["lam"]
    n, w = ex.nonzero_terms(lam)
    for a, q in ((1, 3), (2, 5), (5, 12), (7, 30)):
        exact = ex.S_rational(a, q, sv.psi_residues(lam, q))
        assert abs(exact - ex.S_direct([a / q], n, w)[0]) < 1e-7 * sv.chebyshev_psi(lam)


def test_dirichlet_kernel():
    n0 = 1000
    beta = np.array([0.0, 1e-4, 0.0123, -0.3, 0.5])
    direct = np.array([np.sum(np.exp(2j * np.pi * np.arange(1, n0 + 1) * b)) for b in beta])
    assert np.allclose(ex.dirichlet(beta, n0), direct, atol=1e-8)


# ------------------------------------------------------------------ Ramanujan sums
def test_ramanujan_divisor_formula_vs_definition():
    for q in range(1, 61):
        for n in (0, 1, 2, 6, 12, 30, 35, 210, 999999, 10**6):
            assert ss.ramanujan_divisor(q, n) == pytest.approx(ss.ramanujan_direct(q, n), abs=1e-9)


def test_ramanujan_holder_vs_divisor():
    mu, phi = sv.mobius_phi(500)
    qs = np.arange(1, 501)
    for n in (1, 10**6, 999999, 1000001, 30030, 2**10):
        h = ss.ramanujan_holder(qs, n, mu, phi)
        d = np.array([ss.ramanujan_divisor(int(q), n) for q in qs])
        assert np.array_equal(h, d)


def test_ramanujan_multiplicative_in_q():
    for q in range(1, 40):
        for r in range(1, 40):
            if math.gcd(q, r) != 1:
                continue
            for n in (1, 12, 30, 10**6):
                assert ss.ramanujan_divisor(q * r, n) == ss.ramanujan_divisor(q, n) * ss.ramanujan_divisor(r, n)


# ------------------------------------------------------------------ singular series
@pytest.fixture(scope="module")
def mu_phi_1e6():
    return sv.mobius_phi(10**6)


@pytest.mark.parametrize("n", [10**6, 30030, 2 * 3**4 * 7, 2**20, 999998, 1531530])
def test_SS2_sum_vs_product(n, mu_phi_1e6):
    # Very smooth n (30030, 1531530) converge slowest: the error at Q = 2e5 is still ~1e-4, at Q = 1e6 ~3e-5.
    mu, phi = mu_phi_1e6
    assert ss.SS2_partial(n, 10**6, mu, phi)[-1] == pytest.approx(ss.SS2_product(n), abs=1e-4)


@pytest.mark.parametrize("n", [999999, 1000001, 3 * 5 * 7 * 11 * 13, 10**4 + 1])
def test_SS3_sum_vs_product(n):
    assert ss.SS3_partial(n, 10**4)[-1] == pytest.approx(ss.SS3_product(n), abs=1e-5)


def test_twin_prime_constant():
    assert ss.C2_product(10**6) == pytest.approx(ss.C2_KNOWN, abs=1e-7)


def test_hl_integral_vs_quad():
    from scipy.integrate import quad
    for n in (100, 10**4, 10**6):
        ref = quad(lambda t: 1 / (math.log(t) * math.log(n - t)), 2, n - 2, limit=400)[0]
        assert ss.hl_integral([n])[0] == pytest.approx(ref, rel=1e-10)


# ------------------------------------------------------------------ arcs
def test_disjointness_assertion():
    mu, phi = sv.mobius_phi(200)
    with pytest.raises(AssertionError):
        ar.MajorArcs(N, M, 159, 20, mu, phi)          # 2*20*159^2 = 1,011,240 >= N: rejected up front
    A = ar.MajorArcs(N, M, 150, 20, mu, phi)          # 900,000 < N: accepted
    assert A.n_arcs == int(phi[1:151].sum())
    # every arc contributes the same number of grid points (+-1) and they add up: no point was dropped
    assert abs(len(A.idx) - A.n_arcs * 2 * 20 * M / N) <= A.n_arcs


def test_on_grid_overlap_detector(monkeypatch):
    """Bypass the 2cQ^2 < N pre-check: genuinely overlapping arcs must be caught on the grid."""
    mu, phi = sv.mobius_phi(50)
    monkeypatch.setattr(ar, "admissible", lambda *a: True)
    with pytest.raises(AssertionError, match="overlap on the grid"):
        ar.MajorArcs(1000, 2**14, 30, 10, mu, phi)    # 2*10*30^2 = 18,000 >> 1000: arcs overlap


def test_arc_membership_exact():
    mu, phi = sv.mobius_phi(60)
    A = ar.MajorArcs(N, M, 50, 10, mu, phi)
    k = A.idx
    d = np.abs(k / M - A.a / A.q)
    d = np.minimum(d, 1 - d)
    assert np.all(d <= 10 / N + 1e-15)
    # membership on the grid agrees with the analytic test A.contains
    sample = np.unique(np.concatenate([A.idx[:: 97], (A.idx[:: 97] + 50) % M]))
    assert np.array_equal(np.isin(sample, A.idx), A.contains(sample / M))


def test_major_plus_minor_is_whole(base):
    mu, phi = sv.mobius_phi(30)
    A = ar.MajorArcs(N, M, 20, 10, mu, phi)
    S = base["S"]
    d = ar.arc_integrals(S, A, N)
    Pm, Im = ar.minor_integrals(S, A.mask(), N)
    assert d["I_M"].real + Im.real == pytest.approx(base["R2"], rel=1e-12)
    assert d["P_M"] + Pm == pytest.approx(base["L2"], rel=1e-12)
