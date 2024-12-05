import numpy as np
import ngm
from numpy.testing import assert_array_equal
from numpy.testing import assert_allclose
import pytest


def test_dominant_eigen_simple():
    X = np.array([[1, 2], [2, 1]])
    e = ngm.dominant_eigen(X)
    assert np.isclose(e.value, 3.0)
    assert np.isclose(e.vector, np.array([0.5, 0.5])).all()


def test_dominant_eigen_bigger():
    X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    e = ngm.dominant_eigen(X)
    assert np.isclose(e.value, 16.116843969807043)
    assert np.isclose(e.vector, np.array([0.14719267, 1.0 / 3, 0.51947399])).all()


def test_R_vax():
    R_novax = np.array([[10.0, 0.1], [0.1, 1.0]])
    n = np.array([2.0, 8.0])
    n_vax = np.array([1.0, 0.0])
    ve = 1.0

    current = ngm.reduce_R(R=R_novax, p_vax=n_vax / n, ve=ve)
    expected = np.array(
        [
            [10.0 * 0.5, 0.1 * 0.5],
            [0.1, 1.0],
        ]
    )

    assert_array_equal(current, expected)


def test_simulate():
    # Tests ngm against itself
    n = np.array([200, 200, 100, 500])
    n_total = n.sum()
    beta = np.array(
        [
            [3.0, 0.5, 3.0, 0.5],
            [0.5, 0.5, 0.5, 0.5],
            [3.0, 0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5, 0.5],
        ]
    )
    R_novax = (beta.T * (n / n_total)).T

    n_vax = np.array([0, 0, 0, 0])
    p_severe = np.array([0.02, 0.06, 0.02, 0.02])
    ve = 1.0
    current = ngm.simulate(R_novax=R_novax, n=n, n_vax=n_vax, p_severe=p_severe, ve=ve)

    assert set(current.keys()) == {"Re", "R", "infections", "severe_infection_ratio"}
    assert np.isclose(current["Re"], 0.9213240982914677)
    assert_allclose(
        current["R"],
        np.array(
            [
                [0.6, 0.1, 0.6, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.3, 0.05, 0.05, 0.05],
                [0.25, 0.25, 0.25, 0.25],
            ]
        ),
    )
    assert_allclose(
        current["infections"],
        np.array([0.44507246, 0.10853944, 0.17503951, 0.2713486]),
    )

    # this test has values borrowed from when we had defined this output differently,
    # which is why the calculation looks so wacky
    assert_allclose(
        current["severe_infection_ratio"],
        sum(np.array([0.00820112, 0.006, 0.00322536, 0.005])) / current["Re"],
        rtol=1e-5,
    )


def test_ensure_positive():
    assert_array_equal(
        np.array([1, 2, 3]), ngm._ensure_positive_array(np.array([1, 2, 3]))
    )

    assert_array_equal(
        ngm._ensure_positive_array(np.array([-1, -2, -3])), np.array([1, 2, 3])
    )

    with pytest.raises(RuntimeError, match="all positive"):
        ngm._ensure_positive_array(np.array([1, -1]))


def test_kr():
    beta = np.array([[10, 0.1], [0.1, 1]])
    n = np.array([0.2, 0.8])
    R = (beta.T * n).T

    # see Keeling & Rohani, page 61
    R_p61 = np.array([[2, 0.02], [0.08, 0.8]])
    assert np.isclose(R, R_p61).all()

    r0 = ngm.dominant_eigen(R).value

    assert np.isclose(r0, 2.0013, atol=5e-5)


def test_eigenvectors():
    r = np.array([[3.1, 0.15, 1.7], [0.78, 1.5, 0.1], [0.32, 0.98, 1.1]])

    brute_force = np.linalg.matrix_power(r, 200) @ np.array([1, 0, 0])
    brute_force = brute_force / brute_force.sum()
    assert np.isclose(ngm.dominant_eigen(r).vector, brute_force).all()


def test_distribute_vaccine_even():
    N_i = np.array([1.0, 2.0, 3.0])
    V = 1.0
    n_vax = ngm.distribute_vaccines(V, N_i, strategy="even")
    assert_allclose(n_vax, np.array([1.0 / 6, 2.0 / 6, 3.0 / 6]))


def test_distribute_vaccine():
    N_i = np.array([1.0, 2.0, 3.0])
    V = 2.0
    n_vax = ngm.distribute_vaccines(V, N_i, strategy="0")
    # note that population 0 got filled, and remaining 1.0 doses were
    # distributed according to denominator 5 (i.e., of the remaining population)
    assert_allclose(n_vax, np.array([1.0, 2.0 / 5, 3.0 / 5]))


def test_distribute_vaccine_error():
    """If there is too much vaccine, error out"""
    N_i = np.array([1.0, 2.0, 3.0])

    ngm.distribute_vaccines(V=6.0, N_i=N_i, strategy="even")

    with pytest.raises(AssertionError):
        ngm.distribute_vaccines(V=10.0, N_i=N_i, strategy="even")
