import numpy as np
import ngm
from numpy.testing import assert_array_equal
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


def test_vax_beta():
    beta = np.array([[10.0, 0.1], [0.1, 1.0]])
    n = np.array([200, 800])
    n_vax = np.array([100, 0])
    ve = 1.0

    current = ngm.get_R(beta=beta, n=n, n_vax=n_vax, ve=ve)
    expected = np.array(
        [
            [10.0 * 1 / 9, 0.1 * 8/9],
            [0.1 * 1/9, 1 * 8/9],
        ]
    )

    assert_array_equal(current, expected)


def test_simulate():
    n = np.array([200, 200, 100, 500])
    n_vax = np.array([0, 0, 0, 0])
    beta = np.array([[3.0, 0.5, 3, 0.5],
                     [0.5, 0.5, 0.5, 0.5],
                     [3, 0.5, 0.5, 0.5],
                     [0.5, 0.5, 0.5, 0.5]])
    p_severe = np.array([0.02, 0.06, 0.02, 0.02])
    ve = 1.0
    current = ngm.simulate(n=n, n_vax=n_vax, beta=beta, p_severe=p_severe, ve=ve)

    assert set(current.keys()) == {"Re", "R", "infections", "severe_infections"}
    assert np.isclose(current["Re"], 0.92)
    assert_array_equal(
        current["R"],
        np.array([[0.6, 0.1, 0.6, 0.1],
                  [0.1, 0.1, 0.1, 0.1],
                  [0.3, .05, 0.05, 0.05],
                  [0.25, 0.25, 0.25, 0.25]]),
    )
    assert_array_equal(
        np.round(current["infections"], 6), np.array([0.44507246, 0.10853944, 0.17503951, 0.2713486])
    )
    assert_array_equal(
        np.round(current["severe_infections"], 6),
        np.array([0.00820112, 0.006, 0.00322536, 0.005]),
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
