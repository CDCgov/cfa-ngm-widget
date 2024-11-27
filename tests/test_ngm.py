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


def test_vax_k():
    K = np.array([[1.0, 2.0], [3.0, 4.0]])
    p_vax = np.array([0.1, 0.2])
    ve = 0.3

    current = ngm.vaccinated_K(K=K, p_vax=p_vax, ve=ve)
    expected = np.array(
        [
            [1.0 * (1.0 - 0.1 * 0.3), 2.0 * (1.0 - 0.2 * 0.3)],
            [3.0 * (1.0 - 0.1 * 0.3), 4.0 * (1.0 - 0.2 * 0.3)],
        ]
    )

    assert_array_equal(current, expected)


def test_simulate():
    n = np.array([300, 200, 100])
    n_vax = np.array([0, 11, 22])
    K = np.array([[3.0, 0.5, 0.5], [0.5, 3.0, 0.5], [0.5, 0.5, 2.0]])
    p_severe = np.array([0.0, 0.1, 0.2])
    ve = 0.8
    current = ngm.simulate(n=n, n_vax=n_vax, K=K, p_severe=p_severe, ve=ve)

    assert set(current.keys()) == {"Re", "reduced_K", "infections", "severe_infections"}
    assert np.isclose(current["Re"], 3.630075754903929)
    assert_array_equal(
        current["reduced_K"],
        np.array([[3.0, 0.478, 0.412], [0.5, 2.868, 0.412], [0.5, 0.478, 1.648]]),
    )
    assert_array_equal(
        np.round(current["infections"], 6), np.array([0.419582, 0.382363, 0.198055])
    )
    assert_array_equal(
        np.round(current["severe_infections"], 6),
        np.array([0.0, 0.038236, 0.039611]),
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
