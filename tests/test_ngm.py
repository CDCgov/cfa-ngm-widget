import numpy as np
import ngm
import numpy.testing


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

    numpy.testing.assert_array_equal(current, expected)
