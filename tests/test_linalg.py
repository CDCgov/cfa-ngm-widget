import numpy as np
import pytest
from numpy.testing import assert_array_equal

import ngm.linalg


class TestDominantEigen:
    def test_dominant_eigen_simple(self):
        X = np.array([[1, 2], [2, 1]])
        e = ngm.linalg.dominant_eigen(X)
        assert np.isclose(e.value, 3.0)
        assert np.isclose(e.vector, np.array([0.5, 0.5])).all()

    def test_dominant_eigen_bigger(self):
        X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        e = ngm.linalg.dominant_eigen(X)
        assert np.isclose(e.value, 16.116843969807043)
        assert np.isclose(e.vector, np.array([0.14719267, 1.0 / 3, 0.51947399])).all()

    def test_kr(self):
        beta = np.array([[10, 0.1], [0.1, 1]])
        n = np.array([0.2, 0.8])
        R = (beta.T * n).T

        # see Keeling & Rohani, page 61
        R_p61 = np.array([[2, 0.02], [0.08, 0.8]])
        assert np.isclose(R, R_p61).all()

        r0 = ngm.linalg.dominant_eigen(R).value

        assert np.isclose(r0, 2.0013, atol=5e-5)

    def test_eigenvectors(self):
        r = np.array([[3.1, 0.15, 1.7], [0.78, 1.5, 0.1], [0.32, 0.98, 1.1]])

        brute_force = np.linalg.matrix_power(r, 200) @ np.array([1, 0, 0])
        brute_force = brute_force / brute_force.sum()
        assert np.isclose(ngm.linalg.dominant_eigen(r).vector, brute_force).all()

    def test_eigen_returns_real(self):
        M = np.array([[3.0, 0.0, 0.2], [40.0, 1.0, 500], [0.25, 1.0, 1.5]])
        eigen = ngm.linalg.dominant_eigen(M)
        assert np.isreal(eigen.value)
        assert all(np.isreal(eigen.vector))

    def test_eigen_returns_error(self):
        M = np.array([[0, -1], [1, 0]])
        with pytest.raises(RuntimeError, match="non-negative"):
            ngm.linalg.dominant_eigen(M)

    def test_weird_shape(self):
        M = np.array([[0, 1, 0], [1, 1, 0], [0, 0, 0]])
        eig = ngm.linalg.dominant_eigen(M)
        assert np.isclose(eig.value, 1.6180, rtol=0.0, atol=1e-3)
        assert np.allclose(
            eig.vector, np.array([[0.3820, 0.6180, 0.0]]), rtol=0.0, atol=1e-3
        )


class TestEnsureReal:
    def test_trivial(self):
        e = ngm.linalg.Eigen(value=2.0, vector=np.array([1, 2, 3]))
        assert ngm.linalg._ensure_real_eigen(e) == e

    def test_success(self):
        e_in = ngm.linalg.Eigen(
            value=2.0 + 0j, vector=np.array([1 + 0j, 2 + 0j, 3 + 0j])
        )
        e_out = ngm.linalg._ensure_real_eigen(e_in)

        assert e_out.value == 2.0
        assert_array_equal(e_out.vector, np.array([1, 2, 3]))

    def test_fail_complex_value(self):
        with pytest.raises(RuntimeError, match="Complex-valued"):
            ngm.linalg._ensure_real_eigen(
                ngm.linalg.Eigen(value=0 + 1j, vector=np.array([0 + 0j]))
            )

    def test_fail_complex_vector(self):
        with pytest.raises(RuntimeError, match="Complex-valued"):
            ngm.linalg._ensure_real_eigen(
                ngm.linalg.Eigen(value=1 + 0j, vector=np.array([0 + 1j]))
            )


class TestProbVectorEigen:
    def test_simple(self):
        e = ngm.linalg.Eigen(value=2.0, vector=np.array([1, 2, 3]))
        out = ngm.linalg._ensure_prob_vector_eigen(e)
        assert_array_equal(out.vector, np.array([1 / 6, 2 / 6, 3 / 6]))

    def test_flip(self):
        e = ngm.linalg.Eigen(value=2.0, vector=np.array([-1, -2, -3]))
        out = ngm.linalg._ensure_prob_vector_eigen(e)
        assert_array_equal(out.vector, np.array([1 / 6, 2 / 6, 3 / 6]))


class TestEnsurePositive:
    def test_trivial(self):
        e = ngm.linalg.Eigen(value=1, vector=np.array([1, 2, 3]))
        assert ngm.linalg._ensure_positive_eigen(e) == e

    def test_fail_mixed(self):
        with pytest.raises(RuntimeError, match="mixed"):
            ngm.linalg._ensure_positive_eigen(
                ngm.linalg.Eigen(value=1, vector=np.array([1, -1]))
            )

    def test_fail_negative_value(self):
        with pytest.raises(RuntimeError, match="Negative eigenvalue"):
            ngm.linalg._ensure_positive_eigen(
                ngm.linalg.Eigen(value=-1, vector=np.array([1]))
            )


class TestIsIrreducible:
    def test_simple_true(self):
        assert ngm.linalg.is_irreducible(np.array([[1, 2], [3, 4]]))

    def test_simple_false(self):
        L_matrix = np.array([[1, 0, 0], [1, 0, 0], [1, 1, 1]])
        assert not ngm.linalg.is_irreducible(L_matrix)


class TestIsDiagonalizable:
    def test_simple_true(self):
        assert ngm.linalg.is_diagonalizable(np.array([[1, 2], [3, 4]])) is True

    def test_simple_dontknow(self):
        L_matrix = np.array([[1, 0, 0], [1, 0, 0], [1, 1, 1]])
        assert ngm.linalg.is_diagonalizable(L_matrix) is None
