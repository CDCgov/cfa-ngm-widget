from collections import namedtuple
from typing import Optional

import numpy as np
import numpy.linalg as la

Eigen = namedtuple("Eigen", ["value", "vector"])


def is_irreducible(X: np.ndarray) -> bool:
    """Is a matrix irreducible?

    Args:
        X (np.ndarray): square matrix

    Returns:
        bool: is irreducible?
    """
    n = _square_n(X)
    out = la.matrix_power(np.identity(n=n) + X, n - 1)
    return bool((out > 0.0).all())


def is_diagonalizable(X: np.ndarray) -> Optional[bool]:
    """Is the matrix diagonalizable?

    If an nxn matrix has n unique eigenvalues, it is diagonalizable.
    But this is not iff, so we might return `None` to mean "don't know."

    Args:
        X (np.ndarray): square matrix

    Returns:
        bool: is diagonalizable?
    """
    n = _square_n(X)
    eigen = la.eig(X)
    if len(set(eigen.eigenvalues)) == n:
        return True
    else:
        return None


def _is_nonnegative_vector(x: np.ndarray) -> bool:
    return all(x >= 0.0) or all(x <= 0.0)


def _square_n(X: np.ndarray) -> int:
    assert X.shape[0] == X.shape[1], "Matrix must be square"
    return X.shape[0]


def dominant_eigen(X: np.ndarray) -> Eigen:
    """Dominant eigenvalue and eigenvector of a matrix

    Ensure that:
    - Dominant eigenvalue is real and positive
    - Returned eigenvector is a probability vector

    Args:
        X (np.array): matrix

    Returns:
        namedtuple: with entries `value` and `vector`
    """

    if not (X >= 0.0).all():
        raise RuntimeError("Matrix must be non-negative")

    n = _square_n(X)

    # do the eigenvalue analysis, getting all eigenvalues and eigenvectors
    eigen_all = la.eig(X)

    # find the index of the dominant eigenvalue
    # note that the i-th eigenvector is the i-th column of a matrix; i.e.,
    # eig().eigenvectors is a matrix not a list
    spectral_radius = np.max(np.abs(eigen_all.eigenvalues))
    idx = [
        i
        for i in range(n)
        if eigen_all.eigenvalues[i] == spectral_radius
        and _is_nonnegative_vector(eigen_all.eigenvectors[:, i])
    ]

    assert len(idx) == 1
    idx = idx[0]

    eigen = Eigen(
        value=eigen_all.eigenvalues[idx], vector=eigen_all.eigenvectors[:, idx]
    )

    # ensure the dominant eigenvalue and eigenvector are real and positive
    eigen = _ensure_real_eigen(eigen)
    eigen = _ensure_positive_eigen(eigen)
    # ensure eigenvector is a distribution
    assert eigen is not None
    eigen = _ensure_prob_vector_eigen(eigen)

    return eigen


def _ensure_real_eigen(e: Eigen) -> Eigen:
    """Verify that eigenvalue/vector are real-valued. Then ensure that they
    are also real-typed."""
    is_real_typed = np.isrealobj(e.value) and np.isrealobj(e.vector)
    is_complex_typed = np.iscomplexobj(e.value) and np.iscomplexobj(e.vector)
    is_real_valued = np.isreal(e.value) and all(np.isreal(e.vector))
    is_complex_valued = np.iscomplex(e.value) or any(np.iscomplex(e.vector))

    if is_real_typed:
        # if value and vector are real-typed, nothing to do
        return e
    elif is_complex_typed and is_real_valued:
        # cast from complex to real type
        return Eigen(value=np.real(e.value), vector=np.real(e.vector))
    elif is_complex_typed and is_complex_valued:
        raise RuntimeError("Complex-valued eigenvalue or eigenvector")
    else:
        raise RuntimeError("Unexpected types or values")


def _ensure_positive_eigen(e: Eigen) -> Optional[Eigen]:
    """Ensure eigenvalue and eigenvector have positive sign"""
    assert np.isrealobj(e.value)
    assert np.isrealobj(e.vector)

    positive_value = e.value > 0.0
    one_sign_vector = all(e.vector >= 0.0) or all(e.vector <= 0.0)

    if positive_value and one_sign_vector:
        return e
    elif positive_value and not one_sign_vector:
        raise RuntimeError("Eigenvector has mixed signs")
    elif not positive_value:
        raise RuntimeError("Negative eigenvalue")


def _ensure_prob_vector_eigen(e: Eigen) -> Eigen:
    """Ensure the eigenvector is a probability vector (i.e., entries sum to 1)"""
    return Eigen(value=e.value, vector=e.vector / sum(e.vector))
