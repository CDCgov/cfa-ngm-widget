from collections import namedtuple
import numpy as np


def ngm_sir(n, doses, r, v_e):
    """
    Function to calculate Re and distribution of infections for a 4-group SIR model with equal recovery rates of 1

    Args:
        n (array): Population sizes for each group
        doses (array): Number of vaccine doses administered to each group
        r (2D array): Square matrix with entries representing between and within group beta (when gamma =1 these entires are within and between group R0)
        v_e (float): Vaccine efficacy, all or nothing

    Returns:
        dict: Contains dominant eigenvalue, dominant eigenvector, and adjusted NGM accounting for vaccination
    """
    n_sus = n - doses * v_e
    n_t = n_sus.reshape(-1, 1)  # transpose
    R = r * n_t / sum(n_t)

    eigenvalues, eigenvectors = np.linalg.eig(R)
    dominant_index = np.argmax(np.abs(eigenvalues))
    dominant_eigenvalue = eigenvalues[dominant_index]
    dominant_vector = eigenvectors[:, dominant_index]
    dominant_vector_rescaled = dominant_vector / dominant_vector.sum()

    return {
        "dominant_eigenvalue": dominant_eigenvalue,
        "dominant_eigenvector": dominant_vector_rescaled,
        "ngm_adjusted": R,
    }


def dominant_eigen(X: np.array, norm: str = "L1") -> namedtuple:
    """Dominant eigenvalue and eigenvector of a matrix

    Args:
        X (np.array): matrix
        norm (str, optional): Vector norm. `np.linalg.eig()` returns
          a result with `"L2"` norm. Defaults to "L1", in which case
          the sum of the vector values is 1.

    Returns:
        namedtuple: with entries `value` and `vector`
    """
    # do the eigenvalue analysis
    e = np.linalg.eig(X)
    # which eigenvalue is the dominant one?
    i = np.argmax(e.eigenvalues)

    value = e.eigenvalues[i]
    vector = e.eigenvectors[:, i]

    if norm == "L2":
        pass
    elif norm == "L1":
        vector /= sum(vector)
    else:
        raise RuntimeError(f"Unknown norm '{norm}'")

    return namedtuple("DominantEigen", ["value", "vector"])(value=value, vector=vector)
