from collections import namedtuple
import numpy as np
from typing import Any


def simulate(
    n: np.array, n_vax: np.array, beta: np.array, p_severe: np.array, ve: float
) -> dict[str, Any]:
    """
    Calculate Re and distribution of infections

    Args:
        n (np.array): Population sizes for each group
        n_vax (np.array): Number of people vaccinated in each group
        beta (np.array): Square matrix with entries representing contact between and within groups
        p_severe (np.array): Group-specific probability of severe infection
        ve (float): Vaccine efficacy

    Returns:
        dict: Contains dominant eigenvalue, dominant eigenvector, and adjusted NGM accounting for vaccination
    """
    n_groups = len(n)
    assert len(n_vax) == n_groups
    assert len(p_severe) == n_groups
    assert beta.shape[0] == n_groups
    assert beta.shape[1] == n_groups
    assert all(n >= n_vax), "Vaccinated cannot exceed population size"

    # eigen analysis
    R = get_R(beta = beta, n = n, n_vax = n_vax, ve = ve)
    eigen = dominant_eigen(R, norm = "L1")

    return {
        "R": R,
        "Re": eigen.value,
        "infections": eigen.vector,
        "severe_infections": eigen.value * eigen.vector * p_severe,
    }


def get_R(beta: np.array, n: np.array, n_vax: np.array, ve: float) -> np.array:
    """Adjust a next generation matrix with vaccination

    Matrix element beta_ij is the matrix of who acquires infection from whom and
    captures mixing between different groups. When recovery rate = 1, R (the next
    generation matrix) is calculated by multiplying each row by the relative size
    of the susceptible population of each group.

    The size of the susceptible population is calculated with the population size
    of each group subtracted n_vax * ve when n_vax > 1.

    Args:
        n (np.array): Population sizes for each group
        n_vax (np.array): Number of people vaccinated in each group
        beta (np.array): Square matrix with entries representing contact between and within groups
        ve (float): Vaccine efficacy

    Returns:
        np.array: matrix R from which to calculate R0
    """
    assert len(beta.shape) == 2 and beta.shape[0] == beta.shape[1], "beta must be square"
    assert beta.shape[0] == len(n_vax), "Input dimensions must match"
    assert 0 <= ve <= 1.0

    s_i = n
    s_vax = (n - n_vax * ve) / n
    return beta * (s_i / n.sum()) * s_vax


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
    i = np.argmax(np.abs(e.eigenvalues))

    value = e.eigenvalues[i]
    vector = _ensure_positive_array(e.eigenvectors[:, i])

    if not value > 0:
        raise RuntimeError(f"Negative dominant eigenvalue: {value}")
    if not all(vector >= 0):
        raise RuntimeError(f"Negative dominant eigenvector values: {vector}")

    if norm == "L2":
        pass
    elif norm == "L1":
        vector /= sum(vector)
    else:
        raise RuntimeError(f"Unknown norm '{norm}'")

    return namedtuple("DominantEigen", ["value", "vector"])(value=value, vector=vector)


def _ensure_positive_array(x: np.array) -> np.array:
    """Ensure all entries of an array are positive"""
    if all(x >= 0):
        return x
    elif all(x < 0):
        return -x
    else:
        raise RuntimeError(f"Cannot make vector all positive: {x}")
