from collections import namedtuple
import numpy as np
from typing import Any

DominantEigen = namedtuple("DominantEigen", ["value", "vector"])


def run_ngm(
    M_novax: np.ndarray,
    n: np.ndarray,
    n_vax: np.ndarray,
    ve: float,
) -> dict[str, Any]:
    """
    Calculate Re and distribution of infections

    Args:
        M_novax: Next Generation Matrix in the absence of administering any vaccines
        n (np.array): Population sizes for each group
        n_vax (np.array): Number of people vaccinated in each group
        ve (float): Vaccine efficacy

    Returns:
        dict: Contains dominant eigenvalue, dominant eigenvector, and adjusted NGM accounting for vaccination
    """
    n_groups = len(n)
    assert len(n_vax) == n_groups
    assert M_novax.shape[0] == n_groups
    assert M_novax.shape[1] == n_groups
    assert all(n >= n_vax), "Vaccinated cannot exceed population size"

    # eigen analysis
    M_vax = vaccinate_M(M=M_novax, p_vax=n_vax / n, ve=ve)
    eigen = dominant_eigen(M_vax, norm="L1")

    return {"M": M_vax, "Re": eigen.value, "infection_distribution": eigen.vector}


def severity(eigenvalue: float, eigenvector: np.ndarray, p_severe: np.ndarray, G: int
) -> np.ndarray:

    """
    Calculate cumulative severe infections up to and including the Gth generation.

    The first generation is the generation produced by the index case, so G = 1 includes the index
    infection (generation 0) and one generation of spread.
    All infections in every generation (including the index case) are distributed proportionately
    according to the eigenvector (defining the stationary distribution of infections) into groups,
    and for each group into severe infections according to the provided probability of severe infections.

    Args:
        eigenvalue: eigenvalue, which is number of new infections caused by each infected (Re)
        eigenvector (np.array): eigenvector, representing distribution of infections in each group
        p_severe (np.array): Probability of severe outcome in each group
        G (int): Number of generations of infections which have occurred.

    Returns:
        np.ndarray: contains total number of severe infections in each category.
    """

    return (eigenvalue ** np.arange(G + 1)).sum() * eigenvector * p_severe


def vaccinate_M(M: np.ndarray, p_vax: np.ndarray, ve: float) -> np.ndarray:
    """Adjust a next generation matrix with vaccination"""
    assert len(M.shape) == 2 and M.shape[0] == M.shape[1], "M must be square"
    n_groups = M.shape[0]
    assert len(p_vax) == n_groups, "Input dimensions must match"
    assert (0 <= p_vax).all() and (
        p_vax <= 1.0
    ).all(), "Vaccine coverage must be in [0, 1]"
    assert 0 <= ve <= 1.0

    return (M.T * (1 - p_vax * ve)).T


def dominant_eigen(X: np.ndarray, norm: str = "L1") -> DominantEigen:
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

    return DominantEigen(value=value, vector=vector)


def _ensure_positive_array(x: np.ndarray) -> np.ndarray:
    """Ensure all entries of an array are positive"""
    if all(x >= 0):
        return x
    elif all(x < 0):
        return -x
    else:
        raise RuntimeError(f"Cannot make vector all positive: {x}")


def distribute_vaccines(
    V: float, N_i: np.ndarray, strategy: str = "even"
) -> np.ndarray:
    """
    Distribute vaccines based on the specified strategy.

    Parameters:
    V (int): Number of vaccine doses.
    N_i (np.ndarray): Population sizes for each group.
    strategy (str): If "even", then distribute evenly. If a string representation of
        an integer (or just an integer), then distribute to that group first,
        and divide according to population sizes for the other groups.

    Returns:
    np.ndarray: Array of vaccine doses distributed to each group.
    """

    # Ensure V and N_i are of type float
    V = float(V)
    N_i = N_i.astype(float)

    assert V <= sum(N_i), "Can't vaccinate more people than there are in the population"

    n_groups = len(N_i)
    population_proportions = N_i / np.sum(N_i)

    if strategy == "even":
        # Distribute doses according to the proportion in each group
        n_vax = V * population_proportions
    else:
        target_indices = list(map(int, strategy.split("_")))
        N_prioritized = sum(N_i[i] for i in target_indices)

        if V <= N_prioritized:
            n_vax = np.zeros(n_groups)
            for i in target_indices:
                n_vax[i] = V * (N_i[i] / N_prioritized)
        else:
            # Fill up the selected groups
            n_vax = np.zeros(n_groups)
            for i in target_indices:
                n_vax[i] = N_i[i]
            remaining_doses = V - N_prioritized

            # Exclude the selected indices
            remaining_population = np.sum(np.delete(N_i, target_indices))
            remaining_proportions = np.where(
                np.isin(np.arange(n_groups), target_indices, invert=True),
                N_i / remaining_population,
                0.0
            )

            n_vax += remaining_doses * remaining_proportions

    assert sum(n_vax) == V
    assert len(n_vax) == n_groups

    return n_vax
