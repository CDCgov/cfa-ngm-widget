from collections import namedtuple
from typing import Any, Optional

import numpy as np

Eigen = namedtuple("Eigen", ["value", "vector"])


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


def severity(
    eigenvalue: float, eigenvector: np.ndarray, p_severe: np.ndarray, G: int
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


def dominant_eigen(X: np.ndarray, norm: str = "L1") -> Eigen:
    """Dominant eigenvalue and eigenvector of a matrix

    Args:
        X (np.array): matrix
        norm (str, optional): Vector norm. `np.linalg.eig()` returns
          a result with `"L2"` norm. Defaults to "L1", in which case
          the sum of the vector values is 1.

    Returns:
        namedtuple: with entries `value` and `vector`
    """
    # do the eigenvalue analysis, getting all eigenvalues and eigenvectors
    eigen_all = np.linalg.eig(X)
    # which eigenvalue is the dominant one?
    i = np.argmax(np.abs(eigen_all.eigenvalues))
    # note that the i-th eigenvector is the i-th column of a matrix; i.e.,
    # eig().eigenvectors is a matrix not a list
    eigen = Eigen(value=eigen_all.eigenvalues[i], vector=eigen_all.eigenvectors[:, i])

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
                0.0,
            )

            n_vax += remaining_doses * np.array(remaining_proportions)

    assert sum(n_vax) == V
    assert len(n_vax) == n_groups

    return n_vax


def exp_growth_model_severity(R_e, inf_distribution, p_severe, G) -> np.ndarray:
    """
    Get cumulative infections and severe infections in generations 0, 1, ..., G

    Parameters:
    V (int): Number of vaccine doses.
    N_i (np.ndarray): Population sizes for each group.
    strategy (str): If "even", then distribute evenly. If a string representation of
        an integer (or just an integer), then distribute to that group first,
        and divide according to population sizes for the other groups.

    Returns:
    np.ndarray: array of infections
        [:,0] is the generation
        [:,1] is the number of infections
        [:,2] is the number of severe infections
    """
    gens = np.arange(G + 1)
    infections = np.cumsum(R_e**gens)
    severe = np.outer(infections, inf_distribution * p_severe).sum(axis=1)

    return np.stack((gens, infections, severe), 1)
