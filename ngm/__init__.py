from typing import Any

import numpy as np
import numpy.random

import ngm.linalg


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
    eigen = ngm.linalg.dominant_eigen(M_vax)

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
    assert (0 <= p_vax).all() and (p_vax <= 1.0).all(), (
        "Vaccine coverage must be in [0, 1]"
    )
    assert 0 <= ve <= 1.0

    return (M.T * (1 - p_vax * ve)).T


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


def run_stochastic(
    M: np.ndarray,
    x: np.ndarray,
    max_generations: int,
    rng: None | int | np.random.Generator = None,
) -> np.ndarray:
    """
    Run a stochastic simulation using a given next generation matrix and initial conditions.

    Args:
        M (np.ndarray): Next generation matrix
        x (np.ndarray): Initial number of infections
        max_generations (int): Maximum number of generations to simulate
        rng (None | int | np.random.Generator): Random number generator or seed. If
            None, uses the default random number generator. If an int, uses it as a seed.
            If a Generator, uses it directly. Default is None.

    Returns:
        np.ndarray: Array of infections over time
    """

    n = M.shape[0]
    assert M.shape[1] == n, "M must be square"
    assert len(x) == n, "Input dimensions must match"
    assert (x >= 0).all(), "Initial infections must be non-negative"
    assert x.sum() > 0, "Must start with nonzero number of infections"
    assert max_generations > 0, "Must simulate at least one generation"

    if rng is None:
        rng = np.random.default_rng()
    elif isinstance(rng, int):
        rng = np.random.default_rng(rng)
    elif isinstance(rng, np.random.Generator):
        pass
    else:
        raise ValueError(f"rng {rng} is not one of None, int, or np.random.Generator")

    out = np.zeros((max_generations, n))
    out[0, :] = x

    for g in range(1, max_generations):
        rates = M @ out[g - 1, :]
        x = rng.poisson(lam=rates, size=None)

        if x.sum() == 0:
            break

        out[g, :] = x

    return out
