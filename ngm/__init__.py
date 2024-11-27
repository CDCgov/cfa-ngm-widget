import numpy as np


def ngm_sir(n, v, k, v_e, p_s):
    """
    Function to calculate Re and distribution of infections for a 4-group SIR model

    Args:
        n (array): Population sizes for each group
        v (array): Number of vaccine doses administered to each group
        k (2D array): Square matrix with entries representing between and within group R0
        v_e (float): Vaccine efficacy, all or nothing
        p_s (array): Group specific probability of severe infection

    Returns:
        dict: Contains R-effective, distribution of infections, severe infections, adjusted K matrix
    """
    # check sizes of inputs? do we need this
    assert np.all(n >= n), "Vaccinated cannot exceed population size"
    assert len(k.shape) == 2 and k.shape[0] == k.shape[1]
    assert len(n) == len(v) == len(p_s) == k.shape[0], "Input dimensions must match"

    # Calculate susceptibles for NGM
    s = n - v_e * v
    k_adjusted = k * s / n

    # Eigenvalue computation for NGM
    eigenvalues, eigenvectors = np.linalg.eig(k_adjusted)
    r_effective = np.max(np.abs(eigenvalues))

    # Calculate distribution of infections and severe infections
    dominant_vector = eigenvectors[:, np.argmax(np.abs(eigenvalues))]
    infections = dominant_vector / dominant_vector.sum()
    severe_infections = infections * p_s

    return {
        "R_effective": r_effective,
        "Infections": infections.real,
        "Severe_Infections": severe_infections.real,
        "K_adjusted": k_adjusted,
    }
