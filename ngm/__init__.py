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
    assert np.all(n >= doses), "Vaccinated cannot exceed population size"
    assert len(r.shape) == 2 and r.shape[0] == r.shape[1]
    assert len(n) == len(doses) == r.shape[0], "Input dimensions must match"

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


# ###example from keeling and rohani
# n = np.array([200, 800])
# doses = np.array([0,0])
# v_e = 1
# r = np.array([[10, 0.1], [0.1, 1]])
# p_s = [0.01, 0.01]
# print(ngm_sir(n=n, doses=doses, r=r, v_e=v_e))
