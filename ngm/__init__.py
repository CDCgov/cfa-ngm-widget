import numpy as np


class CoreGroupNgm:
    def __init__(
        self,
        num_groups: int,
        r0_matrix: np.ndarray,
        # Other parameters go here
    ):
        self.k = self.make_k(num_groups, r0_matrix)
        # Other parameters get stored here
        raise NotImplementedError()

    def make_k(self, num_groups: int, r0_matrix: np.ndarray):
        raise NotImplementedError()

    def calculate_r_eff(self) -> float:
        r"""
        Compute $R_e$

        Returns
        -------
        float
            $R_e$ via the spectral radius of the next generation matrix
        """
        # Could consider caching eigenvalues/vectors, but always re-computing feels safer
        eval = np.linalg.eig(self.k).eigenvalues
        # eval at index (via argmax) of maximum absolute value
        # @TODO: do we need to check for imaginary components?
        return eval[np.argmax(np.abs(eval))]

    def calculate_infectious_distribution(self) -> np.ndarray:
        r"""
        Compute distribution of infections

        Returns
        -------
        np.ndarray
            The PMF on infections, I think?
        """
        # evec = np.linalg.eig(self.k).eigenvectors
        raise NotImplementedError()

    def calculate_severe_outcomes(self, inf_dist: np.ndarray):
        r"""
        Compute distribution of severe outcomes from the distribution of infections

        Parameters
        ----------
        inf_dist : np.ndarray
            Output of self.calculate_infectious_distribution()

        Returns
        -------
        np.ndarray
            The PMF on severe outcomes, I think?
        """
        raise NotImplementedError()
