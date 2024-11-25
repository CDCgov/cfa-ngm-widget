from collections import namedtuple

import numpy as np

NgmSummary = namedtuple(
    "NgmSummary", ["r_eff", "infectious_dist", "outcome_dist"]
)


class CoreGroupNgm:
    def __init__(
        self,
        r0_matrix: np.ndarray,
        # Other parameters go here
    ):
        assert (
            len(r0_matrix) == 2 and r0_matrix.shape[0] == r0_matrix.shape[1]
        ), "r0_matrix must be square 2-D matrix"
        self.k = self.make_k(r0_matrix)
        # Other parameters get stored here
        raise NotImplementedError()

    def make_k(self, r0_matrix: np.ndarray):
        raise NotImplementedError()

    def summarize(self) -> NgmSummary:
        raise NotImplementedError

        # eigendecomp = np.linalg.eig(self.k)

        # # Index of dominant eigenvalue determines r_eff and infectious distribution
        # dom = np.argmax(np.abs(eigendecomp.eigenvalues))

        # return NgmSummary(
        #     r_eff=eigendecomp.eigenvalues[dom],
        #     infectious_dist=None, #TBD: is this just eigendecomp.eigenvectors[dom] / eigendecomp.eigenvectors[dom].sum()?
        #     outcome_dist=None, #TBD
        # )
