from collections import namedtuple
import numpy as np


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
