# Linear algebra of next-generation matrices

## BLUF

For many next-generation matrices $\mathbf{R}$, the population-wide reproduction number will approach the spectral radius of $\mathbf{R}$ (i.e., the absolute value of the eigenvector with the largest absolute value) and the distribution of infections will approach the corresponding eigenvector.

## Motivation: Reproduction numbers in a single population

The basic reproduction number $R_0$ is the average number of infections caused by each infected person in an otherwise fully susceptible population. The effective reproduction number $R_\mathrm{eff} \leq R_0$ is the average number of infections caused by each infected person, which is time-varying and depends on the number of remaining susceptibles.

Make two assumptions:

1. The number of infected people is much greater than 1, so that the actual number of secondary infections per infection approaches the theoretical mean $R_\mathrm{eff}$.
2. _Disease-free equilibrium_: The number of infected people is small compared to the population so that $R_\mathrm{eff}$ approaches $R_0$.

Let $I(0)$ be the number of people infected at some point in time. Under these two assumptions, the number of people in the next generation of infections will then be $I(1) = R_0 \times I(0)$. The number of infections grows exponentially with the number of generations $g$:

$$
I(g) = I(0) \times R_0^g
$$

## Next-generation matrices: reproduction numbers in multiple subpopulations

In a structured population, with multiple subpopulations, the effective reproduction number is replaced by the _next-generation matrix_ $\mathbf{R}$ with entries $R_{ij}$, which are the number of infections in subpopulation $i$ caused by an infected person in subpopulation $j$. Note that the sum of entries in column $j$, i.e. $\sum_k R_{kj}$, is the number of infections caused by an infected person in subpopulation $j$.

Make two assumptions similar to the above:

1. $I_k(0) \gg 1$, where $I_k$ is the number infected in subpopulation $k$, for each subpopulation $k$.
2. The total number of infected $\sum_k I_k(0)$ is small compared to the size of every subpopulation.

Let $\mathbf{R}_0$ be the next-generation matrix under the two assumptions above. Then the number of infections in each subpopulation in generation $g$ follows $\vec{I}(g) = \mathbf{R}_0^g \vec{I}(0)$, that is, the matrix $\mathbf{R}_0$ is applied $g$ times to the initial vector of infections $\vec{I}(0)$. Note that this is mathematically equivalent to a deterministic multi-type branching process model: each infection in each type $j$ gives rise to exactly $R_{ij}$ infections in subpopulation $i$ in the next generation.

## Eigen analysis for growth rates under a stable distribution of infections

A vector $\vec{v}$ is an _eigenvector_ of matrix $\mathbf{M}$ with _eigenvalue_ $\lambda$ if $\mathbf{M} \vec{v} = \lambda \vec{v}$. In other words, the application of a matrix to one of its eigenvectors is to simply multiply that eigenvector by its corresponding eigenvalue.

If a next-generation matrix $\mathbf{R}$ has a positive eigenvalue $\lambda$ with a corresponding non-negative eigenvector $\vec{v}$ (i.e., a vector with no negative entries), then we can interpret $\vec{v}$ as a stable distribution of infections across populations and $\lambda$ as the corresponding population-wide reproduction number, since:

$$
\mathbf{R}^g \vec{v} = \lambda^g \vec{v}
$$

## All NGMs have a dominant eigenvalue and eigenvector

NGMs are square matrices with non-negative entries. By an extension to the [Perron-Frobenius theorem](https://en.wikipedia.org/wiki/Perron%E2%80%93Frobenius_theorem), these matrices have these properties:

1. There is a _dominant eigenvalue_, which is a real, non-negative eigenvalue that is greater than or equal to the absolute value of any other eigenvalue (i.e., is equal to the spectral radius).
2. The eigenvector (or one of the eigenvectors) corresponding to the dominant eigenvalue, called the _dominant eigenvector_, is non-negative (i.e., does not have a mix of positive and negative entries).

NGMs can furthermore be _irreducible_, meaning that there is a way for an infection in any subpopulation to eventually cause an infection in every other subpopulation. (This is a sensible requirement for our analyses. If there were separate "blocks" of subpopulations that were epidemiologically independent, we could model them in separate NGMs.) All positive matrices (i.e., with no zeros) are irreducible.

Irreducible NGMs have [these additional properties](https://en.wikipedia.org/wiki/Perron%E2%80%93Frobenius_theorem#Perron%E2%80%93Frobenius_theorem_for_irreducible_non-negative_matrices):

1. The dominant eigenvalue is positive.
2. All other eigenvalues have absolute values smaller than the dominant eigenvalue.
3. The dominant eigenvector has at least one positive entry.
4. There are no other non-negative eigenvectors.

## The dominant eigenvalue is a proxy for population-wide $R_0$

If the numbers of infections $\vec{x}$ across subpopulations were equal to the dominant eigenvector, then $R_0$ is the dominant eigenvalue.

This is also about as fast as the dynamics can get, regardless of the starting vector of infections $\vec{x}$. For all square, non-negative, irreducible NGMs, the Collatz-Wielandt formula also shows that the dominant eigenvalue is equal to:

$$
\max_{\vec{x}} \min_{i, x_i \neq 0} \frac{[\mathbf{R}\vec{x}]_i}{x_i}
$$

The term $[\mathbf{R}\vec{x}]_i / x_i$ is the growth in subpopulation $i$ that will occur after one generation, starting from the vector of infection counts $\vec{x}$. The minimum is over the subpopulations (excluding those that had no infections to start with). The maximum is over all possible vectors of infection counts.

## For diagonalizable NGMs, infections tend to approach the dominant eigenvector

An NGM may also be [_diagonalizable_](https://en.wikipedia.org/wiki/Diagonalizable_matrix), in which case its eigenvectors [form a basis](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Eigendecomposition_of_a_matrix). (All positive matrices are both irreducible and diagonalizable.) In this case, any distribution of infections $\vec{x}$ can be written as a linear combination of the eigenvectors $\vec{v}_i$:

$$
\vec{x} = \sum_i \alpha_i \vec{v}_i
$$

In these cases, after $g$ generations, that vector of infections will become:

$$
\mathbf{R}_0^g \vec{x} = \sum_i \alpha_i \mathbf{R}_0^g \vec{v}_i = \sum_i \alpha_i \lambda_i^g \vec{v}_i
$$

where $\lambda_i$ are the corresponding eigenvalues.

Without loss of generality, let $\lambda_1$ be the dominant eigenvalue. Then, after a sufficiently large number $g$ of generations, the growth in the first eigenvector will outpace the others, that is, $\lambda_1^g \gg \lambda_i^g$ for any other $i \neq 1$. In that limit:

$$
\mathbf{R}_0^g \vec{x} \approx \lambda_1^g \vec{v}_1
$$

<!-- A proof of why this is reasonable to expect would be helpful here. -->
so long as we actually begin with a vector $\vec{x}$ such that $\alpha_1 > 0$. Thus, the population-wide reproduction number will approach $\lambda_1$ and the distribution of infections will approach $\vec{v}_1$.

## Caveats to this interpretation

### Disease-free equilibrium but also approaching the stable distribution

We need the number of generations to be small enough that exponential growth has not depleted a meaningful number of susceptibles, but also large enough that the population-wide reproduction number approaches the dominant eigenvalue of $\mathbf{R}$.

### Counterexample matrices

An example of a non-negative, reducible, non-diagonalizable matrix is the "L" matrix:

$$
\begin{pmatrix}
1 & 0 & 0 \\
1 & 0 & 0 \\
1 & 1 & 1
\end{pmatrix}
$$

These are curious dynamics:

- Each infection in subpopulation 1 produces 1 infection in each of the 3 subpopulations.
- Each infection in subpopulation 2 and 3 produces 1 infection in subpopulation 3.

This matrix has two distinct eigenvalues: 1 and 0. Eigenvalue 1 has two eigenvectors that are multiples of one another: $(0, 0, 1)$ and $(0, 0, -1)$. Eigenvalue 0 has eigenvector $(0, \tfrac{1}{2}, -\tfrac{1}{2})$. As per the Perron-Frobenius theorem, there is a dominant eigenvalue 1 and dominant eigenvector $(0, 0, 1)$. However, because the matrix is reducible, the dominant eigenvalue is not larger than the absolute values of all other eigenvectors.

Because the matrix is non-diagonalizable, we were not guaranteed that the eigenvectors would form a basis. In this particular case, any infection in any subpopulation would indeed lead to dynamics that approach the dominant eigenvalue and eigenvector.

An example of a non-negative and irreducible matrix with multiple eigenvalues equal to the spectral radius is:

$$
\begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix}
$$

This NGM might very roughly model, for example, a sexually-transmitted disease in which each infected man infects one woman, and each infected woman infects one man. It has two eigenvalues: $1$ and $-1$, with corresponding eigenvectors $(\tfrac{1}{2}, \tfrac{1}{2})$ and $(\tfrac{1}{2}, -\tfrac{1}{2})$. Note that, as per the Perron-Frobenius theorem, there is a single non-negative eigenvalue that is greater than the absolute value of all other eigenvalues, and the corresponding eigenvector is non-negative.

## Computational notes

There are standard algorithms for finding matrices' eigenvectors and eigenvalues (e.g., Python's [`numpy.linalg.eig`](https://numpy.org/doc/2.1/reference/generated/numpy.linalg.eig.html) and R's [`eigen`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/eigen.html)), which may yield confusing results.

### Eigenvectors can be rescaled

If $\mathbf{M}$ has eigenvector $\vec{v}$ with corresponding eigenvalue $\lambda$, then:

$$
\mathbf{M} \vec{v} = \lambda \vec{v}
$$

For any scalar $\alpha$, it follows that:

$$
\mathbf{M} (\alpha \vec{v}) = \lambda (\alpha \vec{v})
$$

Thus, eigenvectors can be rescaled (including having their signs all changed). For example, an algorithm might return an eigenvector with all negative entries, but we would consider this a "positive" eigenvector, because we can swap the signs.

The eigenvectors returned by an algorithm are likely L2-normed (i.e., the square root of the sum of squares of the entries sum to 1), to form an orthnormal basis. Because a stable _distribution_ of infections should be a probability vector (i.e., entries sum to 1), you may need to rescale the eigenvector.

### Real-valued vs. real-typed

An algorithm might return a eigenvalues or eigenvectors that are real-valued but have a complex type. Be sure to check for the value, not the type, of the dominant eigenvalue and eigenvector. For example:

```
>>> np.isreal(1+0j)
True
>>> np.isrealobj(1+0j)
False
```

### Determining if a matrix is irreducible

An $n \times n$ non-negative matrix $\mathbf{R}$ is irreducible if all entries of the matrix

$$
(\mathbb{I} + \mathbf{R})^{n-1}
$$

are positive, where $\mathbb{I}$ is the identity matrix.

### Determining if a matrix is diagonalizable

An $n \times n$ matrix is diagonalizable if it has $n$ distinct eigenvalues. This is easy to check during an eigen analysis.

## Further reading

- [_Matrix Analysis_](https://epubs.siam.org/doi/book/10.1137/1.9781611977448), which has a [free pdf](http://matrixanalysis.com/ErrataPdfFiles/Sections8.2_8.3.pdf) of the most relevant section
