# Linear algebra of next-generation matrices

## BLUF

For many next-generation matrices $\mathbf{R}$, the population-wide reproduction number will approach the spectral radius of $\mathbf{R}$ (i.e., the absolute value of the eigenvector with the largest absolute value) and the distribution of infections will approach the corresponding eigenvector.

## Motivation: Reproduction numbers in a single population

The basic reproduction number $R_0$ is the average number of infections caused by each infected person in an otherwise fully susceptible population. The effective reproduction number $R_\mathrm{eff} < R_0$ is the average number of infections caused by each infected person, which is time-varying and depends on the number of remaining susceptibles.

Assume the population is large compared to the number of infections. In this limit, stochastic effects are unimportant and the actual number of infections caused by each infected person approaches $R_0$. This is the _disease-free equilibrium_ approximation: to model early dynamics of an outbreak, we assume that, to a first-order approximation, the number of susceptibles is not declining, so that $R_\mathrm{eff} \approx R_0$.

Let $I(0)$ be the number of people infected at some point in time. The number of people in the next generation of infections will then be $I(1) = R_0 \times I(0)$. The number of infections grows exponentially with the number of generations $g$:

$$
I(g) = I(0) \times R_0^g
$$

## Next-generation matrices: reproduction numbers in multiple population

In a structured population, with multiple subpopulations, the effective reproduction number is replaced by the _next-generation matrix_ $\mathbf{R}$ with entries $R_{ij}$, which are the number of infections in subpopulation $i$ caused by an infected person in subpopulation $j$. Note that the sum of entries in column $j$, i.e. $\sum_k R_{kj}$, is the number of infections caused by an infected person in subpopulation $j$.

Let $\mathbf{R}_0$ be the next-generation matrix at the disease-free equilibrium. In this limit, the number of infections in each subpopulation in generation $g$ follows $\vec{I}(g) = \mathbf{R}_0^g \vec{I}(0)$, that is, the matrix $\mathbf{R}_0$ is applied $g$ times to the initial vector of infections $\vec{I}(0)$. Note that this is mathematically equivalent to a deterministic multi-type branching process model: each infection in each type $j$ gives rise to exactly $R_{ij}$ infections in subpopulation $i$ in the next generation.

## Eigen analysis for growth rates under a stable distribution of infections

A vector $\vec{v}$ is an _eigenvector_ of matrix $\mathbf{M}$ with _eigenvalue_ $\lambda$ if $\mathbf{M} \vec{v} = \lambda \vec{v}$. In other words, the application of a matrix to one of its eigenvectors is to simply multiply that eigenvector by its corresponding eigenvalue.

If a next-generation matrix $\mathbf{R}$ has a nonnegative eigenvector $\vec{v}$ (i.e., with no negative entries) with a positive eigenvalue $\lambda$, then we can interpret $\vec{v}$ as a stable distribution of infections across populations and $\lambda$ as the corresponding population-wide reproduction number $R_0$, since:

$$
\mathbf{R}^g \vec{v} = \lambda^g \vec{v}
$$

## Infections tend to approach the dominant eigenvalue and eigenvector

In certain cases (e.g., if $\mathbf{R}$ is a [positive semi-definite matrix](https://en.wikipedia.org/wiki/Definite_matrix)):

1. the matrix $\mathbf{R}$ has a number of eigenvalues and eigenvectors equal to the number of populations,
2. the eigenvalues are all nonnegative, and
3. the eigenvectors will form a _basis_, such that any distribution $\vec{x}$ of infections can be written as a linear combination of the eigenvectors $\vec{v}_i$:

$$
\vec{x} = \sum_i \alpha_i \vec{v}_i
$$

In these cases, after $g$ generations, that vector of infections will become:

$$
\mathbf{R}_0^g \vec{x} = \sum_i \alpha_i \mathbf{R}_0^g \vec{v}_i = \sum_i \alpha_i \lambda_i^g \vec{v}_i
$$

where $\lambda_i$ are the corresponding eigenvalues.

Without loss of generality, let $\lambda_1$ be the largest eigenvalue. Then, after a sufficiently large number $g$ of generations, the growth in the first eigenvector will outpace the others: $\lambda_1^g \gg \lambda_i^g$ for any other $i \neq 1$. In that limit:

$$
\mathbf{R}_0^g \vec{x} \approx \lambda_1^g \vec{v}_1
$$

so long as $\alpha_1 > 0$. Thus, the population-wide reproduction number will approach $\lambda_1$ and the distribution of infections will approach $\vec{v}_1$.

## Caveats to this interpretation

### Eigenvectors can be rescaled

There are standard algorithms for finding matrices' eigenvectors and eigenvalues (e.g., Python's [`numpy.linalg.eig`](https://numpy.org/doc/2.1/reference/generated/numpy.linalg.eig.html) and R's [`eigen`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/eigen.html)), which may yield confusing results.

If $\mathbf{M}$ has eigenvector $\vec{v}$ with corresponding eigenvalue $\lambda$, then:

$$
\mathbf{M} \vec{v} = \lambda \vec{v}
$$

For any scalar $\alpha$, it follows that:

$$
\mathbf{M} (\alpha \vec{v}) = \lambda (\alpha \vec{v})
$$

Thus, eigenvectors can be rescaled (including having their signs all changed). To be interpretable for an NGM, the dominant eigenvalue must be real and positive, and all entries of the eigenvector must have the same sign. If they are all negative, we can simply swap them all to positive.

The eigenvectors returned by an algorithm are likely L2-normed (i.e., the square root of the sum of squares of the entries sum to 1), to form an orthnormal basis. Because a stable _distribution_ of infections should be a probability vector (i.e., entries sum to 1), you may need to rescale the eigenvector.

### Disease-free equilibrium but also sufficient generations

We need the number of generations to be small enough that exponential growth has not depleted a meaningful number of susceptibles, but also large enough that the population-wide reproduction number approaches the spectral radius of $\mathbf{R}$.

### Matrices that don't work nicely

In the case of diagonal (or block-diagonal) NGM, it is easy to see why the $\alpha_1 > 0$ caveat matters: if you have totally independent subpopulations, and don't start with an infection in the subpopulation with the largest reproduction number, the population-wide dynamics will not approach the dynamics of that subpopulation.

Not all relevant matrices are positive semi-definite. Consider a toy model of a sexually transmitted disease, in which each male infects exactly one female, and each female infects exactly one male:

$$
\mathbf{R} = \begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix}
$$

This matrix has two eigenvectors: $[\tfrac{1}{2}, \tfrac{1}{2}]$ with eigenvalue $1$, and $[\tfrac{1}{2}, -\tfrac{1}{2}]$ with eigenvalue $-1$. Both eigenvalues have the same absolute value, but only the positive eigenvector and eigenvalue are interpretable for us. We should filter only for positive eigenvalues and eigenvectors.

The matrix

$$
\begin{pmatrix}
2 & 1 \\
1 & 1
\end{pmatrix}
$$

has two eigenvalues, but the dominant one has an eigenvector with mixed positive and negative values.

## Eigenvalues and eigenvectors can be complex

[Example of a matrix that has a sensible dominant eigenvalue but then has complex eigenvectors.]
