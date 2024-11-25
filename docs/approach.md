# Using Next Generation Matrices

The NGM for a 4-Group Infectious Disease Model with compartments $S_g$, $I_g$, $R_g$: Susceptible, Infected, and Recovered compartments in the $G$ groups $g = 1, \dots, G$. Where possible we try to be general, otherwise we will take $K = 4$

Dynamics for $I_k$ in each group given by:

$$
\frac{d I_g}{dt} = \sum_{j} \frac{\beta_{jg} S_j I_j}{N_j} - \gamma_g I_g
$$

where:

- $\beta_{jg}$: Transmission rate from group $j$ to group $g$,
- $S_j$: Susceptible population in group $j$,
- $N_j$: Total population in group $j$,
- $\gamma_g$: Recovery rate in group $g$.

The NGM is calculated at the disease free equilibrium (DFE) where

$$
I_g = 0, S_g = N_g \  \text{for all\ } g
$$

---

And then the NGM `K` is given by:

$$
\mathbf{K} = \mathbf{F} \mathbf{V}^{-1}
$$

where $\mathbf{F}$ is the matrix of new infections and $\mathbf{V}$ is the matrix of transitions between compartments, not representing new infections.

The elements of $\mathbf{F}$ are

$$
\mathbf{F}_{ij} = \frac{\beta_{ij} S_j}{N_j}
$$

while $\mathbf{V}$ is a diagonal matrix with $\mathbf{V}_{ii} = \gamma$ where the recovery rate is shared among all groups $g$.

Since $\mathbf{V}$ is diagonal, its inverse is as well, with $(\mathbf{V}^{-1})_{ij} = 1 / \gamma_i$.

Thus, $\mathbf{K}$ is given by

$$
\mathbf{K}_{ij} = \frac{\beta_{ij} S_j}{\gamma_j N_j}
$$

which we can re-write as

$$
\mathbf{K}_{ij} = R_{ij} \frac{S_j}{N_j}
$$

The basic reproductive number $R_0$ is calculated as the dominant eigenvalue of $K$, while $R_0$ multiplied by the fraction of susceptibles yields the effective reproduction number $R_e$.

The distribution of infections is calculated from the dominant eigenvector of $K$.

Severe infections are calculated by multiplying the proportion of infections in each group by a group-specific probability of severe infection.
