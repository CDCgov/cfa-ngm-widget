# Using Next Generation Matrices

The NGM for a 3-Group Infectious Disease Model with compartments $S_i$, $I_i$, $R_i$: Susceptible, Infected, and Recovered compartments.

The dynamics for $I$ in each group given by:

$$
\frac{d I_i}{dt} = \sum_{j} \frac{\beta_{ij} S_i I_j}{N} - \gamma_i I_i
$$

where:

- $\beta_{ij}$: Transmission rate from group $j$ to group $i$,
- $S_j$: Susceptible population in group $j$,
- $N_j$: Total population in group $j$,
- $\gamma_i$: Recovery rate in group $i$.

The NGM is calculated at the disease free equilibrium (DFE) where

$$
I_i = 0, S_i \approx N_i \  \text{for all\ } i
$$

The NGM $\mathbf{R}_{ij}$ is given by:

$$
\mathbf{R}_{ij} = \frac{\beta_{ij} S_{i}}{\gamma N}
$$

Where we assume all $\gamma_i = \gamma$ and SIR dynamics.

The basic reproductive number $R_0$ is calculated as the [spectral radius](https://en.wikipedia.org/wiki/Spectral_radius) (the largest absolute value of the eigenvalues) of $\mathbf{R}$.

This model incorporates vaccination by recalculating the number of susceptible individuals in each group $S_{i}^{\mathrm{vax}}$ (assuming all or nothing vaccination, with vaccine efficacy given by $\mathrm{VE}$ and the proportion of $i$ vacinated is $v_i$):

$$
\mathbf{S_{i}^\mathrm{vax}} = S_{i} - v_{i} * \mathrm{VE}
$$

So that $S_i^\mathrm{vax}$ is the population $i$ that is still susceptible post vaccination administration. Then $\mathbf{R}_{ij}$ is, accounting for vaccination,

$$
\mathbf{R}_{ij}^\mathrm{vax} = \mathbf{R}_{ij} \frac{S_{i}^\mathrm{vax}}{N_i}
$$

The spectral radius of $\mathbf{R}^\mathrm{vax}$ is $R_e$.

The distribution of infections is calculated from the eigenvector that corresponds to the spectral radius of $\mathbf{R}^\mathrm{vax}$ (specifying that the elements of said eigenvector sum to 1).

Severe infections are calculated by multiplying the proportion of infections in each group by a group-specific probability of severe infection.

## Derived quantities

Let $x$ be the L1-normed eigenvector associated with the dominant eigenvalue, and let $p_\mathrm{severe}$ be a vector representing the proportion of infections in each population that are severe. Then:

- The distribution of severe infections is the element-wise product $x \odot p_\mathrm{severe}$.
- The number of severe infections in each group per infection in the prior generation is $R_e (x \odot p_\mathrm{severe})$. The number of severe infections after $G$ generations is $R_e^G (x \odot p_\mathrm{severe})$.
- The population-wide ratio of severe infections to all infections is the dot product $x \cdot p_\mathrm{severe}$.
