# Using Next Generation Matrices

The NGM for a 4-Group Infectious Disease Model with compartments $S_i$, $I_i$, $R_i$: Susceptible, Infected, and Recovered compartments.

The dynamics for $i$ in each group given by:

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

And then the NGM $\mathbf{R}_{ij}$ is given by:

$$
\mathbf{R}_{ij} = \frac{\beta_{ij} S_{i}}{\gamma N}
$$

If all $\gamma_i = \gamma$ and we assume SIR dynamics.

The basic reproductive number $R_0$ is calculated as the dominant eigenvalue of $R$.

This model incorporates vaccination by recalculating the distribution of susceptible individuals in each group $S_{i}^{vax}$ (assuming all or nothing vaccination, with vaccine efficacy given by $ve$ and the proportion of $i$ vacinated is $v_i$):

$$
\mathbf{S_{i}^\mathrm{vax}} = S_{i} - v_{i} * \mathrm{VE}
$$

So that $S_i^\mathrm{vax}$ is the population $i$ that is still susceptible post vaccination administration. Then \mathbf{R}_{ij} with vaccination factored in is given by

$$
\mathbf{R}_{ij}^{vax} = \mathbf{R}_{ij} \frac{S_{i}^\mathrm{vax}}{N_i}
$$


The dominant eigenvalue of $R$ with vaccination is $R_e$.

The distribution of infections is calculated from the dominant eigenvector of $R$ (specifying that its elements add to 1).

Severe infections are calculated by multiplying the proportion of infections in each group by a group-specific probability of severe infection.
