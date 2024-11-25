# Next generation matrix widget

⚠️ This is a work in progress

## Overview

This repo contains code for investigating the potential efficacy of
vaccination allocation for a disease interactively via a
[streamlit](https://streamlit.io/) app using a next generation matrix
approach.


## Model Description

This repo contains code to apply the next-generation method of Diekman et al. (1990) to calculate R0 for an SIR model with 4 risk groups and flexible inputs for varying vaccine allocation to the 4 groups.

### Next generation matrix calculation

We calculate the next-generation matrix for a 4-Group Infectious Disease Model with compartments $S_i$, $I_i$, $R_i$: Susceptible, Infected, and Recovered compartments in group $i$, where $i = 1, 2, 3, 4$. Transmission dynamics for $I_i$ in each group given by:

$$dI_i/dt = Σ (β_ij * S_i * I_j / N_j) - γ * I_i$$

where:
- $β_ij$: Transmission rate from group $j$ to group $i$,
- $S_i$: Susceptible population in group $i$,
- $N_j$: Total population in group $j$,
- $γ$: Recovery rate (same for all groups).

The NGM is calculated at the disease free equilibrium (DFE) where

$$I_1 = I_2 = I_3 = I_4 = 0$$
$$S_i = N_i$$

And then the NGM $K$ is given by:

$$K = F \cdot V^-1$$

where $F$ is the matrix of new infections and V is the matrix of transitions between compartments, not representing new infections:

$$ F = {\left\lbrack \matrix{
β_11 * S_1 / N_1 & β_12 * S_1 / N_2 & β_13 * S_1 / N_3 & β_14 * S_1 / N_4 \cr
β_21 * S_2 / N_1 & β_22 * S_2 / N_2 & β_23 * S_2 / N_3 & β_24 * S_2 / N_4 \cr
β_31 * S_3 / N_1 & β_32 * S_3 / N_2 & β_33 * S_3 / N_3 & β_34 * S_3 / N_4 \cr
β_41 * S_4 / N_1 & β_42 * S_4 / N_2 & β_43 * S_4 / N_3 & β_44 * S_4 / N_4} \right\rbrack} $$

$$ V = {\left\lbrack \matrix{
γ & 0 & 0 & 0 \cr
0 & γ & 0 & 0 \cr
0 & 0 & γ & 0 \cr
0 & 0 & 0 & γ} \right\rbrack} $$

Since $V$ is diagonal, its inverse is:

$$ V^-1 = {\left\lbrack \matrix{
1/γ & 0 & 0 & 0 \cr
0 & 1/γ & 0 & 0 \cr
0 & 0 & 1/γ & 0 \cr
0 & 0 & 0 & 1/γ} \right\rbrack} $$

Multiplying $F$ and $V^-1$:

$$ K = {\left\lbrack \matrix{
    β_{11} * S_1 / (γ * N_1) & β_{12} * S_1 / (γ * N_2) & β_{13} * S_1 / (γ * N_3) & β_{14} * S_1 / (γ * N_4) \cr
    β_{21} * S_2 / (γ * N_1) & β_{22} * S_2 / (γ * N_2) & β_{23} * S_2 / (γ * N_3) & β_{24} * S_2 / (γ * N_4) \cr
    β_{31} * S_3 / (γ * N_1) & β_{32} * S_3 / (γ * N_2) & β_{33} * S_3 / (γ * N_3) & β_{34} * S_3 / (γ * N_4) \cr
    β_{41} * S_4 / (γ * N_1) & β_{42} * S_4 / (γ * N_2) & β_{43} * S_4 / (γ * N_3) & β_{44} * S_4 / (γ * N_4)} \right\rbrack} $$

Note that each entry $K_{ij}$ represents the expected number of secondary infections, $R_{ij}$ in group $i$ caused by an infected individual in group $j$:

$$ K = {\left\lbrack \matrix{
    R_{11} / N_1 & R_{12} / N_2 & R_{13} / N_3 & R_{14} / N_4 \cr
    R_{21} / N_1 & R_{22} / N_2 & R_{23} / N_3 & R_{24} / N_4 \cr
    R_{31} / N_1 & R_{32} / N_2 & R_{33} / N_3 & R_{34} / N_4 \cr
    R_{41} / N_1 & R_{42} / N_2 & R_{43} / N_3 & R_{44} / N_4} \right\rbrack} $$

The effective reproductive number is calculated as the dominant eigenvalue of $K$ (when the population has vaccination?).

The distribution of infections is calculated from the dominant eigenvector of $K$.

Severe infections are calculated by multiplying the proportion of infections in each group by a group-specific probability of severe infection.

### Model assumptions

Vaccination is assumed to be all or nothing  -- each individual's immunity is determined by a coin flip with probability of being immune equal to the vaccine efficacy. In the NGM, each $N_i = N_i * (1 - VE)$.

# Streamlit app

We make the following assumptions about transmission between and within groups for the default parameters:

| Group | Probability of onward, within group transmission | Probability of being infected from outside group | Probability of severe outcome in group |
|----------|----------|----------|----------|
| Core      | High     | High     | Low     |
| Children  | Low      | Low      | High    |
| Travelers | Low      | High     | Low     |
| General   | Low      | Low      | Low     |

Inputs:

* Sizes of the groups, $N_i$
* Vaccination efficacy (VE) and number of doses allocated to each group ($V_i$)
* Within and between group reproduction numbers; the entries to the NGM, $R_{ij}$
* Per-group probability of severe infection

Outputs:

* Effective reproductive number
* Distribution of infections and severe infections

### References

Diekmann O, Heesterbeek JA, Metz JA. On the definition and the computation of the basic reproduction ratio R0 in models for infectious diseases in heterogeneous populations. J Math Biol. 1990;28(4):365-82. doi: 10.1007/BF00178324. PMID: 2117040.

Diekmann O, Heesterbeek JA, Roberts MG. The construction of next-generation matrices for compartmental epidemic models. J R Soc Interface. 2010 Jun 6;7(47):873-85. doi: 10.1098/rsif.2009.0386. Epub 2009 Nov 5. PMID: 19892718; PMCID: PMC2871801.

 van den Driessche P, Watmough J. Reproduction numbers and sub-threshold endemic equilibria for compartmental models of disease transmission. Math Biosci. 2002 Nov-Dec;180:29-48. doi: 10.1016/s0025-5564(02)00108-6. PMID: 12387915.


## Authors

- Paige Miller <yub1@cdc.gov>
- Scott Olesen <ulp7@cdc.gov>
- Andy Magee <rzg0@cdc.gov>

## General Disclaimer

This repository was created for use by CDC programs to collaborate on public
health related projects in support of the
[CDC mission](https://www.cdc.gov/about/organization/mission.htm). GitHub is not
hosted by the CDC, but is a third party website used by CDC and its partners to
share information and collaborate on software. CDC use of GitHub does not imply
an endorsement of any one particular service, product, or enterprise.

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is
in the public domain within the United States, and copyright and related rights
in the work worldwide are waived through the
[CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication.
By submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice

This repository is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or
modify it under the terms of the Apache Software License version 2, or (at your
option) any later version.

This source code in this repository is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software
License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md) and
[Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit
[http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by
[forking](https://help.github.com/articles/fork-a-repo) and submitting a pull
request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law,
including but not limited to the Federal Records Act, and may be archived. Learn
more at
[http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice

This repository is not a source of government records but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).
