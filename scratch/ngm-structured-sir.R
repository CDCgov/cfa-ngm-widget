# function to calculate next generation for 4 group SIR model

ngm_sir <- function(N, V, R_between, R_within, VE, p_severe) {
  stopifnot(all(N >= V))
  stopifnot(length(N) == length(V))
  n_groups <- length(N)

  S <- N - VE * V
  K <- matrix(
    R_between,
    nrow = n_groups, ncol = n_groups
  )
  diag(K) <- R_within

  K <- K * S / N

  eigenvalues <- eigen(K)
  r_effective <- max(eigenvalues$values)
  infections <- eigenvalues$vectors[, which.max(eigenvalues$values)] * N
  severe_infections <- infections * p_severe
  return(list(r_e = r_effective, infections = infections, severe_infections = severe_infections))
}

# inputs

N <- c(100, 100, 10, 790) # pop size: kids, core, travelers, general
V <- c(100, 0, 0, 0) # doses
p_severe <- c(0.9, 0.3, 0.3, 0.3)
R_within <- 3 # secondary infections # should be a matrix?
R_between <- 1 # secondary infections # should be a matrix?
VE <- 0.7 # vaccine efficacy

ngm_sir(N = N, V = V, R_between = R_between, R_within = R_within, VE = VE, p_severe = p_severe)
