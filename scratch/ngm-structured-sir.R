# function to calculate next generation for 4 group SIR model

ngm_sir <- function(N, V, K, VE, p_severe) {
  stopifnot(all(N >= V))
  stopifnot(length(N) == length(V))

  S <- N - VE * V
  K <- K * S / N

  eigenvalues <- eigen(K)
  r_effective <- max(eigenvalues$values)
  infections <- eigenvalues$vectors[, which.max(eigenvalues$values)] * N
  severe_infections <- infections * p_severe
  return(list(r_e = r_effective, infections = infections, severe_infections = severe_infections))
}

# inputs
K <- matrix(c(
  3, 1, 3, 1, # core
  1, 1, 1, 1, # kids
  3, 1, 1, 1, # travelers
  1, 1, 1, 1 # general
), nrow = 4, ncol = 4)



N <- c(100, 100, 10, 790) # pop size: core, kids, travelers, general
V <- c(10, 0, 0, 0) # doses
VE <- 0.7 # vaccine efficacy
p_severe <- c(0.03, 0.09, 0.03, 0.03)

ngm_sir(N = N, V = V, K = K, VE = VE, p_severe = p_severe)
