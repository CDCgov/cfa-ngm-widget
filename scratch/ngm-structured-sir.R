# function to calculate next generation for 4 group SIR model
ngm_sir_4 <- function(N, V, b_mult, ve, gam) {
  low <- 0.1
  high <- b_mult * low

  params <- list(
    bKK = low, bKC = low, bKG = low, bKT = low,
    bCK = b_mult, bCC = b_mult, bCG = b_mult, bCT = b_mult,
    bGK = low, bGC = low, bGG = low, bGT = low,
    bTK = low, bTC = b_mult, bTG = low, bTT = low,
    gamma = gam, ve = ve
  )
  istates <- c("IK", "IC", "IG", "IT")
  flist <- c(
    dIKdt = quote(bKK * SK * IK + bKC * SK * IC + bKG * SK * IG + bKT * SK * IT),
    dICdt = quote(bCK * SC * IK + bCC * SC * IC + bCG * SC * IG + bCT * SC * IT),
    dIGdt = quote(bGK * SG * IK + bGC * SG * IC + bGG * SG * IG + bGT * SG * IT),
    dITdt = quote(bTK * ST * IK + bTC * ST * IC + bTG * ST * IG + bTT * ST * IT)
  )
  V1 <- quote(gamma * IK)
  V2 <- quote(gamma * IC)
  V3 <- quote(gamma * IG)
  V4 <- quote(gamma * IT)
  vlist <- c(V1, V2, V3, V4)
  df <- list(
    SK = N[1] * (1 - V[1] * ve), SC = N[2] * (1 - V[2] * ve),
    SG = N[3] * (1 - V[3] * ve), ST = N[4] * (1 - V[4] * ve),
    IK = 0, IC = 0, IG = 0, IT = 0,
    RK = V[1] * params$ve, RC = V[2] * params$ve, RG = V[3] * params$ve, RT = V[4] * params$ve
  )

  nextgenR0(Istates = istates, Flist = flist, Vlist = vlist, parameters = params, dfe = df)
}

# inputs

N <- c(.100, .100, .800, 0)
V <- c(0, 0, 0, 0)
b_mult <- 3
ve <- 0.7
gam <- 1 / 7

ngm_sir_4(N = N, V = V, b_mult = b_mult, ve = ve, gam = gam)
