# script to calculate next generation for 4 group SIR model


ngm_sir_4 <- function(N, V, params) {
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

N <- c(.1, .1, .8, 0)
V <- rep(0, 4) # prop to vaccinate c(.25, .25, .25, .25)
params <- list(
  bKK = 1, bKC = 1, bKG = 1, bKT = 1,
  bCK = 10, bCC = 10, bCG = 10, bCT = 10,
  bGK = 1, bGC = 1, bGG = 1, bGT = 1,
  bTK = 1, bTC = 10, bTG = 1, bTT = 1,
  gamma = 1, ve = 0.7
)

ngm_sir_4(N, V, params)
