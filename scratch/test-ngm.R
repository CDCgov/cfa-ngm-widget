# script to calculate next generation for some basic models

library(epimdr2)

##### SEIR model ###########
# example from package vignette
istates <- c("E", "I")
flist <- c(dEdt = quote(beta * S * I / N), dIdt = quote(0))
Vm1 <- quote(mu * E + sigma * E)
Vm2 <- quote(mu * I + alpha * I + gamma * I)
Vp1 <- 0
Vp2 <- quote(sigma * E)
V1 <- substitute(a - b, list(a = Vm1, b = Vp1))
V2 <- substitute(a - b, list(a = Vm2, b = Vp2))
vlist <- c(V1, V2)
params <- list(mu = 0, alpha = 0, beta = 5, gamma = .8, sigma = 1.2, N = 1)
df <- list(S = 1, E = 0, I = 0, R = 0)
nextgenR0(
  Istates = istates, Flist = flist,
  Vlist = vlist, parameters = params, dfe = df
)
params$beta / params$gamma

###### SIR Model with two risk groups (high and low) #########
# example from keeling and rohani
istates <- c("IH", "IL")
flist <- c(
  dIHdt = quote(bHH * SH * IH + bHL * SH * IL),
  dILdt = quote(bLH * SL * IH + bLL * SL * IL)
)
Vm1 <- quote(gamma * IH)
Vm2 <- quote(gamma * IL)
Vp1 <- 0
Vp2 <- 0
V1 <- substitute(a - b, list(a = Vm1, b = Vp1))
V2 <- substitute(a - b, list(a = Vm2, b = Vp2))
vlist <- c(V1, V2)
params <- list(
  bHH = 10, bHL = 0.1,
  bLH = 0.1, bLL = 1,
  gamma = 1
)
df <- list(
  SH = 0.2, SL = 0.8,
  IH = 0, IL = 0,
  RH = 0, RL = 0
)
nextgenR0(Istates = istates, Flist = flist, Vlist = vlist, parameters = params, dfe = df)

J <- matrix(c(
  params$bHH * df$SH - params$gamma, params$bHL * df$SH,
  params$bLH * df$SL, params$bLL * df$SL - params$gamma
), nrow = 2)
# when params$gamma = 1, R0=
(params$bHH * df$SH + params$bLH * df$SL) * (1 - (norm(eigen(J)$vectors) - 1)) +
  (params$bHL * df$SH + params$bLL * df$SL) * (norm(eigen(J)$vectors) - 1)
(params$bHH * df$SH + params$bLH * df$SL) * .9376 +
  (params$bHL * df$SH + params$bLL * df$SL) * .0624

####### SIR Model with 3 classes ######
# K, C, G: kids, core, general
istates <- c("IK", "IC", "IG")
flist <- c(
  dIKdt = quote(bKK * SK * IK + bKC * SK * IC + bKG * SK * IG),
  dICdt = quote(bCK * SC * IK + bCC * SC * IC + bCG * SC * IG),
  dIGdt = quote(bGK * SG * IK + bGC * SG * IC + bGG * SG * IG)
)
Vm1 <- quote(gamma * IK)
Vm2 <- quote(gamma * IC)
Vm3 <- quote(gamma * IG)
Vp1 <- 0
Vp2 <- 0
Vp3 <- 0
V1 <- substitute(a - b, list(a = Vm1, b = Vp1))
V2 <- substitute(a - b, list(a = Vm2, b = Vp2))
V3 <- substitute(a - b, list(a = Vm3, b = Vp3))
vlist <- c(V1, V2, V3)
params <- list(
  bKK = 10, bKC = 1, bKG = 1,
  bCK = 1, bCC = 10, bCG = 1,
  bGK = 1, bGC = 1, bGG = 1,
  gamma = 1
)
df <- list(
  SK = 0.1, SC = 0.1, SG = 0.8,
  IK = 0, IC = 0, IG = 0,
  RK = 0, RC = 0, RG = 0
)
nextgenR0(Istates = istates, Flist = flist, Vlist = vlist, parameters = params, dfe = df)
