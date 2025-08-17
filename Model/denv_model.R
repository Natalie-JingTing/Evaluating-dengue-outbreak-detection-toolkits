
##################################################################################
# MOSQUITO MODEL

# --------------------------
# Parameters

kappa <- user(0.6, min = 0)   # biting rate
omega <- user(0.1, min = 0)   # larval mosquito mortality rate
epsilon <- user(0.1, min = 0) # adult mosquito mortality rate
alpha <- user(1/19, min = 0)  # rate at which mosquitoes leave larval stage
b <- user(0.396775, min = 0)  # rate at which adult female mosquitoes produce larvae
eta <- user(1/10, min = 0)    # rate at which mosquitoes leave the incubation period

# carrying capacity model
K0 <- user(1e6, min = 0)
delta_k <- user(0.3)

M <- A + H1_1 + H2_1 + H3_1 + H4_1 + 
  H1_2 + H2_2 + H3_2 + H4_2 + 
  H1_3 + H2_3 + H3_3 + H4_3 + 
  H1_4 + H2_4 + H3_4 + H4_4 + 
  Y1 + Y2 + Y3 + Y4

Y <- Y1 + Y2 + Y3 + Y4

# --------------------------
# Initialisation

initial(L) <- 0 # larvae
initial(A) <- 0 # uninfected adult female mozzies
initial(H1_1) <- 0 #incubation, serotype 1
initial(H2_1) <- 0 
initial(H3_1) <- 0 
initial(H4_1) <- 0 

initial(H1_2) <- 0 #incubation, serotype 2
initial(H2_2) <- 0 
initial(H3_2) <- 0 
initial(H4_2) <- 0 

initial(H1_3) <- 0 #incubation, serotype 3
initial(H2_3) <- 0 
initial(H3_3) <- 0 
initial(H4_3) <- 0 

initial(H1_4) <- 0 #incubation, serotype 4
initial(H2_4) <- 0 
initial(H3_4) <- 0 
initial(H4_4) <- 0 

initial(Y1) <- 1e3 #infected female mozzies
initial(Y2) <- 1e3 
initial(Y3) <- 1e3
initial(Y4) <- 1e3

# --------------------------
# State updates

# forced seasonal carrying capacity
K <- K0*(1 + delta_k*sin(2*pi_val*(step / 365)))

update(L) <- L + b*M - alpha*L - omega*L*(1 + L / K)
update(A) <- A + alpha*L - (psi1 + psi2 + psi3 + psi4)*A - epsilon*A

update(H1_1) <- H1_1 + psi1*A - (4*eta+epsilon)*H1_1 
update(H2_1) <- H2_1 + 4*eta*H1_1 - (4*eta+epsilon)*H2_1
update(H3_1) <- H3_1 + 4*eta*H2_1 - (4*eta+epsilon)*H3_1
update(H4_1) <- H4_1 + 4*eta*H3_1 - (4*eta+epsilon)*H4_1

update(H1_2) <- H1_2 + psi2*A - (4*eta+epsilon)*H1_2
update(H2_2) <- H2_2 + 4*eta*H1_2 - (4*eta+epsilon)*H2_2
update(H3_2) <- H3_2 + 4*eta*H2_2 - (4*eta+epsilon)*H3_2
update(H4_2) <- H4_2 + 4*eta*H3_2 - (4*eta+epsilon)*H4_2

update(H1_3) <- H1_3 + psi3*A - (4*eta+epsilon)*H1_3
update(H2_3) <- H2_3 + 4*eta*H1_3 - (4*eta+epsilon)*H2_3
update(H3_3) <- H3_3 + 4*eta*H2_3 - (4*eta+epsilon)*H3_3
update(H4_3) <- H4_3 + 4*eta*H3_3 - (4*eta+epsilon)*H4_3

update(H1_4) <- H1_4 + psi4*A - (4*eta+epsilon)*H1_4 
update(H2_4) <- H2_4 + 4*eta*H1_4 - (4*eta+epsilon)*H2_4
update(H3_4) <- H3_4 + 4*eta*H2_4 - (4*eta+epsilon)*H3_4
update(H4_4) <- H4_4 + 4*eta*H3_4 - (4*eta+epsilon)*H4_4

update(Y1) <- Y1 + 4*eta*H4_1 - epsilon*Y1
update(Y2) <- Y2 + 4*eta*H4_2 - epsilon*Y2
update(Y3) <- Y3 + 4*eta*H4_3 - epsilon*Y3
update(Y4) <- Y4 + 4*eta*H4_4 - epsilon*Y4

##################################################################################
# HUMAN MODEL

# --------------------------
# Parameters

N <- user(1e6, min = 0)         # total human population size
p_death <- user(0.0136 / 365, min = 0) # death rate
p_birth <- user(0.0136 / 365, min = 0) # birth rate
p_recov <- user(1 / 365, min = 0) # one over duration of cross-protection
extInf <- user(2e-8)            # additional foi on humans due to external sources
beta_hm <- user(1, min = 0)     # per bite transmission probability from humans to mosquitoes
incub <- user(6, min = 0)       # duration of intrinsic incubation period (time before a human becomes infectious)
inf_period <- user(4, min = 0)  # duration of infectious period in humans
beta_mh_1 <- user(1, min = 0)   # per bite transmission probability from mosquitoes to humans (serotype 1)
beta_mh_2 <- user(1, min = 0)   # per bite transmission probability from mosquitoes to humans (serotype 2)
beta_mh_3 <- user(1, min = 0)   # per bite transmission probability from mosquitoes to humans (serotype 3)
beta_mh_4 <- user(1, min = 0)   # per bite transmission probability from mosquitoes to humans (serotype 4)

N_human <- S0 + R1 + R2 + R3 + R4 + S1 + S2 + S3 + S4 +
  R12 + R13 + R14 + R23 + R24 + R34 +
  S12 + S13 + S14 + S23 + S24 + S34 + 
  R123 + R124 + R134 + R234 +
  S123 + S124 + S134 + S234 + 
  R1234

# old code showing how to set transmission parameter from mosquito to human
# based on desired R0 for each serotype
#beta_mh_1 <- R0_1*epsilon*(eta + epsilon)*N_human / (kappa^2*beta_hm*inf_period*eta*M)
#beta_mh_2 <- R0_2*epsilon*(eta + epsilon)*N_human / (kappa^2*beta_hm*inf_period*eta*M)
#beta_mh_3 <- R0_3*epsilon*(eta + epsilon)*N_human / (kappa^2*beta_hm*inf_period*eta*M)
#beta_mh_4 <- R0_4*epsilon*(eta + epsilon)*N_human / (kappa^2*beta_hm*inf_period*eta*M)

# Lambda_i: force of infection on humans due to serotype i
foi1 <- kappa*beta_mh_1*Y1 / N_human + extInf 
foi2 <- kappa*beta_mh_2*Y2 / N_human + extInf 
foi3 <- kappa*beta_mh_3*Y3 / N_human + extInf 
foi4 <- kappa*beta_mh_4*Y4 / N_human + extInf 

# R0 of each serotype, including both human-to-mosquito and mosquito-to-human components
R0_1 <- kappa^2*beta_mh_1*beta_hm*inf_period*eta*M / (epsilon*(eta + epsilon)*N_human)
R0_2 <- kappa^2*beta_mh_2*beta_hm*inf_period*eta*M / (epsilon*(eta + epsilon)*N_human)
R0_3 <- kappa^2*beta_mh_3*beta_hm*inf_period*eta*M / (epsilon*(eta + epsilon)*N_human)
R0_4 <- kappa^2*beta_mh_4*beta_hm*inf_period*eta*M / (epsilon*(eta + epsilon)*N_human)

# R0 corrected for math issues
R0_1b <- kappa^2 * beta_mh_1 * beta_hm * (4*eta / (epsilon + 4*eta))^4 * (1 / epsilon) * (M / N_human) *
  (1 - p_death)^incub * 1 / (1 / inf_period - log(1 - p_death))
R0_2b <- kappa^2 * beta_mh_2 * beta_hm * (4*eta / (epsilon + 4*eta))^4 * (1 / epsilon) * (M / N_human) *
  (1 - p_death)^incub * 1 / (1 / inf_period - log(1 - p_death))
R0_3b <- kappa^2 * beta_mh_3 * beta_hm * (4*eta / (epsilon + 4*eta))^4 * (1 / epsilon) * (M / N_human) *
  (1 - p_death)^incub * 1 / (1 / inf_period - log(1 - p_death))
R0_4b <- kappa^2 * beta_mh_4 * beta_hm * (4*eta / (epsilon + 4*eta))^4 * (1 / epsilon) * (M / N_human) *
  (1 - p_death)^incub * 1 / (1 / inf_period - log(1 - p_death))

# --------------------------
# Initialisation

initial(S0) <- N
initial(S1) <- 0
initial(S2) <- 0
initial(S3) <- 0
initial(S4) <- 0

initial(S12) <- 0
initial(S13) <- 0
initial(S14) <- 0
initial(S23) <- 0
initial(S24) <- 0
initial(S34) <- 0

initial(S123) <- 0
initial(S124) <- 0
initial(S234) <- 0
initial(S134) <- 0

initial(R1) <- 0
initial(R2) <- 0
initial(R3) <- 0
initial(R4) <- 0

initial(R12) <- 0
initial(R13) <- 0
initial(R14) <- 0
initial(R23) <- 0
initial(R24) <- 0
initial(R34) <- 0

initial(R123) <- 0
initial(R124) <- 0
initial(R234) <- 0
initial(R134) <- 0

initial(R1234) <- 0

initial(infectious1) <- 0
initial(exposed1) <- 0
initial(infectious2) <- 0
initial(exposed2) <- 0
initial(infectious3) <- 0
initial(exposed3) <- 0
initial(infectious4) <- 0
initial(exposed4) <- 0

# --------------------------
# State updates

update(S0) <- S0 + n_births - nS0_R1 - nS0_R2 - nS0_R3 - nS0_R4 - nS0_deaths

update(S1) <- S1 + nR1_S1 - nS1_R12 - nS1_R13 - nS1_R14 - nS1_deaths
update(S2) <- S2 + nR2_S2 - nS2_R12 - nS2_R23 - nS2_R24 - nS2_deaths
update(S3) <- S3 + nR3_S3 - nS3_R13 - nS3_R23 - nS3_R34 - nS3_deaths
update(S4) <- S4 + nR4_S4 - nS4_R14 - nS4_R24 - nS4_R34 - nS4_deaths

update(S12) <- S12 + nR12_S12 - nS12_R123 - nS12_R124 - nS12_deaths
update(S13) <- S13 + nR13_S13 - nS13_R123 - nS13_R134 - nS13_deaths
update(S14) <- S14 + nR14_S14 - nS14_R124 - nS14_R134 - nS14_deaths
update(S23) <- S23 + nR23_S23 - nS23_R123 - nS23_R234 - nS23_deaths
update(S24) <- S24 + nR24_S24 - nS24_R124 - nS24_R234 - nS24_deaths
update(S34) <- S34 + nR34_S34 - nS34_R134 - nS34_R234 - nS34_deaths

update(S123) <- S123 + nR123_S123 - nS123_R1234 - nS123_deaths
update(S124) <- S124 + nR124_S124 - nS124_R1234 - nS124_deaths
update(S134) <- S134 + nR134_S134 - nS134_R1234 - nS134_deaths
update(S234) <- S234 + nR234_S234 - nS234_R1234 - nS234_deaths

update(R1) <- R1 + nS0_R1 - nR1_S1 - nR1_deaths
update(R2) <- R2 + nS0_R2 - nR2_S2 - nR2_deaths
update(R3) <- R3 + nS0_R3 - nR3_S3 - nR3_deaths
update(R4) <- R4 + nS0_R4 - nR4_S4 - nR4_deaths

update(R12) <- R12 + nS1_R12 + nS2_R12 - nR12_S12 - nR12_deaths
update(R13) <- R13 + nS1_R13 + nS3_R13 - nR13_S13 - nR13_deaths
update(R14) <- R14 + nS1_R14 + nS4_R14 - nR14_S14 - nR14_deaths
update(R23) <- R23 + nS2_R23 + nS3_R23 - nR23_S23 - nR23_deaths
update(R24) <- R24 + nS2_R24 + nS4_R24 - nR24_S24 - nR24_deaths
update(R34) <- R34 + nS3_R34 + nS4_R34 - nR34_S34 - nR34_deaths

update(R123) <- R123 + nS12_R123 + nS13_R123 + nS23_R123 - nR123_S123 - nR123_deaths
update(R124) <- R124 + nS12_R124 + nS14_R124 + nS24_R124 - nR124_S124 - nR124_deaths
update(R134) <- R134 + nS13_R134 + nS14_R134 + nS34_R134 - nR134_S134 - nR134_deaths
update(R234) <- R234 + nS23_R234 + nS24_R234 + nS34_R234 - nR234_S234 - nR234_deaths

update(R1234) <- R1234 + nS123_R1234 + nS124_R1234 + nS134_R1234 + nS234_R1234 - nR1234_deaths

# --------------------------
# Manual integration to produce Psi (force of infection on mosquitoes)

# calculate the total incidence of infection of each serotype, summed over all states
inc_inf1 <- nS0_R1 + 
  nS2_R12 + nS3_R13 + nS4_R14 + 
  nS23_R123 + nS24_R124 + nS34_R134 + 
  nS234_R1234
inc_inf2 <- nS0_R2 + 
  nS1_R12 + nS3_R23 + nS4_R24 +
  nS13_R123 + nS14_R124 + nS34_R234 + 
  nS134_R1234
inc_inf3 <- nS0_R3 + 
  nS1_R13 + nS2_R23 + nS4_R34 +
  nS12_R123 + nS14_R134 + nS24_R234 + 
  nS124_R1234
inc_inf4 <- nS0_R4 + 
  nS1_R14 + nS2_R24 + nS3_R34 + 
  nS12_R124 + nS13_R134 + nS23_R234 +
  nS123_R1234

update(exposed1) <- exposed1 + inc_inf1 - exposed1/incub
update(infectious1) <- infectious1 - infectious1/inf_period + exposed1/incub

update(exposed2) <- exposed2 + inc_inf2 - exposed2/incub
update(infectious2) <- infectious2 - infectious2/inf_period + exposed2/incub

update(exposed3) <- exposed3 + inc_inf3 - exposed3/incub
update(infectious3) <- infectious3 - infectious3/inf_period + exposed3/incub

update(exposed4) <- exposed4 + inc_inf4 - exposed4/incub
update(infectious4) <- infectious4 - infectious4/inf_period + exposed4/incub

psi1 <- beta_hm*kappa*infectious1 / N_human
psi2 <- beta_hm*kappa*infectious2 / N_human
psi3 <- beta_hm*kappa*infectious3 / N_human
psi4 <- beta_hm*kappa*infectious4 / N_human

# Option for manually setting foi on mosquitoes
# psi1 <- user(0.02, min=0)
# psi2 <- user(0.02, min=0)
# psi3 <- user(0.02, min=0)
# psi4 <- user(0.02, min=0)

# ================================================================================
# UPDATE ELEMENTS

# --------------------------
# SUSCEPTIBLE

# --------------------------
# Zero past infections

## S0
n_births <- rbinom(N_human, p_birth)

sum_pS0 <- foi1 + foi2 + foi3 + foi4 + p_death
prS0[1] <- foi1 / sum_pS0 
prS0[2] <- foi2 / sum_pS0
prS0[3] <- foi3 / sum_pS0 
prS0[4] <- foi4 / sum_pS0
prS0[5] <- p_death / sum_pS0

pS0 <- 1 - exp(-(sum_pS0))
nS0 <- rbinom(S0, pS0)

nS0_event[] <- rmultinom(nS0, prS0)

nS0_R1 <- nS0_event[1]
nS0_R2 <- nS0_event[2] 
nS0_R3 <- nS0_event[3]
nS0_R4 <- nS0_event[4]
nS0_deaths <- nS0_event[5]

# --------------------------
# One past infection

## S1
sum_pS1 <- foi2 + foi3 + foi4 + p_death
pS1 <- 1 - exp(-(sum_pS1))

nS1 <- rbinom(S1, pS1)
nS1_event[] <- rmultinom(nS1, prS1)
prS1[1] <- foi2 / sum_pS1
prS1[2] <- foi3 / sum_pS1
prS1[3] <- foi4 / sum_pS1
prS1[4] <- p_death / sum_pS1

nS1_R12 <- nS1_event[1]
nS1_R13 <- nS1_event[2]
nS1_R14 <- nS1_event[3]
nS1_deaths <- nS1_event[4]

## S2
sum_pS2 <- foi1 + foi3 + foi4 + p_death
pS2 <- 1 - exp(-(sum_pS2))

nS2 <- rbinom(S2, pS2)
nS2_event[] <- rmultinom(nS2, prS2)
prS2[1] <- foi1 / sum_pS2
prS2[2] <- foi3 / sum_pS2
prS2[3] <- foi4 / sum_pS2
prS2[4] <- p_death / sum_pS2

nS2_R12 <- nS2_event[1]
nS2_R23 <- nS2_event[2]
nS2_R24 <- nS2_event[3]
nS2_deaths <- nS2_event[4]

## S3
sum_pS3 <- foi1 + foi2 + foi4 + p_death
pS3 <- 1 - exp(-(sum_pS3))

nS3 <- rbinom(S3, pS3)
nS3_event[] <- rmultinom(nS3, prS3)
prS3[1] <- foi1 / sum_pS3
prS3[2] <- foi2 / sum_pS3
prS3[3] <- foi4 / sum_pS3
prS3[4] <- p_death / sum_pS3

nS3_R13 <- nS3_event[1]
nS3_R23 <- nS3_event[2]
nS3_R34 <- nS3_event[3]
nS3_deaths <- nS3_event[4] 

## S4
sum_pS4 <- foi1 + foi2 + foi3 + p_death
pS4 <- 1 - exp(-(sum_pS4))

nS4 <- rbinom(S4, pS4)
nS4_event[] <- rmultinom(nS4, prS4)
prS4[1] <- foi1 / sum_pS4
prS4[2] <- foi2 / sum_pS4
prS4[3] <- foi3 / sum_pS4
prS4[4] <- p_death / sum_pS4

nS4_R14 <- nS4_event[1]
nS4_R24 <- nS4_event[2]
nS4_R34 <- nS4_event[3]
nS4_deaths <- nS4_event[4]

# --------------------------
# Two past infections

## S12
sum_pS12 <- foi3 + foi4 + p_death
pS12 <- 1 - exp(-(sum_pS12))

nS12 <- rbinom(S12, pS12)
nS12_event[] <- rmultinom(nS12, prS12)
prS12[1] <- foi3 / sum_pS12
prS12[2] <- foi4 / sum_pS12
prS12[3] <- p_death / sum_pS12

nS12_R123 <- nS12_event[1]
nS12_R124 <- nS12_event[2]
nS12_deaths <- nS12_event[3]

## S13
sum_pS13 <- foi2 + foi4 + p_death
pS13 <- 1 - exp(-(sum_pS13))

nS13 <- rbinom(S13, pS13)
nS13_event[] <- rmultinom(nS13, prS13)
prS13[1] <- foi2 / sum_pS13
prS13[2] <- foi4 / sum_pS13
prS13[3] <- p_death / sum_pS13

nS13_R123 <- nS13_event[1]
nS13_R134 <- nS13_event[2]
nS13_deaths <- nS13_event[3]

## S14
sum_pS14 <- foi2 + foi3 + p_death
pS14 <- 1 - exp(-(sum_pS14))

nS14 <- rbinom(S14, pS14)
nS14_event[] <- rmultinom(nS14, prS14)
prS14[1] <- foi2 / sum_pS14
prS14[2] <- foi3 / sum_pS14
prS14[3] <- p_death / sum_pS14

nS14_R124 <- nS14_event[1]
nS14_R134 <- nS14_event[2]
nS14_deaths <- nS14_event[3]

## S23
sum_pS23 <- foi1 + foi4 + p_death
pS23 <- 1 - exp(-(sum_pS23))

nS23 <- rbinom(S23, pS23)
nS23_event[] <- rmultinom(nS23, prS23)
prS23[1] <- foi1 / sum_pS23
prS23[2] <- foi4 / sum_pS23
prS23[3] <- p_death / sum_pS23

nS23_R123 <- nS23_event[1]
nS23_R234 <- nS23_event[2]
nS23_deaths <- nS23_event[3]

## S24
sum_pS24 <- foi1 + foi3 + p_death
pS24 <- 1 - exp(-(sum_pS24))

nS24 <- rbinom(S24, pS24)
nS24_event[] <- rmultinom(nS24, prS24)
prS24[1] <- foi1 / sum_pS24
prS24[2] <- foi3 / sum_pS24
prS24[3] <- p_death / sum_pS24

nS24_R124 <- nS24_event[1]
nS24_R234 <- nS24_event[2]
nS24_deaths <- nS24_event[3]

## S34
sum_pS34 <- foi1 + foi2 + p_death
pS34 <- 1 - exp(-(sum_pS34))

nS34 <- rbinom(S34, pS34)
nS34_event[] <- rmultinom(nS34, prS34)
prS34[1] <- foi1 / sum_pS34
prS34[2] <- foi2 / sum_pS34
prS34[3] <- p_death / sum_pS34

nS34_R134 <- nS34_event[1]
nS34_R234 <- nS34_event[2]
nS34_deaths <- nS34_event[3]

# --------------------------
# Three past infections

## S123
sum_pS123 <- foi4 + p_death
pS123 <- 1 - exp(-(sum_pS123))

nS123 <- rbinom(S123, pS123)
nS123_event[] <- rmultinom(nS123, prS123)
prS123[1] <- foi4 / sum_pS123
prS123[2] <- p_death / sum_pS123

nS123_R1234 <- nS123_event[1]
nS123_deaths <- nS123_event[2]

## S124
sum_pS124 <- foi3 + p_death
pS124 <- 1 - exp(-(sum_pS124))

nS124 <- rbinom(S124, pS124)
nS124_event[] <- rmultinom(nS124, prS124)
prS124[1] <- foi3 / sum_pS124
prS124[2] <- p_death / sum_pS124

nS124_R1234 <- nS124_event[1]
nS124_deaths <- nS124_event[2]

## S234
sum_pS234 <- foi1 + p_death
pS234 <- 1 - exp(-(sum_pS234))

nS234 <- rbinom(S234, pS234)
nS234_event[] <- rmultinom(nS234, prS234)
prS234[1] <- foi1 / sum_pS234
prS234[2] <- p_death / sum_pS234

nS234_R1234 <- nS234_event[1]
nS234_deaths <- nS234_event[2]

## S134
sum_pS134 <- foi2 + p_death
pS134 <- 1 - exp(-(sum_pS134))

nS134 <- rbinom(S134, pS134)
nS134_event[] <- rmultinom(nS134, prS134)
prS134[1] <- foi2 / sum_pS134
prS134[2] <- p_death / sum_pS134

nS134_R1234 <- nS134_event[1]
nS134_deaths <- nS134_event[2]

# --------------------------
# RECOVERED

# --------------------------
# One past infection

## R1
sum_pR1 <- p_recov + p_death
pR1 <- 1 - exp(-(sum_pR1))

nR1 <- rbinom(R1, pR1)
nR1_event[] <- rmultinom(nR1, prR1)
prR1[1] <- p_recov / sum_pR1
prR1[2] <- p_death / sum_pR1

nR1_S1 <- nR1_event[1]
nR1_deaths <- nR1_event[2]

## R2
sum_pR2 <- p_recov + p_death
pR2 <- 1 - exp(-(sum_pR2))

nR2 <- rbinom(R2, pR2)
nR2_event[] <- rmultinom(nR2, prR2)
prR2[1] <- p_recov / sum_pR2
prR2[2] <- p_death / sum_pR2

nR2_S2 <- nR2_event[1]
nR2_deaths <- nR2_event[2]

## R3
sum_pR3 <- p_recov + p_death
pR3 <- 1 - exp(-(sum_pR3))

nR3 <- rbinom(R3, pR3)
nR3_event[] <- rmultinom(nR3, prR3)
prR3[1] <- p_recov / sum_pR3
prR3[2] <- p_death / sum_pR3

nR3_S3 <- nR3_event[1]
nR3_deaths <- nR3_event[2]

## R4
sum_pR4 <- p_recov + p_death
pR4 <- 1 - exp(-(sum_pR4))

nR4 <- rbinom(R4, pR4)
nR4_event[] <- rmultinom(nR4, prR4)
prR4[1] <- p_recov / sum_pR4
prR4[2] <- p_death / sum_pR4

nR4_S4 <- nR4_event[1]
nR4_deaths <- nR4_event[2]

# --------------------------
# Two past infections

## R12
sum_pR12 <- p_recov + p_death
pR12 <- 1 - exp(-(sum_pR12))

nR12 <- rbinom(R12, pR12)
nR12_event[] <- rmultinom(nR12, prR12)
prR12[1] <- p_recov / sum_pR12
prR12[2] <- p_death / sum_pR12

nR12_S12 <- nR12_event[1]
nR12_deaths <- nR12_event[2]

## R13
sum_pR13 <- p_recov + p_death
pR13 <- 1 - exp(-(sum_pR13))

nR13 <- rbinom(R13, pR13)
nR13_event[] <- rmultinom(nR13, prR13)
prR13[1] <- p_recov / sum_pR13
prR13[2] <- p_death / sum_pR13

nR13_S13 <- nR13_event[1]
nR13_deaths <- nR13_event[2]

## R14
sum_pR14 <- p_recov + p_death
pR14 <- 1 - exp(-(sum_pR14))

nR14 <- rbinom(R14, pR14)
nR14_event[] <- rmultinom(nR14, prR14)
prR14[1] <- p_recov / sum_pR14
prR14[2] <- p_death / sum_pR14

nR14_S14 <- nR14_event[1]
nR14_deaths <- nR14_event[2]

## R23
sum_pR23 <- p_recov + p_death
pR23 <- 1 - exp(-(sum_pR23))

nR23 <- rbinom(R23, pR23)
nR23_event[] <- rmultinom(nR23, prR23)
prR23[1] <- p_recov / sum_pR23
prR23[2] <- p_death / sum_pR23

nR23_S23 <- nR23_event[1]
nR23_deaths <- nR23_event[2]

## R24
sum_pR24 <- p_recov + p_death
pR24 <- 1 - exp(-(sum_pR24))

nR24 <- rbinom(R24, pR24)
nR24_event[] <- rmultinom(nR24, prR24)
prR24[1] <- p_recov / sum_pR24
prR24[2] <- p_death / sum_pR24

nR24_S24 <- nR24_event[1]
nR24_deaths <- nR24_event[2]

## R34
sum_pR34 <- p_recov + p_death
pR34 <- 1 - exp(-(sum_pR34))

nR34 <- rbinom(R34, pR34)
nR34_event[] <- rmultinom(nR34, prR34)
prR34[1] <- p_recov / sum_pR34
prR34[2] <- p_death / sum_pR34

nR34_S34 <- nR34_event[1]
nR34_deaths <- nR34_event[2]

# --------------------------
# Three past infections

## R123
sum_pR123 <- p_recov + p_death
pR123 <- 1 - exp(-(sum_pR123))

nR123 <- rbinom(R123, pR123)
nR123_event[] <- rmultinom(nR123, prR123)
prR123[1] <- p_recov / sum_pR123
prR123[2] <- p_death / sum_pR123

nR123_S123 <- nR123_event[1]
nR123_deaths <- nR123_event[2]

## R124
sum_pR124 <- p_recov + p_death
pR124 <- 1 - exp(-(sum_pR124))

nR124 <- rbinom(R124, pR124)
nR124_event[] <- rmultinom(nR124, prR124)
prR124[1] <- p_recov / sum_pR124
prR124[2] <- p_death / sum_pR124

nR124_S124 <- nR124_event[1]
nR124_deaths <- nR124_event[2]

## R234
sum_pR234 <- p_recov + p_death
pR234 <- 1 - exp(-(sum_pR234))

nR234 <- rbinom(R234, pR234)
nR234_event[] <- rmultinom(nR234, prR234)
prR234[1] <- p_recov / sum_pR234
prR234[2] <- p_death / sum_pR234

nR234_S234 <- nR234_event[1]
nR234_deaths <- nR234_event[2]

## R134
sum_pR134 <- p_recov + p_death
pR134 <- 1 - exp(-(sum_pR134))

nR134 <- rbinom(R134, pR134)
nR134_event[] <- rmultinom(nR134, prR134)
prR134[1] <- p_recov / sum_pR134
prR134[2] <- p_death / sum_pR134

nR134_S134 <- nR134_event[1]
nR134_deaths <- nR134_event[2]

# --------------------------
# Four past infections (fully saturated)

## R1234
nR1234_deaths <- rbinom(R1234, p_death)


##################################################################################
# INITIALISING OBJECTS

dim(prS0) <- 5

dim(prR1) <- 2
dim(prR2) <- 2
dim(prR3) <- 2
dim(prR4) <- 2

dim(prS1) <- 4
dim(prS2) <- 4
dim(prS3) <- 4
dim(prS4) <- 4

dim(prR12) <- 2
dim(prR13) <- 2
dim(prR14) <- 2
dim(prR23) <- 2
dim(prR24) <- 2
dim(prR34) <- 2

dim(prS12) <- 3
dim(prS13) <- 3
dim(prS14) <- 3
dim(prS23) <- 3
dim(prS24) <- 3
dim(prS34) <- 3

dim(prR123) <- 2
dim(prR124) <- 2
dim(prR134) <- 2
dim(prR234) <- 2

dim(prS123) <- 2
dim(prS124) <- 2
dim(prS134) <- 2
dim(prS234) <- 2

dim(nS0_event) <- 5

dim(nR1_event) <- 2
dim(nR2_event) <- 2
dim(nR3_event) <- 2
dim(nR4_event) <- 2

dim(nS1_event) <- 4
dim(nS2_event) <- 4
dim(nS3_event) <- 4
dim(nS4_event) <- 4

dim(nR12_event) <- 2
dim(nR13_event) <- 2
dim(nR14_event) <- 2
dim(nR23_event) <- 2
dim(nR24_event) <- 2
dim(nR34_event) <- 2

dim(nS12_event) <- 3
dim(nS13_event) <- 3
dim(nS14_event) <- 3
dim(nS23_event) <- 3
dim(nS24_event) <- 3
dim(nS34_event) <- 3

dim(nR123_event) <- 2
dim(nR124_event) <- 2
dim(nR134_event) <- 2
dim(nR234_event) <- 2

dim(nS123_event) <- 2
dim(nS124_event) <- 2
dim(nS134_event) <- 2
dim(nS234_event) <- 2


##################################################################################
# OUTPUTS

# --------------------------
# Mosquito

output(K) <- K
output(Y) <- Y
output(M) <- M

output(psi1) <- psi1
output(psi2) <- psi2 
output(psi3) <- psi3
output(psi4) <- psi4

# --------------------------
# Human

output(N_human) <- N_human

output(foi1) <- foi1
output(foi2) <- foi2
output(foi3) <- foi3
output(foi4) <- foi4

output(inc_inf1) <- inc_inf1
output(inc_inf2) <- inc_inf2
output(inc_inf3) <- inc_inf3
output(inc_inf4) <- inc_inf4

output(rate_inc_inf1) <- kappa * beta_mh_1 / N_human * (4*eta*H4_1 - epsilon*Y1)
output(rate_inc_inf2) <- kappa * beta_mh_2 / N_human * (4*eta*H4_2 - epsilon*Y2)
output(rate_inc_inf3) <- kappa * beta_mh_3 / N_human * (4*eta*H4_3 - epsilon*Y3)
output(rate_inc_inf4) <- kappa * beta_mh_4 / N_human * (4*eta*H4_4 - epsilon*Y4)

output(Rt_1) <- R0_1b * A/M * (S0 + S2 + S3 + S4 + S23 + S24 + S34 + S234)/N
output(Rt_2) <- R0_2b * A/M * (S0 + S1 + S3 + S4 + S13 + S14 + S34 + S134)/N
output(Rt_3) <- R0_3b * A/M * (S0 + S1 + S2 + S4 + S12 + S14 + S24 + S124)/N
output(Rt_4) <- R0_4b * A/M * (S0 + S1 + S2 + S3 + S12 + S13 + S23 + S123)/N

output(oldRt_1) <- R0_1 * A/M * (S0 + S2 + S3 + S4 + S23 + S24 + S34 + S234)/N
output(oldRt_2) <- R0_2 * A/M * (S0 + S1 + S3 + S4 + S13 + S14 + S34 + S134)/N
output(oldRt_3) <- R0_3 * A/M * (S0 + S1 + S2 + S4 + S12 + S14 + S24 + S124)/N
output(oldRt_4) <- R0_4 * A/M * (S0 + S1 + S2 + S3 + S12 + S13 + S23 + S123)/N

##################################################################################
# MISC

pi_val <- 3.141593
