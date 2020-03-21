
library(emplik, quietly = TRUE, warn.conflicts = FALSE)
data(myeloma)

n <- nrow(myeloma)
time <-  myeloma[ , 1]    ## survival or censoring time
vstatus <- myeloma[ , 2]  ##  0 = alive (i.e., censored)
alive <- vstatus == 0
cens_time <- rep(NA, n)
cens_time[alive] <- time[alive]
cens_time[!alive] <- Inf
time[alive] <- NA
## covariates:
logBUN <- myeloma[ , 3]
HGB <- myeloma[ , 4]
logBUN <- (logBUN - mean(logBUN)) / sd(logBUN)
HGB <- (HGB - mean(HGB)) / sd(HGB)

codeAFT <- nimbleCode({
    for(i in 1:n) {
        x[i] ~ dweib(alpha, lambda[i])
        is_cens[i] ~ dinterval(x[i], c[i])  ## right-censoring
        lambda[i] <- exp(eta[i] + Z[i,1]*delta[1] + Z[i,2]*delta[2])
        eta[i] <- etaTilde[xi[i]]  ## mix over eta; mu = exp(eta)
    }
    xi[1:n] ~ dCRP(conc, size = n) ## CRP for mixture components
    conc ~ dgamma(1, 1)
    for(i in 1:nSub)    ## cap the number of clusters for faster computation
        etaTilde[i] ~ dunif(b0, B0) ## base measure H_b
    alpha ~ dunif(a0, A0)
    for(j in 1:p)
        delta[j] ~ dflat()
})

nSub = 15
constants = list(b0 = -10, B0 = 10, a0 = 0.1, A0 = 10, p = 2, n = n,
                 c = cens_time, Z = cbind(logBUN, HGB), nSub = nSub)
data = list(is_cens = as.numeric(alive), x = time)
xInit <- rep(NA, n)
xInit[alive] <- cens_time[alive] + 10
inits = list(alpha = 1, delta = c(0, 0), conc = 1,
             etaTilde = runif(nSub, constants$b0, constants$B0),
             xi = sample(1:3, n, replace = TRUE), x = xInit)
model <- nimbleModel(codeAFT, constants = constants, data = data, inits = inits)
