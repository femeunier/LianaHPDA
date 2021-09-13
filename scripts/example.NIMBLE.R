
library(nimble)
library(igraph)
library(tidyverse)

Nobs <- 100

pumpConsts <- list(N = Nobs,
                   t = runif(Nobs,0,300))

alpha_true <- 0.8
beta_true <- 1.2

theta_true <- rgamma(Nobs,alpha_true,beta_true)
lambda_true <- theta_true*pumpConsts$t
xtrue <- rpois(Nobs,lambda_true)

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
pumpData <- list(x = xtrue)

pumpCode <- nimbleCode(
  {
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))

pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

# pump$getNodeNames()
# pump$plotGraph()

set.seed(0)
pump$simulate("theta")
pump$theta
pump$calculate()

# Cpump <- compileNimble(pump)
# Cpump$theta
# Cpump$x
# Cpump$calculate()
# Cpump$simulate()

mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                       data = pumpData, inits = pumpInits,
                       monitors=c("alpha","beta","theta","lambda","x"),
                       nchains = 4, niter = 25000,
                       summary = TRUE, WAIC = TRUE,samplesAsCodaMCMC = TRUE)

Samples <- mcmc.out$samples
plot(Samples[,2])

mcmc.out$summary

mcmc.out$WAI


df <- data.frame(mcmc.out$samples$chain1)
df_l <- df %>% select(alpha, beta, theta.1., lambda.1.) %>% gather(key="parameter", value="value")

p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 60)
p + facet_wrap(~parameter, scales = "free") + theme_bw()
