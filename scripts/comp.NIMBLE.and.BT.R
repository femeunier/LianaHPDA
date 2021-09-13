rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)

e1_approx <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))

    A <- log((0.56146 / x + 0.65) * (1 + x))
    B <- x^4 * exp(7.7 * x) * (2 + x)^3.7

    return((A^-7.7 + B)^-0.13)
  })

NIMprospect5 <- nimbleFunction(
  run = function(N = double(0),Cab = double(0),Car = double(0), Cw = double(0), Cm = double(0),
                 dataspec_p5 = double(2),talf = double(1),t12 = double(1),t21 = double(1), Nwl = double(0)) {

    cc <- matrix(NA,nrow = 5,ncol = 1)
    k <- numeric(length = Nwl)

    Cbrown <- 0

    cc[1,1] <- Cab / N
    cc[2,1] <- Car / N
    cc[3,1] <- Cbrown / N
    cc[4,1] <- Cw / N
    cc[5,1] <- Cm / N

    k[] <- dataspec_p5[,] %*% cc[,]

    trans <- (1 - k)*exp(-k) + k^2 *e1_approx(k)
    trans[trans < 0] <- 0

    ralf <- 1 - talf
    r12 <- 1 - t12
    r21 <- 1 - t21

    denom <- 1 - (r21 ^ 2) * (trans ^ 2)
    Ta <- talf * trans * t21 / denom
    Ra <- ralf + r21 * trans * Ta

    tt <- t12 * trans * t21 / denom
    rr <- r12 + r21 * trans * tt

    gt1 <- rr + tt >= 1
    tgt1 <- tt[gt1]

    Tsub <- 0*tt
    Rsub <- 0*tt
    r <- 0*tt
    t <- 0*tt

    Tsub[gt1] <- tgt1 / (tgt1 + (1 - tgt1) * (N - 1))
    Rsub[gt1] <- 1 - Tsub[gt1]

    inf <- rr == 0 | tt == 0
    Tsub[inf] <- 0
    Rsub[inf] <- 0

    r <- rr[!gt1 & !inf]
    t <- tt[!gt1 & !inf]

    D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    r2 <- r ^ 2
    t2 <- t ^ 2
    va <- (1 + r2 - t2 + D) / (2 * r)
    vb <- (1 - r2 + t2 + D) / (2 * t)

    vbNN <- vb ^ (N - 1)
    vbNN2 <- vbNN ^ 2
    va2 <- va ^ 2
    denomx <- va2 * vbNN2 - 1
    Rsub[!gt1 & !inf] <- va * (vbNN2 - 1) / denomx
    Tsub[!gt1 & !inf] <- vbNN * (va2 - 1) / denomx

    denomy <- 1 - Rsub * rr

    reflectance <- Ra + Ta * Rsub * tt / denomy

    returnType(double(1))
    # print(c(N,Cab, Car, Cbrown, Cw, Cm,mean(reflectance)))
    return(reflectance)
  })

run_prospect5 <- nimbleCode({

  reflectance[1:Nwl] <- NIMprospect5(Nmean,Cabmean,Carmean,Cwmean,Cmmean,
                                     dataspec_p5[,], talf[],t12[],t21[], Nwl)

  for (i in 1:Nsamples){

    creflectance[,i] <- reflectance[1:Nwl]

    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(creflectance[j,i], sd = Standard.Dev)
    }
  }

  Nmean ~ dunif(1.,5)
  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)
  Standard.Dev ~ dunif(0,1)

})


Nleaves <- 10
WLa <- 400
WLb <- 2500
Delta_WL <- 20
WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)


Constants <- list(Nsamples = Nleaves,
                  Nwl = Nwl,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = rrtm:::dataspec_p5[pos,1:5])


Ntrue <- 2 ; Cabtrue <- 40 ; Cartrue <- 20 ; Cwtrue <- 0.01 ; Cmtrue <- 0.01

fac <- 4
Nall <- pmax(1.1,rnorm(Nleaves,Ntrue,Ntrue/fac))
Caball <- pmax(10,rnorm(Nleaves,Cabtrue,Cabtrue/fac))
Carall <- pmax(5,rnorm(Nleaves,Cartrue,Cartrue/fac))
Cwall <- pmax(0,rnorm(Nleaves,Cwtrue,Cwtrue/fac))
Cmall <- pmax(0,rnorm(Nleaves,Cmtrue,Cmtrue/fac))

Data <- list(obs_reflectance = matrix(unlist(lapply(1:length(Nall),function(ileaf){
  rrtm::prospect5(N = Nall[ileaf],
                  Cab = Caball[ileaf],
                  Car = Carall[ileaf],
                  Cbrown = 0,
                  Cw = Cwall[ileaf],
                  Cm = Cmall[ileaf])[["reflectance"]][pos]})),ncol = Nleaves))

matplot(WLs,Data$obs_reflectance,type = 'l')

Inits <- list(Nmean = Ntrue,
              Cabmean = Cabtrue,
              Carmean = Cartrue,
              Cwmean = Cwtrue,
              Cmmean = Cmtrue)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl),
                                         creflectance = c(Nwl,Nleaves)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

# compiled_P5model <- compileNimble(P5model,
#                                   showCompilerOutput = TRUE)

Nchains = 2
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("Nmean","Cabmean","Carmean","Cwmean","Cmmean","Standard.Dev","logProb_obs_reflectance"),
                       data = Data,
                       inits = Inits,
                       nburnin = 10000,
                       nchains = Nchains,
                       niter = 100000,
                       summary = TRUE,
                       WAIC = FALSE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples

param <- MCMCsamples[,1:6]
# plot(param[,6])
# pairs(as.matrix(param[,c(1:5)]), pch = '.')

hist(param[[1]][,6])

Nsimu <- 1000

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:6]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:6]))
}

Simu <- matrix(unlist(lapply(1:Nsimu,function(isimu){
  rrtm::prospect5(N = param_all[isimu,5],
                  Cab = param_all[isimu,1],
                  Car = param_all[isimu,2],
                  Cbrown = 0,
                  Cw = param_all[isimu,4],
                  Cm = param_all[isimu,3])[["reflectance"]][pos]})),ncol = Nsimu)

matplot(WLs,Simu,type = 'l',col = "darkgrey",ylim = c(0,max(c(Data$obs_reflectance,Simu))*1.1))
matlines(WLs,rowMeans(Data$obs_reflectance),col = "red",lwd = 2)
matlines(WLs,Data$obs_reflectance,col = "red")


df_l <- data.frame(mcmc.out$samples$chain1) %>% rename(N = Nmean,
                                                       Cab = Cabmean,
                                                       Car = Carmean,
                                                       Cw = Cwmean,
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev"),
               names_to = "param",
               values_to = "value")

params_true <- data.frame(N = Ntrue,
                          Cab = Cabtrue,
                          Car = Cartrue,
                          Cw = Cwtrue,
                          Cm = Cmtrue) %>% pivot_longer(cols = c("N","Cab","Car","Cw","Cm"),
                                                        names_to = "param",
                                                        values_to = "value")

params_all <- data.frame(N = Nall,
                         Cab = Caball,
                          Car = Carall,
                          Cw = Cwall,
                          Cm = Cmall) %>% pivot_longer(cols = c("N","Cab","Car","Cw","Cm"),
                                                       names_to = "param",
                                                       values_to = "value")



ggplot(df_l,aes(value)) +
  geom_histogram(aes(y= ..density..),bins = 60, alpha = 0.4, color = "darkgrey") +
  geom_point(data = params_true,
             aes(x = value,y = 0), color = 'red') +
  geom_point(data = params_all,
             aes(x = value,y = 0), color = 'red',shape = "|",size = 4) +
  facet_wrap(~param, scales = "free") +
  theme_bw()

######################################################################################
# Compare with Bayesian tools

library(BayesianTools)
library(rrtm)
library(Rprospect)
library(pracma)

run_prospect <- function(params,waves){

  # PROSPECT parameters
  Nlayers <- params[1]
  Cab <- params[2]
  Car <- params[3]
  Cw <- params[4]
  Cm <- params[5]


  # Call RTM
  result <- tryCatch(
    Rprospect::prospect5(Nlayers,Cab,Car,Cw,Cm),
    error = function(e) NULL)
  if (is.null(result)) return(-1e20)
  reflectance <- result[["Reflectance"]]
  if (any(!is.finite(reflectance))) return(-1e20)
  if (any(reflectance < 0)) return(-1e20)
  if (any(reflectance > 1)) return(-1e20)

  reflectance_waves <- interp1(x = result[["Wavelength"]],
                               y = reflectance,
                               waves)
  return(reflectance_waves)
}

create_likelihood <- function(observed, waves) {
  function(params) {

    ssigma <- params[6]
    reflectance_waves <- run_prospect(params,waves)

    # Calculate likelihood
    ll <- 0

    for (i in seq(1,ncol(observed))){
      ll <- ll + sum(dnorm(reflectance_waves, observed[,i], ssigma, log = TRUE))
    }

    return(ll)
  }
}

Prospect_param_names <- c("N","Cab","Car","Cw","Cm","Standard.Dev")
pft_lowers <- c(N = 1, Cab = 0 , Car = 0,Cw = 0, Cm = 0, Standard.Dev = 0)
pft_uppers <-  c(N = 5, Cab = 100, Car = 50,Cw = 0.1, Cm = 0.1, Standard.Dev = 1)
params <- (pft_lowers + pft_uppers)/2

prior <- createUniformPrior(pft_lowers, pft_uppers)
likelihood <- create_likelihood(observed=(Data$obs_reflectance),
                                waves = WLs)
settings_MCMC <- list(iterations = 10000, nrChains = Nchains)

# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)

pos.simu_BT <- sample(1:nrow(getSample(samples,start = 1000)),Nsimu)
param_all_BT <- getSample(samples) [pos.simu_BT,1:6]

colnames(param_all_BT) <- Prospect_param_names

Simu_BT <- matrix(unlist(lapply(1:Nsimu,function(isimu){
  rrtm::prospect5(N = param_all_BT[isimu,1],
                  Cab = param_all_BT[isimu,2],
                  Car = param_all_BT[isimu,3],
                  Cbrown = 0,
                  Cw = param_all_BT[isimu,4],
                  Cm = param_all_BT[isimu,5])[["reflectance"]][pos]})),ncol = Nsimu)

matplot(WLs,Simu_BT,type = 'l',col = "darkgrey",ylim = c(0,max(c(Data$obs_reflectance,Simu))*1.1))
matlines(WLs,rowMeans(Data$obs_reflectance),col = "red",lwd = 2)
matlines(WLs,Data$obs_reflectance,col = "red")


df_all <- bind_rows(list(df_l %>% mutate(type = "NIMBLE"),
                         as.data.frame(param_all_BT) %>% pivot_longer(cols = c(Prospect_param_names),
                                                       names_to = "param",
                                                       values_to = "value") %>% mutate(type = "Bayes.Tools")))

ggplot(df_all,aes(value)) +
  geom_histogram(aes(y= ..density.., fill = type),bins = 60, alpha = 0.4) +
  geom_point(data = params_true,
             aes(x = value,y = 0), color = 'red') +
  geom_point(data = params_all,
             aes(x = value,y = 0), color = 'red',shape = "|",size = 4) +
  facet_wrap(~param, scales = "free") +
  theme_bw()




