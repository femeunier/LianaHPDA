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

  for (i in 1:Nsamples){

    N[i] ~ dnorm(Nmean, sd = Nmean/100)
    Cab[i] ~ dnorm(Cabmean, sd = Cabmean/100)
    Car[i] ~ dnorm(Carmean, sd = Carmean/100)
    Cw[i] ~ dnorm(Cwmean, sd = Cwmean/100)
    Cm[i] ~ dnorm(Cmmean, sd = Cmmean/100)

    reflectance[,i] <- NIMprospect5(N[i],Cab[i],Car[i],Cw[i],Cm[i],
                                    dataspec_p5[,], talf[],t12[],t21[], Nwl)


    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(reflectance[j,i], sd = Standard.Dev)
      # obs_reflectance[j,i] ~ dunif(0.99*reflectance[j,i], 1.01*reflectance[j,i])
    }
  }

  Nmean ~ dunif(1,5)
  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
  Cwmean ~ dunif(0.00001,0.1)
  Cmmean ~ dunif(0.00001,0.1)
  Standard.Dev ~ dunif(0,1)


})

Nleaves <- 3
WLa <- 400
WLb <- 2500
Delta_WL <- 5
WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)

Constants <- list(Nsamples = Nleaves,
                  Nwl = Nwl,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = rrtm:::dataspec_p5[pos,1:5])


Data <- list(obs_reflectance = readRDS("./data/All.spectra.RDS") %>% filter(Ind == 1, site == "PNM",GF == "Liana",Species == 6) %>%
  dplyr::select(wv,name,value) %>% pivot_wider(names_from = name,
                                               values_from = value) %>% dplyr::select(-c(wv)) %>% as.matrix())

matplot(WLs,Data$obs_reflectance,type = 'l')

Inits <- list(Nmean = 2,
              Cabmean = 40,
              Carmean = 20,
              Cwmean = 0.01,
              Cmmean = 0.01)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl,Nleaves)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

Nchains = 2
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("Nmean","Cabmean","Carmean","Cwmean","Cmmean","Standard.Dev"),
                       data = Data,
                       inits = Inits,
                       nburnin = 5000,
                       nchains = Nchains,
                       niter = 20000,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples

param <- MCMCsamples[,1:5]
plot(param[,5])
# pairs(as.matrix(param[,c(1:5)]), pch = '.')



Nsimu <- 1000

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:5]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:5]))
}

Simu <- matrix(unlist(lapply(1:Nsimu,function(isimu){
  rrtm::prospect5(N = param_all[isimu,5],
                  Cab = param_all[isimu,1],
                  Car = param_all[isimu,2],
                  Cbrown = 0,
                  Cw = param_all[isimu,4],
                  Cm = param_all[isimu,3])[["reflectance"]][pos]})),ncol = Nsimu)

matplot(WLs,Simu,type = 'l',col = "darkgrey",ylim = c(0,max(c(Data$obs_reflectance,Simu))*1.1))
matlines(WLs,Data$obs_reflectance,col = "red")
lines(WLs,rowMeans(Data$obs_reflectance),col = "red",lwd = 2)

meanSimu <- rowMeans(Simu)
meanData <- rowMeans(Data$obs_reflectance)

hist(meanSimu - meanData)
hist((meanSimu - meanData)/meanData)

plot(meanSimu,(meanSimu - meanData))
abline(h = 0)


df_l <- data.frame(mcmc.out$samples$chain1) %>% rename(N = Nmean,
                                                       Cab = Cabmean,
                                                       Car = Carmean,
                                                       Cw = Cwmean,
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev"),
               names_to = "param",
               values_to = "value")


ggplot(df_l,aes(value)) +
  geom_histogram(aes(y= ..density..),bins = 60, alpha = 0.4, color = "darkgrey") +
  facet_wrap(~param, scales = "free") +
  theme_bw()



