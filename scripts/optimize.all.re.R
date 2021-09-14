rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(reshape2)

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

    if (N < 1.1 | Car < 0 | Cab < 0 | Cm <= 0 | Cw <0) {
      return(-9999*(t12**0))
    }

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

    return(reflectance)
  })

run_prospect5 <- nimbleCode({

  Standard.Dev ~ dunif(0,1)

  Nmean ~ dunif(1.,5)
  Cabmean ~ dunif(0,200)
  Carmean ~ dunif(0,100)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)

  # GF effect
  alpha_N ~ dunif(-1,1)
  alpha_Cab ~ dunif(-100,100)
  alpha_Car ~ dunif(-50,50)
  alpha_Cw ~ dunif(-0.05,0.05)
  alpha_Cm ~ dunif(-0.05,0.05)

  # Site effect
  beta_N ~ dunif(-1,1)
  beta_Cab ~ dunif(-100,100)
  beta_Car ~ dunif(-50,50)
  beta_Cw ~ dunif(-0.05,0.05)
  beta_Cm ~ dunif(-0.05,0.05)

  # Random effects
  sigma_RE_N ~ dunif(-1, 1)
  sigma_RE_Cab ~ dunif(-10, 10)
  sigma_RE_Car ~ dunif(-10, 10)
  sigma_RE_Cm ~ dunif(-0.01, 0.01)
  sigma_RE_Cw ~ dunif(-0.01, 0.01)

  mean_RE_N ~ dunif(0, 1)
  mean_RE_Cab ~ dunif(0, 10)
  mean_RE_Car ~ dunif(0, 10)
  mean_RE_Cm ~ dunif(0, 0.05)
  mean_RE_Cw ~ dunif(0, 0.05)

  for(ispectra in 1:Nspectra){

    gamma_RE_N[ispectra] ~ dnorm(mean_RE_N, sd = sigma_RE_N)
    gamma_RE_Cab[ispectra] ~ dnorm(mean_RE_Cab, sd = sigma_RE_Cab)
    gamma_RE_Car[ispectra] ~ dnorm(mean_RE_Car, sd = sigma_RE_Car)
    gamma_RE_Cm[ispectra] ~ dnorm(mean_RE_Cm, sd = sigma_RE_Cm)
    gamma_RE_Cw[ispectra] ~ dnorm(mean_RE_Cw, sd = sigma_RE_Cw)

    reflectance[1:Nwl,ispectra] <- NIMprospect5(Nmean + is.liana[ispectra] * alpha_N + is.FTS[ispectra] * beta_N + gamma_RE_N[ispectra],
                                                Cabmean + is.liana[ispectra] * alpha_Cab + is.FTS[ispectra] * beta_Cab + gamma_RE_Cab[ispectra],
                                                Carmean + is.liana[ispectra] * alpha_Car + is.FTS[ispectra] * beta_Car + gamma_RE_Car[ispectra],
                                                Cwmean + is.liana[ispectra] * alpha_Cw + is.FTS[ispectra] * beta_Cw + gamma_RE_Cw[ispectra],
                                                Cmmean + is.liana[ispectra] * alpha_Cm + is.FTS[ispectra] * beta_Cm + gamma_RE_Cm[ispectra],
                                                dataspec_p5[,], talf[],t12[],t21[], Nwl)

    for (j in 1:Nwl){
      obs_reflectance[j,ispectra] ~ dnorm(reflectance[j,ispectra], sd = Standard.Dev)
    }
  }
})

WLa <- 400
WLb <- 2500
Delta_WL <- 50

WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)

data.spectra <- readRDS("./data/All.spectra.ID.RDS") %>% group_by(GF,Species,Ind,site) %>% mutate(ind.id = cur_group_id()) %>%
  group_by(GF,Species,Ind,site,name) %>% mutate(leaf.id = cur_group_id()) %>%
  mutate(value2 = value,
         value = case_when(wv %in% c(400:680,800:2500) ~ value,
                           TRUE ~ NA_real_))

data.spectra_T_PNM <- data.spectra %>% filter(GF == "Tree",
                                              site == "PNM") %>%
  filter(wv %in% WLs) %>% arrange(wv) %>% ungroup() %>%
  dplyr::select(value,value2)

data.spectra_T_FTS <- data.spectra %>% filter(GF == "Tree",
                                              site == "FTS") %>%
  filter(wv %in% WLs) %>% arrange(wv) %>% ungroup() %>%
  dplyr::select(value,value2)

data.spectra_L_PNM <- data.spectra %>% filter(GF == "Liana",
                                              site == "PNM") %>%
  filter(wv %in% WLs) %>% arrange(wv) %>% ungroup() %>%
  dplyr::select(value,value2)

data.spectra_L_FTS <- data.spectra %>% filter(GF == "Liana",
                                              site == "FTS") %>%
  filter(wv %in% WLs) %>% arrange(wv) %>% ungroup() %>%
  dplyr::select(value,value2)

N = 5
Nspectra = N*4

N_T_PNM <- ncol((data.spectra_T_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])
N_T_FTS <- ncol((data.spectra_T_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])
N_L_PNM <- ncol((data.spectra_L_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])
N_L_FTS <- ncol((data.spectra_L_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])

Data <- list(obs_reflectance = (cbind((data.spectra_T_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                     (data.spectra_T_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                     (data.spectra_L_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                     (data.spectra_L_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])),
             is.liana = c(rep(0,2*N),rep(1,2*N)),
             is.FTS = c(rep(0,N),rep(1,N),rep(0,N),rep(1,N)))

Data2 <- list(obs_reflectance = cbind((data.spectra_T_PNM %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                      (data.spectra_T_FTS %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                      (data.spectra_L_PNM %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                      (data.spectra_L_FTS %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N]))

par(mfrow = c(1,1))
matplot(WLs,Data2$obs_reflectance[,1:(N_T_PNM + N_T_FTS)],type = 'l',col = "darkgreen")
matlines(WLs,Data2$obs_reflectance[,(N_T_PNM + N_T_FTS + 1) : ncol(Data2$obs_reflectance)],type = 'l',col = "darkblue")

p5_data <- rrtm:::dataspec_p5[pos,1:5]

Constants <- list(Nwl = Nwl,
                  Nspectra = Nspectra,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = p5_data)

Inits <- list(Nmean = 2,
              Cabmean = 40,
              Carmean = 10,
              Cwmean = 0.01,
              Cmmean = 0.01,

              Standard.Dev = 0.05,

              alpha_N = 0,
              alpha_Cab = 0,
              alpha_Car = 0,
              alpha_Cw = 0,
              alpha_Cm = 0,

              beta_N = 0,
              beta_Cab = 0,
              beta_Car = 0,
              beta_Cw = 0,
              beta_Cm = 0,

              sigma_RE_N = 0,
              sigma_RE_Cab = 0,
              sigma_RE_Car = 0,
              sigma_RE_Cw = 0,
              sigma_RE_Cm = 0,

              gamma_RE_Cab = 0,
              gamma_RE_Car = 0,
              gamma_RE_N = 0,
              gamma_RE_Cm = 0,
              gamma_RE_Cw = 0)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl,Nspectra)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

Nchains = 2
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("Nmean","Cabmean","Carmean","Cwmean","Cmmean","Standard.Dev",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm",
                                    "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm",
                                    "gamma_RE_N","gamma_RE_Cab","gamma_RE_Car","gamma_RE_Cw","gamma_RE_Cm",
                                    "mean_RE_N","mean_RE_Cab","mean_RE_Car","mean_RE_Cw","mean_RE_Cm",
                                    "sigma_RE_N","sigma_RE_Cab","sigma_RE_Car","sigma_RE_Cw","sigma_RE_Cm"),
                       data = Data,
                       inits = Inits,
                       nburnin = 2000,
                       nchains = Nchains,
                       niter = 100000,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,1:121]

param_X = 5

plot(param[,c(0,6,11) + param_X])
plot(param[,which(grepl("gamma_RE_N",colnames(param$chain1)))[1:3]])
