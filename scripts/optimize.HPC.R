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

    if (N < 1.1 | Car < 0 | Cab < 0 | Cm <= 0 | Cw <0) {
      return(-9999*(reflectance**0))
    } else{
      return(reflectance)
    }
  })

run_prospect5 <- nimbleCode({

  reflectance_T_PNM[1:Nwl] <- NIMprospect5(Nmean,
                                           Cabmean,
                                           Carmean,
                                           Cwmean,
                                           Cmmean,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  reflectance_T_FTS[1:Nwl] <- NIMprospect5(Nmean + beta_N,
                                           Cabmean + beta_Cab,
                                           Carmean + beta_Car,
                                           Cwmean + beta_Cw,
                                           Cmmean + beta_Cm,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  reflectance_L_PNM[1:Nwl] <- NIMprospect5(Nmean + alpha_N,
                                           Cabmean + alpha_Cab,
                                           Carmean + alpha_Car,
                                           Cwmean + alpha_Cw,
                                           Cmmean + alpha_Cm,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  reflectance_L_FTS[1:Nwl] <- NIMprospect5(Nmean + alpha_N + beta_N + gamma_N,
                                           Cabmean + alpha_Cab + beta_Cab + gamma_Cab,
                                           Carmean + alpha_Car + beta_Car + gamma_Car,
                                           Cwmean + alpha_Cw + beta_Cw + gamma_Cw,
                                           Cmmean + alpha_Cm + beta_Cm + gamma_Cm,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  creflectance[,1] <- reflectance_T_PNM[1:Nwl]
  creflectance[,2] <- reflectance_T_FTS[1:Nwl]
  creflectance[,3] <- reflectance_L_PNM[1:Nwl]
  creflectance[,4] <- reflectance_L_FTS[1:Nwl]

  for (i in 1:N_T_PNM){
    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(creflectance[j,1], sd = Standard.Dev)
    }
  }

  for (i in 1:N_T_FTS){
    for (j in 1:Nwl){
      obs_reflectance[j,i + (N_T_PNM)] ~ dnorm(creflectance[j,2], sd = max(0,Standard.Dev + beta_SD))
    }
  }

  for (i in 1:N_L_PNM){
    for (j in 1:Nwl){
      obs_reflectance[j,i + (N_T_PNM + N_T_FTS)] ~ dnorm(creflectance[j,3], sd = max(0,Standard.Dev + alpha_SD))
    }
  }

  for (i in 1:N_L_FTS){
    for (j in 1:Nwl){
      obs_reflectance[j,i + (N_T_PNM + N_T_FTS + N_L_PNM)] ~ dnorm(creflectance[j,4], sd = max(0,Standard.Dev + alpha_SD + beta_SD + gamma_SD))
    }
  }

  Nmean ~ dunif(1.,5)
  Cabmean ~ dunif(0,200)
  Carmean ~ dunif(0,100)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)
  Standard.Dev ~ dunif(0,1)

  # GF effect
  alpha_N ~ dunif(-1,1)
  alpha_Cab ~ dunif(-400,400)
  alpha_Car ~ dunif(-200,200)
  alpha_Cw ~ dunif(-0.2,0.2)
  alpha_Cm ~ dunif(-0.2,0.2)
  alpha_SD ~ dunif(-0.2,0.2)

  # Site effect
  beta_N ~ dunif(-1,1)
  beta_Cab ~ dunif(-400,400)
  beta_Car ~ dunif(-200,200)
  beta_Cw ~ dunif(-0.2,0.2)
  beta_Cm ~ dunif(-0.2,0.2)
  beta_SD ~ dunif(-0.2,0.2)

  # Interaction effect
  gamma_N ~ dunif(-1,1)
  gamma_Cab ~ dunif(-400,400)
  gamma_Car ~ dunif(-200,200)
  gamma_Cw ~ dunif(-0.2,0.2)
  gamma_Cm ~ dunif(-0.2,0.2)
  gamma_SD ~ dunif(-0.2,0.2)

})

WLa <- 400
WLb <- 2500
Delta_WL <- 20

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

N_T_PNM <- ncol((data.spectra_T_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()))
N_T_FTS <- ncol((data.spectra_T_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()))
N_L_PNM <- ncol((data.spectra_L_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()))
N_L_FTS <- ncol((data.spectra_L_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()))

Data <- list(obs_reflectance = cbind((data.spectra_T_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()),
                                     (data.spectra_T_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()),
                                     (data.spectra_L_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t()),
                                     (data.spectra_L_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())))

Data2 <- list(obs_reflectance = cbind((data.spectra_T_PNM %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t()),
                                      (data.spectra_T_FTS %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t()),
                                      (data.spectra_L_PNM %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t()),
                                      (data.spectra_L_FTS %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())))

p5_data <- rrtm:::dataspec_p5[pos,1:5]

Constants <- list(Nwl = Nwl,
                  N_T_PNM = N_T_PNM,
                  N_T_FTS = N_T_FTS,
                  N_L_PNM = N_L_PNM,
                  N_L_FTS = N_L_FTS,
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
              alpha_SD = 0,

              beta_N = 0,
              beta_Cab = 0,
              beta_Car = 0,
              beta_Cw = 0,
              beta_Cm = 0,
              beta_SD = 0,

              gamma_N = 0,
              gamma_Cab = 0,
              gamma_Car = 0,
              gamma_Cw = 0,
              gamma_Cm = 0,
              gamma_SD = 0)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance_T_PNM = c(Nwl),
                                         reflectance_T_FTS = c(Nwl),
                                         reflectance_L_PNM = c(Nwl),
                                         reflectance_L_FTS = c(Nwl),
                                         creflectance = c(Nwl,4)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

Nchains = 1
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("Nmean","Cabmean","Carmean","Cwmean","Cmmean","Standard.Dev",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                    "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD",
                                    "gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm","gamma_SD"),
                       data = Data,
                       inits = Inits,
                       nburnin = 2000,
                       nchains = Nchains,
                       niter = 10000,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

saveRDS(mcmc.out,"./outputs/mcmc.RDS")

