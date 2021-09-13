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

data.spectra <- readRDS("/home/femeunier/Documents/projects/LianaHPDA/data/All.spectra.ID.RDS") %>% group_by(GF,Species,Ind,site) %>% mutate(ind.id = cur_group_id()) %>%
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

N_T_PNM <- ncol((data.spectra_T_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])
N_T_FTS <- ncol((data.spectra_T_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])
N_L_PNM <- ncol((data.spectra_L_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])
N_L_FTS <- ncol((data.spectra_L_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N])

Data <- list(obs_reflectance = cbind((data.spectra_T_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                     (data.spectra_T_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                     (data.spectra_L_PNM %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                     (data.spectra_L_FTS %>% pull(value) %>% matrix(ncol = length(WLs)) %>% t())[,1:N]))

Data2 <- list(obs_reflectance = cbind((data.spectra_T_PNM %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                      (data.spectra_T_FTS %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                      (data.spectra_L_PNM %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N],
                                      (data.spectra_L_FTS %>% pull(value2) %>% matrix(ncol = length(WLs)) %>% t())[,1:N]))

par(mfrow = c(1,1))
matplot(WLs,Data2$obs_reflectance[,1:(N_T_PNM + N_T_FTS)],type = 'l',col = "darkgreen")
matlines(WLs,Data2$obs_reflectance[,(N_T_PNM + N_T_FTS + 1) : ncol(Data2$obs_reflectance)],type = 'l',col = "darkblue")

p5_data <- rrtm:::dataspec_p5[pos,1:5]
matplot(WLs,p5_data,type = 'l')

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

Nchains = 2
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
                       niter = 100000,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,1:24]

param_X = 5

plot(param[,c(0,6,12,18) + param_X])
pairs(as.matrix(param[,c(0,6,12,18) + param_X]), pch = '.')

Nsimu <- 1000

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:24]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:24]))
}


Simu <- cbind(
  # Tree, PNM
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5],
                    Cab = param_all[ileaf,1],
                    Car = param_all[ileaf,2],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4],
                    Cm = param_all[ileaf,3])[["reflectance"]][pos]})),ncol = Nsimu),
  # Tree, FTS
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,17],
                    Cab = param_all[ileaf,1] + param_all[ileaf,13],
                    Car = param_all[ileaf,2] + param_all[ileaf,14],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,16],
                    Cm = param_all[ileaf,3] + param_all[ileaf,15])[["reflectance"]][pos]})),ncol = Nsimu),
  # Liana, PNM
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,11],
                    Cab = param_all[ileaf,1] + param_all[ileaf,7],
                    Car = param_all[ileaf,2] + param_all[ileaf,8],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,10],
                    Cm = param_all[ileaf,3] + param_all[ileaf,9])[["reflectance"]][pos]})),ncol = Nsimu),
  # Liana, FTS
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,11] + param_all[ileaf,17] + param_all[ileaf,23],
                    Cab = param_all[ileaf,1] + param_all[ileaf,7] + param_all[ileaf,13] + param_all[ileaf,19],
                    Car = param_all[ileaf,2] + param_all[ileaf,8] + param_all[ileaf,14] + param_all[ileaf,20],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,10] + param_all[ileaf,16] + param_all[ileaf,22],
                    Cm = param_all[ileaf,3] + param_all[ileaf,9] + param_all[ileaf,15] + param_all[ileaf,21])[["reflectance"]][pos]})),ncol = Nsimu))

# par(mfrow = c(1,1))
# test1 <- rrtm::prospect5(N = 2,
#                          Cab = 70,
#                          Car = 10,
#                          Cbrown = 0,
#                          Cw = 0.015,
#                          Cm = 0.008)[["reflectance"]]
# test2 <- NIMprospect5(N = 2,
#                       Cab = 70,Car = 10,Cw = 0.015,Cm = 0.008,
#                       dataspec_p5 = rrtm:::dataspec_p5[pos,1:5],talf = rrtm:::p45_talf[pos],
#                       t12 = rrtm:::p45_t12[pos],t21 = rrtm:::p45_t21[pos], Nwl = Nwl)
# plot(400:2500,test1,type = 'l')
# lines(WLs,test2, col = "red")

par(mfrow = c(2,2))
ileaf = 1:1000
hist(param_all[ileaf,0 + param_X])                                            # Tree, PNM
hist(param_all[ileaf,0 + param_X] + param_all[ileaf,0 + param_X + 12])        # Tree, FTS
hist(param_all[ileaf,0 + param_X] + param_all[ileaf,0 + param_X + 6])         # Liana, PNM
hist(param_all[ileaf,0 + param_X] + param_all[ileaf,0 + param_X + 6] +
       param_all[ileaf,0 + param_X + 12] + param_all[ileaf,0 + param_X + 18]) # Liana, FTS

# PNM
par(mfrow = c(1,2))
matplot(WLs,Simu[,1:Nsimu],type = 'l',col = "darkgreen",ylim = c(0,0.6))
matlines(WLs,Simu[,2*Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")
lines(WLs,rowMeans(Data2$obs_reflectance[,1:N_T_PNM]),type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))
lines(WLs,rowMeans(Data2$obs_reflectance[,(N_T_PNM + N_T_PNM) + (1:N_L_PNM)]),type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))

# FTS
matplot(WLs,Simu[,Nsimu + (1:Nsimu)],type = 'l',col = "darkgreen",ylim = c(0,0.6))
matlines(WLs,Simu[,3*Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")
lines(WLs,rowMeans(Data2$obs_reflectance[,N_T_PNM + (1:N_T_FTS)]),type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))
lines(WLs,rowMeans(Data2$obs_reflectance[,(N_T_PNM + N_T_PNM + N_L_PNM) + (1:N_L_FTS)]),type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))


df_l <- data.frame(mcmc.out$samples$chain1) %>% rename(N = Nmean,
                                                       Cab = Cabmean,
                                                       Car = Carmean,
                                                       Cw = Cwmean,
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                                                                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                                                                        "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD",
                                                                                        "gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm","gamma_SD")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev",
                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                        "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD",
                        "gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm","gamma_SD"),
               names_to = "param",
               values_to = "value") %>% mutate(param = factor(param,
                                                              levels = c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                                                         "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                                                         "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD",
                                                                         "gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm","gamma_SD")))

vlines <- bind_rows(list(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm")) %>% group_by(param) %>% summarise(m = mean(value),
                                                                                                                    .groups = "keep"),
                         df_l %>% filter(param %in% c("alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm")) %>% group_by(param) %>% summarise(m = 0),
                         df_l %>% filter(param %in% c("beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm")) %>% group_by(param) %>% summarise(m = 0),
                         df_l %>% filter(param %in% c("gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm","gamma_SD")) %>% group_by(param) %>% summarise(m = 0)))


ggplot(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                    "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD",
                                    "gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm","gamma_SD")),
       aes(value)) +
  geom_histogram(aes(y= ..density..),bins = 60, alpha = 0.4, color = "darkgrey") +
  geom_vline(data = vlines,
             aes(xintercept = m)) +
  facet_wrap(~param, scales = "free",nrow = 4) +
  theme_bw()

param.names <- c("Cab","Car","Cm","Cw","N","Standard.Dev")
Nparams <- length(param.names)

df.param.all <- data.frame()
for (iparam in seq(1,Nparams)){
  df.param <- bind_rows(list(data.frame(GF = "Tree",
                                        site = "PNM",
                                        param = param.names[iparam],
                                        value = param_all[,iparam]),
                             data.frame(GF = "Tree",
                                        site = "FTS",
                                        param = param.names[iparam],
                                        value = param_all[,iparam] + param_all[,iparam + 12]),
                             data.frame(GF = "Liana",
                                        site = "PNM",
                                        param = param.names[iparam],
                                        value = param_all[,iparam] + param_all[,iparam + 6]),
                             data.frame(GF = "Liana",
                                        site = "FTS",
                                        param = param.names[iparam],
                                        value = param_all[,iparam] + param_all[,iparam + 6] + param_all[,iparam + 12] + param_all[,iparam + 18])
  ))

  df.param.all <- bind_rows(list(df.param.all,
                                 df.param %>% mutate(id = rep(1:(nrow(df.param)/4),4))))
}

GFs <- c("Tree","Tree","Liana","Liana")
Sites <- c("PNM","FTS","PNM","FTS")

df.param.spectra <- df.param.all
for (i in seq(1,4)) {
  current.spectra <- Simu[,(i-1)*Nsimu + (1:Nsimu)]
  R750 <- current.spectra[which.min(abs(WLs-750)),]
  R705 <- current.spectra[which.min(abs(WLs-705)),]
  R445 <- current.spectra[which.min(abs(WLs-445)),]
  mND705 <- (R750 - R705)/(R750 + R705 - 2*R445)
  mSR705 <- (R750 - R445)/(R705 - R445)
  rededge <- colMeans(current.spectra[WLs %in% seq(720,750),])


  df.param.spectra <- bind_rows(list(df.param.spectra,
                                     data.frame(GF = GFs[i],
                                                site = Sites[i],
                                                param = "mND705",
                                                value = mND705,
                                                id = 1:length(mND705)),
                                     data.frame(GF = GFs[i],
                                                site = Sites[i],
                                                param = "mSR705",
                                                value = mSR705,
                                                id = 1:length(mND705)),
                                     data.frame(GF = GFs[i],
                                                site = Sites[i],
                                                param = "rededge",
                                                value = rededge,
                                                id = 1:length(mND705))
  ))

}

df.param.all.wide <- df.param.spectra %>%
  pivot_wider(names_from = param,
              values_from = value)

ggplot(data = df.param.all.wide,
       aes(x = Cab, y = mND705, shape = as.factor(GF), color = as.factor(GF))) +
  geom_point() +
  facet_wrap(~ site) +
  stat_smooth(method = "lm") +
  scale_shape_manual(values = c(16,1)) +
  theme_bw()

df.param.all.wide %>% group_by(site,GF) %>% summarise(slope = coef(lm(mND705 ~ Cab))[2])

ggplot(data = df.param.all.wide,
       aes(x = Cab, y = mSR705, shape = as.factor(GF), color = as.factor(GF))) +
  geom_point() +
  facet_wrap(~ site) +
  stat_smooth(method = "lm") +
  scale_shape_manual(values = c(16,1)) +
  theme_bw()

df.param.all.wide %>% group_by(site,GF) %>% summarise(slope = coef(lm(mSR705 ~ Cab))[2])


ggplot(data = df.param.all.wide,
       aes(x = Cab, y = rededge, shape = as.factor(GF), color = as.factor(GF))) +
  geom_point() +
  facet_wrap(~ site, nrow = 2) +
  stat_smooth(method = "lm") +
  scale_shape_manual(values = c(16,1)) +
  theme_bw()


ggplot(data = df.param.all %>% filter(!(param == c("Standard.Dev")))) +
  geom_density_ridges(aes(x = value, y = site,fill = GF),alpha = .4, color = NA) +
  facet_wrap(~ param,scale = "free") +
  scale_fill_manual(values = c("darkgreen","darkblue")) +
  labs(x = "",y = "") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25),
        text = element_text(20))

dataSanchez <- readRDS(file = "./data/DATA_Sanchez_2009_Table2.RDS") %>%
  filter(name %in% c("Cab","Car","Cm","Cw","N")) %>% rename(param = name) %>%
  mutate(GF = factor(GF,levels = c("Tree","Liana")))

N.ref <- df.param.all %>% filter(!(param == c("Standard.Dev"))) %>% ungroup() %>%
  filter(param == "N",GF == "Tree",site == "PNM") %>% pull(value) %>% mean()
param2plot <- df.param.all %>% filter(!(param == c("Standard.Dev"))) %>% ungroup() %>%
  mutate(value = case_when(param == "N" ~ value/N.ref,
                           TRUE ~ value))

param2plot %>% filter(param == "N",GF == "Tree",site == "PNM")

ggplot(data = param2plot) +
  geom_density_ridges(aes(x = value, y = as.factor(site),fill = GF),alpha = 0.4, color = NA) +
  # geom_errorbarh(data = dataSanchez,
  #                aes(xmin = m - sd, xmax = m + sd,y = as.factor(Site), color = GF), height = 0.1) +
  geom_point(data = dataSanchez,
             aes(x = m,y = as.factor(Site), color = GF), size = 2, alpha  = 0.4) +
  facet_wrap(~ param,scale = "free") +
  scale_fill_manual(values = c("darkgreen","darkblue")) +
  scale_color_manual(values = c("darkgreen","darkblue")) +
  labs(x = "",y = "") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25),
        text = element_text(20))

# ggplot(data = df.param.all %>% filter(!(param == c("Standard.Dev")))) +
#   geom_density_ridges(aes(x = value, y = site,fill = GF),alpha = .8, color = NA) +
#   facet_wrap(~ param,scale = "free", nrow = 1) +
#   theme_bw()

# Cw <- param2plot %>% filter(param == "Cw") %>% pull(value)
# Cm <- param2plot %>% filter(param == "Cm") %>% pull(value)
# WC = (1-(1/(1+Cw/Cm)))*100
# Cwbis = -Cm*(1 - 1/(1 - WC/100))
# plot(Cw-Cwbis)
# abline(a = 0, b = 1)

ggplot(data = df_l %>% mutate(type = gsub(".*_","",param),
                              effect = gsub("\\_.*","",param)) %>%
         mutate(effect = case_when(effect %in% c("alpha","beta","gamma") ~ effect,
                                   TRUE ~ "base"),
                type = case_when(type %in% c("N","Cab","Car","Cm","Cw") ~ type,
                                 TRUE ~ "Standard.Dev"))) +
  geom_density_ridges(aes(x = value, y = 0,fill = effect),alpha = .8, color = NA) +
  facet_wrap(~ type,scale = "free") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw()


df.data <- bind_rows(list(
  data.frame(
    GF = "Tree",
    site = "PNM",
    wv = WLs,
    R = rowMeans(Data$obs_reflectance[, 1:N_T_PNM])
  ),
  data.frame(
    GF = "Tree",
    site = "FTS",
    wv = WLs,
    R = rowMeans(Data$obs_reflectance[, N_T_PNM + (1:N_T_FTS)])
  ),
  data.frame(
    GF = "Liana",
    site = "PNM",
    wv = WLs,
    R = rowMeans(Data$obs_reflectance[, N_T_PNM + N_T_FTS + (1:N_L_PNM)])
  ),
  data.frame(
    GF = "Liana",
    site = "FTS",
    wv = WLs,
    R = rowMeans(Data$obs_reflectance[, N_T_PNM + N_T_FTS + N_L_PNM + (1:N_L_FTS)])
  )
)) %>% mutate(GF_site = paste(GF, site, sep = "_")) %>%
  dplyr::select(-c("GF", "site")) %>%
  pivot_wider(values_from = "R",
              names_from = "GF_site") %>%
  mutate(
    site.effect = Tree_FTS - Tree_PNM,
    liana.effect = Liana_PNM - Tree_PNM,
    site.liana.effect = Liana_FTS - Tree_PNM
  ) %>%
  pivot_longer(cols = c("Tree_PNM","Tree_FTS", "Liana_PNM", "Liana_FTS",
                        "site.effect","liana.effect","site.liana.effect")) %>%
  mutate(type = "Data")



meanSimu <- cbind(rowMeans(Simu[,1:1000]),rowMeans(Simu[,1001:2000]),rowMeans(Simu[,2001:3000]),rowMeans(Simu[,3001:4000]))

colnames(meanSimu) <- c("Tree_PNM","Tree_FTS", "Liana_PNM", "Liana_FTS")
rownames(meanSimu) <- WLs

df.compiled <- df.data %>% select(-type) %>% rename(Obs = value) %>%
  filter(name %in% c("Tree_PNM", "Tree_FTS", "Liana_PNM", "Liana_FTS")) %>%
  left_join(melt(meanSimu) %>% rename(wv = Var1,
                                      name = Var2,
                                      Mod = value),
            by = c("wv", "name")) %>%
  mutate(site = gsub(".*_","",name),
         GF = gsub("_.*","",name))

ggplot(data = df.compiled) +
  stat_smooth(aes(x = Obs, y = Mod, group = wv),
              method = "lm", col = "darkgrey",se = FALSE) +
  geom_point(aes(x = Obs, y = Mod, color = GF, shape = site)) +
  scale_color_manual(values = c("darkgreen","darkblue")) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw()

df.sum <- df.compiled %>% filter(!is.na(Obs)) %>% group_by(wv) %>% summarise(slope = coef(lm(Mod ~ Obs))[2],
                                                                             r2 = summary(lm(Mod ~ Obs))$r.squared)

summary(df.sum$r2)


# ggplot(data = df.compiled %>% filter(wv %in% (df.sum %>% filter(r2 < 0.5) %>% pull(wv)))) +
#   stat_smooth(aes(x = Obs, y = Mod, group = wv),
#               method = "lm", col = "darkgrey",se = FALSE) +
#   geom_point(aes(x = Obs, y = Mod, color = GF, shape = site)) +
#   scale_color_manual(values = c("darkgreen","darkblue")) +
#   scale_shape_manual(values = c(16,1)) +
#   theme_bw()


df.Simu <- bind_rows(list(
  data.frame(
    GF = "Tree",
    site = "PNM",
    wv = WLs,
    R = meanSimu[,1]
  ),
  data.frame(
    GF = "Tree",
    site = "FTS",
    wv = WLs,
    R = meanSimu[,2]
  ),
  data.frame(
    GF = "Liana",
    site = "PNM",
    wv = WLs,
    R = meanSimu[,3]
  ),
  data.frame(
    GF = "Liana",
    site = "FTS",
    wv = WLs,
    R = meanSimu[,4]
  )
)) %>% mutate(GF_site = paste(GF, site, sep = "_")) %>%
  dplyr::select(-c("GF", "site")) %>%
  pivot_wider(values_from = "R",
              names_from = "GF_site") %>%
  mutate(
    site.effect = Tree_FTS - Tree_PNM,
    liana.effect = Liana_PNM - Tree_PNM,
    site.liana.effect = Liana_FTS - Tree_PNM
  ) %>%
  pivot_longer(cols = c("Tree_PNM","Tree_FTS", "Liana_PNM", "Liana_FTS",
                        "site.effect","liana.effect","site.liana.effect")) %>%
  mutate(type = "Mod")

df.all <- bind_rows(list(df.data,
                         df.Simu)) %>% pivot_wider(values_from = "value",
                                                   names_from = "type")

ggplot(data = df.all %>% filter(name %in% c("Liana_FTS","Liana_PNM",
                                            "Tree_FTS","Tree_PNM")))+
  geom_point(aes(x = Data,y = Mod)) +
  labs(x = "Observed",y = "Modelled") +
  facet_wrap(~name,scales = "free", nrow = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  theme(text = element_text(size = 20))

summary(lm(data = df.all,formula = Mod ~ Data))

df.all %>% filter(name %in% c("Liana_FTS","Liana_PNM",
                              "Tree_FTS","Tree_PNM")) %>% group_by(name) %>% summarise(m =mean(abs(Data),na.rm = TRUE),
                                                                                       r2 = summary(lm(formula = Data ~ Mod))[["adj.r.squared"]],
                                                                                       RMSE = sqrt(c(crossprod(lm(formula = Data ~ Mod)[["residuals"]]))/length(Data[!is.na(Data)])),
                                                                                       .groups = "keep") %>%
  mutate(RMSE.rel = RMSE/m)

ggplot(data = df.all %>% filter(!(name %in% c("Liana_FTS","Liana_PNM",
                                              "Tree_FTS","Tree_PNM"))))+
  geom_point(aes(x = Data,y = Mod)) +
  facet_wrap(~name,scales = "free", nrow = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw()

df.all %>% filter(!(name %in% c("Liana_FTS","Liana_PNM",
                                "Tree_FTS","Tree_PNM"))) %>% group_by(name) %>% summarise(m =mean(abs(Data),na.rm = TRUE),
                                                                                          r2 = summary(lm(formula = Data ~ Mod))[["adj.r.squared"]],
                                                                                          RMSE = sqrt(c(crossprod(lm(formula = Data ~ Mod)[["residuals"]]))/length(Data[!is.na(Data)])),
                                                                                          .groups = "keep") %>%
  mutate(RMSE.rel = RMSE/m)


df.all %>% filter((name %in% c("Liana_FTS","Liana_PNM",
                               "Tree_FTS","Tree_PNM"))) %>% summarise(m =mean(abs(Data),na.rm = TRUE),
                                                                      r2 = summary(lm(formula = Data ~ Mod))[["adj.r.squared"]],
                                                                      RMSE = sqrt(c(crossprod(lm(formula = Data ~ Mod)[["residuals"]]))/length(Data[!is.na(Data)])),
                                                                      .groups = "keep") %>%
  mutate(RMSE.rel = RMSE/m)

# m    r2   RMSE RMSE.rel
# <dbl> <dbl>  <dbl>    <dbl>
#   1 0.257 0.996 0.0107   0.0417
#   1 0.269 0.987 0.0189   0.0701

# name          m    r2    RMSE RMSE.rel
# <chr>     <dbl> <dbl>   <dbl>    <dbl>
# 1 Liana_FTS 0.257 0.995 0.0113    0.0438
# 2 Liana_PNM 0.272 0.998 0.00752   0.0276
# 3 Tree_FTS  0.262 0.995 0.0123    0.0470
# 4 Tree_PNM  0.285 0.996 0.0108    0.0380
#
# name                   m    r2    RMSE RMSE.rel
# <chr>              <dbl> <dbl>   <dbl>    <dbl>
# 1 liana.effect      0.0154 0.913 0.00416    0.270
# 2 site.effect       0.0268 0.809 0.00942    0.351
# 3 site.liana.effect 0.0300 0.738 0.0103     0.344
