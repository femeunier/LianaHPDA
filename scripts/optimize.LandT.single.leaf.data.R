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

  reflectance_L_FTS[1:Nwl] <- NIMprospect5(Nmean + alpha_N + beta_N,
                                           Cabmean + alpha_Cab + beta_Cab,
                                           Carmean + alpha_Car + beta_Car,
                                           Cwmean + alpha_Cw + beta_Cw,
                                           Cmmean + alpha_Cm + beta_Cm,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  creflectance[,1] <- reflectance_T_PNM[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,1] ~ dnorm(creflectance[j,1], sd = Standard.Dev)
  }

  creflectance[,2] <- reflectance_T_FTS[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,2] ~ dnorm(creflectance[j,2], sd = max(0,Standard.Dev + beta_SD))
  }

  creflectance[,3] <- reflectance_L_PNM[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,3] ~ dnorm(creflectance[j,3], sd = max(0,Standard.Dev + alpha_SD))
  }

  creflectance[,4] <- reflectance_L_FTS[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,4] ~ dnorm(creflectance[j,4], sd = max(0,Standard.Dev + alpha_SD + beta_SD))
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

})

WLa <- 400
WLb <- 2500
Delta_WL <- 5

WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)

data.spectra <- readRDS("./data/All.spectra.ID.RDS") %>% group_by(GF,Species,Ind,site) %>% mutate(ind.id = cur_group_id()) %>%
  group_by(GF,Species,Ind,site,name) %>% mutate(leaf.id = cur_group_id()) %>%
  mutate(value2 = value,
         value = case_when(wv %in% c(500:680,800:2500) ~ value,
                           TRUE ~ NA_real_))

data.spectra_T_PNM <- data.spectra %>% filter(GF == "Tree",
                                              site == "PNM") %>%
  group_by(wv) %>% filter(wv %in% WLs) %>% arrange(wv) %>%
  summarise(value = mean(value),
            value2 = mean(value2),
            .groups = "drop") %>% dplyr::select(value,value2)

data.spectra_T_FTS <- data.spectra %>% filter(GF == "Tree",
                                              site == "FTS") %>%
  group_by(wv) %>% filter(wv %in% WLs) %>% arrange(wv) %>%
  summarise(value = mean(value),
            value2 = mean(value2),
            .groups = "drop") %>% dplyr::select(value,value2)

data.spectra_L_PNM <- data.spectra %>% filter(GF == "Liana",
                                              site == "PNM") %>%
  group_by(wv) %>% filter(wv %in% WLs) %>% arrange(wv) %>%
  summarise(value = mean(value),
            value2 = mean(value2),
            .groups = "drop") %>% dplyr::select(value,value2)

data.spectra_L_FTS <- data.spectra %>% filter(GF == "Liana",
                                              site == "FTS") %>%
  group_by(wv) %>% filter(wv %in% WLs) %>% arrange(wv) %>%
  summarise(value = mean(value),
            value2 = mean(value2),
            .groups = "drop") %>% dplyr::select(value,value2)


Data <- list(obs_reflectance = cbind(data.spectra_T_PNM %>% pull(value),
                                     data.spectra_T_FTS %>% pull(value),
                                     data.spectra_L_PNM %>% pull(value),
                                     data.spectra_L_FTS %>% pull(value)))

Data2 <- list(obs_reflectance = cbind(data.spectra_T_PNM %>% pull(value2),
                                      data.spectra_T_FTS %>% pull(value2),
                                      data.spectra_L_PNM %>% pull(value2),
                                      data.spectra_L_FTS %>% pull(value2)))

par(mfrow = c(1,1))
matplot(WLs,Data2$obs_reflectance[,1:2],type = 'l',col = "darkgreen")
matlines(WLs,Data2$obs_reflectance[,3:4],type = 'l',col = "darkblue")

Constants <- list(Nwl = Nwl,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = rrtm:::dataspec_p5[pos,1:5])

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
              beta_SD = 0)

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
                                    "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD"),
                       data = Data,
                       inits = Inits,
                       nburnin = 2000,
                       nchains = Nchains,
                       niter = 10000,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,1:18]

param_X = 4

plot(param[,c(0,6,12) + param_X])
pairs(as.matrix(param[,c(0,6,12) + param_X]), pch = '.')

Nsimu <- 1000

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:18]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:18]))
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
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,11] + param_all[ileaf,17],
                    Cab = param_all[ileaf,1] + param_all[ileaf,7] + param_all[ileaf,13],
                    Car = param_all[ileaf,2] + param_all[ileaf,8] + param_all[ileaf,14],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,10] + param_all[ileaf,16],
                    Cm = param_all[ileaf,3] + param_all[ileaf,9] + param_all[ileaf,15])[["reflectance"]][pos]})),ncol = Nsimu))

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
       param_all[ileaf,0 + param_X + 12])                                     # Liana, FTS

# PNM
par(mfrow = c(1,2))
matplot(WLs,Simu[,1:Nsimu],type = 'l',col = "darkgreen",ylim = c(0,0.6))
matlines(WLs,Simu[,2*Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")
matlines(WLs,Data2$obs_reflectance[,c(1,3)],type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))

# FTS
matplot(WLs,Simu[,Nsimu + (1:Nsimu)],type = 'l',col = "darkgreen",ylim = c(0,0.6))
matlines(WLs,Simu[,3*Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")
matlines(WLs,Data2$obs_reflectance[,c(2,4)],type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))

df_l <- data.frame(mcmc.out$samples$chain1) %>% rename(N = Nmean,
                                                       Cab = Cabmean,
                                                       Car = Carmean,
                                                       Cw = Cwmean,
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                                                                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                                                                        "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev",
                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                        "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD"),
               names_to = "param",
               values_to = "value") %>% mutate(param = factor(param,
                                                              levels = c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                                                         "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                                                         "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD")))

vlines <- bind_rows(list(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm")) %>% group_by(param) %>% summarise(m = mean(value),
                                                                                                                    .groups = "keep"),
                         df_l %>% filter(param %in% c("alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm")) %>% group_by(param) %>% summarise(m = 0),
                         df_l %>% filter(param %in% c("beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm")) %>% group_by(param) %>% summarise(m = 0)))


ggplot(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD",
                                    "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm","beta_SD")),
       aes(value)) +
  geom_histogram(aes(y= ..density..),bins = 60, alpha = 0.4, color = "darkgrey") +
  geom_vline(data = vlines,
             aes(xintercept = m)) +
  facet_wrap(~param, scales = "free",nrow = 3) +
  theme_bw()

