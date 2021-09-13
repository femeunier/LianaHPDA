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

  reflectance[1:Nwl] <- NIMprospect5(Nmean,
                                     Cabmean,
                                     Carmean,
                                     Cwmean,
                                     Cmmean,
                                     dataspec_p5[,], talf[],t12[],t21[], Nwl)

  creflectance[,1] <- reflectance[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,1] ~ dnorm(creflectance[j,1], sd = Standard.Dev)
  }

  creflectance[,2] <- reflectance[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,2] ~ dnorm(creflectance[j,2], sd = max(0,Standard.Dev))
  }

  creflectance[,3] <- reflectance[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,3] ~ dnorm(creflectance[j,3], sd = max(0,Standard.Dev))
  }

  creflectance[,4] <- reflectance[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,4] ~ dnorm(creflectance[j,4], sd = max(0,Standard.Dev))
  }

  Nmean ~ dunif(1.,5)
  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)
  Standard.Dev ~ dunif(0,1)

})

data.spectra <- readRDS("./data/All.spectra.RDS") %>% group_by(GF,Species,Ind,site) %>% mutate(ind.id = cur_group_id()) %>%
  group_by(GF,Species,Ind,site,name) %>% mutate(leaf.id = cur_group_id()) %>%
  mutate(value2 = value,
         value = case_when(wv %in% c(500:680,800:2500) ~ value,
                           TRUE ~ NA_real_))

data.spectra_T_PNM <- data.spectra %>% filter(GF == "Tree",
                                              site == "PNM") %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)

data.spectra_T_FTS <- data.spectra %>% filter(GF == "Tree",
                                              site == "FTS") %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)

data.spectra_L_PNM <- data.spectra %>% filter(GF == "Liana",
                                              site == "PNM") %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)

data.spectra_L_FTS <- data.spectra %>% filter(GF == "Liana",
                                              site == "FTS") %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)


WLa <- 400
WLb <- 2500
Delta_WL <- 5

WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)

Data <- list(obs_reflectance = cbind(data.spectra_T_PNM,
                                     data.spectra_T_FTS,
                                     data.spectra_L_PNM,
                                     data.spectra_L_FTS))

colnames(Data$obs_reflectance) <- NULL

par(mfrow = c(1,1))

matplot(WLs,Data$obs_reflectance[,1:2],type = 'l',col = "darkgreen") # Trees
matlines(WLs,Data$obs_reflectance[,3:4],type = 'l',col = "darkblue") # Lianas

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
              Standard.Dev = 0.05)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl),
                                         creflectance = c(Nwl,4)),
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
                       niter = 15000,
                       summary = TRUE,
                       WAIC = FALSE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
# matplot((matrix(MCMCsamples$chain1[1,19:1702],ncol = 4)),type = 'l')

param <- MCMCsamples[,1:6]

param_X = 5

plot(param[,c(0) + param_X])
# pairs(as.matrix(param[,c(0) + param_X]), pch = '.')

Nsimu <- 1000

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:6]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:6]))
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
    rrtm::prospect5(N = param_all[ileaf,5],
                    Cab = param_all[ileaf,1],
                    Car = param_all[ileaf,2],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4],
                    Cm = param_all[ileaf,3])[["reflectance"]][pos]})),ncol = Nsimu),
  # Liana, PNM
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5],
                    Cab = param_all[ileaf,1],
                    Car = param_all[ileaf,2],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4],
                    Cm = param_all[ileaf,3])[["reflectance"]][pos]})),ncol = Nsimu),
  # Liana, FTS
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5],
                    Cab = param_all[ileaf,1],
                    Car = param_all[ileaf,2],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4],
                    Cm = param_all[ileaf,3])[["reflectance"]][pos]})),ncol = Nsimu))

par(mfrow = c(2,2))
ileaf = 1:1000
hist(param_all[ileaf,0 + param_X]) # Tree, PNM
hist(param_all[ileaf,0 + param_X]) # Tree, FTS
hist(param_all[ileaf,0 + param_X]) # Liana, PNM
hist(param_all[ileaf,0 + param_X]) # Liana, FTS

# PNM
par(mfrow = c(1,2))
matplot(WLs,Simu[,1:Nsimu],type = 'l',col = "darkgreen",ylim = c(0,0.6))
matlines(WLs,Simu[,2*Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")
matlines(WLs,Data$obs_reflectance[,c(1,3)],type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))

# FTS
matplot(WLs,Simu[,Nsimu + (1:Nsimu)],type = 'l',col = "darkgreen",ylim = c(0,0.6))
matlines(WLs,Simu[,3*Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")
matlines(WLs,Data$obs_reflectance[,c(2,4)],type = 'l',col = c("black","black"),lwd = 2,lty = c(1,2))

df_l <- data.frame(mcmc.out$samples$chain1) %>% rename(N = Nmean,
                                                       Cab = Cabmean,
                                                       Car = Carmean,
                                                       Cw = Cwmean,
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev"),
               names_to = "param",
               values_to = "value")

vlines <- bind_rows(list(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm")) %>% group_by(param) %>% summarise(m = mean(value),
                                                                                                                    .groups = "keep")))

ggplot(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm")),
       aes(value)) +
  geom_histogram(aes(y= ..density..),bins = 60, alpha = 0.4, color = "darkgrey") +
  geom_vline(data = vlines,
             aes(xintercept = m)) +
  facet_wrap(~param, scales = "free",nrow = 3) +
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
                                        value = param_all[,iparam]),
                             data.frame(GF = "Liana",
                                        site = "PNM",
                                        param = param.names[iparam],
                                        value = param_all[,iparam]),
                             data.frame(GF = "Liana",
                                        site = "FTS",
                                        param = param.names[iparam],
                                        value = param_all[,iparam])
  ))

  df.param.all <- bind_rows(list(df.param.all,
                                 df.param))
}

ggplot(data = df.param.all) +
  geom_density_ridges(aes(x = value, y = site,fill = GF),alpha = .8, color = NA) +
  facet_wrap(~ param,scale = "free") +
  theme_bw()

ggplot(data = df.param.all %>% filter(!(param == c("Standard.Dev")))) +
  geom_density_ridges(aes(x = value, y = site,fill = GF),alpha = .8, color = NA) +
  facet_wrap(~ param,scale = "free", nrow = 1) +
  theme_bw()

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
    R = Data$obs_reflectance[, 1]
  ),
  data.frame(
    GF = "Tree",
    site = "FTS",
    wv = WLs,
    R = Data$obs_reflectance[, 2]
  ),
  data.frame(
    GF = "Liana",
    site = "PNM",
    wv = WLs,
    R = Data$obs_reflectance[, 3]
  ),
  data.frame(
    GF = "Liana",
    site = "FTS",
    wv = WLs,
    R = Data$obs_reflectance[, 4]
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
  facet_wrap(~name,scales = "free", nrow = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw()

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
# 1 0.269 0.987 0.0189   0.0701

