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

  reflectance_T_FTS[1:Nwl] <- NIMprospect5(Nmean,
                                           Cabmean,
                                           Carmean,
                                           Cwmean,
                                           Cmmean,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  reflectance_L_PNM[1:Nwl] <- NIMprospect5(Nmean + alpha_N,
                                           Cabmean + alpha_Cab,
                                           Carmean + alpha_Car,
                                           Cwmean + alpha_Cw,
                                           Cmmean + alpha_Cm,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  reflectance_L_FTS[1:Nwl] <- NIMprospect5(Nmean + alpha_N,
                                           Cabmean + alpha_Cab,
                                           Carmean + alpha_Car,
                                           Cwmean + alpha_Cw,
                                           Cmmean + alpha_Cm,
                                           dataspec_p5[,], talf[],t12[],t21[], Nwl)

  creflectance[,1] <- reflectance_T_PNM[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,1] ~ dnorm(creflectance[j,1], sd = Standard.Dev)
  }

  creflectance[,2] <- reflectance_T_FTS[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,2] ~ dnorm(creflectance[j,2], sd = max(0,Standard.Dev))
  }

  creflectance[,3] <- reflectance_L_PNM[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,3] ~ dnorm(creflectance[j,3], sd = max(0,Standard.Dev + alpha_SD))
  }

  creflectance[,4] <- reflectance_L_FTS[1:Nwl]
  for (j in 1:Nwl){
    obs_reflectance[j,4] ~ dnorm(creflectance[j,4], sd = max(0,Standard.Dev + alpha_SD))
  }

  Nmean ~ dunif(1.,5)

  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
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

})

WLa <- 400
WLb <- 2500
Delta_WL <- 50

WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)

data.spectra <- readRDS("./data/All.spectra.ID.RDS") %>%
  group_by(GF,site,Species) %>% mutate(species.id = cur_group_id()) %>% group_by(GF,site) %>% mutate(species.id = species.id - min(species.id) + 1) %>%
  group_by(GF,site,species.id,Ind) %>% mutate(ind.id = cur_group_id()) %>% group_by(GF,site,species.id) %>% mutate(ind.id = ind.id - min(ind.id) + 1) %>%
  group_by(GF,site,species.id,ind.id,name) %>% mutate(leaf.id = cur_group_id()) %>% group_by(GF,site,species.id,ind.id) %>% mutate(leaf.id = leaf.id - min(leaf.id) + 1) %>%
  mutate(value2 = value,
         value = case_when(wv %in% c(400:680,800:2500) ~ value,
                           TRUE ~ NA_real_))

data.spectra_T_PNM <- data.spectra %>% filter(GF == "Tree",
                                              site == "PNM",
                                              wv %in% WLs) %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)

data.spectra_T_FTS <- data.spectra %>% filter(GF == "Tree",
                                              site == "FTS",
                                              wv %in% WLs) %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)

data.spectra_L_PNM <- data.spectra %>% filter(GF == "Liana",
                                              site == "PNM",
                                              wv %in% WLs) %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)

data.spectra_L_FTS <- data.spectra %>% filter(GF == "Liana",
                                              site == "FTS",
                                              wv %in% WLs) %>%
  group_by(wv) %>% summarise(value = mean(value),
                             .groups = "drop") %>% dplyr::select(value)


GFs <- c('Tree','Liana')
sites <- c('PNM','FTS')

maxNspecies <- Inf
maxNind <- Inf

array_obs_reflectance <- array(data = NA, dim = c(Nwl,2,2,30,15,3))
Nspecies <- array(data = NA, dim = c(2,2))
Nind <- array(data = NA, dim = c(2,2,30))

for (iGF in seq(1,length(GFs))){
  for (isite in seq(1,length(sites))){
    cdata <- data.spectra %>% filter(GF == GFs[iGF],
                                     site == sites[isite])
    Nspecies[iGF,isite] <- min(maxNspecies,length(unique(cdata$species.id)))

    for (ispecies in seq(1,Nspecies[iGF,isite])){
      ccdata <- cdata %>% filter(species.id == ispecies)
      Nind[iGF,isite,ispecies] <- min(maxNind,length(unique(ccdata$ind.id)))

      for (iind in seq(1, Nind[iGF,isite,ispecies])){
        cccdata <- ccdata %>% filter(ind.id == iind)

        array_obs_reflectance[,iGF,isite,ispecies,iind,1] <- cccdata %>% filter(name == 1) %>% filter(wv %in% WLs) %>% arrange(wv) %>% pull(value)
        array_obs_reflectance[,iGF,isite,ispecies,iind,2] <- cccdata %>% filter(name == 2) %>% filter(wv %in% WLs) %>% arrange(wv) %>% pull(value)
        array_obs_reflectance[,iGF,isite,ispecies,iind,3] <- cccdata %>% filter(name == 3) %>% filter(wv %in% WLs) %>% arrange(wv) %>% pull(value)
      }
    }
  }
}

maxNspecies <- max(Nspecies,na.rm = TRUE)
maxNind <- max(Nind,na.rm = TRUE)


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
              Standard.Dev = 0.05,

              alpha_N = 0,
              alpha_Cab = 0,
              alpha_Car = 0,
              alpha_Cw = 0,
              alpha_Cm = 0,
              alpha_SD = 0)

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
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD"),
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



param <- MCMCsamples[,1:12]

param_X = 5

plot(param[,c(0,6) + param_X])
pairs(as.matrix(param[,c(0,6) + param_X]), pch = '.')

Nsimu <- 100

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:12]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:12]))[sample(1:(Nchains*Nsimu),Nsimu),]
}

array_mod_reflectance <- array(data = NA, dim = c(dim(array_obs_reflectance)[1:5],Nsimu))
all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <- array(data = NA, dim = c(dim(array_obs_reflectance)[2:5],Nsimu),dimnames = list(c("Tree","Liana"),
                                                                                                                                   c("PNM","FTS"),
                                                                                                                                   seq(1,30),
                                                                                                                                   seq(1:15),
                                                                                                                                   seq(1:Nsimu)))
NGF <- 2
Nsite <- c(2,2)

for(iGF in 1:NGF){
  for (isite in 1:Nsite[iGF]){
    for (ispecies in 1:Nspecies[iGF,isite]){

      print(c(iGF,isite,ispecies))

      for (iind in 1:Nind[iGF,isite,ispecies]){

        cN <- param_all[,"Nmean"] + (iGF - 1)*param_all[,"alpha_N"]
        all_N[iGF,isite,ispecies,iind,] <- cN

        cCab <- param_all[,"Cabmean"] + (iGF - 1)*param_all[,"alpha_Cab"]
        all_Cab[iGF,isite,ispecies,iind,] <- cCab

        cCar <- param_all[,"Carmean"] + (iGF - 1)*param_all[,"alpha_Car"]
        all_Car[iGF,isite,ispecies,iind,] <- cCar

        cCw <- param_all[,"Cwmean"] + (iGF - 1)*param_all[,"alpha_Cw"]
        all_Cw[iGF,isite,ispecies,iind,] <- cCw

        cCm <- param_all[,"Cmmean"] + (iGF - 1)*param_all[,"alpha_Cm"]
        all_Cm[iGF,isite,ispecies,iind,] <- cCm

        tmp <- matrix(unlist(lapply(1:Nsimu,function(ileaf){
          rrtm::prospect5(N = cN[ileaf],
                          Cab = cCab[ileaf],
                          Car = cCar[ileaf],
                          Cbrown = 0,
                          Cw = cCw[ileaf],
                          Cm = cCm[ileaf])[["reflectance"]][pos]})),ncol = Nsimu)



        array_mod_reflectance[,iGF,isite,ispecies,iind,] <- tmp

        # matplot(WLs,tmp[,1:Nsimu],type = 'l',col = "black")
        # matlines(WLs,Data$obs_reflectance[,iGF,isite,ispecies,iind,],col = "red")
        # lines(WLs,apply(Data$obs_reflectance[,iGF,isite,ispecies,iind,],1,mean),col = "red",lwd = 2)

      }
    }
  }
}


# Individual level
Y <- as.vector(apply(array_obs_reflectance,c(1,2,3,4,5),mean))
X <- as.vector(apply(array_mod_reflectance,c(1,2,3,4,5),mean))
plot(X,Y)

abline(a = 0, b = 1, col ='red')
LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

summary(LM)
anovobj<-aov(LM)
allssq<-summary(anovobj)[[1]][,2]

sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

# GF/site level
Y <- as.vector(apply(array_obs_reflectance,c(1,2,3),mean,na.rm = TRUE))
X <- as.vector(apply(array_mod_reflectance,c(1,2,3),mean,na.rm = TRUE))
plot(X,Y)
abline(a = 0, b = 1, col ='red')
LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

summary(LM)

anovobj<-aov(LM)
allssq<-summary(anovobj)[[1]][,2]

sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

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
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,11],
                    Cab = param_all[ileaf,1] + param_all[ileaf,7],
                    Car = param_all[ileaf,2] + param_all[ileaf,8],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,10],
                    Cm = param_all[ileaf,3] + param_all[ileaf,9])[["reflectance"]][pos]})),ncol = Nsimu),
  # Liana, FTS
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,11],
                    Cab = param_all[ileaf,1] + param_all[ileaf,7],
                    Car = param_all[ileaf,2] + param_all[ileaf,8],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,10],
                    Cm = param_all[ileaf,3] + param_all[ileaf,9])[["reflectance"]][pos]})),ncol = Nsimu))

par(mfrow = c(2,2))
ileaf = 1:1000
hist(param_all[ileaf,0 + param_X])                                    # Tree, PNM
hist(param_all[ileaf,0 + param_X])                                    # Tree, FTS
hist(param_all[ileaf,0 + param_X] + param_all[ileaf,0 + param_X + 6]) # Liana, PNM
hist(param_all[ileaf,0 + param_X] + param_all[ileaf,0 + param_X + 6]) # Liana, FTS

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
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                                                                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev",
                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD"),
               names_to = "param",
               values_to = "value")

vlines <- bind_rows(list(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm")) %>% group_by(param) %>% summarise(m = mean(value),
                                                                                                                    .groups = "keep"),
                         df_l %>% filter(param %in% c("alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm")) %>% group_by(param) %>% summarise(m = 0)))

ggplot(df_l %>% filter(param %in% c("N","Cab","Car","Cw","Cm",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm")),
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
                                        value = param_all[,iparam] + param_all[,iparam + 6]),
                             data.frame(GF = "Liana",
                                        site = "FTS",
                                        param = param.names[iparam],
                                        value = param_all[,iparam] + param_all[,iparam + 6])
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

