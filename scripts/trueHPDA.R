rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(reshape2)
library(rrtm)

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

  mean_N ~ dunif(1.,5)
  mean_Cab ~ dunif(0,200)
  mean_Car ~ dunif(0,100)
  mean_Cw ~ dunif(0.,0.1)
  mean_Cm ~ dunif(0.,0.1)

  # GF effect
  alpha_N ~ dunif(-1,1)
  alpha_Cab ~ dunif(-100,100)
  alpha_Car ~ dunif(-50,50)
  alpha_Cw ~ dunif(-0.05,0.05)
  alpha_Cm ~ dunif(-0.05,0.05)


  for(iGF in 1:NGF){

    # Site effect
    beta_N[iGF] ~ dunif(-1,1)
    beta_Cab[iGF] ~ dunif(-100,100)
    beta_Car[iGF] ~ dunif(-50,50)
    beta_Cw[iGF] ~ dunif(-0.05,0.05)
    beta_Cm[iGF] ~ dunif(-0.05,0.05)

    for (isite in 1:Nsite[iGF]){

      for (ispecies in 1:Nspecies[iGF,isite]){

        # Species effect
        gamma_N[iGF,isite,ispecies] ~ dunif(-1,1)
        gamma_Cab[iGF,isite,ispecies] ~ dunif(-100,100)
        gamma_Car[iGF,isite,ispecies] ~ dunif(-50,50)
        gamma_Cw[iGF,isite,ispecies] ~ dunif(-0.05,0.05)
        gamma_Cm[iGF,isite,ispecies] ~ dunif(-0.05,0.05)

        for (iind in 1:Nind[iGF,isite,ispecies]){

          # Ind effect
          delta_N[iGF,isite,ispecies,iind] ~ dunif(-1,1)
          delta_Cab[iGF,isite,ispecies,iind] ~ dunif(-100,100)
          delta_Car[iGF,isite,ispecies,iind] ~ dunif(-50,50)
          delta_Cw[iGF,isite,ispecies,iind] ~ dunif(-0.05,0.05)
          delta_Cm[iGF,isite,ispecies,iind] ~ dunif(-0.05,0.05)

          reflectance[1:Nwl,iGF,isite,ispecies,iind] <-
            NIMprospect5(mean_N + (iGF - 1)*alpha_N + (isite - 1)*beta_N[iGF] + gamma_N[iGF,isite,ispecies] + delta_N[iGF,isite,ispecies,iind],
                         mean_Cab + (iGF - 1)*alpha_Cab + (isite - 1)*beta_Cab[iGF] + gamma_Cab[iGF,isite,ispecies] + delta_Cab[iGF,isite,ispecies,iind],
                         mean_Car + (iGF - 1)*alpha_Car + (isite - 1)*beta_Car[iGF] + gamma_Car[iGF,isite,ispecies] + delta_Car[iGF,isite,ispecies,iind],
                         mean_Cw + (iGF - 1)*alpha_Cw + (isite - 1)*beta_Cw[iGF] + gamma_Cw[iGF,isite,ispecies] + delta_Cw[iGF,isite,ispecies,iind],
                         mean_Cm + (iGF - 1)* alpha_Cm + (isite - 1)*beta_Cm[iGF] + gamma_Cm[iGF,isite,ispecies] + delta_Cm[iGF,isite,ispecies,iind],
                         dataspec_p5[,], talf[],t12[],t21[], Nwl)

          for (ileaf in 1:Nleaves){
            for (j in 1:Nwl){
              obs_reflectance[j,iGF,isite,ispecies,iind,ileaf] ~ dnorm(reflectance[j,iGF,isite,ispecies,iind], sd = Standard.Dev)
            }
          }

        }
      }
    }
  }

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

Inits <- list(mean_N = 2,
              mean_Cab = 40,
              mean_Car = 10,
              mean_Cw = 0.01,
              mean_Cm = 0.01,

              Standard.Dev = 0.05,

              alpha_N = 0,
              alpha_Cab = 0,
              alpha_Car = 0,
              alpha_Cw = 0,
              alpha_Cm = 0,

              beta_N = c(0,0),
              beta_Cab = c(0,0),
              beta_Car = c(0,0),
              beta_Cw = c(0,0),
              beta_Cm = c(0,0))

Data <- list(obs_reflectance = array_obs_reflectance)
p5_data <- rrtm:::dataspec_p5[pos,1:5]

Nleaves <- 3
NGF <- 2
Nsite <- c(2,2)

Constants <- list(Nwl = Nwl,
                  Nleaves = Nleaves,
                  NGF = NGF,
                  Nsite = Nsite,
                  Nspecies = Nspecies,
                  Nind = Nind,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = p5_data)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl,2,2,maxNspecies,maxNind)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

Nchains = 1
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("mean_N","mean_Cab","mean_Car","mean_Cw","mean_Cm","Standard.Dev",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm",
                                    "beta_N","beta_Cab","beta_Car","beta_Cw","beta_Cm",
                                    "gamma_N","gamma_Cab","gamma_Car","gamma_Cw","gamma_Cm",
                                    "delta_N","delta_Cab","delta_Car","delta_Cw","delta_Cm"),
                       data = Data,
                       inits = Inits,
                       nburnin = 2000,
                       thin = 10,
                       nchains = Nchains,
                       niter = 3000,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,]

param_X = 5

plot(param[,c("mean_N","alpha_N","beta_N[1]","beta_N[2]")])
pairs(as.matrix(param[,c("mean_N","alpha_N","beta_N[1]","beta_N[2]")]), pch = '.')

plot(param[,which(grepl("delta_Cw",colnames(param)))[1:4]])


Nsimu <- min(1000,nrow(MCMCsamples))
pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[pos.simu,]))

# iGF <- isite <- ispecies <- iind <- 2

array_mod_reflectance <- array(data = NA, dim = c(dim(array_obs_reflectance)[1:5],Nsimu))
all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <- array(data = NA, dim = c(dim(array_obs_reflectance)[2:5],Nsimu),dimnames = list(c("Tree","Liana"),
                                                                                                                                   c("PNM","FTS"),
                                                                                                                                   seq(1,30),
                                                                                                                                   seq(1:15),
                                                                                                                                   seq(1:Nsimu)))

for(iGF in 1:NGF){
  for (isite in 1:Nsite[iGF]){
    for (ispecies in 1:Nspecies[iGF,isite]){
      for (iind in 1:Nind[iGF,isite,ispecies]){

        cN <- param_all[,"mean_N"] + (iGF - 1)*param_all[,"alpha_N"] + (isite - 1)*param_all[,paste0("beta_N[",iGF,"]")] +
          param_all[,paste0("gamma_N[",iGF,", ",isite,", ",ispecies,"]")] + param_all[,paste0("delta_N[",iGF,", ",isite,", ",ispecies,", ",iind,"]")]
        all_N[iGF,isite,ispecies,iind,] <- cN

        cCab <- param_all[,"mean_Cab"] + (iGF - 1)*param_all[,"alpha_Cab"] + (isite - 1)*param_all[,paste0("beta_Cab[",iGF,"]")] +
          param_all[,paste0("gamma_Cab[",iGF,", ",isite,", ",ispecies,"]")] + param_all[,paste0("delta_Cab[",iGF,", ",isite,", ",ispecies,", ",iind,"]")]
        all_Cab[iGF,isite,ispecies,iind,] <- cCab

        cCar <- param_all[,"mean_Car"] + (iGF - 1)*param_all[,"alpha_Car"] + (isite - 1)*param_all[,paste0("beta_Car[",iGF,"]")] +
          param_all[,paste0("gamma_Car[",iGF,", ",isite,", ",ispecies,"]")] + param_all[,paste0("delta_Car[",iGF,", ",isite,", ",ispecies,", ",iind,"]")]
        all_Car[iGF,isite,ispecies,iind,] <- cCar

        cCw <- param_all[,"mean_Cw"] + (iGF - 1)*param_all[,"alpha_Cw"] + (isite - 1)*param_all[,paste0("beta_Cw[",iGF,"]")] +
          param_all[,paste0("gamma_Cw[",iGF,", ",isite,", ",ispecies,"]")] + param_all[,paste0("delta_Cw[",iGF,", ",isite,", ",ispecies,", ",iind,"]")]
        all_Cw[iGF,isite,ispecies,iind,] <- cCw

        cCm <- param_all[,"mean_Cm"] + (iGF - 1)*param_all[,"alpha_Cm"] + (isite - 1)*param_all[,paste0("beta_Cm[",iGF,"]")] +
          param_all[,paste0("gamma_Cm[",iGF,", ",isite,", ",ispecies,"]")] + param_all[,paste0("delta_Cm[",iGF,", ",isite,", ",ispecies,", ",iind,"]")]
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

all_parameters <- bind_rows(list(melt(all_N) %>% mutate(param = "N"),
                             melt(all_Cab) %>% mutate(param = "Cab"),
                             melt(all_Car) %>% mutate(param = "Car"),
                             melt(all_Cw) %>% mutate(param = "Cw"),
                             melt(all_Cm) %>% mutate(param = "Cm"))) %>% filter(!is.na(value)) %>% rename(GF = Var1,
                                                                                                          site = Var2,
                                                                                                          species = Var3,
                                                                                                          Ind = Var4,
                                                                                                          simu = Var5)

ggplot(data = all_parameters) +
  geom_density_ridges(aes(x = value,y = as.factor(site),fill = as.factor(GF)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()


dataSanchez <- readRDS(file = "./data/DATA_Sanchez_2009_Table2.RDS") %>%
  filter(name %in% c("Cab","Car","Cm","Cw","N")) %>% rename(param = name) %>%
  mutate(GF = factor(GF,levels = c("Tree","Liana")))

N.ref <- all_parameters %>% filter(!(param == c("Standard.Dev"))) %>% ungroup() %>%
  filter(param == "N",GF == "Tree",site == "PNM") %>% pull(value) %>% mean()

param2plot <- all_parameters %>% filter(!(param == c("Standard.Dev"))) %>% ungroup() %>%
  mutate(value = case_when(param == "N" ~ value/N.ref,
                           TRUE ~ value))


ggplot(data = param2plot) +
  geom_density_ridges(aes(x = value, y = as.factor(site),fill = GF),alpha = 0.4, color = NA) +
  geom_point(data = dataSanchez,
             aes(x = m,y = as.factor(Site), color = GF), size = 2, alpha  = 0.4) +
  facet_wrap(~ param,scale = "free") +
  scale_fill_manual(values = c("darkgreen","darkblue")) +
  scale_color_manual(values = c("darkgreen","darkblue")) +
  labs(x = "",y = "") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25),
        text = element_text(20))


merged.params <- param2plot %>% group_by(param,GF,site) %>% summarise(mod = mean(value,na.rm = TRUE),
                                                                      .groups = "keep") %>%
  left_join(dataSanchez %>% dplyr::select(Site,GF,param,m) %>% rename(site = Site,
                                                                      obs = m),
            by = c("site","GF","param"))

ggplot(data = merged.params) +
  geom_point(aes(x = mod, y = obs)) +
  facet_wrap(~param, scales = "free") +
  geom_abline(slope = 1, color = "red", linetype = 2) +
  theme_bw()

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

Y = as.vector(array_obs_reflectance[1,,,,,1])
X = apply(all_N,c(1,2,3,4),mean)
X1 = apply(all_Cab,c(1,2,3,4),mean)
varpart(Y[!is.na(Y)],X[!is.na(X)],X1[!is.na(X1)])
