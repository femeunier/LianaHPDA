rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(LianaHPDA)
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


  sd_species_N ~ dunif(0.,1)
  sd_species_Cab ~ dunif(0.,50)
  sd_species_Car ~ dunif(0.,20)
  sd_species_Cw ~ dunif(0.,0.05)
  sd_species_Cm ~ dunif(0.,0.05)

  for (ispecies in seq(1,Nspecies)){

    sd_leaf_N[ispecies] ~ dunif(0.0,0.5)
    sd_leaf_Cab[ispecies] ~ dunif(0.,30)
    sd_leaf_Car[ispecies] ~ dunif(0.,10)
    sd_leaf_Cw[ispecies] ~ dunif(0.,0.02)
    sd_leaf_Cm[ispecies] ~ dunif(0.,0.02)

    # Species fixed effect

    nu_species_N[ispecies] ~ dnorm(0, sd = sd_species_N)
    nu_species_Cab[ispecies] ~ dnorm(0, sd = sd_species_Cab)
    nu_species_Car[ispecies] ~ dnorm(0, sd = sd_species_Car)
    nu_species_Cw[ispecies] ~ dnorm(0, sd = sd_species_Cw)
    nu_species_Cm[ispecies] ~ dnorm(0, sd = sd_species_Cm)

    for (ileaf in seq(1,Nleaves[ispecies])){

      # leaf random effect
      nu_leaf_N[ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_N[ispecies])
      nu_leaf_Cab[ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Cab[ispecies])
      nu_leaf_Car[ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Car[ispecies])
      nu_leaf_Cw[ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Cw[ispecies])
      nu_leaf_Cm[ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Cm[ispecies])

      reflectance[1:Nwl,ispecies,ileaf] <- NIMprospect5(Nmean + nu_species_N[ispecies] + nu_leaf_N[ispecies,ileaf],
                                                        Cabmean + nu_species_Cab[ispecies] + nu_leaf_Cab[ispecies,ileaf],
                                                        Carmean + nu_species_Car[ispecies] + nu_leaf_Car[ispecies,ileaf],
                                                        Cwmean + nu_species_Cw[ispecies] + nu_leaf_Cw[ispecies,ileaf],
                                                        Cmmean + nu_species_Cm[ispecies] + nu_leaf_Cm[ispecies,ileaf],
                                                        dataspec_p5[,], talf[],t12[],t21[], Nwl)

      for (j in 1:Nwl){
        obs_reflectance[j,ispecies,ileaf] ~ dnorm(reflectance[j,ispecies,ileaf], sd = max(0,Standard.Dev))
      }
    }

  }

  Standard.Dev ~ dunif(0,1)

  Nmean ~ dunif(1.1,5)
  Cabmean ~ dunif(0,250)
  Carmean ~ dunif(0,100)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)

})

WLa <- 400
WLb <- 2500
Delta_WL <- 20

WLs <- seq(WLa,WLb,Delta_WL)
# WLs <- WLs[(WLs < 650 | WLs > 800) & (WLs > 450)]
pos <- which((WLa:WLb %in% WLs))
Nwl <- length(pos)

data.raw <- LianaHPDA::array_obs_reflectance[pos,,,,,]

data.mean <- data.raw[,2,1,3:6,1:4,1:3]
dims <- dim(data.mean)
data.minus1d <- array(data = data.mean, dim = c(dims[1:(length(dims)-2)],dims[length(dims)-1]*dims[length(dims)]))

data.2d.NA <- matrix(data.mean,nrow = Nwl)

par(mfrow = c(1,1))
matplot(WLs,data.2d.NA,type = 'l')

Data <- list(obs_reflectance = data.minus1d)

colnames(Data$obs_reflectance) <- NULL

Nspecies <- dims[2]
Nleaves <- rep(dims[3]*dims[4],Nspecies)
Constants <- list(Nwl = Nwl,
                  Nspecies = Nspecies,
                  Nleaves = Nleaves,
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
                                         reflectance = c(Nwl,max(Nspecies),max(Nleaves))),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

Nchains = 2
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("Nmean","Cabmean","Carmean","Cwmean","Cmmean",
                                    "Standard.Dev",
                                    "nu_species_N","nu_species_Cab","nu_species_Car","nu_species_Cw","nu_species_Cm",
                                    "nu_leaf_N","nu_leaf_Cab","nu_leaf_Car","nu_leaf_Cw","nu_leaf_Cm"),
                       data = Data,
                       inits = Inits,
                       nburnin = 2000,
                       nchains = Nchains,
                       niter = 4000,
                       thin = 10,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)
# system2("rsync",paste("-avz","hpc:/data/gent/vo/000/gvo00074/felicien/R/LianaHPDA/outputs/MCMC.multiplespecies.RDS","./outputs/"))
# mcmc.out <- readRDS("./outputs/MCMC.multiplespecies.RDS")

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,]

Nsimu <- min(1000,nrow(MCMCsamples[[1]]))
pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,]))[sample(1:(Nchains*Nsimu),Nsimu),]

plot(param[,c("Nmean","Cabmean","Carmean","Cwmean","Cmmean")])
plot(param[,"nu_species_N[4]"])
plot(param[,"nu_leaf_N[2, 2]"])

hist(as.vector(as.matrix(param[,c("Nmean")])))
hist(as.vector(as.matrix(param[,which(grepl("nu_species_N",colnames(param$chain1)))])))
hist(as.vector(as.matrix(param[,which(grepl("nu_leaf_N",colnames(param$chain1)))])))

hist(as.vector(as.matrix(param[,c("Nmean")])))

array_mod_reflectance <- array(data = NA,c(dim(data.minus1d),Nsimu))

all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <-
  array(data = NA, c(Nspecies,max(Nleaves),Nsimu))

species_effect_N <- species_effect_Cab <- species_effect_Car <-
  species_effect_Cw <- species_effect_Cm <-
  array(data = NA, dim = c(Nspecies,Nsimu))

leaf_effect_N <- leaf_effect_Cab <- leaf_effect_Car <-
  leaf_effect_Cm <- leaf_effect_Cw <- array(data = NA, dim = c(Nspecies,max(Nleaves),Nsimu))

RMSE <- array(data = NA, dim = c(Nspecies,max(Nleaves)))


for (ispecies in seq(1,Nspecies)){

  print(ispecies)

  species_effect_N[ispecies,] <- param_all[,paste0("nu_species_N[",ispecies,"]")]
  species_effect_Cab[ispecies,] <- param_all[,paste0("nu_species_Cab[",ispecies,"]")]
  species_effect_Car[ispecies,] <- param_all[,paste0("nu_species_Car[",ispecies,"]")]
  species_effect_Cw[ispecies,] <- param_all[,paste0("nu_species_Cw[",ispecies,"]")]
  species_effect_Cm[ispecies,] <- param_all[,paste0("nu_species_Cm[",ispecies,"]")]

  for (ileaf in seq(1,Nleaves[ispecies])){

    leaf_effect_N[ispecies,ileaf,] <- param_all[,paste0("nu_leaf_N[",ispecies,", ",ileaf,"]")]
    leaf_effect_Cab[ispecies,ileaf,] <- param_all[,paste0("nu_leaf_Cab[",ispecies,", ",ileaf,"]")]
    leaf_effect_Car[ispecies,ileaf,] <- param_all[,paste0("nu_leaf_Car[",ispecies,", ",ileaf,"]")]
    leaf_effect_Cm[ispecies,ileaf,] <- param_all[,paste0("nu_leaf_Cm[",ispecies,", ",ileaf,"]")]
    leaf_effect_Cw[ispecies,ileaf,] <- param_all[,paste0("nu_leaf_Cw[",ispecies,", ",ileaf,"]")]

    cN <- param_all[,"Nmean"] + species_effect_N[ispecies,] + leaf_effect_N[ispecies,ileaf,]
    cCab <- param_all[,"Cabmean"] + species_effect_Cab[ispecies,] + leaf_effect_Cab[ispecies,ileaf,]
    cCar <- param_all[,"Carmean"] + species_effect_Car[ispecies,] + leaf_effect_Car[ispecies,ileaf,]
    cCw <- param_all[,"Cwmean"] + species_effect_Cw[ispecies,] + leaf_effect_Cw[ispecies,ileaf,]
    cCm <- param_all[,"Cmmean"] + species_effect_Cm[ispecies,] + leaf_effect_Cm[ispecies,ileaf,]

    all_N[ispecies,ileaf,] <- cN
    all_Cab[ispecies,ileaf,] <- cCab
    all_Car[ispecies,ileaf,] <- cCar
    all_Cw[ispecies,ileaf,] <- cCw
    all_Cm[ispecies,ileaf,] <- cCm

    tmp <- matrix(unlist(lapply(1:Nsimu,function(isimu){
      rrtm::prospect5(N = cN[isimu],
                      Cab = cCab[isimu],
                      Car = cCar[isimu],
                      Cbrown = 0,
                      Cw = cCw[isimu],
                      Cm = cCm[isimu])[["reflectance"]][pos]})),ncol = Nsimu)

    array_mod_reflectance[,ispecies,ileaf,] <- tmp

    X <- apply(tmp,c(1),mean)
    Y <- Data$obs_reflectance[,ispecies,ileaf]
    LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

    RMSE[ispecies,ileaf] <- sqrt(c(crossprod(LM$residuals))/length(LM$residuals))


    # ispecies = 2; ileaf = 8
    # matplot(WLs,array_mod_reflectance[,ispecies,ileaf,],type = 'l')
    # lines(WLs,data.minus1d[,ispecies,ileaf], col = "red",lwd = 2)

  }
}


all_parameters <- bind_rows(list(melt(all_N) %>% mutate(param = "N"),
                                 melt(all_Cab) %>% mutate(param = "Cab"),
                                 melt(all_Car) %>% mutate(param = "Car"),
                                 melt(all_Cw) %>% mutate(param = "Cw"),
                                 melt(all_Cm) %>% mutate(param = "Cm"))) %>% filter(!is.na(value)) %>% rename(species = Var1,
                                                                                                              leaf = Var2,
                                                                                                              simu = Var3)

ggplot(data = all_parameters) +
  geom_density(aes(x = value, fill = as.factor(species)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()



Mean_effects <- bind_rows(list(data.frame(param = "N",value = as.vector(as.matrix(param[,"Nmean"]))),
                               data.frame(param = "Cab",value = as.vector(as.matrix(param[,"Cabmean"]))),
                               data.frame(param = "Car",value = as.vector(as.matrix(param[,"Carmean"]))),
                               data.frame(param = "Cw",value = as.vector(as.matrix(param[,"Cwmean"]))),
                               data.frame(param = "Cm",value = as.vector(as.matrix(param[,"Cmmean"])))))

ggplot(data = Mean_effects) +
  geom_density(aes(x = value, fill = as.factor(param)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()


all_species_effects <- bind_rows(list(melt(species_effect_N) %>% mutate(param = "N"),
                                      melt(species_effect_Cab) %>% mutate(param = "Cab"),
                                      melt(species_effect_Car) %>% mutate(param = "Car"),
                                      melt(species_effect_Cm) %>% mutate(param = "Cm"),
                                      melt(species_effect_Cw) %>% mutate(param = "Cw"))) %>% filter(!is.na(value)) %>% rename(species = Var1,
                                                                                                                              simu = Var2)

ggplot(data = all_species_effects) +
  geom_density(aes(x = value, fill = as.factor(species)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

all_leaves_effects <- bind_rows(list(melt(leaf_effect_N) %>% mutate(param = "N"),
                                     melt(leaf_effect_Cab) %>% mutate(param = "Cab"),
                                     melt(leaf_effect_Car) %>% mutate(param = "Car"),
                                     melt(leaf_effect_Cm) %>% mutate(param = "Cm"),
                                     melt(leaf_effect_Cw) %>% mutate(param = "Cw"))) %>% filter(!is.na(value)) %>% rename(species = Var1,
                                                                                                                          leaf = Var2,
                                                                                                                          simu = Var3)

ggplot(data = all_leaves_effects) +
  geom_density(aes(x = value, fill = as.factor(leaf)), alpha = 0.4) +
  facet_wrap(param ~ species, scales = "free",ncol = 4) +
  geom_vline(xintercept = 0) +
  theme_bw()


plot(as.vector(Data$obs_reflectance[,2,]),as.vector(apply(array_mod_reflectance[,2,,],c(1,2),mean)))
abline(a = 0, b = 1, col ='red')

matplot(WLs,matrix(Data$obs_reflectance[,3,],nrow = Nwl),type = 'l', col = "black")
matlines(WLs,matrix(apply(array_mod_reflectance[,3,,],c(1,2),mean),nrow = Nwl), col = "red")

X <- as.vector(apply(array_mod_reflectance[,,,],c(1,2,3),mean,na.rm = TRUE))
Y <- as.vector(Data$obs_reflectance[,,])
plot(X,Y)

abline(a = 0, b = 1, col ='red')
LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

summary(LM)
sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

par(mfrow = c(1,2))
matplot(WLs,data.2d.NA,type = 'l',col = "black")
matlines(WLs,Data$obs_reflectance[,3,],type = 'l',col = "red")

matplot(WLs,matrix(apply(array_mod_reflectance,c(1,2,3),mean),nrow = Nwl),type = 'l',col = "black")
matlines(WLs,apply(array_mod_reflectance[,3,,],c(1,2),mean,na.rm = TRUE),type = 'l',col = "red")

par(mfrow = c(1,2))
matplot(WLs,data.2d.NA,type = 'l')
matplot(WLs,matrix(apply(array_mod_reflectance,c(1,2,3),mean),nrow = Nwl),type = 'l')

df.merged <- bind_rows(list(melt(Data$obs_reflectance) %>% mutate(type = 'obs'),
                            melt(apply(array_mod_reflectance,c(1,2,3),mean,na.rm = TRUE)) %>% mutate(type = 'mod'))) %>%
  pivot_wider(names_from = "type",
              values_from = "value")

ggplot(data = df.merged,
       aes(x = mod,y = obs,group = Var1, color = Var1)) +
  geom_point(alpha = 0.4, size = 0.1) +
  stat_smooth(se = FALSE,method = "lm", size = 0.1) +
  theme_bw()

par(mfrow=c(1,1))
df.r2 <- df.merged %>% group_by(Var1) %>% summarise(r2 = summary(lm(formula = obs ~ mod))[["r.squared"]])
plot(WLs,df.r2$r2)
abline(h = c(0.5,0.8,0.9), col = 'red')


ggplot(data = df.merged %>% filter(Var1 %in% (df.r2 %>% filter(r2 < 0.8) %>% pull(Var1))),
       aes(x = mod,y = obs,group = Var1, color = as.factor(Var1))) +
  geom_point() +
  stat_smooth(se = FALSE,method = "lm") +
  theme_bw()
print(WLs[df.r2 %>% filter(r2 < 0.8) %>% pull(Var1)])
