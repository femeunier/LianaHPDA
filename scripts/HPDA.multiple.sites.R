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


  sd_site_N ~ dunif(0.0,0.5)
  sd_site_Cab ~ dunif(0.,30)
  sd_site_Car ~ dunif(0.,10)
  sd_site_Cw ~ dunif(0.,0.02)
  sd_site_Cm ~ dunif(0.,0.02)

  for (isite in seq(1,Nsites)){

    # Site fixed effect

    nu_site_N[isite] ~ dnorm(0, sd = sd_site_N)
    nu_site_Cab[isite] ~ dnorm(0, sd = sd_site_Cab)
    nu_site_Car[isite] ~ dnorm(0, sd = sd_site_Car)
    nu_site_Cw[isite] ~ dnorm(0, sd = sd_site_Cw)
    nu_site_Cm[isite] ~ dnorm(0, sd = sd_site_Cm)

    sd_species_N[isite] ~ dunif(0.0,0.5)
    sd_species_Cab[isite] ~ dunif(0.,30)
    sd_species_Car[isite] ~ dunif(0.,10)
    sd_species_Cw[isite] ~ dunif(0.,0.02)
    sd_species_Cm[isite] ~ dunif(0.,0.02)

    for (ispecies in seq(1,Nspecies[isite])){

      # Species fixed effect

      nu_species_N[isite,ispecies] ~ dnorm(0, sd = sd_species_N[isite])
      nu_species_Cab[isite,ispecies] ~ dnorm(0, sd = sd_species_Cab[isite])
      nu_species_Car[isite,ispecies] ~ dnorm(0, sd = sd_species_Car[isite])
      nu_species_Cw[isite,ispecies] ~ dnorm(0, sd = sd_species_Cw[isite])
      nu_species_Cm[isite,ispecies] ~ dnorm(0, sd = sd_species_Cm[isite])

      sd_leaf_N[isite,ispecies] ~ dunif(0.0,0.5)
      sd_leaf_Cab[isite,ispecies] ~ dunif(0.,30)
      sd_leaf_Car[isite,ispecies] ~ dunif(0.,10)
      sd_leaf_Cw[isite,ispecies] ~ dunif(0.,0.02)
      sd_leaf_Cm[isite,ispecies] ~ dunif(0.,0.02)


      for (ileaf in seq(1,Nleaves[isite,ispecies])){

        # leaf random effect
        nu_leaf_N[isite,ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_N[isite,ispecies])
        nu_leaf_Cab[isite,ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Cab[isite,ispecies])
        nu_leaf_Car[isite,ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Car[isite,ispecies])
        nu_leaf_Cw[isite,ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Cw[isite,ispecies])
        nu_leaf_Cm[isite,ispecies,ileaf] ~ dnorm(0, sd = sd_leaf_Cm[isite,ispecies])

        reflectance[1:Nwl,isite,ispecies,ileaf] <- NIMprospect5(Nmean +   nu_site_N[isite] +   nu_species_N[isite,ispecies] +   nu_leaf_N[isite,ispecies,ileaf],
                                                                Cabmean + nu_site_Cab[isite] + nu_species_Cab[isite,ispecies] + nu_leaf_Cab[isite,ispecies,ileaf],
                                                                Carmean + nu_site_Car[isite] + nu_species_Car[isite,ispecies] + nu_leaf_Car[isite,ispecies,ileaf],
                                                                Cwmean +  nu_site_Cw[isite] +  nu_species_Cw[isite,ispecies] +  nu_leaf_Cw[isite,ispecies,ileaf],
                                                                Cmmean +  nu_site_Cm[isite] +  nu_species_Cm[isite,ispecies] +  nu_leaf_Cm[isite,ispecies,ileaf],
                                                          dataspec_p5[,], talf[],t12[],t21[], Nwl)

        for (j in 1:Nwl){
          obs_reflectance[j,isite,ispecies,ileaf] ~ dnorm(reflectance[j,isite,ispecies,ileaf], sd = max(0,Standard.Dev))
        }
      }

    }
  }

  Standard.Dev ~ dunif(0,1)

  Nmean ~ dunif(1.1,5)
  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)

})

WLa <- 400
WLb <- 2500
Delta_WL <- 20

WLs <- seq(WLa,WLb,Delta_WL)
WLs <- WLs[WLs < 680 | WLs > 800]
pos <- which((WLa:WLb %in% WLs))
Nwl <- length(pos)

data.raw <- LianaHPDA::array_obs_reflectance[pos,,,,,]

data.mean <- data.raw[,1,1:2,1:4,1:4,1:3]
dims <- dim(data.mean)
data.minus1d <- array(data = data.mean, dim = c(dims[1:(length(dims)-2)],dims[length(dims)-1]*dims[length(dims)]))
data.2d.NA <- matrix(data.mean,nrow = Nwl)

par(mfrow = c(1,1))
matplot(WLs,data.2d.NA,type = 'l')

Data <- list(obs_reflectance = data.minus1d)

colnames(Data$obs_reflectance) <- NULL

Nsites <- dims[2]
Nspecies <- rep(dims[3],Nsites)
Nleaves <- array(dims[4]*dims[5],dim = c(dims[2],dims[3]))

Constants <- list(Nwl = Nwl,
                  Nsites = Nsites,
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
                                         reflectance = c(Nwl,max(Nsites),max(Nspecies),max(Nleaves))),
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
                                    "nu_site_N","nu_site_Cab","nu_site_Car","nu_site_Cw","nu_site_Cm",
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

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,]

Nsimu <- min(1000,nrow(MCMCsamples[[1]]))
pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,]))[sample(1:(Nchains*Nsimu),Nsimu),]

plot(param[,c("Nmean","Cabmean","Cwmean","Cmmean")])
plot(param[,"nu_species_Cab[2, 1]"])

hist(as.vector(as.matrix(param[,which(grepl("nu_site_N",colnames(param$chain1)))])))
hist(as.vector(as.matrix(param[,which(grepl("nu_species_N",colnames(param$chain1)))])))
hist(as.vector(as.matrix(param[,which(grepl("nu_leaf_Cab",colnames(param$chain1)))])))

array_mod_reflectance <- array(data = NA,c(dim(data.minus1d),Nsimu))

all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <-
  array(data = NA, c(max(Nsites),max(Nspecies),max(Nleaves),Nsimu))

for (isite in seq(1,Nsites)){
  for (ispecies in seq(1,Nspecies[isite])){

    print(c(isite,ispecies))

    for (ileaf in seq(1,Nleaves[isite,ispecies])){

      cN <- param_all[,"Nmean"] +
        param_all[,paste0("nu_site_N[",isite,"]")] +
        param_all[,paste0("nu_species_N[",isite,", ",ispecies,"]")] +
        param_all[,paste0("nu_leaf_N[",isite,", ",ispecies,", ",ileaf,"]")]

      cCab <- param_all[,"Cabmean"] +
        param_all[,paste0("nu_site_Cab[",isite,"]")] +
        param_all[,paste0("nu_species_Cab[",isite,", ",ispecies,"]")] +
        param_all[,paste0("nu_leaf_Cab[",isite,", ",ispecies,", ",ileaf,"]")]

      cCar <- param_all[,"Carmean"] +
        param_all[,paste0("nu_site_Car[",isite,"]")] +
        param_all[,paste0("nu_species_Car[",isite,", ",ispecies,"]")] +
        param_all[,paste0("nu_leaf_Car[",isite,", ",ispecies,", ",ileaf,"]")]

      cCw <- param_all[,"Cwmean"] +
        param_all[,paste0("nu_site_Cw[",isite,"]")] +
        param_all[,paste0("nu_species_Cw[",isite,", ",ispecies,"]")] +
        param_all[,paste0("nu_leaf_Cw[",isite,", ",ispecies,", ",ileaf,"]")]

      cCm <- param_all[,"Cmmean"] +
        param_all[,paste0("nu_site_Cm[",isite,"]")] +
        param_all[,paste0("nu_species_Cm[",isite,", ",ispecies,"]")] +
        param_all[,paste0("nu_leaf_Cm[",isite,", ",ispecies,", ",ileaf,"]")]

      all_N[isite,ispecies,ileaf,] <- cN
      all_Cab[isite,ispecies,ileaf,] <- cCab
      all_Car[isite,ispecies,ileaf,] <- cCar
      all_Cw[isite,ispecies,ileaf,] <- cCw
      all_Cm[isite,ispecies,ileaf,] <- cCm

      tmp <- matrix(unlist(lapply(1:Nsimu,function(isimu){
        rrtm::prospect5(N = cN[isimu],
                        Cab = cCab[isimu],
                        Car = cCar[isimu],
                        Cbrown = 0,
                        Cw = cCw[isimu],
                        Cm = cCm[isimu])[["reflectance"]][pos]})),ncol = Nsimu)

      array_mod_reflectance[,isite,ispecies,ileaf,] <- tmp

      # matplot(WLs,array_mod_reflectance[,ispecies,ileaf,],type = 'l')
      # lines(WLs,data.minus1d[,ispecies,ileaf], col = "red",lwd = 2)

    }
  }
}

all_parameters <- bind_rows(list(melt(all_N) %>% mutate(param = "N"),
                                 melt(all_Cab) %>% mutate(param = "Cab"),
                                 melt(all_Car) %>% mutate(param = "Car"),
                                 melt(all_Cw) %>% mutate(param = "Cw"),
                                 melt(all_Cm) %>% mutate(param = "Cm"))) %>% filter(!is.na(value)) %>% rename(site = Var1,
                                                                                                              species = Var2,
                                                                                                              leaf = Var3,
                                                                                                              simu = Var4)

ggplot(data = all_parameters) +
  geom_density_ridges(aes(x = value, y = as.factor(site), fill = as.factor(species)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

X <- as.vector(apply(array_mod_reflectance,c(1,2,3,4),mean,na.rm = TRUE))
Y <- as.vector(data.minus1d)
plot(X,Y)

abline(a = 0, b = 1, col ='red')
LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

summary(LM)
sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

par(mfrow = c(1,2))
matplot(WLs,data.2d.NA,type = 'l')
matplot(WLs,matrix(apply(array_mod_reflectance,c(1,2,3,4),mean),nrow = Nwl),type = 'l')

df.merged <- bind_rows(list(melt(Data$obs_reflectance) %>% mutate(type = 'obs'),
                            melt(apply(array_mod_reflectance,c(1,2,3,4),mean,na.rm = TRUE)) %>% mutate(type = 'mod'))) %>%
  pivot_wider(names_from = "type",
              values_from = "value")
ggplot(data = df.merged,
       aes(x = mod,y = obs,group = Var1)) +
  geom_point() +
  stat_smooth(se = FALSE,method = "lm") +
  theme_bw()

par(mfrow=c(1,1))
df.r2 <- df.merged %>% group_by(Var1) %>% summarise(r2 = summary(lm(formula = obs ~ mod))[["r.squared"]])
plot(WLs,df.r2$r2)


ggplot(data = df.merged %>% filter(Var1 %in% (df.r2 %>% filter(r2 < 0.8) %>% pull(Var1))),
       aes(x = mod,y = obs,group = Var1)) +
  geom_point() +
  stat_smooth(se = FALSE,method = "lm") +
  theme_bw()

