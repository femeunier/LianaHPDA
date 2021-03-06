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

  for (ileaf in seq(1,Nleaves)){

    nu_leaf_N[ileaf] ~ dnorm(0, sd = sd_leaf_N)
    nu_leaf_Cab[ileaf] ~ dnorm(0, sd = sd_leaf_Cab)
    nu_leaf_Car[ileaf] ~ dnorm(0, sd = sd_leaf_Car)
    nu_leaf_Cw[ileaf] ~ dnorm(0, sd = sd_leaf_Cw)
    nu_leaf_Cm[ileaf] ~ dnorm(0, sd = sd_leaf_Cm)

    reflectance[1:Nwl,ileaf] <- NIMprospect5(Nmean + nu_leaf_N[ileaf],
                                             Cabmean + nu_leaf_Cab[ileaf],
                                             Carmean + nu_leaf_Car[ileaf],
                                             Cwmean + nu_leaf_Cw[ileaf],
                                             Cmmean + nu_leaf_Cm[ileaf],
                                             dataspec_p5[,], talf[],t12[],t21[], Nwl)

    for (j in 1:Nwl){
      obs_reflectance[j,ileaf] ~ dnorm(reflectance[j,ileaf], sd = max(0,Standard.Dev))
    }
  }

  Standard.Dev ~ dunif(0,1)

  Nmean ~ dunif(1.1,5)
  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)

  sd_leaf_N ~ dunif(0.0,0.5)
  sd_leaf_Cab ~ dunif(0.,30)
  sd_leaf_Car ~ dunif(0.,10)
  sd_leaf_Cw ~ dunif(0.,0.02)
  sd_leaf_Cm ~ dunif(0.,0.02)

})

WLa <- 400
WLb <- 2500
Delta_WL <- 20

WLs <- seq(WLa,WLb,Delta_WL)
WLs <- WLs[(WLs > 450 & WLs < 680) | WLs > 800]
pos <- which((WLa:WLb %in% WLs))
Nwl <- length(pos)

data.raw <- LianaHPDA::array_obs_reflectance[pos,,,,,]

data.mean <- data.raw[,1,1,1,,]
dims <- dim(data.mean)
data.2d <- matrix(data = data.mean,nrow = dims[1])
data.2d.NA <- data.2d[,!is.na(data.2d[1,])]

matplot(WLs,data.2d.NA,type = 'l')

Data <- list(obs_reflectance = data.2d.NA)

colnames(Data$obs_reflectance) <- NULL

Nleaves <- ncol(data.2d.NA)

Constants <- list(Nwl = Nwl,
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
                                         reflectance = c(Nwl,Nleaves)),
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
                                    "nu_leaf_N","nu_leaf_Cab","nu_leaf_Car","nu_leaf_Cw","nu_leaf_Cm",
                                    "sd_leaf_N","sd_leaf_Cab","sd_leaf_Car","sd_leaf_Cw","sd_leaf_Cm"),
                       data = Data,
                       inits = Inits,
                       nburnin = 50000,
                       nchains = Nchains,
                       niter = 100000,
                       thin = 50,
                       summary = TRUE,
                       WAIC = TRUE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples
param <- MCMCsamples[,]

Nsimu <- min(1000,nrow(MCMCsamples))
pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,]))[sample(1:(Nchains*Nsimu),Nsimu),]

param.names <- colnames(param_all)
plot(param[,c("Nmean","Cabmean","Cwmean","Cmmean")])
plot(param[,"nu_leaf_Cab[1]"])

hist(as.vector(as.matrix(param[,which(grepl("nu_leaf_Cab",colnames(param$chain1)))[1:10]])))

array_mod_reflectance <- array(data = NA,c(dim(data.2d.NA),Nsimu))

all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <- array(data = NA, c(Nleaves,Nsimu))


for (ileaf in seq(1,Nleaves)){

  print(ileaf)

  cN <- param_all[,"Nmean"] + param_all[,paste0("nu_leaf_N[",ileaf,"]")]
  cCab <- param_all[,"Cabmean"] + param_all[,paste0("nu_leaf_Cab[",ileaf,"]")]
  cCar <- param_all[,"Carmean"] + param_all[,paste0("nu_leaf_Car[",ileaf,"]")]
  cCw <- param_all[,"Cwmean"] + param_all[,paste0("nu_leaf_Cw[",ileaf,"]")]
  cCm <- param_all[,"Cmmean"] + param_all[,paste0("nu_leaf_Cm[",ileaf,"]")]

  all_N[ileaf,] <- cN
  all_Cab[ileaf,] <- cCab
  all_Car[ileaf,] <- cCar
  all_Cw[ileaf,] <- cCw
  all_Cm[ileaf,] <- cCm

  tmp <- matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = cN[ileaf],
                    Cab = cCab[ileaf],
                    Car = cCar[ileaf],
                    Cbrown = 0,
                    Cw = cCw[ileaf],
                    Cm = cCm[ileaf])[["reflectance"]][pos]})),ncol = Nsimu)

  array_mod_reflectance[,ileaf,] <- tmp

}

all_parameters <- bind_rows(list(melt(all_N) %>% mutate(param = "N"),
                                 melt(all_Cab) %>% mutate(param = "Cab"),
                                 melt(all_Car) %>% mutate(param = "Car"),
                                 melt(all_Cw) %>% mutate(param = "Cw"),
                                 melt(all_Cm) %>% mutate(param = "Cm"))) %>% filter(!is.na(value)) %>% rename(leaf = Var1,
                                                                                                              simu = Var2)

ggplot(data = all_parameters) +
  geom_density(aes(x = value), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

X <- as.vector(apply(array_mod_reflectance,c(1,2),mean,na.rm = TRUE))
Y <- as.vector(data.2d.NA)
plot(X,Y)

abline(a = 0, b = 1, col ='red')
LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

summary(LM)
sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

par(mfrow = c(1,2))
matplot(WLs,data.2d.NA,type = 'l')
matplot(WLs,apply(array_mod_reflectance,c(1,2),mean,na.rm = TRUE),type = 'l')

df.merged <- bind_rows(list(melt(data.2d.NA) %>% mutate(type = 'obs'),
                            melt(apply(array_mod_reflectance,c(1,2),mean,na.rm = TRUE)) %>% mutate(type = 'mod'))) %>%
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

