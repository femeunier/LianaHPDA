rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

WLa <- 400
WLb <- 2500
Delta_WL <- 20

WLs <- seq(WLa,WLb,Delta_WL)
WLs <- WLs[(WLs > 450 & WLs < 680) | WLs > 800]
pos <- which((WLa:WLb %in% WLs))
Nwl <- length(pos)

GFs <- c("Liana","Tree")
sites <- c("PNM","FTS")

ispecies = 3
isite = 1
iGF = 2

data.raw <- LianaHPDA::array_obs_reflectance[pos,iGF,isite,ispecies,,]
dims <- dim(data.raw)
data.2d <- matrix(data = data.raw,nrow = dims[1])
data.2d.NA <- data.2d[,!is.na(data.2d[1,])]
Nleaves <- ncol(data.2d.NA)

matplot(WLs,data.2d.NA,type = 'l')

param <- readRDS(file.path("/home/femeunier/Documents/projects/LianaHPDA/out",
                          paste(GFs[iGF],sites[isite],ispecies,sep = "."),
                          paste0("MCMC.single.species.GF",iGF,".site",isite,".species",ispecies,".RDS")))

# param <- readRDS(file.path("/home/femeunier/Documents/projects/LianaHPDA/out/Liana.PNM.1/MCMC.single.species.GF1.site1.species1.RDS"))

Nchains = length(param)
Nsimu <- min(1000,nrow(param[[1]]))
pos.simu <- sample(1:nrow(param[[1]]),round(Nsimu/length(param)))
param_all <- do.call(rbind,lapply(1:Nchains,function(i) param[[i]][pos.simu,]))

param.names <- colnames(param_all)
plot(param$chain2[,c("Nmean","Cabmean","Cwmean","Cmmean")])
plot(param$chain1[,"nu_leaf_Cab[1]"])

hist(as.vector(as.matrix(param$chain2[,which(grepl("nu_leaf_Cab",colnames(param$chain1)))[1:10]])))

array_mod_reflectance <- array(data = NA,c(dim(data.2d.NA),Nsimu))

all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <- array(data = NA, c(Nleaves,Nsimu))

RMSE <- array(data = NA, dim = c(Nleaves))

leaf_effect_N <- leaf_effect_Cab <- leaf_effect_Car <-
  leaf_effect_Cm <- leaf_effect_Cw <- array(data = NA, dim = c(Nleaves,Nsimu))

for (ileaf in seq(1,Nleaves)){

  print(ileaf/Nleaves)

  leaf_effect_N[ileaf,] <- param_all[,paste0("nu_leaf_N[",ileaf,"]")]
  leaf_effect_Cab[ileaf,] <- param_all[,paste0("nu_leaf_Cab[",ileaf,"]")]
  leaf_effect_Car[ileaf,] <- param_all[,paste0("nu_leaf_Car[",ileaf,"]")]
  leaf_effect_Cm[ileaf,] <- param_all[,paste0("nu_leaf_Cm[",ileaf,"]")]
  leaf_effect_Cw[ileaf,] <- param_all[,paste0("nu_leaf_Cw[",ileaf,"]")]

  cN <- param_all[,"Nmean"] + leaf_effect_N[ileaf,]
  cCab <- param_all[,"Cabmean"] + leaf_effect_Cab[ileaf,]
  cCar <- param_all[,"Carmean"] + leaf_effect_Car[ileaf,]
  cCw <- param_all[,"Cwmean"] + leaf_effect_Cw[ileaf,]
  cCm <- param_all[,"Cmmean"] + leaf_effect_Cm[ileaf,]

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

  X <- as.vector(apply(tmp,c(1),mean))
  Y <- as.vector(data.2d.NA[,ileaf])
  LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

  # plot(WLs,X)
  # lines(WLs,Y,col = "red")
  # plot(X,Y)
  # abline(a = 0, b = 1, col ='red')

  RMSE[ileaf] <- sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

}

hist(RMSE)

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

Mean_effects <- bind_rows(list(data.frame(param = "N",value = as.vector(as.matrix(param_all[,"Nmean"]))),
                               data.frame(param = "Cab",value = as.vector(as.matrix(param_all[,"Cabmean"]))),
                               data.frame(param = "Car",value = as.vector(as.matrix(param_all[,"Carmean"]))),
                               data.frame(param = "Cw",value = as.vector(as.matrix(param_all[,"Cwmean"]))),
                               data.frame(param = "Cm",value = as.vector(as.matrix(param_all[,"Cmmean"])))))

ggplot(data = Mean_effects) +
  geom_density(aes(x = value, fill = as.factor(param)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

all_leaves_effects <- bind_rows(list(melt(leaf_effect_N) %>% mutate(param = "N"),
                                     melt(leaf_effect_Cab) %>% mutate(param = "Cab"),
                                     melt(leaf_effect_Car) %>% mutate(param = "Car"),
                                     melt(leaf_effect_Cm) %>% mutate(param = "Cm"),
                                     melt(leaf_effect_Cw) %>% mutate(param = "Cw"))) %>% filter(!is.na(value)) %>% rename(leaf = Var1,
                                                                                                                          simu = Var2)

ggplot(data = all_leaves_effects) +
  geom_density(aes(x = value, fill = as.factor(leaf)), alpha = 0.4) +
  facet_wrap(~ param, scales = "free",ncol = 4) +
  geom_vline(xintercept = 0) +
  theme_bw()

ggplot(data = all_leaves_effects) +
  geom_density(aes(x = value), alpha = 0.4) +
  facet_wrap(~ param, scales = "free",ncol = 4) +
  geom_vline(xintercept = 0) +
  theme_bw()

X <- as.vector(apply(array_mod_reflectance,c(1,2),mean,na.rm = TRUE))
Y <- as.vector(data.2d.NA)
plot(X,Y)

abline(a = 0, b = 1, col ='red')
LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

summary(LM)
sqrt(c(crossprod(LM$residuals))/length(LM$residuals))
