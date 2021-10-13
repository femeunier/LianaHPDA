rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggridges)

Nsimus <- 100

WLa <- 400
WLb <- 2500
Delta_WL <- 20

WLs <- seq(WLa,WLb,Delta_WL)
WLs <- WLs[(WLs > 450 & WLs < 680) | WLs > 800]
pos <- which((WLa:WLb %in% WLs))
Nwl <- length(pos)

GFs <- c("Liana","Tree")
sites <- c("PNM","FTS")

df.RMSE <- df.params <- data.frame()
leaf.count <- 0

for (iGF in seq(1,2)){
  for(isite in seq(1,2)){
    for (ispecies in seq(1,30)){

      data.raw <- LianaHPDA::array_obs_reflectance[pos,iGF,isite,ispecies,,]

      if (!all(is.na(data.raw))){

        print(c(iGF,isite,ispecies))

        dims <- dim(data.raw)
        data.2d <- matrix(data = data.raw,nrow = dims[1])
        data.2d.NA <- data.2d[,!is.na(data.2d[1,])]
        Nleaves <- ncol(data.2d.NA)

        OP.file <- file.path("/home/femeunier/Documents/projects/LianaHPDA/out",
                             paste(GFs[iGF],sites[isite],ispecies,sep = "."),
                             paste0("MCMC.single.species.GF",iGF,".site",isite,".species",ispecies,".RDS"))

        if (file.exists(OP.file)){

          param <- readRDS(OP.file)

          Nchains = length(param)
          Nsimu <- min(Nsimus,nrow(param[[1]]))
          pos.simu <- sample(1:nrow(param[[1]]),round(Nsimu/length(param)))
          param_all <- do.call(rbind,lapply(1:Nchains,function(i) param[[i]][pos.simu,]))

          param.names <- colnames(param_all)

          array_mod_reflectance <- array(data = NA,c(dim(data.2d.NA),Nsimu))

          all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <- array(data = NA, c(Nleaves,Nsimu))

          RMSE <- c()

          leaf_effect_N <- leaf_effect_Cab <- leaf_effect_Car <-
            leaf_effect_Cm <- leaf_effect_Cw <- array(data = NA, dim = c(Nleaves,Nsimu))

          for (ileaf in seq(1,Nleaves)){

            # print(ileaf/Nleaves)

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

            # plot(WLs,X,type = 'l')
            # lines(WLs,Y,col = "red")
            #
            # plot(X,Y)
            # abline(a = 0, b = 1, col ='red')

            RMSE[ileaf] <- sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

          }

          all_parameters <- bind_rows(list(melt(all_N) %>% mutate(param = "N"),
                                           melt(all_Cab) %>% mutate(param = "Cab"),
                                           melt(all_Car) %>% mutate(param = "Car"),
                                           melt(all_Cw) %>% mutate(param = "Cw"),
                                           melt(all_Cm) %>% mutate(param = "Cm"))) %>% filter(!is.na(value)) %>% rename(leaf = Var1,
                                                                                                                        simu = Var2)

          Mean_effects <- bind_rows(list(data.frame(param = "N",value = as.vector(as.matrix(param_all[,"Nmean"]))),
                                         data.frame(param = "Cab",value = as.vector(as.matrix(param_all[,"Cabmean"]))),
                                         data.frame(param = "Car",value = as.vector(as.matrix(param_all[,"Carmean"]))),
                                         data.frame(param = "Cw",value = as.vector(as.matrix(param_all[,"Cwmean"]))),
                                         data.frame(param = "Cm",value = as.vector(as.matrix(param_all[,"Cmmean"])))))


          all_leaves_effects <- bind_rows(list(melt(leaf_effect_N) %>% mutate(param = "N"),
                                               melt(leaf_effect_Cab) %>% mutate(param = "Cab"),
                                               melt(leaf_effect_Car) %>% mutate(param = "Car"),
                                               melt(leaf_effect_Cm) %>% mutate(param = "Cm"),
                                               melt(leaf_effect_Cw) %>% mutate(param = "Cw"))) %>% filter(!is.na(value)) %>% rename(leaf = Var1,
                                                                                                                                    simu = Var2)


          # X <- as.vector(apply(array_mod_reflectance,c(1,2),mean,na.rm = TRUE))
          # Y <- as.vector(data.2d.NA)
          # plot(X,Y)
          #
          # abline(a = 0, b = 1, col ='red')
          # LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)
          #
          # summary(LM)
          # sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

          mod <- melt(apply(array_mod_reflectance,c(1,2),mean,na.rm = TRUE)) %>% rename(wl = Var1,
                                                                                        leaf = Var2) %>%
            mutate(wl = WLs[wl])

          obs <- melt(data.2d.NA) %>% rename(wl = Var1,
                                             leaf = Var2) %>%
            mutate(wl = WLs[wl])

          df.species <- mod %>% rename(mod = value) %>% left_join(obs %>% rename(obs = value),
                                                                  by = c("wl","leaf"))


          # matplot(WLs,data.2d.NA,type = 'l')
          # plot(param$chain2[,c("Nmean","Cabmean","Cwmean","Cmmean")])
          # plot(param$chain1[,"nu_leaf_Cab[1]"])
          #
          # hist(as.vector(as.matrix(param$chain2[,which(grepl("nu_leaf_Cab",colnames(param$chain1)))[1:10]])))
          #
          # hist(RMSE)
          #
          # ggplot(data = all_parameters) +
          #   geom_density(aes(x = value), alpha = 0.4) +
          #   facet_wrap(~ param, scales = "free") +
          #   theme_bw()
          #
          # ggplot(data = Mean_effects) +
          #   geom_density(aes(x = value, fill = as.factor(param)), alpha = 0.4) +
          #   facet_wrap(~ param, scales = "free") +
          #   theme_bw()
          #
          # ggplot(data = all_leaves_effects) +
          #   geom_density(aes(x = value, fill = as.factor(leaf)), alpha = 0.4) +
          #   facet_wrap(~ param, scales = "free",ncol = 4) +
          #   geom_vline(xintercept = 0) +
          #   theme_bw()
          #
          # ggplot(data = all_leaves_effects) +
          #   geom_density(aes(x = value), alpha = 0.4) +
          #   facet_wrap(~ param, scales = "free",ncol = 4) +
          #   geom_vline(xintercept = 0) +
          #   theme_bw()
          #
          # ggplot(data = df.species,
          #        aes(x = mod, y = obs, color = as.factor(leaf))) +
          #   geom_point(alpha = 0.4) +
          #   stat_smooth(method = "lm", se = FALSE) +
          #   theme_bw()
          #
          ggplot(data = df.species,
                 aes(x = mod, y = obs, group = as.factor(wl))) +
            geom_point(alpha = 0.4) +
            stat_smooth(method = "lm", se = FALSE) +
            theme_bw()

          df.RMSE <- bind_rows(list(df.RMSE,
                                    data.frame(GF = GFs[iGF],
                                               site = sites[isite],
                                               species = ispecies,
                                               RMSE = RMSE,
                                               leaf = 1:length(RMSE),
                                               leaf.id = leaf.count + (1:length(RMSE)))))

          df.params <- bind_rows(list(df.params,
                                      all_parameters %>% mutate(data.frame(GF = GFs[iGF],
                                                                           site = sites[isite],
                                                                           species = ispecies,
                                                                           leaf.id = leaf.count + leaf))))

          leaf.count <- leaf.count + length(RMSE)

        } else {
          warning(paste0(OP.file," does not exist"))
        }
      }
    }
  }
}

leaf.selection <- df.RMSE %>% filter(RMSE < 0.015) %>% pull(leaf.id)

ggplot(data = df.RMSE) +
  geom_density(aes(x = RMSE)) +
  facet_wrap(GF ~ site, scales = "free") +
  theme_bw()

ggplot(data = df.params %>% filter(leaf.id %in% leaf.selection)) +
  geom_density_ridges(aes(x = value, y = site, fill = GF),
                      alpha = 0.3) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

df.params %>% filter(leaf.id %in% leaf.selection) %>% group_by(GF,site,param) %>% summarise(value.m = median(value)) %>% arrange(param,site,GF)



