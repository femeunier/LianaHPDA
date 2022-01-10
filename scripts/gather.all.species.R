rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggridges)
library(RColorBrewer)

Nsimus <- 100

WLa <- 400
WLb <- 2500
Delta_WL <- 20

all.WLs <- WLa:WLb

WLs <- seq(WLa,WLb,Delta_WL)
WLs <- WLs[(WLs > 450 & WLs < 680) | WLs > 800]

pos <- which((WLa:WLb %in% WLs))
Nwl <- length(pos)

GFs <- c('Tree','Liana')
sites <- c('PNM','FTS')

basename <- "MCMC.single.species.GF"
basename <- "Rcorrected"

df.RMSE <- df.params <- data.frame()
leaf.count <- 0

array_obs_reflectance_all <- LianaHPDA::array_obs_reflectance
array_mod_reflectance_all <- array(data = NA, dim = c(dim(array_obs_reflectance_all), Nsimus))

Nleaf.per.ind <- 3

for (iGF in seq(1,2)){
  for(isite in seq(1,2)){
    for (ispecies in seq(1,30)){

      data.raw <- array_obs_reflectance_all[pos,iGF,isite,ispecies,,]

      if (!all(is.na(data.raw))){

        print(c(iGF,isite,ispecies))

        dims <- dim(data.raw)
        data.2d <- matrix(data = data.raw,nrow = dims[1])
        data.2d.NA <- data.2d[,!is.na(data.2d[1,])]

        Nleaves <- ncol(data.2d.NA)

        Nactual.ind <- Nleaves/3

        OP.file <- file.path("/home/femeunier/Documents/projects/LianaHPDA/out",
                             paste(GFs[iGF],sites[isite],ispecies,sep = "."),
                             paste0(basename,".GF",iGF,".site",isite,".species",ispecies,".RDS"))

        if (file.exists(OP.file)){

          param <- readRDS(OP.file)

          Nchains = length(param)
          Nsimu <- min(Nsimus,nrow(param[[1]]))
          pos.simu <- sample(1:nrow(param[[1]]),round(Nsimu/length(param)))
          param_all <- do.call(rbind,lapply(1:Nchains,function(i) param[[i]][pos.simu,]))

          param.names <- colnames(param_all)

          array_mod_reflectance <- array(data = NA,c(dim(data.2d.NA),Nsimu))

          all_N <- all_Cab <- all_Car <- all_Cw <- all_Cm <-
            all_mND705 <- all_mSR705 <- array(data = NA, c(Nleaves,Nsimu))

          RMSE <- c()

          leaf_effect_N <- leaf_effect_Cab <- leaf_effect_Car <-
            leaf_effect_Cm <- leaf_effect_Cw <- all_redgedge <-
            array(data = NA, dim = c(Nleaves,Nsimu))

          compt.ind <- 1

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
                              Cm = cCm[ileaf])[["reflectance"]]})),ncol = Nsimu)

            current.spectra <- tmp
            tmp <- tmp[pos,]

            array_mod_reflectance[,ileaf,] <- tmp

            R750 <- current.spectra[which.min(abs(all.WLs-750)),]
            R705 <- current.spectra[which.min(abs(all.WLs-705)),]
            R445 <- current.spectra[which.min(abs(all.WLs-445)),]
            mND705 <- (R750 - R705)/(R750 + R705 - 2*R445)
            mSR705 <- (R750 - R445)/(R705 - R445)
            rededge <- colMeans(current.spectra[all.WLs %in% seq(720,750),])

            all_mND705[ileaf,] <- mND705
            all_mSR705[ileaf,] <- mSR705
            all_redgedge[ileaf,] <- rededge

            X <- as.vector(apply(tmp,c(1),mean))
            Y <- as.vector(data.2d.NA[,ileaf])
            LM <- lm(data.frame(x = X,y = Y),formula = y ~ x)

            # df.all <- bind_rows(list(df.all,
            #                          data.frame(wv = WLs,
            #                                     mod = X,
            #                                     obs = Y,
            #                                     ind = compt.ind,
            #                                     leaf = ceiling(ileaf/Nactual.ind),
            #                                     GF = GFs[iGF],
            #                                     site = sites[isite],
            #                                     species = ispecies)))
            # plot(WLs,X,type = 'l')
            # lines(WLs,Y,col = "red")
            #
            # plot(X,Y)
            # abline(a = 0, b = 1, col ='red')

            RMSE[ileaf] <- sqrt(c(crossprod(LM$residuals))/length(LM$residuals))

            array_mod_reflectance_all[pos,iGF,isite,ispecies,compt.ind,ceiling(ileaf/Nactual.ind),] <- tmp

            compt.ind <- compt.ind + 1
            if (compt.ind > Nactual.ind){
              compt.ind <- 1
            }

          }

          all_parameters <- bind_rows(list(melt(all_N) %>% mutate(param = "N"),
                                           melt(all_Cab) %>% mutate(param = "Cab"),
                                           melt(all_Car) %>% mutate(param = "Car"),
                                           melt(all_Cw) %>% mutate(param = "Cw"),
                                           melt(all_Cm) %>% mutate(param = "Cm"),
                                           melt(all_mND705) %>% mutate(param = "mND705"),
                                           melt(all_mSR705) %>% mutate(param = "mSR705"),
                                           melt(all_redgedge) %>% mutate(param = "red.edge"))) %>% filter(!is.na(value)) %>% rename(leaf = Var1,
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


          # matplot(WLs,data.2d.NA[,(3*(4-1)) + 1:3],type = 'l')
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
          # ggplot(data = df.species,
          #        aes(x = mod, y = obs, group = as.factor(wl))) +
          #   geom_point(alpha = 0.4) +
          #   stat_smooth(method = "lm", se = FALSE) +
          #   theme_bw()

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
          # stop()
          warning(paste0(OP.file," does not exist"))
        }
      }
    }
  }
}

stop()

df.mod <- melt(apply(array_mod_reflectance_all[pos,,,,,,],c(1,2,3,4,5,6),mean,na.rm = TRUE)) %>%
  rename(wv = Var1,
         GF = Var2,
         site = Var3,
         species = Var4,
         ind = Var5,
         leaf = Var6,
         mod = value) %>%
  mutate(wv = WLs[wv],
         GF = GFs[GF],
         site = sites[site]) %>%
  filter(!is.na(mod))

df.obs <- melt(array_obs_reflectance_all[pos,,,,,]) %>%
  rename(wv = Var1,
         GF = Var2,
         site = Var3,
         species = Var4,
         ind = Var5,
         leaf = Var6,
         obs = value) %>%
  filter(!is.na(obs))

df.all <- df.mod %>% left_join(df.obs,
                               by = c("wv","GF","site","species","ind","leaf"))

ggplot(data = df.all,
       aes(x = mod,y = obs)) +
  geom_hex(bins = 200) +
  scale_fill_distiller(palette="OrRd",trans = "reverse") +
  geom_abline(slope = 1, intercept = 0, color = "black",linetype = 3) +
  theme_bw()

ggplot(data = df.all,
       aes(x = mod,y = obs)) +
  geom_point(aes(color = as.factor(species))) +
  scale_fill_distiller(palette="OrRd",trans = "reverse") +
  geom_abline(slope = 1, intercept = 0, color = "black",linetype = 3) +
  facet_wrap(GF ~ site) +
  theme_bw()


# df.all.all <- df.all %>% left_join(df.all2 %>% rename(mod2 = mod,obs2 = obs),
#                                    by = c("wv","GF","site","species","ind","leaf")) %>% mutate(diff.obs = obs2 - obs,
#                                                                                          diff.mod = mod2 - mod)




#############################################################################################################################

df.mod.species <- bind_cols(list(
  melt(apply(
    array_mod_reflectance_all[pos, , , , , , ], c(1, 2, 3, 4), mean, na.rm = TRUE
  )) %>% rename(
    wv = Var1,
    GF = Var2,
    site = Var3,
    species = Var4,
    mod = value
  ) %>%
    mutate(wv = WLs[wv],
           GF = GFs[GF],
           site = sites[site]),
  melt(apply(
    array_mod_reflectance_all[pos, , , , , , ], c(1, 2, 3, 4), sd, na.rm = TRUE)) %>% dplyr::select(value) %>% rename(mod.sd = value))) %>%
  filter(!is.na(mod))


df.obs.species <- bind_cols(list(
  melt(apply(array_obs_reflectance_all[pos,,,,,],c(1,2,3,4),mean,na.rm = TRUE)) %>% rename(wv = Var1,
                                                                                           GF = Var2,
                                                                                           site = Var3,
                                                                                           species = Var4,
                                                                                           obs = value),
  melt(apply(
    array_obs_reflectance_all[pos, , , , , ], c(1, 2, 3, 4), sd, na.rm = TRUE)) %>% dplyr::select(value) %>% rename(obs.sd = value))) %>%
  filter(!is.na(obs))

df.all.species <- df.mod.species %>% left_join(df.obs.species,
                                               by = c("wv","GF","site","species"))


ggplot(data = df.all.species,
       aes(x = mod, y = obs, color = as.factor(species))) +
  geom_point(size = 0.1) +
  geom_errorbar(aes(ymin = obs - obs.sd, ymax = obs + obs.sd), show.legend = FALSE) +
  geom_errorbarh(aes(xmin = mod - mod.sd, xmax = mod + mod.sd), show.legend = FALSE) +
  scale_fill_distiller(palette="OrRd",trans = "reverse") +
  geom_abline(slope = 1, intercept = 0, color = "black",linetype = 3) +
  facet_grid(GF ~ site) +
  theme_bw()

df.all.species %>% mutate(SE = (mod - obs)**2) %>% group_by(GF,site,species) %>% summarise(RMSE = sqrt(sum(SE)/length(SE))) %>% arrange(desc(RMSE))


ggplot() +
  geom_line(data = df.all.species %>% filter(GF == "Liana",site == "FTS",species == 1),
            aes(x = wv, y = mod), color = "black") +
  geom_line(data = df.all.species %>% filter(GF == "Liana",site == "FTS",species == 1),
            aes(x = wv, y = obs), color = "red") +
  theme_bw()




df.RMSE.all <- bind_rows(list(df.RMSE,
                              df.RMSE %>% mutate(GF = "All", site = "")))

ggplot(data = df.RMSE.all) +
  geom_density_ridges(aes(x = RMSE, y = interaction(GF,site), fill = interaction(GF,site))) +
  scale_x_continuous(limits = c(0,1.1*max(df.RMSE$RMSE))) +
  theme_bw()

leaf.selection <- df.RMSE %>% filter(RMSE < Inf) %>% pull(leaf.id)
df.params$GF <- factor(df.params$GF,levels = c("Liana","Tree"))


df.params.wide.long <- df.params %>% pivot_wider(values_from = "value",
                                                 names_from = "param") %>% mutate(WC =  (1-(1/(1+Cw/Cm)))*100) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","WC","mND705","mSR705","red.edge"),
               values_to = "value",
               names_to = "param")

df.params.mod <- df.params.wide.long %>% mutate(value = case_when(param == "Cm" ~ 1/(10*value),
                                                                  TRUE ~ value),
                                      param = case_when(param == "Cm" ~ "SLA",
                                                        TRUE ~ param))

df.params.mod %>% filter(leaf.id %in% leaf.selection) %>%
  group_by(GF,site,param) %>%
  summarise(value.m = median(value),
            .groups = "keep") %>% arrange(param,site,GF)

ggplot(data = df.params.mod %>%
         filter(param %in% c("N","Cab","Car","SLA","WC")) %>%
         filter(param != "SLA" | (param == "SLA" & value <= 30))) +
  geom_density_ridges(aes(x = value, y = site, fill = GF),
                      alpha = 0.3) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()

ggplot(data = df.params.mod %>% filter(leaf.id %in% leaf.selection,
                                       param %in% c("N","Cab","Car","SLA","WC")) %>%
         filter(param != "SLA" | (param == "SLA" & value <= 50)),
       aes(x = site, fill = GF,y = value)) +
  geom_violin(trim=FALSE, position = position_dodge(0.9))+
  geom_boxplot(aes(group = interaction(site,GF)),position = position_dodge(0.9),fill = NA, width = 0.1,
               outlier.shape=NA) +
  facet_wrap(~ param, scales = "free", nrow = 1) +
  theme_bw() +
  labs(x = "", y = "", fill = "Growth form") +
  theme(legend.position = c(0.14,0.85))

df.params.wide <- df.params.mod %>% pivot_wider(values_from = "value",
                                            names_from = "param")

df.params.wide.species <- df.params.wide %>% group_by(GF,site,species) %>%
  summarise(N.mean = mean(N),
            Cab.mean = mean(Cab),
            Car.mean = mean(Car),
            SLA.mean = mean(SLA),
            Cw.mean = mean(Cw),
            N.mean = mean(N),
            mND705.mean = mean(mND705),
            mSR705.mean = mean(mSR705),
            red.edge.mean = mean(red.edge),
            .groups = "keep")

df.params.wide.species %>% group_by(GF,site) %>% summarise(Nspecies = length(N.mean),
                                                           .groups = "keep")

ggplot(data = df.params.wide.species,
       aes(x = Cab.mean, y =  mND705.mean, color = GF)) +
  geom_point(alpha = 0.5) +
  stat_smooth(se = TRUE,method = "lm", alpha = 0.1) +
  facet_wrap(~ site) +
  theme_bw()

ggplot(data = df.params.wide.species,
       aes(x = Cab.mean, y =  mSR705.mean, color = GF,fill = GF)) +
  geom_point(alpha = 0.5) +
  stat_smooth(se = TRUE,method = "lm", alpha = 0.1) +
  facet_wrap(~ site) +
  theme_bw()

ggplot(data = df.params.wide.species,
       aes(x = Cab.mean, y =  red.edge.mean, color = GF,fill = GF)) +
  geom_point(alpha = 0.5) +
  stat_smooth(se = TRUE,method = "lm", alpha = 0.1) +
  facet_wrap(~ site, nrow = 2) +
  theme_bw()



