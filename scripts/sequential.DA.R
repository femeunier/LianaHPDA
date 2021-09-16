rm(list = ls())

library(pracma)
library(dplyr)
library(ggplot2)
library(tidyr)
library(BayesianTools)
library(purrr)
library(ggridges)
library(cowplot)
library(ggpointdensity)
library(Rprospect)

wv.min = 400
wv.max = 2500

Nprospect <- 1000

directory <- "/home/femeunier/Documents/projects/LianaHPDA/data/Guzman_all"
metadata.file <- file.path(directory,"meta_data_all2.csv")
metadata <- read.csv(metadata.file) %>% mutate(GF = case_when(substr(Code,1,1) == "L" ~ "Liana",
                                                              substr(Code,1,1) == "T" ~ "Tree"),
                                               Species = (substring(sub("\\_.*", "", Code),2)),
                                               Ind = as.numeric(sub(".*\\_", "", Code)),
                                               N  = 1 + To.ASD.spectra - From.ASD.spectra,
                                               data.exist = FALSE,
                                               site = case_when(substr(Site.date,1,1) == "P" ~ "PNM",
                                                                substr(Site.date,1,1) == "F" ~ "FTS"))


files <- paste0(as.character(unique(metadata$Site.date)),".csv")

data.files <- list()

for (ifile in seq(files)){
  data.files[[tools::file_path_sans_ext(files[ifile])]] <- read.csv(file.path(directory,files[ifile]))
}

All.cols <- unique(c(colnames(data.files[[1]]),colnames(data.files[[2]]),colnames(data.files[[3]])))

data.all <- list("Liana" = list(),"Tree" = list())

df.data.all <- data.frame()

for (i in seq(1,nrow(metadata))){

  print(i/nrow(metadata))

  temp <- metadata[i,]
  wv <- data.files[[as.character(temp[["Site.date"]])]][["Wavelength"]]
  columns <- paste0(as.character(temp$ASD.code..prefix.),".",sprintf("%03d",seq(temp$From.ASD.spectra,temp$To.ASD.spectra)))
  if (any(columns %in% colnames(data.files[[as.character(temp[["Site.date"]])]]))){

    columns <- columns[which(columns %in% colnames(data.files[[as.character(temp[["Site.date"]])]]))]

    Reflectance <- data.files[[as.character(temp[["Site.date"]])]][,columns]
    metadata[["data.exist"]][i] <- TRUE

    temp.mat <- cbind(wv,Reflectance)
    colnames(temp.mat) <- c("wv",1,2,3)
    df.single <- temp.mat %>% pivot_longer(c("1","2","3")) %>% filter(wv >= wv.min,wv <= wv.max)

    data.all[[metadata[["GF"]][i]]] [[metadata[["Species"]][i]]] [[metadata[["Ind"]][i]]] <- list()
    data.all[[metadata[["GF"]][i]]] [[metadata[["Species"]][i]]] [[metadata[["Ind"]][i]]] [["spectrum"]] <- df.single
    data.all[[metadata[["GF"]][i]]] [[metadata[["Species"]][i]]] [[metadata[["Ind"]][i]]] [["site"]] <- metadata[["site"]][i]

    df.single_small <- df.single %>% filter(wv %in% seq(wv.min,wv.max,1))
    df.data.all <- rbind(df.data.all,
                         df.single_small %>% mutate(GF = metadata[["GF"]][i],
                                                    Species = metadata[["Species"]][i],
                                                    Ind = metadata[["Ind"]][i],
                                                    site = metadata[["site"]][i]))

  }
}

saveRDS(df.data.all,file = "./data/All.spectra.RDS")

summary <- metadata %>% filter(data.exist) %>% ungroup() %>% group_by(GF,site) %>% summarise(Nspectra = sum(N),
                                                                                             Nspecies = length(unique(Species)),
                                                                                             max_Nleaf = length(unique(Ind)))

sum(summary$Nspectra)

df.data.all_sum <- df.data.all %>% group_by(site,GF,wv) %>% summarise(R = mean(value),
                                                                      Rmin = quantile(value,0.025),
                                                                      Rmax = quantile(value,0.975))

ggplot(data = df.data.all_sum) +
  geom_ribbon(aes(x = wv,ymin=Rmin,ymax = Rmax,fill = GF),color = NA,alpha = 0.2) +
  geom_line(aes(x = wv,y=R,color = GF)) +
  facet_wrap(~ site) +
  theme_bw()

run_prospect <- function(params,waves){
  # PROSPECT parameters
  Nlayers <- params[1]
  Cab <- params[2]
  Car <- params[3]
  Cw <- params[4]
  Cm <- params[5]


  # Call RTM
  result <- tryCatch(
    Rprospect::prospect5(Nlayers,Cab,Car,Cw,Cm),
    error = function(e) NULL)
  if (is.null(result)) return(-1e20)
  reflectance <- result[["Reflectance"]]
  if (any(!is.finite(reflectance))) return(-1e20)
  if (any(reflectance < 0)) return(-1e20)
  if (any(reflectance > 1)) return(-1e20)

  reflectance_waves <- interp1(x = result[["Wavelength"]],
                               y = reflectance,
                               waves)
  return(reflectance_waves)
}

# test optimization
create_likelihood <- function(observed, waves) {
  function(params) {

    ssigma <- params[6]
    reflectance_waves <- run_prospect(params,waves)
    # Calculate likelihood
    ll <- sum(dnorm(reflectance_waves, observed, ssigma, log = TRUE))
    return(ll)
  }
}

# Uniform prior
Prospect_param_names <- c("Nlayers","Cab","Car","Cw","Cm","ssigma")
pft_lowers <- c(Nlayers = 1, Cab = 0 , Car = 0,Cw = 0, Cm = 0, ssigma = 0)
pft_uppers <-  c(Nlayers = 5, Cab = 100, Car = 50,Cw = 0.1, Cm = 0.1, ssigma = 1)

df.data.all_mod <- df.data.all
compt <- 1
for (pft in seq(names(data.all))){
  for (species in seq(names(data.all[[pft]]))){
    for (ind in seq(length(data.all[[pft]][[species]]))) {

      print(compt)
      data.temp <- data.all[[pft]][[species]][[ind]]

      if (!is.null(data.temp)){
        temp <- data.temp[["spectrum"]] %>% filter(!(wv > 680 & wv < 750))
        tp <- temp %>% group_by(wv) %>% summarise(max_diff = max(value)-min(value),
                                                  sd = sd(value),
                                                  m = mean(value))

        # If one spectrum is really too different, we remove it
        if (mean(tp$sd/tp$m)>0.10){

          # ggplot(data = data.temp$spectrum) +
          #   geom_line(aes(x = wv,y = value,group = name)) +
          #   theme_bw()

          df.data.all_mod <- df.data.all_mod %>% filter(!(GF == names(data.all)[pft] &
                                                            Species == species &
                                                            Ind == ind))

          print("Changing those data:")
          print(paste("pft = ",pft,", species = ",species,", ind = ",ind))
          tp2 <- temp %>% group_by(wv) %>% summarise(diff_1 = abs(mean(value[name==1]-mean(value[name!=1]))),
                                                     diff_2 = abs(mean(value[name==2]-mean(value[name!=2]))),
                                                     diff_3 = abs(mean(value[name==3]-mean(value[name!=3])))) %>% dplyr::select(!wv)
          temp <- temp %>% dplyr::filter(name != which.max(apply(tp2,2,sum)))
          data.all[[pft]][[species]][[ind]][["spectrum"]] <- data.temp[["spectrum"]] %>% dplyr::filter(name != which.max(apply(tp2,2,sum)))

          df.data.all_mod <- df.data.all_mod %>% filter(!(GF == names(data.all)[pft] &
                                                            Species == species &
                                                            Ind == ind &
                                                            name == which.max(apply(tp2,2,sum))))
        }

        observed <- as.vector(t(temp[,"value"]))
        waves <- as.vector(t(temp[,"wv"]))

        prior <- createUniformPrior(pft_lowers, pft_uppers)
        likelihood <- create_likelihood(observed=observed, waves)
        settings_MCMC <- list(iterations = 10000, nrChains = 1)

        # Run inversion
        setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
        samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
        samples <- BayesianTools::runMCMC(samples, settings = settings_MCMC)

        best_param <- BayesianTools::MAP(samples)$parametersMAP
        names(best_param) <- Prospect_param_names

        best_run <- run_prospect(best_param,waves)

        data.all[[pft]][[species]][[ind]] [["param"]] <- best_param
        data.all[[pft]][[species]][[ind]] [["best_run"]] <- best_run

        LM <- lm(formula = y ~ x, data = data.frame(x = observed, y = best_run))
        data.all[[pft]][[species]][[ind]] [["r2"]] <- summary(LM)$r.squared
        data.all[[pft]][[species]][[ind]] [["rmse"]] <- sqrt(mean(LM$residuals^2))
        data.all[[pft]][[species]][[ind]] [["samples"]] <- samples

        ggplot() +
          geom_line(data = temp, aes(x = wv, y = value, group = name,colour = name)) +
          geom_line(data = data.frame(wv = waves,
                                      name = temp$name,
                                      value = best_run), aes(x = wv, y = value, group = name),color = 'red') +
          theme_bw()

        print(c(summary(LM)$r.squared, 100*sqrt(mean(LM$residuals^2))/mean(observed)))

        compt <- compt +1
      }
    }
  }
}

df.data.all_mod_sum <- df.data.all_mod %>% group_by(site,GF,wv) %>% summarise(R = mean(value),
                                                                              Rmin = quantile(value,0.025),
                                                                              Rmax = quantile(value,0.975))

# saveRDS(file = "Sanchez_allfits.RDS",object = data.all)

wv_select <- seq(wv.min,wv.max,50)

# # Add posterior runs
# compt <- 1
# for (pft in seq(names(data.all))){
#   for (species in seq(names(data.all[[pft]]))){
#     for (ind in seq(length(data.all[[pft]][[species]]))) {
#       print(compt)
#
#       cdata <- data.all[[pft]][[species]][[ind]]
#
#       if (!is.null(cdata)){
#         cspectrum <- cdata[["spectrum"]]
#         observed <- as.vector(t(cspectrum[,"value"]))
#         waves <- as.vector(t(cspectrum[,"wv"]))
#         waves_select <- which(waves %in% wv_select)
#
#         posteriors <- getSample(cdata[["samples"]],numSamples = Nprospect)
#         colnames(posteriors) <- Prospect_param_names
#
#         df_runs <- data.frame()
#         for (irun in seq(1,min(nrow(posteriors),Nprospect))){
#           # print(irun/min(nrow(posteriors),Nprospect))
#           cparams <- c(posteriors[irun,"Nlayers"],posteriors[irun,"Cab"],posteriors[irun,"Car"],posteriors[irun,"Cw"],posteriors[irun,"Cm"])
#           crun <- run_prospect(cparams,waves[waves_select])
#           df_runs <- rbind(df_runs,
#                            data.frame(waves = waves[waves_select],obs = observed[waves_select],sim = crun,run = irun))
#         }
#
#         # leaves <- as.numeric(unique(as.vector(t(cspectrum[,"name"]))))
#         df_runs_sum <- df_runs %>% group_by(waves) %>% summarise(sim_m = quantile(sim,0.5),
#                                                                  sim_low = quantile(sim,0.025),
#                                                                  sim_high = quantile(sim,0.975),
#                                                                  obs_m = mean(obs[run==1]),
#                                                                  obs_low = mean(obs[run==1])-1.96*sd(obs[run==1]/length(obs[run==1])),
#                                                                  obs_high = mean(obs[run==1])+1.96*sd(obs[run==1]/length(obs[run==1])))
#
#         # ggplot(df_runs_sum,aes(x=obs_m,y = sim_m)) +
#         #   geom_point() +
#         #   geom_errorbarh(aes(xmin = obs_low,xmax=obs_high)) +
#         #   geom_errorbar(aes(ymin = sim_low,ymax=sim_high)) +
#         #   geom_abline(slope = 1,linetype = 2) +
#         #   theme_bw()
#
#         data.all[[pft]][[species]][[ind]][["posterior"]] <- df_runs_sum
#         compt <- compt + 1
#
#       }
#     }
#   }
# }

##############################################################
# saveRDS(file = "Sanchez_allfits2.RDS",object = data.all)

# data.all <- readRDS(file = "./Guzman_allfits_filter.RDS")
data.all <- readRDS(file = "./scripts/Sanchez_allfits2.RDS")
##############################################################
GFs <- names(data.all)
Prospect_param_names <- c("Nlayers","Cab","Car","Cw","Cm","ssigma")

# Add the full spectra
df.data.all_sim <- data.frame()
compt <- 1
Ntot = 0
for (pft in seq(names(data.all))){
  for (species in seq(names(data.all[[pft]]))){
    for (ind in seq(length(data.all[[pft]][[species]]))) {

      print(compt)
      data.temp <- data.all[[pft]][[species]][[ind]]

      if (!is.null(data.temp)){
        temp <- data.temp[["spectrum"]]
        waves <- unique(as.vector(t(temp[,"wv"])))

        best_param <- data.all[[pft]][[species]][[ind]][["param"]]
        names(best_param) <- Prospect_param_names

        best_run <- run_prospect(best_param,waves)

        data.all[[pft]][[species]][[ind]] [["best_run_full"]] <- best_run

        temp.mat <- cbind(waves,best_run)
        colnames(temp.mat) <- c("wv","value")

        df.single_small <- as.data.frame(temp.mat) %>% filter(wv >= wv.min,wv <= wv.max) %>% filter(wv %in% seq(wv.min,wv.max,5))
        df.data.all_sim <- rbind(df.data.all_sim,
                                 df.single_small %>% mutate(GF = GFs[pft],
                                                            Species = species,
                                                            Ind = ind,
                                                            site =data.all[[pft]][[species]][[ind]][["site"]]))

        Ntot = Ntot + length(unique(temp$name))

        compt <- compt +1
      }
    }
  }
}

print(Ntot)

df.data.all_sim_sum <- df.data.all_sim %>% group_by(site,GF,wv) %>% summarise(R = mean(value),
                                                                              Rmin = quantile(value,0.025),
                                                                              Rmax = quantile(value,0.975))

sims <- ggplot(data = df.data.all_sim_sum) +
  geom_ribbon(aes(x = wv,ymin=Rmin,ymax = Rmax,fill = GF),color = NA,alpha = 0.2) +
  geom_line(aes(x = wv,y=R,color = GF)) +
  facet_wrap(~ site) +
  theme_bw()

obs <- ggplot(data = df.data.all_mod_sum) +
  geom_ribbon(aes(x = wv,ymin=Rmin,ymax = Rmax,fill = GF),color = NA,alpha = 0.2) +
  geom_line(aes(x = wv,y=R,color = GF)) +
  facet_wrap(~ site) +
  theme_bw()

saveRDS(object = df.data.all_mod_sum %>% filter(site == "PNM") %>% rename(wavelength = wv,
                                                                          Reflectance = R) %>% mutate(pft = case_when(GF == "Liana" ~ "Liana_optical",
                                                                                                                      GF == "Tree" ~ "Tree_optical")) %>%
          ungroup() %>% dplyr::select(-c(Rmin,Rmax,site,GF)) %>% filter(!(wavelength > 680 & wavelength < 750)),
        file = "~/data/RTM/Figure6_sanchez2009_PNM_all2.rds")

saveRDS(object = df.data.all_mod_sum %>% filter(site == "FTS") %>% rename(wavelength = wv,
                                                                          Reflectance = R) %>% mutate(pft = case_when(GF == "Liana" ~ "Liana_optical",
                                                                                                                      GF == "Tree" ~ "Tree_optical")) %>%
          ungroup() %>% dplyr::select(-c(Rmin,Rmax,site,GF)) %>% filter(!(wavelength > 680 & wavelength < 750)),
        file = "~/data/RTM/Figure6_sanchez2009_FTS_all2.rds")

ggplot(data = df.data.all) +
  geom_line(aes(x = wv,y = value,color = GF,group = interaction(GF,Species,Ind,name))) +
  # geom_line(aes(x = wv,y=R,color = GF)) +
  scale_x_continuous(limits = c(400,700)) +
  facet_wrap(~ site) +
  theme_bw()

plot_grid(sims, obs,
          align = "h", axis = "b", nrow = 2, rel_widths = c(1, 1))


temp <- df.data.all_sim_sum %>% mutate(sim = R) %>% dplyr::select(site,GF,wv,sim) %>% left_join(df.data.all_sum%>% mutate(obs = R) %>% dplyr::select(site,GF,wv,obs))

ggplot(data = temp) +
  geom_line(aes(x = wv,y = sim,color = GF)) +
  geom_line(aes(x = wv,y = obs,color = GF),linetype = 2) +
  facet_wrap(~site) +
  theme_bw()
saveRDS(object = temp,file = "~/data/RTM/fit_Sanchez.RDS")

##############################################################

stats.all <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  r2 <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'r2'))
  rmse <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'rmse'))
  pos <- which(!sapply(data.all[[GFs[iGF]]][[i]],is.null))
  site <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'site'))
  if (is.null(r2)){
    return(NULL)
  } else{
    data.frame(r2 = r2, rmse = rmse, site = site,species = i, ind = pos, GF = GFs[iGF])
  }
}))
}))

ggplot(stats.all %>% pivot_longer(c(r2,rmse))) +
  geom_boxplot(aes(x = GF, y = value,fill = GF),alpha = 0.5) +
  facet_grid(name ~ site,scales = "free_y") +
  scale_color_manual(values = c("darkblue","darkgreen")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  theme_bw()
ggsave(plot = last_plot(),dpi = 300, width = 20,height = 20, filename = file.path("./Figures","Quality.of.fit.png"),units = "cm")


params.all <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  params <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'param'))
  site <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'site'))
  if (is.null(params)){
    return(NULL)
  } else{
    data.frame(value = params, params = sub(".*\\.", "", names(params)), site = site, species = i, GF = GFs[iGF])
  }
}))
}))

params.all2 <- params.all %>% filter(!(params == "ssigma"))
params.all2 %>% group_by(GF,params,site) %>% summarise(m = mean(value)) %>% arrange(params,site)
saveRDS(file = "~/data/RTM/params_sanchez.RDS",object = params.all2)

ggplot(data = params.all2,aes(x = value, y = GF, fill = GF)) +
  geom_density_ridges(alpha= 0.5) +
  facet_grid(site ~ params,scales = "free") +
  scale_color_manual(values = c("darkblue","darkgreen")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  theme_bw()

ggplot(data = params.all2,aes(x = value, y = GF, fill = GF)) +
  geom_density_ridges(alpha= 0.5) +
  facet_grid(~ params,scales = "free") +
  scale_color_manual(values = c("darkblue","darkgreen")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  theme_bw()

data.all_formated <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))

  clist <- data.all[[GFs[iGF]]][[i]]
  clist[sapply(clist, is.null)] <- NULL

  waves <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,1]}))
  pos <- which(!(waves > 680 & waves < 750))
  waves = waves[pos]
  num <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,2]}))[pos]
  Reflectance <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,3]}))[pos]
  best_run <- unlist(sapply(clist, '[', 'best_run'))
  site <- unlist(sapply(clist, '[', 'site'))
  site <- site[1]
  names(site) <- NULL

  select <- which(waves %in% wv_select)

  if (is.null(best_run)){
    return(NULL)
  } else{
    data.frame(wv = waves[select], num = num[select], obs = Reflectance[select], site = site,sim = best_run[select], species = i ,GF = GFs[iGF])
  }
}))
}))

ggplot(data.all_formated) +
  geom_abline(slope = 1, color = "black",size = 1.3) +
  geom_point(aes(x=obs,y=sim,color = GF, shape = site),alpha = 0.2,size = 0.5) +
  theme_bw()
#
ggplot(data.all_formated, aes(x=obs, y=sim)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density") +
  facet_wrap(site~GF,scales = "free") +
  geom_abline(slope = 1,color = "black") +
  theme_bw()
#   #+ geom_point(shape = '.')

ggplot(data.all_formated, aes(x=obs, y=sim)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon',alpha = 0.4) +
  scale_fill_viridis_c(name = "density",alpha = 0.4) +
  facet_wrap(~GF,scales = "free") +
  geom_abline(slope = 1,color = "black") +
  theme_bw()

A <- ggplot(data.all_formated %>% filter(GF == "Liana"), aes(x=obs, y=sim)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon',alpha = 0.4) +
  scale_fill_gradient2(low = ("white"),
                       high = ("#1E64C8")) +
  labs(x = "Observed reflectance [-]",
       y = "Simulated reflectance [-]") +
  scale_x_continuous(expand = c(0,0),limits = c(0,0.6)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.6)) +
  coord_fixed(ratio = 1) +
  geom_abline(slope = 1,color = "black") +
  theme_bw() + guides(fill = FALSE)

B <- ggplot(data.all_formated %>% filter(GF == "Tree"), aes(x=obs, y=sim)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon',alpha = 0.4) +
  scale_fill_gradient2(low = ("white"),
                       high = ("#137300")) +
  scale_x_continuous(expand = c(0,0),limits = c(0,0.6)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.6)) +
  labs(x = "Observed reflectance [-]",
       y = "Simulated reflectance [-]") +
  coord_fixed(ratio = 1) +
  geom_abline(slope = 1,color = "black") +
  theme_bw() + guides(fill = FALSE)

C <- ggplot(data.all_formated %>% filter(!(wv %in% seq(680,720))), aes(x=sim-obs,fill = GF)) +
  geom_histogram(alpha = 0.6,color = NA)+
  geom_hline(yintercept = 0) +
  labs(x = "Simulated - Observed",y = "") +
  scale_x_continuous(limits = c(-0.05,0.05)) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  theme_bw() + guides(fill = FALSE)
#   #+ geom_point(shape = '.')

plot_grid(plot_grid(A,B,
                    align = c("hv"),
                    nrow = 1,rel_widths = c(1,1,1)),C,nrow = 2,rel_widths = c(1,1))

obsvssim <- readRDS("~/data/RTM/obsvssim.RDS")


data.all_formated_merged <- bind_rows(list(data.all_formated,
                                           obsvssim %>% filter(ref %in% c("Castro_FTS","Castro_PNM","Guzman","Kalacska")) %>% ungroup() %>% dplyr::select(wavelength,Reflectance_median,posterior_median,pft,ref) %>%
                                             mutate(num = as.factor(1)) %>% rename(wv = wavelength,
                                                                                   obs = Reflectance_median,
                                                                                   sim = posterior_median,
                                                                                   site = ref,
                                                                                   GF = pft) %>% mutate(GF = case_when(GF == "Liana_optical" ~ "Liana",
                                                                                                                       GF == "Tree_optical" ~ "Tree"))))

ggplot(data = data.all_formated_merged) +
  geom_point(aes(x = obs,y = sim,color = GF)) +
  facet_wrap(~ site) +
  theme_bw()

data.all_formated_merged %>% ungroup() %>% summarise(summary(lm(formula = sim ~ obs))$r.squared)
data.all_formated_merged %>% group_by(site,GF) %>% summarise(summary(lm(formula = sim ~ obs))$r.squared)
data.all_formated_merged %>% group_by(GF) %>% summarise(summary(lm(formula = sim ~ obs))$r.squared)
data.all_formated_merged %>% group_by(GF) %>% summarise(RMSE = sqrt(sum((obs-sim)**2,na.rm = TRUE)/length(obs)))

data.all_formated_merged %>% ungroup() %>% summarise(coef(lm(formula = sim ~ obs))[2])
data.all_formated_merged %>% group_by(site,GF) %>% summarise(coef(lm(formula = sim ~ obs))[2])

data.all_formated %>% group_by(GF) %>% summarise(RMSE = sqrt(sum((obs-sim)**2)/length(obs)),
                                                 bias = sum(sim - obs)/length(obs),
                                                 m = mean(obs),
                                                 RMSE_r = RMSE/m,
                                                 N = length(obs))

data.all_formated %>% mutate(Band = case_when(wv < 700 ~ "VIS",
                                              wv < 1400 ~ "NIR",
                                              wv > 1500 & wv < 2500 ~ "SWIR")) %>% filter(Band %in% c("VIS","NIR","SWIR")) %>%
  group_by(GF,Band) %>% summarise(RMSE = sqrt(sum((obs-sim)**2)/length(obs)),
                                  bias = sum(sim - obs)/length(obs),
                                  m = mean(obs),
                                  RMSE_r = RMSE/m,
                                  N = length(obs))

ggsave(plot = last_plot(),dpi = 300, width = 22,height = 15, filename = file.path("./Figures","Quality.of.fit_final.png"),units = "cm")

ggplot(data.all_formated, aes(x=obs, y=sim)) +
  geom_pointdensity() + scale_color_viridis_c() +
  geom_abline(slope = 1, color = "black",size = 1.3) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(site~GF,scales = "free") +
  theme_bw()

ggplot(data.all_formated, aes(x=obs, y=sim)) +
  geom_pointdensity() + scale_color_viridis_c() +
  geom_abline(slope = 1, color = "black",size = 1.3) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw()

data.all_formated %>% group_by(GF,site) %>% summarise(r2=summary(lm(formula = sim ~ obs))$r.squared)
