rm(list = ls())

library(ggplot2)
library(dplyr)
library(pracma)

data.Sanchez <- readRDS(file = "/home/femeunier/Documents/projects/LianaHPDA/data/All.spectra.ID.RDS") %>% group_by(wv,site,GF) %>%
  summarise(R = mean(value),
            .groups = "keep") %>% mutate(ref = "Sanchez")

data.rem <- readRDS("/home/femeunier/Documents/projects/LianaHPDA/data/All_leaf_spectra.rds") %>% filter(!(ref %in% c("Sanchez_FTS","Sanchez_PNM"))) %>%
  dplyr::select(ref,pft,wavelength,Reflectance_median) %>%
  rename(wv = wavelength,
         R = Reflectance_median) %>%
  mutate(site = gsub(".*_","",ref),
         ref =  gsub("\\_.*","",ref),
         GF = case_when(pft == "Tree_optical" ~ "Tree",
                        pft == "Liana_optical" ~ "Liana")) %>% ungroup() %>% dplyr::select(-pft) %>%
  mutate(site = case_when(site == "Guzman" ~ "SRNP",
                          site == "Kalacska" ~ "PNM",
                          TRUE ~ site))

data.all <- bind_rows(list(data.rem,
                           data.Sanchez))


all.wls <- 400:2500

studies <- unique(data.all$ref) ; sites <- unique(data.all$site) ; GFs <- unique(data.all$GF)

data.all.interp <- data.frame()
for (i in seq(1:length(studies))){
  for (j in seq(1:length(GFs))){
    for (k in seq(1:length(sites))){

      data.select <- data.all %>% filter(ref == studies[i],
                                         GF == GFs[j],
                                         site == sites[k])

      if (nrow(data.select)>0){

        c.wls <- data.select$wv
        c.R <- data.select$R

        n.wls <- all.wls[all.wls >= min(c.wls) & all.wls <= max(c.wls)]
        interp.R <- interp1(c.wls,c.R,n.wls)
        all.R <- rep(NA,length(all.wls))
        all.R[all.wls %in% n.wls] <- interp.R

        data.all.interp <- bind_rows(list(data.all.interp,
                                          data.frame(ref = studies[i],
                                                     GF = GFs[j],
                                                     site = sites[k],
                                                     wv = all.wls,
                                                     R = all.R)))
      }

    }
  }
}


ggplot(data = data.all) +
  geom_line(aes(x = wv, y = R, color = GF, group = interaction(GF,ref,site))) +
  theme_bw()

saveRDS(object = data.all.interp,
        file = "/home/femeunier/Documents/projects/LianaHPDA/data/All_leaf_spectra_allstudies.rds")
