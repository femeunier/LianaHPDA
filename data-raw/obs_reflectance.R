## code to prepare the spectra dataset goes here
rm(list = ls())

library(dplyr)
library(purrr)
library(tidyr)

data.spectra <- readRDS("./data/All.spectra_corrected.RDS") %>%
  group_by(GF,site,Species) %>% mutate(species.id = cur_group_id()) %>% group_by(GF,site) %>% mutate(species.id = species.id - min(species.id) + 1) %>%
  group_by(GF,site,species.id,Ind) %>% mutate(ind.id = cur_group_id()) %>% group_by(GF,site,species.id) %>% mutate(ind.id = ind.id - min(ind.id) + 1) %>%
  group_by(GF,site,species.id,ind.id,name) %>% mutate(leaf.id = cur_group_id()) %>% group_by(GF,site,species.id,ind.id) %>% mutate(leaf.id = leaf.id - min(leaf.id) + 1)

GFs <- c('Tree','Liana')
sites <- c('PNM','FTS')

maxNspecies <- Inf
maxNind <- Inf

array_obs_reflectance <- array(data = NA, dim = c(length(400:2500),2,2,30,15,3),dimnames = list(seq(400,2500),
                                                                                                GFs,
                                                                                                sites,
                                                                                                seq(1,30),
                                                                                                seq(1:15),
                                                                                                seq(1:3)))
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

        array_obs_reflectance[,iGF,isite,ispecies,iind,1] <- cccdata %>% filter(name == 1) %>% arrange(wv) %>% pull(value.corrected)
        array_obs_reflectance[,iGF,isite,ispecies,iind,2] <- cccdata %>% filter(name == 2) %>% arrange(wv) %>% pull(value.corrected)
        array_obs_reflectance[,iGF,isite,ispecies,iind,3] <- cccdata %>% filter(name == 3) %>% arrange(wv) %>% pull(value.corrected)
      }
    }
  }
}

# Hard-coded: datasets with problems
# Clean dataset
array_obs_reflectance[,1,1,5,,] <- NA # Weird species
array_obs_reflectance[,2,1,1,4,] <- NA # Weird curve
array_obs_reflectance[,2,1,13,,] <- NA # Weird species
array_obs_reflectance[,2,1,17,,] <- NA # Weird species
array_obs_reflectance[,1,2,6,c(1,3,4,5,6,7,10),] <- NA # Weird curves

# Too much water (weird single cuves)
array_obs_reflectance[,1,2,5,2,1] <- NA
array_obs_reflectance[,1,2,5,5,1] <- NA
array_obs_reflectance[,1,2,5,10,1] <- NA

array_obs_reflectance[,1,1,2,,] <- NA # Species with very large reflectances ...
array_obs_reflectance[,2,2,3,5,3] <- NA # weird curve (zero reflectance)

maxNspecies <- max(Nspecies,na.rm = TRUE)
maxNind <- max(Nind,na.rm = TRUE)

usethis::use_data(data.spectra, overwrite = TRUE)
usethis::use_data(array_obs_reflectance, overwrite = TRUE)
usethis::use_data(Nspecies, overwrite = TRUE)
usethis::use_data(Nind, overwrite = TRUE)

