## code to prepare the spectra dataset goes here
rm(list = ls())

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

data.spectra <- readRDS("./data/All.spectra.RDS") %>%
  group_by(GF,site,Species) %>% mutate(species.id = cur_group_id()) %>% group_by(GF,site) %>% mutate(species.id = species.id - min(species.id) + 1) %>%
  group_by(GF,site,species.id,Ind) %>% mutate(ind.id = cur_group_id()) %>% group_by(GF,site,species.id) %>% mutate(ind.id = ind.id - min(ind.id) + 1) %>%
  group_by(GF,site,species.id,ind.id,name) %>% mutate(leaf.id = cur_group_id()) %>% group_by(GF,site,species.id,ind.id) %>% mutate(leaf.id = leaf.id - min(leaf.id) + 1) %>%
  mutate(Identification = paste(GF,site,species.id,ind.id,leaf.id,sep = "_"))

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

# iGF = 1
# isite = 2
# ispecies = 2
# iind = 10

wls <- c(550,750,1450,2000)
Nwl <- length(wls)

for (iGF in seq(1,length(GFs))){
  for (isite in seq(1,length(sites))){
    cdata <- data.spectra %>% filter(GF == GFs[iGF],
                                     site == sites[isite])
    Nspecies[iGF,isite] <- min(maxNspecies,length(unique(cdata$species.id)))

    for (ispecies in seq(1,Nspecies[iGF,isite])){
      ccdata <- cdata %>% filter(species.id == ispecies)

      outliers <- data.frame()
      for (iwl in seq(1,Nwl)){
        values <- ccdata %>% filter(wv == wls[iwl]) %>% pull(value)
        Q <- quantile(values,c(0.25,0.5,0.75))
        IQR <- Q[3] - Q[1]

        outliers <- bind_rows(list(
          outliers,
          ccdata %>% filter(wv == wls[iwl]) %>% ungroup() %>% filter(value < (Q[1] - 1.5*IQR) | value > (Q[3] + 1.5*IQR))))
      }

      test.all <- names(table(outliers$Identification))[table(outliers$Identification) >= 2]

      if (length(test.all)>0){


        print(test.all[["Identification"]][1])

        print(ggplot(data = ccdata) +
          geom_line(aes(x = wv, y = value,
                        group = interaction(as.factor(ind.id),as.factor(leaf.id))),
                    color = "black") +
          geom_line(data = ccdata %>% filter(Identification %in% test.all),
                    aes(x = wv, y = value,
                        group = interaction(as.factor(ind.id),as.factor(leaf.id))),
                    color = "red") +
          theme_bw())

      }


    }
  }
}


#
# ccdata <- data.spectra %>% filter(GF == "Tree",
#                                   site == "FTS",
#                                   species.id == 1)
#
# outliers <- data.frame()
# for (iwl in seq(1,Nwl)){
#   values <- ccdata %>% filter(wv == wls[iwl]) %>% pull(value)
#   Q <- quantile(values,c(0.25,0.5,0.75))
#   IQR <- Q[3] - Q[1]
#
#   outliers <- bind_rows(list(
#     outliers,
#     ccdata %>% filter(wv == wls[iwl]) %>% ungroup() %>% filter(value < (Q[1] - 1.5*IQR) | value > (Q[3] + 1.5*IQR))))
# }





