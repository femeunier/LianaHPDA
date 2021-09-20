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
  mutate(Identification = paste(GF,site,species.id,ind.id,leaf.id,sep = "_")) %>%
  mutate(value.corrected = value)


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

plot(NA,NA,xlim = c(400,2500),ylim = c(0,0.7),xlab = "",ylab = "R")
compt <- 1

# wv.min = 600

# Correct the jumps
for (iGF in seq(1,length(GFs))){
  for (isite in seq(1,length(sites))){
    cdata <- data.spectra %>% filter(GF == GFs[iGF],
                                     site == sites[isite])
    Nspecies[iGF,isite] <- min(maxNspecies,length(unique(cdata$species.id)))

    for (ispecies in seq(1,Nspecies[iGF,isite])){

      print(c(iGF,isite,ispecies))

      ccdata <- cdata %>% filter(species.id == ispecies)
      Nind[iGF,isite,ispecies] <- min(maxNspecies,length(unique(ccdata$ind.id)))

      for (iind in seq(1,Nind[iGF,isite,ispecies])){
        for (ileaf in seq(1,3)){

          cccdata <- ccdata %>% filter(ind.id == iind,
                                       leaf.id == ileaf)

          Delta_R <- mean(cccdata$value[cccdata$wv>1000 & cccdata$wv<=1001]) -
            mean(cccdata$value[cccdata$wv>998 & cccdata$wv<1000])
          rel_R <- mean(cccdata$value[cccdata$wv>1000 & cccdata$wv<=1001])/
            mean(cccdata$value[cccdata$wv>998 & cccdata$wv<1000])

          Delta_R2 <- mean(cccdata$value[cccdata$wv>1830 & cccdata$wv<1840]) -
            mean(cccdata$value[cccdata$wv>1820 & cccdata$wv<1830])
          rel_R2 <- mean(cccdata$value[cccdata$wv>1831 & cccdata$wv<1833])/
            mean(cccdata$value[cccdata$wv>1828 & cccdata$wv<=1830])


          # if (abs(Delta_R) > 0.01){
          cccdata <- cccdata %>% ungroup() %>% mutate(value.corrected = case_when(wv > 1000 ~ value.corrected,
                                                                                  TRUE ~ value.corrected *rel_R))
          # }

          # if (abs(Delta_R2) > 0.01){
          cccdata <- cccdata %>% ungroup() %>% mutate(value.corrected = case_when(wv <= 1830 ~ value.corrected,
                                                                                  TRUE ~ value.corrected /rel_R2))
          # }


          # lines(cccdata$wv,cccdata$value,type = 'l')
          # lines(cccdata$wv,cccdata$value.corrected,type = 'l',col = "red")
          # stop()

          # First
          if (abs(Delta_R2) > 0.02 & abs(Delta_R) > 0.02){

            # print(paste0("here, ", c(iGF,isite,ispecies)))


            lines(cccdata$wv,cccdata$value,type = 'l')
            lines(cccdata$wv,cccdata$value.corrected,type = 'l',col = "red")

            compt <- compt + 1
            if (compt > 5){
              plot(NA,NA,xlim = c(400,2500),ylim = c(0.,0.7),xlab = "",ylab = "R")
              compt <- 1
            }
          }

          # Correction
          data.spectra[data.spectra$GF == GFs[iGF] &
                         data.spectra$site == sites[isite] &
                         data.spectra$species.id == ispecies &
                         data.spectra$ind.id == iind &
                         data.spectra$leaf.id == ileaf, "value.corrected"] <- cccdata[["value.corrected"]]


        }
      }
    }
  }
}

saveRDS(data.spectra,file = "./data/All.spectra_corrected.RDS")

ggplot(data = data.spectra %>% filter(GF == "Liana",
                                      site == "PNM",
                                      species.id %in% seq(6,10)),
aes(x = wv, y = value, group = interaction(site,GF,as.factor(species.id),as.factor(ind.id),as.factor(leaf.id)))) +
  # geom_line() +
  geom_line(aes(y = value.corrected),color = "red",linetype = 1) +
  # scale_x_continuous(limits = c(950,1050)) +
  # scale_y_continuous(limits = c(0.45,0.55)) +
  theme_bw()


"# Correct the shift for the red-edge/green peak

for (iGF in seq(1,length(GFs))){
  for (isite in seq(1,length(sites))){
    cdata <- data.spectra %>% filter(GF == GFs[iGF],
                                     site == sites[isite])
    Nspecies[iGF,isite] <- min(maxNspecies,length(unique(cdata$species.id)))

    for (ispecies in seq(1,Nspecies[iGF,isite])){

      print(c(iGF,isite,ispecies))

      ccdata <- cdata %>% filter(species.id == ispecies)
      Nind[iGF,isite,ispecies] <- min(maxNspecies,length(unique(ccdata$ind.id)))

      for (iind in seq(1,Nind[iGF,isite,ispecies])){
        for (ileaf in seq(1,3)){

          cccdata <- ccdata %>% filter(ind.id == iind,
                                       leaf.id == ileaf)

          WLs <- cccdata$wv
          R <- cccdata$value.corrected
          gradient.data <-pracma::gradient(R)


          if (!(mean(R[WLs > 525 & WLs < 575]) > mean(R[WLs < 525 & WLs > 500]) &
              mean(R[WLs > 525 & WLs < 575]) > mean(R[WLs < 600 & WLs > 575]))){
            plot(cccdata$wv,cccdata$value.corrected,xlim = c(400,700),ylim = c(0,0.12),type = "l")
          }

        }
      }
    }
  }
}


# Flag the bad-quality (soil-like)

# Flag the very different spectra

# Flag the very low waterered spectra?

# Spectra with R = 0

