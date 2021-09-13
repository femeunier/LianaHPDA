rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)

data.spectra <- readRDS("./data/All.spectra.RDS") %>%
  group_by(GF,Species,site) %>% mutate(species.id = cur_group_id()) %>%
  group_by(GF,Species,Ind,site) %>% mutate(ind.id = cur_group_id()) %>%
  group_by(GF,Species,Ind,site,name) %>% mutate(leaf.id = cur_group_id())

data.spectra %>% filter(wv == 400) %>% nrow()

data.spectra %>% filter(wv == 400) %>% group_by(GF,species.id,site,ind.id,wv) %>% summarise(N = length(value))
data.spectra %>% filter(wv == 400) %>% group_by(GF,species.id,site,wv) %>% summarise(N = length(value))
data.spectra %>% filter(wv == 400) %>% group_by(GF,site,wv) %>% summarise(N = length(value),
                                                                          Nspecies = length(unique(Species)))

data.spectra.L.FTS <- data.spectra %>% filter(GF == "Liana",site == "FTS") %>% group_by(Species) %>%
  mutate(species.id = cur_group_id())
data.spectra.T.FTS <- data.spectra %>% filter(GF == "Tree",site == "FTS") %>% group_by(Species) %>%
  mutate(species.id = cur_group_id())
data.spectra.L.PNM <- data.spectra %>% filter(GF == "Liana",site == "PNM") %>% group_by(Species) %>%
  mutate(species.id = cur_group_id())
data.spectra.T.PNM <- data.spectra %>% filter(GF == "Tree",site == "PNM") %>% group_by(Species) %>%
  mutate(species.id = cur_group_id())

data.spectra.all <- bind_rows(list(data.spectra.L.FTS,data.spectra.T.FTS,data.spectra.L.PNM,data.spectra.T.PNM))
saveRDS(object = data.spectra.all,
        file = "./data/All.spectra.ID.RDS")

Species.sum <- data.spectra.all %>% filter(wv == 400) %>% group_by(GF,species.id,site,wv) %>% summarise(N = length(value),
                                                                                                        .groups = "keep")
GF.sum <- data.spectra.all %>% filter(wv == 400) %>% group_by(GF,site,wv) %>% summarise(N = length(value),
                                                                                        Nspecies = length(unique(Species)),
                                                                                        .groups = "keep")


data.spectra.sum <- data.spectra.all %>% group_by(GF,species.id,site,wv) %>% summarise(value.m = mean(value),
                                                                                value.min = min(value),
                                                                                value.max = max(value),
                                                                                value.sd = sd(value),
                                                                                .groups = "keep")

data.spectra.all <- data.spectra.all %>% mutate(site = as.factor(site))
data.spectra.all$site <- relevel(data.spectra.all$site,"PNM")

data2plot.tree <- data.spectra.sum %>% filter(GF == "Tree") %>% mutate(site = as.factor(site))
data2plot.tree$site <- relevel(data2plot.tree$site,"PNM")

ggplot(data = data2plot.tree,
       aes(x = wv, y = value.m,
           group = (species.id))) +
  geom_ribbon(aes(ymin = value.min, ymax = value.max),
              alpha = 0.4, color = NA, fill = 'darkgreen') +
  geom_line(data = data.spectra.all %>% filter(GF == "Tree"),
            aes(group = interaction(species.id,leaf.id),y = value),size = 0.2) +
  # scale_x_continuous(limits = c(690,750)) +
  # geom_line(color = 'darkgreen',size = 1) +
  facet_wrap(species.id ~ site,nrow = 9) +
  labs(x = 'Wavelength (nm)', y = "Reflectance") +
  theme_bw() +
  theme(strip.text.x = element_blank())

data2plot.liana <- data.spectra.sum %>% filter(GF == "Liana") %>% mutate(site = as.factor(site))
data2plot.liana$site <- relevel(data2plot.liana$site,"PNM")

# %>%
#   mutate(site.2 = factor(case_when(site == "FTS" ~ "FTS",
#                             site == "PNM" & species.id <= 13 ~ "PNM.1",
#                             TRUE ~ "PNM.2"),
#                          levels = c("PNM.1","PNM.2","FTS"))) %>%
#   mutate(species.id.2 = as.factor(case_when(site.2 == "PNM.2" ~ as.numeric(species.id - 13),
#                                             TRUE ~ as.numeric(species.id))))
#
# data2plotall.liana <- data.spectra.all %>% filter(GF == "Liana") %>%
#   mutate(site.2 = factor(case_when(site == "FTS" ~ "FTS",
#                                    site == "PNM" & species.id <= 13 ~ "PNM.1",
#                                    TRUE ~ "PNM.2"),
#                          levels = c("PNM.1","PNM.2","FTS"))) %>%
#   mutate(species.id.2 = as.factor(case_when(site.2 == "PNM.2" ~ as.numeric(species.id - 13),
#                                             TRUE ~ as.numeric(species.id))))

ggplot(data = data2plot.liana,
       aes(x = wv, y = value.m,
           group = (species.id))) +
  geom_ribbon(aes(ymin = value.min, ymax = value.max),
              alpha = 0.4, color = NA, fill = "darkblue") +
  geom_line(data = data.spectra.all %>% filter(GF == "Liana"),
            aes(group = interaction(species.id,leaf.id),y = value),size = 0.2) +
  # geom_line(color = 'red',size = 1) +
  # scale_x_continuous(limits = c(710,730)) +
  facet_wrap(site ~ species.id,nrow = 13,dir="v") +
  labs(x = 'Wavelength (nm)', y = "Reflectance") +
  theme_bw() +
  theme(strip.text.x = element_blank())


# ggplot(data = data2plot.liana %>% filter(site == "FTS"),
#        aes(x = wv, y = value.m,
#            group = (species.id))) +
#   geom_ribbon(aes(ymin = value.min, ymax = value.max),
#               alpha = 0.4, color = NA, fill = "darkblue") +
#   geom_line(data = data.spectra.all %>% filter(GF == "Liana",
#                                                site == "FTS"),
#             aes(group = interaction(species.id,leaf.id),y = value),size = 0.2) +
#   # geom_line(color = 'red',size = 1) +
#   # scale_x_continuous(limits = c(710,730)) +
#   facet_wrap(site ~ species.id,nrow = 13,dir="v") +
#   labs(x = 'Wavelength (nm)', y = "Reflectance") +
#   theme_bw() +
#   theme(strip.text.x = element_blank())
