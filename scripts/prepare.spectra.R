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
library(tcltk2)

wv.min = 400
wv.max = 2500

directory <- "/home/femeunier/Documents/projects/LianaHPDA/data/Guzman_all"
metadata.file <- file.path(directory,"meta_data_all2.csv")
metadata <- read.csv(metadata.file) %>% mutate(GF = case_when(substr(Code,1,1) == "L" ~ "Liana",
                                                              substr(Code,1,1) == "T" ~ "Tree"),
                                               Species.name = Species,
                                               Species = substring(sub("\\_.*", "", Code),2),
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

df.data.all <- data.frame()

total <- nrow(metadata)
pb <- txtProgressBar(min = 0, max = total, style = 3)

for (i in seq(1,total)){

  setTxtProgressBar(pb, i)

  temp <- metadata[i,]
  wv <- data.files[[as.character(temp[["Site.date"]])]][["Wavelength"]]
  columns <- paste0(as.character(temp$ASD.code..prefix.),".",sprintf("%03d",seq(temp$From.ASD.spectra,temp$To.ASD.spectra)))
  if (any(columns %in% colnames(data.files[[as.character(temp[["Site.date"]])]]))){

    columns <- columns[which(columns %in% colnames(data.files[[as.character(temp[["Site.date"]])]]))]

    Reflectance <- data.files[[as.character(temp[["Site.date"]])]][,columns]
    metadata[["data.exist"]][i] <- TRUE

    temp.mat <- cbind(wv,Reflectance)
    colnames(temp.mat) <- c("wv",1,2,3)
    df.single <- temp.mat %>% pivot_longer(c("1","2","3"),
                                           names_to = "rep") %>% filter(wv >= wv.min,wv <= wv.max)

    df.single_small <- df.single %>% filter(wv %in% seq(wv.min,wv.max,5))
    df.data.all <- rbind(df.data.all,
                         df.single_small %>% mutate(GF = metadata[["GF"]][i],
                                                    Species = metadata[["Species.name"]][i],
                                                    Ind = metadata[["Ind"]][i],
                                                    site = metadata[["site"]][i]))

  }
}

df.data.all_sum <- df.data.all %>% group_by(site,GF,wv) %>% summarise(R = mean(value),
                                                                      Rmin = quantile(value,0.025),
                                                                      Rmax = quantile(value,0.975))

df.summary <- df.data.all %>% filter(wv == 400) %>% group_by(GF,site) %>% summarise(N = length(rep),
                                                                                    Nspecies = length(unique(Species)),
                                                                                    NInd = length(unique(Ind)),
                                                                                    Nrep = length(rep)/length(unique(Ind))/length(unique(Species)),
                                                                                    .groups = "keep")

df.data.all %>% filter(site == "PNM",GF == "Tree") %>% pull(Species) %>% unique() %>% as.character() %>% sort()

metadata %>% filter(site == "PNM",GF == "Tree") %>% pull(Species.name) %>% unique() %>% as.character() %>% sort()

ggplot(data = df.data.all_sum) +
  geom_ribbon(aes(x = wv,ymin=Rmin,ymax = Rmax,fill = GF),color = NA,alpha = 0.2) +
  geom_line(aes(x = wv,y=R,color = GF)) +
  facet_wrap(~ site) +
  theme_bw()
