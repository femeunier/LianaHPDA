rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(LianaHPDA)
library(reshape2)

param.simu <- list(WLa = 400,
                   WLb = 2500,
                   Delta_WL = 20,
                   WL.min = 450,
                   WL.trans.min = 680,
                   WL.trans.max = 800)

param.MCMC <- list(Nchains = 2,
                   nburnin = 500000,
                   niter = 1000000,
                   thin = 500)

GFs <- c('Tree','Liana')
sites <- c('PNM','FTS')

main.dir <- "/data/gent/vo/000/gvo00074/felicien/R/LianaHPDA"
main.dir.run <- file.path(main.dir,"run")
main.dir.out <- file.path(main.dir,"out")

if (!dir.exists(main.dir.run)) dir.create(main.dir.run)
if (!dir.exists(main.dir.out)) dir.create(main.dir.out)

data.raw <- LianaHPDA::array_obs_reflectance

# Submission parameters
args <- c('-l walltime=72:00:00',
          '-l nodes=1:ppn=16',
          '-o log/BCI.o$PBS_JOBID',
          '-e log/BCI.e$PBS_JOBID')

run_per_nodes <- 16
folders.all <- c()

for (iGF in seq(1,2)){
  for(isite in seq(1,2)){
    for (ispecies in seq(1,30)){

      if (!all(is.na(data.raw[,iGF,isite,ispecies,,]))){

        id <- list(iGF = iGF,
                   isite = isite,
                   ispecies = ispecies)

        dir.run <- file.path(main.dir.run,paste(GFs[iGF],sites[isite],ispecies,sep = "."))
        dir.out <- file.path(main.dir.out,paste(GFs[iGF],sites[isite],ispecies,sep = "."))

        if (!dir.exists(dir.run)) dir.create(dir.run)
        if (!dir.exists(dir.out)) dir.create(dir.out)

        create_Rscript_PROSPECT5(file.path(dir.run,"script.R"),id,param.MCMC,param.simu,dir.out, basename = "Rcorrected_noleaf",fun = "invert.single.species_noleafeffect")
        folders.all <- c(folders.all,dir.run)

      }
    }
  }
}

Njobs <- ceiling(length(folders.all)/run_per_nodes)
folder.runs <- folders.all[seq(1,length(folders.all),run_per_nodes)]

for (ijob in seq(1,Njobs)){
  cfolder <- folder.runs[ijob]
  simus <- folders.all[((ijob - 1)*run_per_nodes + 1):(ijob*run_per_nodes)]

  job_list_file <- file.path(cfolder,'joblist.txt')
  writeLines(paste0("echo ",'"',"source('script.R')",'"',"| R --save"),con = job_list_file)

  cmd_l <- lapply(simus,function(sim){
    write(sim,file=job_list_file,append=TRUE)})

  launcher_file <- file.path(cfolder,'launcher.sh')
  writeLines("#!/bin/bash",con = launcher_file)
  write("ml purge",file=launcher_file,append=TRUE)
  write("ml R/3.5.1-intel-2018b",file=launcher_file,append=TRUE)
  write(paste("mpirun",file.path(main.dir,"modellauncher","modellauncher"),job_list_file),file=launcher_file,append=TRUE)
}

# Submitting files
cmd <- "qsub"
for (ijob in seq(1,Njobs)){
  cfolder <- folder.runs[ijob]
  launcher_file <- file.path(cfolder,'launcher.sh')
  out <- system2(cmd, c(args,launcher_file), stdout = TRUE, stderr = TRUE)
}
