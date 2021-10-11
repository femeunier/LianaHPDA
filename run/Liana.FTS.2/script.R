rm(list=ls())

library(nimble)
library(coda)
library(dplyr)
library(tidyr)
library(LianaHPDA)
library(reshape2)

param.simu <- list(WLa = 400,WLb = 2500, Delta_WL = 20, WL.min = 450, WL.trans.min = 680, WL.trans.max = 800)

id <- list(iGF = 1,isite = 2, ispecies = 2)

param.MCMC <- list(Nchains = 2,nburnin = 5e+05, niter = 1e+06, thin = 500)

invert.single.species(id,param.MCMC,param.simu,'/data/gent/vo/000/gvo00074/felicien/R/LianaHPDA/out/Liana.FTS.2')
