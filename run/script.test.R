rm(list=ls())

library(nimble)
library(coda)
library(dplyr)
library(tidyr)
library(LianaHPDA)
library(reshape2)

param.simu <- list(WLa = 400,WLb = 2500, Delta_WL = 20, WL.min = 450, WL.trans.min = 680, WL.trans.max = 800)

id <- list(iGF = 2,isite = 2, ispecies = 3)

param.MCMC <- list(Nchains = 2,nburnin = 100, niter = 200, thin = 1)

invert.single.species(id,param.MCMC,param.simu,'./outputs/')
