create_Rscript_PROSPECT5 <- function(file,id,param.MCMC,param.simu,OP.folder) {

  ########################################################################

  iGF <- id[["iGF"]]
  isite <- id[["isite"]]
  ispecies <- id[["ispecies"]]

  Nchains <- param.MCMC[["Nchains"]]
  nburnin <- param.MCMC[["nburnin"]]
  niter <- param.MCMC[["niter"]]
  thin <- param.MCMC[["thin"]]

  WLa <- param.simu[["WLa"]]
  WLb <- param.simu[["WLb"]]
  Delta_WL <- param.simu[["Delta_WL"]]
  WL.min <- param.simu[["WL.min"]]
  WL.trans.min <- param.simu[["WL.trans.min"]]
  WL.trans.max <- param.simu[["WL.trans.max"]]

  ########################################################################

  writeLines("rm(list=ls())",con = file)

  write("",file=file,append=TRUE)
  write("library(nimble)",file=file,append=TRUE)
  write("library(coda)",file=file,append=TRUE)
  write("library(dplyr)",file=file,append=TRUE)
  write("library(tidyr)",file=file,append=TRUE)
  write("library(LianaHPDA)",file=file,append=TRUE)
  write("library(reshape2)",file=file,append=TRUE)

  write("",file=file,append=TRUE)
  write(paste0("param.simu <- list(WLa = ",WLa,",WLb = ",WLb,", Delta_WL = ",Delta_WL,", WL.min = ",WL.min,", WL.trans.min = ",WL.trans.min,", WL.trans.max = ",WL.trans.max,")"),file=file,append=TRUE)

  write("",file=file,append=TRUE)
  write(paste0("id <- list(iGF = ",iGF,",isite = ",isite,", ispecies = ",ispecies,")"),file=file,append=TRUE)

  write("",file=file,append=TRUE)
  write(paste0("param.MCMC <- list(Nchains = ",Nchains,",nburnin = ",nburnin,", niter = ",niter,", thin = ",thin,")"), file=file,append=TRUE)

  write("",file=file,append=TRUE)
  write(paste0("invert.single.species(id,param.MCMC,param.simu,'",OP.folder,"')"),file=file,append=TRUE)

}
