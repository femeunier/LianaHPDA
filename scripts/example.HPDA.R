rm(list = ls())

if (!require("PEcAn.photosynthesis", character.only = TRUE)) {
  devtools::install_github("PecanProject/pecan/modules/photosynthesis")
}
knitr::opts_chunk$set(cache = TRUE)

library(PEcAn.photosynthesis)

## Get list of LI-COR 6400 file names (ASCII not xls)
filenames <- system.file("extdata", paste0("flux-course-",rep(1:6,each=2),c("aci","aq")), package = "PEcAn.photosynthesis")

## Load files to a list
master = lapply(filenames, read_Licor)

dat <- do.call("rbind", master)

## if QC was done, remove both unchecked points and those that fail QC
if("QC" %in% colnames(dat)){
  dat <- dat[-which(dat$QC < 1),]
} else {
  QC <- rep(1, nrow(dat))
  dat <- cbind(dat, QC)
}

if(file.exists("fit.RData")){
  load("fit.RData")
} else{
  fit <- fitA(dat)
  save(fit, file = "fit.RData")
}

writeLines(fit$model)

par(mfrow = c(2, 1))
plot(fit$params, auto.layout = FALSE)    ## MCMC diagnostic plots
summary(fit$params) ## parameter estimates
gelman.plot(fit$params, auto.layout = FALSE)
gelman.diag(fit$params)

## predicted vs observed plot
par(mfrow = c(1, 1))
mstats <- summary(fit$predict)
pmean <- mstats$statistics[grep("pmean", rownames(mstats$statistics)), 1]
plot(pmean, dat$Photo, pch = "+", xlab = "Predicted A", ylab = "Observed A")
abline(0, 1, col = 2, lwd = 2)

##
plot_photo(dat, fit)

A.model <- list(a.fixed = NULL, a.random = "leaf",
                V.fixed = NULL, V.random = "leaf",
                n.iter = 5000, match = "fname")

if(file.exists("fitI.RData")){
  load("fitI.RData")
} else{
  fitI <- fitA(dat,model = A.model)
  save(fitI, file="fitI.RData")
}

par(mfrow = c(2, 1))
plot(fitI$params, auto.layout = FALSE)    ## MCMC diagnostic plotssummary
gelman.plot(fitI$params)
gelman.diag(fitI$params)

## predicted vs observed plot
par(mfrow = c(1, 1))
mstats <- summary(fitI$predict)
pmean <- mstats$statistics[grep("pmean",rownames(mstats$statistics)),1]
plot(pmean, dat$Photo, pch = "+", xlab = "Predicted A", ylab = "Observed A")
abline(0,1,col=2,lwd=2)

##
plot_photo(dat,fitI)
