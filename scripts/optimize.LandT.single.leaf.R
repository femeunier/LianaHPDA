rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)

e1_approx <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))

    A <- log((0.56146 / x + 0.65) * (1 + x))
    B <- x^4 * exp(7.7 * x) * (2 + x)^3.7

    return((A^-7.7 + B)^-0.13)
  })

NIMprospect5 <- nimbleFunction(
  run = function(N = double(0),Cab = double(0),Car = double(0), Cw = double(0), Cm = double(0),
                 dataspec_p5 = double(2),talf = double(1),t12 = double(1),t21 = double(1), Nwl = double(0)) {

    cc <- matrix(NA,nrow = 5,ncol = 1)
    k <- numeric(length = Nwl)

    Cbrown <- 0

    cc[1,1] <- Cab / N
    cc[2,1] <- Car / N
    cc[3,1] <- Cbrown / N
    cc[4,1] <- Cw / N
    cc[5,1] <- Cm / N

    k[] <- dataspec_p5[,] %*% cc[,]

    trans <- (1 - k)*exp(-k) + k^2 *e1_approx(k)
    trans[trans < 0] <- 0

    ralf <- 1 - talf
    r12 <- 1 - t12
    r21 <- 1 - t21

    denom <- 1 - (r21 ^ 2) * (trans ^ 2)
    Ta <- talf * trans * t21 / denom
    Ra <- ralf + r21 * trans * Ta

    tt <- t12 * trans * t21 / denom
    rr <- r12 + r21 * trans * tt

    gt1 <- rr + tt >= 1
    tgt1 <- tt[gt1]

    Tsub <- 0*tt
    Rsub <- 0*tt
    r <- 0*tt
    t <- 0*tt

    Tsub[gt1] <- tgt1 / (tgt1 + (1 - tgt1) * (N - 1))
    Rsub[gt1] <- 1 - Tsub[gt1]

    inf <- rr == 0 | tt == 0
    Tsub[inf] <- 0
    Rsub[inf] <- 0

    r <- rr[!gt1 & !inf]
    t <- tt[!gt1 & !inf]

    D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    r2 <- r ^ 2
    t2 <- t ^ 2
    va <- (1 + r2 - t2 + D) / (2 * r)
    vb <- (1 - r2 + t2 + D) / (2 * t)

    vbNN <- vb ^ (N - 1)
    vbNN2 <- vbNN ^ 2
    va2 <- va ^ 2
    denomx <- va2 * vbNN2 - 1
    Rsub[!gt1 & !inf] <- va * (vbNN2 - 1) / denomx
    Tsub[!gt1 & !inf] <- vbNN * (va2 - 1) / denomx

    denomy <- 1 - Rsub * rr

    reflectance <- Ra + Ta * Rsub * tt / denomy

    returnType(double(1))
    # print(c(N,Cab, Car, Cbrown, Cw, Cm,mean(reflectance)))
    return(reflectance)
  })

run_prospect5 <- nimbleCode({

  reflectance_T[1:Nwl] <- NIMprospect5(Nmean,Cabmean,Carmean,Cwmean,Cmmean,
                                       dataspec_p5[,], talf[],t12[],t21[], Nwl)

  reflectance_L[1:Nwl] <- NIMprospect5(max(1.1,Nmean + alpha_N),
                                       max(10,Cabmean + alpha_Cab),
                                       max(5,Carmean + alpha_Car),
                                       max(0.0001,Cwmean + alpha_Cw),
                                       max(0.0001,Cmmean + alpha_Cm),
                                       dataspec_p5[,], talf[],t12[],t21[], Nwl)

  for (i in 1:Nsamples){
    creflectance[,i] <- reflectance_T[1:Nwl]
    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(creflectance[j,i], sd = Standard.Dev)
    }
  }

  for (i in 1:Nsamples){
    creflectance[,Nsamples + i] <- reflectance_L[1:Nwl]
    for (j in 1:Nwl){
      obs_reflectance[j,Nsamples + i] ~ dnorm(creflectance[j,Nsamples + i], sd = max(0.00001,Standard.Dev + alpha_SD))
    }
  }

  Nmean ~ dunif(1.,5)

  Cabmean ~ dunif(0,100)
  Carmean ~ dunif(0,50)
  Cwmean ~ dunif(0.,0.1)
  Cmmean ~ dunif(0.,0.1)
  Standard.Dev ~ dunif(0,1)

  alpha_N ~ dunif(-1,1)
  alpha_Cab ~ dunif(-30,30)
  alpha_Car ~ dunif(-20,20)
  alpha_Cw ~ dunif(-0.02,0.02)
  alpha_Cm ~ dunif(-0.02,0.02)

  alpha_SD ~ dunif(-0.1,0.1)

})

Nleaves <- 10
WLa <- 400
WLb <- 2500
Delta_WL <- 20
WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)


Constants <- list(Nsamples = Nleaves,
                  Nwl = Nwl,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = rrtm:::dataspec_p5[pos,1:5])


Ntrue <- 2 ; Cabtrue <- 40 ; Cartrue <- 20 ; Cwtrue <- 0.01 ; Cmtrue <- 0.01
Ntrue_L <- 1.5 ; Cabtrue_L <- 30 ; Cartrue_L <- 15 ; Cwtrue_L <- 0.005 ; Cmtrue_L <- 0.005

fac <- 4

Nall <- pmax(1.1,rnorm(Nleaves,Ntrue,Ntrue/fac))
Caball <- pmax(10,rnorm(Nleaves,Cabtrue,Cabtrue/fac))
Carall <- pmax(5,rnorm(Nleaves,Cartrue,Cartrue/fac))
Cwall <- pmax(0,rnorm(Nleaves,Cwtrue,Cwtrue/fac))
Cmall <- pmax(0,rnorm(Nleaves,Cmtrue,Cmtrue/fac))

Nall_L <- pmax(1.1,rnorm(Nleaves,Ntrue_L,Ntrue/fac))
Caball_L <- pmax(10,rnorm(Nleaves,Cabtrue_L,Cabtrue/fac))
Carall_L <- pmax(5,rnorm(Nleaves,Cartrue_L,Cartrue/fac))
Cwall_L <- pmax(0,rnorm(Nleaves,Cwtrue_L,Cwtrue/fac))
Cmall_L <- pmax(0,rnorm(Nleaves,Cmtrue_L,Cmtrue/fac))

Data <- list(obs_reflectance = cbind(matrix(unlist(lapply(1:length(Nall),function(ileaf){
  rrtm::prospect5(N = Nall[ileaf],
                  Cab = Caball[ileaf],
                  Car = Carall[ileaf],
                  Cbrown = 0,
                  Cw = Cwall[ileaf],
                  Cm = Cmall[ileaf])[["reflectance"]][pos]})),ncol = Nleaves),
  matrix(unlist(lapply(1:length(Nall),function(ileaf){
    rrtm::prospect5(N = Nall_L[ileaf],
                    Cab = Caball_L[ileaf],
                    Car = Carall_L[ileaf],
                    Cbrown = 0,
                    Cw = Cwall_L[ileaf],
                    Cm = Cmall_L[ileaf])[["reflectance"]][pos]})),ncol = Nleaves)))

matplot(WLs,Data$obs_reflectance[,1:Nleaves],type = 'l',col = "darkgreen")
matlines(WLs,Data$obs_reflectance[,Nleaves + (1:Nleaves)],type = 'l',col = "darkblue")

Inits <- list(Nmean = Ntrue,
              Cabmean = Cabtrue,
              Carmean = Cartrue,
              Cwmean = Cwtrue,
              Cmmean = Cmtrue,
              alpha_N = 0,
              alpha_Cab = 0,
              alpha_Car = 0,
              alpha_Cw = 0,
              alpha_Cm = 0)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance_T = c(Nwl),
                                         reflectance_L = c(Nwl),
                                         creflectance = c(Nwl,Nleaves*2)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

P5model$initializeInfo()

# compiled_P5model <- compileNimble(P5model,
#                                   showCompilerOutput = TRUE)

Nchains = 2
mcmc.out <- nimbleMCMC(code = P5model,
                       constants = Constants,
                       monitors = c("Nmean","Cabmean","Carmean","Cwmean","Cmmean","Standard.Dev",
                                    "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD"),
                       data = Data,
                       inits = Inits,
                       nburnin = 10000,
                       nchains = Nchains,
                       niter = 100000,
                       summary = TRUE,
                       WAIC = FALSE,
                       samplesAsCodaMCMC = TRUE)

MCMCsamples <- mcmc.out$samples

param <- MCMCsamples[,1:12]
plot(param[,1])
# pairs(as.matrix(param[,c(1:5)]), pch = '.')

hist(param[[1]][,7])

Nsimu <- 1000

if (Nchains == 1){
  pos.simu <- sample(1:nrow(MCMCsamples),Nsimu)
  param_all <- MCMCsamples[pos.simu,1:12]
} else {
  pos.simu <- sample(1:nrow(MCMCsamples[[1]]),Nsimu)
  param_all <- do.call(rbind,lapply(1:Nchains,function(i) MCMCsamples[[i]][pos.simu,1:12]))
}

Simu <- cbind(matrix(unlist(lapply(1:Nsimu,function(ileaf){
  rrtm::prospect5(N = param_all[ileaf,5],
                  Cab = param_all[ileaf,1],
                  Car = param_all[ileaf,2],
                  Cbrown = 0,
                  Cw = param_all[ileaf,4],
                  Cm = param_all[ileaf,3])[["reflectance"]][pos]})),ncol = Nsimu),
  matrix(unlist(lapply(1:Nsimu,function(ileaf){
    rrtm::prospect5(N = param_all[ileaf,5] + param_all[ileaf,11],
                    Cab = param_all[ileaf,1] + param_all[ileaf,7],
                    Car = param_all[ileaf,2] + param_all[ileaf,8],
                    Cbrown = 0,
                    Cw = param_all[ileaf,4] + param_all[ileaf,10],
                    Cm = param_all[ileaf,3] + param_all[ileaf,9])[["reflectance"]][pos]})),ncol = Nsimu))

matplot(WLs,Simu[,1:Nsimu],type = 'l',col = "darkgreen",ylim = c(0,max(c(Data$obs_reflectance,Simu))*1.1))
matlines(WLs,Simu[,Nsimu + (1:Nsimu)],type = 'l',col = "darkblue")

matlines(WLs,rowMeans(Data$obs_reflectance[,1:10]),col = "red",lwd = 2)
matlines(WLs,rowMeans(Data$obs_reflectance[,10 + (1:10)]),col = "red",lwd = 2)


# matlines(WLs,Data$obs_reflectance[,1:Nleaves],type = 'l',col = "darkgreen")
# matlines(WLs,Data$obs_reflectance[,Nleaves + (1:Nleaves)],type = 'l',col = "darkblue")


df_l <- data.frame(mcmc.out$samples$chain1) %>% rename(N = Nmean,
                                                       Cab = Cabmean,
                                                       Car = Carmean,
                                                       Cw = Cwmean,
                                                       Cm = Cmmean) %>% dplyr::select(c("N","Cab","Car","Cw","Cm","Standard.Dev",
                                                                                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD")) %>%
  pivot_longer(cols = c("N","Cab","Car","Cw","Cm","Standard.Dev",
                        "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm","alpha_SD"),
               names_to = "param",
               values_to = "value")

params_true <- data.frame(N = mean(Nall),
                          Cab = mean(Caball),
                          Car = mean(Carall),
                          Cw = mean(Cwall),
                          Cm = mean(Cmall),
                          alpha_N = mean(Nall_L - Nall),
                          alpha_Cab = mean(Caball_L - Caball),
                          alpha_Car = mean(Carall_L - Carall),
                          alpha_Cw = mean(Cwall_L - Cwall),
                          alpha_Cm = mean(Cmall_L - Cmall)) %>% pivot_longer(cols = c("N","Cab","Car","Cw","Cm",
                                                                                      "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm"),
                                                                             names_to = "param",
                                                                             values_to = "value")

params_all <- data.frame(N = Nall,
                         Cab = Caball,
                         Car = Carall,
                         Cw = Cwall,
                         Cm = Cmall,
                         alpha_N = Nall_L - Nall,
                         alpha_Cab = Caball_L - Caball,
                         alpha_Car = Carall_L - Carall,
                         alpha_Cw = Cwall_L - Cwall,
                         alpha_Cm = Cmall_L - Cmall) %>% pivot_longer(cols = c("N","Cab","Car","Cw","Cm",
                                                                               "alpha_N","alpha_Cab","alpha_Car","alpha_Cw","alpha_Cm"),
                                                                      names_to = "param",
                                                                      values_to = "value")

ggplot(df_l,aes(value)) +
  geom_histogram(aes(y= ..density..),bins = 60, alpha = 0.4, color = "darkgrey") +
  geom_point(data = params_true,
             aes(x = value,y = 0), color = 'red') +
  geom_point(data = params_all,
             aes(x = value,y = 0), color = 'red',shape = "|",size = 4) +
  facet_wrap(~param, scales = "free") +
  theme_bw()


