rm(list = ls())

library(nimble)
library(igraph)
library(coda)
library(dplyr)

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
    # cc[,] <- t(c(Cab, Car, Cbrown, Cw, Cm) / N)
    cc[1,1] <- Cab / N
    cc[2,1] <- Car / N
    cc[3,1] <- 0 / N
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
    # print(c(N,mean(reflectance)))

    returnType(double(1))

    print(c(N,Cab, Car, Cbrown, Cw, Cm,mean(reflectance)))
    return(reflectance)
  })

run_prospect5 <- nimbleCode({


  for (i in 1:Nsamples){


    reflectance[,i] <- NIMprospect5(Nmean,Cabmean,Carmean,Cwmean,Cmmean,
                                    dataspec_p5[,], talf[],t12[],t21[], Nwl)

    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(reflectance[j,i], sd = 0.1)
      # obs_reflectance[j,i] ~ dunif(0.99*reflectance[j,i], 1.01*reflectance[j,i])
    }
  }

  Nmean ~ dunif(1.1,5)
  Cabmean ~ dunif(10,50)
  Carmean ~ dunif(10,30)
  Cwmean ~ dunif(0.0001,0.1)
  Cmmean ~ dunif(0.0001,0.1)
})

Nleaves <- 1
WLa <- 400
WLb <- 2500
Delta_WL <- 100
WLs <- seq(WLa,WLb,Delta_WL)
pos <- which(400:2500 %in% WLs)
Nwl <- length(pos)


Constants <- list(Nsamples = Nleaves,
                  Nwl = Nwl,
                  talf = rrtm:::p45_talf[pos],
                  t12 = rrtm:::p45_t12[pos],
                  t21 = rrtm:::p45_t21[pos],
                  dataspec_p5 = rrtm:::dataspec_p5[pos,1:5])


Ntrue <- 2 ; Cabtrue <- 20 ; Cartrue <- 20 ; Cwtrue <- 0.02 ; Cmtrue <- 0.02

dataspec_p5 = rrtm:::dataspec_p5[pos,1:5] ; talf = rrtm:::p45_talf[pos] ; t12 = rrtm:::p45_t12[pos] ; t21 = rrtm:::p45_t21[pos]

fac <- 100
Nall <- rnorm(Nleaves,Ntrue,Ntrue/fac)
Caball <- rnorm(Nleaves,Cabtrue,Cabtrue/fac)
Carall <- rnorm(Nleaves,Cartrue,Cartrue/fac)
Cwall <- rnorm(Nleaves,Cwtrue,Cwtrue/fac)
Cmall <- rnorm(Nleaves,Cmtrue,Cmtrue/fac)

Data <- list(obs_reflectance = matrix(unlist(lapply(1:length(Nall),function(ileaf){
  rrtm::prospect5(N = Nall[ileaf],
                  Cab = Caball[ileaf],
                  Car = Carall[ileaf],
                  Cbrown = 0,
                  Cw = Cwall[ileaf],
                  Cm = Cmall[ileaf])[["reflectance"]][pos]})),ncol = Nleaves))

matplot(WLs,Data$obs_reflectance,type = 'l')

Inits <- list(Nmean = Ntrue,
              Cabmean = Cabtrue,
              Carmean = Cartrue,
              Cwmean = Cwtrue,
              Cmmean = Cmtrue)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl,Nleaves)),
                       data = Data,
                       constants = Constants,
                       debug = FALSE,
                       inits = Inits)

compiled_P5model <- compileNimble(P5model,
                                  showCompilerOutput = TRUE)


compiled_P5model$simulate()


dumb <- NIMprospect5(N = compiled_P5model$Nmean,Cab = compiled_P5model$Cabmean,Car = compiled_P5model$Carmean, Cw = compiled_P5model$Cwmean, Cm = compiled_P5model$Cmmean,
                     dataspec_p5 = rrtm:::dataspec_p5[pos,1:5],talf = rrtm:::p45_talf[pos],t12 = rrtm:::p45_t12[pos],t21 = rrtm:::p45_t21[pos], Nwl = Nwl)
