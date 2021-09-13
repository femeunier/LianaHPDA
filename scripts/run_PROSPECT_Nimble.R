rm(list = ls())

library(nimble)
library(igraph)

e1_approx <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))

    A <- log((0.56146 / x + 0.65) * (1 + x))
    B <- x^4 * exp(7.7 * x) * (2 + x)^3.7

    return((A^-7.7 + B)^-0.13)
  })

prospect5 <- nimbleFunction(
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
    # print(c(N,Cab,Cm,Cw,Car,mean(Ra),mean(Ta),mean(Rsub),mean(tt),mean(denomy),mean(reflectance)))

    returnType(double(1))
    return(reflectance)
  })

run_prospect5 <- nimbleCode({

  for (i in 1:Nsamples){

    N[i] ~ dunif(1.1,5)
    Cab[i] ~ dunif(30,50)
    Car[i] ~ dunif(5,25)
    Cw[i] ~ dunif(0.0001,0.1)
    Cm[i] ~ dunif(0.0001,0.1)

    reflectance[,i] <- prospect5(N[i],Cab[i],Car[i],Cw[i],Cm[i],dataspec_p5[,], talf[],t12[],t21[], Nwl)
  }
})

Nleaves <- 3
Nwl <- 2101
Constants <- list(Nwl = Nwl,
                  Nsamples = Nleaves,
                  talf = rrtm:::p45_talf[1:Nwl],
                  t12 = rrtm:::p45_t12[1:Nwl],
                  t21 = rrtm:::p45_t21[1:Nwl],
                  dataspec_p5 = rrtm:::dataspec_p5[1:Nwl,1:5])

reflectance <- rrtm::prospect5(N = 2,Cab = 50,Car = 10,Cbrown = 0,Cw = 0.002,Cm = 0.001)[["reflectance"]][1:Nwl]
Data <- list(reflectance = matrix(rep(reflectance,Nleaves),ncol = Nleaves) * cbind(rep(runif(1,0.5,1.5),Nwl),
                                                                                   rep(runif(1,0.5,1.5),Nwl),
                                                                                   rep(runif(1,0.5,1.5),Nwl)))

Inits <- list(Cab = rep(40,Nleaves),
              Car = rep(10,Nleaves),
              Cm = rep(0.001,Nleaves),
              Cw = rep(0.002,Nleaves),
              N = rep(2,Nleaves))

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl,Nleaves)),
                       constants = Constants,
                       inits = Inits)

P5model$initializeInfo()

P5model$N <- runif(3,1.1,5)
P5model$Cab <- runif(3,30,50)
P5model$Car <- runif(3,5,25)
P5model$Cw <- runif(3,0.0001,0.1)
P5model$Cm <- runif(3,0.0001,0.1)

P5model$calculate()
P5model$simulate()
P5model$getLogProb()

matplot(P5model$reflectance,type = 'l',col = "black",ylim = c(0,1))
matlines(Data$reflectance,col = 'red',shape = 1,type = "l")

compiled_Pmodel <- compileNimble(P5model,showCompilerOutput = TRUE)

compiled_Pmodel$calculate()
compiled_Pmodel$simulate()
compiled_Pmodel$getLogProb()

matplot(P5model$reflectance,type = 'l')
i = 3
lines(rrtm::prospect5(N = P5model$N[i],
                      Cab = P5model$Cab[i],
                      Car = P5model$Car[i],
                      Cbrown = 0,
                      Cw = P5model$Cw[i],
                      Cm = P5model$Cm[i])[["reflectance"]],
      col = 'red')


all(compiled_Pmodel$simulate() == P5model$simulate())

