rm(list = ls())

library(nimble)
library(igraph)
library(coda)

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
    cc[,] <- t(c(Cab, Car, Cbrown, Cw, Cm) / N)

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
    return(reflectance)
  })


run_prospect5 <- nimbleCode({

  Nmean ~ dunif(1.1,5)

  for (i in 1:Nsamples){

    N[i] ~ dnorm(Nmean,1/(0.1*Nmean)**2)

    reflectance[,i] <- prospect5(N[i],Cab,Car,Cw,Cm,dataspec_p5[,], talf[],t12[],t21[], Nwl)

    for (j in 1:Nwl){
      obs_reflectance[j,i] ~ dnorm(reflectance[j,i],1/(0.001**2))
    }
  }
})

Nleaves <- 10
Nwl <- 2101
Constants <- list(Nsamples = Nleaves,
                  Nwl = Nwl,
                  talf = rrtm:::p45_talf[1:Nwl],
                  t12 = rrtm:::p45_t12[1:Nwl],
                  t21 = rrtm:::p45_t21[1:Nwl],
                  dataspec_p5 = rrtm:::dataspec_p5[1:Nwl,1:5])

Ntrue <- 3 ; Cabtrue <- 40 ; Cartrue <- 20 ; Cwtrue <- 0.002 ; Cmtrue <- 0.001

Nall <- rnorm(Nleaves,Ntrue,Ntrue*0.1)

Data <- list(obs_reflectance = matrix(unlist(reflectance <- lapply(1:length(Nall),function(ileaf){
  rrtm::prospect5(N = Nall[ileaf],Cab = Cabtrue,Car = Cartrue,Cbrown = 0,Cw = Cwtrue,Cm = Cmtrue)[["reflectance"]][1:Nwl]
  })),ncol = Nleaves))

matplot(Data$obs_reflectance,type = 'l')

# N <- 3 ; Cab <- 40 ; Car <- 20 ; Cw <- 0.002 ; Cm <- 0.001 ; Nmean = 5
# talf = rrtm:::p45_talf[1:Nwl] ; t12 = rrtm:::p45_t12[1:Nwl] ; t21 = rrtm:::p45_t21[1:Nwl] ; dataspec_p5 = rrtm:::dataspec_p5[1:Nwl,1:5]

Inits <- list(Cab = Cabtrue,
              Car = Cartrue,
              Cm = Cmtrue,
              Cw = Cwtrue,
              Nmean = 5)

P5model <- nimbleModel(run_prospect5,
                       dimensions = list(dataspec_p5 = c(Nwl,5),
                                         talf = Nwl,
                                         t12 = Nwl,
                                         t21 = Nwl,
                                         reflectance = c(Nwl,Nleaves)),
                       constants = Constants,
                       data = Data,
                       inits = Inits)

P5model$simulate()
P5model$calculate()

compiled_P5model <- compileNimble(P5model,
                                  showCompilerOutput = TRUE)

configMCMC <- configureMCMC(model = compiled_P5model,
                            constants = Constants,
                            data = Data,
                            inits = Inits,
                            nchains = 1, niter = 1000,
                            summary = TRUE,WAIC = TRUE)

P5modelMCMC <- buildMCMC(configMCMC)
Cobj <- compileNimble(P5model, P5modelMCMC)

niter <- 100000
set.seed(0)
Cobj$P5modelMCMC$run(niter)

MCMCsamples <- as.matrix(Cobj$P5modelMCMC$mvSamples)

# pairs(MCMCsamples[,1], pch = '.')
hist(MCMCsamples[,1])

# box = list( list(c('Nmean'), c(1.1, 5)))
#
# prospectMCEM <- buildMCEM(model = compiled_P5model,
#                           latentNodes = c('N[1:3]'),
#                           boxConstraints = box)
# prospectMLE <- prospectMCEM$run(initM = 1001)
# prospectMLE
