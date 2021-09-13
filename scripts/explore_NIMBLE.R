rm(list = ls())

# install.packages("nimble", repos = "http://r-nimble.org", type = "source")

library(nimble)

dyesCode <- nimbleCode({
  # Model
  for (i in 1:BATCHES) {
    for (j in 1:SAMPLES) {
      y[i,j] ~ dnorm(mu[i], sd = sigma.within);
    }
    mu[i] ~ dnorm(theta, sd = sigma.between);
  }

  # Priors
  theta ~ dnorm(0.0, 1.0E-10);
  sigma.within ~ dunif(0, 100)
  sigma.between ~ dunif(0, 100)
})

dyesModel <- nimbleModel(dyesCode,
                         constants = list(BATCHES = 6, SAMPLES = 5))
dyesModel$initializeInfo()

data <- matrix(c(1545, 1540, 1595, 1445, 1595, 1520, 1440, 1555, 1550,
                 1440, 1630, 1455, 1440, 1490, 1605, 1595, 1515, 1450, 1520, 1560,
                 1510, 1465, 1635, 1480, 1580, 1495, 1560, 1545, 1625, 1445), nrow = 6)

dyesModel$setData(list(y = data))
dyesModel$y

dyesModel$theta <- 1500
dyesModel$mu <- rnorm(6, 1500, 50)
dyesModel$sigma.within <- 20
dyesModel$sigma.between <- 20
dyesModel$y[1,]

dyesModel$calculate(c("y[1,2]","theta"))

dyesModel$simulate()
dyesModel$simulate(c('mu'))
dyesModel$mu
dyesModel$y
dyesModel$theta

library(igraph)

plot(dyesModel$getGraph())

compiled_dyesModel <- compileNimble(dyesModel,showCompilerOutput = TRUE)

compiled_dyesModel$theta <- 1450
compiled_dyesModel$calculate()

####################################################################################
rm(list = ls())

library(nimble)

timesTwo <- nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))
    return(2*x)
  })


e1_approx <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))

    A <- log((0.56146 / x + 0.65) * (1 + x))
    B <- x^4 * exp(7.7 * x) * (2 + x)^3.7

    return((A^-7.7 + B)^-0.13)
  })

code <- nimbleCode({
  for(i in 1:2101) {
    mu[i] ~ dgamma(1, 2)
    mu_e1_approx[i] <- e1_approx(mu[i])
  }

  mu2_e1_approx[1:2101] <- e1_approx(mu[1:2101])
  ralf[1:2101] <- 1 - talf[1:2101]

})


model <- nimbleModel(code)

model$setData(list(talf = rrtm:::p45_t12,
                   t12 <- rrtm:::p45_t12,
                   t21 <- rrtm:::p45_t21,
                   dataspec_p5 <- rrtm:::dataspec_p5))


model$simulate()

model$mu
model$mu_e1_approx
model$mu2_e1_approx
rrtm:::e1_approx(x = model$mu)
model$ralf
