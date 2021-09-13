rm(list = ls())

library(rrtm)

# Inputs
N = 2; Cab = 40; Car = 10; Cw = 0.002; Cm = 0.0001

# Constants
Cbrown = 0
talf <- rrtm:::p45_talf
t12 <- rrtm:::p45_t12
t21 <- rrtm:::p45_t21
dataspec_p5 <- rrtm:::dataspec_p5
e1fun = rrtm:::e1_approx

# Computation
cc <- rbind(Cab, Car, Cbrown, Cw, Cm) / N
k <-  dataspec_p5 %*% cc

# GPM
gt0 <- k > 0
k0 <- k[gt0]
trans <- rep(1.0, length(k))
trans[gt0] <- (1 - k0) * exp(-k0) + k0 ^ 2 * e1fun(k0)

tlt0 <- trans < 0
if (any(tlt0)) {
  warning(
    sum(tlt0),
    " transmissivity values less than zero ",
    "were set to zero. ",
    "This happens under very high absorptivity values ",
    "(e.g. very high pigment concentrations) ",
    "and when using the fast exponential integral approximation. ",
    "The transmissivity is on the order of <= 1e-7 ",
    "so this warning can probably be ignored.",
    "However, if you want more precise results, you can try setting ",
    "`e1fun = gsl::expint_E1` to use a more precise approximation."
  )
  trans[tlt0] <- 0
}

# Reflectance/transmittance of a single layer
# Commented quantities are precalculated
## talf <- tav_abs(40, refractive)
ralf <- 1 - talf
## t12 <- tav_abs(90, refractive)
r12 <- 1 - t12
## t21 <- t12 / (refractive ^ 2)
r21 <- 1 - t21

denom <- 1 - r21 ^ 2 * trans ^ 2
Ta <- talf * trans * t21 / denom
Ra <- ralf + r21 * trans * Ta

tt <- t12 * trans * t21 / denom
rr <- r12 + r21 * trans * tt

Tsub <- numeric(2101)
Rsub <- numeric(2101)

gt1 <- rr + tt >= 1
tgt1 <- tt[gt1]
Tsub[gt1] <- tgt1 / (tgt1 + (1 - tgt1) * (N - 1))
Rsub[gt1] <- 1 - Tsub[gt1]

# Extremely high absorptivity leads to zero reflectance and transmittance
inf <- rr == 0 | tt == 0
Tsub[inf] <- 0
Rsub[inf] <- 0

# Reflectance/transmittance of N layers
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
## result <- matrix(NA_real_, nrow = 2101, ncol = 2,
##                  dimnames = list(NULL, c("reflectance", "transmittance")))
## result[, "reflectance"] <- Ra + Ta * Rsub * t / denomy
## result[, "transmittance"] <- Ta * Tsub / denomy

# NOTE: A matrix here might be better, but returning a list is much faster.
result <- list(
  reflectance = Ra + Ta * Rsub * tt / denomy,
  transmittance = Ta * Tsub / denomy
)
plot(result$reflectance,type = 'l')
