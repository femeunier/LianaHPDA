rm(list = ls())

library(rrtm)

gsl::expint_E1


N = 1.1; Cab = 100; Car = 10; Cw = 0.002; Cm = 0.0001

result1 <- rrtm::prospect5(N = N,Cab = Cab,Car = Car,Cbrown = 0,Cw = Cw,Cm = Cm)
result2 <- rrtm::prospect5(N = N,Cab = Cab,Car = Car,Cbrown = 0,Cw = Cw,Cm = Cm,e1fun = gsl::expint_E1)


plot(result1$transmittance,ylim = c(0,1),type = 'l')
lines(result2$transmittance,ylim = c(0,1),col = 'red')

plot((result1$reflectance-result2$reflectance),type = 'l')
plot((result1$transmittance-result2$transmittance),type = 'l')
