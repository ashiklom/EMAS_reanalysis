# Process resutls from Bayesian analysis
library(coda)
load("Ca_Al.Rdata")
load("Ca.Rdata")
load("Al.Rdata")
plot(caal.mcmc)
plot(ca.mcmc)
plot(al.mcmc)
