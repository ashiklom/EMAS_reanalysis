### Hierarchical analysis
library(nimble)
library(coda)
source("format_data.R")

### Ca/Al code
CaAl.code <- nimbleCode({
        mu ~ dlnorm(3, 0.25)
        tau ~ dgamma(0.01, 0.01)
        tau.sample ~ dgamma(0.01, 0.01)
        tau.flux ~ dgamma(0.01, 0.01)
        tau.month ~ dgamma(0.01, 0.01)
        for(s in 1:n.sample) {alpha.sample[s] ~ dnorm(0, tau.sample)}
        for(f in 1:n.flux) {alpha.flux[f] ~ dnorm(0, tau.flux)}
        for(m in 1:n.month) {alpha.month[m] ~ dnorm(0, tau.month)}
        
        
        for(i in 1:n.value){
                Ex[i] <- mu +
                        alpha.flux[flux[i]] +
                        alpha.sample[sample[i]] +
                        alpha.month[Month[i]]
                Ca.Al[i] ~ dnorm(Ex[i], tau)
        }
})

setkey(chem.sub, Sample)
chem.spec <- chem.sub[c("BE", "YP")]
chem.spec <- chem.spec[complete.cases(chem.spec)]
CaAl.constants <- list(Ca.al = chem.spec[,Ca.Al],
                       sample = as.numeric(as.factor(chem.spec[,Sample])),
                       flux = as.numeric(as.factor(chem.spec[,Flux])),
                       Month = chem.spec[,Month],
                       n.value = nrow(chem.spec))
CaAl.constants <- within(CaAl.constants, {
        n.month = max(Month)
        n.sample = max(sample)
        n.flux = max(flux)
})

CaAl <- nimbleModel(CaAl.code, name = "CaAl", constants = CaAl.constants)
CaAl.spec <- configureMCMC(CaAl)
CaAl.spec$addMonitors(c("alpha.month", "alpha.sample", "alpha.flux"))
CaAl.MCMC <- buildMCMC(CaAl.spec, project=CaAl)
CaAl.proj <- compileNimble(CaAl, CaAl.MCMC)

n.iter <- 50000
n.chains <- 5
thin <- 10
burnin <- 40000
bt <- seq(burnin, n.iter, by=thin)
caal.mcmc <- list()
for(i in 1:n.chains){
        CaAl.proj$CaAl.MCMC$run(n.iter)
        caal.mcmc[[i]] <- mcmc(as.matrix(CaAl.proj$CaAl.MCMC$mvSamples)[bt,])
}
caal.mcmc <- mcmc.list(caal.mcmc)
save(caal.mcmc, file="Ca_Al_ncv.Rdata")

print("#############################")
print("COMPLETED CA_AL ANALYSIS. STARTING CA ANALYSIS")
print("#############################")

Ca.constants <- CaAl.constants
Ca.constants$Ca.Al <- chem.spec[,Ca]
Ca.model <- nimbleModel(CaAl.code, name = "Ca_model", constants = Ca.constants)
Ca.spec <- configureMCMC(Ca.model)
Ca.spec$addMonitors(c("alpha.month", "alpha.sample", "alpha.flux"))
Ca.MCMC <- buildMCMC(Ca.spec, project=Ca_model)
Ca.proj <- compileNimble(Ca.model, Ca.MCMC)
ca.mcmc <- list()
for(i in 1:n.chains){
        Ca.proj$Ca.MCMC$run(n.iter)
        ca.mcmc[[i]] <- mcmc(as.matrix(Ca.proj$Ca.MCMC$mvSamples)[bt,])
}
ca.mcmc <- mcmc.list(ca.mcmc)
save(ca.mcmc, file = "Ca_ncv.Rdata")

print("#############################")
print("COMPLETED CA ANALYSIS. STARTING AL ANALYSIS")
print("#############################")

Al.constants <- CaAl.constants
Al.constants$Ca.Al <- chem.spec[,Total.Al]
Al.model <- nimbleModel(CaAl.code, name = "Ca_model", constants = Al.constants)
Al.spec <- configureMCMC(Al.model)
Al.spec$addMonitors(c("alpha.month", "alpha.sample", "alpha.flux"))
Al.MCMC <- buildMCMC(Al.spec, project=Al_model)
Al.proj <- compileNimble(Al.model, Al.MCMC)
al.mcmc <- list()
for(i in 1:n.chains){
        Al.proj$Al.MCMC$run(n.iter)
        al.mcmc[[i]] <- mcmc(as.matrix(Al.proj$Al.MCMC$mvSamples)[bt,])
}
al.mcmc <- mcmc.list(al.mcmc)
save(al.mcmc, file = "Al_ncv.Rdata")

print("#############################")
print("COMPLETED AL ANALYSIS")
print("#############################")

