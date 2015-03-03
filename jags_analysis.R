### Hierarchical analysis
library(rjags)
library(coda)
source("format_data.R")

burnin <- 80000
thin <- 20
n.chain <- 5
n.iter <- 5000*thin/n.chain 

### Ca/Al code
CaAl.code <- "
model{
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
}
"

setkey(chem.sub, Sample)
chem.spec <- chem.sub[c("BE", "YP")]
chem.spec <- chem.spec[complete.cases(chem.spec)]
CaAl.constants <- list(Ca.Al = chem.spec[,Ca.Al],
                       sample = as.numeric(as.factor(chem.spec[,Sample])),
                       flux = as.numeric(as.factor(chem.spec[,Flux])),
                       Month = chem.spec[,Month],
                       n.value = nrow(chem.spec))
CaAl.constants <- within(CaAl.constants, {
        n.month = max(Month)
        n.sample = max(sample)
        n.flux = max(flux)
})

CaAl.model <- jags.model(textConnection(CaAl.code), data=CaAl.constants, n.chains = n.chain)
update(CaAl.model, n.iter=burnin)
monitors <- c("mu", "tau", "tau.sample", "tau.flux", "tau.month",
              "alpha.sample", "alpha.flux", "alpha.month")
caal.mcmc <- coda.samples(model = CaAl.model,
                             variable.names = monitors,
                             n.iter = n.iter,
                             thin = thin)

print("#############################")
print("COMPLETED CA_AL ANALYSIS. STARTING CA ANALYSIS")
print("#############################")

Ca.constants <- CaAl.constants
Ca.constants$Ca.Al <- chem.spec[,Ca]
Ca.model <- jags.model(textConnection(CaAl.code), data=Ca.constants, n.chains = 5)
update(Ca.model, n.iter=burnin)
ca.mcmc <- coda.samples(model = CaAl.model,
                          variable.names = monitors,
                          n.iter = n.iter,
                          thin = thin)

print("#############################")
print("COMPLETED CA ANALYSIS. STARTING AL ANALYSIS")
print("#############################")

Al.constants <- CaAl.constants
Al.constants$Ca.Al <- chem.spec[,Total.Al]
Al.model <- jags.model(textConnection(CaAl.code), data=Al.constants, n.chains = 5)
update(Al.model, n.iter=burnin)
al.mcmc <- coda.samples(model = Al.model,
                        variable.names = monitors,
                        n.iter = n.iter,
                        thin = thin)

print("#############################")
print("COMPLETED AL ANALYSIS")
print("#############################")

