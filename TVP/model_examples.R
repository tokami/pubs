## Example applications of time-variant model extensions of spict
## March 2019


## install and load required spict branch (uncomment if not installed yet)
## ------------------------------------------------------------
## install.packages("ellipse")
## install.packages("devtools")
## devtools::install_github("tokami/spict/spict", ref="seaprodTVP", force=TRUE)
library(spict)


## Model applications
## ------------------------------------------------------------



## 1. Original spict
## ------------------------------------------------------------

## simulate data
## ------------------------------------------------------------
nyears <- 30
logsdb <- log(0.1)
logsdi <- log(0.1)
logsdf <- log(0.1)
logsdc <- log(0.1)
timeIlist <- list(seq(1/8, nyears, by=1),
                  seq(3/8,nyears,by=1),
                  seq(5/8,nyears,by=1),
                  seq(7/8,nyears,by=1))
logphi <- log(c(2.5,1.2,0.34))  ## identical pattern: log(c(1.7,0.8,0.5))
simPhase <- pi * 1.28                              ## timing of seasonal production
simAmp <- 0.73
## not changing in scenarios
dteuler <- 1/8
nseaC <- 4
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
logqi <- log(c(0.03,0.04,0.02,0.03))
logbkfraci <- log(0.95)
nt <- length(seq(0,nyears+0.25,dteuler))
## Albacore
logKi <- log(201.48)   
logmi <- log(c(22.58, 20))
logni <- log(0.69)
psi <- log(0.05)
sdm <- log(0.01)
nt <- length(seq(0,nyears+0.25,dteuler))
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
## set up simulation parameters
inp <- list(nseasons=4,
            seasontype = 1,
            splineorder=3,                    
            dteuler = dteuler)
## time of catches and indices
inp$timeC <- timeC
inp$timeI <- timeIlist
## initial parameters
inp$ini <- list(logK=logKi,
                logm=logmi,
                logn=logni,
                logq=logqi,
                logbkfrac=logbkfraci,
                logsdb=logsdb,
                logsdi=logsdi,
                logsdf=logsdf, 
                logsdc=logsdc,
                logphi=logphi)
## additional settings
inp$msytype <- "d"             ## for n < 1 stochastic have not been proven to be correct        
inp$ini$logF0 <- log(0.01)
inp$Fpattern <- 3
inp$Fmax <- 2 

## seasonal (OFF)
inp$seaprod <- 0
inp$simlogm <- logmi
inp$ampSP <- simAmp
inp$phaseSP <- simPhase

## regime shift
nt7 <- length(seq(0,7,dteuler))    
inp$MSYregime <- as.factor(c(rep(1,nt)))

## simulate
inpori <- spict::sim.spictSPRSTVG(inp)


## fit model
## ------------------------------------------------------------
inp <- inpori

## seaprod
inp$optimiser.control <- list("iter.max" = 10000,
                              "eval.max" = 10000)

inp$msytype <- "d"

inp$priors$logn <- c(0,0,0)
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)

fit <- spict::fit.spict(inp)

plot(fit)




## 2. Seasonal spict
## ------------------------------------------------------------

## simulate data
## ------------------------------------------------------------
nyears <- 30
logsdb <- log(0.1)
logsdi <- log(0.1)
logsdf <- log(0.1)
logsdc <- log(0.1)
timeIlist <- list(seq(1/8, nyears, by=1),
                  seq(3/8,nyears,by=1),
                  seq(5/8,nyears,by=1),
                  seq(7/8,nyears,by=1))
logphi <- log(c(2.5,1.2,0.34))  ## identical pattern: log(c(1.7,0.8,0.5))
simPhase <- pi * 1.28                              ## timing of seasonal production
simAmp <- 0.73
## not changing in scenarios
dteuler <- 1/8
nseaC <- 4
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
logqi <- log(c(0.03,0.04,0.02,0.03))
logbkfraci <- log(0.95)
nt <- length(seq(0,nyears+0.25,dteuler))
## Albacore
logKi <- log(201.48)   
logmi <- log(c(22.58, 20))
logni <- log(0.69)
psi <- log(0.05)
sdm <- log(0.01)
nt <- length(seq(0,nyears+0.25,dteuler))
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
## set up simulation parameters
inp <- list(nseasons=4,
            seasontype = 1,
            splineorder=3,                    
            dteuler = dteuler)
## time of catches and indices
inp$timeC <- timeC
inp$timeI <- timeIlist
## initial parameters
inp$ini <- list(logK=logKi,
                logm=logmi,
                logn=logni,
                logq=logqi,
                logbkfrac=logbkfraci,
                logsdb=logsdb,
                logsdi=logsdi,
                logsdf=logsdf, 
                logsdc=logsdc,
                logphi=logphi)
## additional settings
inp$msytype <- "d"             ## for n < 1 stochastic have not been proven to be correct        
inp$ini$logF0 <- log(0.01)
inp$Fpattern <- 3
inp$Fmax <- 2 

## seasonal 
inp$seaprod <- 3
inp$simlogm <- logmi
inp$ampSP <- simAmp
inp$phaseSP <- simPhase

## regime shift
nt7 <- length(seq(0,7,dteuler))    
inp$MSYregime <- as.factor(c(rep(1,nt)))

## simulate
inpori <- spict::sim.spictSPRSTVG(inp)


## fit model
## ------------------------------------------------------------
inp <- inpori

## seaprod
inp$optimiser.control <- list("iter.max" = 10000,
                              "eval.max" = 10000)

inp$msytype <- "d"

inp$priors$logn <- c(0,0,0)
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)

fit <- spict::fit.spict(inp)

## plot results
plot(fit)

## seasonal graphs
opar <- par(mfrow=c(2,1))
spict::plotspict.season(fit)
spict::plotspict.seaprod(fit)
par(opar)






## 3. Regime-shift spict
## ------------------------------------------------------------

## simulate data
## ------------------------------------------------------------
nyears <- 30
logsdb <- log(0.1)
logsdi <- log(0.1)
logsdf <- log(0.1)
logsdc <- log(0.1)
timeIlist <- list(seq(1/8, nyears, by=1),
                  seq(3/8,nyears,by=1),
                  seq(5/8,nyears,by=1),
                  seq(7/8,nyears,by=1))
logphi <- log(c(2.5,1.2,0.34))  ## identical pattern: log(c(1.7,0.8,0.5))
simPhase <- pi * 1.28                              ## timing of seasonal production
simAmp <- 0.73
## not changing in scenarios
dteuler <- 1/8
nseaC <- 4
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
logqi <- log(c(0.03,0.04,0.02,0.03))
logbkfraci <- log(0.95)
nt <- length(seq(0,nyears+0.25,dteuler))
## Albacore
logKi <- log(201.48)   
logmi <- log(c(22.58, 20))
logni <- log(0.69)
psi <- log(0.05)
sdm <- log(0.01)
nt <- length(seq(0,nyears+0.25,dteuler))
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
## set up simulation parameters
inp <- list(nseasons=4,
            seasontype = 1,
            splineorder=3,                    
            dteuler = dteuler)
## time of catches and indices
inp$timeC <- timeC
inp$timeI <- timeIlist
## initial parameters
inp$ini <- list(logK=logKi,
                logm=logmi,
                logn=logni,
                logq=logqi,
                logbkfrac=logbkfraci,
                logsdb=logsdb,
                logsdi=logsdi,
                logsdf=logsdf, 
                logsdc=logsdc,
                logphi=logphi)
## additional settings
inp$msytype <- "d"             ## for n < 1 stochastic have not been proven to be correct        
inp$ini$logF0 <- log(0.01)
inp$Fpattern <- 3
inp$Fmax <- 2 

## seasonal (OFF)
inp$seaprod <- 0
inp$simlogm <- logmi
inp$ampSP <- simAmp
inp$phaseSP <- simPhase

## regime shift
nt7 <- length(seq(0,7,dteuler))    
inp$MSYregime <- as.factor(c(rep(1,nt7), rep(2,nt-nt7)))

## simulate
inpori <- spict::sim.spictSPRSTVG(inp)


## fit model
## ------------------------------------------------------------
inp <- inpori

## seaprod
inp$optimiser.control <- list("iter.max" = 10000,
                              "eval.max" = 10000)

inp$msytype <- "d"

inp$priors$logn <- c(0,0,0)
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)

fit <- spict::fit.spict(inp)


## plot results
plot(fit)






## 4. Seasonal regime-shift spict
## ------------------------------------------------------------

## simulate data
## ------------------------------------------------------------
nyears <- 30
logsdb <- log(0.1)
logsdi <- log(0.1)
logsdf <- log(0.1)
logsdc <- log(0.1)
timeIlist <- list(seq(1/8, nyears, by=1),
                  seq(3/8,nyears,by=1),
                  seq(5/8,nyears,by=1),
                  seq(7/8,nyears,by=1))
logphi <- log(c(2.5,1.2,0.34))  ## identical pattern: log(c(1.7,0.8,0.5))
simPhase <- pi * 1.28                              ## timing of seasonal production
simAmp <- 0.73
## not changing in scenarios
dteuler <- 1/8
nseaC <- 4
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
logqi <- log(c(0.03,0.04,0.02,0.03))
logbkfraci <- log(0.95)
nt <- length(seq(0,nyears+0.25,dteuler))
## Albacore
logKi <- log(201.48)   
logmi <- log(c(22.58, 20))
logni <- log(0.69)
psi <- log(0.05)
sdm <- log(0.01)
nt <- length(seq(0,nyears+0.25,dteuler))
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
## set up simulation parameters
inp <- list(nseasons=4,
            seasontype = 1,
            splineorder=3,                    
            dteuler = dteuler)
## time of catches and indices
inp$timeC <- timeC
inp$timeI <- timeIlist
## initial parameters
inp$ini <- list(logK=logKi,
                logm=logmi,
                logn=logni,
                logq=logqi,
                logbkfrac=logbkfraci,
                logsdb=logsdb,
                logsdi=logsdi,
                logsdf=logsdf, 
                logsdc=logsdc,
                logphi=logphi)
## additional settings
inp$msytype <- "d"             ## for n < 1 stochastic have not been proven to be correct        
inp$ini$logF0 <- log(0.01)
inp$Fpattern <- 3
inp$Fmax <- 2 

## seasonal 
inp$seaprod <- 3
inp$simlogm <- logmi
inp$ampSP <- simAmp
inp$phaseSP <- simPhase

## regime shift
nt7 <- length(seq(0,7,dteuler))    
inp$MSYregime <- as.factor(c(rep(1,nt7), rep(2,nt-nt7)))

## simulate
inpori <- spict::sim.spictSPRSTVG(inp)


## fit model
## ------------------------------------------------------------
inp <- inpori

## seaprod
inp$optimiser.control <- list("iter.max" = 10000,
                              "eval.max" = 10000)

inp$msytype <- "d"

inp$priors$logn <- c(0,0,0)
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)

fit <- spict::fit.spict(inp)

## plot results
plot(fit)

## seasonal graphs
opar <- par(mfrow=c(2,1))
spict::plotspict.season(fit)
spict::plotspict.seaprod(fit)
par(opar)





## 5. Gradual shift spict
## ------------------------------------------------------------

## simulate data
## ------------------------------------------------------------
nyears <- 30
logsdb <- log(0.1)
logsdi <- log(0.1)
logsdf <- log(0.1)
logsdc <- log(0.1)
timeIlist <- list(seq(1/8, nyears, by=1),
                  seq(3/8,nyears,by=1),
                  seq(5/8,nyears,by=1),
                  seq(7/8,nyears,by=1))
logphi <- log(c(2.5,1.2,0.34))  ## identical pattern: log(c(1.7,0.8,0.5))
simPhase <- pi * 1.28                              ## timing of seasonal production
simAmp <- 0.73
dteuler <- 1/8
nseaC <- 4
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
logqi <- log(c(0.03,0.04,0.02,0.03))
logbkfraci <- log(0.95)
nt <- length(seq(0,nyears+0.25,dteuler))
## Albacore
logKi <- log(201.48)   
logmi <- log(c(22.58, 20))
logni <- log(0.69)
logpsi <- log(0.05)
logsdm <- log(0.01)
nt <- length(seq(0,nyears+0.25,dteuler))
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
## set up simulation parameters
inp <- list(nseasons=4,
            seasontype = 1,
            splineorder = 3,                    
            dteuler = dteuler)
## time of catches and indices
inp$timeC <- timeC
inp$timeI <- timeIlist
## initial parameters
inp$ini <- list(logK=logKi,
                logm=logmi,
                logn=logni,
                logq=logqi,
                logbkfrac=logbkfraci,
                logsdb=logsdb,
                logsdi=logsdi,
                logsdf=logsdf, 
                logsdc=logsdc,
                logphi=logphi)
## additional settings
inp$msytype <- "d"             ## for n < 1 stochastic have not been proven to be correct        
inp$ini$logF0 <- log(0.01)
inp$Fpattern <- 3
inp$Fmax <- 2 

## seasonal (OFF)
inp$seaprod <- 0
inp$simlogm <- logmi
inp$ampSP <- simAmp
inp$phaseSP <- simPhase


## time varying growth/productivity
inp$ini$logpsi <- logpsi
inp$ini$logsdm <- logsdm
inp$timevaryinggrowth <- TRUE
inp$simlogmre0 <- logmi + 0.5
inp$MREpattern <- 0
inp$simlogmre <- seq(-0.5, 0, length.out = nt)


## simulate
inpori <- spict::sim.spictSPRSTVG(inp)


## fit model
## ------------------------------------------------------------
inp <- inpori

## seaprod
inp$optimiser.control <- list("iter.max" = 10000,
                              "eval.max" = 10000)

inp$msytype <- "d"

inp$priors$logn <- c(0,0,0)
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)

fit <- spict::fit.spict(inp)


## plot results
plot(fit)






## 4. Seasonal gradual shift spict
## ------------------------------------------------------------

## simulate data
## ------------------------------------------------------------
nyears <- 30
logsdb <- log(0.1)
logsdi <- log(0.1)
logsdf <- log(0.1)
logsdc <- log(0.1)
timeIlist <- list(seq(1/8, nyears, by=1),
                  seq(3/8,nyears,by=1),
                  seq(5/8,nyears,by=1),
                  seq(7/8,nyears,by=1))
logphi <- log(c(2.5,1.2,0.34))  ## identical pattern: log(c(1.7,0.8,0.5))
simPhase <- pi * 1.28                              ## timing of seasonal production
simAmp <- 0.73
## not changing in scenarios
dteuler <- 1/8
nseaC <- 4
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
logqi <- log(c(0.03,0.04,0.02,0.03))
logbkfraci <- log(0.95)
nt <- length(seq(0,nyears+0.25,dteuler))
## Albacore
logKi <- log(201.48)   
logmi <- log(c(22.58, 20))
logni <- log(0.69)
psi <- log(0.05)
sdm <- log(0.01)
nt <- length(seq(0,nyears+0.25,dteuler))
timeC <- seq(0, nyears - 1/nseaC, by=1/nseaC)
## set up simulation parameters
inp <- list(nseasons=4,
            seasontype = 1,
            splineorder=3,                    
            dteuler = dteuler)
## time of catches and indices
inp$timeC <- timeC
inp$timeI <- timeIlist
## initial parameters
inp$ini <- list(logK=logKi,
                logm=logmi,
                logn=logni,
                logq=logqi,
                logbkfrac=logbkfraci,
                logsdb=logsdb,
                logsdi=logsdi,
                logsdf=logsdf, 
                logsdc=logsdc,
                logphi=logphi)
## additional settings
inp$msytype <- "d"             ## for n < 1 stochastic have not been proven to be correct        
inp$ini$logF0 <- log(0.01)
inp$Fpattern <- 3
inp$Fmax <- 2 

## seasonal 
inp$seaprod <- 3
inp$simlogm <- logmi
inp$ampSP <- simAmp
inp$phaseSP <- simPhase

## time varying growth/productivity
inp$ini$logpsi <- logpsi
inp$ini$logsdm <- logsdm
inp$timevaryinggrowth <- TRUE
inp$simlogmre0 <- logmi + 0.5
inp$MREpattern <- 0
inp$simlogmre <- seq(-0.5, 0, length.out = nt)

## simulate
inpori <- spict::sim.spictSPRSTVG(inp)


## fit model
## ------------------------------------------------------------
inp <- inpori

## seaprod
inp$optimiser.control <- list("iter.max" = 10000,
                              "eval.max" = 10000)

inp$msytype <- "d"

inp$priors$logn <- c(0,0,0)
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)

fit <- spict::fit.spict(inp)

## plot results
plot(fit)

## seasonal graphs
opar <- par(mfrow=c(2,1))
spict::plotspict.season(fit)
spict::plotspict.seaprod(fit)
par(opar)
