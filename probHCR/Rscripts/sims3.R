## stochHCR manuscript
## April 2020
## Tobias Mildenberger <t.k.mildenberger@gmail.com>

## I learnt:
## AC could have huge effect, should make additional secnario with high and low AC otherwise all get intermediate AC
## F process shouldn't be 0.001 this doesn't work well for ling 0.2 seems to work for ling, does it for the other species too?
## which scenario with changing AC? figure out
## send new scripts to server

## input <- as.character(Sys.getenv("INPUT"))

## server <- FALSE
server <- TRUE
checkLab <- ""


## load packages
##--------------------------------------------------------------------------------
## remotes::install_github("tokami/spict/spict", ref="manage4.0") ## version ‘1.3.0’
## remotes::install_github("DLMtool/DLMtool") ## version: ‘5.4.999’
library(spict)
library(DLMtool)


## Scenario settings
##--------------------------------------------------------------------------------
if(server){
    datdir <- "../data/"
    savedir <- "../res/"
    source("stochSimFuncs.R")
    input2 <- input <- as.character(commandArgs(TRUE))
##    input <- as.character(Sys.getenv("INPUT"))  ## scenario number
##    input2 <- unlist(strsplit(input, " "))
    specNum <- as.integer(input2[1])
    scenNum <- as.integer(input2[2])
    acNum <- 1 ## as.integer(input2[3])
    sigmaRNum <- 1 ## as.integer(input2[4])
    seedNum <- as.integer(input2[3])
    nCores <- as.numeric(Sys.getenv('LSB_DJOB_NUMPROC'))
    parallel <- TRUE
}else{
    datdir <- "../data/"
    savedir <- "../res/"
    source("stochSimFuncs.R")
    specNum <- 16
    scenNum <- 1
    acNum <- 1
    sigmaRNum <- 1
    seedNum <- 1
    nCores <- 6
    parallel <- FALSE
}


## All updated stocks
##--------------------------------------------------------------------------------
stockNames <- c("haddock", "herring", "lemonsole", "ling", "lobster",
                "megrim", "plaice", "pollack", "redmullet", "rosefish",
                "sandeels", "turbot", "whiteanglerfishC",
                "whiteanglerfishN", "whiting", "anchovy", "pilchard",
                "ling2","UNKNOWN","UNKNOWN","UNKNOWN")
nstocks <- length(stockNames)
species <- stockNames[specNum]

## Load data
##--------------------------------------------------------------------------------
load(paste0(datdir,"allOM2020May.RData"))
## select spec/stock
OM <- allOM[[specNum]]


## Baseline settings
##--------------------------------------------------------------------------------
## dteuler
dteuler <- 1/16  ## 1 ## 1/4 ## 1/8 ## 1/16
## SD of logn prior
npriorSD <- NA
## number of historic years to use in assessment
nhist <- 50 ## means to use all (only changed in scenario 6)
## assessment interval
amtint <- 1

## Noise levels
lowI <- c(0.29,0.29)     ## based on Wiedenmann 2015
midI <- rep(mean(c(0.29,0.63)),2)
highI <- c(0.63,0.63)    ## based on Wiedenmann 2015
lowC <- c(0.05,0.05)     ## based on Wiedenmann 2015
midC <- rep(mean(c(0.05,0.3)),2)
highC <- c(0.3,0.3)      ## based on Wiedenmann 2015

## Remove all gradients, additional noise
OMup <- tinyErr(OM)
OMup@Cbiascv <- 0
OMup@Ibiascv <- 0
OMup@beta <- c(1,1)
## inter-variability in M
OMup@Msd <- c(0,0)     ## c(0.15,0.15)
OMup@Esd <- c(0.15,0.15) ## c(0.15,0.15) ## Carruthers
OMup@qinc <- c(0,0)
OMup@qcv <- c(0,0)
OMup@Linfsd <- c(0,0)
OMup@Ksd <- c(0,0)

## standard effort too low to get good contrast
OMup@EffYears <- c(0,0.6,0.8,1)
OMup@EffLower <- c(0.2,1.6,1.7,1.2)
OMup@EffUpper <- c(0.4,2,2,1.5)

## check selectivity:
print(OMup@L5)
print(OMup@LFS)
print(OMup@isRel)

OMup@maxF <- 4
OMup@R0 <- 1e7
OMup@D <- c(0.05,0.1) ## c(0.05,0.1)
OMup@seed <- seq(1,1e4,by=51)[seedNum]  ## cannot be subsequent numbers because that is used in DLMtool for run_parallel
OMup@reps <- 1
OMup@proyears <- 20
OMup@nyears <- nyears <- 50
OMup@nsim <- nsim <- 50


## species
if(specNum == 16){
    ## Anchovy
    ## Engraulis encrasicolus
    ## Actinopterygii (ray-finned fishes) > Clupeiformes (Herrings) > Engraulidae (Anchovies) > Engraulinae
    ## Thorson 2014: Clupeiformes
    OMup@Perr <- rep(0.771,2)
    OMup@AC <- rep(0.456,2)
}else if(specNum == 1){
    ## Haddock
    ## Melanogrammus aeglefinus
    ## Actinopterygii (ray-finned fishes) > Gadiformes (Cods) > Gadidae (Cods and haddocks)
    ## Thorson 2014: Gadiformes
    OMup@Perr <- rep(0.747,2)
    OMup@AC <- rep(0.420,2)
}else if(specNum == 4){
    ## Ling
    ## Molva molva
    ## Actinopterygii (ray-finned fishes) > Gadiformes (Cods) > Lotidae (Hakes and burbots)
    ## Thorson 2014: Gadiformes
    OMup@Perr <- rep(0.747,2)
    OMup@AC <- rep(0.420,2)
}


## Scenarios for each species
##--------------------------------------------------------------------------------
if(scenNum == 1){ ## scenario 1: low scientific noise
    OMup@Cobs <- lowC
    OMup@Iobs <- lowI
}else if(scenNum == 2){ ## scenario 2: mid scientific noise
    OMup@Cobs <- midC
    OMup@Iobs <- midI
}else if(scenNum == 3){ ## scenario 3: high scentific noise
    OMup@Cobs <- highC
    OMup@Iobs <- highI
}else if(scenNum == 4){ ## scenario 4: shorter time series
    nhist <- 30
    OMup@Cobs <- midC
    OMup@Iobs <- midI
}else if(scenNum == 5){ ## scenario 5: higher recruitment variability
    OMup@Perr <- 1.2 * OMup@Perr
    OMup@Cobs <- midC
    OMup@Iobs <- midI
}


## ## Run historic simulations (needed for vulnerable biomass index)
## Hist <- runMSE(OMup, Hist=TRUE)
## ##
## OMup2 <- OMup
## # historical vulnerable biomass
## histVB <- Hist@TSdata$VB
## histVB2 <- apply(histVB, 2, function(x) exp(log(x) + rnorm(length(x), 0, sdI) - sdI^2/2))
## V <- Hist@AtAge$Select[,,OMup2@nyears] # selectivity for vuln biomass in last historical year
## ##
## Data <- new("Data")
## Data@AddInd <- array(histVB2, dim=c(OMup2@nsim, 1, OMup2@nyears))
## Data@CV_AddInd <- array(0, dim=c(OMup2@nsim, 1, OMup2@nyears))
## Data@AddIndV <- array(V, dim=c(OMup2@nsim, 1, OMup2@nyears))
## OMup2@cpars$Data <- Data
## OMup2@cpars$dummy = rep(0, OMup2@nsim)


## AR
##--------------------------------------------------------------------------------

ar <- c(
    "FMSYref2",
    "MSYHSref",
    ## HSX
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 3),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 2.5),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 2),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 1.8),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 1.6),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 1.4),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 1.2),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 1),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.9),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.8),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.7),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.6),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.5), ## HS - ICES default
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.4),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.3),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.2),
    get.MP(npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0.1),
    ## MSY rule
    get.MP(npriorSD = npriorSD,
           breakpointB = 0,
           nhist=nhist),
    ## msy rules fractiles on C (P* approach)
    get.MP(fractileC = 0.45,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.4,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.35,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.3,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.25,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.15,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.05,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.01,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0),
    get.MP(fractileC = 0.001,npriorSD = npriorSD,
           nhist=nhist, breakpointB = 0))

gc()

## print to console
print(species)
print(scenNum)
print(acNum)
print(sigmaRNum)
print(seedNum)

## Parallel processing
##--------------------------------------------------------------------------------
if(parallel){
    setup(cpus = nCores)
    sfLibrary(spict)
    sfExportAll()
}

## MSE <- runMSE(OMup, MPs = ar, Hist = TRUE)
## x=3
## opar <- par(mfrow=c(2,1))
## plot(MSE@TSdata$B[x,]/MSE@Ref$BMSY[x],ty='l',ylim=c(0,3))
## abline(h=1,lty=2)
## plot(MSE@TSdata$F[x,]/MSE@Ref$FMSY[x],ty='l',ylim=c(0,10))
## abline(h=1,lty=2)
##options(error=recover, warn=2)


## MSE
##--------------------------------------------------------------------------------
MSE <- runMSE(OMup, MPs = ar, timelimit = 200, CheckMPs = FALSE,
              parallel = parallel, PPD = TRUE)   ## PPD saves data for the last year


## Stop parallel processing
##--------------------------------------------------------------------------------
if(parallel) sfStop()

## Save output
##--------------------------------------------------------------------------------
save(MSE,file= paste0(savedir,"resMSE6_spec",specNum,"_scen",scenNum,"_nsim",nsim,
                      "_seed",seedNum,checkLab,".RData"))


if(FALSE){

    ## when debugging in HCR for creating spict example plots:
    pdf("spict_spec4_scen1.pdf", width=7, height=6)
    opar <- par(mfrow=c(2,2),mar=c(2,1,2,2),oma=c(3,3,2,0))
    plotspict.bbmsy(rep)
    plotspict.ffmsy(rep)
    plotspict.production(rep)
    plotspict.fb(rep)
    par(opar)
    dev.off()


    ## Fleets
    fleets <- avail("Fleet")
    class?Fleet

    opar <- par(mfrow=c(4,4),mar=c(2,1,2,2),oma=c(3,3,2,0))
    for(i in 1:length(fleets)){
        plotFleet(get(fleets[i]), xlab="", ylab="")
        mtext(fleets[i])
    }
    mtext("Effort", 2, 1.5, outer = TRUE)
    mtext("Time", 1, 1, outer = TRUE)
    par(opar)


    pdf("fig_effort.pdf", width=7, height=6)
    opar <- par(mar=c(5,5,4,2))
    plotFleet(OMup, xlab="", ylab="",yaxt="n", ylim=c(0,2), xaxt="n")
    axis(1, axis.cex=1.4,at = c(0,0.25,0.5,0.75,1), labels = c(0,10,20,30,40))
    axis(2, axis.cex=1.4)
    mtext("Effort", 2, 3,cex=1.3)
    mtext("Year", 1, 3, cex=1.3)
    par(opar)
    dev.off()

}
