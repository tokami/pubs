## Setting everything up - final development script for probHCR
## August 2020 - updated February 2021
## T.K. Mildenberger

## TODO: error messages when only one time step, or one rep!
## TODO: implement SSBfinal and use that for recruitment instead of FAA last
## season (selectivity can be very different between season 1 and 4)!!

## Script settings
## ------------------------
server <- TRUE
dbg <- 0
addLab <- ""


## Load packages
## ------------------------
## remotes::install_github("fishfollower/SAM/stockassessment")
## remotes::install_github("tokami/iamse", ref = "contAge")
require(iamse)
## remotes::install_github("tokami/spict/spict", ref="pubFF")
require(spict)


## Scenario settings
##--------------------------------------------------------------------------------
if(server){
    savedir <- "../res/"
    funcdir <- ""
    datdir <- ""
    input2 <- input <- as.character(commandArgs(TRUE))
    specNum <- as.integer(input2[1])
    scenNum <- as.integer(input2[2])
    seedNum <- as.integer(input2[3])
    hcrNum <- as.integer(input2[4])
    ncores <- as.numeric(Sys.getenv('LSB_DJOB_NUMPROC'))
    ncores <- ifelse(ncores > 5, ncores - 5, ncores)  ## for more memory
    nyhist <- 35
    nysim <- 35
    if(dbg > 0) nrep <- 250 else nrep <- 250
}else{
    savedir <- "res/"
    funcdir <- ""
    datdir <- ""
    specNum <- 2
    scenNum <- 1
    seedNum <- 1
    hcrNum <- 1
    ncores <- 1
    nyhist <- 35
    nysim <- 5
    nrep <- 1
}


## set seed
## ------------------------
seed <- seq(2,1e4,by=251)[seedNum]
set.seed(seed)


maxF <- 20 ## 100

## Data
## ------------------------
load(paste0(datdir,"stocklist6.RData"))
if(specNum == 1){
    spec <- stocklist[["anchovy"]]
    ns <- 2
    surveyBeforeAssessment <- c(TRUE, FALSE)
}else if(specNum == 2){
    spec <- stocklist[["haddock"]]
    ns <- 1
    surveyBeforeAssessment <- c(FALSE, FALSE)
}else if(specNum == 3){
    spec <- stocklist[["halibut"]]
    ns <- 1
    surveyBeforeAssessment <- c(FALSE, FALSE)
}else if(specNum == 4){
    spec <- stocklist[["anchovy2"]]
    ns <- 2
    surveyBeforeAssessment <- c(TRUE, FALSE)
}else if(specNum == 5){
    spec <- stocklist[["haddock2"]]
    ns <- 1
    surveyBeforeAssessment <- c(FALSE, FALSE)
}else if(specNum == 6){
    spec <- stocklist[["halibut2"]]
    ns <- 1
    surveyBeforeAssessment <- c(FALSE, FALSE)
}else if(specNum == 7){
    spec <- stocklist[["anchovy3"]]
    ns <- 2
    surveyBeforeAssessment <- c(TRUE, FALSE)
}else if(specNum == 8){
    spec <- stocklist[["haddock3"]]
    ns <- 1
    surveyBeforeAssessment <- c(FALSE, FALSE)
}else if(specNum == 9){
    spec <- stocklist[["halibut3"]]
    ns <- 1
    surveyBeforeAssessment <- c(FALSE, FALSE)
}


## else if(specNum == 7){
##     spec <- stocklist[["anchovyX1"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 8){
##     spec <- stocklist[["anchovyX2"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 9){
##     spec <- stocklist[["anchovyX3"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 10){
##     spec <- stocklist[["anchovyY1"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 11){
##     spec <- stocklist[["anchovyY2"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 12){
##     spec <- stocklist[["anchovyY3"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 13){
##     spec <- stocklist[["anchovyZ1"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 14){
##     spec <- stocklist[["anchovyZ2"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 15){
##     spec <- stocklist[["anchovyZ3"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 16){
##     spec <- stocklist[["anchovyH1"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 17){
##     spec <- stocklist[["anchovyH2"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }else if(specNum == 18){
##     spec <- stocklist[["anchovyH3"]]
##     ns <- 2
##     maxF <- 100
##     surveyBeforeAssessment <- c(TRUE, FALSE)
## }



spec$nseasons <- ns
spec$surveyBeforeAssessment <- surveyBeforeAssessment
dat <- check.dat(spec)




## Baseline settings
## ------------------------
set <- check.set()
set$maxF <- maxF
set$seed <- seed
set$nyhist <- nyhist
set$nysim <- nysim
set$nrep <- nrep
set$noiseR <- c(dat$sigmaR, dat$rhoR, 1)
set$noiseF <- c(0.15, 0, 1)
set$assessmentInterval <- 1
## Carruthers 2016: sdC: 0.2-0.4 and 0.3-0.6 and sdI: 0.1-0.3 and 0.2-0.6
## Wiedemann 2017: sdC: 0.15 and sdI: 0.29 and 0.63 (but 80 years of data)
set$noiseC <- c(0.3,0,0)
set$noiseI <- c(0.3,0,0)
set$noiseE <- c(0.3,0,0)
schaefer <- FALSE
## SPiCT prior: c(log(2),2,1)
## priorlogn <- c(log(2),2,1)
## Thorson's prior: c(log(1.478),sqrt(0.849^2/1.478^2),1)
priorlogn <- c(log(1.478),sqrt(0.849^2/1.478^2),1)
## priorlogn <- c(log(1.478),2,1)  ## wide Thorson's prior
priorlogsdf <- c(3,1,0)
priorlogsdc <- c(log(0.2),2,0)
priorlogalpha <- c(log(1),2,1)
priorlogbeta <- c(log(1),2,1)
priorlogbkfrac <- c(log(0.8),2,0)
dteuler <- 1/4
stab <- TRUE
clType <- "observed"
clyears <- 1 ## 3
clfac <- 100
lower <- 1/clfac
upper <- clfac



## Scenarios (1-X)
## ------------------------
## Baseline 1 => no changes
if(scenNum == 1){
    hcrID <- "base"
}else if(scenNum %in% 2:5){
    hcrID <- "noise"
}else if(scenNum %in% 6:13){
    hcrID <- "sensitivity"
}else if(scenNum %in% 14:30){
    hcrID <- "pseudo"
    ## assume assessment every season for pseudo analysis
    if(specNum %in% c(1,4,7)){  ## HERE: added 7
        set$assessmentTiming <- c(1,2)
        set$assessmentInterval <- 1
    }
}else if(scenNum %in% 31:40){
    hcrID <- "noAuto"
}


## Observation noise levels
if(scenNum == 2){
    set$noiseC <- 0.5 * c(0.3,0,0)
    set$noiseI <- 0.5 * c(0.3,0,0)
}else if(scenNum == 3){
    set$noiseC <- 2 * c(0.3,0,0)
    set$noiseI <- 2 * c(0.3,0,0)
}



## Recruitment noise levels
if(scenNum == 4){
    set$noiseR <- c(0.5 * dat$sigmaR, dat$rhoR, 1)
}else if(scenNum == 5){
    set$noiseR <- c(1.5 * dat$sigmaR, dat$rhoR, 1)
}



## Shorter time series
if(scenNum == 6){
    set$nyhist <- 20
    dat$nyC <- 20
    dat$nyI <- 20
    dat$selI[[2]] <- NULL
    dat$surveyTimes <- 7/12
}


## Sensitivity runs
## ------------------------------
## dteuler
if(scenNum == 7){
    dteuler <- 1/8
}

## intermediate year + constant F
## intermediate year + constant catch
manstartdY <- ifelse(scenNum %in% c(8,9), 1, 0)
intC <- ifelse(scenNum == 9, 1, NA)


## additional data-limited scenario: no survey, but effort
if(scenNum == 10){
    set$nyhist <- 35
    dat$nyC <- 35
    dat$nyI <- NA
    dat$nyE <- 35
    dat$selI[[1]] <- NULL
    dat$selI[[2]] <- NULL
    dat$surveyTimes <- NA
    priorlogalpha <- c(log(1),2,0)
}


## implementation uncertainty
## 0.1 assumed by Walsh (2018)
## 0.2 assumed by Nieland (2008)
## others: 0.1 assumed by Fisher (2020) AND 0.18 (?) by Dichmont (2006)
if(scenNum == 11){
    set$noiseImp <- c(0.15,0,1)
}

## default priors
if(scenNum == 12){
    priorlogn <- c(log(2),2,1)
    priorlogalpha <- c(log(1),2,1)
    priorlogbeta <- c(log(1),2,1)
}
## no priors
if(scenNum == 13){
    priorlogn <- c(0,0,0)
    priorlogalpha <- c(log(1),2,0)
    priorlogbeta <- c(log(1),2,0)
}


## Pseudo-assessment scenarios (tacSD based on tacs$cp.sd: 0.0-0.6)
if(scenNum == 14){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0
    set$tacSD <- 0
}
if(scenNum == 15){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0
    set$tacSD <- 0.3
}
if(scenNum == 16){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0
    set$tacSD <- 0.6
}
## biases (based on MB: 0.01-4)
if(scenNum == 17){  ## pos. bias in Fmsy implies neg. bias in F/Fmsy
    set$bbmsySD <- 0
    set$bbmsyBias <- 0
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0.5
    set$tacSD <- 0.3
}
if(scenNum == 18){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- -0.5
    set$tacSD <- 0.3
}
if(scenNum == 19){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0.5
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0
    set$tacSD <- 0.3
}
if(scenNum == 20){
    set$bbmsySD <- 0
    set$bbmsyBias <- -0.5
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0
    set$tacSD <- 0.3
}
if(scenNum == 21){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0.5
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0.5
    set$tacSD <- 0.3
}
if(scenNum == 22){
    set$bbmsySD <- 0
    set$bbmsyBias <- -0.5
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- -0.5
    set$tacSD <- 0.3
}
if(scenNum == 23){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0.5
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- -0.5
    set$tacSD <- 0.3
}
if(scenNum == 24){
    set$bbmsySD <- 0
    set$bbmsyBias <- -0.5
    set$ffmsySD <- 0
    set$ffmsyBias <- 0.5
    set$fmsyBias <- 0.5
    set$tacSD <- 0.3
}
if(scenNum == 25){
    set$bbmsySD <- 0
    set$bbmsyBias <- 0
    set$ffmsySD <- 0
    set$ffmsyBias <- 0
    set$fmsyBias <- 0
    set$tacSD <- 0.9
}

## HCRs
## ------------------------
source(paste0(funcdir,"loadHCRs.R"))  ## uses hcrNum (or dbg) to set HCRs
set$hcr <- hcrsAll
print(set$hcr)


## MSE (reference levels validation)
## ------------------------
t1 <- Sys.time()
resMSE <- run.mse(dat, set, ncores = ncores, verbose = FALSE)
t2 <- Sys.time()
print(t2 - t1)


## RData
## ------------------------
filepath <- paste0(savedir,
                   "resMSE_spec",specNum,"_scen",scenNum,"_seed",
                   seedNum,"_hcr",hcrNum, addLab,".RData")
save(resMSE, file = filepath)




## Print some checks
## ------------------------
print("Convergence rate - replicates")
print(lapply(resMSE, function(x) mean(sapply(x, function(x) all(x$tacs$conv)))))

print("Convergence rate - assessments")
print(lapply(resMSE, function(x) mean(sapply(x, function(x) x$tacs$conv))))

print("Average prob. B[last] < Blim")
print(lapply(resMSE, function(x) mean(sapply(x, function(x) x$TSBfinal[35]) < dat$ref$Blim[35])))




## OLDER
if(FALSE){

    saveFig <- FALSE
    ## Plot
    ## ------------------------
    ## quant over time
    filepath <- paste0(savedir, "tsplot_spec",specNum,"_scen",scenNum,"_seed",seedNum,"_hcr",hcrNum,
                       addLab,".pdf")
    if(saveFig) pdf(filepath, width = 14, height = 12)
    hcrs <- set$hcr
    opar <- par(mfrow=c(3,1),mar=c(3,4,2,2),oma=c(2,1,1,0))
    plotmse.b(dat, set, resMSE, hcrs=hcrs)##, ylim=c(0,3e10))
    plotmse.f(dat, set, resMSE, hcrs=hcrs)##, ylim=c(0,3))
    plotmse.cw(dat, set, resMSE, hcrs=hcrs)##, ylim=c(0,4e9))
    par(opar)
    if(saveFig) dev.off()

    mets <- estMets(resMSE, dat, mets=c("PBBlim","CMSY"))
    mets[[2]]
}
