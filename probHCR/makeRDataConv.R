## Stoch MSE paper
## May 2020
## Tobias Mildenberger <t.k.mildenberger@gmail.com>


## Settings
##--------------------------------------------------------------------------------
## input data
date_MSElistAll <- "2020-05-17"
## only consider evalyears for removing non converged reps and for total stats
evalyears <- 3:20
server <- FALSE



## load packages
##--------------------------------------------------------------------------------
##  devtools::install_github("tokami/spict/spict", ref="wklife9")
## install.packages("DLMtool")
library(spict)
library(DLMtool)
library(ggplot2)
library(lattice)
library(scales)
library(xtable)
library(plotrix)
library(reshape2)  # to draw pattern in boxplot


## set environment + load RData (from makeRData.R)
## ---------------------------
if(server){
    resdir <- "../res/"
    source("funcs3.R")
    load(paste0("MSElistAll_server.RData"))
}else{
    source("funcs3.R")
    load(paste0("MSElistAll_",date_MSElistAll,".RData"))
}

## ---------------------------
## list of 3 (species), list of 8 (scenarios), list of 15 (quantities)
nspec <- length(MSElistAll)
nscen <- length(MSElistAll[[1]])


## quantities (might change):
## --------------------------
## - bbmsy [rep, hcr, years]
## - bmsyb0 [rep]
## - ffmsy [rep, hcr, years]
## - catch [rep, hcr, years]
## - refcatch [rep]
## - tac [rep, hcr, years]
## - conv [rep, hcr, years-1]
## - spictTAC [rep, hcr, years-1]
## - spictBBMSY [rep, hcr, years-1]
## - spictFFMSY [rep, hcr, years-1]
## - spictBMSY [rep, hcr, years-1]
## - spictFMSY [rep, hcr, years-1]
## - BMSY [rep]
## - FMSY [rep]
## - MSY [rep]



## Remove non-converged replicates
## ---------------------------
## all HCRs
hcrsAll = c("Ref","Ref",
            ## HSx
            "HS300","HS250","HS200","HS180","HS160","HS140","HS120",
            "HS100","HS90","HS80","HS70","HS60","HS50","HS40","HS30","HS20","HS10",
            ## MSY-Cx
            "MSY-Med",
            "MSY-C45","MSY-C40","MSY-C35","MSY-C30","MSY-C25","MSY-C15","MSY-C05",
            "MSY-C01","MSY-C001",
            ## Ref again
            "Ref","Ref",
            ## HS-Cx
            "HS-C45","HS-C40","HS-C35","HS-C30","HS-C25","HS-C15","HS-C05","HS-C01","HS-C001",
            ## MSY-Ax
            "MSY-A45","MSY-A40","MSY-A35","MSY-A30","MSY-A25","MSY-A15","MSY-A05",
            ## HS-Ax
            "HS-A45","HS-A40","HS-A35","HS-A30","HS-A25","HS-A15","HS-A05")
## new alt HS100-x:
## c("HS100-Med","HS100-A45","HS100-A40","HS100-A35","HS100-A30","HS100-A25",
## "HS100-A15","HS100-A05",
## "HS100-C45","HS100-C40","HS100-C35","HS100-C30","HS100-C25",
## "HS100-C15","HS100-C05","HS100-C01","HS100-C001")


## use all spict HCRs + Refs
keepExtra <- 1:2
indexConv <- c(3:29,32:54)


## get MSElist with only converged reps
MSElistConv <- list()
for(spec in 1:nspec){
    ## container
    MSElistConv[[spec]] <- list()
    for(scen in 1:nscen){
        MSElistConv[[spec]][[scen]] <- get.conv(MSElistAll[[spec]][[scen]],
                                                indexConv = indexConv,
                                                convInfo = NULL,
                                                keepExtra = keepExtra,
                                                all = FALSE, evalyears = evalyears)
    }
}



## number of converged reps
print(lapply(MSElistConv, function(x) lapply(x, function(y) length(y$refcatch))))


## Derived quantities
## ---------------------------
for(spec in 1:nspec){
    for(scen in 1:nscen){

        ## BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyYear <- get.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyTot <- get.bbmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob B in interval 10% BBMSY
        ## MSElistConv[[spec]][[scen]]$bbmsy10 <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$bbmsy10Tot <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyYear <- get.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyTot <- get.ffmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob F in interval 10% FFMSY
        ## MSElistConv[[spec]][[scen]]$ffmsy10 <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$ffmsy10Tot <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## Risk (Prop of B < Blim)
        MSElistConv[[spec]][[scen]]$pblim <- get.pblim(MSElistConv[[spec]][[scen]], ref = 0.3)
        MSElistConv[[spec]][[scen]]$pblimTot <- get.pblim(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          ref = 0.3, evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$pblimY6_20 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                               ref = 0.3,
                                                               years = FALSE,
                                                               evalyears = 6:20)
        MSElistConv[[spec]][[scen]]$pblimY3_5 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                              ref = 0.3,
                                                              years = FALSE,
                                                              evalyears = 3:5)

        ## Rel. yield
        MSElistConv[[spec]][[scen]]$yield <- get.yield(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$yieldTot <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$yieldY3_5 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                           evalyears = 3:5)
        MSElistConv[[spec]][[scen]]$yieldY6_20 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                              evalyears = 6:20)

        ## AAV
        MSElistConv[[spec]][[scen]]$yieldDiff <- get.aav(MSElistConv[[spec]][[scen]],
                                                                evalyears = evalyears)

        ## spict est ref levels - BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyRE <- get.re.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyRETot <- get.re.bbmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$bbmsySD <- get.bbmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsySDTot <- get.bbmsySD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## spict est ref levels - FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyRE <- get.re.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyRETot <- get.re.ffmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$ffmsySD <- get.ffmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsySDTot <- get.ffmsySD(MSElistConv[[spec]][[scen]],
                                                              years = FALSE)

        ## "rel err in TAC" (= diff between TAC of HCR and Ref-HCR)
        MSElistConv[[spec]][[scen]]$tacRE <- get.re.tac(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$tacRETot <- get.re.tac(MSElistConv[[spec]][[scen]],
                                                         evalyears = evalyears, years = FALSE)
        ## SD of spict Cp
        MSElistConv[[spec]][[scen]]$cpSD <- get.cpSD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$cpSDTot <- get.cpSD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## Sample size
        MSElistConv[[spec]][[scen]]$sampleSize <- length(MSElistConv[[spec]][[scen]]$refcatch)

        ## Time to recovery
        MSElistConv[[spec]][[scen]]$timeRecTot <- get.timeRec(MSElistConv[[spec]][[scen]], ref = 0.3,
                                                              maxyear = 20)
    }
}


## save file
save(MSElistConv, file = paste0("MSElistConv_Yrs",
                                paste(range(evalyears),collapse = "-"),
                                "_",Sys.Date(),".RData"))
rm(MSElistConv)
gc()




## converged runs between 1 - 3
## ----------------------------------------------------------------------------------------------------
scenariosOverlapp <- 1:3

convInfo <- list()
for(spec in 1:nspec){
    convInfo[[spec]] <- list()
    for(scen in 1:nscen){
        convInfo[[spec]][[scen]] <- get.convInd(MSElistAll[[spec]][[scen]]$conv,
                                                index = indexConv,
                                                evalyears = evalyears)
    }
}
## index of overlapping replicates between scenarios
for(spec in 1:nspec){
    tmp <- convInfo[[spec]][[scenariosOverlapp[1]]]
    for(i in scenariosOverlapp[-1]){
        tmp <- intersect(tmp, convInfo[[spec]][[i]])
    }
    for(i in scenariosOverlapp){
        convInfo[[spec]][[i]] <- tmp
    }
}
print(lapply(convInfo, function(x) lapply(x,length)))


MSElistConv <- list()
for(spec in 1:nspec){
    ## container
    MSElistConv[[spec]] <- list()
    for(scen in 1:nscen){
        convTmp <- convInfo[[spec]][[scen]]
        MSElistConv[[spec]][[scen]] <- get.conv(MSElistAll[[spec]][[scen]],
                                                indexConv = indexConv,
                                                convInfo = convTmp,
                                                keepExtra = keepExtra,
                                                all = FALSE, evalyears = evalyears)
    }
}



## Derived quantities
## ---------------------------
for(spec in 1:nspec){
    for(scen in 1:nscen){

        ## BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyYear <- get.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyTot <- get.bbmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob B in interval 10% BBMSY
        ## MSElistConv[[spec]][[scen]]$bbmsy10 <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$bbmsy10Tot <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyYear <- get.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyTot <- get.ffmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob F in interval 10% FFMSY
        ## MSElistConv[[spec]][[scen]]$ffmsy10 <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$ffmsy10Tot <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## Risk (Prop of B < Blim)
        MSElistConv[[spec]][[scen]]$pblim <- get.pblim(MSElistConv[[spec]][[scen]], ref = 0.3)
        MSElistConv[[spec]][[scen]]$pblimTot <- get.pblim(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          ref = 0.3, evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$pblimY6_20 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                               ref = 0.3,
                                                               years = FALSE,
                                                               evalyears = 6:20)
        MSElistConv[[spec]][[scen]]$pblimY3_5 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                              ref = 0.3,
                                                              years = FALSE,
                                                              evalyears = 3:5)

        ## Rel. yield
        MSElistConv[[spec]][[scen]]$yield <- get.yield(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$yieldTot <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$yieldY3_5 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                           evalyears = 3:5)
        MSElistConv[[spec]][[scen]]$yieldY6_20 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                              evalyears = 6:20)

        ## AAV
        MSElistConv[[spec]][[scen]]$yieldDiff <- get.aav(MSElistConv[[spec]][[scen]],
                                                                evalyears = evalyears)

        ## spict est ref levels - BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyRE <- get.re.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyRETot <- get.re.bbmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$bbmsySD <- get.bbmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsySDTot <- get.bbmsySD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## spict est ref levels - FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyRE <- get.re.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyRETot <- get.re.ffmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$ffmsySD <- get.ffmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsySDTot <- get.ffmsySD(MSElistConv[[spec]][[scen]],
                                                              years = FALSE)

        ## "rel err in TAC" (= diff between TAC of HCR and Ref-HCR)
        MSElistConv[[spec]][[scen]]$tacRE <- get.re.tac(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$tacRETot <- get.re.tac(MSElistConv[[spec]][[scen]],
                                                         evalyears = evalyears, years = FALSE)
        ## SD of spict Cp
        MSElistConv[[spec]][[scen]]$cpSD <- get.cpSD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$cpSDTot <- get.cpSD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## Sample size
        MSElistConv[[spec]][[scen]]$sampleSize <- length(MSElistConv[[spec]][[scen]]$refcatch)

        ## Time to recovery
        MSElistConv[[spec]][[scen]]$timeRecTot <- get.timeRec(MSElistConv[[spec]][[scen]], ref = 0.3,
                                                              maxyear = 20)
    }
}

MSElistConv1_3 <- MSElistConv


## save file
save(MSElistConv1_3, file = paste0("MSElistConv1_3_Yrs",
                                paste(range(evalyears),collapse = "-"),
                                "_",Sys.Date(),".RData"))
rm(MSElistConv,MSElistConv1_3)
gc()




## converged runs between 2 & 4
## ----------------------------------------------------------------------------------------------------
scenariosOverlapp <- c(2,4)

convInfo <- list()
for(spec in 1:nspec){
    convInfo[[spec]] <- list()
    for(scen in 1:nscen){
        convInfo[[spec]][[scen]] <- get.convInd(MSElistAll[[spec]][[scen]]$conv,
                                                index = indexConv,
                                                evalyears = evalyears)
    }
}
## index of overlapping replicates between scenarios
for(spec in 1:nspec){
    tmp <- convInfo[[spec]][[scenariosOverlapp[1]]]
    for(i in scenariosOverlapp[-1]){
        tmp <- intersect(tmp, convInfo[[spec]][[i]])
    }
    for(i in scenariosOverlapp){
        convInfo[[spec]][[i]] <- tmp
    }
}
print(lapply(convInfo, function(x) lapply(x,length)))


MSElistConv <- list()
for(spec in 1:nspec){
    ## container
    MSElistConv[[spec]] <- list()
    for(scen in 1:nscen){
        convTmp <- convInfo[[spec]][[scen]]
        MSElistConv[[spec]][[scen]] <- get.conv(MSElistAll[[spec]][[scen]],
                                                indexConv = indexConv,
                                                convInfo = convTmp,
                                                keepExtra = keepExtra,
                                                all = FALSE, evalyears = evalyears)
    }
}



## Derived quantities
## ---------------------------
for(spec in 1:nspec){
    for(scen in 1:nscen){

        ## BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyYear <- get.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyTot <- get.bbmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob B in interval 10% BBMSY
        ## MSElistConv[[spec]][[scen]]$bbmsy10 <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$bbmsy10Tot <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyYear <- get.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyTot <- get.ffmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob F in interval 10% FFMSY
        ## MSElistConv[[spec]][[scen]]$ffmsy10 <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$ffmsy10Tot <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## Risk (Prop of B < Blim)
        MSElistConv[[spec]][[scen]]$pblim <- get.pblim(MSElistConv[[spec]][[scen]], ref = 0.3)
        MSElistConv[[spec]][[scen]]$pblimTot <- get.pblim(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          ref = 0.3, evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$pblimY6_20 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                               ref = 0.3,
                                                               years = FALSE,
                                                               evalyears = 6:20)
        MSElistConv[[spec]][[scen]]$pblimY3_5 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                              ref = 0.3,
                                                              years = FALSE,
                                                              evalyears = 3:5)

        ## Rel. yield
        MSElistConv[[spec]][[scen]]$yield <- get.yield(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$yieldTot <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$yieldY3_5 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                           evalyears = 3:5)
        MSElistConv[[spec]][[scen]]$yieldY6_20 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                              evalyears = 6:20)

        ## AAV
        MSElistConv[[spec]][[scen]]$yieldDiff <- get.aav(MSElistConv[[spec]][[scen]],
                                                                evalyears = evalyears)

        ## spict est ref levels - BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyRE <- get.re.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyRETot <- get.re.bbmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$bbmsySD <- get.bbmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsySDTot <- get.bbmsySD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## spict est ref levels - FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyRE <- get.re.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyRETot <- get.re.ffmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$ffmsySD <- get.ffmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsySDTot <- get.ffmsySD(MSElistConv[[spec]][[scen]],
                                                              years = FALSE)

        ## "rel err in TAC" (= diff between TAC of HCR and Ref-HCR)
        MSElistConv[[spec]][[scen]]$tacRE <- get.re.tac(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$tacRETot <- get.re.tac(MSElistConv[[spec]][[scen]],
                                                         evalyears = evalyears, years = FALSE)
        ## SD of spict Cp
        MSElistConv[[spec]][[scen]]$cpSD <- get.cpSD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$cpSDTot <- get.cpSD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## Sample size
        MSElistConv[[spec]][[scen]]$sampleSize <- length(MSElistConv[[spec]][[scen]]$refcatch)

        ## Time to recovery
        MSElistConv[[spec]][[scen]]$timeRecTot <- get.timeRec(MSElistConv[[spec]][[scen]], ref = 0.3,
                                                              maxyear = 20)
    }
}

MSElistConv2_4 <- MSElistConv


## save file
save(MSElistConv2_4, file = paste0("MSElistConv2_4_Yrs",
                                paste(range(evalyears),collapse = "-"),
                                "_",Sys.Date(),".RData"))
rm(MSElistConv,MSElistConv2_4)
gc()


## converged runs between 2 & 5
## ----------------------------------------------------------------------------------------------------
scenariosOverlapp <- c(2,5)

convInfo <- list()
for(spec in 1:nspec){
    convInfo[[spec]] <- list()
    for(scen in 1:nscen){
        convInfo[[spec]][[scen]] <- get.convInd(MSElistAll[[spec]][[scen]]$conv,
                                                index = indexConv,
                                                evalyears = evalyears)
    }
}
## index of overlapping replicates between scenarios
for(spec in 1:nspec){
    tmp <- convInfo[[spec]][[scenariosOverlapp[1]]]
    for(i in scenariosOverlapp[-1]){
        tmp <- intersect(tmp, convInfo[[spec]][[i]])
    }
    for(i in scenariosOverlapp){
        convInfo[[spec]][[i]] <- tmp
    }
}
print(lapply(convInfo, function(x) lapply(x,length)))


MSElistConv <- list()
for(spec in 1:nspec){
    ## container
    MSElistConv[[spec]] <- list()
    for(scen in 1:nscen){
        convTmp <- convInfo[[spec]][[scen]]
        MSElistConv[[spec]][[scen]] <- get.conv(MSElistAll[[spec]][[scen]],
                                                indexConv = indexConv,
                                                convInfo = convTmp,
                                                keepExtra = keepExtra,
                                                all = FALSE, evalyears = evalyears)
    }
}



## Derived quantities
## ---------------------------
for(spec in 1:nspec){
    for(scen in 1:nscen){

        ## BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyYear <- get.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyTot <- get.bbmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob B in interval 10% BBMSY
        ## MSElistConv[[spec]][[scen]]$bbmsy10 <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$bbmsy10Tot <- get.bbmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyYear <- get.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyTot <- get.ffmsy(MSElistConv[[spec]][[scen]],
                                                          years = FALSE,
                                                          evalyears = evalyears)

        ## Prob F in interval 10% FFMSY
        ## MSElistConv[[spec]][[scen]]$ffmsy10 <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                  cv = 0.1, interval = TRUE)
        ## MSElistConv[[spec]][[scen]]$ffmsy10Tot <- get.ffmsy(MSElistConv[[spec]][[scen]],
        ##                                                     years = FALSE, interval = TRUE,
        ##                                                     cv = 0.1, evalyears = evalyears)

        ## Risk (Prop of B < Blim)
        MSElistConv[[spec]][[scen]]$pblim <- get.pblim(MSElistConv[[spec]][[scen]], ref = 0.3)
        MSElistConv[[spec]][[scen]]$pblimTot <- get.pblim(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          ref = 0.3, evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$pblimY6_20 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                               ref = 0.3,
                                                               years = FALSE,
                                                               evalyears = 6:20)
        MSElistConv[[spec]][[scen]]$pblimY3_5 <- get.pblim(MSElistConv[[spec]][[scen]],
                                                              ref = 0.3,
                                                              years = FALSE,
                                                              evalyears = 3:5)

        ## Rel. yield
        MSElistConv[[spec]][[scen]]$yield <- get.yield(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$yieldTot <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                          evalyears = evalyears)
        MSElistConv[[spec]][[scen]]$yieldY3_5 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                           evalyears = 3:5)
        MSElistConv[[spec]][[scen]]$yieldY6_20 <- get.yield(MSElistConv[[spec]][[scen]], years = FALSE,
                                                              evalyears = 6:20)

        ## AAV
        MSElistConv[[spec]][[scen]]$yieldDiff <- get.aav(MSElistConv[[spec]][[scen]],
                                                                evalyears = evalyears)

        ## spict est ref levels - BBMSY
        MSElistConv[[spec]][[scen]]$bbmsyRE <- get.re.bbmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsyRETot <- get.re.bbmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$bbmsySD <- get.bbmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$bbmsySDTot <- get.bbmsySD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## spict est ref levels - FFMSY
        MSElistConv[[spec]][[scen]]$ffmsyRE <- get.re.ffmsy(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsyRETot <- get.re.ffmsy(MSElistConv[[spec]][[scen]],
                                                               evalyears = evalyears, years = FALSE)
        MSElistConv[[spec]][[scen]]$ffmsySD <- get.ffmsySD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$ffmsySDTot <- get.ffmsySD(MSElistConv[[spec]][[scen]],
                                                              years = FALSE)

        ## "rel err in TAC" (= diff between TAC of HCR and Ref-HCR)
        MSElistConv[[spec]][[scen]]$tacRE <- get.re.tac(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$tacRETot <- get.re.tac(MSElistConv[[spec]][[scen]],
                                                         evalyears = evalyears, years = FALSE)
        ## SD of spict Cp
        MSElistConv[[spec]][[scen]]$cpSD <- get.cpSD(MSElistConv[[spec]][[scen]])
        MSElistConv[[spec]][[scen]]$cpSDTot <- get.cpSD(MSElistConv[[spec]][[scen]], years = FALSE)

        ## Sample size
        MSElistConv[[spec]][[scen]]$sampleSize <- length(MSElistConv[[spec]][[scen]]$refcatch)

        ## Time to recovery
        MSElistConv[[spec]][[scen]]$timeRecTot <- get.timeRec(MSElistConv[[spec]][[scen]], ref = 0.3,
                                                              maxyear = 20)
    }
}

MSElistConv2_5 <- MSElistConv


## save file
save(MSElistConv2_5, file = paste0("MSElistConv2_5_Yrs",
                                   paste(range(evalyears),collapse = "-"),
                                   "_",Sys.Date(),".RData"))
rm(MSElistConv,MSElistConv2_5)
gc()
