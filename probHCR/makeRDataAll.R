## Stoch MSE paper
## April 2020
## Tobias Mildenberger <t.k.mildenberger@gmail.com>

server <- FALSE
checkLab <- ""
scenNum <- c(1,2,3,4,5)
seedvec <- c(1:10)


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
library(abind)


## set environment
## ---------------------------
## activate topic
## ---------------------------
addlabel <- "" ##"_altSel"   ## for file saving                          ## HERE:
checkLab <- ""   ## for file loading
if(server){
    resdir <- "../res/"
    source("stochSimFuncs.R")
}else{
    source("stochSimFuncs.R")
    resdir <- "res/2020-05-14/"                                          ## HERE:
}


speciesAll <- c("anchovy","haddock", "ling")##,c("anchovy","haddock","ling")
specNumAll <-  c(16,1,4) ##  c(18,19,20) ##                             ## HERE:
addLabAll <- c("")


## Scenario settings
dteuler <- 1/16
amtint <- 1
proyears <- 20
nsim <- 50     ## per seed
nhist <- 50 ## /30
pnpriorSD <- NA
nscen <- 8
nhcrs <- 19


## automatic
MSElistAll <- vector("list",length(speciesAll))
for(ii in 1:length(speciesAll)){

    species <- speciesAll[ii]
    specNum <- specNumAll[ii]
    addLab <- addLabAll[ii]

    ## Load data
    ##--------------------------------------------------------------------------------
    files <- dir(paste0(resdir))
    files0 <- strsplit(files, "_")

    ## all seeds
    seedi <- which(unlist(lapply(files0, function(x) strsplit(x[5],".RData"))) == "seed1")
    scenarios <- paste0("scen",scenNum)
    scenariosAll <- paste0(rep("scen",length(scenNum)),scenNum)
    ## combine
    expgri <- expand.grid(specNum,scenarios)
    combis <- apply(expand.grid(species,scenarios),1,paste0,collapse = "_")
    combisAll <- apply(expand.grid(species,scenariosAll),1,paste0,collapse = "_")

    ## build MSElists
    MSElist <- vector("list",length(combisAll))
    for(i in 1:length(combisAll)){
        MSElistS <- vector("list",length(seedvec))
        for(j in 1:length(seedvec)){
            seedNum <- seedvec[j]
            fileName <- paste0(resdir,"resMSE6_spec",expgri[i,1],"_",             ## resMSE4_spec
                               expgri[i,2],"_nsim",nsim,"_seed",seedNum,checkLab,".RData")
            writeLines(fileName)
            tmp <- try(load(fileName),silent=FALSE)
            fileNameAdd <- paste0(resdir,"resMSE6Add_spec",expgri[i,1],"_",             ## resMSE4_spec
                               expgri[i,2],"_nsim",nsim,"_seed",seedNum,checkLab,".RData")
            writeLines(fileNameAdd)
            tmpAdd <- try(load(fileNameAdd),silent=FALSE)
            writeLines(paste0(all(round(MSEAdd@OM$SSB0,3) == round(MSE@OM$SSB0,3))))
            MSE@SSB <- abind(MSE@SSB, MSEAdd@SSB, along = 2)
            MSE@B_BMSY <- abind(MSE@B_BMSY, MSEAdd@B_BMSY, along = 2)
            MSE@F_FMSY <- abind(MSE@F_FMSY, MSEAdd@F_FMSY, along = 2)
            MSE@C <- abind(MSE@C, MSEAdd@C, along = 2)
            MSE@TAC <- abind(MSE@TAC, MSEAdd@TAC, along = 2)
            MSE@Misc$Data <- c(MSE@Misc$Data,MSEAdd@Misc$Data) ## spict related quantities
            if(class(tmp) != "try-error") MSElistS[[j]] <- MSE
        }
        ## join MSE with different seed
        if(length(seedvec)>1){
            tmp <- try(joinMSE(MSElistS),silent=FALSE)
            resMSE <- list()
            if(class(tmp) != "try-error"){

                ## for Biomass
                resMSE$bbmsy <- tmp@B_BMSY
                resMSE$bmsyb0 <- tmp@OM$BMSY_B0 ## to check spict priors against

                ## for SSB (and recovery definition dependent on ssb)
                resMSE$ssb <- tmp@SSB
                resMSE$ssb0 <- tmp@OM$SSB0

                ## for F
                resMSE$ffmsy <- tmp@F_FMSY

                ## for Yield
                resMSE$catch <- tmp@C
                resMSE$refcatch <- tmp@OM$RefY

                ## for TAC
                resMSE$tac <- tmp@TAC

                ## for convergence
                resMSE$conv <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spict)))
                for(k in 1:dim(resMSE$conv)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$conv[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$conv[,k,] <- matrix(NA, ncol=ncol(resMSE$conv[,k,]), nrow=nrow(resMSE$conv[,k,]))
                    }
                }

                ## comparing TACs
                resMSE$spictTAC <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$TAC)))
                for(k in 1:dim(resMSE$spictTAC)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictTAC[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictTAC[,k,] <- matrix(NA, ncol=ncol(resMSE$spictTAC[,k,]), nrow=nrow(resMSE$spictTAC[,k,]))
                    }
                }

                ## spict BBmsy - est
                resMSE$spictBBMSY <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictBest)))
                for(k in 1:dim(resMSE$spictBBMSY)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictBBMSY[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictBBMSY[,k,] <- matrix(NA, ncol=ncol(resMSE$spictBBMSY[,k,]),
                                                         nrow=nrow(resMSE$spictBBMSY[,k,]))
                    }
                }
                ## spict FFmsy - est
                resMSE$spictFFMSY <- array(NA, c(dim(tmp@F_FMSY)[1:2],dim(tmp@F_FMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictFest)))
                for(k in 1:dim(resMSE$spictFFMSY)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictFFMSY[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictFFMSY[,k,] <- matrix(NA, ncol=ncol(resMSE$spictFFMSY[,k,]),
                                                         nrow=nrow(resMSE$spictFFMSY[,k,]))
                    }
                }
                ## spict BBmsy - sd
                resMSE$spictBsd <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictBsd)))
                for(k in 1:dim(resMSE$spictBsd)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictBsd[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictBsd[,k,] <- matrix(NA, ncol=ncol(resMSE$spictBsd[,k,]),
                                                       nrow=nrow(resMSE$spictBsd[,k,]))
                    }
                }
                ## spict FFmsy - sd
                resMSE$spictFsd <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictFsd)))
                for(k in 1:dim(resMSE$spictFsd)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictFsd[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictFsd[,k,] <- matrix(NA, ncol=ncol(resMSE$spictFsd[,k,]),
                                                       nrow=nrow(resMSE$spictFsd[,k,]))
                    }
                }
                ## spict Cp - est
                resMSE$spictCest <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictCest)))
                for(k in 1:dim(resMSE$spictCest)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictCest[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictCest[,k,] <- matrix(NA, ncol=ncol(resMSE$spictCest[,k,]),
                                                       nrow=nrow(resMSE$spictCest[,k,]))
                    }
                }
                ## spict Cp - sd
                resMSE$spictCsd <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictCsd)))
                for(k in 1:dim(resMSE$spictCsd)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictCsd[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictCsd[,k,] <- matrix(NA, ncol=ncol(resMSE$spictCsd[,k,]),
                                                       nrow=nrow(resMSE$spictCsd[,k,]))
                    }
                }
                ## spict MSY
                resMSE$spictMSY <- array(NA, c(dim(tmp@B_BMSY)[1:2],dim(tmp@B_BMSY)[3]-1))
                tmp2 <- lapply(tmp@Misc$Data, function(x) do.call(rbind, lapply(x@Misc, function(y) y$spictMSY)))
                for(k in 1:dim(resMSE$spictMSY)[2]){
                    if(!is.null(tmp2[[k]])){
                        resMSE$spictMSY[,k,] <- tmp2[[k]]
                    }else{
                        resMSE$spictMSY[,k,] <- matrix(NA, ncol=ncol(resMSE$spictMSY[,k,]),
                                                       nrow=nrow(resMSE$spictMSY[,k,]))
                    }
                }

                resMSE$BMSY <- tmp@OM$BMSY
                resMSE$FMSY <- tmp@OM$FMSY
                resMSE$MSY <- tmp@OM$MSY

            }
            MSElist[[i]] <- resMSE
        }else{
            MSElist[[i]] <- MSElistS[[1]]
        }
    }

    MSElistAll[[ii]] <- MSElist

    ## garbage collection
    rm(MSElistS, MSElist)
    gc()
}

## save
if(server){
    save(MSElistAll, file = paste0("MSElistAll_server.RData"))
}else{
    save(MSElistAll, file = paste0("MSElistAll_",Sys.Date(),".RData"))
}
