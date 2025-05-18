##  --------------------------------------------------------------------------
## Copyright (c) 2023, Tobias K. Mildenberger <tobm@aqua.dtu.dk>, Marc H. Taylor
## <marc.taylor@thuenen.de>
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in the
##     documentation and/or other materials provided with the distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL TOBIAS MILDENBERGER OR MARC TAYLOR BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
## --------------------------------------------------------------------------

## This scripts converts SAM output to spict input data.


## Directories
## -----------------------
datdir <- "../spmprod/data"
resdir <- "../output"
figdir <- "../res/figs/general"


dir.create(file.path("../res/figs","general"), showWarnings = FALSE)


## Load Packages
## -----------------------
## remotes::install_github("flr/FLCore")
library(FLCore)
## remotes::install_github('fishfollower/SAM/stockassessment')
## remotes::install_github('flr/FLSAM')
library(FLSAM)
## remotes::install_github('flr/ggplotFL')
library(ggplotFL)
## remotes::install_github('ices-tools-prod/icesSAG')
require(icesSAG)
library(mgcv)
library(reshape2)
library(data.table)
## install.packages("unikn")
require(unikn)
require(RColorBrewer)


source("funcs.R")


## Variables
## -----------------------
final.year <- 2021


## Load data
## -----------------------
stockInfo <- read.csv("../data/stockInfo.csv")
stockInfo <- stockInfo[which(stockInfo$file_name2 != ""),]
## To get ices keys for stocks:
## all.info <- icesSAG::getListStocks(year = 2021)
## all.stocks <- all.info$StockKeyLabel

stocks <- stockInfo$stock
print(stocks)
nstocks <- length(stocks)

## Loop through stocks
## -----------------------
stockData <- vector("list", nstocks)
names(stockData) <- stocks
for(stocki in seq(stockData)){

    print(stocks[stocki])

    ## load
    ## -----------------------
    objname <- load(file.path(datdir, "WGNSSK_2022","est_catch_versions",
                              stockInfo$file_name2[stocki]),
                    verbose = T)

    ## fixes
    ## -----------------------
    if(!(stocks[stocki] %in% c("PLE-EC", "TUR", "SOL-NS"))){
        stock <- get(objname)
    }
    ## if(stocks[stocki] == "HAD"){
    ##     stock@catch <- computeCatch(stock)
    ## }
    if(stocks[stocki] == "TUR"){
        stock <- TUR
        catch(stock) <- landings(stock)
        discards(stock) <- landings(stock)-catch(stock)
        ## add numbers
        ## catch.n(stock)
        ## landings.n(stock)
        discards.n(stock)[] <- 0
        ## catch.wt(stock)
        ## landings.wt(stock)
        discards.wt(stock)[is.na(discards.wt(stock))] <- 0
        ## add stock numbers
        stock.n(stock) <- TUR.sam@stock.n
        ## stock@stock.wt
        stock(stock) <- computeStock(stock)
        ## add harvest
        harvest(stock) <- TUR.sam@harvest
    }
    if(stocks[stocki] == "WIT"){
        ## Exclude time series before 1983 for witch flounder
        stock <- window(stock, 1983, final.year)
    }
    if(stocks[stocki] == "PLE-EC"){
        stock <- ass.stock
        harvest(stock)@units <- "f"
    }
    if(stocks[stocki] == "SOL-NS"){
        stock <- sol.27.4[["estimated"]]
    }

    ## r from stocks, but eggs per biomass unknown for most stocks
    print(mean(r(m(stock), mat(stock) * stock@stock.wt * 600)))

    ## Get info from SAG
    ## -----------------------
    sag.info0 <- icesSAG::getSAG(stock = stockInfo$ices_code[stocki],
                                 year = final.year)
    sag.info <- subset(sag.info0, Purpose == "Advice" & Year <= final.year)


    ## General
    ## -----------------------
    ageRange <- range(as.numeric(dimnames(stock@stock.n)$age))
    minAge <- ageRange[1]
    maxAge <- ageRange[2]
    ages <- minAge:maxAge
    nages <- length(ages)

    ## Cut stock
    minYear <- min(as.numeric(dimnames(stock@stock.n)$year))
    stock <- window(stock, minYear, final.year)

    yearRange <- range(as.numeric(dimnames(stock@stock.n)$year))
    minYear <- yearRange[1]
    maxYear <- yearRange[2]
    years <- minYear:maxYear


    ## Main data table
    dat <- data.frame(year = years)


    ## Catch in weight
    ## -----------------------
    C <- as.data.frame(catch(stock))

    ## dis <- as.data.frame(stock@discards.n[,44,1,1,1,1])$data
    ## land <- as.data.frame(stock@landings.n[,44,1,1,1,1])$data

    ## ns <- as.data.frame(stock@stock.n[,44,1,1,1,1])$data
    ## fs <- as.data.frame(stock@harvest[,44,1,1,1,1])$data
    ## ms <- as.data.frame(stock@m[,44,1,1,1,1])$data

    ## fs/(fs+ms) * ns * (1 - exp(-(fs+ms)))
    ## land + dis

    ## Fix for WIT
    if(stocks[stocki] == "WIT"){
        off.land <- tail(c(2379,2568,2554,2263,1590,2027,1936,1991,2678,2333,2563,
                           2093,1419,1523,1760,1356,1137,1125,1174,891,597,843,908,
                           1494,1138,1841,1496,1618,1664,1572,1883,1933,3155,3606,
                           3903,3979,3579,3700,3290,3841,3862,3641,3164,2673,2696,
                           2810,2790,3494,3786,4024,4422,4206,3640,3281,3029,2813,
                           2303,2236,1953),length(1983:2008))
        C$data[which(C$year %in% 1983:2008)] <- off.land
    }

    dat$C <- C$data[match(C$year,dat$year)]



    ## Fishing mortality
    ## -----------------------
    Fmean <- as.data.frame(apply(stock@harvest, 2, mean))
    ## For Fbar: scale F by percentage of total catch weight of each age:
    ## cwp <- apply(stock@catch.n * stock@catch.wt, 2, function(x) x/sum(x))
    ## Fbar <- as.data.frame(apply(stock@harvest * cwp, 2, mean))
    ## Better: use def of Fbar for each stock
    barAges <- strsplit(stockInfo$fbarDef[stocki],"-")[[1]]
    Fbar <- as.data.frame(apply(stock@harvest[barAges[1]:barAges[2],,,,,], 2, mean))

    dat$Fmean <- Fmean$data[match(Fmean$year,dat$year)]
    dat$Fbar <- Fbar$data[match(Fbar$year,dat$year)]



    ## TSB
    ## -----------------------
    TSB <- as.data.frame(apply(stock@stock.n * stock@stock.wt, 2, sum))
    dat$TSB <- TSB$data[match(TSB$year,dat$year)]



    ## SSB
    ## -----------------------
    SSB <- as.data.frame(apply(stock@stock.n *
                               exp(-stock@harvest*stock@harvest.spwn -
                                   stock@m*stock@m.spwn) *
                               stock@stock.wt * stock@mat, 2, sum))
    dat$SSB <- SSB$data[match(SSB$year,dat$year)]



    ## ESB
    ## -----------------------
    tmp <- as.data.frame(stock)
    mod.df <- reshape2::dcast(subset(tmp, slot %in% c("stock.wt", "catch.wt","catch.n",
                                                      "stock.n", "harvest", "m")),
                              formula = age + year ~ slot,
                              fun.aggregate = mean, value.var = "data")
    mod.df$biomass <- mod.df$stock.n * mod.df$stock.wt
    mod.df$age <- an(mod.df$age)
    mod.df$yearFac <- factor(mod.df$year, ordered = TRUE)
    mod.df <- as.data.table(mod.df)
    mod.df[, scHarvest := harvest/max(harvest, na.rm = T), by = list(year)]

    if(FALSE){

    if(stocks[stocki] %in% c("WIT","TUR")){
        fitSel <- gam(harvest ~ yearFac + s(age, k = nages-1) +
                          s(age, k = nages-1, by = yearFac),
                      data = mod.df, family = gaussian(link = "log"),
                      weights = stock.n)
    }else{
        fitSel <- gam(harvest ~ yearFac + s(age, k = nages-1) +
                          s(age, k = nages-1, by = yearFac),
                      data = mod.df, family = gaussian(link = "log"),
                      weights = catch.n)
    }
    mod.df$pred <- predict(object = fitSel, newdata = mod.df, type = "response")

    mod.df[, sel := pred/max(pred, na.rm = T), by = list(year)]


    ## Stock weights
    mod.df$ESB <- mod.df$stock.n * mod.df$stock.wt * mod.df$sel
    dat$ESB_wStock <- aggregate(ESB ~ year, mod.df, sum)$ESB

    ## Catch weights
    mod.df$sel2 <- mod.df$sel
    mod.df$ESB2 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
        mod.df$catch.wt * mod.df$sel2
    dat$ESB_wCatch <- aggregate(ESB2 ~ year, mod.df, sum)$ESB2

    ## No year-age interaction in GAM
    if(stocks[stocki] %in% c("WIT","TUR")){
        fitSel2 <- gam(harvest ~ yearFac + s(age, k = nages-1),
                       data = mod.df, family = gaussian(link = "log"),
                       weights = stock.n)
    }else{
        fitSel2 <- gam(harvest ~ yearFac + s(age, k = nages-1),
                       data = mod.df, family = gaussian(link = "log"),
                       weights = catch.n)
    }
    mod.df$pred <- predict(object = fitSel2, newdata = mod.df, type = "response")
    mod.df[, sel3 := pred/max(pred, na.rm = T), by = list(year)]


    ## Stock weights + no interaction
    mod.df$ESB3 <- mod.df$stock.n * mod.df$stock.wt * mod.df$sel3
    dat$ESB_wStock_av <- aggregate(ESB3 ~ year, mod.df, sum)$ESB3

    ## Catch weights + no interaction
    mod.df$sel4 <- mod.df$sel3
    mod.df$ESB4 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
        mod.df$catch.wt * mod.df$sel4
    dat$ESB_wCatch_av <- aggregate(ESB4 ~ year, mod.df, sum)$ESB4


    ## No year effect in GAM
    if(stocks[stocki] %in% c("WIT","TUR")){
        fitSel3 <- gam(harvest ~ s(age, k = nages-1),
                       data = mod.df, family = gaussian(link = "log"),
                       weights = stock.n)
    }else{
        fitSel3 <- gam(harvest ~ s(age, k = nages-1),
                       data = mod.df, family = gaussian(link = "log"),
                       weights = catch.n)
    }
    mod.df$pred <- predict(object = fitSel3, newdata = mod.df, type = "response")
    mod.df[, sel5 := pred/max(pred, na.rm = T), by = list(year)]

    ## Stock weights + no interaction
    mod.df$ESB5 <- mod.df$stock.n * mod.df$stock.wt * mod.df$sel5
    dat$ESB_wStock_ov <- aggregate(ESB5 ~ year, mod.df, sum)$ESB5

    ## Catch weights + no interaction
    mod.df$sel6 <- mod.df$sel5
    mod.df$ESB6 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
        mod.df$catch.wt * mod.df$sel6
        dat$ESB_wCatch_ov <- aggregate(ESB6 ~ year, mod.df, sum)$ESB6

    mods <- list(fitSel, fitSel2)

        ## Raw
        mod.df[, sel7 := harvest/max(harvest, na.rm = T), by = list(year)]
        mod.df$ESB7 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
            mod.df$catch.wt * mod.df$sel7
        dat$ESB_raw <- aggregate(ESB7 ~ year, mod.df, sum)$ESB7
        ## sel <- apply(stock@harvest, 2, function(x){x/max(x)})
        ## ESB <- as.data.frame(apply(stock@stock.n * stock@stock.wt * sel, 2, sum))
        ## dat$ESB_raw <- ESB$data[match(ESB$year,dat$year)]

        ## Raw - average
        mod.df[, sel8 := harvest/max(harvest, na.rm = T), by = list(year)]
        mod.df[, sel8 := sel8/mean(sel8, na.rm = T), by = list(age)]
        mod.df[, sel8 := sel8/max(sel8, na.rm = T), by = list(year)]
        mod.df$ESB8 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
            mod.df$catch.wt * mod.df$sel8
        dat$ESB_raw_av <- aggregate(ESB8 ~ year, mod.df, sum)$ESB8
        ## sel <- as.numeric(apply(apply(stock@harvest, 2, function(x){x/max(x)}), 1, mean))
        ## ESB <- as.data.frame(apply(stock@stock.n * stock@stock.wt * sel, 2, sum))
        ## dat$ESB_raw_av <- ESB$data[match(ESB$year,dat$year)]

        ## Raw - overall
        mod.df[, sel9 := harvest/max(harvest, na.rm = T), by = list(year)]
        mod.df[, sel9 := sel9/max(sel9, na.rm = T), by = list(age)]
        mod.df$ESB9 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
            mod.df$catch.wt * mod.df$sel9
        dat$ESB_raw_ov <- aggregate(ESB9 ~ year, mod.df, sum)$ESB9

        ## Raw - overall - stock weight
        mod.df[, sel10 := harvest/max(harvest, na.rm = T), by = list(year)]
        mod.df[, sel10 := sel10/max(sel10, na.rm = T), by = list(age)]
        mod.df$ESB10 <- mod.df$stock.n * mod.df$stock.wt * mod.df$sel10
        dat$ESB_raw_ov_sWeight <- aggregate(ESB10 ~ year, mod.df, sum)$ESB10

    }

    ## Raw new - USE!
    mod.df$sel0 <- NA
    years2 <- unique(mod.df$year)
    for(i in 1:length(years2)){
        tmp <- mod.df$harvest[mod.df$year == years2[i]]
        mod.df$sel0[mod.df$year == years2[i]] <- tmp / max(tmp)
    }
    mod.df$sel11 <- NA
    ages2 <- unique(mod.df$age)
    for(i in 1:length(ages2)){
        tmp <- mod.df$sel0[mod.df$age == ages2[i]]
        mod.df$sel11[mod.df$age == ages2[i]] <- mean(tmp)
    }
    mod.df$ESB11 <- mod.df$stock.n * (exp(-(mod.df$harvest + mod.df$m)/2)) *
        mod.df$catch.wt * mod.df$sel11
    dat$ESB_raw_ov_new <- aggregate(ESB11 ~ year, mod.df, sum)$ESB11

    mod.df$ESB12 <- mod.df$stock.n * mod.df$stock.wt * mod.df$sel11
    dat$ESB_raw_ov_sWeight_new <- aggregate(ESB12 ~ year, mod.df, sum)$ESB12

    mod.df$ESB13 <- mod.df$stock.n * mod.df$stock.wt * mod.df$sel0
    dat$ESB_raw_new <- aggregate(ESB13 ~ year, mod.df, sum)$ESB13

    ## plot(mod.df$age, mod.df$sel11)
    ## points(mod.df$age, mod.df$sel9)

    ## plot(dat$ESB_raw_ov,ty='b', ylim = c(0,6e5))
    ## lines(dat$ESB_raw_ov_new,ty='b',col=4)


    ## Alternative F (C/ESB)
    dat$Fesb <- dat$C / dat$ESB_raw_ov_sWeight_new ## ESB_wCatch


    ## Recruits
    ## -----------------------
    rec <- as.data.frame(rec(stock))
    ## NOTE: rec() not always recruits but N of min age class (for saithe: age 3)
    print(paste0("Minimum age in stock: ",min(as.data.frame(stock@stock.n)$age)))

    dat$rec <- rec$data[match(rec$year,dat$year)]
    dat$recSSB <- dat$rec/dat$SSB


    ## CV of SSB (from icesSAG)
    ## -----------------------
    est <- sag.info$SSB
    ul <- sag.info$high_SSB
    cv <- ((ul - est) / 1.96) / est
    cv <- cv[sag.info$Year %in% SSB$year]
    if(min(sag.info$Year) > min(SSB$year)){
        cv <- c(cv,rep(cv[1], min(sag.info$Year) -  min(SSB$year)))
    }

    dat$cv_ssb <- cv[match(SSB$year,dat$year)]


    ## Age composition
    ## -----------------------
    ageComp0 <- as.data.frame(stock@stock.n * stock@stock.wt)
    ageComp <- NULL
    for(i in ages) ageComp <- cbind(ageComp, ageComp0$data[ageComp0$age == i])
    ageComp <- as.data.frame(ageComp)
    rownames(ageComp) <- years
    colnames(ageComp) <- paste0("age", ages)


    ## Natural mortality
    ## -----------------------
    natMort0 <- as.data.frame(stock@m)
    natMort <- NULL
    for(i in ages) natMort <- cbind(natMort, natMort0$data[natMort0$age == i])
    natMort <- as.data.frame(natMort)
    rownames(natMort) <- years
    colnames(natMort) <- paste0("age", ages)


    ## Stock weight
    ## -----------------------
    wStock0 <- as.data.frame(stock@stock.wt)
    wStock <- NULL
    for(i in ages) wStock <- cbind(wStock, wStock0$data[wStock0$age == i])
    wStock <- as.data.frame(wStock)
    rownames(wStock) <- years
    colnames(wStock) <- paste0("age", ages)


    ## Catch weight
    ## -----------------------
    wCatch0 <- as.data.frame(stock@catch.wt)
    wCatch <- NULL
    for(i in ages) wCatch <- cbind(wCatch, wCatch0$data[wCatch0$age == i])
    wCatch <- as.data.frame(wCatch)
    rownames(wCatch) <- years
    colnames(wCatch) <- paste0("age", ages)


    ## Selectivity
    ## -----------------------
    sel0 <- as.data.frame(apply(stock@harvest, 2, function(x){x/max(x)}))
    sel <- NULL
    for(i in ages) sel <- cbind(sel, sel0$data[sel0$age == i])
    rownames(sel) <- years
    colnames(sel) <- paste0("age", ages)


    ## Maturity
    ## -----------------------
    mat0 <- as.data.frame(stock@mat)
    mat <- NULL
    for(i in ages) mat <- cbind(mat, mat0$data[mat0$age == i])
    rownames(mat) <- years
    colnames(mat) <- paste0("age", ages)


    ## Results
    ## -----------------------
    stockData[[stocki]] <- list(data = dat,
                                sel = sel, mat = mat,
                                natMort = natMort, wStock = wStock,
                                wCatch = wCatch,
                                ageComp = ageComp,
##                                mods = mods,
                                mod.df = mod.df)

    rm(stock, objname)

}



# Save
## -----------------------
save(stockData, file = file.path("../data", "stockData.RData"))

## load(file.path("../data", "stockData.RData"))



## ESB based refs and stock status in stockInfo (refs from make_refpts.R script)
## -----------------------
last.year <- as.data.frame(t(sapply(stockData,function(x) x$data[x$data$year == final.year,])))
stockInfo <- read.csv(file.path("../data", "stockInfo.csv"))
refs <- read.csv("../spmprod/output/BRPs2.csv")

stockInfo0 <- stockInfo

## refs
## stockInfo[,c("EQ_FmsyESB","EQ_ESBmsy","EQ_MSY")] <- refs[match(refs$X, stockInfo$stock),c("EQ_FmsyESB","EQ_ESBmsy","EQ_MSY")]
stockInfo[,c("EQ_FmsyESB","EQ_ESBmsy","EQ_MSY")] <- refs[match(refs$X, stockInfo$stock),c("EQ_Fmsy","EQ_Bmsy","EQ_MSY")]
## stock status
stockInfo$EQ_FFmsyESB <- as.numeric(last.year$Fesb[match(rownames(last.year), stockInfo$stock)])/as.numeric(stockInfo$EQ_FmsyESB)
stockInfo$EQ_ESBESBmsy <- as.numeric(last.year$ESB_raw_ov_new[match(rownames(last.year), stockInfo$stock)])/as.numeric(stockInfo$EQ_ESBmsy)
stockInfo$EQ_CMSY <- as.numeric(last.year$C[match(rownames(last.year), stockInfo$stock)])/as.numeric(stockInfo$EQ_MSY)

stockInfo$EQ_FmsyESB
stockInfo0$EQ_FmsyESB

## write.csv(stockInfo, file.path("../data", "stockInfo.csv"))
