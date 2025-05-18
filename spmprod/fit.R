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

## This scripts fits spict to demersal fish stocks.

## Packages
## -----------------
require(R.devices)
require(corrplot)
## remotes::install_github("tokami/spict/spict", ref="mChangeK2")
require(spict)
## install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.0.5.tar.gz",
##                  repo=NULL, type="source")
## install.packages("rfishbase")
require(rfishbase)
## remotes::install_github("kaskr/TMB_contrib_R/TMBhelper")
## remotes::install_github("james-thorson/FishLife")
## remotes::install_github("James-Thorson-NOAA/FishLife") v.3.1.0
require(FishLife)
## remotes::install_github("henning-winker/SPMpriors")
require(SPMpriors)

sessionInfo()
## other attached packages:
## [1] SPMpriors_1.0.4    FishLife_3.1.0     rfishbase_3.1.9.99 spict_1.3.7
## [5] TMB_1.9.10         corrplot_0.92      R.devices_2.17.2

## Functions
## -----------------
source("funcs.R")


## Variables
## -----------------
ndt <- 8
dteuler <- 1/ndt
server <- FALSE
if(server){
    mc.cores <- as.numeric(Sys.getenv('LSB_DJOB_NUMPROC'))
    input2 <- input <- as.character(commandArgs(TRUE))
    n.scen <- 12
    print(n.scen)
    r.scen <- as.integer(input2[2])
    print(r.scen)
    q.scen <- as.integer(input2[3])
    print(q.scen)
    sdi.scen <- as.integer(input2[4])
    sdc.scen <- as.integer(input2[4])
    print(sdi.scen)
    esb.type <- 0
    print(esb.type)
    prior.scen <- as.integer(input2[6])
    print(prior.scen)
    sdwn <- as.numeric(input2[7])
    print(sdwn)
    init.scen <- 0
    dbg <- 0
    RhpcBLASctl::blas_set_num_threads(mc.cores)
    datdir <- "../data"
    resdir <- "/work3/tobm/spmprod"
}else{
    require(unix)
    rlimit_as(1e12)
    sdwn <- 0.0
    n.scen <- 12
    r.scen <- 0
    q.scen <- 0
    sdi.scen <- 0
    sdc.scen <- 0
    esb.type <- 0
    add.scen <- 5
    init.scen <- 0
    dbg <- 0
    datdir <- "../data"
    resdir <- "/media/tobm/tobm_ext/outsourcing/spmprod/res/raw"
}


## Data
## -----------------
load(file.path(datdir,"stockData.RData"))
stocks <- names(stockData)
stocks <- stocks[stocks %in% c("COD-NS","HAD","POK","WHG-NS","PLE-NS","SOL-NS","TUR")]
nstocks <- length(stocks)
stockInfo <- read.csv(file.path(datdir,"stockInfo.csv"))
stockInfo <- stockInfo[match(stocks, stockInfo$stock),]


## Label
## -----------------
label <- paste0("sd_",sdwn, "_nScen_", n.scen, "_rScen_", r.scen,
                "_qScen_", q.scen, "_sdiScen_", sdi.scen, "_esb_",esb.type,
                "_add_",add.scen,"_ndt_",1/dteuler,"_ini_",init.scen)
print(label)


## Directories
## -----------------
dir.create(file.path(resdir), showWarnings = FALSE)
dir.create(file.path(resdir,"supp"), showWarnings = FALSE)
dir.create(file.path(resdir,"supp",label), showWarnings = FALSE)
dir.create(file.path(resdir,"fits"), showWarnings = FALSE)
dir.create(file.path(resdir,"fits",label), showWarnings = FALSE)
dir.create(file.path(resdir,"sum"), showWarnings = FALSE)
dir.create(file.path(resdir,"sum",label), showWarnings = FALSE)


## Models
## -----------------
models <- c("CM","KM","RM","PKRM","LKRM","IKRM")

get.priors <- FALSE
if(get.priors){
    rpriorsL <- vector("list",nstocks)
    for(i in 1:nstocks){
        stock.label <- paste0("fit_stock_",stocks[i])
        print(stock.label)
        rpriorsL[[i]] <- get.fishlife.prior(stockInfo$species[stockInfo$stock == stocks[i]])
    }
    rpriors <- data.frame(stocks, do.call(rbind, rpriorsL))
    rpriors$r_mu <- signif(exp(rpriors[,2]),2)
    rpriors[,c(1,4)]
    paste(rpriors[,1],collapse=", ")
    paste(rpriors[,4],collapse=", ")
}
rpriors <- data.frame(stocks = c("COD-NS", "HAD",
                                 "PLE-NS", "POK",
                                 "SOL-NS", "TUR",
                                 "WHG-NS"),
                      mu = c(0.35, 0.29, 0.25, 0.26,
                             0.27, 0.55, 0.38))


## Fitting
## -----------------
for(i in 1:nstocks){

    stock.label <- paste0("fit_stock_",stocks[i])
    print(stock.label)

    stocki <- which(names(stockData) == stocks[i])

    ## Data
    dat0 <- list()

    dat0$timeC <- stockData[[stocki]]$data$year
    dat0$obsC <- stockData[[stocki]]$data$C

    ## Abundance index
    if(esb.type == 0){ ## Default (stock weight)
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$ESB_raw_ov_sWeight_new
    }else if(esb.type == 1){ ## GAM + stock.wt
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$ESB_wStock
    }else if(esb.type == 2){  ## no interaction in GAM
        dat0$timeI <- stockData[[stocki]]$data$year + 0.5
        dat0$obsI <- stockData[[stocki]]$data$ESB_wCatch_av
    }else if(esb.type == 3){  ## no interaction in GAM + stock.wt
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$ESB_wStock_av
    }else if(esb.type == 4){  ## no year effect in GAM
        dat0$timeI <- stockData[[stocki]]$data$year + 0.5
        dat0$obsI <- stockData[[stocki]]$data$ESB_wCatch_ov
    }else if(esb.type == 5){  ## no year effect in GAM + stock.wt
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$ESB_wStock_ov
    }else if(esb.type == 6){  ## RAW
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$ESB_raw_new
    }else if(esb.type == 7){  ## RAW mean
        dat0$timeI <- stockData[[stocki]]$data$year + 0.5
        dat0$obsI <- stockData[[stocki]]$data$ESB_raw_av
    }else if(esb.type == 8){  ## GAM
        dat0$timeI <- stockData[[stocki]]$data$year + 0.5
        dat0$obsI <- stockData[[stocki]]$data$ESB_wCatch
    }else if(esb.type == 9){  ## SSB
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$SSB
    }else if(esb.type == 10){  ## TSB
        dat0$timeI <- stockData[[stocki]]$data$year
        dat0$obsI <- stockData[[stocki]]$data$TSB
    }else if(esb.type == 11){  ## default but catch weight
        dat0$timeI <- stockData[[stocki]]$data$year + 0.5
        dat0$obsI <- stockData[[stocki]]$data$ESB_raw_ov_new
    }

    names(stockData[[stocki]]$data)

    ind <- which(dat0$timeC == 2022)
    if(length(ind) > 0){
        dat0$timeC <- dat0$timeC[-ind]
        dat0$obsC <- dat0$obsC[-ind]
    }
    ind <- which(dat0$timeI == 2022)
    if(length(ind) > 0){
        dat0$timeI <- dat0$timeI[-ind]
        dat0$obsI <- dat0$obsI[-ind]
    }

    dat <- dat0

    set.seed(1234)

    if(any(dat$obsI < 0)){
        stop("Neg obs - CHECK")
    }

    ## Default spict
    inp0 <- dat

    ## Dteuler
    inp0$dteuler <- dteuler

    ## Prios
    ## ---------------
    if(n.scen == 0){
        inp0$priors$logn <- c(0,0,0)
    }else if(n.scen == 1){
        inp0$priors$logn <- c(log(2),2,1)
    }else if(n.scen == 2){
        inp0$priors$logn <- c(0,0,0)
        inp0$ini$logn <- log(2)
        inp0$phases$logn <- -1
    }else if(n.scen == 3){
        inp0 <- add.thorson.prior(inp0, stockInfo$order[stocki])
    }else if(n.scen == 4){
        inp0 <- add.winker.prior(inp0, genus[stocki], species[stocki], n = TRUE, r = FALSE)
        if(inp0$priors$logr[2] < 0.5) inp0$priors$logn[2] <- 0.5
    }else if(n.scen == 5){
        inp0 <- add.winker.prior(inp0, genus[stocki], species[stocki], n = TRUE, r = FALSE)
    }else if(n.scen == 6){
        inp0$priors$logn <- c(log(2),0.5,1)
    }else if(n.scen == 7){
        inp0$priors$logn <- c(0,0,0)
        inp0$ini$logn <- add.thorson.prior(inp0, stockInfo$order[stocki])$priors$logn[1]
        inp0$phases$logn <- -1
    }else if(n.scen == 8){
        inp0$priors$logn <- c(0,0,0)
        inp0$ini$logn <- log(1.001)
        inp0$phases$logn <- -1
    }else if(n.scen == 9){
        inp0$priors$logn <- c(0,0,0)
        inp0$ini$logn <- log(1.5)
        inp0$phases$logn <- -1
    }else if(n.scen == 10){
        inp0$priors$logn <- c(log(5),2,1)
        inp0$ini$logn <- log(5)
    }else if(n.scen == 11){
        inp0$priors$logn <- c(log(5),2,1)
        inp0$ini$logn <- log(5)
    }else if(n.scen == 12){
        inp0$priors$logn <- c(0,0,0)
        inp0$ini$logn <- log(5)
        inp0$phases$logn <- -1
    }

    ## r prior
    if(r.scen == 0){
        ## no prior
        inp0$priors$logr <- c(0,0,0)
    }else if(r.scen == 1){
        ## FishLife prior
        inp0 <- add.fishlife.prior(inp0, genus[stocki], species[stocki])
        ## sd is quite small
        if(inp0$priors$logr[2] < 0.5) inp0$priors$logr[2] <- 0.5
    }else if(r.scen == 2){
        ## FishLife prior
        inp0 <- add.fishlife.prior(inp0, genus[stocki], species[stocki])
    }else if(r.scen == 3){
        ## Add Henning's prior
        inp0 <- add.winker.prior(inp0, genus[stocki], species[stocki], n = FALSE, r = TRUE)
        ## sd is quite small
        if(inp0$priors$logr[2] < 0.5) inp0$priors$logr[2] <- 0.5
    }else if(r.scen == 4){
        ## Add Henning's prior
        inp0 <- add.winker.prior(inp0, genus[stocki], species[stocki], n = FALSE, r = TRUE)
    }else if(r.scen == 5){
        ## Customised FishBase prior
        inp0 <- add.fishlife.prior(inp0, genus[stocki], species[stocki])
        inp0$priors$logr[1] <- inp0$priors$logr[1] + inp0$ini$logn - log(2)
        inp0$priors$logr[2] <- 0.5
    }else if(r.scen == 6){
        ## Customised FishBase prior
        inp0 <- add.fishlife.prior(inp0, genus[stocki], species[stocki])
        inp0$priors$logr[1] <- inp0$priors$logr[1] + inp0$ini$logn - log(2)
        inp0$priors$logr[2] <- 1
    }



    ## additional priors
    if(add.scen == 1){
        inp0$priors$logbkfrac <- c(log(0.8),1,1)
    }else if(add.scen == 2){
        inp0$priors$logsdb <- c(log(0.15),0.5,1)
    }else if(add.scen == 3){
        inp0$priors$logsdb <- c(log(0.15),1,1)
    }else if(add.scen == 4){
        inp0$stdevfacC <- rep(1, length(inp0$timeC))
        inp0$stdevfacC[inp0$timeC < 1990] <- 2
    }else if(add.scen == 5){
        sd <- 1
        inp0$priors$logsdb <- c(log(0.15)-0.5*sd^2,sd,1)
    }

    inp0$priors$logalpha <- c(0,0,0)
    inp0$priors$logbeta <- c(0,0,0)
    inp0$priors$logsdm <- c(0,0,0)
    inp0$priors$logpsi <- c(0,0,0)
    inp0$priors$logsdK <- c(0,0,0)
    inp0$priors$logpsiK <- c(0,0,0)


    ## Fix parameters
    ## ----------------
    ## catchability scenarios
    if(q.scen == 0){
        inp0$ini$logq <- log(1); inp0$phases$logq <- -1
    }else if(q.scen == 1){
        inp0$priors$logq <- c(log(1),2,1)
    }else if(q.scen == 2){
    }

    ## Fix sdi
    if(sdi.scen == 0){
        inp0$ini$logsdi <- log(sdwn + 1e-3)
        inp0$phases$logsdi <- -1
    }else if(sdi.scen == 1){
        inp0$priors$logsdi <- c(log(sdwn + 1e-3),2,1)
    }else if(sdi.scen == 2){
    }

    ## Fix sdc
    if(sdc.scen == 0){
        inp0$ini$logsdc <- log(sdwn + 1e-3)
        inp0$phases$logsdc <- -1
    }else if(sdc.scen == 1){
        inp0$priors$logsdc <- c(log(sdwn + 1e-3),2,1)
    }else if(sdc.scen == 2){
    }


    ## Better initial values
    ## ----------------
    if(init.scen == 0){
        init.vals <- c(0.2,0.2)
    }else if(init.scen == 1){
        init.vals <- c(0.05,0.05)
    }else if(init.scen == 2){
        init.vals <- c(0.5,0.5)
    }else if(init.scen == 3){
        init.vals <- c(0.05,0.5)
    }else if(init.scen == 4){
        init.vals <- c(0.5,0.05)
    }
    ## for AR1 on m
    inp0$ini$logsdm <- log(init.vals[1])
    inp0$ini$logpsi <- log(init.vals[2])
    ## for AR1 on K
    inp0$ini$logsdK <- log(init.vals[1])
    inp0$ini$logpsiK <- log(init.vals[2])
    ## for AR1 on q
    inp0$ini$logsdq <- log(init.vals[1])
    inp0$ini$logpsiq <- log(init.vals[2])

    print(init.vals)

    inp0$msytype <- "d"



    fit <- list()
    for(modi in 1:length(models)){

        mod.lab <- models[modi]
        print(mod.lab)

        if(mod.lab == "CM"){

            ## Constant productivity
            ## -----------------
            inp <- check.inp(inp0)

        }else if(mod.lab == "KM"){

            ## Time-variant m + K (constant r)
            ## -----------------
            inp <- inp0
            inp$tvKConsR <- TRUE
            inp <- check.inp(inp)

        }else if(mod.lab == "RM"){

            ## Time-variant m + r (constant K)
            ## -----------------
            inp <- inp0
            inp$timevaryinggrowth <- TRUE
            inp <- check.inp(inp)

        }else if(mod.lab == "PKRM"){

            ## Time-variant m + K + r (constant slope between r and K)
            ## -----------------
            inp <- inp0
            inp$tvPropChange <- TRUE
            inp <- check.inp(inp)
        }else if(mod.lab == "LKRM"){

            ## Time-variant m + r + K (coupled AR(1)s + LM on log scale)
            ## -----------------
            inp <- inp0
            inp$tvNonPropChange <- TRUE
            inp <- check.inp(inp)

        }else if(mod.lab == "IKRM"){

            ## Time-variant m + r + K (uncoupled AR(1)s)
            ## -----------------
            inp <- inp0
            inp$timevaryinggrowth <- TRUE
            inp$timevaryingK <- TRUE
            inp <- check.inp(inp)
        }

        if(n.scen == 11){
            if(mod.lab == "CM" && stock.label != "fit_stock_WIT"){
                inp$priors$logn <- c(log(5),2,1)
                inp$ini$logn <- log(5)
            }else{
                inp$priors$logn <- c(log(2),2,1)
                inp$ini$logn <- log(2)
            }
        }

        ## fit spict
        set.seed(1234)
        fit <- try(fit.spict(inp), silent = TRUE)

        if(inherits(fit, "try-error") || fit$opt$convergence != 0 ||
               any(is.infinite(fit$sd)) || any(is.na(fit$sd)) ||
               any(is.infinite(exp(sqrt(diag(fit$cov.fixed))))) ||
           any(is.infinite(get.par("logBBmsy", fit, exp = TRUE)[,2]))){
            print("Not converged according to hard criteria!")
        }else if(any(is.na(fit$sd)) ||
           any(is.na(exp(sqrt(diag(fit$cov.fixed))))) ||
           any.invalid(get.par("logBBmsy", fit, exp = TRUE)) ||
           any.invalid(get.par("logFFmsy", fit, exp = TRUE)) ||
           any.invalid(get.par("logFnotS", fit, exp = TRUE)) ||
           any.invalid(get.par("logB", fit, exp = TRUE)) ||
           any.invalid(get.par("logrvec", fit, exp = TRUE)) ||
           any.invalid(get.par("logKvec", fit, exp = TRUE)) ||
           any.invalid(get.par("logmvec", fit, exp = TRUE)) ||
           any.invalid(get.par("logFmsyvec", fit, exp = TRUE)) ||
           any.invalid(get.par("logBmsyvec", fit, exp = TRUE)) ||
           any.invalid(get.par("logFmsy", fit, exp = TRUE)) ||
           any.invalid(get.par("logBmsy", fit, exp = TRUE)) ||
           any(is.infinite(fit$sd)) ||
           ## meaningful CVs
           any(is.infinite(exp(sqrt(diag(fit$cov.fixed))))) ||
           any.invalid2(get.par("logBBmsy", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logFFmsy", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logFnotS", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logB", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logrvec", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logKvec", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logmvec", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logFmsyvec", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logBmsyvec", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logFmsy", fit, exp = TRUE),Inf) ||
           any.invalid2(get.par("logBmsy", fit, exp = TRUE),Inf) ||
           ## max(get.par("logFFmsy", fit, exp = TRUE)[,2],
           ##     na.rm = TRUE) < 0.1 ||
           ## Order of magnitude
           any.invalid3(get.par("logBBmsy", fit, exp = TRUE)) ||
           any.invalid3(get.par("logFFmsy", fit, exp = TRUE)) ||
           any.invalid3(get.par("logFnotS", fit, exp = TRUE)) ||
           any.invalid3(get.par("logB", fit, exp = TRUE)) ||
           any.invalid3(get.par("logrvec", fit, exp = TRUE)) ||
           any.invalid3(get.par("logKvec", fit, exp = TRUE)) ||
           any.invalid3(get.par("logmvec", fit, exp = TRUE)) ||
           any.invalid3(get.par("logFmsyvec", fit, exp = TRUE)) ||
           any.invalid3(get.par("logBmsyvec", fit, exp = TRUE)) ||
           any.invalid3(get.par("logFmsy", fit, exp = TRUE)) ||
           any.invalid3(get.par("logBmsy", fit, exp = TRUE)) ||
           any.invalid3(get.par("logsdb", fit, exp = TRUE)) ||
           any.invalid3(get.par("logsdf", fit, exp = TRUE)) ||
           any.invalid3(get.par("logpsiK", fit, exp = TRUE)) ||
           any.invalid3(get.par("logsdK", fit, exp = TRUE)) ||
           any.invalid3(get.par("logpsim", fit, exp = TRUE)) ||
           any.invalid3(get.par("logsdm", fit, exp = TRUE)) ||
           sd(get.par("logBBmsy", fit, exp = TRUE)[get.par("logBBmsy", fit, exp = TRUE)[,2] > quantile(get.par("logBBmsy", fit, exp = TRUE)[,2],0.025) & get.par("logBBmsy", fit, exp = TRUE)[,2] < quantile(get.par("logBBmsy", fit, exp = TRUE)[,2],0.975), 2]) < 0.06
           ){
            print("Not converged according to super hard criteria!")
        }

        if(!inherits(fit, "try-error")){

            print(paste0("AIC: ", round(AIC(fit))))
            print(paste0("m = ", round(get.par("logm", fit, exp = TRUE)[,2]/1e3)))
            print(paste0("K = ", round(get.par("logK", fit, exp = TRUE)[,2]/1e3)))

            fit <- try(calc.osa.resid(fit))

            if(inherits(fit, "try-error")){

                print(paste0("Model not coverged for: ", stocks[i], " - ", models[modi]))

            }else{

                write.table(capture.output(fit),
                            file = file.path(resdir,"sum",label,paste0(stock.label,"_def.txt")))

                ## Default fit plot
                png(file.path(resdir,"supp",label,
                              paste0(stock.label,"_fit_",mod.lab,".png")),
                    width = 1000, height = 1000, res = 120)
                try(plot(fit),silent=TRUE)
                dev.off()

                ## MacCall plot
                png(file.path(resdir,"supp",label,
                              paste0(stock.label,"_maccall_",mod.lab,".png")),
                    width = 800, height = 1000, res = 120)
                try(plotextra.maccall(fit),silent=TRUE)
                dev.off()

                ## TV plot
                png(file.path(resdir,"supp",label,
                              paste0(stock.label,"_tv_",mod.lab,".png")),
                    width = 800, height = 1200, res = 120)
                try(plotextra.tv(fit),silent=TRUE)
                dev.off()

                ## Observation residuals
                png(file.path(resdir,"supp",label,
                              paste0(stock.label,"_resid_",mod.lab,".png")),
                    width = 1000, height = 1200, res = 120)
                try(plotspict.diagnostic(fit),silent=TRUE)
                dev.off()

                ## Correlation matrix
                ind <- rownames(fit$cov.fixed) %in% c("logm","logK","logn")
                refcor <- cov2cor(fit$cov.fixed[ind,ind])
                png(file.path(resdir,"supp",label,
                              paste0(stock.label,"_corr_",mod.lab,".png")),
                    width = 700, height = 800, res = 120)
                corrplot(refcor, method = c("number"), type = "upper")
                dev.off()
            }
        }

        ## Save
        save(fit, file = file.path(resdir,"fits",label,
                                   paste0(stock.label,"_",mod.lab,".RData")),
             version = 2)
    }
}
