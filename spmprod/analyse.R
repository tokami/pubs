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


## This scripts analyses fits and creates tables and figures.


## Packages --------------------------
## remotes::install_github("tokami/spict/spict", ref="mChangeK2", force = TRUE)
## remotes::install_version("Matrix", version = "1.6-1.1")
## or: remove.packages("Matrix")
## and: install.packages("TMB", dependencies = TRUE)
require(spict)
require(car)
require(MASS)

## Additional functions
source("spmprod_0_funcs.R")


## Paths --------------------------
datdir <- "../data/"
resdir <- "../res"
fitdir <- "/media/tobm/tobm_ext/outsourcing/spmprod/res/raw"


## Variables --------------------------
final.year <- 2021
compile.rdata <- FALSE
save.tables <- TRUE
save.figs <- TRUE
sdigits <- 2
scaleMetric <- 1e3
ndt <- 8
dteuler <- 1/ndt
n.scen <- 2
esb.type <- 0
init.scen <- 0

sdwn <- 0
r.scen <- 0
q.scen <- 0
sdi.scen <- 0
sdc.scen <- 0
add.scen <- 5
label <- paste0("sd_",sdwn,"_nScen_", n.scen, "_rScen_", r.scen,
                "_qScen_", q.scen,
                "_sdiScen_", sdi.scen, "_esb_",esb.type,"_add_",add.scen,
                "_ndt_",1/dteuler,"_ini_",init.scen)
dir(file.path(fitdir,"fits"))
print(label)

if(n.scen %in% c(1,2,6)){
    nPar <- 2
}else if(n.scen %in% c(3,7,9)){
    nPar <- 1.5
}else if(n.scen %in% c(8)){
    nPar <- 1.000001
}else if(n.scen %in% c(10,11,12)){
    nPar <- 5
}


load(file.path(datdir,"stockData.RData")) ## load stockData
stocks <- names(stockData)
stocks.true <- stocks[stocks %in% c("COD-NS","HAD","POK","WHG-NS","PLE-NS","SOL-NS","TUR")]
print(stocks.true)
stocks <- sapply(strsplit(stocks.true, "-"), function(x) x[1])
stockInfo <- read.csv(file.path(datdir,"stockInfo.csv"))
stockInfo <- stockInfo[match(stocks.true, stockInfo$stock),]


envData <- read.csv(file.path(datdir,"env_data.csv"))

dir.create(file.path(resdir,"figs"), showWarnings = FALSE)
dir.create(file.path(resdir,"tabs"), showWarnings = FALSE)
dir.create(file.path(resdir,"rdata"), showWarnings = FALSE)
dir.create(file.path(resdir,"figs",label), showWarnings = FALSE)
dir.create(file.path(resdir,"tabs",label), showWarnings = FALSE)
dir.create(file.path(resdir,"figs",label,"fits"), showWarnings = FALSE)
dir.create(file.path(resdir,"rdata",label), showWarnings = FALSE)

if(compile.rdata){
    files0 <- dir(file.path(fitdir,"fits",label))[grep(".RData", dir(file.path(fitdir,"fits",label)))]
    labels <- sapply(strsplit(unlist(strsplit(files0, ".RData")), "fit_stock_"),"[[",2)
    stocks0 <- c("COD-NS", "HAD", "POK", "WHG-NS",
                 "PLE-NS", "SOL-NS", "TUR")
    common.names <- c("Cod","Haddock","Saithe","Whiting",
                      "Plaice","Sole","Turbot")
    ## "PLE-EC", "SOL-EC", "WIT"
    stocks.true <- unique(sapply(strsplit(labels, "_"), "[[", 1))
    stocks.true <- stocks.true[match(stocks0,stocks.true)]
    stocks <- sapply(strsplit(stocks.true, "-"), function(x) x[1])
    nstocks <- length(stocks.true)
    stockData <- stockData[match(stocks.true, names(stockData))]
    stockInfo <- stockInfo[match(stocks.true, stockInfo$stock),]
    mods <- unique(sapply(strsplit(labels, "_"), "[[", 2))
    modLabels <- c("CM","RM", "KM","PKRM","LKRM","IKRM")
    mods <- mods[match(modLabels, mods)]
    stocks <- common.names
}else{
    stocks.true <- c("COD-NS", "HAD", "POK", "WHG-NS",
                 "PLE-NS", "SOL-NS", "TUR")
    common.names <- c("Cod","Haddock","Saithe","Whiting",
                      "Plaice","Sole","Turbot")
    nstocks <- length(stocks.true)
    stocks <- sapply(strsplit(stocks.true, "-"), function(x) x[1])
    stockData <- stockData[match(stocks.true, names(stockData))]
    stockInfo <- stockInfo[match(stocks.true, stockInfo$stock),]
    mods <- modLabels <- c("CM","RM", "KM","PKRM","LKRM","IKRM")
    mods <- mods[match(modLabels, mods)]
    stocks <- common.names
}


pars <- c("time","logK","logm","logr","logKvec","logmvec","logrvec","logBBmsy",
          "logFFmsynotS","logB","logFnotS","logBmsy","logFmsy",
          "logMSYvec","logFmsyvec","logBmsyvec","kobe",
          "logsdm","logpsi","logsdK","logpsiK","mkb",
          "logsdb","logsdf","cmsy","logn")

max(sapply(stockData, function(x) min(x$data$year)))

if(!all(names(stockData) == stockInfo$stock)) stop("ices data and stockData does not match!")
if(compile.rdata){

    df <- df2 <- nonConv2 <- nonConv1 <- NULL
    presList <- resList <- lrList <- vector("list", length(stocks.true))

    for(i in 1:length(stocks.true)){
        presList[[i]] <- resList[[i]] <- vector("list", length(mods))
        lrList[[i]] <- vector("list", 3)
        fitList <- vector("list", length(mods))
        for(j in 1:length(mods)){
            resList[[i]][[j]] <- vector("list", length(pars))
            tmp <- try(load(file.path(fitdir,"fits",label,files0[grep(paste0(stocks.true[i],"_",mods[j],"\\b"), files0)])))

            fitList[[j]] <- fit

            if(inherits(tmp, "try-error") ||
               inherits(fit, "try-error") ||
               fit$opt$convergence != 0 ||
               any(is.na(fit$sd)) ||
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
               (!is.na(any.invalid2(get.par("logsdK", fit, exp = TRUE),Inf)) &&
                any.invalid2(get.par("logsdK", fit, exp = TRUE),Inf)) ||
               (!is.na(any.invalid2(get.par("logpsiK", fit, exp = TRUE),Inf)) &&
                any.invalid2(get.par("logpsiK", fit, exp = TRUE),Inf)) ||
               (!is.na(any.invalid2(get.par("logsdm", fit, exp = TRUE),Inf)) &&
                any.invalid2(get.par("logsdm", fit, exp = TRUE),Inf)) ||
               (!is.na(any.invalid2(get.par("logpsi", fit, exp = TRUE),Inf)) &&
                any.invalid2(get.par("logpsi", fit, exp = TRUE),Inf)) ||
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
               any.invalid2(get.par("logBmsy", fit, exp = TRUE),Inf)
               ){
                for(k in 1:length(pars)) resList[[i]][[j]][[k]] <- matrix(NA, 1, 5)
                nonConv1 <- c(nonConv1, paste0(stocks.true[i],"_",modLabels[j]))
                df <- rbind(df, c(stocks[i],modLabels[j], rep(NA, 8)))
                df2 <- rbind(df2, c(stocks[i],modLabels[j], rep(NA, 15)))
            }else{
                if(
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
                   any.invalid3(get.par("logsdm", fit, exp = TRUE))
                   ){
                    nonConv2 <- c(nonConv2, paste0(stocks.true[i],"_",modLabels[j]))
                }

                png(file.path(resdir,"figs",label,"fits",paste0("fit_main_",stocks.true[i],"_",modLabels[j],".png")),width = 1000, height = 1100, res = 120)
                plot(fit)
                dev.off()
                png(file.path(resdir,"figs",label,"fits",paste0("fit_main2_",stocks.true[i],"_",modLabels[j],".png")), width = 1000, height = 800, res = 120)
                plot2(fit)
                dev.off()
                png(file.path(resdir,"figs",label,"fits",paste0("fit_diag_",stocks.true[i],"_",modLabels[j],".png")), width = 1000, height = 1000, res = 120)
                plotspict.diagnostic(fit)
                dev.off()
                fit2 <- try(calc.process.resid(fit))
                if(!inherits(fit2, "try-error")){
                    png(file.path(resdir,"figs",label,"fits",paste0("fit_diag_proc_",stocks.true[i],"_",modLabels[j],".png")), width = 900, height = 1000, res = 120)
                    plotspict.diagnostic.process(fit2)
                    dev.off()
                    presList[[i]][[j]] <- fit2$process.resid
                }
                ## dteuler
                dt <- fit$inp$dteuler
                ## Results table
                df <- rbind(df, c(stocks[i], modLabels[j],
                                  AIC(fit),
                                  AICc(fit),
                                  fit$opt$objective,
                                  paste0(
                                      signif(get.par("logBBmsy", fit, exp = TRUE)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),2], sdigits)," (",
                                      paste(signif(diff(get.par("logBBmsy", fit, exp = TRUE, CI = 0.75)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)]),sdigits),collapse = "-"),")"),
                                  paste0(
                                      signif(get.par("logFFmsynotS", fit, exp = TRUE)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),2], sdigits)," (",
                                      paste(signif(diff(get.par("logFFmsynotS", fit, exp = TRUE, CI = 0.75)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)]),sdigits),collapse = "-"),")"),
                                  signif(get.TAC(fit)/scaleMetric,sdigits),
                                  paste0(
                                      signif(get.par("logBBmsy", fit, exp = TRUE)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),2], sdigits)," (",
                                      paste(signif(get.par("logBBmsy", fit, exp = TRUE, CI = 0.95)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)],sdigits),collapse = "-"),")"),
                                  paste0(
                                      signif(get.par("logFFmsynotS", fit, exp = TRUE)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),2], sdigits)," (",
                                      paste(signif(get.par("logFFmsynotS", fit, exp = TRUE, CI = 0.95)[
                                          which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)],sdigits),collapse = "-"),")")
                                  )
                            )
                df2 <- rbind(df2, c(stocks[i], modLabels[j],
                                    paste0(
                                        signif(get.par("logB", fit, exp = TRUE)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),2]/scaleMetric,sdigits)," (",
                                        paste(signif(diff(get.par("logB", fit, exp = TRUE, CI = 0.75)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)]/scaleMetric),sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logFnotS", fit, exp = TRUE)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),2],sdigits)," (",
                                        paste(signif(diff(get.par("logFnotS", fit, exp = TRUE, CI = 0.75)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)]),sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logBmsy", fit, exp = TRUE)[,2]/scaleMetric,sdigits)," (",
                                        paste(signif(diff(get.par("logBmsy", fit, exp = TRUE, CI = 0.75)[,c(1,3)]/scaleMetric),sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logFmsy", fit, exp = TRUE)[,2],sdigits)," (",
                                        paste(signif(diff(get.par("logFmsy", fit, exp = TRUE, CI = 0.75)[,c(1,3)]),sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logB", fit, exp = TRUE)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),2]/scaleMetric,sdigits)," (",
                                        paste(signif(get.par("logB", fit, exp = TRUE, CI = 0.95)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)]/scaleMetric,sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logFnotS", fit, exp = TRUE)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),2],sdigits)," (",
                                        paste(signif(get.par("logFnotS", fit, exp = TRUE, CI = 0.95)[
                                            which(fit$inp$time == fit$inp$timerangeObs[2]),c(1,3)],sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logBmsy", fit, exp = TRUE)[,2]/scaleMetric,sdigits)," (",
                                        paste(signif(get.par("logBmsy", fit, exp = TRUE, CI = 0.95)[,c(1,3)]/scaleMetric,sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logFmsy", fit, exp = TRUE)[,2],sdigits)," (",
                                        paste(signif(get.par("logFmsy", fit, exp = TRUE, CI = 0.95)[,c(1,3)],sdigits),collapse = "-"),")"),
                                    paste0(
                                        signif(get.par("logMSY", fit, exp = TRUE)[,2]/scaleMetric,sdigits)," (",
                                        paste(signif(get.par("logMSY", fit, exp = TRUE, CI = 0.95)[,c(1,3)]/scaleMetric,sdigits),collapse = "-"),")"),

                                    round((get.mean.par("logrvec",fit, 2021) - get.mean.par("logrvec",fit, 1983)) /
                                          get.mean.par("logrvec",fit, 1983) * 100),

                                    round((get.mean.par("logKvec",fit, 2021) - get.mean.par("logKvec",fit, 1983)) /
                                          get.mean.par("logKvec",fit, 1983) * 100),

                                    round((get.mean.par("logmvec",fit, 2021) - get.mean.par("logmvec",fit, 1983)) /
                                          get.mean.par("logmvec",fit, 1983) * 100),

                                    round(get.mean.par("logrvec",fit, 2021) /
                                          get.mean.par("logrvec",fit, 1983) * 100),

                                    round(get.mean.par("logKvec",fit, 2021) /
                                          get.mean.par("logKvec",fit, 1983) * 100),

                                    round(get.mean.par("logmvec",fit, 2021) /
                                          get.mean.par("logmvec",fit, 1983) * 100)
                                    ))

                ## For plots
                for(k in 1:length(pars)){
                    if(pars[k] == "kobe"){
                        ## cl <- try(spict:::make.rpellipse(fit))
                        cl <- try(make.rpellipse2(fit))
                        if(inherits(cl, "try-error")){
                            cl <- matrix(NA, 300, 2)
                        }
                        resList[[i]][[j]][[k]] <- cl
                    }else if(pars[k] == "time"){
                        resList[[i]][[j]][[k]] <- fit$inp$time
                    }else if(pars[k] %in% c("logm","logmvec")){
                        resList[[i]][[j]][[k]] <- get.par(pars[k], fit, exp = TRUE) ## * dt
                    }else if(pars[k] %in% c("mkb")){
                        resList[[i]][[j]][[k]] <- get.par(pars[k], fit, exp = FALSE)
                    }else if(pars[k] %in% c("cmsy")){
                        resList[[i]][[j]][[k]] <- matrix(NA, length(fit$inp$indest), 5)
                        resList[[i]][[j]][[k]][,2] <- fit$report$Cpredsub[fit$inp$indest] / (get.par("logmvec",fit, exp = TRUE)[fit$inp$indest,2] * fit$inp$dteuler)
                    }else{
                        resList[[i]][[j]][[k]] <- get.par(pars[k], fit, exp = TRUE)
                    }
                }
            }
            names(resList[[i]][[j]]) <- pars

        }## end model loop
        names(resList[[i]]) <- modLabels
        names(fitList) <- modLabels
        df <- rbind(df, c(stocks[i], "ICES",
                          NA, NA, NA,
                          signif(stockInfo$BBmsy[i], sdigits),
                          signif(stockInfo$FFmsy[i], sdigits),
                          NA,
                          signif(stockInfo$BBmsy[i], sdigits),
                          signif(stockInfo$FFmsy[i], sdigits)))
        df2 <- rbind(df2, c(stocks[i], "ICES",
                            signif(stockInfo$ssb[i] / scaleMetric, sdigits),
                            signif(stockInfo$Fbar[i], sdigits),
                            signif(stockInfo$EQ_Bmsy[i] / scaleMetric, sdigits),
                            signif(stockInfo$Fmsy[i], sdigits),
                            signif(stockInfo$ssb[i] / scaleMetric, sdigits),
                            signif(stockInfo$Fbar[i], sdigits),
                            signif(stockInfo$EQ_Bmsy[i] / scaleMetric, sdigits),
                            signif(stockInfo$Fmsy[i], sdigits),
                            signif(stockInfo$EQ_MSY[i] / scaleMetric, sdigits),
                            NA, NA, NA, NA, NA, NA))

        lrList[[i]][[1]] <- try(signif(lrtest.spict(fitList[1], fitList[-1]),
                                       sdigits), silent = TRUE)
        lrList[[i]][[2]] <- try(signif(lrtest.spict(fitList[5], fitList[-c(5,6)]),
                                       sdigits), silent = TRUE)
        lrList[[i]][[3]] <- try(signif(lrtest.spict(fitList[6], fitList[-6]),
                                       sdigits), silent = TRUE)

        fini <- length(which(sapply(fitList, function(x) !all(is.na(x)) &&
                                             !inherits(x, "try-error"))))
        widthi <- fini * (1700 / 6)
        png(file.path(resdir,"figs",label,"fits",
                      paste0("maccall_",stocks.true[i],".png")),
            width = widthi, height = 600, res = 120)
        layout(matrix(1:(fini*2), 2, fini))
        par(mar = c(1,1.5,0,1), oma = c(3.5,3,1.5,1))
        for(k in 1:length(mods)){
            if(k == 1) labels <- TRUE else labels <- FALSE
            tmp <- try(plotextra.maccall(fitList[[k]], bcor = 1e3, prodcor = 1e3,
                                         add = TRUE, labels = labels,
                                         legend = ifelse(k == fini,
                                                         TRUE, FALSE),
                                         lab2 = mods[k]),
                       silent=TRUE)
        }
        dev.off()

        xlimList <- ylimList1 <- ylimList2 <- vector("list",length(mods))
        for(k in 1:length(mods)){
            tmp <- try(get.maccall.axes(fitList[[k]], scaled = TRUE))
            if(!inherits(tmp, "try-error")){
                xlimList[[k]] <- tmp$xlim
                ylimList1[[k]] <- tmp$ylim1
                ylimList2[[k]] <- tmp$ylim2
            }
        }
        xlim <- range(c(0,unlist(xlimList)))
        ylim1 <- range(c(0,unlist(ylimList1)))
        ylim2 <- range(c(0,unlist(ylimList2)))

        fini <- max(which(sapply(fitList, function(x) !all(is.na(x)) &&
                                          !inherits(x, "try-error"))))
        widthi <- fini * (1700 / 6)
        png(file.path(resdir,"figs",label,"fits",
                      paste0("maccall2_",stocks.true[i],".png")),
            width = widthi, height = 600, res = 120)
        layout(matrix(1:12, 2, 6))
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(4.5,4,1,1))
        for(k in 1:length(mods)){
            labels <- FALSE
            if(k == 1) axes <- TRUE else axes <- FALSE
            tmp <- try(plotextra.maccall(fitList[[k]],
                                         add = TRUE, labels = labels,
                                         axes = axes, xlim = xlim, ylim1 = ylim1,
                                         ylim2 = ylim2, legend = ifelse(k == fini,
                                                                        TRUE, FALSE),
                                         lab2 = mods[k],
                                         scaled = TRUE),
                       silent=TRUE)
            if(inherits(tmp, "try-error")){
                plot.new()
                plot.new()
            }
            if(k == 1){
                mtext(expression("rel. dB" ~ dt^{-1}), 2, 1.7,
                      outer = TRUE, adj = 0.88)
                mtext(expression("rel. dB" ~ dt^{-1} ~ B^{-1}), 2, 2.5)
                mtext("B/K", 1, 2.5, outer = TRUE)
            }
        }
        dev.off()

        save(fitList,
             file = file.path(resdir, "rdata", label,
                              paste0("fitList_",label,"_",stocks.true[i],".RData")),
             version = 2)

    } ## end stock loop
    names(resList) <- stocks
    names(lrList) <- stocks

    df <- as.data.frame(df)
    colnames(df) <- c("stock", "model","AIC","AICc",
                      "LogLik","BBmsy","FFmsy","TAC","BBmsy2", "FFmsy2")
    rownames(df) <- NULL
    seli <- colnames(df) %in% c("AIC","AICc","LogLik")
    df <- cbind(df[,!seli], apply(df[,seli], 2, as.numeric))


    df2 <- as.data.frame(df2)
    colnames(df2) <- c("stock", "model","B","F","Bmsy","Fmsy","B2","F2", "Bmsy2", "Fmsy2",
                       "MSY","rChange","KChange","mChange","rDiff","KDiff","mDiff")
    rownames(df2) <- NULL
    seli <- colnames(df2) %in% c("rChange","KChange","mChange","rDiff","KDiff","mDiff")
    df2 <- cbind(df2[,!seli], apply(df2[,seli], 2, as.numeric))

    if(save.tables){
        write.csv(df, file.path(resdir, "tabs", label, paste0("results_",label,".csv")))
        write.csv(df2, file.path(resdir, "tabs", label, paste0("results_abs_",label,".csv")))
    }


    ## save all results
    save(resList, presList, lrList, nonConv2, nonConv1,
         file = file.path(resdir, "rdata", label, paste0("resList_",label,".RData")),
         version = 2)

}else{

    load(file.path(resdir, "rdata", label, paste0("resList_",label,".RData")))
    df <- read.csv(file.path(resdir, "tabs", label,
                             paste0("results_",label,".csv")))
    df2 <- read.csv(file.path(resdir, "tabs", label,
                              paste0("results_abs_",label,".csv")))

}


## Selected models
selStocks <- seq(stocks)
allMods <- lapply(seq(nstocks), function(x) mods)
selMods <- lapply(seq(nstocks), function(x) mods)  ## 6
names(selMods) <- stocks
rMods <- lapply(seq(nstocks), function(x) c("CM","RM"))

## Best models (AICc)
bestMods <- vector("list",length(stocks))
resList2 <- vector("list", length(stocks))
for(i in 1:length(stocks)){
    tmp <- df[which(df$stock == stocks[i]),]
    indi <- which(!is.na(tmp$AICc) & tmp$AICc == min(tmp$AICc, na.rm = TRUE))
    if(length(indi) > 1) browser()
    if(any(!is.na(indi))){
        bestMods[[i]] <- tmp$model[indi]
        resList2[[i]] <- resList[[i]][[indi]]
    }
}
cbind(stocks,
      unlist(bestMods))


nmods <- length(mods)
npars <- c(5,7,7,7,8,9)

df3 <- c("Stock","Model","LogLik","AICc","deltaAICc")

for(i in 1:length(stocks)){
    tmp <- df[which(df$stock == stocks[i]),]
    tmp <- tmp[which(tmp$model != "ICES"),]
    tmpi <- rep(NA, nmods)
    tmpi[modLabels %in% bestMods[[i]]] <- "x"
    df3 <- rbind(df3, cbind(rep(stocks[i], nmods),
                            modLabels,
                            signif(tmp[,"LogLik"],3),
                            signif(tmp[,"AICc"],3),
                            signif(tmp[,"AICc"] -
                                   min(tmp[,"AICc"], na.rm = TRUE),3)
                            ))
}
rownames(df3) <- NULL
colnames(df3) <- NULL

df3[is.na(df3)] <- ""

df3

tmp <- df3[-1,]
tmp <- tmp[which(tmp[,4] != ""),]
colnames(tmp) <- c("Stock", "Model", "logLik", "AIC\\textsubscript{c}",
                   "${\\Delta}$AIC\\textsubscript{c}")

tmp[,2] <- paste0(tmp[,2], ifelse(tmp[,5]  == "0", "\\textsuperscript{*}", ""))

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Log likelihood value, AIC\\textsubscript{c}, and ${\\Delta}$AIC\\textsubscript{c} for all stocks and models. Asterisks indicate the best model for each stock according to AIC\\textsubscript{c}.",
                                                     label = "tab:aic",
                                                     align="ll|lccc"),
                                      file = file.path(resdir, "tabs", label,
                                                       paste0("table_aic_",label,".tex")),
                                      include.rownames = FALSE, include.colnames = TRUE,
                                      sanitize.text.function = identity,
                                      hline.after = c(-1,nrow(tmp)),
                                      caption.placement = "top",
                                      floating=TRUE,
                                      booktabs=TRUE,
                                      size="\\fontsize{8pt}{8pt}\\selectfont"
                                      )


tmp <- as.data.frame(df3[-1,])
colnames(tmp) <- df3[1,]

## Convergence
(nstocks * nmods) - sum(tmp$LogLik == "")
(nstocks * nmods)
round((sum(tmp$LogLik != "")) / (nstocks * nmods)*100)
length(nonConv1)
nonConv1

length(nonConv2)
nonConv2


## Non-converged models
tmpi <- lapply(mods, function(x) stocks[tmp$LogLik[tmp$Model == x] == ""])
names(tmpi) <- mods
tmpi

nonConv1
length(nonConv2)
nonConv2

## Remaining
nstocks * nmods - sum(tmp$LogLik == "") - length(nonConv2)
round((sum(tmp$LogLik != "") - length(nonConv2))/ (nstocks * nmods)*100)

## By stock
sapply(stocks, function(x) sum(tmp$LogLik[tmp$Stock == x] != ""))
sapply(stocks, function(x) tmp$Model[tmp$Stock == x][tmp$LogLik[tmp$Stock == x] != ""])

tmpi

cbind(stocks, unlist(bestMods))

range(sapply(stocks, function(x) as.numeric(tmp$deltaAICc[tmp$Stock == x][1])))
sapply(stocks, function(x) as.numeric(tmp$deltaAICc[tmp$Stock == x][1]))

signif(mean(sapply(stocks, function(x) as.numeric(tmp$deltaAICc[tmp$Stock == x][-1])),
            na.rm = TRUE),2)

sapply(stocks, function(x) as.numeric(tmp$deltaAICc[tmp$Stock == x][-1]))


max(as.numeric(tmp$deltaAICc[tmp$Model %in% mods[-1]]), na.rm = TRUE)
signif(mean(as.numeric(tmp$deltaAICc[tmp$Model %in% mods[-1]]), na.rm = TRUE),2)
signif(mean(as.numeric(tmp$deltaAICc[tmp$Model %in% mods[1]]), na.rm = TRUE),2)

nmods <- length(mods)

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("proc_resid_b_",label,".png")),
        width = 1800, height = 1500, res = 120)
par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,3,3))
layout(matrix(1:(nmods*nstocks), ncol = nmods, nrow = nstocks))
min(sapply(presList, function(x) sapply(x, function(x) min(x$time))))
xlimi <- c(1957,2022)
ylimi <- c(-1.4,1.4)
for(j in 1:nmods){
for(i in 1:nstocks){
        if(!all(is.null(presList[[i]][[j]]$time))){
            xaxt <- ifelse(i == nstocks, "s", "n")
            yaxt <- ifelse(j == 1, "s", "n")
            plot(presList[[i]][[j]]$time,
                 presList[[i]][[j]]$B,
                 ty = 'n',
                 xlim = xlimi, ylim = ylimi,
                 xaxt = xaxt, yaxt = yaxt,
                 xlab = "", ylab = "")
            if(modLabels[j] == bestMods[[i]]){
                rect(par("usr")[1], par("usr")[3],
                     par("usr")[2], par("usr")[4],
                     col = adjustcolor("goldenrod2",0.2))
            }
            if(paste0(stocks.true[i],"_",mods[j]) %in% nonConv2){
                adji <- 0.5
                colii <- gray(0.6)
                flagi <- TRUE
            }else{
                adji <- 1
                colii <- "black"
                flagi <- FALSE
            }
            coli <- adjustcolor(c("chartreuse3","firebrick3")[(presList[[i]][[j]]$B < 0) + 1], adji)
            abline(h = 0, col = adjustcolor(gray(0.7),adji))
            lines(presList[[i]][[j]]$time,
                  presList[[i]][[j]]$B,
                  ty = 'b', pch = NA, col = adjustcolor(gray(0.7),adji),
                  lwd = 1.5)
            points(presList[[i]][[j]]$time,
                   presList[[i]][[j]]$B,
                   col = coli,
                   lwd = 1.5)
            box(lwd = 1.5, col = colii)
            if(flagi) points(2022, 1.3, pch = 8, xpd = TRUE)
        }else{
            plot.new()
        }
        if(i == 1){
            mtext(mods[j], 3, 0.5, font = 2)
        }
        if(j == nmods){
            mtext(stocks[i], 4, 1, font = 2)
        }
    }
}
mtext("Year", 1, 3, outer = TRUE)
mtext("Residuals", 2, 3, outer = TRUE)
if(save.figs) dev.off()

pval.list <- vector("list", nstocks)
for(i in 1:nstocks){
    for(j in 1:nmods){
        pvals <- c()
        resid <- presList[[i]][[j]]$B
        if(!is.null(resid)){
            pvals <- c(
                Box.test(resid, lag=4, fitdf=1)$p.value ## ACF
                       )
        }else{
            pvals <- c(NA)
        }
        pval.list[[i]] <- rbind(pval.list[[i]], c(stocks[i], mods[j],
                                                  signif(pvals,2)))
    }
}

pval.tab <- do.call(rbind, pval.list)
pval.tab <- data.frame(pval.tab[,1:2],
                       pval = sapply(pval.tab[,-c(1,2)], function(x) as.numeric(x)))
pval.tab <- data.frame(pval.tab[,1:2],
                       pval = sapply(pval.tab[,-c(1,2)], function(x) if(is.na(x)){
                                        x
                                    }else if(x < 0.01){
                                        "<0.01"
                                    }else if(x < 0.05){
                                        "<0.05"
                                    }else x))
pval.tab <- pval.tab[!sapply(pval.tab[,-c(1,2)],function(x) all(is.na(x))),]

rownames(pval.tab) <- NULL
colnames(pval.tab) <- c("Stock","Model","BoxTest")

pval.tab[is.na(pval.tab)] <- ""

pval.tab


tmp <- pval.tab
colnames(tmp) <- c("Stock", "Model",
                   "Box-test\\textsubscript{B}")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
    paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

capi <- "P-values of the auto-correlation test (Ljung-Box test) for the biomass process residuals."

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = capi,
                                                     label = "tab:si-resid",
                                                     align="ll|lc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("table_resid_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{10pt}{10pt}\\selectfont"
                      )



nmods <- length(mods)
vars <- c("logrvec","logKvec","logmvec","logFmsyvec","logBmsyvec","logFnotS","logB")
df4 <- as.data.frame(matrix(NA, nrow = nstocks * length(mods), ncol = 2 + length(vars)))
for(i in 1:nstocks){
    for(j in 1:length(mods)){
        tmp <- c(stocks[i], mods[j])
        for(k in 1:length(vars)){
            if(j > 1 && !all(is.na(resList[[i]][[j]]$time))){
                ind <- resList[[i]][[j]]$time < 2022
            }else{
                ind <- 1
            }
            if(j == 1 && vars[k] == "logFmsyvec"){
                vari <- "logFmsy"
            }else if(j == 1 && vars[k] == "logBmsyvec"){
                vari <- "logBmsy"
            }else{
                vari <- vars[k]
            }
            tmp <- c(tmp, round(mean(resList[[i]][[j]][[vari]][ind,5], na.rm = TRUE) *
                                100, 0))
        }
    df4[(i-1)*length(mods)+j,] <- tmp
    }
}
rownames(df4) <- NULL
colnames(df4) <- c("Stock","Model","r","K","MSY","Fmsy","Bmsy","F","B")


df4[df4 == "NaN"] <- ""

for(i in 3:ncol(df4)){
    df4[as.numeric(df4[,i]) > 100 & !is.na(as.numeric(df4[,i])),i] <- ">100"
}

df4 <- df4[-which(apply(df4[,3:ncol(df4)],1, function(x) all(x == ""))),]

df4

tmp <- df4 ## [-1,]
colnames(tmp) <- c("Stock", "Model", "r", "K", "MSY", "F_\\textsubscript{MSY}","B_\\textsubscript{MSY}","F","B")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Average coefficient of variation (CV) in percent (\\%) of main model parameters and reference points.",
                                                     label = "tab:si-cv",
                                                     align="ll|lccccccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("SI_table_cv_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE
                      )


tmp <- df4

indi.above100 <- which(apply(tmp[,3:ncol(tmp)],1, function(x) any(x == ">100")))
tmp[indi.above100,]

sdigits <- 2
nmods <- length(mods)

vars <- c("logrvec","logKvec", "logsdb","logsdf")

df5 <- c("Stock","Model","av_r","av_K","sdb","sdf")
for(i in 1:nstocks){
    for(j in 1:length(mods)){
        tmp <- c(stocks[i], mods[j])
        for(k in 1:length(vars)){
            if(mods[j] == "CM" &&
               vars[k] == "logFmsyvec"){
                varii <- "logFmsy"
            }else if(mods[j] == "CM" &&
               vars[k] == "logBmsyvec"){
                varii <- "logBmsy"
            }else{
                varii <- vars[k]
            }
            est <- resList[[i]][[j]][[varii]][,2]
            lo <- resList[[i]][[j]][[varii]][,1]
            up <- resList[[i]][[j]][[varii]][,3]
            if(vars[k] %in% c("logmvec","logKvec","logBmsyvec")){
                cori <- 1e3
            }else if(vars[k] %in% c("logsdb","logsdf")){
                cori <- 1/1e2
            }else{
                cori <- 1
            }
            tmp <- c(tmp, paste0(
                              signif(mean(est)/cori, sdigits)," (",
                              signif(mean(lo)/cori, sdigits),"-",
                              signif(mean(up)/cori, sdigits),")"))
        }
    df5 <- rbind(df5, tmp)
    }
}
rownames(df5) <- NULL
colnames(df5) <- NULL


df5[df5 == "NA (NA-NA)"] <- NA

indili <- !apply(is.na(df5[,3:ncol(df5)]),1,all)
df5 <- cbind(df5[indili,1:2], df5[indili,3:ncol(df5)])

df5[is.na(df5)] <- ""

df5

tmp <- df5[-1,]
colnames(tmp) <- c("Stock", "Model", "$\\overline{r}$ [$yr\\textsuperscript{-1}$]", "$\\overline{K}$ ['000t]", "$\\sigma_B$ [$1e\\textsuperscript{-2}$]", "$\\sigma_F$ [$1e\\textsuperscript{-2}$]")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Constant and average time-varying parameters with lower and upper 95\\% confidence intervals.",
                                                     label = "tab:pars1",
                                                     align="ll|lcccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("table_modPars1_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{8pt}{8pt}\\selectfont"
                      )


sdigits <- 2
nmods <- length(mods)

vars <- c("logpsi","logsdm","logpsiK","logsdK","mkb")

df5a <- c("Stock","Model","psi_m", "sd_m","psi_K","sd_K","lambda")
for(i in 1:nstocks){
    for(j in 1:length(mods)){
        tmp <- c(stocks[i], mods[j])
        for(k in 1:length(vars)){
            if(mods[j] == "CM" &&
               vars[k] == "logFmsyvec"){
                varii <- "logFmsy"
            }else if(mods[j] == "CM" &&
               vars[k] == "logBmsyvec"){
                varii <- "logBmsy"
            }else{
                varii <- vars[k]
            }
            est <- resList[[i]][[j]][[varii]][,2]
            lo <- resList[[i]][[j]][[varii]][,1]
            up <- resList[[i]][[j]][[varii]][,3]
            if(vars[k] %in% c("logmvec","logKvec","logBmsyvec")){
                cori <- 1e3
            }else if(vars[k] %in% c("logsdb","logsdf")){
                cori <- 1/1e2
            }else{
                cori <- 1
            }
            tmp <- c(tmp, paste0(
                              signif(mean(est)/cori, sdigits)," (",
                              signif(mean(lo)/cori, sdigits),"-",
                              signif(mean(up)/cori, sdigits),")"))
        }
    df5a <- rbind(df5a, tmp)
    }
}
rownames(df5a) <- NULL
colnames(df5a) <- NULL


df5a <- df5a[-which(df5a[,2] == "CM"),]
df5a[df5a == "NA (NA-NA)"] <- NA
df5a[df5a[,2] %in% c("RM"),c(5,6)] <- NA
df5a[df5a[,2] %in% c("KM","PKRM","LKRM"),c(3,4)] <- NA
df5a[df5a[,2] %in% c("KM","RM","PKRM","IKRM"),c(7)] <- NA

df5a <- df5a[-which(apply(df5a[,3:ncol(df5a)],1,function(x) all(is.na(x)))),]

indili <- !apply(is.na(df5a[,3:ncol(df5a)]),1,all)
df5a <- cbind(df5a[indili,1:2], df5a[indili,3:ncol(df5a)])

df5a[is.na(df5a)] <- ""

df5a

tmp <- df5a[-1,]
colnames(tmp) <- c("Stock", "Model", "$\\psi_m$", "$\\sigma_m$", "$\\psi_K$", "$\\sigma_K$", "$\\lambda$")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Estimated parameters of the time-varying processes with lower and upper 95\\% confidence intervals.",
                                                     label = "tab:pars2",
                                                     align="ll|lccccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("table_modPars2_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{8pt}{8pt}\\selectfont"
                      )

sdigits <- 2
nmods <- length(mods)

vars <- c("logrvec","logKvec", "logsdb","logsdf",
          "logpsi","logsdm","logpsiK","logsdK","mkb")

df5b <- c("Stock","Model","av_r","av_K","sdb","sdf",
          "psi_m", "sd_m","psi_K","sd_K","lambda")
for(i in 1:nstocks){
    for(j in 1:length(mods)){
        tmp <- c(stocks[i], mods[j])
        for(k in 1:length(vars)){
            if(mods[j] == "CM" &&
               vars[k] == "logFmsyvec"){
                varii <- "logFmsy"
            }else if(mods[j] == "CM" &&
               vars[k] == "logBmsyvec"){
                varii <- "logBmsy"
            }else{
                varii <- vars[k]
            }
            est <- resList[[i]][[j]][[varii]][,2]
            lo <- resList[[i]][[j]][[varii]][,1]
            up <- resList[[i]][[j]][[varii]][,3]
            if(vars[k] %in% c("logmvec","logKvec","logBmsyvec")){
                cori <- 1e3
            }else if(vars[k] %in% c("logsdb","logsdf")){
                cori <- 1/1e2
            }else{
                cori <- 1
            }
            tmp <- c(tmp, paste0(
                              signif(mean(est)/cori, sdigits)," (",
                              signif(mean(lo)/cori, sdigits),"-",
                              signif(mean(up)/cori, sdigits),")"))
        }
    df5b <- rbind(df5b, tmp)
    }
}
rownames(df5b) <- NULL
colnames(df5b) <- NULL

df5b[df5b == "NA (NA-NA)"] <- NA
df5b[df5b[,2] %in% c("CM","RM"),c(5,6)+4] <- NA
df5b[df5b[,2] %in% c("CM","KM","PKRM","LKRM"),c(3,4)+4] <- NA
df5b[df5b[,2] %in% c("CM","KM","RM","PKRM","IKRM"),c(7)+4] <- NA

df5b <- df5b[-which(apply(df5b[,3:ncol(df5b)],1,function(x) all(is.na(x)))),]


indili <- !apply(is.na(df5b[,3:ncol(df5b)]),1,all)
df5b <- cbind(df5b[indili,1:2], df5b[indili,3:ncol(df5b)])

df5b[is.na(df5b)] <- ""

df5b

tmp <- df5b[-1,]
colnames(tmp) <- c("Stock", "Model", "$\\overline{r}$ [$yr\\textsuperscript{-1}$]", "$\\overline{K}$ ['000t]", "$\\sigma_B$ [$1e\\textsuperscript{-2}$]", "$\\sigma_F$ [$1e\\textsuperscript{-2}$]", "$\\psi_m$", "$\\sigma_m$", "$\\psi_K$", "$\\sigma_K$", "$\\lambda$")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Constant and average time-varying parameters with lower and upper 95\\% confidence intervals.",
                                                     label = "tab:si-pars",
                                                     align="ll|lccccccccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("table_modPars_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{6pt}{6pt}\\selectfont"
                      )

tmp <- as.data.frame(df5b[-1,])
colnames(tmp) <- df5b[1,]

tmpi <- sapply(stocks, function(x) as.numeric(sapply(tmp$av_r[tmp$Stock == x],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))
signif(sapply(tmpi,mean),3)
signif(range(unlist(tmpi), na.rm = TRUE),2)
signif(mean(unlist(tmpi), na.rm = TRUE),2)

lapply(tmpi, function(x) sd(x) / mean(x))

tmpi <- sapply(stocks, function(x) as.numeric(sapply(tmp$av_K[tmp$Stock == x],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))
signif(sapply(tmpi,mean),3)
mean(unlist(tmpi))
range(unlist(tmpi))

lapply(tmpi, function(x) sd(x) / mean(x))


tmpi <- sapply(stocks, function(x) as.numeric(sapply(tmp$av_K[tmp$Stock == x &
                                                              tmp$Model == "CM"],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))
tmpi2 <- sapply(stocks, function(x) as.numeric(sapply(tmp$av_K[tmp$Stock == x &
                                                               !tmp$Model %in% c("CM",
                                                                                 "KM")],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))
tmpi3 <- lapply(1:nstocks, function(x) (tmpi[[x]] - tmpi2[[x]]) / tmpi2[[x]]*100)
lapply(tmpi3,mean,na.rm = TRUE)
tmpi

tmpi <- sapply(stocks, function(x) as.numeric(sapply(tmp$sdb[tmp$Stock == x],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))
signif(mean(unlist(tmpi)*0.01),2)
range(unlist(tmpi)*0.01)

lapply(tmpi, function(x) sd(x) / mean(x))

lapply(tmpi,function(x) signif(mean(x * 0.01),1))


tmpi <- sapply(stocks, function(x) as.numeric(sapply(tmp$sdf[tmp$Stock == x],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))
signif(mean(unlist(tmpi)*0.01),2)
range(unlist(tmpi))

lapply(tmpi, function(x) sd(x) / mean(x))

lapply(tmpi,function(x) signif(mean(x * 0.01),2))

tmpi1 <- sapply(stocks, function(x) as.numeric(sapply(tmp$psi_m[tmp$Stock == x],
                                                      function(x) strsplit(x, " \\(")[[1]][1])))
tmpi2 <- sapply(stocks, function(x) as.numeric(sapply(tmp$sd_m[tmp$Stock == x],
                                                      function(x) strsplit(x, " \\(")[[1]][1])))
tmpi3 <- sapply(stocks, function(x) as.numeric(sapply(tmp$psi_K[tmp$Stock == x],
                                                      function(x) strsplit(x, " \\(")[[1]][1])))
tmpi4 <- sapply(stocks, function(x) as.numeric(sapply(tmp$sd_K[tmp$Stock == x],
                                                     function(x) strsplit(x, " \\(")[[1]][1])))

tmpi11 <- sapply(tmpi1,mean,na.rm = TRUE)
signif((mean(tmpi11[1:4],na.rm = TRUE)-mean(tmpi11[5:7], na.rm = TRUE)) /
       mean(tmpi11[5:7],na.rm = TRUE)*100,2)
tmpi33 <- sapply(tmpi3,mean,na.rm = TRUE)
signif((mean(tmpi33[1:4],na.rm=TRUE)-mean(tmpi33[5:7],na.rm=TRUE)) / mean(tmpi33[5:7],na.rm=TRUE)*100,2)

tmpi22 <- sapply(tmpi2,mean,na.rm = TRUE)
signif((mean(tmpi22[1:4],na.rm = TRUE)-mean(tmpi22[5:7], na.rm = TRUE)) /
       mean(tmpi22[5:7],na.rm = TRUE)*100,2)
tmpi44 <- sapply(tmpi4,mean,na.rm = TRUE)
signif((mean(tmpi44[1:4],na.rm=TRUE)-mean(tmpi44[5:7],na.rm=TRUE)) / mean(tmpi44[5:7],na.rm=TRUE)*100,2)

tmpi2 <- sapply(tmpi2,mean,na.rm = TRUE)
signif(tmpi2,2)
signif((mean(tmpi2[1:4])-mean(tmpi2[5:7])) / mean(tmpi2[5:7])*100,2)
tmpi4 <- sapply(tmpi4,mean,na.rm = TRUE)
signif(tmpi4,2)
signif((mean(tmpi4[1:4],na.rm=TRUE)-mean(tmpi4[5:7],na.rm=TRUE)) / mean(tmpi4[5:7],na.rm=TRUE)*100,2)

signif(sapply(tmpi2,mean,na.rm = TRUE),2)
signif(sapply(tmpi4,mean,na.rm = TRUE),2)


head(tmp)
signif(mean(unlist(tmpi)*0.01),2)
range(unlist(tmpi))

## lambda
tmpi <- sapply(stocks, function(x) as.numeric(sapply(tmp$lambda[tmp$Stock == x],
                                                      function(x) strsplit(x, " \\(")[[1]][1])))

range(tmpi, na.rm = TRUE)

ri <- resList[[3]][[5]]$logrvec[,2]
ri <- ri/min(ri)
ki <- resList[[3]][[5]]$logKvec[,2]
ki <- ki/min(ki)

res <- 120
width <- 1000
height <- 1200


tmp0 <- sapply(resList, function(x) sapply(x, function(y) round(y$logn[2],1)))
tmp <- cbind(rownames(tmp0), tmp0)
colnames(tmp) <- c("Model",colnames(tmp0))

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_all_",label,".png")),
        width = 1200, height = 1600, res = res)
plot.time.all2(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = NULL,
               remove.scenario = TRUE)
if(save.figs) dev.off()



if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("si_ts_all_",label,".png")),
        width = 1200, height = 1600, res = res)
plot.time.all2(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = NULL,
               remove.scenario = 2)
if(save.figs) dev.off()

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_r_",label,".png")),
        width = 1600, height = 1400, res = res)
plot.time.all3(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = 6,
               remove.scenario = FALSE, quant = 1,
               nonConv2 = nonConv2)
if(save.figs) dev.off()

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_k_",label,".png")),
        width = 1650, height = 1400, res = res)
plot.time.all3(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = 6,
               remove.scenario = 3, quant = 2,
               nonConv2 = nonConv2)
if(save.figs) dev.off()

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_m_",label,".png")),
        width = 1650, height = 1400, res = res)
plot.time.all3(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = 6,
               remove.scenario = 3, quant = 3,
               nonConv2 = nonConv2)
if(save.figs) dev.off()


if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_b_",label,".png")),
        width = 1625, height = 1400, res = res)
plot.time.all4(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = 6,
               remove.scenario = 3, quant = 4,
               nonConv2 = nonConv2)
if(save.figs) dev.off()


if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_f_",label,".png")),
        width = 1625, height = 1400, res = res)
plot.time.all4(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = 6,
               remove.scenario = 3, quant = 5,
               nonConv2 = nonConv2)
if(save.figs) dev.off()

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("ts_c_",label,".png")),
        width = 1625, height = 1400, res = res)
plot.time.all4(resList, stockData, stockInfo, selMods, 2021, nPar = nPar,
               bestMods = bestMods, drop.mods = 6,
               remove.scenario = 3, quant = 6,
               nonConv2 = nonConv2)
if(save.figs) dev.off()



resList.sel <- lapply(1:nstocks, function(x) resList[[x]][names(resList[[x]]) == bestMods[[x]]][[1]])

tmpi <- sapply(1:nstocks, function(x) round((mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 2021 & resList.sel[[x]]$time < 2022)]) -
mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 1981 & resList.sel[[x]]$time < 1982)])) / mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 1981 & resList.sel[[x]]$time < 1982)]) * 100))

round(mean(tmpi))

## roundfish
tmpi <- sapply(c(1:4), function(x) round((mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 2021 & resList.sel[[x]]$time < 2022)]) -
mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 1981 & resList.sel[[x]]$time < 1982)])) / mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 1981 & resList.sel[[x]]$time < 1982)]) * 100))

round(mean(tmpi))

## flatfish
tmpi <- sapply(c(5:7), function(x) round((mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 2021 & resList.sel[[x]]$time < 2022)]) -
mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 1981 & resList.sel[[x]]$time < 1982)])) / mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= 1981 & resList.sel[[x]]$time < 1982)]) * 100))

round(mean(tmpi))




sdigits <- 2
nmods <- length(mods)

vars <- c("logFmsyvec","logBmsyvec","logmvec",
          "logFFmsynotS","logBBmsy")

df5c <- c("Stock","Model","aic","av_fmsy","av_bmsy","av_m",
          "av_ffmsy","av_bbmsy")  ## average in 2021
for(i in 1:nstocks){
    for(j in 1:length(mods)){
        tmp <- c(stocks[i], mods[j])
        aic <- signif(df[which(df$stock == stocks[i] &
                                      df$model == mods[j]),"AICc"], 5)
        tmp <- c(tmp, aic)
        for(k in 1:length(vars)){

            if(!all(is.na(resList[[i]][[j]]$time))){
                indi <- which(resList[[i]][[j]]$time == final.year)
            }else indi <- 1
            if(mods[j] == "CM" &&
               vars[k] == "logFmsyvec"){
                varii <- "logFmsy"
                indi <- 1
            }else if(mods[j] == "CM" &&
               vars[k] == "logBmsyvec"){
                varii <- "logBmsy"
                indi <- 1
            }else{
                varii <- vars[k]
            }
            if(all(is.na(indi))){
                indi <- 1
            }
            est <- resList[[i]][[j]][[varii]][indi,2]
            lo <- resList[[i]][[j]][[varii]][indi,1]
            up <- resList[[i]][[j]][[varii]][indi,3]
            if(vars[k] %in% c("logmvec","logKvec","logBmsyvec")){
                cori <- 1e3
            }else if(vars[k] %in% c("logsdb","logsdf")){
                cori <- 1/1e2
            }else{
                cori <- 1
            }
            tmp <- c(tmp, paste0(
                              signif(mean(est)/cori, sdigits)," (",
                              signif(mean(lo)/cori, sdigits),"-",
                              signif(mean(up)/cori, sdigits),")"))
        }
    df5c <- rbind(df5c, tmp)
    }
    tmp <- c(stocks[i], "EqSim","")
    tmp <- c(tmp, signif(as.numeric(stockInfo[stockInfo$stock %in%
                                       stocks.true[i],
                                       c("EQ_FmsyESB","EQ_ESBmsy","EQ_MSY")]) /
             c(1,1e3,1e3),sdigits))
    tmp <- c(tmp, signif(as.numeric(stockInfo[stockInfo$stock %in%
                                       stocks.true[i],
                                       c("EQ_FFmsyESB","EQ_ESBESBmsy")]),sdigits))
    df5c <- rbind(df5c, tmp)
}
rownames(df5c) <- NULL
colnames(df5c) <- NULL


df5c[df5c == "NA (NA-NA)"] <- NA
df5c[df5c == "NA"] <- NA

indili <- !apply(is.na(df5c[,3:ncol(df5c)]),1,all)
df5c <- cbind(df5c[indili,1:2], df5c[indili,3:ncol(df5c)])

df5c[is.na(df5c)] <- ""

df5c


tmp <- df5c[-1,]

for(i in 1:length(stocks)){
    tmpi <- as.numeric(tmp[tmp[,1] == stocks[i],3])
    indi <- which.min(tmpi)
    tmp[tmp[,1] == stocks[i],3] <- signif(tmpi - tmpi[indi],3)
}

tmp

colnames(tmp) <- c("Stock", "Model", "$\\Delta$AIC\\textsubscript{c}", "$F$\\textsubscript{MSY,$t$} [$yr\\textsuperscript{-1}$]", "$B$\\textsubscript{MSY,$t$} ['000 t]", "MSY\\textsubscript{$t$} ['000 t]", "$F_t/F$\\textsubscript{MSY,$t$}", "$B_t/B$\\textsubscript{MSY,$t$}")

indi <- which(paste0(stocks.true[match(tmp[,1], stocks)],"_",tmp[,2]) %in% nonConv2)
if(length(indi) > 0){
    tmp[,2][indi] <- paste0(tmp[,2][indi],"\\textsuperscript{†}")
}

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
    paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")


tmp[tmp[,3] == "0" & !is.na(tmp[,3]),2:ncol(tmp)] <- sapply(tmp[tmp[,3] == "0" & !is.na(tmp[,3]),2:ncol(tmp)], function(x) paste0("\\textbf{",x,"}"))


if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Model, $\\Delta$AIC\\textsubscript{c}, and terminal year (2021) reference points and stock status ratios, with lower and upper 95\\% confidence intervals in parentheses. B\\textsubscript{MSY,$t$} and MSY$_t$ are given in thousand tonnes. The best model ($\\Delta$AIC\\textsubscript{c}=0) and associated values are highlighted in bold font. Models with the † symbol showed large uncertainty in some model parameters and/or derived quantites and need to be interpreted with caution.",
                                                     label = "tab:refs",
                                                     align="ll|lcccccc"),
                                      file = file.path(resdir, "tabs", label,
                                                       paste0("SI_table_refs_CIs_",label,".tex")),
                                      include.rownames = FALSE, include.colnames = TRUE,
                                      sanitize.text.function = identity,
                                      hline.after = c(-1,nrow(tmp)),
                                      caption.placement = "top",
                                      floating=TRUE,
                                      booktabs=TRUE,
                                      size="\\fontsize{7pt}{7pt}\\selectfont"
                                      )

sdigits <- 2
nmods <- length(mods)

vars <- c("logFmsyvec","logBmsyvec","logmvec",
          "logFFmsynotS","logBBmsy")

df5c <- c("Stock","Model","aic","av_fmsy","av_bmsy","av_m",
          "av_ffmsy","av_bbmsy")  ## average in 2021
for(i in 1:nstocks){
    for(j in 1:length(mods)){
        tmp <- c(stocks[i], mods[j])
        aic <- signif(df[which(df$stock == stocks[i] &
                                      df$model == mods[j]),"AICc"], 5)
        tmp <- c(tmp, aic)
        for(k in 1:length(vars)){
            if(!all(is.na(resList[[i]][[j]]$time))){
                indi <- which(resList[[i]][[j]]$time == final.year)
            }else indi <- 1
            if(mods[j] == "CM" &&
               vars[k] == "logFmsyvec"){
                varii <- "logFmsy"
                indi <- 1
            }else if(mods[j] == "CM" &&
               vars[k] == "logBmsyvec"){
                varii <- "logBmsy"
                indi <- 1
            }else{
                varii <- vars[k]
            }
            if(all(is.na(indi))){
                indi <- 1
            }
            est <- resList[[i]][[j]][[varii]][indi,2]
            lo <- resList[[i]][[j]][[varii]][indi,1]
            up <- resList[[i]][[j]][[varii]][indi,3]
            if(vars[k] %in% c("logmvec","logKvec","logBmsyvec")){
                cori <- 1e3
            }else if(vars[k] %in% c("logsdb","logsdf")){
                cori <- 1/1e2
            }else{
                cori <- 1
            }
            tmp <- c(tmp, paste0(
                              signif(mean(est)/cori, sdigits)))
        }
    df5c <- rbind(df5c, tmp)
    }
    tmp <- c(stocks[i], "EqSim","")
    tmp <- c(tmp, signif(as.numeric(stockInfo[stockInfo$stock %in%
                                       stocks.true[i],
                                       c("EQ_FmsyESB","EQ_ESBmsy","EQ_MSY")]) /
             c(1,1e3,1e3),sdigits))
    tmp <- c(tmp, signif(as.numeric(stockInfo[stockInfo$stock %in%
                                       stocks.true[i],
                                       c("EQ_FFmsyESB","EQ_ESBESBmsy")]),sdigits))
    df5c <- rbind(df5c, tmp)
}
rownames(df5c) <- NULL
colnames(df5c) <- NULL


df5c[df5c == "NA (NA-NA)"] <- NA
df5c[df5c == "NA"] <- NA

indili <- !apply(is.na(df5c[,3:ncol(df5c)]),1,all)
df5c <- cbind(df5c[indili,1:2], df5c[indili,3:ncol(df5c)])

df5c[is.na(df5c)] <- ""

df5c

tmp <- df5c[-1,]


for(i in 1:length(stocks)){
    tmpi <- as.numeric(tmp[tmp[,1] == stocks[i],3])
    indi <- which.min(tmpi)
    tmp[tmp[,1] == stocks[i],3] <- signif(tmpi - tmpi[indi],3)
}

tmp

colnames(tmp) <- c("Stock", "Model", "$\\Delta$AIC\\textsubscript{c}", "$F$\\textsubscript{MSY,$t$} [$yr\\textsuperscript{-1}$]", "$B$\\textsubscript{MSY,$t$} ['000 t]", "MSY\\textsubscript{$t$} ['000 t]", "$F_t/F$\\textsubscript{MSY,$t$}", "$B_t/B$\\textsubscript{MSY,$t$}")

indi <- which(paste0(stocks.true[match(tmp[,1], stocks)],"_",tmp[,2]) %in% nonConv2)
if(length(indi) > 0){
    tmp[,2][indi] <- paste0(tmp[,2][indi],"\\textsuperscript{†}")
}

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
    paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")


tmp[tmp[,3] == "0" & !is.na(tmp[,3]),2:ncol(tmp)] <- sapply(tmp[tmp[,3] == "0" & !is.na(tmp[,3]),2:ncol(tmp)], function(x) paste0("\\textbf{",x,"}"))


if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Model, $\\Delta$AIC\\textsubscript{c}, and terminal year (2021) reference points and stock status ratios. B\\textsubscript{MSY,$t$} and MSY$_t$ are given in thousand tonnes. The best model ($\\Delta$AIC\\textsubscript{c}=0) and associated values are highlighted in bold font. Models with the † symbol showed large uncertainty in some model parameters and/or derived quantites and need to be interpreted with caution.",
                                                     label = "tab:refs",
                                                     align="ll|lcccccc"),
                                      file = file.path(resdir, "tabs", label,
                                                       paste0("table_refs_",label,".tex")),
                                      include.rownames = FALSE, include.colnames = TRUE,
                                      sanitize.text.function = identity,
                                      hline.after = c(-1,nrow(tmp)),
                                      caption.placement = "top",
                                      floating=TRUE,
                                      booktabs=TRUE,
                                      size="\\fontsize{7pt}{7pt}\\selectfont"
                                      )



tmp <- as.data.frame(df5c[-1,])
colnames(tmp) <- df5c[1,]

tmp2 <- tmp[!tmp$Model %in% c("EqSim"),]
tmp3 <- tmp[tmp$Model %in% c("EqSim"),]


resList.sel <- lapply(1:nstocks, function(x) resList[[x]][names(resList[[x]]) == bestMods[[x]]][[1]])

histi <- c(1981,1981+5)
reci <- c(2022-5,2022)

tmpi1 <- sapply(1:nstocks, function(x) round((mean(resList.sel[[x]]$logFmsyvec[which(resList.sel[[x]]$time >= reci[1] & resList.sel[[x]]$time < reci[2])]) -
                                              mean(resList.sel[[x]]$logFmsyvec[which(resList.sel[[x]]$time >= histi[1] & resList.sel[[x]]$time < histi[2])])) / mean(resList.sel[[x]]$logFmsyvec[which(resList.sel[[x]]$time >= histi[1] & resList.sel[[x]]$time < histi[2])]) * 100))

tmpi2 <- sapply(1:nstocks, function(x) round((mean(resList.sel[[x]]$logBmsyvec[which(resList.sel[[x]]$time >= reci[1] & resList.sel[[x]]$time < reci[2])]) -
                                              mean(resList.sel[[x]]$logBmsyvec[which(resList.sel[[x]]$time >= histi[1] & resList.sel[[x]]$time < histi[2])])) / mean(resList.sel[[x]]$logBmsyvec[which(resList.sel[[x]]$time >= histi[1] & resList.sel[[x]]$time < histi[2])]) * 100))

tmpi3 <- sapply(1:nstocks, function(x) round((mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= reci[1] & resList.sel[[x]]$time < reci[2])]) -
                                             mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= histi[1] & resList.sel[[x]]$time < histi[2])])) / mean(resList.sel[[x]]$logMSYvec[which(resList.sel[[x]]$time >= histi[1] & resList.sel[[x]]$time < histi[2])]) * 100))

data.frame(stock = stocks,
           fmsy = tmpi1,
           bmsy = tmpi2,
           msy = tmpi3)

(mean(resList.sel[[5]]$logFmsyvec[,2]) - as.numeric(tmp3[5,4])) / as.numeric(tmp3[5,4])



## Fmsy
apply(tmp3, 2, range)
range(as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)))

tmpii <- cbind(tmp2$Stock, as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)))
tmpii[tmpii[,2] > 0.5,]

## Bmsy
apply(tmp3, 2, range)
range(as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)))

as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1))

range(as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)))
range(as.numeric(sapply(strsplit(tmp2$av_m," \\("),"[[",1)))


as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1))

by(as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE))
mean(by(as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE)),na.rm = TRUE) * 100

by(as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE))
mean(by(as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE)),na.rm = TRUE) * 100

by(as.numeric(sapply(strsplit(tmp2$av_m," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE))
mean(by(as.numeric(sapply(strsplit(tmp2$av_m," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE)),na.rm = TRUE) * 100


tmp2 <- tmp[!tmp$Model %in% c("CM","EqSim"),]

by(as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE))

by(as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE))

by(as.numeric(sapply(strsplit(tmp2$av_m," \\("),"[[",1)),
   tmp2$Stock, function(x) sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE))


tmp2 <- tmp[!tmp$Model %in% c("EqSim"),]
tmp3 <- tmp[tmp$Model == "EqSim",]

## Fmsy
tmpii <- data.frame(stock = tmp2$Stock, fmsy = as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, fmsy2 = as.numeric(tmp3$av_fmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
by(round((tmpiii$fmsy - tmpiii$fmsy2) / tmpiii$fmsy2 * 100), tmpiii$stock, mean)
signif(mean(by((tmpiii$fmsy - tmpiii$fmsy2) / tmpiii$fmsy2 * 100, tmpiii$stock, mean)),2)

## Bmsy
tmpii <- data.frame(stock = tmp2$Stock, bmsy = as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, bmsy2 = as.numeric(tmp3$av_bmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
by(round((tmpiii$bmsy - tmpiii$bmsy2) / tmpiii$bmsy2 * 100), tmpiii$stock, mean)
signif(mean(by((tmpiii$bmsy - tmpiii$bmsy2) / tmpiii$bmsy2 * 100, tmpiii$stock, mean)),2)

## MSY
tmpii <- data.frame(stock = tmp2$Stock, msy = as.numeric(sapply(strsplit(tmp2$av_m," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, msy2 = as.numeric(tmp3$av_m))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
by(round((tmpiii$msy - tmpiii$msy2) / tmpiii$msy2 * 100), tmpiii$stock, mean)
signif(mean(by((tmpiii$msy - tmpiii$msy2) / tmpiii$msy2 * 100, tmpiii$stock, mean)),2)


## Differences between EqSim and KBPMs per model
tmpii <- data.frame(model = tmp2$Model,
                    stock = tmp2$Stock,
                    fmsy = as.numeric(sapply(strsplit(tmp2$av_fmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, fmsy2 = as.numeric(tmp3$av_fmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
tmpiii$re <- (tmpiii$fmsy - tmpiii$fmsy2) / tmpiii$fmsy2 * 100
aggregate(tmpiii$re, by=list(model = tmpiii$model), mean)


tmpii <- data.frame(model = tmp2$Model,
                    stock = tmp2$Stock,
                    bmsy = as.numeric(sapply(strsplit(tmp2$av_bmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, bmsy2 = as.numeric(tmp3$av_bmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
tmpiii$re <- (tmpiii$bmsy - tmpiii$bmsy2) / tmpiii$bmsy2 * 100
aggregate(tmpiii$re, by=list(model = tmpiii$model), mean)


tmpii <- data.frame(model = tmp2$Model,
                    stock = tmp2$Stock,
                    msy = as.numeric(sapply(strsplit(tmp2$av_m," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, msy2 = as.numeric(tmp3$av_m))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
tmpiii$re <- (tmpiii$msy - tmpiii$msy2) / tmpiii$msy2 * 100
aggregate(tmpiii$re, by=list(model = tmpiii$model), mean)



tmp2 <- tmp[!tmp$Model %in% c("EqSim"),]

range(as.numeric(sapply(strsplit(tmp2$av_ffmsy," \\("),"[[",1)))
range(as.numeric(sapply(strsplit(tmp2$av_bbmsy," \\("),"[[",1)))

by(as.numeric(sapply(strsplit(tmp2$av_ffmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) signif(mean(x,na.rm = TRUE),2))

by(as.numeric(sapply(strsplit(tmp2$av_bbmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) signif(mean(x,na.rm = TRUE),2))

by(as.numeric(sapply(strsplit(tmp2$av_ffmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) signif(sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE)*100,2))

by(as.numeric(sapply(strsplit(tmp2$av_bbmsy," \\("),"[[",1)),
   tmp2$Stock, function(x) signif(sd(x, na.rm = TRUE)/mean(x,na.rm = TRUE)*100,2))


tmp2 <- tmp[!tmp$Model %in% c("CM","EqSim"),]
tmp3 <- tmp[tmp$Model == "EqSim",]

tmp2

## F/Fmsy
tmpii <- data.frame(stock = tmp2$Stock,
                    fmsy = as.numeric(sapply(strsplit(tmp2$av_ffmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock,
                     fmsy2 = as.numeric(tmp3$av_ffmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
by(round((tmpiii$fmsy2 - tmpiii$fmsy) / tmpiii$fmsy * 100), tmpiii$stock, mean)
mean(by((tmpiii$fmsy2 - tmpiii$fmsy) / tmpiii$fmsy * 100, tmpiii$stock, mean))

cbind(tmp2$Stock,tmp2$Model, round((tmpiii$fmsy2 - tmpiii$fmsy) / tmpiii$fmsy * 100))


## B/Bmsy
tmpii <- data.frame(stock = tmp2$Stock,
                    fmsy = as.numeric(sapply(strsplit(tmp2$av_bbmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock,
                     fmsy2 = as.numeric(tmp3$av_bbmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")

## Quiet variable between mdoels!
cbind(tmp2$Stock,tmp2$Model, round((tmpiii$fmsy2 - tmpiii$fmsy) / tmpiii$fmsy * 100))

by(round((tmpiii$fmsy2 - tmpiii$fmsy) / tmpiii$fmsy * 100), tmpiii$stock, mean)
mean(by((tmpiii$fmsy2 - tmpiii$fmsy) / tmpiii$fmsy * 100, tmpiii$stock, mean))


tmp2 <- tmp[!tmp$Model %in% c("EqSim"),]
tmp3 <- tmp[tmp$Model == "EqSim",]

tmpii <- data.frame(model = tmp2$Model,
                    stock = tmp2$Stock,
                    msy = as.numeric(sapply(strsplit(tmp2$av_ffmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, msy2 = as.numeric(tmp3$av_ffmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
tmpiii$re <- (tmpiii$msy - tmpiii$msy2) / tmpiii$msy2 * 100
aggregate(tmpiii$re, by=list(model = tmpiii$model), function(x) signif(mean(x),2))


tmpii <- data.frame(model = tmp2$Model,
                    stock = tmp2$Stock,
                    msy = as.numeric(sapply(strsplit(tmp2$av_bbmsy," \\("),"[[",1)))
tmpiii <- data.frame(stock = tmp3$Stock, msy2 = as.numeric(tmp3$av_bbmsy))
tmpiii <- plyr::join(tmpii, tmpiii, by = "stock")
tmpiii$re <- (tmpiii$msy - tmpiii$msy2) / tmpiii$msy2 * 100
aggregate(tmpiii$re, by=list(model = tmpiii$model), function(x) signif(mean(x),2))


width <- 1000
height <- 1300

corsEnv0 <- get.cors(resList, stockData, stockInfo,
                     bestMods,
                      envData = envData,
                      maxYear = 2021,
                      var = c("m","temp.bot","sal.bot"),
                      fit.poly = TRUE,
                     lag = 0)

sapply(resList, function(x) min(x[[4]]$time))

indi1 <- which(envData$year >= 1981 &
               envData$year < 1991)
indi2 <- which(envData$year >= 2011 &
               envData$year < 2021)

signif(mean(envData$temp.bot[indi1]),2)
signif(mean(envData$temp.bot[indi2]),2)

signif(mean(envData$sal.bot[min(indi1):max(indi2)]),3)

sapply(corsEnv0, function(x){
    tmp <- predict(x$mod, newdata = data.frame(expl1 = c(
                                                   signif(mean(envData$temp.bot[indi1]),2),
                                                   signif(mean(envData$temp.bot[indi2]),2)),
                                               expl2 = signif(mean(envData$sal.bot[min(indi1):max(indi2)]),3)), type = "response")
    signif((tmp[2] - tmp[1]) / tmp[1] * 100,2)
})

signif(mean(envData$sal.bot[indi1]),3)
signif(mean(envData$sal.bot[indi2]),3)
signif(mean(envData$temp.bot[min(indi1):max(indi2)]),2)

sapply(corsEnv0, function(x){
    tmp <- predict(x$mod, newdata = data.frame(expl2 = c(
                                                   signif(signif(mean(envData$sal.bot[indi1]),3),3),
                                                   signif(signif(mean(envData$sal.bot[indi2]),3))),
                                               expl1 = signif(mean(envData$temp.bot[min(indi1):max(indi2)]),2)), type = "response")
    signif((tmp[2] - tmp[1]) / tmp[1] * 100,2)
})


## For text: model terms
lapply(corsEnv0, function(x) x$mod)


## For plot + caption: min and max observered salinity
signif(quantile(envData$sal.bot[envData$year >= 1981 &
                envData$year <= 2021],c(0.1,0.9)),3)


if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("cor_temp_sal_lag0_",label,".png")),
        width = width, height = height, res = res)
plot.cor2(resList, stockData, stockInfo, selMods = bestMods,
          corsEnv0,
          leg.ncol = 2,
          mfrow = c(4,2), byrow = FALSE,
          legend.extra = TRUE)
if(save.figs) dev.off()


sdigits <- 2
nlags <- 11
cori <- vector("list", nlags)
for(i in 1:nlags){
    cori[[i]] <- get.cors(resList, stockData, stockInfo,
                     bestMods,
                     envData = envData,
                     maxYear = 2021,
                     var = c("m","temp.bot","sal.bot"),
                     fit.poly = TRUE,
                     lag = i-1)
}
df6 <- c("Stock","Lag","Intercept","Temp","Temp2","Sal","R2","AICc")
daic.all <- NULL
for(i in 1:nstocks){
    daic <- NULL
    for(j in 1:nlags){
        tmp <- c(stocks[i], j-1)
        tmpi <- signif(coefficients(cori[[j]][[i]]$mod),sdigits)
        alli <- rep(NA, 4)
        alli[match(names(tmpi),
                   c("(Intercept)",
                     "poly(expl1, degree = degree)1",
                     "poly(expl1, degree = degree)2",
                     "expl2"))] <- tmpi
        alli <- c(alli,signif(with(summary(cori[[j]][[i]]$mod),
                                   1 - deviance/null.deviance),sdigits),
                  round(AIC(cori[[j]][[i]]$mod),1))
        daic <- c(daic,AIC(cori[[j]][[i]]$mod))
        tmp <- c(tmp, alli)
        df6 <- rbind(df6, tmp)
    }
    daic.all <- c(daic.all,daic - min(daic))
}
rownames(df6) <- NULL
colnames(df6) <- NULL

df6[is.na(df6)] <- ""


tmp <- df6[-1,]
colnames(tmp) <- c("Stock", "Lag", "Intercept",
                   "Temp", "Temp\\textsuperscript{2}",
                   "Salinity", "R\\textsuperscript{2}","AIC")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")


tmp[,2] <- paste0(tmp[,2], ifelse(daic.all  == 0, "\\textsuperscript{*}", ""))


if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Parameters, McFadden's R\\textsuperscript{2}, and AIC for the best model for each stock and time lags between 0 and 8 years. Asterisks indicate the best model for each stock according to AIC.",
                                                     label = "tab:si-cor-lag",
                                                     align="ll|ccccccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("SI_table_cor_lags_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)), ## cumsum(rle(tmp[,1])$lengths),
                      caption.placement = "top", #"top", NULL
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{8pt}{8pt}\\selectfont"
                      )



corsEnv.best <- get.cors(resList, stockData, stockInfo,
                     bestMods,
                      envData = envData,
                      maxYear = 2021,
                      var = c("m","temp.bot","sal.bot"),
                      fit.poly = TRUE,
                     lag = as.numeric(sapply(strsplit(tmp[,2][daic.all == 0],"*"),"[[",1)))


if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("cor_temp_sal_lagBest_",label,".png")),
        width = width, height = height, res = res)
plot.cor2(resList, stockData, stockInfo, selMods = bestMods,
          corsEnv.best,
          leg.ncol = 2,
          mfrow = c(4,2), byrow = FALSE,
          legend.extra = TRUE)
if(save.figs) dev.off()


sdigits <- 2
nlags <- 6
cori <- vector("list", nlags)
for(i in 1:nlags){
    besti <- as.list(rep(c("CM","RM","KM","PKRM","LKRM","IKRM")[i],
                         length(resList)))
    if(besti[[1]] == "CM"){
        mod.link <- "identity"
        pres <- presList
    }else{
        mod.link <- "log"
        pres <- NULL
    }
    cori[[i]] <- get.cors(resList, stockData, stockInfo,
                     besti,
                     envData = envData,
                     maxYear = 2021,
                     var = c("m","temp.bot","sal.bot"),
                     fit.poly = TRUE,
                     mod.link = mod.link,
                     lag = 0, pres = pres)
}

df7 <- c("Stock","Model","Intercept","Temp","Temp2","Sal","R2","AIC")
daic.all <- NULL
for(i in 1:nstocks){
    daic <- NULL
    for(j in 1:nlags){
        tmp <- c(stocks[i], c("CM","RM","KM","PKRM","LKRM","IKRM")[j])
        alli <- rep(NA, 4)
        if(!is.null(cori[[j]][[i]]$mod)){
            tmpi <- signif(coefficients(cori[[j]][[i]]$mod),sdigits)
            alli[match(names(tmpi),
                       c("(Intercept)",
                         "poly(expl1, degree = degree)1",
                         "poly(expl1, degree = degree)2",
                         "expl2"))] <- tmpi
            alli <- c(alli,signif(with(summary(cori[[j]][[i]]$mod),
                                       1 - deviance/null.deviance),sdigits),
                  round(AIC(cori[[j]][[i]]$mod),1))
        daic <- c(daic,AIC(cori[[j]][[i]]$mod))
        }else{
            alli <- c(alli, NA,NA)
        }
        tmp <- c(tmp, alli)
        df7 <- rbind(df7, tmp)
    }
    daic.all <- c(daic.all,daic - min(daic))
}
rownames(df7) <- NULL
colnames(df7) <- NULL

df7 <- df7[-which(apply(df7[,3:ncol(df7)],1,function(x) all(is.na(x)))),]

df7[is.na(df7)] <- ""

## remove AIC
df7 <- df7[,-8]

tmp <- df7[-1,]
colnames(tmp) <- c("Stock", "Model", "Intercept",
                   "Temp", "Temp\\textsuperscript{2}",
                   "Salinity", "R\\textsuperscript{2}")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
    paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")

if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Parameters and McFadden's R\\textsuperscript{2} for all converged models and a lag of 0.",
                                                     label = "tab:si-cor-mod",
                                                     align="ll|cccccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("SI_table_cor_mods_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{8pt}{8pt}\\selectfont"
                      )


here <- 6

sdigits <- 2
nlags <- 9
cori <- vector("list", nlags)
besti <- as.list(rep(c("CM","RM","KM","PKRM","LKRM","IKRM")[here],
                     length(resList)))
for(i in 1:nlags){
    if(besti[[1]] == "CM"){
        mod.link <- "identity"
        pres <- presList
    }else{
        mod.link <- "log"
        pres <- NULL
    }
    cori[[i]] <- get.cors(resList, stockData, stockInfo,
                     besti,
                     envData = envData,
                     maxYear = 2021,
                     var = c("m","temp.bot","sal.bot"),
                     fit.poly = TRUE,
                     mod.link = mod.link,
                     lag = i-1,
                     pres = pres)
}
df6 <- c("Stock","Lag","Intercept","Temp","Temp2","Sal","R2","AIC")
daic.all <- NULL
for(i in 1:nstocks){
    daic <- NULL
    for(j in 1:nlags){
        tmp <- c(stocks[i], j-1)
        alli <- rep(NA, 4)
        if(!is.null(coefficients(cori[[j]][[i]]$mod))){
        tmpi <- signif(coefficients(cori[[j]][[i]]$mod),sdigits)
        alli[match(names(tmpi),
                   c("(Intercept)",
                     "poly(expl1, degree = degree)1",
                     "poly(expl1, degree = degree)2",
                     "expl2"))] <- tmpi
        alli <- c(alli,signif(with(summary(cori[[j]][[i]]$mod),
                                   1 - deviance/null.deviance),sdigits),
                  round(AIC(cori[[j]][[i]]$mod),1))
        daic <- c(daic,AIC(cori[[j]][[i]]$mod))
        }else{
            alli <- c(alli, NA, NA)
            daic <- c(daic, NA)
        }
        tmp <- c(tmp, alli)
        df6 <- rbind(df6, tmp)
    }
    daic.all <- c(daic.all,daic - min(daic,na.rm = TRUE))
}
rownames(df6) <- NULL
colnames(df6) <- NULL


df6

df6[-1,2][daic.all == 0 & !is.na(daic.all)]


corsEnv.lag0.cm <- get.cors(resList, stockData, stockInfo,
                        as.list(rep("CM",7)),
                        envData = envData,
                        maxYear = 2021,
                        var = c("m","temp.bot","sal.bot"),
                        fit.poly = TRUE,
                        mod.link = "identity",
                        lag = 0, pres = presList)

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("cor_temp_sal_cm_",label,".png")),
        width = width, height = height, res = res)
plot.cor2(resList, stockData, stockInfo, selMods = bestMods,
          corsEnv.lag0.cm,
          leg.ncol = 2,
          mfrow = c(4,2), byrow = FALSE,
          legend.extra = TRUE)
if(save.figs) dev.off()




sdigits <- 2
scenarios <- c("Baseline",
               "Oxygen",
               "Surface",
               "Linear",
               "tvr",
               "tvK")
nlags <- length(scenarios)
cori <- vector("list", nlags)
for(i in 1:nlags){
    if(i == 1){
        fit.poly <- TRUE
        vars <- c("m","temp.bot","sal.bot")
        besti <- bestMods
    }else if(i == 2){
        fit.poly <- TRUE
        vars <- c("m","oxy.bot","sal.bot")
        besti <- bestMods
    }else if(i == 3){
        fit.poly <- TRUE
        vars <- c("m","temp.surf","sal.surf")
        besti <- bestMods
    }else if(i == 4){
        fit.poly <- FALSE
        vars <- c("m","temp.bot","sal.bot")
        besti <- bestMods
    }else if(i == 5){
        fit.poly <- TRUE
        vars <- c("r","temp.bot","sal.bot")
        besti <- as.list(rep("RM", nstocks))
    }else if(i == 6){
        fit.poly <- TRUE
        vars <- c("K","temp.bot","sal.bot")
        besti <- as.list(rep("KM", nstocks))
    }
    cori[[i]] <- get.cors(resList, stockData, stockInfo,
                     besti,
                     envData = envData,
                     maxYear = 2021,
                     var = vars,
                     fit.poly = fit.poly,
                     lag = 0)
}

df8 <- c("Stock","Scenario","Intercept","Temp/Oxy","Temp2/Oxy2","Sal","R2","AIC")
daic.all <- NULL
for(i in 1:nstocks){
    daic <- NULL
    for(j in 1:nlags){
        tmp <- c(stocks[i], scenarios[j])
        alli <- rep(NA, 4)
        if(!is.null(coefficients(cori[[j]][[i]]$mod))){
            tmpi <- signif(coefficients(cori[[j]][[i]]$mod),sdigits)
            if(scenarios[j] == "Linear"){
                alli[match(names(tmpi),
                           c("(Intercept)",
                             "expl1",
                             "",
                             "expl2"))] <- tmpi
            }else{
                alli[match(names(tmpi),
                           c("(Intercept)",
                             "poly(expl1, degree = degree)1",
                             "poly(expl1, degree = degree)2",
                             "expl2"))] <- tmpi
            }
            alli <- c(alli,signif(with(summary(cori[[j]][[i]]$mod),
                                       1 - deviance/null.deviance),sdigits),
                      round(AIC(cori[[j]][[i]]$mod),1))
            daic <- c(daic,AIC(cori[[j]][[i]]$mod))
        }else{
            alli <- c(alli, NA,NA)
            daic <- c(daic, NA)
        }
        tmp <- c(tmp, alli)
        df8 <- rbind(df8, tmp)
    }
    daic.all <- c(daic.all,daic - min(daic, na.rm = TRUE))
}
rownames(df8) <- NULL
colnames(df8) <- NULL

df8[is.na(df8)] <- ""

## remove AIC
df8 <- df8[,-8]

df8

tmp <- df8[-1,]
colnames(tmp) <- c("Stock", "Scenario", "Intercept",
                   "Temp/Oxy", "Temp\\textsuperscript{2}/Oxy\\textsuperscript{2}",
                   "Salinity", "R\\textsuperscript{2}")

rle.lengths <- rle(tmp[,1])$lengths
first <- !duplicated(tmp[,1])
tmp[,1][!first] <- ""
tmp[,1][first] <-
   paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", tmp[,1][first], "}}")


if(save.tables) xtable:::print.xtable(xtable::xtable(tmp, type = "latex",
                                                     caption = "Parameters and McFadden's R\\textsuperscript{2} for the best model for each stock, a time lag of 0 and various scenarios. Scenario 'Oxygen' replaces the polynomial term for bottom temperature with bottom oxygen in mol m\\textsuperscript{-3}, scenario 'Surface' uses the sea surface temperature and salinity values instead of sea bottom values, scenario 'Linear' replaces the polynomial term for sea bottom temperature with a linear term for sea bottom temperature, scenarion 'tvr' correlates the intrinsic growth rate (r) predicted by RM with sea bottom temperature and salinity, and scenario 'tvK' correlates the carrying capacity (K) predicted by KM with sea bottom temperature and salinity. All predictions assume a time lag of 0 between response and explanatory variables.",
                                                     label = "tab:si-cor-other",
                                                     align="ll|cccccc"),
                      file = file.path(resdir, "tabs", label,
                                       paste0("SI_table_cor_other_",label,".tex")),
                      include.rownames = FALSE, include.colnames = TRUE,
                      sanitize.text.function = identity,
                      hline.after = c(-1,nrow(tmp)),
                      caption.placement = "top",
                      floating=TRUE,
                      booktabs=TRUE,
                      size="\\fontsize{8pt}{8pt}\\selectfont"
                      )





indi <- envData$year >= 1981 &
                 envData$year <= 2021

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("env_cors_",label,".png")),
        width = 700, height = 1000, res = res)
par(mfrow = c(2,1), mar = c(2,2,1,1), oma = c(3,3,1,1))
plot(envData$temp.bot[indi],
     envData$oxy.bot[indi],
     xlab = "", ylab = "")
mtext("Dissolved sea bottom oxygen", 2, 3)
plot(envData$temp.bot[indi],
     envData$sal.bot[indi],
     xlab = "", ylab = "")
mtext("Sea bottom salinity", 2, 3)
mtext("Sea bottom temperature", 1, 3)
if(save.figs) dev.off()



cor.test(envData$temp.bot[indi],
     envData$oxy.bot[indi])

corsF0 <- get.cors(resList, stockData, stockInfo,
                  bestMods,
                  envData = envData,
                  maxYear = 2021,
                  var = c("m","F2"),
                  fit.poly = TRUE,
                  lag = 0)

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("cor_f_lag0_",label,".png")),
        width = width, height = height, res = res)
plot.cor2(resList, stockData, stockInfo, selMods = bestMods,
          corsF0,
          leg.ncol = 2,
          xlim = c(0.09,1.5),
          mfrow = c(4,2), byrow = FALSE,
          legend.extra = TRUE)
if(save.figs) dev.off()


corsFFmsy0 <- get.cors(resList, stockData, stockInfo,
                  bestMods,
                  envData = envData,
                  maxYear = 2021,
                  var = c("m","FFmsy2"),
                  fit.poly = TRUE,
                  lag = 0)

if(save.figs)
    png(file.path(resdir, "figs", label,
                  paste0("cor_ffmsy_lag0_",label,".png")),
        width = width, height = height, res = res)
plot.cor2(resList, stockData, stockInfo, selMods = bestMods,
          corsFFmsy0,
          leg.ncol = 2,
          xlim = c(0.09,5),
          mfrow = c(4,2), byrow = FALSE,
          legend.extra = TRUE)
if(save.figs) dev.off()
