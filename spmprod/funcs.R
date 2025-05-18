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

## This scripts includes additional functions.

plotspict.diagnostic.process <- function(rep, lag.max=4, qlegend=TRUE, plot.data=TRUE, mfcol=FALSE,
                                         add.loess = FALSE, span = 0.75, stamp=get.version()){

    spict:::check.rep(rep)
    inp <- rep$inp

    ## Check that process residuals calculated
    if(!any(names(rep) == "process.resid")){
        stop("No process residuals found. Please run the function 'process.resid' on your fitted spict object.")
    }

    testsB <- res.diagn(rep$process.resid$B, "B", "B")
    testsF <- res.diagn(rep$process.resid$F, "F", "F")


    fun <- function(time, res, add=FALSE, add.legend=FALSE, col=1, pch=1,
                    add.vline.at=NULL, ...){
        nrem <- length(time) - length(res)
        if (nrem > 0){
            time <- time[-nrem]
            warning('length of residual vector and length of corresponding time vector are not equal!')
        }
        spict:::plot.col(time, res, pch=pch, add=add, add.legend=add.legend, typ='p', xlab='Time',
                 add.vline.at=add.vline.at, ...)
        dum <- rep(NA, length(res))
        dum[is.na(res)] <- 0
        text(time, dum, labels='NA', cex=0.8, col=col)
    }

    time.full <- inp$time
    dt <- diff(rep$process.resid$time)[1]
    time.agg <- time.full[which(time.full %% (dt) == 0)]
    time.rep <- rep(time.agg, each = 1/(inp$dteuler/dt))
    nta <- length(time.agg)
    logB <- get.par("logB",rep)[,2]
    logF <- get.par("logF",rep)[,2]
    logFs <- get.par("logFs",rep)[,2]
    P <- rep$obj$report()$P
    logB.agg <- logF.agg <- rep(0, nta)
    P.agg <- logFs.agg <- rep(0, nta)
    for (i in 1:nta){
        inds <- which(time.agg[i]==time.rep)
        logB.agg[i] <- mean(logB[inds])
        logF.agg[i] <- mean(logF[inds])
        logFs.agg[i] <- mean(logFs[inds])
        P.agg[i] <- sum(P[inds])
    }

    mar <- c(4.7, 4.1, 2.5, 2)
    mfrow <- c(3 + as.numeric(plot.data), 2)
    if (mfcol){
        opar <- par(mfcol=rev(mfrow), mar=mar)
    } else {
        opar <- par(mfrow=mfrow, mar=mar)
    }
    on.exit(par(opar))

    ## Plot data
    if(plot.data){
        spict:::plot.col(time.agg, logB.agg, ylab='logB',
                 main='Biomass', xlab='Time')
        spict:::plot.col(time.agg, logF.agg, ylab='logF',
                 main='Fishing mortality', xlab='Time')
    }

    pval <- round(testsB$biasB.p, 4)
    colmain <- ifelse(pval < 0.05, 'red', 'forestgreen')
    fun(rep$process.resid$time, rep$process.resid$B,
        add.legend=qlegend, ylab='B residuals',
        main=paste0('Bias p-val: ', pval),
        col.main=colmain)
    abline(h=0, lty=3)
    if(add.loess){
        mod <- loess(B ~ time, data=rep$process.resid, span = span)
        j <- order(rep$process.resid$time)
        lines(rep$process.resid$time[j],mod$fitted[j], lwd=1.5)
    }
    pval <- round(testsF$biasF.p, 4)
    colmain <- ifelse(pval < 0.05, 'red', 'forestgreen')
    fun(rep$process.resid$time, rep$process.resid$F,
        add.legend=FALSE, ylab='F residuals',
        main=paste0('Bias p-val: ', pval),
        col.main=colmain)
    abline(h=0, lty=3)
    if(add.loess){
        mod <- loess(F ~ time, data=rep$process.resid, span = span)
        j <- order(rep$process.resid$time)
        lines(rep$process.resid$time[j], mod$fitted[j], lwd=1.5)
    }

    pvalacfB <- round(testsB$LBoxB.p, 4)
    spict:::osar.acf.plot(rep$process.resid$B, lag.max, pvalacfB, ylab='B ACF')
    pvalacfF <- round(testsF$LBoxF.p, 4)
    spict:::osar.acf.plot(rep$process.resid$F, lag.max, pvalacfF, ylab='F ACF')

    pvalB <- round(testsB$shapiroB.p, 4)
    spict:::osar.qq.plot(rep$process.resid$B, pvalB)
    pvalF <- round(testsF$shapiroF.p, 4)
    spict:::osar.qq.plot(rep$process.resid$F, pvalF)

    spict:::txt.stamp(stamp, do.flag=TRUE)
}


calc.process.resid <- function(rep, dt = NULL){

    repin <- rep

    spict:::check.rep(rep)
    inp <- rep$inp

    tmp <- 2^(0:10)
    if(is.null(dt)){
        dt <- (1/inp$nseasons)
    }else if(!(dt %in% (1/(tmp)[which(tmp <= 1/inp$dteuler)]))){
        stop(paste0("Don't know how to handle dt=",dt,". dt has to be between dteuler and 1, and a divisor of inp$dteuler."))
    }

    ## helper function to draw multivariate normal samples
    rmvnorm <- function(n = 1, mu, Sigma){
        p <- length(mu)
        if(!all(dim(Sigma) == c(p, p))){
            stop("incompatible arguments")
        }
        idx <- diag(Sigma) > .Machine$double.xmin
        L <- matrix(0,p,p)
        if(any(idx)){
            L[idx,idx] <- chol(Sigma[idx,idx])
        }
        X <- matrix(rnorm(p * n), n)
        X <- drop(mu) + t(X %*% L)
        if(n == 1){
            drop(X)
        }else{
            t(X)
        }
    }

    if('sderr' %in% names(repin))
        cat('WARNING: Could not calculate standard deviations. The optimum found may be invalid. Proceed with caution.\n')

    rep.co <- list()
    rep.co$obj <- make.obj(make.datin(repin$inp), repin$pl, repin$inp, phase=1)
    rep.co$obj$env$data$residFlag <- 1
    rep.co$obj$retape()
    class(rep.co) <- "spictcls"
    sdrep <- sdreport(rep.co$obj, repin$opt$par)

        ## Use same seed if not specified in global env
    if (exists(".Random.seed")){
        oldseed <- get(".Random.seed", .GlobalEnv)
        oldRNGkind <- RNGkind()
    }
    set.seed(123456)

    ## Biomass
    idx <- which(names(sdrep$value) == "residB")
    residB <- rmvnorm(1, mu = sdrep$value[idx], Sigma = sdrep$cov[idx,idx])

    ## Fishing mortality
    idx <- which(names(sdrep$value) == "residF")
    residF <- rmvnorm(1, mu = sdrep$value[idx], Sigma = sdrep$cov[idx,idx])


    if (exists("oldseed")){
        do.call("RNGkind",as.list(oldRNGkind))
        assign(".Random.seed", oldseed, .GlobalEnv)
    }

    time.full <- inp$time[-length(inp$time)]
    time.agg <- time.full[which(time.full %% (dt) == 0)]
    time.rep <- rep(time.agg, each = 1/(inp$dteuler/dt))
    nta <- length(time.agg)
    residB.agg <- residF.agg <- rep(0, nta)
    for (i in 1:nta){
        inds <- which(time.agg[i]==time.rep)
        residB.agg[i] <- mean(residB[inds])
        residF.agg[i] <- mean(residF[inds])
    }

    res <- repin
    res$process.resid <- as.data.frame(cbind(
        time = time.agg,
        B = residB.agg,
        F = residF.agg
    ))

    if(!'diagn' %in% names(res)) res$diagn <- list()
    diagnF <- res.diagn(res$process.resid$F, "F", "F")
    for (nm in names(diagnF)){
        res$diagn[[nm]] <- diagnF[[nm]]
    }
    diagnB <- res.diagn(res$process.resid$B, "B", "B")
    for (nm in names(diagnB)){
        res$diagn[[nm]] <- diagnB[[nm]]
    }

    return(res)
}


## AIC
AIC.spictcls <- function(x) 2*x$opt$objective + 2*length(x$opt$par)


## AICc
AICc <- function(x){
    k <- length(x$opt$par)
    n <- length(x$inp$obssrt)
    2*x$opt$objective + 2*k + (2*k^2 + 2*k) / (n - k - 1)
}


## Colours and characters
spmprod.cols <- function(){
    cols <- c("grey70","#e6194b", "#3cb44b", "#0082c8", "#f58231",
              "#911eb4", "#46f0f0",
              "#f032e6", "#d2f53c", "#008080", "#aa6e28","darkgreen", "#8c591e",
              "#800000","#808080", "#000000")
    cols <- c("grey70","#0072bd", "#d95319","#edb120", "#7e2f8e", "#77ac30")
    cols <- c("grey70",  ## "firebrick3",
              ##           "chartreuse4",
              ## "darkolivegreen4",
              "dodgerblue",
              "#41a400", ## "forestgreen",
              "chocolate4", ##"plum3", ## "purple3",
              ## "goldenrod1",
              "#e7ab13", ## "#ffbd12",
              "tomato2", ## "firebrick2",
              "darkorchid3",
              "darkorange2",
              "dodgerblue4",
              "grey40",
              "orchid2")
    cols <- c("grey50",  ## "firebrick3",
              "dodgerblue1",
              "#41a400", ## "forestgreen",
              "tomato2", ## "firebrick2",
              "#e7ab13", ## "#ffbd12",
              "purple2",
              "darkorange2", ## "purple1",
              "darkorchid3",
              "darkorange2",
              "dodgerblue4",
              "grey40",
              "orchid2")
    return(cols)
}
spmprod.pchs <- function(){
    return(c(3,21:25))
}
spmprod.ltys <- function(){
    return(rep(1,7))
}


any.invalid <- function(x){
    if(all(is.na(x))) return(FALSE)
    any(is.na(x[,c(1:3,5)]))
}

any.invalid2 <- function(x, cv.limit = 5){
    if(all(is.na(x))) return(FALSE)
    any(is.infinite(x[,c(1:3,5)])) || any(x[,5] > cv.limit)
}

any.invalid3 <- function(x, om.limit = 2){
    if(all(is.na(x))) return(FALSE)
    any(abs(floor(log(abs(x[,3]),10)) - floor(log(abs(x[,1]),10))) > om.limit)
}


## kobe plot
plot.kobe <- function(resList, stockData, stockInfo, selMods = NULL,
                      selYear = 2021, add.ices = TRUE, xlim = NULL,
                      ylim = NULL, xlab = expression("B/"*B[MSY]),
                      ylab = expression("F/"*F[MSY]), cex = 1.5,
                      yfixed = TRUE, xfixed = TRUE, mfrow = c(2,5),
                      age.ref.period = FALSE, add.traj = FALSE,
                      byrow = FALSE){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
    }else{
        resListSel <- resList
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    if(length(add.traj) < max(sapply(selMods,length))){
        add.traj <- rep(add.traj[1], max(sapply(selMods,length)))
    }

    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,2,0.5,1), oma = c(5,2,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        par(mar = c(3,3,2,2), oma = c(3,2,2,4))
    }
    for(i in 1:length(resListSel)){
        if(is.null(xlim0)) xlim <- c(0,5)
        if(is.null(ylim0)) ylim <- c(0,5)
        ## if(!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow))
        ##     xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        if(!yfixed || i %in% sapply(seq(mfrow[1]),function(x) (x-1) * (mfrow[2]) + 1))
            yaxt = "s" else yaxt = "n"
        plot(1, 1,
             ty = 'n',
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = yaxt,
             xlab = "",
             ylab = ""
             )
        abline(h = 1, col = "grey80")
        abline(v = 1, col = "grey80")
        for(j in 1:length(resListSel[[i]])){
            years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            if(age.ref.period){
                indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
            }else{
                indi <- years >= selYear & years < selYear + 1
            }
            xdat <- mean(resListSel[[i]][[j]]$logBBmsy[indi,2])
            ydat <- mean(resListSel[[i]][[j]]$logFFmsynotS[indi,2])

            if(!age.ref.period){
            polygon(mean(resListSel[[i]][[j]]$logB[indi,2]) /
                    exp(resListSel[[i]][[j]]$kobe[,1]),
                    mean(resListSel[[i]][[j]]$logFnotS[indi,2]) /
                    exp(resListSel[[i]][[j]]$kobe[,2]),
                    col = adjustcolor(cols[indiList[[i]][j]], 0.2),
                    border = NA)
            }
            points(xdat, ydat, col = cols[indiList[[i]][j]],
                   bg = cols[indiList[[i]][j]],
                   pch = pchs[indiList[[i]][j]], cex = cex)
            if(add.traj[j]){
                indi <- years-(years%%1)
                indi <- indi[indi <= selYear]
                xdat <- as.vector(by(resListSel[[i]][[j]]$logBBmsy[1:length(indi),2], indi, mean))
                ydat <- as.vector(by(resListSel[[i]][[j]]$logFFmsynotS[1:length(indi),2], indi, mean))
                lines(xdat, ydat, col = adjustcolor(cols[indiList[[i]][j]],0.4),
                      lty = 1.5)
            }
        }
        if(add.ices){
            years <- stockData[[i]]$data$year
            if(age.ref.period){
                indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
            }else{
                indi <- which.min((years - selYear)^2)
            }
            xdat <- stockData[[i]]$data$SSB[indi] /
                stockInfo[i,"EQ_Bmsy"]
            ydat <- stockData[[i]]$data$Fbar[indi] / stockInfo[i,"Fmsy"]
            points(mean(xdat), mean(ydat), col = 1, pch = 1, cex = 1.1)
            xdat <- stockData[[i]]$data$ESB_wCatch[indi] /
                stockInfo[i,"EQ_ESBmsy"]
            ydat <- stockData[[i]]$data$Fesb[indi] / stockInfo[i,"EQ_FmsyESB"]
            points(mean(xdat), mean(ydat), col = 1, pch = 8, cex = 1.1)
        }
        legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)
        if(i == mfrow[2]){
            leg <- unique(unlist(lapply(resListSel, function(x) names(x))))
            colleg <- unique(unlist(lapply(indiList, function(x) cols[x])))
            pchleg <- unique(unlist(lapply(indiList, function(x) pchs[x])))
            if(add.ices){
                leg <- c(leg, "APM","APM (ESB)")
                colleg <- c(colleg, "black","black")
                pchleg <- c(pchleg, 8, 1)
            }
            legend("topright", legend = leg,
                   pch = pchleg,
                   col = colleg,
                   pt.bg = colleg,
                   pt.cex = 1.3,
                   bg = "white",
                   ncol = 2,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
        box(lwd = 1.5)
    }
    mtext(xlab,1,3,outer = TRUE)
    mtext(ylab,2,3,outer = TRUE)
}


get.ices.data <- function(stockData, stockInfo, var, i, maxYear){
    vars.ices <- c("logFmsyvec","logBmsyvec","logMSYvec",
                   "logFFmsynotS","logBBmsy","cmsy",
                   "logFnotS","logB","logC")
    if(any(var %in% vars.ices)){
        years <- c(stockData[[i]]$data$year,
                   max(stockData[[i]]$data$year) + 0.9)
        ny <- length(years)
        if(any(var %in% vars.ices[1:6])){
            indi <- which(years <= maxYear & years >= stockInfo$EQ_start[i])
        }else{
            indi <- which(years <= maxYear)
        }
        if(any(var == "logBBmsy")){
            ydat <- stockData[[i]]$data$SSB / stockInfo[i,"EQ_Bmsy"]
        }else if(any(var == "logFFmsynotS")){
            ydat <- stockData[[i]]$data$Fbar / stockInfo[i,"EQ_Fmsy"]
        }else if(any(var == "logFmsyvec")){
            ydat <- rep(stockInfo[i,"EQ_FmsyESB"], ny)
        }else if(any(var == "logBmsyvec")){
            ydat <- rep(stockInfo[i,"EQ_ESBmsy"], ny)
        }else if(any(var == "logMSYvec")){
            ydat <- rep(stockInfo[i,"EQ_MSY"], ny)
        }else if(any(var == "logFnotS")){
            ydat <- stockData[[i]]$data$Fesb
        }else if(any(var == "logB")){
            ydat <- stockData[[i]]$data$ESB_raw_ov_sWeight_new
        }else if(any(var == "cmsy")){
            ydat <- stockData[[i]]$data$C / stockInfo[i,"EQ_MSY"]
        }else if(any(var == "logC")){
            ydat <- stockData[[i]]$data$C
        }
        res <- list(years = years[indi],
                    ydat = ydat[indi])
    }else{
        res <- NULL
    }
    return(res)
}



## time plot single
plot.time.single <- function(resList, stockData, stockInfo, selMods = NULL,
                             maxYear = 2021, var = "logrvec", scale = 1,
                             rel = 1,
                             add.ices = TRUE, xlim = NULL,
                             ylim = NULL, xlab = "Time",
                             ylab = "Instrinsic growth rate",
                             xfixed = TRUE, yfixed = FALSE, plot.ci = TRUE,
                             plot.obs = FALSE, mfrow = c(2,5),
                             plot.ref = FALSE,
                             remove.cm.roundfish = FALSE,
                             byrow = FALSE,
                             legend.extra = FALSE){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()


    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
       names(resListSel) <- names(resList)
    }else{
        resListSel <- resList
        names(resListSel) <- names(resList)
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed && length(var) > 1){
        par(mar = c(0.5,3,0.5,2), oma = c(4.2,2.6,2.5,3.6))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,2,0.5,1), oma = c(4.2,2.6,2.5,0.6)) ##  oma = c(5,3.5,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        par(mar = c(3,3,2,2), oma = c(3,2,2,4))
    }
    for(i in 1:length(resListSel)){
        if(is.null(xlim0)) xlim <- c(1958,2021)
        if(is.null(ylim0) && plot.ci){
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish %in% c(1,2)){
                ylim <- range(unlist(lapply(resListSel[[i]][names(resListSel[[i]]) != "CM"], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
                if(remove.cm.roundfish == 2){
                    ylim2 <- range(unlist(lapply(resListSel[[i]], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
                    ylim <- range(ylim,ylim2)
                }
            }else{
                ylim <- range(unlist(lapply(resListSel[[i]],
                                            function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
            }
        }else if(is.null(ylim0)){
            ylim <- range(unlist(lapply(resListSel[[i]], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
        }
        if(add.ices || plot.obs){
            ices.dat <- get.ices.data(stockData, stockInfo, var, i, maxYear)
            ylim <- range(ylim, ices.dat$ydat / scale[1] / rel)
        }
        if(any(is.infinite(ylim))) ylim <- c(0,1)
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        ylim <- c(1,1.05) * ylim
        if(!yfixed && i == mfrow[2]) ylim <- c(1,1.05) * ylim
        plot(1, 1,
             ty = "n",
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = "n",
             xlab = "",
             ylab = ""
             )
        ati <- pretty(c(0.9, 1.1) * ylim)
        if(plot.ref){
            abline(h = 1, col = adjustcolor("grey70",1), lwd = 1, lty = 2)
        }else{
            abline(h = ati, col = adjustcolor("grey70",0.3), lwd = 1, lty = 2)
        }
        if(!yfixed || i %in% sapply(seq(mfrow[1]),function(x) (x-1) * (mfrow[2]) + 1))
            axis(2, at = ati, labels = ati)
        if(!yfixed && length(var) > 1)
            axis(4, at = ati, labels = round(ati * rel,2))
        for(j in 1:length(resListSel[[i]])){
            years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            indi <- years < maxYear
            years <- years[indi]
            ydat <- resListSel[[i]][[j]][[var[1]]][which(indi),1:3] / scale[1]
            if(plot.ci){
                if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                   remove.cm.roundfish %in% c(1,2) && names(resListSel[[i]])[j] == "CM"){
                }else{
                    polygon(c(years, rev(years)),
                            c(ydat[,1],rev(ydat[,3])),
                            col = adjustcolor(cols[indiList[[i]][j]], 0.2),
                            border = NA)
                }
            }
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish == 1 && names(resListSel[[i]])[j] == "CM"){
            }else{
                lines(years, ydat[,2], col = cols[indiList[[i]][j]],
                      lwd = 2)
            }
        }
        if(add.ices || plot.obs){
            if(!is.null(ices.dat)){
                if(plot.obs){
                    points(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                           col = "grey70")
                }else if(add.ices){
                    lines(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                          col = 1, lwd = 2, lty = 2)
                }
            }
        }
        legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)
        if(i == mfrow[2] && !legend.extra){
            texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
            legend("topright", legend = texti,
                   lty = unique(unlist(lapply(indiList, function(x) ltys[x]))),
                   col = unique(unlist(lapply(indiList, function(x) cols[x]))),
                   lwd = 2,
                   bg = "white",
                   ncol = 3,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
        box(lwd = 1.5)
        if(i == 1) mtext("Roundfish",3,0.5, font = 2,cex = 0.94)
        if(i == mfrow[1]+1) mtext("Flatfish",3,0.5, font = 2,cex = 0.94)
    }
    mtext(xlab,1,2.6,outer = TRUE)
    mtext(ylab[[1]],2,0.1,outer = TRUE)
    if(length(var) > 1) mtext(ylab[[2]],4,1.6,outer = TRUE)
    if(legend.extra){
        plot.new()
        texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
        legend("center", legend = texti,
               title = "Models", title.cex = 1.2, title.font = 1,
               lty = unique(unlist(lapply(indiList, function(x) ltys[x]))),
               col = unique(unlist(lapply(indiList, function(x) cols[x]))),
               lwd = 2,
               bg = "white",
               ncol = 2,
               x.intersp = 1.1,
               y.intersp = 1.2,
               cex = 1.2)
    }
}


## time plot single
plot.time.single2 <- function(resList, stockData, stockInfo, selMods = NULL,
                             maxYear = 2021, var = "logrvec", scale = 1,
                             rel = 1,
                             add.ices = TRUE, xlim = NULL,
                             ylim = NULL, xlab = "Time",
                             ylab = "Instrinsic growth rate",
                             xfixed = TRUE, yfixed = FALSE, plot.ci = TRUE,
                             plot.obs = FALSE, mfrow = c(2,5),
                             plot.ref = FALSE,
                             remove.cm.roundfish = FALSE,
                             byrow = FALSE,
                             legend.extra = FALSE,
                             plot.legend = TRUE,
                             adj = c(0,0),
                             title = NULL,
                             title.line = 2,
                             stock.in.plot = TRUE,
                             stock.in.plot.line = 3,
                             bestMods = NULL,
                             leg.text = NULL,
                             drop.mods = NULL, gap = TRUE){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()


    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
       names(resListSel) <- names(resList)
    }else{
        resListSel <- resList
        names(resListSel) <- names(resList)
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    for(i in 1:length(resListSel)){
        if(is.null(xlim0)) xlim <- c(1958,2021)
        if(is.null(ylim0) && plot.ci){
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish %in% c(1,2)){
                ylim <- range(unlist(lapply(resListSel[[i]][names(resListSel[[i]]) != "CM"], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
                if(remove.cm.roundfish == 2){
                if(var[1] == "logMSYvec"){
                    varii <- "logmvec"
                }else{
                    varii <- var[1]
                }
                    ylim2 <- range(unlist(lapply(resListSel[[i]], function(x) x[[varii]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
                    ylim <- range(ylim,ylim2)
                }
            }else{
                if(var[1] == "logMSYvec"){
                    varii <- "logmvec"
                }else{
                    varii <- var[1]
                }
                ylim <- range(unlist(lapply(resListSel[[i]],
                                            function(x) x[[varii]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
            }
        }else if(is.null(ylim0)){
            ylim <- range(unlist(lapply(resListSel[[i]], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
        }
        if(add.ices || plot.obs){
            ices.dat <- get.ices.data(stockData, stockInfo, var, i, maxYear)
            ylim <- range(ylim, ices.dat$ydat / scale[1] / rel)
        }
        if(any(is.infinite(ylim))) ylim <- c(0,1)
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        ylim <- c(1,1.05) * ylim
        if(!yfixed && i == mfrow[2]) ylim <- c(1,1.05) * ylim
        plot(1, 1,
             ty = "n",
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = "n",
             xlab = "",
             ylab = ""
             )
        ati <- pretty(c(0.9, 1.1) * ylim)
        if(plot.ref){
            abline(h = 1, col = adjustcolor("grey70",1), lwd = 1, lty = 2)
        }else{
            abline(h = ati, col = adjustcolor("grey70",0.3), lwd = 1, lty = 2)
        }
        if(!yfixed || i %in% sapply(seq(mfrow[1]),function(x) (x-1) * (mfrow[2]) + 1))
            axis(2, at = ati, labels = ati)
        if(!yfixed && length(var) > 1)
            axis(4, at = ati, labels = round(ati * rel,2))
        if(!is.null(bestMods)){
            jbest <- which(names(resListSel[[i]]) == bestMods[[i]])
        }else{
            jbest <- NULL
        }
        for(j in 1:length(resListSel[[i]])){
            years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            indi <- years < maxYear
            years <- years[indi]
            if(var[1] == "logMSYvec"){
                varii <- "logmvec"
            }else{
                varii <- var[1]
            }
            ydat <- resListSel[[i]][[j]][[varii]][which(indi),1:3] / scale[1]
            if(is.null(jbest) || jbest != j){
                if(plot.ci){
                    if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                       remove.cm.roundfish %in% c(1,2) && names(resListSel[[i]])[j] == "CM"){
                    }else{
                        polygon(c(years, rev(years)),
                                c(ydat[,1],rev(ydat[,3])),
                                col = adjustcolor(cols[indiList[[i]][j]], 0.15),
                                border = NA)
                    }
                }
                if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                   remove.cm.roundfish == 1 && names(resListSel[[i]])[j] == "CM"){
                }else{
                    lines(years, ydat[,2],
                          col = cols[indiList[[i]][j]],
                          ## lty = c(2,1,2,1,2,2)[j],
                          lwd = 2)
                }
            }
        }
        if(add.ices || plot.obs){
            if(!is.null(ices.dat)){
                if(plot.obs){
                    points(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                           col = "grey70")
                }else if(add.ices){
                    lines(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                          col = 1, lwd = 1, lty = 3)
                }
            }
        }
        if(!is.null(bestMods)){
            j <- which(names(resListSel[[i]]) == bestMods[[i]])
            if(length(j) > 0){
                years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
                indi <- years < maxYear
                years <- years[indi]
                if(var[1] == "logMSYvec"){
                    varii <- "logmvec"
                }else{
                    varii <- var[1]
                }
                ydat <- resListSel[[i]][[j]][[varii]][which(indi),1:3] / scale[1]
                if(plot.ci){
                    if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                       remove.cm.roundfish %in% c(1,2) && names(resListSel[[i]])[j] == "CM"){
                    }else{
                        polygon(c(years, rev(years)),
                                c(ydat[,1],rev(ydat[,3])),
                                col = adjustcolor(cols[indiList[[i]][j]], 0.15),
                                border = NA)
                    }
                }
                if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                   remove.cm.roundfish == 1 && names(resListSel[[i]])[j] == "CM"){
                }else{
                    lines(years, ydat[,2],
                          col = cols[indiList[[i]][j]],
                          lwd = 2)
                }
                points(tail(years,1), tail(ydat[,2],1), pch = "*",
                       col = cols[indiList[[i]][j]],
                       cex = 3)
            }
        }
        if(as.integer(stock.in.plot) == 1)
            legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)
        if(i == 1 && !legend.extra && plot.legend){
            if(is.null(leg.text)){
                texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
            }else{
                texti <- leg.text
            }
            if(!is.null(drop.mods)) texti <- texti[-drop.mods]
            ltyi <- unique(unlist(lapply(indiList, function(x) ltys[x])))
            if(!is.null(drop.mods)) ltyi <- ltyi[-drop.mods]
            coli <- unique(unlist(lapply(indiList, function(x) cols[x])))
            if(!is.null(drop.mods)) coli <- coli[-drop.mods]
            legend("topright", legend = texti,
                  lty = ltyi,
                   ## lty = c(2,1,2,1,1,2),
                   col = coli,
                   lwd = 2,
                   bg = "white",
                   ncol = 1,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
        box(lwd = 1.5)
        if(i == 1 && !is.null(title)) mtext(title[[1]], 3, title.line, font = 1)
        if(i == 1 && !is.null(title) && length(title) > 1)
            mtext(title[[2]], 3, 0.5, font = 1)
        if(as.integer(stock.in.plot) == 2)
            mtext(names(resList)[i], 4, stock.in.plot.line, font = 1)
        if(gap && i == 4) plot.new()
    }
    mtext(xlab,1,2.6,outer = TRUE)
    mtext(ylab[[1]],2,0,outer = TRUE, padj = adj[1])
    if(length(var) > 1 && length(ylab) > 1) mtext(ylab[[2]],4,0,outer = TRUE, padj = adj[2])
    if(legend.extra){
        plot.new()
        texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
        if(!is.null(drop.mods)) texti <- texti[-drop.mods]
        ltyi <- unique(unlist(lapply(indiList, function(x) ltys[x])))
        if(!is.null(drop.mods)) ltyi <- ltyi[-drop.mods]
        coli <- unique(unlist(lapply(indiList, function(x) cols[x])))
        if(!is.null(drop.mods)) coli <- coli[-drop.mods]
        legend("center", legend = texti,
               title = "Models", title.cex = 1.2, title.font = 1,
               lty = ltyi,
               ## lty = c(2,1,2,1,1,2),
               col = coli,
               lwd = 2,
               bg = "white",
               ncol = 2,
               x.intersp = 1.1,
               y.intersp = 1.2,
               cex = 1.2)
    }
}


## Get spict variables
var2vars <- function(var){
    nvar <- length(var)

    res <- list()
    res$var.ori <- res$var <- res$lab <- res$scale <- NULL
    for(i in 1:nvar){
        res$var.ori[i] <- var[i]
        if(var[i] == "r"){
            res$var[i] <- "logrvec"
            res$lab[[i]] <- expression("Intrinsic growth rate ["*yr^-1*"]")
            res$scale[i] <- 1
        }else if(var[i] == "K"){
            res$var[i] <- "logKvec"
            res$lab[[i]] <- expression("Carrying capacity ["*"'000 t"*"]")
            res$scale[i] <- 1e3
        }else if(var[i] == "m"){
            res$var[i] <- "logmvec"
            res$lab[[i]] <- expression("Maximum net productivity ['000 t"~yr^-1*"]")
            res$scale[i] <- 1e3
        }else if(var[i] == "Fmsy"){
            res$var[i] <- "logFmsyvec"
            res$lab[[i]] <- expression(F[MSY] * " ["*yr^-1*"]")
            res$scale[i] <- 1
        }else if(var[i] == "Bmsy"){
            res$var[i] <- "logBmsyvec"
            res$lab[[i]] <- expression(B[MSY] * " ['000 t]")
            res$scale[i] <- 1e3
        }else if(var[i] == "MSY"){
            res$var[i] <- "logMSYvec"
            res$lab[[i]] <- expression("MSY ['000 t]")
            res$scale[i] <- 1e3
        }else if(var[i] == "FFmsy"){
            res$var[i] <- "logFFmsynotS"
            res$lab[[i]] <- expression(F/F[MSY])
            res$scale[i] <- 1
        }else if(var[i] == "BBmsy"){
            res$var[i] <- "logBBmsy"
            res$lab[[i]] <- expression(B/B[MSY])
            res$scale[i] <- 1
        }else if(var[i] == "CMSY"){
            res$var[i] <- "cmsy"
            res$lab[[i]] <- expression(C/MSY)
            res$scale[i] <- 1
        }else if(var[i] == "C"){
            res$var[i] <- "logC"
            res$lab[[i]] <- expression("Catch ['000 t]")
            res$scale[i] <- 1e3
        }else if(var[i] == "F"){
            res$var[i] <- "logFnotS"
            res$lab[[i]] <- expression("Fishing mortality ["*yr^-1*"]")
            res$scale[i] <- 1
        }else if(var[i] == "F2"){
            res$var[i] <- "F2"
            res$lab[[i]] <- expression("Fishing mortality ["*yr^-1*"]")
            res$scale[i] <- 1
        }else if(var[i] == "FFmsy2"){
            res$var[i] <- "FFmsy2"
            res$lab[[i]] <- expression(F/F[MSY])
            res$scale[i] <- 1
        }else if(var[i] == "B"){
            res$var[i] <- "logB"
            res$lab[[i]] <- expression("Biomass ['000 t]")
            res$scale[i] <- 1e3
        }else if(var[i] == "temp.bot"){
            res$var[i] <- "temp.bot"
            res$lab[[i]] <- expression("Temperature [\u00B0C]")
            res$scale[i] <- 1
        }else if(var[i] == "sal.bot"){
            res$var[i] <- "sal.bot"
            res$lab[[i]] <- expression("Salinity [PSU]")
            res$scale[i] <- 1
        }else if(var[i] == "oxy.bot"){
            res$var[i] <- "oxy.bot"
            res$lab[[i]] <- expression("Oxygen [mol " * m^{-3}*"]")
            res$scale[i] <- 1
        }else if(var[i] == "temp.surf"){
            res$var[i] <- "temp.surf"
            res$lab[[i]] <- expression("Temperature [\u00B0C]")
            res$scale[i] <- 1
        }else if(var[i] == "sal.surf"){
            res$var[i] <- "sal.surf"
            res$lab[[i]] <- expression("Salinity [PSU]")
            res$scale[i] <- 1
        }else if(var[i] == "natMort"){
            res$var[i] <- "natMort"
            res$lab[[i]] <- expression("Natural mortality ["*yr^-1*"]")
            res$scale[i] <- 1
        }else if(var[i] == "mat"){
            res$var[i] <- "mat"
            res$lab[[i]] <- expression("Maturity")
            res$scale[i] <- 1
        }else if(var[i] == "wStock"){
            res$var[i] <- "wStock"
            res$lab[[i]] <- expression("Mean weight-at-age (stock) [kg]")
            res$scale[i] <- 1
        }else if(var[i] == "wCatch"){
            res$var[i] <- "wCatch"
            res$lab[[i]] <- expression("Mean weight-at-age (catch) [kg]")
            res$scale[i] <- 1
        }else if(var[i] == "ageComp"){
            res$var[i] <- "ageComp"
            res$lab[[i]] <- expression("Proportion of individuals of max. age")
            res$scale[i] <- 1
        }else if(var[i] == "selYoung"){
            res$var[i] <- "selYoung"
            res$lab[[i]] <- expression("Prob. of selection of min. age [%]")
            res$scale[i] <- 1/100
        }else if(var[i] == "selOld"){
            res$var[i] <- "selOld"
            res$lab[[i]] <- expression("Prob. of selection of max. age [%]")
            res$scale[i] <- 1/100
        }else if(var[i] == "rec"){
            res$var[i] <- "rec"
            res$lab[[i]] <- expression("Number of recruits ['000]")
            res$scale[i] <- 1e3
        }else if(var[i] == "SSB"){
            res$var[i] <- "SSB"
            res$lab[[i]] <- expression("Spawning stock biomass ['000 t]")
            res$scale[i] <- 1e3
        }else if(var[i] == "recSSB"){
            res$var[i] <- "recSSB"
            res$lab[[i]] <- expression("Number of recruits / SSB ["*t^-1*"]")
            res$scale[i] <- 1
        }else if(var[i] == "comb"){
            res$var[i] <- "comb"
            res$lab[[i]] <- ""
            res$scale[i] <- 1
        }else stop("Don't know var!")
    }

    return(res)
}


## time plot
plot.time <- function(resList, stockData, stockInfo, selMods = NULL,
                      maxYear = 2021, var = "r", rel = 1,
                      add.ices = TRUE, xlim = NULL,
                      ylim = NULL, xlab = "Time",
                      ylab = NULL,
                      xfixed = TRUE, yfixed = FALSE,
                      plot.ci = TRUE, plot.obs = FALSE, mfrow = c(2,5),
                      plot.ref = FALSE,
                      remove.cm.roundfish = FALSE,
                             byrow = FALSE,
                             legend.extra = FALSE){

    vars <- var2vars(var)

    plot.time.single(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = if(any(is.null(ylab))) vars$lab else ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra)

}


plot.time.all <- function(resList, stockData, stockInfo, selMods = NULL,
                          maxYear = 2021, nPar = 2,
                          adj = list(c(0,0),c(0,0),c(0,0))){

    var <- c(1,2)

    mfrow = c(7,3)
    byrow = FALSE
    xfixed = TRUE
    yfixed = FALSE
    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed && length(var) > 1){
        par(mar = c(0.5,6,0.5,5), oma = c(4.2,0,1,0))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,2,0.5,1), oma = c(4.2,2.6,2.5,0.6)) ##  oma = c(5,3.5,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        par(mar = c(3,3,2,2), oma = c(3,2,2,4))
    }

    ## r
    vars <- var2vars(c("r","Fmsy"))
    add.ices = TRUE
    xlim = NULL
    ylim = NULL
    xlab = "Time"
    ylab = NULL
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = FALSE
    legend.extra = FALSE
    rel = 1/nPar
    plot.legend = FALSE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = if(any(is.null(ylab))) vars$lab else ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[1]])
    ## K
    vars <- var2vars(c("K","Bmsy"))
    add.ices = TRUE
    xlim = NULL
    ylim = NULL
    xlab = ""
    ylab = NULL
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = 2
    legend.extra = FALSE
    rel = c((1/nPar)^(1/(nPar-1)))
    plot.legend = FALSE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = if(any(is.null(ylab))) vars$lab else ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[2]])
    ## m
    vars <- var2vars(c("m","MSY"))
    add.ices = TRUE
    xlim = NULL
    ylim = NULL
    xlab = ""
    ylab = NULL
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = 2
    rel = 1
    byrow = FALSE
    legend.extra = FALSE
    plot.legend = TRUE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = if(any(is.null(ylab))) vars$lab else ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[3]])

}

plot.time.all2 <- function(resList, stockData, stockInfo, selMods = NULL,
                          maxYear = 2021, nPar = 2,
                          adj = list(c(0,0),c(0,0),c(0,0)),
                          bestMods = NULL, add.ices = TRUE,
                          leg.text = NULL, drop.mods = NULL,
                          remove.scenario = TRUE, gap = TRUE){

    var <- c(1,2)

    if(gap){
        mfrow = c(8,3)
    }else{
        mfrow = c(7,3)
    }
    byrow = FALSE
    xfixed = TRUE
    yfixed = FALSE
    if(gap){
        layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow),
               heights = c(rep(1,4),0.15,rep(1,3)))
    }else{
        layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    }
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed && length(var) > 1){
        par(mar = c(0.5,3,0.5,2), oma = c(4.2,2,4,3))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,2,0.5,1), oma = c(4.2,2.6,2.5,0.6)) ##  oma = c(5,3.5,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        par(mar = c(3,3,2,2), oma = c(3,2,2,4))
    }

    if(gap) par(oma = c(4.2,2,4,6))

    ## r
    vars <- var2vars(c("r","Fmsy"))
    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = "Time"
    ylab = ""
    title = list(parse(text = paste0(vars$lab[[1]], '* " &"')),
                 vars$lab[[2]])
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = FALSE
    legend.extra = FALSE
    rel = 1/nPar
    plot.legend = FALSE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = if(any(is.null(ylab))) vars$lab else ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[1]],
                     title = title,
                     stock.in.plot = FALSE,
                     bestMods = bestMods,
                     leg.text = leg.text,
                     drop.mods = drop.mods, gap = gap)
    ## K
    vars <- var2vars(c("K","Bmsy"))
    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = ""
    ylab = NULL
    title = list(parse(text = paste0(vars$lab[[1]], '* " &"')),
                 vars$lab[[2]])
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = remove.scenario
    legend.extra = FALSE
    rel = c((1/nPar)^(1/(nPar-1)))
    plot.legend = FALSE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[2]],
                     title = title,
                     stock.in.plot = FALSE,
                     bestMods = bestMods,
                     leg.text = leg.text,
                     drop.mods = drop.mods, gap = gap)
    ## m
    vars <- var2vars(c("m","MSY"))
    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = ""
    ylab = NULL
    title = list(parse(text = paste0(vars$lab[[1]], '* " &"')),
                 vars$lab[[2]])
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = remove.scenario
    rel = 1
    byrow = FALSE
    legend.extra = FALSE
    plot.legend = TRUE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[3]],
                     title = title,
                     stock.in.plot = 2,
                     stock.in.plot.line = 3,
                     bestMods = bestMods,
                     leg.text = leg.text,
                     drop.mods = drop.mods, gap = gap)

    ti <- 2
    par(xpd = NA)
    y_top <- par("usr")[4] + 21  ## top
    y_bottom <- par("usr")[3] + 3 ## bottom
    x_pos <- par("usr")[2] + 18    ## right
    segments(x_pos, y_bottom, x_pos, y_top)
    segments(x_pos - ti, y_top, x_pos, y_top)
    segments(x_pos - ti, y_bottom, x_pos, y_bottom)
    segments(x_pos, y_bottom + (y_top - y_bottom)/2,
             x_pos + ti, y_bottom + (y_top - y_bottom)/2)

    mtext("Flatfish", 4, 3.8, outer = TRUE, font = 2, adj = 0.195)


    y_top <- par("usr")[4] + 71  ## top
    y_bottom <- par("usr")[3] + 40 ## bottom
    x_pos <- par("usr")[2] + 18    ## right
    segments(x_pos, y_bottom, x_pos, y_top)
    segments(x_pos - ti, y_top, x_pos, y_top)
    segments(x_pos - ti, y_bottom, x_pos, y_bottom)
    segments(x_pos, y_bottom + (y_top - y_bottom)/2,
             x_pos + ti, y_bottom + (y_top - y_bottom)/2)

    mtext("Roundfish", 4, 3.8, outer = TRUE, font = 2, adj = 0.735)

    mtext("Respective value", 2, 0, outer = TRUE, adj = 0.176)
    mtext("Respective value", 2, 0, outer = TRUE, adj = 0.736)

}


plot.time.all.noref <- function(resList, stockData, stockInfo, selMods = NULL,
                          maxYear = 2021, nPar = 2,
                          adj = list(c(0,0),c(0,0),c(0,0)),
                          bestMods = NULL, add.ices = TRUE){

    var <- c(1)

    mfrow = c(7,3)
    byrow = FALSE
    xfixed = TRUE
    yfixed = FALSE
    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed && length(var) > 1){
        par(mar = c(0.5,3,0.5,2), oma = c(4.2,2,4,3))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,3,0.5,0), oma = c(4.2,2.6,3,3)) ##  oma = c(5,3.5,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        par(mar = c(3,3,2,2), oma = c(3,2,2,4))
    }

    ## r
    vars <- var2vars(c("r"))
    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = "Time"
    ylab = "Respective value"
    title = list(parse(text = paste0(vars$lab[[1]])))
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = FALSE
    legend.extra = FALSE
    rel = 1
    plot.legend = FALSE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = if(any(is.null(ylab))) vars$lab else ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[1]],
                     title = title,
                     title.line = 1,
                     stock.in.plot = FALSE,
                     bestMods = bestMods)
    ## K
    vars <- var2vars(c("K"))
    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = ""
    ylab = NULL
    title = list(parse(text = paste0(vars$lab[[1]])))
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = FALSE
    legend.extra = FALSE
    rel = 1
    plot.legend = FALSE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[2]],
                     title = title,
                     title.line = 1,
                     stock.in.plot = FALSE,
                     bestMods = bestMods)
    ## m
    vars <- var2vars(c("m"))
    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = ""
    ylab = NULL
    title = list(parse(text = paste0(vars$lab[[1]])))
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = FALSE
    rel = 1
    byrow = FALSE
    legend.extra = FALSE
    plot.legend = TRUE
    plot.time.single2(resList, stockData, stockInfo, selMods = selMods,
                     maxYear = maxYear, var = vars$var, scale = vars$scale,
                     rel = rel,
                     add.ices = add.ices, xlim = xlim,
                     ylim = ylim, xlab = xlab,
                     ylab = ylab,
                     xfixed = xfixed, yfixed = yfixed,
                     plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                     plot.ref = plot.ref,
                     remove.cm.roundfish = remove.cm.roundfish,
                     byrow = byrow,
                     legend.extra = legend.extra,
                     plot.legend = plot.legend,
                     adj = adj[[3]],
                     title = title,
                     title.line = 1,
                     stock.in.plot = 2,
                     stock.in.plot.line = 1,
                     bestMods = bestMods)

}

## get correlations
get.cors <- function(resList, stockData, stockInfo, selMods = NULL,
                     envData = NULL, maxYear = 2021, var = NULL,
                     mod.link = "log", fit.poly = TRUE, use.aicc = FALSE,
                     lag = 0, pres = NULL){

    if(!is.null(selMods)){
        presListSel <- resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
            if(!is.null(pres)){
                presListSel[[i]] <- pres[[i]][indi]
            }
        }
    }else{
        resListSel <- resList
        indiList <- lapply(resList, function(x) seq(length(x)))
        presListSel <- pres
    }


    vars <- var2vars(var)
    nvar <- length(vars$var)
    expl.is.env <- length(grep("temp", vars$var[2])) > 0 || length(grep("sal", vars$var[2])) > 0 || length(grep("oxy", vars$var[2])) > 0

    if(length(lag) == 1) lag <- rep(lag, length(resListSel))

    res <- vector("list", length(resListSel))
    for(i in 1:length(resListSel)){
        ## Response variable
        if(all(!is.na(resListSel[[i]][[1]]$logB)) ||
           !is.null(pres)){

            if(!is.null(pres)){
                yearsResp <- presListSel[[i]][[1]]$time[presListSel[[i]][[1]]$time <= maxYear]
                resp <- presListSel[[i]][[1]]$B[presListSel[[i]][[1]]$time <= maxYear]
            }else{
                tmp <- spict::annual(as.numeric(rownames(resListSel[[i]][[1]]$logB)),
                                     resListSel[[i]][[1]][[vars$var[1]]][,2], mean)
                yearsResp <- tmp$anntime[tmp$anntime <= maxYear]
                resp <- tmp$annvec[tmp$anntime <= maxYear]
            }

            ## Explanatory variable
            expl <- vector("list", nvar-1)
            if(expl.is.env){
                yearsExpl <- envData$year
                for(j in 2:nvar){
                    expl[[j-1]] <- envData[yearsExpl <= (max(yearsResp)-lag[i]) &
                                           yearsExpl >= (min(yearsResp)-lag[i]),
                                           vars$var[j]]
                }
            }else{
                yearsExpl <- stockData[[i]]$data$year
                expl <- list()
                if(vars$var[2] == "rec"){
                    expl[[1]] <- stockData[[i]]$data$rec[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                    vars$lab[[2]] <- var2vars("rec")$lab[[1]]
                }else if(vars$var[2] == "recSSB"){
                    expl[[1]] <- stockData[[i]]$data$recSSB[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                    vars$lab[[2]] <- var2vars("recSSB")$lab[[1]]
                }else if(vars$var[2] == "SSB"){
                    expl[[1]] <- stockData[[i]]$data$SSB[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                    vars$lab[[2]] <- var2vars("SSB")$lab[[1]]
                }else if(vars$var[2] == "natMort" || (vars$var[2] == "comb" && i %in% c(2))){
                    expl[[1]] <- rowMeans(stockData[[i]]$natMort[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                    vars$lab[[2]] <- var2vars("natMort")$lab[[1]]
                }else if(vars$var[2] == "mat" || (vars$var[2] == "comb" && i %in% c(100))){
                    expl[[1]] <- rowMeans(stockData[[i]]$mat[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                    vars$lab[[2]] <- var2vars("mat")$lab[[1]]
                }else if(vars$var[2] == "ageComp" || (vars$var[2] == "comb" && i %in% c(4,6,9))){
                    expl[[1]] <- stockData[[i]]$ageComp[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),ncol(stockData[[i]]$ageComp)] / rowSums(stockData[[i]]$ageComp[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                    vars$lab[[2]] <- var2vars("ageComp")$lab[[1]]
                }else if(vars$var[2] == "selYoung" || (vars$var[2] == "comb" && i %in% c(1,7,10))){
                    expl[[1]] <- stockData[[i]]$sel[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),1]
                    vars$lab[[2]] <- var2vars("selYoung")$lab[[1]]
                }else if(vars$var[2] == "selOld" || (vars$var[2] == "comb" && i %in% c(100))){
                    expl[[1]] <- stockData[[i]]$sel[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),ncol(stockData[[i]]$sel)]
                    vars$lab[[2]] <- var2vars("selOld")$lab[[1]]
                }else if(vars$var[2] == "wCatch" || (vars$var[2] == "comb" && i %in% c(5,8))){
                    expl[[1]] <- rowMeans(stockData[[i]]$wCatch[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                    vars$lab[[2]] <- var2vars("wCatch")$lab[[1]]
                }else if(vars$var[2] == "wStock" || (vars$var[2] == "comb" && i %in% c(3))){
                    expl[[1]] <- rowMeans(stockData[[i]]$wStock[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                    vars$lab[[2]] <- var2vars("wStock")$lab[[1]]
                }else if(vars$var[2] == "F2"){
                    expl[[1]] <- stockData[[i]]$data$Fesb[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                }else if(vars$var[2] == "FFmsy2"){
                    expl[[1]] <- stockData[[i]]$data$Fesb[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)] / stockInfo$EQ_FmsyESB[i]
                }
            }

            mod0 <- glm(resp ~ 1, family = gaussian(link = mod.link))

            if(i == 1 && vars$var.ori[2] == "natMort"){
                resp <- resp[1:48]
                expl[[1]] <- expl[[1]][1:48]
            }

            ## if(lag > 0){
            ##     resp <- resp[1:(length(resp)-lag)]
            ##     for(j in 1:length(expl)){
            ##         expl[[j]] <- expl[[j]][(lag+1):length(expl[[j]])]
            ##     }
            ## }else if(is.na(lag)){
            ##     stop("'lag' not numeric!")
            ## }


            if(all(diff(resp) < 1e-8) || all(diff(expl[[1]]) < 1e-8)){
                modSel <- mod0
            }else{
                degree <- 2
                if(length(unique(diff(expl[[1]]))) < 4) degree <- 1
                expl1 <- expl[[1]]
                if(nvar > 2){
                    expl2 <- expl[[2]]
                    if(fit.poly){
                        modFull <- glm(resp ~ poly(expl1, degree = degree) + expl2,
                                       family = gaussian(link = mod.link))
                    }else{
                        modFull <- glm(resp ~ expl1 + expl2,
                                       family = gaussian(link = mod.link))
                    }
                }else{
                    if(fit.poly){
                        modFull <- glm(resp ~ poly(expl1, degree = degree),
                                       family = gaussian(link = mod.link))
                    }else{
                        modFull <- glm(resp ~ expl1,
                                       family = gaussian(link = mod.link))
                    }
                }


                if(!use.aicc){
                    modSel <- MASS::stepAIC(modFull, direction="backward",
                                            test = "Chisq",
                                            scope=formula(mod0), trace=0)

                }else{

                    if(nvar > 2){
                        modSel1 <- modFull
                        aic1 <- spind::aic.calc(modFull, "gaussian",
                                                cbind(resp, expl1, expl2), modSel1$fitted)$AICc

                        if(fit.poly){
                            modSel2 <- glm(resp ~ expl1 + expl2,
                                           family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                    cbind(resp, expl1, expl2), modSel2$fitted)$AICc

                            if(aic2 < aic1){
                                modSel1 <- modSel2
                                aic1 <- aic2
                                modSel2 <- glm(resp ~ expl1, family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                        cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                    aic1 <- aic2
                                    modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                    aic2 <- spind::aic.calc(modSel2, "gaussian", cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                    if(aic2 < aic1){
                                        modSel1 <- modSel2
                                    }
                                }else{
                                    modSel1 <- modSel2
                                    aic1 <- aic2
                                    modSel2 <- glm(resp ~ expl2, family = gaussian(link = mod.link))
                                    aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                            cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                    if(aic2 < aic1){
                                        modSel1 <- modSel2
                                        aic1 <- aic2
                                        modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                        aic2 <- spind::aic.calc(modSel2, "gaussian", cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                        if(aic2 < aic1){
                                            modSel1 <- modSel2
                                        }
                                    }
                                }
                            }else{
                                modSel2 <- glm(resp ~ poly(expl1, degree = degree),
                                               family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                        cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                    aic1 <- aic2
                                    modSel2 <- glm(resp ~ expl1, family = gaussian(link = mod.link))
                                    aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                            cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                    if(aic2 < aic1){
                                        modSel1 <- modSel2
                                        aic1 <- aic2
                                        modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                        aic2 <- spind::aic.calc(modSel2, "gaussian", cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                        if(aic2 < aic1){
                                            modSel1 <- modSel2
                                        }
                                    }
                                }
                            }
                        }else{
                            modSel2 <- glm(resp ~ expl1,
                                           family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                    cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                            if(aic2 < aic1){
                                modSel1 <- modSel2
                                aic1 <- aic2
                                modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                        cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                }
                            }else{
                                modSel1 <- modSel2
                                aic1 <- aic2
                                modSel2 <- glm(resp ~ exp2, family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                        cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                    aic1 <- aic2
                                    modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                    aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                            cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                    if(aic2 < aic1){
                                        modSel1 <- modSel2
                                    }
                                }
                            }
                        }
                    }else{
                        modSel1 <- modFull
                        aic1 <- spind::aic.calc(modFull, "gaussian",
                                                cbind(resp, expl1), modSel1$fitted)$AICc
                        if(fit.poly){
                            modSel2 <- glm(resp ~ expl1,
                                           family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                    cbind(resp, expl1), modSel2$fitted)$AICc
                            if(aic2 < aic1){
                                modSel1 <- modSel2
                                aic1 <- aic2
                                modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                        cbind(resp, expl1), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                }
                            }else{
                                modSel2 <- glm(resp ~ 1,
                                               family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                        cbind(resp, expl1), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                }
                            }
                        }
                    }

                    modSel <- modSel1

                }

                if(i == 1 && vars$var.ori[2] == "natMort"){
                    hihi <<- predict(modSel,
                                     newdata = data.frame(expl1 = seq(min(expl1), max(expl1),
                                                                      length.out = 1e2)),
                                     type = "link", se.fit = TRUE)
                }

            }
            res[[i]] <- list(mod = modSel,
                             link = mod.link,
                             var = vars,
                             year = yearsResp,
                             resp = resp,
                             expl = expl)
        }
    }
    names(res) <- names(resList)

    return(res)
}

get.cor <- function(resList, envData,
                    maxYear = 2021,
                    env.var = "temp.bot",
                    do.log = FALSE){
    res <- vector("list", length(resList))
    for(i in 1:length(resList)){
        res[[i]] <- vector("list", length(resList[[i]]))
        for(j in 1:length(resList[[i]])){
            if(all(!is.na(resList[[i]][[j]]$logB))){
                tmp <- spict::annual(as.numeric(rownames(resList[[i]][[j]]$logB)),
                                     resList[[i]][[j]][["logmvec"]][,2], mean)
                yearsResp <- tmp$anntime[tmp$anntime <= maxYear]
                resp <- tmp$annvec[tmp$anntime <= maxYear]
                yearsExpl <- envData$year
                expl <- envData[yearsExpl <= max(yearsResp) &
                                yearsExpl >= min(yearsResp),env.var]
                if(do.log){
                    res[[i]][[j]] <- cor(log(resp), log(expl))
                }else{
                    res[[i]][[j]] <- cor(resp, expl)
                }
            }
        }
        names(res[[i]]) <- names(resList[[i]])
    }
    names(res) <- names(resList)

    all.mods <- c("CM","KM","RM","PKRM","LKRM","IKRM")

    tmpi <- lapply(res, function(x) do.call(rbind, x))
    indi <- lapply(tmpi, function(x) match(rownames(x), all.mods))

    rep <- matrix(NA, length(tmpi), length(all.mods))
    for(i in 1:length(tmpi)){
        rep[i,indi[[i]]] <- tmpi[[i]]
    }

    colnames(rep) <- all.mods
    rownames(rep) <- names(tmpi)

    return(rep)
}


## get correlations
get.cors2 <- function(resList, stockData, stockInfo, selMods = NULL,
                      envData = NULL, maxYear = 2021, var = NULL,
                      mod.link = "log", fit.poly = TRUE, use.aicc = FALSE,
                      lag = 0){

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
    }else{
        resListSel <- resList
        indiList <- lapply(resList, function(x) seq(length(x)))
    }


    vars <- var2vars(var)
    nvar <- length(vars$var)
    expl.is.env <- length(grep("temp", vars$var[2])) > 0 || length(grep("sal", vars$var[2])) > 0 || length(grep("oxy", vars$var[2])) > 0


    res <- vector("list", length(resListSel))
    for(i in 1:length(resListSel)){
    res[[i]] <- vector("list", length(resListSel[[i]]))
        for(ii in 1:length(resListSel[[i]])){
            ## Response variable
            if(all(!is.na(resList[[i]][[ii]]$logB))){
        tmp <- spict::annual(as.numeric(rownames(resListSel[[i]][[ii]]$logB)),
                             resListSel[[i]][[ii]][[vars$var[1]]][,2], mean)
        yearsResp <- tmp$anntime[tmp$anntime <= maxYear]
        resp <- tmp$annvec[tmp$anntime <= maxYear]

        ## Explanatory variable
        expl <- vector("list", nvar-1)
        if(expl.is.env){
            yearsExpl <- envData$year
            for(j in 2:nvar){
                expl[[j-1]] <- envData[yearsExpl <= max(yearsResp) &
                                       yearsExpl >= min(yearsResp),vars$var[j]]
            }
        }else{
            yearsExpl <- stockData[[i]]$data$year
            expl <- list()
            if(vars$var[2] == "rec"){
                expl[[1]] <- stockData[[i]]$data$rec[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                vars$lab[[2]] <- var2vars("rec")$lab[[1]]
            }else if(vars$var[2] == "recSSB"){
                expl[[1]] <- stockData[[i]]$data$recSSB[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                vars$lab[[2]] <- var2vars("recSSB")$lab[[1]]
            }else if(vars$var[2] == "SSB"){
                expl[[1]] <- stockData[[i]]$data$SSB[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp)]
                vars$lab[[2]] <- var2vars("SSB")$lab[[1]]
            }else if(vars$var[2] == "natMort" || (vars$var[2] == "comb" && i %in% c(2))){
                expl[[1]] <- rowMeans(stockData[[i]]$natMort[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                vars$lab[[2]] <- var2vars("natMort")$lab[[1]]
            }else if(vars$var[2] == "mat" || (vars$var[2] == "comb" && i %in% c(100))){
                expl[[1]] <- rowMeans(stockData[[i]]$mat[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                vars$lab[[2]] <- var2vars("mat")$lab[[1]]
            }else if(vars$var[2] == "ageComp" || (vars$var[2] == "comb" && i %in% c(4,6,9))){
                expl[[1]] <- stockData[[i]]$ageComp[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),ncol(stockData[[i]]$ageComp)] / rowSums(stockData[[i]]$ageComp[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                vars$lab[[2]] <- var2vars("ageComp")$lab[[1]]
            }else if(vars$var[2] == "selYoung" || (vars$var[2] == "comb" && i %in% c(1,7,10))){
                expl[[1]] <- stockData[[i]]$sel[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),1]
                vars$lab[[2]] <- var2vars("selYoung")$lab[[1]]
            }else if(vars$var[2] == "selOld" || (vars$var[2] == "comb" && i %in% c(100))){
                expl[[1]] <- stockData[[i]]$sel[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),ncol(stockData[[i]]$sel)]
                vars$lab[[2]] <- var2vars("selOld")$lab[[1]]
            }else if(vars$var[2] == "wCatch" || (vars$var[2] == "comb" && i %in% c(5,8))){
                expl[[1]] <- rowMeans(stockData[[i]]$wCatch[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                vars$lab[[2]] <- var2vars("wCatch")$lab[[1]]
            }else if(vars$var[2] == "wStock" || (vars$var[2] == "comb" && i %in% c(3))){
                expl[[1]] <- rowMeans(stockData[[i]]$wStock[yearsExpl <= max(yearsResp) & yearsExpl >= min(yearsResp),])
                vars$lab[[2]] <- var2vars("wStock")$lab[[1]]
            }
        }

        mod0 <- glm(resp ~ 1, family = gaussian(link = mod.link))

        if(i == 1 && vars$var.ori[2] == "natMort"){
            resp <- resp[1:48]
            expl[[1]] <- expl[[1]][1:48]
        }

        if(lag > 0){
            resp <- resp[1:(length(resp)-lag)]
            for(j in 1:length(expl)){
                expl[[j]] <- expl[[j]][(lag+1):length(expl[[j]])]
            }
        }else if(is.na(lag)){
            stop("'lag' not numeric!")
        }


        if(all(diff(resp) < 1e-8) || all(diff(expl[[1]]) < 1e-8)){
            modSel <- mod0
        }else{
            degree <- 2
            if(length(unique(diff(expl[[1]]))) < 4) degree <- 1
            expl1 <- expl[[1]]
            if(nvar > 2){
                expl2 <- expl[[2]]
                if(fit.poly){
                    modFull <- glm(resp ~ poly(expl1, degree = degree) + expl2,
                                   family = gaussian(link = mod.link))
                }else{
                    modFull <- glm(resp ~ expl1 + expl2,
                                   family = gaussian(link = mod.link))
                }
            }else{
                if(fit.poly){
                    modFull <- glm(resp ~ poly(expl1, degree = degree),
                                   family = gaussian(link = mod.link))
                }else{
                    modFull <- glm(resp ~ expl1,
                                   family = gaussian(link = mod.link))
                }
            }


            if(!use.aicc){
            modSel <- MASS::stepAIC(modFull, direction="backward", test = "Chisq",
                                    scope=formula(mod0), trace=0)

            }else{

            if(nvar > 2){
                modSel1 <- modFull
                aic1 <- spind::aic.calc(modFull, "gaussian",
                                        cbind(resp, expl1, expl2), modSel1$fitted)$AICc

                if(fit.poly){
                    modSel2 <- glm(resp ~ expl1 + expl2,
                                   family = gaussian(link = mod.link))
                    aic2 <- spind::aic.calc(modSel2, "gaussian",
                                            cbind(resp, expl1, expl2), modSel2$fitted)$AICc

                    if(aic2 < aic1){
                        modSel1 <- modSel2
                        aic1 <- aic2
                        modSel2 <- glm(resp ~ expl1, family = gaussian(link = mod.link))
                        aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                        if(aic2 < aic1){
                            modSel1 <- modSel2
                            aic1 <- aic2
                            modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian", cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                            if(aic2 < aic1){
                                modSel1 <- modSel2
                            }
                        }else{
                            modSel1 <- modSel2
                            aic1 <- aic2
                            modSel2 <- glm(resp ~ expl2, family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                    cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                            if(aic2 < aic1){
                                modSel1 <- modSel2
                                aic1 <- aic2
                                modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian", cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                }
                            }
                        }
                    }else{
                        modSel2 <- glm(resp ~ poly(expl1, degree = degree),
                                       family = gaussian(link = mod.link))
                        aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                        if(aic2 < aic1){
                            modSel1 <- modSel2
                            aic1 <- aic2
                            modSel2 <- glm(resp ~ expl1, family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                    cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                            if(aic2 < aic1){
                                modSel1 <- modSel2
                                aic1 <- aic2
                                modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                                aic2 <- spind::aic.calc(modSel2, "gaussian", cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                                if(aic2 < aic1){
                                    modSel1 <- modSel2
                                }
                            }
                        }
                    }
                }else{
                    modSel2 <- glm(resp ~ expl1,
                                   family = gaussian(link = mod.link))
                    aic2 <- spind::aic.calc(modSel2, "gaussian",
                                            cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                    if(aic2 < aic1){
                        modSel1 <- modSel2
                        aic1 <- aic2
                        modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                        aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                        if(aic2 < aic1){
                            modSel1 <- modSel2
                        }
                    }else{
                        modSel1 <- modSel2
                        aic1 <- aic2
                        modSel2 <- glm(resp ~ exp2, family = gaussian(link = mod.link))
                        aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                        if(aic2 < aic1){
                            modSel1 <- modSel2
                            aic1 <- aic2
                            modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                            aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                    cbind(resp, expl1, expl2), modSel2$fitted)$AICc
                            if(aic2 < aic1){
                                modSel1 <- modSel2
                            }
                        }
                    }
                }
            }else{
                modSel1 <- modFull
                aic1 <- spind::aic.calc(modFull, "gaussian",
                                        cbind(resp, expl1), modSel1$fitted)$AICc
                if(fit.poly){
                    modSel2 <- glm(resp ~ expl1,
                                   family = gaussian(link = mod.link))
                    aic2 <- spind::aic.calc(modSel2, "gaussian",
                                            cbind(resp, expl1), modSel2$fitted)$AICc
                    if(aic2 < aic1){
                        modSel1 <- modSel2
                        aic1 <- aic2
                        modSel2 <- glm(resp ~ 1, family = gaussian(link = mod.link))
                        aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                cbind(resp, expl1), modSel2$fitted)$AICc
                        if(aic2 < aic1){
                            modSel1 <- modSel2
                        }
                    }else{
                        modSel2 <- glm(resp ~ 1,
                                       family = gaussian(link = mod.link))
                        aic2 <- spind::aic.calc(modSel2, "gaussian",
                                                cbind(resp, expl1), modSel2$fitted)$AICc
                        if(aic2 < aic1){
                            modSel1 <- modSel2
                        }
                    }
                }
            }

            modSel <- modSel1

            }

            if(i == 1 && vars$var.ori[2] == "natMort"){
                hihi <<- predict(modSel,
                                 newdata = data.frame(expl1 = seq(min(expl1), max(expl1),
                                                                  length.out = 1e2)),
                                 type = "link", se.fit = TRUE)
            }

        }
            res[[i]][[ii]] <- list(mod = modSel,
                                   link = mod.link,
                                   var = vars,
                                   year = yearsResp,
                                   resp = resp,
                                   expl = expl)
            }
        }
    names(res[[i]]) <- names(resList[[i]])
    }
    names(res) <- names(resList)

    return(res)
}


## correlation plot
plot.cor <- function(resList, stockData, stockInfo, selMods = NULL,
                     cors, maxYear = 2021, pred.expl2 = c(34.75,34.95),
                     xlim = NULL,
                     ylim = NULL, xlab = NULL,
                     ylab = NULL, cex = 1.3,
                     yfixed = FALSE, xfixed = TRUE, plot.ci = TRUE,
                     mfrow = c(2,5), leg.ncol = 1, every.ylab = FALSE,
                     byrow = TRUE, legend.extra = FALSE){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
    }else{
        resListSel <- resList
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    is.comb <- any(unlist(lapply(cors, function(x) x$var$var)) == "comb")
    if(is.comb){
        xfixed <- FALSE
    }


    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,2,0.5,1), oma = c(5,3.5,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        if(is.comb && every.ylab){
            par(mar = c(4,1.5,2,1.5), oma = c(1,4,0.5,1))
        }else{
            par(mar = c(1.5,1.5,1.5,1.5), oma = c(4,4,1,1))
        }
    }
    for(i in 1:length(resListSel)){
        nvar <- length(cors[[i]]$var$var)

        obs.y <- cors[[i]]$resp / cors[[i]]$var$scale[1]
        obs.x <- cors[[i]]$expl[[1]]

        explPred <- list()
        explPred[[1]] <- seq(min(cors[[i]]$expl[[1]]), max(cors[[i]]$expl[[1]]),
                             length.out = 1e2)
        if(nvar > 2){
            explPred[[2]] <- seq(pred.expl2[1], pred.expl2[2], length.out = 2)
            predData <- expand.grid(expl1 = explPred[[1]],
                                    expl2 = explPred[[2]])
        }else{
            predData <- expand.grid(expl1 = explPred[[1]])
        }

        if(i == 1 && cors[[i]]$var$var.ori[2] == "natMort"){
            preds <- hihi
        }else{
            preds <- predict(cors[[i]]$mod,
                         newdata = predData,
                         type = "link", se.fit = TRUE)
        }

        critval <- 1.96 ## approx 95% CI
        if(cors[[i]]$link == "log"){
            upr <- exp(preds$fit + (critval * preds$se.fit))
            lwr <- exp(preds$fit - (critval * preds$se.fit))
            fit <- exp(preds$fit)
        }else{
            upr <- preds$fit + (critval * preds$se.fit)
            lwr <- preds$fit - (critval * preds$se.fit)
            fit <- preds$fit
        }
        newdata <- cbind(predData,
                         fit = fit / cors[[i]]$var$scale[1],
                         lwr = lwr / cors[[i]]$var$scale[1],
                         upr = upr / cors[[i]]$var$scale[1])

        if(is.null(xlim0)) xlim <- range(newdata$expl1)
        if(plot.ci){
            if(is.null(ylim0)) ylim <- range(obs.y, newdata[,c("fit","lwr","upr")])
        }else{
            if(is.null(ylim0)) ylim <- range(obs.y, newdata[,"fit"])
        }
        ## if(!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow))
        ##     xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        if(!yfixed || i %in% sapply(seq(mfrow[1]),function(x) (x-1) * (mfrow[2]) + 1))
            yaxt = "s" else yaxt = "n"
        plot(1, 1,
             ty = 'n',
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = yaxt,
             xlab = "",
             ylab = ""
             )

        nj <- as.integer(nvar>2)+1

        for(j in 1:nj){
            if(nj > 1){
                tmp <- newdata[newdata$expl2 == pred.expl2[j],]
            }else{
                tmp <- newdata
            }

            if(plot.ci){
                polygon(c(tmp$expl1,rev(tmp$expl1)),
                        c(tmp$lwr, rev(tmp$upr)),
                        border = NA,
                        col = adjustcolor(cols[indiList[[i]]], 0.2))
            }
            lines(tmp$expl1,
                  tmp$fit,
                  lwd = 2,
                  col = cols[indiList[[i]]],
                  lty = j)
        }

        points(obs.x, obs.y,
               col = adjustcolor("grey20",0.4))

        legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)

        if(i == 1 &&
           any(sapply(cors, function(x) any(names(coefficients(x$mod)) == "expl2")))){

            legend("topright",
                   title = cors[[1]]$var$lab[[3]],
                   legend = pred.expl2,
                   lwd = 2,
                   lty = 1:2,
                   col = "grey60",
                   bg = "white",
                   ncol = 1, ##leg.ncol,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   text.width = 0.3,
                   cex = 1.2)
        }
        if(i == mfrow[2] && !legend.extra){
            texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
            legend("topright", legend = texti,
                   lty = unique(unlist(lapply(indiList, function(x) ltys[x]))),
                   col = unique(unlist(lapply(indiList, function(x) cols[x]))),
                   lwd = 2,
                   bg = "white",
                   ncol = 3,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
        box(lwd = 1.5)
        if(is.comb && (i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && !every.ylab){
            mtext(cors[[i]]$var$lab[[2]],1,4)
        }else if(is.comb && every.ylab){
            mtext(cors[[i]]$var$lab[[2]],1,3.5)
        }
    }
    if(legend.extra){
        plot.new()
        texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
        legend("center", legend = texti,
               title = "Models", title.cex = 1.2, title.font = 1,
               lty = unique(unlist(lapply(indiList, function(x) ltys[x]))),
               col = unique(unlist(lapply(indiList, function(x) cols[x]))),
               lwd = 2,
               bg = "white",
               ncol = 1,
               x.intersp = 1.1,
               y.intersp = 1.2,
               cex = 1.3)
    }
    if(!is.comb){
        mtext(cors[[1]]$var$lab[[2]],1,3,outer = TRUE)
        mtext(cors[[1]]$var$lab[[1]],2,1,outer = TRUE)
    }else{
        mtext(cors[[1]]$var$lab[[1]],2,2,outer = TRUE)
    }
}


## correlation plot
plot.cor2 <- function(resList, stockData, stockInfo, selMods = NULL,
                      cors, maxYear = 2021, pred.expl2 = c(34.7,34.9),
                                                ## c(34.75,34.95),
                     xlim = NULL,
                     ylim = NULL, xlab = NULL,
                     ylab = NULL, cex = 1.3,
                     yfixed = FALSE, xfixed = TRUE, plot.ci = TRUE,
                     mfrow = c(2,5), leg.ncol = 1, every.ylab = FALSE,
                     byrow = TRUE, legend.extra = FALSE,
                     pred.all = FALSE){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- rep("black",10) ## spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
    }else{
        resListSel <- resList
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    is.comb <- any(unlist(lapply(cors, function(x) x$var$var)) == "comb")
    if(is.comb){
        xfixed <- FALSE
    }


    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    if(xfixed && yfixed){
        par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    }else if(xfixed && !yfixed){
        par(mar = c(0.5,2,0.5,1), oma = c(5,3.5,2,2))
    }else if(!xfixed && yfixed){
        par(mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    }else{
        if(is.comb && every.ylab){
            par(mar = c(4,1.5,2,1.5), oma = c(1,4,0.5,1))
        }else{
            par(mar = c(1.5,1.5,1.5,1.5), oma = c(4,4,1,1))
        }
    }
    if(pred.all){
        xlim0 <- range(unlist(lapply(cors, function(x) x$expl[[1]])))
    }
    for(i in 1:length(resListSel)){
        nvar <- length(cors[[i]]$var$var)

        obs.y <- cors[[i]]$resp / cors[[i]]$var$scale[1]
        obs.x <- cors[[i]]$expl[[1]]

        explPred <- list()
        if(!is.null(cors[[i]])){
        explPred[[1]] <- seq(min(cors[[i]]$expl[[1]]), max(cors[[i]]$expl[[1]]),
                             length.out = 1e2)
        if(nvar > 2){
            explPred[[2]] <- seq(pred.expl2[1], pred.expl2[2], length.out = 2)
            predData <- expand.grid(expl1 = explPred[[1]],
                                    expl2 = explPred[[2]])
        }else{
            predData <- expand.grid(expl1 = explPred[[1]])
        }

        if(i == 1 && cors[[i]]$var$var.ori[2] == "natMort"){
            preds <- hihi
        }else{
            preds <- predict(cors[[i]]$mod,
                         newdata = predData,
                         type = "link", se.fit = TRUE)
        }

        critval <- 1.96 ## approx 95% CI
        if(cors[[i]]$link == "log"){
            upr <- exp(preds$fit + (critval * preds$se.fit))
            lwr <- exp(preds$fit - (critval * preds$se.fit))
            fit <- exp(preds$fit)
        }else{
            upr <- preds$fit + (critval * preds$se.fit)
            lwr <- preds$fit - (critval * preds$se.fit)
            fit <- preds$fit
        }
        newdata <- cbind(predData,
                         fit = fit / cors[[i]]$var$scale[1],
                         lwr = lwr / cors[[i]]$var$scale[1],
                         upr = upr / cors[[i]]$var$scale[1])

        if(is.null(xlim0)) xlim <- range(newdata$expl1)
        if(plot.ci){
            if(is.null(ylim0)) ylim <- range(obs.y, newdata[,c("fit","lwr","upr")])
        }else{
            if(is.null(ylim0)) ylim <- range(obs.y, newdata[,"fit"])
        }
        ## if(!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow))
        ##     xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        if(!yfixed || i %in% sapply(seq(mfrow[1]),function(x) (x-1) * (mfrow[2]) + 1))
            yaxt = "s" else yaxt = "n"
        plot(1, 1,
             ty = 'n',
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = yaxt,
             xlab = "",
             ylab = ""
             )

        nj <- as.integer(nvar>2)+1

        for(j in 1:nj){
            if(nj > 1){
                tmp <- newdata[newdata$expl2 == pred.expl2[j],]
            }else{
                tmp <- newdata
            }

            if(plot.ci){
                polygon(c(tmp$expl1,rev(tmp$expl1)),
                        c(tmp$lwr, rev(tmp$upr)),
                        border = NA,
                        col = adjustcolor(cols[indiList[[i]]], 0.2))
            }
            lines(tmp$expl1,
                  tmp$fit,
                  lwd = 2,
                  col = cols[indiList[[i]]],
                  lty = j)
        }

        points(obs.x, obs.y,
               col = adjustcolor("grey20",0.4))

        legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)

        if(i == mfrow[2] && !legend.extra){
            texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
            legend("topright", legend = texti,
                   lty = unique(unlist(lapply(indiList, function(x) ltys[x]))),
                   col = unique(unlist(lapply(indiList, function(x) cols[x]))),
                   lwd = 2,
                   bg = "white",
                   ncol = 3,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
        legend("topright",
               legend = bquote(R^2 ~ "=" ~
                                   .(signif(with(summary(cors[[i]]$mod),
                                                 1 - deviance/null.deviance),2))),
               pch = NA, x.intersp = 0.1)
        box(lwd = 1.5)
        if(is.comb && (i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && !every.ylab){
            mtext(cors[[i]]$var$lab[[2]],1,4)
        }else if(is.comb && every.ylab){
            mtext(cors[[i]]$var$lab[[2]],1,3.5)
        }
        }
    }
    ## if(legend.extra){
    ##     plot.new()
    ##     texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
    ##     legend("center", legend = texti,
    ##            title = "Models", title.cex = 1.2, title.font = 1,
    ##            lty = unique(unlist(lapply(indiList, function(x) ltys[x]))),
    ##            col = unique(unlist(lapply(indiList, function(x) cols[x]))),
    ##            lwd = 2,
    ##            bg = "white",
    ##            ncol = 1,
    ##            x.intersp = 1.1,
    ##            y.intersp = 1.2,
    ##            cex = 1.3)
    ## }
       if(any(sapply(cors, function(x) any(names(coefficients(x$mod)) == "expl2")))){
    plot.new()
    legend("center",
           title = cors[[1]]$var$lab[[3]],
           legend = pred.expl2,
           lwd = 2,
           lty = 1:2,
           ## col = "grey60",
           bg = "white",
           ncol = 1, ##leg.ncol,
           x.intersp = 1.1,
           y.intersp = 1.0,
           text.width = 0.3,
           cex = 1.2)
    }
    if(!is.comb){
        mtext(cors[[1]]$var$lab[[2]],1,3,outer = TRUE)
        mtext(cors[[1]]$var$lab[[1]],2,1,outer = TRUE)
    }else{
        mtext(cors[[1]]$var$lab[[1]],2,2,outer = TRUE)
    }
}



plot.schematic.maccall <- function(r = 0.2, K = 1, n = 2,
                                   nrep = 4,
                                   a = 5, b1 = 3, b2 = 0.7,
                                   xlim = NULL, ylim1 = NULL,
                                   ylim2 = NULL,
                                   xlab = expression("Density"),
                                   ylab1 = expression("Surplus production (" * dB%.%dt^{-1} * ")"),
                                   ylab2 = expression("Per capita growth rate (" * dB%.%dt^{-1}%.%B^{-1} * ")")
                                   ){

    ## ylab1 <-  expression(dB%.%dt^{-1}~'[wt'%.%yr^{-1}*"]")

    xlim0 <- xlim
    ylim0 <- range(ylim1,ylim2)

    cols <- lapply(spmprod.cols()[2:6], function(x)
        sapply(seq(0.7,0.1,length.out = nrep),
               function(y) adjustcolor(x, offset = rep(y,4))))
    lwds <- rep(2, nrep)
    ltys <- rep(1, nrep)

    calcProd <- function(r = 0.2, K = 1, n = 2.0){
        df <- data.frame(t = seq(0,200,0.1), B = 0.01,
                         dBdt = NaN, dBdt_B = NaN)
        for(i in seq(nrow(df))[-1]){
            Bt <- df$B[i-1]
            dt <- df$t[i] - df$t[i-1]
            dBdt <- r/(n-1)*Bt*(1-(Bt/K)^(n-1))
            df$dBdt[i-1] <- dBdt
            df$B[i] <- Bt + dBdt * dt
        }
        df$dBdt_B <- df$dBdt / df$B
        return(df)
    }

    calcProd2 <- function(r = 0.2, K = 1, n = 2.0, m = NULL){
        if(is.null(r)){
            r <- m / K * n^(n/(n-1))
        }else if(is.null(K)){
            K <- m / r * n^(n/(n-1))
        }else if(is.null(m)){
            m <- K * r / n^(n/(n-1))
        }
        print(m)
        gamma <- n^(n/(n-1)) / (n-1)
        B <- seq(1, K, length.out = 2e2)
        P <- gamma * m * (B/K) - gamma * m * (B/K)^n
        df <- data.frame(B = B, dBdt = P, dBdt_B = P/B)
        return(df)
    }

    op <- par(no.readonly = T)
    ## par(mar = c(3.5,3.5,0.5,0.5), oma = c(0,0,2,0),
    ##     mfcol = c(2,5), mgp = c(2,0.5,0), ps = 9)
    par(mfcol = c(2,5), mar = c(0.5,0.5,0.5,0.5), oma = c(4.5,4.5,2,1))

    # KM
    vars <- seq(K - 0.2, K + 0.2, length.out = nrep)
    scs <- list()
    for(i in 1:nrep){
        scs[[i]] <- calcProd(r = r, K = vars[i], n = n)
    }
    ## B ~ dB/dt
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim1 <- range(unlist(lapply(scs, function(x) x$dBdt)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim1,
         xaxs = "i", yaxs = "i",
         xaxt = "n",
         col = cols[[1]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep) lines(dBdt ~ B, scs[[i]], col = cols[[1]][i], lwd = lwds[i], lty = ltys[i])
    mtext("Varying K, fixed r (KM)", side = 3, line = 0.25, adj = 0, cex = 0.8)
    mtext(ylab1, 2, 3)
    box(lwd = 1.5)
    ## B ~ dB/dt/B
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim2 <- range(unlist(lapply(scs, function(x) x$dBdt_B)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt_B ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim2, xaxs = "i", yaxs = "i",
         col = cols[[1]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep){
        lines(dBdt_B ~ B, scs[[i]],
              col = cols[[1]][i], lwd = lwds[i], lty = ltys[i])
        points(dBdt_B ~ B, tail(scs[[i]], 2)[1,], pch = 16, col = cols[[1]][1])
        tmp <- tail(scs[[i]], 2)[1,]
        text(tmp$B,tmp$dBdt_B,
             labels = bquote(K[.(i)]), pos = 3, col = cols[[1]][1])
    }
    points(dBdt_B ~ B, head(scs[[1]], 1), pch = 16, col = cols[[1]][1])
    text(dBdt_B ~ B, head(scs[[1]], 1),
         label = expression(r), pos = 4, col = cols[[1]][1])
    mtext(ylab2, 2, 3)
    box(lwd = 1.5)

    # RM
    vars <- seq(r - 0.05, r + 0.05, length.out = nrep)
    scs <- list()
    for(i in 1:nrep){
        scs[[i]] <- calcProd(r = vars[i], K = K, n = n)
    }
    ## B ~ dB/dt
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim1 <- range(unlist(lapply(scs, function(x) x$dBdt)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim1, xaxs = "i", yaxs = "i",
         xaxt = "n", yaxt = "n",
         col = cols[[2]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep) lines(dBdt ~ B, scs[[i]], col = cols[[2]][i], lwd = lwds[i], lty = ltys[i])
    mtext("Varying r, fixed K (RM)", side = 3, line = 0.25, adj = 0, cex = 0.8)
    box(lwd = 1.5)
    ## B ~ dB/dt/B
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim2 <- range(unlist(lapply(scs, function(x) x$dBdt_B)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt_B ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim2, xaxs = "i", yaxs = "i",
         yaxt = "n",
         col = cols[[2]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep){
        lines(dBdt_B ~ B, scs[[i]],
              col = cols[[2]][i], lwd = lwds[i], lty = ltys[i])
        points(dBdt_B ~ B, head(scs[[i]], 1), pch = 16, col = cols[[2]][1])
        tmp <- head(scs[[i]], 1)
        text(tmp$B, tmp$dBdt_B,
             label = bquote(r[.(i)]), pos = 4, col = cols[[2]][1])
    }
    points(dBdt_B ~ B, tail(scs[[1]], 2)[1,], pch = 16, col = cols[[2]][1])
    text(dBdt_B ~ B, tail(scs[[1]], 2)[1,],
         label = expression(K), pos = 3, col = cols[[2]][1])
    box(lwd = 1.5)

    # Proportional (PKRM)
    vars <- seq(r-0.05,r+0.05, length.out = nrep)
    scs <- list()
    for(i in 1:nrep){
        scs[[i]] <- calcProd(r = vars[i], K = exp(log(vars[i]) + log(a)), n = n)
    }
    ## B ~ dB/dt
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim1 <- range(unlist(lapply(scs, function(x) x$dBdt)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim1, xaxs = "i", yaxs = "i",
         xaxt = "n", yaxt = "n",
         col = cols[[3]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep) lines(dBdt ~ B, scs[[i]], col = cols[[3]][i], lwd = lwds[i], lty = ltys[i])
    mtext("Varying K and r (PKRM/LKRM/IKRM)", side = 3, line = 0.25, adj = 0, cex = 0.8)
    box(lwd = 1.5)
    ## B ~ dB/dt/B
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim2 <- range(unlist(lapply(scs, function(x) x$dBdt_B)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt_B ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim2, xaxs = "i", yaxs = "i",
         yaxt = "n",
         col = cols[[3]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep){
        lines(dBdt_B ~ B, scs[[i]],
              col = cols[[3]][i], lwd = lwds[i], lty = ltys[i])
        points(dBdt_B ~ B, head(scs[[i]], 1), pch = 16, col = cols[[3]][1])
        tmp <- head(scs[[i]], 1)
        text(tmp$B, tmp$dBdt_B,
             label = bquote(r[.(i)]), pos = 4, col = cols[[3]][1])
        points(dBdt_B ~ B, tail(scs[[i]], 2)[1,], pch = 16, col = cols[[3]][1])
        tmp <- tail(scs[[i]], 2)[1,]
        text(tmp$B, tmp$dBdt_B,
             label = bquote(K[.(i)]), pos = 3, col = cols[[3]][1])
    }
    box(lwd = 1.5)

    # Linear relationship (LKRM)
    a1 <- exp(-log(b1)*log(r)) # to get coefficient value resulting in K = 1.0 for r = 0.2
    vars <- seq(r-0.05,r+0.05, length.out = nrep)
    scs <- list()
    for(i in 1:nrep){
        scs[[i]] <- calcProd(r = vars[i], K = exp(log(a1)+log(b1)*log(vars[i])), n = n)
    }
    ## B ~ dB/dt
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim1 <- range(unlist(lapply(scs, function(x) x$dBdt)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim1, xaxs = "i", yaxs = "i",
         xaxt = "n", yaxt = "n",
         col = cols[[4]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep) lines(dBdt ~ B, scs[[i]], col = cols[[4]][i], lwd = lwds[i], lty = ltys[i])

    mtext("Varying K and r (LKRM/IKRM)", side = 3, line = 0.25, adj = 0, cex = 0.8)
    box(lwd = 1.5)
    ## B ~ dB/dt/B
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim2 <- range(unlist(lapply(scs, function(x) x$dBdt_B)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt_B ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim2, xaxs = "i", yaxs = "i",
         yaxt = "n",
         col = cols[[4]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep){
        lines(dBdt_B ~ B, scs[[i]],
              col = cols[[4]][i], lwd = lwds[i], lty = ltys[i])
        points(dBdt_B ~ B, head(scs[[i]], 1), pch = 16, col = cols[[4]][1])
        tmp <- head(scs[[i]], 1)
        text(tmp$B, tmp$dBdt_B,
             label = bquote(r[.(i)]), pos = 4, col = cols[[4]][1])
        points(dBdt_B ~ B, tail(scs[[i]], 2)[1,], pch = 16, col = cols[[4]][1])
        tmp <- tail(scs[[i]], 2)[1,]
        text(tmp$B, tmp$dBdt_B,
             label = bquote(K[.(i)]), pos = 3, col = cols[[4]][1])
    }
    box(lwd = 1.5)


    # Independent (IKRM)
    a2 <- exp(-log(b2)*log(r)) # to get coefficient value resulting in K = 1.0 for r = 0.2
    vars <- seq(r-0.05,r+0.05, length.out = nrep)
    scs <- list()
    for(i in 1:nrep){
        scs[[i]] <- calcProd(r = vars[i], K = exp(log(a2)+log(b2)*log(vars[i])), n = n)
    }
    ## B ~ dB/dt
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim1 <- range(unlist(lapply(scs, function(x) x$dBdt)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim1, xaxs = "i", yaxs = "i",
         xaxt = "n", yaxt = "n",
         col = cols[[5]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep) lines(dBdt ~ B, scs[[i]], col = cols[[5]][i], lwd = lwds[i], lty = ltys[i])
    mtext("Varying K and r (LKRM/IKRM)", side = 3, line = 0.25, adj = 0, cex = 0.8)
    box(lwd = 1.5)
    ## B ~ dB/dt/B
    if(is.null(xlim0)) xlim <- range(unlist(lapply(scs, function(x) x$B)), na.rm = TRUE) * c(0, 1.05)
    if(is.null(ylim0)) ylim2 <- range(unlist(lapply(scs, function(x) x$dBdt_B)), na.rm = TRUE) * c(0, 1.05)
    plot(dBdt_B ~ B, scs[[i]], t="n",
         xlim = xlim, ylim = ylim2, xaxs = "i", yaxs = "i",
         yaxt = "n",
         col = cols[[5]][1], lwd = lwds[1], lty = ltys[1],
         xlab = "", ylab = "")
    for(i in 1:nrep){
        lines(dBdt_B ~ B, scs[[i]],
              col = cols[[5]][i], lwd = lwds[i], lty = ltys[i])
        points(dBdt_B ~ B, head(scs[[i]], 1), pch = 16, col = cols[[5]][1])
        tmp <- head(scs[[i]], 1)
        text(tmp$B, tmp$dBdt_B,
             label = bquote(r[.(i)]), pos = 4, col = cols[[5]][1])
        points(dBdt_B ~ B, tail(scs[[i]], 2)[1,], pch = 16, col = cols[[5]][1])
        tmp <- tail(scs[[i]], 2)[1,]
        text(tmp$B, tmp$dBdt_B,
             labels = bquote(K[.(i)]), pos = 3, col = cols[[5]][1])
    }
    box(lwd = 1.5)


    mtext(xlab, 1, 2.5, outer = TRUE)

    par(op)
}


## Thorson's RAM n prior
add.thorson.prior <- function(inp, order){
    if(order=="Pleuronectiformes") { n.est <- 1.353; sdn <- 0.739 }
    else if(order=="Gadiformes") { n.est <- 1.729; sdn <- 0.937 }
    else if(order=="Perciformes") { n.est <- 1.064; sdn <- 0.590 }
    else if(order=="Clupeiformes") { n.est <- 0.599; sdn <- 0.342 }
    else if(order=="Scorpaeniformes") { n.est <- 1.970; sdn <- 1.074 }
    else if(order=="Other") { n.est <- 1.431; sdn <- 0.805 }
    else stop("Order not found")
    logn.est <- log(n.est)
    var.logn <- sdn^2 / n.est^2
    inp$priors$logn <- c(logn.est, sqrt(var.logn), 1)
    return(inp)
}

## FishLife r prior
get.fishlife.prior <- function(species){
    tmp <- strsplit(species,"_")[[1]]
    FLpars <- R.devices::suppressGraphics({
        FishLife::Plot_taxa(FishLife::Search_species(Genus = tmp[1],
                                                     Species = tmp[2],
                                                     add_ancestors =
                                                         FALSE)$match_taxonomy)
    })
    logmu <- FLpars[[1]]$Mean_pred[names(FLpars[[1]]$Mean_pred) == "ln_r"]
    sd <- diag(FLpars[[1]]$Cov_pred)[names(diag(FLpars[[1]]$Cov_pred)) == "ln_r"]
    return(c(logmu,sd))
}

add.fishlife.prior <- function(inp, species){
    tmp <- get.fishlife.prior(species)
    logmu <- tmp[1]
    sd <- tmp[2]
    inp$priors$logr <- c(logmu - 0.5 * sd^2, sd, 1)
    return(inp)
}

## Henning's r SPMpriors
add.winker.prior <- function(inp, genus, species, n = TRUE, r = TRUE){
    stk <- SPMpriors::flmvn_traits(Genus = genus, Species = species, Plot = FALSE)
    tmp <- SPMpriors::fl2asem(stk, mc = 1000)
    if(n) inp$priors$logn <- c(log(tmp[1,"shape"]), tmp[2,"shape"], 1)
    if(r) inp$priors$logr <- c(log(tmp[1,"r"]), tmp[2,"r"], 1)
    return(inp)
}


## get production from spict fit
get.prod <- function(rep, n.plotyears=40, main='Production curve',
                     stamp=get.version(), CI = 0.95){
    spict:::check.rep(rep)
    res <- list()
    if (!'sderr' %in% names(rep)){
        inp <- rep$inp
        tvgflag <- try(rep$inp$timevaryinggrowth || rep$inp$logmcovflag || rep$inp$tvPropChange || rep$inp$tvNonPropChange || rep$inp$tvKConsR, silent=TRUE)
        if(inherits(tvgflag,"try-error")) tvgflag <- FALSE
        tvKflag <- try(rep$inp$timevaryingK || rep$inp$logKcovflag || rep$inp$tvKConsR || rep$inp$tvPropChange || rep$inp$tvNonPropChange, silent = TRUE)
        if(inherits(tvKflag,"try-error")) tvKflag <- FALSE
        Kest <- get.par('logK', rep, exp=TRUE, CI = CI)
        mest <- get.par('logm', rep, exp=TRUE, CI = CI)
        nr <- dim(mest)[1]
        nK <- dim(Kest)[1]
        ntv <- max(nr,nK)
        if(nr > 1){
            im <- 1:nr
        }else im <- rep(1,ntv)
        if(nK > 1){
            iK <- 1:nK
            Kmax <- max(Kest[,2])
        }else{
            iK <- rep(1,ntv)
            Kmax <- Kest[2]
        }
        gamma <- get.par('gamma', rep, CI = CI)
        n <- get.par('logn', rep, exp=TRUE, CI = CI)
        if (tvgflag){
            yscal <- get.par('logMSYvec', rep, exp=TRUE, CI = CI)[, 2]
        } else {
            yscal <- 1 ## rep(1, length(binds))
        }
        nBplot <- 200
        Bplot <- seq(0.5*1e-8, Kmax, length=nBplot)
        # Calculate production curve (Pst)
        pfun <- function(gamma, m, K, n, B) gamma*m/K*B*(1 - (B/K)^(n-1))
        Pst <- list()
        for (i in 1:ntv){
            Pst[[i]] <- pfun(gamma[2], mest[im[i],2], Kest[iK[i],2], n[2], Bplot)
        }
        if (inp$reportall){
            Best <- get.par('logB', rep, exp=TRUE, CI = CI)
            Bplot <- seq(0.5*min(c(1e-8, Best[, 2])), 1*max(c(Kest[2], Best[, 2])), length=nBplot)
            for (i in 1:ntv){
                Pst[[i]] <- pfun(gamma[2], mest[im[i],2], Kest[iK[i],2], n[2], Bplot)
            }
            Bvec <- Best[, 2]
        }
        dt <- inp$dt[-1]
        inde <- inp$indest[-length(inp$indest)]
        indp <- inp$indpred[-1]-1
        x <- Bplot/Kmax
        y <- Pst[[1]]
        yNorm <- Pst[[1]]/max(unlist(Pst))
        ## lines(Bvec/Kmax, Pest[, 2]/yscal, col=4, lwd=1.5)
        ## points(Bvec/Kmax, Pest[, 2]/yscal, col=4, pch=20, cex=0.7)
        res$x <- x
        res$y <- y
        res$yNorm <- yNorm
    }

    return(res)
}


## get production from spict fit
get.prod.year <- function(rep, n.plotyears=40, main='Production curve',
                          stamp=get.version(), CI = 0.95){
    spict:::check.rep(rep)
    res <- list()
    if (!'sderr' %in% names(rep)){
        inp <- rep$inp

        tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag || rep$inp$tvPropChange || rep$inp$tvNonPropChange || rep$inp$tvKConsR
        tvKflag <- rep$inp$timevaryingK | rep$inp$logKcovflag || rep$inp$tvKConsR || rep$inp$tvPropChange || rep$inp$tvNonPropChange
        if(tvKflag){
            Kest <- get.par('logKvec', rep, exp=TRUE, CI = CI)
        }else{
            Kest <- get.par('logK', rep, exp=TRUE, CI = CI)
        }
        if(tvgflag){
            mest <- get.par('logmvec', rep, exp=TRUE, CI = CI)
        }else{
            mest <- get.par('logm', rep, exp=TRUE, CI = CI)
        }
        Best <- get.par('logB', rep, exp=TRUE, CI = CI)
        gamma <- get.par('gamma', rep, CI = CI)
        n <- get.par('logn', rep, exp=TRUE, CI = CI)
        nBplot <- 200
        pfun <- function(gamma, m, K, n, B) gamma*m/K*B*(1 - (B/K)^(n-1))
        years <- floor(rep$inp$time[rep$inp$indest])
        uni.years <- unique(years)
        nyears <- length(uni.years)
        x <- list()
        xNorm1 <- list()
        for(i in 1:nyears){
            if(tvKflag){
                K.an <- mean(Kest[years == uni.years[i],2])
            }else{
                K.an <- Kest[,2]
            }
            Bplot <- seq(0.5*1e-8, K.an, length=nBplot)
            x[[i]] <- Bplot
            xNorm1[[i]] <- Bplot / K.an
        }
        xNorm2 <- lapply(x, function(x) x/max(unlist(x)))
        y <- list()
        yNorm1 <- list()
        for(i in 1:nyears){
            if(tvgflag){
                m.an <- mean(mest[years == uni.years[i],2])
            }else{
                m.an <- mest[,2]
            }
            if(tvKflag){
                K.an <- mean(Kest[years == uni.years[i],2])
            }else{
                K.an <- Kest[,2]
            }
            Bplot <- seq(0.5*1e-8, K.an, length=nBplot)
            y[[i]] <- pfun(gamma[2], m.an, K.an, n[2], Bplot)
            yNorm1[[i]] <- y[[i]] / m.an
        }
        yNorm2 <- lapply(y, function(x) x/max(unlist(y)))
        res$x <- x
        names(res$x) <- uni.years
        res$xNorm1 <- xNorm1
        res$xNorm2 <- xNorm2
        res$y <- y
        names(res$y) <- uni.years
        res$yNorm1 <- yNorm1
        res$yNorm2 <- yNorm2
    }

    return(res)
}


make.x <- function(x){sapply(x, function(x) ifelse(x, "x", ""))}


make.rpellipse2 <- function(rep){


    indi <- (rep$inp$dtpredcinds[1]-1-1/rep$inp$dt[1]):(rep$inp$dtpredcinds[1]-1)

    ## logFmsyvec(i) = log(mvec(i) / Kvec(i) * pow(1.0/n, 1.0/(n-1.0)));
    ## logBmsyvec(i) = log(Kvec(i) * pow(1.0/n, 1.0/(n-1.0)));
    ## logMSYvec(i) = log(mvec(i));
    logmveci <- as.numeric(rep$value[which(names(rep$value)=='logmvec')[indi]])
    ni <- get.par("logn",rep,exp = TRUE)[,2]
    logbmsyveci <-  log(exp(as.numeric(rep$value[which(names(rep$value)=='logKvec')[indi]])) * ((1/ni) ^ (1/(ni-1))))
    logfmsyveci <- log(exp(logmveci)/(exp(as.numeric(rep$value[which(names(rep$value)=='logKvec')[indi]])) * ((1/ni) ^ (1/(ni-1)))))

    ## sd the same for Kvec and Bmsyvec (for CM Bmsyvec, Fmsyvec not reported....
    sds <- c(mean(rep$sd[which(names(rep$value)=='logKvec')[indi]]), ## logBmsyvec
             mean(rep$sd[which(names(rep$value)=='logrvec')[indi]])) ## logFmsyvec
    cova <- (rep$sd[which(names(rep$value)=='logBmsyPluslogFmsy')]^2 -
             sds[1]^2 - sds[2]^2)/2
    covBF <- matrix(c(sds[1]^2, cova, cova, sds[2]^2),2,2,byrow=TRUE)
    corBF <- cov2cor(covBF)
    ## parBF <- c(mean(rep$value[which(names(rep$value)=='logBmsyvec')[indi]]),
    ##            mean(rep$value[which(names(rep$value)=='logFmsyvec')[indi]]))
    parBF <- c(mean(logbmsyveci),
               mean(logfmsyveci))
    return(ellipse::ellipse(corBF[1,2], scale=sqrt(diag(covBF)), centre=parBF, npoints=300))
}


get.indi <- function(fit, year){
    which(fit$inp$time >= year & fit$inp$time < (year+1))
}

get.mean.par <- function(par, fit, year = NULL){
    if(is.null(year)){
        mean(get.par(par, fit, exp = TRUE)[,2])
    }else{
        mean(get.par(par, fit, exp = TRUE)[get.indi(fit,year),2])
    }
}


make.cor.table <- function(cors, vars){
tmpList <- vector("list", nstocks)
for(i in 1:nstocks){
    for(j in 1:2){
        if(j == 1){
            ani <- try(summary(cors[[i]]$mod))
            if(!inherits(ani, "try-error")){
                tmp <- as.data.frame(ani$coefficients)
            }else{
                tmp <- matrix(NA,1,7)
            }
        }else{
            ani <- try(car::Anova(cors[[i]]$mod), silent = TRUE)
            if(!inherits(ani, "try-error")){
                tmp <- as.data.frame(ani)
            }else{
                tmp <- matrix(NA,1,3)
            }
        }
        tmp2 <- rownames(tmp)
        tmp2[tmp2 == "(Intercept)"] <- "Intercept"
        tmp2[tmp2 == "expl1"] <- vars[1]
        tmp2[tmp2 == "poly(expl1, degree = degree)"] <- vars[1]
        tmp2[tmp2 == "poly(expl1, degree = degree)1"] <- paste0(vars[1],"1")
        tmp2[tmp2 == "poly(expl1, degree = degree)2"] <- paste0(vars[1],"2")
        tmp2[tmp2 == "expl2"] <- vars[2]
        tmp2[tmp2 == "poly(expl1, degree = degree)1:expl2"] <- paste0(vars[1],"1:",vars[2])
        tmp2[tmp2 == "poly(expl1, degree = degree)2:expl2"] <- paste0(vars[1],"2:",vars[2])
        if(j == 1){
            tmpi <- cbind(rep(stocks[i],nrow(tmp)), tmp2, tmp)
        }else{
            if(any(rownames(tmp) == "poly(expl1, degree = degree):expl2")){
                tmpi <- cbind(tmpi, rbind(rep(NA, ncol(tmp)),
                                          rep(NA, ncol(tmp)),tmp,
                                          rep(NA, ncol(tmp))))
            }else if(any(rownames(tmp) == "poly(expl1, degree = degree)")){
                tmpi <- cbind(tmpi, rbind(rep(NA, ncol(tmp)),
                                          rep(NA, ncol(tmp)),tmp))
            }else if(length(tmp2) == 0){
                tmpi <- cbind(tmpi, tmp)
            }else{
                tmpi <- cbind(tmpi, rbind(rep(NA, ncol(tmp)),tmp))
            }
        }
    }
    colnames(tmpi) <- 1:ncol(tmpi)
    rownames(tmpi) <- NULL
    tmpList[[i]] <- tmpi
}
tmpdf <- do.call(rbind,tmpList)
tmpdf <- cbind(tmpdf[,1:2], apply(tmpdf[,c(3:9)],2, function(x) signif(x,2)))
tmpdf[which(as.numeric(tmpdf[,6]) < 0.01),6] <- "<0.01"
tmpdf[which(as.numeric(tmpdf[,9]) < 0.01),9] <- "<0.01"


tmp2 <- rbind(c("Stock","Term","Estimate","SD","tval","Pr(>|t|)","LR(Chisq)","Df","p(Chisq)"),tmpdf)
rownames(tmp2) <- NULL
colnames(tmp2) <- NULL


tmp2[is.na(tmp2)] <- ""


return(tmp2)

}



## hcr plot
plot.hcr <- function(resList, stockData, stockInfo, selMods = NULL,
                     selYear = 2021, add.ices = TRUE, xlim = NULL,
                     ylim = NULL, xlab = expression("Biomass ['000t]"),
                     ylab = expression("Fishing mortality ["*yr^-1*"]"),
                     cex = 1.5,
                     yfixed = TRUE, xfixed = TRUE, mfrow = c(2,5),
                     age.ref.period = FALSE, xscale = 1e3,
                     blim.fac = 0.3, btrigger.fac = 0.5,
                     remove.cm.roundfish = FALSE,
                     byrow = FALSE){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
        names(resListSel) <- names(resList)
    }else{
        resListSel <- resList
        names(resListSel) <- names(resList)
        indiList <- lapply(resList, function(x) seq(length(x)))
    }


    ## if(xfixed && yfixed){
    ##     par(mfrow = mfrow, mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,2,2))
    ## }else if(xfixed && !yfixed){
    ##     par(mfrow = mfrow, mar = c(0.5,2,0.5,1), oma = c(5,2,2,2))
    ## }else if(!xfixed && yfixed){
    ##     par(mfrow = mfrow, mar = c(2,0.5,1,0.5), oma = c(3,5,2,2))
    ## }else{
        ## par(mfrow = mfrow, mar = c(3,3,2,2), oma = c(3,2,2,4))
    ## }
    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    par(mar = c(2,2,1,1), oma = c(3,3.5,2,1))
    for(i in 1:length(resListSel)){
        if(is.null(xlim0)){
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               as.integer(remove.cm.roundfish) %in% c(1,2)){
                xlim <- range(0,sapply(resListSel[[i]][names(resListSel[[i]]) != "CM"],
                                       function(x){
                    years <- as.numeric(rownames(x$logB))
                    if(age.ref.period){
                        indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
                    }else{
                        indi <- years >= selYear & years < selYear + 1
                    }
                    if(all(is.na(x$logBmsyvec))){
                        x$logBmsy[,2]
                    }else{
                        mean(x$logBmsyvec[indi,2])
                    }
                }), na.rm = TRUE)
            }else{
                xlim <- range(0,sapply(resListSel[[i]], function(x){
                    years <- as.numeric(rownames(x$logB))
                    if(age.ref.period){
                        indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
                    }else{
                        indi <- years >= selYear & years < selYear + 1
                    }
                    if(all(is.na(x$logBmsyvec))){
                        x$logBmsy[,2]
                    }else{
                        mean(x$logBmsyvec[indi,2])
                    }
                }), na.rm = TRUE)
            }
            if(add.ices){
                xlim <- range(xlim, stockInfo[i,"EQ_Bmsy"],
                              stockInfo[i,"EQ_ESBmsy"])
            }
            xlim <- xlim / xscale
        }
        if(is.null(ylim0)){
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               as.integer(remove.cm.roundfish) %in% c(1,2)){
                ylim <- range(0,sapply(resListSel[[i]][names(resListSel[[i]]) != "CM"],
                                       function(x){
                    years <- as.numeric(rownames(x$logB))
                    if(age.ref.period){
                        indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
                    }else{
                        indi <- years >= selYear & years < selYear + 1
                    }
                    if(all(is.na(x$logFmsyvec))){
                        x$logFmsy[,2]
                    }else{
                        mean(x$logFmsyvec[indi,2])
                    }
                }), na.rm = TRUE)
            }else{
            ylim <- range(0,sapply(resListSel[[i]], function(x){
                years <- as.numeric(rownames(x$logB))
                if(age.ref.period){
                    indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
                }else{
                    indi <- years >= selYear & years < selYear + 1
                }
                if(all(is.na(x$logFmsyvec))){
                    x$logFmsy[,2]
                }else{
                    mean(x$logFmsyvec[indi,2])
                }
            }), na.rm = TRUE)
            }
            if(add.ices){
                ylim <- range(ylim, stockInfo[i,"Fmsy"],
                              stockInfo[i,"EQ_FmsyESB"])
            }
            ylim <- c(0,1.3) * ylim
        }
        ## if(!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow))
        ##     xaxt = "s" else xaxt = "n"
        ## if(!yfixed || i %in% sapply(seq(mfrow[1]),function(x) (x-1) * (mfrow[2]) + 1))
        ##     yaxt = "s" else yaxt = "n"
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        if(i == 2) ylim <- c(0,1.3) * ylim
        plot(1, 1,
             ty = 'n',
             xlim = xlim,
             ylim = ylim,
             xlab = "",
             ylab = ""
             )
        ## abline(h = 1, col = "grey80")
        ## abline(v = 1, col = "grey80")
        for(j in 1:length(resListSel[[i]])){
            years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            if(age.ref.period){
                indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
            }else{
                indi <- years >= selYear & years < selYear + 1
            }
            if(all(is.na(resListSel[[i]][[j]]$logBmsyvec))){
                ## xdat <- resListSel[[i]][[j]]$logBmsyvec[,2]
                xdat <- resListSel[[i]][[j]]$logBmsy[,2]
            }else{
                xdat <- mean(resListSel[[i]][[j]]$logBmsyvec[indi,2])
            }
            if(all(is.na(resListSel[[i]][[j]]$logFmsyvec))){
                ## ydat <- resListSel[[i]][[j]]$logFmsy[,2]
                ydat <- resListSel[[i]][[j]]$logFmsy[,2]
            }else{
                ydat <- mean(resListSel[[i]][[j]]$logFmsyvec[indi,2])
            }

            xdat <- xdat / xscale

            lty <- ifelse(names(resListSel[[i]])[j] == "IKRM", 2,1)

            segments(0,0, xdat * blim.fac, 0,
                     col = cols[indiList[[i]][j]],
                     lty = lty,
                     lwd = 1.5)
            segments(xdat * blim.fac, 0,
                     xdat * btrigger.fac, ydat,
                     col = cols[indiList[[i]][j]],
                     lty = lty,
                     lwd = 1.5)
            segments(xdat * btrigger.fac, ydat,
                     ifelse(xlim[2] > xdat * btrigger.fac, xlim[2],
                            xdat * btrigger.fac), ydat,
                     col = cols[indiList[[i]][j]],
                     lty = lty,
                     lwd = 1.5)
        }
        if(add.ices){
            years <- stockData[[i]]$data$year
            if(age.ref.period){
                indi <- years >= stockInfo[i,"EQ_start"] & years < stockInfo[i,"EQ_end"] + 1
            }else{
                indi <- which.min((years - selYear)^2)
            }
            xdat <- stockInfo[i,"Btrigger"]
            xdat <- xdat / xscale
            ydat <- stockInfo[i,"Fmsy"]
            segments(0,0, xdat * blim.fac, 0,
                     col = 1,
                     lwd = 1.5, lty = 2)
            segments(xdat * blim.fac, 0,
                     xdat, ydat,
                     col = 1,
                     lwd = 1.5, lty = 2)
            segments(xdat, ydat,
                     xlim[2], ydat,
                     col = 1,
                     lwd = 1.5, lty = 2)
        }
        legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)
        if(i == mfrow[2]){
            leg <- unique(unlist(lapply(resListSel, function(x) names(x))))
            colleg <- unique(unlist(lapply(indiList, function(x) cols[x])))
            lty <- unique(unlist(lapply(resListSel, function(x) names(x))))
            lty[lty != "IKRM"] <- 1
            lty[lty == "IKRM"] <- 2
            lty <- as.numeric(lty)
            if(add.ices){
                leg <- c(leg, "APM")
                lty <- c(lty, 2)
                colleg <- c(colleg, "black")
            }
            legend("topright", legend = leg,
                   col = colleg,
                   lwd = 1.5,
                   lty = lty,
                   bg = "white",
                   ncol = 3,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
        box(lwd = 1.5)
    }
    mtext(xlab,1,1.5,outer = TRUE)
    mtext(ylab,2,1.5,outer = TRUE)
}


get.maccall.axes <- function(fit, bcor = 1, prodcor = 1,
                             scaled = FALSE,
                             plot.obs = TRUE){

    tmp <- get.prod.year(fit)
    tmp$x <- lapply(tmp$x, function(x) x / bcor)
    tmp$y <- lapply(tmp$y, function(x) x / prodcor)

    Pest <- get.par('P', fit)[,2] / bcor
    Best <- get.par('logB', fit, exp=TRUE)[fit$inp$ic[1:length(Pest)],2] / prodcor

    if(scaled){
        xlim <- c(0,1)
        ylim1 <- c(0,1)
        ylim2 <- c(0,1)
    }else{
        xlim <- range(tmp$x)
        ylim1 <- range(unlist(tmp$y))
        ylim2 <- range(unlist(tmp$y) / unlist(tmp$x))
    }

    if(scaled){
        xplot <- lapply(tmp$x, function(x) x/max(x))
        yplot <- lapply(tmp$y, function(x) x/max(x))
        y22 <- lapply(as.list(1:length(tmp$y)),
                         function(i) tmp$y[[i]] / tmp$x[[i]])
        yplot2 <- lapply(y22, function(x) x/max(x))
        PBest <- (Pest/Best) / sapply(y22,max)
        Best <- Best / sapply(tmp$x,max)
        Pest <- Pest / sapply(tmp$y,max)
    }else{
        xplot <- tmp$x
        yplot <- tmp$y
        yplot2 <- lapply(as.list(1:length(tmp$y)),
                         function(i) tmp$y[[i]] / tmp$x[[i]])
        PBest <- Pest/Best
    }

    if(plot.obs){
        xlim <- range(c(xlim, Best))
        ylim1 <- range(c(ylim1, Pest))
        ylim2 <- range(c(ylim2, PBest))
    }


    return(list(xlim = xlim,
                ylim1 = ylim1,
                ylim2 = ylim2))

}

plotextra.maccall <- function(fit, bcor = 1, prodcor = 1,
                              add = FALSE, labels = TRUE,
                              models = NULL, axes = TRUE,
                              xlim = NULL, ylim1 = NULL, ylim2 = NULL,
                              plot.obs = TRUE, legend = TRUE,
                              lab2 = NULL, scaled = FALSE,
                              sel = c(1,2), font.models = 1){
    tmp <- get.prod.year(fit)
    tmp$x <- lapply(tmp$x, function(x) x / bcor)
    tmp$y <- lapply(tmp$y, function(x) x / prodcor)

   if(plot.obs){
        Pest <- get.par('P', fit)[,2] / bcor
        Best <- get.par('logB', fit, exp=TRUE)[fit$inp$ic[1:length(Pest)],2] / prodcor
   }


    cols <- rep(c("#f7ff58", "#f7f655", "#f8ed52", "#f8e54f", "#f8dc4c", "#f8d349",
                  "#f9ca46", "#f9c143", "#f9b940", "#f9b03d", "#faa73a", "#fa9e37",
                  "#fa9534", "#fb8d31", "#fb842e", "#fb7b2a", "#fb7227", "#fc6a24",
                  "#fc6121", "#fc581e", "#fd4f1b", "#fd4618", "#fd3e15", "#fd3512",
                  "#fe2c0f", "#fe230c", "#fe1a09", "#fe1206", "#ff0903", "#ff0000"),each = 3)


    cols <- rep(c("#ffd000", "#f7ca07", "#efc40e", "#e7be15", "#dfb81c",
                  "#d7b222", "#cfac29", "#c7a630", "#c0a037", "#b89a3e",
                  "#b09445", "#a88e4c", "#a08853", "#98825a", "#907c61",
                  "#887667", "#80706e", "#786a75", "#70647c", "#685e83",
                  "#60588a", "#585291", "#514c98", "#49469f", "#4140a6",
                  "#393aac", "#3134b3", "#292eba", "#2128c1", "#1922c8"), each = 3)

    cols <- colorRampPalette(c("dodgerblue1","grey80","darkred"))(length(tmp[[1]]))


    if(!add) par(mfrow = c(2,1), mar = c(1,2,0.5,1), oma = c(3,2,1,1))
    ## prod curve
    if(axes){
        xaxt <- "s"
        yaxt <- "s"
    }else{
        xaxt <- "n"
        yaxt <- "n"
    }
    if(is.null(xlim)){
        if(scaled){
            xlim <- c(0,1)
        }else{
            xlim <- range(tmp$x)
        }
    }else{
        xlim <- xlim / bcor
    }
    if(is.null(ylim1)){
        if(scaled){
            ylim1 <- c(0,1)
        }else{
            ylim1 <- range(unlist(tmp$y))
        }
    }else{
        ylim1 <- ylim1 / prodcor
    }
    if(is.null(ylim2)){
        if(scaled){
            ylim2 <- c(0,1)
        }else{
            ylim2 <- range(unlist(tmp$y) / unlist(tmp$x))
        }
    }
    if(scaled){
        xplot <- lapply(tmp$x, function(x) x/max(x))
        yplot <- lapply(tmp$y, function(x) x/max(x))
        y22 <- lapply(as.list(1:length(tmp$y)),
                         function(i) tmp$y[[i]] / tmp$x[[i]])
        yplot2 <- lapply(y22, function(x) x/max(x))
        PBest <- (Pest/Best) / sapply(y22,max)
        Best <- Best / sapply(tmp$x,max)
        Pest <- Pest / sapply(tmp$y,max)
    }else{
        xplot <- tmp$x
        yplot <- tmp$y
        yplot2 <- lapply(as.list(1:length(tmp$y)),
                         function(i) tmp$y[[i]] / tmp$x[[i]])
        PBest <- Pest/Best
    }
    if(plot.obs){
        xlim <- range(c(xlim, Best))
        ylim1 <- range(c(ylim1, Pest))
        ylim2 <- range(c(ylim2, PBest))
    }
    if(1 %in% sel){
    plot(xplot[[1]], yplot[[1]], ty = 'n',
         xlab = "", ylab = "",
         xaxt = ifelse(length(sel) == 1, xaxt, "n"),
         yaxt = yaxt,
         xlim = xlim,
         ylim = ylim1)
    for(i in 1:length(yplot)){
        lines(xplot[[i]], yplot[[i]],
              col = adjustcolor(cols[i], 0.3),
              lwd = 1.5)
    }
    if(plot.obs){
        points(Best, Pest, col = cols[1:length(Best)])
    }
    if(!is.null(lab2)){
        legend("topleft",
               legend = lab2,
               x.intersp = -0.5,
               cex = 0.8,
               bg = "white")
    }
    box(lwd = 1.5)
    if(labels) mtext(expression("dB" ~ dt^{-1} * " ['000 t" ~ yr^{-1}*"]"), 2, 2.5)
    if(!is.null(models)) mtext(models[1], 3, 0.5, font = font.models)
    }
    ## dBdt / B
    if(2 %in% sel){
    plot(xplot[[1]], yplot2[[1]], ty = 'n',
         yaxt = yaxt,
         xlab = "", ylab = "",
         xlim = xlim,
         ylim = ylim2)
    for(i in 1:length(yplot)){
        lines(xplot[[i]], yplot2[[i]],
              col = adjustcolor(cols[i], 0.3),
              lwd = 1.5)
    }
    if(plot.obs){
        points(Best, PBest, col = cols[1:length(Best)])
    }
    if(legend){
        indi <- seq(1,length(tmp$x), length.out = 4)
        legend("topright",
               legend = names(tmp$x)[indi],
               col = cols[indi],
               lwd = 1.5)
    }
    box(lwd = 1.5)
    if(labels) mtext(expression("dB" ~ dt^{-1} ~ B^{-1}), 2, 2.5)
    if(labels) mtext("Biomass ['000 t]", 1, 1.5, outer = TRUE)
    }


}



plotextra.maccall.prod <- function(fit, bcor = 1, prodcor = 1,
                                   years = 1963:2022,
                              add = FALSE, labels = TRUE,
                              models = NULL, axes = TRUE,
                              xlim = NULL, ylim1 = NULL, ylim2 = NULL,
                              plot.obs = TRUE, legend = TRUE,
                              lab2 = NULL, scaled = FALSE,
                              black.line = FALSE,
                              title = NULL, best = FALSE,
                              sel = c(1,2), font.models = 1){
    tmp <- get.prod.year(fit)
    tmp$x <- lapply(tmp$x, function(x) x / bcor)
    tmp$y <- lapply(tmp$y, function(x) x / prodcor)

   if(plot.obs){
        Pest <- get.par('P', fit)[,2] / bcor
        Best <- get.par('logB', fit, exp=TRUE)[fit$inp$ic[1:length(Pest)],2] / prodcor
   }


    cols <- rep(c("#f7ff58", "#f7f655", "#f8ed52", "#f8e54f", "#f8dc4c", "#f8d349",
                  "#f9ca46", "#f9c143", "#f9b940", "#f9b03d", "#faa73a", "#fa9e37",
                  "#fa9534", "#fb8d31", "#fb842e", "#fb7b2a", "#fb7227", "#fc6a24",
                  "#fc6121", "#fc581e", "#fd4f1b", "#fd4618", "#fd3e15", "#fd3512",
                  "#fe2c0f", "#fe230c", "#fe1a09", "#fe1206", "#ff0903", "#ff0000"),each = 3)


    cols <- rep(c("#ffd000", "#f7ca07", "#efc40e", "#e7be15", "#dfb81c",
                  "#d7b222", "#cfac29", "#c7a630", "#c0a037", "#b89a3e",
                  "#b09445", "#a88e4c", "#a08853", "#98825a", "#907c61",
                  "#887667", "#80706e", "#786a75", "#70647c", "#685e83",
                  "#60588a", "#585291", "#514c98", "#49469f", "#4140a6",
                  "#393aac", "#3134b3", "#292eba", "#2128c1", "#1922c8"), each = 3)

    cols <- colorRampPalette(c("dodgerblue1","grey80","darkred"))(length(years))

    years.local <- as.numeric(names(tmp[[1]]))
    indi.y <- match(years.local, years)


    if(!add) par(mfrow = c(2,1), mar = c(1,2,0.5,1), oma = c(3,2,1,1))
    ## prod curve
    if(axes){
        xaxt <- "s"
        yaxt <- "s"
    }else{
        xaxt <- "n"
        yaxt <- "n"
    }
    if(is.null(xlim)){
        if(scaled){
            xlim <- c(0,1)
        }else{
            xlim <- range(tmp$x)
        }
    }else{
        xlim <- xlim / bcor
    }
    if(is.null(ylim1)){
        if(scaled){
            ylim1 <- c(0,1)
        }else{
            ylim1 <- range(unlist(tmp$y))
        }
    }else{
        ylim1 <- ylim1 / prodcor
    }
    if(is.null(ylim2)){
        if(scaled){
            ylim2 <- c(0,1)
        }else{
            ylim2 <- range(unlist(tmp$y) / unlist(tmp$x))
        }
    }
    if(scaled){
        xplot <- lapply(tmp$x, function(x) x/max(x))
        yplot <- lapply(tmp$y, function(x) x/max(x))
        y22 <- lapply(as.list(1:length(tmp$y)),
                         function(i) tmp$y[[i]] / tmp$x[[i]])
        yplot2 <- lapply(y22, function(x) x/max(x))
        PBest <- (Pest/Best) / sapply(y22,max)
        Best <- Best / sapply(tmp$x,max)
        Pest <- Pest / sapply(tmp$y,max)
    }else{
        xplot <- tmp$x
        yplot <- tmp$y
        yplot2 <- lapply(as.list(1:length(tmp$y)),
                         function(i) tmp$y[[i]] / tmp$x[[i]])
        PBest <- Pest/Best
    }
    if(plot.obs){
        xlim <- range(c(xlim, Best))
        ylim1 <- range(c(ylim1, Pest))
        ylim2 <- range(c(ylim2, PBest))
    }
    if(1 %in% sel){
    plot(xplot[[1]], yplot[[1]], ty = 'n',
         xlab = "", ylab = "",
         xaxt = xaxt,
         yaxt = yaxt,
         xlim = xlim,
         ylim = c(0,1.2) * ylim1)
    if(best){
        rect(par("usr")[1], par("usr")[3],
             par("usr")[2], par("usr")[4],
             col = adjustcolor("goldenrod2",0.2))
    }
    if(black.line){
        lines(xplot[[1]], yplot[[1]],
              col = "black",
              lwd = 1.5)
    }else{
        for(i in 1:length(yplot)){
            lines(xplot[[i]], yplot[[i]],
                  col = adjustcolor(cols[indi.y[i]], 0.3),
                  lwd = 1.5)
        }
    }
    if(plot.obs){
        points(Best, Pest, col = cols[indi.y])
    }
    if(!is.null(lab2)){
        legend("topleft",
               legend = lab2,
               x.intersp = -0.5,
               cex = 0.8,
               bg = "white")
    }
    if(!is.null(title)){
        mtext(title, 3, 0.5, font = 2)
    }
    if(legend){

        embedPlot2(expr = expression({
            imageScale2(x, col = y, axis.pos = 2)
            ## mtext(text = "Year", side = 1, line = 2)
        }),
        ## at = c(0.5,0.98,0.15,0.2),
        at = c(0.9,0.95,0.05,0.95),
            envir = list(x = as.numeric(names(tmp$x)),
                         y = cols))

        ## indi <- seq(1,length(tmp$x), length.out = 5)
        ## legend("topright",
        ##        legend = c("const", names(tmp$x)[indi]),
        ##        col = c(1,cols[indi]),
        ##        ncol = 1,
        ##        pch = c(NA,rep(1,6)),
        ##        lwd = 1.5)
    }
    box(lwd = 1.5)
    ## if(labels) mtext(expression("dB" ~ dt^{-1} * " ['000 t" ~ yr^{-1}*"]"), 2, 2.5)
    ## if(!is.null(models)) mtext(models[1], 3, 0.5, font = font.models)
    }

}


embedPlot2 <- function (expr = expression({
    plot(1, t = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, font = 3, labels = "embedPlot error:\nexpression\nnot defined")
}), at = c(0.5, 0.95, 0.6, 0.95),
envir = ls())
{
    space_convert <- function(vec1, vec2) {
        vec1[1:2] <- vec1[1:2] * diff(vec2)[1] + vec2[1]
        vec1[3:4] <- vec1[3:4] * diff(vec2)[3] + vec2[3]
        vec1
    }
    plt <- par("plt")
    plt_space <- space_convert(at, plt)
    par(plt = plt_space, new = TRUE)
    eval(expr = expr, envir = envir)
    par(plt = plt)
}


imageScale2 <- function (z, zlim, col = hcl.colors(12, "YlOrRd", rev = TRUE),
                         breaks, axis.pos = 1, add.axis = TRUE, xlim = NULL,
                         ylim = NULL,
                         at = NULL,
                         ...){
    if (!missing(breaks)) {
        if (length(breaks) != (length(col) + 1)) {
            stop("must have one more break than colour")
        }
    }
    if (missing(breaks) & !missing(zlim)) {
        breaks <- seq(zlim[1], zlim[2], length.out = (length(col) +
            1))
    }
    if (missing(breaks) & missing(zlim)) {
        zlim <- range(z, na.rm = TRUE)
        breaks <- seq(zlim[1], zlim[2], length.out = (length(col) +
            1))
    }
    poly <- vector(mode = "list", length(col))
    for (i in seq(poly)) {
        poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1],
            breaks[i])
    }
    if (axis.pos %in% c(1, 3)) {
        YLIM <- c(0, 1)
        XLIM <- range(breaks)
    }
    if (axis.pos %in% c(2, 4)) {
        YLIM <- range(breaks)
        XLIM <- c(0, 1)
    }
    if (!missing(ylim)) {
        YLIM <- ylim
    }
    if (!missing(xlim)) {
        XLIM <- xlim
    }
    plot(1, 1, t = "n", ylim = YLIM, xlim = XLIM, axes = FALSE,
        xlab = "", ylab = "", xaxs = "i", yaxs = "i", ...)
    for (i in seq(poly)) {
        if (axis.pos %in% c(1, 3)) {
            polygon(poly[[i]], c(0, 0, 1, 1), col = col[i], border = col[i],
                lwd = 0.01)
        }
        if (axis.pos %in% c(2, 4)) {
            polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = col[i],
                lwd = 0.01)
        }
    }
    box()
    if (add.axis) {
        axis(axis.pos, at = at)
    }
}



plot.time.single3 <- function(resList, stockData, stockInfo, selMods = NULL,
                             maxYear = 2021, var = "logrvec", scale = 1,
                             rel = 1,
                             add.ices = TRUE, xlim = NULL,
                             ylim = NULL, xlab = "Time",
                             ylab = "Instrinsic growth rate",
                             xfixed = TRUE, yfixed = FALSE, plot.ci = TRUE,
                             plot.obs = FALSE, mfrow = c(2,5),
                             plot.ref = FALSE,
                             remove.cm.roundfish = FALSE,
                             byrow = FALSE,
                             legend.extra = FALSE,
                             plot.legend = TRUE,
                             adj = c(0,0),
                             title = NULL,
                             title.line = 0.5,
                             stock.in.plot = TRUE,
                             stock.in.plot.line = 3,
                             bestMods = NULL,
                             leg.text = NULL,
                             drop.mods = NULL,
                             yaxt = FALSE,
                             yaxt2 = FALSE,
                             ylab2 = "",
                             nonConv2 = NULL,
                             add.lab.right = FALSE
                             ){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
       names(resListSel) <- names(resList)
    }else{
        resListSel <- resList
        names(resListSel) <- names(resList)
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    for(i in 1:length(resListSel)){
        if(is.null(xlim0)) xlim <- c(1957,2021)
        if(is.null(ylim0) && plot.ci){
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish %in% c(1:3)){
                ylim <- range(unlist(lapply(resList[[i]][names(resList[[i]]) != "CM"], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
                if(remove.cm.roundfish == 2){
                    if(var[1] == "logMSYvec"){
                        varii <- "logmvec"
                    }else{
                        varii <- var[1]
                    }
                    ylim2 <- range(unlist(lapply(resList[[i]], function(x) x[[varii]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
                    ylim <- range(ylim,ylim2)
                }
            }else{
                if(var[1] == "logMSYvec"){
                    varii <- "logmvec"
                }else{
                    varii <- var[1]
                }
                ylim <- range(unlist(lapply(resList[[i]],
                                            function(x) x[[varii]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
            }
        }else if(is.null(ylim0)){
            ylim <- range(unlist(lapply(resList[[i]], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
        }
        if(add.ices || plot.obs){
            ices.dat <- get.ices.data(stockData, stockInfo, var, i, maxYear)
            ylim <- range(ylim, ices.dat$ydat / scale[1] / rel)
        }
        if(any(is.infinite(ylim))) ylim <- c(0,1)
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        ylim <- c(1,1.05) * ylim
        if(!yfixed && i == mfrow[2]) ylim <- c(1,1.05) * ylim
        if(length(resListSel[[i]]) == 1 &&
           all(is.na(resListSel[[i]][[1]][[vari]]))){
            plot.new()
        if(i == 1 && !is.null(title)) mtext(title[[1]], 3, title.line, font = 2)
        }else{
            if(paste0(stocks.true[i],"_",selMods[[1]]) %in% nonConv2){
                adji <- 0.5
                colii <- gray(0.6)
                flagi <- TRUE
            }else{
                adji <- 1
                colii <- "black"
                flagi <- FALSE
            }
        plot(1, 1,
             ty = "n",
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = "n",
             xlab = "",
             ylab = ""
             )
        ati <- pretty(c(0.9, 1.1) * ylim)
        if(plot.ref){
            abline(h = 1, col = adjustcolor("grey70",1), lwd = 1, lty = 2)
        }else{
            abline(h = ati, col = adjustcolor("grey70",0.3), lwd = 1, lty = 2)
        }
        ## if(!yfixed || i %in% sapply(seq(mfrow[1]),
        ##                             function(x) (x-1) * (mfrow[2]) + 1))
        ##     axis(2, at = ati, labels = ati)
        ## if(!yfixed && length(var) > 1)
        ##     axis(4, at = ati, labels = round(ati * rel,2))
        if(yaxt[i]){
            axis(2, at = ati, labels = ati)
        }
        if(yaxt2[i]){
            axis(4, at = ati, labels = round(ati * rel,2))
        }
        for(j in 1:length(resListSel[[i]])){
            years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            indi <- years < maxYear
            years <- years[indi]
            if(var[1] == "logMSYvec"){
                varii <- "logmvec"
            }else{
                varii <- var[1]
            }
            ydat <- resListSel[[i]][[j]][[varii]][which(indi),1:3] / scale[1]
            if(plot.ci){
                if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                   remove.cm.roundfish %in% c(1:3) &&
                   names(resListSel[[i]])[j] == "CM"){
                }else{
                    polygon(c(years, rev(years)),
                            c(ydat[,1],rev(ydat[,3])),
                            col = adjustcolor(cols[indiList[[i]][j]], 0.2 * adji),
                            border = NA)
                }
            }
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish == 1 && names(resListSel[[i]])[j] == "CM"){
            }else{
                lines(years, ydat[,2],
                      col = adjustcolor(cols[indiList[[i]][j]], adji),
                      lwd = 2)
            }
        }
        if(add.ices || plot.obs){
            if(!is.null(ices.dat)){
                if(plot.obs){
                    points(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                           col = adjustcolor("grey70",adji))
                }else if(add.ices){
                    lines(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                          col = adjustcolor(1,adji),
                          lwd = 1.5, lty = 2)
                }
            }
        }
        if(!is.null(bestMods)){
            j = which(names(resListSel[[i]]) == bestMods[[i]])
            if(length(j) > 0){
            ##     years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            ##     indi <- years < maxYear
            ##     years <- years[indi]
            ##     if(var[1] == "logMSYvec"){
            ##         varii <- "logmvec"
            ##     }else{
            ##         varii <- var[1]
            ##     }
            ##     ydat <- resListSel[[i]][[j]][[varii]][which(indi),1:3] / scale[1]
            ##     points(tail(years,1), tail(ydat[,2],1), pch = "*",
            ##            col = cols[indiList[[i]][j]],
                ##            cex = 3)
                rect(par("usr")[1], par("usr")[3],
                     par("usr")[2], par("usr")[4],
                     col = adjustcolor("goldenrod2",0.2))
            }
        }
        if(as.integer(stock.in.plot) == 1)
            legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)
        if(i == 1 && !legend.extra && plot.legend){
            if(is.null(leg.text)){
                texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
            }else{
                texti <- leg.text
            }
            if(!is.null(drop.mods)) texti <- texti[-drop.mods]
            ltyi <- unique(unlist(lapply(indiList, function(x) ltys[x])))
            if(!is.null(drop.mods)) ltyi <- ltyi[-drop.mods]
            coli <- unique(unlist(lapply(indiList, function(x) cols[x])))
            if(!is.null(drop.mods)) coli <- coli[-drop.mods]
            legend("topright", legend = texti,
                   lty = ltyi,
                   col = coli,
                   lwd = 2,
                   bg = "white",
                   ncol = 1,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
            if(flagi) points(2021, ylim[2]*0.98, pch = 8, xpd = TRUE)
        box(lwd = 1.5, col = colii)
        if(i == 1 && !is.null(title)) mtext(title[[1]], 3, title.line, font = 2)
        if(i == 1 && !is.null(title) && length(title) > 1)
            mtext(title[[2]], 3, 0.5, font = 1)
        if(as.integer(stock.in.plot) == 2)
            mtext(names(resList)[i], 4, stock.in.plot.line, font = 1)
        }
        if(add.lab.right){
            mtext(names(resList)[i],4,5.5, font = 2)
        }
    }
    mtext(xlab,1,2.6,outer = TRUE)
    mtext(ylab,2,2.5,outer = TRUE, padj = adj[1])
    mtext(ylab2,4,3,outer = TRUE, padj = adj[2])
    if(legend.extra){
        plot.new()
        texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
        if(!is.null(drop.mods)) texti <- texti[-drop.mods]
        ltyi <- unique(unlist(lapply(indiList, function(x) ltys[x])))
        if(!is.null(drop.mods)) ltyi <- ltyi[-drop.mods]
        coli <- unique(unlist(lapply(indiList, function(x) cols[x])))
        if(!is.null(drop.mods)) coli <- coli[-drop.mods]
        legend("center", legend = texti,
               title = "Models", title.cex = 1.2, title.font = 1,
               lty = ltyi,
               col = coli,
               lwd = 2,
               bg = "white",
               ncol = 2,
               x.intersp = 1.1,
               y.intersp = 1.2,
               cex = 1.2)
    }
}


plot.time.all3 <- function(resList, stockData, stockInfo, selMods = NULL,
                          maxYear = 2021, nPar = 2,
                          adj = list(c(0,0),c(0,0),c(0,0)),
                          bestMods = NULL, add.ices = TRUE,
                          leg.text = NULL, drop.mods = NULL,
                          remove.scenario = TRUE,
                          quant = 1,
                          nonConv2 = NULL){

    var <- c(1,2)

    nstocks <- length(resList)
    nmods <- length(resList[[1]])

    mfrow = c(nstocks, nmods)
    byrow = FALSE
    xfixed = TRUE
    yfixed = FALSE
    widths <- rep(1, nmods)
    if(as.integer(remove.scenario) %in% c(1:3)){
        widths[1] <- 1.3
        widths <- widths / mean(widths)
    }
    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow),
           widths = widths)
    par(oma = c(5,5,3,7.5))

    if(quant == 1){
        vars <- var2vars(c("r","Fmsy"))
    }else if(quant == 2){
        vars <- var2vars(c("K","Bmsy"))
    }else if(quant == 3){
        vars <- var2vars(c("m","MSY"))
    }

    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = "Time"
    ylab = vars$lab[[1]]
    ylab2 <- vars$lab[[2]]
    plot.ci = TRUE
    plot.obs = FALSE
    plot.ref = FALSE
    remove.cm.roundfish = remove.scenario
    legend.extra = FALSE
    rel = 1/nPar
    plot.legend = FALSE


    for(i in 1:nmods){
        if(as.integer(remove.scenario) %in% c(1:3)){
            if(i == 1){
                par(mar = c(0.5,0.5,0.5,5))
            }else{
                par(mar = c(0.5,0.5,0.5,0.5))
            }
        }else{
            par(mar = c(0.5,0.5,0.5,0.5))
        }

        title <- mods[i]
        selMods <- as.list(rep(title, nstocks))
        yaxt <- rep(ifelse(i == 1, TRUE, FALSE), nstocks)
        yaxt2 <- rep(ifelse(i == nmods, TRUE, FALSE), nstocks)
        if(i == nmods - 1){
            yaxt2 <- sapply(resList, function(x) all(is.na(x[[nmods]][["logB"]])))
        }
        if(as.integer(remove.scenario) %in% c(1:3) && i == 1){
            yaxt2 <- rep(TRUE, nstocks)
        }
        if(as.integer(remove.scenario) %in% c(1:3) && i == 2){
            yaxt <- rep(TRUE, nstocks)
        }
        if(as.integer(remove.scenario) == 3){
            remove.cm.roundfish <- ifelse(i == 1, 2, 3)
        }
        plot.time.single3(resList, stockData, stockInfo, selMods = selMods,
                          maxYear = maxYear, var = vars$var, scale = vars$scale,
                          rel = rel,
                          add.ices = add.ices, xlim = xlim,
                          ylim = ylim, xlab = ifelse(i == 1, xlab, ""),
                          ylab = ifelse(i == 1, ylab, ""),
                          xfixed = xfixed, yfixed = yfixed,
                          plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                          plot.ref = plot.ref,
                          remove.cm.roundfish = remove.cm.roundfish,
                          byrow = byrow,
                          legend.extra = legend.extra,
                          plot.legend = plot.legend,
                          adj = adj[[1]],
                          title = title,
                          stock.in.plot = FALSE,
                          bestMods = bestMods,
                          leg.text = leg.text,
                          drop.mods = drop.mods,
                          yaxt = yaxt,
                          yaxt2 = yaxt2,
                          ylab2 = ifelse(i == nmods, ylab2, ""),
                          nonConv2 = nonConv2,
                          add.lab.right = ifelse(i == nmods, TRUE, FALSE)
                          )
    }

}


plot.time.single4 <- function(resList, stockData, stockInfo, selMods = NULL,
                             maxYear = 2021, var = "logrvec", scale = 1,
                             rel = 1,
                             add.ices = TRUE, xlim = NULL,
                             ylim = NULL, xlab = "Time",
                             ylab = "Instrinsic growth rate",
                             xfixed = TRUE, yfixed = FALSE, plot.ci = TRUE,
                             plot.obs = FALSE, mfrow = c(2,5),
                             plot.ref = FALSE,
                             remove.cm.roundfish = FALSE,
                             byrow = FALSE,
                             legend.extra = FALSE,
                             plot.legend = TRUE,
                             adj = c(0,0),
                             title = NULL,
                             title.line = 0.5,
                             stock.in.plot = TRUE,
                             stock.in.plot.line = 3,
                             bestMods = NULL,
                             leg.text = NULL,
                             drop.mods = NULL,
                             yaxt = FALSE,
                             yaxt2 = FALSE,
                             ylab2 = "",
                             nonConv2 = NULL,
                             add.lab.right = FALSE
                             ){

    xlim0 <- xlim
    ylim0 <- ylim
    cols <- spmprod.cols()
    pchs <- spmprod.pchs()
    ltys <- spmprod.ltys()

    if(!is.null(selMods)){
        resListSel <- indiList <- vector("list",length(resList))
        for(i in 1:length(resList)){
            indi <- which(names(resList[[i]]) %in% selMods[[i]])
            resListSel[[i]] <- resList[[i]][indi]
            names(resListSel[[i]]) <- names(resList[[i]])[indi]
            indiList[[i]] <- indi
        }
       names(resListSel) <- names(resList)
    }else{
        resListSel <- resList
        names(resListSel) <- names(resList)
        indiList <- lapply(resList, function(x) seq(length(x)))
    }

    for(i in 1:length(resListSel)){
        if(is.null(xlim0)) xlim <- c(1957,2021)
        if(is.null(ylim0) && plot.ci){
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish %in% c(1:3)){
                ylim <- range(unlist(lapply(resList[[i]][names(resList[[i]]) != "CM"], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
                if(remove.cm.roundfish == 2){
                    if(var[1] == "logMSYvec"){
                        varii <- "logmvec"
                    }else{
                        varii <- var[1]
                    }
                    ylim2 <- range(unlist(lapply(resList[[i]], function(x) x[[varii]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
                    ylim <- range(ylim,ylim2)
                }
            }else{
                if(var[1] == "logMSYvec"){
                    varii <- "logmvec"
                }else{
                    varii <- var[1]
                }
                ylim <- range(unlist(lapply(resList[[i]],
                                            function(x) x[[varii]][which(as.numeric(rownames(x$logB)) < maxYear),1:3] / scale[1])))
            }
        }else if(is.null(ylim0)){
            ylim <- range(unlist(lapply(resList[[i]], function(x) x[[var[1]]][which(as.numeric(rownames(x$logB)) < maxYear),2] / scale[1])))
        }
        if(var[1] == "logC"){
            ylim <- range(unlist(lapply(resList[[i]], function(x){
                tmp <- x[["logB"]][which(as.numeric(rownames(x$logB)) < maxYear),2] * x[["logFnotS"]][which(as.numeric(rownames(x$logB)) < maxYear),2] * diff(x[["time"]])[1]
                tmp <- spict::annual(as.numeric(rownames(x$logB))[which(as.numeric(rownames(x$logB)) < maxYear)], tmp, sum)
                tmp$annvec / scale[1]
            })))
        }
        if(add.ices || plot.obs){
            ices.dat <- get.ices.data(stockData, stockInfo, var, i, maxYear)
            ylim <- range(ylim, ices.dat$ydat / scale[1] / rel)
        }
        if(any(is.infinite(ylim))) ylim <- c(0,1)
        if((!xfixed || i %in% (prod(mfrow)-mfrow[2]+1):prod(mfrow)) && byrow)
            xaxt = "s" else xaxt = "n"
        if((!xfixed || i %in% seq(mfrow[1],prod(mfrow),mfrow[1]) ||
            i == length(resListSel)) && !byrow)
            xaxt = "s" else xaxt = "n"
        ylim <- c(1,1.05) * ylim
        if(!yfixed && i == mfrow[2]) ylim <- c(1,1.05) * ylim
        if(length(resListSel[[i]]) == 1 &&
           all(is.na(resListSel[[i]][[1]][[vari]]))){
            plot.new()
            if(i == 1 && !is.null(title)) mtext(title[[1]], 3, title.line, font = 2)
        }else{
            if(paste0(stocks.true[i],"_",selMods[[1]]) %in% nonConv2){
                adji <- 0.5
                colii <- gray(0.6)
                flagi <- TRUE
            }else{
                adji <- 1
                colii <- "black"
                flagi <- FALSE
            }
        plot(1, 1,
             ty = "n",
             xlim = xlim,
             ylim = ylim,
             xaxt = xaxt,
             yaxt = "n",
             xlab = "",
             ylab = ""
             )
        ati <- pretty(c(0.9, 1.1) * ylim)
        if(plot.ref){
            abline(h = 1, col = adjustcolor("grey70",1), lwd = 1, lty = 2)
        }else{
            abline(h = ati, col = adjustcolor("grey70",0.3), lwd = 1, lty = 2)
        }
        ## if(!yfixed || i %in% sapply(seq(mfrow[1]),
        ##                             function(x) (x-1) * (mfrow[2]) + 1))
        ##     axis(2, at = ati, labels = ati)
        ## if(!yfixed && length(var) > 1)
        ##     axis(4, at = ati, labels = round(ati * rel,2))
        if(yaxt[i]){
            axis(2, at = ati, labels = ati)
        }
        if(yaxt2[i]){
            axis(4, at = ati, labels = round(ati * rel,2))
        }
            for(j in 1:length(resListSel[[i]])){
            years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            indi <- years < maxYear
            years <- years[indi]
            if(var[1] == "logMSYvec"){
                varii <- "logmvec"
            }else{
                varii <- var[1]
            }
                if(varii == "logC"){
                    ydat <- resListSel[[i]][[j]][["logB"]][which(indi),1:3] * resListSel[[i]][[j]][["logFnotS"]][which(indi),1:3] * diff(resListSel[[i]][[j]]$time)[1]
                    tmp1 <- annual(as.numeric(rownames(resListSel[[i]][[j]][["logB"]][which(indi),])), ydat[,1], sum)
                    tmp2 <- annual(as.numeric(rownames(resListSel[[i]][[j]][["logB"]][which(indi),])), ydat[,2], sum)
                    tmp3 <- annual(as.numeric(rownames(resListSel[[i]][[j]][["logB"]][which(indi),])), ydat[,3], sum)
                    years <- tmp1$anntime
                    ydat <- data.frame(ll = tmp1$annvec,
                                       est = tmp2$annvec,
                                       ul = tmp3$annvec)
                    ydat <- ydat / scale[1]
                }else{
                    ydat <- resListSel[[i]][[j]][[varii]][which(indi),1:3] / scale[1]
                }

            if(plot.ci){
                if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
                   remove.cm.roundfish %in% c(1:3) &&
                   names(resListSel[[i]])[j] == "CM"){
                }else{
                    polygon(c(years, rev(years)),
                            c(ydat[,1],rev(ydat[,3])),
                            col = adjustcolor(cols[indiList[[i]][j]], 0.2 * adji),
                            border = NA)
                }
            }
            if((names(resListSel)[i] %in% c("Cod","Haddock","Saithe","Whiting")) &&
               remove.cm.roundfish == 1 && names(resListSel[[i]])[j] == "CM"){
            }else{
                lines(years, ydat[,2],
                      col = adjustcolor(cols[indiList[[i]][j]], adji),
                      lwd = 2)
            }
        }
        if(add.ices || plot.obs){
            if(!is.null(ices.dat)){
                names(ices.dat)
                ices.dat$ydat
                if(plot.obs){
                    points(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                           col = "black")
                           ## col = adjustcolor("grey70",adji))
                }else if(add.ices){
                    lines(ices.dat$years, ices.dat$ydat / scale[1] / rel,
                          col = adjustcolor(1,adji),
                          lwd = 1.5, lty = 2)
                }
            }
        }
        if(!is.null(bestMods)){
            j = which(names(resListSel[[i]]) == bestMods[[i]])
            if(length(j) > 0){
            ##     years <- as.numeric(rownames(resListSel[[i]][[j]]$logB))
            ##     indi <- years < maxYear
            ##     years <- years[indi]
            ##     if(var[1] == "logMSYvec"){
            ##         varii <- "logmvec"
            ##     }else{
            ##         varii <- var[1]
            ##     }
            ##     ydat <- resListSel[[i]][[j]][[varii]][which(indi),1:3] / scale[1]
            ##     points(tail(years,1), tail(ydat[,2],1), pch = "*",
            ##            col = cols[indiList[[i]][j]],
                ##            cex = 3)
                rect(par("usr")[1], par("usr")[3],
                     par("usr")[2], par("usr")[4],
                     col = adjustcolor("goldenrod2",0.2))
            }
        }
        if(as.integer(stock.in.plot) == 1)
            legend("topleft", legend = names(resList)[i],
               pch = NA, bg = "white", text.font = 2,
               x.intersp = 0, cex = 1.2)
        if(i == 1 && !legend.extra && plot.legend){
            if(is.null(leg.text)){
                texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
            }else{
                texti <- leg.text
            }
            if(!is.null(drop.mods)) texti <- texti[-drop.mods]
            ltyi <- unique(unlist(lapply(indiList, function(x) ltys[x])))
            if(!is.null(drop.mods)) ltyi <- ltyi[-drop.mods]
            coli <- unique(unlist(lapply(indiList, function(x) cols[x])))
            if(!is.null(drop.mods)) coli <- coli[-drop.mods]
            legend("topright", legend = texti,
                   lty = ltyi,
                   col = coli,
                   lwd = 2,
                   bg = "white",
                   ncol = 1,
                   x.intersp = 1.1,
                   y.intersp = 1.0,
                   cex = 1.1)
        }
            if(flagi) points(2021, ylim[2]*0.98, pch = 8, xpd = TRUE)
        box(lwd = 1.5, col = colii)
        if(i == 1 && !is.null(title)) mtext(title[[1]], 3, title.line, font = 2)
        if(i == 1 && !is.null(title) && length(title) > 1)
            mtext(title[[2]], 3, 0.5, font = 1)
        if(as.integer(stock.in.plot) == 2)
            mtext(names(resList)[i], 4, stock.in.plot.line, font = 1)
        }
        if(add.lab.right){
            mtext(names(resList)[i],4,2, font = 2)
        }
    }
    mtext(xlab,1,2.6,outer = TRUE)
    mtext(ylab,2,2.5,outer = TRUE, padj = adj[1])
    ## mtext(ylab2,4,3,outer = TRUE, padj = adj[2])
    if(legend.extra){
        plot.new()
        texti <- unique(unlist(lapply(resListSel, function(x) names(x))))
        if(!is.null(drop.mods)) texti <- texti[-drop.mods]
        ltyi <- unique(unlist(lapply(indiList, function(x) ltys[x])))
        if(!is.null(drop.mods)) ltyi <- ltyi[-drop.mods]
        coli <- unique(unlist(lapply(indiList, function(x) cols[x])))
        if(!is.null(drop.mods)) coli <- coli[-drop.mods]
        legend("center", legend = texti,
               title = "Models", title.cex = 1.2, title.font = 1,
               lty = ltyi,
               col = coli,
               lwd = 2,
               bg = "white",
               ncol = 2,
               x.intersp = 1.1,
               y.intersp = 1.2,
               cex = 1.2)
    }
}

plot.time.all4 <- function(resList, stockData, stockInfo, selMods = NULL,
                          maxYear = 2021, nPar = 2,
                          adj = list(c(0,0),c(0,0),c(0,0)),
                          bestMods = NULL, add.ices = TRUE,
                          leg.text = NULL, drop.mods = NULL,
                          remove.scenario = TRUE,
                          quant = 1,
                          nonConv2 = NULL){

    var <- c(1,2)

    nstocks <- length(resList)
    nmods <- length(resList[[1]])

    mfrow = c(nstocks, nmods)
    byrow = FALSE
    xfixed = TRUE
    yfixed = FALSE
    widths <- rep(1, nmods)
    layout(matrix(1:prod(mfrow),mfrow[1],mfrow[2],byrow=byrow))
    par(oma = c(5,5,3,4))

    if(quant == 1){
        vars <- var2vars(c("r","Fmsy"))
    }else if(quant == 2){
        vars <- var2vars(c("K","Bmsy"))
    }else if(quant == 3){
        vars <- var2vars(c("m","MSY"))
    }else if(quant == 4){
        vars <- var2vars("B")
    }else if(quant == 5){
        vars <- var2vars("F")
    }else if(quant == 6){
        vars <- var2vars("C")
    }

    add.ices = add.ices
    xlim = NULL
    ylim = NULL
    xlab = "Time"
    ylab = vars$lab[[1]]
    ylab2 <- NULL ## vars$lab[[2]]
    plot.ci = TRUE
    plot.obs = TRUE
    plot.ref = FALSE
    remove.cm.roundfish = remove.scenario
    legend.extra = FALSE
    rel = 1 ## /nPar
    plot.legend = FALSE


    for(i in 1:nmods){
        if(as.integer(remove.scenario) %in% c(1:3)){
            ## if(i == 1){
            ##     par(mar = c(0.5,0.5,0.5,5))
            ## }else{
            ## }
                par(mar = c(0.5,0.5,0.5,0.5))
        }else{
            par(mar = c(0.5,0.5,0.5,0.5))
        }

        title <- mods[i]
        selMods <- as.list(rep(title, nstocks))
        yaxt <- rep(ifelse(i == 1, TRUE, FALSE), nstocks)
        yaxt2 <- rep(FALSE, nstocks)
        ## if(i == nmods - 1){
        ##     yaxt2 <- sapply(resList, function(x) all(is.na(x[[nmods]][["logB"]])))
        ## }
        ## if(as.integer(remove.scenario) %in% c(1:3) && i == 1){
        ##     yaxt2 <- rep(TRUE, nstocks)
        ## }
        ## if(as.integer(remove.scenario) %in% c(1:3) && i == 2){
        ##     yaxt <- rep(TRUE, nstocks)
        ## }
        if(as.integer(remove.scenario) == 3){
            remove.cm.roundfish <- ifelse(i == 1, 2, 3)
        }
        plot.time.single4(resList, stockData, stockInfo, selMods = selMods,
                          maxYear = maxYear, var = vars$var, scale = vars$scale,
                          rel = rel,
                          add.ices = add.ices, xlim = xlim,
                          ylim = ylim, xlab = ifelse(i == 1, xlab, ""),
                          ylab = ifelse(i == 1, ylab, ""),
                          xfixed = xfixed, yfixed = yfixed,
                          plot.ci = plot.ci, plot.obs = plot.obs, mfrow = mfrow,
                          plot.ref = plot.ref,
                          remove.cm.roundfish = remove.cm.roundfish,
                          byrow = byrow,
                          legend.extra = legend.extra,
                          plot.legend = plot.legend,
                          adj = adj[[1]],
                          title = title,
                          stock.in.plot = FALSE,
                          bestMods = bestMods,
                          leg.text = leg.text,
                          drop.mods = drop.mods,
                          yaxt = yaxt,
                          yaxt2 = yaxt2,
                          ylab2 = ifelse(i == nmods, ylab2, ""),
                          nonConv2 = nonConv2,
                          add.lab.right = ifelse(i == nmods, TRUE, FALSE)
                          )
    }

}
