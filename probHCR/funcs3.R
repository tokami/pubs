


## get indices of converged replicates
get.convInd <- function(conv, index = 1:dim(conv)[2], evalyears = 1:dim(conv)[3]){
    if(max(evalyears) > dim(conv)[3]) evalyears <- min(evalyears):dim(conv)[3] ## conv only length 19
    sel <- conv[,index,evalyears]
    conv <- apply(sel, c(1,2), function(x) all(x == 2))
    convInd <- which(apply(conv, 1, function(x) all(x == TRUE)) == TRUE)
}

## get quantities of converged replicates only
get.conv <- function(x, indexConv = 1:dim(x$bbmsy)[2],
                     keepExtra = NULL, all = FALSE,
                     evalyears = NULL,
                     convInfo = NULL){

    ## evalyears
    if(is.null(evalyears)) evalyears <- 1:dim(x$bbmsy)[3]

    ## converged replicates
    if(is.null(convInfo)){
        if(!all){
            convInd <- get.convInd(x$conv, indexConv, evalyears)
        }else convInd <- 1:dim(x$bbmsy)[1]
    }else{
        convInd <- convInfo
    }

    ## hcrs
    keepInd <- sort(unique(c(keepExtra, indexConv)))

    reslist <- list()

    ## bbmsy
    reslist$bbmsy <- x$bbmsy[convInd,keepInd,]

    ## ffmsy
    reslist$ffmsy <- x$ffmsy[convInd,keepInd,]

    ## catch
    reslist$catch <- x$catch[convInd,keepInd,]

    ## refcatch
    reslist$refcatch <- x$refcatch[convInd]

    ## tac
    reslist$tac <- x$tac[convInd,keepInd,]

    ## spict bbmsy
    reslist$spictBBMSY <- x$spictBBMSY[convInd,keepInd,]

    ## spict ffmsy
    reslist$spictFFMSY <- x$spictFFMSY[convInd,keepInd,]

    ## spict SD of bbmsy
    reslist$spictBsd <- x$spictBsd[convInd,keepInd,]

    ## spict SD of ffmsy
    reslist$spictFsd <- x$spictFsd[convInd,keepInd,]

    ## spict SD of Cp
    reslist$spictCsd <- x$spictCsd[convInd,keepInd,]

    ## spict MSY
    reslist$spictMSY <- x$spictMSY[convInd,keepInd,]

    ## spict Cp est
    reslist$spictCest <- x$spictCest[convInd,keepInd,]

    ## bmsy
    reslist$bmsy <- x$BMSY[convInd]

    ## fmsy
    reslist$fmsy <- x$FMSY[convInd]

    ## return
    return(reslist)
}

## check sufficiency of sample size!
check.sample.size <- function(x, ref = 0.3, evalyears = NULL){

    ## pblim Y3_5
    tmp <- apply(x$bbmsy, c(1,2,3), function(x) x < ref)
    dims <- dim(tmp)
    risks <- vector("list", dims[1])
    for(i in 2:dims[1]){
        risks[[i]] <- apply(tmp[1:i,,], c(2), function(x){
            x <- as.numeric(x[,3:5])
            prop.test(sum(x), n = length(x),
                      conf.level = 0.95, correct = FALSE)$estimate
        }
        )
    }
    pblimY3_5 <- do.call(rbind,risks)

    ## pblim Y6_20
    tmp <- apply(x$bbmsy, c(1,2,3), function(x) x < ref)
    dims <- dim(tmp)
    risks <- vector("list", dims[1])
    for(i in 2:dims[1]){
        risks[[i]] <- apply(tmp[1:i,,], c(2), function(x){
            x <- as.numeric(x[,6:20])
            prop.test(sum(x), n = length(x),
                      conf.level = 0.95, correct = FALSE)$estimate
        }
        )
    }
    pblimY6_20 <- do.call(rbind,risks)

    ## yield 3-5
    tmp <- x$catch
    for(i in 1:length(x$refcatch)){
        tmp[i,,] <- tmp[i,,] / x$refcatch[i]
    }
    dims <- dim(tmp)
    yields <- vector("list", dims[1])
    for(i in 2:dims[1]){
        yields[[i]] <- apply(tmp[1:i,,], c(2), function(x){
            median(as.numeric(x[,3:5]))
        }
        )
    }
    yieldY3_5 <- do.call(rbind,yields)

    ## yield 6-20
    tmp <- x$catch
    for(i in 1:length(x$refcatch)){
        tmp[i,,] <- tmp[i,,] / x$refcatch[i]
    }
    dims <- dim(tmp)
    yields <- vector("list", dims[1])
    for(i in 2:dims[1]){
        yields[[i]] <- apply(tmp[1:i,,], c(2), function(x){
            median(as.numeric(x[,6:20]))
        }
        )
    }
    yieldY6_20 <- do.call(rbind,yields)

    ## AAV
    tmp <- x$catch
    for(i in 1:length(x$refcatch)){
        tmp[i,,] <- tmp[i,,] / x$refcatch[i]
    }
    if(is.null(evalyears)) evalyears <- 1:dim(tmp)[3]
    if(min(evalyears) == 1){
        evalyears1 <- evalyears
        evalyears2 <- evalyears[-1]
    }else{
        evalyears1 <- c(min(evalyears)-1,evalyears)
        evalyears2 <- evalyears
    }
    tmp2 <- apply(tmp[,,evalyears1], c(1,2), function(x) sum(abs(diff(x))))
    tmp3 <- apply(tmp[,,evalyears2], c(1,2), sum)
    tmp4 <- tmp2/tmp3
    dims <- dim(tmp4)
    aavs <- vector("list", dims[1])
    for(i in 2:dims[1]){
        aavs[[i]] <- apply(tmp4[1:i,], c(2), function(x){
            median(as.numeric(x))
        }
        )
    }
    aav <- do.call(rbind,aavs)

    return(list(pblimY3_5 = pblimY3_5, pblimY6_20 = pblimY6_20,
                yieldY3_5 = yieldY3_5, yieldY6_20 = yieldY6_20,
                aav = aav))

}

## get RE between true and spict est. rel. states
get.re.bbmsy <- function(x, evalyears = NULL, years = TRUE){
    spictbbmsy <- x$spictBBMSY
    dims <- dim(x$bbmsy)
    ny <- dims[3]
    nhcr <- dims[2]
    nrep <- dims[1]
    bbmsy <- x$bbmsy[,,-ny]

    if(is.null(evalyears)) evalyears <- 1:dim(bbmsy)[3]
    if(max(evalyears) > dim(bbmsy)[3]) evalyears <- min(evalyears):dim(bbmsy)[3]

    RE <- array(NA, dim = c(nrep, nhcr, length(evalyears)))
    for(hcr in 1:nhcr){
        for(rep in 1:nrep){
            RE[rep,hcr,] <- abs((bbmsy[rep,hcr,evalyears] - spictbbmsy[rep,hcr,evalyears]) / bbmsy[rep,hcr,evalyears])
        }
    }

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(RE, dimKeep, function(x){
        x <- as.numeric(x)
        if(any(!is.na(x))){
            res <- wilcox.test(x,
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }else{
            c(med = NA, lo = NA, up = NA, qcv = NA)
        }
    }
    )

    return(test)
}

get.re.ffmsy <- function(x, evalyears = NULL, years = TRUE){
    spictffmsy <- x$spictFFMSY
    dims <- dim(x$ffmsy)
    ny <- dims[3]
    nhcr <- dims[2]
    nrep <- dims[1]
    ffmsy <- x$ffmsy[,,-ny]

    if(is.null(evalyears)) evalyears <- 1:dim(ffmsy)[3]
    if(max(evalyears) > dim(ffmsy)[3]) evalyears <- min(evalyears):dim(ffmsy)[3]

    RE <- array(NA, dim = c(nrep, nhcr, length(evalyears)))
    for(hcr in 1:nhcr){
        for(rep in 1:nrep){
            RE[rep,hcr,] <- abs((ffmsy[rep,hcr,evalyears] - spictffmsy[rep,hcr,evalyears]) / ffmsy[rep,hcr,evalyears])
        }
    }

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(RE, dimKeep, function(x){
        x <- as.numeric(x)
        if(any(!is.na(x))){
            res <- wilcox.test(x,
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }else{
            c(med = NA, lo = NA, up = NA, qcv = NA)
        }
    }
    )

    return(test)
}


get.re.tac <- function(x, evalyears = NULL, years = TRUE){

    tac <- x$tac
    dims <- dim(tac)
    nrep <- dims[1]
    nhcr <- dims[2]

    if(is.null(evalyears)) evalyears <- 1:dims[3]
    if(max(evalyears) > dims[3]) evalyears <- min(evalyears):dims[3]

    ind <- c(3:19)
    ind2 <- c(20:36)
    tacRE <- array(NA, dim = c(dims[1],nhcr,length(evalyears)))
    for(i in 1:17){
        tacRE[,2+i,] <- abs((tac[,1,evalyears] - tac[,ind[i],evalyears]) / tac[,1,evalyears])
        tacRE[,2+i+17,] <- abs((tac[,2,evalyears] - tac[,ind2[i],evalyears]) / tac[,2,evalyears])
    }

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(tacRE, dimKeep, function(x){
        x <- as.numeric(x)
        if(any(!is.na(x))){
            res <- wilcox.test(x,
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }else{
            c(med = NA, lo = NA, up = NA, qcv = NA)
        }
    }
    )

    return(test)
}


## get SD of BBmsy dist in spict
get.bbmsySD <- function(x, evalyears = NULL, years = TRUE){
    spictbsd <- x$spictBsd
    dims <- dim(x$bbmsy)
    ny <- dims[3]
    nhcr <- dims[2]
    nrep <- dims[1]
    bbmsy <- x$bbmsy[,,-ny]

    if(is.null(evalyears)) evalyears <- 1:dim(bbmsy)[3]
    if(max(evalyears) > dim(bbmsy)[3]) evalyears <- min(evalyears):dim(bbmsy)[3]
    spictbsd <- spictbsd[,,evalyears]

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(spictbsd, dimKeep, function(x){
        x <- as.numeric(x)
        if(any(!is.na(x))){
            res <- wilcox.test(x,
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }else{
            c(med = NA, lo = NA, up = NA, qcv = NA)
        }
    }
    )

    return(test)
}

## get SD of FFmsy dist in spict
get.ffmsySD <- function(x, evalyears = NULL, years = TRUE){
    spictfsd <- x$spictFsd
    dims <- dim(x$ffmsy)
    ny <- dims[3]
    nhcr <- dims[2]
    nrep <- dims[1]
    ffmsy <- x$ffmsy[,,-ny]

    if(is.null(evalyears)) evalyears <- 1:dim(ffmsy)[3]
    if(max(evalyears) > dim(ffmsy)[3]) evalyears <- min(evalyears):dim(ffmsy)[3]
    spictfsd <- spictfsd[,,evalyears]

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(spictfsd, dimKeep, function(x){
        x <- as.numeric(x)
        if(any(!is.na(x))){
            res <- wilcox.test(x,
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }else{
            c(med = NA, lo = NA, up = NA, qcv = NA)
        }
    }
    )

    return(test)
}

## get SD of Cp dist in spict
get.cpSD <- function(x, evalyears = NULL, years = TRUE){
    spictcsd <- x$spictCsd
    dims <- dim(x$ffmsy)
    ny <- dims[3]
    nhcr <- dims[2]
    nrep <- dims[1]
    ffmsy <- x$ffmsy[,,-ny]

    if(is.null(evalyears)) evalyears <- 1:dim(ffmsy)[3]
    if(max(evalyears) > dim(ffmsy)[3]) evalyears <- min(evalyears):dim(ffmsy)[3]
    spictcsd <- spictcsd[,,evalyears]

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(spictcsd, dimKeep, function(x){
        x <- as.numeric(x)
        if(any(!is.na(x))){
            res <- wilcox.test(x,
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }else{
            c(med = NA, lo = NA, up = NA, qcv = NA)
        }
    }
    )

    return(test)
}



get.bbmsy <- function(x, years = TRUE, interval = FALSE, cv = 0.1,evalyears = NULL){
    quant <- x$bbmsy

    if(years) dimKeep <- c(2,3) else dimKeep <- 2
    if(is.null(evalyears)) evalyears <- 1:dim(quant)[3]

    quant <- quant[,,evalyears]

    if(interval){
        low <- 1-cv
        upp <- 1+cv
        tmp <- apply(quant, c(1,2,3), function(x) low <= x & x <= upp)
        test <- apply(tmp, dimKeep, function(x){
            x <- as.numeric(x)
            res <- prop.test(sum(x), n = length(x), conf.level = 0.95, correct = FALSE)
            med <- res$estimate
            lo <- as.numeric(res$conf.int)[1]
            up <- as.numeric(res$conf.int)[2]
            round(c(med,lo,up),3)
        }
        )
    }else{
        test <- apply(quant, dimKeep, function(x){
            res <- wilcox.test(as.numeric(x),
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }
        )
    }
    return(test)
}



get.ffmsy <- function(x, years = TRUE, interval = FALSE, cv = 0.1,evalyears = NULL){
    quant <- x$ffmsy

    if(years) dimKeep <- c(2,3) else dimKeep <- 2
    if(is.null(evalyears)) evalyears <- 1:dim(quant)[3]

    quant <- quant[,,evalyears]

    if(interval){
        low <- 1-cv
        upp <- 1+cv
        tmp <- apply(quant, c(1,2,3), function(x) low <= x & x <= upp)
        test <- apply(tmp, dimKeep, function(x){
            x <- as.numeric(x)
            res <- prop.test(sum(x), n = length(x), conf.level = 0.95, correct = FALSE)
            med <- res$estimate
            lo <- as.numeric(res$conf.int)[1]
            up <- as.numeric(res$conf.int)[2]
            round(c(med,lo,up),3)
        }
        )
    }else{
        test <- apply(quant, dimKeep, function(x){
            res <- wilcox.test(as.numeric(x),
                               alternative="two.sided",
                               correct=TRUE,
                               conf.int=TRUE,
                               conf.level=0.95)
            med=res$estimate
            lo=res$conf.int[1]
            up=res$conf.int[2]
            ## "Quartile based Coefficient of Variation" [or QCV]
            qcv = ((up - lo) / med) * 100
            round(c(med=med,lo=lo,up=up, qcv = qcv),3)
        }
        )
    }
    return(test)
}



## get probability of overfishing
## for now defined as 30% of BMSY (citation: WKBUT 2013)
## == 50% of productivity (for Schaefer) ideally adjusted
## prop.test with correct = FALSE uses the Wilson score interval method
## for Bionomial porportion confidence intervals (alternative: binom.test for
## "exact" method = Clopper-Pearson interval
## check if intervals are correct if p == 0 or p == 1,
## if not there is "the rule of 3" for p == 0: (0, 3/n), for p == 1: (1-3/n,1)
get.pblim <- function(x, ref = 0.3, evalyears = NULL, years = TRUE){
    quant <- x$bbmsy
    tmp <- apply(quant, c(1,2,3), function(x) x < ref)
    if(is.null(evalyears)) evalyears <- 1:dim(tmp)[3]

    tmp <- tmp[,,evalyears]

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    test <- apply(tmp, dimKeep, function(x){
        x <- as.numeric(x)
        res <- prop.test(sum(x), n = length(x), conf.level = 0.95, correct = FALSE)
        round(c(res$estimate,lo=as.numeric(res$conf.int)[1],up=as.numeric(res$conf.int)[2]),3)
    }
    )
    return(test)
}


modifiedCox <- function(x){
    n <- length(x)
    y <- log(x)
    y.m <- mean(y)
    y.var <- var(y)

    my.t <- qt(0.975, df = n-1) # 95% Confidence interval

    my.mean <- mean(x)
    up <- y.m + y.var/2 + my.t*sqrt(y.var/n + y.var^2/(2*(n - 1)))
    lo <- y.m + y.var/2 - my.t*sqrt(y.var/n + y.var^2/(2*(n - 1)))

    return(list(up = exp(up), mean = my.mean, lo = exp(lo)))

}

get.timeRec <- function(x, evalyears = NULL, ref = 0.3, maxyear=20){
    tmp <- x$bbmsy
    if(is.null(evalyears)) evalyears <- 1:dim(tmp)[3]
    tmp <- tmp[,,evalyears]

    ## median
    test1 <- apply(tmp, c(1,2), function(x){
        tmp <- which(x > ref)
        if(length(tmp) == 0) tmp <- maxyear
        min(tmp)
    })

    ## median
    test <- apply(test1, 2, function(x){
        res <- wilcox.test(as.numeric(x),
                           alternative="two.sided",
                           correct=TRUE,
                           conf.int=TRUE,
                           conf.level=0.95)
        med=res$estimate
        lo=res$conf.int[1]
        up=res$conf.int[2]
        ## "Quartile based Coefficient of Variation" [or QCV]
        qcv = ((up - lo) / med) * 100
        round(c(med=med,lo=lo,up=up, qcv = qcv),3)
    }
    )
    return(test)
}


## get yield
## Yield not normally distributed => modified cox method for CIs of mean
## mean and 95% CIs based on the modified cox method (Olsson 2005: )
## get stats over reps for each year
get.yield <- function(x, evalyears = NULL, years = TRUE){
    tmp <- x$catch
    for(i in 1:length(x$refcatch)){
        tmp[i,,] <- tmp[i,,] / x$refcatch[i]
    }
    if(is.null(evalyears)) evalyears <- 1:dim(tmp)[3]

    if(years) dimKeep <- c(2,3) else dimKeep <- 2

    tmp <- tmp[,,evalyears]

    ## median
    test <- apply(tmp, dimKeep, function(x){
        res <- wilcox.test(as.numeric(x),
                           alternative="two.sided",
                           correct=TRUE,
                           conf.int=TRUE,
                           conf.level=0.95)
        med=res$estimate
        lo=res$conf.int[1]
        up=res$conf.int[2]
        ## "Quartile based Coefficient of Variation" [or QCV]
        qcv = ((up - lo) / med) * 100
        round(c(med=med,lo=lo,up=up, qcv = qcv),3)
    }
    )
    return(test)
}


## get AAV
get.aav <- function(x, evalyears = NULL){
    tmp <- x$catch
    for(i in 1:length(x$refcatch)){
        tmp[i,,] <- tmp[i,,] / x$refcatch[i]
    }

    if(is.null(evalyears)) evalyears <- 1:dim(tmp)[3]

    if(min(evalyears) == 1){
        evalyears1 <- evalyears
        evalyears2 <- evalyears[-1]
    }else{
        evalyears1 <- c(min(evalyears)-1,evalyears)
        evalyears2 <- evalyears
    }


    tmp2 <- apply(tmp[,,evalyears1], c(1,2), function(x) sum(abs(diff(x))))
    tmp3 <- apply(tmp[,,evalyears2], c(1,2), sum)

    tmp4 <- tmp2/tmp3

    test <- apply(tmp4, c(2), function(x){
        res <- wilcox.test(as.numeric(x),
                           alternative="two.sided",
                           correct=TRUE,
                           conf.int=TRUE,
                           conf.level=0.95)
        med=res$estimate
        lo=res$conf.int[1]
        up=res$conf.int[2]
        ## "Quartile based Coefficient of Variation" [or QCV]
        qcv = ((up - lo) / med) * 100
        round(c(med=med,lo=lo,up=up, qcv = qcv),3)
    }
    )
    return(test)
}


## Plots
##############################################################################


plot.fig1 <- function(x, scenario = 1,
                      hcrBoth = NULL,
                      hcrref = 1,
                      hcrrefKobe = c(1,2), yearKobe = 20,
                      hcrs1 = NULL, hcrs2 = NULL,
                      cols1 = NULL, cols2 = NULL){

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)
    cex = 1.9
    cols1P <- sapply(cols1, function(x) rgb(t(col2rgb(x))/255,alpha=0.2))
    cols2P <- sapply(cols2, function(x) rgb(t(col2rgb(x))/255,alpha=0.2))
    legendwidth <- 0.4

    ## Kobe lims
    hcrsUse <- c(hcrBoth, hcrs1, hcrs2)
    limBKobe = c(0,1) * range(unlist(lapply(x, function(y)
        c(y[[scenario]]$bbmsyYear[1:3,hcrsUse,yearKobe],
          y[[scenario]]$bbmsyYear[1,hcrrefKobe,yearKobe]))))
    limBKobe[2] <- max(limBKobe[2], 1)
    limFKobe = c(0,1) * range(unlist(lapply(x, function(y)
        c(y[[scenario]]$ffmsyYear[1:3,hcrsUse,yearKobe],
          y[[scenario]]$ffmsyYear[1,hcrrefKobe,yearKobe]))))
    limFKobe[2] <- max(limFKobe[2], 1)

    ## Y-R long-term lims
    hcrsUse <- c(hcrref,hcrBoth,hcrs1,hcrs2)
    limRlong <- c(0,1.05) * range(unlist(lapply(x, function(y) y[[scenario]]$pblimY6_20[1:3,hcrsUse])))
    limYlong <- c(0,1.05) * range(unlist(lapply(x, function(y) y[[scenario]]$yieldY6_20[1:3,hcrsUse])))

    ## Y-R short-term lims
    hcrsUse <- c(hcrref,hcrBoth,hcrs1,hcrs2)
    limRshort <- c(0,1.05) * range(unlist(lapply(x, function(y) y[[scenario]]$pblimY3_5[1:3,hcrsUse])))
    limYshort <- c(0,1.05) * range(unlist(lapply(x, function(y) y[[scenario]]$yieldY3_5[1:3,hcrsUse])))

    ## AAV lims
    hcrsUse <- c(hcrBoth,hcrs1,hcrs2)
    limA <- c(0,1.05) * range(unlist(lapply(x, function(y)
        c(y[[scenario]]$yieldDiff[1:3,hcrsUse],
          y[[scenario]]$yieldDiff[1,hcrref]))))

    ly <- layout(matrix(c(1:12,rep(13,4)),
                        byrow = FALSE,ncol=4,nrow=4),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,4))
    opar <- par(oma=c(0,6,2,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        tmp <- x[[spec]][[scenario]]

        ## Kobe plot
        bbmsy1 <- tmp$bbmsyYear[1:3,hcrs1,yearKobe]
        ffmsy1 <- tmp$ffmsyYear[1:3,hcrs1,yearKobe]
        bbmsy2 <- tmp$bbmsyYear[1:3,hcrs2,yearKobe]
        ffmsy2 <- tmp$ffmsyYear[1:3,hcrs2,yearKobe]
        bbmsyref <- tmp$bbmsyYear[1,hcrrefKobe,yearKobe]
        ffmsyref <- tmp$ffmsyYear[1,hcrrefKobe,yearKobe]
        bbmsyBo <- tmp$bbmsyYear[1:3,hcrBoth,yearKobe]
        ffmsyBo <- tmp$ffmsyYear[1:3,hcrBoth,yearKobe]
        if(spec == 3) par(mar=c(4,0,2,1)) else par(mar=c(4,0,2,0))
        plot(bbmsy1[1,],ffmsy1[1,],
             ylim = limFKobe, xlim = limBKobe, ty= 'n',
             ylab = "", xlab = "",
             xaxt = "n", yaxt = "n")
        ## Area between two refs
        polygon(c(range(bbmsyref),rev(range(bbmsyref))),
                c(rep(min(ffmsyref),2),rep(max(ffmsyref),2)),
                border = NA, col = rgb(t(col2rgb("grey40"))/255,alpha=0.15))
        abline(h=1, col="grey60",lwd=2)
        abline(v=1, col="grey60",lwd=2)
        ## points
        for(i in 1:nhcrs1){
            polygon(bbmsy1[c(2,3,3,2),i],ffmsy1[c(2,2,3,3),i],
                    col=cols1P[i], border = NA)
        }
        for(i in 1:nhcrs2){
            polygon(bbmsy2[c(2,3,3,2),i],ffmsy2[c(2,2,3,3),i],
                    col=cols2P[i], border = NA)
        }
        polygon(bbmsyBo[c(2,3,3,2)],ffmsyBo[c(2,2,3,3)],
                col=rgb(t(col2rgb("black"))/255,alpha=0.2), border = NA)
        points(bbmsy1[1,],ffmsy1[1,], cex = cex, pch = 15, col = cols1)
        points(bbmsy2[1,],ffmsy2[1,], cex = cex, pch = 17, col = cols2)
        points(bbmsyBo[1],ffmsyBo[1], cex = cex, pch = 16, col = "black")
        axis(1, cex.axis = 1.1)
        mtext(bquote(B/B[MSY]), 1, 2.9)
        if(spec == 1){
            axis(2, cex.axis = 1.1)
            mtext(bquote(F/F[MSY]), 2, 3.2)
        }
        box()
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, outer=TRUE, line=-1, font=2,
              adj = c(0.12,0.42,0.73)[spec])

        ## Risk-yield short-term
        pblim <- tmp$pblimY3_5[1:3,]
        yield <- tmp$yieldY3_5[1:3,]
        plot(pblim[1,], yield[1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             xlab="",ylab="",
             xlim = limRshort,
             ylim = limYshort)
        abline(v = seq(0,2,0.05), col = "grey80", lty = 2)
        abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        ats <- seq(0,0.25,0.05)
        axis(1, at = ats, cex.axis = 1.1)
        mtext(bquote("Risk [years 3-5]"), 1, 2.8)
        if(spec == 1){
            axis(2, cex.axis = 1.1)
            mtext(bquote("Rel. yield [years 3-5]"), 2, 3.6)
        }
        if(spec == 1) axis(2, cex.axis = 1.1)
        ## 95% CI
        polygon(pblim[c(2,3,3,2),hcrref],yield[c(2,2,3,3),hcrref],
                col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
        for(i in 1:nhcrs1){
            polygon(pblim[c(2,3,3,2),hcrs1[i]],yield[c(2,2,3,3),hcrs1[i]],
                    col=cols1P[i], border = NA)
        }
        for(i in 1:nhcrs2){
            polygon(pblim[c(2,3,3,2),hcrs2[i]],yield[c(2,2,3,3),hcrs2[i]],
                    col=cols2P[i], border = NA)
        }
        polygon(pblim[c(2,3,3,2),hcrBoth],yield[c(2,2,3,3),hcrBoth],
                col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)

        ## HCR type connectors
        lines(pblim[1,c(hcrBoth,hcrs1)], yield[1,c(hcrBoth,hcrs1)], col="grey10", lty=3)
        lines(pblim[1,c(hcrBoth,hcrs2)], yield[1,c(hcrBoth,hcrs2)], col="grey10", lty=3)
        ## points
        ##            points(pblim[1,1], yield[1,1], pch = 8, col=1, cex=cex)
        points(pblim[1,hcrs1], yield[1,hcrs1], pch = 15, col = cols1, cex=cex)
        points(pblim[1,hcrs2], yield[1,hcrs2], pch = 17, col = cols2, cex=cex)
        points(pblim[1,hcrBoth], yield[1,hcrBoth], pch = 16, col=1, cex=cex)
        points(pblim[1,hcrref], yield[1,hcrref], pch = 8, col=1, cex=cex)
        box()

        ## Risk-yield long-term
        pblim <- tmp$pblimY6_20[1:3,]
        yield <- tmp$yieldY6_20[1:3,]
        plot(pblim[1,], yield[1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             xlab="",ylab="",
             xlim = limRlong,
             ylim = limYlong)
        abline(v = seq(0,2,0.05), col = "grey80", lty = 2)
        abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        ats <- seq(0,0.25,0.05)
        axis(1, at = ats, cex.axis = 1.1)
        mtext(bquote("Risk [years 6-20]"), 1, 2.8)
        if(spec == 1){
            axis(2, cex.axis = 1.1)
            mtext(bquote("Rel. yield [years 6-20]"), 2, 3.6)
        }
        if(spec == 1) axis(2, cex.axis = 1.1)
        ## 95% CI
        polygon(pblim[c(2,3,3,2),hcrref],yield[c(2,2,3,3),hcrref],
                col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
        for(i in 1:nhcrs1){
            polygon(pblim[c(2,3,3,2),hcrs1[i]],yield[c(2,2,3,3),hcrs1[i]],
                    col=cols1P[i], border = NA)
        }
        for(i in 1:nhcrs2){
            polygon(pblim[c(2,3,3,2),hcrs2[i]],yield[c(2,2,3,3),hcrs2[i]],
                    col=cols2P[i], border = NA)
        }
        polygon(pblim[c(2,3,3,2),hcrBoth],yield[c(2,2,3,3),hcrBoth],
                col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)

        ## HCR type connectors
        lines(pblim[1,c(hcrBoth,hcrs1)], yield[1,c(hcrBoth,hcrs1)], col="grey10", lty=3)
        lines(pblim[1,c(hcrBoth,hcrs2)], yield[1,c(hcrBoth,hcrs2)], col="grey10", lty=3)
        ## points
        ##            points(pblim[1,1], yield[1,1], pch = 8, col=1, cex=cex)
        points(pblim[1,hcrs1], yield[1,hcrs1], pch = 15, col = cols1, cex=cex)
        points(pblim[1,hcrs2], yield[1,hcrs2], pch = 17, col = cols2, cex=cex)
        points(pblim[1,hcrBoth], yield[1,hcrBoth], pch = 16, col=1, cex=cex)
        points(pblim[1,hcrref], yield[1,hcrref], pch = 8, col=1, cex=cex)
        box()

        ## AAV
        yieldcv <- tmp$yieldDiff[1:3,]
        hcrsUse <- c(hcrBoth,hcrs1,hcrs2)
        xplot <- seq(c(hcrBoth, hcrs1, hcrs2))
        xplotB <- xplot[hcrsUse %in% hcrBoth]
        xplot1 <- xplot[hcrsUse %in% hcrs1]
        xplot2 <- xplot[hcrsUse %in% hcrs2]
        plot(xplot, yieldcv[1,hcrsUse], ty = 'n',
             xaxt="n", yaxt= "n",
             xlab="",ylab="",ylim = limA)
        ##            xplot <- seq(hcrs) + c(-0.075,-0.025,0.025,0.075)[scen]
        abline(h = yieldcv[1,hcrref],lty=2,col="grey80",lwd=2)
        abline(h = yieldcv[1,hcrBoth],lty=1,col="grey80",lwd=2)
        lines(xplot1, yieldcv[1,hcrs1], col = tail(cols1,1), lty=2)
        lines(xplot2, yieldcv[1,hcrs2], col = tail(cols2,1), lty=2)
        segments(xplotB,yieldcv[2,hcrBoth],xplotB,yieldcv[3,hcrBoth],
                 col = 1, lwd = 2)
        segments(xplot1,yieldcv[2,hcrs1],xplot1,yieldcv[3,hcrs1],
                 col = cols1, lwd = 2)
        segments(xplot2,yieldcv[2,hcrs2],xplot2,yieldcv[3,hcrs2],
                 col = cols2, lwd = 2)
        points(xplotB, yieldcv[1,hcrBoth], col = 1,
               pch = 16, cex = cex)
        points(xplot1, yieldcv[1,hcrs1], col = cols1,
               pch = 15, cex = cex)
        points(xplot2, yieldcv[1,hcrs2], col = cols2,
               pch = 17, cex = cex)
##n        axis(1, at = xplot, labels = hcrsAll[hcrsUse], cex.axis = 1.1, las=2)
        mtext(bquote("HCR"), 1, 1)
        if(spec == 1){
            axis(2, cex.axis = 1.1)
            mtext(bquote("AAV"), 2, 3.5)
        }
    }
    ## legend
    par(mar=c(0,0,0,0))
    plot.new()
    hcrsLegend <- c(hcrref,hcrBoth,hcrs1,hcrs2)
    pchLegend <- c(8,16,rep(15,nhcrs1),rep(17,nhcrs2))
    colLegend <- c(1,1,cols1,cols2)
    legend("center", legend=hcrsAll[hcrsLegend], title.adj=0.1,
           pch=pchLegend,
           title = "HCR",
           col=colLegend,
           y.intersp=1.2,
           ncol = 1, bty="n", cex=1.2)
}


## plotting functions
plot.risk.yield <- function(x,
                            showCI = FALSE, blackwhite = FALSE, hcrref = c(1),
                            hcrBoth = 3,
                            hcrs1 = c(46:37,70:64),
                            hcrs2 = c(11:19), showref = TRUE,
                            cols1 = NULL,
                            cols2 = NULL){## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)

    ## plot character color
    if(is.null(cols1) && is.null(cols2)){
        if(blackwhite){
            cols1 <- paste0(rep("grey",nhcrs1),round(seq(10,75,length.out = nhcrs1)))
            cols2 <- paste0(rep("grey",nhcrs2),round(seq(10,75,length.out = nhcrs2)))
        }else{
            ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
            colours()
            cols1 <- colorRampPalette(c("darkorange4","darkorange1"))(nhcrs1)
            cols2 <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(nhcrs2)
        }
        ## cols1P <- rgb(t(col2rgb("darkorange3"))/255,alpha=0.2)
        ## cols2P <- rgb(t(col2rgb("dodgerblue3"))/255,alpha=0.2)
    }
    cols1P <- sapply(cols1, function(x) rgb(t(col2rgb(x))/255,alpha=0.2))
    cols2P <- sapply(cols2, function(x) rgb(t(col2rgb(x))/255,alpha=0.2))
    ## plot character size
    cex <- 1.7
    legendwidth <- 0.4

    hcrsUse <- c(hcrBoth,hcrs1,hcrs2)
    if(showref) hcrsUse <- c(hcrref,hcrsUse)

    xlim <- c(0,1.05) * range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1:3,hcrsUse]))))
    ylim <- c(0,1.05) * range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldTot[1:3,hcrsUse]))))


    ly <- layout(matrix(c(1:12,rep(13,4)),
                        byrow = FALSE,ncol=4,nrow=4),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        nscen <- length(x[[spec]])
        xlim <- c(0,1.05) * range(unlist(lapply(x[[spec]], function(z) z$pblimTot[1:3,hcrsUse])))
        for(scen in 1:nscen){
            pblim <- x[[spec]][[scen]]$pblimTot[1:3,]
            yield <- x[[spec]][[scen]]$yieldTot[1:3,]
            samplesize <- x[[spec]][[scen]]$sampleSize
            if(spec == 3) par(mar=c(0,0,0,1)) else par(mar=c(0,0,0,0))
            plot(pblim[1,], yield[1,],
                 type = "n",
                 xaxt = "n", yaxt = "n",
                 xlim = xlim,
                 ylim = ylim)
            abline(v = seq(0,2,0.1), col = "grey80", lty = 2)
            abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
            if(scen %in% c(4)) axis(1, cex.axis = 1.1)
            if(spec == 1 &&  scen %in% c(1:4)) axis(2, cex.axis = 1.1)
            ##
            ##xlab="Risk", ylab="log Yield")
            ## 95% CI
            if(!is.null(showCI) && showCI == "segments"){
                ## segments(pblim[2,1],yield[1,1],pblim[3,1],yield[1,1], col=1, lty=2)
                segments(pblim[2,indA[1]],yield[1,indA[1]],pblim[3,indA[1]],yield[1,indA[1]], col=1, lty=2)
                segments(pblim[2,indA[-1]],yield[1,indA[-1]],pblim[3,indA[-1]],yield[1,indA[-1]], col=cols1, lty=2)
                segments(pblim[2,indC[-1]],yield[1,indC[-1]],pblim[3,indC[-1]],yield[1,indC[-1]], col=cols2, lty=2)
                ## segments(pblim[1,1],yield[2,1],pblim[1,1],yield[3,1], col=1, lty=2)
                segments(pblim[1,indA[1]],yield[2,indA[1]],pblim[1,indA[1]],yield[3,indA[1]], col=1, lty=2)
                segments(pblim[1,indA[-1]],yield[2,indA[-1]],pblim[1,indA[-1]],yield[3,indA[-1]], col=cols1, lty=2)
                segments(pblim[1,indC[-1]],yield[2,indC[-1]],pblim[1,indC[-1]],yield[3,indC[-1]], col=cols2, lty=2)
            }
            if(!is.null(showCI) && showCI == "polygon"){
                ## polygon(pblim[c(2,3,3,2),1], yield[c(2,2,3,3),1],
                ##         col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
                if(showref) polygon(pblim[c(2,3,3,2),hcrref],yield[c(2,2,3,3),hcrref],
                                    col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
                for(i in 1:nhcrs1){
                    polygon(pblim[c(2,3,3,2),hcrs1[i]],yield[c(2,2,3,3),hcrs1[i]],
                            col=cols1P[i], border = NA)
                }
                for(i in 1:nhcrs2){
                    polygon(pblim[c(2,3,3,2),hcrs2[i]],yield[c(2,2,3,3),hcrs2[i]],
                            col=cols2P[i], border = NA)
                }
                polygon(pblim[c(2,3,3,2),hcrBoth],yield[c(2,2,3,3),hcrBoth],
                        col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
            }
            ## HCR type connectors
            lines(pblim[1,c(hcrBoth,hcrs1)], yield[1,c(hcrBoth,hcrs1)], col="grey10", lty=3)
            lines(pblim[1,c(hcrBoth,hcrs2)], yield[1,c(hcrBoth,hcrs2)], col="grey10", lty=3)
            ## points
            ##            points(pblim[1,1], yield[1,1], pch = 8, col=1, cex=cex)
            points(pblim[1,hcrs1], yield[1,hcrs1], pch = 15, col = cols1, cex=cex)
            points(pblim[1,hcrs2], yield[1,hcrs2], pch = 17, col = cols2, cex=cex)
            points(pblim[1,hcrBoth], yield[1,hcrBoth], pch = 16, col=1, cex=cex)
            if(showref) points(pblim[1,hcrref], yield[1,hcrref], pch = 8, col=1, cex=cex)
            usr <- par("usr")
            txt <- bquote("N"~"="~.(samplesize))
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.01
            ## sample size
            xt <- 0.86 * usr[2]
            yt <- -0.7 * usr[3]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz - 0.01,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz + 0.01,
                col = "white"
            )
            text(xt, yt, txt, font=1, cex=1.1)
            ## plot title
            txt <- paste0(c("A","H","L")[spec],scen)
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.01
            xt <- -1 * usr[1]
            yt <- 0.94 * usr[4]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz - 0.01,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz + 0.01,
                col = "grey40"
            )
            text(xt, yt, txt, font=2, cex=1.2, col = "white")
            box()
        }
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, outer=TRUE, line=1.2, font=2,
              adj = c(0.136,0.466,0.788)[spec])
    }
    ## Axes labels
    mtext("Rel. yield", 2, outer=TRUE, line=4)
    mtext("Risk", 1, outer=TRUE, line=3.5, adj = 0.48)
    ## legend
    par(mar=c(0,0,0,0))
    plot.new()
    hcrsLegend <- c(hcrBoth,hcrs1,hcrs2)
    pchLegend <- c(16,rep(15,nhcrs1),rep(17,nhcrs2))
    colLegend <- c(1,cols1,cols2)
    if(showref){
        hcrsLegend <- c(hcrref,hcrsLegend)
        pchLegend <- c(8,pchLegend)
        colLegend <- c(1,colLegend)
    }
    legend("center", legend=hcrsAll[hcrsLegend], title.adj=0.1,
           pch=pchLegend,
           title = "HCR",
           col=colLegend,
           y.intersp=1.2,
           ncol = 1, bty="n", cex=1.2)

}

## plotting functions
plot.risk.yield.single <- function(x, species = 1, scenario = 1,
                                   showCI = FALSE, blackwhite = FALSE, hcrref = c(1),
                                   hcrBoth = 3,
                                   hcrs1 = c(46:37,70:64),
                                   hcrs2 = c(11:19), showref = TRUE){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)

    ## plot character color
    if(blackwhite){
        col1 <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = nhcrs1)))
        col2 <- paste0(rep("grey",length(indC[-1])),round(seq(10,75,length.out = nhcrs2)))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        col1 <- colorRampPalette(c("darkorange4","darkorange1"))(nhcrs1)
        col2 <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(nhcrs2)
    }
    col1P <- rgb(t(col2rgb("darkorange3"))/255,alpha=0.2)
    col2P <- rgb(t(col2rgb("dodgerblue3"))/255,alpha=0.2)
    ## plot character size
    cex <- 1.7
    legendwidth <- 0.55
    spec = species
    scen = scenario

    hcrsUse <- c(hcrBoth,hcrs1,hcrs2)
    if(showref) hcrsUse <- c(hcrref,hcrsUse)

    xlim <- c(0,1.05) * range(x[[spec]][[scen]]$pblimTot[1:3,hcrsUse])
    ylim <- c(0,1.05) * range(x[[spec]][[scen]]$yieldTot[1:3,hcrsUse])

    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,2))
    nspec <- length(x)
    pblim <- x[[spec]][[scen]]$pblimTot[1:3,]
    yield <- x[[spec]][[scen]]$yieldTot[1:3,]
    samplesize <- x[[spec]][[scen]]$sampleSize
    ly <- layout(matrix(1:2,
                        byrow = FALSE,ncol=2,nrow=1),
                 widths=c(rep(1,1),legendwidth),
                 heights=rep(1,1))
    opar <- par(mar=c(0,0,0,1),oma=c(5,5.5,3.5,0))
    plot(pblim[1,], yield[1,],
         type = "n",
         xaxt = "n", yaxt = "n",
         xlim = xlim,
         ylim = ylim)
    abline(v = seq(0,2,0.05), col = "grey80", lty = 2)
    abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
    axis(1, cex.axis = 1.1)
    axis(2, cex.axis = 1.1)
    ##xlab="Risk", ylab="log Yield")
    ## 95% CI
    if(!is.null(showCI) && showCI == "segments"){
        ## segments(pblim[2,1],yield[1,1],pblim[3,1],yield[1,1], col=1, lty=2)
        segments(pblim[2,indA[1]],yield[1,indA[1]],pblim[3,indA[1]],yield[1,indA[1]], col=1, lty=2)
        segments(pblim[2,indA[-1]],yield[1,indA[-1]],pblim[3,indA[-1]],yield[1,indA[-1]], col=colA, lty=2)
        segments(pblim[2,indC[-1]],yield[1,indC[-1]],pblim[3,indC[-1]],yield[1,indC[-1]], col=colC, lty=2)
        ## segments(pblim[1,1],yield[2,1],pblim[1,1],yield[3,1], col=1, lty=2)
        segments(pblim[1,indA[1]],yield[2,indA[1]],pblim[1,indA[1]],yield[3,indA[1]], col=1, lty=2)
        segments(pblim[1,indA[-1]],yield[2,indA[-1]],pblim[1,indA[-1]],yield[3,indA[-1]], col=colA, lty=2)
        segments(pblim[1,indC[-1]],yield[2,indC[-1]],pblim[1,indC[-1]],yield[3,indC[-1]], col=colC, lty=2)
    }
    if(!is.null(showCI) && showCI == "polygon"){
        ## polygon(pblim[c(2,3,3,2),1], yield[c(2,2,3,3),1],
        ##         col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
        if(showref) polygon(pblim[c(2,3,3,2),hcrref],yield[c(2,2,3,3),hcrref],
                            col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
        for(i in 1:nhcrs1){
            polygon(pblim[c(2,3,3,2),hcrs1[i]],yield[c(2,2,3,3),hcrs1[i]],
                    col=col1P, border = NA)
        }
        for(i in 1:nhcrs2){
            polygon(pblim[c(2,3,3,2),hcrs2[i]],yield[c(2,2,3,3),hcrs2[i]],
                    col=col2P, border = NA)
        }
        polygon(pblim[c(2,3,3,2),hcrBoth],yield[c(2,2,3,3),hcrBoth],
                            col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
    }
    ## HCR type connectors
    lines(pblim[1,c(hcrBoth,hcrs1)], yield[1,c(hcrBoth,hcrs1)], col="grey10", lty=3)
    lines(pblim[1,c(hcrBoth,hcrs2)], yield[1,c(hcrBoth,hcrs2)], col="grey10", lty=3)
    ## points
    ##            points(pblim[1,1], yield[1,1], pch = 8, col=1, cex=cex)
    points(pblim[1,hcrs1], yield[1,hcrs1], pch = 15, col = col1, cex=cex)
    points(pblim[1,hcrs2], yield[1,hcrs2], pch = 17, col = col2, cex=cex)
    points(pblim[1,hcrBoth], yield[1,hcrBoth], pch = 16, col=1, cex=cex)
    if(showref) points(pblim[1,hcrref], yield[1,hcrref], pch = 8, col=1, cex=cex)
    usr <- par("usr")
    txt <- bquote("N"~"="~.(samplesize))
    sw   <- strwidth(txt)
    sh   <- strheight(txt)
    frsz <- 0.01
    ## sample size
    xt <- 0.86 * usr[2]
    yt <- -0.7 * usr[3]
    rect(
        xt - sw/2 - frsz + 0.008,
        yt - sh/2 - frsz - 0.01,
        xt + sw/2 + frsz - 0.008,
        yt + sh/2 + frsz + 0.01,
        col = "white"
    )
    text(xt, yt, txt, font=1, cex=1.1)
    ## plot title
    txt <- paste0(c("A","H","L")[spec],scen)
    sw   <- strwidth(txt)
    sh   <- strheight(txt)
    frsz <- 0.01
    xt <- -1 * usr[1]
    yt <- 0.94 * usr[4]
    rect(
        xt - sw/2 - frsz + 0.008,
        yt - sh/2 - frsz - 0.01,
        xt + sw/2 + frsz - 0.008,
        yt + sh/2 + frsz + 0.01,
        col = "grey40"
    )
    text(xt, yt, txt, font=2, cex=1.2, col = "white")
    box()
    mtext(c("Anchovy","Haddock","Ling")[species], 3, outer=FALSE, line=1.2, font=2)

    ## Axes labels
    mtext("Rel. yield", 2, outer=TRUE, line=3.5)
    mtext("Risk", 1, outer=TRUE, line=3, adj = 0.33)
    ## legend
    plot.new()
    hcrsLegend <- c(hcrBoth,hcrs1,hcrs2)
    pchLegend <- c(16,rep(15,nhcrs1),rep(17,nhcrs2))
    colLegend <- c(1,col1,col2)
    if(showref){
        hcrsLegend <- c(hcrref,hcrsLegend)
        pchLegend <- c(8,pchLegend)
        colLegend <- c(1,colLegend)
    }
    legend("center", legend=hcrsAll[hcrsLegend], title.adj=0.1,
           pch=pchLegend,
           title = "HCR",
           col=colLegend,
           y.intersp=1.1,
           ncol = 2, bty="n", cex=0.8)

}



## AAV in yield
plot.yield.cv <- function(x, showCI = FALSE, blackwhite = FALSE, outline = TRUE,
                          hcrref = 1, hcrBoth = 3,
                          hcrs1 = c(46:37,70:64),
                          hcrs2 = c(11:19), showref=TRUE, hcrs3=NULL){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)
    nhcrs3 <- length(hcrs3)

    ## plot character color
    if(blackwhite){
        col1 <- paste0(rep("grey",nhcrs1),round(seq(10,75,length.out = nhcrs1)))
        col2 <- paste0(rep("grey",nhcrs2),round(seq(10,75,length.out = nhcrs2)))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        col1 <- colorRampPalette(c("darkorange4","darkorange1"))(nhcrs1)
        col2 <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(nhcrs2)
        col3 <- colorRampPalette(c("purple4","purple1"))(nhcrs3)
    }
    col1P <- rgb(t(col2rgb("darkorange3"))/255,alpha=0.2)
    col2P <- rgb(t(col2rgb("dodgerblue3"))/255,alpha=0.2)
    col3P <- rgb(t(col2rgb("purple3"))/255,alpha=0.2)
    ## plot character size
    cex <- 1.9

    hcrsUse <- c(hcrBoth,hcrs1,hcrs2,hcrs3)
    colsUse <- c(1,col1,col2,col3)
    if(showref){
        hcrsUse <- c(hcrref,hcrsUse)
        colsUse <- c("grey80",colsUse)
    }
    ly <- layout(matrix(c(1:12),
                        byrow = FALSE,ncol=3,nrow=4),
                 widths=c(rep(1,3)),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(7.5,7,3.5,2))
    nscen <- length(x[[1]])
    ylim <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        z$yieldDiff[1:3,hcrsUse]))))

    for(spec in 1:nspec){
        yieldcv <- lapply(x[[spec]], function(y) y$yieldDiff[1:3,])
        for(scen in 1:nscen){
            ##            ylim <- range(yieldcv[[scen]])
            xplot <- seq(hcrsUse)
            xplotB <- xplot[hcrsUse %in% hcrBoth]
            xplot1 <- xplot[hcrsUse %in% hcrs1]
            xplot2 <- xplot[hcrsUse %in% hcrs2]
            xplot3 <- xplot[hcrsUse %in% hcrs3]
            plot(xplot, yieldcv[[scen]][1,hcrsUse], ty = 'n',
                 col = colsUse, xaxt="n", yaxt= "n",
                 xlab="",ylab="",ylim = ylim)
            ##            xplot <- seq(hcrs) + c(-0.075,-0.025,0.025,0.075)[scen]
            if(showref) abline(h = yieldcv[[scen]][1,hcrref],lty=2,col="grey80",lwd=2)
            abline(h = yieldcv[[scen]][1,hcrBoth],lty=1,col="grey80",lwd=2)
            lines(xplot1, yieldcv[[scen]][1,hcrs1], col = tail(col1,1), lty=2)
            lines(xplot2, yieldcv[[scen]][1,hcrs2], col = tail(col2,1), lty=2)
            if(!is.null(hcrs3)) lines(xplot3, yieldcv[[scen]][1,hcrs3], col = tail(col3,1), lty=2)
            segments(xplotB,yieldcv[[scen]][2,hcrBoth],xplotB,yieldcv[[scen]][3,hcrBoth],
                     col = 1, lwd = 2)
            segments(xplot1,yieldcv[[scen]][2,hcrs1],xplot1,yieldcv[[scen]][3,hcrs1],
                     col = col1, lwd = 2)
            segments(xplot2,yieldcv[[scen]][2,hcrs2],xplot2,yieldcv[[scen]][3,hcrs2],
                     col = col2, lwd = 2)
            if(!is.null(hcrs3)) segments(xplot3,yieldcv[[scen]][2,hcrs3],
                                         xplot3,yieldcv[[scen]][3,hcrs3],
                                         col = col3, lwd = 2)
            points(xplotB, yieldcv[[scen]][1,hcrBoth], col = 1,
                   pch = c(15,16,17,18)[scen], cex = cex)
            points(xplot1, yieldcv[[scen]][1,hcrs1], col = col1,
                   pch = c(15,16,17,18)[scen], cex = cex)
            points(xplot2, yieldcv[[scen]][1,hcrs2], col = col2,
                   pch = c(15,16,17,18)[scen], cex = cex)
            if(!is.null(hcrs3)) points(xplot3, yieldcv[[scen]][1,hcrs3], col = col3,
                                       pch = c(15,16,17,18)[scen], cex = cex)
            usr <- par("usr")
            ## plot title
            txt <- paste0(c("A","H","L")[spec],scen)
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.45
            xt <- 1.3 ##6 * usr[1]
            yt <- 0.93 * usr[4]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz + 0.44,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz - 0.44,
                col = "grey40"
            )
            text(xt, yt, txt, font=2, cex=1.2, col = "white")
            if(spec == 1) axis(2)
            if(scen == 4) axis(1,at=xplot,labels = hcrsAll[hcrsUse], las=2)
            if(scen == 1) mtext(c("Anchovy","Haddock","Ling")[spec],font=2,line=0.5)
            box()
        }
    }

    ## Axes labels
    mtext("AAV in yield", 2, outer=TRUE, line=4)
    mtext("HCRs", 1, outer=TRUE, line=6)
}


plot.diff.scen <- function(x, scenarios = c(1,3),
                           hcrBoth = NULL,
                           hcrref = 1, showref=TRUE,
                           hcrs1 = NULL,
                           hcrs2 = NULL,
                           cols1 = NULL, cols2 = NULL){

    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)

    cols1P <- sapply(cols1, function(x) rgb(t(col2rgb(x))/255,alpha=0.2))
    cols2P <- sapply(cols2, function(x) rgb(t(col2rgb(x))/255,alpha=0.2))
    cex <- 1.9

    ## plot character size
    cex <- 1.7
    legendwidth <- 0.5

    hcrsUse <- c(hcrBoth,hcrs1,hcrs2)
    if(showref) hcrsUse <- c(hcrref,hcrsUse)

    ly <- layout(matrix(c(1:15,rep(16,5)),
                        byrow = FALSE,ncol=4,nrow=5),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,5))
    opar <- par(mar=c(0,0,0,0),oma=c(7,7,4.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        for(q in 1:5){
            quant <- c("pblimY3_5","pblimY6_20",
                       "yieldY3_5","yieldY6_20","yieldDiff")[q]
            tmp1 <- x[[spec]][[scenarios[1]]][[quant]]
            tmp2 <- x[[spec]][[scenarios[2]]][[quant]]
            lim <- c(0.1,1.1) * range(tmp1[1:3,hcrsUse],tmp2[1:3,hcrsUse])
            if(spec == 3) par(mar=c(0,0,0,0.7)) else par(mar=c(0,0,0,0))
            plot(tmp1[1,hcrs1],
                 tmp2[1,hcrs1],
                 ty='n', ylim=lim, xlim = lim,
                 ylab = "", xlab = "", xaxt = "n", yaxt = "n")
            abline(0,1,lty=2, lwd = 1.5)
            mtext(c("Anchovy","Haddock","Ling")[spec], 3, outer=TRUE, line=1, font=2,
                  adj = c(0.12,0.42,0.73)[spec])
            if(spec == 1){
                if(q == 1){
                    mtext(bquote("Risk [years 3-5]"), 2, 1.5)
                }else if(q == 2){
                    mtext(bquote("Risk [years 6-20]"), 2, 1.5)
                }else if(q == 3){
                    mtext(bquote("Rel. yield [years 3-5]"), 2, 1.5)
                }else if(q == 4){
                    mtext(bquote("Rel. yield [years 6-20]"), 2, 1.5)
                }else if(q == 5){
                    mtext(bquote("AAV"), 2, 1.5)
                }
            }
            ## ref
            if(showref) polygon(tmp1[c(2,3,3,2),hcrref],tmp2[c(2,2,3,3),hcrref],
                                col=rgb(t(col2rgb("black"))/255,alpha=0.15),
                                border = NA)
            ## both
            polygon(tmp1[c(2,3,3,2),hcrBoth],tmp2[c(2,2,3,3),hcrBoth],
                    col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
            for(i in 1:nhcrs1){
                polygon(tmp1[c(2,3,3,2),hcrs1[i]],tmp2[c(2,2,3,3),hcrs1[i]],
                        col=cols1P[i], border = NA)
            }
            for(i in 1:nhcrs1){
                polygon(tmp1[c(2,3,3,2),hcrs2[i]],tmp2[c(2,2,3,3),hcrs2[i]],
                        col=cols2P[i], border = NA)
            }
            if(showref) points(tmp1[1,hcrref],tmp2[1,hcrref],
                               pch = 8, cex=cex)
            points(tmp1[1,hcrBoth],tmp2[1,hcrBoth],
                   pch = 16, cex=cex)
            points(tmp1[1,hcrs1],tmp2[1,hcrs1],
                   pch = 15, cex=cex, col = cols1)
            points(tmp1[1,hcrs2],tmp2[1,hcrs2],
                   pch = 17, cex=cex, col = cols2)
            lines(tmp1[1,hcrs1],tmp2[1,hcrs1], lty = 3, lwd = 1,
                  col = cols1[1])
            lines(tmp1[1,hcrs2],tmp2[1,hcrs2], lty = 3, lwd = 1,
                  col = cols2[1])
            box()
            if(q == 5) mtext(bquote("Respective scaled quantity"), 1, 1.5)
        }
    }
    mtext(paste0("Scenario ", scenarios[1]), 1, 4, outer = TRUE, font = 2, adj=0.45)
    mtext(paste0("Scenario ", scenarios[2]), 2, 4.5, outer = TRUE, font = 2)
    ## legend
    par(mar=c(0,0,0,0))
    plot.new()
    hcrsLegend <- c(hcrref,hcrBoth,hcrs1,hcrs2)
    pchLegend <- c(8,16,rep(15,nhcrs1),rep(17,nhcrs2))
    colLegend <- c(1,1,cols1,cols2)
    legend("center", legend=hcrsAll[hcrsLegend], title.adj=0.1,
           pch=pchLegend,
           title = "HCR",
           col=colLegend,
           y.intersp=1.2,
           ncol = 1, bty="n", cex=1.2)
}



plot.diff.scen.old <- function(x,
                           showCI = FALSE, blackwhite = FALSE,
                           hcrref = 1, hcrBoth = 3,
                           hcrs1 = c(46:37,70:64),
                           hcrs2 = c(11:19), showref=TRUE, hcrs3=NULL){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)
    nhcrs3 <- length(hcrs3)

    ## plot character color
    if(blackwhite){
        col1 <- paste0(rep("grey",nhcrs1),round(seq(10,75,length.out = nhcrs1)))
        col2 <- paste0(rep("grey",nhcrs2),round(seq(10,75,length.out = nhcrs2)))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        col1 <- colorRampPalette(c("darkorange4","darkorange1"))(nhcrs1)
        col2 <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(nhcrs2)
        col3 <- colorRampPalette(c("purple4","purple1"))(nhcrs3)
    }
    col1P <- rgb(t(col2rgb("darkorange3"))/255,alpha=0.2)
    col2P <- rgb(t(col2rgb("dodgerblue3"))/255,alpha=0.2)
    col3P <- rgb(t(col2rgb("purple3"))/255,alpha=0.2)
    ## plot character size
    cex <- 1.9

    hcrsUse <- c(hcrBoth,hcrs1,hcrs2,hcrs3)
    colsUse <- c(1,col1,col2,col3)
    if(showref){
        hcrsUse <- c(hcrref,hcrsUse)
        colsUse <- c("grey80",colsUse)
    }

    nscen <- length(x[[1]])

    risk <- lapply(x, function(y) do.call(cbind,lapply(y[1:3], function(z) z$pblimTot[1,])))
    yield <- lapply(x, function(y) do.call(cbind,lapply(y[1:3], function(z) z$yieldTot[1,])))
    aav <- lapply(x, function(y) do.call(cbind,lapply(y[1:3], function(z) z$yieldDiff[1,])))

    ylimR <- range(unlist(lapply(risk,function(y) y[hcrsUse,])))
    ylimY <- range(unlist(lapply(yield,function(y) y[hcrsUse,])))
    ylimA <- range(unlist(lapply(aav,function(y) y[hcrsUse,])))



    ly <- layout(matrix(c(1:12),
                        byrow = FALSE,ncol=3,nrow=4),
                 widths=c(rep(1,3)),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(7.5,7,3.5,2))
    for(spec in 1:nspec){
        xplot <- 1:3
        plot(xplot, risk[[spec]][1,], ty = 'n',
             col = colsUse, xaxt="n", yaxt= "n",
             xlab="",ylab="",ylim = ylimR)
        if(showref) lines(xplot, risk[[spec]][hcrref,], col = "grey60", lty=2,lwd=2)
        for(i in 1:nhcrs1){
            lines(xplot, risk[[spec]][hcrs1[i],], col = col1[i], lty=1,lwd=2)
        }
        for(i in 1:nhcrs2){
            lines(xplot, risk[[spec]][hcrs2[i],], col = col2[i], lty=1,lwd=2)
        }
        c
        lines(xplot1, yieldcv[[scen]][1,hcrs1], col = tail(col1,1), lty=2)
        lines(xplot2, yieldcv[[scen]][1,hcrs2], col = tail(col2,1), lty=2)
        if(!is.null(hcrs3)) lines(xplot3, yieldcv[[scen]][1,hcrs3], col = tail(col3,1), lty=2)
        segments(xplotB,yieldcv[[scen]][2,hcrBoth],xplotB,yieldcv[[scen]][3,hcrBoth],
                 col = 1, lwd = 2)
        segments(xplot1,yieldcv[[scen]][2,hcrs1],xplot1,yieldcv[[scen]][3,hcrs1],
                 col = col1, lwd = 2)
        segments(xplot2,yieldcv[[scen]][2,hcrs2],xplot2,yieldcv[[scen]][3,hcrs2],
                 col = col2, lwd = 2)
        if(!is.null(hcrs3)) segments(xplot3,yieldcv[[scen]][2,hcrs3],
                                     xplot3,yieldcv[[scen]][3,hcrs3],
                                     col = col3, lwd = 2)
        points(xplotB, yieldcv[[scen]][1,hcrBoth], col = 1,
               pch = c(15,16,17,18)[scen], cex = cex)
        points(xplot1, yieldcv[[scen]][1,hcrs1], col = col1,
               pch = c(15,16,17,18)[scen], cex = cex)
        points(xplot2, yieldcv[[scen]][1,hcrs2], col = col2,
               pch = c(15,16,17,18)[scen], cex = cex)
        if(!is.null(hcrs3)) points(xplot3, yieldcv[[scen]][1,hcrs3], col = col3,
                                   pch = c(15,16,17,18)[scen], cex = cex)
        usr <- par("usr")
        ## plot title
        txt <- paste0(c("A","H","L")[spec],scen)
        sw   <- strwidth(txt)
        sh   <- strheight(txt)
        frsz <- 0.45
        xt <- 1.3 ##6 * usr[1]
        yt <- 0.93 * usr[4]
        rect(
            xt - sw/2 - frsz,
            yt - sh/2 - frsz + 0.44,
            xt + sw/2 + frsz,
            yt + sh/2 + frsz - 0.44,
            col = "grey40"
        )
        text(xt, yt, txt, font=2, cex=1.2, col = "white")
        if(spec == 1) axis(2)
        if(scen == 4) axis(1,at=xplot,labels = hcrsAll[hcrsUse], las=2)
        if(scen == 1) mtext(c("Anchovy","Haddock","Ling")[spec],font=2,line=0.5)

    }

    ## Axes labels
    mtext("AAV in yield", 2, outer=TRUE, line=4)
    mtext("HCRs", 1, outer=TRUE, line=6)
}



plot.diff.perc <- function(x, scen = 1, hcrs1 = 20:29, hcrs2 = NULL,
                           cols1 = colorRampPalette(c("dodgerblue4","dodgerblue1"))(9),
                           cols2 = NULL){
    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)
    nhcrs <- nhcrs1-1 + ifelse(is.null(hcrs2), 0, nhcrs2-1)
    hcrsUse <- c(hcrs1[-1],hcrs2[-1])

    ylimR <- range(unlist(lapply(x, function(y) c(
                                                    abs(diff(y[[scen]]$pblimTot[1,hcrs1])),
                                                    abs(diff(y[[scen]]$pblimTot[1,hcrs2]))
))))
    ylimY <- range(unlist(lapply(x, function(y) c(
                                                    abs(diff(y[[scen]]$yieldTot[1,hcrs1])),
                                                    abs(diff(y[[scen]]$yieldTot[1,hcrs2]))
))))
    ylimA <- range(unlist(lapply(x, function(y) c(
                                                    abs(diff(y[[scen]]$yieldDiff[1,hcrs1])),
                                                    abs(diff(y[[scen]]$yieldDiff[1,hcrs2]))
                                                  ))))

    cols <- c(cols1,cols2)


    ly <- layout(matrix(c(1:9),
                        byrow = FALSE,ncol=3,nrow=3),
                 widths=c(rep(1,3)),
                 heights=rep(1,3))
    opar <- par(mar=c(0,0,0,0),oma=c(7.5,7,3.5,2))
    for(spec in 1:nspec){
        risk <- c(abs(diff(x[[spec]][[scen]]$pblimTot[1,hcrs1])),
                  abs(diff(x[[spec]][[scen]]$pblimTot[1,hcrs2])))
        yield <- c(abs(diff(x[[spec]][[scen]]$yieldTot[1,hcrs1])),
                   abs(diff(x[[spec]][[scen]]$yieldTot[1,hcrs2])))
        aav <- c(abs(diff(x[[spec]][[scen]]$yieldDiff[1,hcrs1])),
                 abs(diff(x[[spec]][[scen]]$yieldDiff[1,hcrs2])))
        xplot <- seq(nhcrs)
        plot(xplot, risk, ty='n',
             xaxt = "n", yaxt = "n",
             xlab = "", ylab = "",
             xlim = range(xplot) + c(-0.2,0.2),
             ylim = ylimR)
        for(i in 1:length(xplot)){
            rect(xplot[i]-0.4,
                 -1,
                 xplot[i]+0.4,
                 risk[i],
                 col = cols[i]
                 )
        }
        if(spec == 1){
            axis(2)
            mtext(bquote(Delta ~ "Risk"), 2, 3)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec],3,1,font=2)
        box()
        ## yield
        plot(xplot, yield, ty='n',
             xaxt = "n", yaxt = "n",
             xlab = "", ylab = "",
             xlim = range(xplot) + c(-0.2,0.2),
             ylim = ylimY)
        for(i in 1:length(xplot)){
            rect(xplot[i]-0.4,
                 -1,
                 xplot[i]+0.4,
                 yield[i],
                 col = cols[i]
                 )
        }
        if(spec == 1){
            axis(2)
            mtext(bquote(Delta ~ "Yield"), 2, 3)
        }
        box()
        ## aav
        plot(xplot, aav, ty='n',
             xaxt = "n", yaxt = "n",
             xlab = "", ylab = "",
             xlim = range(xplot) + c(-0.2,0.2),
             ylim = ylimA)
        for(i in 1:length(xplot)){
            rect(xplot[i]-0.4,
                 -1,
                 xplot[i]+0.4,
                 aav[i],
                 col = cols[i]
                 )
        }
        axis(1, at = xplot, labels = hcrsAll[hcrsUse])
        if(spec == 1){
            axis(2)
            mtext(bquote(Delta ~ "AAV"), 2, 3)
        }
        box()
    }
    mtext("HCR", 1, 3, outer = TRUE)
}




## plotting functions
plot.risk.yield.q2 <- function(x,
                                showCI = FALSE, blackwhite = FALSE, hcrref = c(1),
                               hcrs1 = c(3,11:19), hcrs2 = c(20,28:36), hcrs3 = c(47,55:63),
                               showref = TRUE){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)
    nhcrs3 <- length(hcrs3)

    ## plot character color
    if(blackwhite){
        col1 <- paste0(rep("grey",nhcrs1),round(seq(10,75,length.out = nhcrs1)))
        col2 <- paste0(rep("grey",nhcrs2),round(seq(10,75,length.out = nhcrs2)))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        col1 <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(nhcrs1)
        col2 <- colorRampPalette(c("darkgreen","seagreen1"))(nhcrs2)
        col3 <- colorRampPalette(c("purple4","purple1"))(nhcrs3)
    }
    col1P <- rgb(t(col2rgb("dodgerblue3"))/255,alpha=0.2)
    col2P <- rgb(t(col2rgb("darkgreen"))/255,alpha=0.2)
    col3P <- rgb(t(col2rgb("purple3"))/255,alpha=0.2)
    ## plot character size
    cex <- 1.7
    legendwidth <- 0.4

    hcrsUse <- c(hcrs1,hcrs2,hcrs3)
    if(showref) hcrsUse <- c(hcrref,hcrsUse)

    xlim <- c(0,1.05) * range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1:3,hcrsUse]))))
    ylim <- c(0,1.05) * range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldTot[1:3,hcrsUse]))))


    ly <- layout(matrix(c(1:12,rep(13,4)),
                        byrow = FALSE,ncol=4,nrow=4),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,2))
    nspec <- length(x)
    for(spec in 1:nspec){
        nscen <- length(x[[spec]])
        xlim <- c(0,1.05) * range(unlist(lapply(x[[spec]], function(z) z$pblimTot[1:3,hcrsUse])))
        for(scen in 1:nscen){
            pblim <- x[[spec]][[scen]]$pblimTot[1:3,]
            yield <- x[[spec]][[scen]]$yieldTot[1:3,]
            samplesize <- x[[spec]][[scen]]$sampleSize
            plot(pblim[1,], yield[1,],
                 type = "n",
                 xaxt = "n", yaxt = "n",
                 xlim = xlim,
                 ylim = ylim)
            abline(v = seq(0,2,0.1), col = "grey80", lty = 2)
            abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
            if(scen %in% c(4)) axis(1, cex.axis = 1.1)
            if(spec == 1 &&  scen %in% c(1:4)) axis(2, cex.axis = 1.1)
            ##
            ##xlab="Risk", ylab="log Yield")
            ## 95% CI
            if(!is.null(showCI) && showCI == "segments"){
                ## segments(pblim[2,1],yield[1,1],pblim[3,1],yield[1,1], col=1, lty=2)
                segments(pblim[2,indA[1]],yield[1,indA[1]],pblim[3,indA[1]],yield[1,indA[1]], col=1, lty=2)
                segments(pblim[2,indA[-1]],yield[1,indA[-1]],pblim[3,indA[-1]],yield[1,indA[-1]], col=colA, lty=2)
                segments(pblim[2,indC[-1]],yield[1,indC[-1]],pblim[3,indC[-1]],yield[1,indC[-1]], col=colC, lty=2)
                ## segments(pblim[1,1],yield[2,1],pblim[1,1],yield[3,1], col=1, lty=2)
                segments(pblim[1,indA[1]],yield[2,indA[1]],pblim[1,indA[1]],yield[3,indA[1]], col=1, lty=2)
                segments(pblim[1,indA[-1]],yield[2,indA[-1]],pblim[1,indA[-1]],yield[3,indA[-1]], col=colA, lty=2)
                segments(pblim[1,indC[-1]],yield[2,indC[-1]],pblim[1,indC[-1]],yield[3,indC[-1]], col=colC, lty=2)
            }
            if(!is.null(showCI) && showCI == "polygon"){
                ## polygon(pblim[c(2,3,3,2),1], yield[c(2,2,3,3),1],
                ##         col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
                if(showref) polygon(pblim[c(2,3,3,2),hcrref],yield[c(2,2,3,3),hcrref],
                                    col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
                for(i in 1:nhcrs1){
                    polygon(pblim[c(2,3,3,2),hcrs1[i]],yield[c(2,2,3,3),hcrs1[i]],
                            col=col1P, border = NA)
                }
                for(i in 1:nhcrs2){
                    polygon(pblim[c(2,3,3,2),hcrs2[i]],yield[c(2,2,3,3),hcrs2[i]],
                            col=col2P, border = NA)
                }
                for(i in 1:nhcrs3){
                    polygon(pblim[c(2,3,3,2),hcrs3[i]],yield[c(2,2,3,3),hcrs3[i]],
                            col=col3P, border = NA)
                }
            }
            ## HCR type connectors
            lines(pblim[1,hcrs1], yield[1,hcrs1], col="grey10", lty=3)
            lines(pblim[1,hcrs2], yield[1,hcrs2], col="grey10", lty=3)
            lines(pblim[1,hcrs3], yield[1,hcrs3], col="grey10", lty=3)
            ## points
            ##            points(pblim[1,1], yield[1,1], pch = 8, col=1, cex=cex)
            if(showref) points(pblim[1,hcrref], yield[1,hcrref], pch = 8, col=1, cex=cex)
            points(pblim[1,hcrs1], yield[1,hcrs1], pch = 15, col = col1, cex=cex)
            points(pblim[1,hcrs2], yield[1,hcrs2], pch = 17, col = col2, cex=cex)
            points(pblim[1,hcrs3], yield[1,hcrs3], pch = 16, col = col3, cex=cex)
            usr <- par("usr")
            txt <- bquote("N"~"="~.(samplesize))
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.01
            ## sample size
            xt <- 0.86 * usr[2]
            yt <- -0.7 * usr[3]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz - 0.01,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz + 0.01,
                col = "white"
            )
            text(xt, yt, txt, font=1, cex=1.1)
            ## plot title
            txt <- paste0(c("A","H","L")[spec],scen)
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.01
            xt <- -1 * usr[1]
            yt <- 0.94 * usr[4]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz - 0.01,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz + 0.01,
                col = "grey40"
            )
            text(xt, yt, txt, font=2, cex=1.2, col = "white")
            box()
        }
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, outer=TRUE, line=1.2, font=2,
              adj = c(0.136,0.466,0.788)[spec])
    }
    ## Axes labels
    mtext("Rel. yield", 2, outer=TRUE, line=4)
    mtext("Risk", 1, outer=TRUE, line=3.5, adj = 0.48)
    ## legend
    plot.new()
    hcrsLegend <- c(hcrs1,hcrs2,hcrs3)
    pchLegend <- c(rep(15,nhcrs1),rep(17,nhcrs2),rep(16,nhcrs3))
    colLegend <- c(col1,col2,col3)
    if(showref){
        hcrsLegend <- c(hcrref,hcrsLegend)
        pchLegend <- c(8,pchLegend)
        colLegend <- c(1,colLegend)
    }
    legend("center", legend=hcrsAll[hcrsLegend], title.adj=0.1,
           pch=pchLegend,
           title = "HCR",
           col=colLegend,
           y.intersp=1.2,
           ncol = 1, bty="n", cex=1.2)

}



## AAV in yield as boxplots (but considerable variation between scenarios)
plot.yield.cv.box <- function(x, showCI = FALSE, blackwhite = FALSE, outline = TRUE, hcrs = c(1,3:19)){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    ## plot character color
    if(blackwhite){
        colA <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        colC <- paste0(rep("grey",length(indC)),round(seq(10,75,length.out = length(indC))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        colA <- colorRampPalette(c("darkorange4","darkorange1"))(length(indA[-1]))
        colC <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(indC))
    }
    cols <- c("grey60","grey30",colA,colC)
    ## plot character size
    cex <- 1.7

    opar <- par(mfrow=c(3,1), mar=c(2,2,2,2),oma=c(6.5,4,2,0.5))
    nscen <- length(x[[1]])
    ylim <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$yieldDiff[1,hcrs]))))

    for(spec in 1:nspec){
        yieldcv <- do.call(rbind, lapply(x[[spec]], function(y) y$yieldDiff[1,hcrs]))
        boxplot(yieldcv, col = cols, xaxt="n", outline = outline)
        if(spec == 3) axis(1,at=1:18,labels = hcrsAll[hcrs], las=2)
        mtext(c("Anchovy","Haddock","Ling")[spec],font=2,line=0.5)
    }

    ## Axes labels
    mtext("AAV in yield", 2, outer=TRUE, line=1.5)
    mtext("HCRs", 1, outer=TRUE, line=4.5)
}


## Plot rel error in spict est. quantities and SD of quantities
plot.rel.err.sd <- function(x, hcrs = c(3,4)){

    ylimB <- c(0.93,1.07) * range(unlist(lapply(x,function(y) lapply(y, function(z) z$bbmsyRETot[1,hcrs]))))
    ylimBsd <- c(0.95,1.07) * range(unlist(lapply(x,function(y) lapply(y, function(z) z$bbmsySDTot[1,hcrs]))))
    ylimF <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z) z$ffmsyRETot[1,hcrs]))))
    ylimFsd <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z) z$ffmsySDTot[1,hcrs]))))
    ylimT <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z) z$tacRETot[1,hcrs]))))
    ylimCsd <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z) z$cpSDTot[1,hcrs]))))

##    ylimB <- ylimBsd <- ylimF <- ylimFsd <- ylimT <- ylimCsd <- c(0,3.6)

    scens <- 1:4
    nhcr <- length(hcrs)

    col1 <- "dodgerblue2"
##    col1 <- rgb(t(col2rgb("red"))/255,alpha = 0.5)
    col2 <- "white"

    par(mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:15),
                        byrow = FALSE,ncol=3,nrow=5),
                 widths=c(rep(1,3)),
                 heights=rep(1,5))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,7))
    nspec <- length(x)
    for(spec in 1:nspec){
        bb <- sapply(x[[spec]], function(x) x$bbmsyRETot[1,hcrs])
        bbsd <- sapply(x[[spec]], function(x) x$bbmsySDTot[1,hcrs])
        ff <- sapply(x[[spec]], function(x) x$ffmsyRETot[1,hcrs])
        ffsd <- sapply(x[[spec]], function(x) x$ffmsySDTot[1,hcrs])
        tac <- sapply(x[[spec]], function(x) x$tacRETot[1,hcrs])
        csd <- sapply(x[[spec]], function(x) x$cpSDTot[1,hcrs])
        nx <- 4
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimB,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = seq(0,2,0.1), col = "grey80", lty = 2)
        abline(h = median(bb[,2]), col = col1, lwd = 1.5, lty = 1)
        vioplot(bb[,1],bb[,2],bb[,3],bb[,4], col = col2,
                xaxt="n", yaxt="n", axes = FALSE, add = TRUE)
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ B/B[MSY]),2, line = 3)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec],3,line=1.5, font=2)
        box()
        ## FF
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimF,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = seq(0,2,0.25), col = "grey80", lty = 2)
        abline(h = median(ff[,2]), col = col1, lwd = 1.5, lty = 1)
        vioplot(ff[,1],ff[,2],ff[,3],ff[,4], col = col2,
                xaxt="n", yaxt="n", axes = FALSE, add = TRUE)
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Bsd
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimBsd,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = median(bbsd[,2]), col = col1, lwd = 1.5, lty = 1)
        vioplot(bbsd[,1],bbsd[,2],bbsd[,3],bbsd[,4], col = col2,
                xaxt="n", yaxt="n", axes = FALSE, add = TRUE)
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ B/B[MSY]),2, line = 3)
        }
        box()
        ## Fsd
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimFsd,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = median(ffsd[,2]), col = col1, lwd = 1.5, lty = 1)
        vioplot(ffsd[,1],ffsd[,2],ffsd[,3],ffsd[,4], col = col2,
                xaxt="n", yaxt="n", axes = FALSE, add = TRUE)
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Cp
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimCsd,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = median(csd[,2]), col = col1, lwd = 1.5, lty = 1)
        vioplot(csd[,1],csd[,2],csd[,3],csd[,4], col = col2,
                xaxt="n", yaxt="n", axes = FALSE, add = TRUE)
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ C[p]),2, line = 3)
        }
        axis(1, at = scens,labels=paste0("S",scens))
        box()
    }
    mtext("Scenarios",1,line = 3, outer = TRUE)
}



plot.rel.err.sd.box <- function(x, outline = TRUE, hcrs = c(3,4)){


    ylimB <- c(0.93,1.07) * range(unlist(lapply(x,function(y) lapply(y, function(z)
        boxplot.stats(z$bbmsyRETot[1,hcrs])$stats[c(1,5)]))))
    ylimBsd <- c(0.95,1.07) * range(unlist(lapply(x,function(y) lapply(y, function(z)
        boxplot.stats(z$bbmsySDTot[1,hcrs])$stats[c(1,5)]))))
    ylimF <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z)
        boxplot.stats(z$ffmsyRETot[1,hcrs])$stats[c(1,5)]))))
    ylimFsd <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z)
        boxplot.stats(z$ffmsySDTot[1,hcrs])$stats[c(1,5)]))))
    ylimT <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z)
        boxplot.stats(z$tacRETot[1,hcrs])$stats[c(1,5)]))))
    ylimCsd <- c(0.95,1.05) * range(unlist(lapply(x,function(y) lapply(y, function(z)
        boxplot.stats(z$cpSDTot[1,hcrs])$stats[c(1,5)]))))

##    ylimB <- ylimBsd <- ylimF <- ylimFsd <- ylimT <- ylimCsd <- c(0,3.6)

    scens <- 1:4
    nhcr <- length(hcrs)

    col1 <- "dodgerblue2"
##    col1 <- rgb(t(col2rgb("red"))/255,alpha = 0.5)
    col2 <- "white"

    par(mfrow = c(2,3),mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:15),
                        byrow = FALSE,ncol=3,nrow=5),
                 widths=c(rep(1,3)),
                 heights=rep(1,5))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,7))
    nspec <- length(x)
    for(spec in 1:nspec){
        bb <- sapply(x[[spec]], function(x) x$bbmsyRETot[1,hcrs])
        bbsd <- sapply(x[[spec]], function(x) x$bbmsySDTot[1,hcrs])
        ff <- sapply(x[[spec]], function(x) x$ffmsyRETot[1,hcrs])
        ffsd <- sapply(x[[spec]], function(x) x$ffmsySDTot[1,hcrs])
        tac <- sapply(x[[spec]], function(x) x$tacRETot[1,hcrs])
        csd <- sapply(x[[spec]], function(x) x$cpSDTot[1,hcrs])
        nx <- 4
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimB,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = at[seq(1,10,2)], col = "grey80", lty = 2)
        bx <- boxplot(bb, xaxt="n", yaxt="n", axes = FALSE, plot = FALSE,
                      outline = outline, add = TRUE, col = col2)
        ##      lines(scens, bx$stats[3,], lty=2)
        abline(h = bx$stats[3,2], col = col1, lwd = 1.5, lty = 1)
        bx <- boxplot(bb, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, add = TRUE, col = col2)
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ B/B[MSY]),2, line = 3)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec],3,line=1.5, font=2)
        box()
        ## FF
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimF,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = at[seq(1,10,2)], col = "grey80", lty = 2)
        bx <- boxplot(ff, xaxt="n", yaxt="n", axes = FALSE,plot=FALSE,
                      outline = outline, add = TRUE, col = col2)
        ##      lines(scens, bx$stats[3,], lty=2)
        abline(h = bx$stats[3,2], col = col1, lwd = 2, lty = 1)
        bx <- boxplot(ff, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, add = TRUE, col = col2)
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Bsd
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimBsd,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = at[seq(1,10,2)], col = "grey80", lty = 2)
        bx <- boxplot(bbsd, xaxt="n", yaxt="n", axes = FALSE, plot =FALSE,
                      outline = outline, add = TRUE, col = col2)
        ##      lines(scens, bx$stats[3,], lty=2)
        abline(h = bx$stats[3,2], col = col1, lwd = 2, lty = 1)
        bx <- boxplot(bbsd, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, add = TRUE, col = col2)
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ B/B[MSY]),2, line = 3)
        }
        box()
        ## Fsd
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimFsd,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = at[seq(1,10,2)], col = "grey80", lty = 2)
        bx <- boxplot(ffsd, xaxt="n", yaxt="n", axes = FALSE,plot=FALSE,
                      outline = outline, add = TRUE, col = col2)
        ##      lines(scens, bx$stats[3,], lty=2)
        abline(h = bx$stats[3,2], col = col1, lwd = 2, lty = 1)
        bx <- boxplot(ffsd, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, add = TRUE, col = col2)
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Cp
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimCsd,
             xlim = c(-0.5,0.5) + range(scens))
        at <- axTicks(2)
        abline(h = at[seq(1,10,2)], col = "grey80", lty = 2)
        bx <- boxplot(csd, xaxt="n", yaxt="n", axes = FALSE,plot=FALSE,
                      outline = outline, add = TRUE, col = col2)
        ##      lines(scens, bx$stats[3,], lty=2)
        abline(h = bx$stats[3,2], col = col1, lwd = 2, lty = 1)
        bx <- boxplot(csd, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, add = TRUE, col = col2)
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ C[p]),2, line = 3)
        }
        axis(1, at = scens,labels=paste0("S",scens))
        box()
    }
    mtext("Scenarios",1,line = 3, outer = TRUE)
}


plot.hs.vs.msy <- function(x, blackwhite = FALSE, hcrs1 = c(3,11:19), hcrs2 = c(47,55:63)){

    nhcrs1 <- length(hcrs1)
    nhcrs2 <- length(hcrs2)

    if(nhcrs1 != nhcrs2) stop("Two HCR vector do not have the same length!")

    ## plot character color
    if(blackwhite){
        col1 <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        col2 <- paste0(rep("grey",length(indC)),round(seq(10,75,length.out = length(indC))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        col2 <- colorRampPalette(c("darkorange4","darkorange1"))(nhcrs2)
        col1 <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(nhcrs1)
    }
    cols <- c(col1,col2)
    ## plot character size
    cex <- 1.7
    legendwidth = 0.4

    hcrsUse <- c(hcrs1,hcrs2)

    ## lims
    limR <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1,hcrsUse]))))
    limY <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldTot[1,hcrsUse]))))
    limA <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldDiff[1,hcrsUse]))))
    limB <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$bbmsyTot[1,hcrsUse]))))
    limF <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$ffmsyTot[1,hcrsUse]))))

    ly <- layout(matrix(c(1:15,rep(16,5)),
                        byrow = FALSE,ncol=4,nrow=5),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,5))
##    opar <- par(mar=c(2,0,1,0),oma=c(6,8,4.5,2))
    opar <- par(mar=c(4,3,2,1),oma=c(5,6,2.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        nscen <- length(x[[spec]])
        ## Risk
        limR <- c(0.95,1.05) * range(unlist(lapply(x[[spec]], function(z) z$pblimTot[1,hcrsUse])))
        plot(1,1,ty='n', axes = FALSE,xlim = limR, ylim = limR,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$pblimTot[1:3,hcrs1]
            yplot <- x[[spec]][[scen]]$pblimTot[1:3,hcrs2]
            lines(xplot[1,], yplot[1,], lty = 3,
                  col = tail(col1,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c(col1))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("Risk",1,3)
        if(spec == 1){
            mtext("Risk",2,3)
        }
        box()
        mtext(c("Anchovy","Haddock","Ling")[spec],3,1.5,font=2)
        ## Yield
        limY <- c(0.95,1.05) * range(unlist(lapply(x[[spec]], function(z) z$yieldTot[1,hcrsUse])))
        plot(1,1,ty='n', axes = FALSE,xlim = limY, ylim = limY,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$yieldTot[1:3,hcrs1]
            yplot <- x[[spec]][[scen]]$yieldTot[1:3,hcrs2]
            lines(xplot[1,], yplot[1,], lty = 3,
                  col = tail(col1,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c(col1))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("Rel. yield",1,3)
        if(spec == 1){
            mtext("Rel. yield",2,3)
        }
        box()
        ## AAV
        limA <- c(0.95,1.05) * range(unlist(lapply(x[[spec]], function(z) z$yieldDiff[1,hcrsUse])))
        plot(1,1,ty='n', axes = FALSE,xlim = limA, ylim = limA,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$yieldDiff[1:3,hcrs1]
            yplot <- x[[spec]][[scen]]$yieldDiff[1:3,hcrs2]
            lines(xplot[1,], yplot[1,], lty = 3,
                  col = tail(col1,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c(col1))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("AAV",1,3)
        if(spec == 1){
            mtext("AAV",2,3)
        }
        box()
        ## BBmsy
        limB <- c(0.95,1.05) * range(unlist(lapply(x[[spec]], function(z) z$bbmsyTot[1,hcrsUse])))
        plot(1,1,ty='n', axes = FALSE,xlim = limB, ylim = limB,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$bbmsyTot[1:3,hcrs1]
            yplot <- x[[spec]][[scen]]$bbmsyTot[1:3,hcrs2]
            lines(xplot[1,], yplot[1,], lty = 3,
                  col = tail(col1,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c(col1))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("BBmsy",1,3)
        if(spec == 1){
            mtext("BBmsy",2,3)
        }
        box()
        ## FFmsy
        limF <- c(0.95,1.05) * range(unlist(lapply(x[[spec]], function(z) z$ffmsyTot[1,hcrsUse])))
        plot(1,1,ty='n', axes = FALSE,xlim = limF, ylim = limF,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$ffmsyTot[1:3,hcrs1]
            yplot <- x[[spec]][[scen]]$ffmsyTot[1:3,hcrs2]
            lines(xplot[1,], yplot[1,], lty = 3,
                  col = tail(col1,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c(col1))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("FFmsy",1,3)
        if(spec == 1){
            mtext("FFmsy",2,3)
        }
        box()
    }
    mtext("HS", 2, 3.5, outer = TRUE, font = 2, cex = 1.2)
    mtext("MSY", 1, 3.5, outer = TRUE, font = 2, adj = 0.45, cex = 1.2)
    ## legend
    plot.new()
    legend(0.2,0.9, legend=c("Med","C45","C40","C35","C30","C25","C15","C05","C01","C001"),
           title.adj=0.1,
           pch=16,
           title = "HCR",
           col=c(col1),
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)
    legend(0.2,0.3, legend=paste0("S",1:nscen),
           title.adj=0.1,
           pch=c(15,16,17,18),
           title = "Scenario",
           col="black",
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)

}








## plot differences in risk between scenarios for all HCRs
plot.risk.diff.scen <- function(x, blackwhite = FALSE){

    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    ## plot character color
    if(blackwhite){
        colA <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        colC <- paste0(rep("grey",length(indC[-1])),round(seq(10,75,length.out = length(indC[-1]))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        colA <- c("grey60",colorRampPalette(c("darkorange4","darkorange1"))(length(indA[-1])))
        colC <- c("grey60",colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(indC[-1])))
    }
    ## plot character size
    cex <- 1.7
    legendwidth <- 0.3

    hcrs <- 3:36
    hcrsList <- list(hcrsHSa = c(20:27), ## HS - A
                    hcrsMSYa = c(3:10), ## MSY - A
                    hcrsHSc = c(20,28:36), ## HS - C
                    hcrsMSYc = c(3,11:19)) ## MSY - C
    ##
    pstarsA <- c(0.5,0.45,0.4,0.35,0.3,0.25,0.15,0.05)
    pstarsC <- c(0.5,0.45,0.4,0.35,0.3,0.25,0.15,0.05,0.01,0.001)


##    lim <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1,hcrs]))))

    ##
    par(mfrow = c(2,3),mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:9,rep(10,3)),
                        byrow = FALSE,ncol=4,nrow=3),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,3))
    opar <- par(mar=c(4,3,2,1),oma=c(2,3,2.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        ## Diff between S2 and S1
        lim <- range(unlist(lapply(x[[spec]][1:2], function(z) z$pblimTot[1,hcrs])))
        plot(1, 1, ty= 'n', xlim = lim, ylim = lim,
             xlab="",xaxt="n", yaxt = "n",ylab='')
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, 2, font = 2)
        axis(1)
        axis(2)
        if(spec == 2) mtext("Risk S1",1,3)
        if(spec == 1){
            mtext(paste0("Risk S2"),2,4)
        }
        abline(0,1,col="grey70",lty=1)
        for(i in 1:4){
            pblim1 <- x[[spec]][[1]]$pblimTot[1:3,hcrsList[[i]]]
            pblim2 <- x[[spec]][[2]]$pblimTot[1:3,hcrsList[[i]]]
            if(i %in% c(1,2)){
                cols = colA
                pstars <- pstarsA
            }else{
                cols = colC
                pstars <- pstarsC
            }
            if(i %in% c(1,3)){
                pch = 16
            }else{
                pch = 17
            }
            lines(pblim1[1,],pblim2[1,], lwd=1.5, col = tail(cols,1), lty=3)
            points(pblim1[1,], pblim2[1,], cex=cex, pch=pch, col = cols)
        }
        box()
        ## Diff between S3 and S2
        lim <- range(unlist(lapply(x[[spec]][2:3], function(z) z$pblimTot[1,hcrs])))
        plot(1, 1, ty= 'n', xlim = lim, ylim = lim,
             xlab="",xaxt="n", yaxt = "n",ylab='')
        axis(1)
        axis(2)
        if(spec == 2) mtext("Risk S2",1,3)
        if(spec == 1){
            mtext(paste0("Risk S3"),2,4)
        }
        abline(0,1,col="grey70",lty=1)
        for(i in 1:4){
            pblim1 <- x[[spec]][[2]]$pblimTot[1:3,hcrsList[[i]]]
            pblim2 <- x[[spec]][[3]]$pblimTot[1:3,hcrsList[[i]]]
            if(i %in% c(1,2)){
                cols = colA
                pstars <- pstarsA
            }else{
                cols = colC
                pstars <- pstarsC
            }
            if(i %in% c(1,3)){
                pch = 16
            }else{
                pch = 17
            }
            lines(pblim1[1,],pblim2[1,], lwd=1.5, col = tail(cols,1), lty=3)
            points(pblim1[1,], pblim2[1,], cex=cex, pch=pch, col = cols)
        }
        box()
        ## Diff between S4 and S2
        lim <- range(unlist(lapply(x[[spec]][c(2,4)], function(z) z$pblimTot[1,hcrs])))
        plot(1, 1, ty= 'n', xlim = lim, ylim = lim,
             xlab="",xaxt="n", yaxt = "n",ylab='')
        axis(1)
        axis(2)
        if(spec == 2) mtext("Risk S2",1,3)
        if(spec == 1){
            mtext(paste0("Risk S4"),2,4)
        }
        abline(0,1,col="grey70",lty=1)
        for(i in 1:4){
            pblim1 <- x[[spec]][[2]]$pblimTot[1:3,hcrsList[[i]]]
            pblim2 <- x[[spec]][[4]]$pblimTot[1:3,hcrsList[[i]]]
            if(i %in% c(1,2)){
                cols = colA
                pstars <- pstarsA
            }else{
                cols = colC
                pstars <- pstarsC
            }
            if(i %in% c(1,3)){
                pch = 16
            }else{
                pch = 17
            }
            lines(pblim1[1,],pblim2[1,], lwd=1.5, col = tail(cols,1), lty=3)
            points(pblim1[1,], pblim2[1,], cex=cex, pch=pch, col = cols)
        }
        box()
    }
    ## legend
    plot.new()
    legend(0.2,0.6, legend=c("HS-A","HS-C","MSY-A","MSY-C"),
           title.adj=0.1,
           pch=c(16,16,17,17),
           title = "HCR type",
           col=c(tail(colA,1),tail(colC,1),tail(colA,1),tail(colC,1)),
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)

}


## plot differences in yield between scenarios for all HCRs
plot.yield.diff.scen <- function(x, blackwhite = FALSE){

    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")

    ## plot character color
    if(blackwhite){
        colA <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        colC <- paste0(rep("grey",length(indC[-1])),round(seq(10,75,length.out = length(indC[-1]))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        colA <- c("grey60",colorRampPalette(c("darkorange4","darkorange1"))(length(indA[-1])))
        colC <- c("grey60",colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(indC[-1])))
    }
    ## plot character size
    cex <- 1.7
    legendwidth <- 0.3

    hcrs <- 3:36
    hcrsList <- list(hcrsHSa = c(20:27), ## HS - A
                    hcrsMSYa = c(3:10), ## MSY - A
                    hcrsHSc = c(20,28:36), ## HS - C
                    hcrsMSYc = c(3,11:19)) ## MSY - C
    ##
    pstarsA <- c(0.5,0.45,0.4,0.35,0.3,0.25,0.15,0.05)
    pstarsC <- c(0.5,0.45,0.4,0.35,0.3,0.25,0.15,0.05,0.01,0.001)


##    lim <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1,hcrs]))))

    ##
    par(mfrow = c(2,3),mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:9,rep(10,3)),
                        byrow = FALSE,ncol=4,nrow=3),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,3))
    opar <- par(mar=c(4,3,2,1),oma=c(2,3,2.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        ## Diff between S2 and S1
        lim <- range(unlist(lapply(x[[spec]][1:2], function(z) z$yieldTot[1,hcrs])))
        plot(1, 1, ty= 'n', xlim = lim, ylim = lim,
             xlab="",xaxt="n", yaxt = "n",ylab='')
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, 2, font = 2)
        axis(1)
        axis(2)
        if(spec == 2) mtext("Yield S1",1,3)
        if(spec == 1){
            mtext(paste0("Yield S2"),2,4)
        }
        abline(0,1,col="grey70",lty=1)
        for(i in 1:4){
            pblim1 <- x[[spec]][[1]]$yieldTot[1:3,hcrsList[[i]]]
            pblim2 <- x[[spec]][[2]]$yieldTot[1:3,hcrsList[[i]]]
            if(i %in% c(1,2)){
                cols = colA
                pstars <- pstarsA
            }else{
                cols = colC
                pstars <- pstarsC
            }
            if(i %in% c(1,3)){
                pch = 16
            }else{
                pch = 17
            }
            lines(pblim1[1,],pblim2[1,], lwd=1.5, col = tail(cols,1), lty=3)
            points(pblim1[1,], pblim2[1,], cex=cex, pch=pch, col = cols)
        }
        box()
        ## Diff between S3 and S2
        lim <- range(unlist(lapply(x[[spec]][2:3], function(z) z$yieldTot[1,hcrs])))
        plot(1, 1, ty= 'n', xlim = lim, ylim = lim,
             xlab="",xaxt="n", yaxt = "n",ylab='')
        axis(1)
        axis(2)
        if(spec == 2) mtext("Yield S2",1,3)
        if(spec == 1){
            mtext(paste0("Yield S3"),2,4)
        }
        abline(0,1,col="grey70",lty=1)
        for(i in 1:4){
            pblim1 <- x[[spec]][[2]]$yieldTot[1:3,hcrsList[[i]]]
            pblim2 <- x[[spec]][[3]]$yieldTot[1:3,hcrsList[[i]]]
            if(i %in% c(1,2)){
                cols = colA
                pstars <- pstarsA
            }else{
                cols = colC
                pstars <- pstarsC
            }
            if(i %in% c(1,3)){
                pch = 16
            }else{
                pch = 17
            }
            lines(pblim1[1,],pblim2[1,], lwd=1.5, col = tail(cols,1), lty=3)
            points(pblim1[1,], pblim2[1,], cex=cex, pch=pch, col = cols)
        }
        box()
        ## Diff between S4 and S2
        lim <- range(unlist(lapply(x[[spec]][c(2,4)], function(z) z$yieldTot[1,hcrs])))
        plot(1, 1, ty= 'n', xlim = lim, ylim = lim,
             xlab="",xaxt="n", yaxt = "n",ylab='')
        axis(1)
        axis(2)
        if(spec == 2) mtext("Yield S2",1,3)
        if(spec == 1){
            mtext(paste0("Yield S4"),2,4)
        }
        abline(0,1,col="grey70",lty=1)
        for(i in 1:4){
            pblim1 <- x[[spec]][[2]]$yieldTot[1:3,hcrsList[[i]]]
            pblim2 <- x[[spec]][[4]]$yieldTot[1:3,hcrsList[[i]]]
            if(i %in% c(1,2)){
                cols = colA
                pstars <- pstarsA
            }else{
                cols = colC
                pstars <- pstarsC
            }
            if(i %in% c(1,3)){
                pch = 16
            }else{
                pch = 17
            }
            lines(pblim1[1,],pblim2[1,], lwd=1.5, col = tail(cols,1), lty=3)
            points(pblim1[1,], pblim2[1,], cex=cex, pch=pch, col = cols)
        }
        box()
    }
    ## legend
    plot.new()
    legend(0.2,0.6, legend=c("HS-A","HS-C","MSY-A","MSY-C"),
           title.adj=0.1,
           pch=c(16,16,17,17),
           title = "HCR type",
           col=c(tail(colA,1),tail(colC,1),tail(colA,1),tail(colC,1)),
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)

}



## risk and yield relative to reference HCR
plot.delta.risk.yield <- function(x, showCI = FALSE, blackwhite = FALSE, outline = TRUE,
                                  hcrType = "MSY", ylim = c(-1,1)){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    ## plot character color
    if(blackwhite){
        colA <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        colC <- paste0(rep("grey",length(indC[-1])),round(seq(10,75,length.out = length(indC[-1]))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        colA <- colorRampPalette(c("darkorange4","darkorange1"))(length(indA[-1]))
        colC <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(indC[-1]))
    }
##    cols <- c("grey30","white",colA,rep("white",length(colA)),colC,rep("white",length(colC)))
    cols <- c("grey30",colA,colC)
    ## plot character size
    cex <- 1.7
    legendwith = 0.3

    if(hcrType == "MSY"){
        hcrsSel <- c(3:19)
        hcrRef <- 1
    }else{
        hcrsSel <- c(20:36)
        hcrRef <- 2
    }
    hcrs <- c(hcrRef,hcrsSel)

    ly <- layout(matrix(c(1:3,rep(4,3)),
                        byrow = FALSE,ncol=2,nrow=3),
                 widths=c(1,legendwith),
                 heights=rep(1,3))
    opar <- par(mar=c(2,2,2,2),oma=c(6.5,4,2,0.5))
    nscen <- length(x[[1]])
    for(spec in 1:nspec){
        ## risk
        pblim <- do.call(rbind,lapply(x[[spec]], function(y) y$pblimTot[1,hcrs]))
        diffRisk <- t(apply(pblim, 1, function(y) y[1] - y))[,-1]
        ## yield
        yield <- do.call(rbind,lapply(x[[spec]], function(y) y$yieldTot[1,hcrs]))
        diffYield <- t(apply(yield, 1, function(y) y - y[1]))[,-1]
        ## AAV
        aav <- do.call(rbind,lapply(x[[spec]], function(y) y$yieldDiff[1,hcrs]))
        diffAAV <- t(apply(aav, 1, function(y) y[1] - y))[,-1]
        ## combine
        comb <- do.call(cbind,lapply(as.list(1:ncol(diffYield)), function(y)
            cbind(diffRisk[,y],diffYield[,y],diffAAV[,y])))
        cols2 <- rep("white",ncol(comb))
        cols2[seq(1,ncol(comb),3)] <- cols
        bods2 <- rep("black",ncol(comb))
        bods2[seq(3,ncol(comb),3)] <- cols
        bx0 <- boxplot(comb, col = cols2, xaxt="n", outline = outline, plot = FALSE)
        plot(1,1, ty="n", axes = FALSE,
             xlab = "", ylab = "",
            xlim = c(0,ncol(comb)), ylim = ylim) ## range(bx0$stats))
        abline(h = 0, col = "grey30",lwd=1.5)
        lines(seq(1,ncol(comb),3), bx0$stats[3, seq(1,ncol(comb),3)], lty=3)
        bx <- boxplot(comb, col = cols2, xaxt="n", border = bods2,
                      outline = outline, add = TRUE)
        rect(seq(2,ncol(comb),3) - 0.4,
             bx$stats[2, seq(2,ncol(comb),3)],
             seq(2,ncol(comb),3) + 0.4,
             bx$stats[4,seq(2,ncol(comb),3)],
             density=30, angle=45, col = cols)
        if(spec == 3) axis(1,at=seq(1,ncol(comb),2) + 0.5,labels = hcrs[-1], las=2,cex =1.2)
        box()
        mtext(c("Anchovy","Haddock","Ling")[spec],font=2,line=0.5)
    }
    ## Axes labels
    mtext("Risk and yield relative to Ref", 2, outer=TRUE, line=1.5)
    mtext("HCRs", 1, outer=TRUE, line=4.5)
    ## legend
    plot.new()
    legend(0.2,0.6, legend = c("Risk","Yield","AAV"),
           title.adj=0.1,
           fill = c("grey20","grey20","white"),
           border = c("grey20","grey20","grey20"),
           title = "Metric",
           density = c(NA,20,NA),
           y.intersp=1.3, pt.cex = 4,
           ncol = 1, bty="n", cex=1.5)
}



## quants over time
plot.time <- function(x, hcrs = c(3,20,15),hcrref = 1,scenario = 1,
                      cols = c("dodgerblue2","darkorange","darkgreen","darkred")
                      ){

    ylimBB <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$bbmsyYear[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimFF <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$ffmsyYear[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimRisk <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$pblim[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimYield <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$yield[1,hcrs,],2,median,na.rm=TRUE)))))

    legendwidth <- 0.3
    years <- 1:20
    col1 <- rgb(t(col2rgb("grey60"))/255,alpha=0.6)

    ylimB <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        y[[scenario]]$bbmsyYear[1,c(hcrref,hcrs),])))
    ylimF <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        y[[scenario]]$ffmsyYear[1,c(hcrref,hcrs),])))
    ylimR <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        y[[scenario]]$pblim[1,c(hcrref,hcrs),])))
    ylimY <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        y[[scenario]]$yield[1,c(hcrref,hcrs),])))

    ly <- layout(matrix(c(1:12,rep(13,4)),
                        byrow = FALSE,ncol=4,nrow=4),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,2))
    nspec <- length(x)
    for(spec in 1:nspec){
        risk <- x[[spec]][[scenario]]$pblim[1,hcrs,]
        riskref <- x[[spec]][[scenario]]$pblim[1:3,hcrref,]
        yield <- x[[spec]][[scenario]]$yield[1,hcrs,]
        yieldref <- x[[spec]][[scenario]]$yield[1:3,hcrref,]
        bbmsy <- x[[spec]][[scenario]]$bbmsyYear[1,hcrs,]
        bbmsyref <- x[[spec]][[scenario]]$bbmsyYear[1:3,hcrref,]
        ffmsy <- x[[spec]][[scenario]]$ffmsyYear[1,hcrs,]
        ffmsyref <- x[[spec]][[scenario]]$ffmsyYear[1:3,hcrref,]
        ## Risk
        par(mar=c(0,0,0,0))
        if(spec == 3) par(mar=c(0,0,0,1))
        plot(years, risk[1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimR)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        ## ref
        polygon(c(years,rev(years)), c(riskref[2,],rev(riskref[3,])),
                border=NA, col=col1)
        lines(years, riskref[1,], lwd=1.5, col=col1)
        ## HCRs
        for(i in 1:nrow(risk)){
            lines(years, risk[i,], lwd=1.5, col = cols[i])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Risk"),2, line = 3)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, 0.7, font=2,cex=1)
        box()
        ## Yield
        plot(years, yield[1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimY)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = 1, col = "grey60", lty = 2)
        ## ref
        polygon(c(years,rev(years)), c(yieldref[2,],rev(yieldref[3,])),
                border=NA, col=col1)
        lines(years, yieldref[1,], lwd=1.5, col=col1)
        ## HCRs
        for(i in 1:nrow(yield)){
            lines(years, yield[i,], lwd=1.5, col = cols[i])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. yield"),2, line = 3)
        }
        box()
        ## BBmsy
        plot(years, bbmsy[1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimB)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = 1, col = "grey60", lty = 2)
        ## ref
        polygon(c(years,rev(years)), c(bbmsyref[2,],rev(bbmsyref[3,])),
                border=NA, col=col1)
        lines(years, bbmsyref[1,], lwd=1.5, col=col1)
        ## HCRs
        for(i in 1:nrow(bbmsy)){
            lines(years, bbmsy[i,], lwd=1.5, col = cols[i])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("BBmsy"),2, line = 3)
        }
        box()
        ## FFmsy
        plot(years, ffmsy[1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimF)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = 1, col = "grey60", lty = 2)
        ## ref
        polygon(c(years,rev(years)), c(ffmsyref[2,],rev(ffmsyref[3,])),
                border=NA, col=col1)
        lines(years, ffmsyref[1,], lwd=1.5, col=col1)
        ## HCRs
        for(i in 1:nrow(ffmsy)){
            lines(years, ffmsy[i,], lwd=1.5, col = cols[i])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("FFmsy"),2, line = 3)
        }
        axis(1)
        box()
    }
    mtext("Projection period [yr]", 1, 3, outer=TRUE)
    par(mar=c(5,0,4,0.1))
    plot(c(0,1),type="n", axes=F, xlab="", ylab="")
    legend("center", legend = hcrsAll[c(hcrref,hcrs)],
           col=c(col1,cols),
           title= "HCRs",
           bty='n',
           cex=1.2, lwd=1.5)
}



## quants over time
plot.time.focus <- function(x, hcrs = c(40,14),hcrref = 1,scenario = 1,
                      cols = c("darkorange4","dodgerblue4"), traj = 5
                      ){

    ylimBB <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$bbmsyYear[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimFF <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$ffmsyYear[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimRisk <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$pblim[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimYield <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$yield[1,hcrs,],2,median,na.rm=TRUE)))))

    legendwidth <- 0.3
    years <- 1:20
    col1 <- rgb(t(col2rgb("grey60"))/255,alpha=0.6)
    col1a <- "grey60"
    colsP <- sapply(cols, function(x) rgb(t(col2rgb(x))/255,alpha=0.4))

    ylimB <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        list(y[[scenario]]$bbmsyYear[1,c(hcrref,hcrs),],
             y[[scenario]]$bbmsy[traj,c(hcrref,hcrs),]))))
    ylimF <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        list(y[[scenario]]$ffmsyYear[1,c(hcrref,hcrs),],
             y[[scenario]]$ffmsy[traj,c(hcrref,hcrs),]))))
    ylimR <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        y[[scenario]]$pblim[1,c(hcrref,hcrs),])))
    ylimY <- c(0.95,1.05) * range(unlist(lapply(x,function(y)
        list(y[[scenario]]$yield[1,c(hcrref,hcrs),],
             y[[scenario]]$catch[traj,c(hcrref,hcrs),]/y[[scenario]]$refcatch[traj]))))


    ## ly <- layout(matrix(c(1:12,rep(13,4)),
    ##                     byrow = FALSE,ncol=4,nrow=4),
    ##              widths=c(rep(1,3),legendwidth),
    ##              heights=rep(1,4))
    ly <- layout(matrix(c(1:12),
                        byrow = FALSE,ncol=3,nrow=4),
                 widths=c(rep(1,3)),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(5.5,6,4.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        risk <- x[[spec]][[scenario]]$pblim[1:3,hcrs,]
        riskref <- x[[spec]][[scenario]]$pblim[1:3,hcrref,]
        yield <- x[[spec]][[scenario]]$yield[1:3,hcrs,]
        yieldref <- x[[spec]][[scenario]]$yield[1:3,hcrref,]
        yieldtraj <- x[[spec]][[scenario]]$catch[traj, hcrs, ] / x[[spec]][[scenario]]$refcatch[traj]
        yieldreftraj <- x[[spec]][[scenario]]$catch[traj, hcrref, ] / x[[spec]][[scenario]]$refcatch[traj]
        bbmsy <- x[[spec]][[scenario]]$bbmsyYear[1:3,hcrs,]
        bbmsyref <- x[[spec]][[scenario]]$bbmsyYear[1:3,hcrref,]
        bbmsytraj <- x[[spec]][[scenario]]$bbmsy[traj, hcrs, ]
        bbmsyreftraj <- x[[spec]][[scenario]]$bbmsy[traj,hcrref,]
        ffmsy <- x[[spec]][[scenario]]$ffmsyYear[1:3,hcrs,]
        ffmsyref <- x[[spec]][[scenario]]$ffmsyYear[1:3,hcrref,]
        ffmsyreftraj <- x[[spec]][[scenario]]$ffmsy[traj,hcrref,]
        ffmsytraj <- x[[spec]][[scenario]]$ffmsy[traj, hcrs, ]
        ## Risk
        par(mar=c(0,0,0,0))
        if(spec == 3) par(mar=c(0,0,0,1))
        plot(years, risk[1,1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimR)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        ## ref
        polygon(c(years,rev(years)), c(riskref[2,],rev(riskref[3,])),
                border=NA, col=col1)
        lines(years, riskref[1,], lwd=2, col=col1)
        ## HCRs
        for(i in 1:nrow(risk[1,,])){
            polygon(c(years,rev(years)), c(risk[2,i,],rev(risk[3,i,])),
                    border = NA, col = colsP[i])
            lines(years, risk[1,i,], lwd=2, col = cols[i])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Risk"),2, line = 3)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, 0.7, font=2,cex=1)
        if(spec==3) legend("topright", legend = hcrsAll[c(hcrref,hcrs)],
           col=c(col1,cols),
           bty='n',
           cex=1.1, lwd=1.5)
        box()
        ## Yield
        plot(years, yield[1,1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimY)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = 1, col = "grey60", lty = 2, lwd=2)
        ## ref
        polygon(c(years,rev(years)), c(yieldref[2,],rev(yieldref[3,])),
                border=NA, col=col1)
        lines(years, yieldref[1,], lwd=2, col=col1a)
        lines(years, yieldreftraj, lwd=2, col = col1a, lty=3)
        ## HCRs
        for(i in 1:nrow(yield[1,,])){
            polygon(c(years,rev(years)), c(yield[2,i,],rev(yield[3,i,])),
                    border = NA, col = colsP[i])
            lines(years, yield[1,i,], lwd=2, col = cols[i])
            lines(years, yieldtraj[i,], lwd=2, col = cols[i], lty=3)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. yield"),2, line = 3)
        }
        box()
        ## BBmsy
        plot(years, bbmsy[1,1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimB)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = 1, col = "grey60", lty = 2, lwd=2)
        ## ref
        polygon(c(years,rev(years)), c(bbmsyref[2,],rev(bbmsyref[3,])),
                border=NA, col=col1)
        lines(years, bbmsyref[1,], lwd=2, col=col1a)
        lines(years, bbmsyreftraj, lwd=2, col = col1a, lty=3)
        ## HCRs
        for(i in 1:nrow(bbmsy[1,,])){
            polygon(c(years,rev(years)), c(bbmsy[2,i,],rev(bbmsy[3,i,])),
                    border = NA, col = colsP[i])
            lines(years, bbmsy[1,i,], lwd=2, col = cols[i])
            lines(years, bbmsytraj[i,], lwd=2, col = cols[i], lty=3)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("BBmsy"),2, line = 3)
        }
        box()
        ## FFmsy
        plot(years, ffmsy[1,1,],
             type = "n",
             xaxt = "n", yaxt = "n",
             ylim = ylimF)
##        abline(v = seq(0,20,1), col = "grey80", lty = 2)
        ##abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
        abline(h = 1, col = "grey60", lty = 2, lwd=2)
        ## ref
        polygon(c(years,rev(years)), c(ffmsyref[2,],rev(ffmsyref[3,])),
                border=NA, col=col1)
        lines(years, ffmsyref[1,], lwd=2, col=col1a)
        lines(years, ffmsyreftraj, lwd=2, col = col1a, lty=3)
        ## HCRs
        for(i in 1:nrow(ffmsy[1,,])){
            polygon(c(years,rev(years)), c(ffmsy[2,i,],rev(ffmsy[3,i,])),
                    border = NA, col = colsP[i])
            lines(years, ffmsy[1,i,], lwd=2, col = cols[i])
            lines(years, ffmsytraj[i,], lwd=2, col = cols[i], lty=3)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("FFmsy"),2, line = 3)
        }
        axis(1)
        box()
    }
    mtext("Projection period [yr]", 1, 3, outer=TRUE)
    ## par(mar=c(5,0,4,0.1))
    ## plot(c(0,1),type="n", axes=F, xlab="", ylab="")
    ## legend("center", legend = hcrsAll[c(hcrref,hcrs)],
    ##        col=c(col1,cols),
    ##        title= "HCR",
    ##        bty='n',
    ##        cex=1, lwd=1.5)
}




## plot sufficiency of sample size
plot.suff.samp <- function(x, quant = "pblimY3_5"){
    if(quant == "pblimY3_5"){
        ylab <- "Risk [years 3-5]"
    }else if(quant == "pblimY6_20"){
        ylab <- "Risk [years 6-20]"
    }else if(quant == "yieldY3_5"){
        ylab <- "Rel. yield [years 3-5]"
    }else if(quant == "yieldY6_20"){
        ylab <- "Rel. yield [years 6-20]"
    }else if(quant == "aav"){
        ylab <- "AAV"
    }else ylab <- quant

    ly <- layout(matrix(c(1:12),
                        byrow = FALSE,ncol=3,nrow=4),
                 widths=c(rep(1,3)),
                 heights=rep(1,4))
    opar <- par(mar=c(2,2,2,2),oma=c(4,4,1,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        nscen <- length(x[[spec]])
        for(scen in 1:nscen){
            matplot(suffSamp[[spec]][[scen]][[quant]],ty='l')
            usr <- par("usr")
            ## plot title
            txt <- paste0(c("A","H","L")[spec],scen)
            mtext(txt, 3,0.1,font=2, cex=0.8, col = "black")
            box()
        }
    }
    mtext(ylab, 2, 1.5, outer = TRUE, cex = 1)
    mtext("Sample size", 1, 1.5, outer = TRUE, cex = 1)
}


## plot spict convergence by year
plot.conv.year <- function(x, cutoff = 3){
    nn <- ncol(x)
    layout(matrix(c(1,2),nrow=1), width=c(4,0.8))
    par(mar=c(5,5,4,0))
    matplot(x,type="b",pch = 1:7, col = 1:6, lty = 1:5,
            lwd = 1.5, cex = 1.2,
            ylab = "Convergence [%]", xlab="Projection period [yrs]")
    abline(v = cutoff, lty = 2, lwd= 1.5)
    par(mar=c(5,0,4,0.6))
    plot(c(0,1),type="n", axes=F, xlab="", ylab="")
    legend("center", colnames(x),col=rep(1:6,10),
           title= "Scenario",
           cex=1,pch = rep(1:7,10),, lty=rep(1:5,10))
}



## OLD:
#################################################################################

## errors and sd over time
plot.time.old <- function(x, hcrs = c(3:36)){
    ylimB <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$bbmsyRE[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimBsd <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$bbmsySD[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimF <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$ffmsyRE[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimFsd <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$ffmsySD[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimCsd <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$cpSD[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimRisk <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$pblim[1,hcrs,],2,median,na.rm=TRUE)))))
    ylimYield <- range(unlist(lapply(x,function(y) lapply(y, function(z)
        apply(z$yield[1,hcrs,],2,median,na.rm=TRUE)))))

    xplot <- 1:20
    xlim <- range(xplot)
    xplot2 <- 1:19
    legendwidth <- 0.3
    cols <- c("dodgerblue2","darkorange","darkgreen","darkred")

    par(mfrow = c(2,3),mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:21,rep(22,7)),
                        byrow = FALSE,ncol=4,nrow=7),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,7))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,2))
    nspec <- length(x)
    for(spec in 1:nspec){
        bb <- sapply(x[[spec]], function(x) apply(x$bbmsyRE[1,hcrs,],2,median,na.rm=TRUE))
        bbsd <- sapply(x[[spec]], function(x) apply(x$bbmsySD[1,hcrs,],2,median,na.rm=TRUE))
        ff <- sapply(x[[spec]], function(x) apply(x$ffmsyRE[1,hcrs,],2,median,na.rm=TRUE))
        ffsd <- sapply(x[[spec]], function(x) apply(x$ffmsySD[1,hcrs,],2,median,na.rm=TRUE))
        csd <- sapply(x[[spec]], function(x) apply(x$cpSD[1,hcrs,],2,median,na.rm=TRUE))
        risk <- sapply(x[[spec]], function(x) apply(x$pblim[1,hcrs,],2,median,na.rm=TRUE))
        yield <- sapply(x[[spec]], function(x) apply(x$yield[1,hcrs,],2,median,na.rm=TRUE))
        nscen <- ncol(risk)
        ## risk
        par(mar=c(0,0,0,0))
        if(spec == 3) par(mar=c(0,0,0,1))
        plot(xplot, risk[,1], ty= 'n',
             ylim = ylimRisk,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot, risk[,scen], col = cols[scen], lwd = 1.5)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, 0.7, font=2,cex=1)
        if(spec == 1){
            axis(2)
            mtext(bquote("Risk"),2, line = 3)
        }
        box()
        ## yield
        plot(xplot, yield[,1], ty= 'n',
             ylim = ylimYield,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot, yield[,scen], col = cols[scen], lwd = 1.5)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. yield"),2, line = 3)
        }
        box()
        ## BBmsy
        plot(xplot2, bb[,1], ty= 'n', xlim = xlim, ylim = ylimB,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot2, bb[,scen], col = cols[scen], lwd = 1.5)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ B/B[MSY]),2, line = 3)
        }
        box()
        ## FFmsy
        plot(xplot2, ff[,1], ty= 'n', xlim = xlim, ylim = ylimF,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot2, ff[,scen], col = cols[scen], lwd = 1.5)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## BB SD
        plot(xplot2, bbsd[,1], ty= 'n', xlim = xlim, ylim = ylimBsd,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot2, bbsd[,scen], col = cols[scen], lty = 1, lwd = 1.5)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ B/B[MSY]),2, line = 3)
        }
        box()
        ## FF SD
        plot(xplot2, ffsd[,1], ty= 'n', xlim = xlim, ylim = ylimFsd,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot2, ffsd[,scen], col = cols[scen], lty = 1, lwd = 1.5)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Cp
        plot(xplot2, csd[,1], ty= 'n', xlim = xlim, ylim = ylimCsd,
             xaxt= 'n', yaxt = 'n', xlab = '', ylab = '')
        for(scen in 1:nscen){
            lines(xplot2, csd[,scen], col = cols[scen], lty = 1, lwd = 1.5)
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ C[p]),2, line = 3)
        }
        axis(1)
        box()
    }
    mtext("Projection period [yr]", 1, 3, outer=TRUE)
    par(mar=c(5,0,4,0.1))
    plot(c(0,1),type="n", axes=F, xlab="", ylab="")
    legend("center", legend = paste0("S",1:nscen),
           col=cols[1:nscen],
           title= "Scenario",
           bty='n',
           cex=1.2, lwd=1.5)
}


## Plot rel error in spict est. quantities and SD of quantities
plot.rel.err.sd.single <- function(x, outline = TRUE, hcrs = c(3,4)){

    ylimB <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$bbmsyRETot[1:3,hcrs]))))
    ylimBsd <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$bbmsySDTot[1:3,hcrs]))))
    ylimF <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$ffmsyRETot[1:3,hcrs]))))
    ylimFsd <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$ffmsySDTot[1:3,hcrs]))))
    ylimT <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$tacRETot[1:3,hcrs]))))
    ylimCsd <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$cpSDTot[1:3,hcrs]))))

    ylimB <- ylimBsd <- ylimF <- ylimFsd <- ylimT <- ylimCsd <- c(0,2.3)

    scens <- 1:4
    nhcr <- length(hcrs)


    col1 <- "darkgreen"

    par(mfrow = c(2,3),mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:15),
                        byrow = FALSE,ncol=3,nrow=5),
                 widths=c(rep(1,3)),
                 heights=rep(1,5))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,7))
    nspec <- length(x)
    for(spec in 1:nspec){
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimB)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        for(hcr in 1:nhcr){
            bb <- sapply(x[[spec]], function(x) x$bbmsyRETot[1:3,hcrs[hcr]])
            points(scens, bb[1,], pch = hcr)
            lines(scens, bb[1,], lty=2)
            segments(scens, bb[2,],scens,bb[3,])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ B/B[MSY]),2, line = 3)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec],3,line=1.5, font=2)
        box()
        ## FF
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimF)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        for(hcr in 1:nhcr){
            ff <- sapply(x[[spec]], function(x) x$ffmsyRETot[1:3,hcrs[hcr]])
            points(scens, ff[1,], pch = hcr)
            lines(scens, ff[1,], lty=2)
            segments(scens, ff[2,],scens,ff[3,])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Bsd
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimBsd)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        for(hcr in 1:nhcr){
            bbsd <- sapply(x[[spec]], function(x) x$bbmsySDTot[1:3,hcrs])
            points(scens, bbsd[1,], pch = hcr)
            lines(scens, bbsd[1,], lty=2)
            segments(scens, bbsd[2,],scens,bbsd[3,])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ B/B[MSY]),2, line = 3)
        }
        box()
        ## Fsd
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimFsd)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        for(hcr in 1:nhcr){
            ffsd <- sapply(x[[spec]], function(x) x$ffmsySDTot[1:3,hcrs])
            points(scens, ffsd[1,], pch = hcr)
            lines(scens, ffsd[1,], lty=2)
            segments(scens, ffsd[2,],scens,ffsd[3,])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ F/F[MSY]),2, line = 3)
        }
        box()
        ## Cp
        plot(scens,rep(1,4), ty="n", axes = FALSE,
             xlab = "", ylab = "", ylim = ylimCsd)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        for(hcr in 1:nhcr){
            csd <- sapply(x[[spec]], function(x) x$cpSDTot[1:3,hcrs])
            points(scens, csd[1,], pch = hcr)
            lines(scens, csd[1,], lty=2)
            segments(scens, csd[2,],scens,csd[3,])
        }
        if(spec == 1){
            axis(2)
            mtext(bquote("SD" ~ C[p]),2, line = 3)
        }
        axis(1, at = scens,labels=paste0("S",scens))
        box()
    }
    mtext("Scenarios",1,line = 3, outer = TRUE)
}


## Plot rel error in spict est. quantities and SD of quantities
plot.rel.err.sd.old <- function(x, outline = TRUE, hcrs = c(3:19)){

    ylimB <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$bbmsyRETot[1:3,hcrs]))))
    ylimBsd <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$bbmsySDTot[1:3,hcrs]))))
    ylimF <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$ffmsyRETot[1:3,hcrs]))))
    ylimFsd <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$ffmsySDTot[1:3,hcrs]))))
    ylimT <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$tacRETot[1:3,hcrs]))))
    ylimCsd <- range(unlist(lapply(x,function(y) lapply(y, function(z) z$cpSDTot[1:3,hcrs]))))

    ylimB <- ylimBsd <- ylimF <- ylimFsd <- ylimT <- ylimCsd <- c(0,2.5)


    col1 <- "darkgreen"

    par(mfrow = c(2,3),mar=c(1,5,0,2),oma=c(5,0,4,0))
    ly <- layout(matrix(c(1:9),
                        byrow = FALSE,ncol=3,nrow=3),
                 widths=c(rep(1,3)),
                 heights=rep(1,2))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,7))
    nspec <- length(x)
    for(spec in 1:nspec){
        bb <- sapply(x[[spec]], function(x) x$bbmsyRETot[1,hcrs])
        bbsd <- sapply(x[[spec]], function(x) x$bbmsySDTot[1,hcrs])
        ff <- sapply(x[[spec]], function(x) x$ffmsyRETot[1,hcrs])
        ffsd <- sapply(x[[spec]], function(x) x$ffmsySDTot[1,hcrs])
        tac <- sapply(x[[spec]], function(x) x$tacRETot[1,hcrs])
        csd <- sapply(x[[spec]], function(x) x$cpSDTot[1,hcrs])
        nx <- 2*ncol(bb)
        bx <- boxplot(bbsd, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, plot = FALSE)
        plot(1,1, ty="n", axes = FALSE,
             xlab = "", ylab = "",
             xlim = c(0.5,nx+0.5), ylim = ylimB)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        boxplot(bb, xaxt="n", yaxt="n", axes = FALSE,
                outline = outline, add = TRUE, at = seq(1,nx,2))
        boxplot(bbsd, xaxt="n", yaxt="n", axes = FALSE, col = col1,
                outline = outline, add = TRUE, at = seq(2,nx,2))
        rect(seq(2,nx,2) - 0.4,
             bx$stats[2,],
             seq(2,nx,2) + 0.4,
             bx$stats[4,],
             density=30, angle=45, col = "black")
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ B/B[MSY]),2, line = 3)
        }
        if(spec == 3){
            at <- axTicks(4)
            axis(4, col.ticks = col1, col = col1, labels = FALSE)
            mtext(side = 4, text = at, at = at, col = col1, line = 1,cex=0.8)
            mtext(bquote("SD" ~ B/B[MSY]),4, line = 3, col = col1)
        }
        mtext(c("Anchovy","Haddock","Ling")[spec],3,line=1.5, font=2)
        box()
        ## FF
        bx <- boxplot(ffsd, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, plot = FALSE)
        plot(1,1, ty="n", axes = FALSE,
             xlab = "", ylab = "",
             xlim = c(0.5,nx+0.5), ylim = ylimF)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        boxplot(ff, xaxt="n", yaxt="n", axes = FALSE,
                outline = outline, add = TRUE, at = seq(1,nx,2))
        boxplot(ffsd, xaxt="n", yaxt="n", axes = FALSE,col = col1,
                outline = outline, add = TRUE, at = seq(2,nx,2))
        rect(seq(2,nx,2) - 0.4,
             bx$stats[2,],
             seq(2,nx,2) + 0.4,
             bx$stats[4,],
             density=30, angle=45, col = "black")
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ F/F[MSY]),2, line = 3)
        }
        if(spec == 3){
            at <- axTicks(4)
            axis(4, col.ticks = col1, col = col1, labels = FALSE)
            mtext(side = 4, text = at, at = at, col = col1, line = 1,cex=0.8)
            mtext(bquote("SD" ~ F/F[MSY]),4, line = 3, col = col1)
        }
        box()
        ## Csd
        bx <- boxplot(csd, xaxt="n", yaxt="n", axes = FALSE,
                      outline = outline, plot = FALSE)
        plot(1,1, ty="n", axes = FALSE,
             xlab = "", ylab = "",
             xlim = c(0.5,nx+0.5), ylim = ylimT)
        abline(h = seq(0,5,0.5), col = "grey80", lty = 2)
        boxplot(tac, xaxt="n", yaxt="n", axes = FALSE,
                outline = outline, add = TRUE, at = seq(1,nx,2))
        boxplot(csd, xaxt="n", yaxt="n", axes = FALSE,col = col1,
                outline = outline, add = TRUE, at = seq(2,nx,2))
        rect(seq(2,nx,2) - 0.4,
             bx$stats[2,],
             seq(2,nx,2) + 0.4,
             bx$stats[4,],
             density=30, angle=45, col = "black")
        box()
        if(spec == 1){
            axis(2)
            mtext(bquote("Rel. Err." ~ TAC),2, line = 3)
        }
        if(spec == 3){
            at <- axTicks(4)
            axis(4, col.ticks = col1, col = col1, labels = FALSE)
            mtext(side = 4, text = at, at = at, col = col1, line = 1,cex=0.8)
            mtext(bquote("SD" ~ C[p]),4, line = 3, col = col1)
        }
        axis(1, at = seq(1.5,nx+0.5,2),labels=paste0("S",1:4))
    }
    mtext("Scenarios",1,line = 3, outer = TRUE)
}



## plotting functions
plot.risk.yield.old <- function(x, showCI = FALSE, blackwhite = FALSE, hcrs = c(3:19)){
    ## checks
    quantAvail <- names(x[[1]][[1]])
    if(!any("pblimTot" == quantAvail)) stop("pblimTot is missing!")
    if(!any("yieldTot" == quantAvail)) stop("yieldTot is missing!")
    if(!any("sampleSize" == quantAvail)) stop("sampleSize is missing!")

    ## plot character color
    if(blackwhite){
        colA <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        colC <- paste0(rep("grey",length(indC[-1])),round(seq(10,75,length.out = length(indC[-1]))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        colA <- colorRampPalette(c("darkorange4","darkorange1"))(length(indA[-1]))
        colC <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(indC[-1]))
    }
    ## plot character size
    cex <- 1.7
    legendwidth <- 0.4

    xlim <- c(0,1.05) * range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1:3,hcrs]))))
    ylim <- c(0,1.05) * range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldTot[1:3,hcrs]))))

    ly <- layout(matrix(c(1:12,rep(13,4)),
                        byrow = FALSE,ncol=4,nrow=4),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,4))
    opar <- par(mar=c(0,0,0,0),oma=c(6,7,4.5,2))
    nspec <- length(x)
    for(spec in 1:nspec){
        nscen <- length(x[[spec]])
        xlim <- c(0,1.05) * range(unlist(lapply(x[[spec]], function(z) z$pblimTot[1:3,hcrs])))
        for(scen in 1:nscen){
            pblim <- x[[spec]][[scen]]$pblimTot[1:3,hcrs]
            yield <- x[[spec]][[scen]]$yieldTot[1:3,hcrs]
            samplesize <- x[[spec]][[scen]]$sampleSize
            plot(pblim[1,], yield[1,],
                 type = "n",
                 xaxt = "n", yaxt = "n",
                 xlim = xlim,
                 ylim = ylim)
            abline(v = seq(0,2,0.1), col = "grey80", lty = 2)
            abline(h = seq(0,2,0.2), col = "grey80", lty = 2)
            if(scen %in% c(4)) axis(1, cex.axis = 1.1)
            if(spec == 1 &&  scen %in% c(1:4)) axis(2, cex.axis = 1.1)
            ##xlab="Risk", ylab="log Yield")
            ## 95% CI
            if(!is.null(showCI) && showCI == "segments"){
                ## segments(pblim[2,1],yield[1,1],pblim[3,1],yield[1,1], col=1, lty=2)
                segments(pblim[2,indA[1]],yield[1,indA[1]],pblim[3,indA[1]],yield[1,indA[1]], col=1, lty=2)
                segments(pblim[2,indA[-1]],yield[1,indA[-1]],pblim[3,indA[-1]],yield[1,indA[-1]], col=colA, lty=2)
                segments(pblim[2,indC[-1]],yield[1,indC[-1]],pblim[3,indC[-1]],yield[1,indC[-1]], col=colC, lty=2)
                ## segments(pblim[1,1],yield[2,1],pblim[1,1],yield[3,1], col=1, lty=2)
                segments(pblim[1,indA[1]],yield[2,indA[1]],pblim[1,indA[1]],yield[3,indA[1]], col=1, lty=2)
                segments(pblim[1,indA[-1]],yield[2,indA[-1]],pblim[1,indA[-1]],yield[3,indA[-1]], col=colA, lty=2)
                segments(pblim[1,indC[-1]],yield[2,indC[-1]],pblim[1,indC[-1]],yield[3,indC[-1]], col=colC, lty=2)
            }
            if(!is.null(showCI) && showCI == "polygon"){
                ## polygon(pblim[c(2,3,3,2),1], yield[c(2,2,3,3),1],
                ##         col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
                polygon(pblim[c(2,3,3,2),1],yield[c(2,2,3,3),1],
                        col=rgb(t(col2rgb("black"))/255,alpha=0.15), border = NA)
                for(i in 2:8){
                    polygon(pblim[c(2,3,3,2),i],yield[c(2,2,3,3),i],
                            col=rgb(t(col2rgb("darkorange3"))/255,alpha=0.2), border = NA)
                }
                for(i in c(2,9:17)){
                    polygon(pblim[c(2,3,3,2),i],yield[c(2,2,3,3),i],
                            col=rgb(t(col2rgb("dodgerblue3"))/255,alpha=0.2), border = NA)
                }
            }
            ## HCR type connectors
            lines(pblim[1,1:8], yield[1,1:8], col="grey10", lty=3)
            lines(pblim[1,c(1,9:17)], yield[1,c(1,9:17)], col="grey10", lty=3)
            ## points
##            points(pblim[1,1], yield[1,1], pch = 8, col=1, cex=cex)
            points(pblim[1,1], yield[1,1], pch = 16, col=1, cex=cex)
            points(pblim[1,2:8], yield[1,2:8], pch = 15, col = colA, cex=cex)
            points(pblim[1,9:17], yield[1,9:17], pch = 17, col = colC, cex=cex)
            usr <- par("usr")
            txt <- bquote("N"~"="~.(samplesize))
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.01
            ## sample size
            xt <- 0.86 * usr[2]
            yt <- -0.7 * usr[3]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz - 0.01,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz + 0.01,
                col = "white"
            )
            text(xt, yt, txt, font=1, cex=1.1)
            ## plot title
            txt <- paste0(c("A","H","L")[spec],scen)
            sw   <- strwidth(txt)
            sh   <- strheight(txt)
            frsz <- 0.01
            xt <- -1 * usr[1]
            yt <- 0.94 * usr[4]
            rect(
                xt - sw/2 - frsz,
                yt - sh/2 - frsz - 0.01,
                xt + sw/2 + frsz,
                yt + sh/2 + frsz + 0.01,
                col = "grey40"
            )
            text(xt, yt, txt, font=2, cex=1.2, col = "white")
            box()
        }
        mtext(c("Anchovy","Haddock","Ling")[spec], 3, outer=TRUE, line=1.2, font=2,
              adj = c(0.136,0.466,0.788)[spec])
    }
    ## Axes labels
    mtext("Rel. yield", 2, outer=TRUE, line=4)
    mtext("Risk", 1, outer=TRUE, line=3.5, adj = 0.48)
    ## legend
    plot.new()
    legend("center", legend=hcrsAll[hcrs[-1]], title.adj=0.1,
           pch=c(16,rep(15,length(indA[-1])),rep(17,length(indC[-1]))),
           title = "HCR",
           col=c(1,colA, colC),
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)

}



plot.hs.vs.msy.old <- function(x, blackwhite = FALSE){

    ## plot character color
    if(blackwhite){
        colA <- paste0(rep("grey",length(indA[-1])),round(seq(10,75,length.out = length(indA[-1]))))
        colC <- paste0(rep("grey",length(indC)),round(seq(10,75,length.out = length(indC))))
    }else{
        ##        colA <- colorRampPalette(c("gold4","gold1"))(length(indA[-1]))
        colours()
        colA <- colorRampPalette(c("darkorange4","darkorange1"))(length(indA[-1]))
        colC <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(indC))
    }
    cols <- c("grey60","grey30",colA,colC)
    ## plot character size
    cex <- 1.7
    legendwidth = 0.4
    hcrsMSY = 3:19
    hcrsHS = 20:36

    ## lims
    limR <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$pblimTot[1,c(hcrsMSY,hcrsHS)]))))
    limY <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldTot[1,c(hcrsMSY,hcrsHS)]))))
    limA <- range(unlist(lapply(x, function(y) lapply(y, function(z) z$yieldDiff[1,c(hcrsMSY,hcrsHS)]))))

    ly <- layout(matrix(c(1:9,rep(10,3)),
                        byrow = FALSE,ncol=4,nrow=3),
                 widths=c(rep(1,3),legendwidth),
                 heights=rep(1,3))
##    opar <- par(mar=c(2,0,1,0),oma=c(6,8,4.5,2))
    opar <- par(mar=c(4,3,2,1),oma=c(5,6,2.5,1))
    nspec <- length(x)
    for(spec in 1:nspec){
        nscen <- length(x[[spec]])
        ## Risk
        limR <- range(unlist(lapply(x[[spec]], function(z) z$pblimTot[1,c(hcrsMSY,hcrsHS)])))
        plot(1,1,ty='n', axes = FALSE,xlim = limR, ylim = limR,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$pblimTot[1:3,hcrsMSY]
            yplot <- x[[spec]][[scen]]$pblimTot[1:3,hcrsHS]
            lines(xplot[1,c(1:8)], yplot[1,c(1:8)], lty = 3,
                  col = tail(colA,1))
            lines(xplot[1,c(1,9:17)], yplot[1,c(1,9:17)], lty = 3,
                  col = tail(colC,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c("grey60",colA,colC))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("Risk",1,3)
        if(spec == 1){
            mtext("Risk",2,3)
        }
        box()
        mtext(c("Anchovy","Haddock","Ling")[spec],3,1.5,font=2)
        ## Yield
        limY <- range(unlist(lapply(x[[spec]], function(z) z$yieldTot[1,c(hcrsMSY,hcrsHS)])))
        plot(1,1,ty='n', axes = FALSE,xlim = limY, ylim = limY,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$yieldTot[1:3,hcrsMSY]
            yplot <- x[[spec]][[scen]]$yieldTot[1:3,hcrsHS]
            lines(xplot[1,c(1:8)], yplot[1,c(1:8)], lty = 3,
                  col = tail(colA,1))
            lines(xplot[1,c(1,9:17)], yplot[1,c(1,9:17)], lty = 3,
                  col = tail(colC,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c("grey60",colA,colC))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("Rel. yield",1,3)
        if(spec == 1){
            mtext("Rel. yield",2,3)
        }
        box()
        ## AAV
        limA <- range(unlist(lapply(x[[spec]], function(z) z$yieldDiff[1,c(hcrsMSY,hcrsHS)])))
        plot(1,1,ty='n', axes = FALSE,xlim = limA, ylim = limA,
             ylab='',xlab='')
        abline(0,1,col="grey70",lty=1)
        for(scen in 1:nscen){
            xplot <- x[[spec]][[scen]]$yieldDiff[1:3,hcrsMSY]
            yplot <- x[[spec]][[scen]]$yieldDiff[1:3,hcrsHS]
            lines(xplot[1,c(1:8)], yplot[1,c(1:8)], lty = 3,
                  col = tail(colA,1))
            lines(xplot[1,c(1,9:17)], yplot[1,c(1,9:17)], lty = 3,
                  col = tail(colC,1))
            points(xplot[1,], yplot[1,], cex = cex,
                   pch = c(15,16,17,18)[scen], col = c("grey60",colA,colC))
        }
        axis(1)
        axis(2)
        if(spec == 2) mtext("AAV",1,3)
        if(spec == 1){
            mtext("AAV",2,3)
        }
        box()
    }
    mtext("HS", 2, 3.5, outer = TRUE, font = 2, cex = 1.2)
    mtext("MSY", 1, 3.5, outer = TRUE, font = 2, adj = 0.45, cex = 1.2)
    ## legend
    plot.new()
    legend(0.2,0.9, legend=c("Med","A45","A40","A35","A30","A25","A15","A05",
                             "C45","C40","C35","C30","C25","C15","C05","C01","C001"),
           title.adj=0.1,
           pch=16,
           title = "HCR",
           col=c("grey60",colA, colC),
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)
    legend(0.2,0.3, legend=paste0("S",1:nscen),
           title.adj=0.1,
           pch=c(15,16,17,18),
           title = "Scenario",
           col="black",
           y.intersp=1.3,
           ncol = 1, bty="n", cex=1.3)

}
