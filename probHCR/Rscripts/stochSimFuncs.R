## Define functions for WKLife 8
##--------------------------------------------------------------------------------

require(RColorBrewer) ## display.brewer.all()

#' @name get.MP
#' @title Get the management procedure (harvest control rule, advice rule)
#' @details This function creates harvest control rules (HCRs) which can be incorporated into a
#' management strategy evaluation framework (DLMtool package). HCRs are saved with a
#' generic name to the global environment and the names of the HCRs are returned if results of the
#' function are assigned to an object. HCR runs a SPiCT assessment using catch and
#' relative biomass index observations. Stock status estimates are used to set the TAC
#' for the next year. TAC can be based on the distribution of predicted catches (percentileC)
#' and/or the distribution of the Fmsy reference level (percentileFmsy).
#' Additionally, a cap can be applied to account for low biomass levels (below Bmsy).
#' Arguments of returned function are 'x' - the position in a data-limited mehods data object,
#' 'Data' - the data-limited methods data object (see DLMtool), and 'reps' - the number of
#' stochastic samples of the TAC recommendation (not used for this HCR).
#' One or several arguments of the function can be provided as vectors to generate several
#' HCRs at once (several vectors have to have same length).
#'
#' @param fractileC The fractile of the catch distribution to be used for setting TAC. Default
#'   is median (0.5).
#' @param fractileFFmsy The fractile of the distribution of F/Fmsy. Default is 0.5 (median).
#' @param fractileBBmsy The fractile of the distribution of B/Bmsy. Default is 0.5 (median).
#' @param pa Logical; indicating if the precautionary approach should be applied (reduce F if P(B<Blim) < prob). Default is FALSE.
#' @param prob Probability for the precautionary approach (see argument 'pa', default is 0.95).
#' @param bbmsyfrac
#' @param stabilityClause Logical; If true TAC is bound between two values set in lower and upper. Default: FALSE.
#' @param lower lower bound of the stability clause. Default is 0.8, used if stabilityClause = TRUE.
#' @param upper upper bound of the stability clause. Default is 1.2, used if stabilityClause = TRUE.
#' @param amtint Assessment interval. Default is 1, which indicates annual assessments.
#' @param npriorSD standard deviation of logn prior (Default: 2). If NA, the logn prior is removed
#' @param nhist number of historic years to use for assessment (default = NA, which means to use all available years)
#' @param env environment where the harvest control rule function(s) are assigned to.
#' @param package
#' @return A function which can estimate TAC recommendations based on SPiCT assessment,
#'   taking assessment uncertainty into account.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Put together an operating model from the available DLM toolkit examples
#' StockEx <- Herring
#' FleetEx <- Generic_IncE
#' ObsEx <- Precise_Unbiased
#' ## Remove changes in life history parameters
#' StockEx@Mgrad <- c(0,0)
#' StockEx@Kgrad <- c(0,0)
#' StockEx@Linfgrad <- c(0,0)
#' StockEx@Size_area_1 <- c(0.5,0.5)
#' StockEx@Frac_area_1 <- c(0.5,0.5)
#' StockEx@Prob_staying <- c(0.5,0.5)
#' ## Set the depletion level
#' StockEx@D <- c(0.3, 0.4)
#' ## create Operation Model
#' OMex <- new("OM", Stock = StockEx, Fleet = FleetEx,
#'                 Obs = ObsEx)
#' ## Set simulation options
#' OMex@nsim <- 10
#' OMex@nyears <- 25
#' OMex@proyears <- 5
#' ## Get SPiCT HCRs
#' MPname <- c(spict2DLMtool(),
#'             spict2DLMtool(pa=TRUE))
#' ## run MSE
#' MSEex <- runMSE(OMex,
#'                 MPs = MPname,
#'                 timelimit = 150,
#'                 CheckMPs = FALSE)
#' ## example plot of results
#' Pplot2(MSEex, traj="quant", quants=c(0.2, 0.8))
#' }
#'
get.MP <- function(fractileC = 0.5,
                   fractileFFmsy = 0.5,
                   fractileBBmsy = 0.5,
                   breakpointB=0.5,
                   amtint = 1,
                   dteuler = 1/16,
                   npriorSD = 2,
                   nhist = NA,
                   env = globalenv(),
                   package="dlmtool"){
    ## allowing for multiple generation of MPs
    argList <- list(fracc=fractileC, fracf=fractileFFmsy,
                    fracb=fractileBBmsy, breakpointB=breakpointB,
                    dteuler=dteuler,
                    npriorSD=npriorSD)
    template  <- expression(paste0(
        'structure(function(x, Data, reps = 1,
          fractileC = ',a,',
          fractileFFmsy=',b,',
          fractileBBmsy=',c,',
          breakpointB=',d,',
          npriorSD=',npriorSD,',
          nhist=',nhist,'){
            dependencies <- "Data@Year, Data@Cat, Data@VInd, Data@Misc, Data@MPrec, Data@LHYear"
            ## limit data to use for assessment (with argument nhist)
            LH <- Data@LHYear   ## last historical year
            if(is.na(nhist)) nhist <- LH else nhist <- nhist     ## number of historical data years to use
            curr.yr <- length(Data@Year)
            fst.yr <- max(LH-nhist+1, 1)  ## 1 ### why 1??
            yr.ind <-  fst.yr:curr.yr
            n.yr <- length(yr.ind)
            time <- yr.ind
            Catch <- Data@Cat[x,yr.ind]
            Index <- Data@VInd[x,yr.ind]
##            Index <- Data@AddInd[x,1,yr.ind]
            ## last years TAC (in first year = Clast)
            taclast <- Data@MPrec[x]
            ## compile inp list
            inp <- list(timeC=time, obsC=Catch,
                        timeI=time, obsI=Index,
                        dteuler = ',e,')
            ## inp$do.sd.report <- TRUE
            ## inp$getReportCovariance <- FALSE
            inp$reportmode <- 1
            ## assessment + advice year
            manstart = inp$timeC[length(inp$timeC)] + 2
            inp$maninterval <- c(manstart, manstart + 1)
            ## check inp
            inp <- check.inp(inp, verbose = FALSE)
            inp$optimiser.control <- list(iter.max = 1e3, eval.max = 1e3)
            ## stronger prior
            if(!is.null(npriorSD) && is.finite(npriorSD)){
               inp$priors$logn <- c(log(2),npriorSD,1)
            }else{
               ## inp$priors$logn <- c(0,0,0)
               ## inp$priors$logalpha <- c(0,0,0)
               ## inp$priors$logbeta <- c(0,0,0)
            }
            ## use Thorsons prior
            ## inp$priors$logn <- c(log(1.478),sqrt(0.849^2/1.478^2),1)
            ## use more informative n prior
            inp$priors$logn <- c(log(2),1,1)
            ## fit spict
            rep <- try(spict::fit.spict(inp, verbose = FALSE),silent=TRUE)
            ## save reference levels
            bbmsy = try(get.par("logBmBmsy",rep,exp=TRUE),silent=TRUE)
            if(is.numeric(bbmsy)){
                best = as.numeric(bbmsy[,2])
                bsd = as.numeric(bbmsy[,4])
            }else{
                best = NA
                bsd = NA
            }
            ffmsy = try(get.par("logFmFmsynotS",rep,exp=TRUE),silent=TRUE)
            if(is.numeric(ffmsy)){
                fest = as.numeric(ffmsy[,2])
                fsd = as.numeric(ffmsy[,4])
            }else{
                fest = NA
                fsd = NA
            }
            cp = try(get.par("logCp",rep,exp=TRUE),silent=TRUE)
            if(is.numeric(cp)){
                cest = as.numeric(cp[,2])
                csd = as.numeric(cp[,4])
            }else{
                cest = NA
                csd = NA
            }
            msy = try(rep$report$MSY,silent=TRUE)
            msy = ifelse(is.numeric(msy), msy, NA)
            ## stop if not converged
            if(is.null(rep) || is(rep, "try-error") || rep$opt$convergence != 0 ||
                   any(is.infinite(rep$sd))){
               taci <- NA
               conv <- 0
            }else{
               TAC <- try(spict:::get.TAC(rep=rep,
                                          fractiles = list(catch=',a,', ffmsy=',b,',bbmsy=',c,'),
                                          breakpointB = ',d,', intermediatePeriodCatch = taclast, verbose = FALSE),
                          silent=TRUE)
               if(is.null(TAC) || is(TAC, "try-error")){
                  taci <- NA
                  conv <- 1
               }else{
                  if(is.na(TAC)){
                      taci <- NA
                      conv <- 3
                  }else{
                      taci <- TAC
                      conv <- 2
                  }
              }
            }
            Rec <- new("Rec")
            Rec@TAC <- as.numeric(taci)
            Rec@Misc <- list()
            Rec@Misc$spict <- as.numeric(c(unlist(Data@Misc[[x]]$spict), conv))
            Rec@Misc$TAC <- as.numeric(c(unlist(Data@Misc[[x]]$TAC), as.numeric(taci)))
            Rec@Misc$spictBest <- as.numeric(c(unlist(Data@Misc[[x]]$spictBest), best))
            Rec@Misc$spictBsd <- as.numeric(c(unlist(Data@Misc[[x]]$spictBsd), bsd))
            Rec@Misc$spictFest <- as.numeric(c(unlist(Data@Misc[[x]]$spictFest), fest))
            Rec@Misc$spictFsd <- as.numeric(c(unlist(Data@Misc[[x]]$spictFsd), fsd))
            Rec@Misc$spictCest <- as.numeric(c(unlist(Data@Misc[[x]]$spictCest), cest))
            Rec@Misc$spictCsd <- as.numeric(c(unlist(Data@Misc[[x]]$spictCsd), csd))
            Rec@Misc$spictMSY <- as.numeric(c(unlist(Data@Misc[[x]]$spictMSY), msy))
            return(Rec)
        },
        class="MP")'))

    ## create MPs as functions
    nami <- NA
    I=1
    argListCor <- argList
    subList <- lapply(argListCor, "[[", I)
    names(subList) <- letters[1:6]
    templati <- eval(parse(text=paste(parse(text = eval(template, subList)),collapse=" ")))

    ## save names of MPs
    c0 <- "MSY"
    if(argListCor[[1]][I] == 0.5){
        c1 <- ""
    }else{
        c1 <- paste0("_C",argListCor[[1]][I])
    }
    if(argListCor[[2]][I] == 0.5){
        c2 <- ""
    }else{
        c2 <- paste0("_FFmsy",argListCor[[2]][I])
    }
    if(argListCor[[3]][I] == 0.5){
        c3 <- ""
    }else{
        c3 <- paste0("_BBmsy",argListCor[[3]][I])
    }
    if(argListCor[[4]][I] == 0){
        c4 <- ""
    }else{
        c4 <- paste0("_bpB",argListCor[[4]][I])
    }
    if(argListCor[[5]][I] == 1){
        c5 <- ""
    }else{
        c5 <- paste0("_dt",round(argListCor[[5]][I],2))
    }
    c6 <- paste0("_nsd",round(npriorSD,2))

    ## put everythin together
    nami[I] <- paste0("spict",c0,c1,c2,c3,c4,c5,c6)
    assign(value=templati, x=nami[I], envir=env)

    ## allow for assigning names
    invisible(nami)
}


## helper and wrapper functions
## ----------------------------------------------------------------------------------------
## get limits for y and x axes
getLimits <- function (MSEobj, PM, remove23=FALSE, risk=1, yrs=NULL){
    if(remove23){
        MPs <- 3:(MSEobj@nMPs-2)
        MSEobj <- Sub(MSEobj, MPs)
    }
    if(is.null(yrs)){
        return(eval(call(PM, MSEobj))@Mean)
    }else{
        return(eval(call(PM, MSEobj,Yrs=yrs))@Mean)
    }
}


## subset all simulation with certain fraction of NA
subNAsims <- function(MSE, frac=0.1){
    tmp <- apply(MSE@TAC, c(1,2), function(x) sum(is.na(x)))
    tmp2 <- tmp/MSE@proyears
    simsNotNA <- apply(tmp2, 1, function(x) all(x < frac))
    MSE2 <- Sub(MSEobj = MSE, sims = simsNotNA)
    return(MSE2)
}




run_parallelTKM <- function (i, itsim, OM, MPs, CheckMPs, timelimit, Hist, ntrials,
                             fracD, CalcBlow, HZN, Bfrac, AnnualMSY, silent, PPD, control,
                             parallel = FALSE){
    OM@nsim <- itsim[i]
    OM@seed <- OM@seed + i
    sfCat(paste("Seed ", OM@seed), sep="\n")
    mse <- DLMtool:::runMSE_int(OM, MPs, CheckMPs, timelimit, Hist, ntrials,
                                fracD, CalcBlow, HZN, Bfrac, AnnualMSY, silent, PPD = PPD,
                                control = control, parallel = parallel)
    return(mse)
}

runMSETKM <- function (OM = DLMtool::testOM, MPs = c("AvC", "DCAC", "FMSYref",
                                                     "curE", "matlenlim", "MRreal"), CheckMPs = FALSE, timelimit = 1,
                       Hist = FALSE, ntrials = 50, fracD = 0.05, CalcBlow = TRUE,
                       HZN = 2, Bfrac = 0.5, AnnualMSY = TRUE, silent = FALSE, PPD = FALSE,
                       parallel = FALSE, save_name = NULL, checks = FALSE, control = NULL){
    if (class(OM) != "OM")
        stop("OM is not class 'OM'", call. = FALSE)
    if (Hist & parallel) {
        message("Sorry! Historical simulations currently can't use parallel.")
        parallel <- FALSE
    }
    rm(list = ls(DLMenv), envir = DLMenv)
    tt <- suppressWarnings(try(lsf.str(envir = globalenv()),
                               silent = TRUE))
    if (class(tt) != "try-error") {
        gl.funs <- as.vector(tt)
        pkg.funs <- as.vector(ls.str("package:DLMtool"))
        if ("package:MSEtool" %in% search())
            pkg.funs <- c(pkg.funs, as.vector(ls.str("package:MSEtool")))
        if (length(gl.funs) > 0) {
            gl.clss <- unlist(lapply(lapply(gl.funs, get), class))
            gl.MP <- gl.funs[gl.clss %in% "MP"]
            if (length(gl.MP) > 0) {
                inc.gl <- gl.MP[gl.MP %in% MPs]
                if (length(inc.gl) > 0) {
                    dup.MPs <- inc.gl[inc.gl %in% pkg.funs]
                    if (length(dup.MPs) > 0) {
                        stop("Custom MP names already in DLMtool: ",
                             paste0(dup.MPs, " "), "\nRename Custom MPs")
                    }
                }
            }
        }
    }
    if (!all(is.na(MPs))) {
        for (mm in MPs) {
            chkMP <- try(get(mm), silent = TRUE)
            if (class(chkMP) != "MP")
                stop(mm, " is not a valid MP", call. = FALSE)
        }
    }
    if (parallel) {
        if (!snowfall::sfIsRunning()) {
            message("Parallel processing hasn't been initialized. Calling 'setup()' now")
            setup(cpu = nCores, slaveOutfile=paste0("log/log_",species, "_scen",scenNum,"_nsim",
                                                    nsim, "_dtE",dteuler,checkLab,".txt"))
            snowfall::sfExport(list=list("sample_index", "run_parallelTKM"))
        }
        if (all(is.na(MPs)))
            MPs <- avail("MP")
        cMPs <- MPs[!MPs %in% pkg.funs]
        if (length(cMPs) > 0)
            snowfall::sfExport(list = cMPs)
        ncpu <- snowfall::sfCpus()
        print(ncpu)
        if (OM@nsim < 48)
            stop("nsim must be >=48")
        nits <- ceiling(OM@nsim/48)
        itsim <- rep(48, nits)
        if (nits < ncpu) {
            nits <- ncpu
            itsim <- rep(ceiling(OM@nsim/ncpu), ncpu)
        }
        cnt <- 1
        while (sum(itsim) != OM@nsim | any(itsim < 2)) {
            diff <- OM@nsim - sum(itsim)
            if (diff > 0) {
                itsim[cnt] <- itsim[cnt] + 1
            }
            if (diff < 0) {
                itsim[cnt] <- itsim[cnt] - 1

            }
            cnt <- cnt + 1
            if (cnt > length(itsim))
                cnt <- 1
        }
        if (!silent)
            message("Running MSE in parallel on ", ncpu, " processors")
        temp <- snowfall::sfClusterApplyLB(1:nits, run_parallelTKM,
                                           itsim = itsim, OM = OM, MPs = MPs, CheckMPs = CheckMPs,
                                           timelimit = timelimit, Hist = Hist, ntrials = ntrials,
                                           fracD = fracD, CalcBlow = CalcBlow, HZN = HZN, Bfrac = Bfrac,
                                           AnnualMSY = AnnualMSY, silent = TRUE, PPD = PPD,
                                           control = control, parallel = parallel)
        if (!is.null(save_name) && is.character(save_name))
            saveRDS(temp, paste0(save_name, ".rdata"))
        MSE1 <- joinMSE(temp)
        if (class(MSE1) == "MSE") {
            if (!silent)
                message("MSE completed")
        }
        else {
            warning("MSE completed but could not join MSE objects. Re-run with `save_name ='MyName'` to debug")
        }
    }
    if (!parallel) {
        if (OM@nsim > 48 & !silent & !Hist)
            message("Suggest using 'parallel = TRUE' for large number of simulations")
        MSE1 <- DLMtool:::runMSE_int(OM, MPs, CheckMPs, timelimit, Hist,
                                     ntrials, fracD, CalcBlow, HZN, Bfrac, AnnualMSY,
                                     silent, PPD, checks = checks, control = control)
    }
    return(MSE1)
}



## Extract Stock from OMs
extractStockFromOM <- function(OM, name){
    newStock <- new("Stock")
    stockslots <- slotNames("Stock")
    for(i in 1:length(stockslots)){
        tmp <- try(slot(OM,stockslots[i]), silent=TRUE)
        if(class(tmp) != "try-error"){
            slot(newStock,stockslots[i]) <- tmp
        }
    }
    assign(paste0(name,"Stock"), newStock)
    return(newStock)
}



## customized update OM
updateOM <- function(OMin, pars){
    OM <- Replace(OMin, get(pars$fleetmod))
    OM <- Replace(OM, get(pars$obsmod))
    OM <- Replace(OM, get(pars$impmod))
    ## Assessment interval
    OM@interval <- pars$amtint
    ## Projection years
    OM@proyears <- pars$proyears
    ## Number of simulations
    OM@nsim <- pars$nsim
    ## Number of historical years
    OM@nyears <- pars$nyears
    ## Number of samples of the MP
    OM@reps <- pars$reps
    ## Seed value
    OM@seed <- pars$seed
    ## Depletion level
    OM@D <- pars$d
    ## Make sure no grad changes and 1 area
    OM@Mgrad <- c(0,0)
    OM@Kgrad <- c(0,0)
    OM@Linfgrad <- c(0,0)
    OM@Size_area_1 <- c(0.5,0.5) ## The size of area 1 relative to area 2.
    OM@Frac_area_1 <- c(0.5,0.5) ## The fraction of the unfished biomass in stock 1.
    OM@Prob_staying <- c(0.5,0.5) ## The probability of inviduals in area 1 remaining in area 1 over the course of one year.
    return(OM)
}



getSPB <- function(OMlist, nsim=100, proyears=10){
    bklist <- vector("list",length(OMlist))
    splist <- vector("list",length(OMlist))
    for(i in 1:length(OMlist)){
        OM <- OMlist[[i]]
        OM@nsim <- nsim
        OM@proyears <- proyears
        MSE <- runMSE(OM, MPs="NFref")
        ##
        b <- MSE@SSB_hist
        fm <- MSE@FM_hist
        ## approximate yield
        c <- b * fm
        ## summarize over ages and areas
        bmat <- apply(b, c(1,3), sum, na.rm=TRUE)
        cmat <- apply(c, c(1,3), sum, na.rm=TRUE)
        ##
        spmat <- matrix(NA, nrow=dim(bmat)[1], ncol=dim(bmat)[2])
        for(j in 1:dim(spmat)[2]){
            if(j==dim(spmat)[2]){
                spmat[,j] <- NA
            }else{
                spmat[,j] <- bmat[,(j+1)] - bmat[,j] + cmat[,j]
            }
        }
        tmp <- apply(bmat, 1, max, na.rm=TRUE)
        bkmat <- bmat/tmp
        ##
        bklist[[i]] <- apply(bkmat,2,mean,na.rm=TRUE)
        splist[[i]] <- apply(spmat,2,mean,na.rm=TRUE)
    }
    splist <- lapply(1:length(OMlist),function(x) splist[[x]]/lapply(splist,max,na.rm=TRUE)[[x]])
    return(list(splist=splist,bklist=bklist))
}


## remove certain number of years from MSE (which are influenced by initial values)
removeYears <- function(MSEin, yearsRemove=1:5){
    MSE <- MSEin
    MSE@proyears <- MSEin@proyears - max(yearsRemove)
    elements <- c("B_BMSY", "F_FMSY", "B", "SSB", "VB", "FM", "C", "TAC", "Effort")
    for(i in 1:length(elements)){
        slot(MSE, elements[i]) <- getElement(MSEin, elements[i])[,,-yearsRemove]
    }
    return(MSE)
}





## custom plotting functions
## ----------------------------------------------------------------------------------------
## plot effort of fleet objects
plotFleet <- function(fleet,ylim=c(0,1), xlab="time", ylab="effort", yaxt=NULL, xaxt=NULL){
    y <- fleet@EffYears
    efflo <- fleet@EffLower
    effup <- fleet@EffUpper
    plot(y,rep(1,length(y)), ty='n',ylim=ylim, xlab=xlab, ylab=ylab, yaxt = yaxt, xaxt=xaxt)
    polygon(c(y,rev(y)), c(efflo,rev(effup)), border=NA, col="grey60")
}



## plot surplus production vs. biomass for determination of n
plotSurplus <- function(OMlist){
    bklist <- vector("list",3)
    splist <- vector("list",3)
    for(i in 1:3){
        MSE <- runMSE(OMlist[[i]], MPs="AvC")
        ##
        b <- MSE@SSB_hist
        fm <- MSE@FM_hist
        ## approximate yield
        c <- b * fm
        ## summarize over ages and areas
        bmat <- apply(b, c(1,3), sum, na.rm=TRUE)
        cmat <- apply(c, c(1,3), sum, na.rm=TRUE)
        ##
        spmat <- matrix(NA, nrow=dim(bmat)[1], ncol=dim(bmat)[2])
        for(j in 1:dim(spmat)[2]){
            if(j==dim(spmat)[2]){
                spmat[,j] <- NA
            }else{
                spmat[,j] <- bmat[,(j+1)] - bmat[,j] + cmat[,j]
            }
        }
        tmp <- apply(bmat, 1, max, na.rm=TRUE)
        bkmat <- bmat/tmp
        ##
        bklist[[i]] <- apply(bkmat,2,mean,na.rm=TRUE)
        splist[[i]] <- apply(spmat,2,mean,na.rm=TRUE)
    }
    splist <- lapply(1:3,function(x) splist[[x]]/lapply(splist,max,na.rm=TRUE)[[x]])
    plot(bklist[[i]],splist[[i]], ty='n',lwd=2,
         ylab = "Surplus production", xlab = "B/K",
         xlim=c(0,1),ylim=range(unlist(splist), na.rm=TRUE))
    for(i in 1:3){
        lines(bklist[[i]],splist[[i]],lwd=2, lty=i)
    }
    abline(v=0.5, col="grey70")
}





## plot hist and mse SSB
## trajectories (or quantiles) of SSB over historical and projection years
plotSSB <- function(MSE, MPnum = 1,blim=FALSE,xlab="Year",ylab="SSB/SBBmsy",
                    xaxt=NULL,yaxt=NULL, trendline=TRUE,
                    quant = FALSE, ylim=NA, pro=TRUE){
    nyears <- MSE@nyears
    proyears <- MSE@proyears
    nsim <- MSE@nsim
    ## SSB
    ssbhist <- MSE@SSB_hist ## nsim, nages, nyears, nareas.
    ssb <- MSE@SSB          ## nsim, nMPs, proyears
    ssbhistA <- apply(ssbhist, c(1,3), sum)  ## nsim, nyears
    ## Bmsy (for each sim)
    bmsy <- MSE@OM$SSBMSY
    ssbbmsy <- ssb / bmsy
    ssbhistbmsy <- ssbhistA / bmsy
    if(pro) xlim=c(1,nyears+proyears) else xlim=c(1,nyears)
    ## plot
    if(is.na(ylim[1])) ylim <- range(as.numeric(ssbhistbmsy[,]))
    plot(ssbhistbmsy[1,],ty='n',
         ylab=ylab,xlab=xlab,xaxt=xaxt,yaxt=yaxt,
         xlim=xlim,
         ylim=ylim)
##    if(pro) abline(v=nyears+0.5, col="grey50",lwd=1.5, lty=2)
    if(quant){
        statshist <- apply(ssbhistbmsy[, , drop = FALSE], 2,
                           quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        stats <- apply(ssbbmsy[, MPnum, , drop = FALSE], 3,
                       quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        polygon(x = c(1:(nyears+1), (nyears+1):1),
                y = c(c(statshist[1, ],stats[1,1]),
                      rev(c(statshist[3,],stats[3,1]))),
                col = rgb(t(col2rgb("grey50"))/255,alpha=0.5), border = FALSE)
        lines(1:(nyears+1), c(statshist[2, ],stats[2, 1]), lwd = 3)
        if(pro){
            polygon(x = c((nyears+1):(nyears+proyears), rev((nyears+1):(nyears+proyears))),
                    y = c(stats[1, ], rev(stats[3,])),
                    col = rgb(t(col2rgb("darkgreen"))/255,alpha=0.5),
                    border = FALSE)
            lines((nyears+1):(nyears+proyears), stats[2, ], lwd = 3)
        }
    }else{
        for(i in 1:nsim){
            lines(ssbhistbmsy[i,],lwd=1.3,
                  col=rgb(t(col2rgb("grey50"))/255,alpha=0.7))
            if(pro){
                lines((nyears+1):(nyears+proyears),ssbbmsy[i,MPnum,], lwd=1.3,
                      col=rgb(t(col2rgb("darkgreen"))/255,alpha=0.7))
            }
        }
    }
    abline(h=1,col=1)
    if(trendline) lines(1:(nyears+proyears),c(ssbhistbmsy[3,],ssbbmsy[3,MPnum,]),
                        lwd=1.3,lty=2,
                        col=1)

    if(blim) abline(h=0.5, col="black", lwd=1.5, lty=2)
    box()
}


## plot hist and MSE FM
## trajectories (or quantiles) of FM over historical and projection years
plotFM <- function(MSE, MPnum=1, quant=FALSE, ylim=NA, trendline=TRUE,
                   ylab="F/Fmsy",xlab="Years",xaxt="s",yaxt="s"){
    nyears <- MSE@nyears
    proyears <- MSE@proyears
    nsim <- MSE@nsim
    ## FM
    fmhist <- MSE@FM_hist ## nsim, nages, nyears, nareas.
    fm <- MSE@FM          ## nsim, nMPs, proyears
    fmhistA <- apply(fmhist, c(1,3), mean)  ## nsim, nyears
    ## Fmsy (for each sim)
    fmsy <- MSE@OM$FMSY
    fmfmsy <- fm / fmsy
    fmhistfmsy <- fmhistA / fmsy
    ## plot
    if(is.na(ylim[1])) ylim <- range(as.numeric(fmhistfmsy[,]))
    plot(fmhistfmsy[1,],ty='n',
         ylab=ylab,xlab=xlab,
         xlim=c(1,nyears+proyears),
         ylim=ylim,xaxt=xaxt,yaxt=yaxt)
    if(quant){
        statshist <- apply(fmhistfmsy[, , drop = FALSE], 2,
                           quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        stats <- apply(fmfmsy[, MPnum, , drop = FALSE], 3,
                       quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        polygon(x = c(1:(nyears), (nyears):1),
                y = c(c(statshist[1, ]), rev(c(statshist[3,]))),
                col = rgb(t(col2rgb("grey50"))/255,alpha=0.5), border = FALSE)
        lines(1:(nyears+1), c(statshist[2, ], stats[2,1]), lwd = 3)
        polygon(x = c((nyears+1):(nyears+proyears), rev((nyears+1):(nyears+proyears))),
                y = c(stats[1, ], rev(stats[3,])),
                col = rgb(t(col2rgb("dodgerblue4"))/255,alpha=0.5),
                border = FALSE)
        lines((nyears+1):(nyears+proyears), stats[2, ], lwd = 3)
        abline(v=nyears, col="grey60",lty=2,lwd=1.5)
    }else{
        for(i in 1:nsim){
            lines(fmhistfmsy[i,],lwd=1.3,
                  col=rgb(t(col2rgb("grey50"))/255,alpha=0.7))
            lines((nyears+1):(nyears+proyears),fmfmsy[i,MPnum,], lwd=1.3,
                  col=rgb(t(col2rgb("darkred"))/255,alpha=0.7))
        }
    }
    abline(h=1, col="black", lwd=1.5)
    if(trendline) lines(0:(nyears+proyears-1),
                        c(fmhistfmsy[3,,drop = FALSE],fmfmsy[3,MPnum,]),
                        lwd=1.3,lty=2,
                        col="black")
    box()
}


plotY <- function(MSE, MPnum=1, quant=FALSE, ylim=NA, trendline=TRUE,
                  ylab="Rel. catch",xlab="Years",xaxt="s",yaxt="s"){
    nyears <- MSE@nyears
    proyears <- MSE@proyears
    nsim <- MSE@nsim
    ## Yield
    cahist <- MSE@CB_hist ## nsim, nages, nyears, nareas.
    ca <- MSE@C          ## nsim, nMPs, proyears
    cahistA <- apply(cahist, c(1,3), sum)  ## nsim, nyears
    ## Catch (for each sim)
    refy <- MSE@OM$RefY
    carefy <- ca / refy
    cahistrefy <- cahistA / refy
    ## plot
    if(is.na(ylim[1])) ylim <- range(as.numeric(cahistrefy[,]),as.numeric(carefy[,,]))
    plot(cahistrefy[1,],ty='n',
         ylab=ylab,xlab=xlab,
         xlim=c(1,nyears+proyears),
         ylim=ylim,xaxt=xaxt,yaxt=yaxt)
    if(quant){
        statshist <- apply(cahistrefy[, , drop = FALSE], 2,
                           quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        stats <- apply(carefy[, MPnum, , drop = FALSE], 3,
                       quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        polygon(x = c(1:(nyears+1), (nyears+1):1),
                y = c(c(statshist[1, ],stats[1,1]), rev(c(statshist[3,],stats[3,1]))),
                col = rgb(t(col2rgb("grey50"))/255,alpha=0.5), border = FALSE)
        lines(1:(nyears+1), c(statshist[2, ],stats[2,1]), lwd = 3)
        polygon(x = c((nyears+1):(nyears+proyears), rev((nyears+1):(nyears+proyears))),
                y = c(stats[1, ], rev(stats[3,])), col = rgb(t(col2rgb("dodgerblue4"))/255,alpha=0.5),
                border = FALSE)
        lines((nyears+1):(nyears+proyears), stats[2, ], lwd = 3)
    }else{
        for(i in 1:nsim){
            lines(cahistrefy[i,],lwd=1.3,
                  col=rgb(t(col2rgb("grey50"))/255,alpha=0.7))
            lines((nyears+1):(nyears+proyears),carefy[i,MPnum,], lwd=1.3,
                  col=rgb(t(col2rgb("dodgerblue4"))/255,alpha=0.7))
        }
    }
    abline(h=1, col="black", lwd=1)
    if(trendline) lines(1:(nyears+proyears),c(cahistrefy[3,,drop = FALSE],carefy[3,MPnum,]),
                        lwd=1.3,lty=2,
                        col=1)
    box()
}



plotTrade <- function (MSEobj, PMlist = c("PBBlim","YieldTKM"), risk=1,
                       Title = NULL,
                       Yrs = NULL, Refs = NULL,
                       xlim=NULL, ylim=NULL,
                       xcap=NULL,ycap=NULL){

    if (class(MSEobj) != "MSE")
        stop("Object must be class `MSE`", call. = FALSE)
    if (is.null(PMlist)) {
        PMlist <- unlist(list(...))
    }else{
        PMlist <- unlist(PMlist)
    }
    if (length(PMlist) == 0)
        PMlist <- c("P50", "YieldTKM")
    if (class(PMlist) != "character")
        stop("Must provide names of PM methods")
    for (X in seq_along(PMlist)) if (!PMlist[X] %in% avail("PM"))
                                     stop(PMlist[X], " is not a valid PM method")
    if (length(PMlist) < 2)
        stop("Must provided more than 1 PM method")
    nPMs <- length(PMlist)
    if (nPMs%%2 != 0) {
        message("Odd number of PMs. Recycling first PM")
        PMlist <- c(PMlist, PMlist[1])
        nPMs <- length(PMlist)
    }
    runPM <- vector("list", length(PMlist))
    for (X in 1:length(PMlist)) {
        ref <- Refs[[PMlist[X]]]
        yrs <- Yrs##[[PMlist[X]]]
        if (is.null(ref)) {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs))
                }
            }
        }else {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk,
                                            Yrs = yrs))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref,
                                            Yrs = yrs))
                }
            }
        }
    }
    nplots <- nPMs/2
    n.col <- ceiling(sqrt(nplots))
    n.row <- ceiling(nplots/n.col)
    ##    m <- matrix(1:(n.col * n.row), ncol = n.col, nrow = n.row,
    ##        byrow = FALSE)
    xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
    plots <- listout <- list()
    xInd <- seq(1, by = 2, length.out = nplots)
    yInd <- xInd + 1
    for (pp in 1:nplots) {
        yPM <- PMlist[yInd[pp]]
        yvals <- runPM[[match(yPM, PMlist)]]@Mean
        ycap <- runPM[[match(yPM, PMlist)]]@Caption
        yname <- runPM[[match(yPM, PMlist)]]@Name
        xPM <- PMlist[xInd[pp]]
        xvals <- runPM[[match(xPM, PMlist)]]@Mean
        xcap <- runPM[[match(xPM, PMlist)]]@Caption
        xname <- runPM[[match(xPM, PMlist)]]@Name
        if(is.null(xlim)) xlim <- range(xvals)
        if(is.null(ylim)) ylim <- range(yvals)
        if(yPM == "AAVY"){
            ycap <- paste0("P(AAVY<0.2) [years ",yrs[1],"-",yrs[2],"]")
        }
        ##        MPType <- MPtype(MSEobj@MPs)
        ##        Class <- MPType[match(MSEobj@MPs, MPType[, 1]), 2]
        df <- data.frame(x = xvals, y = yvals, label = MSEobj@MPs,
                         ##            Class = Class,
                         fontface = "plain", xPM = xPM, yPM = yPM)
        df$fontface <- as.character(df$fontface)
        df$fontface <- factor(df$fontface)
        listout[[pp]] <- df
        plots[[pp]] <- ggplot2::ggplot() +
            ggplot2::geom_point(data = df, ggplot2::aes(x, y),
                                size = 2) + ggrepel::geom_text_repel(data = df,
                                                                     ggplot2::aes(x, y, label = label,
                                                                                  fontface = 1), show.legend = FALSE, size = 4,
                                                                     na.rm = TRUE) + ggplot2::theme_bw() + ggplot2::xlab(xcap) + ggplot2::ylab(ycap) +
                                                     ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
                                                     ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                                                                    axis.title = ggplot2::element_text(size=15),
                                                                    plot.margin = unit(c(1,1,1,1), "cm"),
                                                                    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 25,
                                                                                                                                  r = 0, b = 0, l = 0)),
                                                                    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0,
                                                                                                                                  r = 25, b = 0, l = 0)),
                                                                    plot.title = ggplot2::element_text(hjust = 0.5, face="bold"))+
                                                     ggplot2::ggtitle(Title)
    }
    invisible(plots[[pp]])
}



plotTrade2 <- function (MSEobj, PMlist = c("PBBlim","YieldTKM"), risk=1,
                        numMSY = 6, numPA = 10,
                        Title = NULL,MPs=NA,
                        Yrs = NULL, Refs = NULL,
                        xlim=NULL, ylim=NULL,
                        xcap=NULL,ycap=NULL, remove23=FALSE){
    if(remove23) MPs <- 3:(MSEobj@nMPs-2)
    if (!all(is.na(MPs)))
        MSEobj <- Sub(MSEobj, MPs = MPs)
    if (class(MSEobj) != "MSE")
        stop("Object must be class `MSE`", call. = FALSE)
    if (is.null(PMlist)) {
        PMlist <- unlist(list(...))
    }else{
        PMlist <- unlist(PMlist)
    }
    if (length(PMlist) == 0)
        PMlist <- c("P50", "YieldTKM")
    if (class(PMlist) != "character")
        stop("Must provide names of PM methods")
    for (X in seq_along(PMlist)) if (!PMlist[X] %in% avail("PM"))
                                     stop(PMlist[X], " is not a valid PM method")
    if (length(PMlist) < 2)
        stop("Must provided more than 1 PM method")
    nPMs <- length(PMlist)
    if (nPMs%%2 != 0) {
        message("Odd number of PMs. Recycling first PM")
        PMlist <- c(PMlist, PMlist[1])
        nPMs <- length(PMlist)
    }
    runPM <- vector("list", length(PMlist))
    for (X in 1:length(PMlist)) {
        ref <- Refs[[PMlist[X]]]
        yrs <- Yrs##[[PMlist[X]]]
        if (is.null(ref)) {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs))
                }
            }
        }else {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk,
                                            Yrs = yrs))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref,
                                            Yrs = yrs))
                }
            }
        }
    }
    nplots <- nPMs/2
    n.col <- ceiling(sqrt(nplots))
    n.row <- ceiling(nplots/n.col)
    ##    m <- matrix(1:(n.col * n.row), ncol = n.col, nrow = n.row,
    ##        byrow = FALSE)
    xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
    plots <- listout <- list()
    xInd <- seq(1, by = 2, length.out = nplots)
    yInd <- xInd + 1
    ##
    yPM <- PMlist[yInd[1]]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    if(is.null(ycap)) ycap <- runPM[[match(yPM, PMlist)]]@Caption else ycap <- ycap
    yname <- runPM[[match(yPM, PMlist)]]@Name
    xPM <- PMlist[xInd[1]]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    if(is.null(xcap)) xcap <- runPM[[match(xPM, PMlist)]]@Caption else xcap <- xcap
    xname <- runPM[[match(xPM, PMlist)]]@Name
    if(is.null(xlim)) xlim <- range(xvals)
    if(is.null(ylim)) ylim <- range(yvals)
    if(yPM == "AAVY"){
        ycap <- paste0("P(AAVY<0.2) [years ",yrs[1],"-",yrs[2],"]")
    }
    ##    if(yPM == "CVyield"){
    ##        ycap <- paste0(" [years ",yrs[1],"-",yrs[2],"]")
    ##    }
    df <- data.frame(x = xvals, y = yvals, label = MSEobj@MPs,
                     fontface = "plain", xPM = xPM, yPM = yPM)
    df$fontface <- as.character(df$fontface)
    df$fontface <- factor(df$fontface)
    ##
    ##    numMSY <- 6
    ##    numPA <- 10
    ##
##    browser()

    if(remove23){
        pchs <- c(17, 2, 6, 15)
        argroup <- c(rep(1,numMSY), 2,3, rep(4,numPA))
        arcols <- c(paste0(rep("grey",numMSY),round(seq(75,10,length.out = numMSY))),1,1,
                    paste0(rep("grey",numPA),round(seq(75,10,length.out = numPA))))
        aralphas <- c(seq(0.3,1,length.out = numMSY),1,1, seq(0.3,1,length.out = numPA))
        arcex <- c(rep(1.2,numMSY),rep(1.25,2),rep(1.2,numPA),1.2)
        ind1 <- 1:(numMSY)
        ind2 <- (numMSY+3):(length(argroup))
        ind3 <- c(6,4,7)  ## with different stab clauses
    }else{
        pchs <- c(3,4, 17,2,6, 15,1,10)
        argroup <- c(1,2,rep(3,numMSY),4,5, rep(6,numPA),7,8)
        arcols <- c(1,1,paste0(rep("grey",numMSY),round(seq(75,10,length.out = numMSY))),1,1,
                    paste0(rep("grey",numPA),round(seq(75,10,length.out = numPA))),1,1)
        aralphas <- c(rep(1,2), seq(0.3,1,length.out = numMSY),1,1,  seq(0.3,1,length.out = numPA),1,1)
        arcex <- c(rep(1.2,2),rep(1.2,numMSY),rep(1.25,2),rep(1.2,numPA),1.2,1.2)
        ind1 <- 3:(numMSY+2)
        ind2 <- (numMSY+5):(length(argroup)-2)
        ind3 <- c(8,6,9)  ## with different stab clauses
    }


    ##    print(MSEobj@MPs)
    plot(xvals, yvals, ty='n',main=Title,
         xlab = xcap, ylab = ycap, xlim=xlim, ylim=ylim)
    if(xPM == "PBBlim") abline(v=0.05, col="grey70", lty=2, lwd=1.4)
    ##    abline(v=xvals[length(xvals)], col="grey70", lty=3, lwd=1.4)
    ##    text(xvals,yvals,labels=MSEobj@MPs)
    points(xvals, yvals,lwd=1.3,
           pch=pchs[argroup],
           col=arcols,##rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           cex=arcex)
    points(xvals[ind1], yvals[ind1], ty='b', pch=NA, col="grey60", lwd=1.4)
    points(xvals[ind2], yvals[ind2], ty='b', pch=NA, col="grey60", lwd=1.4)
    points(xvals[ind3], yvals[ind3], ty='b', pch=NA, col="grey60", lwd=1.4,
           lty=3)
    box()

    if(FALSE){

        pchs <- c(3,4, 17,15)
        argroup <- c(1,2,3,rep(4,10))
        aralphas <- c(rep(1,2), 1, seq(0.01,1,length.out = 10))
        ind2 <- (4):length(argroup)
        ##    print(MSEobj@MPs)
        plot(xvals, yvals, ty='n',main=Title,
             xlab = xcap, ylab = ycap)
        abline(v=0.05, col="grey70", lty=2, lwd=1.4)
        ##    text(xvals,yvals,labels=MSEobj@MPs)
        points(xvals, yvals,
               pch=pchs[argroup],
               col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
               cex=1.2)
        ##    points(xvals[ind1], yvals[ind1], ty='b', pch=NA, col="grey60", lwd=1.4)
        points(xvals[ind2], yvals[ind2], ty='b', pch=NA, col="grey60", lwd=1.4)
        dev.print(pdf,"PAapproachingMSY05,pdf")
        ##    getwd()


        ## checkEuler

        pchs <- c(17,15)
        argroup <- c(rep(1,6),rep(2,3))
        aralphas <- rep(1,9)##c(seq(0.01,1,length.out = 10))
        ##        ind2 <- (4):length(argroup)
        ##    print(MSEobj@MPs)
        plot(xvals, yvals, ty='n',main=Title,
             xlab = xcap, ylab = ycap)
        abline(v=0.05, col="grey70", lty=2, lwd=1.4)
        text(xvals+0.002,yvals+0.0015,labels=MSEobj@MPs)
        points(xvals, yvals,
               pch=pchs[argroup],
               col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
               cex=1.2)
        ##    points(xvals[ind1], yvals[ind1], ty='b', pch=NA, col="grey60", lwd=1.4)
        ##        points(xvals[ind2], yvals[ind2], ty='b', pch=NA, col="grey60", lwd=1.4)
        ##        dev.print(pdf,"PAapproachingMSY05,pdf")
        ##    getwd()

    }

}


plotTrade2E <- function (MSEobj, PMlist = c("PBBlim","YieldTKM"), risk=1,
                         numMSY = 6, numPA = 10,
                         Title = NULL,MPs=NA,
                         Yrs = NULL, Refs = NULL,
                         xlim=NULL, ylim=NULL,
                         xcap=NULL,ycap=NULL, remove23=FALSE){
    if(remove23) MPs <- 3:MSEobj@nMPs
    if (!all(is.na(MPs)))
        MSEobj <- Sub(MSEobj, MPs = MPs)
    if (class(MSEobj) != "MSE")
        stop("Object must be class `MSE`", call. = FALSE)
    if (is.null(PMlist)) {
        PMlist <- unlist(list(...))
    }else{
        PMlist <- unlist(PMlist)
    }
    if (length(PMlist) == 0)
        PMlist <- c("P50", "YieldTKM")
    if (class(PMlist) != "character")
        stop("Must provide names of PM methods")
    for (X in seq_along(PMlist)) if (!PMlist[X] %in% avail("PM"))
                                     stop(PMlist[X], " is not a valid PM method")
    if (length(PMlist) < 2)
        stop("Must provided more than 1 PM method")
    nPMs <- length(PMlist)
    if (nPMs%%2 != 0) {
        message("Odd number of PMs. Recycling first PM")
        PMlist <- c(PMlist, PMlist[1])
        nPMs <- length(PMlist)
    }
    runPM <- vector("list", length(PMlist))
    for (X in 1:length(PMlist)) {
        ref <- Refs[[PMlist[X]]]
        yrs <- Yrs##[[PMlist[X]]]
        if (is.null(ref)) {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs))
                }
            }
        }else {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk,
                                            Yrs = yrs))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref,
                                            Yrs = yrs))
                }
            }
        }
    }
    nplots <- nPMs/2
    n.col <- ceiling(sqrt(nplots))
    n.row <- ceiling(nplots/n.col)
    xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
    plots <- listout <- list()
    xInd <- seq(1, by = 2, length.out = nplots)
    yInd <- xInd + 1
    ##
    yPM <- PMlist[yInd[1]]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    if(is.null(ycap)) ycap <- runPM[[match(yPM, PMlist)]]@Caption else ycap <- ycap
    yname <- runPM[[match(yPM, PMlist)]]@Name
    xPM <- PMlist[xInd[1]]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    if(is.null(xcap)) xcap <- runPM[[match(xPM, PMlist)]]@Caption else xcap <- xcap
    xname <- runPM[[match(xPM, PMlist)]]@Name
    if(is.null(xlim)) xlim <- range(xvals)
    if(is.null(ylim)) ylim <- range(yvals)
    df <- data.frame(x = xvals, y = yvals, label = MSEobj@MPs,
                     fontface = "plain", xPM = xPM, yPM = yPM)
    df$fontface <- as.character(df$fontface)
    df$fontface <- factor(df$fontface)
    ## checkEuler
    pchs <- c(17,16,15)
    argroup <- c(rep(1,4),rep(2,4),rep(3,4))
    aralphas <- rep(seq(0.3,1,length.out = 3),4)
    ##        ind2 <- (4):length(argroup)
    ##    print(MSEobj@MPs)
    plot(xvals, yvals, ty='n',main=Title,
         xlab = xcap, ylab = ycap, ylim=ylim)
    abline(v=0.05, col="grey70", lty=2, lwd=1.4)
    ##    text(xvals+0.002,yvals+0.0015,labels=MSEobj@MPs)
    points(xvals, yvals,
           pch=pchs[argroup],
           col=rep(c("grey70","grey50","grey30","grey10"),3),
           cex=1.2)
    points(xvals[1:4], yvals[1:4], ty='b', pch=NA, col="grey60", lwd=1.4)
    points(xvals[5:8], yvals[5:8], ty='b', pch=NA, col="grey60", lwd=1.4)
    points(xvals[9:13], yvals[9:13], ty='b', pch=NA, col="grey60", lwd=1.4)
}



## for logn analysis
plotTrade2N <- function (MSEobj, PMlist = c("PBBlim","YieldTKM"), risk=1,
                         numMSY = 6, numPA = 10,
                         Title = NULL,MPs=NA,
                         Yrs = NULL, Refs = NULL,
                         xlim=NULL, ylim=NULL,
                         xcap=NULL,ycap=NULL, remove23=FALSE){
    if(remove23) MPs <- 3:(MSEobj@nMPs-2)
    if (!all(is.na(MPs)))
        MSEobj <- Sub(MSEobj, MPs = MPs)
    if (class(MSEobj) != "MSE")
        stop("Object must be class `MSE`", call. = FALSE)
    if (is.null(PMlist)) {
        PMlist <- unlist(list(...))
    }else{
        PMlist <- unlist(PMlist)
    }
    if (length(PMlist) == 0)
        PMlist <- c("P50", "YieldTKM")
    if (class(PMlist) != "character")
        stop("Must provide names of PM methods")
    for (X in seq_along(PMlist)) if (!PMlist[X] %in% avail("PM"))
                                     stop(PMlist[X], " is not a valid PM method")
    if (length(PMlist) < 2)
        stop("Must provided more than 1 PM method")
    nPMs <- length(PMlist)
    if (nPMs%%2 != 0) {
        message("Odd number of PMs. Recycling first PM")
        PMlist <- c(PMlist, PMlist[1])
        nPMs <- length(PMlist)
    }

    PMlist <- c(PMlist, "SDyield") ## "SDPBBlim",
    runPM <- vector("list", length(PMlist))
    for (X in 1:length(PMlist)) {
        ref <- Refs[[PMlist[X]]]
        yrs <- Yrs##[[PMlist[X]]]
        if (is.null(ref)) {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy","SDPBBlim")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy","SDPBBlim")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs))
                }
            }
        }else {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy","SDPBBlim")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy","SDPBBlim")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk,
                                            Yrs = yrs))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref,
                                            Yrs = yrs))
                }
            }
        }
    }



    nplots <- nPMs/2
    n.col <- ceiling(sqrt(nplots))
    n.row <- ceiling(nplots/n.col)
    ##    m <- matrix(1:(n.col * n.row), ncol = n.col, nrow = n.row,
    ##        byrow = FALSE)
    xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
    plots <- listout <- list()
    xInd <- seq(1, by = 2, length.out = nplots)
    yInd <- xInd + 1
    ##
    yPM <- PMlist[yInd[1]]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    ycap <- runPM[[match(yPM, PMlist)]]@Caption
    yname <- runPM[[match(yPM, PMlist)]]@Name
    xPM <- PMlist[xInd[1]]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    xcap <- runPM[[match(xPM, PMlist)]]@Caption
    xname <- runPM[[match(xPM, PMlist)]]@Name
    ## 95% CIs
    if(FALSE){
        y2PM <- PMlist[3]
        ysds <- runPM[[match(y2PM, PMlist)]]@Mean
        ycis <- matrix(c(yvals+2*ysds,yvals-2*ysds),nrow=2,byrow=TRUE)
        x2PM <- PMlist[4]
        xsds <- runPM[[match(x2PM, PMlist)]]@Mean
        xcis <- matrix(c(xvals+2*xsds,xvals-2*xsds),nrow=2,byrow=TRUE)
    }
    if(is.null(xlim)) xlim <- range(xvals)
    if(is.null(ylim)) ylim <- range(yvals)
    if(yPM == "AAVY"){
        ycap <- paste0("P(AAVY<0.2) [years ",yrs[1],"-",yrs[2],"]")
    }
    df <- data.frame(x = xvals, y = yvals, label = MSEobj@MPs,
                     fontface = "plain", xPM = xPM, yPM = yPM)
    df$fontface <- as.character(df$fontface)
    df$fontface <- factor(df$fontface)
    if(remove23){
        pchs <- c(17, 15)
        argroup <- c(rep(1,numMSY), rep(2,numPA))
        arcols <- c(paste0(rep("grey",numMSY),round(seq(75,10,length.out = numMSY))),
                    paste0(rep("grey",numPA),round(seq(75,10,length.out = numPA))))
        aralphas <- c(seq(0.3,1,length.out = numMSY),seq(0.3,1,length.out = numPA))
        arcex <- c(rep(1.2,numMSY),rep(1.2,numPA))
        ind1 <- 1:(numMSY)
        ind2 <- (numMSY+1):(length(argroup))
    }else{
        pchs <- c(3,4, 17,15,1,10)
        argroup <- c(1,2,rep(3,numMSY),rep(4,numPA),5,6)
        arcols <- c(1,1,paste0(rep("grey",numMSY),round(seq(75,10,length.out = numMSY))),
                    paste0(rep("grey",numPA),round(seq(75,10,length.out = numPA))),1,1)
        aralphas <- c(rep(1,2), seq(0.3,1,length.out = numMSY),seq(0.3,1,length.out = numPA),1,1)
        arcex <- c(rep(1.2,2),rep(1.2,numMSY),rep(1.2,numPA),1.2,1.2)
        ind1 <- 3:(numMSY+2)
        ind2 <- (numMSY+3):(length(argroup)-2)
    }


    ##    print(MSEobj@MPs)
    plot(xvals, yvals, ty='n',main=Title,
         xlab = xcap, ylab = ycap, xlim=xlim, ylim=ylim)
    ##    if(xPM == "PBBlim") abline(v=0.05, col="grey70", lty=2, lwd=1.4)
    abline(v=xvals[length(xvals)], col="grey70", lty=3, lwd=1.4)
    ##    text(xvals,yvals,labels=MSEobj@MPs)
    ##    for(i in 1:length(xvals)){
    ##        segments(x0=xcis[2,i],x1=xcis[1,i],y0=yvals[i],y1=yvals[i], col="grey70")
    ##        segments(x0=xvals[i],x1=xvals[i],y0=ycis[2,i],y1=ycis[1,i], col="grey70")
    ##    }
    points(xvals, yvals,lwd=1.3,
           pch=pchs[argroup],
           col=arcols,##rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           cex=arcex)
    points(xvals[ind1], yvals[ind1], ty='b', pch=NA, col="grey40", lwd=1.4)
    points(xvals[ind2], yvals[ind2], ty='b', pch=NA, col="grey40", lwd=1.4)
    box()
}





## tade plot for average over all stocks
plotTradeMean <- function (MSEobj,
                           meanrisk, meanyield,
                           Title = NULL,
                           Yrs = NULL, Refs = NULL,
                           xlim=NULL, ylim=NULL,
                           xcap=NULL,ycap=NULL){
    xcap <- paste0("P(B<Blim) [years ",Yrs[1], "-", Yrs[2], "]")
    ycap <- paste0("Mean Relative Yield [years ", Yrs[1],"-", Yrs[2], "]")
    if(is.null(xlim)) xlim <- range(meanrisk)
    if(is.null(ylim)) ylim <- range(meanyield)
    df <- data.frame(x = meanrisk, y = meanyield, label = MSEobj@MPs,
                     fontface = "plain")
    df$fontface <- as.character(df$fontface)
    df$fontface <- factor(df$fontface)
    plots <- ggplot2::ggplot() +
        ggplot2::geom_point(data = df, ggplot2::aes(x, y),
                            size = 2) + ggrepel::geom_text_repel(data = df,
                                                                 ggplot2::aes(x, y, label = label,
                                                                              fontface = 1), show.legend = FALSE, size = 4,
                                                                 na.rm = TRUE) + ggplot2::theme_bw() + ggplot2::xlab(xcap) + ggplot2::ylab(ycap) +
                                                 ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
                                                 ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                                                                axis.title = ggplot2::element_text(size=15),
                                                                plot.margin = unit(c(1,1,1,1), "cm"),
                                                                axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 25,
                                                                                                                              r = 0, b = 0, l = 0)),
                                                                axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0,
                                                                                                                              r = 25, b = 0, l = 0)),
                                                                plot.title = ggplot2::element_text(hjust = 0.5, face="bold"))+
                                                 ggplot2::ggtitle(Title)
    return(invisible(plots))
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

                                        # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

                                        # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
                                        # Make the panel
                                        # ncol: Number of columns of plots
                                        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
                                        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                                        # Make each plot, in the correct location
        for (i in 1:numPlots) {
                                        # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}


## Kplot alternative
KplotTKM <- function (MSEobj, maxsim = 60, MPs = NA, sims = NULL, maxMP = 9,
                      nam = NA, cex.leg = 1.5){
    if (!is.null(sims) & all(is.na(MPs)))
        MSEobj <- Sub(MSEobj, sims = sims)
    if (!is.null(sims) & all(!is.na(MPs)))
        MSEobj <- Sub(MSEobj, sims = sims, MPs = MPs)
    if (is.null(sims) & !all(is.na(MPs)))
        MSEobj <- Sub(MSEobj, MPs = MPs)
    nMPs <- MSEobj@nMPs
    nsim <- MSEobj@nsim
    if (is.null(sims) & nsim > maxsim)
        MSEobj <- Sub(MSEobj, sims = 1:maxsim)
    if (nMPs > maxMP) {
        message("MSE object has more than ", maxMP, " MPs. Plotting the first ",
                maxMP)
        MSEobj <- Sub(MSEobj, MPs = 1:maxMP)
        nMPs <- MSEobj@nMPs
    }
    nr <- floor((MSEobj@nMPs)^0.5)
    nc <- ceiling((MSEobj@nMPs)/nr)
    Cex <- 1.5
    TitleCex <- 1.5
    FMSYr <- quantile(MSEobj@F_FMSY, c(0.001, 0.9), na.rm = T)
    BMSYr <- quantile(MSEobj@B_BMSY, c(0.001, 0.975), na.rm = T)
    if (is.na(nam))
        op <- par(mfrow = c(nr, nc), mar = c(2, 2, 3, 1), oma = c(3,
                                                                  3.5, 1.2, 0))
    if (!is.na(nam))
        op <- par(mfrow = c(nr, nc), mar = c(2, 2, 3, 1), oma = c(3,
                                                                  3.5, 3, 0))
    colsse <- rainbow(MSEobj@proyears, start = 0.63, end = 0.95)[1:MSEobj@proyears]
    colsse <- makeTransparent(colsse, 95)
    XLim <- c(0, 3)
    YLim <- c(0, 2.5)
    pmat <- matrix(NA, nrow = nc, ncol = nr, byrow = FALSE)
    pmat[1:nMPs] <- 1:nMPs
    pmat <- t(pmat)
    for (mm in 1:MSEobj@nMPs) {
        plot(MSEobj@B_BMSY[1, mm, 1], MSEobj@F_FMSY[1, mm, 1],
             xlim = XLim, ylim = YLim, col = colsse[1], bty = "n",
             axes = FALSE)
        if (nrow(pmat) > 1) {
            if (mm %in% pmat[, 1]) {
                axis(side = 2, labels = TRUE, las = 1)
            }
            else axis(side = 2, labels = FALSE)
            if (mm %in% pmat[nr, ]) {
                axis(side = 1, labels = TRUE)
            }
            else axis(side = 1, labels = FALSE)
            nas <- apply(pmat, 1, sum)
            rr <- which.max(nas)
            nas2 <- apply(pmat, 2, sum)
            cc <- which(is.na(nas2))
            if (mm %in% pmat[rr, cc])
                axis(side = 1, labels = TRUE)
        }
        else {
            if (mm == 1)
                axis(side = 2, labels = TRUE, las = 1)
            if (mm != 1)
                axis(side = 2, labels = FALSE)
            axis(side = 1, labels = TRUE)
        }
        OO <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] <
                        1 & MSEobj@F_FMSY[, mm, MSEobj@proyears] > 1, na.rm = T)/MSEobj@nsim *
                    100, 1)
        OU <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] >
                        1 & MSEobj@F_FMSY[, mm, MSEobj@proyears] > 1, na.rm = T)/MSEobj@nsim *
                    100, 1)
        UO <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] <
                        1 & MSEobj@F_FMSY[, mm, MSEobj@proyears] < 1, na.rm = T)/MSEobj@nsim *
                    100, 1)
        UU <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] >
                        1 & MSEobj@F_FMSY[, mm, MSEobj@proyears] < 1, na.rm = T)/MSEobj@nsim *
                    100, 1)
        abline(h = 1, col = "grey", lwd = 3)
        abline(v = 1, col = "grey", lwd = 3)
        y <- 1:(MSEobj@proyears - 1)
        rng <- 1:min(maxsim, MSEobj@nsim)
        points(MSEobj@B_BMSY[rng, mm, 1], MSEobj@F_FMSY[rng,
                                                        mm, 1], pch = 19, cex = 0.8, col = colsse[1])
        points(MSEobj@B_BMSY[rng, mm, MSEobj@proyears], MSEobj@F_FMSY[rng,
                                                                      mm, MSEobj@proyears], pch = 19, cex = 0.8, col = colsse[MSEobj@proyears])
        if (mm == 1)
            legend("right", c("Start", "End"), bty = "n", text.col = c(colsse[1],
                                                                       colsse[MSEobj@proyears]), pch = 19, col = c(colsse[1],
                                                                                                                   colsse[MSEobj@proyears]))
        legend("topleft", paste(OO, "%", sep = ""), bty = "n",
               text.font = 2, cex = cex.leg)
        legend("topright", paste(OU, "%", sep = ""), bty = "n",
               text.font = 2, cex = cex.leg)
        legend("bottomleft", paste(UO, "%", sep = ""), bty = "n",
               text.font = 2, cex = cex.leg)
        legend("bottomright", paste(UU, "%", sep = ""), bty = "n",
               text.font = 2, cex = cex.leg)
        mtext(MSEobj@MPs[mm], 3, line = 0.6, cex = TitleCex)
    }
    mtext(expression(B/B[MSY]), 1, outer = T, line = 2, cex = Cex)
    mtext(expression(F/F[MSY]), 2, outer = T, line = 1.2, cex = Cex)
    if (!is.na(nam))
        mtext(nam, 3, outer = TRUE, line = 0.25, font = 2, cex = TitleCex)
    par(op)
}


## Customized Performance metrics
calcProbTKM <- function (PM, MSEobj, risk = 1){   ## nsim, nMPs, proyears
    if (MSEobj@nMPs > 1) {   ## nsim, nMPs, proyears
        mar <- 2:3
    }else mar <- 1:2
    ## Risk 1: mean over nsim then mean over years
    if(risk%in%c(1,3)){
        apply(PM, mar, mean, na.rm = TRUE)
    }else if(risk==2){
        apply(PM, mar, function(x) as.numeric(any(x)))
    }
}

calcProbMedian <- function (PM, MSEobj){
    if (MSEobj@nMPs > 1) {
        mar <- 2
    }
    else mar <- 1
    mar <- 1:mar
    apply(PM, mar, median, na.rm = TRUE)
}

calcMedian <- function (Prob){
    if (class(Prob) == "matrix")
        return(apply(Prob, 2, median, na.rm = TRUE))
    if (class(Prob) == "numeric")
        return(median(Prob, na.rm = TRUE))
}

calcMeanTKM <- function (Prob,risk=1){
    if(risk%in%c(1,2)){                             ## Risk 1: average risk over annual risks
        ## Risk 2: prob at least once during ny years
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, mean, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(mean(Prob, na.rm = TRUE))
    }else if(risk==3){                              ## Risk 3: Max risk of annual risks
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, max, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(max(Prob, na.rm = TRUE))
    }
}

calcMedianTKM <- function (Prob,risk=1){
    if(risk%in%c(1,2)){                             ## Risk 1: average risk over annual risks
        ## Risk 2: prob at least once during ny years
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, median, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(median(Prob, na.rm = TRUE))
    }else if(risk==3){                              ## Risk 3: Max risk of annual risks
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, max, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(max(Prob, na.rm = TRUE))
    }
}

calcSDTKM <- function (Prob,risk=1){
    if(risk%in%c(1,2)){                             ## Risk 1: average risk over annual risks
        ## Risk 2: prob at least once during ny years
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, sd, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(sd(Prob, na.rm = TRUE))
    }else if(risk==3){                              ## Risk 3: Max risk of annual risks
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, max, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(max(Prob, na.rm = TRUE))
    }
}

calcSDTKM <- function (Prob,risk=1){
    if(risk%in%c(1,2)){                             ## Risk 1: average risk over annual risks
        ## Risk 2: prob at least once during ny years
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, sd, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(sd(Prob, na.rm = TRUE))
    }else if(risk==3){                              ## Risk 3: Max risk of annual risks
        if (class(Prob) == "matrix")
            return(apply(Prob, 1, max, na.rm = TRUE))
        if (class(Prob) == "numeric")
            return(max(Prob, na.rm = TRUE))
    }
}

calcSD <- function (Prob){
    if (class(Prob) == "matrix")
        return(apply(Prob, 2, sd, na.rm = TRUE))
    if (class(Prob) == "numeric")
        return(sd(Prob, na.rm = TRUE))
}

calcTimeTKM <- function (Prob){
    if (class(Prob) == "matrix"){
        tmp <- apply(Prob, 1, function(x) min(which(x <= 0.05),na.rm=TRUE))
        tmp[is.infinite(tmp)] <- length(Prob[1,])
        return(tmp)
    }
    if (class(Prob) == "numeric"){
        tmp <- min(which(Prob <= 0.05), na.rm = TRUE)
        tmp[is.infinite(tmp)] <- length(Prob)
        return(tmp)
    }
}


PBBlim <- function (MSEobj = NULL, risk = 1, Ref = 0.3, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- "B rel to Blim"
    if (Ref != 1) {
        PMobj@Caption <- paste0("P(B<Blim) [years ",
                                Yrs[1], "-", Yrs[2], "]")
    }else {
        PMobj@Caption <- paste0("P(B<Blim) [years ", Yrs[1],
                                "-", Yrs[2], "]")
    }
    PMobj@Ref <- Ref
    PMobj@Stat <- MSEobj@B_BMSY[, , Yrs[1]:Yrs[2]]
    PMobj@Prob <- calcProbTKM(PMobj@Stat < PMobj@Ref, MSEobj, risk = risk)
    PMobj@Mean <- calcMeanTKM(PMobj@Prob, risk=risk)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(PBBlim) <- "PM"

PBBmsy <- function (MSEobj = NULL, risk=1, Ref = 1, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- "B rel to Bmsy"
    if (Ref != 1) {
        PMobj@Caption <- paste0("P(B<Bmsy) [years ",
                                Yrs[1], "-", Yrs[2], "]")
    }else {
        PMobj@Caption <- paste0("P(B<Bmsy) [years ", Yrs[1],
                                "-", Yrs[2], "]")
    }
    PMobj@Ref <- Ref
    PMobj@Stat <- MSEobj@B_BMSY[, , Yrs[1]:Yrs[2]]
    PMobj@Prob <- calcProbTKM(PMobj@Stat < PMobj@Ref, MSEobj,risk=risk)
    PMobj@Mean <- calcMeanTKM(PMobj@Prob,risk=risk)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(PBBmsy) <- "PM"


YieldTKM <- function (MSEobj = NULL, Ref = 1, Yrs = NULL)
{
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- paste0("Yield relative to Reference Yield (Years ",
                         Yrs[1], "-", Yrs[2], ")")
    PMobj@Caption <- paste0("Mean Relative Yield [years ", Yrs[1],
                            "-", Yrs[2], "]")
    RefYd <- array(MSEobj@OM$RefY, dim = dim(MSEobj@C[, , Yrs[1]:Yrs[2]]))
    PMobj@Stat <- MSEobj@C[, , Yrs[1]:Yrs[2]]/RefYd
    PMobj@Ref <- Ref
    PMobj@Prob <- calcProb(PMobj@Stat, MSEobj)
    PMobj@Mean <- calcMean(PMobj@Prob)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(YieldTKM) <- "PM"


CVyield <- function (MSEobj = NULL, Ref = 1, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- paste0("Yield relative to Reference Yield (Years ",
                         Yrs[1], "-", Yrs[2], ")")
    PMobj@Caption <- paste0("CV of yield [years ", Yrs[1],
                            "-", Yrs[2], "]")
    CmeanSim <- apply(MSEobj@C[, , Yrs[1]:Yrs[2]],c(2,3), mean, na.rm=TRUE)  ## nsim, nMPs, proyears
    Cmean <- apply(CmeanSim,1, mean, na.rm=TRUE)  ## nMPs, proyears
    Csd <- apply(CmeanSim,1, sd, na.rm=TRUE)  ## nMPs, proyears
    PMobj@Stat <- CmeanSim
    PMobj@Ref <- Ref
    PMobj@Prob <- CmeanSim ##calcProb(PMobj@Stat, MSEobj)
    PMobj@Mean <- Csd/(Cmean+1.e-8)##calcMean(PMobj@Prob)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(CVyield) <- "PM"


SDyield <- function (MSEobj = NULL, Ref = 1, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- paste0("SD of yield relative to Reference Yield (Years ",
                         Yrs[1], "-", Yrs[2], ")")
    PMobj@Caption <- paste0("SD of yield [years ", Yrs[1],
                            "-", Yrs[2], "]")
    RefYd <- array(MSEobj@OM$RefY, dim = dim(MSEobj@C[, , Yrs[1]:Yrs[2]]))
    CmeanSim <- apply(MSEobj@C[, , Yrs[1]:Yrs[2]],c(2,3), mean, na.rm=TRUE)  ## nsim, nMPs, proyears
    ##    Cmean <- apply(CmeanSim,1, mean, na.rm=TRUE)  ## nMPs, proyears
    ##    Csd <- apply(CmeanSim,1, sd, na.rm=TRUE)  ## nMPs, proyears
    PMobj@Stat <- MSEobj@C[, , Yrs[1]:Yrs[2]]/RefYd
    PMobj@Ref <- Ref
    PMobj@Prob <- calcProb(PMobj@Stat, MSEobj)    ## mean over year or simulations?
    PMobj@Mean <- calcSD(PMobj@Prob)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(SDyield) <- "PM"


SDPBBlim <- function (MSEobj = NULL, risk = 1, Ref = 0.3, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- "SD B rel to Blim"
    if (Ref != 1) {
        PMobj@Caption <- paste0("SD P(B<Blim) [years ",
                                Yrs[1], "-", Yrs[2], "]")
    }else {
        PMobj@Caption <- paste0("SD P(B<Blim) [years ", Yrs[1],
                                "-", Yrs[2], "]")
    }
    PMobj@Ref <- Ref
    PMobj@Stat <- MSEobj@B_BMSY[, , Yrs[1]:Yrs[2]]
    PMobj@Prob <- calcProb(PMobj@Stat < PMobj@Ref, MSEobj)##, risk = risk)
    PMobj@Mean <- calcSD(PMobj@Prob) ##, risk=risk)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(SDPBBlim) <- "PM"


AAVYTKM <- function (MSEobj = NULL, Ref = 1, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- paste0("Average Annual Variability in Yield (Years ",
                         Yrs[1], "-", Yrs[2], ")")
    PMobj@Caption <- paste0("AAVY [years ",
                            Yrs[1], "-", Yrs[2], "]")
    y1 <- Yrs[1]:(Yrs[2] - 1)
    y2 <- (Yrs[1] + 1):Yrs[2]
    if (MSEobj@nMPs > 1) {
        AAVY <- apply(((((MSEobj@C[, , y1] - MSEobj@C[, , y2])/MSEobj@C[,, y2])^2)^0.5), c(1, 2),
                      median, na.rm=TRUE)
    }
    else {
        AAVY <- array(apply(((((MSEobj@C[, 1, y1] - MSEobj@C[, 1, y2])/
                               MSEobj@C[, 1, y2])^2)^0.5), c(1), median, na.rm=TRUE))
    }
    PMobj@Stat <- AAVY
    PMobj@Ref <- Ref
    PMobj@Prob <- calcProbMedian(PMobj@Stat, MSEobj)
    PMobj@Mean <- calcMedian(PMobj@Prob)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(AAVYTKM) <- "PM"



timetoPBBlim <- function (MSEobj = NULL, risk = 1, Ref = 0.3, Yrs = NULL){
    Yrs <- ChkYrs(Yrs, MSEobj)
    PMobj <- new("PMobj")
    PMobj@Name <- "Time ti B rel to Blim"
    if (Ref != 1) {
        PMobj@Caption <- paste0("Time to P(B<Blim) [years ",
                                Yrs[1], "-", Yrs[2], "]")
    }else {
        PMobj@Caption <- paste0("Time to P(B<Blim) [years ", Yrs[1],
                                "-", Yrs[2], "]")
    }
    PMobj@Ref <- Ref
    PMobj@Stat <- MSEobj@B_BMSY[, , Yrs[1]:Yrs[2]]
    PMobj@Prob <- calcProbTKM(PMobj@Stat < PMobj@Ref, MSEobj, risk = risk)
    PMobj@Mean <- calcTimeTKM(PMobj@Prob)
    PMobj@MPs <- MSEobj@MPs
    PMobj
}
class(timetoPBBlim) <- "PM"





ConvergeTKM <- function (MSEobj, PMs = list(Yield, P10, AAVY), maxMP = 15, thresh = 0.5,
                         ref.it = 20, inc.leg = FALSE, all.its = FALSE, nrow = NULL,
                         ncol = NULL)
{
    if (class(MSEobj) != "MSE")
        stop("MSEobj must be object of class 'MSE'", call. = FALSE)
    if (MSEobj@nMPs < 2)
        stop("Converge requires more than 1 MP in the MSE object",
             call. = FALSE)
    nPMs <- length(PMs)
    if (is.null(ncol))
        ncol <- floor(sqrt(nPMs))
    if (is.null(nrow))
        nrow <- ceiling(nPMs)/ncol
    if (ncol * nrow < nPMs)
        stop("ncol x nrow must be > length(PMs)")
    if (MSEobj@nMPs > maxMP) {
        nplots <- ceiling(MSEobj@nMPs/maxMP)
    }
    else {
        nplots <- 1
    }
    nsim <- MSEobj@nsim
    if (nsim - ref.it < 1) {
        ref.it.new <- nsim - 1
        message("nsim (", nsim, ") - ref.it (", ref.it, ") < 1. Setting ref.it to ",
                ref.it.new, "\n")
        ref.it <- ref.it.new
    }
    message("Checking if order of MPs is changing in last ",
            ref.it, " iterations\n")
    message("Checking average difference in PM over last ", ref.it,
            " iterations is > ", thresh, "\n")
    SwitchOrd <- vector(mode = "list", length = nPMs)
    NonCon <- vector(mode = "list", length = nPMs)
    PMName <- vector("character", length = nPMs)
    getPalette <- colorRampPalette(c("#00007F", "blue", "#007FFF",
                                     "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    st <- 1
    end <- min(maxMP, MSEobj@nMPs)
    it <- MP <- NULL
    for (nplot in 1:nplots) {
        message("Plotting MPs ", st, " - ", end)
        subMSE <- Sub(MSEobj, MPs = MSEobj@MPs[st:end])
        nMPs <- subMSE@nMPs
        values <- rep(c("solid", "dashed", "dotted"), nMPs)
        values <- values[1:nMPs]
        plist <- list()
        for (xx in 1:nPMs) {
            PMval <- PMs[[xx]](subMSE)
            PMName[xx] <- PMval@Name
            PMval@Prob[!is.finite(PMval@Prob)] <- 0
            cum_mean <- apply(PMval@Prob, 2, cumsum)/apply(PMval@Prob,
                                                           2, seq_along) * 100
            vals <- as.vector(cum_mean)
            mp <- rep(subMSE@MPs, each = subMSE@nsim)
            df <- data.frame(it = 1:subMSE@nsim, vals, MP = mp,
                             name = PMval@Name)
            if (!all.its)
                df <- subset(df, it %in% (nsim - ref.it + 1):nsim)
            p <- ggplot2::ggplot(df, ggplot2::aes(x = it, y = vals,
                                                  color = MP, linetype = MP)) + ggplot2::scale_linetype_manual(values = values) +
                          ggplot2::scale_color_manual(values = getPalette(nMPs)) +
                          ggplot2::geom_line() + ggplot2::theme_classic() +
                          ggplot2::labs(x = "# Simulations", y = PMval@Caption) +
                          ggplot2::ggtitle(PMval@Name) + ggplot2::geom_vline(xintercept = MSEobj@nsim -
                                                                                 ref.it + 1, color = "darkgray", linetype = "dashed") +
                          ggplot2::geom_vline(xintercept = MSEobj@nsim,
                                              color = "darkgray", linetype = "dashed") +
                          ggplot2::coord_cartesian(xlim = c(min(df$it),
                                                            max(df$it)))
            if (!inc.leg)
                p <- p + ggplot2::theme(legend.position = "none")
            ords <- apply(cum_mean[(MSEobj@nsim - ref.it + 1):MSEobj@nsim,
                                   ], 1, order, decreasing = FALSE)
            rownames(ords) <- subMSE@MPs
            mat <- matrix(subMSE@MPs[ords], nrow = subMSE@nMPs,
                          ncol = ref.it)
            tab <- table(unlist(apply(mat, 1, unique)))
            SwitchOrd[[xx]] <- append(SwitchOrd[[xx]], rownames(tab)[which(tab >
                                                                           1)])
            NonCon[[xx]] <- append(NonCon[[xx]], subMSE@MPs[apply(cum_mean,
                                                                  2, DLMtool:::Chk, MSEobj = MSEobj, thresh = thresh, ref.it)])
            noncoverg <- unique(c(SwitchOrd[[xx]], NonCon[[xx]]))
            if (length(noncoverg) > 0) {
                df2 <- subset(df, MP %in% noncoverg)
                if (dim(df2)[1] > 0) {
                    ##                  p <- p + ggrepel::geom_text_repel(data = subset(df2,
                    ##                    it == max(it)), ggplot2::aes(label = MP),
                    ##                    size = 4, segment.color = "grey", direction = "both",
                    ##                    show.legend = FALSE)
                }
            }
            plist[[xx]] <- p
        }
        if (inc.leg)
            join_plots(plist, nrow = nrow, ncol = ncol, position = "right")
        if (!inc.leg)
            gridExtra::grid.arrange(grobs = plist, nrow = nrow,
                                    ncol = ncol)
        st <- st + maxMP
        end <- min(end + maxMP, MSEobj@nMPs)
    }
    for (x in 1:nPMs) {
        SwitchOrds <- SwitchOrd[[x]]
        NonCons <- NonCon[[x]]
        if (length(SwitchOrds) > 0 | length(NonCons) > 0) {
            message("\n", PMName[x], "\n")
            if (length(SwitchOrds) > 0) {
                message("Order over last ", ref.it, " iterations is not consistent for:\n ",
                        paste(SwitchOrds, ""), " \n")
            }
            if (length(NonCons) > 0) {
                message("Mean difference over last ", ref.it,
                        " iterations is > ", thresh, " for:\n", paste(NonCons,
                                                                      ""), "\n")
            }
        }
    }
}




## plot  surplus production vs biomass over K
plotSPB <- function(splist, bklist){
    plot(bklist[[1]],splist[[1]], ty='n',lwd=2,
         ylab = "Surplus production", xlab = "B/K",
         xlim=c(0,1),ylim=range(unlist(splist), na.rm=TRUE))
    abline(v=0.5, col="grey70")
    for(i in 1:length(splist)){
        lines(bklist[[i]],splist[[i]],lwd=2, lty=1, ty="b", pch = NA)
        text(bklist[[i]],splist[[i]],label=i)
    }
}






## custom management procedures:
## ------------------------------------------------------------------------------------------
## custom hockey-stick MSY rule
MSYHSref <- function (x, Data, reps = 1){
    # runs with no obs error
##    B_BMSY <- Data@SpAbun[x]/Data@OM$SSBMSY[x] # spawning biomass/sb_msy
    # or:
##    B_BMSY <- Data@Abun[x]/Data@OM$BMSY[x] # total biomass/b_msy
    # or, updated true spawning abundance ignores obs model:
    B_BMSY <- Data@OM$Asp[x]/Data@OM$SSBMSY[x] # spawning biomass/sb_msy
    B_Btrigger <- B_BMSY * 2

    rec <- new("Rec")
    TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
    # linear decline in TAC
    if (B_Btrigger < 1) TAC <- TAC * B_Btrigger

    if (TAC <= 0) TAC <- tiny
    if(is.null(TAC) || !is.numeric(TAC)) TAC <- NA
    rec@TAC <- as.numeric(TAC)
    rec
}
class(MSYHSref) <- "MP"

## custom FMSYref (no TACfilter)
FMSYref2 <- function (x, Data, reps = 1){
    rec <- new("Rec")
    TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])),
                   0.01)
    if(is.null(TAC) || !is.numeric(TAC)) TAC <- NA
    rec@TAC <- as.numeric(TAC)
    rec
}
class(FMSYref2) <- "MP"

## Template function for 2 over 3 rule
make2o3rules <- function(reps = 1, uncertCap = TRUE,
                         lower = 0.8, upper = 1.2,
                         paBuffer = TRUE,
                         paInt = NA,
                         nhist = NA,
                         num = 2, denom = 3,
                         env = globalenv()){

    ## allowing for multiple generation of MPs
    argList <- list(uncertCap, lower, upper, paBuffer, paInt,  num, denom)
    argLengths <- sapply(argList, length)
    maxi <- max(argLengths)
    maxl  <- which(argLengths == maxi)
    if(maxi>1){
        if(max(argLengths[(1:7)[-maxl]]) > 1)
            stop("Specified arguments have different lengths, they should have the same length or length = 1.")
    }
    argListCor <- argList
    argListCor[(1:7)[-maxl]] <- lapply(argList[(1:7)[-maxl]], function(x) rep(unlist(x), maxi))
    template  <- expression(paste0(
        'structure(function(x, Data, reps = 1,
          uncertCap=',a,',
          lower=',b,',
          upper=',c,',
          paBuffer=',d,',
          paInt=',e,',
          num=',f,',
          denom=',g,'){
            dependencies <- "Data@Year, Data@Cat, Data@Ind, Data@Misc, Data@MPrec"
            ## limit data to use for assessment (with argument nhist)
            LH <- Data@LHYear   ## last historical year
            if(is.na(nhist)) nhist <- LH else nhist <- nhist
            curr.yr <- length(Data@Year)
            fst.yr <- max(LH-nhist+1, 1)
            yr.ind <-  fst.yr:curr.yr
            n.yr <- length(yr.ind)
            time <- yr.ind
            Catch <- Data@Cat[x,yr.ind]
       ##     Index <- Data@AddInd[x,1,yr.ind]
            Index <- Data@Ind[x,yr.ind]
            ## get indices
            I.num <- Index[(length(Index)-(num-1)):length(Index)]
            I.den <- Index[(length(Index)-(num+denom-1)):(length(Index)-num)]
            ## catch multiplication factor
            cmult <- mean(I.num, na.rm = TRUE)/mean(I.den, na.rm = TRUE)
            ## uncertainty cap
            if(uncertCap){
                cmult[cmult >= upper] <- upper
                cmult[cmult <= lower] <- lower
            }
            ## TAC based on cmult
            ## Ccurrent = most recent advice provided (in first year = Clast)
            Ccurr = Data@MPrec[x]
            TAC <- Ccurr * cmult
            ## PA buffer (0.2 red every int years)
            if(paBuffer){
                if(is.null(Data@Misc[[x]])){
                    pacounter <- paInt ## initalize
                }else{
                    pacounter <- Data@Misc[[x]]$pacounter
                }
                pacounter <- pacounter + 1
                if(pacounter > (paInt-1)){
                   TAC <- TAC * 0.8
                   pacounter <- 0
                }
            }else pacounter <- NULL
            ## return TAC
            Rec <- new("Rec")
            if(is.null(TAC) || !is.numeric(TAC) || (TAC < 0)) taci <- NA else taci <- TAC
            Rec@TAC <- as.numeric(taci)
            Rec@Misc <- list(pacounter = pacounter, TAC = taci)
            return(Rec)
        },
        class="MP")'))
    nami <- rep(NA,maxi)
    for(I in 1:maxi){
        ## create MPs as functions
        subList <- lapply(argListCor, "[[", I)
        names(subList) <- letters[1:7]
        templati <- eval(parse(text=paste(parse(text = eval(template, subList)),collapse=" ")))

        c0 <- paste0(argListCor[[6]][I],"over",argListCor[[7]][I])
        ## save names of MPs
        if(argListCor[[1]][I] == FALSE){
            c1 <- ""
        }else{
            c1 <- paste0("_uC")
        }
        if(argListCor[[2]][I] == 0.8){
            c2 <- ""
        }else{
            c2 <- paste0("_lo",argListCor[[2]][I])
        }
        if(argListCor[[3]][I] == 1.2){
            c3 <- ""
        }else{
            c3 <- paste0("_up",argListCor[[3]][I])
        }
        if(argListCor[[4]][I] == FALSE){
            c4 <- ""
        }else{
            c4 <- paste0("_paB")
        }
        if(is.na(argListCor[[5]][I])){
            c5 <- ""
        }else{
            c5 <- paste0("_paI",round(argListCor[[5]][I],2))
        }

        ## put everythin together
        nami[I] <- paste0("ind",c0,c1,c2,c3,c4,c5)
        assign(value=templati, x=nami[I], envir=env)
    }
    ## allow for assigning names
    invisible(nami)
}




### NEW 2019
plotdYdR <- function (MSEobj, PMlist = c("PBBlim","YieldTKM"), risk=1, pchAAVY = TRUE,
                      pchs=NA,argroup=NA,arcols=NA,aralphas=aralphas,
                      connect=NA, xaxt = "s", yaxt="s", relref=FALSE,
                        Title = NULL,MPs=NA,
                        Yrs = NULL, Refs = NULL,
                        xlim=NULL, ylim=NULL,
                      xcap=NULL,ycap=NULL){

    nPMs <- length(PMlist)
    runPM <- vector("list", length(PMlist))
    for (X in 1:length(PMlist)) {
        ref <- Refs[[PMlist[X]]]
        yrs <- Yrs##[[PMlist[X]]]
        if (is.null(ref)) {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs = yrs))
                }
            }
        }else {
            if (is.null(yrs)) {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref))
                }
            }else {
                if(PMlist[X]%in%c("PBBlim","PBBmsy")){
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref, risk = risk,
                                            Yrs = yrs))
                }else{
                    runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref = ref,
                                            Yrs = yrs))
                }
            }
        }
    }

    nplots <- nPMs/2
    n.col <- ceiling(sqrt(nplots))
    n.row <- ceiling(nplots/n.col)
    ##    m <- matrix(1:(n.col * n.row), ncol = n.col, nrow = n.row,
    ##        byrow = FALSE)
    xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
    plots <- listout <- list()
    xInd <- seq(1, by = 2, length.out = nplots)
    yInd <- xInd + 1
    ##
    yPM <- PMlist[yInd[1]]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    if(is.null(ycap)) ycap <- runPM[[match(yPM, PMlist)]]@Caption else ycap <- ycap
    yname <- runPM[[match(yPM, PMlist)]]@Name
    xPM <- PMlist[xInd[1]]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    if(is.null(xcap)) xcap <- runPM[[match(xPM, PMlist)]]@Caption else xcap <- xcap
    xname <- runPM[[match(xPM, PMlist)]]@Name

    if(pchAAVY){
        runPMz <- eval(call("AAVYTKM", MSEobj, Yrs = yrs))
        zvals <- runPMz@Mean
        ix1 <- which(zvals < 0.1)
        ix2 <- which(zvals < 0.3 & zvals >= 0.1)
        ix3 <- which(zvals >= 0.3)
        arcex <- zvals
        arcex[ix1] <- 1.5
        arcex[ix2] <- 2.25
        arcex[ix3] <- 3
##        arcex <- rescale(zvals, to = c(1.5,3)) ##scale(zvals, min(zvals), 1/2)+1.5
##        print(range(arcex))
    }else{
        arcex <- rep(1.5,nrow(dfnew))
    }

    if(yPM == "AAVY"){
        ycap <- paste0("P(AAVY<0.2) [years ",yrs[1],"-",yrs[2],"]")
    }
    ##    if(yPM == "CVyield"){
    ##        ycap <- paste0(" [years ",yrs[1],"-",yrs[2],"]")
    ##    }
    df <- data.frame(x = xvals, y = yvals, label = MSEobj@MPs,
                     fontface = "plain", xPM = xPM, yPM = yPM)

    ref <- df[which(df$label == "FMSYref2"),]
    dfnew <- df##[-which(df$label %in% c("NFref","FMSYref2")),]
    dfnew$dx <- dfnew$x * 100 ##- ref$x
    dfnew$dy <- dfnew$y
    df$fontface <- as.character(df$fontface)
    df$fontface <- factor(df$fontface)

    if(is.null(xlim)) xlim <- range(c(0,dfnew$dx))
    if(is.null(ylim)) ylim <- range(c(1,dfnew$dy))

    if(is.na(pchs[1])) pchs <- c(3,4, 17,15,16)
    if(is.na(argroup[1])) argroup <- c(1,2,rep(3,4),4,5)
    if(is.na(arcols[1])) arcols <- c(1,1,paste0(rep("grey",4),round(seq(75,10,length.out = 4))),1,1)
    if(is.na(aralphas[1])) aralphas <- c(rep(1,3), seq(0.3,1,length.out = 3),1)

    if(is.na(connect[1])) connect <- 1

    if(relref){
        ind <- which("Ref" == df$label)
        if(length(ind) != 1) stop("Problem with relref!")
        refi <- dfnew[ind,]
        dfnew <- dfnew[-ind,]
        if(refi$dx == 0) refi$dx <- 1e-8
        dfnew$dx <- dfnew$dx / refi$dx
        dfnew$dy <- dfnew$dy / refi$dy
        argroup <- argroup[-ind]
        arcols <- arcols[-ind]
        arcex <- arcex[-ind]
    }

    ##    print(MSEobj@MPs)
    plot(dfnew$dx, dfnew$dy, ty='n',main=Title,
         xlab = expression(Delta~Risk), xaxt="n",yaxt="n",
         ylab = expression(Delta~Yield), xlim=xlim,
         ylim=ylim)
    if(xaxt == "s" && relref){
        axis(1,at=0:100, labels=0:100,cex.axis=1.5)
    }else if(xaxt == "s"){
        axis(1,cex.axis=1.5)
    }
    if(yaxt == "s") axis(2,at=seq(0,2,0.5), labels=seq(0,2,0.5),cex.axis=1.5)
    ## segments for reference rule
    if(!relref) abline(h=dfnew$dy[1], lty=2, col="grey50", lwd=1.2)
    if(!relref) abline(v=dfnew$dx[1], lty=2, col="grey50", lwd=1.2)
    if(relref) abline(h=1, lty=2, col="grey60", lwd=1.5)
    if(relref) abline(v=1, lty=2, col="grey60", lwd=1.5)
    ## segments(x0=-1,y0=dfnew$dy[1],x1=dfnew$dx[1],y1=dfnew$dy[1],lty=2,
    ##          col="grey40")
    ## segments(x0=dfnew$dx[1],y0=-1,x1=dfnew$dx[1],y1=dfnew$dy[1],lty=2,
    ##          col="grey40")
##    abline(v=0.05,lty=2,col="grey80")
    ## abline(h=1,col="grey70")
    ## abline(v=0,col="grey70")
    points(dfnew$dx, dfnew$dy, lwd=1.5,
           pch=pchs[argroup], col=arcols,cex=arcex)
    points(dfnew$dx[connect[[1]]], dfnew$dy[connect[[1]]], ty='b', pch=NA,
           col="grey40", lwd=2, lty=3)
    if(class(connect) == "list" && length(connect) > 1)
        points(dfnew$dx[connect[[2]]], dfnew$dy[connect[[2]]], ty='b', pch=NA,
               col="grey40", lwd=2, lty=3)
    if(class(connect) == "list" && length(connect) > 2)
        points(dfnew$dx[connect[[3]]], dfnew$dy[connect[[3]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    if(class(connect) == "list" && length(connect) > 3)
        points(dfnew$dx[connect[[4]]], dfnew$dy[connect[[4]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    if(class(connect) == "list" && length(connect) > 4)
        points(dfnew$dx[connect[[5]]], dfnew$dy[connect[[5]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    if(class(connect) == "list" && length(connect) > 5)
        points(dfnew$dx[connect[[6]]], dfnew$dy[connect[[6]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    if(class(connect) == "list" && length(connect) > 6)
        points(dfnew$dx[connect[[7]]], dfnew$dy[connect[[7]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    if(class(connect) == "list" && length(connect) > 7)
        points(dfnew$dx[connect[[8]]], dfnew$dy[connect[[8]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    if(class(connect) == "list" && length(connect) > 8)
        points(dfnew$dx[connect[[9]]], dfnew$dy[connect[[9]]], ty='b', pch=NA,
               col="grey60", lwd=1.5, lty=2)
    ## points(xvals[ind2], yvals[ind2], ty='b', pch=NA, col="grey60", lwd=1.4)
    ## points(xvals[ind3], yvals[ind3], ty='b', pch=NA, col="grey60", lwd=1.4,
    ##        lty=3)
    box()

}



## non covergence in MSElist[[1]]@Misc$spict
## previously: MSElist[[2]]@Misc$Data[[12]]@Misc
nonconv <- function(MSEobj){
    ## number of non-conv assessments
    ## non conv assess
    nca <- lapply(MSEobj@Misc$Data, function(x) unlist(lapply(x@Misc, function(y){
        if(length(y$spict) > 0){
            length(which(as.numeric(y$spict) != 2))
        }else 0
    })))

    ## nca <- lapply(MSEobj@Misc$spict, function(y) unlist(lapply(y, function(x){
    ##     if(length(x) > 0){
    ##         length(which(as.numeric(x) != 2))
    ##     }else 0
    ## })))

    ## any non-converged assessment (for exclusion)
    ## non conv sims

    ncs <- lapply(MSEobj@Misc$Data, function(x) unlist(lapply(x@Misc, function(y){
        if(length(y$spict) > 0){
            any(as.numeric(y$spict) != 2)  ## x[5:19]  ## to not use first 4 years
        }else FALSE
    })))

    ## ncs <- lapply(MSEobj@Misc$spict, function(y) unlist(lapply(y, function(x){
    ##     if(length(x) > 0){
    ##         any(as.numeric(x) != 2)  ## x[5:19]  ## to not use first 4 years
    ##     }else FALSE
    ## })))
    ## return
    return(list(nca,ncs))
}


## calculating production curve for OM

calcProdCurve <- function(MSEobj, plot = TRUE){##, varname="K", roundi=2){

    dims <- dim(MSEobj@TSdata$Catch)

    splist <- list()
    for(i in 1:dims[1]){
        sp <- rep(NA,(dims[2]-1))
        for(j in 2:dims[2]){
            sp[j-1] <- MSEobj@TSdata$B[i,j] - MSEobj@TSdata$B[i,(j-1)] + MSEobj@TSdata$C[i,(j-1)]
        }
        splist[[i]] <- sp
    }

    bklist <- list()
    for(i in 1:dims[1]){
        bklist[[i]] <- MSEobj@TSdata$B[i,-1] / MSEobj@Ref$B0[i]
    }

    spmsylist <- list()
    for(i in 1:dims[1]){
        spmsylist[[i]] <- splist[[i]] / MSEobj@Ref$MSY[i]
    }

    ## maxima
    bmsyk <- rep(NA,dims[1])
    for(i in 1:dims[1]){
        bmsyk[i] <- bklist[[i]][which.max(spmsylist[[i]])]
    }
    bmsykstats <- c(quantile(bmsyk,0.025),median(bmsyk),quantile(bmsyk,0.975))

    ## vari <- unlist(MSEobj@OM[varname])
    ## idx <- sort(unique(round(vari,roundi)))
    ## tabi <- cbind(idx,brewer.pal(length(idx),"Blues"))

    ## browser()

    ## idx <- which(bmsyk < 0.43)

    ## colMeans(MSEobj@OM[idx,])

    ## colMeans(MSEobj@OM[-idx,])

    res <- list(spmsylist,bklist,bmsykstats)


    if(plot){
        ylim <- range(unlist(spmsylist))
        plot(bklist[[i]],spmsylist[[1]], ty='n',
             xlim=c(0,1),ylim=ylim,
             ylab="SP/MSY", xlab="B/K")
        for(i in 1:dims[1]){
            lines(bklist[[i]],spmsylist[[i]], col="grey60")##tabi[match(round(vari,2)[i], tabi),2])
        }
        polygon(x=c(bmsykstats[1],bmsykstats[3],bmsykstats[3],bmsykstats[1]),
                y= c(-2,-2,2*c(ylim[2],ylim[2])),
                border=NA, col=rgb(t(col2rgb("blue")/255),alpha=0.4))
        abline(v=bmsykstats[2],col=4,lwd=2)
    }

    return(res)

}


## select only converged simulations
getConvSims <- function(MSElist, indexMPs=NULL){
    idx <- if(is.null(indexMPs)) 1:MSElist[[1]]@nMPs else indexMPs
    nscenarios <- length(MSElist)
    MSElistConv <- list()
    for(i in 1:nscenarios){
        nonconvi <- nonconv(MSElist[[i]])
        tmp <- do.call(cbind,nonconvi[[2]])[,idx]
        if(length(idx)>1){
            convSims <- which(apply(tmp,1, function(x) all(x==FALSE)))
        }else{
            convSims <- which(tmp == FALSE)
        }
        MSElistConv[[i]] <- Sub(MSElist[[i]], sims=convSims)
    }
    return(MSElistConv)
}



## specific plots

## haddock plot
plotTradeHaddock <- function(MSElist, evalyears,
                             xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                             ptit1="",ptit2="",ptit3="",ptit4=""){
    if(is.null(xlimUpper)) xlimUpper <- rep(50,4)
    if(is.null(xlimLower)) xlimLower <- rep(0,4)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,4)
    if(is.null(ylimLower)) ylimLower <- rep(0,4)
    ly <- layout(matrix(c(1,2,5,3,4,5),byrow = TRUE,ncol=3,nrow=2),
                 widths=c(1,1,0.4,1,1,0.4), heights=c(1,1))
    opar <- par(mar=c(3,2,3,2),oma=c(5,4,1,0))
    legends=c("ref","2/3","2/3uC","2/3uC_PA3",
          "1/2","1/2uC","1/2uC_PA3",
          "MSY50","MSY45","MSY35","MSY25",
          "MSY50-S1","MSY45-S1","MSY35-S1","MSY25-S1",
          "MSY50-S05","MSY45-S05","MSY35-S05","MSY25-S05",
          "MSY-C45","MSY-C35","MSY-C25",
          "MSY-C45-S05","MSY-C35-S05","MSY-C25-S05")
    pchs <- c(1,3,4,16,15,22,17,24)
    argroup <- c(1,2,2,2,3,3,3,
                 rep(4,4),
                 rep(5,4),
                 rep(6,4),
                 rep(7,3),
                 rep(8,3))
    arcols <- c(1,
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4)))[-1],
                paste0(rep("grey",4),round(seq(75,10,length.out = 4)))[-1])
    aralphas <- c(1,seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4)[-1],
                  seq(0.3,1,length.out = 4)[-1])
    connect <- list(c(2:4),c(5:7),c(8:11),c(12:15),c(16:19),c(8,20:22),c(16,23:25))
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[1],ylimUpper[1]),
             xlim=c(xlimLower[1],xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit1,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[2]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[2],ylimUpper[2]),
             xlim=c(xlimLower[2],xlimUpper[2]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit2,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[3]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[3],ylimUpper[3]),
             xlim=c(xlimLower[3],xlimUpper[3]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit3,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[4]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[4],ylimUpper[4]),
             xlim=c(xlimLower[4],xlimUpper[4]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit4,
          font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 0.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.2, outer=TRUE, cex=1.2)
    plot.new()
    legend(0,1, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.6)
    legend(0,0.25, title = "MIAVY", title.adj=0,
           legend=c("<0.25","0.25-0.75"," >=0.75"), cex=1.6,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}

plotTradeHaddock_spict <- function(MSElist, evalyears,
                             xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                             ptit1="",ptit2="",ptit3="",ptit4=""){
    if(is.null(xlimUpper)) xlimUpper <- rep(50,4)
    if(is.null(xlimLower)) xlimLower <- rep(0,4)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,4)
    if(is.null(ylimLower)) ylimLower <- rep(0,4)
    ly <- layout(matrix(c(1,2,5,3,4,5),byrow = TRUE,ncol=3,nrow=2),
                 widths=c(1,1,0.4,1,1,0.4), heights=c(1,1))
    opar <- par(mar=c(3,2,3,2),oma=c(5,4,1,0))
    legends=c(
##        "ref","2/3","2/3uC","2/3uC_PA3",
##          "1/2","1/2uC","1/2uC_PA3",
          "MSY50","MSY45","MSY35","MSY25",
##          "MSY50-S1","MSY45-S1","MSY35-S1","MSY25-S1",
          "MSY50-S05","MSY45-S05","MSY35-S05","MSY25-S05",
          "MSY-C45","MSY-C35","MSY-C25",
          "MSY-C45-S05","MSY-C35-S05","MSY-C25-S05")
    pchs <- c(16,15,21,22)
    argroup <- c(##1,2,2,2,3,3,3,
                 rep(1,4),
                 rep(2,4),
                 rep(3,3),
                 rep(4,3))
    arcols <- c(##1,
##                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
##                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4)))[-1],
                paste0(rep("grey",4),round(seq(75,10,length.out = 4)))[-1])
    aralphas <- c(##1,seq(0.3,1,length.out = 3),
##                  seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4)[-1],
                  seq(0.3,1,length.out = 4)[-1])
    connect <- list(##c(2:4),c(5:7),
        c(1:4),c(5:8),c(1,9:11),c(5,12:14))
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[1],ylimUpper[1]),
             xlim=c(xlimLower[1],xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit1,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[2]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[2],ylimUpper[2]),
             xlim=c(xlimLower[2],xlimUpper[2]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit2,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[3]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[3],ylimUpper[3]),
             xlim=c(xlimLower[3],xlimUpper[3]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit3,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[4]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[4],ylimUpper[4]),
             xlim=c(xlimLower[4],xlimUpper[4]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(ptit4,
          font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 0.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.2, outer=TRUE, cex=1.2)
    plot.new()
    legend(0,1, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.6)
    legend(0,0.25, title = "MIAVY", title.adj=0,
           legend=c("<0.25","0.25-0.75"," >=0.75"), cex=1.6,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}




## plot SSB/SSBMSY
plotSSBSSBMSY <- function(MSElist, savedir, species, scenarios, addLab = "", ylim=NULL){
    if(is.null(ylim)) ylim = c(0,2)
    for(i in 1:length(MSElist)){
        MSEobj <- MSElist[[i]]
        pdf(paste0(savedir,paste0("/SSBBMSYPlot_",species,"_scen",scenarios[i],"_",addLab,".pdf")),
            width=12,height=10)
        opar <- par(mfrow=c(5,4), mar=c(3,2,3,2),oma=c(5,4,1,0))
        for(j in 1:MSEobj@nMPs){
            if(j %in% c(2)){plot.new()
                plot.new()
                plot.new()}
            if(j %in% c(5,8)) rep(plot.new(),1)
            plotSSB(MSEobj, j, quant = TRUE, ylim=ylim)
            mtext(MSEobj@MPs[j])
        }
        mtext("Years",side=1,line=1,outer=TRUE)
        mtext("SSB/SSBmsy",side=2,line=1,outer=TRUE)
        mtext(paste0(species, "  Scenario = ",scenarios[i],"  ",addLab),side=3,line=-1,outer=TRUE)
        par(opar)
        dev.off()
    }
}


## FFMSY
plotFFMSY <- function(MSElist, savedir, species, scenarios, addLab = "", ylim = NULL){
    if(is.null(ylim)) ylim = c(0,4)
    for(i in 1:length(MSElist)){
        MSEobj <- MSElist[[i]]
        pdf(paste0(savedir,paste0("/FFMSYPlot_",species,"_scen",scenarios[i],"_",addLab,".pdf")),
            width=12,height=10)
        opar <- par(mfrow=c(5,4), mar=c(3,2,3,2),oma=c(5,4,1,0))
        for(j in 1:MSEobj@nMPs){
            if(j %in% c(2)){plot.new()
                plot.new()
                plot.new()}
            if(j %in% c(5,8)) rep(plot.new(),1)
            plotFM(MSEobj, j, quant = TRUE, ylim=ylim)
            mtext(MSEobj@MPs[j])
        }
        mtext("Years",side=1,line=1,outer=TRUE)
        mtext("F/Fmsy",side=2,line=1,outer=TRUE)
        mtext(paste0(species,"  Scenario = ",scenarios[i],"  ",addLab),side=3,line=-1,outer=TRUE)
        par(opar)
        dev.off()
    }
}





## P(B<Blim) (Risk 1)
plotPBBlim <- function(MSElist, savedir, species, scenarios, addLab = "", ylim = NULL){
    if(is.null(ylim)) ylim = c(-0.1,1)
    for(i in 1:length(MSElist)){
        MSEobj <- MSElist[[i]]
        tmp <- PBBlim(MSEobj, risk = 1)@Prob
        pdf(paste0(savedir,paste0("/PBBlimPlot_",species,"_scen",scenarios[i],"_",addLab,".pdf")),
            width=12,height=10)
        opar <- par(mfrow=c(5,4), mar=c(3,2,3,2),oma=c(5,4,1,0))
        for(j in 1:MSEobj@nMPs){
            if(j %in% c(2)){plot.new()
                plot.new()
                plot.new()}
            if(j %in% c(5,8)) rep(plot.new(),1)
            plot(tmp[j,], ty='n',
                 ylim=ylim,
                 ylab="",xlab="")
            abline(h=seq(0.2,1,0.2),lwd=1,lty=2,col="grey80")
            abline(h=0.05, lwd=1.5, col="grey60")
            lines(tmp[j,], lwd=1.5)
            mtext(MSEobj@MPs[j])
        }
        mtext("Years",side=1,line=1,outer=TRUE)
        mtext("P(B<Blim)",side=2,line=1,outer=TRUE)
        mtext(paste0(species,"  Scenario = ",scenarios[i],"  ",addLab),side=3,line=-1,outer=TRUE)
        par(opar)
        dev.off()
    }
}


plotYield <- function(MSElist, savedir, species, scenarios, addLab = "", ylim = NULL){
    if(is.null(ylim)) ylim = c(0,2)
    for(i in 1:length(MSElist)){
        MSEobj <- MSElist[[i]]
        pdf(paste0(savedir,paste0("/relYieldPlot_",species,"_scen",scenarios[i],"_",addLab,".pdf")),
            width=12,height=10)
        opar <- par(mfrow=c(5,4), mar=c(3,2,3,2),oma=c(5,4,1,0))
        for(j in 1:MSEobj@nMPs){
            if(j %in% c(2)){plot.new()
                plot.new()
                plot.new()}
            if(j %in% c(5,8)) rep(plot.new(),1)
            plotY(MSEobj, j, quant = TRUE, ylim=ylim)
            mtext(MSEobj@MPs[j])
        }
        mtext("Years",side=1,line=1,outer=TRUE)
        mtext("Rel. Catch",side=2,line=1,outer=TRUE)
        mtext(paste0(species,"  Scenario = ",scenarios[i],"  ",addLab),side=3,line=-1,outer=TRUE)
        par(opar)
        dev.off()
    }
}






## ling anchovy plot
plotTradeLA <- function(MSElist, evalyears, xlimUpper=NULL, plottitle1="", plottitle2="",ylims=NULL){
    if(is.null(xlimUpper)) xlimUpper <- rep(50,4)
    if(is.null(ylims)) ylims <- c(0,1.05)
    ly <- layout(matrix(c(1,2,3),byrow = TRUE,ncol=3,nrow=1),
                 widths=c(1,1,0.7), heights=c(1,1))
    opar <- par(mar=c(3,2,3,2),oma=c(5,4,1,0))
legends=c("ref","2/3","2/3uC","2/3uC_PA3",
          "1/2","1/2uC","1/2uC_PA3",
          "MSY50","MSY45","MSY35","MSY25",
          "MSY50-S1","MSY45-S1","MSY35-S1","MSY25-S1",
          "MSY50-S05","MSY45-S05","MSY35-S05","MSY25-S05",
          "MSY-C45-S1","MSY-C35-S1","MSY-C25-S1",
          "MSY-C45-S05","MSY-C35-S05","MSY-C25-S05")
    pchs <- c(8,3,4,17,15,16,2,1)
    argroup <- c(1,2,2,2,3,3,3,rep(4,4),rep(5,4),rep(6,4),rep(7,3),rep(8,3))
    arcols <- c(1,
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4)))[2:4],
                paste0(rep("grey",4),round(seq(75,10,length.out = 4)))[2:4])
    aralphas <- c(1,seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4)[2:4],
                  seq(0.3,1,length.out = 4)[2:4])
    connect <- list(c(2:4),c(5:7),c(8:11),c(12:15),c(16:19),c(8,20:22),c(16,23:25))
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = ylims, xlim=c(0,xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext(plottitle1,
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[2]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = ylims, xlim=c(0,xlimUpper[2]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
##    abline(v=12.380952,lty=3)
    mtext(plottitle2,
          font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 0.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.2, outer=TRUE, cex=1.2)
    plot.new()
    legend(-0.03,1.05, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.6)
    legend(0.6,0.66, title = "MIAVY", title.adj=0,
           legend=c("<0.25","0.25-0.75"," >=0.75"), cex=1.6,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}



## ling anchovy plot
plotTrade1 <- function(MSElist, evalyears, xlimUpper=NULL){
    if(is.null(xlimUpper)) xlimUpper <- rep(50,4)
    ly <- layout(matrix(c(1,2),byrow = TRUE,ncol=2,nrow=1),
                 widths=c(1,0.4), heights=c(1))
    opar <- par(mar=c(3,2,3,2),oma=c(5,4,1,0))
    legends=c("ref","2/3","2/3uC","2/3uC_PA3",
              "1/2","1/2uC","1/2uC_PA3",
              "MSY50","MSY45","MSY35","MSY25",
              "MSY50-S05","MSY45-S05","MSY35-S05","MSY25-S05",
              "MSY50-S01","MSY45-S01","MSY35-S01","MSY25-S01")
    pchs <- c(1,3,4,17,16,15)
    argroup <- c(1,2,2,2,3,3,3,rep(4,4),rep(5,4),rep(6,4))
    arcols <- c(1,
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))))
    aralphas <- c(1,seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4),
                  seq(0.3,1,length.out = 4))
    connect <- list(c(2:4),c(5:7),c(8:11),c(12:15),c(16:19))
    ylims <- c(0,1.05)
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = ylims, xlim=c(0,xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext("S9",##bquote(sigma["R"] %in% "[0.5;0.7],"~~sigma["C"]~"="~sigma["I"]~"= 0.35,"~~"n = 20,"~~"D ~ 0.5"~B["MSY"]),
          font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 0.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.2, outer=TRUE, cex=1.2)
    plot.new()
    legend(0,1, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1)
    legend(0,0.25, title = "MIAVY", title.adj=0,
           legend=c("<0.25","0.25-0.75"," >=0.75"), cex=1,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))
}



## ling anchovy plot
plotTrade1Schaefer <- function(MSElist, evalyears, xlimUpper=NULL){
    if(is.null(xlimUpper)) xlimUpper <- rep(50,4)
    ly <- layout(matrix(c(1,2),byrow = TRUE,ncol=2,nrow=1),
                 widths=c(1,0.4), heights=c(1))
    opar <- par(mar=c(3,2,3,2),oma=c(5,4,1,0))
    legends=c("ref","2/3","2/3uC","2/3uC_PA3",
              "1/2","1/2uC","1/2uC_PA3",
              "MSY-S","MSY-45-S","MSY-35-S","MSY-25-S")
    pchs <- c(1,3,4,15)
    argroup <- c(1,2,2,2,3,3,3,rep(4,4))
    arcols <- c(1,
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",3),round(seq(75,10,length.out = 3))),
                paste0(rep("grey",4),round(seq(75,10,length.out = 4))))
    aralphas <- c(1,seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 3),
                  seq(0.3,1,length.out = 4))
    connect <- list(c(2:4),c(5:7),c(8:11))
    ylims <- c(0,1.05)
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = ylims, xlim=c(0,xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    abline(v=5,lty=2)
    mtext("S9", ##bquote(sigma["R"] %in% "[0.8;0.9],"~~sigma["C"]~"="~sigma["I"] %in% "[0.9;1],"~~"n = 20,"~~"D ~ 0.5"~B["MSY"]),
          font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 0.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.2, outer=TRUE, cex=1.2)
    plot.new()
    legend(0,1, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1)
    legend(0,0.25, title = "MIAVY", title.adj=0,
           legend=c("<0.25","0.25-0.75"," >=0.75"), cex=1,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))
}






## testing
## runMSE <- function (OM = DLMtool::testOM, MPs = c("AvC", "DCAC", "FMSYref",
##     "curE", "matlenlim", "MRreal"), CheckMPs = FALSE, timelimit = 1,
##     Hist = FALSE, ntrials = 100, fracD = 0.05, CalcBlow = TRUE,
##     HZN = 2, Bfrac = 0.5, AnnualMSY = TRUE, silent = FALSE, PPD = TRUE,
##     parallel = FALSE, save_name = NULL, checks = FALSE, control = NULL)
## {
##     if (class(OM) != "OM")
##         stop("OM is not class 'OM'", call. = FALSE)
##     rm(list = ls(DLMenv), envir = DLMenv)
##     tt <- suppressWarnings(try(lsf.str(envir = globalenv()),
##         silent = TRUE))
##     if (class(tt) != "try-error") {
##         gl.funs <- as.vector(tt)
##         pkg.funs <- as.vector(ls.str("package:DLMtool"))
##         if ("package:MSEtool" %in% search())
##             pkg.funs <- c(pkg.funs, as.vector(ls.str("package:MSEtool")))
##         if (length(gl.funs) > 0) {
##             gl.clss <- unlist(lapply(lapply(gl.funs, get), class))
##             gl.MP <- gl.funs[gl.clss %in% "MP"]
##             if (length(gl.MP) > 0) {
##                 inc.gl <- gl.MP[gl.MP %in% MPs]
##                 if (length(inc.gl) > 0) {
##                   dup.MPs <- inc.gl[inc.gl %in% pkg.funs]
##                   if (length(dup.MPs) > 0) {
##                     stop("Custom MP names already in DLMtool: ",
##                       paste0(dup.MPs, " "), "\nRename Custom MPs")
##                   }
##                 }
##             }
##         }
##     }
##     if (!all(is.na(MPs))) {
##         for (mm in MPs) {
##             chkMP <- try(get(mm), silent = TRUE)
##             if (class(chkMP) != "MP")
##                 stop(mm, " is not a valid MP", call. = FALSE)
##         }
##     }
##     if (parallel) {
##         if (OM@nsim < 48)
##             stop("nsim must be >=48 for parallel processing",
##                 call. = FALSE)
##         if (!snowfall::sfIsRunning()) {
##             message("Parallel processing hasn't been initialized. Calling 'setup()' now")
##             setup()
##         }
##         if (all(is.na(MPs)))
##             MPs <- avail("MP")
##         cMPs <- MPs[!MPs %in% pkg.funs]
##         globalMP <- NULL
##         extra_package <- NULL
##         for (mm in seq_along(cMPs)) {
##             nmspace <- utils::find(cMPs[mm])
##             if (nmspace == ".GlobalEnv") {
##                 globalMP <- c(globalMP, cMPs[mm])
##             }
##             else {
##                 extra_package <- c(extra_package, strsplit(nmspace,
##                   ":")[[1]][2])
##             }
##             extra_package <- unique(extra_package)
##         }
##         if (!is.null(globalMP)) {
##             message("Exporting custom MPs in global environment")
##             snowfall::sfExport(list = globalMP)
##         }
##         if (!is.null(extra_package)) {
##             message("Exporting additional packages with MPs")
##             for (pk in extra_package) sfLibrary(pk, character.only = TRUE,
##                 verbose = FALSE)
##         }
##         ncpu <- snowfall::sfCpus()
##         nits <- ceiling(OM@nsim/5)
##         itsim <- rep(5, nits)
##         if (nits < ncpu) {
##             if (nits < 4) {
##                 nits <- 4
##                 itsim <- rep(ceiling(OM@nsim/4), 4)
##             }
##             else {
##                 nits <- ncpu
##                 itsim <- rep(ceiling(OM@nsim/ncpu), ncpu)
##             }
##         }
##         cnt <- 1
##         while (sum(itsim) != OM@nsim | any(itsim < 2)) {
##             diff <- OM@nsim - sum(itsim)
##             if (diff > 0) {
##                 itsim[cnt] <- itsim[cnt] + 1
##             }
##             if (diff < 0) {
##                 itsim[cnt] <- itsim[cnt] - 1
##             }
##             cnt <- cnt + 1
##             if (cnt > length(itsim))
##                 cnt <- 1
##         }
##         if (!silent & !Hist)
##             message("Running MSE in parallel on ", ncpu, " processors")
##         if (!silent & Hist)
##             message("Running historical simulations in parallel on ",
##                 ncpu, " processors")
##         temp <- snowfall::sfClusterApplyLB(1:nits, DLMtool:::run_parallel,
##             itsim = itsim, OM = OM, MPs = MPs, CheckMPs = CheckMPs,
##             timelimit = timelimit, Hist = Hist, ntrials = ntrials,
##             fracD = fracD, CalcBlow = CalcBlow, HZN = HZN, Bfrac = Bfrac,
##             AnnualMSY = AnnualMSY, silent = TRUE, PPD = PPD,
##             control = control, parallel = parallel)
##         if (!is.null(save_name) && is.character(save_name))
##             saveRDS(temp, paste0(save_name, ".rdata"))
##         MSE1 <- joinMSE(temp)
##         if (class(MSE1) == "MSE") {
##             if (!silent)
##                 message("MSE completed")
##         }
##         else if (class(MSE1) == "Hist") {
##             if (!silent)
##                 message("Historical simulations completed")
##         }
##         else {
##             warning("MSE completed but could not join MSE objects. Re-run with `save_name ='MyName'` to debug")
##         }
##     }
##     if (!parallel) {
##         if (OM@nsim > 48 & !silent & !Hist)
##             message("Suggest using 'parallel = TRUE' for large number of simulations")
##         MSE1 <- DLMtool:::runMSE_int(OM, MPs, CheckMPs, timelimit, Hist,
##             ntrials, fracD, CalcBlow, HZN, Bfrac, AnnualMSY,
##             silent, PPD, checks = checks, control = control)
##     }
##     if (class(MSE1) == "MSE") {
##         if (class(MSE1@Misc$TryMP) == "list") {
##             ok <- unlist(MSE1@Misc$TryMP) == "Okay"
##             fail <- unlist(MSE1@Misc$TryMP)
##         }
##         if (class(MSE1@Misc$TryMP) == "matrix") {
##             ok <- colSums(MSE1@Misc$TryMP == "Okay") == nrow(MSE1@Misc$TryMP)
##             fail <- t(MSE1@Misc$TryMP)
##             if (any(grepl("could not find function", unique(fail[!ok,
##                 ])))) {
##                 warning("MPs may have been dropped because of non-exported functions in parallel mode. \nUse `setup(); snowfall::sfExport('FUNCTION1', 'FUNCTION2')` to export functions to cores")
##             }
##         }
##         if (any(!ok)) {
##             failedMPs <- MSE1@MPs[!ok]
##             warning("Dropping failed MPs: ", paste(failedMPs,
##                 collapse = ", "), "\n\nSee MSE@Misc$TryMP for error messages\n\n")
##             if (length(failedMPs) == MSE1@nMPs) {
##                 warning("All MPs failed.")
##                 return(MSE1)
##             }
##             MSE1 <- Sub(MSE1, MPs = MSE1@MPs[!MSE1@MPs %in% failedMPs])
##         }
##     }
##     return(MSE1)
## }



runMSE <- function (OM = DLMtool::testOM, MPs = c("AvC", "DCAC", "FMSYref",
    "curE", "matlenlim", "MRreal"), CheckMPs = FALSE, timelimit = 1,
    Hist = FALSE, ntrials = 100, fracD = 0.05, CalcBlow = TRUE,
    HZN = 2, Bfrac = 0.5, AnnualMSY = TRUE, silent = FALSE, PPD = TRUE,
    parallel = FALSE, save_name = NULL, checks = FALSE, control = NULL)
{
    if (class(OM) != "OM")
        stop("OM is not class 'OM'", call. = FALSE)
    rm(list = ls(DLMenv), envir = DLMenv)
    tt <- suppressWarnings(try(lsf.str(envir = globalenv()),
        silent = TRUE))
    if (class(tt) != "try-error") {
        gl.funs <- as.vector(tt)
        pkg.funs <- as.vector(ls.str("package:DLMtool"))
        if ("package:MSEtool" %in% search())
            pkg.funs <- c(pkg.funs, as.vector(ls.str("package:MSEtool")))
        if (length(gl.funs) > 0) {
            gl.clss <- unlist(lapply(lapply(gl.funs, get), class))
            gl.MP <- gl.funs[gl.clss %in% "MP"]
            if (length(gl.MP) > 0) {
                inc.gl <- gl.MP[gl.MP %in% MPs]
                if (length(inc.gl) > 0) {
                  dup.MPs <- inc.gl[inc.gl %in% pkg.funs]
                  if (length(dup.MPs) > 0) {
                    stop("Custom MP names already in DLMtool: ",
                      paste0(dup.MPs, " "), "\nRename Custom MPs")
                  }
                }
            }
        }
    }
    if (!all(is.na(MPs))) {
        for (mm in MPs) {
            chkMP <- try(get(mm), silent = TRUE)
            if (class(chkMP) != "MP")
                stop(mm, " is not a valid MP", call. = FALSE)
        }
    }
    if (parallel) {
        if (OM@nsim < 48)
            stop("nsim must be >=48 for parallel processing",
                call. = FALSE)
        if (!snowfall::sfIsRunning()) {
            message("Parallel processing hasn't been initialized. Calling 'setup()' now")
            setup()
        }
        if (all(is.na(MPs)))
            MPs <- DLMtool::avail("MP")
        cMPs <- MPs[!MPs %in% pkg.funs]
        globalMP <- NULL
        extra_package <- NULL
        for (mm in seq_along(cMPs)) {
            nmspace <- utils::find(cMPs[mm])
            if (nmspace == ".GlobalEnv") {
                globalMP <- c(globalMP, cMPs[mm])
            }
            else {
                extra_package <- c(extra_package, strsplit(nmspace,
                  ":")[[1]][2])
            }
            extra_package <- unique(extra_package)
        }
        if (!is.null(globalMP)) {
            message("Exporting custom MPs in global environment")
            snowfall::sfExport(list = globalMP)
        }
        if (!is.null(extra_package)) {
            message("Exporting additional packages with MPs")
            for (pk in extra_package) sfLibrary(pk, character.only = TRUE,
                verbose = FALSE)
        }
        ncpu <- snowfall::sfCpus()
        nits <- ceiling(OM@nsim/5)
        itsim <- rep(5, nits)
        if (nits < ncpu) {
            if (nits < 4) {
                nits <- 4
                itsim <- rep(ceiling(OM@nsim/4), 4)
            }
            else {
                nits <- ncpu
                itsim <- rep(ceiling(OM@nsim/ncpu), ncpu)
            }
        }
        cnt <- 1
        while (sum(itsim) != OM@nsim | any(itsim < 2)) {
            diff <- OM@nsim - sum(itsim)
            if (diff > 0) {
                itsim[cnt] <- itsim[cnt] + 1
            }
            if (diff < 0) {
                itsim[cnt] <- itsim[cnt] - 1
            }
            cnt <- cnt + 1
            if (cnt > length(itsim))
                cnt <- 1
        }
        if (!silent & !Hist)
            message("Running MSE in parallel on ", ncpu, " processors")
        if (!silent & Hist)
            message("Running historical simulations in parallel on ",
                ncpu, " processors")
        temp <- snowfall::sfClusterApplyLB(1:nits, DLMtool:::run_parallel,
            itsim = itsim, OM = OM, MPs = MPs, CheckMPs = CheckMPs,
            timelimit = timelimit, Hist = Hist, ntrials = ntrials,
            fracD = fracD, CalcBlow = CalcBlow, HZN = HZN, Bfrac = Bfrac,
            AnnualMSY = AnnualMSY, silent = TRUE, PPD = PPD,
            control = control, parallel = parallel)
        if (!is.null(save_name) && is.character(save_name))
            saveRDS(temp, paste0(save_name, ".rdata"))
        MSE1 <- DLMtool::joinMSE(temp)
        if (class(MSE1) == "MSE") {
            if (!silent)
                message("MSE completed")
        }
        else if (class(MSE1) == "Hist") {
            if (!silent)
                message("Historical simulations completed")
        }
        else {
            warning("MSE completed but could not join MSE objects. Re-run with `save_name ='MyName'` to debug")
        }
    }
    if (!parallel) {
        if (OM@nsim > 48 & !silent & !Hist)
            message("Suggest using 'parallel = TRUE' for large number of simulations")
        MSE1 <- DLMtool:::runMSE_int(OM, MPs, CheckMPs, timelimit, Hist,
            ntrials, fracD, CalcBlow, HZN, Bfrac, AnnualMSY,
            silent, PPD, checks = checks, control = control)
    }
    if (class(MSE1) == "MSE") {
        if (class(MSE1@Misc$TryMP) == "list") {
            ok <- unlist(MSE1@Misc$TryMP) == "Okay"
            fail <- unlist(MSE1@Misc$TryMP)
        }
        if (class(MSE1@Misc$TryMP) == "matrix") {
            ok <- colSums(MSE1@Misc$TryMP == "Okay") == nrow(MSE1@Misc$TryMP)
            fail <- t(MSE1@Misc$TryMP)
            if (any(grepl("could not find function", unique(fail[!ok,
                ])))) {
                warning("MPs may have been dropped because of non-exported functions in parallel mode. \nUse `setup(); snowfall::sfExport('FUNCTION1', 'FUNCTION2')` to export functions to cores")
            }
        }
        if (any(!ok)) {
            failedMPs <- MSE1@MPs[!ok]
            warning("Dropping failed MPs: ", paste(failedMPs,
                collapse = ", "), "\n\nSee MSE@Misc$TryMP for error messages\n\n")
            if (length(failedMPs) == MSE1@nMPs) {
                warning("All MPs failed.")
                return(MSE1)
            }
            MSE1 <- DLMtool::Sub(MSE1, MPs = MSE1@MPs[!MSE1@MPs %in% failedMPs])
        }
    }
    return(MSE1)
}





## December 2019
plotAll <- function(MSElist, evalyears, legends = NULL,
                             xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                    ptitles=rep("",length(MSElist))){
    nscen = length(MSElist)
    if(is.null(xlimUpper)) xlimUpper <- rep(50,10)
    if(is.null(xlimLower)) xlimLower <- rep(0,10)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,10)
    if(is.null(ylimLower)) ylimLower <- rep(0,10)
    ly <- layout(matrix(c(1,2,9,3,4,9,5,6,9,7,8,9),
                        byrow = TRUE,ncol=3,nrow=4),
                 widths=c(1,1,0.44,1,1,0.44,1,1,0.44,1,1,0.44),
                 heights=rep(1,5))
    opar <- par(mar=c(3,2,3,2),oma=c(5,5.5,1,0.5))
    legends=legends
    pchs <- c(3,17,15,16)
    argroup <- c(1,2,
                 rep(3,9),
                 rep(4,9))
    arcols <- c(1,1,
                paste0(rep("grey",9),round(seq(10,75,length.out = 9))),
                paste0(rep("grey",9),round(seq(10,75,length.out = 9))))
    aralphas <- c(1,1,
                  seq(1,0.3,length.out = 10)[-1],
                  seq(1,0.3,length.out = 10)[-1])
    connect <- list(c(2:11),c(2,12:21))
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[1],ylimUpper[1]),
             xlim=c(xlimLower[1],xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[1],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[2]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[2],ylimUpper[2]),
             xlim=c(xlimLower[2],xlimUpper[2]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[2],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[3]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[3],ylimUpper[3]),
             xlim=c(xlimLower[3],xlimUpper[3]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[3],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[4]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[4],ylimUpper[4]), xlim=c(xlimLower[4],xlimUpper[4]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[4],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[5]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[5],ylimUpper[5]), xlim=c(xlimLower[5],xlimUpper[5]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[5],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[6]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[6],ylimUpper[6]), xlim=c(xlimLower[6],xlimUpper[6]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
##    abline(v=12.380952,lty=3)
    mtext(ptitles[6],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[7]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[7],ylimUpper[7]), xlim=c(xlimLower[7],xlimUpper[7]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[7],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[8]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[8],ylimUpper[8]), xlim=c(xlimLower[8],xlimUpper[8]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
##    abline(v=5,lty=1,col="grey80")
##    abline(v=12.380952,lty=3)
    mtext(ptitles[8],
          font=2,line=0.5,cex=1.3)
    ## additional
    ## plotdYdR(MSElist[[9]], Yrs=evalyears, risk = 1,
    ##          connect=connect,
    ##          ylim = c(ylimLower[9],ylimUpper[9]), xlim=c(xlimLower[9],xlimUpper[9]),
    ##          pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=2)
##     mtext(ptitles[9],
##           font=2,line=0.5,cex=1.3)
##     plotdYdR(MSElist[[10]], Yrs=evalyears, risk = 1,
##              connect=connect,
##              ylim = c(ylimLower[10],ylimUpper[10]), xlim=c(xlimLower[10],xlimUpper[10]),
##              pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
##     abline(v=5,lty=2)
## ##   abline(v=12.380952,lty=3)
##     mtext(ptitles[10],
##           font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 1.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.9, outer=TRUE, cex=1.2)
    plot.new()
    legend(-0.05,0.68, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.5)
    legend(-0.05,0.32, title = "MIAVY", title.adj=0,
           legend=c("<0.1","0.1-0.3"," >=0.3"), cex=1.5,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}




require(dplyr)
if(FALSE){
joinMSE <- function (MSEobjs = NULL)
{
    if (class(MSEobjs) != "list")
        stop("MSEobjs must be a list")
    if (length(MSEobjs) < 2)
        stop("MSEobjs list doesn't contain multiple MSE objects")
    lapply(MSEobjs, checkMSE)
    ishist <- all(lapply(MSEobjs, slotNames) %>% unlist() %>%
        unique() %in% slotNames("Hist"))
    if (ishist) {
        out <- new("Hist")
        sls <- slotNames("Hist")
        nsim <- MSEobjs[[1]]@Ref$B0 %>% length()
        for (sl in sls) {
            obj <- lapply(MSEobjs, slot, name = sl)
            if (sl == "Data") {
                out@Data <- joinData(obj)
            }
            else {
                if (class(obj[[1]]) == "data.frame") {
                  slot(out, sl) <- do.call("rbind", obj)
                }
                if (class(obj[[1]]) == "list") {
                  out.list <- list()
                  for (nm in names(obj[[1]])) {
                    obj2 <- lapply(obj, "[[", nm)
                    ind <- which(dim(obj2[[1]]) == nsim)
                    if (length(ind) > 1)
                      ind <- ind[1]
                    if (length(ind) > 0) {
                      if (class(obj2[[1]]) == "array") {
                        tempVal <- lapply(obj2, dim)
                        tdf <- lapply(obj2, dim) %>% unlist() %>%
                          matrix(nrow = length(obj2), ncol = length(dim(obj2[[1]])),
                            byrow = TRUE)
                        nBins <- tdf[, 2]
                        Max <- max(nBins)
                        nyrs <- max(tdf[, 3])
                        nsims <- sapply(tempVal, function(x) x[1])
                        if (!mean(nBins) == max(nBins)) {
                          index <- which(nBins < Max)
                          for (kk in index) {
                            dif <- Max - dim(obj2[[kk]])[2]
                            obj2[[kk]] <- abind::abind(obj2[[kk]],
                              array(0, dim = c(nsims[kk], dif,
                                nyrs)), along = 2)
                          }
                        }
                      }
                      out.list[[nm]] <- abind::abind(obj2, along = ind)
                    }
                    else {
                      out.list[[nm]] <- unlist(obj2)
                    }
                  }
                  slot(out, sl) <- out.list
                }
            }
        }
        return(out)
    }
    MPNames <- lapply(MSEobjs, getElement, name = "MPs")
    allsame <- length(unique(lapply(MPNames, unique))) == 1
    if (!allsame) {
        mpnames <- unlist(MPNames)
        npack <- length(MSEobjs)
        tab <- table(mpnames)
        ind <- tab == npack
        commonMPs <- names(tab)[ind]
        if (length(commonMPs) < 1)
            stop("No common MPs in MSE objects", call. = FALSE)
        MSEobjs <- lapply(MSEobjs, Sub, MPs = commonMPs)
        message("MPs not in all MSE objects:")
        message(paste(names(tab)[!ind], ""))
        message("Dropped from final MSE object.")
    }
    Nobjs <- length(MSEobjs)
    for (X in 1:Nobjs) {
        tt <- MSEobjs[[X]]
        assign(paste0("obj", X), tt)
        if (X > 1) {
            tt <- MSEobjs[[X]]
            tt2 <- MSEobjs[[X - 1]]
            if (!all(slotNames(tt) == slotNames(tt2)))
                stop("The MSE objects don't have the same slots")
            if (any(tt@MPs != tt2@MPs))
                stop("MPs must be the same for all MSE objects")
        }
    }
    chkmat <- matrix(NA, nrow = Nobjs, ncol = 2)
    nms <- NULL
    for (X in 1:Nobjs) {
        tt <- get(paste0("obj", X))
        chkmat[X, ] <- c(tt@nyears, tt@proyears)
        if (X > 1)
            if (!any(grepl(tt@Name, nms)))
                stop("MSE objects have different names")
        nms <- append(nms, tt@Name)
    }
    chk <- all(colSums(chkmat) == chkmat[1, ] * Nobjs)
    if (!chk)
        stop("The MSE objects have different number of nyears or proyears")
    Allobjs <- mget(paste0("obj", 1:Nobjs))
    sns <- slotNames(Allobjs[[1]])
    sns <- sns[sns != "Misc"]
    outlist <- vector("list", length(sns))
    for (sn in 1:length(sns)) {
        templs <- lapply(Allobjs, slot, name = sns[sn])
        if (class(templs[[1]]) == "character") {
            outlist[[sn]] <- templs[[1]]
        }
        if (class(templs[[1]]) == "numeric" | class(templs[[1]]) ==
            "integer") {
            if (sns[sn] == "CALbins") {
                tempInd <- which.max(unlist(lapply(templs, length)))
                CALbins <- templs[[tempInd]]
            }
            else {
                outlist[[sn]] <- do.call(c, templs)
            }
        }
        if (class(templs[[1]]) == "matrix" | class(templs[[1]]) ==
            "data.frame") {
            outlist[[sn]] <- do.call(rbind, templs)
        }
        if (class(templs[[1]]) == "array") {
            if (sns[sn] == "CAL") {
                tempVal <- lapply(templs, dim)
                if (all(unlist(lapply(tempVal, length)) == 3)) {
                  nBins <- sapply(tempVal, function(x) x[3])
                  nsims <- sapply(tempVal, function(x) x[1])
                  nMPs <- sapply(tempVal, function(x) x[2])
                  if (!mean(nBins) == max(nBins)) {
                    Max <- max(nBins)
                    index <- which(nBins < Max)
                    for (kk in index) {
                      dif <- Max - dim(templs[[kk]])[3]
                      templs[[kk]] <- abind::abind(templs[[kk]],
                        array(0, dim = c(nsims[kk], nMPs[kk],
                          dif)), along = 3)
                    }
                  }
                  outlist[[sn]] <- abind::abind(templs, along = 1)
                }
                else {
                  outlist[[sn]] <- templs[[1]]
                }
            }
            else {
                outlist[[sn]] <- abind::abind(templs, along = 1)
            }
        }
    }
    names(outlist) <- sns
    Misc <- list()
    if (length(MSEobjs[[1]]@Misc) > 0) {
        if (!is.null(MSEobjs[[1]]@Misc$Data)) {
            Misc$Data <- list()
            for (i in 1:length(MSEobjs[[1]]@Misc$Data)) Misc$Data[[i]] <- joinData(lapply(MSEobjs,
                function(x) slot(x, "Misc")$Data[[i]]))
        }
        if (!is.null(MSEobjs[[1]]@Misc$RInd.stats)) {
            nms <- unique(MSEobjs[[1]]@Misc$RInd.stats$Index) %>%
                as.character()
            temp <- list()
            for (nm in seq_along(nms)) {
                temp1 <- list()
                for (i in 1:length(MSEobjs)) {
                  temp1[[i]] <- MSEobjs[[i]]@Misc$RInd.stats %>%
                    dplyr::filter(Index == nms[nm])
                }
                temp[[nm]] <- do.call("rbind", temp1)
            }
            Misc$RInd.stats <- do.call("rbind", temp)
        }
        if (!is.null(MSEobjs[[1]]@Misc$TryMP)) {
            temp1 <- list()
            for (i in 1:length(MSEobjs)) {
                temp1[[i]] <- MSEobjs[[i]]@Misc$TryMP
            }
            Misc$TryMP <- do.call("rbind", temp1)
        }
        if (!is.null(MSEobjs[[1]]@Misc$Unfished)) {
            Misc$Unfished <- list()
            temp1 <- temp2 <- list()
            for (i in 1:length(MSEobjs)) {
                temp1[[i]] <- MSEobjs[[i]]@Misc$Unfished$Refs
                temp2[[i]] <- MSEobjs[[i]]@Misc$Unfished$ByYear
            }
            Misc$Unfished$Refs <- do.call("cbind", temp1)
            for (nm in names(temp2[[1]])) {
                tt = lapply(temp2, "[[", nm)
                tt <- do.call("rbind", tt)
                Misc$Unfished$ByYear[[nm]] <- tt
            }
        }
        if (!is.null(MSEobjs[[1]]@Misc$MSYRefs)) {
            Misc$MSYRefs <- list()
            temp1 <- temp2 <- list()
            for (i in 1:length(MSEobjs)) {
                temp1[[i]] <- MSEobjs[[i]]@Misc$MSYRefs$Refs
                temp2[[i]] <- MSEobjs[[i]]@Misc$MSYRefs$ByYear
            }
            Misc$MSYRefs$Refs <- do.call("rbind", temp1)
            for (nm in names(temp2[[1]])) {
                tt = lapply(temp2, "[[", nm)
                tt <- do.call("rbind", tt)
                Misc$MSYRefs$ByYear[[nm]] <- tt
            }
        }
        temp <- list()
        nsim <- ncol(Misc$Unfished$Ref)
        dims <- dim(MSEobjs[[1]]@Misc$LatEffort)
        Misc$LatEffort <- array(NA, dim = c(nsim, dims[2], dims[3]))
        Misc$Revenue <- array(NA, dim = c(nsim, dims[2], dims[3]))
        Misc$Cost <- array(NA, dim = c(nsim, dims[2], dims[3]))
        Misc$TAE <- array(NA, dim = c(nsim, dims[2], dims[3]))
        st <- 1
        for (i in 1:length(MSEobjs)) {
            dims <- dim(MSEobjs[[i]]@Misc$LatEffort)
            indvec <- st:(st + dims[1] - 1)
            st <- indvec[length(indvec)] + 1
            Misc$LatEffort[indvec, , ] <- MSEobjs[[i]]@Misc$LatEffort
            Misc$Cost[indvec, , ] <- MSEobjs[[i]]@Misc$Cost
            Misc$Revenue[indvec, , ] <- MSEobjs[[i]]@Misc$Revenue
            Misc$TAE[indvec, , ] <- MSEobjs[[i]]@Misc$TAE
        }

        ## spict specific (here: for number of non-converged assessments)
        resList <- vector("list",MSEobjs[[1]]@nMPs)
        for(i in 1:MSEobjs[[1]]@nMPs){
            tmpList <- list()
            for(j in 1:length(MSEobjs)){
                tmpList <- c(tmpList, MSEobjs[[j]]@Misc$spict[[i]])
                             ##lapply(MSEobjs[[j]]@Misc[[4]][[i]]@Misc, "[[", "spict"))
            }
            resList[[i]] <- tmpList
        }
        Misc$spict <- resList
    }
    newMSE <- new("MSE", Name = outlist$Name, nyears = unique(outlist$nyears),
        proyears = unique(outlist$proyears), nMP = unique(outlist$nMP),
        MPs = unique(outlist$MPs), nsim = sum(outlist$nsim),
        OM = outlist$OM, Obs = outlist$Obs, B_BMSY = outlist$B_BMSY,
        F_FMSY = outlist$F_FMSY, outlist$B, outlist$SSB, outlist$VB,
        outlist$FM, outlist$C, outlist$TAC, outlist$SSB_hist,
        outlist$CB_hist, outlist$FM_hist, outlist$Effort, outlist$PAA,
        outlist$CAA, outlist$CAL, CALbins, Misc = Misc)
    newMSE
}
}


plotAllAll <- function(MSElistin, evalyears, legends = NULL, relref =FALSE,
                       xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                       ptitles=rep("",length(MSElist))){
    nscen = length(MSElistin)
    if(is.null(xlimUpper)) xlimUpper <- rep(50,10)
    if(is.null(xlimLower)) xlimLower <- rep(0,10)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,10)
    if(is.null(ylimLower)) ylimLower <- rep(0,10)
    ly <- layout(matrix(c(c(c(1,3,5,7,2,4,6,8),rep(9,4)),
                          c(c(10,12,14,16,11,13,15,17),rep(18,4)),
                          c(c(19,21,23,25,20,22,24,26),rep(27,4))),
                        byrow = FALSE, ncol=9, nrow=4),
                 widths=c(rep(1,2),0.25,rep(1,2),0.25,rep(1,2),0.8),
                 heights=rep(1,4))
    ##    opar <- par(mar=c(3,2,3,2),oma=c(5,5.5,1,0.5))
    opar <- par(mar=c(0,0,0,0),oma=c(9,9,8,5))
    legends=legends
    pchs <- c(3,17,15,16)
    argroup <- c(1,2,
                 rep(3,9),
                 rep(4,9))
    arcols <- c(1,1,
                paste0(rep("grey",9),round(seq(10,75,length.out = 9))),
                paste0(rep("grey",9),round(seq(10,75,length.out = 9))))
    aralphas <- c(1,1,
                  seq(1,0.3,length.out = 10)[-1],
                  seq(1,0.3,length.out = 10)[-1])
    connect <- list(c(2:11),c(2,12:21))
    if(relref){
        aralphas <- c(1,
                      seq(1,0.3,length.out = 10)[-1],
                      seq(1,0.3,length.out = 10)[-1])
        connect <- list(c(1:10),c(1,11:20))
    }
    for(j in 1:length(MSElistin)){
        MSElist <- MSElistin[[j]]
        ## first column
        for(k in 1:8){
            indi <- k + ((j-1) * 8)
            xaxt <- "n" ## ifelse(k %in% c(7,8), "s", "n")
            yaxt <- "n" ## ifelse(k %in% c(1,3,5,7), "s", "n")
            plotdYdR(MSElist[[k]], Yrs=evalyears, risk = 1, relref = relref,
                     connect=connect, ylim = c(ylimLower[indi],ylimUpper[indi]),
                     xlim=c(xlimLower[indi],xlimUpper[indi]), xaxt=xaxt, yaxt=yaxt,
                     pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
            text(par("usr")[1]+0.02*xlimUpper[indi],
                 par("usr")[4]-0.06*par("usr")[4],
                 srt=0, adj = 0, cex=1.6, font=2,
                 labels = ptitles[k])
            if(k == 1){
                mtext(speciesAll[j], side =3, font=2, adj=c(1.6,1.6,1.15)[j],
                      line=3, cex=1.5)
            }
        }
        if(j < length(MSElistin)) plot.new()
    }
    if(!relref) mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
                     side = 1, line = 4.5, outer=TRUE, adj = 0.45, cex=1.3)
    if(relref) mtext(bquote("Rel." ~ risk ~ "[years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 4.5, outer=TRUE, adj = 0.45, cex=1.3)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 4.5, outer=TRUE, cex=1.3)
    plot.new()
    if(relref){
        legends <- legends[-1]
        argroup <- argroup[-1]
        aralphas <- aralphas[-1]
    }
    legend(0.2,0.9, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.5)
    legend(0.2,0.24, title = "MIAVY", title.adj=0,
           legend=c("<0.1","0.1-0.3"," >=0.3"), cex=1.5,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}



plotFM2 <- function(MSE, hcrsNum=1,cols="dodgerblue4",
                    ylim=NA,
                    ylab="F/Fmsy",xlab="Years",
                    xaxt="s",yaxt="s"){
    nyears <- MSE@nyears
    proyears <- MSE@proyears
    nsim <- MSE@nsim
    ## FM
    fmhist <- MSE@FM_hist ## nsim, nages, nyears, nareas.
    fm <- MSE@FM          ## nsim, nMPs, proyears
    fmhistA <- apply(fmhist, c(1,3), mean)  ## nsim, nyears
    ## Fmsy (for each sim)
    fmsy <- MSE@OM$FMSY
    fmfmsy <- fm / fmsy
    fmhistfmsy <- fmhistA / fmsy
    statshist <- apply(fmhistfmsy[, , drop = FALSE], 2,
                       quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
    ## plot hist
    if(is.na(ylim[1])) ylim <- range(as.numeric(fmhistfmsy[,]))
    plot(fmhistfmsy[1,],ty='n',
         ylab=ylab,xlab=xlab,
         xlim=c(1,nyears+proyears),
         ylim=ylim,xaxt=xaxt,yaxt=yaxt)
    polygon(x = c(1:(nyears), (nyears):1),
            y = c(c(statshist[1, ]), rev(c(statshist[3,]))),
            col = rgb(t(col2rgb("grey50"))/255,alpha=0.5), border = FALSE)
    lines(1:(nyears), c(statshist[2, ]), lwd = 2)
    ## hcrs
    for(i in 1:length(hcrsNum)){
        MPnum <- hcrsNum[i]
        stats <- apply(fmfmsy[, MPnum, , drop = FALSE], 3,
                       quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        polygon(x = c((nyears+1):(nyears+proyears), rev((nyears+1):(nyears+proyears))),
                y = c(stats[1, ], rev(stats[3,])),
                col = rgb(t(col2rgb(cols[i]))/255,alpha=0.5),
                border = FALSE)
    }
    for(i in 1:length(hcrsNum)){
        MPnum <- hcrsNum[i]
        stats <- apply(fmfmsy[, MPnum, , drop = FALSE], 3,
                       quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
        lines((nyears+1):(nyears+proyears), stats[2, ], lwd = 2,col=cols[i])
    }
    abline(v=nyears, col="grey60",lty=2,lwd=1.5)
    abline(h=1, col="grey60", lwd=1.5,lty=2)
    box()
}


plotAllAll2 <- function(MSElistin, evalyears, legends = NULL, relref =FALSE,
                       xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                       ptitles=rep("",length(MSElist))){
    evalyearsin <- evalyears
    nscen = length(MSElistin)
    if(is.null(xlimUpper)) xlimUpper <- rep(50,10)
    if(is.null(xlimLower)) xlimLower <- rep(0,10)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,10)
    if(is.null(ylimLower)) ylimLower <- rep(0,10)
    ly <- layout(matrix(c(c(c(1,3,5,7,2,4,6,8),rep(9,4)),
                          c(c(10,12,14,16,11,13,15,17),rep(18,4)),
                          c(c(19,21,23,25,20,22,24,26),rep(27,4))),
                        byrow = FALSE, ncol=9, nrow=4),
                 widths=c(rep(1,2),0.25,rep(1,2),0.25,rep(1,2),0.8),
                 heights=rep(1,4))
    ##    opar <- par(mar=c(3,2,3,2),oma=c(5,5.5,1,0.5))
    opar <- par(mar=c(0,0,0,0),oma=c(7,7,8,4))
    legends=legends
    pchs <- c(3,17,15,16)
    argroup <- c(1,2,
                 rep(3,9),
                 rep(4,9))
    arcols <- c(1,1,
                colorRampPalette(c("gold4","gold1"))(9), ##paste0(rep("goldenrod",9),round(seq(10,75,length.out = 9))),
                colorRampPalette(c("dodgerblue4","dodgerblue1"))(9)) ##paste0(rep("dodgerblue",9),round(seq(10,75,length.out = 9))))
    aralphas <- c(1,1,
                  seq(1,0.3,length.out = 10)[-1],
                  seq(1,0.3,length.out = 10)[-1])
    connect <- list(c(2:11),c(2,12:21))
    if(relref){
        aralphas <- c(1,
                      seq(1,0.3,length.out = 10)[-1],
                      seq(1,0.3,length.out = 10)[-1])
        connect <- list(c(1:10),c(1,11:20))
    }
    for(j in 1:length(MSElistin)){
        if(evalyearsin[1] == "flex"){
            evalyears <- list(c(2,6),c(3,9),c(6,20))[[j]]
        }else if(evalyearsin[1] == "flex1"){
            evalyears <- list(c(1,6),c(1,9),c(1,20))[[j]]
        }
        MSElist <- MSElistin[[j]]
        ## first column
        for(k in 1:8){
            indi <- k + ((j-1) * 8)
            xaxt <- "n" ## ifelse(k %in% c(7,8), "s", "n")
            yaxt <- "n" ## ifelse(k %in% c(1,3,5,7), "s", "n")
            plotdYdR(MSElist[[k]], Yrs=evalyears, risk = 1, relref = relref,
                     connect=connect, ylim = c(ylimLower[indi],ylimUpper[indi]),
                     xlim=c(xlimLower[indi],xlimUpper[indi]), xaxt=xaxt, yaxt=yaxt,
                     pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
            plotrix::boxed.labels(par("usr")[1]+0.07*xlimUpper[indi],
                                  par("usr")[4]-0.06*par("usr")[4],
                                  ptitles[k], cex = 1.1,
                                  border = NA,
                                  bg ="white", font=2,
                                  xpad = 1.3, ypad = 1.3)
            box()
            ## text(par("usr")[1]+0.02*xlimUpper[indi],
            ##      par("usr")[4]-0.07*par("usr")[4],
            ##      srt=0, adj = 0, cex=1.6, font=2,
            ##      labels = ptitles[k])
            if(k == 1){
                mtext(speciesAll[j], side =3, font=2, adj=c(1.55,1.55,1.18)[j],
                      line=3, cex=1.5)
            }
        }
        if(j < length(MSElistin)) plot.new()
    }
    if(!relref && !evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
                     side = 1, line = 3.5, outer=TRUE, adj = 0.45, cex=1.3)
    }else if(relref && !evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote("Rel." ~ risk ~ "[years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 3.5, outer=TRUE, adj = 0.43, cex=1.3)
    }else{
        mtext("Rel. risk",
          side = 1, line = 3.5, outer=TRUE, adj = 0.43, cex=1.3)
    }
    if(!evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
              side = 2, line = 3, outer=TRUE, cex=1.3)
    }else{
        mtext("Rel. yield",
              side = 2, line = 3, outer=TRUE, cex=1.3)
    }
    plot.new()
    if(relref){
        legends <- legends[-1]
        argroup <- argroup[-1]
        aralphas <- aralphas[-1]
    }
    legend(0.2,0.9, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=c(1,colorRampPalette(c("gold4","gold1"))(9),
                 colorRampPalette(c("dodgerblue4","dodgerblue1"))(9)), ##rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.5)
    legend(0.2,0.26, title = "MIAVY", title.adj=0,
           legend=c("<0.1","0.1-0.3"," >=0.3"), cex=1.5,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}



## final: vertical arrangement of species absolute axes
plotAllAll3 <- function(MSElistin, evalyears, legends = NULL, relref =FALSE,
                       xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                       ptitles=rep("",length(MSElist))){
    evalyearsin <- evalyears
    nscen = length(MSElistin)
    if(is.null(xlimUpper)) xlimUpper <- rep(50,10)
    if(is.null(xlimLower)) xlimLower <- rep(0,10)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,10)
    if(is.null(ylimLower)) ylimLower <- rep(0,10)
    mat = matrix(c(1:4,27,5:8,27,
                   rep(9,4),27,
                   10:13,27,14:17,27,
                   rep(18,4),27,
                   19:22,27,23:26,27),
                 byrow = TRUE, ncol=5, nrow=8)
    ly <- layout(mat,
                 widths=c(1,1,1,1,0.82),
                 heights=rep(c(1,1,0.4),3))
    ##    opar <- par(mar=c(3,2,3,2),oma=c(5,5.5,1,0.5))
    opar <- par(mar=c(2,0,1,0),oma=c(5,7.3,5.5,4))
    legends=legends
    pchs <- c(3,17,15,16)
    argroup <- c(1,2,
                 rep(3,9),
                 rep(4,9))
    arcols <- c(1,1,
                colorRampPalette(c("gold4","gold1"))(9), ##paste0(rep("goldenrod",9),round(seq(10,75,length.out = 9))),
                colorRampPalette(c("dodgerblue4","dodgerblue1"))(9)) ##paste0(rep("dodgerblue",9),round(seq(10,75,length.out = 9))))
    aralphas <- c(1,1,
                  seq(1,0.3,length.out = 10)[-1],
                  seq(1,0.3,length.out = 10)[-1])
    connect <- list(c(2:11),c(2,12:21))
    if(relref){
        aralphas <- c(1,
                      seq(1,0.3,length.out = 10)[-1],
                      seq(1,0.3,length.out = 10)[-1])
        connect <- list(c(1:10),c(1,11:20))
    }

    for(j in 1:length(MSElistin)){
        if(evalyearsin[1] == "flex"){
            evalyears <- list(c(2,6),c(3,9),c(6,20))[[j]]
        }else if(evalyearsin[1] == "flex1"){
            evalyears <- list(c(1,6),c(1,9),c(1,20))[[j]]
        }
        MSElist <- MSElistin[[j]]
        ## first column
        for(k in 1:8){
            indi <- k + ((j-1) * 8)
            xaxt <- "n" ## ifelse(k %in% c(7,8), "s", "n")
            yaxt <- "n" ## ifelse(k %in% c(1,3,5,7), "s", "n")
            plotdYdR(MSElist[[k]], Yrs=evalyears, risk = 1, relref = relref,
                     connect=connect, ylim = c(ylimLower[indi],ylimUpper[indi]),
                     xlim=c(xlimLower[indi],xlimUpper[indi]), xaxt=xaxt, yaxt=yaxt,
                     pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
            axis(1, axis.cex=1.2)
            if(k %in% c(1,5)) axis(2, axis.cex=1.2)
            plotrix::boxed.labels(par("usr")[1]+0.08*xlimUpper[indi],
                                  par("usr")[4]-0.07*par("usr")[4],
                                  ptitles[k], cex = 0.9,
                                  border = NA,
                                  bg ="white", font=2,
                                  xpad = 1.3, ypad = 1.3)
            box()
            ## text(par("usr")[1]+0.02*xlimUpper[indi],
            ##      par("usr")[4]-0.07*par("usr")[4],
            ##      srt=0, adj = 0, cex=1.6, font=2,
            ##      labels = ptitles[k])
            if(k == 2){
                mtext(speciesAll[j], side =3, font=2, adj=c(1.55,1.55,1.14)[j],
                      line=1.5, cex=1.2)
            }
        }
        if(j < length(MSElistin)) plot.new()
    }

    if(!relref && !evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
                     side = 1, line = 2.5, outer=TRUE, adj = 0.42, cex=1.2)
    }else if(relref && !evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote("Rel." ~ risk ~ "[years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 2.5, outer=TRUE, adj = 0.4, cex=1.2)
    }else{
        mtext("Rel. risk",
          side = 1, line = 2.5, outer=TRUE, adj = 0.4, cex=1.2)
    }
    if(!evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
              side = 2, line = 3.6, outer=TRUE, cex=1.2)
    }else{
        mtext("Rel. yield",
              side = 2, line = 3.6, outer=TRUE, cex=1.2)
    }
    plot.new()
    if(relref){
        legends <- legends[-1]
        argroup <- argroup[-1]
        aralphas <- aralphas[-1]
    }
    legend(0.2,0.715, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=c(1,colorRampPalette(c("gold4","gold1"))(9),
                 colorRampPalette(c("dodgerblue4","dodgerblue1"))(9)), ##rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.32)
    legend(0.2,0.358, title = "MIAVY", title.adj=0,
           legend=c("<0.1","0.1-0.3"," >=0.3"), cex=1.32,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}


## final: vertical arrangement of species absolute axes
plotAllAll3Abs <- function(MSElistin, evalyears, legends = NULL, relref =FALSE,
                       xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                       ptitles=rep("",length(MSElist))){
    evalyearsin <- evalyears
    nscen = length(MSElistin)
    if(is.null(xlimUpper)) xlimUpper <- rep(50,10)
    if(is.null(xlimLower)) xlimLower <- rep(0,10)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,10)
    if(is.null(ylimLower)) ylimLower <- rep(0,10)
    mat = matrix(c(1:4,27,5:8,27,
                   rep(9,4),27,
                   10:13,27,14:17,27,
                   rep(18,4),27,
                   19:22,27,23:26,27),
                 byrow = TRUE, ncol=5, nrow=8)
    ly <- layout(mat,
                 widths=c(1,1,1,1,0.85),
                 heights=rep(c(1,1,0.4),3))
    ##    opar <- par(mar=c(3,2,3,2),oma=c(5,5.5,1,0.5))
    opar <- par(mar=c(2,3,1,0),oma=c(5,5,5.5,4))
    legends=legends
    pchs <- c(3,17,15,16)
    argroup <- c(1,2,
                 rep(3,9),
                 rep(4,9))
    arcols <- c(1,1,
                colorRampPalette(c("gold4","gold1"))(9), ##paste0(rep("goldenrod",9),round(seq(10,75,length.out = 9))),
                colorRampPalette(c("dodgerblue4","dodgerblue1"))(9)) ##paste0(rep("dodgerblue",9),round(seq(10,75,length.out = 9))))
    aralphas <- c(1,1,
                  seq(1,0.3,length.out = 10)[-1],
                  seq(1,0.3,length.out = 10)[-1])
    connect <- list(c(2:11),c(2,12:21))
    if(relref){
        aralphas <- c(1,
                      seq(1,0.3,length.out = 10)[-1],
                      seq(1,0.3,length.out = 10)[-1])
        connect <- list(c(1:10),c(1,11:20))
    }

    for(j in 1:length(MSElistin)){
        if(evalyearsin[1] == "flex"){
            evalyears <- list(c(2,6),c(3,9),c(6,20))[[j]]
        }else if(evalyearsin[1] == "flex1"){
            evalyears <- list(c(1,6),c(1,9),c(1,20))[[j]]
        }
        MSElist <- MSElistin[[j]]
        ## first column
        for(k in 1:8){
            indi <- k + ((j-1) * 8)
            xaxt <- "n" ## ifelse(k %in% c(7,8), "s", "n")
            yaxt <- "n" ## ifelse(k %in% c(1,3,5,7), "s", "n")
            plotdYdR(MSElist[[k]], Yrs=evalyears, risk = 1, relref = relref,
                     connect=connect, ylim = c(ylimLower[indi],ylimUpper[indi]),
                     xlim=c(xlimLower[indi],xlimUpper[indi]), xaxt=xaxt, yaxt=yaxt,
                     pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
            axis(1, axis.cex=1.2)
            axis(2, axis.cex=1.2)
            plotrix::boxed.labels(par("usr")[1]+0.08*xlimUpper[indi],
                                  par("usr")[4]-0.07*par("usr")[4],
                                  ptitles[k], cex = 0.9,
                                  border = NA,
                                  bg ="white", font=2,
                                  xpad = 1.3, ypad = 1.3)
            box()
            ## text(par("usr")[1]+0.02*xlimUpper[indi],
            ##      par("usr")[4]-0.07*par("usr")[4],
            ##      srt=0, adj = 0, cex=1.6, font=2,
            ##      labels = ptitles[k])
            if(k == 2){
                mtext(speciesAll[j], side =3, font=2, adj=c(1.65,1.65,1.3)[j],
                      line=2, cex=1.2)
            }
        }
        if(j < length(MSElistin)) plot.new()
    }

    if(!relref && !evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
                     side = 1, line = 2.5, outer=TRUE, adj = 0.47, cex=1.2)
    }else if(relref && !evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote("Rel." ~ risk ~ "[years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 2.5, outer=TRUE, adj = 0.44, cex=1.2)
    }else{
        mtext("Abs. risk [%]",
          side = 1, line = 2.5, outer=TRUE, adj = 0.44, cex=1.2)
    }
    if(!evalyearsin[1] %in% c("flex","flex1")){
        mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
              side = 2, line = 2, outer=TRUE, cex=1.2)
    }else{
        mtext("Rel. yield",
              side = 2, line = 2, outer=TRUE, cex=1.2)
    }
    plot.new()
    if(relref){
        legends <- legends[-1]
        argroup <- argroup[-1]
        aralphas <- aralphas[-1]
    }
    legend(0.2,0.73, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=c(1,colorRampPalette(c("gold4","gold1"))(9),
                 colorRampPalette(c("dodgerblue4","dodgerblue1"))(9)), ##rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.32)
    legend(0.2,0.35, title = "MIAVY", title.adj=0,
           legend=c("<0.1","0.1-0.3"," >=0.3"), cex=1.32,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}




plotAllFeb2020 <- function(MSElist, evalyears, legends = NULL,
                             xlimUpper=NULL, xlimLower=NULL, ylimUpper=NULL, ylimLower=NULL,
                    ptitles=rep("",length(MSElist))){
    nscen = length(MSElist)
    if(is.null(xlimUpper)) xlimUpper <- rep(50,10)
    if(is.null(xlimLower)) xlimLower <- rep(0,10)
    if(is.null(ylimUpper)) ylimUpper <- rep(1.05,10)
    if(is.null(ylimLower)) ylimLower <- rep(0,10)
    ly <- layout(matrix(c(1,2,9,3,4,9,5,6,9,7,8,9),
                        byrow = TRUE,ncol=3,nrow=4),
                 widths=c(1,1,0.44,1,1,0.44,1,1,0.44,1,1,0.44),
                 heights=rep(1,5))
    opar <- par(mar=c(3,2,3,2),oma=c(5,5.5,1,0.5))
    legends=legends
    nmp1 <- 7
    nmp2 <- 9
    pchs <- c(3,17,15,16)
    argroup <- c(1,2,
                 rep(3,nmp1),
                 rep(4,nmp2))
    arcols <- c(1,1,
                paste0(rep("grey",nmp1),round(seq(10,75,length.out = nmp1))),
                paste0(rep("grey",nmp2),round(seq(10,75,length.out = nmp2))))
    aralphas <- c(1,1,
                  seq(1,0.3,length.out = nmp1+2)[-1],
                  seq(1,0.3,length.out = nmp2+1)[-1])
    connect <- list(c(2:(nmp1+2)),c(2,(nmp1+3):(nmp1+nmp2)))
    plotdYdR(MSElist[[1]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[1],ylimUpper[1]),
             xlim=c(xlimLower[1],xlimUpper[1]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[1],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[2]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[2],ylimUpper[2]),
             xlim=c(xlimLower[2],xlimUpper[2]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[2],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[3]], Yrs=evalyears, risk = 1,
             connect=connect, ylim = c(ylimLower[3],ylimUpper[3]),
             xlim=c(xlimLower[3],xlimUpper[3]), ## xlim=c(-0.2,0.35), ##
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[3],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[4]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[4],ylimUpper[4]), xlim=c(xlimLower[4],xlimUpper[4]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[4],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[5]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[5],ylimUpper[5]), xlim=c(xlimLower[5],xlimUpper[5]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[5],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[6]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[6],ylimUpper[6]), xlim=c(xlimLower[6],xlimUpper[6]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
##    abline(v=12.380952,lty=3)
    mtext(ptitles[6],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[7]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[7],ylimUpper[7]), xlim=c(xlimLower[7],xlimUpper[7]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=1,col="grey80")
    mtext(ptitles[7],
          font=2,line=0.5,cex=1.3)
    plotdYdR(MSElist[[8]], Yrs=evalyears, risk = 1,
             connect=connect,
             ylim = c(ylimLower[8],ylimUpper[8]), xlim=c(xlimLower[8],xlimUpper[8]),
             pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
##    abline(v=5,lty=1,col="grey80")
##    abline(v=12.380952,lty=3)
    mtext(ptitles[8],
          font=2,line=0.5,cex=1.3)
    ## additional
    ## plotdYdR(MSElist[[9]], Yrs=evalyears, risk = 1,
    ##          connect=connect,
    ##          ylim = c(ylimLower[9],ylimUpper[9]), xlim=c(xlimLower[9],xlimUpper[9]),
    ##          pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
    ## abline(v=5,lty=2)
##     mtext(ptitles[9],
##           font=2,line=0.5,cex=1.3)
##     plotdYdR(MSElist[[10]], Yrs=evalyears, risk = 1,
##              connect=connect,
##              ylim = c(ylimLower[10],ylimUpper[10]), xlim=c(xlimLower[10],xlimUpper[10]),
##              pchs=pchs,argroup=argroup,arcols=arcols,aralphas=aralphas)
##     abline(v=5,lty=2)
## ##   abline(v=12.380952,lty=3)
##     mtext(ptitles[10],
##           font=2,line=0.5,cex=1.3)
    mtext(bquote(Risk ~ "[% years"~.(evalyears[1]) * "-" * .(evalyears[2])*"]"),
          side = 1, line = 1.5, outer=TRUE, adj = 0.38, cex=1.2)
    mtext(bquote("Rel." ~ yield ~ "[years"~.(evalyears[1])*"-"*.(evalyears[2])*"]"),
          side = 2, line = 1.9, outer=TRUE, cex=1.2)
    plot.new()
    legend(-0.05,0.68, legend=legends, title.adj=0,
           pch=pchs[argroup], title = "HCR",
           col=rgb(t(col2rgb("grey10"))/255,alpha=aralphas),
           ncol = 1, bty="n", cex=1.5)
    legend(-0.05,0.32, title = "MIAVY", title.adj=0,
           legend=c("<0.1","0.1-0.3"," >=0.3"), cex=1.5,
           pch=16,col=1,bty="n", pt.cex=c(1.5,1.5+(3-1.5)/2,3))

}





SubTKM <- function (MSEobj, MPs = NULL, sims = NULL, years = NULL){
    checkMSE(MSEobj)
    Class <- class(MPs)
    if (Class == "NULL")
        subMPs <- MSEobj@MPs
    if (Class == "integer" | Class == "numeric")
        subMPs <- MSEobj@MPs[as.integer(MPs)]
    if (Class == "character")
        subMPs <- MPs
    if (Class == "factor")
        subMPs <- as.character(MPs)
    subMPs <- subMPs[!is.na(subMPs)]
    subMPs <- unique(subMPs)
    SubMPs <- match(subMPs, MSEobj@MPs)
    if (any(is.na(SubMPs))) {
        missing <- subMPs[is.na(SubMPs)]
        stop(paste0(missing, collapse = ","), " not found in MSE object",
            call. = FALSE)
    }
    not <- (subMPs %in% MSEobj@MPs)
    ind <- which(not == FALSE)
    newMPs <- MSEobj@MPs[SubMPs]
    if (length(SubMPs) < 1)
        stop("MPs not found")
    if (length(ind > 0)) {
        message("Warning: MPs not found - ", paste0(subMPs[ind],
            " "))
        message("Subsetting by MPs: ", paste0(newMPs, " "))
    }
    ClassSims <- class(sims)
    if (ClassSims == "NULL")
        SubIts <- 1:MSEobj@nsim
    if (ClassSims == "integer" | ClassSims == "numeric") {
        SubIts <- as.integer(sims)
    }
    if (ClassSims == "logical")
        SubIts <- which(sims)
    nsim <- length(SubIts)
    ClassYrs <- class(years)
    AllNYears <- MSEobj@proyears
    if (ClassYrs == "NULL")
        Years <- 1:AllNYears
    if (ClassYrs == "integer" | ClassYrs == "numeric")
        Years <- years
    if (max(Years) > AllNYears)
        stop("years exceeds number of years in MSE")
    if (min(Years) <= 0)
        stop("years must be positive")
    if (min(Years) != 1) {
        message("Not starting from first year. Are you sure you want to do this?")
        message("Probably a bad idea!")
    }
    if (!all(diff(Years) == 1))
        stop("years are not consecutive")
    if (length(Years) <= 1)
        stop("You are going to want more than 1 projection year")
    MSEobj@proyears <- max(Years)
    SubF <- MSEobj@F_FMSY[SubIts, SubMPs, Years, drop = FALSE]
    SubB <- MSEobj@B_BMSY[SubIts, SubMPs, Years, drop = FALSE]
    SubC <- MSEobj@C[SubIts, SubMPs, Years, drop = FALSE]
    SubBa <- MSEobj@B[SubIts, SubMPs, Years, drop = FALSE]
    SubFMa <- MSEobj@FM[SubIts, SubMPs, Years, drop = FALSE]
    SubTACa <- MSEobj@TAC[SubIts, SubMPs, Years, drop = FALSE]
    OutOM <- MSEobj@OM[SubIts, ]
    tt <- try(slot(MSEobj, "Effort"), silent = TRUE)
    if (class(tt) == "try-error")
        slot(MSEobj, "Effort") <- array(NA)
    if (all(is.na(tt)) || all(tt == 0))
        slot(MSEobj, "Effort") <- array(NA)
    if (all(is.na(MSEobj@Effort))) {
        SubEffort <- array(NA)
    }
    else {
        SubEffort <- MSEobj@Effort[SubIts, SubMPs, Years, drop = FALSE]
    }
    tt <- try(slot(MSEobj, "SSB"), silent = TRUE)
    if (class(tt) == "try-error")
        slot(MSEobj, "SSB") <- array(NA)
    if (all(is.na(tt)) || all(tt == 0))
        slot(MSEobj, "SSB") <- array(NA)
    if (all(is.na(MSEobj@SSB))) {
        SubSSB <- array(NA)
    }
    else {
        SubSSB <- MSEobj@SSB[SubIts, SubMPs, Years, drop = FALSE]
    }
    tt <- try(slot(MSEobj, "VB"), silent = TRUE)
    if (class(tt) == "try-error")
        slot(MSEobj, "VB") <- array(NA)
    if (all(is.na(tt)) || all(tt == 0))
        slot(MSEobj, "VB") <- array(NA)
    if (all(is.na(MSEobj@VB))) {
        SubVB <- array(NA)
    }
    else {
        SubVB <- MSEobj@VB[SubIts, SubMPs, Years, drop = FALSE]
    }
    tt <- try(slot(MSEobj, "PAA"), silent = TRUE)
    if (class(tt) == "try-error")
        slot(MSEobj, "PAA") <- array(NA)
    if (all(is.na(tt)) || all(tt == 0))
        slot(MSEobj, "PAA") <- array(NA)
    if (all(is.na(MSEobj@PAA))) {
        SubPAA <- array(NA)
    }
    else {
        SubPAA <- MSEobj@PAA[SubIts, SubMPs, , drop = FALSE]
    }
    tt <- try(slot(MSEobj, "CAL"), silent = TRUE)
    if (class(tt) == "try-error")
        slot(MSEobj, "CAL") <- array(NA)
    if (all(is.na(tt)) || all(tt == 0))
        slot(MSEobj, "CAL") <- array(NA)
    if (all(is.na(MSEobj@CAL))) {
        SubCAL <- array(NA)
    }
    else {
        SubCAL <- MSEobj@CAL[SubIts, SubMPs, , drop = FALSE]
    }
    tt <- try(slot(MSEobj, "CAA"), silent = TRUE)
    if (class(tt) == "try-error")
        slot(MSEobj, "CAA") <- array(NA)
    if (all(is.na(tt)) || all(tt == 0))
        slot(MSEobj, "CAA") <- array(NA)
    if (all(is.na(MSEobj@CAA))) {
        SubCAA <- array(NA)
    }
    else {
        SubCAA <- MSEobj@CAA[SubIts, SubMPs, , drop = FALSE]
    }
    CALbins <- MSEobj@CALbins
    mpind <- match(newMPs, MSEobj@MPs)
    MSEobj@Misc$Unfished$Refs <- MSEobj@Misc$Unfished$Refs[,
        SubIts]
    for (i in 1:length(MSEobj@Misc$Unfished$ByYear)) {
        MSEobj@Misc$Unfished$ByYear[[i]] <- MSEobj@Misc$Unfished$ByYear[[i]][SubIts,
            Years]
    }
    MSEobj@Misc$MSYRefs$Refs <- MSEobj@Misc$MSYRefs$Refs[SubIts,]
    for (i in 1:length(MSEobj@Misc$MSYRefs$ByYear)) {
        MSEobj@Misc$MSYRefs$ByYear[[i]] <- MSEobj@Misc$MSYRefs$ByYear[[i]][SubIts,mpind, Years]
    }
    MSEobj@Misc$TryMP <- MSEobj@Misc$TryMP[mpind]
    MSEobj@Misc$LatEffort <- MSEobj@Misc$LatEffort[SubIts, mpind,
        ]
    MSEobj@Misc$Revenue <- MSEobj@Misc$Revenue[SubIts, mpind,
        ]
    MSEobj@Misc$Cost <- MSEobj@Misc$Cost[SubIts, mpind, ]
    MSEobj@Misc$TAE <- MSEobj@Misc$TAE[SubIts, mpind, ]
    MSEobj@Misc$ErrList$Cbiasa <- MSEobj@Misc$ErrList$Cbiasa[SubIts,
        Years]
    MSEobj@Misc$ErrList$Cerr <- MSEobj@Misc$ErrList$Cerr[SubIts,
        Years]
    MSEobj@Misc$ErrList$Ierr <- MSEobj@Misc$ErrList$Ierr[SubIts,
        Years]
    MSEobj@Misc$ErrList$SpIerr <- MSEobj@Misc$ErrList$SpIerr[SubIts,
        Years]
    MSEobj@Misc$ErrList$VIerr <- MSEobj@Misc$ErrList$VIerr[SubIts,
        Years]
    MSEobj@Misc$ErrList$Recerr <- MSEobj@Misc$ErrList$Recerr[SubIts,
        Years]
    if ("Data" %in% names(MSEobj@Misc)) {
        MSEobj@Misc$Data <- MSEobj@Misc$Data[match(newMPs, MSEobj@MPs)]
    }
    SubResults <- new("MSE", Name = MSEobj@Name, nyears = MSEobj@nyears,
        proyears = MSEobj@proyears, nMPs = length(SubMPs), MPs = newMPs,
        nsim = length(SubIts), OM = OutOM, Obs = MSEobj@Obs[SubIts,
            , drop = FALSE], B_BMSY = SubB, F_FMSY = SubF, B = SubBa,
        SSB = SubSSB, VB = SubVB, FM = SubFMa, SubC, TAC = SubTACa,
        SSB_hist = MSEobj@SSB_hist[SubIts, , , , drop = FALSE],
        CB_hist = MSEobj@CB_hist[SubIts, , , , drop = FALSE],
        FM_hist = MSEobj@FM_hist[SubIts, , , , drop = FALSE],
        Effort = SubEffort, PAA = SubPAA, CAL = SubCAL, CAA = SubCAA,
        CALbins = CALbins, Misc = MSEobj@Misc)
    return(SubResults)
}
