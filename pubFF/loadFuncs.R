## probhcr - additional functions
## July 2020
## Tobias Mildenberger <t.k.mildenberger@gmail.com>

## estimate Blim as Bimass where SP = 50% of MSY
est.blim <- function(dat, set=NULL,
                     msy = NULL,
                     ref = 0.5,
                     fmin = 0.001,
                     fmax = 10,
                     ncores = parallel::detectCores()-1,
                     plot = FALSE){

    ## Checks
    if(is.null(set)) set <- check.set()
    dist <- NULL
    if(!(set$refMethod %in% c("mean","median"))){
        stop("'set$refMethod' not known! Has to be 'mean' or 'median'!")
    }

    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    amax <- dat$amax + 1
    asmax <- amax * ns
    nyref <- set$refYears
    nrep <- set$refN
    nyrefmsy <- set$refYearsMSY
    tvflag <- FALSE

    ## Time-variant processes
    ## natural mortality
    mtv <- length(unique(as.numeric(dat$M)))  ## CHECK: if M is matrix
    ms <- unique(dat$M)
    mind <- match(dat$M, ms)
    if(length(dat$Msel) > 1){
        msel <- dat$Msel[!duplicated(dat$Msel)]
        mseltv <- length(msel)
    }else{
        msel <- dat$Msel[1]
        mseltv <- 1
    }
    if(mseltv > 1 && mseltv != mtv) stop("Both natural mortality over time (dat$M) and over age (dat$Msel) are time-variant, but do not have the same dimensions. This is not yet implemented, please let both vary equally or keep one of them constant.")
    alltv <- max(c(mtv, mseltv))
    ## selectivity
    if(length(dat$sel) > 1){
        sel <- dat$sel[!duplicated(dat$sel)]
        seltv <- length(sel)
    }else{
        sel <- dat$sel[1]
        seltv <- 1
    }
    if(seltv > 1 && alltv > 1 && seltv != alltv) stop("Both gear selectivity (dat$sel) and natural mortality (dat$M or dat$Msel) are time-variant, but do not have the same dimensions. This is not yet implemented, please let both vary equally or keep one of them constant.")
    alltv <- max(c(alltv,seltv))

    ##
    ## errors (have to be re-used for estimation of Bmsy)
    errs <- vector("list", nrep)
    for(i in 1:nrep){
        errs[[i]] <- vector("list", 9)
        errs[[i]]$eF <- gen.noise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
##        errs[[i]]$eF <- rep(1.0, nyref)
        errs[[i]]$eR <- gen.noise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        errs[[i]]$eM <- gen.noise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        errs[[i]]$eH <- gen.noise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        errs[[i]]$eW <- gen.noise(nyref, set$noiseW[1], set$noiseW[2], set$noiseW[3])
        errs[[i]]$eR0 <- gen.noise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        errs[[i]]$eMat <- gen.noise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        errs[[i]]$eSel <- gen.noise(nyref, set$noiseSel[1], set$noiseSel[2], set$noiseSel[3])
        errs[[i]]$eImp <- gen.noise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
    }
    ##
    datx <- dat

    ##
    ## datx$yvec <- rep(1:nyref, each = ns)
    ## datx$svec <- rep(1:ns, each = nyref)
    ## datx$s1vec <- seq(1, nyref * ns, ns)
    ## datx$as2a <- rep(1:amax, each = ns)
    ## datx$as2s <- rep(1:ns, amax)
    ## datx$inds <- seq(1,asmax,ns)

    ## browser()

    ## setx$refMethod <- "median"
    ## setx$refYearsMSY <- 10
    ## f <- 1
    ## x = 1
    ## setx <- c(set, errs[[x]])
    ## (ref*msy - unlist(simpop(log(f), datx, setx, out=1)))^2


    ## Fmsy
    res <- parallel::mclapply(as.list(1:nrep), function(x){
        setx <- c(set, errs[[x]])
        tmp <- rep(NA, alltv)
        for(i in 1:alltv){
            ## datx$M <- as.matrix(ms[mtv[i],])
            ## datx$Msel <- msels[mseltv[i]]
            ## datx$sel <- sels[seltv[i]]
            setx$tvm <- 1
            setx$tvmsel <- 1
            setx$tvsel <- 1
            ## datx$M <- rep(dat$M[i], nyref)
            ## ind <- (i-1)*ns+1
            ## datx$Ms <- rep(dat$Ms[ind:(ind+ns)], nyref)
            ## ind2 <- ifelse(ntv2 > 1, i, 1)
            ## datx$Msels <- msels[ind2]
            ## datx$Msel <- lapply(datx$Msels, rowMeans)
            opt <- optimise(function(x) (ref * msy - unlist(simpop(x, datx, setx, out=1)))^2,
                            log(c(fmin,fmax)), maximum = FALSE)
            tmp[i] <- exp(opt$minimum)
        }
        return(tmp)
    }, mc.cores = ncores)
    fmsys <- do.call(rbind, res)


    ## MSY and Biomass reference points
    res <- parallel::mclapply(as.list(1:nrep), function(x){
        setx <- c(set, errs[[x]])
        tmp <- vector("list", alltv)
        for(i in 1:alltv){
            ## datx$M <- rep(dat$M[i], nyref)
            ## ind <- (i-1)*ns+1
            ## datx$Ms <- rep(dat$Ms[ind:(ind+ns)], nyref)
            ## ind2 <- ifelse(ntv2 > 1, i, 1)
            ## datx$Msels <- msels[ind2]
            ## datx$Msel <- lapply(datx$Msels, rowMeans)
            setx$tvm <- 1
            setx$tvmsel <- 1
            setx$tvsel <- 1
            tmp0 <- simpop(log(fmsys[x,i]), datx, setx, out=0)
            if(set$refMethod == "mean"){
                tmp[[i]] <- c(mean(tail(tmp0$CW,nyrefmsy)), mean(tail(tmp0$TSB,nyrefmsy)),
                              mean(tail(tmp0$ESB,nyrefmsy)), mean(tail(tmp0$SSB,nyrefmsy)))
            }else if(set$refMethod == "median"){
                tmp[[i]] <- c(median(tail(tmp0$CW,nyrefmsy)), median(tail(tmp0$TSB,nyrefmsy)),
                              median(tail(tmp0$ESB,nyrefmsy)), median(tail(tmp0$SSB,nyrefmsy)))
            }
        }
        return(tmp)
    }, mc.cores = ncores)

    ## sort in list by reference point with matrix (nrep, ntv)
    brefs <- vector("list", 4) ## c("MSY","Bmsy","ESBmsy","SSBmsy")
    for(i in 1:4){
        brefs[[i]] <- do.call(rbind,
                              lapply(as.list(1:nrep),
                                     function(x) sapply(1:alltv, function(j) res[[x]][[j]][[i]])))
    }


    ## all refs in one list
    refs <- c(list(fmsys), brefs)

    ## remove runs where long-term SP is smaller or equal to 0
    for(i in 1:alltv){
        ind <- which(brefs[[1]][,i] <= 0)
        if(length(ind) > 0){
            for(j in 1:6) refs[[j]][,i] <- refs[[j]][-ind,i]
        }
    }

    ## overall refs
    if(set$refMethod == "mean"){
        meds <- lapply(refs, function(x){
            tmp <- rep(NA, alltv)
            for(i in 1:alltv) tmp[i] <- mean(x[,i])
            return(tmp)
        })
    }else if(set$refMethod == "median"){
        meds <- lapply(refs, function(x){
            tmp <- rep(NA, alltv)
            for(i in 1:alltv) tmp[i] <- median(x[,i])
            return(tmp)
        })
    }


    refdist <- refs
    names(refdist) <- c("Fmsy","MSY", "Bmsy", "ESBmsy", "SSBmsy")
    dat$refdist <- refdist

    refmed <- matrix(NA, ncol=length(refs), nrow=nt)
    for(i in 1:5){
        refmed[,i] <- meds[[i]][mind]
    }
    colnames(refmed) <- c("Fmsy","MSY", "Bmsy", "ESBmsy", "SSBmsy")
    dat$ref <- as.data.frame(refmed)

    if(plot){
        if(alltv < 4){
            cols <- c("dodgerblue2","darkorange","darkgreen","purple")
            nr <- floor(length(refdist)/2)
            par(mfrow=c(nr,2))
            for(i in 1:length(refdist)){
                hist(refdist[[i]][,1], main = names(refdist)[i],
                     breaks=20, freq = TRUE, xlim = range(refdist[[i]]),
                     xlab = "", col = rgb(t(col2rgb(cols[1]))/255,alpha=0.4))
                if(alltv > 1){
                    for(j in 2:alltv) hist(refdist[[i]][,j],
                                         breaks=20, freq = TRUE,
                                         add = TRUE, col = rgb(t(col2rgb(cols[j]))/255,alpha=0.4))
                }
                ## abline(v=mean(dist[,i]), lty=1, lwd=1.5, col=4)
                ## abline(v=median(dist[,i]), lty=2, lwd=1.5, col=4)
                ## if(i == 1) legend("topright", legend = c("mean","median"),
                ##                   col=4, lty=c(1,2),lwd=1.5)
            }
        }else{
            nr <- floor(length(refdist)/2)
            par(mfrow=c(nr,2))
            for(i in 1:length(refdist)){
                plot(refmed[,i], main = colnames(refmed)[i],
                     ty = 'l', lwd=1.5, xlab = "Time", ylab = colnames(refmed)[i])
            }
        }
    }

    ## return
    return(dat)
}



## colours
colfuncBlue <- colorRampPalette(c("#7BAFF3", "#2C65AF"))
colfuncOrange <- colorRampPalette(c("#E8DC3D", "#DFB529"))
colfuncGreen <- colorRampPalette(c("#78D51A", "#529A0A"))
colfuncPurple <- colorRampPalette(c("#E991F6", "#810794"))  ## old light: D559E8  ## old dark: 9E09B5
colfuncGrey <- colorRampPalette(c("#BAB4B4","#504D4D"))
colfuncRed <- colorRampPalette(c("#FC6666", "#9F1212"))
colfuncRed <- colorRampPalette(c("#EF9631", "#D07309"))
colfuncTurquoise <- colorRampPalette(c("#4BE3EE","#0E8F98"))


## lables
change.hcrs.labels <- function(hcrs){

    plyr::revalue(hcrs,
                  c("MSY" = expression("Median"),

                    "MSY-C45" = expression(italic(f)^{C}*"=0.45"),
                    "MSY-C35" = expression(italic(f)^{C}*"=0.35"),
                    "MSY-C25" = expression(italic(f)^{C}*"=0.25"),
                    "MSY-C15" = expression(italic(f)^{C}*"=0.15"),
                    "MSY-C05" = expression(italic(f)^{C}*"=0.05"),
                    "MSY-C01" = expression(italic(f)^{C}*"=0.01"),

                    "BT50" = expression(italic(B)[T]*"=0.5"),
                    "BT100" = expression(italic(B)[T]*"=1"),
                    "BT200" = expression(italic(B)[T]*"=2"),
                    "BT300" = expression(italic(B)[T]*"=3"),
                    "BT400" = expression(italic(B)[T]*"=4"),

                    "BT50-BL30" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5"),
                    "BT100-BL30" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=1"),
                    "BT200-BL30" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=2"),
                    "BT300-BL30" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=3"),
                    "BT400-BL30" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=4"),

                    "BT100-BL50" = expression(italic(B)[L]*"=0.5,"*italic(B)[T]*"=1"),
                    "BT200-BL50" = expression(italic(B)[L]*"=0.5,"*italic(B)[T]*"=2"),
                    "BT300-BL50" = expression(italic(B)[L]*"=0.5,"*italic(B)[T]*"=3"),
                    "BT400-BL50" = expression(italic(B)[L]*"=0.5,"*italic(B)[T]*"=4"),

                    "BT50-BL30-C45" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{C}*"=0.45"),
                    "BT50-BL30-C35" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{C}*"=0.35"),
                    "BT50-BL30-C25" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{C}*"=0.25"),
                    "BT50-BL30-C15" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{C}*"=0.15"),

                    "BT50-BL30-FB45" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"B,F"}*"=0.45"),
                    "BT50-BL30-FB35" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"B,F"}*"=0.35"),
                    "BT50-BL30-FB25" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"B,F"}*"=0.25"),
                    "BT50-BL30-FB15" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"B,F"}*"=0.15"),

                    "BT50-BL30-CFB45" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"C,B,F"}*"=0.45"),
                    "BT50-BL30-CFB35" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"C,B,F"}*"=0.35"),
                    "BT50-BL30-CFB25" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"C,B,F"}*"=0.25"),
                    "BT50-BL30-CFB15" =
                        expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=0.5,"*italic(f)^{"C,B,F"}*"=0.15"),

                    "BT100-BL30-C45" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=1,"*italic(f)^{C}*"=0.45"),
                    "BT100-BL30-C35" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=1,"*italic(f)^{C}*"=0.35"),
                    "BT100-BL30-C25" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=1,"*italic(f)^{C}*"=0.25"),
                    "BT100-BL30-C15" = expression(italic(B)[L]*"=0.3,"*italic(B)[T]*"=1,"*italic(f)^{C}*"=0.15"),

                    "BT50-BL50" = expression(italic(B)[L]*"="*italic(B)[T]*"=0.5"),
                    "BT100-BL100" = expression(italic(B)[L]*"="*italic(B)[T]*"=1")


                    )
                  )

}


change.hcrs.tex <- function(hcrs){

    plyr::revalue(hcrs,
                  c("MSY" = "Median",

                    "MSY-C45" = "$f^{\\text{C}}=0.45$",
                    "MSY-C35" = "$f^{\\text{C}}=0.35$",
                    "MSY-C25" = "$f^{\\text{C}}=0.25$",
                    "MSY-C15" = "$f^{\\text{C}}=0.15$",
                    "MSY-C05" = "$f^{\\text{C}}=0.05$",
                    "MSY-C01" = "$f^{\\text{C}}=0.01$",

                    "MSY-FB45" = "$f^{\\text{B,F}}=0.45$",
                    "MSY-FB35" = "$f^{\\text{B,F}}=0.35$",
                    "MSY-FB25" = "$f^{\\text{B,F}}=0.25$",
                    "MSY-FB15" = "$f^{\\text{B,F}}=0.15$",
                    "MSY-FB05" = "$f^{\\text{B,F}}=0.05$",
                    "MSY-FB01" = "$f^{\\text{B,F}}=0.01$",

                    "MSY-CFB45" = "$f^{\\text{C,B,F}}=0.45$",
                    "MSY-CFB35" = "$f^{\\text{C,B,F}}=0.35$",
                    "MSY-CFB25" = "$f^{\\text{C,B,F}}=0.25$",
                    "MSY-CFB15" = "$f^{\\text{C,B,F}}=0.15$",
                    "MSY-CFB05" = "$f^{\\text{C,B,F}}=0.05$",
                    "MSY-CFB01" = "$f^{\\text{C,B,F}}=0.01$",

                    "BT50" = "$\\text{B}_{\\text{T}}=0.5$",
                    "BT100" = "$\\text{B}_{\\text{T}}=1$",
                    "BT200" = "$\\text{B}_{\\text{T}}=2$",
                    "BT300" = "$\\text{B}_{\\text{T}}=3$",
                    "BT400" = "$\\text{B}_{\\text{T}}=4$",

                    "BT50-BL30" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5$",
                    "BT100-BL30" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1$",
                    "BT200-BL30" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=2$",
                    "BT300-BL30" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=3$",
                    "BT400-BL30" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=4$",

                    "BT50-BL50" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5$",
                    "BT100-BL50" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=1$",
                    "BT200-BL50" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=2$",
                    "BT300-BL50" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=3$",
                    "BT400-BL50" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=4$",

                    "BT90-BL10" = "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9$",
                    "BT80-BL20" = "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8$",

                    "BT50-BL50" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5$",
                    "BT100-BL100" = "$\\text{B}_{\\text{L}}=1, \\text{B}_{\\text{T}}=1$",

                    "BT50-C45" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.45$",
                    "BT50-C35" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.35$",
                    "BT50-C25" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.25$",
                    "BT50-C15" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.15$",
                    "BT50-C05" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.05$",

                    "BT50-FB45" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.45$",
                    "BT50-FB35" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.35$",
                    "BT50-FB25" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.25$",
                    "BT50-FB15" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.15$",
                    "BT50-FB05" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.05$",

                    "BT50-CFB45" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.45$",
                    "BT50-CFB35" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.35$",
                    "BT50-CFB25" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.25$",
                    "BT50-CFB15" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.15$",
                    "BT50-CFB05" = "$\\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.05$",

                    "BT100-C45" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C}}=0.45$",
                    "BT100-C35" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C}}=0.35$",
                    "BT100-C25" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C}}=0.25$",
                    "BT100-C15" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C}}=0.15$",
                    "BT100-C05" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C}}=0.05$",

                    "BT100-FB45" = "$\\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.45$",
                    "BT100-FB35" = "$\\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.35$",
                    "BT100-FB25" = "$\\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.25$",
                    "BT100-FB15" = "$\\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.15$",
                    "BT100-FB05" = "$\\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.05$",

                    "BT100-CFB45" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.45$",
                    "BT100-CFB35" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.35$",
                    "BT100-CFB25" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.25$",
                    "BT100-CFB15" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.15$",
                    "BT100-CFB05" = "$\\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.05$",

                    "BT50-BL30-C45" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.45$",
                    "BT50-BL30-C35" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.35$",
                    "BT50-BL30-C25" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.25$",
                    "BT50-BL30-C15" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.15$",
                    "BT50-BL30-C05" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.05$",

                    "BT50-BL30-FB45" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.45$",
                    "BT50-BL30-FB35" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.35$",
                    "BT50-BL30-FB25" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.25$",
                    "BT50-BL30-FB15" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.15$",
                    "BT50-BL30-FB05" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.05$",

                    "BT50-BL30-CFB45" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.45$",
                    "BT50-BL30-CFB35" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.35$",
                    "BT50-BL30-CFB25" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.25$",
                    "BT50-BL30-CFB15" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.15$",
                    "BT50-BL30-CFB05" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.05$",

                    "BT100-BL30-C45" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.45$",
                    "BT100-BL30-C35" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.35$",
                    "BT100-BL30-C25" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.25$",
                    "BT100-BL30-C15" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.15$",
                    "BT100-BL30-C05" = "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.05$",

                    "BT100-BL30-FB45" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.45$",
                    "BT100-BL30-FB35" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.35$",
                    "BT100-BL30-FB25" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.25$",
                    "BT100-BL30-FB15" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.15$",
                    "BT100-BL30-FB05" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{B,F}}=0.05$",

                    "BT100-BL30-CFB45" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.45$",
                    "BT100-BL30-CFB35" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.35$",
                    "BT100-BL30-CFB25" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.25$",
                    "BT100-BL30-CFB15" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.15$",
                    "BT100-BL30-CFB05" =
                        "$\\text{B}_{\\text{L}}=0.3, \\text{B}_{\\text{T}}=1, f^{\\text{C,B,F}}=0.05$",

                    "BT80-BL20-C45" = "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C}}=0.45$",
                    "BT80-BL20-C35" = "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C}}=0.35$",
                    "BT80-BL20-C25" = "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C}}=0.25$",
                    "BT80-BL20-C15" = "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C}}=0.15$",
                    "BT80-BL20-C05" = "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C}}=0.05$",

                    "BT80-BL20-FB45" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{B,F}}=0.45$",
                    "BT80-BL20-FB35" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{B,F}}=0.35$",
                    "BT80-BL20-FB25" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{B,F}}=0.25$",
                    "BT80-BL20-FB15" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{B,F}}=0.15$",
                    "BT80-BL20-FB05" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{B,F}}=0.05$",

                    "BT80-BL20-CFB45" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C,B,F}}=0.45$",
                    "BT80-BL20-CFB35" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C,B,F}}=0.35$",
                    "BT80-BL20-CFB25" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C,B,F}}=0.25$",
                    "BT80-BL20-CFB15" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C,B,F}}=0.15$",
                    "BT80-BL20-CFB05" =
                        "$\\text{B}_{\\text{L}}=0.2, \\text{B}_{\\text{T}}=0.8, f^{\\text{C,B,F}}=0.05$",

                    "BT90-BL10-C45" = "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C}}=0.45$",
                    "BT90-BL10-C35" = "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C}}=0.35$",
                    "BT90-BL10-C25" = "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C}}=0.15$",
                    "BT90-BL10-C15" = "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C}}=0.15$",
                    "BT90-BL10-C05" = "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C}}=0.05$",

                    "BT90-BL10-FB45" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{B,F}}=0.45$",
                    "BT90-BL10-FB35" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{B,F}}=0.35$",
                    "BT90-BL10-FB25" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{B,F}}=0.15$",
                    "BT90-BL10-FB15" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{B,F}}=0.15$",
                    "BT90-BL10-FB05" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{B,F}}=0.05$",

                    "BT90-BL10-CFB45" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C,B,F}}=0.45$",
                    "BT90-BL10-CFB35" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C,B,F}}=0.35$",
                    "BT90-BL10-CFB25" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C,B,F}}=0.15$",
                    "BT90-BL10-CFB15" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C,B,F}}=0.15$",
                    "BT90-BL10-CFB05" =
                        "$\\text{B}_{\\text{L}}=0.1, \\text{B}_{\\text{T}}=0.9, f^{\\text{C,B,F}}=0.05$",

                    "BT50-BL50-C45" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.45$",
                    "BT50-BL50-C35" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.35$",
                    "BT50-BL50-C25" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.25$",
                    "BT50-BL50-C15" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.15$",
                    "BT50-BL50-C05" = "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C}}=0.05$",

                    "BT50-BL50-FB45" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.45$",
                    "BT50-BL50-FB35" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.35$",
                    "BT50-BL50-FB25" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.25$",
                    "BT50-BL50-FB15" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.15$",
                    "BT50-BL50-FB05" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{B,F}}=0.05$",

                    "BT50-BL50-CFB45" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.45$",
                    "BT50-BL50-CFB35" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.35$",
                    "BT50-BL50-CFB25" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.25$",
                    "BT50-BL50-CFB15" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.15$",
                    "BT50-BL50-CFB05" =
                        "$\\text{B}_{\\text{L}}=0.5, \\text{B}_{\\text{T}}=0.5, f^{\\text{C,B,F}}=0.05$"




                    )
                  )

}




arrowLine <- function(x,y,...){
    lines(x,y,...)
    for(i in 2:length(x)){
        Ax=seq(x[i-1],x[i],length=3)
        Ay=seq(y[i-1],y[i],length=3)
        arrows(Ax[1],Ay[1],Ax[2],Ay[2],length = 0.1, ...)
    }
}


arrowLine2 <- function(x0,y0,x1,y1,nArrow=1,...)
{
  lines(c(x0,x1),c(y0,y1),...)
  Ax=seq(x0,x1,length=nArrow+1)
    Ay=seq(y0,y1,length=nArrow+1)
  for (i in 1:nArrow)
  {
    arrows(Ax[i],Ay[i],Ax[i+1],Ay[i+1],...)
  }
}


estCI.ratioMeans <- function(ys1, ys2, largeSample =FALSE){

    ## mean: mu
    mu1 <- mean(ys1)
    mu2 <- mean(ys2)

    ## var
    var1 <- var(ys1)
    var2 <- var(ys2)

    ## var
    n1 <- length(ys1)
    n2 <- length(ys2)

    ## ratio
    t <- mu1/mu2

    ## variance
    var <- var1/(n1 * mu1^2) + var2/(n2 * mu2^2)

    ## df
    df <- (var1/(mu1^2*n1) + var2/(mu2^2*n2))^2 / (var1^2/(mu1^4*n1^2*(n1-1)) + var2^2/(mu2^4*n2^2*(n2-1)))

    ## ci
    lo <- exp(log(t) - qt(1-0.025, df) * sqrt(var))
    up <- exp(log(t) + qt(1-0.025, df) * sqrt(var))

    return(c(lo,t,up))
}



estCI.ratioMedians <- function(ys1, ys2, largeSample =FALSE){

    # median: eta
    eta1 <- median(ys1)
    eta2 <- median(ys2)

    ## ratio
    t <- eta1/eta2

    # log median: etaStar
    etaStar1 <- log(eta1)
    etaStar2 <- log(eta2)

    # sample size: n
    n1 <- length(ys1)
    n2 <- length(ys2)

    ## variable: c (nearest nonzero integer)
    c1 <- round((n1+1)/2 - sqrt(n1))
    c1 <- ifelse(c1 == 0, 1, c1)
    c2 <- round((n2+1)/2 - sqrt(n2))
    c2 <- ifelse(c2 == 0, 1, c2)

    ## variable: p
    p1 <- 2*0.5^n1 * sum(sapply(0:(c1-1), function(i) factorial(n1)/(factorial(i) * factorial(n1-i))))
    p2 <- 2*0.5^n2 * sum(sapply(0:(c2-1), function(i) factorial(n2)/(factorial(i) * factorial(n2-i))))

    ## variable: z
    z1 <- qnorm(1-p1/2)
    z2 <- qnorm(1-p2/2)

    ## set z = 2 if n large...
    if(n1 > 1e3 || largeSample) z1 <- 2
    if(n2 > 1e3 || largeSample) z2 <- 2

    ## variance
    var1 <- ((log(ys1[n1-c1+1]) - log(ys1[c1])) / (2*z1))^2
    var2 <- ((log(ys2[n1-c1+1]) - log(ys2[c1])) / (2*z2))^2

    ## ci
    lo <- (eta1/eta2) * exp(-qnorm(1-0.025) * sqrt(var1 + var2))
    up <- (eta1/eta2) * exp(qnorm(1-0.025) * sqrt(var1 + var2))

    return(c(lo,t,up))
}


estCI.subMeans <- function(ys1, ys2, largeSample =FALSE){

    ## mean: mu
    mu1 <- mean(ys1)
    mu2 <- mean(ys2)

    ## var
    var1 <- var(ys1)
    var2 <- var(ys2)

    ## var
    n1 <- length(ys1)
    n2 <- length(ys2)

    ## ratio
    t <- mu1 - mu2

    ## variance
    var <- var1/n1 + var2/n2

    ## ci
    lo <- t - qnorm(1-0.025) * sqrt(var)
    up <- t + qnorm(1-0.025) * sqrt(var)

    return(c(lo,t,up))
}




## y1=y2 vectors with TRUE,FALSE
estCI.ratioRisks <- function(y1, y2){

    x <- sum(y1)
    m <- length(y1)
    y <- sum(y2)
    n <- length(y2)

    ## ratio
    t <- (x/m) / (y/n)

    ## variance
    var <- 1/x - 1/m + 1/y - 1/n

    ## ci
    lo <- t * exp(-qnorm(1-0.025) * sqrt(var))
    up <- t * exp(qnorm(1-0.025) * sqrt(var))

    return(c(lo,t,up))

}


## y1=y2 vectors with TRUE,FALSE
estCI.subRisks <- function(y1, y2){

    x <- sum(y1)
    m <- length(y1)
    y <- sum(y2)
    n <- length(y2)
    pt1 <- prop.test(x, n = m, conf.level = 0.95, correct = FALSE)
    pt2 <- prop.test(y, n = n, conf.level = 0.95, correct = FALSE)
    p1 <- pt1$estimate ## or x/m
    p2 <- pt2$estimate ## or y/n

    ## risk difference
    t <- p1 - p2

    ## variables
    l1 <- pt1$conf.int[1]
    u1 <- pt1$conf.int[2]
    l2 <- pt2$conf.int[1]
    u2 <- pt2$conf.int[2]

    ## variance
    var1 <- (p1 - l1)^2 + (u2 - p2)^2
    var2 <- (p2 - l2)^2 + (u1 - p1)^2

    ## ci
    lo <- t - qnorm(1-0.025) * sqrt(var1)
    up <- t + qnorm(1-0.025) * sqrt(var2)

    return(c(lo,t,up))

}

estMetsDiff <- function(x, dat, xRefFmsy, proy = 41:80){

    nrep <- length(x[[1]]) ## should be same as x[[2]]

    if(!all(names(x[[1]]) == names(x[[2]]))) stop("Something went wrong with extracting converged reps for scenario and reference scenario. Replicate names do not match... (Check code for getting mselist2)")

    ## TODO: this should be SSBfinal
    ## riskR <- estCI.subRisks(unlist(lapply(x[[2]], function(x) x$SSB[proy,ncol(x$SSB)]/dat$ref$SSBlim < 1)),
    ##                         unlist(lapply(x[[1]], function(x) x$SSB[proy,ncol(x$SSB)]/dat$ref$SSBlim < 1)))

    riskR <- estCI.subRisks(unlist(lapply(x[[2]], function(x) x$TSBfinal[proy]/dat$ref$Blim < 1)),
                            unlist(lapply(x[[1]], function(x) x$TSBfinal[proy]/dat$ref$Blim < 1)))



    cmsyR <- estCI.subMeans(unlist(lapply(1:nrep,
                                              function(i){
                                                  indi <- which(names(xRefFmsy[[2]]) == names(x[[2]])[i])
                                                  refyield <- apply(xRefFmsy[[2]][[indi]]$CW,1,sum)
                                                  median(apply(x[[2]][[i]]$CW,1,sum)[proy]/refyield[proy])
                                              })),
                                unlist(lapply(1:nrep,
                                              function(i){
                                                  indi <- which(names(xRefFmsy[[1]]) == names(x[[1]])[i])
                                                  refyield <- apply(xRefFmsy[[1]][[indi]]$CW,1,sum)
                                                  median(apply(x[[1]][[i]]$CW,1,sum)[proy]/refyield[proy])
                                              })))

    res <- rbind(riskR, cmsyR)
    colnames(res) <- c("lo","mu","up")
    rownames(res) <- c("risk","cmsy")
    return(res)
}




getConvsScens <- function(mse, convyears = "all", convhcrs = "all",
                          convscens = "all", verbose = FALSE, forMetDiffs=FALSE){

    nysim <- nrow(mse[[1]][[1]][[1]]$tacs)
    dims <- dim(mse[[1]][[1]][[1]]$TSB)
    ny <- dims[1] - nysim
    ns <- dims[2]

    if(!is.na(convyears[1]) && convyears[1] == "all") convyears <- 1:ny
    if(!is.na(convhcrs[1]) && convhcrs[1] == "all") convhcrs <- 1:length(mse[[1]])
    if(!is.na(convscens[1]) && convscens[1] == "all") convscens <- 1:length(mse)

    if(!is.numeric(convscens)){
        nscen <- length(mse[convscens])
    }else if(is.numeric(convscens)){
        nscen <- length(mse)
    }
    nhcr <- length(mse[[1]])
    hcrs <- names(mse[[1]])
    nrep <- length(mse[[1]][[1]])

    if(!forMetDiffs){
        res <- vector("list",nscen)
        for(hcr in 1:nhcr){
            nami <- names(mse[[1]][[hcr]])
            for(scen in 2:nscen){
                nami <- nami[which(nami %in% names(mse[[convscens[scen]]][[hcr]]))]
            }
            for(scen in 1:nscen){
                if(hcr == 1) res[[scen]] <- vector("list",nhcr)
                inds <- which(names(mse[[convscens[scen]]][[hcr]]) %in% nami)
                if(hcrs[hcr] != "refFmsy"){
                    res[[scen]][[hcr]] <- mse[[convscens[scen]]][[hcr]][inds]
                }else res[[scen]][[hcr]] <- mse[[convscens[scen]]][[hcr]]
                names(res[[scen]]) <- hcrs
            }
        }
    }else{
    if(is.numeric(convscens)){
        if(length(convscens) == 1){
            res <- vector("list",nhcr)
            for(hcr in 1:nhcr){
                namiRef <- names(mse[[convscens]][[hcr]])
                res[[hcr]] <- vector("list",nscen)   ## convscens has to have length=1
                for(scen in 1:nscen){
                    res[[hcr]][[scen]] <- vector("list",2)
                    nami <- namiRef[which(namiRef %in% names(mse[[scen]][[hcr]]))]
                    inds <- which(names(mse[[convscens]][[hcr]]) %in% nami)
                    res[[hcr]][[scen]][[1]] <- mse[[convscens]][[hcr]][inds]
                    inds <- which(names(mse[[scen]][[hcr]]) %in% nami)
                    res[[hcr]][[scen]][[2]] <- mse[[scen]][[hcr]][inds]
                }
            }
        }else{
            res <- vector("list",nhcr)
            for(hcr in 1:nhcr){
                nami <- names(mse[[1]][[hcr]])
                for(scen in 2:nscen){
                    nami <- nami[which(nami %in% names(mse[[convscens[scen]]][[hcr]]))]
                }
                res[[hcr]] <- vector("list",nscen)   ## convscens has to have length=1
                for(scen in 1:nscen){
                    res[[hcr]][[scen]] <- vector("list",2)
                    if(hcrs[hcr] != "refFmsy"){
                        inds <- which(names(mse[[1]][[hcr]]) %in% nami)
                    }else inds <- 1:length(mse[[1]][[hcr]])
                    res[[hcr]][[scen]][[1]] <- mse[[1]][[hcr]][inds]
                    if(hcrs[hcr] != "refFmsy"){
                        inds <- which(names(mse[[convscens[scen]]][[hcr]]) %in% nami)
                    }else inds <- 1:length(mse[[convscens[scen]]][[hcr]])
                    res[[hcr]][[scen]][[2]] <- mse[[convscens[scen]]][[hcr]][inds]
                }
            }
        }
    }
    }

    ## return
    return(res)

}


subsetMSE <- function(mse, hcrs = "all"){

    nspec <- length(mse)
    nscen <- length(mse[[1]])
    if(hcrs[1] == "all") hcrs <- 1:length(mse[[1]][[1]])
    nhcr <- length(hcrs)

    res <- vector("list",nspec)
    for(spec in 1:nspec){
        res[[spec]] <- vector("list",nscen)
        for(scen in 1:nscen){
            res[[spec]][[scen]] <- mse[[spec]][[scen]][hcrs]
        }
    }

    ## return
    return(res)

}



## rename HCRs (For tables)
renameHCRs <- function(x){

    plyr::revalue(x,
                  c("MSY" = "$f_{0.5}$",
                      "HS25" = "$\\text{BT}_{0.25}$",
              "HS50" = "$\\text{BT}_{0.50}$",
              "HS75" = "$\\text{BT}_{0.75}$",
              "HS100" = "$\\text{BT}_{1.00}$",
              "HS-100" = "$\\text{BT}_{1.00}$",
              "HS125" = "$\\text{BT}_{1.25}$",
              "HS150" = "$\\text{BT}_{1.50}$",
              "HS175" = "$\\text{BT}_{1.75}$",
              "HS200" = "$\\text{BT}_{2.00}$",
              "HS300" = "$\\text{BT}_{3.00}$",
              "MSY-C05" = "$f^{\\text{C}}_{0.005}$",
              "MSY-C5" = "$f^{\\text{C}}_{0.05}$",
              "MSY-C15" = "$f^{\\text{C}}_{0.15}$",
              "MSY-C25" = "$f^{\\text{C}}_{0.25}$",
              "MSY-C35" = "$f^{\\text{C}}_{0.35}$",
              "MSY-C45" = "$f^{\\text{C}}_{0.45}$",
              "MSY-F05" = "$f^{\\text{F}}_{0.005}$",
              "MSY-F5" = "$f^{\\text{F}}_{0.05}$",
              "MSY-F15" = "$f^{\\text{F}}_{0.15}$",
              "MSY-F25" = "$f^{\\text{F}}_{0.25}$",
              "MSY-F35" = "$f^{\\text{F}}_{0.35}$",
              "MSY-F45" = "$f^{\\text{F}}_{0.45}$",
              "MSY-A05" = "$f^{\\text{CF}}_{0.005}$",
              "MSY-A5" = "$f^{\\text{CF}}_{0.05}$",
              "MSY-A15" = "$f^{\\text{CF}}_{0.15}$",
              "MSY-A25" = "$f^{\\text{CF}}_{0.25}$",
              "MSY-A35" = "$f^{\\text{CF}}_{0.35}$",
              "MSY-A45" = "$f^{\\text{CF}}_{0.45}$",
              "HS50-C05" = "$\\text{BT}_{50}\\text{-f}^{\\text{C}}_{0.005}$",
              "HS50-C5" = "$\\text{BT}_{50}f^{\\text{C}}_{0.05}$",
              "HS50-C15" = "$\\text{BT}_{50}f^{\\text{C}}_{0.15}$",
              "HS50-C25" = "$\\text{BT}_{50}f^{\\text{C}}_{0.25}$",
              "HS50-C35" = "$\\text{BT}_{50}f^{\\text{C}}_{0.35}$",
              "HS50-C45" = "$\\text{BT}_{50}f^{\\text{C}}_{0.45}$",
              "HS50-F05" = "$\\text{BT}_{50}f^{\\text{F}}_{0.005}$",
              "HS50-F5" = "$\\text{BT}_{50}f^{\\text{F}}_{0.05}$",
              "HS50-F15" = "$\\text{BT}_{50}f^{\\text{F}}_{0.15}$",
              "HS50-F25" = "$\\text{BT}_{50}f^{\\text{F}}_{0.25}$",
              "HS50-F35" = "$\\text{BT}_{50}f^{\\text{F}}_{0.35}$",
              "HS50-F45" = "$\\text{BT}_{50}f^{\\text{F}}_{0.45}$",
              "HS50-B05" = "$\\text{BT}_{50}f^{\\text{B}}_{0.005}$",
              "HS50-B5" = "$\\text{BT}_{50}f^{\\text{B}}_{0.05}$",
              "HS50-B15" = "$\\text{BT}_{50}f^{\\text{B}}_{0.15}$",
              "HS50-B25" = "$\\text{BT}_{50}f^{\\text{B}}_{0.25}$",
              "HS50-B35" = "$\\text{BT}_{50}f^{\\text{B}}_{0.35}$",
              "HS50-B45" = "$\\text{BT}_{50}f^{\\text{B}}_{0.45}$",
              "HS50-FB05" = "$\\text{BT}_{50}f^{\\text{FB}}_{0.005}$",
              "HS50-FB5" = "$\\text{BT}_{50}f^{\\text{FB}}_{0.05}$",
              "HS50-FB15" = "$\\text{BT}_{50}f^{\\text{FB}}_{0.15}$",
              "HS50-FB25" = "$\\text{BT}_{50}f^{\\text{FB}}_{0.25}$",
              "HS50-FB35" = "$\\text{BT}_{50}f^{\\text{FB}}_{0.35}$",
              "HS50-FB45" = "$\\text{BT}_{50}f^{\\text{FB}}_{0.45}$",
              "HS50-A05" = "$\\text{BT}_{50}f^{\\text{CFB}}_{0.005}$",
              "HS50-A5" = "$\\text{BT}_{50}f^{\\text{CFB}}_{0.05}$",
              "HS50-A15" = "$\\text{BT}_{50}f^{\\text{CFB}}_{0.15}$",
              "HS50-A25" = "$\\text{BT}_{50}f^{\\text{CFB}}_{0.25}$",
              "HS50-A35" = "$\\text{BT}_{50}f^{\\text{CFB}}_{0.35}$",
              "HS50-A45" = "$\\text{BT}_{50}f^{\\text{CFB}}_{0.45}$"
              ))

}


## for plotting header, labels, etc.
renameHCRs2 <- function(x){

    plyr::revalue(x,
                  c("MSY" = expression(italic(f)[0.5]),

                    "HS25" = expression(BT[0.25]),
                    "HS50" = expression(BT[0.5]),
                    "HS75" = expression(BT[0.75]),
                    "HS100" = expression(BT[1]),
                    "HS-100" = expression(BT[1]),
                    "HS125" = expression(BT[1.25]),
                    "HS150" = expression(BT[1.5]),
                    "HS175" = expression(BT[1.75]),
                    "HS200" = expression(BT[2]),
                    "HS300" = expression(BT[3]),
                    "MSY-C05" = expression(italic(f)[0.005]^{C}),
                    "MSY-C5" = expression(italic(f)[0.5]^{C}),
                    "MSY-C15" = expression(italic(f)[0.15]^{C}),
                    "MSY-C25" = expression(italic(f)[0.25]^{C}),
                    "MSY-C35" = expression(italic(f)[0.35]^{C}),
                    "MSY-C45" = expression(italic(f)[0.45]^{C}),

                    "MSY-F05" = expression(italic(f)[0.005]^{F}),
                    "MSY-F5" = expression(italic(f)[0.05]^{F}),
                    "MSY-F15" = expression(italic(f)[0.15]^{F}),
                    "MSY-F25" = expression(italic(f)[0.25]^{F}),
                    "MSY-F35" = expression(italic(f)[0.35]^{F}),
                    "MSY-F45" = expression(italic(f)[0.45]^{F}),

                    "MSY-A05" = expression(italic(f)[0.005]^{CF}),
                    "MSY-A5" = expression(italic(f)[0.05]^{CF}),
                    "MSY-A15" = expression(italic(f)[0.15]^{CF}),
                    "MSY-A25" = expression(italic(f)[0.25]^{CF}),
                    "MSY-A35" = expression(italic(f)[0.35]^{CF}),
                    "MSY-A45" = expression(italic(f)[0.45]^{CF}),

                    "HS50-C05" = expression(BT[0.5]*italic(f)[0.005]^{C}),
                    "HS50-C5" = expression(BT[0.5]*italic(f)[0.05]^{C}),
                    "HS50-C15" = expression(BT[0.5]*italic(f)[0.15]^{C}),
                    "HS50-C25" = expression(BT[0.5]*italic(f)[0.25]^{C}),
                    "HS50-C35" = expression(BT[0.5]*italic(f)[0.35]^{C}),
                    "HS50-C45" = expression(BT[0.5]*italic(f)[0.45]^{C}),

                    "HS50-F05" = expression(BT[0.5]*italic(f)[0.005]^{F}),
                    "HS50-F5" = expression(BT[0.5]*italic(f)[0.05]^{F}),
                    "HS50-F15" = expression(BT[0.5]*italic(f)[0.15]^{F}),
                    "HS50-F25" = expression(BT[0.5]*italic(f)[0.25]^{F}),
                    "HS50-F35" = expression(BT[0.5]*italic(f)[0.35]^{F}),
                    "HS50-F45" = expression(BT[0.5]*italic(f)[0.45]^{F}),

                    "HS50-B05" = expression(BT[0.5]*italic(f)[0.005]^{B}),
                    "HS50-B5" = expression(BT[0.5]*italic(f)[0.05]^{B}),
                    "HS50-B15" = expression(BT[0.5]*italic(f)[0.15]^{B}),
                    "HS50-B25" = expression(BT[0.5]*italic(f)[0.25]^{B}),
                    "HS50-B35" = expression(BT[0.5]*italic(f)[0.35]^{B}),
                    "HS50-B45" = expression(BT[0.5]*italic(f)[0.45]^{B}),

                    "HS50-FB05" = expression(BT[0.5]*italic(f)[0.005]^{FB}),
                    "HS50-FB5" = expression(BT[0.5]*italic(f)[0.05]^{FB}),
                    "HS50-FB15" = expression(BT[0.5]*italic(f)[0.15]^{FB}),
                    "HS50-FB25" = expression(BT[0.5]*italic(f)[0.25]^{FB}),
                    "HS50-FB35" = expression(BT[0.5]*italic(f)[0.35]^{FB}),
                    "HS50-FB45" = expression(BT[0.5]*italic(f)[0.45]^{FB}),

                    "HS50-A05" = expression(BT[0.5]*italic(f)[0.005]^{CFB}),
                    "HS50-A5" = expression(BT[0.5]*italic(f)[0.05]^{CFB}),
                    "HS50-A15" = expression(BT[0.5]*italic(f)[0.15]^{CFB}),
                    "HS50-A25" = expression(BT[0.5]*italic(f)[0.25]^{CFB}),
                    "HS50-A35" = expression(BT[0.5]*italic(f)[0.35]^{CFB}),
                    "HS50-A45" = expression(BT[0.5]*italic(f)[0.45]^{CFB})
                    ))

}


## rename HCRs additional (For tables)
renameHCRs3 <- function(x){

    plyr::revalue(x,
                  c("MSY" = "$f_{0.5}$",
                    "HSx" = "$\\text{BT}$",
                    "MSY-C" = "$f^{\\text{C}}$",
                    "MSY-F" = "$f^{\\text{F}}$",
                    "MSY-A" = "$f^{\\text{FB}}$",
                    "HS50-C" = "$\\text{BT}_{0.5}f^{\\text{C}}$",
                    "HS50-F" = "$\\text{BT}_{0.5}f^{\\text{F}}$",
                    "HS50-B" = "$\\text{BT}_{0.5}f^{\\text{B}}$",
                    "HS50-FB" = "$\\text{BT}_{0.5}f^{\\text{FB}}$",
                    "HS50-A" = "$\\text{BT}_{0.5}f^{\\text{CFB}}$"
                    ))
}



## plotting


getPchCol <- function(hcr){
    hcr2 <- unlist(strsplit(as.character(hcr), "-"))[1]
    if(hcr %in% c("noF")){
        pch <- 3
        col <- "grey50"
        cex <- 1.2
    }else if(hcr2 %in% c("refFmsy")){
        pch <- 3
        col <- "grey50"
        cex <- 1.2
    }else if(hcr %in% "MSY"){
        pch <- 18
        col <- "black"
        cex <- 2.5
    }else if(hcr %in% paste0("HS50")){
        pch <- 15
        col <- "black"
        cex <- 2
    }else if(hcr %in% paste0("HS",c(25,seq(75,300,25)))){
        pch <- 17
        col <- "chartreuse4"
        cex <- 2
    }else if(hcr %in% paste0("HS-",c(25,seq(75,300,25)))){
        pch <- 17
        col <- "chartreuse4"
        cex <- 2
    }else if(hcr %in% c("MSY-C05",paste0("MSY-C",seq(5,45,5)))){
        pch <- 16
        col <- "#083b8f" ## "deepskyblue4"
        cex <- 2
    }else if(hcr %in% c("MSY-F05",paste0("MSY-F",seq(5,45,5)))){
        pch <- 14 # 16
        col <- "#3073e3" ## "deepskyblue2"
        cex <- 2
    }else if(hcr %in% c("MSY-A05",paste0("MSY-A",seq(5,45,5)))){
        pch <- 13 # 16
        col <- "#9ebdf1" ## "cadetblue2"
        cex <- 2
    }else if(hcr %in% c("HS50-C05",paste0("HS50-C",seq(5,45,5)))){
        pch <- 12 #15
        col <- "#f2ac41" ## "#f4dc8a" ## "gold1"
        cex <- 2
    }else if(hcr %in% c("HS50-F05",paste0("HS50-F",seq(5,45,5)))){
        pch <- 9 # 15
        col <- "#f1650a" ## "darkorange1"
        cex <- 2
    }else if(hcr %in% c("HS50-B05",paste0("HS50-B",seq(5,45,5)))){
        pch <- 8 #15
        col <- "#c2420b"  ## "firebrick2"
        cex <- 2
    }else if(hcr %in% c("HS50-FB05",paste0("HS50-FB",seq(5,45,5))) || hcr %in% c("HS50-F35-B35")){
        pch <- 5 # 15
        col <- "#ed2215" ## "firebrick4"
        cex <- 2
    }else if(hcr %in% c("HS50-A05",paste0("HS50-A",seq(5,45,5)))){
        pch <- 2 #15
        col <- "#7c0404" ##"brown"
        cex <- 2
    }
    return(c(pch,col,cex))
}



getPchCol2 <- function(hcr){
    hcr2 <- unlist(strsplit(as.character(hcr), "-"))[1]
    if(hcr %in% c("noF")){
        pch <- 3
        col <- "grey50"
        cex <- 1.2
    }else if(hcr2 %in% c("refFmsy")){
        pch <- 3
        col <- "grey50"
        cex <- 1.2
    }else if(hcr %in% "MSY"){
        pch <- 18
        col <- "black"
        cex <- 2.5
    }else if(hcr %in% paste0("HS50")){
        pch <- 15
        col <- "black"
        cex <- 2
    }else if(hcr %in% paste0("HS",c(25,seq(75,300,25)))){
        pch <- 17
        col <- "chartreuse4"
        cex <- 2
    }else if(hcr %in% c("MSY-C05",paste0("MSY-C",seq(5,45,5)))){
        pch <- 1
        col <- "dodgerblue1"
        cex <- 2
    }else if(hcr %in% c("MSY-F05",paste0("MSY-F",seq(5,45,5)))){
        pch <- 14 # 16
        col <- "dodgerblue1"
        cex <- 2
    }else if(hcr %in% c("MSY-A05",paste0("MSY-A",seq(5,45,5)))){
        pch <- 13 # 16
        col <- "dodgerblue1"
        cex <- 2
    }else if(hcr %in% c("HS50-C05",paste0("HS50-C",seq(5,45,5)))){
        pch <- 12 #15
        col <- "darkorange"
        cex <- 2
    }else if(hcr %in% c("HS50-F05",paste0("HS50-F",seq(5,45,5)))){
        pch <- 9 # 15
        col <- "darkorange"
        cex <- 2
    }else if(hcr %in% c("HS50-B05",paste0("HS50-B",seq(5,45,5)))){
        pch <- 8 #15
        col <- "darkorange"
        cex <- 2
    }else if(hcr %in% c("HS50-FB05",paste0("HS50-FB",seq(5,45,5))) || hcr %in% c("HS50-F35-B35")){
        pch <- 5 # 15
        col <- "darkorange"
        cex <- 2
    }else if(hcr %in% c("HS50-A05",paste0("HS50-A",seq(5,45,5)))){
        pch <- 2 #15
        col <- "darkorange"
        cex <- 2
    }
    return(c(pch,col,cex))
}


getPchCol3 <- function(hcr){
    hcr2 <- unlist(strsplit(as.character(hcr), "-"))[1]
    cols <- brewer.pal(8,"Dark2")
    if(hcr %in% c("noF")){
        pch <- 3
        col <- "grey50"
        cex <- 1.2
    }else if(hcr2 %in% c("refFmsy")){
        pch <- 3
        col <- "grey50"
        cex <- 1.2
    }else if(hcr %in% "MSY"){
        pch <- 16
        col <- "black"
        cex <- 2
    }else if(hcr %in% paste0("HS50")){
        pch <- 15
        col <- "black"
        cex <- 2
    }else if(hcr %in% paste0("HS",c(25,seq(75,300,25)))){
        pch <- 17
        col <- "chartreuse4"
        cex <- 2.5
    }else if(hcr %in% c("MSY-C05",paste0("MSY-C",seq(5,45,5)))){
        pch <- 16
        col <- cols[1] ##"dodgerblue1"
        cex <- 2
    }else if(hcr %in% c("MSY-F05",paste0("MSY-F",seq(5,45,5)))){
        pch <- 16
        col <- cols[2] ## "firebrick"
        cex <- 2
    }else if(hcr %in% c("MSY-A05",paste0("MSY-A",seq(5,45,5)))){
        pch <- 16
        col <- cols[3] ##"darkorange"
        cex <- 2
    }else if(hcr %in% c("HS50-C05",paste0("HS50-C",seq(5,45,5)))){
        pch <- 15
        col <- cols[4] ##"chartreuse4"
        cex <- 2
    }else if(hcr %in% c("HS50-F05",paste0("HS50-F",seq(5,45,5)))){
        pch <- 15
        col <- cols[5] ##"purple"
        cex <- 2
    }else if(hcr %in% c("HS50-B05",paste0("HS50-B",seq(5,45,5)))){
        pch <- 15
        col <- cols[6] ##"goldenrod3"
        cex <- 2
    }else if(hcr %in% c("HS50-FB05",paste0("HS50-FB",seq(5,45,5))) || hcr %in% c("HS50-F35-B35")){
        pch <- 15
        col <- cols[7] ##"brown"
        cex <- 2
    }else if(hcr %in% c("HS50-A05",paste0("HS50-A",seq(5,45,5)))){
        pch <- 15
        col <- cols[8] ##"grey50"
        cex <- 2
    }
    return(c(pch,col,cex))
}
