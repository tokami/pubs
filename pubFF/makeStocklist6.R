## ---------------##
## Define species ##
## ---------------##


## updates after probHCR review:
## - updated IMASE
## - species-dependent number of intra-annual time steps
## - species-dependent spawning event
## - parameterisation based on assessment reports (e.g. not Gislason for M)

## install.packages(c("TMB","Rcpp","doBy"))
## remotes::install_github("fishfollower/SAM/stockassessment")
## remotes::install_github("DTUAqua/spict/spict")

## load package
## remotes::install_github("tokami/iamse", ref = "pubFF")
require(iamse)
try(source("loadFuncs.R"))
try(source("../res/loadFuncs.R"))


## Variables
saveFig <- TRUE
saveDat <- TRUE
estDepl <- TRUE
estProd <- FALSE


if(saveFig) pdf("stocklist6.pdf")

## settings
ny <- 35
nrepDepl <- 1e3
deplmethod <- "percentile"
ncores <- 10
fmax <- 100
nf <- 1e3
maxF <- 10

## no noise
set0 <- check.set()
set0$burnin <- 500
set0$refN <- 1e3
set0$refYears <- 300
set0$refYearsMSY <- 100
set0$refMethod <- "median"
set0$spType <- 0 ## TSB


## with noise
set <- set0
set$noiseF <- c(0.15,0,1)
set$noiseM <- c(0,0,0)
set$noiseH <- c(0,0,0)
set$noiseMat <- c(0,0,0)
set$noiseR0 <- c(0,0,0)
set$noiseImp <- c(0,0,0)
set$noiseC <- c(0.3,0,0)
set$noiseI <- c(0.3,0,0)


## Anchovy in Subarea 8 (Bay of Biscay)
## --------------------
## Sources:

## 1) actually from Andres, waiting on his response how to cite them

## 1) ICES. 2019. Working Group on Southern Horse Mackerel, Anchovy and Sardine
## (WGHANSA). ICES Scientific Reports. 1:34. 653 pp.
## http://doi.org/10.17895/ices.pub.4983

## 4) Froese, R. and D. Pauly. Editors. 2019. FishBase. World Wide Web electronic
## publication. www.fishbase.org, ( 12/2019 )

## 5) Thorson, J. T., Jensen, O. P., & Zipkin, E. F. (2014). How variable is
## recruitment for exploited marine fishes? A hierarchical model for testing
## life history theory. Canadian Journal of Fisheries and Aquatic Sciences,
## 71(7), 973-983.

## number of seasons
ns <- 2

## maximum F for historic F vector (per season)
fmax2 <- 1.4

anchovy <- list(
    stock = "ane.27.8",
    species = "Anchovy", ## European Anchovy
    latin = "Engraulis encrasicolus",
    family = "Engraulidae",   ## Clupeiformes
    icesWG = "WGHANSA",
    ny = ny,
    ns = ns,
    amax = 4,           ## (1)  Fishbase = 5, Andres = 6, HANSA no data for older than 4
    R0 = 1e5,
    SR = "bevholt",
    h = 0.75, ## before: 0.74       ## older: 0.75, ## (1), alternatively: 0.75, 0.9 ## assumed by Andres et al.
    bp = 0,
    Linf = 18.69,
    K = 0.89,
    t0 = -0.02,
    lwa = 0.004799048,
    lwb = 3.134380952,
    fecun = 1,
    q = 0.05,
    depl.quant = "Blim",
    depl = 1,
    depl.prob = 0.5,
    sigmaR = 0.766,  ## (5)
    rhoR = 0.435,    ## (5)
    biascorR = 1,     ## (5)
    age0 = 0
)
anchovy <- check.dat(anchovy)
set$noiseR <- c(anchovy$sigmaR, anchovy$rhoR, anchovy$biascorR)
set$maxF <- maxF

## ICES HANSA 2020 p. 35
MAA <- c(1.2, 0.8, rep(1.2, anchovy$amax-1)) ## age0 p.97
anchovy$M <- matrix(rep(max(MAA)/ns, ny * ns), nrow=ny, ncol=ns)
anchovy$Msel[[1]] <- matrix(rep(MAA/max(MAA),ns), nrow=anchovy$amax+1, ncol=ns)

## Andres
anchovy$mat <- matrix(c(0,1,
                        rep(rep(1,anchovy$amax),each=ns)), nrow=anchovy$amax+1, ncol=ns, byrow=TRUE)

## Andres
anchovy$sel[[1]] <- matrix(c(0.01,1,
                             rep(rep(1,anchovy$amax),each=ns)), nrow=anchovy$amax+1, ncol=ns, byrow=TRUE)




## ICES HANSA 2020 p.17 (values for 2019) ## no value provided for age 0 quarter 1 -> assumed based on VBGF + ALK
## anchovy$weight <- matrix(c(0.4, 20.2, 27.4, 32.2, 27.7,
##                            11.9, 21, 26, 33.6, 33.6),
##                          nrow=anchovy$amax+1, ncol=ns)
## anchovy$weightF <- matrix(c(0.4, 20.2, 27.4, 32.2, 27.7,
##                            11.9, 21, 26, 33.6, 33.6),
##                          nrow=anchovy$amax+1, ncol=ns)
## USE VBGF parameters for weight at age, VBGF parameters from Andres + Sonia

FM <- c(seq(0.01,fmax2, length.out=ny-5),rep(fmax2,5))
anchovy$FM <- matrix(rep(FM/ns, ns), nrow=ny, ncol=ns)

## Previous biological studies based on commercial samples of GoC anchovy (9a S (ES))
## indicate that the species’ spawning season extends from late winter to early autumn
## with a peak spawning time for the whole population occurring from June to August
## (Millán, 1999)
if(ns == 2) anchovy$spawning <- c(0,1)

## survey selectivities and timings
anchovy$selI <- list()
anchovy$selI[[1]] <- anchovy$sel[[1]]
anchovy$selI[[2]] <- anchovy$sel[[1]]

anchovy$surveyTimes <- c(1/12,7/12)

anchovy$nyC <- 35
anchovy$nyI <- c(35, 27)

anchovy$catchSeasons <- 1

anchovy$surveyBeforeAssessment <- c(TRUE, FALSE)




## Reference points
par(mfrow=c(1,1))
plot.new()
title("Anchovy")


anchovy <- est.ref.levels.stochastic(anchovy, set, fmax = fmax,
                                     ncores=ncores, plot = TRUE)

fmsy <- anchovy$ref$Fmsy[1]

print(paste0("Fmsy = ", round(anchovy$ref$Fmsy[1],3)))
print(paste0("Bmsy = ", round(anchovy$ref$Bmsy[1],2)))
print(paste0("MSY = ", round(anchovy$ref$MSY[1],2)))
print(paste0("B0 = ", round(anchovy$ref$B0[1],2)))
print(paste0("ESBmsy = ", round(anchovy$ref$ESBmsy[1],2)))
print(paste0("SSBmsy = ", round(anchovy$ref$SSBmsy[1],2)))

print(paste0("Bmsy/B0 = ", round(anchovy$ref$Bmsy[1]/anchovy$ref$B0[1]*100,2)))

fmax  <- 10

## Blim
tmp <- est.blim(anchovy, set,
                msy = anchovy$ref$MSY[1],
                ref = 0.5,
                fmin = anchovy$ref$Fmsy[1],
                fmax = fmax)
anchovy$ref$Blim <- tmp$ref$Bmsy


## Stock-recruitment relationship + Blim
ssbs <- seq(0, max(anchovy$ref$B0), length.out = 5e4)
recs <- recfunc(anchovy$h, get.ssb0(as.numeric(t(anchovy$M[1,]*anchovy$Msel[[1]])), as.numeric(t(anchovy$mat)),
                                    as.numeric(t(anchovy$weight)), 1,
                                    (anchovy$amax+1) * anchovy$ns, anchovy$ns,
                                    anchovy$spawning, anchovy$R0,anchovy$indage0,
                                    season = 1)/anchovy$R0,
                ssbs, anchovy$R0, anchovy$SR, NULL, NULL, NULL)

r0 <- max(recs)  ## or: r0 <- dat$R0
## ssblim <- 0.2 * anchovy$ref$B0 (Dichmont 2017)
## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
## blim <- 0.3 * anchovy$ref$Bmsy[1]
## anchovy$ref$Blim <- blim
print(paste0("Blim = ", round(anchovy$ref$Blim[1],2)))
print(paste0("Blim/Bmsy = ",  round(anchovy$ref$Blim[1] / anchovy$ref$Bmsy[1],2)))

par(mfrow=c(1,2))
plot(ssbs, recs, ty='l',
     ylim = range(0,1.2*anchovy$R0), lwd=1.5)
abline(h = anchovy$R0, lty=2, lwd = 1.5)
## abline(v = anchovy$ref$SSBlim, lwd=1.5, lty=2)


## plot recruitment deviations
n <- 5
tmp <- list()
for(i in 1:n) tmp[[i]] <-  anchovy$R0 * gen.noise(35, set$noiseR[1], set$noiseR[2], set$noiseR[3])
plot(1:35, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
for(i in 1:n){
    lines(1:35, tmp[[i]], ty='b', col=i)
}


## Set FM based on depletion level
if(estDepl) anchovy <- est.depletion(anchovy, set, fmax = fmax, nrep = nrepDepl,
                                     method = deplmethod)


## Production curve
if(estProd){
    par(mfrow=c(1,1))
    prod <- est.productivity.stochastic(anchovy, set, 100, nf, ncores = ncores, plot = TRUE)
    abline(h=0)
    abline(v = anchovy$ref$Blim[1],lty=2)
    abline(v = round(anchovy$ref$Bmsy[1],2),lty=3)
}






## Haddock in Celtic Seas (divisions 7.b,c,e-k)
## --------------------------------------
## Sources:

## 1) ICES 2020. WGCSE

## sd <- 0.16
## tmp <- rnorm(1e4, log(0.74), sd) - sd^2/2
## round(quantile(exp(tmp), c(0.2,0.8)),2)

## number of seasons
ns <- 1

## maximum F for historic F vector (per season)
fmax2 <- 1

haddock <- list(
    stock = "had",
    species = "Haddock",
    latin = "Melanogrammus aeglefinus",
    family = "Gadidae",
    icesWG = "WGCSE",
    ny = ny,
    ns = ns,
    amax = 8,      ## potentially older (15), but plus group and parameters available up until age 8
    R0 = 1e6,
    SR = "bevholt",
    h = 0.75, ## 0.63, ## thorson 2019 Fish&Fisheries ## older: 0.74,       ## (3)
    bp = 0,
    fecun = 1,
    q = 0.05,
    depl.quant = "Blim",
    depl = 1,
    depl.prob = 0.5,
    sigmaR = 0.748, ## (5)
    rhoR = 0.404,   ## (5)
    biascorR = 1,    ## (5)
    sigmaH = 0.16, ## Myers et. 1999
    rhoH = 0,   ## Myers et. 1999
    biascorH = 1    ## Myers et. 1999
)
haddock <- check.dat(haddock)
set$noiseR <- c(haddock$sigmaR, haddock$rhoR, haddock$biascorR)
## set$noiseH <- c(haddock$sigmaH, haddock$rhoH, haddock$biascorH)
set$maxF <- maxF

## ICES WGCSE 2020 p.321
MAA <- c(1.087,0.721,0.575,0.483,0.44,0.406,0.398,0.385,0.36)
haddock$M <- matrix(max(MAA), nrow=ny, ncol=ns)
haddock$Msel[[1]] <- matrix(rep(MAA/max(MAA)/ns,ns), nrow=haddock$amax+1, ncol=ns)

## ICES WGCSE 2020 p.321
haddock$mat <- matrix(rep(c(0,0.039,0.911,0.970,0.98,
                            rep(1,haddock$amax-4)),ns), nrow=haddock$amax+1, ncol=ns)

## ICES WGCSE 2020 (stockassessment.org: mean stock weights)
WAA <- c(0.065,0.175,0.368,0.693,1.118,1.220,1.826,2.008,2.209)
haddock$weight <- matrix(rep(WAA,ns), nrow=haddock$amax+1, ncol=ns)
##WAA <- c(0,0.184,0.428,0.700,1.024,1.204,1.721,1.623,1.840)
haddock$weightF <- matrix(rep(WAA,ns), nrow=haddock$amax+1, ncol=ns)

## ICES WGCSE 2020 p. 337 (values for 2019)
## SAM estimated F-at-age:
FAA <- c(0.01,0.116,0.33,0.371,0.403,0.444,0.502,0.693,0.693) ## age0 = 0 (see stock annex)
FM <- matrix(c(seq(0.01,fmax2, length.out=ny-5),rep(fmax2,5)), nrow=ny, ncol=ns)
haddock$FM <- FM/max(FM) * max(FAA) / ns
haddock$sel[[1]] <- matrix(rep(FAA/max(FAA),ns), nrow=haddock$amax+1, ncol=ns)

## Date of peak spawning: Day 85 (26 March), independent of location or age group Calculated from Field survey data , 29 March ± 13 April 1992: see Appendix 1  (Heath & Gallego, 1998), normal dist with sigma = 16 days
if(ns==4) haddock$spawning <- c(0,1,0,0)

## survey selectivities and timings
haddock$selI <- list()
haddock$selI[[1]] <- haddock$sel[[1]]
haddock$selI[[2]] <- haddock$sel[[1]]

haddock$surveyTimes <- c(1/12,7/12)

haddock$nyC <- 35
haddock$nyI <- c(35, 27)

haddock$catchSeasons <- 1

haddock$surveyBeforeAssessment <- c(FALSE, FALSE)


## Reference points
par(mfrow=c(1,1))
plot.new()
title("Haddock")

haddock <- est.ref.levels.stochastic(haddock, set, fmax = fmax, ncores=ncores, plot = TRUE)

fmsy <- haddock$ref$Fmsy[1]

## ICES MSY reference levels p.337: Fmsy = 0.354, Blim = 9227
print(paste0("Fmsy = ", round(haddock$ref$Fmsy[1],3)))
print(paste0("Bmsy = ", round(haddock$ref$Bmsy[1],2)))
print(paste0("MSY = ", round(haddock$ref$MSY[1],2)))
print(paste0("B0 = ", round(haddock$ref$B0[1],2)))
print(paste0("ESBmsy = ", round(haddock$ref$ESBmsy[1],2)))
print(paste0("SSBmsy = ", round(haddock$ref$SSBmsy[1],2)))

print(paste0("Bmsy/B0 = ", round(haddock$ref$Bmsy[1]/haddock$ref$B0[1]*100,2)))

## Blim
tmp <- est.blim(haddock, set,
                msy = haddock$ref$MSY[1],
                ref = 0.5,
                fmin = haddock$ref$Fmsy[1],
                fmax = fmax)
haddock$ref$Blim <- tmp$ref$Bmsy

## Stock-recruitment relationship + Blim
ssbs <- seq(0, max(haddock$ref$B0), length.out = 5e4)
ssb0 <- get.ssb0(as.numeric(t(haddock$M[1,]*haddock$Msel[[1]])),
                 as.numeric(t(haddock$mat)),
                 as.numeric(t(haddock$weight)), 1,
                 (haddock$amax+1) * haddock$ns, haddock$ns,
                 haddock$spawning, haddock$R0,haddock$indage0,
                 season = 1)
recs <- recfunc(haddock$h, ssb0/haddock$R0, ssbs, haddock$R0, haddock$SR, NULL, NULL, NULL)

r0 <- max(recs)  ## or: r0 <- dat$R0
## ssblim <- 0.2 * haddock$ref$B0 (Dichmont 2017)
## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
## blim <- 0.3 * haddock$ref$Bmsy[1]
## haddock$ref$Blim <- blim   ## in terms of SSB now (compare to SSBfinal)
print(paste0("Blim/Bmsy = ", round(haddock$ref$Blim[1] / haddock$ref$Bmsy[1],2)))
print(paste0("Blim = ", round(haddock$ref$Blim[1],2)))

par(mfrow=c(1,2))
plot(ssbs, recs, ty='l',
     ylim = range(0,1.2*haddock$R0), lwd=1.5)
abline(h = haddock$R0, lty=2, lwd = 1.5)
## abline(v = haddock$ref$SSBlim, lwd=1.5, lty=2)


## plot recruitment deviations
n <- 5
tmp <- list()
for(i in 1:n) tmp[[i]] <-  haddock$R0 * gen.noise(35, set$noiseR[1], set$noiseR[2], set$noiseR[3])
plot(1:35, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
for(i in 1:n){
    lines(1:35, tmp[[i]], ty='b', col=i)
}

par(mfrow=c(1,1))
## pdf("recruitment_haddock.pdf",width=8, height=6)
years <- 1993:2019
icesRecs <- c(138345, 399725, 467464, 164566, 59414,90819, 355678, 351066, 478013, 964657, 233720, 351353,
              254303, 198880, 715742, 411306, 2283436, 203340, 90043, 68005, 647545, 220385, 509828,
              111283, 176611, 946209, 493214)
n <- 3
tmp <- list()
for(i in 1:n) tmp[[i]] <-  haddock$R0 * gen.noise(27, set$noiseR[1], set$noiseR[2], set$noiseR[3])
## plot(years, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
plot(years, icesRecs, ty='n', ylim = range(c(icesRecs,unlist(tmp))),
     ylab="Recruits", xlab = "")
for(i in 1:n){
    lines(years, tmp[[i]], ty='b', col="grey70", lty=2, lwd=1.5)
}
lines(years, icesRecs, ty='b', col=1, lwd=2)
legend("topright", legend=c("ICES","simulated"),
       lty= c(1,2), col = c(1, "grey70"),
       pch=1, bg="white",
       lwd=c(2,1.5), cex=1.2)
box(lwd=1.5)
## dev.off()


## Set FM based on depletion level
if(estDepl) haddock <- est.depletion(haddock, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)


## Production curve
if(estProd){
    par(mfrow=c(1,1))
    prod <- est.productivity.stochastic(haddock, set, 100, nf, ncores = ncores, plot = TRUE)
    abline(h=0)
    abline(v = haddock$ref$Blim[1],lty=2)
    abline(v = round(haddock$ref$Bmsy[1],2),lty=3)
}



## Greenland halibut
## ----------------------------------------

## Sources:

## 1) ICES. 2018. Report of the North-Western Working Group (NWWG), 26 April–3
## May, 2018, ICES HQ, Copenhagen, Denmark. ICES CM 2018/ACOM:09. 733 pp.

## 2) Jardim, E., Azevedo, M., & Brites, N. M. (2015). Harvest control rules for
## data limited stocks using length-based reference points and survey biomass
## indices. Fisheries research, 171, 12-19.

## 3) Myers, R. A., Bowen, K. G., & Barrowman, N. J. (1999). Maximum
## reproductive rate of fish at low population sizes. Canadian Journal of
## Fisheries and Aquatic Sciences, 56(12), 2404-2419.

## 4) Froese, R. and D. Pauly. Editors. 2019.FishBase. World Wide Web electronic
## publication. www.fishbase.org, ( 12/2019 )

## 5) Rickman, S.J., N.K. Dulvy, S. Jennings and J.D. Reynolds, 2000.
## Recruitment variation related to fecundity in marine fishes. Can. J. Fish.
## Aquat. Sci. 57:116-124.

## 6) ICES. 2013. Report of the Benchmark Workshop on Greenland Halibut Stocks
## (WKBUT), 26–29 November 2013, Copenhagen, Denmark. ICES CM 2013/ ACOM:44. 367
## pp.

## 7) Thorson, J. T., Jensen, O. P., & Zipkin, E. F. (2014). How variable is
## recruitment for exploited marine fishes? A hierarchical model for testing
## life history theory. Canadian Journal of Fisheries and Aquatic Sciences,
## 71(7), 973-983.

## 8) Thorson, J. T., Predicting recruitment density dependence and intrinsic
## growth rate for all fishes worldwide using a data‐integrated life‐history
## model (2019)

## sd <- 0.3
## tmp <- rnorm(1e4, log(0.79), sd) - sd^2/2
## round(quantile(exp(tmp), c(0.2,0.8)),2)

ns <- 1

fmax2 <- 0.5

halibut <- list(
    stock = "ghl-dis",
    species = "Greenland halibut",
    latin = "Reinhardtius hippoglossoides",
    family = "Pleuronectidae",
    icesWG = "NWWG",
    ny = ny,
    ns = ns,
    amax = 27,     ## (4) (FishBase: 30 years but plus group)
    lwa = 0.00333, ## (2)
    lwb = 3.249,   ## (2)
    a0 = -0.1,     ## (2)
    K = 0.073,     ## (2)
    Linf = 120,    ## (2)
    Lm50 = 71.2,   ## (5)
    Lm95 = 81.2,   ## (5)
    Ls50 = 51,      ## (6)
    Ls95 = 58.23,   ## (6)
    R0 = 1e3,
    SR = "bevholt",
    h = 0.75, ## 0.76 ## (8)   ## older: ## 0.79,      ## (3)
    bp = 0,
    fecun = 1,
    q = 0.05,
    binwidth = 1,
    CVlen = 0.1,
    depl.quant = "Blim",
    depl = 1,
    depl.prob = 0.5,
    sigmaR = 0.636,  ## (7)
    rhoR = 0.437,    ## (7)
    biascorR = 1,    ## (7)
    sigmaH = 0.3,    ## Myers et al. 1999
    rhoRH = 0,       ## Myers et al. 1999
    biascorH = 1     ## Myers et al. 1999
)
halibut <- check.dat(halibut)
set$noiseR <- c(halibut$sigmaR, halibut$rhoR, halibut$biascorR)
## set$noiseH <- c(haddock$sigmaH, haddock$rhoH, haddock$biascorH)
set$maxF <- maxF

## M as assumed in Gadget models WKBUT 2013
halibut$M <- matrix(0.1, nrow = ny, ncol = ns)
halibut$Msel[[1]] <- matrix(1, nrow = halibut$amax + 1, ncol = ns)

## FM
FM <- c(seq(0.01,fmax2, length.out=ny-5),rep(fmax2,5))
halibut$FM <- matrix(rep(FM, ns), nrow=ny, ncol=ns)

halibut$spawning <- 1

## survey selectivities and timings
halibut$selI <- list()
halibut$selI[[1]] <- halibut$sel[[1]]
halibut$selI[[2]] <- halibut$sel[[1]]

halibut$surveyTimes <- c(1/12,7/12)

halibut$nyC <- 35
halibut$nyI <- c(35, 27)

halibut$catchSeasons <- 1

halibut$surveyBeforeAssessment <- c(FALSE, FALSE)

## Figure out selectivity
if(FALSE){
    ## selectivity (1)
    l50 <- 51
    ptarget <- 0.08 ## read from graph in (6)
    fn <- function(l95){
        pl45 <- 1 /(1 + exp(-log(19)*(45 - l50)/(l95 - l50)))
        (pl45 - ptarget)^2
    }
    optimise(fn, c(0,100))

}


## Reference points
par(mfrow=c(1,1))
plot.new()
title("Halibut")


halibut <- est.ref.levels.stochastic(halibut, set, fmax = fmax,
                                     ncores=ncores, plot = TRUE)

fmsy <- halibut$ref$Fmsy[1]

print(paste0("Fmsy = ", round(halibut$ref$Fmsy[1],3)))
print(paste0("Bmsy = ", round(halibut$ref$Bmsy[1],2)))
print(paste0("MSY = ", round(halibut$ref$MSY[1],2)))
print(paste0("B0 = ", round(halibut$ref$B0[1],2)))
print(paste0("ESBmsy = ", round(halibut$ref$ESBmsy[1],2)))
print(paste0("SSBmsy = ", round(halibut$ref$SSBmsy[1],2)))

print(paste0("Bmsy/B0 = ", round(halibut$ref$Bmsy[1]/halibut$ref$B0[1]*100,2)))

## Blim
tmp <- est.blim(halibut, set,
                msy = halibut$ref$MSY[1],
                ref = 0.5,
                fmin = halibut$ref$Fmsy[1],
                fmax = fmax)
halibut$ref$Blim <- tmp$ref$Bmsy


## Stock-recruitment relationship + Blim
ssbs <- seq(0, max(halibut$ref$B0), length.out = 5e4)
recs <- recfunc(halibut$h, get.ssb0(as.numeric(t(halibut$M[1,]*halibut$Msel[[1]])),
                                    as.numeric(t(halibut$mat)),
                                    as.numeric(t(halibut$weight)), 1,
                                    (halibut$amax+1) * halibut$ns, halibut$ns,
                                    halibut$spawning, halibut$R0,halibut$indage0,
                                    season = 1)/halibut$R0,
                ssbs, halibut$R0, halibut$SR, NULL, NULL, NULL)

r0 <- max(recs)  ## or: r0 <- dat$R0
## ssblim <- 0.2 * halibut$ref$B0 (Dichmont 2017)
## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
## blim <- 0.3 * halibut$ref$Bmsy[1]
## halibut$ref$Blim <- blim   ## in terms of SSB now (compare to SSBfinal)
print(paste0("Blim/Bmsy = ", round(halibut$ref$Blim[1] / halibut$ref$Bmsy[1],2)))
print(paste0("Blim = ", round(halibut$ref$Blim[1],2)))

par(mfrow=c(1,2))
plot(ssbs, recs, ty='l',
     ylim = range(0,1.2*halibut$R0), lwd=1.5)
abline(h = halibut$R0, lty=2, lwd = 1.5)
## abline(v = halibut$ref$SSBlim, lwd=1.5, lty=2)

## plot recruitment deviations
n <- 5
tmp <- list()
for(i in 1:n) tmp[[i]] <-  halibut$R0 * gen.noise(35, halibut$sigmaR, halibut$rhoR, halibut$biascorR)
plot(1:35, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
for(i in 1:n){
    lines(1:35, tmp[[i]], ty='b', col=i)
}


## Define Blim as 0.15 B0 according to ICES fisheries management reference
## points for category 1 and 2 stocks (2017)
## Blim = 0.2 Blim (Dichmont 2017)
## halibut$ref$Blim <- 0.2 * halibut$ref$B0
## print(paste0("Blim = ", round(halibut$ref$Blim,2)))


## Set FM based on depletion level
if(estDepl) halibut <- est.depletion(halibut, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)


## Production curve
if(estProd){
    par(mfrow=c(1,1))
    prod <- est.productivity.stochastic(halibut, set, 100, nf, ncores = ncores, plot = TRUE)
    abline(h=0)
    abline(v = halibut$ref$Blim[1],lty=2)
    abline(v = round(halibut$ref$Bmsy[1],2),lty=3)
}






## Underfished stocks
## --------------------------------------
print("Underfished stocks")

## Anchovy
## --------------
ns <- 2
anchovy2 <- anchovy
anchovy2$depl <- 2
anchovy2$depl.quant <- "Bmsy"
anchovy2$depl.prob <- 0.5
fmax2 <- 1.4
FM <- c(seq(0.01,fmax2, length.out=ny-5),rep(fmax2,5))
## FM <- rev(seq(0.01,fmax2, length.out=ny))
anchovy2$FM <- matrix(rep(FM/ns, ns), nrow=ny, ncol=anchovy2$ns)
set$maxF <- maxF
set$noiseR <- c(anchovy2$sigmaR, anchovy2$rhoR, anchovy2$biascorR)
if(estDepl) anchovy2 <- est.depletion(anchovy2, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)

## Haddock
## --------------
ns <- 1
haddock2 <- haddock
haddock2$depl <- 2
haddock2$depl.quant <- "Bmsy"
haddock2$depl.prob <- 0.5
fmax2 <- 1
FM <- c(seq(0.01,fmax2, length.out=ny-5),rep(fmax2,5))
## FM <- rev(seq(0.01,fmax2, length.out=ny))
haddock2$FM <- matrix(rep(FM, ns), nrow=ny, ncol=haddock2$ns)
set$maxF <- maxF
set$noiseR <- c(haddock2$sigmaR, haddock2$rhoR, haddock2$biascorR)
if(estDepl) haddock2 <- est.depletion(haddock2, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)

## Halibut
## --------------
ns <- 1
halibut2 <- halibut
halibut2$depl <- 2
halibut2$depl.quant <- "Bmsy"
halibut2$depl.prob <- 0.5
fmax2 <- 0.5
FM <- c(seq(0.01,fmax2, length.out=ny-5),rep(fmax2,5))
## FM <- rev(seq(0.01,fmax2, length.out=ny))
halibut2$FM <- matrix(rep(FM, ns), nrow=ny, ncol=halibut2$ns)
set$maxF <- maxF
set$noiseR <- c(halibut2$sigmaR, halibut2$rhoR, halibut2$biascorR)
if(estDepl) halibut2 <- est.depletion(halibut2, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)




## Alternative steepness (h)
## --------------------------------------
print("Alternative steepness (h)")

## Anchovy
## --------------
ns <- 2
anchovy3 <- anchovy
anchovy3$h <- 0.9
set$maxF <- maxF
set$noiseR <- c(anchovy3$sigmaR, anchovy3$rhoR, anchovy3$biascorR)

## Reference points
par(mfrow=c(1,1))
plot.new()
title("Anchovy (h=0.9)")

anchovy3 <- est.ref.levels.stochastic(anchovy3, set, fmax = fmax,
                                     ncores=ncores, plot = TRUE)

fmsy <- anchovy3$ref$Fmsy[1]

print(paste0("Fmsy = ", round(anchovy3$ref$Fmsy[1],3)))
print(paste0("Bmsy = ", round(anchovy3$ref$Bmsy[1],2)))
print(paste0("MSY = ", round(anchovy3$ref$MSY[1],2)))
print(paste0("B0 = ", round(anchovy3$ref$B0[1],2)))
print(paste0("ESBmsy = ", round(anchovy3$ref$ESBmsy[1],2)))
print(paste0("SSBmsy = ", round(anchovy3$ref$SSBmsy[1],2)))

print(paste0("Bmsy/B0 = ", round(anchovy3$ref$Bmsy[1]/anchovy3$ref$B0[1]*100,2)))


## Blim
tmp <- est.blim(anchovy3, set,
                msy = anchovy3$ref$MSY[1],
                ref = 0.5,
                fmin = anchovy3$ref$Fmsy[1],
                fmax = fmax)
anchovy3$ref$Blim <- tmp$ref$Bmsy


## Stock-recruitment relationship + Blim
ssbs <- seq(0, max(anchovy3$ref$B0), length.out = 5e4)
recs <- recfunc(anchovy3$h, get.ssb0(as.numeric(t(anchovy3$M[1,]*anchovy3$Msel[[1]])), as.numeric(t(anchovy3$mat)),
                                    as.numeric(t(anchovy3$weight)), 1,
                                    (anchovy3$amax+1) * anchovy3$ns, anchovy3$ns,
                                    anchovy3$spawning, anchovy3$R0,anchovy3$indage0,
                                    season = 1)/anchovy3$R0,
                ssbs, anchovy3$R0, anchovy3$SR, NULL, NULL, NULL)

r0 <- max(recs)  ## or: r0 <- dat$R0
## ssblim <- 0.2 * anchovy3$ref$B0 (Dichmont 2017)
## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
## blim <- 0.3 * anchovy3$ref$Bmsy[1]
## anchovy3$ref$Blim <- blim
print(paste0("Blim = ", round(anchovy3$ref$Blim[1],2)))
print(paste0("Blim/Bmsy = ",  round(anchovy3$ref$Blim[1] / anchovy3$ref$Bmsy[1],2)))

par(mfrow=c(1,2))
plot(ssbs, recs, ty='l',
     ylim = range(0,1.2*anchovy3$R0), lwd=1.5)
abline(h = anchovy3$R0, lty=2, lwd = 1.5)
## abline(v = anchovy3$ref$SSBlim, lwd=1.5, lty=2)


## plot recruitment deviations
n <- 5
tmp <- list()
for(i in 1:n) tmp[[i]] <-  anchovy3$R0 * gen.noise(35, set$noiseR[1], set$noiseR[2], set$noiseR[3])
plot(1:35, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
for(i in 1:n){
    lines(1:35, tmp[[i]], ty='b', col=i)
}


## Set FM based on depletion level
if(estDepl) anchovy3 <- est.depletion(anchovy3, set, fmax = fmax, nrep = nrepDepl,
                                     method = deplmethod)

## Production curve
if(estProd){
    par(mfrow=c(1,1))
    prod <- est.productivity.stochastic(anchovy3, set, 100, nf, ncores = ncores, plot = TRUE)
    abline(h=0)
    abline(v = anchovy3$ref$Blim[1],lty=2)
    abline(v = round(anchovy3$ref$Bmsy[1],2),lty=3)
}


## Haddock
## ---------------------
ns <- 1
haddock3 <- haddock
haddock3$h <- 0.9
set$maxF <- maxF
set$noiseR <- c(haddock3$sigmaR, haddock3$rhoR, haddock3$biascorR)

## Reference points
par(mfrow=c(1,1))
plot.new()
title("Haddock (h=0.9)")

haddock3 <- est.ref.levels.stochastic(haddock3, set, fmax = fmax, ncores=ncores, plot = TRUE)

fmsy <- haddock3$ref$Fmsy[1]

## ICES MSY reference levels p.337: Fmsy = 0.354, Blim = 9227
print(paste0("Fmsy = ", round(haddock3$ref$Fmsy[1],3)))
print(paste0("Bmsy = ", round(haddock3$ref$Bmsy[1],2)))
print(paste0("MSY = ", round(haddock3$ref$MSY[1],2)))
print(paste0("B0 = ", round(haddock3$ref$B0[1],2)))
print(paste0("ESBmsy = ", round(haddock3$ref$ESBmsy[1],2)))
print(paste0("SSBmsy = ", round(haddock3$ref$SSBmsy[1],2)))

print(paste0("Bmsy/B0 = ", round(haddock3$ref$Bmsy[1]/haddock3$ref$B0[1]*100,2)))

## Blim
tmp <- est.blim(haddock3, set,
                msy = haddock3$ref$MSY[1],
                ref = 0.5,
                fmin = haddock3$ref$Fmsy[1],
                fmax = fmax)
haddock3$ref$Blim <- tmp$ref$Bmsy

## Stock-recruitment relationship + Blim
ssbs <- seq(0, max(haddock3$ref$B0), length.out = 5e4)
ssb0 <- get.ssb0(as.numeric(t(haddock3$M[1,]*haddock3$Msel[[1]])),
                 as.numeric(t(haddock3$mat)),
                 as.numeric(t(haddock3$weight)), 1,
                 (haddock3$amax+1) * haddock3$ns, haddock3$ns,
                 haddock3$spawning, haddock3$R0,haddock3$indage0,
                 season = 1)
recs <- recfunc(haddock3$h, ssb0/haddock3$R0, ssbs, haddock3$R0, haddock3$SR, NULL, NULL, NULL)

r0 <- max(recs)  ## or: r0 <- dat$R0
## ssblim <- 0.2 * haddock3$ref$B0 (Dichmont 2017)
## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
## blim <- 0.3 * haddock3$ref$Bmsy[1]
## haddock3$ref$Blim <- blim   ## in terms of SSB now (compare to SSBfinal)
print(paste0("Blim/Bmsy = ", round(haddock3$ref$Blim[1] / haddock3$ref$Bmsy[1],2)))
print(paste0("Blim = ", round(haddock3$ref$Blim[1],2)))

par(mfrow=c(1,2))
plot(ssbs, recs, ty='l',
     ylim = range(0,1.2*haddock3$R0), lwd=1.5)
abline(h = haddock3$R0, lty=2, lwd = 1.5)
## abline(v = haddock3$ref$SSBlim, lwd=1.5, lty=2)


## plot recruitment deviations
n <- 5
tmp <- list()
for(i in 1:n) tmp[[i]] <-  haddock3$R0 * gen.noise(35, set$noiseR[1], set$noiseR[2], set$noiseR[3])
plot(1:35, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
for(i in 1:n){
    lines(1:35, tmp[[i]], ty='b', col=i)
}

par(mfrow=c(1,1))
## pdf("recruitment_haddock3.pdf",width=8, height=6)
years <- 1993:2019
icesRecs <- c(138345, 399725, 467464, 164566, 59414,90819, 355678, 351066, 478013, 964657, 233720, 351353,
              254303, 198880, 715742, 411306, 2283436, 203340, 90043, 68005, 647545, 220385, 509828,
              111283, 176611, 946209, 493214)
n <- 3
tmp <- list()
for(i in 1:n) tmp[[i]] <-  haddock3$R0 * gen.noise(27, set$noiseR[1], set$noiseR[2], set$noiseR[3])
## plot(years, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
plot(years, icesRecs, ty='n', ylim = range(c(icesRecs,unlist(tmp))),
     ylab="Recruits", xlab = "")
for(i in 1:n){
    lines(years, tmp[[i]], ty='b', col="grey70", lty=2, lwd=1.5)
}
lines(years, icesRecs, ty='b', col=1, lwd=2)
legend("topright", legend=c("ICES","simulated"),
       lty= c(1,2), col = c(1, "grey70"),
       pch=1, bg="white",
       lwd=c(2,1.5), cex=1.2)
box(lwd=1.5)
## dev.off()


## Set FM based on depletion level
if(estDepl) haddock3 <- est.depletion(haddock3, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)


## Production curve
if(estProd){
    par(mfrow=c(1,1))
    prod <- est.productivity.stochastic(haddock3, set, 100, nf, ncores = ncores, plot = TRUE)
    abline(h=0)
    abline(v = haddock3$ref$Blim[1],lty=2)
    abline(v = round(haddock3$ref$Bmsy[1],2),lty=3)
}


## Halibut
## ---------------------
ns <- 1
halibut3 <- halibut
halibut3$h <- 0.9
set$maxF <- maxF
set$noiseR <- c(halibut3$sigmaR, halibut3$rhoR, halibut3$biascorR)

## Reference points
par(mfrow=c(1,1))
plot.new()
title("Halibut (h = 0.9)")

halibut3 <- est.ref.levels.stochastic(halibut3, set, fmax = fmax,
                                     ncores=ncores, plot = TRUE)

fmsy <- halibut3$ref$Fmsy[1]

print(paste0("Fmsy = ", round(halibut3$ref$Fmsy[1],3)))
print(paste0("Bmsy = ", round(halibut3$ref$Bmsy[1],2)))
print(paste0("MSY = ", round(halibut3$ref$MSY[1],2)))
print(paste0("B0 = ", round(halibut3$ref$B0[1],2)))
print(paste0("ESBmsy = ", round(halibut3$ref$ESBmsy[1],2)))
print(paste0("SSBmsy = ", round(halibut3$ref$SSBmsy[1],2)))

print(paste0("Bmsy/B0 = ", round(halibut3$ref$Bmsy[1]/halibut3$ref$B0[1]*100,2)))

## Blim
tmp <- est.blim(halibut3, set,
                msy = halibut3$ref$MSY[1],
                ref = 0.5,
                fmin = halibut3$ref$Fmsy[1],
                fmax = fmax)
halibut3$ref$Blim <- tmp$ref$Bmsy


## Stock-recruitment relationship + Blim
ssbs <- seq(0, max(halibut3$ref$B0), length.out = 5e4)
recs <- recfunc(halibut3$h, get.ssb0(as.numeric(t(halibut3$M[1,]*halibut3$Msel[[1]])),
                                    as.numeric(t(halibut3$mat)),
                                    as.numeric(t(halibut3$weight)), 1,
                                    (halibut3$amax+1) * halibut3$ns, halibut3$ns,
                                    halibut3$spawning, halibut3$R0,halibut3$indage0,
                                    season = 1)/halibut3$R0,
                ssbs, halibut3$R0, halibut3$SR, NULL, NULL, NULL)

r0 <- max(recs)  ## or: r0 <- dat$R0
## ssblim <- 0.2 * halibut3$ref$B0 (Dichmont 2017)
## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
## blim <- 0.3 * halibut3$ref$Bmsy[1]
## halibut3$ref$Blim <- blim   ## in terms of SSB now (compare to SSBfinal)
print(paste0("Blim/Bmsy = ", round(halibut3$ref$Blim[1] / halibut3$ref$Bmsy[1],2)))
print(paste0("Blim = ", round(halibut3$ref$Blim[1],2)))

par(mfrow=c(1,2))
plot(ssbs, recs, ty='l',
     ylim = range(0,1.2*halibut3$R0), lwd=1.5)
abline(h = halibut3$R0, lty=2, lwd = 1.5)
## abline(v = halibut3$ref$SSBlim, lwd=1.5, lty=2)

## plot recruitment deviations
n <- 5
tmp <- list()
for(i in 1:n) tmp[[i]] <-  halibut3$R0 * gen.noise(35, halibut3$sigmaR, halibut3$rhoR, halibut3$biascorR)
plot(1:35, tmp[[1]], ty='n', ylim = range(unlist(tmp)))
for(i in 1:n){
    lines(1:35, tmp[[i]], ty='b', col=i)
}


## Define Blim as 0.15 B0 according to ICES fisheries management reference
## points for category 1 and 2 stocks (2017)
## Blim = 0.2 Blim (Dichmont 2017)
## halibut3$ref$Blim <- 0.2 * halibut3$ref$B0
## print(paste0("Blim = ", round(halibut3$ref$Blim,2)))


## Set FM based on depletion level
if(estDepl) halibut3 <- est.depletion(halibut3, set, fmax = fmax, nrep = nrepDepl, method = deplmethod)


## Production curve
if(estProd){
    par(mfrow=c(1,1))
    prod <- est.productivity.stochastic(halibut3, set, 100, nf, ncores = ncores, plot = TRUE)
    abline(h=0)
    abline(v = halibut3$ref$Blim[1],lty=2)
    abline(v = round(halibut3$ref$Bmsy[1],2),lty=3)
}




if(saveFig) dev.off()



## Save all species
## -----------------------------
stocklist <- list(
    "anchovy" = anchovy,
    "haddock" = haddock,
    "halibut" = halibut,
    "anchovy2" = anchovy2,
    "haddock2" = haddock2,
    "halibut2" = halibut2,
    "anchovy3" = anchovy3,
    "haddock3" = haddock3,
    "halibut3" = halibut3)

filename <- paste0("stocklist6.RData")
if(saveDat) save(stocklist, file = filename, version = 2)

sapply(stocklist, function(x) x$species)













if(FALSE){
    stocklist[["anchovyX1"]] <- anchovyX1
    stocklist[["anchovyX2"]] <- anchovyX2
    stocklist[["anchovyX3"]] <- anchovyX3

    stocklist[["anchovyY1"]] <- anchovyY1
    stocklist[["anchovyY2"]] <- anchovyY2
    stocklist[["anchovyY3"]] <- anchovyY3

    stocklist[["anchovyZ1"]] <- anchovyZ1
    stocklist[["anchovyZ2"]] <- anchovyZ2
    stocklist[["anchovyZ3"]] <- anchovyZ3

    stocklist[["anchovyH1"]] <- anchovyH1
    stocklist[["anchovyH2"]] <- anchovyH2
    stocklist[["anchovyH3"]] <- anchovyH3
}







## Plot for reviewers comments

if(FALSE){

## NEW: anchovy: without autocorrelation:
if(FALSE){

    anchovyX <- stocklist[["anchovy"]]

    anchovyX$rhoR <- 0
    anchovyX$sigmaR <- 0.766
    anchovyX <- check.dat(anchovyX)
    set$noiseR <- c(anchovyX$sigmaR, anchovyX$rhoR, anchovyX$biascorR)

    anchovyX <- est.ref.levels.stochastic(anchovyX, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    title("sigmaR")

    fmsy <- anchovyX$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyX$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyX$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyX$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyX$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyX$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyX$ref$SSBmsy[1],2)))

    print(round(anchovyX$ref$Bmsy[1]/anchovyX$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyX$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyX$h, get.ssb0(as.numeric(t(anchovyX$M[1,]*anchovyX$Msel[[1]])), as.numeric(t(anchovyX$mat)),
                                         as.numeric(t(anchovyX$weight)), 1,
                                         (anchovyX$amax+1) * anchovyX$ns, anchovyX$ns,
                                         anchovyX$spawning, anchovyX$R0,anchovyX$indage0,
                                         season = 1)/anchovyX$R0,
                    ssbs, anchovyX$R0, anchovyX$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyX$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyX$ref$Bmsy[1]
    anchovyX$ref$Blim <- blim
    print(blim / anchovyX$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyX$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyX <- est.depletion(anchovyX, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyX1 <- anchovyX
    rm(anchovyX)

    ##

    anchovyX <- stocklist[["anchovy"]]
    anchovyX$rhoR <- 0
    anchovyX$sigmaR <- 0.5  * 0.766
    anchovyX <- check.dat(anchovyX)
    set$noiseR <- c(anchovyX$sigmaR, anchovyX$rhoR, anchovyX$biascorR)

    dev.new()
    anchovyX <- est.ref.levels.stochastic(anchovyX, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    title("0.5 * sigmaR")

    fmsy <- anchovyX$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyX$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyX$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyX$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyX$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyX$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyX$ref$SSBmsy[1],2)))

    print(round(anchovyX$ref$Bmsy[1]/anchovyX$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyX$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyX$h, get.ssb0(as.numeric(t(anchovyX$M[1,]*anchovyX$Msel[[1]])), as.numeric(t(anchovyX$mat)),
                                         as.numeric(t(anchovyX$weight)), 1,
                                         (anchovyX$amax+1) * anchovyX$ns, anchovyX$ns,
                                         anchovyX$spawning, anchovyX$R0,anchovyX$indage0,
                                         season = 1)/anchovyX$R0,
                    ssbs, anchovyX$R0, anchovyX$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyX$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyX$ref$Bmsy[1]
    anchovyX$ref$Blim <- blim
    print(blim / anchovyX$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyX$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyX <- est.depletion(anchovyX, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyX2 <- anchovyX
    rm(anchovyX)

    anchovyX <- stocklist[["anchovy"]]
    anchovyX$rhoR <- 0
    anchovyX$sigmaR <- 1.5 * 0.766
    anchovyX <- check.dat(anchovyX)
    set$noiseR <- c(anchovyX$sigmaR, anchovyX$rhoR, anchovyX$biascorR)

    dev.new()
    anchovyX <- est.ref.levels.stochastic(anchovyX, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    title("1.5 * sigmaR")

    fmsy <- anchovyX$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyX$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyX$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyX$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyX$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyX$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyX$ref$SSBmsy[1],2)))

    print(round(anchovyX$ref$Bmsy[1]/anchovyX$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyX$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyX$h, get.ssb0(as.numeric(t(anchovyX$M[1,]*anchovyX$Msel[[1]])), as.numeric(t(anchovyX$mat)),
                                         as.numeric(t(anchovyX$weight)), 1,
                                         (anchovyX$amax+1) * anchovyX$ns, anchovyX$ns,
                                         anchovyX$spawning, anchovyX$R0,anchovyX$indage0,
                                         season = 1)/anchovyX$R0,
                    ssbs, anchovyX$R0, anchovyX$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyX$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyX$ref$Bmsy[1]
    anchovyX$ref$Blim <- blim
    print(blim / anchovyX$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyX$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyX <- est.depletion(anchovyX, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyX3 <- anchovyX
    rm(anchovyX)


    ## with autocorrelation (just to doublecheck)
    ####################################33

    anchovyY <- stocklist[["anchovy"]]

    anchovyY$sigmaR <- 0.766
    anchovyY <- check.dat(anchovyY)
    set$noiseR <- c(anchovyY$sigmaR, anchovyY$rhoR, anchovyY$biascorR)

    anchovyY <- est.ref.levels.stochastic(anchovyY, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyY$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyY$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyY$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyY$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyY$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyY$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyY$ref$SSBmsy[1],2)))

    print(round(anchovyY$ref$Bmsy[1]/anchovyY$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyY$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyY$h, get.ssb0(as.numeric(t(anchovyY$M[1,]*anchovyY$Msel[[1]])), as.numeric(t(anchovyY$mat)),
                                         as.numeric(t(anchovyY$weight)), 1,
                                         (anchovyY$amax+1) * anchovyY$ns, anchovyY$ns,
                                         anchovyY$spawning, anchovyY$R0,anchovyY$indage0,
                                         season = 1)/anchovyY$R0,
                    ssbs, anchovyY$R0, anchovyY$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyY$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyY$ref$Bmsy[1]
    anchovyY$ref$Blim <- blim
    print(blim / anchovyY$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyY$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyY <- est.depletion(anchovyY, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyY1 <- anchovyY
    rm(anchovyY)

    ##

    anchovyY <- stocklist[["anchovy"]]
    anchovyY$sigmaR <- 0.5  * 0.766
    anchovyY <- check.dat(anchovyY)
    set$noiseR <- c(anchovyY$sigmaR, anchovyY$rhoR, anchovyY$biascorR)

    dev.new()
    anchovyY <- est.ref.levels.stochastic(anchovyY, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("0.5 * sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyY$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyY$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyY$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyY$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyY$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyY$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyY$ref$SSBmsy[1],2)))

    print(round(anchovyY$ref$Bmsy[1]/anchovyY$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyY$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyY$h, get.ssb0(as.numeric(t(anchovyY$M[1,]*anchovyY$Msel[[1]])), as.numeric(t(anchovyY$mat)),
                                         as.numeric(t(anchovyY$weight)), 1,
                                         (anchovyY$amax+1) * anchovyY$ns, anchovyY$ns,
                                         anchovyY$spawning, anchovyY$R0,anchovyY$indage0,
                                         season = 1)/anchovyY$R0,
                    ssbs, anchovyY$R0, anchovyY$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyY$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyY$ref$Bmsy[1]
    anchovyY$ref$Blim <- blim
    print(blim / anchovyY$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyY$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyY <- est.depletion(anchovyY, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyY2 <- anchovyY
    rm(anchovyY)

    anchovyY <- stocklist[["anchovy"]]
    anchovyY$sigmaR <- 1.5 * 0.766
    anchovyY <- check.dat(anchovyY)
    set$noiseR <- c(anchovyY$sigmaR, anchovyY$rhoR, anchovyY$biascorR)

    dev.new()
    anchovyY <- est.ref.levels.stochastic(anchovyY, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("1.5 * sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyY$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyY$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyY$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyY$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyY$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyY$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyY$ref$SSBmsy[1],2)))

    print(round(anchovyY$ref$Bmsy[1]/anchovyY$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyY$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyY$h, get.ssb0(as.numeric(t(anchovyY$M[1,]*anchovyY$Msel[[1]])), as.numeric(t(anchovyY$mat)),
                                         as.numeric(t(anchovyY$weight)), 1,
                                         (anchovyY$amax+1) * anchovyY$ns, anchovyY$ns,
                                         anchovyY$spawning, anchovyY$R0,anchovyY$indage0,
                                         season = 1)/anchovyY$R0,
                    ssbs, anchovyY$R0, anchovyY$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyY$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyY$ref$Bmsy[1]
    anchovyY$ref$Blim <- blim
    print(blim / anchovyY$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyY$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyY <- est.depletion(anchovyY, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyY3 <- anchovyY
    rm(anchovyY)



    ## with autocorrelation and BevHolt SR (h=0.75)
    ####################################

    anchovyZ <- stocklist[["anchovy"]]

    anchovyZ$sigmaR <- 0.766
    anchovyZ$SR <- "bevholt"
    anchovyZ$h <- 0.75
    anchovyZ <- check.dat(anchovyZ)
    set$noiseR <- c(anchovyZ$sigmaR, anchovyZ$rhoR, anchovyZ$biascorR)

    anchovyZ <- est.ref.levels.stochastic(anchovyZ, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyZ$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyZ$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyZ$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyZ$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyZ$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyZ$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyZ$ref$SSBmsy[1],2)))

    print(round(anchovyZ$ref$Bmsy[1]/anchovyZ$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyZ$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyZ$h, get.ssb0(as.numeric(t(anchovyZ$M[1,]*anchovyZ$Msel[[1]])), as.numeric(t(anchovyZ$mat)),
                                         as.numeric(t(anchovyZ$weight)), 1,
                                         (anchovyZ$amax+1) * anchovyZ$ns, anchovyZ$ns,
                                         anchovyZ$spawning, anchovyZ$R0,anchovyZ$indage0,
                                         season = 1)/anchovyZ$R0,
                    ssbs, anchovyZ$R0, anchovyZ$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyZ$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyZ$ref$Bmsy[1]
    anchovyZ$ref$Blim <- blim
    print(blim / anchovyZ$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyZ$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyZ <- est.depletion(anchovyZ, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyZ1 <- anchovyZ
    rm(anchovyZ)

    ##

    anchovyZ <- stocklist[["anchovy"]]
    anchovyZ$sigmaR <- 0.5  * 0.766
    anchovyZ$SR <- "bevholt"
    anchovyZ$h <- 0.75
    anchovyZ <- check.dat(anchovyZ)
    set$noiseR <- c(anchovyZ$sigmaR, anchovyZ$rhoR, anchovyZ$biascorR)

    dev.new()
    anchovyZ <- est.ref.levels.stochastic(anchovyZ, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("0.5 * sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyZ$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyZ$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyZ$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyZ$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyZ$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyZ$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyZ$ref$SSBmsy[1],2)))

    print(round(anchovyZ$ref$Bmsy[1]/anchovyZ$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyZ$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyZ$h, get.ssb0(as.numeric(t(anchovyZ$M[1,]*anchovyZ$Msel[[1]])), as.numeric(t(anchovyZ$mat)),
                                         as.numeric(t(anchovyZ$weight)), 1,
                                         (anchovyZ$amax+1) * anchovyZ$ns, anchovyZ$ns,
                                         anchovyZ$spawning, anchovyZ$R0,anchovyZ$indage0,
                                         season = 1)/anchovyZ$R0,
                    ssbs, anchovyZ$R0, anchovyZ$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyZ$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyZ$ref$Bmsy[1]
    anchovyZ$ref$Blim <- blim
    print(blim / anchovyZ$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyZ$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyZ <- est.depletion(anchovyZ, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyZ2 <- anchovyZ
    rm(anchovyZ)

    anchovyZ <- stocklist[["anchovy"]]
    anchovyZ$sigmaR <- 1.5 * 0.766
    anchovyZ$SR <- "bevholt"
    anchovyZ$h <- 0.75
    anchovyZ <- check.dat(anchovyZ)
    set$noiseR <- c(anchovyZ$sigmaR, anchovyZ$rhoR, anchovyZ$biascorR)

    dev.new()
    anchovyZ <- est.ref.levels.stochastic(anchovyZ, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("1.5 * sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyZ$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyZ$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyZ$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyZ$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyZ$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyZ$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyZ$ref$SSBmsy[1],2)))

    print(round(anchovyZ$ref$Bmsy[1]/anchovyZ$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyZ$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyZ$h, get.ssb0(as.numeric(t(anchovyZ$M[1,]*anchovyZ$Msel[[1]])), as.numeric(t(anchovyZ$mat)),
                                         as.numeric(t(anchovyZ$weight)), 1,
                                         (anchovyZ$amax+1) * anchovyZ$ns, anchovyZ$ns,
                                         anchovyZ$spawning, anchovyZ$R0,anchovyZ$indage0,
                                         season = 1)/anchovyZ$R0,
                    ssbs, anchovyZ$R0, anchovyZ$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyZ$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyZ$ref$Bmsy[1]
    anchovyZ$ref$Blim <- blim
    print(blim / anchovyZ$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyZ$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyZ <- est.depletion(anchovyZ, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyZ3 <- anchovyZ
    rm(anchovyZ)



    ## with autocorrelation and BevHolt SR (h=0.9)
    ####################################

    anchovyH <- stocklist[["anchovy"]]

    anchovyH$sigmaR <- 0.766
    anchovyH$SR <- "bevholt"
    anchovyH$h <- 0.9
    anchovyH <- check.dat(anchovyH)
    set$noiseR <- c(anchovyH$sigmaR, anchovyH$rhoR, anchovyH$biascorR)

    anchovyH <- est.ref.levels.stochastic(anchovyH, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyH$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyH$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyH$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyH$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyH$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyH$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyH$ref$SSBmsy[1],2)))

    print(round(anchovyH$ref$Bmsy[1]/anchovyH$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyH$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyH$h, get.ssb0(as.numeric(t(anchovyH$M[1,]*anchovyH$Msel[[1]])), as.numeric(t(anchovyH$mat)),
                                         as.numeric(t(anchovyH$weight)), 1,
                                         (anchovyH$amax+1) * anchovyH$ns, anchovyH$ns,
                                         anchovyH$spawning, anchovyH$R0,anchovyH$indage0,
                                         season = 1)/anchovyH$R0,
                    ssbs, anchovyH$R0, anchovyH$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyH$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyH$ref$Bmsy[1]
    anchovyH$ref$Blim <- blim
    print(blim / anchovyH$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyH$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyH <- est.depletion(anchovyH, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyH1 <- anchovyH
    rm(anchovyH)

    ##

    anchovyH <- stocklist[["anchovy"]]
    anchovyH$sigmaR <- 0.5  * 0.766
    anchovyH$SR <- "bevholt"
    anchovyH$h <- 0.9
    anchovyH <- check.dat(anchovyH)
    set$noiseR <- c(anchovyH$sigmaR, anchovyH$rhoR, anchovyH$biascorR)

    dev.new()
    anchovyH <- est.ref.levels.stochastic(anchovyH, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("0.5 * sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyH$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyH$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyH$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyH$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyH$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyH$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyH$ref$SSBmsy[1],2)))

    print(round(anchovyH$ref$Bmsy[1]/anchovyH$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyH$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyH$h, get.ssb0(as.numeric(t(anchovyH$M[1,]*anchovyH$Msel[[1]])), as.numeric(t(anchovyH$mat)),
                                         as.numeric(t(anchovyH$weight)), 1,
                                         (anchovyH$amax+1) * anchovyH$ns, anchovyH$ns,
                                         anchovyH$spawning, anchovyH$R0,anchovyH$indage0,
                                         season = 1)/anchovyH$R0,
                    ssbs, anchovyH$R0, anchovyH$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyH$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyH$ref$Bmsy[1]
    anchovyH$ref$Blim <- blim
    print(blim / anchovyH$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyH$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyH <- est.depletion(anchovyH, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyH2 <- anchovyH
    rm(anchovyH)

    anchovyH <- stocklist[["anchovy"]]
    anchovyH$sigmaR <- 1.5 * 0.766
    anchovyH$SR <- "bevholt"
    anchovyH$h <- 0.9
    anchovyH <- check.dat(anchovyH)
    set$noiseR <- c(anchovyH$sigmaR, anchovyH$rhoR, anchovyH$biascorR)

    dev.new()
    anchovyH <- est.ref.levels.stochastic(anchovyH, set, fmax = fmax,
                                          ncores=ncores, plot = TRUE)
    mtext("1.5 * sigmaR",3,-2,outer=TRUE)

    fmsy <- anchovyH$ref$Fmsy[1]

    print(paste0("Fmsy = ", round(anchovyH$ref$Fmsy[1],3)))
    print(paste0("Bmsy = ", round(anchovyH$ref$Bmsy[1],2)))
    print(paste0("MSY = ", round(anchovyH$ref$MSY[1],2)))
    print(paste0("B0 = ", round(anchovyH$ref$B0[1],2)))
    print(paste0("ESBmsy = ", round(anchovyH$ref$ESBmsy[1],2)))
    print(paste0("SSBmsy = ", round(anchovyH$ref$SSBmsy[1],2)))

    print(round(anchovyH$ref$Bmsy[1]/anchovyH$ref$B0[1]*100,2))

    ## Stock-recruitment relationship + Blim
    ssbs <- seq(0, max(anchovyH$ref$B0), length.out = 5e4)
    recs <- recfunc(anchovyH$h, get.ssb0(as.numeric(t(anchovyH$M[1,]*anchovyH$Msel[[1]])), as.numeric(t(anchovyH$mat)),
                                         as.numeric(t(anchovyH$weight)), 1,
                                         (anchovyH$amax+1) * anchovyH$ns, anchovyH$ns,
                                         anchovyH$spawning, anchovyH$R0,anchovyH$indage0,
                                         season = 1)/anchovyH$R0,
                    ssbs, anchovyH$R0, anchovyH$SR, NULL, NULL, NULL)

    r0 <- max(recs)  ## or: r0 <- dat$R0
    ## ssblim <- 0.2 * anchovyH$ref$B0 (Dichmont 2017)
    ## ssblim <- ssbs[which.min(abs(recs - 0.7 * r0))]
    blim <- 0.3 * anchovyH$ref$Bmsy[1]
    anchovyH$ref$Blim <- blim
    print(blim / anchovyH$ref$Bmsy[1])
    print(paste0("Blim = ", round(anchovyH$ref$Blim[1],2)))

    ## Set FM based on depletion level
    ## if(estDepl) anchovyH <- est.depletion(anchovyH, set, fmax = fmax, nrep = nrepDepl,
    ##                                       method = deplmethod, do.opt = FALSE)

    anchovyH3 <- anchovyH
    rm(anchovyH)

}



    anchovy2 <- anchovy
    anchovy2$ny <- 35
    anchovy2 <- check.dat(anchovy2)
    anchovy2$FM <- matrix(0.00001,nrow=anchovy2$ny,ncol=anchovy2$ns)

    set$noiseR <- c(1,0,1)

    popA <- list()
    for(i in 1:5){
        set.seed(432 + i)
        popA[[i]] <- initPop(anchovy2, set)
    }

    ## plot(popA$rec[,2], ty='b')
    ## hist(popA$rec[,2], breaks = 30)
    ## hist(popA$SSB[,2], breaks = 30)
    ## plot(popA$SSB[,2], popA$rec[,2])

    haddock2 <- haddock
    haddock2$ny <- 35
    haddock2 <- check.dat(haddock2)
    haddock2$FM <- matrix(0.00001,nrow=haddock2$ny,ncol=haddock2$ns)

    set$noiseR <- c(1,0,1)

    popH <- list()
    for(i in 1:5){
        set.seed(432 + i)
        popH[[i]] <- initPop(haddock2, set)
    }

    ## plot(pop$rec[,2], ty='b')
    ## hist(pop$rec[,2], breaks = 30)
    ## hist(pop$SSB[,2], breaks = 30)
    ## plot(pop$SSB[,2], pop$rec[,2])
    ## plot(pop$SSB[,2])

    ## inp <- pop$obs[c("timeC","obsC","timeI","obsI")]
    ## plotspict.data(inp)
    ## fit <- fit.spict(inp)
    ## plot(fit)





par(mfrow=c(2,1))
plot(popA[[1]]$rec[,2], ty='n',
     ylim = range(sapply(popA,function(x) x$rec[,2])))
for(i in 1:5){
    lines(popA[[i]]$rec[,2], col=i)
}
mtext("Anchovy",3,0.5,font=2)
plot(popH[[1]]$rec[,2], ty='n',
     ylim = range(sapply(popA,function(x) x$rec[,2])))
for(i in 1:5){
    lines(popH[[i]]$rec[,2], col=i)
}
mtext("Haddock",3,0.5,font=2)

dev.print(pdf, "exampleRecs.pdf")


}














## Validate depletion levels, production curves and reference points
if(FALSE){

    require(mse)

    load("stocklist.RData")

    set <- check.set()
    set$sigmaR <- 0.736
    set$rhoR <- 0.451
    set$sigmaF <- 0.15
    set$nyhist <- 35
    set$nysim <- 35
    set$nrep <- 100
    set$burnin <- 100

    set$sigmaR <- 0.0#736
    set$rhoR <- 0.0##451
    set$sigmaF <- 0.0#15

    i=1

    dat <- dab
    dat <- anchovy

    for(i in 1:3){

        dat <- stocklist[[i]]

        ## high depletion
        dat$depl.quant <- "Blim"
        dat$depl <- 1
        dat <- est.depletion(dat, fmax = 10)

        ## low depletion
        dat$depl.quant <- "Bmsy"
        dat$depl <- 1.5
        dat$FM <- rev(dat$FM)
        dat$Fs <- dat$FM / dat$ns
        dat <- est.depletion(dat, fmax = 10)

        ## production
        prod <- estProd(dat, set, 1e4, plot = TRUE)


        ## Reference points
        dat$ny <- 100
        dat$FM <- rep(dat$ref$Fmsy, dat$ny)
        dat$Fs <- dat$FM / dat$ns
        pops <- parallel::mclapply(1:5e2,
                                   function(x) initPop(dat, set),
                                   mc.cores=5)
        tsbs <- do.call(rbind,lapply(pops, function(x) x$TSB[,1]))
        tsbsR <- apply(tsbs, 2, quantile, prob=c(0.025,0.5,0.975))
        fms <- do.call(rbind,lapply(pops, function(x) rowSums(x$FM)))
        fmsR <- apply(fms, 2, quantile, prob=c(0.025,0.5,0.975))
        ## plot
        par(mfrow=c(2,1))
        plot(tsbs[1,], ty='n', ylim = c(0,dat$ref$B0*1.1))
        polygon(c(1:dat$ny,dat$ny:1),
                c(tsbsR[1,],rev(tsbsR[3,])),
                border = NA, col = "grey80")
        lines(1:dat$ny, tsbsR[2,], lwd=1.5)
        abline(h=dat$ref$B0, lty=2, lwd=1.5,col=11)
        abline(h=dat$ref$Bmsy, lty=2, lwd=1.5,col=11)
        abline(h=dat$ref$Blim, lty=2, lwd=1.5,col=11)
        plot(fms[1,], ty='n', ylim=c(0.8,1.2)*range(fmsR, dat$ref$Fmsy))
        polygon(c(1:dat$ny,dat$ny:1),
                c(fmsR[1,],rev(fmsR[3,])),
                border = NA, col = "grey80")
        lines(1:dat$ny, fmsR[2,], lwd=1.5)
        abline(h=dat$ref$Fmsy, lty=2, lwd=1.5, col=11)

    }

    plot(rowSums(pops[[1]]$FM))
    abline(h=dat$ref$Fmsy)



    a



    ### Figure out better depletion levels and exploitation patterns
    ## problem:same exploitation pattern and depletion level in final historical year
    ## implies different contrast and different risk level in final year
    ## e.g. halibut overfished for 35 years to meet depletion level with set exploitation pattern
    ## better:
    ## define risk level in last historical year (over 5%, ref by ICES), e.g. 40%
    ## use increasing exploitation pattern which can vary slightly between species
    ## define depletion level based on B0 (that is how it is defined)


    dat$FM <- c(seq(0.4,0.1, length.out=ny-10),rep(0.9,10))

    dat$ref
    dat$FM <- c(seq(1.246, length.out=ny))

    dat$FM <- rep(1.246,ny)
    dat$Fs <- dat$FM/dat$ns
    dat$depl.quant <- "Blim"
    dat$depl <- 1
    dat <- est.depletion(dat, fmax = 10)


    set <- check.set()
    set$sigmaR <- 0
    set$rhoR <- 0
    set$sigmaF <- 0
    tail(initPop(dat,set)$TSBfinal,1) / dat$ref$Blim
    ## if there were no noise, than median (or average) of final Depl level would be at Blim!

    pop <- initPop(dat,set)
    dat$FM

    plot(rowSums(pop$FM))
    plot(pop$TSB[,1])


    dat$ns <- 4
    dat <- check.dat(dat)


    require(mse)
    load("stocklist.RData")
    dat <- stocklist[["halibut"]]


    set <- check.set()
    set$sigmaF <- 0.15
    set$sigmaR <- 0.736
    set$rhoR <- 0.451
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0
    set$sigmaImp <- 0
    set$burnin <- 200
    set$refYears <- 300
    set$refN <- 1e2


    set$sigmaF <- 0
    set$sigmaR <- 0
    set$rhoR <- 0

    dat$FM <- rep(1.246387,dat$ny) ## dab

    dat$ny <- 200
    ##     dat$FM <- rep(0.1004029,dat$ny) ## halibut
    dat$FM <- rep(dat$ref$Fmsy,dat$ny)

    dat$FM <- rep(0.09,dat$ny)
    dat$Fs <- dat$FM/dat$ns


    dat$ns <- 4
    dat <- check.dat(dat)
    ## bevholt works!
    dat$h <- 0.79
    dat$SR <- "bevholt"
    ## hockey-stick bmsy not correct!
    ## bent hyperbola (Blim == alpha)
    dat$SR <- "average"
    dat <- estRef(dat, set, ref = "B0")
    dat$recAlpha <- 0.15 * dat$ref$B0
    dat$SR <- "bent-hyperbola"
    dat$recBeta <- 0.9
    dat$recGamma <- 0

    dat <- estRef(dat, set, fvec = seq(0,0.6,0.001))

    set$sigmaR <- 0.7
    dat$ny <- 200
    dat$FM <- rep(dat$ref$Fmsy,dat$ny)
    dat$Fs <- dat$FM/dat$ns
    tmp <- initPop(dat, set)
    ##plot(tmp$TSB[,1], ylim =range(tmp$TSB[,1],5523574))
    plot(tmp$TSB[,1], ty='l', col=4)
    abline(h=dat$ref$Bmsy)
    plot(rowSums(tmp$FM), col=4, ty='l')
    abline(h=dat$ref$Fmsy)

    dat$ref


    ## sum

    depls <- parallel::mclapply(1:1e3,
                                function(x) initPop(dat, set, out.opt = 2, depl.quant = "Blim"),
                                mc.cores=5)
    mean(unlist(depls) < 1)


    pop <- simpop(1.246, dat, set)

    tail(pop$TSB,20)

    pop2 <- initPop(dat,set)
    pop2$TSB[,1]

    dat$ref


















    set2 <- check.set()
    set2$sigmaR <- 0.736
    set2$rhoR <- 0.451
    set2$sigmaF <- 0.15
    set2$burnin <- 0

    ## set2 <- set

    dat <- haddock
    ## Reference points
    dat$ny <- 2e2
    dat$FM <- rep(dat$ref$Fmsy, dat$ny)
    ##dat$FM <- rep(0.0, dat$ny)
    dat$Fs <- dat$FM / dat$ns
    pops <- parallel::mclapply(1:5e2,
                               function(x) initPop(dat, set2),
                               mc.cores=5)
    tsbs <- do.call(rbind,lapply(pops, function(x) x$TSB[,1]))
    tsbsR <- apply(tsbs, 2, quantile, prob=c(0.025,0.5,0.975))
    fms <- do.call(rbind,lapply(pops, function(x) rowSums(x$FM)))
    fmsR <- apply(fms, 2, quantile, prob=c(0.025,0.5,0.975))
    recs <- do.call(rbind,lapply(pops, function(x) x$rec * x$errs$eR))
    recsR <- apply(recs, 2, quantile, prob=c(0.025,0.5,0.975))


    ## plot
    par(mfrow=c(3,1))
    ## R
    plot(recs[1,], ty='n', ylim = range(recsR))
    polygon(c(1:dat$ny,dat$ny:1),
            c(recsR[1,],rev(recsR[3,])),
            border = NA, col = "grey80")
    for(i in c(2,5,6,8)) lines(1:dat$ny, recs[i,], col=i)
    lines(1:dat$ny, recsR[2,], lwd=1.5)
    ## B
    plot(tsbs[1,], ty='n', ylim = c(0,dat$ref$B0*1.1))
    polygon(c(1:dat$ny,dat$ny:1),
            c(tsbsR[1,],rev(tsbsR[3,])),
            border = NA, col = "grey80")
    for(i in c(2,5,6,8)) lines(1:dat$ny, tsbs[i,], col=i)
    abline(h=dat$ref$B0, lty=2, lwd=1.5,col=11)
    abline(h=dat$ref$Bmsy, lty=2, lwd=1.5,col=11)
    abline(h=dat$ref$Blim, lty=2, lwd=1.5,col=12)
    lines(1:dat$ny, tsbsR[2,], lwd=1.5)
    ## FM
    plot(fms[1,], ty='n', ylim=c(0.8,1.2)*range(fmsR, dat$ref$Fmsy))
    polygon(c(1:dat$ny,dat$ny:1),
            c(fmsR[1,],rev(fmsR[3,])),
            border = NA, col = "grey80")
    for(i in c(2,5,6,8)) lines(1:dat$ny, fms[i,], col=i)
    abline(h=dat$ref$Fmsy, lty=2, lwd=1.5, col=11)
    lines(1:dat$ny, fmsR[2,], lwd=1.5)




    tsbs <- do.call(rbind,lapply(pops, function(x) x$TSB[,1]))
    tsbsR <- apply(tsbs, 2, mean)
    fms <- do.call(rbind,lapply(pops, function(x) rowSums(x$FM)))
    fmsR <- apply(fms, 2, mean)
    recs <- do.call(rbind,lapply(pops, function(x) x$rec * x$errs$eR))
    recsR <- apply(recs, 2, mean)




    ## plot
    par(mfrow=c(3,1))
    ## R
    plot(recs[1,], ty='n', ylim = range(recsR))
    for(i in c(2,5,6,8)) lines(1:dat$ny, recs[i,], col=i)
    lines(1:dat$ny, recsR, lwd=1.5)
    ## B
    plot(tsbs[1,], ty='n', ylim = c(0,dat$ref$B0*1.1))
    for(i in c(2,5,6,8)) lines(1:dat$ny, tsbs[i,], col=i)
    abline(h=dat$ref$B0, lty=2, lwd=1.5,col=11)
    abline(h=dat$ref$Bmsy, lty=2, lwd=1.5,col=11)
    abline(h=dat$ref$Blim, lty=2, lwd=1.5,col=12)
    lines(1:dat$ny, tsbsR, lwd=1.5)
    ## FM
    plot(fms[1,], ty='n', ylim=c(0.8,1.2)*range(fmsR, dat$ref$Fmsy))
    for(i in c(2,5,6,8)) lines(1:dat$ny, fms[i,], col=i)
    abline(h=dat$ref$Fmsy, lty=2, lwd=1.5, col=11)
    lines(1:dat$ny, fmsR, lwd=1.5)



    ## find historic FM to get decent risk in last historic year

    require(mse)

    load("stocklist.RData")


    dat <- stocklist[["anchovy"]]

    set <- check.set()
    set$sigmaF <- 0.15
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0
    set$sigmaImp <- 0
    set$burnin <- 100
    set$sigmaR <- dat$sigmaR
    set$rhoR <- dat$rhoR


    set.seed(15)
    pop <- initPop(dat,set,out.opt = 2)
    pop

    dat$depl.quant
    plot(pop$TSBfinal/dat$ref$Blim, ty='b')
    abline(h=1)

    depl


}




## OLDER SPECIES PARAMETERISATIONS

if(FALSE){
anchovy <- list(
    stock = "ane.27.8",
    species = "Anchovy", ## European Anchovy
    latin = "Engraulis encrasicolus",
    family = "Engraulidae",   ## Clupeiformes
    icesWG = "WGHANSA",
    ny = ny,
    ns = ns,
    amax = 4,           ## (1)  Fishbase = 5, Andres = 6, HANSA no data for older than 4
    lwa = 0.004799048,  ## (1)
    lwb = 3.134380952,  ## (1)
    a0 = -0.02,         ## (1)
    K = 0.89,           ## (1)
    Linf = 18.69,       ## (1)
    Lm50 = 6.48, ## 6.08,  ## (1) derived from mat at age from Andres
    Lm95 = 6.49, ## 6.19,  ## (1) derived from mat at age from Andres
    Ls50 = 5.59, ## 6.15, ## 3.31, ## (1) derived from sel at age from Andres
    Ls95 = 6.35, ##6.61, ## 4.52, ## (1) derived from sel at age from Andres
    R0 = 1e6,
##    SR = "bevholt",
    SR = "average",
    h = 0.74,        ## older: 0.75, ## (1), alternatively: 0.75, 0.9 ## assumed by Andres et al.
    bp = 0,
    fecun = 1,
    binwidth = 1,
    CVlen = 0.1,
    q = 0.05,
    depl.quant = "SSBlim",
    depl = 1,
    depl.prob = 0.25,
    FM = NA,
    sigmaR = 0.766,  ## (5)
    rhoR = 0.435,    ## (5)
    biascorR = 1     ## (5)
)
anchovy <- check.dat(anchovy)
set$noiseR <- c(anchovy$sigmaR, anchovy$rhoR, anchovy$biascorR)









## Haddock in Division 7a (Irish Sea)
## --------------------------------------

## Sources:

## 1) ICES. 2019. Working Group for the Celtic Seas Ecoregion (WGCSE). ICES
## Scientific Reports. 1:29. 1604 pp. http://doi.org/10.17895/ices.pub.4982
## page 275 pp.

## 2) ICES. 2017. Report of the Second Workshop on the Impact of Ecosystem and
## Envi- ronmental Drivers on Irish Sea Fisheries Management (WKIrish2), 26–29
## September 2016, Belfast, Northern Ireland. ICES CM 2016/BSG:02. 199 pp.

## 3) Myers, R. A., Bowen, K. G., & Barrowman, N. J. (1999). Maximum
## reproductive rate of fish at low population sizes. Canadian Journal of
## Fisheries and Aquatic Sciences, 56(12), 2404-2419.

## 4) Froese, R. and D. Pauly. Editors. 2019.FishBase. World Wide Web electronic
## publication. www.fishbase.org, ( 12/2019 )

## 5) Thorson, J. T., Jensen, O. P., & Zipkin, E. F. (2014). How variable is
## recruitment for exploited marine fishes? A hierarchical model for testing
## life history theory. Canadian Journal of Fisheries and Aquatic Sciences,
## 71(7), 973-983.

fmax <- 1.3

haddock <- list(
    stock = "har-iris",
    species = "Haddock",
    latin = "Melanogrammus aeglefinus",
    family = "Gadidae",
    icesWG = "WGCSE",
    ny = ny,
    ns = ns,
    amax = 17,      ## (4)  (20 years in FishBase, but plus group in model)
    lwa = 0.006453, ## (1)
    lwb = 3.108,    ## (1)
    a0 = -0.092,    ## (1)
    K = 0.197,      ## (2)
    Linf = 76.9,    ## (2)
    Lm50 = 22,      ## (2)
    Lm95 = 29,      ## (2)
    L50 = 13.48,    ## (1)
    L95 = 25.17,    ## (1)
    SR = "bevholt",
    h = 0.63, ## thorson ## older: 0.74,       ## (3)
    bp = 0,
    fecun = 1,
    binwidth = 1,
    CVlen = 0.1,
    q = 0.05,
    depl.quant = "SSBlim",
    depl = 1,
    depl.prob = 0.25,
    FM = c(seq(0.1,fmax, length.out=ny-10),rep(fmax,10)),
    sigmaR = 0.748, ## (5)
    rhoR = 0.404,   ## (5)
    biascorR = 1    ## (5)
)
haddock <- check.dat(haddock)
set$noiseR <- c(haddock$sigmaR, haddock$rhoR, haddock$biascorR)




## Greenland halibut
## ----------------------------------------

## Sources:

## 1) ICES. 2018. Report of the North-Western Working Group (NWWG), 26 April–3
## May, 2018, ICES HQ, Copenhagen, Denmark. ICES CM 2018/ACOM:09. 733 pp.

## 2) Jardim, E., Azevedo, M., & Brites, N. M. (2015). Harvest control rules for
## data limited stocks using length-based reference points and survey biomass
## indices. Fisheries research, 171, 12-19.

## 3) Myers, R. A., Bowen, K. G., & Barrowman, N. J. (1999). Maximum
## reproductive rate of fish at low population sizes. Canadian Journal of
## Fisheries and Aquatic Sciences, 56(12), 2404-2419.

## 4) Froese, R. and D. Pauly. Editors. 2019.FishBase. World Wide Web electronic
## publication. www.fishbase.org, ( 12/2019 )

## 5) Rickman, S.J., N.K. Dulvy, S. Jennings and J.D. Reynolds, 2000.
## Recruitment variation related to fecundity in marine fishes. Can. J. Fish.
## Aquat. Sci. 57:116-124.

## 6) ICES. 2013. Report of the Benchmark Workshop on Greenland Halibut Stocks
## (WKBUT), 26–29 November 2013, Copenhagen, Denmark. ICES CM 2013/ ACOM:44. 367
## pp.

## 7) Thorson, J. T., Jensen, O. P., & Zipkin, E. F. (2014). How variable is
## recruitment for exploited marine fishes? A hierarchical model for testing
## life history theory. Canadian Journal of Fisheries and Aquatic Sciences,
## 71(7), 973-983.

## 8) Thorson, J. T., Predicting recruitment density dependence and intrinsic
## growth rate for all fishes worldwide using a data‐integrated life‐history
## model


fmax <- 1.3

halibut <- list(
    stock = "ghl-dis",
    species = "Greenland halibut",
    latin = "Reinhardtius hippoglossoides",
    family = "Pleuronectidae",
    icesWG = "NWWG",
    ny = ny,
    ns = ns,
    amax = 27,     ## (4) (FishBase: 30 years but plus group)
    lwa = 0.00333, ## (2)
    lwb = 3.249,   ## (2)
    a0 = -0.1,     ## (2)
    K = 0.073,     ## (2)
    Linf = 120,    ## (2)
    Lm50 = 71.2,   ## (5)
    Lm95 = 81.2,   ## (5)
    L50 = 51,      ## (6)
    L95 = 58.23,   ## (6)
    SR = "bevholt",
    h = 0.76, ## (8)   ## older: ## 0.79,      ## (3)
    bp = 0,
    fecun = 1,
    q = 0.05,
    binwidth = 1,
    CVlen = 0.1,
    depl.quant = "SSBlim",
    depl = 1,
    depl.prob = 0.25,
    FM = c(seq(0.1,fmax, length.out=ny-10),rep(fmax,10)),
    sigmaR = 0.636,  ## (7)
    rhoR = 0.437,    ## (7)
    biascorR = 1     ## (7)
)
halibut <- check.dat(halibut)
set$noiseR <- c(halibut$sigmaR, halibut$rhoR, halibut$biascorR)

## Figure out selectivity
if(FALSE){
    ## selectivity (1)
    l50 <- 51
    ptarget <- 0.08 ## read from graph in (6)
    fn <- function(l95){
        pl45 <- 1 /(1 + exp(-log(19)*(45 - l50)/(l95 - l50)))
        (pl45 - ptarget)^2
    }
    optimise(fn, c(0,100))

}

}
