## Stoch MSE paper
## May 2020
## Tobias Mildenberger <t.k.mildenberger@gmail.com>


## Settings
##--------------------------------------------------------------------------------
savePlots <- TRUE
date_MSElistConv <- "2020-05-17"
evalyears <- 3:20
convType <- TRUE
seedvec <- c(1:10)
evalType <- paste(range(evalyears),collapse = "-")
saveLabel <- "" ## "_scenario2_4"  ## for comparing scenario 2 and 4 (2,3,4)


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
library(vioplot)


## set environment
## ---------------------------
savedir <- file.path("~/Documents/DrTokami/research/MSE/stochhcr/res",Sys.Date())
source("funcs3.R")


## load RData (from makeRDataConv.R)
## ---------------------------
## list of 3 (species), list of 8 (scenarios), list of 15 (quantities)
load(paste0("MSElistConv_Yrs",
            paste(range(evalyears),collapse = "-"),
            "_",date_MSElistConv,".RData"))
nspec <- length(MSElistConv)
nscen <- length(MSElistConv[[1]])


## check sufficicency of sample size
suffSamp <- list()
for(spec in 1:nspec){
    suffSamp[[spec]] <- list()
    for(scen in 1:nscen){
        suffSamp[[spec]][[scen]] <- check.sample.size(MSElistConv[[spec]][[scen]])
    }
}



## Derived quantities
## --------------------------
## all HCRs
hcrsSim = c("Ref","Ref",
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
##
hcrsAll <- hcrsSim[sort(unique(c(keepExtra,indexConv)))]
## MSY-P*:
msyC <- c(20,21:29)
msyA <- c(20,39:45)
## HS-P*:
hs50C <- c(15,30:38)
hs50A <- c(15,46:52)
## HS:
hsX <- rev(c(3:19))
## all in list
hcrsL <- list(hsX,msyC,hs50C,msyA,hs50A)



##
load(paste0("MSElistAll_",date_MSElistConv,".RData"))
convInfo <- list()
for(spec in 1:nspec){
    convInfo[[spec]] <- list()
    for(scen in 1:nscen){
        convInfo[[spec]][[scen]] <- MSElistAll[[spec]][[scen]]$conv
    }
}
rm(MSElistAll)
gc()



## make directory if does not exist
##---------------------------------
dir.create(file.path(savedir), showWarnings = FALSE)


## Graphs & Tables
##--------------------------------------------------------------------------------
scenarios <- paste0("S",seq(nscen))
species <- c("Anchovy","Haddock","Ling")
nreps <- length(seedvec) * 50
nproyearsSpict <- 19

## fish and fisheries:
## small figures: 80 mm (width)
## large figures: 180 mm (width)
## pdf in inches: 1 inch == 25.4 mm
toInch <- 1/25.4
## DIN A4 format aspect ratio: 1 : 1.414
180 / 1.414

## unified colour code:
## -----------------------------
## HSx : orange
colsHSx <- colorRampPalette(c("darkorange4","darkorange1"))(length(hsX))
## MSY-C : blue
colsMSYC <- colorRampPalette(c("dodgerblue4","dodgerblue1"))(length(msyC)-1)
## HS50-C : green
colsHS50C <- colorRampPalette(c("darkgreen","darkolivegreen3"))(length(hs50C)-1)
## MSY-A : aquamarine
colsMSYA <- colorRampPalette(c("aquamarine4","aquamarine1"))(length(msyA)-1)
## HS50-A : orchid
colsHS50A <- colorRampPalette(c("darkorchid4","darkorchid1"))(length(hs50A)-1)
## Ref - grey


## Question 1 (HS vs P*)
## -------------------------

## Main plot for HS vs MSY-C
filename <- paste0(savedir,"/fig_main_HSx_MSY_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 2,
          hcrBoth = 20, hcrref = 1,
          hcrs1 = hsX, hcrs2 = msyC[-1],
          cols1 = colsHSx, cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{180mm,155mm}' ", filename))


## Time plot (HS50 vs. MSY-C35 vs. MSY-C05)
filename <- paste0(savedir,"/fig_time_HSx_MSY_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
cols <- c(colsHSx[which(hsX %in% 15)],
          colsMSYC[which(msyC %in% c(23,27))])
plot.time.focus(MSElistConv, scenario = 2, hcrref=1,
                cols = cols,
                hcrs = c(15,23,27))
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))


load(paste0("MSElistConv1_3_Yrs",
            paste(range(evalyears),collapse = "-"),
            "_",date_MSElistConv,".RData"))

## Performance vs scenario (scen diff S1 - S3)
filename <- paste0(savedir,"/fig_scenDiff_HSx_MSY_scen1_3.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
plot.diff.scen(MSElistConv1_3,
               scenarios = c(1,3),
               hcrref = 1,
               hcrBoth = 20,
               hcrs1 = hsX,
               hcrs2 = msyC[-1],
               cols1 = colsHSx,
               cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))




## Question 2 (HS + P*)
## -------------------------

## HS-C vs MSY-C
filename <- paste0(savedir,"/fig_main_HS_MSY_C_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 2,
          hcrBoth = NULL, hcrref = 1,
          hcrs1 = hs50C, hcrs2 = msyC,
          cols1 = colsHS50C, cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{180mm,155mm}' ", filename))

## Time plot in SI



## Question 3 (best HCR)
## -------------------------

## HS-A vs MSY-A
filename <- paste0(savedir,"/fig_main_HS_A_C_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 2,
          hcrBoth = 15, hcrref = 1,
          hcrs1 = hs50C[-1], hcrs2 = hs50A[-1],
          cols1 = colsHS50C, cols2 = colsHS50A)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{180mm,155mm}' ", filename))


## Table with HS-X HCR that has highest yield
hcrs <- hsX
tmp <- lapply(MSElistConv, function(x) do.call(rbind,lapply(x, function(y) y$yieldTot[1,hcrs])))
tmp <- do.call(cbind,lapply(tmp, function(x) apply(x, 1, function(y) hcrsAll[hcrs][which.max(y)])))
colnames(tmp) <- c("Anchovy","Haddock","Ling")
rownames(tmp) <- c(paste0("S",1:4))
tmp
filename <- paste0(savedir,"/table_best_HSx")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

## gain/lose by percentiles C+A (HS)
filename <- paste0(savedir,"/fig_diffPerc_HS_C_A_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.diff.perc(MSElistConv,
               scen = 2,
               hcrs1 = hs50C, hcrs2 = hs50A, cols1 = colsHS50C,
               cols2 = colsHS50A)
if(savePlots) dev.off()





## Figures and Tables for SI:
## -----------------------------------------------------------------------------------
## main plot for other scenarios
filename <- paste0(savedir,"/fig_SI_main_HSx_MSY_scen1.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 1,
          hcrBoth = 20, hcrref = 1,
          hcrs1 = hsX, hcrs2 = msyC[-1],
          cols1 = colsHSx, cols2 = colsMSYC)
if(savePlots) dev.off()

filename <- paste0(savedir,"/fig_SI_main_HSx_MSY_scen3.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 3,
          hcrBoth = 20, hcrref = 1,
          hcrs1 = hsX, hcrs2 = msyC[-1],
          cols1 = colsHSx, cols2 = colsMSYC)
if(savePlots) dev.off()

filename <- paste0(savedir,"/fig_SI_main_HSx_MSY_scen4.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 4,
          hcrBoth = 20, hcrref = 1,
          hcrs1 = hsX, hcrs2 = msyC[-1],
          cols1 = colsHSx, cols2 = colsMSYC)
if(savePlots) dev.off()

filename <- paste0(savedir,"/fig_SI_main_HSx_MSY_scen5.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 5,
          hcrBoth = 20, hcrref = 1,
          hcrs1 = hsX, hcrs2 = msyC[-1],
          cols1 = colsHSx, cols2 = colsMSYC)
if(savePlots) dev.off()


## HS-A vs MSY-A
filename <- paste0(savedir,"/fig_main_HS_MSY_A_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.fig1(MSElistConv, scenario = 2,
          hcrBoth = NULL, hcrref = 1,
          hcrs1 = hs50A, hcrs2 = msyA,
          cols1 = colsHS50A, cols2 = colsMSYA)
if(savePlots) dev.off()



## Additional scenDiff plots
## -------------------------
filename <- paste0(savedir,"/fig_scenDiff_HSx_MSY_scen1_2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
## in SI (scen diff S1 - S2)
plot.diff.scen(MSElistConv1_3,
               scenarios = c(1,2),
               hcrref = 1,
               hcrBoth = 20,
               hcrs1 = hsX,
               hcrs2 = msyC[-1],
               cols1 = colsHSx,
               cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))

filename <- paste0(savedir,"/fig_scenDiff_HSx_MSY_scen2_3.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
## in SI (scen diff S2 - S3)
plot.diff.scen(MSElistConv1_3,
               scenarios = c(1,2),
               hcrref = 1,
               hcrBoth = 20,
               hcrs1 = hsX,
               hcrs2 = msyC[-1],
               cols1 = colsHSx,
               cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))

filename <- paste0(savedir,"/fig_scenDiff_HS50C_MSYC_scen1_3.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
## in SI (scen diff S1 - S3 : HS50C vs MSYC)
plot.diff.scen(MSElistConv1_3,
               scenarios = c(1,3),
               hcrref = 1,
               hcrBoth = NULL,
               hcrs1 = hs50C[-1],
               hcrs2 = msyC[-1],
               cols1 = colsHS50C,
               cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))

filename <- paste0(savedir,"/fig_scenDiff_HS50A_HS50C_scen1_3.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
## in SI (scen diff S1 - S3 : HS50C vs MSYC)
plot.diff.scen(MSElistConv1_3,
               scenarios = c(1,3),
               hcrref = 1,
               hcrBoth = 15,
               hcrs1 = hs50A[-1],
               hcrs2 = hs50C[-1],
               cols1 = colsHS50A,
               cols2 = colsHS50C)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))



load(paste0("MSElistConv2_4_Yrs",
            paste(range(evalyears),collapse = "-"),
            "_",date_MSElistConv,".RData"))

## Performance vs scenario (scen diff S2 - S4)
filename <- paste0(savedir,"/fig_scenDiff_HSx_MSY_scen2_4.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
plot.diff.scen(MSElistConv2_4,
               scenarios = c(2,4),
               hcrref = 1,
               hcrBoth = 20,
               hcrs1 = hsX,
               hcrs2 = msyC[-1],
               cols1 = colsHSx,
               cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))


load(paste0("MSElistConv2_5_Yrs",
            paste(range(evalyears),collapse = "-"),
            "_",date_MSElistConv,".RData"))

## Performance vs scenario (scen diff S2 - S5)
filename <- paste0(savedir,"/fig_scenDiff_HSx_MSY_scen2_5.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 10, width = 9)
plot.diff.scen(MSElistConv2_5,
               scenarios = c(2,5),
               hcrref = 1,
               hcrBoth = 20,
               hcrs1 = hsX,
               hcrs2 = msyC[-1],
               cols1 = colsHSx,
               cols2 = colsMSYC)
if(savePlots) dev.off()
system(paste0("pdfjam --outfile ", filename, " --papersize '{80mm,90mm}' ", filename))
## system(paste0("pdfinfo ", filename))



## Addtional percDiff plots
## -------------------------
## gain/lose by percentiles C + A (MSY)
filename <- paste0(savedir,"/fig_SI_diffPerc_MSY_C_A_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.diff.perc(MSElistConv, scen = 2,
               hcrs1 = msyC, hcrs2 = msyA,
               cols1 = colsMSYC,
               cols2 = colsMSYA)
if(savePlots) dev.off()

## gain/lose by percentiles HS
filename <- paste0(savedir,"/fig_SI_diffPerc_HSx_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9.2, width = 11)
plot.diff.perc(MSElistConv,
               scen = 2,
               hcrs1 = hsX, cols1 = colsHSx)
if(savePlots) dev.off()



## Additional time plots
## -------------------------
## Time plot (HS50-C35 vs. MSY-C35)
hcrsAll[c(23,32)]
filename <- paste0(savedir,"/fig_SI_time_HS_MSY_C_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
cols <- c(colsHS50C[which(hs50C %in% 32)],
          colsMSYC[which(msyC %in% 23)])
plot.time.focus(MSElistConv, scenario = 2, hcrref=1,
                cols = cols,
                hcrs = c(32,23))
if(savePlots) dev.off()

## Time plot (HS-C35 vs. HS-A35)
hcrsAll[c(32,48)]
filename <- paste0(savedir,"/fig_SI_time_HS_A_C_scen2.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
cols <- c(colsHS50C[which(hs50C %in% 32)],
          colsHS50A[which(hs50A %in% 48)])
plot.time.focus(MSElistConv, scenario = 2,
                hcrref = 1,
                cols = cols,
                hcrs = c(32,48))
if(savePlots) dev.off()

## time plots for all perc of on HCR type (HSx)
filename <- paste0(savedir,"/fig_SI_time_all_HSx.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
plot.time(MSElistConv, scenario = 2, hcrref=1,
          cols = colsHSx,
          hcrs = hsX)
if(savePlots) dev.off()

## time plots for all perc of on HCR type (HS-A)
filename <- paste0(savedir,"/fig_SI_time_all_HS_A.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
plot.time(MSElistConv, scenario = 2, hcrref=1,
          cols = colsHS50A,
          hcrs = hs50A)
if(savePlots) dev.off()

## time plots for all perc of on HCR type (HS-C)
filename <- paste0(savedir,"/fig_SI_time_all_HS_C.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
plot.time(MSElistConv, scenario = 2, hcrref=1,
          cols = colsHS50C,
          hcrs = hs50C)
if(savePlots) dev.off()


## time plots for all perc of on HCR type (MSY-C)
filename <- paste0(savedir,"/fig_SI_time_all_MSY_C.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
plot.time(MSElistConv, scenario = 2, hcrref=1,
          cols = colsMSYC,
          hcrs = msyC)
if(savePlots) dev.off()

## time plots for all perc of on HCR type (MSY-A)
filename <- paste0(savedir,"/fig_SI_time_all_MSY_A.pdf")
if(savePlots) cairo_pdf(filename,
                        height = 9, width = 8)
plot.time(MSElistConv, scenario = 2, hcrref=1,
          cols = colsMSYA,
          hcrs = msyA)
if(savePlots) dev.off()





## Non convergence
## -------------------------
## non-convergence for plots (across hcr)
tmp <- do.call(cbind,lapply(MSElistConv, function(x) unlist(lapply(x, function(y) y$sampleSize))))
colnames(tmp) <- species
rownames(tmp) <- scenarios
tmp
filename <- paste0(savedir,"/table_SI_sampSize_plot")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

## non-convergence for plots (across hcr) ## percent
tmp <- do.call(cbind,lapply(MSElistConv, function(x) unlist(lapply(x, function(y) round(y$sampleSize/nreps * 100,1)))))
colnames(tmp) <- species
rownames(tmp) <- scenarios
tmp
filename <- paste0(savedir,"/table_SI_sampSize_plot_percent")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

## non-convergence
conv <- lapply(convInfo, function(x) do.call(cbind,lapply(x, function(y){
    sel <- y[,indexConv,]
    apply(sel, c(2), function(x) sum(x == 2) / (nreps * nproyearsSpict) * 100)}
    )))
tmp <- do.call(cbind,lapply(conv,function(x) round(apply(x, 2, median),1)))
colnames(tmp) <- species
rownames(tmp) <- scenarios
tmp
filename <- paste0(savedir,"/table_sampSize_percent")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)


## non-convergence by hcr
tmp <- round(do.call(rbind,lapply(conv,t)),1)
colnames(tmp) <- hcrsAll[-c(1,2)]
rownames(tmp) <- c(paste0("A",1:4),paste0("H",1:4),paste0("L",1:4))
tmp
filename <- paste0(savedir,"/table_SI_sampSize_hcr_percent")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)


## non-convergence by year
conv <- lapply(convInfo, function(x) do.call(cbind,lapply(x, function(y){
    sel <- y[,indexConv,]
    apply(sel, c(3), function(x) sum(x == 2) / (nreps * length(hcrsAll[-c(1,2)])) * 100)}
    )))
tmp <- round(do.call(rbind,lapply(conv,t)),1)
colnames(tmp) <- 1:ncol(tmp)
rownames(tmp) <- c(paste0("A",1:4),paste0("H",1:4),paste0("L",1:4))
tmp
convByYear <- t(tmp)
filename <- paste0(savedir,"/table_SI_sampSize_year_percent")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)


## Rel err and SD of spict estimations over all HCRs
##--------------------------
if(savePlots) pdf(paste0(savedir,"/fig_SI_spict_sd_relErr.pdf"),
                  height = 8, width = 8)
plot.rel.err.sd(MSElistConv, hcrs = c(3:52))
if(savePlots) dev.off()


## Suffiency of sample size
##--------------------------
## Risk Y5
if(savePlots) pdf(paste0(savedir,"/fig_SI_suff_pblimY3_5.pdf"),
                  height = 8 , width = 8 * 28/18)
plot.suff.samp(suffSamp, "pblimY3_5")
if(savePlots) dev.off()
## Risk Y20
if(savePlots) pdf(paste0(savedir,"/fig_SI_suff_pblimY6_20.pdf"),
                  height = 8 , width = 8 * 28/18)
plot.suff.samp(suffSamp, "pblimY6_20")
if(savePlots) dev.off()
## Yield 3-5
if(savePlots) pdf(paste0(savedir,"/fig_SI_suff_yieldY3_5.pdf"),
                  height = 8 , width = 8 * 28/18)
plot.suff.samp(suffSamp, "yieldY3_5")
if(savePlots) dev.off()
## Yield 3-20
if(savePlots) pdf(paste0(savedir,"/fig_SI_suff_yieldY6_20.pdf"),
                  height = 8 , width = 8 * 28/18)
plot.suff.samp(suffSamp, "yieldY6_20")
if(savePlots) dev.off()
## AAV
if(savePlots) pdf(paste0(savedir,"/fig_SI_suff_aav.pdf"),
                  height = 8 , width = 8 * 28/18)
plot.suff.samp(suffSamp, "aav")
if(savePlots) dev.off()


## SPiCT convergence by year
##--------------------------
if(savePlots) pdf(paste0(savedir,"/fig_SI_spictConvYear.pdf"),
                  height = 6, width = 9)
plot.conv.year(convByYear, cutoff = min(evalyears))
if(savePlots) dev.off()


## time to recovery table in SI
##--------------------------
tmp <- t(do.call(rbind,
               lapply(MSElistConv,
                      function(y) do.call(rbind,lapply(y, function(z) round(z$timeRec[1,],2))))))
colnames(tmp) <- c("A1","A2","A3","A4",
                   "H1","H2","H3","H4",
                   "L1","L2","L3","L4")
rownames(tmp) <- hcrsAll
tmp
filename <- paste0(savedir,"/table_SI_timeRec")
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)





































## REMOVE: at some point

if(FALSE){



    ## SI:
    ## could be used for A vs C in SI: same percentiles -> different effect if used
    ## on all quantities
    plot.hs.vs.msy(MSElistConv, hcrs1 = hcrsL[[2]][1:8], hcrs2 = hcrsL[[4]])


    ## one A rule vs one C rule time plot? (with similar risk levels), e.g. A35 vs C15
    cols <- c("darkorange3","dodgerblue3")
    plot.time.focus(MSElistConv, scenario = 2, hcrref=1,
                    cols = cols,
                    hcrs = c(6,16))
    ## very similar => SI?

    ## addiitional graphs to highlight: (?)
    ## compare different species in addition
    ## compare different scenarios

    ## Scenario vs scenario?

    spec = 3
    lim <- c(0,0.2)



    quant <- "pblimY6_20"




    hcrType <- hsX
    cols <- colsHSx
    plot(MSElistConv[[spec]][[1]][[quant]][1,hcrType],
         MSElistConv[[spec]][[2]][[quant]][1,hcrType],
         ty='n', pch = 16, col = cols,
         ylim=lim, xlim = lim)
    abline(0,1,lty=2)
    ## s2 - s1
    lines(MSElistConv[[spec]][[1]][[quant]][1,hcrType],
         MSElistConv[[spec]][[2]][[quant]][1,hcrType],
         ty='b', pch = 16, col = cols)
    ## s3 - s2
    lines(MSElistConv[[spec]][[2]][[quant]][1,hcrType],
         MSElistConv[[spec]][[3]][[quant]][1,hcrType],
         ty='b', pch = 15, col = cols)
    hcrType <- msyC
    cols <- colsMSYC
    ## s2 - s1
    lines(MSElistConv[[spec]][[1]][[quant]][1,hcrType],
         MSElistConv[[spec]][[2]][[quant]][1,hcrType],
         ty='b', pch = 16, col = cols)
    ## s3 - s2
    lines(MSElistConv[[spec]][[2]][[quant]][1,hcrType],
         MSElistConv[[spec]][[3]][[quant]][1,hcrType],
         ty='b', pch = 15, col = cols)



    plot(MSElistConv[[1]][[1]]$pblimY20[1,msyC],
         MSElistConv[[1]][[2]]$pblimY20[1,msyC],
         ty='b', pch = 16, col = colsMSYC)
    abline(0,1,lty=2)

    plot(MSElistConv[[1]][[1]]$pblimY20[1,msyA],
         MSElistConv[[1]][[2]]$pblimY20[1,msyA],
         ty='b', pch = 16, col = colsMSYA)
    abline(0,1,lty=2)

    plot(MSElistConv[[1]][[1]]$pblimTot[1,msyA],
         MSElistConv[[1]][[2]]$pblimTot[1,msyA],
         ty='b', pch = 16, col = colsMSYA)
    abline(0,1,lty=2)


    plot(MSElistConv[[1]][[2]]$pblimTot[1,msyA],
         MSElistConv[[1]][[3]]$pblimTot[1,msyA],
         ty='b', pch = 16, col = colsMSYA)
    abline(0,1,lty=2)


    ## in SI
    ##--------------------------

    MSElistConv[[1]][[1]]$pblimTot == MSElistConv[[1]][[1]]$pblimY5

    MSElistConv[[1]][[1]]$pblimY5
    MSElistConv[[1]][[2]]$pblimY5

    MSElistConv[[1]][[1]]$pblimY20
a









    ## OLDER STUFF

    ## Risk vs yield trade-off intro to plot (one spec and scenario)
    ##--------------------------
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_HS_risk_yield_single.pdf"),
                      height = 10 * 28/18, width = 10)

    plot.risk.yield.single(MSElistConv, showCI = "polygon", hcrs = c(20:36))

    if(savePlots) dev.off()


    plot.risk.yield(MSElistConv, showCI = "polygon")

    plot.risk.yield.single(MSElistConv, species = 2, scenario = 1,
                           showCI = "polygon")

    ## Risk vs yield trade-off all scenarios
    ##--------------------------
    ## HS
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_HS_risk_yield.pdf"),
                      height = 10 * 28/18, width = 10)
    plot.risk.yield(MSElistConv, showCI = "polygon", hcrs = c(20:36))
    if(savePlots) dev.off()
    ## MSY
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_MSY_risk_yield.pdf"),
                      height = 10 * 28/18, width = 10)
    plot.risk.yield(MSElistConv, showCI = "polygon", hcrs = c(3:19))
    if(savePlots) dev.off()


    ## variability in yield
    ##--------------------------
    ## HS
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_HS_yield_cv.pdf"),
                      height = 10, width = 6)
    plot.yield.cv(MSElistConv, outline = FALSE, hcrs = c(20:36))
    if(savePlots) dev.off()
    ## MSY
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_MSY_yield_cv.pdf"),
                      height = 10, width = 6)
    plot.yield.cv(MSElistConv, outline = FALSE, hcrs = c(3:19))
    if(savePlots) dev.off()

    plot.yield.cv(MSElistConv, outline = FALSE)


    ## Time plot
    ##--------------------------

    cols <- c("dodgerblue2","darkorange","darkgreen","darkred")

    cols <- c(colorRampPalette(c("darkorange4","darkorange1"))(8),
              colorRampPalette(c("dodgerblue4","dodgerblue1"))(10))
    plot.time(MSElistConv, scenario = 2, hcrref=1,
              cols = cols,
              hcrs = c(hcrsL[[6]],hcrsL[[7]]))


    ## compare all HS-x rules => difference especially in the beginning
    cols <- c(colorRampPalette(c("darkred","red"))(11))
    plot.time(MSElistConv, scenario = 2, hcrref=1,
              cols = cols,
              hcrs = hcrsL[[1]])


    cols <- c(colorRampPalette(c("darkorange4","darkorange1"))(2),
              colorRampPalette(c("dodgerblue4","dodgerblue1"))(2),
              colorRampPalette(c("darkgreen","seagreen"))(2))

    plot.time(MSElistConv, scenario = 2, hcrref=1,
              cols = cols,
              hcrs = c(3,6,20,23,47,50))





    ## Differences between Scenarios - uncertainty
    ##--------------------------
    ## risk
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_scenDiff_risk.pdf"),
                      height = 8, width = 8)
    plot.risk.diff.scen(MSElistConv)
    if(savePlots) dev.off()
    ## yield
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_scenDiff_yield.pdf"),
                      height = 8, width = 8)
    plot.yield.diff.scen(MSElistConv)
    if(savePlots) dev.off()






    ## All metrics over time
    ##--------------------------
    ## HS
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_HS_time.pdf"),
                      height = 8 * 28/18, width = 8)
    plot.time(MSElistConv, hcrs = c(20:36))
    if(savePlots) dev.off()
    ## MSY
    if(savePlots) pdf(paste0(savedir,"/plot_",convType,"_Yrs",evalType,"_MSY_time.pdf"),
                      height = 8 * 28/18, width = 8)
    plot.time(MSElistConv, hcrs = c(3:19))
    if(savePlots) dev.off()

## Table with percentage canges in risk, yield, aav from scenario 1-2 and 2-3
riskY5Scen <- list()
riskY20Scen <- list()
yieldY3_5Scen <- list()
yieldY6_20Scen <- list()
aavScen <- list()
for(spec in 1:nspec){
    riskY5 <- do.call(cbind,lapply(MSElistConv[[spec]][1:3], function(x) x$pblimY3_5[1,]))
    riskY20 <- do.call(cbind,lapply(MSElistConv[[spec]][1:3], function(x) x$pblimY5_20[1,]))
    yieldY3_5 <- do.call(cbind,lapply(MSElistConv[[spec]][1:3], function(x) x$yieldY3_5[1,]))
    yieldY5_20 <- do.call(cbind,lapply(MSElistConv[[spec]][1:3], function(x) x$yieldY5_20[1,]))
    aav <- do.call(cbind,lapply(MSElistConv[[spec]][1:3], function(x) x$yieldDiff[1,]))
    ##
    riskY5.2 <- matrix(NA, nrow(riskY5), ncol(riskY5)-1)
    riskY20.2 <- matrix(NA, nrow(riskY5), ncol(riskY5)-1)
    yieldY3_5.2 <- matrix(NA, nrow(riskY5), ncol(riskY5)-1)
    yieldY5_20.2 <- matrix(NA, nrow(riskY5), ncol(riskY5)-1)
    aav.2 <- matrix(NA, nrow(riskY5), ncol(riskY5)-1)
    for(scen in 2:3){
        riskY5.2[,scen-1] <- round((riskY5[,scen] - riskY5[,scen-1]) / riskY5[,scen-1] * 100,1)
        riskY20.2[,scen-1] <- round((riskY20[,scen] - riskY20[,scen-1]) / riskY20[,scen-1] * 100,1)
        yieldY3_5.2[,scen-1] <- round((yieldY3_5[,scen] - yieldY3_5[,scen-1]) / yieldY3_5[,scen-1] * 100,1)
        yieldY5_20.2[,scen-1] <- round((yieldY5_20[,scen] - yieldY5_20[,scen-1]) / yieldY5_20[,scen-1] * 100,1)
        aav.2[,scen-1] <- round((aav[,scen] - aav[,scen-1]) / aav[,scen-1] * 100,1)
    }
    riskY5Scen[[spec]] <- riskY5.2
    riskY20Scen[[spec]] <- riskY20.2
    yieldY3_5Scen[[spec]] <- yieldY3_5.2
    yieldY5_20Scen[[spec]] <- yieldY5_20.2
    aavScen[[spec]] <- aav.2
}

## HCRs in main text:
hcrs <- c(1, 20, 23, 41, 10, 15, 32, 48)

quant <- "RiskY5"
tmp <- cbind(riskY5Scen[[1]][hcrs,],
             riskY5Scen[[2]][hcrs,],
             riskY5Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_scenDiff_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "RiskY20"
tmp <- cbind(riskY20Scen[[1]][hcrs,],
             riskY20Scen[[2]][hcrs,],
             riskY20Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_scenDiff_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "YieldY3_5"
tmp <- cbind(yieldY3_5Scen[[1]][hcrs,],
             yieldY3_5Scen[[2]][hcrs,],
             yieldY3_5Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_scenDiff_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "YieldY5_20"
tmp <- cbind(yieldY5_20Scen[[1]][hcrs,],
             yieldY5_20Scen[[2]][hcrs,],
             yieldY5_20Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant,"_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_scenDiff_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "AAV"
tmp <- cbind(aavScen[[1]][hcrs,],
             aavScen[[2]][hcrs,],
             aavScen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_scenDiff_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

## diff scenarios for all HCRs
hcrs <- seq(hcrsAll)

quant <- "RiskY5"
tmp <- cbind(riskY5Scen[[1]][hcrs,],
             riskY5Scen[[2]][hcrs,],
             riskY5Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_SI_scenDiff_all_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "RiskY20"
tmp <- cbind(riskY20Scen[[1]][hcrs,],
             riskY20Scen[[2]][hcrs,],
             riskY20Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_SI_scenDiff_all_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "YieldY3_5"
tmp <- cbind(yieldY3_5Scen[[1]][hcrs,],
             yieldY3_5Scen[[2]][hcrs,],
             yieldY3_5Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_SI_scenDiff_all_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "YieldY5_20-"
tmp <- cbind(yieldY5_20Scen[[1]][hcrs,],
             yieldY5_20Scen[[2]][hcrs,],
             yieldY5_20Scen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant,"_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_SI_scenDiff_all_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)

quant <- "AAV"
tmp <- cbind(aavScen[[1]][hcrs,],
             aavScen[[2]][hcrs,],
             aavScen[[3]][hcrs,])
rownames(tmp) <- hcrsAll[hcrs]
colnames(tmp) <- paste0(quant, "_", rep(c("S1-S2","S2-S3"),3))
tmp
filename <- paste0(savedir,"/table_SI_scenDiff_all_",quant)
capture.output(xtable(tmp),
               file=paste0(filename,".tex"),
               append=FALSE)




}
