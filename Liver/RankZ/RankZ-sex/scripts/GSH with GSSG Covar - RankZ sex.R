# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - RankZ TRANSFORMATION AND DATA PREP

#Make a folder under users (for me, my user is "becca") and title it based on your project. For mine, it's R01_GSH_DO_mapping. Then make a data, results, scripts, and docs folder.

setwd("/users/becca/R01_GSH_DO_mapping_Liver/data")

#I can use "~" instead of "/users/becca/" everytime as it represents my home base

#load the command line tools - see https://github.com/Rdatatable/data.table/wiki/Installation for more information - must do every time you open up the Rproject!
library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (pander)
library (writexl)
library (RSQLite)

#making GSSG a covariate along with sex and sexgen
GSSGcovar = model.matrix (~ sex + zLiverGSSG, data = pheno)[,-1]

#digging more into the GSH/GSSG peak
qtlscan_LiverGSH_GSSGcovar<- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = GSSGcovar, cores=10)
perm_LiverGSH_GSSGcovar <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = GSSGcovar, n_perm = 1000, cores=10)

#set working directory
pdf(file = "GSH GSSG Covar QTL Results - RankZ sex.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH_GSSGcovar = summary(perm_LiverGSH_GSSGcovar, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSH_GSSGcovar, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH with GSSG Covar", ylim = c(0,11))
  abline(h = threshold_LiverGSH_GSSGcovar, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverGSH_GSSGcovar, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_GSSGcovar, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSH_GSSGCovar <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGcovar, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_GSSGcovar, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSH_GSSGCovar)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverGSH_GSSGcovar, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_GSSGcovar, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSH_GSSGCovar <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGcovar, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_GSSGcovar, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSH_GSSGCovar)
  
  
  write_xlsx(list("GSH GSSG Covar gmap (cM)" = gmap_peaksGSH_GSSGCovar,
                  "GSH GSSG Covar pmap (Mbp)" = peaksGSH_GSSGCovar),
             "GSH GSSG Covar Peaks - RankZ sex.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver GSH Plot with GSSG Covar --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_LiverGSH_GSSGCovar_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = GSSGcovar, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGCovar_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH with GSSG Covar BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,65)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGCovar_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH with GSSG Covar BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver GSH Plot with GSSG Covar --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 16
  coef_blup_LiverGSH_GSSGCovar_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = GSSGcovar, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGCovar_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH with GSSG Covar BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,11)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGCovar_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH with GSSG Covar BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


dev.off()





