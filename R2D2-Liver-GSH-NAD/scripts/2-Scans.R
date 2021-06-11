# R2D2 Project
# Created by Becca Gould
# Updated May 2021

#### Phenotypes ####
# Chr 2 (WSB): GSH, GSSG, Total GSH (omitted), Eh, NADP, NADP/NADPH
# Chr 14 (NOD): NADPH, NADP/NADPH

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


#Load in R2D2-Liver-GSH-NAD-RankZ.Rdata
#Run 1-R2D2-Setup.R prior to this script


pdf(file = "Pheno-Scans.pdf")


########################################################################################################
## 
## Liver GSH (WSB - Chromosome 2)
##
########################################################################################################

  # QTL scans - sex
  #qtlscan_LiverGSH_sex <- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sex, cores=2)
  #perm_LiverGSH_sex <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = sex, n_perm = 1000, cores=10)
  
    threshold_LiverGSH_sex = summary(perm_LiverGSH_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex", ylim = c(0,11))
    abline(h = threshold_LiverGSH_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver GSH - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksGSH_sex <- find_peaks(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksGSH_sex <- find_peaks(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  
# QTL scans - sex + gen
  #qtlscan_LiverGSH_sexgen <- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  #perm_LiverGSH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_LiverGSH_sexgen = summary(perm_LiverGSH_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex + gen", ylim = c(0,11))
    abline(h = threshold_LiverGSH_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver GSH - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)


########################################################################################################
## 
## Liver GSSG (WSB - Chromosome 2)
##
########################################################################################################

  # QTL scans - sex
    #qtlscan_LiverGSSG_sex <- scan1(genoprobs = probs, pheno = pheno["zLiverGSSG"], kinship = kinship_loco, addcovar = sex, cores=2)
    #perm_LiverGSSG_sex <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSSG"], addcovar = sex, n_perm = 1000, cores=10)

    threshold_LiverGSSG_sex = summary(perm_LiverGSSG_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSSG - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSSG - sex", ylim = c(0,11))
    abline(h = threshold_LiverGSSG_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver GSSG - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksGSSG_sex <- find_peaks(scan1_output = qtlscan_LiverGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSSG_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksGSSG_sex <- find_peaks(scan1_output = qtlscan_LiverGSSG_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSSG_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
    # QTL scans - sex + gen
    #qtlscan_LiverGSSG_sexgen <- scan1(genoprobs = probs, pheno = pheno["zLiverGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=2)
    #perm_LiverGSSG_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)
    
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_LiverGSSG_sexgen = summary(perm_LiverGSSG_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSSG - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSSG - sex + gen", ylim = c(0,11))
    abline(h = threshold_LiverGSSG_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver GSSG - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksGSSG_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSSG_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksGSSG_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSSG_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSSG_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

     
########################################################################################################
## 
## Liver Eh (WSB - Chromosome 2)
##
########################################################################################################
    
    # QTL scans - sex
    #qtlscan_LiverEh_sex <- scan1(genoprobs = probs, pheno = pheno["zLiverEh"], kinship = kinship_loco, addcovar = sex, cores=2)
    #perm_LiverEh_sex <- scan1perm(genoprobs = probs, pheno = pheno["zLiverEh"], addcovar = sex, n_perm = 1000, cores=10)
    
    threshold_LiverEh_sex = summary(perm_LiverEh_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverEh_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver Eh - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverEh_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver Eh - sex", ylim = c(0,11))
    abline(h = threshold_LiverEh_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverEh_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver Eh - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksEh_sex <- find_peaks(scan1_output = qtlscan_LiverEh_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverEh_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksEh_sex <- find_peaks(scan1_output = qtlscan_LiverEh_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverEh_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
    # QTL scans - sex + gen
    #qtlscan_LiverEh_sexgen <- scan1(genoprobs = probs, pheno = pheno["zLiverEh"], kinship = kinship_loco, addcovar = sexgen, cores=2)
    #perm_LiverEh_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zLiverEh"], addcovar = sexgen, n_perm = 1000, cores=10)
    
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_LiverEh_sexgen = summary(perm_LiverEh_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverEh_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver Eh - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverEh_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver Eh - sex + gen", ylim = c(0,11))
    abline(h = threshold_LiverEh_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverEh_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver Eh - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksEh_sexgen <- find_peaks(scan1_output = qtlscan_LiverEh_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverEh_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksEh_sexgen <- find_peaks(scan1_output = qtlscan_LiverEh_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverEh_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)


########################################################################################################
## 
## Liver NADP (WSB - Chromosome 2)
##
########################################################################################################
    
    # QTL scans - sex
    #qtlscan_LiverNADP_sex <- scan1(genoprobs = probs, pheno = pheno["zNADP"], kinship = kinship_loco, addcovar = sex, cores=2)
    #perm_NADP_sex <- scan1perm(genoprobs = probs, pheno = pheno["zNADP"], addcovar = sex, n_perm = 1000, cores=10)
    
    threshold_NADP_sex = summary(perm_NADP_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP - sex", ylim = c(0,11))
    abline(h = threshold_NADP_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver NADP - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADP_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADP_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADP_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADP_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
    # QTL scans - sex + gen
    #qtlscan_LiverNADP_sexgen <- scan1(genoprobs = probs, pheno = pheno["zNADP"], kinship = kinship_loco, addcovar = sexgen, cores=2)
    #perm_NADP_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zNADP"], addcovar = sexgen, n_perm = 1000, cores=10)
    
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_NADP_sexgen = summary(perm_NADP_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP - sex + gen", ylim = c(0,11))
    abline(h = threshold_NADP_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver NADP - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADP_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADP_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADP_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADP_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
########################################################################################################
## 
## Liver NADP/NADPH (WSB - Chromosome 2)
##
########################################################################################################
    
    # QTL scans - sex
    #qtlscan_LiverNADP_NADPHRatio_sex <- scan1(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, cores=2)
    #perm_NADP_NADPHRatio_sex <- scan1perm(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], addcovar = sex, n_perm = 1000, cores=10)
    
    threshold_NADP_NADPHRatio_sex = summary(perm_NADP_NADPHRatio_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex", ylim = c(0,11))
    abline(h = threshold_NADP_NADPHRatio_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver NADP_NADPHRatio - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADP_NADPHRatio_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADP_NADPHRatio_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    print(gmap_peaksNADP_NADPHRatio_sex)
    #using pmap (Mbp)
    peaksNADP_NADPHRatio_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADP_NADPHRatio_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    print(peaksNADP_NADPHRatio_sex)
    
    
    # QTL scans - sex + gen
    #qtlscan_LiverNADP_NADPHRatio_sexgen <- scan1(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sexgen, cores=2)
    #perm_NADP_NADPHRatio_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
    
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_NADP_NADPHRatio_sexgen = summary(perm_NADP_NADPHRatio_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex + gen", ylim = c(0,11))
    abline(h = threshold_NADP_NADPHRatio_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver NADP_NADPHRatio - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADP_NADPHRatio_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADP_NADPHRatio_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADP_NADPHRatio_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADP_NADPHRatio_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
########################################################################################################
## 
## Liver NADPH (NOD - Chromosome 14)
##
########################################################################################################
    
    # QTL scans - sex
    #qtlscan_LiverNADPH_sex <- scan1(genoprobs = probs, pheno = pheno["zNADPH"], kinship = kinship_loco, addcovar = sex, cores=2)
    #perm_NADPH_sex <- scan1perm(genoprobs = probs, pheno = pheno["zNADPH"], addcovar = sex, n_perm = 1000, cores=10)
    
    threshold_NADPH_sex = summary(perm_NADPH_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADPH - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADPH - sex", ylim = c(0,11))
    abline(h = threshold_NADPH_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 14, main = "Genome Scan for Liver NADPH - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADPH_sex <- find_peaks(scan1_output = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADPH_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADPH_sex <- find_peaks(scan1_output = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADPH_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
    # QTL scans - sex + gen
    #qtlscan_LiverNADPH_sexgen <- scan1(genoprobs = probs, pheno = pheno["zNADPH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
    #perm_NADPH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zNADPH"], addcovar = sexgen, n_perm = 1000, cores=10)
    
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_NADPH_sexgen = summary(perm_NADPH_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADPH - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADPH - sex + gen", ylim = c(0,11))
    abline(h = threshold_NADPH_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 14, main = "Genome Scan for Liver NADPH - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADPH_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADPH_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADPH_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADPH_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

    
########################################################################################################
## 
## Liver NADP/NADPH (NOD - Chromosome 14)
##
########################################################################################################
    
    # QTL scans - sex
    #qtlscan_LiverNADP_NADPHRatio_sex <- scan1(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, cores=2)
    #perm_NADP_NADPHRatio_sex <- scan1perm(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], addcovar = sex, n_perm = 1000, cores=10)
    
    threshold_NADP_NADPHRatio_sex = summary(perm_NADP_NADPHRatio_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex", ylim = c(0,11))
    abline(h = threshold_NADP_NADPHRatio_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 14, main = "Genome Scan for Liver NADP_NADPHRatio - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADP_NADPHRatio_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADP_NADPHRatio_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADP_NADPHRatio_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADP_NADPHRatio_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    
    # QTL scans - sex + gen
    #qtlscan_LiverNADP_NADPHRatio_sexgen <- scan1(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sexgen, cores=2)
    #perm_NADP_NADPHRatio_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zNADP_NADPHRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
    
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_NADP_NADPHRatio_sexgen = summary(perm_NADP_NADPHRatio_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver NADP_NADPHRatio - sex + gen", ylim = c(0,11))
    abline(h = threshold_NADP_NADPHRatio_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 14, main = "Genome Scan for Liver NADP_NADPHRatio - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksNADP_NADPHRatio_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_NADP_NADPHRatio_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    #using pmap (Mbp)
    peaksNADP_NADPHRatio_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_NADP_NADPHRatio_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
    
    
dev.off()

