# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - NADH

#Load in Liver-GSH-NAD-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**

setwd("/users/becca/R01_GSH_DO_mapping_Liver/data")

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

####################################################
## Plot Genome Scans with Permutation Tests
####################################################

qtlscan_LiverNADH <- scan1(genoprobs = probs, pheno = pheno["zLiverNADH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverNADH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADH"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "NADH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADH = summary(perm_LiverNADH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADH", ylim = c(0,11))
  abline(h = threshold_LiverNADH, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksNADH <- find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksNADH)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksNADH <- find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksNADH)
  
  
  write_xlsx(list("NADH gmap (cM)" = gmap_peaksNADH,
                  "NADH pmap (Mbp)" = peaksNADH),
             "NADH Peaks - RankZ sex.xlsx")

  
####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADH --- Chromosome 14 
par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 14
  coef_blup_LiverNADH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADH, main = "Liver NADH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,60)
  plot_coefCC(x = coef_blup_LiverNADH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADH, main = "Liver NADH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "NADH GWAS - RankZ sexgen.pdf")
out_gwas_LiverNADH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADH$lod, out_gwas_LiverNADH$snpinfo, altcol="green4", gap=0, main = "Liver NADH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADH_sexgen <- est_herit(pheno["zLiverNADH"], kinship_lmm, sexgen, cores = 10)
herit_LiverNADH_sex <- est_herit(pheno["zLiverNADH"], kinship_lmm, sex, cores = 10)





