# R01 ALT DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - ALT

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

qtlscan_ALT <- scan1(genoprobs = probs, pheno = pheno["zALT"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_ALT <- scan1perm(genoprobs = probs, pheno = pheno["zALT"], addcovar = sexgen, n_perm = 1000, cores=10)
  
#set working directory
pdf(file = "ALT QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_ALT = summary(perm_ALT, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_ALT, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for ALT", ylim = c(0,11))
  abline(h = threshold_ALT, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_ALT, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_ALT, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksALT <- find_peaks(scan1_output = qtlscan_ALT, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_ALT, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksALT)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_ALT, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_ALT, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksALT <- find_peaks(scan1_output = qtlscan_ALT, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_ALT, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksALT)
  
  
  write_xlsx(list("ALT gmap (cM)" = gmap_peaksALT,
                  "ALT pmap (Mbp)" = peaksALT),
             "ALT Peaks - RankZ sexgen.xlsx")
  
dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "ALT GWAS - RankZ sexgen.pdf")
out_gwas_ALT <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zALT"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_ALT$lod, out_gwas_ALT$snpinfo, altcol="green4", gap=0, main = "ALT GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

  herit_ALT_sex <- est_herit(pheno["zALT"], kinship_lmm, sex, cores = 10)
  herit_ALT_sexgen <- est_herit(pheno["zALT"], kinship_lmm, sexgen, cores = 10)



