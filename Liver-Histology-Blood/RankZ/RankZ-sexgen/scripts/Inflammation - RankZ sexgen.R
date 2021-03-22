# R01 DO Mapping Code 
# Updated March 2021  
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Inflammation

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_Inflammation_DO_mapping_Liver/data")

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

qtlscan_Inflammation <- scan1(genoprobs = probs, pheno = pheno["zInflammation"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_Inflammation <- scan1perm(genoprobs = probs, pheno = pheno["zInflammation"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Inflammation QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Inflammation = summary(perm_Inflammation, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Inflammation", ylim = c(0,11))
  abline(h = threshold_Inflammation, col = c("purple", "red", "blue"), lwd = 2)

#using gmap (cM)
  find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksInflammation <- find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksInflammation)
  
#using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksInflammation <- find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksInflammation)
  

write_xlsx(list("Inflammation gmap (cM)" = gmap_peaksInflammation,
                "Inflammation pmap (Mbp)" = peaksInflammation),
                "Inflammation Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Inflammation GWAS - RankZ sexgen.pdf")
out_gwas_Inflammation <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zInflammation"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_Inflammation$lod, out_gwas_Inflammation$snpinfo, altcol="green4", gap=0, main = "Inflammation GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_Inflammation_sex <- est_herit(pheno["zInflammation"], kinship_lmm, sex, cores = 2)
herit_Inflammation_sexgen <- est_herit(pheno["zInflammation"], kinship_lmm, sexgen, cores = 2)

