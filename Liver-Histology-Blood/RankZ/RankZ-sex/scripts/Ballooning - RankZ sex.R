# R01 Ballooning DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Ballooning

#Load in Liver-Histology-Blood-RankZ-Sex.Rdata
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

qtlscan_Ballooning <- scan1(genoprobs = probs, pheno = pheno["zBallooning"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_Ballooning <- scan1perm(genoprobs = probs, pheno = pheno["zBallooning"], addcovar = sex, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Ballooning QTL Results - RankZ sex.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Ballooning = summary(perm_Ballooning, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Ballooning", ylim = c(0,11))
  abline(h = threshold_Ballooning, col = c("purple", "red", "blue"), lwd = 2)

#using gmap (cM)
  find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksBallooning <- find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksBallooning)
  
#using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksBallooning <- find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksBallooning)
  

write_xlsx(list("Ballooning gmap (cM)" = gmap_peaksBallooning,
                "Ballooning pmap (Mbp)" = peaksBallooning),
                "Ballooning Peaks - RankZ sex.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#no chromosomes of interest 

dev.off()
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Ballooning GWAS - RankZ sex.pdf")
out_gwas_Ballooning <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zBallooning"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_Ballooning$lod, out_gwas_Ballooning$snpinfo, altcol="green4", gap=0, main = "Ballooning GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_Ballooning_sex <- est_herit(pheno["zBallooning"], kinship_lmm, sex, cores = 2)
herit_Ballooning_sexgen <- est_herit(pheno["zBallooning"], kinship_lmm, sexgen, cores = 2)

