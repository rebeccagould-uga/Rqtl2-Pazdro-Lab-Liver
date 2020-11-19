# R01 BUN DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - BUN

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

  qtlscan_BUN <- scan1(genoprobs = probs, pheno = pheno["zBUN"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_BUN <- scan1perm(genoprobs = probs, pheno = pheno["zBUN"], addcovar = sexgen, n_perm = 1000, cores=10)
  
#set working directory
pdf(file = "BUN QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_BUN = summary(perm_BUN, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_BUN, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for BUN", ylim = c(0,11))
  abline(h = threshold_BUN, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_BUN, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_BUN, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksBUN <- find_peaks(scan1_output = qtlscan_BUN, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_BUN, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksBUN)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_BUN, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_BUN, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksBUN <- find_peaks(scan1_output = qtlscan_BUN, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_BUN, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksBUN)
  
  
  write_xlsx(list("BUN gmap (cM)" = gmap_peaksBUN,
                  "BUN pmap (Mbp)" = peaksBUN),
             "BUN Peaks - RankZ sexgen.xlsx")
  
dev.off()  

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "BUN GWAS - RankZ sexgen.pdf")
out_gwas_BUN <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zBUN"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_BUN$lod, out_gwas_BUN$snpinfo, altcol="green4", gap=0, main = "BUN GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

  herit_BUN_sex <- est_herit(pheno["zBUN"], kinship_lmm, sex, cores = 10)
  herit_BUN_sexgen <- est_herit(pheno["zBUN"], kinship_lmm, sexgen, cores = 10)



