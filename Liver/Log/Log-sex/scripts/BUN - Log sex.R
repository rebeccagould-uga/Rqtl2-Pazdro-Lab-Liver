# R01 BUN DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - BUN

#Load in Liver QTL Mapping - Log - Sex.Rdata
#Run Log Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_GSH_DO_mapping_Liver/data")

library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(RSQLite)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (RSQLite)

####################################################
## Plot Genome Scans with Permutation Tests
####################################################

qtlscan_BUN <- scan1(genoprobs = probs, pheno = pheno["logBUN"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_BUN <- scan1perm(genoprobs = probs, pheno = pheno["logBUN"], addcovar = sex, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_BUN = summary(perm_BUN, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_BUN, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for BUN", ylim = c(0,11))
abline(h = threshold_BUN, col = c("purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_BUN, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_BUN, alpha = 0.1), prob = 0.95, expand2markers = FALSE)
#using pmap (cM)
find_peaks(scan1_output = qtlscan_BUN, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_BUN, alpha = 0.1), prob = 0.95, expand2markers = FALSE)

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_BUN <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logBUN"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_BUN$lod, out_gwas_BUN$snpinfo, altcol="green4", gap=0, main = "BUN GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_BUN_sex <- est_herit(pheno["logBUN"], kinship_lmm, sex, cores = 10)
herit_BUN_sexgen <- est_herit(pheno["logBUN"], kinship_lmm, sexgen, cores = 10)



