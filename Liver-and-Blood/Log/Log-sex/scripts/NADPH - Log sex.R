# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - NADPH

#Load in Liver QTL Mapping - Log - sex.Rdata
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

qtlscan_LiverNADPH <- scan1(genoprobs = probs, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_LiverNADPH <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADPH"], addcovar = sex, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverNADPH = summary(perm_LiverNADPH, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADPH", ylim = c(0,11))
abline(h = threshold_LiverNADPH, col = c("purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADPH, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADPH --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 2
  coef_blup_LiverNADPH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADPH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,60)
  plot_coefCC(x = coef_blup_LiverNADPH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

  #using pmap (Mbp)
  chr = 2
  peaksNADPH <- find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, prob = 0.95)
  start = peaksNADPH[peaksNADPH$chr ==  chr,"ci_lo"]
  end = peaksNADPH[peaksNADPH$chr == chr, "ci_hi"] 
  
  variants_LiverNADPH_chr2 <- query_variants(chr, 79, 84)
  out_snps_LiverNADPH_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                                  chr = chr, start = 79, end = 84, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADPH_chr2$lod, out_snps_LiverNADPH_chr2$snpinfo, main = "Liver NADPH SNPs")
  
  LiverNADPH_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = 79, end = 84)
  plot(out_snps_LiverNADPH_chr2$lod, out_snps_LiverNADPH_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADPH_Genes_MGI_chr2, main = "Liver NADPH Genes MGI")
  
#For Liver NADPH --- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 14
  coef_blup_LiverNADPH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADPH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(30,50)
  plot_coefCC(x = coef_blup_LiverNADPH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_LiverNADPH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADPH$lod, out_gwas_LiverNADPH$snpinfo, altcol="green4", gap=0, main = "Liver NADPH GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADPH_sex <- est_herit(pheno["logLiverNADPH"], kinship_lmm, sex, cores = 10)  
herit_LiverNADPH_sexgen <- est_herit(pheno["logLiverNADPH"], kinship_lmm, sexgen, cores = 10)  


