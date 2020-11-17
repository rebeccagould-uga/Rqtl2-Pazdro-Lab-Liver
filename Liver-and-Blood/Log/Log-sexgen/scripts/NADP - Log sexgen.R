# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - NADP

#Load in Liver QTL Mapping - Log - SexGen.Rdata
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

qtlscan_LiverNADP <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverNADP <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP"], addcovar = sexgen, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverNADP = summary(perm_LiverNADP, alpha = c(0.35, 0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP", ylim = c(0,11))
abline(h = threshold_LiverNADP, col = c("orange", "purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADP, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADP --- Chromosome 3
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 3
  coef_blup_LiverNADP_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42,55)
  plot_coefCC(x = coef_blup_LiverNADP_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 3
  peaksNADP <- find_peaks(scan1_output = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, drop = 1.5)
  start = peaksNADP[peaksNADP$chr ==  chr,"ci_lo"]
  end = peaksNADP[peaksNADP$chr == chr, "ci_hi"] 
  
  variants_LiverNADP_chr3 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverNADP_chr3 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                        chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_chr3$lod, out_snps_LiverNADP_chr3$snpinfo, main = "Liver NADP SNPs")
  
  LiverNADP_Genes_MGI_chr3 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverNADP_chr3$lod, out_snps_LiverNADP_chr3$snpinfo, drop_hilit=1.5, genes = LiverNADP_Genes_MGI_chr3, main = "Liver NADP Genes MGI")

#For Liver NADP --- Chromosome 8
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 8
  coef_blup_LiverNADP_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(25,35)
  plot_coefCC(x = coef_blup_LiverNADP_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 8
  peaksNADP <- find_peaks(scan1_output = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, drop = 1.5)
  start = peaksNADP[peaksNADP$chr ==  chr,"ci_lo"]
  end = peaksNADP[peaksNADP$chr == chr, "ci_hi"] 
  
  variants_LiverNADP_chr8 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverNADP_chr8 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_chr8$lod, out_snps_LiverNADP_chr8$snpinfo, main = "Liver NADP SNPs")
  
  LiverNADP_Genes_MGI_chr8 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverNADP_chr8$lod, out_snps_LiverNADP_chr8$snpinfo, drop_hilit=1.5, genes = LiverNADP_Genes_MGI_chr8, main = "Liver NADP Genes MGI")

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_LiverNADP <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADP$lod, out_gwas_LiverNADP$snpinfo, altcol="green4", gap=0, main = "Liver NADP GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADP_sexgen <- est_herit(pheno["logLiverNADP"], kinship_lmm, sexgen, cores = 10)
herit_LiverNADP_sex <- est_herit(pheno["logLiverNADP"], kinship_lmm, sex, cores = 10)





