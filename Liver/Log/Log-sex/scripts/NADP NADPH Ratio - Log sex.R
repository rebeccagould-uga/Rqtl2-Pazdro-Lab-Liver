# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - NADP/NADPH Ratio

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

qtlscan_LiverNADP_NADPHRatio <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_LiverNADP_NADPHRatio <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], addcovar = sex, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverNADP_NADPHRatio = summary(perm_LiverNADP_NADPHRatio, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP/NADPH Ratio", ylim = c(0,11))
abline(h = threshold_LiverNADP_NADPHRatio, col = c("purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.1), prob = 0.95, expand2markers = FALSE)
#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADP/NADPH --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_LiverNADP_NADPHRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,65)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  peaksNADP_NADPHRatio <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, prob = 0.95)
  start = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr ==  chr,"ci_lo"]
  end = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr == chr, "ci_hi"] 
  
  variants_LiverNADP_NADPHRatio_chr2 <- query_variants(chr, 90, 99)
  out_snps_LiverNADP_NADPHRatio_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                                  chr = chr, start = 90, end = 99, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPHRatio_chr2$lod, out_snps_LiverNADP_NADPHRatio_chr2$snpinfo, main = "Liver NADP/NADPH SNPs")
  
  LiverNADP_NADPHRatio_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = 90, end = 99)
  plot(out_snps_LiverNADP_NADPHRatio_chr2$lod, out_snps_LiverNADP_NADPHRatio_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPHRatio_Genes_MGI_chr2, main = "Liver NADP/NADPH Genes MGI")
  

#For Liver NADP/NADPH --- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 14
  coef_blup_LiverNADP_NADPHRatio_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(35,45)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

#For Liver NADP/NADPH Ratio
out_gwas_LiverNADP_NADPHRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADP_NADPHRatio$lod, out_gwas_LiverNADP_NADPHRatio$snpinfo, altcol="green4", gap=0, main = "Liver NADP/NADPH Ratio GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADP_NADPHRatio_sex <- est_herit(pheno["logLiverNADP_NADPHRatio"], kinship_lmm, sex, cores = 10)  
herit_LiverNADP_NADPHRatio_sexgen <- est_herit(pheno["logLiverNADP_NADPHRatio"], kinship_lmm, sexgen, cores = 10) 

