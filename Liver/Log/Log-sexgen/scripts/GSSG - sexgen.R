# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - GSSG
#Run Log Transformation and Data Prep R Script before doing this**

#Load in Liver QTL Mapping - Log - SexGen.Rdata

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

qtlscan_LiverGSSG <- scan1(genoprobs = probs, pheno = pheno["logLiverGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverGSSG <- scan1perm(genoprobs = probs, pheno = pheno["logLiverGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverGSSG = summary(perm_LiverGSSG, alpha = c(0.55, 0.35, 0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSSG", ylim = c(0,11))
abline(h = threshold_LiverGSSG, col = c("yellow", "orange", "purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSSG, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
#using pmap (cM)
find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSSG, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver GSSG --- Chromosome 1
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 1
  coef_blup_LiverGSSG_chr1 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 1
  peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, drop = 1.5)
  start = peaksGSSG[peaksGSSG$chr ==  chr,"ci_lo"]
  end = peaksGSSG[peaksGSSG$chr == chr, "ci_hi"] 
    
  variants_LiverGSSG_chr1 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverGSSG_chr1 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                         chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSSG_chr1$lod, out_snps_LiverGSSG_chr1$snpinfo, main = "Liver GSSG SNPs")
    
  LiverGSSG_Genes_MGI_chr1 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverGSSG_chr1$lod, out_snps_LiverGSSG_chr1$snpinfo, drop_hilit=1.5, genes = LiverGSSG_Genes_MGI_chr1, main = "Liver GSSG Genes MGI")


#For Liver GSSG --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_LiverGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, drop = 1.5)
  start = peaksGSSG[peaksGSSG$chr ==  chr,"ci_lo"]
  end = peaksGSSG[peaksGSSG$chr == chr, "ci_hi"] 
  
  variants_LiverGSSG_chr2 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverGSSG_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSSG_chr2$lod, out_snps_LiverGSSG_chr2$snpinfo, main = "Liver GSSG SNPs")
  
  LiverGSSG_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverGSSG_chr2$lod, out_snps_LiverGSSG_chr2$snpinfo, drop_hilit=1.5, genes = LiverGSSG_Genes_MGI_chr2, main = "Liver GSSG Genes MGI") 

#For Liver GSSG --- Chromosome 18
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 18
  coef_blup_LiverGSSG_chr18 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr18, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,55)
  plot_coefCC(x = coef_blup_LiverGSSG_chr18, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 18
  peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, drop = 1.5)
  start = peaksGSSG[peaksGSSG$chr ==  chr,"ci_lo"]
  end = peaksGSSG[peaksGSSG$chr == chr, "ci_hi"] 
  
  variants_LiverGSSG_chr18 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverGSSG_chr18 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSSG_chr18$lod, out_snps_LiverGSSG_chr18$snpinfo, main = "Liver GSSG SNPs")
  
  LiverGSSG_Genes_MGI_chr18 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverGSSG_chr18$lod, out_snps_LiverGSSG_chr18$snpinfo, drop_hilit=1.5, genes = LiverGSSG_Genes_MGI_chr18, main = "Liver GSSG Genes MGI") 
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_LiverGSSG <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverGSSG$lod, out_gwas_LiverGSSG$snpinfo, altcol="green4", gap=0, main = "Liver GSSG GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverGSSG_sexgen <- est_herit(pheno["logLiverGSSG"], kinship_lmm, sexgen, cores = 10)
herit_LiverGSSG_sex <- est_herit(pheno["logLiverGSSG"], kinship_lmm, sex, cores = 10)


