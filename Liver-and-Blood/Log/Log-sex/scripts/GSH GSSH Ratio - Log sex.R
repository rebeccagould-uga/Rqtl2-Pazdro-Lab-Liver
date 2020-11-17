# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - GSH GSSG Ratio

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

#for Liver GSH/GSSG Ratio
qtlscan_LiverGSH_GSSGRatio<- scan1(genoprobs = probs, pheno = pheno["logLiverGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_LiverGSH_GSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["logLiverGSH_GSSGRatio"], addcovar = sex, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverGSH_GSSGRatio = summary(perm_LiverGSH_GSSGRatio, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH/GSSG Ratio", ylim = c(0,11))
abline(h = threshold_LiverGSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver GSH/GSSG --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
  chr = 16
  coef_blup_LiverGSH_GSSGRatio_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,12)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
#using pmap (Mbp)
  chr = 16
  peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, prob = 0.95)
  start = peaksGSH_GSSGRatio[peaksGSH_GSSGRatio$chr ==  chr,"ci_lo"]
  end = peaksGSH_GSSGRatio[peaksGSH_GSSGRatio$chr == chr, "ci_hi"] 

  variants_LiverGSH_GSSGRatio_chr16 <- query_variants(chr, 5, 12)
  out_snps_LiverGSH_GSSGRatio_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                       chr = chr, start = 5, end = 12, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_GSSGRatio_chr16$lod, out_snps_LiverGSH_GSSGRatio_chr16$snpinfo, main = "Liver GSH/GSSG SNPs")
  
  LiverGSH_GSSGRatio_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 5, end = 12)
  plot(out_snps_LiverGSH_GSSGRatio_chr16$lod, out_snps_LiverGSH_GSSGRatio_chr16$snpinfo, drop_hilit=1.5, genes = LiverGSH_GSSGRatio_Genes_MGI_chr16, main = "Liver GSH/GSSG Genes MGI")

#For Liver GSH/GSSG --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 11
  coef_blup_LiverGSH_GSSGRatio_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,8)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 11
  peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold=6, drop = 1.5)
  start = peaksGSH_GSSGRatio[peaksGSH_GSSGRatio$chr ==  chr,"ci_lo"]
  end = peaksGSH_GSSGRatio[peaksGSH_GSSGRatio$chr == chr, "ci_hi"] 
  
  variants_LiverGSH_GSSGRatio_chr11 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverGSH_GSSGRatio_chr11 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                                 chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_GSSGRatio_chr11$lod, out_snps_LiverGSH_GSSGRatio_chr11$snpinfo, main = "Liver GSH/GSSG SNPs")
  
  LiverGSH_GSSGRatio_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverGSH_GSSGRatio_chr11$lod, out_snps_LiverGSH_GSSGRatio_chr11$snpinfo, drop_hilit=1.5, genes = LiverGSH_GSSGRatio_Genes_MGI_chr11, main = "Liver GSH/GSSG Genes MGI")  

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_LiverGSH_GSSGRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverGSH_GSSGRatio$lod, out_gwas_LiverGSH_GSSGRatio$snpinfo, altcol="green4", gap=0, main = "Liver GSH/GSSG Ratio GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverGSH_GSSGRatio_sex <- est_herit(pheno["logLiverGSH_GSSGRatio"], kinship_lmm, sex, cores = 10)
herit_LiverGSH_GSSGRatio_sexgen <- est_herit(pheno["logLiverGSH_GSSGRatio"], kinship_lmm, sexgen, cores = 10)

