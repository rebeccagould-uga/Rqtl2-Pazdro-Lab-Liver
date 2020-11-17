# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - NADH

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

qtlscan_LiverNADH <- scan1(genoprobs = probs, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_LiverNADH <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADH"], addcovar = sex, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverNADH = summary(perm_LiverNADH, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADH", ylim = c(0,11))
abline(h = threshold_LiverNADH, col = c("purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADH, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADH --- Chromosome 14 
par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
chr = 14
coef_blup_LiverNADH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
plot_coefCC(x = coef_blup_LiverNADH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADH, main = "Liver NADH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(45,60)
plot_coefCC(x = coef_blup_LiverNADH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADH, main = "Liver NADH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 14
peaksNADH <- find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, prob = 0.95)
start = peaksNADH[peaksNADH$chr ==  chr,"ci_lo"]
end = peaksNADH[peaksNADH$chr == chr, "ci_hi"] 

variants_LiverNADH_chr14 <- query_variants(chr, start - 1, end + 1)
out_snps_LiverNADH_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADH"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                     chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
plot_snpasso(out_snps_LiverNADH_chr14$lod, out_snps_LiverNADH_chr14$snpinfo, main = "Liver NADH SNPs")

LiverNADH_Genes_MGI_chr14 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
plot(out_snps_LiverNADH_chr14$lod, out_snps_LiverNADH_chr14$snpinfo, drop_hilit=1.5, genes = LiverNADH_Genes_MGI_chr14, main = "Liver NADH Genes MGI")

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_LiverNADH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADH$lod, out_gwas_LiverNADH$snpinfo, altcol="green4", gap=0, main = "Liver NADH GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADH_sex <- est_herit(pheno["logLiverNADH"], kinship_lmm, sex, cores = 10)
herit_LiverNADH_sexgen <- est_herit(pheno["logLiverNADH"], kinship_lmm, sexgen, cores = 10)


