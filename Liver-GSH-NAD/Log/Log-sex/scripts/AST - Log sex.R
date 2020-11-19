# R01 AST DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - AST

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

qtlscan_AST <- scan1(genoprobs = probs, pheno = pheno["logAST"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_AST <- scan1perm(genoprobs = probs, pheno = pheno["logAST"], addcovar = sex, n_perm = 1000, cores=10)

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_AST = summary(perm_AST, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_AST, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for AST", ylim = c(0,11))
abline(h = threshold_AST, col = c("purple", "red", "blue"), lwd = 2)
#using gmap (cM)
find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_AST, alpha = 0.1), prob = 0.95, expand2markers = FALSE)
#using pmap (cM)
find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_AST, alpha = 0.1), prob = 0.95, expand2markers = FALSE)

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For AST --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 16
  coef_blup_AST_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logAST"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_AST_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(25,45)
  plot_coefCC(x = coef_blup_AST_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 16
  peaksAST <- find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$pmap, threshold = 6.5, prob = 0.95)
  start = peaksAST[peaksAST$chr ==  chr,"ci_lo"]
  end = peaksAST[peaksAST$chr == chr, "ci_hi"] 
    
  variants_AST_chr16 <- query_variants(chr, start, end)
  out_snps_AST_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logAST"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                         chr = chr, start = start, end = end, keep_all_snps = TRUE)
  plot_snpasso(out_snps_AST_chr16$lod, out_snps_AST_chr16$snpinfo, main = "AST SNPs")
    
  AST_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = start, end = end)
  plot(out_snps_AST_chr16$lod, out_snps_AST_chr16$snpinfo, drop_hilit=1.5, genes = AST_Genes_MGI_chr16, main = "AST Genes MGI")

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_AST <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logAST"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_AST$lod, out_gwas_AST$snpinfo, altcol="green4", gap=0, main = "AST GWAS", ylim = c(0,6))

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_AST_sex <- est_herit(pheno["logAST"], kinship_lmm, sex, cores = 10)
herit_AST_sexgen <- est_herit(pheno["logAST"], kinship_lmm, sexgen, cores = 10)



