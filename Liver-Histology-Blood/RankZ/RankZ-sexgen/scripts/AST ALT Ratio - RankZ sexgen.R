# R01 GSH DO Mapping Code 
# Updated January 2021
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD QTL MAPPING - AST/ALT Ratio

#Load in Liver-Histology-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


#setwd

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

  qtlscan_ASTALTRatio <- scan1(genoprobs = probs, pheno = pheno["zASTALTRatio"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_ASTALTRatio <- scan1perm(genoprobs = probs, pheno = pheno["zASTALTRatio"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "AST-ALT-Ratio-QTL-Results-RankZ-sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_ASTALTRatio = summary(perm_ASTALTRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_ASTALTRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for AST/ALT Ratio", ylim = c(0,11))
  abline(h = threshold_ASTALTRatio, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_ASTALTRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for AST/ALT Ratio", ylim = c(0,11))
  abline(h = threshold_ASTALTRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_ASTALTRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_ASTALTRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksASTALTRatio <- find_peaks(scan1_output = qtlscan_ASTALTRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_ASTALTRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksASTALTRatio)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_ASTALTRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_ASTALTRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksASTALTRatio <- find_peaks(scan1_output = qtlscan_ASTALTRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_ASTALTRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksASTALTRatio)
  
 
  write_xlsx(list("AST ALT Ratio gmap (cM)" = gmap_peaksASTALTRatio,
                  "AST ALT Ratio pmap (Mbp)" = peaksASTALTRatio),
             "AST-ALT-Ratio-Peaks-RankZ-sexgen.xlsx")
  
####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "AST-ALT-Ratio-GWAS-RankZ-sexgen.pdf")
out_gwas_ASTALTRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zASTALTRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=2)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_ASTALTRatio$lod, out_gwas_ASTALTRatio$snpinfo, altcol="green4", gap=0, main = "AST/ALT GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_ASTALTRatio_sex <- est_herit(pheno["zASTALTRatio"], kinship_lmm, sex, cores=2)
herit_ASTALTRatio_sexgen <- est_herit(pheno["zASTALTRatio"], kinship_lmm, sexgen, cores=2)



