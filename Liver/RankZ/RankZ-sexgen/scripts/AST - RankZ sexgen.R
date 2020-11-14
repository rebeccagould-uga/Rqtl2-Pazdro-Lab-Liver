# R01 AST DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - AST

#Load in Liver QTL Mapping - RankZ - sexgen.Rdata
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

  qtlscan_AST <- scan1(genoprobs = probs, pheno = pheno["zAST"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_AST <- scan1perm(genoprobs = probs, pheno = pheno["zAST"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "AST QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_AST = summary(perm_AST, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_AST, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for AST", ylim = c(0,11))
  abline(h = threshold_AST, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_AST, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksAST <- find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_AST, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksAST)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_AST, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksAST <- find_peaks(scan1_output = qtlscan_AST, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_AST, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksAST)
  
  
  write_xlsx(list("AST gmap (cM)" = gmap_peaksAST,
                  "AST pmap (Mbp)" = peaksAST),
             "AST Peaks - RankZ sexgen.xlsx")
  
####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For AST --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_AST_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_AST_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,20)
  plot_coefCC(x = coef_blup_AST_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  #could use ci_lo and ci_hi, but in this case I want a specific chr 2 position
  #start = peaksAST[peaksAST$chr ==  chr,"ci_lo"]
  #end = peaksAST[peaksAST$chr == chr, "ci_hi"] 
  
  pander(peaksAST)
  #based on peaksAST, peak of interest is ~12 Mbp
  variants_AST_chr2 <- query_variants(chr, 10, 13.5)
  out_snps_AST_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                  chr = chr, start = 10, end = 13.5, keep_all_snps = TRUE)
  plot_snpasso(out_snps_AST_chr2$lod, out_snps_AST_chr2$snpinfo, main = "AST SNPs")
  
  AST_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = 10, end = 13.5)
  plot(out_snps_AST_chr2$lod, out_snps_AST_chr2$snpinfo, drop_hilit=1.5, genes = AST_Genes_MGI_chr2, main = "AST Genes MGI")
  
#For AST --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 16
  coef_blup_AST_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_AST_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(25,45)
  plot_coefCC(x = coef_blup_AST_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 16
  #could also specify a Mbp location on a specific chromosome, but will use ci_lo and ci_hi
  start = peaksAST[peaksAST$chr ==  chr,"ci_lo"]
  end = peaksAST[peaksAST$chr == chr, "ci_hi"] 
  
  pander(peaksAST)
  #based on peaksAST, peak of interest is ~57 Mbp
  variants_AST_chr16 <- query_variants(chr, start -1, end + 1)
  out_snps_AST_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                         chr = chr, start = start - 1, end = end + 1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_AST_chr16$lod, out_snps_AST_chr16$snpinfo, main = "AST SNPs")
    
  AST_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = start - 1, end = end + 1)
  plot(out_snps_AST_chr16$lod, out_snps_AST_chr16$snpinfo, drop_hilit=1.5, genes = AST_Genes_MGI_chr16, main = "AST Genes MGI")

dev.off()

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "AST GWAS - RankZ sexgen.pdf")
out_gwas_AST <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zAST"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_AST$lod, out_gwas_AST$snpinfo, altcol="green4", gap=0, main = "AST GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_AST_sex <- est_herit(pheno["zAST"], kinship_lmm, sex, cores = 10)
herit_AST_sexgen <- est_herit(pheno["zAST"], kinship_lmm, sexgen, cores = 10)



