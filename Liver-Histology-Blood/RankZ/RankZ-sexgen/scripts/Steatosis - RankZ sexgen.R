# R01 Steatosis DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Steatosis

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_Steatosis_DO_mapping_Liver/data")

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

qtlscan_Steatosis <- scan1(genoprobs = probs, pheno = pheno["zSteatosis"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_Steatosis <- scan1perm(genoprobs = probs, pheno = pheno["zSteatosis"], addcovar = sexgen, n_perm = 1000, cores=10)

#X permutation test
Xcovar = get_x_covar(R01_GSH_DO_QTLdata)
perm_strata = mat2strata(Xcovar)
perm_X_Steatosis <- scan1perm(genoprobs = probs, pheno = pheno["zSteatosis"], addcovar = sexgen, n_perm = 1000, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(R01_GSH_DO_QTLdata$gmap), cores=10)

summary(perm_X_Steatosis, alpha=c(0.2, 0.1, 0.05))
#Autosome LOD thresholds (1000 permutations)
#zSteatosis
#0.2        7.32
#0.1        7.82
#0.05       8.44

#X chromosome LOD thresholds (18090 permutations)
#zSteatosis
#0.2        7.73
#0.1        8.41
#0.05       8.83

#set working directory
pdf(file = "Steatosis QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Steatosis = summary(perm_Steatosis, alpha = c(0.2, 0.1, 0.05))
  threshold_X_Steatosis = summary(perm_X_Steatosis, alpha = c(0.2, 0.1, 0.05))
  
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Steatosis", ylim = c(0,11))
  abline(h = threshold_Steatosis, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Steatosis", ylim = c(0,11))
  abline(h = threshold_Steatosis, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Steatosis (X Chrom)", ylim = c(0,11))
  abline(h = c(7.73, 8.41, 8.83), col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Steatosis (X Chrom)", ylim = c(0,11))
  abline(h = c(7.73, 8.41, 8.83), col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #plotting separate autosome versus X axis significance thresholds
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Steatosis (Autosome vs X)", ylim = c(0,11))
  segments(x0 = 0, y0 = threshold_X_Steatosis$A, x1 = 1695, y1 =   threshold_X_Steatosis$A, col = c("purple", "red", "blue"), lwd = 2)
  segments(x0 = 1695, y0 = threshold_X_Steatosis$X, x1 = 2000, y1 = threshold_X_Steatosis$X, col = c("purple", "red", "blue"), lwd = 2)
  
  #plotting separate autosome versus X axis significance thresholds
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Steatosis (Autosome vs X)", ylim = c(0,11))
  segments(x0 = 0, y0 = threshold_X_Steatosis$A, x1 = 1695, y1 =   threshold_X_Steatosis$A, col = c("purple", "red", "blue"), "dashed", lwd = 2)
  segments(x0 = 1695, y0 = threshold_X_Steatosis$X, x1 = 2000, y1 = threshold_X_Steatosis$X, col = c("purple", "red", "blue"), "dashed", lwd = 2)
  
  
#using gmap (cM)
  find_peaks(scan1_output = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Steatosis, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksSteatosis <- find_peaks(scan1_output = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Steatosis, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksSteatosis)
  
#using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Steatosis, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksSteatosis <- find_peaks(scan1_output = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Steatosis, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksSteatosis)
  

write_xlsx(list("Steatosis gmap (cM)" = gmap_peaksSteatosis,
                "Steatosis pmap (Mbp)" = peaksSteatosis),
                "Steatosis Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Steatosis --- Chromosome 18
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 18
  coef_blup_Steatosis_chr18 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_Steatosis_chr18, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,20)
  plot_coefCC(x = coef_blup_Steatosis_chr18, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 18
  #could use ci_lo and ci_hi, but in this case I want a specific chr 2 position
  #start = peaksSteatosis[peaksSteatosis$chr ==  chr,"ci_lo"]
  #end = peaksSteatosis[peaksSteatosis$chr == chr, "ci_hi"] 
  
  pander(peaksSteatosis)
  #based on peaksAST, peak of interest is ~12 Mbp
  variants_Steatosis_chr18 <- query_variants(chr, 16, 18)
  out_snps_Steatosis_chr18 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                 chr = chr, start = 16, end = 18, keep_all_snps = TRUE)
  plot_snpasso(out_snps_Steatosis_chr18$lod, out_snps_Steatosis_chr18$snpinfo, main = "Steatosis SNPs")
  
  Steatosis_Genes_MGI_chr18 <- query_genes_mgi(chr = chr, start = 16, end = 18)
  plot(out_snps_Steatosis_chr18$lod, out_snps_Steatosis_chr18$snpinfo, drop_hilit=1.5, genes = Steatosis_Genes_MGI_chr18, main = "Steatosis Genes MGI")
  
  plot_genes(Steatosis_Genes_MGI_chr18, main = "Steatosis Genes MGI")
  
  
#For Steatosis --- Chromosome X
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "X"
  coef_blup_Steatosis_chrX <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_Steatosis_chrX, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(10,35)
  plot_coefCC(x = coef_blup_Steatosis_chrX, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "X"
  #could use ci_lo and ci_hi, but in this case I want a specific chr X position
  #start = peaksSteatosis[peaksSteatosis$chr ==  chr,"ci_lo"]
  #end = peaksSteatosis[peaksSteatosis$chr == chr, "ci_hi"] 
  
  pander(peaksSteatosis)
  variants_Steatosis_chrX <- query_variants(chr, 46, 53)
  out_snps_Steatosis_chrX <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                        chr = "X", start = 46, end = 53, keep_all_snps = TRUE)
  plot_snpasso(out_snps_Steatosis_chrX$lod, out_snps_Steatosis_chrX$snpinfo, main = "Steatosis SNPs")
  
  Steatosis_Genes_MGI_chrX <- query_genes_mgi(chr = chr, start = 46, end = 53)
  plot(out_snps_Steatosis_chrX$lod, out_snps_Steatosis_chrX$snpinfo, drop_hilit=1.5, genes = Steatosis_Genes_MGI_chrX, main = "Steatosis Genes MGI")
  
  plot_genes(Steatosis_Genes_MGI_chrX, main = "Steatosis Genes MGI")
  
dev.off()
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Steatosis GWAS - RankZ sexgen.pdf")
out_gwas_Steatosis <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zSteatosis"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_Steatosis$lod, out_gwas_Steatosis$snpinfo, altcol="green4", gap=0, main = "Steatosis GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_Steatosis_sex <- est_herit(pheno["zSteatosis"], kinship_lmm, sex, cores = 2)
herit_Steatosis_sexgen <- est_herit(pheno["zSteatosis"], kinship_lmm, sexgen, cores = 2)

