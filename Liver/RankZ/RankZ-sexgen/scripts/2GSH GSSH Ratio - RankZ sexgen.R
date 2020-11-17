# R01 GSH DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER QTL MAPPING - 2GSH/GSSG Ratio

#Load in Liver QTL Mapping - RankZ 1000 perm - sexgen.Rdata
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

#for Liver 2GSH/GSSG Ratio
qtlscan_Liver2GSH_GSSGRatio<- scan1(genoprobs = probs, pheno = pheno["zLiver2GSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_Liver2GSH_GSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["zLiver2GSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "2GSH GSSG Ratio QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Liver2GSH_GSSGRatio = summary(perm_Liver2GSH_GSSGRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Liver2GSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver 2GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_Liver2GSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_Liver2GSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Liver2GSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaks2GSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_Liver2GSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Liver2GSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaks2GSH_GSSGRatio)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_Liver2GSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Liver2GSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaks2GSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_Liver2GSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Liver2GSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaks2GSH_GSSGRatio)
  
  
  write_xlsx(list("2GSH_GSSGRatio gmap (cM)" = gmap_peaks2GSH_GSSGRatio,
                  "2GSH_GSSGRatio pmap (Mbp)" = peaks2GSH_GSSGRatio),
             "2GSH GSSG Ratio Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver 2GSH/GSSG --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 16
  coef_blup_Liver2GSH_GSSGRatio_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiver2GSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_Liver2GSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Liver2GSH_GSSGRatio, main = "Liver 2GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,12)
  plot_coefCC(x = coef_blup_Liver2GSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Liver2GSH_GSSGRatio, main = "Liver 2GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 16
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 16 peak
  #start = peaks2GSH_GSSGRatio[peaks2GSH_GSSGRatio$chr ==  chr,"ci_lo"]
  #end = peaks2GSH_GSSGRatio[peaks2GSH_GSSGRatio$chr == chr, "ci_hi"] 
  
  pander(peaks2GSH_GSSGRatio)
  #based on peaks2GSH_GSSGRatio, peak of interest is ~9 Mbp
  variants_Liver2GSH_GSSGRatio_chr16 <- query_variants(chr, 7, 11)
  out_snps_Liver2GSH_GSSGRatio_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiver2GSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                 chr = chr, start = 7, end = 11, keep_all_snps = TRUE)
  plot_snpasso(out_snps_Liver2GSH_GSSGRatio_chr16$lod, out_snps_Liver2GSH_GSSGRatio_chr16$snpinfo, main = "Liver 2GSH/GSSG SNPs")
  
  Liver2GSH_GSSGRatio_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 7, end = 11)
  plot(out_snps_Liver2GSH_GSSGRatio_chr16$lod, out_snps_Liver2GSH_GSSGRatio_chr16$snpinfo, drop_hilit=1.5, genes = Liver2GSH_GSSGRatio_Genes_MGI_chr16, main = "Liver 2GSH/GSSG Genes MGI")
  
dev.off()

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "2GSH GSSG Ratio GWAS - RankZ sexgen.pdf")
out_gwas_Liver2GSH_GSSGRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiver2GSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_Liver2GSH_GSSGRatio$lod, out_gwas_Liver2GSH_GSSGRatio$snpinfo, altcol="green4", gap=0, main = "Liver 2GSH/GSSG Ratio GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_Liver2GSH_GSSGRatio_sex <- est_herit(pheno["zLiver2GSH_GSSGRatio"], kinship_lmm, sex, cores = 10)
herit_Liver2GSH_GSSGRatio_sexgen <- est_herit(pheno["zLiver2GSH_GSSGRatio"], kinship_lmm, sexgen, cores = 10)

