# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - NADPH

#Load in Liver QTL Mapping - RankZ 1000 perm - SexGen.Rdata
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

qtlscan_LiverNADPH <- scan1(genoprobs = probs, pheno = pheno["zLiverNADPH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverNADPH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADPH"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "NADPH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADPH = summary(perm_LiverNADPH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADPH", ylim = c(0,11))
  abline(h = threshold_LiverNADPH, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksNADPH <- find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSH)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADPH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksNADPH <- find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADPH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksNADPH)
  
  
  write_xlsx(list("NADPH gmap (cM)" = gmap_peaksNADPH,
                  "NADPH pmap (Mbp)" = peaksNADPH),
             "NADPH Peaks - RankZ sex.xlsx")

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADPH --- Chromosome 12
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 12
  coef_blup_LiverNADPH_chr12 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADPH_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(5,15)
  plot_coefCC(x = coef_blup_LiverNADPH_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 12
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksGSH[peaksGSH$chr ==  chr,"ci_lo"]
  #end = peaksGSH[peaksGSH$chr == chr, "ci_hi"] 
  
  pander(peaksNADPH)
  #based on peaksGSH, peak of interest is ~28 Mbp
  variants_LiverNADPH_chr12 <- query_variants(chr, 27, 30.5)
  out_snps_LiverNADPH_chr12 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                        chr = chr, start = 27, end = 30.5, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADPH_chr12$lod, out_snps_LiverNADPH_chr12$snpinfo, main = "Liver NADPH SNPs")
  
  LiverNADPH_Genes_MGI_chr12 <- query_genes_mgi(chr = chr, start = 27, end = 30.5)
  plot(out_snps_LiverNADPH_chr12$lod, out_snps_LiverNADPH_chr12$snpinfo, drop_hilit=1.5, genes = LiverNADPH_Genes_MGI_chr12, main = "Liver NADPH Genes MGI")

dev.off()  
   

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "NADPH GWAS - RankZ sexgen.pdf")
out_gwas_LiverNADPH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADPH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADPH$lod, out_gwas_LiverNADPH$snpinfo, altcol="green4", gap=0, main = "Liver NADPH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADPH_sexgen <- est_herit(pheno["zLiverNADPH"], kinship_lmm, sexgen, cores = 10)  
herit_LiverNADPH_sex <- est_herit(pheno["zLiverNADPH"], kinship_lmm, sex, cores = 10)  



