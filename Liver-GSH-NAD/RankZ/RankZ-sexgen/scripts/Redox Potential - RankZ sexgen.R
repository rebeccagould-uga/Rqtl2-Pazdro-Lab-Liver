# R01 GSH DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - Redox Potential GSSG/2GSH

#Load in Liver-GSH-NAD-RankZ-SexGen.Rdata
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

#for Liver Redox Potential GSSG 2GSH
qtlscan_LiverRedoxPotentialGSSG2GSH<- scan1(genoprobs = probs, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverRedoxPotentialGSSG2GSH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Redox Potential GSSG 2GSH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverRedoxPotentialGSSG2GSH = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Redox Potential GSSG/2GSH", ylim = c(0,11))
  abline(h = threshold_LiverRedoxPotentialGSSG2GSH, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksLiverRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksLiverRedoxPotentialGSSG2GSH)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksLiverRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksLiverRedoxPotentialGSSG2GSH)
  
  
  write_xlsx(list("LiverRedoxPotential gmap (cM)" = gmap_peaksLiverRedoxPotentialGSSG2GSH,
                  "LiverRedoxPotential pmap (Mbp)" = peaksLiverRedoxPotentialGSSG2GSH),
             "Redox Potential GSSG 2GSH Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver Redox Potential GSSG/2GSH --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 16
  coef_blup_LiverRedoxPotentialGSSG2GSH_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,12)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 16
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 16 peak
  #start = peaksLiverRedoxPotentialGSSG2GSH[peaksLiverRedoxPotentialGSSG2GSH$chr ==  chr,"ci_lo"]
  #end = peaksLiverRedoxPotentialGSSG2GSH[peaksLiverRedoxPotentialGSSG2GSH$chr == chr, "ci_hi"] 
  
  pander(peaksLiverRedoxPotentialGSSG2GSH)
  #based on peaksLiverRedoxPotentialGSSG2GSH, peak of interest is ~9 Mbp
  variants_LiverRedoxPotentialGSSG2GSH_chr16 <- query_variants(chr, 7, 11)
  out_snps_LiverRedoxPotentialGSSG2GSH_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                 chr = chr, start = 7, end = 11, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverRedoxPotentialGSSG2GSH_chr16$lod, out_snps_LiverRedoxPotentialGSSG2GSH_chr16$snpinfo, main = "Liver Redox Potential GSSG/2GSH SNPs")
  
  LiverRedoxPotentialGSSG2GSH_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 7, end = 11)
  plot(out_snps_LiverRedoxPotentialGSSG2GSH_chr16$lod, out_snps_LiverRedoxPotentialGSSG2GSH_chr16$snpinfo, drop_hilit=1.5, genes = LiverRedoxPotentialGSSG2GSH_Genes_MGI_chr16, main = "Liver Redox Potential GSSG/2GSH Genes MGI")
  
dev.off()

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Redox Potential GSSG 2GSH GWAS - RankZ sexgen.pdf")
out_gwas_LiverRedoxPotentialGSSG2GSH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverRedoxPotentialGSSG2GSH$lod, out_gwas_LiverRedoxPotentialGSSG2GSH$snpinfo, altcol="green4", gap=0, main = "Liver Redox Potential GSSG/2GSH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverRedoxPotentialGSSG2GSH_sex <- est_herit(pheno["zLiverRedoxPotentialGSSG2GSH"], kinship_lmm, sex, cores = 10)
herit_LiverRedoxPotentialGSSG2GSH_sexgen <- est_herit(pheno["zLiverRedoxPotentialGSSG2GSH"], kinship_lmm, sexgen, cores = 10)

