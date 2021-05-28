# R01 DO Mapping Code 
# Updated March 2021  
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Inflammation

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_Inflammation_DO_mapping_Liver/data")

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

qtlscan_Inflammation <- scan1(genoprobs = probs, pheno = pheno["zInflammation"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_Inflammation <- scan1perm(genoprobs = probs, pheno = pheno["zInflammation"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Inflammation QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Inflammation = summary(perm_Inflammation, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Inflammation", ylim = c(0,11))
  abline(h = threshold_Inflammation, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Inflammation", ylim = c(0,11))
  abline(h = threshold_Inflammation, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
#using gmap (cM)
  find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksInflammation <- find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksInflammation)
  
#using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksInflammation <- find_peaks(scan1_output = qtlscan_Inflammation, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksInflammation)
  

write_xlsx(list("Inflammation gmap (cM)" = gmap_peaksInflammation,
                "Inflammation pmap (Mbp)" = peaksInflammation),
                "Inflammation Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Inflammation GWAS - RankZ sexgen.pdf")
out_gwas_Inflammation <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zInflammation"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_Inflammation$lod, out_gwas_Inflammation$snpinfo, altcol="green4", gap=0, main = "Inflammation GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_Inflammation_sex <- est_herit(pheno["zInflammation"], kinship_lmm, sex, cores = 2)
herit_Inflammation_sexgen <- est_herit(pheno["zInflammation"], kinship_lmm, sexgen, cores = 2)


##################################################################
## Checking other NAFLD genes 
##################################################################

#set working directory
pdf(file = "Inflammation-NAFLD-Genes-RankZ-sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Pnpla3 - chr 15 39.78 cM
chr = 15
#coef_blup_Inflammation_chr15 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr15, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(38.5,41.5)
plot_coefCC(x = coef_blup_Inflammation_chr15, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Pnpla3 Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Tm6sf2 - chr 8 34.15 cM
chr = 8
#coef_blup_Inflammation_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(32.5,35.5)
plot_coefCC(x = coef_blup_Inflammation_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Tm6sf2 and Ncan Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Mboat7 - chr 7 2.12 cM
chr = 7
#coef_blup_Inflammation_chr7 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr7, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(1.8,3.5)
plot_coefCC(x = coef_blup_Inflammation_chr7, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Mboat7 Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Gckr - chr 5 17.27 cM
chr = 5
#coef_blup_Inflammation_chr5 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(15,18.5)
plot_coefCC(x = coef_blup_Inflammation_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Gckr Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Hsd17b13 - chr 5 50.46 cM
chr = 5
#coef_blup_Inflammation_chr5 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(48.5, 51.5)
plot_coefCC(x = coef_blup_Inflammation_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Hsd17b13 Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Ncan - chr 8 34.15 cM - SAME AS TM6SF2 GENE!
#chr = 8
#coef_blup_Inflammation_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_Inflammation_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
#xlim <- c(48.5, 51.5)
#plot_coefCC(x = coef_blup_Inflammation_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Pnpla3 Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Lyplal1 - chr 1 89.81 cM
chr = 1
#coef_blup_Inflammation_chr1 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(48.5, 51.5)
plot_coefCC(x = coef_blup_Inflammation_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Lyplal1 Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Ppp1r3b - chr 8 21.16 cM 
chr = 8
#coef_blup_Inflammation_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zInflammation"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Inflammation_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19.5, 23.5)
plot_coefCC(x = coef_blup_Inflammation_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Inflammation, main = "Ppp1r3b Position -- Inflammation BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)



dev.off()

