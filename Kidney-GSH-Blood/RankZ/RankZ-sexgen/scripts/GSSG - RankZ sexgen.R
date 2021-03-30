# R01 GSH DO Mapping Code 
# Updated December 2020
# Becca Gould 

#KIDNEY GLUTATHIONE + BLOOD (BUN) MAPPING - GSSG

#Load in Kidney-QTL-Mapping-RankZ-sexgen.Rdata
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

qtlscan_KidneyGSSG <- scan1(genoprobs = probs, pheno = pheno["zKidneyGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_KidneyGSSG <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "GSSG QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyGSSG = summary(perm_KidneyGSSG, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyGSSG, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney GSSG", ylim = c(0,11))
  abline(h = threshold_KidneyGSSG, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_KidneyGSSG, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney GSSG", ylim = c(0,11))
  abline(h = threshold_KidneyGSSG, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_KidneyGSSG, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSSG <- find_peaks(scan1_output = qtlscan_KidneyGSSG, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSSG)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_KidneyGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSSG <- find_peaks(scan1_output = qtlscan_KidneyGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSSG)
  
  
  write_xlsx(list("GSSG gmap (cM)" = gmap_peaksGSSG,
                  "GSSG pmap (Mbp)" = peaksGSSG),
             "GSSG Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()  
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSSG-GWAS-RankZ-sexgen.pdf")
out_gwas_KidneyGSSG <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyGSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_KidneyGSSG$lod, out_gwas_KidneyGSSG$snpinfo, altcol="green4", gap=0, main = "Kidney GSSG GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_KidneyGSSG_sex <- est_herit(pheno["zKidneyGSSG"], kinship_lmm, sex, cores = 2)
herit_KidneyGSSG_sexgen <- est_herit(pheno["zKidneyGSSG"], kinship_lmm, sexgen, cores = 2)


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSSG Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_KidneyGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gpx1 Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_KidneyGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gclc Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_KidneyGSSG_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSSG_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_KidneyGSSG_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gclm Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_KidneyGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_KidneyGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gss Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_KidneyGSSG_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSSG_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_KidneyGSSG_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gsr Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################


