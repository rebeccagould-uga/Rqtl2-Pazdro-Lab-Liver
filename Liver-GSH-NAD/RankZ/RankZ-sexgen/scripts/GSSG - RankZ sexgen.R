# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - GSSG

#Load in Liver-GSH-NAD-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**

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

qtlscan_LiverGSSG <- scan1(genoprobs = probs, pheno = pheno["zLiverGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverGSSG <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "GSSG QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSSG = summary(perm_LiverGSSG, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSSG", ylim = c(0,11))
  abline(h = threshold_LiverGSSG, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSSG", ylim = c(0,11))
  abline(h = threshold_LiverGSSG, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSSG)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSSG)
  
  
  write_xlsx(list("GSSG gmap (cM)" = gmap_peaksGSSG,
                  "GSSG pmap (Mbp)" = peaksGSSG),
             "GSSG Peaks - RankZ sexgen.xlsx")

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver GSSG --- Chromosome 1
  #don't need to include in final report, chromosome 1 peak is p < 0.50
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 1
  coef_blup_LiverGSSG_chr1 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  
#For Liver GSSG --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_LiverGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  
#For Liver GSSG --- Chromosome 18
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 18
  coef_blup_LiverGSSG_chr18 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSSG_chr18, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,55)
  plot_coefCC(x = coef_blup_LiverGSSG_chr18, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
dev.off()
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSSG GWAS - RankZ sexgen.pdf")
out_gwas_LiverGSSG <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverGSSG$lod, out_gwas_LiverGSSG$snpinfo, altcol="green4", gap=0, main = "Liver GSSG GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverGSSG_sexgen <- est_herit(pheno["zLiverGSSG"], kinship_lmm, sexgen, cores = 10)
herit_LiverGSSG_sex <- est_herit(pheno["zLiverGSSG"], kinship_lmm, sex, cores = 10)



##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSSG Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_LiverGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gpx1 Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_LiverGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gclc Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_LiverGSSG_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSSG_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_LiverGSSG_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gclm Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_LiverGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gss Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_LiverGSSG_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSSG_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_LiverGSSG_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gsr Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################

#Looking at significant peak for Liver GSH/GSSG --- Chromosome 16
par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
chr = 16
coef_blup_LiverGSSG_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSSG_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(1,9)
plot_coefCC(x = coef_blup_LiverGSSG_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)









