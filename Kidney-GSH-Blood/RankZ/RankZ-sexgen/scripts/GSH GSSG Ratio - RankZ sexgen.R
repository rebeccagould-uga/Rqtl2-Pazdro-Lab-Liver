# R01 GSH DO Mapping Code 
# Updated December 2020
# Becca Gould 

#KIDNEY GLUTATHIONE + BLOOD (BUN) MAPPING - GSH/GSSG Ratio

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

#for Kidney GSH/GSSG Ratio
qtlscan_KidneyGSH_GSSGRatio<- scan1(genoprobs = probs, pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_KidneyGSH_GSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "GSH GSSG Ratio QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyGSH_GSSGRatio = summary(perm_KidneyGSH_GSSGRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_KidneyGSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_KidneyGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_KidneyGSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_KidneyGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_KidneyGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSH_GSSGRatio)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_KidneyGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_KidneyGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSH_GSSGRatio)
  
  
  write_xlsx(list("GSH_GSSGRatio gmap (cM)" = gmap_peaksGSH_GSSGRatio,
                  "GSH_GSSGRatio pmap (Mbp)" = peaksGSH_GSSGRatio),
             "GSH GSSG Ratio Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()

  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSH-GSSG-Ratio-GWAS-RankZ-sexgen.pdf")
out_gwas_KidneyGSH_GSSGRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_KidneyGSH_GSSGRatio$lod, out_gwas_KidneyGSH_GSSGRatio$snpinfo, altcol="green4", gap=0, main = "Kidney GSH/GSSG Ratio GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_KidneyGSH_GSSGRatio_sex <- est_herit(pheno["zKidneyGSH_GSSGRatio"], kinship_lmm, sex, cores = 2)
herit_KidneyGSH_GSSGRatio_sexgen <- est_herit(pheno["zKidneyGSH_GSSGRatio"], kinship_lmm, sexgen, cores = 2)


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSH GSSG Ratio Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_KidneyGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gpx1 Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_KidneyGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gclc Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_KidneyGSH_GSSGRatio_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gclm Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_KidneyGSH_GSSGRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gss Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_KidneyGSH_GSSGRatio_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gsr Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################


