# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - GSH

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

qtlscan_LiverGSH <- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverGSH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "GSH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH = summary(perm_LiverGSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH", ylim = c(0,11))
  abline(h = threshold_LiverGSH, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH", ylim = c(0,11))
  abline(h = threshold_LiverGSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #segments(x0, y0, x1 = x0, y1 = y0)
  #x0, y0 coordinates of points from which to draw.
  #x1, y1 coordinates of points to which to draw. At least one must be supplied.
  
  #plotting separate autosome versus X axis significance thresholds
  #plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH", ylim = c(0,11))
  #segments(x0 = 0, y0 = threshold_LiverGSH, x1 = 1690, y1 = threshold_LiverGSH, col = c("purple", "red", "blue"), "dashed")
      #segments(x0 = 1690, y0 = 6, x1 = 2000, y1 = 6, col = "purple", "dashed")
     #or
      #segments(x0 = 1690, y0 = threshold_X_LiverGSH, x1 = 2000, y1 = threshold_X_LiverGSH, col = c("purple", "red", "blue"))
 
  ############################## 
  
  #comparing plots - Liver GSH with  Kidney Eh - cM
  pdf(file = "KidneyEh-LiverGlutathione-Overlay-cM.pdf")
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap, ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSH", "Renal Eh"), bg="gray90")
  
    plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap, chr = "14", ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", chr = "14", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSH", "Renal Eh"), bg="gray90")
  
    plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap,  ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Renal Eh"), bg="gray90")
    
    plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap, chr = "14", ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", chr = "14", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Renal Eh"), bg="gray90")
    
    dev.off()
    
  #comparing plots - Liver GSH with  Kidney Eh - Mbp
  pdf(file = "KidneyEh-LiverGlutathione-Overlay-Mbp.pdf")
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$pmap, ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSH", "Renal Eh"), bg="gray90")
    
    plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$pmap, chr = "14", ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", chr = "14", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSH", "Renal Eh"), bg="gray90")
    
    plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap,  ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Renal Eh"), bg="gray90")
    
    plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap, chr = "14", ylim = c(0,11))
    plot_scan1(x = qtlscan_KidneyEh, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", chr = "14", add = TRUE)
    legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Renal Eh"), bg="gray90")
    
    dev.off()
    
  ############################## 
    
  #comparing plots - Liver Glutathione with Steatosis - cM
  pdf(file = "Steatosis-LiverGlutathione-Overlay-cM.pdf")
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap, ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSSG", "Steatosis"), bg="gray90")
    
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap, chr = "18",  ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", chr = "18", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSSG", "Steatosis"), bg="gray90")
  
  plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap, ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Steatosis"), bg="gray90")
  
  plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap, chr = "18", ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$gmap, col = "#009999", chr = "18", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Steatosis"), bg="gray90")
  
  dev.off()
  
  #comparing plots - Liver Glutathione with Steatosis - Mbp
  pdf(file = "Steatosis-LiverGlutathione-Overlay-Mbp.pdf")
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSSG", "Steatosis"), bg="gray90")
  
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, chr = "18",  ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", chr = "18", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic GSSG", "Steatosis"), bg="gray90")
  
  plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap, ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Steatosis"), bg="gray90")
  
  plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap, chr = "18", ylim = c(0,11))
  plot_scan1(x = qtlscan_Steatosis, map = R01_GSH_DO_QTLdata$pmap, col = "#009999", chr = "18", add = TRUE)
  legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Hepatic Total Glutathione", "Steatosis"), bg="gray90")
  
  dev.off()
  
  ############################## 
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSH <- find_peaks(scan1_output = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSH)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSH <- find_peaks(scan1_output = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSH)
  
 
  
  write_xlsx(list("GSH gmap (cM)" = gmap_peaksGSH,
                  "GSH pmap (Mbp)" = peaksGSH),
             "GSH Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver GSH --- Chromosome 14 
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 14
  coef_blup_LiverGSH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,20)
  plot_coefCC(x = coef_blup_LiverGSH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

  #using pmap (Mbp)
  chr = 14
  #could use ci_lo or ci_hi, but in this case, I want a specific chromosome 14 peak
  #start = peaksGSH[peaksGSH$chr ==  chr,"ci_lo"]
  #end = peaksGSH[peaksGSH$chr == chr, "ci_hi"] 
  
  pander(peaksGSH)
  #based on peaksGSH, peak of interest is ~22 Mbp
  variants_LiverGSH_chr14 <- query_variants(chr, 21, 25)
  out_snps_LiverGSH_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                         chr = chr, start = 21, end = 25, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_chr14$lod, out_snps_LiverGSH_chr14$snpinfo, main = "Liver GSH SNPs")
    
  LiverGSH_Genes_MGI_chr14 <- query_genes_mgi(chr = chr, start = 21, end = 25)
  plot(out_snps_LiverGSH_chr14$lod, out_snps_LiverGSH_chr14$snpinfo, drop_hilit=1.5, genes = LiverGSH_Genes_MGI_chr14, main = "Liver GSH Genes MGI")

  plot_genes (LiverGSH_Genes_MGI_chr14, main = "Liver GSH Genes MGI")
  
dev.off()

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSH GWAS - RankZ sexgen.pdf")
out_gwas_LiverGSH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverGSH$lod, out_gwas_LiverGSH$snpinfo, altcol="green4", gap=0, main = "Liver GSH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverGSH_sexgen <- est_herit(pheno["zLiverGSH"], kinship_lmm, sexgen, cores = 10)
herit_LiverGSH_sex <- est_herit(pheno["zLiverGSH"], kinship_lmm, sex, cores = 10)



##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSH Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_LiverGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gpx1 Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_LiverGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gclc Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_LiverGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_LiverGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gclm Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_LiverGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_LiverGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gss Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_LiverGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_LiverGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gsr Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################

#Looking at significant peak for Liver GSH/GSSG --- Chromosome 16
par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
chr = 16
coef_blup_LiverGSH_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(1,9)
plot_coefCC(x = coef_blup_LiverGSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)





