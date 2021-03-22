# R01 GSH DO Mapping Code 
# Updated January 2021
# Becca Gould 

#KIDNEY GLUTATHIONE + BLOOD (BUN) MAPPING - GSH

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

qtlscan_KidneyGSH <- scan1(genoprobs = probs, pheno = pheno["zKidneyGSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_KidneyGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

Xcovar = get_x_covar(R01_GSH_DO_QTLdata)
perm_strata = mat2strata(Xcovar)
perm_X_KidneyGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH"], addcovar = sexgen,n_perm = 1000, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(R01_GSH_DO_QTLdata$gmap), cores=10)

summary(perm_X_KidneyGSH, alpha=c(0.2, 0.1, 0.05))
#Autosome LOD thresholds (1000 permutations)
#zKidneyGSH
#0.2        6.97
#0.1        7.33
#0.05       7.80

#X chromosome LOD thresholds (18090 permutations)
#zKidneyGSH
#0.2        6.50
#0.1        6.89
#0.05       7.40

#set working directory
pdf(file = "GSH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyGSH = summary(perm_KidneyGSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney GSH", ylim = c(0,11))
  abline(h = threshold_KidneyGSH, col = c("purple", "red", "blue"), lwd = 2)

  plot_scan1(x = qtlscan_KidneyGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney GSH (X Chrom)", ylim = c(0,11))
  #perm_X_KidneyGSH_only <- perm_X_KidneyGSH[["X"]]
  #threshold_X_KidneyGSH = summary(perm_X_KidneyGSH_only, alpha = c(0.2, 0.1, 0.05))
  #abline(h = threshold_X_KidneyGSH, col = c("purple", "red", "blue"), lwd = 2)
  abline(h = c(6.50, 6.89, 7.40), col = c("purple", "red", "blue"), lwd = 2)
  
#using gmap (cM)
  find_peaks(scan1_output = qtlscan_KidneyGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSH <- find_peaks(scan1_output = qtlscan_KidneyGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSH)
  
#using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_KidneyGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSH <- find_peaks(scan1_output = qtlscan_KidneyGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSH)
  

write_xlsx(list("GSH gmap (cM)" = gmap_peaksGSH,
                "GSH pmap (Mbp)" = peaksGSH),
                "GSH Peaks - RankZ sexgen.xlsx")

####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Kidney GSH --- Chromosome X
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
#using gmap (cM)
  chr = "X"
  coef_blup_KidneyGSH_chrX <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chrX, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,40)
  plot_coefCC(x = coef_blup_KidneyGSH_chrX, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
#using pmap (Mbp)
  chr = "X"
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome X peak
  #start = peaksGSH[peaksGSH$chr ==  chr,"ci_lo"]
  #end = peaksGSH[peaksGSH$chr == chr, "ci_hi"] 
  
  pander(peaksGSH)
  #based on peaksGSH, peak of interest is ~109 Mbp
  variants_KidneyGSH_chrX <- query_variants(chr, 48.2, 52.9)
  out_snps_KidneyGSH_chrX <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                            chr = chr, start = 48.2, end = 52.9, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyGSH_chrX$lod, out_snps_KidneyGSH_chrX$snpinfo, main = "Kidney GSH SNPs")
  
  KidneyGSH_Genes_MGI_chrX <- query_genes_mgi(chr = chr, start = 48.2, end = 52.9)
  plot(out_snps_KidneyGSH_chrX$lod, out_snps_KidneyGSH_chrX$snpinfo, drop_hilit=1.5, genes = KidneyGSH_Genes_MGI_chrX, main = "Kidney GSH Genes MGI")
  
#For Kidney GSH --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "11"
  coef_blup_KidneyGSH_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,75)
  plot_coefCC(x = coef_blup_KidneyGSH_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "11"
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome X peak
  #start = peaksGSH[peaksGSH$chr ==  chr,"ci_lo"]
  #end = peaksGSH[peaksGSH$chr == chr, "ci_hi"] 
  
  pander(peaksGSH)
  #based on peaksGSH, peak of interest is ~100 Mbp
  variants_KidneyGSH_chr11 <- query_variants(chr, 99, 102)
  out_snps_KidneyGSH_chr11 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = 99, end = 102, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyGSH_chr11$lod, out_snps_KidneyGSH_chr11$snpinfo, main = "Kidney GSH SNPs")
  
  KidneyGSH_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 99, end = 102)
  plot(out_snps_KidneyGSH_chr11$lod, out_snps_KidneyGSH_chr11$snpinfo, drop_hilit=1.5, genes = KidneyGSH_Genes_MGI_chr11, main = "Kidney GSH Genes MGI")
  
  
dev.off()

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSH-GWAS-RankZ-sexgen.pdf")
out_gwas_KidneyGSH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_KidneyGSH$lod, out_gwas_KidneyGSH$snpinfo, altcol="green4", gap=0, main = "Kidney GSH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_KidneyGSH_sex <- est_herit(pheno["zKidneyGSH"], kinship_lmm, sex, cores = 2)
herit_KidneyGSH_sexgen <- est_herit(pheno["zKidneyGSH"], kinship_lmm, sexgen, cores = 2)


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSH Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_KidneyGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gpx1 Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_KidneyGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gclc Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_KidneyGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_KidneyGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gclm Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_KidneyGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_KidneyGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gss Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_KidneyGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_KidneyGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gsr Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################

