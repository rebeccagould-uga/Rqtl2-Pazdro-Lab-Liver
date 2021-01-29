# R01 GSH DO Mapping Code 
# Updated November 2020
# Becca Gould 

#KIDNEY GLUTATHIONE + BLOOD (BUN) MAPPING - Total GSH

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

qtlscan_KidneyTotalGSH<- scan1(genoprobs = probs, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_KidneyTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyTotalGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

Xcovar = get_x_covar(R01_GSH_DO_QTLdata)
perm_strata = mat2strata(Xcovar)
perm_X_KidneyTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyTotalGSH"], addcovar = sexgen, Xcovar = Xcovar, n_perm = 1, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(R01_GSH_DO_QTLdata$gmap), cores=2)

#set working directory
pdf(file = "Total GSH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_KidneyTotalGSH = summary(perm_KidneyTotalGSH, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_KidneyTotalGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Kidney Total GSH", ylim = c(0,11))
abline(h = threshold_KidneyTotalGSH, col = c("purple", "red", "blue"), lwd = 2)

#using gmap (cM)
find_peaks(scan1_output = qtlscan_KidneyTotalGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
gmap_peaksTotalGSH <- find_peaks(scan1_output = qtlscan_KidneyTotalGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_KidneyTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(gmap_peaksTotalGSH)

#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_KidneyTotalGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyTotalGSH, alpha = 0.2), peakdrop = 1.0, prob = 0.95, expand2markers = FALSE)
peaksTotalGSH <- find_peaks(scan1_output = qtlscan_KidneyTotalGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_KidneyTotalGSH, alpha = 0.2), peakdrop = 1.0, prob = 0.95, expand2markers = FALSE)
print(peaksTotalGSH)


write_xlsx(list("Total GSH gmap (cM)" = gmap_peaksTotalGSH,
                "Total GSH pmap (Mbp)" = peaksTotalGSH),
           "Total GSH Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Kidney Total GSH --- Chromosome X
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = X
  coef_blup_KidneyTotalGSH_chrX <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chrX, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,65)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chrX, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = X
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome X peak
  #start = peaksTotalGSH[peaksTotalGSH$chr ==  chr,"ci_lo"]
  #end = peaksTotalGSH[peaksTotalGSH$chr == chr, "ci_hi"] 
  
  pander(peaksTotalGSH)
  #based on peaksTotalGSH, peak of interest is ~109 Mbp
  variants_KidneyTotalGSH_chrX <- query_variants(chr, 107, 111)
  out_snps_KidneyTotalGSH_chrX <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                           chr = chr, start = 107, end = 111, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyTotalGSH_chrX$lod, out_snps_KidneyTotalGSH_chrX$snpinfo, main = "Kidney Total GSH SNPs")
  
  KidneyTotalGSH_Genes_MGI_chrX <- query_genes_mgi(chr = chr, start = 107, end = 111)
  plot(out_snps_KidneyTotalGSH_chrX$lod, out_snps_KidneyTotalGSH_chrX$snpinfo, drop_hilit=1.5, genes = KidneyTotalGSH_Genes_MGI_chrX, main = "Kidney Total GSH Genes MGI")
  
dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Total-GSH-GWAS-RankZ-sexgen.pdf")
out_gwas_KidneyTotalGSH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_KidneyTotalGSH$lod, out_gwas_KidneyTotalGSH$snpinfo, altcol="green4", gap=0, main = "Kidney Total GSH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_KidneyTotalGSH_sex <- est_herit(pheno["zKidneyTotalGSH"], kinship_lmm, sex, cores = 2)
herit_KidneyTotalGSH_sexgen <- est_herit(pheno["zKidneyTotalGSH"], kinship_lmm, sexgen, cores = 2)


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "Total GSH Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_KidneyTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gpx1 Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_KidneyTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gclc Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_KidneyTotalGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gclm Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_KidneyTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gss Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_KidneyTotalGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_KidneyTotalGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gsr Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################



