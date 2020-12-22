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

#For Kidney Total GSH --- Chromosome 2 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_KidneyTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,65)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksTotalGSH[peaksTotalGSH$chr ==  chr,"ci_lo"]
  #end = peaksTotalGSH[peaksTotalGSH$chr == chr, "ci_hi"] 
  
  pander(peaksTotalGSH)
  #based on peaksTotalGSH, peak of interest is ~109 Mbp
  variants_KidneyTotalGSH_chr2 <- query_variants(chr, 107, 111)
  out_snps_KidneyTotalGSH_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                           chr = chr, start = 107, end = 111, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyTotalGSH_chr2$lod, out_snps_KidneyTotalGSH_chr2$snpinfo, main = "Kidney Total GSH SNPs")
  
  KidneyTotalGSH_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = 107, end = 111)
  plot(out_snps_KidneyTotalGSH_chr2$lod, out_snps_KidneyTotalGSH_chr2$snpinfo, drop_hilit=1.5, genes = KidneyTotalGSH_Genes_MGI_chr2, main = "Kidney Total GSH Genes MGI")
  
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


