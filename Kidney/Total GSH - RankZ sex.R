# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - Total GSH

#Load in Liver QTL Mapping - RankZ 1000 perm - sex.Rdata
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

qtlscan_LiverTotalGSH<- scan1(genoprobs = probs, pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_LiverTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverTotalGSH"], addcovar = sex, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Total GSH QTL Results - RankZ sex.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverTotalGSH = summary(perm_LiverTotalGSH, alpha = c(0.2, 0.1, 0.05))
plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Total GSH", ylim = c(0,11))
abline(h = threshold_LiverTotalGSH, col = c("purple", "red", "blue"), lwd = 2)

#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
gmap_peaksTotalGSH <- find_peaks(scan1_output = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(gmap_peaksTotalGSH)

#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverTotalGSH, alpha = 0.2), peakdrop = 1.0, prob = 0.95, expand2markers = FALSE)
peaksTotalGSH <- find_peaks(scan1_output = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverTotalGSH, alpha = 0.2), peakdrop = 1.0, prob = 0.95, expand2markers = FALSE)
print(peaksTotalGSH)


write_xlsx(list("Total GSH gmap (cM)" = gmap_peaksTotalGSH,
                "Total GSH pmap (Mbp)" = peaksTotalGSH),
           "Total GSH Peaks - RankZ sex.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver Total GSH --- Chromosome 2 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_LiverTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,65)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksTotalGSH[peaksTotalGSH$chr ==  chr,"ci_lo"]
  #end = peaksTotalGSH[peaksTotalGSH$chr == chr, "ci_hi"] 
  
  pander(peaksTotalGSH)
  #based on peaksTotalGSH, peak of interest is ~109 Mbp
  variants_LiverTotalGSH_chr2 <- query_variants(chr, 107, 111)
  out_snps_LiverTotalGSH_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                           chr = chr, start = 107, end = 111, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverTotalGSH_chr2$lod, out_snps_LiverTotalGSH_chr2$snpinfo, main = "Liver Total GSH SNPs")
  
  LiverTotalGSH_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = 107, end = 111)
  plot(out_snps_LiverTotalGSH_chr2$lod, out_snps_LiverTotalGSH_chr2$snpinfo, drop_hilit=1.5, genes = LiverTotalGSH_Genes_MGI_chr2, main = "Liver Total GSH Genes MGI")
  
dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Total GSH GWAS - RankZ sex.pdf")
out_gwas_LiverTotalGSH <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverTotalGSH$lod, out_gwas_LiverTotalGSH$snpinfo, altcol="green4", gap=0, main = "Liver Total GSH GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverTotalGSH_sex <- est_herit(pheno["zLiverTotalGSH"], kinship_lmm, sex, cores = 10)
herit_LiverTotalGSH_sexgen <- est_herit(pheno["zLiverTotalGSH"], kinship_lmm, sexgen, cores = 10)


