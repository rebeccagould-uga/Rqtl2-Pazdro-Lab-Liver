# R01 Ballooning DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Ballooning

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_Ballooning_DO_mapping_Liver/data")

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

qtlscan_Ballooning <- scan1(genoprobs = probs, pheno = pheno["zBallooning"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_Ballooning <- scan1perm(genoprobs = probs, pheno = pheno["zBallooning"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Ballooning QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Ballooning = summary(perm_Ballooning, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Ballooning", ylim = c(0,11))
  abline(h = threshold_Ballooning, col = c("purple", "red", "blue"), lwd = 2)
  plot_scan1(x = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Ballooning", ylim = c(0,11))
  abline(h = threshold_Ballooning, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
#using gmap (cM)
  find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksBallooning <- find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksBallooning)
  
#using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksBallooning <- find_peaks(scan1_output = qtlscan_Ballooning, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksBallooning)
  

write_xlsx(list("Ballooning gmap (cM)" = gmap_peaksBallooning,
                "Ballooning pmap (Mbp)" = peaksBallooning),
                "Ballooning Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Ballooning --- Chromosome 6
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 6
  coef_blup_Ballooning_chr6 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_Ballooning_chr6, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(5,30)
  plot_coefCC(x = coef_blup_Ballooning_chr6, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  xlim <- c(60,80)
  plot_coefCC(x = coef_blup_Ballooning_chr6, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 6
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksBallooning[peaksBallooning$chr ==  chr,"ci_lo"]
  #end = peaksBallooning[peaksBallooning$chr == chr, "ci_hi"] 
  
  pander(peaksBallooning)
  #based on peaksBallooning, peak of interest is ~109 Mbp
  variants_Ballooning_chr6 <- query_variants(chr, 35, 37)
  out_snps_Ballooning_chr6 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                      chr = chr, start = 35, end = 37, keep_all_snps = TRUE)
  plot_snpasso(out_snps_Ballooning_chr6$lod, out_snps_Ballooning_chr6$snpinfo, main = "Ballooning SNPs")
  
  Ballooning_Genes_MGI_chr6 <- query_genes_mgi(chr = chr, start = 35, end = 37)
  plot(out_snps_Ballooning_chr6$lod, out_snps_Ballooning_chr6$snpinfo, drop_hilit=1.5, genes = Ballooning_Genes_MGI_chr6, main = "Ballooning Genes MGI")
  
dev.off()
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Ballooning GWAS - RankZ sexgen.pdf")
out_gwas_Ballooning <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zBallooning"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_Ballooning$lod, out_gwas_Ballooning$snpinfo, altcol="green4", gap=0, main = "Ballooning GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_Ballooning_sex <- est_herit(pheno["zBallooning"], kinship_lmm, sex, cores = 2)
herit_Ballooning_sexgen <- est_herit(pheno["zBallooning"], kinship_lmm, sexgen, cores = 2)


##################################################################
## Checking other NAFLD genes 
##################################################################

#set working directory
pdf(file = "Ballooning-NAFLD-Genes-RankZ-sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Pnpla3 - chr 15 39.78 cM
chr = 15
#coef_blup_Ballooning_chr15 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr15, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(38.5,41.5)
plot_coefCC(x = coef_blup_Ballooning_chr15, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Pnpla3 Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Tm6sf2 - chr 8 34.15 cM
chr = 8
#coef_blup_Ballooning_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(32.5,35.5)
plot_coefCC(x = coef_blup_Ballooning_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Tm6sf2 and Ncan Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Mboat7 - chr 7 2.12 cM
chr = 7
#coef_blup_Ballooning_chr7 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr7, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(1.8,3.5)
plot_coefCC(x = coef_blup_Ballooning_chr7, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Mboat7 Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Gckr - chr 5 17.27 cM
chr = 5
#coef_blup_Ballooning_chr5 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(15,18.5)
plot_coefCC(x = coef_blup_Ballooning_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Gckr Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Hsd17b13 - chr 5 50.46 cM
chr = 5
#coef_blup_Ballooning_chr5 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(48.5, 51.5)
plot_coefCC(x = coef_blup_Ballooning_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Hsd17b13 Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Ncan - chr 8 34.15 cM - SAME AS TM6SF2 GENE!
#chr = 8
#coef_blup_Ballooning_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_Ballooning_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
#xlim <- c(48.5, 51.5)
#plot_coefCC(x = coef_blup_Ballooning_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Pnpla3 Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Lyplal1 - chr 1 89.81 cM
chr = 1
#coef_blup_Ballooning_chr1 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(48.5, 51.5)
plot_coefCC(x = coef_blup_Ballooning_chr1, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Lyplal1 Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Ppp1r3b - chr 8 21.16 cM 
chr = 8
#coef_blup_Ballooning_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zBallooning"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_Ballooning_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19.5, 23.5)
plot_coefCC(x = coef_blup_Ballooning_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_Ballooning, main = "Ppp1r3b Position -- Ballooning BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)



dev.off()

