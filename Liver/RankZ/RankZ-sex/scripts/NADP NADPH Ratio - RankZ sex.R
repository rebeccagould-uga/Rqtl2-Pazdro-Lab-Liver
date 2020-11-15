# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER QTL MAPPING - NADP/NADPH Ratio

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

qtlscan_LiverNADP_NADPHRatio <- scan1(genoprobs = probs, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, cores=10)
perm_LiverNADP_NADPHRatio <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADP_NADPHRatio"], addcovar = sex, n_perm = 1000, cores=10)

#set working directory
pdf(file = "NADP NADPH Ratio QTL Results - RankZ sex.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_NADPHRatio = summary(perm_LiverNADP_NADPHRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP/NADPH Ratio", ylim = c(0,11))
  abline(h = threshold_LiverNADP_NADPHRatio, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksNADP_NADPHRatio <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksNADP_NADPHRatio)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksNADP_NADPHRatio <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksNADP_NADPHRatio)
  
  
  write_xlsx(list("NADP NADPH Ratio gmap (cM)" = gmap_peaksNADP_NADPHRatio,
                  "NADP NADPH Ratio pmap (Mbp)" = peaksNADP_NADPHRatio),
             "NADP NADPH Ratio Peaks - RankZ sex.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver NADP/NADPH --- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_LiverNADP_NADPHRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,60)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

  chr = 2
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr ==  chr,"ci_lo"]
  #end = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr == chr, "ci_hi"] 
  
  pander(peaksNADP_NADPHRatio)
  #based on peaksNADP_NADPHRatio, peak of interest is ~109 Mbp
  variants_LiverNADP_NADPHRatio_chr2 <- query_variants(chr, 95, 105)
  out_snps_LiverNADP_NADPHRatio_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                                    chr = chr, start = 95, end = 105, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPHRatio_chr2$lod, out_snps_LiverNADP_NADPHRatio_chr2$snpinfo, main = "Liver NADP/NADPH SNPs")
  
  LiverNADP_NADPHRatio_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = 95, end = 105)
  plot(out_snps_LiverNADP_NADPHRatio_chr2$lod, out_snps_LiverNADP_NADPHRatio_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPHRatio_Genes_MGI_chr2, main = "Liver NADP/NADPH Genes MGI")
  
#For Liver NADP/NADPH --- Chromosome 12
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 12
  coef_blup_LiverNADP_NADPHRatio_chr12 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(20,40)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 12
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr ==  chr,"ci_lo"]
  #end = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr == chr, "ci_hi"] 

  pander(peaksNADP_NADPHRatio)
  #based on peaksNADP_NADPHRatio, peak of interest is ~28 Mbp
  variants_LiverNADP_NADPHRatio_chr12a <- query_variants(chr, 27, 30.5)
  out_snps_LiverNADP_NADPHRatio_chr12a <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                       chr = chr, start = 27, end = 30.5, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPHRatio_chr12a$lod, out_snps_LiverNADP_NADPHRatio_chr12a$snpinfo, main = "Liver NADP/NADPH SNPs")
  
  LiverNADP_NADPHRatio_Genes_MGI_chr12a <- query_genes_mgi(chr = chr, start = 27, end = 30.5)
  plot(out_snps_LiverNADP_NADPHRatio_chr12a$lod, out_snps_LiverNADP_NADPHRatio_chr12a$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPHRatio_Genes_MGI_chr12a, main = "Liver NADP/NADPH Genes MGI")
  
  #using pmap (Mbp)
  chr = 12
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 2 peak
  #start = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr ==  chr,"ci_lo"]
  #end = peaksNADP_NADPHRatio[peaksNADP_NADPHRatio$chr == chr, "ci_hi"] 
  
  pander(peaksNADP_NADPHRatio)
  #based on peaksNADP_NADPHRatio, peak of interest is ~69 Mbp
  variants_LiverNADP_NADPHRatio_chr12b <- query_variants(chr, 60, 75)
  out_snps_LiverNADP_NADPHRatio_chr12b <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, query_func = query_variants,
                                                    chr = chr, start = 60, end = 75, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPHRatio_chr12b$lod, out_snps_LiverNADP_NADPHRatio_chr12b$snpinfo, main = "Liver NADP/NADPH SNPs")
  
  LiverNADP_NADPHRatio_Genes_MGI_chr12b <- query_genes_mgi(chr = chr, start = 60, end = 75)
  plot(out_snps_LiverNADP_NADPHRatio_chr12b$lod, out_snps_LiverNADP_NADPHRatio_chr12b$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPHRatio_Genes_MGI_chr12b, main = "Liver NADP/NADPH Genes MGI")
  
#For Liver NADP/NADPH --- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 14
  coef_blup_LiverNADP_NADPHRatio_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(35,50)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

  
dev.off()  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "NADP NADPH Ratio GWAS - RankZ sex.pdf")
out_gwas_LiverNADP_NADPHRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverNADP_NADPHRatio$lod, out_gwas_LiverNADP_NADPHRatio$snpinfo, altcol="green4", gap=0, main = "Liver NADP/NADPH Ratio GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverNADP_NADPHRatio_sex <- est_herit(pheno["zLiverNADP_NADPHRatio"], kinship_lmm, sex, cores = 10)  
herit_LiverNADP_NADPHRatio_sexgen <- est_herit(pheno["zLiverNADP_NADPHRatio"], kinship_lmm, sexgen, cores = 10) 

