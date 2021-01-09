# R01 GSH DO Mapping Code 
# Updated August 2020
# Becca Gould 

#LIVER GLUTATHIONE + NAD MAPPING - GSH GSSG Ratio

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

#for Liver GSH/GSSG Ratio
qtlscan_LiverGSH_GSSGRatio<- scan1(genoprobs = probs, pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_LiverGSH_GSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "GSH GSSG Ratio QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH_GSSGRatio = summary(perm_LiverGSH_GSSGRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_LiverGSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2)
  
  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksGSH_GSSGRatio)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksGSH_GSSGRatio)
  
  
  write_xlsx(list("GSH_GSSGRatio gmap (cM)" = gmap_peaksGSH_GSSGRatio,
                  "GSH_GSSGRatio pmap (Mbp)" = peaksGSH_GSSGRatio),
             "GSH GSSG Ratio Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver GSH/GSSG --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
  chr = 16
  coef_blup_LiverGSH_GSSGRatio_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,9)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
#using pmap (Mbp)
  chr = 16
  #could use ci_lo or ci_hi, but for this case, I want a specific chromosome 16 peak
  #start = peaksGSH_GSSGRatio[peaksGSH_GSSGRatio$chr ==  chr,"ci_lo"]
  #end = peaksGSH_GSSGRatio[peaksGSH_GSSGRatio$chr == chr, "ci_hi"] 

  pander(peaksGSH_GSSGRatio)
  #based on peaksGSH_GSSGRatio, peak of interest is ~9 Mbp
  variants_LiverGSH_GSSGRatio_chr16 <- query_variants(chr, 7, 11)
  out_snps_LiverGSH_GSSGRatio_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = 7, end = 11, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_GSSGRatio_chr16$lod, out_snps_LiverGSH_GSSGRatio_chr16$snpinfo, main = "Liver GSH/GSSG SNPs")
  
  LiverGSH_GSSGRatio_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 7, end = 11)
  plot(out_snps_LiverGSH_GSSGRatio_chr16$lod, out_snps_LiverGSH_GSSGRatio_chr16$snpinfo, drop_hilit=1.5, genes = LiverGSH_GSSGRatio_Genes_MGI_chr16, main = "Liver GSH/GSSG Genes MGI")

#For Liver GSH/GSSG --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 11
  coef_blup_LiverGSH_GSSGRatio_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,8)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSH GSSG Ratio GWAS - RankZ sexgen.pdf")
out_gwas_LiverGSH_GSSGRatio <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverGSH_GSSGRatio$lod, out_gwas_LiverGSH_GSSGRatio$snpinfo, altcol="green4", gap=0, main = "Liver GSH/GSSG Ratio GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverGSH_GSSGRatio_sexgen <- est_herit(pheno["zLiverGSH_GSSGRatio"], kinship_lmm, sexgen, cores = 10)
herit_LiverGSH_GSSGRatio_sex <- est_herit(pheno["zLiverGSH_GSSGRatio"], kinship_lmm, sex, cores = 10)



##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSH_GSSG Ratio Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_LiverGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gpx1 Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_LiverGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gclc Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_LiverGSH_GSSGRatio_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gclm Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_LiverGSH_GSSGRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gss Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_LiverGSH_GSSGRatio_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gsr Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()

######################

#set working directory
pdf(file = "Chr16 GSH-GSSG-Ratio Peak - GSH-and-GSSG - RankZ sexgen.pdf")
#saved to GSH GSSG Ratio results
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##


chr = 16
xlim <- c(1,9)
plot_coefCC(x = coef_blup_LiverGSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
plot_coefCC(x = coef_blup_LiverGSSG_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
plot_coefCC(x = coef_blup_LiverTotalGSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()



