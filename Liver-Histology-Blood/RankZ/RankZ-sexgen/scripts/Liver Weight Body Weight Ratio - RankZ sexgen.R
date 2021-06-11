# R01 LiverWeightBodyWeight DO Mapping Code 
# Updated November 2020
# Becca Gould 

#LIVER HISTOLOGY AND BLOOD (AST AND ALT) MAPPING - Liver Weight/Body Weight

#Load in Liver-Histology-Blood-RankZ-SexGen.Rdata
#Run RankZ Transformation and Data Prep R Script before doing this**


setwd("/users/becca/R01_LiverWeightBodyWeight_DO_mapping_Liver/data")

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

qtlscan_LiverWeightBodyWeight <- scan1(genoprobs = probs, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco, addcovar = sexgen, cores=2)
perm_LiverWeightBodyWeight <- scan1perm(genoprobs = probs, pheno = pheno["zLiverWeightBodyWeight"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Liver Weight Body Weight Ratio QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_LiverWeightBodyWeight = summary(perm_LiverWeightBodyWeight, alpha = c(0.2, 0.1, 0.05))

plot_scan1(x = qtlscan_LiverWeightBodyWeight, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Weight/Body Weight", ylim = c(0,11))
abline(h = threshold_LiverWeightBodyWeight, col = c("purple", "red", "blue"), lwd = 2)
plot_scan1(x = qtlscan_LiverWeightBodyWeight, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Weight/Body Weight", ylim = c(0,11))
abline(h = threshold_LiverWeightBodyWeight, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")

#using gmap (cM)
find_peaks(scan1_output = qtlscan_LiverWeightBodyWeight, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverWeightBodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
gmap_peaksLiverWeightBodyWeight <- find_peaks(scan1_output = qtlscan_LiverWeightBodyWeight, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverWeightBodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(gmap_peaksLiverWeightBodyWeight)

#using pmap (Mbp)
find_peaks(scan1_output = qtlscan_LiverWeightBodyWeight, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverWeightBodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
peaksLiverWeightBodyWeight <- find_peaks(scan1_output = qtlscan_LiverWeightBodyWeight, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverWeightBodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
print(peaksLiverWeightBodyWeight)


write_xlsx(list("LiverWeightBodyWeight gmap (cM)" = gmap_peaksLiverWeightBodyWeight,
                "LiverWeightBodyWeight pmap (Mbp)" = peaksLiverWeightBodyWeight),
           "LiverWeightBodyWeight Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Liver Weight/Body Weight --- Chromosome 11
par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
chr = 11
coef_blup_LiverWeightBodyWeight_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(15,40)
plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr11, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 11
#could use ci_lo and ci_hi, but in this case I want a specific chr 2 position
#start = peaksLiverWeightBodyWeight[peaksLiverWeightBodyWeight$chr ==  chr,"ci_lo"]
#end = peaksLiverWeightBodyWeight[peaksLiverWeightBodyWeight$chr == chr, "ci_hi"] 

pander(peaksLiverWeightBodyWeight)
#based on peaksAST, peak of interest is ~56 Mbp
variants_LiverWeightBodyWeight_chr11 <- query_variants(chr, 54.5, 58)
out_snps_LiverWeightBodyWeight_chr11 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                      chr = chr, start = 54.5, end = 58, keep_all_snps = TRUE)
plot_snpasso(out_snps_LiverWeightBodyWeight_chr11$lod, out_snps_LiverWeightBodyWeight_chr11$snpinfo, main = "Liver Weight/Body Weight SNPs")

LiverWeightBodyWeight_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 54.5, end = 58)
plot(out_snps_LiverWeightBodyWeight_chr11$lod, out_snps_LiverWeightBodyWeight_chr11$snpinfo, drop_hilit=1.5, genes = LiverWeightBodyWeight_Genes_MGI_chr11, main = "Liver Weight/Body Weight Genes MGI")
plot_genes(LiverWeightBodyWeight_Genes_MGI_chr11, main = "Liver Weight/Body Weight Genes MGI")


#For Liver Weight/Body Weight --- Chromosome 12
par(mar=c(4.1, 4.1, 2.6, 2.6))

#using gmap (cM)
chr = 12
coef_blup_LiverWeightBodyWeight_chr12 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(50,70)
plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 12
#could use ci_lo and ci_hi, but in this case I want a specific chr 2 position
#start = peaksLiverWeightBodyWeight[peaksLiverWeightBodyWeight$chr ==  chr,"ci_lo"]
#end = peaksLiverWeightBodyWeight[peaksLiverWeightBodyWeight$chr == chr, "ci_hi"] 

pander(peaksLiverWeightBodyWeight)
#based on peaksAST, peak of interest is ~40 Mbp
variants_LiverWeightBodyWeight_chr12 <- query_variants(chr, 111.5, 114)
out_snps_LiverWeightBodyWeight_chr12 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                  chr = chr, start = 111.5, end = 114, keep_all_snps = TRUE)
plot_snpasso(out_snps_LiverWeightBodyWeight_chr12$lod, out_snps_LiverWeightBodyWeight_chr12$snpinfo, main = "Liver Weight/Body Weight SNPs")

LiverWeightBodyWeight_Genes_MGI_chr12 <- query_genes_mgi(chr = chr, start = 111.5, end = 114)
plot(out_snps_LiverWeightBodyWeight_chr12$lod, out_snps_LiverWeightBodyWeight_chr12$snpinfo, drop_hilit=1.5, genes = LiverWeightBodyWeight_Genes_MGI_chr12, main = "Liver Weight/Body Weight Genes MGI")
plot_genes(LiverWeightBodyWeight_Genes_MGI_chr12, main = "Liver Weight/Body Weight Genes MGI")

dev.off()


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Liver Weight Body Weight GWAS - RankZ sexgen.pdf")
out_gwas_LiverWeightBodyWeight <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_LiverWeightBodyWeight$lod, out_gwas_LiverWeightBodyWeight$snpinfo, altcol="green4", gap=0, main = "LiverWeightBodyWeight GWAS", ylim = c(0,6))
dev.off()

####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

herit_LiverWeightBodyWeight_sex <- est_herit(pheno["zLiverWeightBodyWeight"], kinship_lmm, sex, cores = 2)
herit_LiverWeightBodyWeight_sexgen <- est_herit(pheno["zLiverWeightBodyWeight"], kinship_lmm, sexgen, cores = 2)

