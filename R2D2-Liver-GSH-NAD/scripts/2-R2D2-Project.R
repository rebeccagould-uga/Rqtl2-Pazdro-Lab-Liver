# R2D2 Project
# Created by Becca Gould
# Updated May 2021

#### Phenotypes ####
# Chr 2 (WSB): GSH, GSSG, Total GSH, Eh, NADP, NADP/NADPH
# Chr 14 (NOD): NADPH, NADP/NADPH

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


#Load in R2D2-Liver-GSH-NAD-RankZ.Rdata
#Run 1-R2D2-Setup.R prior to this script




########################################################################################################
##
##
##
########################################################################################################

  # pull genotype frequencies 
  pull <- pull_genoprobpos(probs, map = R01_GSH_DO_QTLdata$gmap, 2, 47)
  head(pull)






########################################################################################################
## Testing regression matrix - is it ill-conditioned? 
## why the covariates affected the scan
########################################################################################################

  # instructions from Gary Churchill
  
  # kappa() computes the condition number of a matrix
  # Construct a matrix with the genoprobs at the maker locus.  It should have 8 columns A-H.  Compute kappa.
  # Repeat this but drop column H and replace it with a column of all ones – kappa should be the same.
  # Now add a column with 1 for females and 0 for males…compute kappa.
  # Next add several column with dummy variables for cohort…compute kappa.
  # Just to get a feel for things, swap in genotypes from loci where there is no problem.
  # You can compute the Sex and Cohort columns using the model.matrix() function,
  # e.g. model.matrix( ~Sex, data=my.dataframe) or model.matrix(~Sex+Cohort, data=my.dataframe)
  
  kappa()






########################################################################################################
## 
## Fit1 model
##
########################################################################################################

  #sex
  #BLUPs
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sex, model = "normal", se = TRUE, blup = TRUE, cores = 2)
  #default
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sex, model = "normal", se = TRUE, blup = FALSE, cores = 2)
  #sexgen  
  #BLUPs
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, model = "normal", se = TRUE, blup = TRUE, cores = 2)
  #default
  fit1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, model = "normal", se = TRUE, blup = FALSE, cores = 2)




########################################################################################################
## 
## Liver GSH (WSB - Chromosome 2)
##
########################################################################################################

# QTL scans - sex
  #qtlscan_LiverGSH_sex <- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sex, cores=2)
  #perm_LiverGSH_sex <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = sex, n_perm = 1000, cores=10)
  
  pdf(file = "GSH-QTL-Results-RankZ.pdf")

    threshold_LiverGSH_sex = summary(perm_LiverGSH_sex, alpha = 0.05)
    par(mar=c(4.1, 4.1, 2.6, 2.6))
    #plot
    plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex", ylim = c(0,11))
    abline(h = threshold_LiverGSH_sex, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver GSH - sex", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksGSH_sex <- find_peaks(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    print(gmap_peaksGSH_sex)
    #using pmap (Mbp)
    peaksGSH_sex <- find_peaks(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_sex, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    print(peaksGSH_sex)

  
# QTL scans - sex + gen
  #qtlscan_LiverGSH_sexgen <- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  #perm_LiverGSH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

    par(mar=c(4.1, 4.1, 2.6, 2.6))
    threshold_LiverGSH_sexgen = summary(perm_LiverGSH_sexgen, alpha = 0.05)
    #plot
    plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex + gen", ylim = c(0,11))
    #plot with sig thresholds
    plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, main = "Genome Scan for Liver GSH - sex + gen", ylim = c(0,11))
    abline(h = threshold_LiverGSH_sexgen, col = "blue", lwd = 2, lty = "dashed")
    #plot chr2 only
    plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, chr = 2, main = "Genome Scan for Liver GSH - sex + gen", ylim = c(0,11))
    
    #using gmap (cM)
    gmap_peaksGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    print(gmap_peaksGSH_sexgen)
    #using pmap (Mbp)
    peaksGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_sexgen, alpha = 0.05), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
    print(peaksGSH_sexgen)

  dev.off()
  

# founder allele effects 
  chr = 2
  
  #setwd
  pdf(file = "founder-allele-effects.pdf") #gmap

  #sex
    #BLUPs
    coef_blup_LiverGSH_chr2_sex <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sex, main = "Liver GSH BLUPs - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sex, main = "Liver GSH BLUPs - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverGSH_chr2_sex <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 2)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver GSH default - sex ", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sex, map = R01_GSH_DO_QTLdata$gmap[chr], columns = 1:8, main = "Liver GSH default - sex ", legend = "bottomleft", bgcolor="gray95", xlim=xlim)
    
  #sexgen
    #BLUPs
    coef_blup_LiverGSH_chr2_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_blup_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH BLUPs - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
    
    #default
    coef_def_LiverGSH_chr2_sexgen <- scan1coef(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH default - sex + gen", legend = "bottomleft", bgcolor="gray95")
    xlim <- c(45,65)
    plot_coefCC(x = coef_def_LiverGSH_chr2_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sexgen, main = "Liver GSH default - sex + gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  dev.off()
    
    
  
  
  
  
  
########################################################################################################
##
##Allele Plots
##
########################################################################################################
    
# Identify QTL scan names
  ls(pattern = "qtl")
    
# use cbind to combine all of the qtlscans + take those 2D tables and combining them with another table over and over again
# scans is an R object containing your genome scans from scan1() that are loaded in to the R environment 
  scans_sex <- cbind(qtlscan_LiverGSH_sex, qtlscan_LiverGSSG_sex, qtlscan_LiverTotalGSH_sex, qtlscan_LiverEh_sex, qtlscan_LiverNADP_sex, qtlscan_LiverNADPH_sex, qtlscan_LiverNADP_NADPHRatio_sex)
  scans_sexgen <- cbind(qtlscan_LiverGSH_sexgen, qtlscan_LiverGSSG_sexgen, qtlscan_LiverTotalGSH_sexgen, qtlscan_LiverEh_sexgen, qtlscan_LiverNADP_sexgen, qtlscan_LiverNADPH_sexgen, qtlscan_LiverNADP_NADPHRatio_sexgen)
  
# Probability plotting function (created by Greg Keele)
    
    prob_plot <- function(pheno_vec,
                          pheno_name = NULL,
                          genoprobs,
                          qtl_chr,
                          qtl_marker,
                          cols = gray(10000:1/10000),
                          label_col = as.character(qtl2::CCcolors),
                          founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                          main = "") {
      sorted_pheno <- sort(pheno_vec)
      image(genoprobs[[qtl_chr]][names(sorted_pheno), rev(LETTERS[1:8]), qtl_marker] * 2,
            yaxt = "n", xaxt = "n", col = cols)
      axis(2, at = seq(0, 8, 1 + 1/8)/8, labels = FALSE,
           lty = 0, srt = 90, las = 2)
      mtext(text = main, side = 3, padj = -1, cex = 1.25)
      mtext(text = rev(founders), side = 2, col = rev(label_col), at = seq(0, 1, length.out = 8),
            las = 1, cex = 1.25, adj = 1.25, font = 2)
      mtext(text = paste("lower", "<--", ifelse(is.null(pheno_name), "phenotype", pheno_name), "-->", "higher"), side = 1, padj = 1.25, cex = 1.25)
    }
    
    
# Alter the pheno file from a data frame into a matrix
  
  #first need to identify what specifically to make part of the phenotype matrix (only need transformed data!)
  names(pheno)

  #pheno_mat is the matrix of outcomes (phenotypes)
  pheno_mat <- as.matrix(pheno[c(36:42)])
    
  #check rownames to make sure they are already set as the write row names
  rownames(pheno[c(36:42)])
    
    
# Gather QTL peaks from all scans with LOD scores > 6

  qtl_gmap_sex <- find_peaks(scans_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_gmap_sex
  
  qtl_gmap_sexgen <- find_peaks(scans_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_gmap_sexgen
    
  qtl_pmap_sex <- find_peaks(scans_sex, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_pmap_sex
  
  qtl_pmap_sexgen <- find_peaks(scans_sexgen, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_pmap_sexgen
    
# Add marker information
  qtl_gmap_sex$marker.id <- find_marker(map = R01_GSH_DO_QTLdata$gmap, chr = qtl_gmap_sex$chr, pos = qtl_gmap_sex$pos)
  qtl_gmap_sexgen$marker.id <- find_marker(map = R01_GSH_DO_QTLdata$gmap, chr = qtl_gmap_sexgen$chr, pos = qtl_gmap_sexgen$pos)
  qtl_pmap_sex$marker.id <- find_marker(map = R01_GSH_DO_QTLdata$pmap, chr = qtl_pmap_sex$chr, pos = qtl_pmap_sex$pos)
  qtl_pmap_sexgen$marker.id <- find_marker(map = R01_GSH_DO_QTLdata$pmap, chr = qtl_pmap_sexgen$chr, pos = qtl_pmap_sexgen$pos)
 
  #setwd
  write_xlsx(list("QTL List RankZ Sex - cM" = qtl_gmap_sex,
                  "QTL List RankZ Sex - Mbp" = qtl_pmap_sex,
                  "QTL List RankZ SexGeb=n - cM" = qtl_gmap_sexgen,
                  "QTL List RankZ SexGen - Mbp" = qtl_pmap_sexgen),
               "QTL List LOD Score 6.xlsx")
 
     
# Allele plots loop - Sex
  
  #setwd
  pdf(file = "allele-plots-cM-RankZ-sex.pdf") #gmap
    for (i in 1:dim(qtl_gmap_sex)[1]) {
      prob_plot(pheno_vec = pheno_mat[,qtl_gmap_sex$lodcolumn[i]],
                genoprobs = probs,
                qtl_chr = qtl_gmap_sex$chr[i],
                qtl_marker = qtl_gmap_sex$marker.id[i],
                main = paste("lodindex", qtl_gmap_sex$lodindex[i], "Chr", qtl_gmap_sex$chr[i], qtl_gmap_sex$marker.id[i], qtl_gmap_sex$pos[i], qtl_gmap_sex$lodcolumn[i]))
    }
    dev.off()  
    
    pdf(file = "allele-plots-Mbp-RankZ-sex.pdf") #pmap
    for (i in 1:dim(qtl_pmap_sex)[1]) {
      prob_plot(pheno_vec = pheno_mat[,qtl_pmap_sex$lodcolumn[i]],
                genoprobs = probs,
                qtl_chr = qtl_pmap_sex$chr[i],
                qtl_marker = qtl_pmap_sex$marker.id[i],
                main = paste(qtl_pmap_sex$lodindex[i], "Chr", qtl_pmap_sex$chr[i], qtl_pmap_sex$marker.id[i], qtl_pmap_sex$pos[i], qtl_pmap_sex$lodcolumn[i]))
    }
    dev.off()  

    
    
# Allele plots - SexGen
    
    #setwd 
    pdf(file = "allele-plots-cM-RankZ-sexgen.pdf") #gmap
    for (i in 1:dim(qtl_gmap_sexgen)[1]) {
      prob_plot(pheno_vec = pheno_mat[,qtl_gmap_sexgen$lodcolumn[i]],
                genoprobs = probs,
                qtl_chr = qtl_gmap_sexgen$chr[i],
                qtl_marker = qtl_gmap_sexgen$marker.id[i],
                main = paste("lodindex", qtl_gmap_sexgen$lodindex[i], "Chr", qtl_gmap_sexgen$chr[i], qtl_gmap_sexgen$marker.id[i], qtl_gmap_sexgen$pos[i], qtl_gmap_sexgen$lodcolumn[i]))
    }
    dev.off()  
    
    pdf(file = "allele-plots-Mbp-RankZ-sexgen.pdf") #pmap
    for (i in 1:dim(qtl_pmap_sexgen)[1]) {
      prob_plot(pheno_vec = pheno_mat[,qtl_pmap_sexgen$lodcolumn[i]],
                genoprobs = probs,
                qtl_chr = qtl_pmap_sexgen$chr[i],
                qtl_marker = qtl_pmap_sexgen$marker.id[i],
                main = paste(qtl_pmap_sexgen$lodindex[i], "Chr", qtl_pmap_sexgen$chr[i], qtl_pmap_sexgen$marker.id[i], qtl_pmap_sexgen$pos[i], qtl_pmap_sexgen$lodcolumn[i]))
    }
    dev.off()  
    
