
#R01 GSH DO Mapping - Liver RankZ

#Using miQTL package from Greg Keele to identify artificial peaks from missing allele data


options(stringsAsFactors = FALSE)
library(tidyverse)

## qtl2
library(qtl2)
#devtools::install_github("gkeele/miqtl")
library(miqtl)

#set wd
setwd("~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping")

genoprobs36 <- readRDS("~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/data/Pazdro_GigaMUGA_genoprobs_qced_36state_sorted.rds")
## Convert to miQTL
physical_map_df <- read.csv("~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/data/R01_GSH_DO_pmap_sorted.csv", header = TRUE) #physical map with sorted chromosome order
physical_map <- physical_map_df %>%
  qtl2convert::map_df_to_list(pos_column = "pos")
genetic_map_df <- read.csv("~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping//data/R01_GSH_DO_gmap_sorted.csv", header = TRUE) #genetic map with sorted chromosome order
genetic_map <- genetic_map_df %>%
  qtl2convert::map_df_to_list(pos_column = "pos")
holder_cross_object <- list(gmap = genetic_map,
                            pmap = physical_map)

# set wd

## Make genomecache for miQTL
## 13 GB (basically the 36-state and 8-state genoprobs combined)
convert.qtl2.to.HAPPY(qtl2.object = genoprobs36, 
                      cross.object = holder_cross_object,
                      HAPPY.output.path = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache", 
                      diplotype.order = "qtl2")

## Mapping in qtl2
## Compare with miQTL
genoprobs8 <- readRDS("data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")
pheno_dat <- read.csv("data/R01_GSH_DO_pheno.csv", header = TRUE)
pheno <- pheno_dat %>%
  mutate(zLiverGSH = rankZ(Liver_GSH), # add whichever phenotypes you want here
         zLiverAdjGSSG = rankZ(Liver_Adj_GSSG),
         zLiverTotalGSH = rankZ(Liver_Adj_Total_GSH),
         zLiverGSHGSSGRatio = rankZ(Liver_Adj_GSH_GSSG_Ratio),
         zLiverNADP = rankZ(Liver_NADP),
         zLiverNADPH = rankZ(Liver_NADPH),
         zLiverNADPHRatio = rankZ(Liver_NADP_NADPH_Ratio),
         zAST = rankZ(AST),
         zALT = rankZ(ALT),
         zBUN = rankZ(BUN)
         ) %>%
  column_to_rownames("id") %>%
  as.matrix
covar_dat <- read.csv("data/R01_GSH_DO_covar.csv", header = TRUE) %>%
  mutate(generation = as.factor(generation)) %>%
  column_to_rownames("id")
covar_mat <- cbind(model.matrix(~sex, data = covar_dat)[,-1, drop = FALSE],
                   model.matrix(~generation, data = covar_dat)[,-1, drop = FALSE])
rownames(covar_mat) <- rownames(covar_dat)

K <- calc_kinship(genoprobs8, type = "loco")
map_df <- read.csv("~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/data/R01_GSH_DO_pmap_sorted.csv", header = TRUE) 
map <- map_df %>%
  qtl2convert::map_df_to_list(pos_column = "pos")

qtl2_scan <- scan1(genoprobs = genoprobs8, 
                   pheno = pheno[,25:34], # change numbers to reflect phenotypes you want to scan
                   addcovar = covar_mat[,1,drop=FALSE],
                   kinship = K,
                   cores = 10)

plot(qtl2_scan, map)
find_peaks(qtl2_scan, map = map, threshold = 6, sort_by = "pos")

pdf(file = "QTL2 Plots - Pheno by Chr - RankZ sex.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
    plot(qtl2_scan, map, chr = 2, # change chromosome number as needed
         main = dimnames(qtl2_scan)[[2]][i],
         lodcolumn = i)
  
  }
  
  for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
    plot(qtl2_scan, map, chr = 8, # change chromosome number as needed
         main = dimnames(qtl2_scan)[[2]][i],
         lodcolumn = i)
  }
  
  for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
    plot(qtl2_scan, map, chr = 12, # change chromosome number as needed
         main = dimnames(qtl2_scan)[[2]][i],
         lodcolumn = i)
  }
  
  for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
    plot(qtl2_scan, map, chr = 14, # change chromosome number as needed
         main = dimnames(qtl2_scan)[[2]][i],
         lodcolumn = i)
  }
  
  for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
    plot(qtl2_scan, map, chr = 16, # change chromosome number as needed
         main = dimnames(qtl2_scan)[[2]][i],
         lodcolumn = i)
  }

dev.off()

## Mapping in miqtl
library(miqtl)

miqtl_pheno_dat <- pheno %>%
  as.data.frame %>%
  rownames_to_column("SUBJECT.NAME") %>%
  left_join(covar_mat %>%
              as.data.frame %>%
              rownames_to_column("SUBJECT.NAME"))

# to include generation as a covariate
# miqtl_pheno_dat <- pheno %>% 
#   as.data.frame %>%
#   rownames_to_column("SUBJECT.NAME") %>%
#   left_join(covar_mat %>%
#               as.data.frame %>%
#               rownames_to_column("SUBJECT.NAME") %>%
#               mutate(generation = as.factor(generation)))


#################################################################
## for RankZ sex peaks
#################################################################

########################

#set working directory
setwd("~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL")

pdf(file = "miQTL Plot GSH Chr2 - RankZ sex.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#GSH chr 2 scan
miqtl_GSHchr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                         data = miqtl_pheno_dat, 
                         formula = rankZ(Liver_GSH) ~ 1 + sexM, # change phenotype here
                         K = K[["2"]], # change chromosome here
                         model = "additive", 
                         use.multi.impute = FALSE, 
                         chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_GSHchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver GSH - miQTL (run like qtl2)")


miqtl_GSHchr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                         data = miqtl_pheno_dat, 
                         formula = rankZ(Liver_GSH) ~ 1 + sexM, # change phenotype here
                         K = K[["2"]],  # change chromosome here
                         model = "additive", 
                         use.multi.impute = TRUE, 
                         num.imp = 11,
                         chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_GSHchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver GSH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot GSSG Chr2 - RankZ sex.pdf")

#GSSG chr 2 scan
miqtl_GSSGchr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                           data = miqtl_pheno_dat, 
                                           formula = rankZ(Liver_Adj_GSSG) ~ 1 + sexM, # change phenotype here
                                           K = K[["2"]], # change chromosome here
                                           model = "additive", 
                                           use.multi.impute = FALSE, 
                                           chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_GSSGchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver GSSG - miQTL (run like qtl2)")


miqtl_GSSGchr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                            data = miqtl_pheno_dat, 
                                            formula = rankZ(Liver_Adj_GSSG) ~ 1 + sexM, # change phenotype here
                                            K = K[["2"]],  # change chromosome here
                                            model = "additive", 
                                            use.multi.impute = TRUE, 
                                            num.imp = 11,
                                            chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_GSSGchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver GSSG - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot Total GSH Chr2 - RankZ sex.pdf")

#Total GSH chr 2 scan
miqtl_TotalGSHchr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                           data = miqtl_pheno_dat, 
                                           formula = rankZ(Liver_Adj_Total_GSH) ~ 1 + sexM, # change phenotype here
                                           K = K[["2"]], # change chromosome here
                                           model = "additive", 
                                           use.multi.impute = FALSE, 
                                           chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_TotalGSHchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Total Liver GSH - miQTL (run like qtl2)")


miqtl_TotalGSHchr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                            data = miqtl_pheno_dat, 
                                            formula = rankZ(Liver_Adj_Total_GSH) ~ 1 + sexM, # change phenotype here
                                            K = K[["2"]],  # change chromosome here
                                            model = "additive", 
                                            use.multi.impute = TRUE, 
                                            num.imp = 11,
                                            chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_TotalGSHchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Total Liver GSH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot GSH-GSSG-Ratio Chr16 - RankZ sex.pdf")

#GSH/GSSG chr 16 scan
miqtl_GSHGSSGRatiochr16_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                      data = miqtl_pheno_dat, 
                                      formula = rankZ(Liver_Adj_GSH_GSSG_Ratio) ~ 1 + sexM, # change phenotype here
                                      K = K[["16"]], # change chromosome here
                                      model = "additive", 
                                      use.multi.impute = FALSE, 
                                      chr = 16) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_GSHGSSGRatiochr16_rqtl_scan, chr = 16, use.lod = TRUE, main = "Liver GSH/GSSG - miQTL (run like qtl2)")


miqtl_GSHGSSGRatiochr16_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                       data = miqtl_pheno_dat, 
                                       formula = rankZ(Liver_Adj_GSH_GSSG_Ratio) ~ 1 + sexM, # change phenotype here
                                       K = K[["16"]],  # change chromosome here
                                       model = "additive", 
                                       use.multi.impute = TRUE, 
                                       num.imp = 11,
                                       chr = 16) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_GSHGSSGRatiochr16_imput_scan, chr = 16, use.lod = TRUE, main = "Liver GSH/GSSG - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADP Chr2 - RankZ sex.pdf")

#NADP chr 2 scan
miqtl_NADPchr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                       data = miqtl_pheno_dat, 
                                       formula = rankZ(Liver_NADP) ~ 1 + sexM, # change phenotype here
                                       K = K[["2"]], # change chromosome here
                                       model = "additive", 
                                       use.multi.impute = FALSE, 
                                       chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver NADP - miQTL (run like qtl2)")


miqtl_NADPchr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(Liver_NADP) ~ 1 + sexM, # change phenotype here
                                        K = K[["2"]],  # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = TRUE, 
                                        num.imp = 11,
                                        chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver NADP - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADP Chr8 - RankZ sex.pdf")

#NADP chr 8 scan
miqtl_NADPchr8_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                       data = miqtl_pheno_dat, 
                                       formula = rankZ(Liver_NADP) ~ 1 + sexM, # change phenotype here
                                       K = K[["8"]], # change chromosome here
                                       model = "additive", 
                                       use.multi.impute = FALSE, 
                                       chr = 8) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPchr8_rqtl_scan, chr = 8, use.lod = TRUE, main = "Liver NADP - miQTL (run like qtl2)")


miqtl_NADPchr8_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(Liver_NADP) ~ 1 + sexM, # change phenotype here
                                        K = K[["8"]],  # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = TRUE, 
                                        num.imp = 11,
                                        chr = 8) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPchr8_imput_scan, chr = 8, use.lod = TRUE, main = "Liver NADP - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADPH Chr2 - RankZ sex.pdf")

#NADPH chr 2 scan
miqtl_NADPHchr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                       data = miqtl_pheno_dat, 
                                       formula = rankZ(Liver_NADPH) ~ 1 + sexM, # change phenotype here
                                       K = K[["2"]], # change chromosome here
                                       model = "additive", 
                                       use.multi.impute = FALSE, 
                                       chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPHchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver NADPH - miQTL (run like qtl2)")


miqtl_NADPHchr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(Liver_NADPH) ~ 1 + sexM, # change phenotype here
                                        K = K[["2"]],  # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = TRUE, 
                                        num.imp = 11,
                                        chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPHchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver NADPH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADPH Chr12 - RankZ sex.pdf")

#NADPH chr 12 scan
miqtl_NADPHchr12_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(Liver_NADPH) ~ 1 + sexM, # change phenotype here
                                        K = K[["12"]], # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = FALSE, 
                                        chr = 12) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPHchr12_rqtl_scan, chr = 12, use.lod = TRUE, main = "Liver NADPH - miQTL (run like qtl2)")


miqtl_NADPHchr12_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                         data = miqtl_pheno_dat, 
                                         formula = rankZ(Liver_NADPH) ~ 1 + sexM, # change phenotype here
                                         K = K[["12"]],  # change chromosome here
                                         model = "additive", 
                                         use.multi.impute = TRUE, 
                                         num.imp = 11,
                                         chr = 12) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPHchr12_imput_scan, chr = 12, use.lod = TRUE, main = "Liver NADPH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADPH Chr14 - RankZ sex.pdf")

#NADPH chr 14 scan
miqtl_NADPHchr14_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                       data = miqtl_pheno_dat, 
                                       formula = rankZ(Liver_NADPH) ~ 1 + sexM, # change phenotype here
                                       K = K[["14"]], # change chromosome here
                                       model = "additive", 
                                       use.multi.impute = FALSE, 
                                       chr = 14) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPHchr14_rqtl_scan, chr = 14, use.lod = TRUE, main = "Liver NADPH - miQTL (run like qtl2)")


miqtl_NADPHchr14_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(Liver_NADPH) ~ 1 + sexM, # change phenotype here
                                        K = K[["14"]],  # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = TRUE, 
                                        num.imp = 11,
                                        chr = 14) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPHchr14_imput_scan, chr = 14, use.lod = TRUE, main = "Liver NADPH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADP NADPH Ratio Chr2 - RankZ sex.pdf")

#NADPH chr 2 scan
miqtl_NADPNADPHRatiochr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(Liver_NADP_NADPH_Ratio) ~ 1 + sexM, # change phenotype here
                                        K = K[["2"]], # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = FALSE, 
                                        chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPNADPHRatiochr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (run like qtl2)")


miqtl_NADPNADPHRatiochr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                         data = miqtl_pheno_dat, 
                                         formula = rankZ(Liver_NADP_NADPH_Ratio) ~ 1 + sexM, # change phenotype here
                                         K = K[["2"]],  # change chromosome here
                                         model = "additive", 
                                         use.multi.impute = TRUE, 
                                         num.imp = 11,
                                         chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPNADPHRatiochr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADP NADPH Ratio Chr12 - RankZ sex.pdf")

#NADPH chr 12 scan
miqtl_NADPNADPHRatiochr12_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                                  data = miqtl_pheno_dat, 
                                                  formula = rankZ(Liver_NADP_NADPH_Ratio) ~ 1 + sexM, # change phenotype here
                                                  K = K[["12"]], # change chromosome here
                                                  model = "additive", 
                                                  use.multi.impute = FALSE, 
                                                  chr = 12) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPNADPHRatiochr12_rqtl_scan, chr = 12, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (run like qtl2)")


miqtl_NADPNADPHRatiochr12_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                                   data = miqtl_pheno_dat, 
                                                   formula = rankZ(Liver_NADP_NADPH_Ratio) ~ 1 + sexM, # change phenotype here
                                                   K = K[["12"]],  # change chromosome here
                                                   model = "additive", 
                                                   use.multi.impute = TRUE, 
                                                   num.imp = 11,
                                                   chr = 12) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPNADPHRatiochr12_imput_scan, chr = 12, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot NADP NADPH Ratio Chr14 - RankZ sex.pdf")

#NADPH chr 14 scan
miqtl_NADPNADPHRatiochr14_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                                  data = miqtl_pheno_dat, 
                                                  formula = rankZ(Liver_NADP_NADPH_Ratio) ~ 1 + sexM, # change phenotype here
                                                  K = K[["14"]], # change chromosome here
                                                  model = "additive", 
                                                  use.multi.impute = FALSE, 
                                                  chr = 14) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPNADPHRatiochr14_rqtl_scan, chr = 14, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (run like qtl2)")


miqtl_NADPNADPHRatiochr14_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                                   data = miqtl_pheno_dat, 
                                                   formula = rankZ(Liver_NADP_NADPH_Ratio) ~ 1 + sexM, # change phenotype here
                                                   K = K[["14"]],  # change chromosome here
                                                   model = "additive", 
                                                   use.multi.impute = TRUE, 
                                                   num.imp = 11,
                                                   chr = 14) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_NADPNADPHRatiochr14_imput_scan, chr = 14, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot AST Chr2 - RankZ sex.pdf")

#AST chr 2 scan
miqtl_ASTchr2_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                       data = miqtl_pheno_dat, 
                                       formula = rankZ(AST) ~ 1 + sexM, # change phenotype here
                                       K = K[["2"]], # change chromosome here
                                       model = "additive", 
                                       use.multi.impute = FALSE, 
                                       chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_ASTchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver AST - miQTL (run like qtl2)")


miqtl_ASTchr2_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                        data = miqtl_pheno_dat, 
                                        formula = rankZ(AST) ~ 1 + sexM, # change phenotype here
                                        K = K[["2"]],  # change chromosome here
                                        model = "additive", 
                                        use.multi.impute = TRUE, 
                                        num.imp = 11,
                                        chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_ASTchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver AST - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

#################################################################

pdf(file = "miQTL Plot AST Chr16 - RankZ sex.pdf")

#AST chr 16 scan
miqtl_ASTchr16_rqtl_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                                data = miqtl_pheno_dat, 
                                                formula = rankZ(AST) ~ 1 + sexM, # change phenotype here
                                                K = K[["16"]], # change chromosome here
                                                model = "additive", 
                                                use.multi.impute = FALSE, 
                                                chr = 16) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_ASTchr16_rqtl_scan, chr = 16, use.lod = TRUE, main = "Liver AST - miQTL (run like qtl2)")


miqtl_ASTchr16_imput_scan <- scan.h2lmm(genomecache = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/QTL Mapping/Liver/RankZ/RankZ - sex/miQTL/Pazdro_genomecache/", 
                                                 data = miqtl_pheno_dat, 
                                                 formula = rankZ(AST) ~ 1 + sexM, # change phenotype here
                                                 K = K[["16"]],  # change chromosome here
                                                 model = "additive", 
                                                 use.multi.impute = TRUE, 
                                                 num.imp = 11,
                                                 chr = 16) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_ASTchr16_imput_scan, chr = 16, use.lod = TRUE, main = "Liver AST - miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)

dev.off()

############################
#saving plots as JPEGs

#set WD

jpeg(file="MQTL RQTL Results.jpeg")

genome.plotter.chr(miqtl_GSHchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver GSH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_GSHchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver GSH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_GSSGchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver GSSG - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_GSSGchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver GSSG - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_TotalGSHchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Total Liver GSH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_TotalGSHchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Total Liver GSH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_GSHGSSGRatiochr16_rqtl_scan, chr = 16, use.lod = TRUE, main = "Liver GSH/GSSG - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_GSHGSSGRatiochr16_imput_scan, chr = 16, use.lod = TRUE, main = "Liver GSH/GSSG - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver NADP - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver NADP - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPchr8_rqtl_scan, chr = 8, use.lod = TRUE, main = "Liver NADP - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPchr8_imput_scan, chr = 8, use.lod = TRUE, main = "Liver NADP - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPHchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver NADPH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPHchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver NADPH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPHchr12_rqtl_scan, chr = 12, use.lod = TRUE, main = "Liver NADPH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPHchr12_imput_scan, chr = 12, use.lod = TRUE, main = "Liver NADPH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPHchr14_rqtl_scan, chr = 14, use.lod = TRUE, main = "Liver NADPH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPHchr14_imput_scan, chr = 14, use.lod = TRUE, main = "Liver NADPH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPNADPHRatiochr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPNADPHRatiochr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPNADPHRatiochr12_rqtl_scan, chr = 12, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPNADPHRatiochr12_imput_scan, chr = 12, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_NADPNADPHRatiochr14_rqtl_scan, chr = 14, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_NADPNADPHRatiochr14_imput_scan, chr = 14, use.lod = TRUE, main = "Liver NADP/NADPH - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_ASTchr2_rqtl_scan, chr = 2, use.lod = TRUE, main = "Liver AST - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_ASTchr2_imput_scan, chr = 2, use.lod = TRUE, main = "Liver AST - miQTL (11 imputations)", y.max.manual = 10)
genome.plotter.chr(miqtl_ASTchr16_rqtl_scan, chr = 16, use.lod = TRUE, main = "Liver AST - miQTL (run like qtl2)")
genome.plotter.chr(miqtl_ASTchr16_imput_scan, chr = 16, use.lod = TRUE, main = "Liver AST - miQTL (11 imputations)", y.max.manual = 10)

dev.off()


