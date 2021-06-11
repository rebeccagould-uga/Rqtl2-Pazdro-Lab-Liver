#transferring files from Liver-GSH-NAD-RankZ-Sex.Rdata and Liver-GSH-NAD-RankZ-SexGen.Rdata

#Liver-GSH-NAD

#RankZ sex
saveRDS(qtlscan_LiverGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSH_sex.rds")
saveRDS(qtlscan_LiverGSSG, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSSG_sex.rds")
saveRDS(qtlscan_LiverTotalGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverTotalGSH_sex.rds")
saveRDS(qtlscan_LiverRedoxPotentialGSSG2GSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverEh_sex.rds")
saveRDS(qtlscan_LiverNADP, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_sex.rds")
saveRDS(qtlscan_LiverNADPH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADPH_sex.rds")
saveRDS(qtlscan_LiverNADP_NADPHRatio, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_NADPHRatio_sex.rds")

saveRDS(perm_LiverGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSH_sex.rds")
saveRDS(perm_LiverGSSG, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSSG_sex.rds")
saveRDS(perm_LiverTotalGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverTotalGSH_sex.rds")
saveRDS(perm_LiverRedoxPotentialGSSG2GSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverEh_sex.rds")
saveRDS(perm_LiverNADP, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_sex.rds")
saveRDS(perm_LiverNADPH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADPH_sex.rds")
saveRDS(perm_LiverNADP_NADPHRatio, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_NADPHRatio_sex.rds")

#########

qtlscan_LiverGSH_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSH_sex.rds")
qtlscan_LiverGSSG_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSSG_sex.rds")
qtlscan_LiverTotalGSH_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverTotalGSH_sex.rds")
qtlscan_LiverEh_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverEh_sex.rds")
qtlscan_NADP_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_sex.rds")
qtlscan_NADPH_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADPH_sex.rds")
qtlscan_NADP_NADPHRatio_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_NADPHRatio_sex.rds")

perm_LiverGSH_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSH_sex.rds")
perm_LiverGSSG_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSSG_sex.rds")
perm_LiverTotalGSH_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverTotalGSH_sex.rds")
perm_LiverEh_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverEh_sex.rds")
perm_NADP_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_sex.rds")
perm_NADPH_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADPH_sex.rds")
perm_NADP_NADPHRatio_sex <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_NADPHRatio_sex.rds")


################################


#RankZ sexgen
saveRDS(qtlscan_LiverGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSH_sexgen.rds")
saveRDS(qtlscan_LiverGSSG, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSSG_sexgen.rds")
saveRDS(qtlscan_LiverTotalGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverTotalGSH_sexgen.rds")
saveRDS(qtlscan_LiverRedoxPotentialGSSG2GSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverEh_sexgen.rds")
saveRDS(qtlscan_LiverNADP, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_sexgen.rds")
saveRDS(qtlscan_LiverNADPH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADPH_sexgen.rds")
saveRDS(qtlscan_LiverNADP_NADPHRatio, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_NADPHRatio_sexgen.rds")

saveRDS(perm_LiverGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSH_sexgen.rds")
saveRDS(perm_LiverGSSG, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSSG_sexgen.rds")
saveRDS(perm_LiverTotalGSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverTotalGSH_sexgen.rds")
saveRDS(perm_LiverRedoxPotentialGSSG2GSH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverEh_sexgen.rds")
saveRDS(perm_LiverNADP, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_sexgen.rds")
saveRDS(perm_LiverNADPH, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADPH_sexgen.rds")
saveRDS(perm_LiverNADP_NADPHRatio, file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_NADPHRatio_sexgen.rds")

########

qtlscan_LiverGSH_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSH_sexgen.rds")
qtlscan_LiverGSSG_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverGSSG_sexgen.rds")
qtlscan_LiverTotalGSH_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverTotalGSH_sexgen.rds")
qtlscan_LiverEh_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_LiverEh_sexgen.rds")
qtlscan_NADP_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_sexgen.rds")
qtlscan_NADPH_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADPH_sexgen.rds")
qtlscan_NADP_NADPHRatio_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/qtlscan_NADP_NADPHRatio_sexgen.rds")

perm_LiverGSH_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSH_sexgen.rds")
perm_LiverGSSG_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverGSSG_sexgen.rds")
perm_LiverTotalGSH_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverTotalGSH_sexgen.rds")
perm_LiverEh_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_LiverEh_sexgen.rds")
perm_NADP_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_sexgen.rds")
perm_NADPH_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADPH_sexgen.rds")
perm_NADP_NADPHRatio_sexgen <- readRDS(file = "~/Rqtl2-Glutathione-Genetics/R2D2-Liver-GSH-NAD/files/perm_NADP_NADPHRatio_sexgen.rds")






