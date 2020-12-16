#Plotting Spearmans Rank (rho) correlations 
#R01 GSH Liver Data


#Removed all outliers
#Saved data in CSV format
#Glutathione data for the Liver and using Adjusted GSHGSSGRatio values 


library(ggplot2)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library(ggpubr)
#ggpubr info: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r



data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/correlations/data.csv")
head(data, 6)


#setwd
pdf(file = "GSHGSSGRatio-correlation-plots.pdf") #create a file called GSHGSSGRatio-correlation-plots

p1 <- ggscatter(data, x = "GSHGSSGRatio", y = "LiverWeight", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "Liver Weight (g)")

p2 <- ggscatter(data, x = "GSHGSSGRatio", y = "LiverWeight", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "Liver Weight (g)")

#####

p3 <- ggscatter(data, x = "GSHGSSGRatio", y = "AST", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "AST (U/L)")

p4 <- ggscatter(data, x = "GSHGSSGRatio", y = "AST", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "AST (U/L)")

####

p5 <- ggscatter(data, x = "GSHGSSGRatio", y = "ALT", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "ALT (U/L)")

p6 <- ggscatter(data, x = "GSHGSSGRatio", y = "ALT", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "ALT (U/L)")

####

p7 <- ggscatter(data, x = "GSHGSSGRatio", y = "BUN", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "BUN (mg/dL)")

p8 <- ggscatter(data, x = "GSHGSSGRatio", y = "BUN", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "BUN (mg/dL)")

####

p9 <- ggscatter(data, x = "GSHGSSGRatio", y = "Glucose", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "Glucose (mg/dL)")

p10 <- ggscatter(data, x = "GSHGSSGRatio", y = "Glucose", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "Glucose (mg/dL)")

####

p11 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "NADH (pmol/µg)")

p12 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "NADH (pmol/µg)")

####

p13 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADP", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "NADP+ (pmol/µg)")

p14 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADP", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "NADP+ (pmol/µg)")

####

p15 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADPH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "NADPH (pmol/µg)")

p16 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADPH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "NADPH (pmol/µg)")

####

p17 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADPNADPHRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "NADP+/NADPH")

p18 <- ggscatter(data, x = "GSHGSSGRatio", y = "NADPNADPHRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "NADP+/NADPH")

####

p19 <- ggscatter(data, x = "GSHGSSGRatio", y = "GSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "GSH (nmol/mg)")

p20 <- ggscatter(data, x = "GSHGSSGRatio", y = "GSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "GSH (nmol/mg)")

####

p21 <- ggscatter(data, x = "GSHGSSGRatio", y = "GSSG", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "GSSG (nmol/mg)")

p22 <- ggscatter(data, x = "GSHGSSGRatio", y = "GSSG", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "GSSG (nmol/mg)")

####

p23 <- ggscatter(data, x = "GSHGSSGRatio", y = "TotalGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "GSH/GSSG", 
          ylab = "Total GSH (nmol/mg)")

p24 <- ggscatter(data, x = "GSHGSSGRatio", y = "TotalGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "GSH/GSSG", 
          ylab = "Total GSH (nmol/mg)")


####

p25 <- ggscatter(data, x = "GSHGSSGRatio", y = "EhGSSG2GSH", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "GSH/GSSG", 
                 ylab = "Eh GSSG/2GSH")

p26 <- ggscatter(data, x = "GSHGSSGRatio", y = "EhGSSG2GSH", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "GSH/GSSG", 
                 ylab = "Eh GSSG/2GS")


ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, ncol = 1, nrow = 2)


# be sure to turn the graphics output off at the end!
dev.off()

pdf(file = "overall-GSHGSSGRatio-correlation-plots.pdf")
ggarrange(p1, p3, p5, p7, p9, p11, p13, p15, p17, p19, p21, p23, p25, ncol = 1, nrow = 2)
dev.off()



