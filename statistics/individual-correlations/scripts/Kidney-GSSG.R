#Plotting Spearmans Rank (rho) correlations 
#R01 GSH Renal Data


#Removed all outliers
#Saved data in CSV format


library(ggplot2)

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)
#ggpubr info: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r



data <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data-Kidney-QTL-paper-fixednames.csv")
head(data)


#setwd
pdf(file = "Kidney-GSSG-correlation-plots.pdf") 

p1 <- ggscatter(data, x = "KidneyGSSG", y = "LiverWeight", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Liver Weight (g)")

p2 <- ggscatter(data, x = "KidneyGSSG", y = "LiverWeight", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Liver Weight (g)")

#####

p3 <- ggscatter(data, x = "KidneyGSSG", y = "AST", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "AST (U/L)")

p4 <- ggscatter(data, x = "KidneyGSSG", y = "AST", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "AST (U/L)")

####

p5 <- ggscatter(data, x = "KidneyGSSG", y = "ALT", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "ALT (U/L)")

p6 <- ggscatter(data, x = "KidneyGSSG", y = "ALT", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "ALT (U/L)")

####

p7 <- ggscatter(data, x = "KidneyGSSG", y = "BUN", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "BUN (mg/dL)")

p8 <- ggscatter(data, x = "KidneyGSSG", y = "BUN", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "BUN (mg/dL)")

####

p9 <- ggscatter(data, x = "KidneyGSSG", y = "Glucose", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Glucose (mg/dL)")

p10 <- ggscatter(data, x = "KidneyGSSG", y = "Glucose", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Glucose (mg/dL)")

####

p11 <- ggscatter(data, x = "KidneyGSSG", y = "NADH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADH (pmol/µg)")

p12 <- ggscatter(data, x = "KidneyGSSG", y = "NADH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADH (pmol/µg)")

####

p13 <- ggscatter(data, x = "KidneyGSSG", y = "NADP", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADP+ (pmol/µg)")

p14 <- ggscatter(data, x = "KidneyGSSG", y = "NADP", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADP+ (pmol/µg)")

####

p15 <- ggscatter(data, x = "KidneyGSSG", y = "NADPH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADPH (pmol/µg)")

p16 <- ggscatter(data, x = "KidneyGSSG", y = "NADPH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADPH (pmol/µg)")

####

p17 <- ggscatter(data, x = "KidneyGSSG", y = "NADPNADPHRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADP+/NADPH")

p18 <- ggscatter(data, x = "KidneyGSSG", y = "NADPNADPHRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "NADP+/NADPH")

####

p19 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Renal GSH (nmol/mg)")

p20 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Renal GSH (nmol/mg)")

####

p21 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyTotalGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Renal Total GSH (nmol/mg)")

p22 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyTotalGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Renal Total GSH (nmol/mg)")

####

p23 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyGSHGSSGRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Renal GSH/GSSG")

p24 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyGSHGSSGRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "Renal GSSG (nmol/mg)", 
          ylab = "Renal GSH/GSSG")


####

p25 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyEh", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Renal Eh (mV)")

p26 <- ggscatter(data, x = "KidneyGSSG", y = "KidneyEh", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Renal Eh (mV)")

####

p27 <- ggscatter(data, x = "KidneyGSSG", y = "LiverGSH", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic GSH (nmol/mg)")

p28 <- ggscatter(data, x = "KidneyGSSG", y = "LiverGSH", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic GSH (nmol/mg)")

####

p29 <- ggscatter(data, x = "KidneyGSSG", y = "LiverGSSG", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic GSSG (nmol/mg)")

p30 <- ggscatter(data, x = "KidneyGSSG", y = "LiverGSSG", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic GSSG (nmol/mg)")

####

p31 <- ggscatter(data, x = "KidneyGSSG", y = "LiverTotalGSH", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic Total GSH (nmol/mg)")

p32 <- ggscatter(data, x = "KidneyGSSG", y = "LiverTotalGSH", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic Total GSH (nmol/mg)")

####

p33 <- ggscatter(data, x = "KidneyGSSG", y = "LiverGSHGSSGRatio", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic GSH/GSSG")

p34 <- ggscatter(data, x = "KidneyGSSG", y = "LiverGSHGSSGRatio", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic GSH/GSSG")


####

p35 <- ggscatter(data, x = "KidneyGSSG", y = "LiverEh", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic Eh (mV)")

p36 <- ggscatter(data, x = "KidneyGSSG", y = "LiverEh", 
                 color = "black", #color of points
                 add = "reg.line", #add regression line
                 add.params = list(color = "black", fill = "darkgray"), #customize regression line
                 conf.int = TRUE, #add confidence interval
                 cor.coef = TRUE, #add correlation coefficient
                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
                 cor.method = "spearman",
                 facet.by = "Sex",
                 panel.labs = list(Sex = c("F", "M")),
                 xlab = "Renal GSSG (nmol/mg)", 
                 ylab = "Hepatic Eh (mV)")

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, ncol = 1, nrow = 2)


# be sure to turn the graphics output off at the end!
dev.off()

pdf(file = "overall-KidneyGSSG-correlation-plots.pdf")
ggarrange(p1, p3, p5, p7, p9, p11, p13, p15, p17, p19, p21, p23, p25, p27, p29, p31, p33, p35, ncol = 1, nrow = 1)
dev.off()

