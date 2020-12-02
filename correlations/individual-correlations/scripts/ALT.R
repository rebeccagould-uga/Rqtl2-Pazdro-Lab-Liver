#Plotting Spearmans Rank (rho) correlations 
#R01 GSH Liver Data


#Removed all outliers
#Saved data in CSV format
#Glutathione data for the Liver and using Adjusted GSSG values 


library(ggplot2)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
#ggpubr info: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r



data <- read.csv(file = "~/OneDrive - University of Georgia/Pazdro Lab/R01 Redox/Analysis and Results/Statistics/Correlation Plots/Raw data.csv")
head(data, 6)

#to look up a specific correlation
cor.test(x = data$LiverWeight, y = data$AST, method = "spearman", use = "complete.obs")


#setwd
pdf(file = "ALT-correlation-plots.pdf") #create a file called ALT-correlation-plots

p1 <- ggscatter(data, x = "ALT", y = "LiverWeight", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "Liver Weight (g)")

p2 <- ggscatter(data, x = "ALT", y = "LiverWeight", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "Liver Weight (g)")

#####

p3 <- ggscatter(data, x = "ALT", y = "AST", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "AST (U/L)")

p4 <- ggscatter(data, x = "ALT", y = "AST", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "AST (U/L)")

####

p5 <- ggscatter(data, x = "ALT", y = "BUN", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "BUN (mg/dL)")

p6 <- ggscatter(data, x = "ALT", y = "BUN", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "BUN (mg/dL)")

####

p7 <- ggscatter(data, x = "ALT", y = "Glucose", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "Glucose (mg/dL)")

p8 <- ggscatter(data, x = "ALT", y = "Glucose", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "Glucose (mg/dL)")

####

p9 <- ggscatter(data, x = "ALT", y = "NADH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "NADH (pmol/µg)")

p10 <- ggscatter(data, x = "ALT", y = "NADH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "NADH (pmol/µg)")

####

p11 <- ggscatter(data, x = "ALT", y = "NADP", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "NADP+ (pmol/µg)")

p12 <- ggscatter(data, x = "ALT", y = "NADP", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "NADP+ (pmol/µg)")

####

p13 <- ggscatter(data, x = "ALT", y = "NADPH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "NADPH (pmol/µg)")

p14 <- ggscatter(data, x = "ALT", y = "NADPH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "NADPH (pmol/µg)")

####

p15 <- ggscatter(data, x = "ALT", y = "NADPNADPHRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "NADP+/NADPH")

p16 <- ggscatter(data, x = "ALT", y = "NADPNADPHRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "NADP+/NADPH")

####

p17 <- ggscatter(data, x = "ALT", y = "GSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "GSH (nmol/mg)")

p18 <- ggscatter(data, x = "ALT", y = "GSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "GSH (nmol/mg)")

####

p19 <- ggscatter(data, x = "ALT", y = "GSSG", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "GSSG (nmol/mg)")

p20 <- ggscatter(data, x = "ALT", y = "GSSG", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "GSSG (nmol/mg)")

####

p21 <- ggscatter(data, x = "ALT", y = "TotalGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "Total GSH (nmol/mg)")

p22 <- ggscatter(data, x = "ALT", y = "TotalGSH", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "Total GSH (nmol/mg)")

####

p23 <- ggscatter(data, x = "ALT", y = "GSHGSSGRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          xlab = "ALT (U/L)", 
          ylab = "GSH/GSSG")

p24 <- ggscatter(data, x = "ALT", y = "GSHGSSGRatio", 
          color = "black", #color of points
          add = "reg.line", #add regression line
          add.params = list(color = "black", fill = "darkgray"), #customize regression line
          conf.int = TRUE, #add confidence interval
          cor.coef = TRUE, #add correlation coefficient
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.sep = "\n"),
          cor.method = "spearman",
          facet.by = "Sex",
          panel.labs = list(Sex = c("F", "M")),
          xlab = "ALT (U/L)", 
          ylab = "GSH/GSSG")

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, ncol = 1, nrow = 2)


# be sure to turn the graphics output off at the end!
dev.off()

pdf(file = "overall-ALT-correlation-plots.pdf")
ggarrange(p1, p3, p5, p7, p9, p11, p13, p15, p17, p19, p21, p23, ncol = 1, nrow = 2)
dev.off()

