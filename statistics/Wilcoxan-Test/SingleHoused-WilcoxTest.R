#Wilcox Test to see if single-housed mice have different phenotypes
#R01 GSH Data


#Removed all outliers
#Saved data in CSV format
#Glutathione data for the Liver using adjusted GSSG values 
#Glutathione data for the Kidney using unadjusted GSSG values

#load my data into the session using read.csv function
rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/Wilcoxan-Test/data-Wilcoxan-Test.csv", check.names = FALSE)

head(rawdata, 6)

rownames(rawdata) <- rawdata$Mouse
#checking pheno file
rawdata[1:10,]

#remove column 1-3 (Mouse ID, sex, and generation)
data = subset(rawdata, select = -c(1,3,4))
data[1:10,]

#single housed --> 1 = single-housed, 0 = not single housed

housing <- data$`Single-Cage`

wilcox.test(data$`Kidney GSH (nmol/mg)`~housing)
wilcox.test(data$`Kidney GSSG (nmol/mg)`~housing)
wilcox.test(data$`Kidney Total Glutathione (nmol/mg)`~housing)
wilcox.test(data$`Kidney GSH/GSSG`~housing)
wilcox.test(data$`Kidney Eh (mV)`~housing)
wilcox.test(data$`BUN`~housing)


boxplot(data$`Kidney GSH (nmol/mg)`~housing)
boxplot(data$`Kidney GSSG (nmol/mg)`~housing)
boxplot(data$`Kidney Total Glutathione (nmol/mg)`~housing)
boxplot(data$`Kidney GSH/GSSG`~housing)
boxplot(data$`Kidney Eh (mV)`~housing)
boxplot(data$`BUN`~housing)


####################################################
## exporting data
####################################################
#set working directory


library(writexl)
write_xlsx(list("Correlations" = correlations, 
                "P Values" = pvalues),
                "correlations and p values.xlsx")


