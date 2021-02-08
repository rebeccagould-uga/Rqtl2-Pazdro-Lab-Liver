#Plotting Spearmans Rank (rho) correlations in Correlation Matrix Heat Map
#R01 GSH Liver Data


#Removed all outliers
#Saved data in CSV format
#Glutathione data for the Liver  using adjusted GSSG values 
#Glutathione data for the Kidney using unadjusted GSSG values

#plot info reference: https://www.displayr.com/how-to-create-a-correlation-matrix-in-r/

#load my data into the session using read.csv function
  #rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data.csv")
  #rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data-Liver-GSH-QTL-paper.csv", check.names = FALSE)
  rawdata <- read.csv(file = "~/Rqtl2-Glutathione-Genetics/statistics/data-Kidney-QTL-paper.csv", check.names = FALSE)

head(rawdata, 6)

rownames(rawdata) <- rawdata$Mouse
#checking pheno file
rawdata[1:10,]

#remove column 1-3 (Mouse ID, sex, and generation)
data = subset(rawdata, select = -c(1,2,3))

#run correlations for all observations
#use complete.obs tells R to handle missing values by case-wise deletion
data.cor = cor(data, method = c("spearman"), use = "complete.obs")

  #to look up a specific correlation
  #cor.test(x = data$LiverWeight, y = data$AST, method = "spearman", use = "complete.obs")


#install.packages("Hmisc")
library("Hmisc")

#Use the following code to run the correlation matrix with p-values
#Note that the data has to be fed to the rcorr function as a matrix
data_correlations = rcorr(as.matrix(data), type = c("spearman"))
data_correlations

#This generates one table of correlation coefficients (the correlation matrix) and another table of the p-values. 
#By default, the correlations and p-values are stored in an object of class type rcorr. 
#To extract the values from this object into a useable data structure, you can use the following syntax:

data_coefficients = data_correlations$r
data_p = data_correlations$P


#This section provides a simple function for formatting a correlation matrix into a table with 4 columns containing :
  #Column 1 : row names (variable 1 for the correlation test)
  #Column 2 : column names (variable 2 for the correlation test)
  #Column 3 : the correlation coefficients
  #Column 4 : the p-values of the correlations
#The custom function below can be used :

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flattenCorrMatrix(data_coefficients, data_p)


#visualize the correlation matrix
#install.packages("corrplot")
library(corrplot)


#helpful link: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software

#setwd to save files

#to save as PDF, do this instead:
#pdf(file = "correlation-matrix-plots.pdf") #create a file called correlation-matrix-plots

png(height=1500, width=1500, pointsize=25, file="correlation-matrix-number.png") #create a PNG file called correlation-matrix-number
  corrplot(data_coefficients, 
         method = "number", #could do circle, number, or pie as well
         type = "lower", #full correlation matrix, upper or lower
         #order = "hclust", #heriarchical clustering order
         #title = "Correlation Matrix",
         #add coefficient of correlation
         #addCoef.col = "black",  #can't do for this type of plot
         #addCoefasPercent =TRUE, #can print coefficients as percents for space-saving
         #customize text labels 
         tl.col = "black", #text label color
         tl.srt = 45, #text string rotation
         #combine with significance
         p.mat = data_p, 
         sig.level = 0.05, 
         insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)
  dev.off()

  
png(height=1500, width=1500, pointsize=25, file="correlation-matrix-circle.png") #create a PNG file called correlation-matrix-circle
  corrplot(data_coefficients, 
           method = "circle", #could do circle, number, or pie as well
           type = "lower", #full correlation matrix, upper or lower
           #order = "hclust", #heriarchical clustering order
           #title = "Correlation Matrix",
           #add coefficient of correlation
           addCoef.col = "black", #can remove to just show circles
           #addCoefasPercent =TRUE, #can print coefficients as percents for space-saving
           #customize text labels 
           tl.col = "black", #text label color
           tl.srt = 45, #text string rotation
           #combine with significance
           p.mat = data_p, 
           sig.level = 0.05, 
           insig = "blank",
           # hide correlation coefficient on the principal diagonal
           diag=FALSE)
  dev.off()

  
png(height=1500, width=1500, pointsize=25, file="correlation-matrix-color.png") #create a PNG file called correlation-matrix-color
  corrplot(data_coefficients, 
           method = "color", #could do circle, number, or pie as well
           type = "lower", #full correlation matrix, upper or lower
           #order = "hclust", #heriarchical clustering order
           #title = "Correlation Matrix",
           #add coefficient of correlation
           addCoef.col = "black", 
           #addCoefasPercent =TRUE, #can print coefficients as percents for space-saving
           #customize text labels 
           tl.col = "black", #text label color
           tl.srt = 45, #text string rotation
           #combine with significance
           p.mat = data_p, 
           sig.level = 0.05, 
           insig = "blank",
           # hide correlation coefficient on the principal diagonal
           diag=FALSE)
  dev.off()
  
  #In the above plot, correlations with p-value > 0.05 are considered as insignificant. 
  #In this case the correlation coefficient values are left blank.
  #could also add crosses instead

  #A default correlation matrix plot (called a Correlogram) is generated. 
  #Positive correlations are displayed in a blue scale while negative correlations are displayed in a red scale.

  
#######
  #We can also generate a Heatmap object again using our correlation coefficients as input to the Heatmap. 
  #Because the default Heatmap color scheme is quite unsightly, we can first specify a color palette to use in the Heatmap. 
  #The value at the end of the function specifies the amount of variation in the color scale. Typically no more than 20 is needed here. 
  #We then use the heatmap function to create the output:
  
png(height=2500, width=2500, pointsize=25, file="correlation-matrix-heatmap.png") #create a PNG file called correlation-matrix-heatmap
  palette = colorRampPalette(c("green", "white", "red")) (20)
  heatmap(x = data.cor, col = palette, symm = TRUE)

dev.off()



####################################################
## exporting data
####################################################
#set working directory

pvalues <- as.data.frame(data_p)
correlations <- as.data.frame(data.cor)

library(writexl)
write_xlsx(list("Correlations" = correlations, 
                "P Values" = pvalues),
                "correlations and p values.xlsx")


