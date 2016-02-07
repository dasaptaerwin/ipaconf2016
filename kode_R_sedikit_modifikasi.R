#####################################################
#####################################################
#####################################################
#####################################################

# Paper: IPA conf 2016
# Data: Conti et al. (2000)
# Team leader: Dasapta Erwin Irawan
# Code and analysis: Dasapta Erwin Irawan and Prana Ugi

#####################################################
#####################################################
#####################################################
#####################################################

# Note: remove the '#' symbol only if you haven't done it previously


# Set up working folder

# install packages
## from CRAN server
#install.packages("cluster") 
#install.packages("rpart")
#install.packages("vegan")
#install.packages("gam")

#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")

## from Bioconductor server
#pcaMethods package from Bioconductor server
#source("http://bioconductor.org/biocLite.R") 
#biocLite("pcaMethods")

# Loading library
library(pcaMethods) # for pcaMethods package
library(cluster) # for cluster analysis
#library(readxl) # for opening data directly from xls format
library(rpart)
library(vegan)
library(gam)
library(psych) # univar description



# Load data
df <- read.csv("dataPoplain.csv", header=T) # load data 
df
row.names(df) <- df$sites # setting row names
df <- df[4:17]
df
dim(df)
is.na(df) # checking NA if any
#df <- na.omit(df) # omitting NA if any


# Exploratory using pairs() function
# Assesing data patterns
cor.tab <- cor(df)
write.csv(cor.tab, "cortab.csv")
univar.tab <- describe(df) # from psych package
univar.tab
write.csv(univar.tab, "univar.csv")

pairs(df,
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)


# Run PCA 
# Run PCA (using pcamethods package)
## svdImpute = standard pca, with imputation, standardised, method univariate (uv)
pca <- pca(df, 
           method = "svdImpute", 
           scale = "uv",
           center = T,
           nPcs = 5,
           evalPcs = 1:5)
summary(pca)



## Evaluating results
loadings(pca) # loadings of each variables
write.csv(loadings(pca), "loadings.csv")
plotPcs(pca, Pcs = 1:5)
plot(pca)
slplot(pca) # default function in pcamethods but not big enough


# PCA using'vegan' package
df <- scale(df)  #standardize the variables
vegan.pca <- rda(df)
vegan.pca
plot(vegan.pca)
biplot(vegan.pca, scaling = -1, display = 'species')
biplot(vegan.pca, scaling = -1, display = 'sites')
dev.off()

ordilabel(plot(vegan.pca), display="species", font=1, col="gray70") # Add some frame on the label
orditkplot(plot(vegan.pca)) # you see that you can move the labels.

# reading "cleanplot.pca" function to make better visualisation
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:cleanplot.pca?do=export_code&codeblock=0')
cleanplot.pca (vegan.pca)

# reading "evplot" function for better visualisation
source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=0")
# select the data frame with eigenvalues of particular axes:
ev <- vegan.pca$CA$eig
evplot(ev) # calculate axis-importance and draw the barplots



# Cluster analysis (using 'vegan' package)
dis <- vegdist(df, method = 'bray') # percentage cover data are transformed by square root
cluster.single <- hclust (dis, method = 'single')
cluster.complete <- hclust (dis, 'complete')
cluster.average <- hclust (dis, 'average')
par (mfrow = c (1,3)) # will draw all dendrogram into one figure
plot (cluster.single, main = 'Single linkage')
plot (cluster.complete, main = 'Complete linkage')
plot (cluster.average, main = 'Average linkage')
class(df) # it's matrix type
df <- as.data.frame(df) # re-conversion from matrix to data frame
tree.fit <- rpart(Alk ~ ., data = df) # df must be in data frame class
class(df) # now it's data frame type
printcp(tree.fit) # display the results
plotcp(tree.fit) # visualize cross-validation results
summary(tree.fit) # detailed summary of splits
plot(tree.fit, uniform=TRUE,
     main="Tree classification of water samples: Po Plain, Italy")
text(tree.fit, use.n=TRUE, all=TRUE, cex=.8)

## using GAM pkg
mod_lm <- gam(Alk ~ ., data = df)
summary(mod_lm)

# Ref:
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/CA_lab.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/PCA_lab1.pdf
## http://bioconductor.wustl.edu/bioc/vignettes/pcaMethods
## https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
## https://www3.nd.edu/~mclark19/learn/GAMS.pdf





