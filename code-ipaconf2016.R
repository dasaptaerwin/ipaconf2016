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
install.packages("cluster") 
install.packages("rpart")
install.packages("vegan")
install.packages("gam")
install.packages("psych") # univar description
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")

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
library(psych) # univar description
library(ggplot2)
library(mgcv)

# Load data
df <- read.csv("dataPoplain.csv", header=T) # load data 
class(df)
dim(df)
is.na(df) # checking NA if any
#df <- na.omit(df) # omitting NA if any


# Exploratory using pairs() function
# Assesing data patterns
cor.tab <- cor(df[4:15])
cor.tab
write.csv(cor.tab, "cortab.csv")
univar.tab <- describe(df) # from psych package
univar.tab
write.csv(univar.tab, "univar.csv")

install.packages("PerformanceAnalytics")    
library(PerformanceAnalytics)
chart.Correlation(df[4:17], histogram=TRUE, pch=19) # visual PA

install.packages("ggcorrplot")
library(ggcorrplot)
correl <- round(cor(df[4:17]), 1)   # rounding correl matrix
p.mat <- cor_pmat(df[4:17])         # compute p-values
ggcorrplot(correl)              # making heatmap

# Run PCA 
installed.packages("FactoMineR")
installed.packages("factoextra")
library("FactoMineR")
library("factoextra")
res.pca <- PCA(df[4:17], graph = FALSE)
fviz_pca(res.pca, choix = "var", habillage=df$type) + labs(title="Biplot PCA Po Plain")

# Run Cluster 
distance <- dist(scale(df[4:17]), method = "euclidean")
cluster <- hclust(distance, method = "complete")
plot(cluster, cex = 0.6, hang = -1, main = "CA Po Plain") 
rect.hclust(cluster, k = 3, border = 2:4) 


################################################################
################ WE DON'T USE THIS PART ########################
################################################################

# Run PCA (using pcamethods package)
## svdImpute = standard pca, with imputation, standardised, method univariate (uv)
pca <- pca(df[4:15], 
           scale = "uv",
           center = T,
           nPcs = 5,
           evalPcs = 1:5)
pca
write.csv((summary(pca)), "pca.csv")


## Evaluating results
loadings(pca) # loadings of each variables
write.csv(loadings(pca), "loadings.csv")
plotPcs(pca, Pcs = 1:5)
plot(pca)
slplot(pca) # default function in pcamethods but not big enough


# PCA using'vegan' package
vegan.pca <- rda(df[4:17], scale=TRUE)
vegan.pca
plot(vegan.pca)
biplot(vegan.pca, scaling = -1, display = 'species')
biplot(vegan.pca, scaling = -1, display = 'sites')
plot(vegan.pca, labels=rownames(df$area))

dev.off()

ordilabel(plot(vegan.pca), display="species", font=1, col="gray70") # Add some frame on the label
orditkplot(plot(vegan.pca, display = 'sites')) # you see that you can move the labels.

# reading "cleanplot.pca" function to make better visualisation
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:cleanplot.pca?do=export_code&codeblock=0')
cleanplot.pca(vegan.pca)

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

# Regression tree using rpart pkg
# rpart() must be in data frame class
class(df) # it's matrix class
df <- as.data.frame(df) # re-conversion from matrix to data frame
class(df) # now it's data frame type
tree.fit <- rpart(Alk ~ ., data = df) # df must be in data frame class
printcp(tree.fit) # display the results
plotcp(tree.fit) # visualize cross-validation results
summary(tree.fit) # detailed summary of splits
plot(tree.fit, uniform=TRUE,
     main="Tree classification of water samples: Po Plain, Italy")
text(tree.fit, use.n=TRUE, all=TRUE, cex=.8)

## using GAM pkg
mod_lm <- gam(Alk ~ Ca + Na + Mg + 
                K + Cl + SO + H + 
                O + C + Sr + Br + 
                I + Fe, data = df)
summary(mod_lm)
gam.check(mod_lm, type=c("deviance","pearson","response"))

mod_lm1 <- gam(Cl ~ Alk + Ca + Na + Mg + 
                 K + SO + H + 
                 O + C + Sr + Br + 
                 I + Fe, data = df)
summary(mod_lm1)
gam.check(mod_lm1, type=c("deviance","pearson","response"))

AIC(mod_lm)
AIC(mod_lm1)



# Regression tree using rpart pkg
# rpart() must be in data frame class
class(df) # it's matrix class
df <- as.data.frame(df) # re-conversion from matrix to data frame
class(df) # now it's data frame type
tree.fit <- rpart(Alk ~ ., data = df) # df must be in data frame class
printcp(tree.fit) # display the results
plotcp(tree.fit) # visualize cross-validation results
summary(tree.fit) # detailed summary of splits
plot(tree.fit, uniform=TRUE,
     main="Tree classification of water samples: Po Plain, Italy")
text(tree.fit, use.n=TRUE, all=TRUE, cex=.8)

tree.fit2 <- rpart(Cl ~ ., data = df) # df must be in data frame class
printcp(tree.fit2) # display the results
plotcp(tree.fit2) # visualize cross-validation results
summary(tree.fit2) # detailed summary of splits
plot(tree.fit2, uniform=TRUE,
     main="Tree classification of water samples: Po Plain, Italy")
text(tree.fit2, use.n=TRUE, all=TRUE, cex=.8)
par(mfrow)
dev.off()


# Ref:
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/CA_lab.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/PCA_lab1.pdf
## http://bioconductor.wustl.edu/bioc/vignettes/pcaMethods
## https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
## https://www3.nd.edu/~mclark19/learn/GAMS.pdf





