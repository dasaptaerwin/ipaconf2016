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
#install.packages("readxl") 
#install.packages("rpart")
#install.packages("party")
#install.packages("vegan")
#install.packages("mgcv")
#install.packages("gam")
#install.packages("psych")


## from Bioconductor server
#pcaMethods package from Bioconductor server
#source("http://bioconductor.org/biocLite.R") 
#biocLite("pcaMethods")

# Loading library
library(pcaMethods) # for pcaMethods package
library(cluster) # for cluster analysis
#library(readxl) # for opening data directly from xls format
library(rpart)
library(party)
library(vegan)
library(gam)
library(mgcv)
library(psych) # univar description

# Load data
df <- read.csv("dataPoplain.csv", header=T) # load data 
row.names(df) <- df$sites # setting row names
df <- df[4:17]
dim(df)
is.na(df) # checking NA if any
#df <- na.omit(df) # omitting NA if any

# Exploratory using pairs() function
# Assesing data patterns
cor.tab <- cor(df)
write.csv(cor.tab, "cortab.csv")
univar.tab <- describe(df)
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
plot(pca, type="lines")
slplot(pca) # default function in pcamethods but not big enough


## Plotting results (loadings and scores)
plot.new()
par(mfrow=c(1,2))
plot(loadings(pca), 
     pch = 20,
     main = "Variable loadings",
     sub = "Po Plain Italy")
text(loadings(pca), 
     row.names(loadings(pca)),
     cex=0.6, 
     pos=1, 
     col="blue")
abline(v=0, h=0, 
       col = "gray70")

plot(scores(pca), 
     main = "Case scores",
     sub = "Po Plain Italy")
abline(v=0, h=0, 
       col = "gray70")
text(scores(pca), row.names(df),
     cex=0.6, 
     pos=2, 
     col="blue", offset=0)
dev.off()

# PCA using'vegan' package
df <- scale(df)  #standardize the variables
vegan.pca <- rda(df)
vegan.pca
plot(vegan.pca)
biplot(vegan.pca, scaling = -1, display = 'species')
biplot(vegan.pca, scaling = -1, display = 'sites')

ordilabel(plot(vegan.pca), display="species", font=1, col="gray70") # Add some frame on the label
# reading "cleanplot.pca" function
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:cleanplot.pca?do=export_code&codeblock=0')
cleanplot.pca (vegan.pca)

# reading "evplot" function first:
source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=0")
# select the data frame with eigenvalues of particular axes:
ev <- vegan.pca$CA$eig

# calculate axis-importance and draw the barplots:
evplot(ev)
orditkplot(plot(vegan.pca)) # you see that you can move the labels.


# Cluster analysis (using 'vegan' package)
dis <- vegdist(df, method = 'bray') # percentage cover data are transformed by square root
cluster.single <- hclust (dis, method = 'single')
cluster.complete <- hclust (dis, 'complete')
cluster.average <- hclust (dis, 'average')
par (mfrow = c (1,3)) # will draw all dendrogram into one figure
plot (cluster.single, main = 'Single linkage')
plot (cluster.complete, main = 'Complete linkage')
plot (cluster.average, main = 'Average linkage')

# Regression tree
## using rpart pkh
tree.fit <- rpart(Alk ~ ., data = df) # error -> 'data' must be a data.frame
class(df) # identifying df class
df2 <- as.data.frame(df) # conversion from matrix to data frame
tree.fit <- rpart(Alk ~ ., data = df2) 
printcp(tree.fit) # display the results
plotcp(tree.fit) # visualize cross-validation results
summary(tree.fit) # detailed summary of splits
plot(tree.fit, uniform=TRUE,
     main="Tree classification of water samples: Po Plain, Italy")
text(tree.fit, use.n=TRUE, all=TRUE, cex=.8)

## using party pkg
install.packages("rpart.plot")
tree.fit <- ctree(Alk ~ ., data = df2) 
plot(tree.fit)

## using GAM pkg
mod_lm <- gam(Alk ~ ., data = df2)
summary(mod_lm)

# Ref:
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/CA_lab.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/PCA_lab1.pdf
## http://bioconductor.wustl.edu/bioc/vignettes/pcaMethods
## https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
## https://www3.nd.edu/~mclark19/learn/GAMS.pdf



