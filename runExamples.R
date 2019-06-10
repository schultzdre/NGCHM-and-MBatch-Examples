#Clear
rm(list=ls())

#Load libraries
library(NGCHM)
library(MBatch)
library(stringr)

#Load data generation function
source("simulateData.R")

#Simulate data
set.seed(713)
dat <- simulateData(500,20,100,400)

#Initialize variables
groupvec <- c()
repvec <- c()
batchvec <- c()

#Tailor each batch

for (i in 1:2) {
  #Tailor names
  rownames(dat[[i]]) <- paste0(rownames(dat[[i]]),".",dat[[i]][,"batch"])
  
  #Get group data
  vec <- as.character(dat[[i]][,"group"])
  names(vec) <- rownames(dat[[i]])
  groupvec <- c(groupvec,vec)
  
  #Get replicate data
  vec <- as.character(dat[[i]][,"rep"])
  names(vec) <- rownames(dat[[i]])
  repvec <- c(repvec,vec)
  
  #Get batch data
  vec <- as.character(dat[[i]][,"batch"])
  names(vec) <- rownames(dat[[i]])
  batchvec <- c(batchvec,vec)
}

#Make covariate bars
#Group
colcmap <- chmNewColorMap(c("group1","group2"),colors = c("lawngreen","tomato"),missing.color="white")
groupcovar <- chmNewCovariate('Group',values=groupvec,
                            value.properties=colcmap,type="discrete")
#Replicate
colcmap <- chmNewColorMap(c("No","Yes"),colors = c("white","black"),missing.color="white")
repcovar <- chmNewCovariate('Replicate',values=repvec,
                              value.properties=colcmap,type="discrete")
#Batch
colcmap <- chmNewColorMap(c("batch1","batch2"),colors = c("blue4","orange4"),missing.color="white")
batchcovar <- chmNewCovariate('Batch',values=batchvec,
                            value.properties=colcmap,type="discrete")

#Combine data
mat <- cbind(t(dat[[1]][,-1:-4]),t(dat[[2]][,-1:-4]))

#Cluster
rowclust <- hclust(as.dist(1-cor(t(mat),use="pairwise.complete.obs")),method="ward.D2")
colclust <- hclust(as.dist(1-cor(mat,use="pairwise.complete.obs")),method="ward.D2")

#Median center
mcmat <- t(apply(mat,1,function(x){x-median(x,na.rm = TRUE)}))

#NGCHM
rwbmap <- chmNewColorMap(c(-8,0,8),colors = c('blue','white','red'),
                         missing.color = "gray70")
layer1 <- chmNewDataLayer ('ColorMap',mcmat,rwbmap)
chm <- chmNew("uncorected",layer1,
              rowOrder = as.dendrogram(rowclust),
              colOrder = as.dendrogram(colclust))
chm <- chmAddCovariateBar(chm,'column',batchcovar)
chm <- chmAddCovariateBar(chm,'column',groupcovar)
chm <- chmAddCovariateBar(chm,'column',repcovar)

#Export
chmExportToFile(chm, "uncorrected.ngchm", overwrite = TRUE)

#Remove batch from column names
for (i in 1:2) {rownames(dat[[i]]) <- str_replace(rownames(dat[[i]]),"\\.batch[1-2]$","")}

#Correct the data using EBN
y <- EBNPlus_Correction_Structures(t(as.matrix(dat[[1]][,-1:-4])),
                                   t(as.matrix(dat[[2]][,-1:-4])),"batch1","batch2",
                                   theEBNP_BatchWithZero="both",
                                   theEBNP_FixDataSet=as.numeric(NA),
                                   theEBNP_CorrectForZero=FALSE,
                                   theEBNP_ParametricPriorsFlag=TRUE)

#Cluster
rowclust <- hclust(as.dist(1-cor(t(y),use="pairwise.complete.obs")),method="ward.D2")
colclust <- hclust(as.dist(1-cor(y,use="pairwise.complete.obs")),method="ward.D2")

#Median center
mcmat <- t(apply(y,1,function(x){x-median(x,na.rm = TRUE)}))

#NGCHM
rwbmap <- chmNewColorMap(c(-5,0,5),colors = c('blue','white','red'),
                         missing.color = "gray70")
layer1 <- chmNewDataLayer ('ColorMap',mcmat,rwbmap)
chm <- chmNew("EBN_corrected",layer1,
              rowOrder = as.dendrogram(rowclust),
              colOrder = as.dendrogram(colclust))
chm <- chmAddCovariateBar(chm,'column',batchcovar)
chm <- chmAddCovariateBar(chm,'column',groupcovar)
chm <- chmAddCovariateBar(chm,'column',repcovar)

#Export
chmExportToFile(chm, "EBN_corrected.ngchm", overwrite = TRUE)

#Correct Data
y <- RBN_Replicates(t(as.matrix(dat[[1]][,-1:-4])),
                    t(as.matrix(dat[[2]][,-1:-4])),"batch1","batch2")
colnames(y) <- str_replace(colnames(y),"-",".")

#Cluster
rowclust <- hclust(as.dist(1-cor(t(y),use="pairwise.complete.obs")),method="ward.D2")
colclust <- hclust(as.dist(1-cor(y,use="pairwise.complete.obs")),method="ward.D2")

#Median center
mcmat <- t(apply(y,1,function(x){x-median(x,na.rm = TRUE)}))

#NGCHM
rwbmap <- chmNewColorMap(c(-5,0,5),colors = c('blue','white','red'),
                         missing.color = "gray70")
layer1 <- chmNewDataLayer ('ColorMap',mcmat,rwbmap)
chm <- chmNew("RBN_corrected",layer1,
              rowOrder = as.dendrogram(rowclust),
              colOrder = as.dendrogram(colclust))
chm <- chmAddCovariateBar(chm,'column',batchcovar)
chm <- chmAddCovariateBar(chm,'column',groupcovar)
chm <- chmAddCovariateBar(chm,'column',repcovar)

#Export
chmExportToFile(chm, "RBN_corrected.ngchm", overwrite = TRUE)
