simulateData <- function(sampleNum, replicateNum, batchImbalance, geneNum, 
                         setTheSeed = NULL, groupOneMeanSD = c(0,1.5), groupTwoSD = 1, 
                         sampledValueSD = 1, noiseLevel = 0.25,
                         batchEffectLocation = c(1,2), batchEffectScale = c(6,4,4,5)) {

  ####Define suppot variables
  sampleNum.5 <- sampleNum/2
  sigd <- nchar(sampleNum + 2*replicateNum)
  
  ####Load library
  library(MCMCpack)
  library(stringr)
  
  #####Set the Seed
  set.seed(setTheSeed)
  
  ####Generate Groups
  #Group 1
  s1geneMean = rnorm(geneNum, groupOneMeanSD[1], groupOneMeanSD[2]) #Generate mean expression of each gene
  mat1 = sapply(s1geneMean, function(x) rnorm(sampleNum.5, x, sampledValueSD) ) #Gene expression values
  colnames(mat1) = paste("Gene", str_pad(as.character(c(1:geneNum)),nchar(geneNum),"left","0"), sep="") #Name
  rownames(mat1) = paste("Sample", str_pad(as.character(c(1:(sampleNum.5))),sigd,"left","0"), sep="") #Name
  #Group 2
  s2geneMean = sapply(s1geneMean, function(x)  rnorm(1, 1, groupTwoSD)*x) #Generate mean expression of each gene
  mat2 = sapply(s2geneMean, function(x) rnorm(sampleNum.5, x, sampledValueSD) ) #Gene expression values
  colnames(mat2) = paste("Gene", str_pad(as.character(c(1:geneNum)),nchar(geneNum),"left","0"), sep="") #Name
  rownames(mat2) = paste("Sample", str_pad(as.character(c((sampleNum.5+1):sampleNum)),sigd,"left","0"), sep="") #Name
  
  #Bind
  oriData = rbind(mat1, mat2)
  
  oriData = data.frame(group=c(rep("group1", sampleNum.5), rep("group2", sampleNum.5)) , oriData)
  oriData = data.frame(rownames(oriData), oriData)
  colnames(oriData)[1] = "Sample.ID"
  
  ####Define replicates
  repvec <- c((sampleNum.5 - replicateNum + 1):sampleNum.5,
              (sampleNum   - replicateNum + 1):sampleNum) #Index of last replicateNum samples in each group
  
  #Grab replicates
  groupRep = oriData[repvec, ] #Take replicates
  groupRep = rbind(groupRep, groupRep) #Double them
  
  groupRep = data.frame(batch = c(rep("batch1", 2*replicateNum), rep("batch2", 2*replicateNum)) , groupRep) #Add batch label
  groupRep = data.frame(rep = "Yes", groupRep) #Add Replicate label
  
  #Remove replicates from remainder of data
  oriData = oriData[(1:sampleNum)[-repvec], ]
  
  if (nrow(oriData) > 0) {
    ####Implement imbalance
    if (batchImbalance > 0) {
      b1vec <- c(1:batchImbalance,(sampleNum.5 -replicateNum + batchImbalance +1):(sampleNum - 2*replicateNum)) #Samples in batch 1
      vec <- 1:(sampleNum - 2*replicateNum)
      b2vec <- vec[!(vec %in% b1vec)] #Samples in batch 2
      
      oriData = oriData[c(b1vec,b2vec), ] #Reorder data
    }
    oriData = data.frame(batch = c(rep("batch1", sampleNum.5-replicateNum), rep("batch2",sampleNum.5-replicateNum)) , oriData) #Add batch label
    oriData = data.frame(rep = "No", oriData) #Add replicate label
    
    oriData = rbind (oriData, groupRep) #Combine replicate and non replicates
    oriData = oriData[order(oriData[ ,"batch"]), ] #Reorder by batch
  } else {
    oriData = groupRep
  }
  
  ####Add noise level
  merrMat <- matrix(rnorm(nrow(oriData)*geneNum,1,noiseLevel),nrow = nrow(oriData))
  oriData[,5:ncol(oriData)] <- oriData[,5:ncol(oriData)]*merrMat
  
  ####Implement batch effects
  batchData = oriData #Copy data
  
  batchMean = rnorm(geneNum, batchEffectLocation[1], batchEffectLocation[2] )  #Location batch effect values
  batchScale1 = rinvgamma(geneNum, batchEffectScale[1], batchEffectScale[2]) #Scale batch effects for batch ones
  batchScale2 = rinvgamma(geneNum, batchEffectScale[3], batchEffectScale[4]) #Scale batch effects for batch two
  
  y <- batchData[batchData[,'batch'] == 'batch1', 5:(ncol(batchData))]
  for (i in 1:ncol(y)) {tmp <- mean(y[,i]);y[,i] <- (tmp + batchMean[i]) + ((y[,i]-tmp) * batchScale1[i])}
  batchData[batchData[,'batch'] == 'batch1', 5:(ncol(batchData))] = y
  
  y <- batchData[batchData[,'batch'] == 'batch2', 5:(ncol(batchData))]
  for (i in 1:ncol(y)) {tmp <- mean(y[,i]);y[,i] <- (tmp - batchMean[i]) + ((y[,i]-tmp) * batchScale2[i])}
  batchData[batchData[,'batch'] == 'batch2', 5:(ncol(batchData))] = y
  
  #Separate and save
  batchData1 = batchData[batchData[,'batch'] == 'batch1', ]
  rownames(batchData1) <- as.character(batchData1[,'Sample.ID'])
  batchData2 = batchData[batchData[,'batch'] == 'batch2', ]
  rownames(batchData2) <- as.character(batchData2[,'Sample.ID'])
  
  #Sort
  batchData1 <- batchData1[order(batchData1[,'Sample.ID']),]
  batchData2 <- batchData2[order(batchData2[,'Sample.ID']),]
  
  #Return
  return(list(batch1 = batchData1,batch2 = batchData2))
}