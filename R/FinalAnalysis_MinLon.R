## Authors: 
# Christina Bredtmann, 
# Alice Balard

library(reshape2)
library(ggplot2)

# source preprocessed data 
source("FinalPreprocessing_MinLon.R")

# Arguments:
# mygroup = "species" by default, change to choose which groups you want to compare
# ID = "sample_ID" by default, but if we compare OTUs we need ID = ID_OTU e.g.
# folder = "1.species/" default folder for storing plots

SuperCyaFun <- function(mygroup = "species",
                              ID = "sample_ID",
                              folder = "1.species/"){
  
  #post processing: Only comparison MIN and LON
  #load data
  load("./avgSpectra.RData")
  load("./avgTina.info.RData")
  load("./peaks.RData")
  
  # Bin the peaks (make similar peak mass values identical)
  peaks <- binPeaks(peaks)
  
  # save original peaks, before filtering
  peaks_beforeFiltering <- peaks
  
  # get filtered peaks depending on threshold (remove the less frequent ones)
  getFilteredPeaks <- function(threshold){
    peaks <- filterPeaks(peaks_beforeFiltering, 
                         minFrequency = c(threshold), # we keep peaks present in at least e.g. 25% of the peaks within one group
                         labels = avgTina.info[[mygroup]], ## choose which groups you want to compare
                         mergeWhitelists = F) ## filter criteria are applied groupwise
    return(peaks)
  }
  
  # get peaks lenght vector depending on peaks
  getPeaksLength <- function(x){
    myPeaksLength <- 0
    for (i in 1:length(snrPeaks(x))){
      myPeaksLength <- c(myPeaksLength, length(snrPeaks(x)[[i]]))
    }
    return(myPeaksLength[-1])
  }
  
  # Decide on a threshold
  steps = seq(0,1,.01)
  
  # make a data frame with continuous threshold
  myDF <- sapply(steps, function(x) {getPeaksLength(getFilteredPeaks(x))})
  myDF <- data.frame(myDF)
  myDF <- cbind(sampleID = avgTina.info$sample_ID, myDF)
  myDF <- melt(myDF)
  
  myDF$variable <- factor(myDF$variable,
                          levels = levels(factor(myDF$variable)),
                          labels = steps)
  # OUTPUT FIGURE 0
  pdf(file =  paste0("./figures/" , folder , "threshold.pdf"))
  
  ggplot(myDF, aes(x = variable, y = value)) +
    geom_line(aes(group = sampleID))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept=26, col = "red")
   
  dev.off()

  # OUTPUT FIGURE 1
  # Peaks that are shared by 100% (per group) of data --> bio markers
  pdf(file =  paste0("./figures/" , folder , "peakPatterns100.pdf")) 
  peakPatterns(getFilteredPeaks(1), cex.axis = .8)
  dev.off()
  
  # OUTPUT FIGURE 2
  # Peaks that are shared by 25% (per group) of data
  pdf(file =  paste0("./figures/" , folder , "peakPatterns25.pdf")) 
  peakPatterns(getFilteredPeaks(.25), cex.axis = .8)
  dev.off()

  ############## -> we take 25% as a threshold to keep as much information as possible
  filteredpeaks <- getFilteredPeaks(.25)
  
  # Create featureMatrix and label the rows with the corresponding worms ID
  # featurematrix: gives intensity of filteredpeaks at m/z value for each avgspectrum
  featureMatrix <- intensityMatrix(filteredpeaks, avgSpectra)
  
  rownames(featureMatrix) <- avgTina.info[[mygroup]] ## <-- choose which groups you want to compare
  
  ##Clustering: Hierarchical clustering analyis with bootstrapping
  # AU (Approximately Unbiased) p-value and BP (Bootstrap Probability)
  library("pvclust")
  pv <- pvclust (t(featureMatrix),
                 method.hclust="ward.D2",
                 method.dist="euclidean") 
  
  pv[["hclust"]][["labels"]] <- avgTina.info[[ID]]
  
  # OUTPUT FIGURE 3
  # Hierarchical clustering
  pdf(file =  paste0("./figures/", folder , "HierClust.pdf"))
  plot(pv, print.num=FALSE, main="Hierarchical Clustering")
  dev.off()
  
  ## round mass data
  colnames(featureMatrix) <- sprintf("%.3f", as.double(colnames(featureMatrix)))
  
  library(sda)
  Xtrain <- featureMatrix
  Ytrain <- avgTina.info[[mygroup]] ## <-- choose which groups you want to compare

  # Linear discriminant Analysis
  ldar <- sda.ranking(Xtrain=featureMatrix, L=Ytrain, fdr=F, diagonal=F)
  
  # OUTPUT FIGURE 4
  pdf(file =  paste0("./figures/", folder , "Ldar.pdf"))
  FIG4 <- plot(ldar)
  dev.off()
  
  #### Variable selection using Cross validation
  library(crossval)
  
  # Create a prediction function for the cross validation
  # to find out how many and which peaks we need to discriminate
  
  predfun <- function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=F){
    
    #estimate ranking and determine the best numVars variables
    ra <- sda.ranking(Xtrain, Ytrain, verbose=F,diagonal=diagonal,fdr=F)
    selVars <- ra[,'idx'][1:numVars]
    
    #fit and predict
    sda.out <- sda(Xtrain [,selVars, drop=F],Ytrain, diagonal = diagonal, verbose=F)
    ynew2 <- predict(sda.out, Xtest[,selVars, drop=F], verbose=F)$class
    
    #compute accuracy
    acc <- mean(Ytest==ynew2)

    return(acc)
  }
  
  # We want to repeat the cross-validation 20 times and use 5 folds
  K<-5 # number of folds
  B<-20 # number of repetitions
  
  # set.seed to reproductible analysis
  set.seed(1234)
  
  #for top N features (peaks): 
  cv.lda.fun <- function(numVars){
    CV <- crossval(predfun, 
                   X=featureMatrix, 
                   Y=avgTina.info[[mygroup]], 
                   K=K, B=B, 
                   numVars = numVars, 
                   diagonal = T, 
                   verbose=F)
    print(paste0("The accuracy represents the performance of the top ", 
                 numVars,
                 " markers (how well they predict the groups)"))
    print(CV$stat)
  }
  
  #for top 1 features (peaks): 
  cv.lda1 <- cv.lda.fun(1)
  #for top 10 features (peaks): 
  cv.lda10 <- cv.lda.fun(10)
  #for top 20  features (peaks): 
  cv.lda20 <- cv.lda.fun(20)
  #for top 30 features (peaks): 
  cv.lda30 <- cv.lda.fun(30)
  #for top 40 features (peaks): 
  cv.lda40 <- cv.lda.fun(40)
  
    #look for optimal number of peaks (in the top 40)
  npeaks <- c(1:40, ncol(featureMatrix))
  
  #estimate accuracy for LDA (diagonal = F)
  set.seed(1234)
  cvsim.lda <- sapply(npeaks, function(i) {
    cv <- crossval(predfun, 
                   X=featureMatrix, 
                   Y=avgTina.info[[mygroup]],
                   K=K, B=B, 
                   numVars=i, 
                   diagonal=F,
                   verbose=F)
    return(cv$stat)})
  
  # Put the results into a table
  result.sim <- data.frame(cbind(npeaks=npeaks,
                                 "LDA-ACC"=cvsim.lda))
  
  result.sim # shows table with top peaks and probability for discrimination
  
  # Plot linear discriminant Analysis with accuracy values
  
  ## Return a list
  TinaFullAnalysis <- list(pv = pv,
                           result.sim = result.sim)
  
  return(TinaFullAnalysis)
}
 
resultsTina_species <- SuperCyaFun()

ggplot(resultsTina_species$result.sim, aes(x = npeaks, y = LDA.ACC)) +
  geom_point() +
  geom_line(col = "red") +
  # geom_smooth(col = "red", se = F) +
  geom_vline(xintercept = 40) +
  theme_bw()

# Return results according to the different categories
resultsTina_species_sex <- SuperCyaFun(mygroup = "species_sex", folder = "2.species_sex/")
resultsTina_OTU <- SuperCyaFun(mygroup = "OTU", ID = "ID_OTU", folder = "3.OTU/")
resultsTina_OTU_sex <- SuperCyaFun(mygroup = "OTU_sex", ID = "ID_OTU", folder = "4.OTU_sex/")

resultsTina_species$result.sim
resultsTina_species_sex$result.sim
resultsTina_OTU$result.sim
resultsTina_OTU_sex$result.sim

plot(resultsTina_OTU$pv, print.num=F)