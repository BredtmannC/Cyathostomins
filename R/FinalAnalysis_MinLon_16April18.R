## Let's consider a function giving Tina all that she wants

# source preprocessed data OR load library (if on slow local machine...)
source("FinalPreprocessing_MinLon_16April18.R")
## OR load libraries
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(reshape2)
library(ggplot2)

# Arguments:
# mergeWhitelists = F by default, = T if small group size (not enough discriminating data)
# mygroup = "species" by default, change to choose which groups you want to compare
# ID = "sample_ID" by default, but if we compare OTUs we need ID = ID_OTU e.g.
# folder = "1.species/" default folder for storing plots

# TinaSuperFunction <- function(mergeWhitelists = F,
#                               mygroup = "species",
#                               ID = "sample_ID",
#                               folder = "1.species/"){ 

mergeWhitelists = F
mygroup = "species"
ID = "sample_ID"
folder = "1.species/"
  


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
                         minFrequency = c(threshold), # we keep peaks present in at least 25% of the peaks within one group
                         labels = avgTina.info[[mygroup]], ## <-- choose which groups you want to compare ;)
                         mergeWhitelists = mergeWhitelists) ##if F filter criteria are applied groupwise, for smaller groups i would choose T
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
  visualizeThreshold <- function(steps = seq(0,1,.01)){
    
    # make a data frame with continuous threshold
    myDF <- sapply(steps, function(x) {getPeaksLength(getFilteredPeaks(x))})
    myDF <- data.frame(myDF)
    myDF <- cbind(sampleID = avgTina.info$sample_ID, myDF)
    myDF <- melt(myDF)
    
    myDF$variable <- factor(myDF$variable,
                            levels = levels(factor(myDF$variable)),
                            labels = steps)
    
    myplot <- ggplot(myDF, aes(x = variable, y = value)) +
      geom_line(aes(group = sampleID))+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(list(myDF, myplot))
  }
  
visualizeThreshold()
  
## To go from there


  peakPatterns(getFilteredPeaks(1), cex.axis = .8) 

  # Create featureMatrix and label the rows with the corresponding worms ID
  # featurematrix: gives intensity of peaks at m/z value for each avgspectrum
  featureMatrix <- intensityMatrix(peaks, avgSpectra)
  
  rownames(featureMatrix) <- avgTina.info[[mygroup]] ## <-- choose which groups you want to compare
  
  ##Clustering: Hierarchical clustering analyis with bootstrapping
  # AU (Approximately Unbiased) p-value and BP (Bootstrap Probability)
  library("pvclust")
  pv <- pvclust (t(featureMatrix),
                 method.hclust="ward.D2",
                 method.dist="euclidean") ## TODO check different methods of distances
  pv[["hclust"]][["labels"]] <- avgTina.info[[ID]]
  
  pdf(file =  paste0("./figures/", folder , "HierClust.pdf"))
  plot(pv, print.num=FALSE, main="Hierarchical Clustering")
  dev.off()
  
  
  ############## <- not checked after this point
  
  
  
  ## round mass data
  colnames(featureMatrix) <- sprintf("%.3f", as.double(colnames(featureMatrix)))

  ## Diagonal discriminant analysis
  library(sda)
  Xtrain <- featureMatrix
  Ytrain <- avgTina.info[[mygroup]] ## <-- choose which groups you want to compare
  ddar <- sda.ranking(Xtrain = featureMatrix, L = Ytrain, fdr = F, diagonal = T)

  pdf(file =  paste0("./figures/", folder , "Ddar.pdf"))
  plot(ddar)
  dev.off()
  
  # Linear discriminant Analysis
  ldar <- sda.ranking(Xtrain=featureMatrix, L=Ytrain, fdr=F,diagonal=F)
  
  pdf(file =  paste0("./figures/", folder , "Ldar.pdf"))
  plot(ldar)
  dev.off()
  
  #### Variable selection using Cross validation
  library(crossval)
  
  # Create a prediction function for the cross validation
  #to find out , how many and which peaks we need to discriminate
  
  predfun.dda <- function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=F){
    #estimate ranking and determine the best numVars variables
    ra <- sda.ranking(Xtrain, Ytrain, verbose=F,diagonal=diagonal,fdr=F)
    selVars <- ra[,'idx'][1:numVars]
    #fit and predict
    sda.out <- sda(Xtrain [,selVars, drop=F],Ytrain, diagonal = diagonal, verbose=F)
    ynew2 <- predict(sda.out,Xtest[,selVars, drop=F], verbose=F)$class
    #compute accuracy
    acc<- mean(Ytest==ynew2)
    return(acc)
  }
  
  # We want to repeat the cross-validation 20 times and use 5 folds
  K<-5 # number of folds
  B<-20 # number of repetitions
  
  # set.seed to reproductible analysis
  set.seed(1234)
  
  #for top 10 features (peaks): 
  cv.dda10 <- crossval(predfun.dda, X=featureMatrix, Y=avgTina.info[[mygroup]], 
                       K=K, B=B, numVars = 10, diagonal = F, verbose=F)
  cv.dda10$stat
  
  #look for optimal number of peaks (in the top 20)
  npeaks <- c(1:20, ncol(featureMatrix))
  
  #estimate accuracy for DDA (diagonal = T)
  set.seed(1234)
  cvsim.dda <- sapply(npeaks, function(i) {
    cv <- crossval(predfun.dda, 
                   X=featureMatrix, Y=avgTina.info[[mygroup]],
                   K=K, B=B, numVars=i, diagonal=T,
                   verbose=F)
    return(cv$stat)})

  #estimate accuracy for LDA (diagonal = F)
  set.seed(1234)
  cvsim.lda <- sapply(npeaks, function(i) {
    cv <- crossval(predfun.dda, 
                   X=featureMatrix, Y=avgTina.info[[mygroup]],
                   K=K, B=B, numVars=i, diagonal=F,
                   verbose=F)
    return(cv$stat)})
  
  # Combine the results and put them into a table
  result.sim <-cbind(npeaks=npeaks,
                     "DDA-ACC"=cvsim.dda,
                     "LDA-ACC"=cvsim.lda)
  
#  result.sim # shows table with top peaks and probability for discrimination

  # How to talk a bit (we like to chit chat)
  
  print("Hey Tina, you rule!!")
  
  ## Return a list
  TinaFullAnalysis <- list(pv = pv,
                           result.sim = result.sim)
  return(TinaFullAnalysis)
#}

resultsTina_species <- TinaSuperFunction()
resultsTina_species_sex <- TinaSuperFunction(mergeWhitelists = T, mygroup = "species_sex", folder = "2.species_sex/")
resultsTina_OTU <- TinaSuperFunction(mergeWhitelists = T, mygroup = "OTU", ID = "ID_OTU", folder = "3.OTU/")
resultsTina_OUT_sex <- TinaSuperFunction(mergeWhitelists = T, mygroup = "OTU_sex", ID = "ID_OTU", folder = "4.OTU_sex/")

resultsTina_species$result.sim
resultsTina_species_sex$result.sim
resultsTina_OTU$result.sim
resultsTina_OUT_sex$result.sim

plot(resultsTina_species$pv)
