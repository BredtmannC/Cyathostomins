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

TinaSuperFunction <- function(mergeWhitelists = F,
                              mygroup = "species",
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
  pdf(file =  paste0("./figures/5.General/", "threshold.pdf"))
  
  ggplot(myDF, aes(x = variable, y = value)) +
    geom_line(aes(group = sampleID))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept=26, col = "red")
   
  dev.off()

  # OUTPUT FIGURE 1
  # Peaks that are shared by 100% (per group) of data --> bio markers
  pdf(file =  paste0("./figures/5.General/" , "peakPatterns100.pdf"))
  peakPatterns(getFilteredPeaks(1), cex.axis = .8)
  dev.off()
  
  # OUTPUT FIGURE 2
  # Peaks that are shared by 25% (per group) of data
  pdf(file =  paste0("./figures/5.General/" , "peakPatterns25.pdf"))
  peakPatterns(getFilteredPeaks(.25), cex.axis = .8)
  dev.off()

  ############## -> we take 25% as a threshold
  filteredpeaks <- getFilteredPeaks(.25)
  
  # Create featureMatrix and label the rows with the corresponding worms ID
  # featurematrix: gives intensity of filteredpeaks at m/z value for each avgspectrum
  featureMatrix <- intensityMatrix(filteredpeaks, avgSpectra)
  
  rownames(featureMatrix) <- avgTina.info[[mygroup]] ## <-- choose which groups you want to compare
  
  ## TODOOOOO
  
  ##Clustering: Hierarchical clustering analyis with bootstrapping
  # AU (Approximately Unbiased) p-value and BP (Bootstrap Probability)
  library("pvclust")
  pv <- pvclust (t(featureMatrix),
                 method.hclust="ward.D2",
                 method.dist="euclidean") ## TODO check different methods of distances
  # Ward's method (and centroid, and so called "median" methods) 
  # are involved in computing geometrical centroids in euclidean space. 
  # They do it in a way that requires squared euclidean distances (?)
  # https://stats.stackexchange.com/questions/109597/how-to-choose-the-right-distance-matrix-for-clustering
  
  
  pv[["hclust"]][["labels"]] <- avgTina.info[[ID]]
  
  # OUTPUT FIGURE 3
  # Hierarchical clustering
  pdf(file =  paste0("./figures/", folder , "HierClust.pdf"))
  plot(pv, print.num=FALSE, main="Hierarchical Clustering")
  dev.off()
  
  ## round mass data
  colnames(featureMatrix) <- sprintf("%.3f", as.double(colnames(featureMatrix)))
  
  ## Diagonal discriminant analysis
  library(sda)
  Xtrain <- featureMatrix
  Ytrain <- avgTina.info[[mygroup]] ## <-- choose which groups you want to compare
  ddar <- sda.ranking(Xtrain = featureMatrix, L = Ytrain, fdr = F, diagonal = T)
  
  # OUTPUT FIGURE 4
  pdf(file =  paste0("./figures/", folder , "Ddar.pdf"))
  FIG4 <- plot(ddar) # TODO : INCREASE THE SIZE OF LABELS!!!
  dev.off()
  
  # Linear discriminant Analysis
  ldar <- sda.ranking(Xtrain=featureMatrix, L=Ytrain, fdr=F,diagonal=F)
  
  # OUTPUT FIGURE 5
  pdf(file =  paste0("./figures/", folder , "Ldar.pdf"))
  FIG5 <- plot(ldar)
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
  
  #for top 40 features (peaks): 
  cv.dda40 <- crossval(predfun.dda, X=featureMatrix, Y=avgTina.info[[mygroup]], 
                       K=K, B=B, numVars = 40, diagonal = F, verbose=F)
  cv.dda40$stat
  
  #look for optimal number of peaks (in the top 40)
  npeaks <- c(1:40, ncol(featureMatrix))
  
  # estimate accuracy for DDA (diagonal = T)
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
  result.sim <- data.frame(cbind(npeaks=npeaks,
                                 "DDA-ACC"=cvsim.dda,
                                 "LDA-ACC"=cvsim.lda))
  
  result.sim # shows table with top peaks and probability for discrimination
  
  # Plot linear discriminant Analysis with accuracy values
  
  ## TODO FIX IT
  
  # labels.lda <- paste0(colnames(featureMatrix)[1:40], " accuracy ", 
  #                      as.character(round(result.sim$LDA.ACC[-41] * 100, 1)), "%")
  # length(labels.lda)
  # 
  # colnames(featureMatrix)
  # plot(ldar) #, ylab = 1:418)
  # 
  # colnames(featureMatrix) <- as.character(1:417)
  # 
  # # dev.off()
  # 
  # ?plot.sda.ranking
  # 
  # How to talk a bit (we like to chit chat)
  
  print("Hey Tina, you rule!!")
  
  ## Return a list
  TinaFullAnalysis <- list(pv = pv,
                           result.sim = result.sim,
                           visualizeThreshold = visualizeThreshold(),
                           fig1 = FIG1)
  
  return(TinaFullAnalysis)
}
 
resultsTina_species <- TinaSuperFunction()

resultsTina_species_sex <- TinaSuperFunction(mygroup = "species_sex", folder = "2.species_sex/")
resultsTina_OTU <- TinaSuperFunction(mygroup = "OTU", ID = "ID_OTU", folder = "3.OTU/")
resultsTina_OUT_sex <- TinaSuperFunction(mygroup = "OTU_sex", ID = "ID_OTU", folder = "4.OTU_sex/")

resultsTina_species$result.sim
resultsTina_species_sex$result.sim
resultsTina_OTU$result.sim
resultsTina_OUT_sex$result.sim

plot(resultsTina_species$pv)
