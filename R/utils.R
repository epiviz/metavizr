.generatePseudoGUID = function(size) {
  chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
  ret = c()
  indices = sample(1:nchar(chars), size, replace=TRUE)

  for (i in 1:size) {
    ret = c(ret, substr(chars, indices[i], indices[i]))
  }

  return(paste(ret, collapse=""))
}

SelectionType = list(NONE=0, LEAVES=1, NODE=2)


#' Method to select and set aggregation type to nodes in FacetZoom
#'
#' @param feature_names Selected Features
#' @param aggregation_level Level in the hierarchy
#' @param selection_type Expanded, aggregated, or removed
#' @param feature_order Order of features at that level
#' @return A selection object for a metavizControl object to accept
#' @export
#'
generateSelection = function(feature_names, aggregation_level, selection_type, feature_order=NULL){
  fSelection <- list()
  fSelection$featureNames <- feature_names
  fSelection$featureOrder <- feature_order
  fSelection$featureLevel <- aggregation_level
  fSelection$selectionType <- selection_type
  return(fSelection)
}

timeSeriesToMRexperiment<-function(obj,sampleNames=NULL,
                          sampleDescription="timepoints",
                          taxonomyLevels=NULL,
                          taxonomyHierarchyRoot="bacteria",
                          taxonomyDescription="taxonomy",
                          featuresOfInterest = NULL,
                          feature_data=NULL){
  if(is.null(obj)){
    stop("Matrix cannot be null")
  }
  if(is.null(sampleNames)){
    numSamples <- dim(obj[[1]]$fit)[1]
    sampleNames <- paste("Timepoint", 1:numSamples, sep="_")
  }
  
  if(is.null(featuresOfInterest)){
    hasFit <- lapply(1:(length(obj)-1), function(i) which(!is.null(obj[[i]]$fit)))
    featuresOfInterest <- which(hasFit == 1)
    hasFit <- (hasFit == 1)
    hasFit <- !is.na(hasFit)
    temp <- 1:length(hasFit)
    temp[!hasFit] <- 0
    hasFit <- temp
  }
  
  if(is.null(taxonomyLevels)){
    #numLevels <- length(featuresOfInterest)
    #taxonomyLevels <- names(obj)[featuresOfInterest]
    
    numLevels <- 1:length(hasFit)
    taxonomyLevels <- names(obj)[1:length(hasFit)]
  }
  
  numSamples <- length(sampleNames)
  numLevels <- length(taxonomyLevels)
  numFeaturesOfInterest <- length(featuresOfInterest)
  
  rangeSamples <- 1:numSamples
  rangeFeaturesOfInterest <- 1:numFeaturesOfInterest
  print(hasFit)
  
  results <- do.call(rbind, lapply(hasFit,function(i){ if (i != 0) t(obj[[i]]$fit)[1,] else rep(NA, numSamples) }))
  
  dfSamples <- data.frame(x=rangeSamples,row.names=sampleNames)
  metaDataSamples <-data.frame(labelDescription=sampleDescription)
  annotatedDFSamples <- AnnotatedDataFrame()
  pData(annotatedDFSamples) <- dfSamples
  varMetadata(annotatedDFSamples) <- metaDataSamples
  validObject(annotatedDFSamples)
  
  if(is.null(feature_data)){
    dfFeatures <- data.frame(taxonomy1=rep(taxonomyHierarchyRoot, numLevels),taxonomy2=taxonomyLevels)
    metaDataFeatures <-data.frame(labelDescription=paste(taxonomyDescription, 1:2, sep=""))
    annotatedDFFeatures <- AnnotatedDataFrame()
    pData(annotatedDFFeatures) <- dfFeatures
    varMetadata(annotatedDFFeatures) <- metaDataFeatures
    validObject(annotatedDFFeatures)
  }
  else{
    annotatedDFFeatures <- feature_data
  }
  
  fitTimeSeriesMRexp <- newMRexperiment(counts=results,
                                        phenoData=annotatedDFSamples,
                                        featureData=annotatedDFFeatures)
  return(fitTimeSeriesMRexp)
}