#' validate \code{\link[metagenomeSeq]{MRexperiment-class}} object
#' 
#' @param object an object of class \code{\link[metagenomeSeq]{MRexperiment-class}}
#' 
#' @return TRUE or FALSE
#' @export
validateObject <- function(object) {
  
  if(class(object) != "MRexperiment") {
    stop("Object is not an MRexperiment")
  }
  
  # validate feature data
  featureCheck <- TRUE
  fdata <- featureData(object)
  
  if(!validObject(fdata)) {
    featureCheck <- FALSE
    message("Not a valid Annotated Data frame - feature Data")
  }
  
  fdata <- pData(fdata)
  featureLength <- nrow(fdata)
  
  if (ncol(fdata) == 0) {
    featureCheck <- FALSE
    message("feature Data is empty")
  }
  
  # validate sample data
  sampleCheck <- TRUE
  sampleData <- phenoData(object)
  
  if(!validObject(sampleData)) {
    sampleCheck <- FALSE
    message("Not a valid Annotated Data frame - Sample/pheno Data")
  }
  
  sampleData <- pData(sampleData)
  sampleLength <- nrow(sampleData)
  
  if(ncol(sampleData) == 0) {
    sampleCheck <- FALSE
    message("Sample Data is empty")
  }
  
  # validate count/assay data
  assayCheck <- TRUE
  assayData <- assayData(object)
  
  if(!grepl("counts", names(assayData))) {
    assayCheck <- FALSE
    message("count data does not exist")
  }
  else {
    dimCountMatrix = dim(assayData[["counts"]])
    
    if(dimCountMatrix[1] != featureLength) {
      assayCheck <- FALSE
      message("count Matrix in assayData does not match feature length")
    }  
    
    if(dimCountMatrix[2] != sampleLength) {
      assayCheck <- FALSE
      message("count Matrix in assayData does not match sample length")
    }  
  }
  
  featureCheck && sampleCheck && assayCheck
}