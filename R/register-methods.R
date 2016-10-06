#' Generic method to register data to the epiviz data server
#' 
#' @param object The object to register to data server
#' @param columns Name of columns containing data to register
#' @param ... Additonal arguments passed to object constructors
#' @return An \code{\link{EpivizMetagenomicsData-class}} object 
#' @import metagenomeSeq
#' @importMethodsFrom epivizrData register
#' 
setMethod("register", "MRexperiment", function(object, columns=NULL, ...) {
  return(EpivizMetagenomicsData$new(object=object, columns=columns, ...))
})
