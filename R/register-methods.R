#' Generic method to register data to the epiviz data server
#' 
#' @param object The object to register to data server
#' @param ... Additonal arguments passed to object constructors
#' @return Object inheriting from \code{\link{EpivizData}} class
#' @export


#' @describeIn register Register a \code{\link{MRexperiment}} object
#' @import metagenomeSeq
#' @param object Which object to register
#' 
setMethod("register", "MRexperiment", function(object, ...) {
  return(EpivizMetagenomicsData$new(object=object, ...))
})
