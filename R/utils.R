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
#' @examples
#' generateSelection("Bacteroidales", 1L, 2L)
#' 
generateSelection = function(feature_names, aggregation_level, selection_type, feature_order=NULL){
  fSelection <- list()
  fSelection$featureNames <- feature_names
  fSelection$featureOrder <- feature_order
  fSelection$featureLevel <- aggregation_level
  fSelection$selectionType <- selection_type
  return(fSelection)
}
