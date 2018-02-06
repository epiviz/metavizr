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

#' Method to replace NA or null feature labels with Not_Annotated_hierarchy-level
#'
#' @param replacing_na_obj_fData fData from MRexperiment object to replace NA or null
#' @param feature_order Order of features
#' @return hierarchy with NA or null feature labels replaced
#' @export

replaceNAFeatures = function(replacing_na_obj_fData, feature_order) {
  
  for(i in seq(1, length(feature_order))){
    na_indices <- which(is.na(replacing_na_obj_fData[,feature_order[i]]))
    for(j in seq(1, length(na_indices))){
      if(i > 1) {
        replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][na_indices[j]], sep="_")
      } else {
        replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
      }
    }
    na_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NA")
    for(j in seq(1, length(na_indices))){
      if(i > 1){ 
        replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][na_indices[j]], sep="_")
      } else{
        replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
      }
    }
    null_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NULL")
    for(j in seq(1, length(null_indices))){
      if(i > 1){ 
        replacing_na_obj_fData[,feature_order[i]][null_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][null_indices[j]], sep="_")
      } else{
        replacing_na_obj_fData[,feature_order[i]][null_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
      }
    }
  }
  
  replacing_na_obj_fData
  return(replacing_na_obj_fData)
}
