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
#' 
#' @examples
#' library(metagenomeSeq)
#' data(mouseData)
#' feature_order <- colnames(fData(mouseData))
#' replaceNAFeatures(fData(mouseData),feature_order)

replaceNAFeatures = function(replacing_na_obj_fData, feature_order) {
  
  message("replacing feature labels")

  tax <- replacing_na_obj_fData
  feature_order <- colnames(tax)
  for(i in seq(2, length(feature_order))) {
    sub_table <- tax[, feature_order[1:i]]

    lineages <- apply(sub_table, 1, paste, collapse="_")
    lineages <- as.data.frame(lineages)
    lineages$group <- sub_table[, feature_order[i]]

    # grp <- sub_table[, by = feature_order[i]]
    # grouped <- aggregate(sub_table, by=list(sub_table[, feature_order[i]]), FUN=length)
    grouped <- aggregate(lineages$lineages, by=list(lineages$group), FUN=function(x) length(unique(x)))
    grouped <- grouped[grouped$x > 1,]

    list_non_unique <- aggregate(lineages$lineages, by=list(lineages$group), FUN=function(x) unique(x))

    if (nrow(grouped) > 0) {
      for(j in seq(1, nrow(grouped))) {
        trowgrp <- grouped[j,]
        tlist <- list_non_unique[list_non_unique[["Group.1"]] == trowgrp[["Group.1"]],]
        tkcount <- 1
        for(k in tlist[["x"]][[1]]) {
          tkrow <- rownames(lineages[lineages$lineages == k,])
          tax[tkrow, feature_order[i]] <- paste(tax[tkrow, feature_order[i]], tkcount, sep="_")
          tkcount <- tkcount + 1
        }
      }
    }
  }

  for( i in seq(ncol(tax))){
    tax[,i] = as.character(tax[,i])
  }

  replacing_na_obj_fData <- as.data.frame(tax, stringAsFactors=FALSE)
  
  
  message("replacing NA feature labels")
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
  
  # replacing_na_obj_fData
  return(replacing_na_obj_fData)
}

#' Method to fix feature labels across lineages
#'
#' @param tax fData from MRexperiment object
#' @return hierarchy with feature labels replaced
#' @export
#' 
replaceFeatureLabels <- function(tax) {
  message("replacing feature labels")
  feature_order <- colnames(tax)
  for(i in seq(2, length(feature_order))) {
    sub_table <- tax[, feature_order[1:i]]
    
    lineages <- apply(sub_table, 1, paste, collapse="_")
    lineages <- as.data.frame(lineages)
    lineages$group <- sub_table[, feature_order[i]]
    
    # grp <- sub_table[, by = feature_order[i]]
    # grouped <- aggregate(sub_table, by=list(sub_table[, feature_order[i]]), FUN=length)
    grouped <- aggregate(lineages$lineages, by=list(lineages$group), FUN=function(x) length(unique(x)))
    grouped <- grouped[grouped$x > 1,]
    
    list_non_unique <- aggregate(lineages$lineages, by=list(lineages$group), FUN=function(x) unique(x))
    
    if (nrow(grouped) > 0) {
      for(j in seq(1, nrow(grouped))) {
        trowgrp <- grouped[j,]
        tlist <- list_non_unique[list_non_unique[["Group.1"]] == trowgrp[["Group.1"]],]
        tkcount <- 1
        for(k in tlist[["x"]][[1]]) {
          tkrow <- rownames(lineages[lineages$lineages == k,])
          tax[tkrow, feature_order[i]] <- paste(tax[tkrow, feature_order[i]], tkcount, sep="_")
          tkcount <- tkcount + 1
        }
      }
    }
  }
  
  for( i in seq(ncol(tax))){
    tax[,i] = as.character(tax[,i]) 
  }
  
  return(as.data.frame(tax, stringAsFactors=FALSE))
}
