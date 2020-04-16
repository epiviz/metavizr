#' Generic method to register data to the epiviz data server
#' 
#' @param object The object to register to data server
#' @param columns Name of columns containing data to register
#' @param type leafCounts, if data objects has counts at leaf level
#'        or innerNodeCounts, if data object has counts at inner nodes
#' @param ... Additonal arguments passed to object constructors
#' @return An \code{\link{EpivizMetagenomicsData-class}} object 
#' @import metagenomeSeq
#' @importMethodsFrom epivizrData register
#' 
setMethod("register", "MRexperiment", function(object, type="LeafCounts", columns=NULL, ...) {
  if(type == "LeafCounts"){
    return(EpivizMetagenomicsData$new(object=object, columns=columns, ...))
  } else if(type == "innerNodeCounts") {
    return(EpivizMetagenomicsDataInnerNodes$new(object=object, columns=columns, ...))
  } else if(type == "TimeSeries"){ 
    return(EpivizMetagenomicsDataTimeSeries$new(object=object, columns=columns, ...)) 
  } 
})

#' Generic method to register data to the epiviz data server
#' 
#' @param object The object to register to data server
#' @param type leafCounts, if data objects has counts at leaf level
#'        or innerNodeCounts, if data object has counts at inner nodes
#' @param ... Additonal arguments passed to object constructors
#' @return An \code{\link{phyloseq-class}} object 
#' @import metagenomeSeq
#' @importClassesFrom phyloseq phyloseq
#' @importFrom phyloseq phyloseq_to_metagenomeSeq
#' @importMethodsFrom epivizrData register
#' 
setMethod("register", "phyloseq", function(object, type="LeafCounts", ...) {
  phy_obj <- phyloseq_to_metagenomeSeq(physeq = object, ...)
  if(type == "LeafCounts"){
    return(EpivizMetagenomicsData$new(object=phy_obj, ...))
  } else if(type == "innerNodeCounts") {
    return(EpivizMetagenomicsDataInnerNodes$new(object=phy_obj, ...))
  } else if(type == "TimeSeries"){ 
    return(EpivizMetagenomicsDataTimeSeries$new(object=phy_obj, ...)) 
  }
})


#' Generic method to register data to the epiviz data server
#' 
#' @param object The object to register to data server
#' @param columns Name of columns containing data to register
#' @param ... Additonal arguments passed to object constructors
#' @return An \code{\link{EpivizMetagenomicsData-class}} object 
#' @import TreeSummarizedExperiment
#' @importMethodsFrom epivizrData register
#' 
setMethod("register", "TreeSummarizedExperiment", function(object, tree="row", columns=NULL, ...) {

  # counts <- t(assays(object)$counts)
  # tree <- colData(object)
  # rownames(tree) <- rownames(counts)
  # annotations <- rowData(object)
  # rownames(annotations) <- colnames(counts)
  # 
  # sExp <- newMRexperiment(counts, featureData = AnnotatedDataFrame(tree), phenoData = AnnotatedDataFrame(annotations))

  return(EpivizTreeData$new(object=object, tree=tree, columns=columns, ...))
})
