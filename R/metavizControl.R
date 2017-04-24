#' metavizr settings
#'
#' Default settings for the various plotting functions in metavizr.
#'
#' @param aggregateAtDepth Level of the tree to aggregate counts at by default.
#' @param aggregateFun Function to aggregate counts by at the aggregateAtDepth level.
#' @param valuesAnnotationFuns Function for error bars.
#' @param maxDepth Level of the tree to display by default in icicle view.
#' @param maxHistory Value for caching.
#' @param maxValue Maximum value to display.
#' @param minValue Minimum value to display.
#' @param title title.
#' @param n Number of OTUs to include in ranking.
#' @param rankFun Ranking function - single vector function.
#' @param norm Normalize MRexperiment object.
#' @param log Log tranformation of MRexperiment object.
#' @param featureSelection List of features to set as nodeSelections
#' @return List of setting parameters.
#' @examples
#' settings = metavizControl()
#'
#' @export
metavizControl<-function(aggregateAtDepth=3,
                         aggregateFun=function(x) colSums(x),
                         valuesAnnotationFuns=NULL,
                         maxDepth=4,
                         maxHistory=3,
                         maxValue=NULL,
                         minValue=NULL,
                         title="",
                         n=10000,
                         rankFun=stats::sd,
                         norm=TRUE,
                         log=FALSE,
                         featureSelection = NULL){
  
  list(aggregateAtDepth=aggregateAtDepth,
       aggregateFun=aggregateFun,
       valuesAnnotationFuns=valuesAnnotationFuns,
       maxDepth=maxDepth,
       maxHistory=maxHistory,
       maxValue=maxValue,
       minValue=minValue,
       title=title,
       n=n,
       rankFun=rankFun,
       norm=norm,
       log=log,
       featureSelection = featureSelection)
}
