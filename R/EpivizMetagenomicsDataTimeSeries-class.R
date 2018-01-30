#' Data container for MRexperiment objects
#' 
#' Used to serve metagenomic data (used in e.g., icicle plots and heatmaps). Wraps
#' \code{\link[metagenomeSeq]{MRexperiment-class}} objects.
#' @importClassesFrom epivizrData EpivizData
#' @importFrom vegan diversity
#' @import data.table
#' @import digest
#' @import methods
#' @import httr
#' @export
#' @exportClass EpivizMetagenomicsDataTimeSeries
#' @examples
#' 
#' library(metagenomeSeq)
#' data(mouseData)
#' obj <- metavizr:::EpivizMetagenomicsData$new(mouseData, feature_order = colnames(fData(mouseData)))
#' 
EpivizMetagenomicsDataTimeSeries <- setRefClass("EpivizMetagenomicsDataTimeSeries",
  contains = "EpivizMetagenomicsData",
  fields = list(
    .alpha = "ANY",
    .feature_order = "character",
    .original_mr_exp = "ANY",
    .formula = "ANY",
    .time_series_class = "ANY",
    .spline_id = "character",
    .time = "ANY",
    .lvl = "ANY",
    .metavizControl = "ANY",
    .C = "ANY",
    .B = "ANY",
    .seed = "ANY",
    .fitThreshold = "ANY"
  ),
  methods = list(
    initialize=function(object, columns=NULL, control=metavizControl(), feature_order = NULL, formula = NULL, class = NULL, id = NULL, time = NULL, lvl = NULL,
                        C = NULL, B = NULL, seed = NULL, alpha = NULL, runFitTimeSeries = FALSE, fitThreshold = NULL, ...) {
      
      .self$.original_mr_exp <- object
      if(is.null(alpha)){
        .self$.alpha = 1.4
      }
      else{
        .self$.alpha <- alpha
      }
      if(is.null(feature_order)){
        .self$.feature_order <- ""
      }else{
        .self$.feature_order <- feature_order
      }
      .self$.formula <- formula
      .self$.time_series_class <- class
      if(is.null(id)){
        .self$.spline_id <- ""
      }else{
        .self$.spline_id <- id
      }
      .self$.metavizControl = control
      .self$.time <- time
      .self$.lvl <- lvl
      .self$.C <- C
      .self$.B <- B
      .self$.seed <- seed
      .self$.fitThreshold <- fitThreshold
      
      if (runFitTimeSeries == TRUE){
          spline_obj <- metagenomeSeq::fitMultipleTimeSeries(obj = object, formula = .self$.formula, class=.self$.time_series_class, id =.self$.spline_id, time = .self$.time, 
                                                             lvl = .self$.lvl, featureOrder = .self$.feature_order, C = .self$.C, B = .self$.B, seed = .self$.seed, alpha = .self$.alpha)          

        
        feature_order_agg <- .self$.feature_order[1:which(.self$.feature_order == .self$.lvl)]
        
        new_mrexp <- ts2MRexperiment(spline_obj, featureData = featureData(metagenomeSeq::aggregateByTaxonomy(object, lvl=.self$.lvl, featureOrder = feature_order_agg)), sampleNames = spline_obj[[2]]$fit$timePoints)
        
        .self$.feature_order <- feature_order_agg
        
        if(is.null(fitThreshold)){
          callSuper(object = new_mrexp, columns = columns, control = control, feature_order = feature_order_agg, ...)
          
        }else {
          splines_to_plot <- sapply(1:nrow(MRcounts(new_mrexp)), function(i) {max(abs(MRcounts(new_mrexp[i,]))) >= .self$.fitThreshold})
          
          splines_to_plot_indices <- which(splines_to_plot == TRUE)
          
          callSuper(object = new_mrexp[splines_to_plot_indices,], columns = columns, control = control, feature_order = feature_order_agg, ...)
        }
      } else{
          callSuper(object = object, columns = columns, control = control, feature_order = feature_order, ...)
      }

    },
    
    updateSplineAlpha = function(alpha = NULL){
      .self$.alpha <- alpha
      .self$initialize(object = .self$.original_mr_exp, columns = .self$.columns, control = .self$.metavizControl, feature_order = .self$.feature_order, formula = .self$.formula, class = .self$.time_series_class, id = .self$.spline_id, time = .self$.time, lvl = .self$.lvl,
                       C = .self$.C, B = .self$.B, seed = .self$.seed, alpha = .self$.alpha, runFitTimeSeries = TRUE, fitThreshold = .self$.fitThreshold)
    }
  )
)