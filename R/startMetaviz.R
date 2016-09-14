.register_all_metaviz_things <- function(app) {
  
  app$server$register_action("registerChartTypes", function(request_data) {
    app$chart_mgr$.register_available_chart_types(request_data$data)
  })
  
  app$server$register_action("getHierarchy", function(request_data) {
    datasource = request_data$datasourceGroup # TODO: Change to datasource
    nodeId = request_data$nodeId
    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getHierarchy(nodeId)
  })
  
  app$server$register_action("propagateHierarchyChanges", function(request_data) {
    datasource = request_data$datasourceGroup # TODO: Change to datasource
    selection = request_data$selection
    order = request_data$order
    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$propagateHierarchyChanges(selection, order)
  })
  
  app$server$register_action("getRows", function(request_data) {
    datasource = request_data$datasource
    seqName = request_data$seqName
    start = request_data$start
    end = request_data$end
    metadata = request_data$metadata
    
    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getRows(seqName, start, end, metadata)
  })
  
  app$server$register_action("getValues", function(request_data) {
    datasource = request_data$datasource
    measurement = request_data$measurement
    seqName = request_data$seqName
    start = request_data$start
    end = request_data$end
    
    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getValues(measurement, seqName, start, end)
  })
  
  app$server$register_action("getSeqInfos", function(request_data) {
    return(list(
      list("metavizr", 0, .Machine$integer.max)
    ))
  })
}


startMetaviz <- function(register_function = .register_all_metaviz_things, host="http://metaviz.cbcb.umd.edu", ...) {
  app = startEpiviz(register_function = register_function, host = host, ...)
  app
}
