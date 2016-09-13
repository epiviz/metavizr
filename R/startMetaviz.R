.register_all_metaviz_things <- function(app) {
  
  app$server$register_action("registerChartTypes", function(request_data) {
    app$chart_mgr$.register_available_chart_types(request_data$data)
  })
  
  app$server$register_action("getHierarchy", function(mgr, msgData, ...) {
    datasource = msgData$datasourceGroup # TODO: Change to datasource
    nodeId = msgData$nodeId
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getHierarchy(nodeId)
  })
  
  app$server$register_action("propagateHierarchyChanges", function(mgr, msgData, ...) {
    datasource = msgData$datasourceGroup # TODO: Change to datasource
    selection = msgData$selection
    order = msgData$order
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$propagateHierarchyChanges(selection, order)
  })
  
  app$server$register_action("getRows", function(mgr, msgData, ...) {
    datasource = msgData$datasource
    seqName = msgData$seqName
    start = msgData$start
    end = msgData$end
    metadata = msgData$metadata
    
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getRows(seqName, start, end, metadata)
  })
  
  app$server$register_action("getValues", function(mgr, msgData, ...) {
    datasource = msgData$datasource
    measurement = msgData$measurement
    seqName = msgData$seqName
    start = msgData$start
    end = msgData$end
    
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getValues(measurement, seqName, start, end)
  })
  
  app$server$register_action("getSeqInfos", function(mgr, msgData, ...) {
    return(list(
      list("metavizr", 0, .Machine$integer.max)
    ))
  })
}


startMetaviz <- function(register_function = .register_all_metaviz_things, host="http://metaviz.cbcb.umd.edu", ...) {
  app = startEpiviz(register_function = register_function, host = host, ...)
  app
}
