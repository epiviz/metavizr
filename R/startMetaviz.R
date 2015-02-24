startMetaviz <- function(...) {
  mgr = startEpiviz(...)
  mgr$registerType("metagenomics",list(class="EpivizMetagenomicsData", description="Metagenomics data", input_class="MRexperiment"))

  mgr$registerAction("getHierarchy", function(mgr, msgData, ...) {
    datasource = msgData$datasourceGroup # TODO: Change to datasource
    nodeId = msgData$nodeId
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getHierarchy(nodeId)
  })


  mgr$registerAction("propagateHierarchyChanges", function(mgr, msgData, ...) {
    datasource = msgData$datasourceGroup # TODO: Change to datasource
    selection = msgData$selection
    order = msgData$order
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$propagateHierarchyChanges(selection, order)
  })

  mgr$registerAction("getRows", function(mgr, msgData, ...) {
    datasource = msgData$datasource
    seqName = msgData$seqName
    start = msgData$start
    end = msgData$end
    metadata = msgData$metadata

    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getRows(start, end, metadata)
  })

  mgr$registerAction("getValues", function(mgr, msgData, ...) {
    datasource = msgData$datasource
    measurement = msgData$measurement
    seqName = msgData$seqName
    start = msgData$start
    end = msgData$end

    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getValues(measurement, start, end)
  })

  mgr
}
