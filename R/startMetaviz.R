startMetaviz <- function(...) {
  mgr = startEpiviz(...)

  mgr$registerChartType("heatmap", "epiviz.plugins.charts.HeatmapPlot")
  mgr$registerChartType("scatterplot", "epiviz.plugins.charts.ScatterPlot")
  mgr$registerChartType("blocks", "epiviz.plugins.charts.BlocksTrack")
  mgr$registerChartType("line", "epiviz.plugins.charts.LineTrack")
  mgr$registerChartType("lineplot", "epiviz.plugins.charts.LinePlot")
  mgr$registerChartType("stackedplot", "epiviz.plugins.charts.StackedLinePlot")
  mgr$registerChartType("icicle", "epiviz.ui.charts.tree.Icicle")

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
    obj$getRows(seqName, start, end, metadata)
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
    obj$getValues(measurement, seqName, start, end)
  })

  mgr$registerAction("getSeqInfos", function(mgr, msgData, ...) {
    return(list(
      list("metavizr", 0, .Machine$integer.max)
    ))
  })

  mgr
}
