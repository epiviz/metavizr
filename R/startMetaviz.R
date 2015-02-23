startMetaviz <- function(...) {
  mgr = startEpiviz(...)
  mgr$registerType("metagenomics",list(class="EpivizMetagenomicsData", description="Metagenomics data", input_class="MRexperiment"))
  mgr$registerAction("getHierarchy", function(mgr, msgData, ...) {
    datasource = msgData$datasourceGroup # TODO: Change to datasource
    obj <- mgr$.findDatasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getHierarchy()
  })

  mgr
}
