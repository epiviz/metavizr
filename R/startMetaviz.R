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
  
  app$server$register_action("getCombined", function(request_data) {
    datasource = request_data$datasource
    measurements = request_data$measurements
    seqName = request_data$seqName
    start = request_data$start
    end = request_data$end
    order = request_data$order
    nodeSelection = request_data$selection
    selectedLevels = request_data$selectedLevels
    
    
    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getCombined(measurements, seqName, start, end, order, nodeSelection, selectedLevels)
  })
  
  app$server$register_action("getSeqInfos", function(request_data) {
    return(list(
      list("metavizr", 0, 100000)
    ))
  })
  
  app$server$register_action("partitions", function(request_data) {
    return(list(
      list("metavizr", 0, 100000)
    ))
  })
}

#' Start metaviz app and create \code{\link{EpivizApp}} object to manage connection.
#' 
#' @param host (character) host address to launch.
#' @param debug (logical) start metaviz app in debug mode.
#' @param workspace (character) a workspace id to load in the metaviz app on startup.
#' @param scripts (character) URLs for JavaScript plugin scripts to be imported when metaviz is loaded (see \url{http://epiviz.cbcb.umd.edu/help} for details).
#' @param gists (character) Ids for github gists (\url{http://gist.github.com}) containing JavaScript plugin scripts to
#'  be imported when metaviz is loaded (see \url{http://epiviz.cbcb.umd.edu/help} for details).
#' @param use_cookie (logical) use cookies within the epiviz app.
#' @param register_function (function) function used to register actions and charts on the metaviz app.
#' @param open_browser (logical) browse to the metaviz URL before exiting function.
#' @param server (EpivizServer) if not \code{NULL} use this object as underlying WebSocket and HTTP server
#' @param browser_fun (function) function used to browse to URL (\code{browseURL} by default)
#' @param ws_host (character) host address to use for websocket connection ("localhost" by default)
#' @param ... additional parameters passed to \code{\link[epivizrServer]{createServer}}.
#' 
#' @return An object of class \code{\link{EpivizApp}}
#' 
#' @import epivizr
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import epivizrData
#' @examples
#' # see package vignette for example usage
#' app <- startMetaviz(non_interactive=TRUE, open_browser=TRUE)
#' app$stop_app()
#' 
#' @export
startMetaviz <- function(host="http://metaviz.cbcb.umd.edu", register_function = .register_all_metaviz_things, ...) {
  chr="metavizr"
  start=1
  end=1000
  app = startEpiviz(host = host, register_function = register_function, chr=chr, start=start, end=end, ...)
  app
}

#' Start metaviz app in standalone (locally) and create \code{\link{EpivizApp}} object to manage connection.
#' 
#' @param branch (character) branch to pull from metaviz github repo to run standalone.
#' @param debug (logical) start metaviz app in debug mode.
#' @param workspace (character) a workspace id to load in the metaviz app on startup.
#' @param scripts (character) URLs for JavaScript plugin scripts to be imported when metaviz is loaded (see \url{http://epiviz.cbcb.umd.edu/help} for details).
#' @param gists (character) Ids for github gists (\url{http://gist.github.com}) containing JavaScript plugin scripts to
#'  be imported when metaviz is loaded (see \url{http://epiviz.cbcb.umd.edu/help} for details).
#' @param use_cookie (logical) use cookies within the epiviz app.
#' @param register_function (function) function used to register actions and charts on the metaviz app.
#' @param open_browser (logical) browse to the metaviz URL before exiting function.
#' @param server (EpivizServer) if not \code{NULL} use this object as underlying WebSocket and HTTP server
#' @param browser_fun (function) function used to browse to URL (\code{browseURL} by default)
#' @param ws_host (character) host address to use for websocket connection ("localhost" by default)
#' @param ... additional parameters passed to \code{\link[epivizrServer]{createServer}}.
#' 
#' @return An object of class \code{\link{EpivizApp}}
#' 
#' @import epivizrStandalone
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import epivizrData
#' @examples
#' # see package vignette for example usage
#' app <- startMetaviz(non_interactive=TRUE, open_browser=TRUE)
#' app$stop_app()
#' 
#' @export
startMetavizStandalone <- function(branch="metaviz-4.1", register_function = .register_all_metaviz_things, ...) {
  chr="metavizr"
  start=1
  end=1000
  seq <- Seqinfo(seqnames=chr,
                 seqlengths=100000,
                 isCircular=FALSE,
                 genome="metavizr")
  setStandalone(branch=branch)
  app = startStandalone(seqinfo=seq, register_function=register_function, chr=chr, start=start, end=end, ...)
  app
}
