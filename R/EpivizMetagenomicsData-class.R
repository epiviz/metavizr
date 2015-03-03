EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains="EpivizData",
  fields=list(
    .taxonomy="MetavizTree",
    .levels="character",
    .maxDepth="numeric",
    .aggregateAtDepth="numeric",
    .lastRootId="ANY",

    .counts="ANY",
    .sampleAnnotation="ANY",

    .lastRequestRanges="list",
    .lastLeafInfos="list",
    .maxHistory="numeric"
  ),
  methods=list(
    initialize=function(object, maxDepth=3, aggregateAtDepth=3, maxHistory=3, ...) {
      # TODO: Some type checking
      .taxonomy <<- buildMetavizTree(object)
      .levels <<- .taxonomy$levels()
      .maxDepth <<- maxDepth
      .aggregateAtDepth <<- aggregateAtDepth
      .lastRootId <<- .taxonomy$root()$id()

      .counts <<- MRcounts(object, norm=TRUE, log=FALSE)
      .sampleAnnotation <<- pData(object)

      .maxHistory <<- maxHistory
      .lastRequestRanges <<- list()
      .lastLeafInfos <<- list()

      if (.aggregateAtDepth >= 0) {
        nodesAtDepth = .taxonomy$nodesAtDepth(.aggregateAtDepth)
        selection = list()
        for (node in nodesAtDepth) {
          selection[[node$id()]] = SelectionType$NODE
        }
        .taxonomy$updateSelection(selection)
      }

      callSuper(object=object, ...)
    },

    .taxonomyLevels=function(exp) {
      # TODO: Joe
      return(colnames(fData(exp))[c(3:9,1)])
    },

    .getSelectedLeaves=function(start, end) {
      ret = NULL
      for (i in rev(seq_along(.lastRequestRanges))) {
        if (.lastRequestRanges[[i]]$start == start || .lastRequestRanges[[i]]$end == end) {
          ret = .lastLeafInfos[[i]]
          break
        }
      }
      if (is.null(ret)) {
        requestRange = list(start=start, end=end)
        ret = Ptr$new(.taxonomy$selectedLeaves(start, end))
        .lastLeafInfos <<- c(.lastLeafInfos, ret)
        .lastRequestRanges <<- c(.lastRequestRanges, list(requestRange))

        if (length(.lastRequestRanges) > .maxHistory) {
          .lastRequestRanges <<- .lastRequestRanges[2:(.maxHistory+1)]
          .lastLeafInfos <<- .lastLeafInfos[2:(.maxHistory+1)]
        }
      }

      return(ret$.)
    },

    update=function(newObject, ...) {
      # TODO
      callSuper(newObject, ...)
    },
    plot=function(...) {
      # TODO
    }
  )
)

EpivizMetagenomicsData$methods(
  getMeasurements=function() {
    out <- lapply(colnames(.counts), function(sample) {
      list(id=sample,
           name=sample,
           type="feature",
           datasourceId=id,
           datasourceGroup=id,
           defaultChartType="heatmap",
           annotation=as.list(.sampleAnnotation[sample,]),
           minValue=0, # TODO
           maxValue=15, # TODO
           metadata=c(rev(.levels), "colLabel", "ancestors", "hierarchy-path"))
    })
    out
  },
  getHierarchy=function(nodeId) {
    root = NULL
    if (missing(nodeId) || is.null(nodeId)) { root = .taxonomy$root() }
    else {
      root = .taxonomy$parent(.taxonomy$node(nodeId))
      if (is.null(root)) { root = .taxonomy$root() }
    }
    .lastRootId <<- nodeId

    return(root$raw(maxDepth=.maxDepth))
  },
  propagateHierarchyChanges=function(selection, order) {
    if (missing(selection) && missing(order)) { return(getHierarchy(.lastRootId)) }

    if (!missing(selection)) {
      .taxonomy$updateSelection(selection)
    }

    if (!missing(order)) {
      .taxonomy$updateOrder(order)
    }

    .lastRequestRanges <<- list()
    .lastLeafInfos <<- list()

    getHierarchy(.lastRootId)
  },
  getRows=function(seqName, start, end, metadata) {
    leafInfos = .self$.getSelectedLeaves(start, end)
    leafAncestors = lapply(leafInfos, function(info) { .taxonomy$ancestors(info$node) })

    leafTaxonomies = list()
    for (i in seq_along(.levels)) {
      leafTaxonomies[[.levels[[i]]]] = list()
      for (j in seq_along(leafInfos)) {
        if (i > length(leafAncestors[[j]])) { next }
        leafTaxonomies[[.levels[[i]]]][[j]] = leafAncestors[[j]][[i]]
      }
    }

    ret = list(
      id = sapply(leafInfos, function(info) { info$realNodesBefore }),
      start=sapply(leafInfos, function(info) { info$start }),
      end=sapply(leafInfos, function(info) { info$start + info$node$nleaves() - 1 }),
      metadata = c(list(
        colLabel = sapply(leafInfos, function(info) { info$node$name() }),
        ancestors = sapply(leafAncestors, function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$name() }), collapse=",") }), # TODO: Use tree .ancestryByDepth
        "hierarchy-path" = sapply(leafAncestors, function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$id() }), collapse=",") })
      ), sapply(rev(.levels), function(level) { list(lapply(leafTaxonomies[[level]], function(node) { if(is.null(node)) { return("<NA>") }; node$name() })) }))
    )

    globalStartIndex = NULL
    if (length(leafInfos) > 0) { globalStartIndex = leafInfos[[1]]$realNodesBefore }
    if (length(leafInfos) == 1) {
      ret$id = list(ret$id)
      ret$start = list(ret$start)
      ret$end = list(ret$end)
      ret$metadata$colLabel = list(ret$metadata$colLabel)
      ret$metadata$ancestors = list(ret$metadata$ancestors)
      ret$metadata[["hierarchy-path"]] = list(ret$metadata[["hierarchy-path"]])
    }

    return(list(
      globalStartIndex = globalStartIndex,
      values = ret
    ))
  },
  getValues=function(measurement, seqName, start, end) {
    leafInfos = .self$.getSelectedLeaves(start, end)
    globalStartIndex = NULL
    if (length(leafInfos) > 0) { globalStartIndex = leafInfos[[1]]$realNodesBefore }
    ret = list(
      globalStartIndex = globalStartIndex,
      values = unname(lapply(leafInfos, function(info) {
        if (info$node$isLeaf()) { return(log2(.counts[info$node$leafIndex()+1, measurement] + 1)) }

        # TODO Joe: Currently, we compute mean of counts. What should we do instead?
        return(log2(1 + mean(.counts[(info$node$leafIndex()+1):(info$node$leafIndex()+info$node$nleaves()), measurement])))
      }))
    )
    return(ret)
  }
)
