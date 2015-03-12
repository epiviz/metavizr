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

    .minValue="ANY",
    .maxValue="ANY",
    .aggregateFun="ANY",
    .valuesAnnotationFuns="ANY",

    .lastRequestRanges="list",
    .lastLeafInfos="list",
    .lastValues="list",
    .maxHistory="numeric"
  ),
  methods=list(
    initialize=function(object, maxDepth=3, aggregateAtDepth=3, maxHistory=3, minValue=NULL, maxValue=NULL, aggregateFun=function(t) log2(1 + colSums(t)), valuesAnnotationFuns=NULL, ...) {
      # TODO: Some type checking
      .taxonomy <<- buildMetavizTree(object)
      .levels <<- .taxonomy$levels()
      .maxDepth <<- maxDepth
      .aggregateAtDepth <<- aggregateAtDepth
      .lastRootId <<- .taxonomy$root()$id()

      t = .taxonomy$taxonomyTable()
      .counts <<- MRcounts(object[rownames(t),],norm=TRUE,log=FALSE)


      # TODO: Make this consistent with the aggregateFun
      if (is.null(minValue)) {
        minValue = log2(min(.counts) + 1)
      }
      .minValue <<- minValue

      if (is.null(maxValue)) {
        maxValue = log2(max(.counts) + 1)
      }
      .maxValue <<- maxValue
      .aggregateFun <<- aggregateFun
      .valuesAnnotationFuns <<- valuesAnnotationFuns

      .sampleAnnotation <<- pData(object)

      .maxHistory <<- maxHistory
      .lastRequestRanges <<- list()
      .lastLeafInfos <<- list()
      .lastValues <<- list()

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
      if (length(.lastRequestRanges) > 0) {
        for (i in rev(seq_along(.lastRequestRanges))) {
          if (.lastRequestRanges[[i]]$start == start || .lastRequestRanges[[i]]$end == end) {
            ret = .lastLeafInfos[[i]]
            break
          }
        }
      }
      if (is.null(ret)) {
        requestRange = list(start=start, end=end)
        ret = Ptr$new(.taxonomy$selectedLeaves(start, end))

        if (.maxHistory > 0) {
          .lastLeafInfos <<- c(.lastLeafInfos, ret)
          .lastRequestRanges <<- c(.lastRequestRanges, list(requestRange))
          .lastValues <<- c(.lastValues, list(NULL))

          if (length(.lastRequestRanges) > .maxHistory) {
            .lastRequestRanges <<- .lastRequestRanges[2:(.maxHistory+1)]
            .lastLeafInfos <<- .lastLeafInfos[2:(.maxHistory+1)]
            .lastValues <<- .lastValues[2:(.maxHistory+1)] # Not yet a value
          }
        }
      }

      return(ret$.)
    },

    .getSelectedValues=function(measurement, start, end) {
      leafInfos = .getSelectedLeaves(start, end)
      index = -1
      values = NULL
      if (length(.lastRequestRanges) > 0) {
        for (i in rev(seq_along(.lastRequestRanges))) {
          if (.lastRequestRanges[[i]]$start == start || .lastRequestRanges[[i]]$end == end) {
            if (!is.null(.lastValues[[i]])) {
              values = .lastValues[[i]]
            } else {
              index = i
            }
            break
          }
        }
      }

      if (is.null(values)) {
        values = Ptr$new(list(values=unname(lapply(leafInfos, function(info) {
          return(.aggregateFun(.counts[(info$node$leafIndex()+1):(info$node$leafIndex()+info$node$nleaves()),, drop=FALSE]))
        }))))
        if (!is.null(.valuesAnnotationFuns)) {
          for (anno in names(.valuesAnnotationFuns)) {
            fun = .valuesAnnotationFuns[[anno]]
            values$.[[anno]] = unname(lapply(leafInfos, function(info) {
              return(fun(.counts[(info$node$leafIndex()+1):(info$node$leafIndex()+info$node$nleaves()),, drop=FALSE]))
            }))
          }
        }
      }
      if (index > 0) {
        .lastValues[[index]] <<- values
      }

      lapply(values$., function(vals) {
        lapply(vals, function(v) {
          if (length(dim(v)) > 0) { return(v[, measurement]) }
          return(v[[measurement]])
        })
      })
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

# Data analysis features
EpivizMetagenomicsData$methods(
  taxonomyTable=function() { .taxonomy$taxonomyTable() },
  calcNodeId=function(rowIndex, colIndex) { .taxonomy$calcNodeId(rowIndex, colIndex) },
  node=function(nodeId) { .taxonomy$node(nodeId) },
  parent=function(node) { .taxonomy$parent(node) },
  siblings=function(node) { .taxonomy$siblings(node) },

  changeAggregation=function(mgr, nodeId, aggregationType) {
    selection = list()
    selection[[nodeId]] = aggregationType
    .taxonomy$updateSelection(selection)
    mgr$.clearDatasourceGroupCache(.self, TRUE)
  },
  changeAggregationAtDepth=function(mgr, depth, aggregationType) {
    if (depth < 0) { return() }
    nodesAtDepth = .taxonomy$nodesAtDepth(depth)
    selection = list()
    for (node in nodesAtDepth) {
      selection[[node$id()]] = aggregationType
    }
    .taxonomy$updateSelection(selection)
    mgr$.clearDatasourceGroupCache(.self, TRUE)
  }
)

# Epiviz Websockets Protocol
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
           minValue=.minValue,
           maxValue=.maxValue,
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
      ), sapply(rev(.levels), function(level) {
        r = list(lapply(leafTaxonomies[[level]], function(node) { if(is.null(node)) { return("<NA>") }; node$name() }))
        if (length(r[[1]]) == 0) {
          r[[1]] = lapply(seq_along(leafInfos), function(i) { "<NA>" })
        }
        return(r)
      }))
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
      values = .self$.getSelectedValues(measurement, start, end)
    )
    return(ret)
  }
)
