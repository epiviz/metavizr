EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains="EpivizData",
  fields=list(
    taxonomy="MetavizTree",
    levels="character",
    maxDepth="numeric",
    lastRootId="ANY",

    counts="ANY",
    sampleAnnotation="ANY",

    selectedNodes="ANY",
    selectedValues="ANY",
    selectedNodesRanges="ANY",
    selectedNodesAncestors="ANY",
    selectedNodesTaxonomies="ANY",

    lastRequestRange="ANY",
    lastLeafInfos="Ptr"
  ),
  methods=list(
    initialize=function(object, ...) {
      # TODO: Some type checking
      taxonomy <<- buildMetavizTree(object) # TODO: Are the counts needed in the tree or in the node?
      levels <<- taxonomy$levels()
      maxDepth <<- 3 # TODO Make it customizable
      lastRootId <<- taxonomy$root()$id

      counts <<- MRcounts(object, norm=TRUE, log=TRUE)
      sampleAnnotation <<- pData(object)

#       # Change the counts row names to the ids of the leaves corresponding to them
#       idsByNames = list()
#       leaves = taxonomy$leaves()
#       for (leaf in leaves) {
#         idsByNames[[leaf$name]] = leaf$id
#       }
#       counts <<- counts[names(idsByNames),]
#       rownames(counts) <<- lapply(rownames(counts), function(name) { idsByNames[[name]] })

#       .self$.updateSelection()

      lastRequestRange <<- NULL
      lastLeafInfos <<- Ptr$new(NULL)

      callSuper(object=object, ...)
    },

    .taxonomyLevels=function(exp) {
      # TODO: Joe
      return(colnames(fData(exp))[c(3:9,1)])
    },

    .getSelectedLeaves=function(start, end) {
      if (is.null(lastRequestRange) || lastRequestRange$start != start || lastRequestRange$end != end) {
        lastRequestRange <<- list(start=start, end=end)
        lastLeafInfos <<- Ptr$new(taxonomy$selectedLeaves(start, end))
      }

      return(lastLeafInfos$.)
    },
#     .updateSelection=function() {
#       selectedNodes <<- taxonomy$selectedLeaves()
#
#       i = 0
#       nodesRanges = list()
#       for (node in selectedNodes) {
#         nodesRanges = c(nodesRanges, list(c(i, node$nleaves)))
#         i = i + node$nleaves
#       }
#       selectedNodesRanges <<- nodesRanges
#
#       selectedNodesAncestors <<- lapply(selectedNodes, function(node) { taxonomy$ancestors(node) })
#       nodesTaxonomies = list()
#       for (i in seq_along(levels)) {
#         nodesTaxonomies[[levels[[i]]]] = list()
#         for (j in seq_along(selectedNodes)) {
#           browser(expr=(j > length(selectedNodesAncestors)))
#           if (i > length(selectedNodesAncestors[[j]])) { next }
#           nodesTaxonomies[[levels[[i]]]][[j]] = selectedNodesAncestors[[j]][[i]]
#         }
#       }
#       selectedNodesTaxonomies <<- nodesTaxonomies
#
#       selection = .getSelectedValues()
#       selectedValues <<- sapply(selection, function(v) { v$. })
#     },

    .getSelectedValues=function(node=taxonomy$root()) {
      if (node$selectionType == SelectionType$NONE) {
        #return(list(rows=list(), values=list()))
        return(list())
      }
      if (node$nchildren == 0) {
        ret = list()
        ret[[node$id]] = Ptr$new(counts[node$id, ])
        return(ret)
        #return(list(rows=list(node$id), values=list(Ptr$new(counts[node$id, ]))))
      }

      #result = list(rows=list(), values=list())
      result = list()

      if (node$selectionType == SelectionType$LEAVES) {
        for (child in taxonomy$children(node)) {
          childResult = .getSelectedValues(child)
          #result$rows = c(result$rows, childResult$rows)
          #result$values = c(result$values, childResult$values)
          result = c(result, childResult)
        }
        return(result)
      }

      if (node$selectionType == SelectionType$NODE) {
        #result$rows = list(node$id)
        values = list()
        for (child in taxonomy$children(node)) {
          childResult = .getSelectedValues(child)
          values = c(values, childResult)
        }
        #result$values = list(Reduce(function(v1, v2) { return(Ptr$new(v1$. + v2$.)) }, values)$. / length(values))
        result[[node$id]] = Ptr$new(Reduce(function(v1, v2) { return(Ptr$new(v1$. + v2$.)) }, values)$. / length(values))
        return(result)
      }
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
    out <- lapply(colnames(counts), function(sample) {
      list(id=sample,
           name=sample,
           type="feature",
           datasourceId=id,
           datasourceGroup=id,
           defaultChartType="heatmap",
           annotation=as.list(sampleAnnotation[sample,]),
           minValue=0, # TODO
           maxValue=15, # TODO
           metadata=c(rev(levels), "colLabel", "ancestors", "hierarchy-path"))
    })
    out
  },
  getHierarchy=function(nodeId) {
    root = NULL
    if (missing(nodeId) || is.null(nodeId)) { root = taxonomy$root() }
    else {
      root = taxonomy$parent(taxonomy$node(nodeId))
      if (is.null(root)) { root = taxonomy$root() }
    }
    lastRootId <<- nodeId

    browser(expr=is.null(root$nleaves))
    return(root$raw(maxDepth=maxDepth))
  },
  propagateHierarchyChanges=function(selection, order) {
    if (missing(selection) && missing(order)) { return(getHierarchy(lastRootId)) }

    if (!missing(selection)) {
      taxonomy$updateSelection(selection)
    }

    if (!missing(order)) {
      taxonomy$updateOrder(order)

      # TODO!!!
      # Update order in counts data as well
      # But if we update the order in the data, ids will have to change too!
      # So big problem now: how do we figure out how many leaves there are total, and how many have not been selected left and right respectively
    }

    #.self$.updateSelection()

    getHierarchy(lastRootId)
  },
  getRows=function(seqName, start, end, metadata) {
#     leafInfos = taxonomy$selectedLeaves(0, taxonomy$root()$nleaves)
#     startIndex = NULL
#     endIndex = NULL
#     for (i in seq_along(leafInfos)) {
#       info = leafInfos[[i]]
#       #range = selectedNodesRanges[[i]]
#       #if (range[1] <= end && range[1] + range[2] > start) {
#       if (info$start <= end && info$start + info$node$nleaves > start) {
#         if (is.null(startIndex)) { startIndex = i }
#         endIndex = i
#       }
#     }
#     if (is.null(startIndex) || is.null(endIndex)) {
#       return(list(
#         globalStartIndex = NULL,
#         values = list(
#           id = list(),
#           start = list(),
#           end = list(),
#           metadata = c(list(
#             colLabel = list(),
#             ancestors = list(),
#             "hierarchy-path" = list()
#           ), sapply(levels, function(level) { list() }))
#         )
#       ))
#     }
    # TODO: Aici
    leafInfos = .self$.getSelectedLeaves(start, end) #taxonomy$selectedLeaves(start, end)
    leafAncestors = lapply(leafInfos, function(info) { taxonomy$ancestors(info$node) })

    leafTaxonomies = list()
    for (i in seq_along(levels)) {
      leafTaxonomies[[levels[[i]]]] = list()
      for (j in seq_along(leafInfos)) {
        if (i > length(leafAncestors[[j]])) { next }
        leafTaxonomies[[levels[[i]]]][[j]] = leafAncestors[[j]][[i]]
      }
    }

    ret = list(
      #id = 0 : (length(leafInfos) - 1),
      #id = startIndex:endIndex,
      id = sapply(leafInfos, function(info) { info$realNodesBefore }),
      #start = sapply(selectedNodesRanges[startIndex:endIndex], function(range) { range[1] }),
      #end = sapply(selectedNodesRanges[startIndex:endIndex], function(range) { range[1] + range[2] - 1 }),
      start=sapply(leafInfos, function(info) { info$start }),
      end=sapply(leafInfos, function(info) { info$start + info$node$nleaves - 1 }),
      metadata = c(list(
        colLabel = sapply(leafInfos, function(info) { info$node$name }),
        ancestors = sapply(leafAncestors, function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$name }), collapse=",") }),
        "hierarchy-path" = sapply(leafAncestors, function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$id }), collapse=",") })
      ), sapply(rev(levels), function(level) { list(lapply(leafTaxonomies[[level]], function(node) { node$name })) }))
    )

#     if (endIndex - startIndex == 0) {
#       ret$id = list(ret$id)
#       ret$start = list(ret$start)
#       ret$end = list(ret$end)
#       for (key in names(ret$metadata)) {
#         ret$metadata[[key]] = list(ret$metadata[[key]])
#       }
#     }
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
#     startIndex = NULL
#     endIndex = NULL
#     for (i in seq_along(selectedNodes)) {
#       range = selectedNodesRanges[[i]]
#       if (range[1] <= end && range[1] + range[2] > start) {
#         if (is.null(startIndex)) { startIndex = i }
#         endIndex = i
#       }
#     }
#     if (is.null(startIndex) || is.null(endIndex)) {
#       return(list(
#         globalStartIndex = NULL,
#         values = list()
#       ))
#     }
    leafInfos = .self$.getSelectedLeaves(start, end) #taxonomy$selectedLeaves(start, end)
    globalStartIndex = NULL
    if (length(leafInfos) > 0) { globalStartIndex = leafInfos[[1]]$realNodesBefore }
    ret = list(
      globalStartIndex = globalStartIndex,
      values = unname(lapply(leafInfos, function(info) {
        if (info$node$isLeaf()) { return(counts[info$node$leafIndex+1, measurement]) }
        return(mean(counts[(info$node$leafIndex+1):(info$node$leafIndex+info$node$nleaves), measurement]))
      }))
    )
#     if (length(leafInfos) == 1) { ret$values = list(ret$values) }
    return(ret)
  }
)
