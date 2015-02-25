EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains="EpivizData",
  fields=list(
    taxonomy="EpivizTree",
    levels="character",
    maxDepth="numeric",
    lastSubtree="ANY",

    counts="ANY",
    sampleAnnotation="ANY",

    selectedNodes="ANY",
    selectedValues="ANY",
    selectedNodesRanges="ANY",
    selectedNodesAncestors="ANY",
    selectedNodesTaxonomies="ANY"
  ),
  methods=list(
    initialize=function(object, ...) {
      # TODO: Some type checking
      levels <<- .self$.extractTaxonomyLevels(object)
      taxonomy <<- EpivizTree$new(.self$.MRexperimentToTree(object))
      maxDepth <<- 4 # TODO Make it customizable
      lastSubtree <<- taxonomy$node()

      counts <<- MRcounts(object, norm=TRUE, log=TRUE)
      sampleAnnotation <<- pData(object)

      # Change the counts row names to the ids of the leaves corresponding to them
      idsByNames = list()
      leaves = taxonomy$leaves()
      for (leaf in leaves) {
        idsByNames[[leaf$name]] = leaf$id
      }
      counts <<- counts[names(idsByNames),]
      rownames(counts) <<- lapply(rownames(counts), function(name) { idsByNames[[name]] })

      selectedNodes <<- taxonomy$selectedLeaves()
      selectedValues <<- lapply(selectedNodes, function(node) { counts[node$id, ] })

      i = 0
      nodesRanges = list()
      for (node in selectedNodes) {
        nodesRanges = c(nodesRanges, list(c(i, node$nleaves)))
        i = i + node$nleaves
      }
      selectedNodesRanges <<- nodesRanges

      selectedNodesAncestors <<- lapply(selectedNodes, function(node) { taxonomy$ancestors(node) })
      nodesTaxonomies = list()
      for (i in 1:length(levels)) {
        nodesTaxonomies[[levels[[i]]]] = list()
        for (j in 1:length(selectedNodes)) {
          nodesTaxonomies[[levels[[i]]]][[j]] = selectedNodesAncestors[[j]][[i]]
        }
      }
      selectedNodesTaxonomies <<- nodesTaxonomies

      callSuper(object=object, ...)
    },
    .MRexperimentToTree=function(exp) {
      # TODO Joe
      filteredExp = filterData(exp,depth=6000,present=25)
      taxLvls = colnames(fData(filteredExp))[c(3:9,1)]
      tax = fData(filteredExp)[,taxLvls]
      tax[,1] = "Bacteria"
      table = tax

      return(.tableToTree(table))
    },
    .tableToTree=function(t) {
      tableToListTree <- function(t, colIndex=1) {
        if (colIndex > dim(t)[2]) { return(NULL) }
        groups = by(t, t[, colIndex], list, simplify=F)
        nodes = lapply(groups, function(group){
          tableToListTree(group[[1]], colIndex + 1)
        })
        return(nodes)
      }
      listTreeToJSONTree <- function(node, globalDepth=0) {
        nodeId = .generatePseudoGUID(6)
        ret = list(
          name=names(node),
          id=nodeId,
          parentId=NULL,
          depth=globalDepth,
          globalDepth=globalDepth,
          taxonomy=colnames(t)[globalDepth+1],
          nchildren=length(node[[1]]),
          size=1,
          selectionType=SelectionType$LEAVES
        )
        nleaves = 0
        if (length(node[[1]]) > 0) {
          children = c()
          for (i in 1:length(node[[1]])) {
            child = listTreeToJSONTree(node[[1]][i], globalDepth + 1)
            child$order = i - 1
            child$parentId = nodeId
            children = c(children, list(child))
            nleaves = nleaves + child$nleaves
          }
          ret$children = children
        } else {
          nleaves = 1
        }
        ret$nleaves = nleaves
        return(ret)
      }

      listTree = tableToListTree(t)
      tree = listTreeToJSONTree(listTree)
      return(tree)
    },
    .extractTaxonomyLevels=function(exp) {
      # TODO: Joe
      return(colnames(fData(exp))[c(3:9,1)])
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
    if (missing(nodeId) || is.null(nodeId)) { root = taxonomy$node() }
    else {
      root = taxonomy$parent(id = nodeId)
      if (is.null(root)) { root = taxonomy$node() }
    }

    ret = taxonomy$build(function(node) {
      if (is.null(node) || node$globalDepth - root$globalDepth >= maxDepth) { return(NULL) }
      node
    }, root)
    lastSubtree <<- ret
    ret
  },
  propagateHierarchyChanges=function(selection, order) {
    if (missing(selection) && missing(order)) { return(lastSubtree) }
    ret = lastSubtree
    if (!missing(selection)) {
      taxonomy <<- taxonomy$updateSelection(selection)
    }

    if (!missing(order)) {
      taxonomy <<- taxonomy$updateOrder(order)
    }

    ret = getHierarchy(lastSubtree$id)
    lastSubtree <<- ret
    ret
  },
  getRows=function(seqName, start, end, metadata) {
    startIndex = NULL
    endIndex = NULL
    for (i in 1:length(selectedNodes)) {
      range = selectedNodesRanges[[i]]
      if (range[1] <= end && range[1] + range[2] > start) {
        if (is.null(startIndex)) { startIndex = i }
        endIndex = i
      }
    }
    if (is.null(startIndex) || is.null(endIndex)) {
      return(list(
        globalStartIndex = NULL,
        values = list(
          id = list(),
          start = list(),
          end = list(),
          metadata = c(list(
            colLabel = list(),
            ancestors = list(),
            "hierarchy-path" = list()
          ), sapply(levels, function(level) { list() }))
        )
      ))
    }
    ret = list(
      id = start : (start + endIndex - startIndex),
      start = sapply(selectedNodesRanges[startIndex:endIndex], function(range) { range[1] }),
      end = sapply(selectedNodesRanges[startIndex:endIndex], function(range) { range[1] + range[2] - 1 }),
      metadata = c(list(
        colLabel = sapply(selectedNodes[startIndex:endIndex], function(node) { node$name }),
        ancestors = sapply(selectedNodesAncestors[startIndex:endIndex], function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$name }), collapse=",") }),
        "hierarchy-path" = sapply(selectedNodesAncestors[startIndex:endIndex], function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$id }), collapse=",") })
      ), sapply(rev(levels), function(level) { list(lapply(selectedNodesTaxonomies[[level]][startIndex:endIndex], function(node) { node$name })) }))
    )

    if (length(startIndex:endIndex) == 1) {
      ret$id = list(ret$id)
      ret$start = list(ret$start)
      ret$end = list(ret$end)
      for (key in names(ret$metadata)) {
        ret$metadata[[key]] = list(ret$metadata[[key]])
      }
    }

    return(list(
      globalStartIndex = start,
      values = ret
    ))
  },
  getValues=function(measurement, seqName, start, end) {
    # TODO
  }
)
