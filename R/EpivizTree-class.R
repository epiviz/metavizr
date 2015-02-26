EpivizTree <- setRefClass("EpivizTree",
  fields=list(
    container="environment",  # Using environment to contain the root, to prevent RStudio from crashing
    nodesById="ANY"           # nullable environment
  ),
  methods=list(
    initialize=function(root, ...) {
      container <<- new.env()
      if (missing(root)) { container$root <<- NULL; return() }
      container$root <<- root
      nodesById <<- NULL
    },
    root=function() { container$root },

    raw=function() {
      root = root()
      if (!is.null(root)) { return(NULL) }
      root$raw()
    },

    traverse=function(callback, node=root()) {
      callback(node)
      if (is.null(node)) { return() }
      if (length(node$children) == 0) { return() }
      for (child in node$children) {
        traverse(callback, child)
      }
    },

    node=function(id) {
      if (missing(id)) { id = root()$id }
      if (is.null(nodesById)) {
        nodesById <<- new.env()
        .self$traverse(function(node) { nodesById[[node$id]] <<- node })
      }
      return(nodesById[[id]])
    },

    parent=function(node=root()) {
      if (is.null(node)) { return(NULL) }
      if (!is.null(node$parentId)) {
        return(.self$node(node$parentId))
      }
      return(NULL)
    },

    leaves=function(node=root()) {
      if (is.null(node)) { return(NULL) }
      if (node$nchildren == 0) { return(list(node)) }
      ret = list()
      for (i in 1:node$nchildren) {
        ret = c(ret, leaves(node$children[[i]]))
      }
      return(ret)
    },

    # selection: @list {nodeId -> selectionType}
    updateSelection=function(selection) {
      if (missing(selection)) { return() }

      for (nodeId in names(selection)) {
        if (is.null(nodesById[[nodeId]])) { next }
        nodesById[[nodeId]]$selectionType <<- selection[[nodeId]]
      }
    },

    # order: @list {nodeId -> order index}
    updateOrder=function(order) {
      if (missing(order)) { return() }

      modifiedParents = list()

      # First, reassign order
      for (nodeId in names(selection)) {
        if (is.null(nodesById[[nodeId]])) { next }
        nodesById[[nodeId]]$order <<- order[[nodeId]]
        parent = .self$parent(nodesById[[nodeId]])
        if (!is.null(parent)) {
          modifiedParents[[parent$id]] = parent
        }
      }

      # Second, sort
      for (node in modifiedParents) {
        o = sapply(node$children, function(child) { child$order })
        node$children = node$children[order(o)]
      }
    },

    # Get a list of the leaves of the given subtree, respecting the "selectionType" field
    selectedLeaves=function(node=root()) {
      if (is.null(node)) { return(NULL) }
      if (node$selectionType == SelectionType$NONE) { return(NULL) }
      if (length(node$children) == 0 || node$selectionType == SelectionType$NODE) { return(list(node)) }
      ret = list()
      for (i in 1:length(node$children)) {
        ret = c(ret, selectedLeaves(node$children[[i]]))
      }
      return(ret)
    },

    # Gets a list of ancestor nodes in the tree for the given node
    ancestors=function(node=root(), inclusive=TRUE) {
      if (is.null(node)) { return(NULL) }
      ret = list()
      if (inclusive) {
        ret[[node$depth + 1]] = node
      }

      node = parent(node)
      while (!is.null(node)) {
        ret[[node$globalDepth + 1]] = node
        node = parent(node=node)
      }
      return(ret)
    },

    filter=function(node=root(), filterCallback) {
      if (is.null(node)) { return(NULL) }
      if (!filterCallback(node)) { return(NULL) }

      copy = node$copy(recursive=FALSE)
      if (length(node$children) == 0) { return(copy) }
      for (child in node$children) {
        copy$children = c(copy$children, filter(child, filterCallback))
      }
      return(copy)
    }
  )
)

# Used by jsonlite::toJSON
# if (!is.null(getGeneric("asJSON"))) {
#   setMethod(getGeneric("asJSON"), "EpivizTree", function(x, ...){
#     return(epivizr::toJSON(x$root(), ...))
#   })
# }

setGeneric("buildEpivizTree", signature=c("object"),
           function(object, ...) standardGeneric("buildEpivizTree"))

setMethod("buildEpivizTree", "MRexperiment", function(object, ...) {
  # TODO Joe
  filteredExp = filterData(object,depth=6000,present=25)
  taxLvls = colnames(fData(filteredExp))[c(3:9,1)]
  tax = fData(filteredExp)[,taxLvls]
  tax[,1] = "Bacteria"
  table = tax

  .addChild=function(parent, child) {
    if (!is.null(parent$childrenById[[child$id]])) { return(parent$childrenById[[child$id]]) }
    if (!is.null(parent$childrenByName[[child$name]])) { return(parent$childrenByName[[child$name]]) }

    parent$childrenById[[child$id]] = child
    parent$childrenByName[[child$name]] = child

    parent$children = c(parent$children, child)
    parent$nchildren = parent$nchildren + 1
    return(child)
  }

  .computeNLeaves=function(node) {
    if (length(node$children) == 0) { node$nleaves = 1; return() }
    for (child in node$children) {
      .computeNLeaves(child)
      node$nleaves = node$nleaves + child$nleaves
    }
  }

  .tableToTree <- function(t, colIndex=1) {
    root = EpivizNode$new(name=t[1,1], taxonomy=colnames(t)[1])

    for (j in 1:dim(t)[1]) {
      node = root
      for (i in 2:dim(t)[2]) {
        nodeName = as.character(t[j, i])
        if (is.na(nodeName)) { nodeName = "NA" }

        child = EpivizNode$new(
          name=nodeName,
          taxonomy=colnames(t)[i],
          parentId=node$id,
          depth=node$depth + 1,
          order=node$nchildren
        )

        child = .addChild(node, child)
        node = child
      }
    }

    .computeNLeaves(root)

    return(EpivizTree$new(root))
  }

  return(.tableToTree(tax))

})
