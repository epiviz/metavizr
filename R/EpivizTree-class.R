EpivizTree <- setRefClass("EpivizTree",
  fields=list(
    nodes="ANY",
    nodesById="ANY",           # nullable environment
    nodesByDepthAndName="ANY"
  ),
  methods=list(
    initialize=function(nodes, ...) {
      if (missing(nodes)) { nodes <<- NULL; return() }
      nodes <<- nodes
      nodesById <<- list()
      nodesByDepthAndName <<- list()

      for (node in nodes) {
        nodesById[[node$id]] <<- node
        if (length(nodesByDepthAndName) < node$depth + 1 || is.null(nodesByDepthAndName[[node$depth + 1]])) {
          nodesByDepthAndName[[node$depth + 1]] <<- list()
        }
        nodesByDepthAndName[[node$depth + 1]][[node$name]] <<- node
      }
    },
    root=function() { nodesByDepthAndName[[1]][[1]] },

    children=function(node) {

      lapply(node$childrenIds, function(id) { nodesById[[id]] })
    },

    traverse=function(callback, node=root()) {
      callback(node)
      if (is.null(node)) { return() }
      if (length(node$childrenIds) == 0) { return() }
      for (child in children(node)) {
        traverse(callback, child)
      }
    },

     node=function(id) {
      if (missing(id)) { id = root()$id }
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
      if (length(node$childrenIds) == 0) { return(list(node)) }
      ret = list()
      for (child in children(node)) {
        ret = c(ret, leaves(child))
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
      for (nodeId in names(order)) {
        if (is.null(nodesById[[nodeId]])) { next }
        nodesById[[nodeId]]$order <<- order[[nodeId]]
        parent = .self$parent(nodesById[[nodeId]])
        if (!is.null(parent)) {
          modifiedParents[[parent$id]] = parent
        }
      }

      # Second, sort
      for (node in modifiedParents) {
        o = sapply(children(node), function(child) { child$order })
        node$childrenIds = node$childrenIds[order(o)]
      }
    },

    # Get a list of the leaves of the given subtree, respecting the "selectionType" field
    selectedLeaves=function(node=root()) {
      if (is.null(node)) { return(NULL) }
      if (node$selectionType == SelectionType$NONE) { return(NULL) }
      if (length(node$childrenIds) == 0 || node$selectionType == SelectionType$NODE) { return(list(node)) }
      ret = list()
      for (child in children(node)) {
        ret = c(ret, selectedLeaves(child))
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
        ret[[node$depth + 1]] = node
        node = parent(node=node)
      }
      return(ret)
    },


    build=function(node=root(), filter=function(node) { TRUE }) {
      recurse=function(node, filter) {
        if (is.null(node)) { return(NULL) }
        if (!filter(node)) { return(NULL) }
        ret = node$raw()
        if (length(node$childrenIds) == 0) { return(list(ret)) }

        childrenList = list()
        for (child in children(node)) {
          childrenList = c(childrenList, recurse(child, filter))
        }
        if (length(childrenList) > 0) { ret$children = childrenList }
        return(list(ret))
      }
      ret = recurse(node, filter)
      if (is.null(ret)) { return(ret) }

      return(ret[[1]])
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

  .tableToTree <- function(t, colIndex=1, parent=NULL) {
    if (colIndex > dim(t)[2]) {
      if (!is.null(parent)) { parent$nleaves = 1 }
      return(NULL)
    }

    parentId = NULL
    if (!is.null(parent)) { parentId = parent$id }

    groups = by(t, t[, colIndex], list, simplify=F)
    names = names(groups)
    nodes = lapply(seq_along(names), function(i) { EpivizNode$new(name=names[i], parentId=parentId, depth=colIndex-1, taxonomy=colnames(t)[colIndex], order=i) })
    if (!is.null(parent)) {
      parent$nchildren = length(nodes)
      parent$childrenIds = lapply(nodes, function(node) { node$id })
    }

    descendants = unlist(lapply(seq_along(groups), function(i) {
      node = nodes[[i]]
      children = .tableToTree(groups[[i]][[1]], colIndex+1, node)
      if (!is.null(parent)) {
        parent$nleaves = parent$nleaves + node$nleaves
      }
      children
    }))
    ret = c(nodes, descendants)

    return(ret)
  }

  return(EpivizTree$new(.tableToTree(tax)))
})
