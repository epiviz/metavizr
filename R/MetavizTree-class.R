MetavizTree <- setRefClass("MetavizTree",
  fields=list(
    nodesById="ANY",           # nullable environment
    nodesByDepthAndName="ANY",
    taxonomyDF="ANY"
  ),
  methods=list(
    initialize=function(taxonomyDF, ...) {
      if (missing(taxonomyDF)) { taxonomyDF <<- NULL }

      taxonomyDF <<- taxonomyDF

      nodesById <<- list()
      nodesByDepthAndName <<- list()

      childrenGroups = by(taxonomyDF, taxonomyDF[, 2], list, simplify=F)
      root = EpivizNode$new(name=taxonomyDF[1,1], taxonomy=colnames(taxonomyDF)[1], nchildren=length(childrenGroups), nleaves=dim(taxonomyDF)[1])
      nodesById[[root$id]] <<- root
      nodesByDepthAndName[[root$depth + 1]] <<- list()
      nodesByDepthAndName[[root$depth + 1]][[root$name]] <<- root
    },
    root=function() { nodesByDepthAndName[[1]][[1]] },

    children=function(node) {
      if (is.null(node)) { return(NULL) }
      if (node$nchildren > 0 && length(node$childrenIds) == 0) {
        # Generate the children
        nodeDF = taxonomyDF[node$leafIndex:(node$leafIndex + node$nleaves),]
        childrenGroups = by(nodeDF, nodeDF[, node$depth + 2], list, simplify=F)
        names = names(childrenGroups)
        env = new.env()
        env$lastLeafIndex = node$leafIndex
        children = lapply(seq_along(names), function(i) {
          depth = node$depth + 1
          leafIndex = env$lastLeafIndex
          if (depth + 2 > dim(nodeDF)[2]) {
            # This is a leaf
            nchildren = 0
            nleaves = 1
          } else {
            childGroup = childrenGroups[[i]]
            groups = by(childGroup, childGroup[, depth + 1], list, simplify=F)
            nchildren = length(groups)
            nleaves = dim(childGroup)[1]
          }
          env$lastLeafIndex = leafIndex + nleaves
          child = EpivizNode$new(name=names[i], parentId=node$id, depth=depth, taxonomy=colnames(nodeDF)[depth + 1], nchildren=nchildren, nleaves=nleaves, order=i, leafIndex=leafIndex)
          if (nchildren == 0) { child$id = child$name }
          nodesById[[child$id]] <<- child
          if (length(nodesByDepthAndName) < depth + 1 || is.null(nodesByDepthAndName[[depth + 1]])) {
            nodesByDepthAndName[[depth + 1]] <<- list()
          }
          nodesByDepthAndName[[depth + 1]][[child$name]] <<- child
          return(child)
        })

        node$childrenIds = lapply(children, function(child) { child$id })
      }
      lapply(node$childrenIds, function(id) { nodesById[[id]] })
    },

    traverse=function(callback, node=root()) {
      doBreak = callback(node)
      if (doBreak) { return() }
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
#   setMethod(getGeneric("asJSON"), "MetavizTree", function(x, ...){
#     return(epivizr::toJSON(x$root(), ...))
#   })
# }

setGeneric("buildMetavizTree", signature=c("object"),
           function(object, ...) standardGeneric("buildMetavizTree"))

setMethod("buildMetavizTree", "MRexperiment", function(object, ...) {
  taxLvls = colnames(fData(object))[c(3:9,1)]
  tax = fData(object)[,taxLvls]
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

  return(MetavizTree$new(.tableToTree(tax)))
})
