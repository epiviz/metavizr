MetavizTree <- setRefClass("MetavizTree",
  fields=list(
    taxonomyTablePtr="Ptr",
    selectionTypes="Ptr",
    orders="Ptr",
    parents="Ptr",
    leavesCounts="Ptr",
    realLeavesCounts="Ptr",
    depthStrSize="numeric",
    leafIndexStrSize="numeric"
  ),
  methods=list(
    initialize=function(taxonomyTablePtr=Ptr$new(as.data.frame(NULL)), selectionTypes=Ptr$new(list()), orders=Ptr$new(list()),
                        parents=Ptr$new(list()), leavesCounts=Ptr$new(list()), realLeavesCounts=Ptr$new(list()), ...) {
      t = taxonomyTablePtr$.
      t[which(is.na(t), arr.ind=T)] = "NA"
      t = t[do.call(order, t),]
      taxonomyTablePtr <<- Ptr$new(t)
      selectionTypes <<- selectionTypes
      orders <<- orders
      parents <<- parents
      leavesCounts <<- leavesCounts
      realLeavesCounts <<- realLeavesCounts

      depthStrSize <<- ceiling(log(dim(taxonomyTablePtr$.)[2], base=16))
      leafIndexStrSize <<- ceiling(log(dim(taxonomyTablePtr$.)[1], base=16))
    },
    root=function() {
      return(MetavizNode$new(taxonomyTablePtr=taxonomyTablePtr, selectionTypes=selectionTypes, orders=orders, parents=parents,
                             leavesCounts=leavesCounts, realLeavesCounts=realLeavesCounts, depthStrSize=depthStrSize, leafIndexStrSize=leafIndexStrSize))
    },

#     children=function(node) {
#       if (is.null(node)) { return(NULL) }
#       if (node$nchildren > 0 && length(node$childrenIds) == 0) {
#         # Generate the children
#         nodeDF = taxonomyDF[node$leafIndex:(node$leafIndex + node$nleaves),]
#         childrenGroups = by(nodeDF, nodeDF[, node$depth + 2], list, simplify=F)
#         names = names(childrenGroups)
#         env = new.env()
#         env$lastLeafIndex = node$leafIndex
#         children = lapply(seq_along(names), function(i) {
#           depth = node$depth + 1
#           leafIndex = env$lastLeafIndex
#           if (depth + 2 > dim(nodeDF)[2]) {
#             # This is a leaf
#             nchildren = 0
#             nleaves = 1
#           } else {
#             childGroup = childrenGroups[[i]]
#             groups = by(childGroup, childGroup[, depth + 1], list, simplify=F)
#             nchildren = length(groups)
#             nleaves = dim(childGroup)[1]
#           }
#           env$lastLeafIndex = leafIndex + nleaves
#           child = EpivizNode$new(name=names[i], parentId=node$id, depth=depth, taxonomy=colnames(nodeDF)[depth + 1], nchildren=nchildren, nleaves=nleaves, order=i, leafIndex=leafIndex)
#           if (nchildren == 0) { child$id = child$name }
#           nodesById[[child$id]] <<- child
#           if (length(nodesByDepthAndName) < depth + 1 || is.null(nodesByDepthAndName[[depth + 1]])) {
#             nodesByDepthAndName[[depth + 1]] <<- list()
#           }
#           nodesByDepthAndName[[depth + 1]][[child$name]] <<- child
#           return(child)
#         })
#
#         node$childrenIds = lapply(children, function(child) { child$id })
#       }
#       lapply(node$childrenIds, function(id) { nodesById[[id]] })
#     },

#     traverse=function(callback, node=root()) {
#       doBreak = callback(node)
#       if (doBreak) { return() }
#       if (is.null(node)) { return() }
#       if (length(node$childrenIds) == 0) { return() }
#       for (child in children(node)) {
#         traverse(callback, child)
#       }
#     },

    node=function(id) {
      if (missing(id)) { return(root()) }
      info = .fromMetavizNodeId(id, depthStrSize=depthStrSize)
      return(MetavizNode$new(taxonomyTablePtr=taxonomyTablePtr, depth=info$depth, leafIndex=info$leafIndex,
                             selectionTypes=selectionTypes, orders=orders, parents=parents, leavesCounts=leavesCounts,
                             realLeavesCounts=realLeavesCounts, depthStrSize=depthStrSize, leafIndexStrSize=leafIndexStrSize))
    },

    parent=function(child=root()) {
      if (is.null(child)) { return(NULL) }

#       tryCatch({ child$parentId },
#         error=function(e) {
#           browser()
#         }
#       )
      if (!is.null(child$parentId)) {
        return(.self$node(child$parentId))
      }
      return(NULL)
    },

    levels=function() { colnames(taxonomyTablePtr$.) },

    # selection: @list {nodeId -> selectionType}
    updateSelection=function(selection) {
      if (missing(selection)) { return() }

      for (nodeId in names(selection)) {
        selectionTypes$.[[nodeId]] <<- selection[[nodeId]]

        realNleaves = realLeavesCounts$.[[nodeId]]
        nleaves = leavesCounts$.[[nodeId]]
        n = .self$node(nodeId)
        if (is.null(nleaves)) {
          nleaves = n$nleaves
        }
        if (is.null(realNleaves)) { realNleaves = nleaves }

        newNleaves = 0
        if (selection[[nodeId]] == SelectionType$NODE) { newNleaves = 1 }
        else if (selection[[nodeId]] == SelectionType$LEAVES) { newNleaves = nleaves }

        realLeavesCounts$.[[nodeId]] <<- newNleaves
        n = parent(n)
        while (!is.null(n)) {
          if (n$selectionType() != SelectionType$LEAVES) { break }
          l = realLeavesCounts$.[[n$id]]
          if (is.null(l)) { l = n$nleaves }
          realLeavesCounts$.[[n$id]] <<- l - realNleaves + newNleaves
          n = parent(n)
        }
      }
    },

    # order: @list {nodeId -> order index}
    updateOrder=function(order) {
      if (missing(order)) { return() }

      for (nodeId in names(order)) {
        orders$.[[nodeId]] <<- order[[nodeId]]
      }
    },

    # Get a list of the leaves of the given subtree, respecting the "selectionType" field
    # leafStartIndex is 0-based, and leafEndIndex is non-inclusive
    selectedLeaves=function(leafStartIndex, leafEndIndex) {
      if (leafStartIndex >= dim(taxonomyTablePtr$.)[1] || leafEndIndex <= 0) { return(NULL) }

      .iterate=function(n, s, realNodesBefore) {
        if (is.null(n)) { return(NULL) }
        if (n$selectionType() == SelectionType$NONE) { return(NULL) }
        if (s >= leafEndIndex || s + n$nleaves <= leafStartIndex) { return(NULL) }
        if (n$selectionType() == SelectionType$NODE || n$isLeaf()) {
          return(list(list(node=n, start=s, realNodesBefore=realNodesBefore)))
        }

        ret = list()
        children = n$children()
        for (child in children) {
          ret = c(ret, .iterate(child, s, realNodesBefore))
          s = s + child$nleaves
          realNodesBefore = realNodesBefore + child$realNleaves()
        }
        return(ret)
      }

      return(.iterate(root(), 0, 0))
    },

    # leafStartIndex is 0-based, and leafEndIndex is non-inclusive
#     .forSelectedLeaves(leafStartIndex, leafEndIndex, fun) {
#       .first=function() {
#         if (startIndex >= dim(taxonomyTablePtr$.)[1] || endIndex <= 0) { return(NULL) }
#         node = root()
#         if (node$selectionType() == SelectionType$NONE) { return(NULL) }
#
#         nodeStart = 0
#         nodeEnd = node$nleaves
#
#         while (TRUE) {
#           if (is.null(node)) { return(NULL) }
#           if (node$selectionType() == SelectionType$NODE) { return(node) }
#           if (node$isLeaf()) { return(node) }
#           nextNode = NULL
#           children = node$children()
#           s = nodeStart
#           for (child in children) {
#             if (s < endIndex && (s + child$nleaves > startIndex) && child$selectionType() != SelectionType$NONE) {
#               nextNode = child
#               nodeStart = s
#               nodeEnd = s + child$nleaves
#               break
#             }
#             s = s + child$nleaves
#           }
#           if ()
#           if (nextNode$selectionType() == SelectionType$NODE) { return(node) }
#           if (node$isLeaf()) { return(node) }
#
#           node = nextNode
#         }
#       }
#       .next=function(node) {
#
#       }
#
#       node = .first()
#       while (!is.null(node)) {
#         fun(node)
#         node = .next(node)
#       }
#     },
#
    # Gets a list of ancestor nodes in the tree for the given node
    ancestors=function(n=root(), inclusive=TRUE) {
      if (is.null(n)) { return(NULL) }
      ret = list()
      if (inclusive) {
        ret[[n$depth + 1]] = n
      }

      n = parent(n)
      while (!is.null(n)) {
        ret[[n$depth + 1]] = n
        n = parent(n)
      }
      return(ret)
    }
#
#
#     build=function(node=root(), filter=function(node) { TRUE }) {
#       recurse=function(node, filter) {
#         if (is.null(node)) { return(NULL) }
#         if (!filter(node)) { return(NULL) }
#         ret = node$raw()
#         if (length(node$childrenIds) == 0) { return(list(ret)) }
#
#         childrenList = list()
#         for (child in children(node)) {
#           childrenList = c(childrenList, recurse(child, filter))
#         }
#         if (length(childrenList) > 0) { ret$children = childrenList }
#         return(list(ret))
#       }
#       ret = recurse(node, filter)
#       if (is.null(ret)) { return(ret) }
#
#       return(ret[[1]])
#     }
  )
)

setGeneric("buildMetavizTree", signature=c("object"),
           function(object, ...) standardGeneric("buildMetavizTree"))

setMethod("buildMetavizTree", "MRexperiment", function(object, ...) {
  tax = colnames(fData(object))[c(3:9,1)]
  taxonomy = fData(object)[,tax]
  taxonomy[,1] = "Bacteria"
  #taxonomy = fData(object)[,tax]
  MetavizTree$new(Ptr$new(taxonomy))
})

#
# setMethod("buildMetavizTree", "MRexperiment", function(object, ...) {
#   taxLvls = colnames(fData(object))[c(3:9,1)]
#   tax = fData(object)[,taxLvls]
#   tax[,1] = "Bacteria"
#   table = tax
#
#   .tableToTree <- function(t, colIndex=1, parent=NULL) {
#     if (colIndex > dim(t)[2]) {
#       if (!is.null(parent)) { parent$nleaves = 1 }
#       return(NULL)
#     }
#
#     parentId = NULL
#     if (!is.null(parent)) { parentId = parent$id }
#
#     groups = by(t, t[, colIndex], list, simplify=F)
#     names = names(groups)
#     nodes = lapply(seq_along(names), function(i) { EpivizNode$new(name=names[i], parentId=parentId, depth=colIndex-1, taxonomy=colnames(t)[colIndex], order=i) })
#     if (!is.null(parent)) {
#       parent$nchildren = length(nodes)
#       parent$childrenIds = lapply(nodes, function(node) { node$id })
#     }
#
#     descendants = unlist(lapply(seq_along(groups), function(i) {
#       node = nodes[[i]]
#       children = .tableToTree(groups[[i]][[1]], colIndex+1, node)
#       if (!is.null(parent)) {
#         parent$nleaves = parent$nleaves + node$nleaves
#       }
#       children
#     }))
#     ret = c(nodes, descendants)
#
#     return(ret)
#   }
#
#   return(MetavizTree$new(.tableToTree(tax)))
# })
