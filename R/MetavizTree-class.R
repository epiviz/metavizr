#' Taxonomy tree structure wrapper
#'
#' Used to manage tree-cut queries from the
#' Metaviz app UI. 
#' 
#' 
MetavizTree <- setRefClass("MetavizTree",
  fields=list(
    .taxonomyTablePtr="Ptr",
    .selectionTypes="Ptr",
    .orders="Ptr",
    .parents="Ptr",
    .leavesCounts="Ptr",
    .realLeavesCounts="Ptr",
    .ancestryByDepth="Ptr",
    .depthStrSize="numeric",
    .leafIndexStrSize="numeric",
    .nodeMap="Ptr",
    .childrenMap="Ptr",
    .idMap="Ptr"
  ),
  methods=list(
    initialize=function(taxonomyTablePtr=Ptr$new(as.data.frame(NULL)), 
                        selectionTypes=Ptr$new(list()), 
                        orders=Ptr$new(list()),
                        parents=Ptr$new(list()), 
                        leavesCounts=Ptr$new(list()), 
                        realLeavesCounts=Ptr$new(list()), ...) {
      t <- taxonomyTablePtr$.
      t[which(is.na(t), arr.ind=TRUE)] = "NA"
      for(i in seq_along(t)) { t[,i] = as.character(t[,i]) }
      t <- t[do.call(order, t),]
      .self$.taxonomyTablePtr <- Ptr$new(t)

      .self$.ancestryByDepth <- Ptr$new(NULL)
      if (nrow(t) > 2) {
        .self$.ancestryByDepth$. <- lapply(2:(ncol(t)), function(depth) { 
          do.call(paste, c(t[,1:depth], sep=",")) 
        })
      }

      .self$.selectionTypes <- selectionTypes
      .self$.orders <- orders
      .self$.parents <- parents
      .self$.leavesCounts <- leavesCounts
      .self$.realLeavesCounts <- realLeavesCounts

      .self$.depthStrSize <- ceiling(log(dim(.self$.taxonomyTablePtr$.)[2], base=16))
      .self$.leafIndexStrSize <- ceiling(log(dim(.self$.taxonomyTablePtr$.)[1], base=16))
      .self$.nodeMap <- Ptr$new(new.env())
      .self$.childrenMap <- Ptr$new(new.env())
      .self$.idMap <- Ptr$new(matrix(nrow=nrow(t), ncol=ncol(t)))
    },

    .createNode=function(depth=0, leafIndex=0) {
      MetavizNode$new(taxonomyTablePtr=.self$.taxonomyTablePtr, 
                      depth=depth, 
                      leafIndex=leafIndex,
                      selectionTypes=.self$.selectionTypes, 
                      orders=.self$.orders, 
                      parents=.self$.parents, 
                      leavesCounts=.self$.leavesCounts,
                      realLeavesCounts=.self$.realLeavesCounts, 
                      ancestryByDepth=.self$.ancestryByDepth,
                      depthStrSize=.self$.depthStrSize, 
                      leafIndexStrSize=.self$.leafIndexStrSize, 
                      nodeMap=.self$.nodeMap, 
                      childrenMap=.self$.childrenMap, 
                      idMap=.self$.idMap)
    },

    taxonomyTable=function() { .self$.taxonomyTablePtr$. },

    root=function() {
      return(.self$.createNode())
    },

    node=function(id) {
      if (missing(id)) { return(root()) }
      ret <- .self$.nodeMap$.[[id]]
      if (!is.null(ret)) { return(ret) }
      info <- .fromMetavizNodeId(id, depthStrSize=.self$.depthStrSize)
      ret <- .createNode(info$depth, info$leafIndex)
      .self$.nodeMap$.[[id]] <- ret
      return(ret)
    },

    parent=function(child=root()) {
      if (is.null(child)) { return(NULL) }
      if (!is.null(child$parentId())) {
        return(.self$node(child$parentId()))
      }
      NULL
    },

    nodesAtDepth=function(depth=0) {
      if (depth == 0) { return(list(root())) }
      if (depth == ncol(.self$.taxonomyTablePtr$.) - 1) {
        return(lapply(1:nrow(.self$.taxonomyTablePtr$.), function(i) {
          return(.self$.createNode(depth, i-1))
        }))
      }
      if (depth > ncol(.self$.taxonomyTablePtr$.) - 1) {
        return(NULL)
      }

      groups <- by(.self$.taxonomyTablePtr$., 
                   .self$.ancestryByDepth$.[[depth]], 
                   list, simplify=FALSE)
      env = new.env()
      env$leafIndex <- 0
      res <- lapply(groups, function(group) {
        leafIndex <- env$leafIndex
        nleaves <- nrow(group[[1]])
        childId <- .generateMetavizNodeId(depth, leafIndex, 
                                          depthStrSize=.self$.depthStrSize, 
                                          leafIndexStrSize=.self$.leafIndexStrSize, 
                                          idMap=.self$.idMap)
        
        .self$.leavesCounts$.[[childId]] <- nleaves
        env$leafIndex <- env$leafIndex + nleaves
        
        .self$.createNode(depth, leafIndex)
      })
      res
    },

    calcNodeId=function(rowIndex, colIndex) {
      depth <- colIndex - 1
      leafIndex <- -1
      if (depth == 0) { 
        res <- .generateMetavizNodeId(depth, 0, 
                                      depthStrSize=.self$.depthStrSize, 
                                      leafIndexStrSize=.self$.leafIndexStrSize, 
                                      idMap=.self$.idMap)
        return(res) 
      }
      if (depth == 1) {
        t <- .self$.taxonomyTablePtr
        leafIndex = min(which(t$.[,2] == t$.[rowIndex, colIndex])) - 1
      } else if (depth == ncol(.self$.taxonomyTablePtr$.) - 1) {
        leafIndex <- rowIndex - 1
      } else {
        is_candidate <- .ancestryByDepth$.[[depth]] == 
                        .ancestryByDepth$.[[depth]][rowIndex]
        leafIndex <- min(which(is_candidate)) - 1
      }
      id <- .generateMetavizNodeId(depth, leafIndex, 
                                   depthStrSize=.self$.depthStrSize, 
                                   leafIndexStrSize=.self$.leafIndexStrSize, 
                                   idMap=.self$.idMap)
      id
    },

    siblings=function(node) {
      ret <- list()
      p <- .self$parent(node)
      if (is.null(p)) { return(ret) }
      for (child in p$children()) {
        if (child$id() != node$id()) {
          ret <- c(ret, child)
        }
      }
      ret
    },

    levels=function() { colnames(.self$.taxonomyTablePtr$.) },

    # selection: @list {nodeId -> selectionType}
    updateSelection=function(selection) {
      if (missing(selection)) { return(NULL) }

      for (nodeId in names(selection)) {
        .self$.selectionTypes$.[[nodeId]] <- selection[[nodeId]]

        n <- .self$node(nodeId)
        realNleaves <- n$realNleaves()
        nleaves <- n$nleaves()

        newNleaves = 0
        if (selection[[nodeId]] == SelectionType$NODE) { 
          newNleaves <- 1 
        } else if (selection[[nodeId]] == SelectionType$LEAVES) { 
          newNleaves <- nleaves 
        }

        .self$.realLeavesCounts$.[[nodeId]] <- newNleaves
        n <- .self$parent(n)
        while (!is.null(n)) {
          if (n$selectionType() != SelectionType$LEAVES) { break }
          l <- n$realNleaves()
          .self$.realLeavesCounts$.[[n$id()]] <- l - realNleaves + newNleaves
          n <- parent(n)
        }
      }
    },

    # order: @list {nodeId -> order index}
    updateOrder=function(order) {
      if (missing(order)) { return(NULL) }

      for (nodeId in names(order)) {
        .self$.orders$.[[nodeId]] <- order[[nodeId]]
      }
    },

    # Get a list of the leaves of the given subtree, respecting the "selectionType" field
    # leafStartIndex is 0-based, and leafEndIndex is non-inclusive
    selectedLeaves=function(leafStartIndex, leafEndIndex) {
      if (leafStartIndex >= dim(.taxonomyTablePtr$.)[1] || leafEndIndex <= 0) { 
        return(NULL) 
      }

      .iterate=function(node, s, realNodesBefore) {
        if (is.null(node)) { return(NULL) }
        if (node$selectionType() == SelectionType$NONE) { return(NULL) }
        if (s >= leafEndIndex || s + node$nleaves() <= leafStartIndex) { 
          return(NULL) 
        }
        if (node$selectionType() == SelectionType$NODE || node$isLeaf()) {
          return(list(list(node=node, start=s, realNodesBefore=realNodesBefore)))
        }

        ret <- list()
        children <- node$children()
        for (child in children) {
          ret <- c(ret, .iterate(child, s, realNodesBefore))
          s <- s + child$nleaves()
          realNodesBefore <- realNodesBefore + child$realNleaves()
        }
        ret
      }

      .iterate(root(), 0, 0)
    },

    # Gets a list of ancestor nodes in the tree for the given node
    ancestors=function(n=.self$root(), inclusive=TRUE) {
      if (is.null(n)) { return(NULL) }
      ret = list()
      if (inclusive) {
        ret[[n$depth() + 1]] <- n
      }

      n <- .self$parent(n)
      while (!is.null(n)) {
        ret[[n$depth() + 1]] <- n
        n <- .self$parent(n)
      }
      ret
    }
  )
)

#' Build a MetavizTree object from another object
#' 
#' @param object The object from which taxonomy data is extracted
#' @param ... Additional arguments
#' @return a \code{\link{MetavizTree}} object
setGeneric("buildMetavizTree", signature=c("object"),
           function(object, ...) standardGeneric("buildMetavizTree"))

#' @describeIn buildMetavizTree Build tree from a \code{\link[metagenomeSeq]{MRexperiment-class}} object
#' @importFrom Biobase fData
#' @param feature_order Ordering of leaves (features) in taxonomy tree
setMethod("buildMetavizTree", "MRexperiment", function(object, feature_order, ...) {
  taxonomy <- Biobase::fData(object)
  if(!is.null(feature_order)) {
    fOrder <- match(feature_order, colnames(taxonomy))
    taxonomy <- taxonomy[fOrder]
  }
  MetavizTree$new(Ptr$new(taxonomy))
})
