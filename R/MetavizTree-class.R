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
    .leafIndexStrSize="numeric"
  ),
  methods=list(
    initialize=function(taxonomyTablePtr=Ptr$new(as.data.frame(NULL)), selectionTypes=Ptr$new(list()), orders=Ptr$new(list()),
                        parents=Ptr$new(list()), leavesCounts=Ptr$new(list()), realLeavesCounts=Ptr$new(list()), ...) {
      t = taxonomyTablePtr$.
      t[which(is.na(t), arr.ind=T)] = "NA"
      t = t[do.call(order, t),]
      .taxonomyTablePtr <<- Ptr$new(t)

      .ancestryByDepth <<- Ptr$new(NULL)
      if (nrow(t) > 2) {
        .ancestryByDepth$. <<- lapply(2:(ncol(t)-1), function(depth) { do.call(paste, c(t[,1:depth],sep=",")) })
      }

      .selectionTypes <<- selectionTypes
      .orders <<- orders
      .parents <<- parents
      .leavesCounts <<- leavesCounts
      .realLeavesCounts <<- realLeavesCounts

      .depthStrSize <<- ceiling(log(dim(.taxonomyTablePtr$.)[2], base=16))
      .leafIndexStrSize <<- ceiling(log(dim(.taxonomyTablePtr$.)[1], base=16))
    },

    .createNode=function(depth=0, leafIndex=0) {
      return(MetavizNode$new(taxonomyTablePtr=.taxonomyTablePtr, depth=depth, leafIndex=leafIndex,
                             selectionTypes=.selectionTypes, orders=.orders, parents=.parents, leavesCounts=.leavesCounts,
                             realLeavesCounts=.realLeavesCounts, ancestryByDepth=.ancestryByDepth,
                             depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize))
    },

    root=function() {
      return(.createNode())
    },

    node=function(id) {
      if (missing(id)) { return(root()) }
      info = .fromMetavizNodeId(id, depthStrSize=.depthStrSize)
      return(.createNode(info$depth, info$leafIndex))
    },

    parent=function(child=root()) {
      if (is.null(child)) { return(NULL) }
      if (!is.null(child$parentId())) {
        return(.self$node(child$parentId()))
      }
      return(NULL)
    },

    nodesAtDepth=function(depth=0) {
      if (depth == 0) { return(list(root())) }
      if (depth == ncol(.taxonomyTablePtr$.) - 1) {
        return(lapply(1:nrow(.taxonomyTablePtr$.), function(i) {
          return(.createNode(depth, i-1))
        }))
      }
      if (depth > ncol(.taxonomyTablePtr$.) - 1) {
        return(NULL)
      }

      groups = by(.taxonomyTablePtr$., .ancestryByDepth$.[[depth]], list, simplify=F)
      env = new.env()
      env$leafIndex = 0
      return(lapply(groups, function(group) {
        leafIndex = env$leafIndex
        nleaves = nrow(group[[1]])
        childId = .generateMetavizNodeId(depth, leafIndex, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize)

        .leavesCounts$.[[childId]] <<- nleaves
        env$leafIndex = env$leafIndex + nleaves

        return(.createNode(depth, leafIndex))
      }))
    },

    levels=function() { colnames(.taxonomyTablePtr$.) },

    # selection: @list {nodeId -> selectionType}
    updateSelection=function(selection) {
      if (missing(selection)) { return() }

      for (nodeId in names(selection)) {
        .selectionTypes$.[[nodeId]] <<- selection[[nodeId]]

        n = .self$node(nodeId)
        realNleaves = n$realNleaves()
        nleaves = n$nleaves()

        newNleaves = 0
        if (selection[[nodeId]] == SelectionType$NODE) { newNleaves = 1 }
        else if (selection[[nodeId]] == SelectionType$LEAVES) { newNleaves = nleaves }

        .realLeavesCounts$.[[nodeId]] <<- newNleaves
        n = parent(n)
        while (!is.null(n)) {
          if (n$selectionType() != SelectionType$LEAVES) { break }
          l = n$realNleaves()
          .realLeavesCounts$.[[n$id()]] <<- l - realNleaves + newNleaves
          n = parent(n)
        }
      }
    },

    # order: @list {nodeId -> order index}
    updateOrder=function(order) {
      if (missing(order)) { return() }

      for (nodeId in names(order)) {
        .orders$.[[nodeId]] <<- order[[nodeId]]
      }
    },

    # Get a list of the leaves of the given subtree, respecting the "selectionType" field
    # leafStartIndex is 0-based, and leafEndIndex is non-inclusive
    selectedLeaves=function(leafStartIndex, leafEndIndex) {
      if (leafStartIndex >= dim(.taxonomyTablePtr$.)[1] || leafEndIndex <= 0) { return(NULL) }

      .iterate=function(n, s, realNodesBefore) {
        if (is.null(n)) { return(NULL) }
        if (n$selectionType() == SelectionType$NONE) { return(NULL) }
        if (s >= leafEndIndex || s + n$nleaves() <= leafStartIndex) { return(NULL) }
        if (n$selectionType() == SelectionType$NODE || n$isLeaf()) {
          return(list(list(node=n, start=s, realNodesBefore=realNodesBefore)))
        }

        ret = list()
        children = n$children()
        for (child in children) {
          ret = c(ret, .iterate(child, s, realNodesBefore))
          s = s + child$nleaves()
          realNodesBefore = realNodesBefore + child$realNleaves()
        }
        return(ret)
      }

      return(.iterate(root(), 0, 0))
    },

    # Gets a list of ancestor nodes in the tree for the given node
    ancestors=function(n=root(), inclusive=TRUE) {
      if (is.null(n)) { return(NULL) }
      ret = list()
      if (inclusive) {
        ret[[n$depth() + 1]] = n
      }

      n = parent(n)
      while (!is.null(n)) {
        ret[[n$depth() + 1]] = n
        n = parent(n)
      }
      return(ret)
    }
  )
)

setGeneric("buildMetavizTree", signature=c("object"),
           function(object, ...) standardGeneric("buildMetavizTree"))

setMethod("buildMetavizTree", "MRexperiment", function(object, ...) {
  # TODO Joe
  tax = colnames(fData(object))[c(3:9,1)]
  taxonomy = fData(object)[,tax]
  taxonomy[,1] = "Bacteria"
#   taxonomy = fData(object)
  MetavizTree$new(Ptr$new(taxonomy))
})
