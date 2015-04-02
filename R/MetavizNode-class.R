MetavizNode <- setRefClass("MetavizNode", # TODO: Rename to MRexperimentNode
  fields=list(
    .taxonomyTablePtr="Ptr",
    .depth="numeric",
    .leafIndex="numeric",     # the index of the first leaf in the taxonomy table
    .nleaves="ANY",
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
    .idMap="Ptr",

    .id="character",
    .parentId="ANY",
    .name="character"
  ),
  methods=list(
    initialize=function(taxonomyTablePtr=Ptr$new(as.data.frame(NULL)), depth=0, leafIndex=0, selectionTypes=Ptr$new(list()), orders=Ptr$new(list()),
                        parents=Ptr$new(list()), leavesCounts=Ptr$new(list()), realLeavesCounts=Ptr$new(list()), ancestryByDepth=Ptr$new(NULL), depthStrSize=0,
                        leafIndexStrSize=0, nodeMap=Ptr$new(list()), childrenMap=Ptr$new(list()), idMap=Ptr$new(NULL), ...) {
      .taxonomyTablePtr <<- taxonomyTablePtr
      .depth <<- depth
      .leafIndex <<- leafIndex
      .selectionTypes <<- selectionTypes
      .orders <<- orders
      .parents <<- parents
      .ancestryByDepth <<- ancestryByDepth
      .nodeMap <<- nodeMap
      .childrenMap <<- childrenMap
      .idMap <<- idMap

      .depthStrSize <<- depthStrSize
      if (.depthStrSize == 0) { .depthStrSize <<- ceiling(log(dim(.taxonomyTablePtr$.)[2], base=16)) }

      .leafIndexStrSize <<- leafIndexStrSize
      if (.leafIndexStrSize == 0) { .leafIndexStrSize <<- ceiling(log(dim(.taxonomyTablePtr$.)[1], base=16)) }

      .id <<- .generateMetavizNodeId(.depth, .leafIndex, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, idMap=.idMap)
      .name <<- as.character(.taxonomyTablePtr$.[.leafIndex+1, .depth+1])
      .realLeavesCounts <<- realLeavesCounts
      .leavesCounts <<- leavesCounts

      .parentId <<- NULL

      .nleaves <<- .leavesCounts$.[[.id]]
      if (is.null(.nleaves)) {
        if (.depth == 0) {
          .nleaves <<- nrow(.taxonomyTablePtr$.)
        } else if (.depth == ncol(.taxonomyTablePtr$.) - 1) {
          .nleaves <<- 1
        } else {
          .nleaves <<- nrow(.taxonomyTablePtr$.[.ancestryByDepth$.[[.depth]] == .ancestryByDepth$.[[.depth]][.leafIndex+1],])
        }
        .leavesCounts$.[[.id]] <<- .nleaves
      }
    },

    id=function() { .id },

    name=function() { .name },

    parentId=function() {
      if (!is.null(.parentId)) { return(.parentId) }
        .parentId <<- .parents$.[[.id]]
        if (is.null(.parentId) && .depth > 0) {
          if (.depth == 1) { .parentId <<- .generateMetavizNodeId(0, 0, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, idMap=.idMap)  }
          else {
            parentLeafIndex = min(which(.ancestryByDepth$.[[.depth-1]] == .ancestryByDepth$.[[.depth-1]][.leafIndex+1])) - 1
            .parentId <<- .generateMetavizNodeId(.depth-1, parentLeafIndex, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, idMap=.idMap)
          }
          .parents$.[[.id]] <<- .parentId
        }
      return(.parentId)
      #.parentId
    },

    depth=function() { .depth },

    leafIndex=function() { .leafIndex },

    nleaves=function() { .nleaves },

    isLeaf=function() {
      return(.depth == dim(.taxonomyTablePtr$.)[2] - 1)
    },

    selectionType=function() {
      selectionType = .selectionTypes$.[[.id]]
      if (!is.null(selectionType)) { return(selectionType) }

      return(SelectionType$LEAVES)
    },

    nodeOrder=function() {
      .orders$.[[.id]]
    },

    realNleaves=function() {
      ret = .realLeavesCounts$.[[.id]]
      if (is.null(ret)) { ret = .nleaves }
      return(ret)
    },

    children=function() {
      nlevels = ncol(.taxonomyTablePtr$.)
      if (.depth + 1 >= nlevels) { return(NULL) }

      ret = .childrenMap$.[[.id]]
      if (is.null(ret)) {

        if (.depth == nlevels - 2) {
          env = new.env()
          env$order = 0
          ret = lapply(seq(.leafIndex, .leafIndex + .nleaves - 1), function(childLeafIndex) {
            childId = .generateMetavizNodeId(.depth+1, childLeafIndex, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, idMap=.idMap)
            child = .nodeMap$.[[childId]]
            if (is.null(child)) {
              child = createNode(.depth+1, childLeafIndex)
              .nodeMap$.[[childId]] <<- child
            }
            if (is.null(.orders$.[[child$.id]])) {
              .orders$.[[child$.id]] <<- env$order
            }
            env$order = env$order + 1
            return(child)
          })
        } else {
          childrenAncestors = unique(.ancestryByDepth$.[[.depth + 1]][seq(.leafIndex + 1, .leafIndex + .nleaves)])
          env = new.env()
          env$order = 0
          ret = lapply(childrenAncestors, function(a) {
            childLeafIndex = min(which(.ancestryByDepth$.[[.depth + 1]] == a)) - 1
            childId = .generateMetavizNodeId(.depth+1, childLeafIndex, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, idMap=.idMap)
            child = .nodeMap$.[[childId]]
            if (is.null(child)) {
              child = createNode(.depth+1, childLeafIndex)
              .nodeMap$.[[childId]] <<- child
            }
            if (is.null(.orders$.[[child$.id]])) {
              .orders$.[[child$.id]] <<- env$order
            }
            env$order = env$order + 1
            return(child)
          })
        }

        .childrenMap$.[[.id]] <<- ret
      }

      o = sapply(ret, function(node) { node$nodeOrder() })
      ret = ret[order(o)]

      return(ret)
    },

    createNode=function(depth, leafIndex) {
      MetavizNode$new(taxonomyTablePtr=.taxonomyTablePtr, depth=depth, leafIndex=leafIndex, selectionTypes=.selectionTypes, orders=.orders,
                      parents=.parents, leavesCounts=.leavesCounts, realLeavesCounts=.realLeavesCounts, ancestryByDepth=.ancestryByDepth,
                      depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, nodeMap=.nodeMap, childrenMap=.childrenMap, idMap=.idMap)
    },

    raw = function(t=.taxonomyTablePtr$.[(.leafIndex+1):(.leafIndex+.nleaves),], maxDepth=NULL, colIndex=.depth+1, rowIndex=.leafIndex+1) {
      if ((!is.null(maxDepth) && maxDepth < 0) || colIndex > dim(t)[2]) { return(NULL) }

      nodeId = .generateMetavizNodeId(colIndex-1, rowIndex-1, depthStrSize=.depthStrSize, leafIndexStrSize=.leafIndexStrSize, idMap=.idMap)
      selectionType = SelectionType$LEAVES
      if (!is.null(.selectionTypes$.[[nodeId]])) { selectionType = .selectionTypes$.[[nodeId]] }

      ord = NULL
      if (!is.null(.orders$.[[nodeId]])) { ord = .orders$.[[nodeId]] }

      groups = NULL
      if (colIndex < dim(t)[2]) {
        groups = by(t, t[, colIndex + 1], list, simplify=F)
      }

      nodes = NULL
      if (maxDepth > 0) {
        childrenMaxDepth = NULL
        if (!is.null(maxDepth)) { childrenMaxDepth = maxDepth - 1 }
        env = new.env()
        env$order = 0
        nodes = unname(lapply(groups, function(group){
          childRowIndex = min(which(t[,colIndex+1]==group[[1]][1, colIndex+1], arr.ind=TRUE)) + rowIndex - 1
          child = raw(group[[1]], childrenMaxDepth, colIndex + 1, childRowIndex)
          if (is.null(child)) { return(NULL) }
          if (is.null(child$order)) {
            child$order = env$order
            .orders$.[[child$id]] <<- env$order
          }
          child$parentId = nodeId
          .parents$.[[child$id]] <<- nodeId

          env$order = env$order + 1
          return(child)
        }))
      }

      o = sapply(nodes, function(node) { node$order })
      nodes = nodes[order(o)]

      pId = .parents$.[[nodeId]]

      ret = list(
        name=t[1, colIndex],
        id=nodeId,
        globalDepth=colIndex-1,
        depth=colIndex-1,
        taxonomy=colnames(t)[colIndex],
        nchildren=length(groups),
        size=1,
        selectionType=selectionType,
        nleaves=dim(t)[1]
      )
      if (!is.null(ord)) { ret$order = ord }
      if (length(nodes) > 0) { ret$children = nodes }
      if (!is.null(pId)) { ret$parentId = pId }
      .leavesCounts$.[[nodeId]] <<- ret$nleaves
      return(ret)
    }
  )
)

.generateMetavizNodeId=function(depth, leafIndex, taxonomyTablePtr=NULL, depthStrSize=NULL, leafIndexStrSize=NULL, idMap=NULL) {
  if (!missing(idMap) && !is.null(idMap)) {
    id = idMap$.[leafIndex+1, depth+1]
    if (!is.na(id)) { return(id) }
  }
  if (missing(depthStrSize)) { depthStrSize = ceiling(log(dim(taxonomyTablePtr$.)[2], base=16)) }
  if (missing(leafIndexStrSize)) { leafIndexStrSize = ceiling(log(dim(taxonomyTablePtr$.)[1], base=16)) }

  id = paste(str_pad(as.hexmode(depth), depthStrSize, pad="0"), str_pad(as.hexmode(leafIndex), leafIndexStrSize, pad="0"), sep="-")
  if (!missing(idMap) && !is.null(idMap)) {
    idMap$.[leafIndex+1, depth+1] = id
  }
  return(id)
}

.fromMetavizNodeId=function(id, taxonomyTablePtr=NULL, depthStrSize=NULL, leafIndexStrSize=NULL) {
  if (missing(depthStrSize)) {
    if (missing(leafIndexStrSize)) {
      depthStrSize <<- ceiling(log(dim(taxonomyTablePtr$.)[2], base=16))
    } else {
      depthStrSize = nchar(id) - leafIndexStrSize - 1
    }
  }

  depth = as.integer(as.hexmode(substring(id, 1, depthStrSize)))
  leafIndex = as.integer(as.hexmode(substring(id, depthStrSize + 2)))
  return(list(depth=depth, leafIndex=leafIndex))
}
