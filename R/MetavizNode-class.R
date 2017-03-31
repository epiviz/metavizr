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
    initialize=function(taxonomyTablePtr=Ptr$new(as.data.frame(NULL)), 
                        depth=0, 
                        leafIndex=0, 
                        selectionTypes=Ptr$new(list()), 
                        orders=Ptr$new(list()),
                        parents=Ptr$new(list()), 
                        leavesCounts=Ptr$new(list()), 
                        realLeavesCounts=Ptr$new(list()), 
                        ancestryByDepth=Ptr$new(NULL), 
                        depthStrSize=0,
                        leafIndexStrSize=0, 
                        nodeMap=Ptr$new(list()), 
                        childrenMap=Ptr$new(list()), 
                        idMap=Ptr$new(NULL), 
                        ...) {
      .self$.taxonomyTablePtr <- taxonomyTablePtr
      .self$.depth <- depth
      .self$.leafIndex <- leafIndex
      .self$.selectionTypes <- selectionTypes
      .self$.orders <- orders
      .self$.parents <- parents
      .self$.ancestryByDepth <- ancestryByDepth
      .self$.nodeMap <- nodeMap
      .self$.childrenMap <- childrenMap
      .self$.idMap <- idMap

      .self$.depthStrSize <- depthStrSize
      if (.self$.depthStrSize == 0) { 
        .self$.depthStrSize <- ceiling(log(dim(.self$.taxonomyTablePtr$.)[2], base=16)) 
      }

      .self$.leafIndexStrSize <- leafIndexStrSize
      if (.self$.leafIndexStrSize == 0) { 
        .self$.leafIndexStrSize <- ceiling(log(dim(.taxonomyTablePtr$.)[1], base=16)) 
      }

      .self$.id <- .generateMetavizNodeId(.self$.depth, 
                      .self$.leafIndex, 
                      depthStrSize=.self$.depthStrSize, 
                      leafIndexStrSize=.self$.leafIndexStrSize, 
                      idMap=.self$.idMap)
      
      .self$.name <- as.character(.self$.taxonomyTablePtr$.[.leafIndex+1, .depth+1])
      .self$.realLeavesCounts <- realLeavesCounts
      .self$.leavesCounts <- leavesCounts

      .self$.parentId <- NULL

      .self$.nleaves <- .self$.leavesCounts$.[[.self$.id]]
      if (is.null(.self$.nleaves)) {
        if (.self$.depth == 0) {
          .self$.nleaves <- nrow(.self$.taxonomyTablePtr$.)
        } else if (.self$.depth == ncol(.self$.taxonomyTablePtr$.) - 1) {
          .self$.nleaves <- 1
        } else {
          row_index <- .self$.ancestryByDepth$.[[.self$.depth]] == 
                       .self$.ancestryByDepth$.[[.self$.depth]][.self$.leafIndex+1]
          .self$.nleaves <- nrow(.self$.taxonomyTablePtr$.[row_index,])
        }
        .self$.leavesCounts$.[[.id]] <- .self$.nleaves
      }
    },

    id=function() { .self$.id },

    name=function() { .self$.name },

    parentId=function() {
      if (!is.null(.self$.parentId)) { return(.self$.parentId) }
      
      .self$.parentId <- .self$.parents$.[[.self$.id]]
      if (is.null(.self$.parentId) && .self$.depth > 0) {
        if (.self$.depth == 1) { 
          .self$.parentId <- .generateMetavizNodeId(0, 0, 
                                                    depthStrSize=.self$.depthStrSize, 
                                                    leafIndexStrSize=.self$.leafIndexStrSize, 
                                                    idMap=.self$.idMap)  }
        else {
          is_possible_index <- .self$.ancestryByDepth$.[[.self$.depth-1]] == 
            .self$.ancestryByDepth$.[[.self$.depth-1]][.self$.leafIndex+1]
          parentLeafIndex <- min(which(is_possible_index)) - 1
          .self$.parentId <- .generateMetavizNodeId(.self$.depth-1, 
                                                    parentLeafIndex, 
                                                    depthStrSize=.self$.depthStrSize, 
                                                    leafIndexStrSize=.self$.leafIndexStrSize, 
                                                    idMap=.idMap)
        }
        .self$.parents$.[[.self$.id]] <- .self$.parentId
      }
      .self$.parentId
    },

    depth=function() { .self$.depth },

    leafIndex=function() { .self$.leafIndex },

    nleaves=function() { .self$.nleaves },

    isLeaf=function() {
      .self$.depth == (dim(.self$.taxonomyTablePtr$.)[2] - 1)
    },

    selectionType=function() {
      selectionType <- .self$.selectionTypes$.[[.self$.id]]
      if (!is.null(selectionType)) { return(selectionType) }

      SelectionType$LEAVES
    },

    nodeOrder=function() {
      .self$.orders$.[[.self$.id]]
    },

    realNleaves=function() {
      ret <- .self$.realLeavesCounts$.[[.self$.id]]
      if (is.null(ret)) { ret <- .self$.nleaves }
      ret
    },

    children=function() {
      nlevels <- ncol(.self$.taxonomyTablePtr$.)
      if (.self$.depth + 1 >= nlevels) { return(NULL) }

      ret = .self$.childrenMap$.[[.self$.id]]
      if (is.null(ret)) {

        if (.self$.depth == nlevels - 2) {
          env = new.env()
          env$order = 0
          ret = lapply(seq(.leafIndex, .leafIndex + .nleaves - 1), function(childLeafIndex) {
            childId <- .generateMetavizNodeId(.self$.depth+1, 
                                              childLeafIndex, 
                                              depthStrSize=.self$.depthStrSize, 
                                              leafIndexStrSize=.self$.leafIndexStrSize, 
                                              idMap=.self$.idMap)
            child <- .self$.nodeMap$.[[childId]]
            if (is.null(child)) {
              child <- createNode(.self$.depth+1, childLeafIndex)
              .self$.nodeMap$.[[childId]] <- child
            }
            if (is.null(.self$.orders$.[[child$.id]])) {
              .self$.orders$.[[child$.id]] <- env$order
            }
            env$order <- env$order + 1
            return(child)
          })
        } else {
          indices <- seq(.self$.leafIndex + 1, .self$.leafIndex + .self$.nleaves)
          childrenAncestors <- unique(.self$.ancestryByDepth$.[[.self$.depth + 1]][indices])
          env = new.env()
          env$order <- 0
          ret = lapply(childrenAncestors, function(a) {
            childLeafIndex <- min(which(.self$.ancestryByDepth$.[[.self$.depth + 1]] == a)) - 1
            childId <- .generateMetavizNodeId(.self$.depth+1, 
                                              childLeafIndex, 
                                              depthStrSize=.self$.depthStrSize, 
                                              leafIndexStrSize=.self$.leafIndexStrSize, 
                                              idMap=.self$.idMap)
            child <- .self$.nodeMap$.[[childId]]
            if (is.null(child)) {
              child <- createNode(.self$.depth+1, childLeafIndex)
              .self$.nodeMap$.[[childId]] <- child
            }
            
            if (is.null(.self$.orders$.[[child$.id]])) {
              .self$.orders$.[[child$.id]] <- env$order
            }
            
            env$order <- env$order + 1
            child
          })
        }

        .self$.childrenMap$.[[.self$.id]] <- ret
      }

      o <- sapply(ret, function(node) { node$nodeOrder() })
      ret <- ret[order(o)]
      ret
    },

    createNode=function(depth, leafIndex) {
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

    raw = function(t=NULL, 
                   maxDepth=NULL, 
                   colIndex=.self$.depth+1, 
                   rowIndex=.self$.leafIndex+1) {
      if (missing(t) || is.null(t)) {
        t <- .self$.taxonomyTablePtr$.[(.self$.leafIndex+1):(.self$.leafIndex+.nleaves),]
      }      
      
      if ((!is.null(maxDepth) && maxDepth < 0) || colIndex > ncol(t)) { return(NULL) }

      nodeId <- .generateMetavizNodeId(colIndex-1, rowIndex-1, 
                                       depthStrSize=.self$.depthStrSize, 
                                       leafIndexStrSize=.self$.leafIndexStrSize, 
                                       idMap=.self$.idMap)
      selectionType = SelectionType$LEAVES
      if (!is.null(.self$.selectionTypes$.[[nodeId]])) { 
        selectionType = .self$.selectionTypes$.[[nodeId]] 
      }

      ord <- NULL
      if (!is.null(.self$.orders$.[[nodeId]])) { ord <- .self$.orders$.[[nodeId]] }

      groups <- NULL
      if (colIndex < ncol(t)) {
        groups <- by(t, t[, colIndex + 1], list, simplify=FALSE)
      }

      nodes <- NULL
      if (maxDepth > 0) {
        childrenMaxDepth <- NULL
        if (!is.null(maxDepth)) { childrenMaxDepth = maxDepth - 1 }
        env = new.env()
        env$order <- 0
        nodes <- unname(lapply(groups, function(group) {
          childRowIndex <- min(which(t[,colIndex+1]==group[[1]][1, colIndex+1], arr.ind=TRUE)) + rowIndex - 1
          child <- raw(group[[1]], childrenMaxDepth, colIndex + 1, childRowIndex)
          if (is.null(child)) { return(NULL) }
          if (is.null(child$order)) {
            child$order <- env$order
            .self$.orders$.[[child$id]] <- env$order
          }
          child$parentId <- nodeId
          .self$.parents$.[[child$id]] <- nodeId

          env$order <- env$order + 1
          return(child)
        }))
      }

      o <- sapply(nodes, function(node) { node$order })
      nodes <- nodes[order(o)]

      pId <- .self$.parents$.[[nodeId]]

      ret = list(
        name=t[1, colIndex],
        id=nodeId,
        globalDepth=colIndex-1,
        depth=colIndex-1,
        taxonomy=colnames(t)[colIndex],
        nchildren=length(groups),
        size=1,
        selectionType=selectionType,
        nleaves=nrow(t),
        start=rowIndex,
        end=rowIndex+nrow(t)
      )
      if (!is.null(ord)) { ret$order <- ord }
      if (length(nodes) > 0) { ret$children <- nodes }
      if (!is.null(pId)) { ret$parentId <- pId }
      .self$.leavesCounts$.[[nodeId]] <- ret$nleaves
      ret
    }
  )
)

.generateMetavizNodeId=function(depth, leafIndex,
                                taxonomyTablePtr=NULL, 
                                depthStrSize=NULL, 
                                leafIndexStrSize=NULL, 
                                idMap=NULL) {
  if (!missing(idMap) && !is.null(idMap)) {
    id <- idMap$.[leafIndex+1, depth+1]
    if (!is.na(id)) { return(id) }
  }
#   if (missing(depthStrSize)) { depthStrSize = ceiling(log(dim(taxonomyTablePtr$.)[2], base=16)) }
#   if (missing(leafIndexStrSize)) { leafIndexStrSize = ceiling(log(dim(taxonomyTablePtr$.)[1], base=16)) }
#
#   id = paste(str_pad(as.hexmode(depth), depthStrSize, pad="0"), str_pad(as.hexmode(leafIndex), leafIndexStrSize, pad="0"), sep="-")
#   if (!missing(idMap) && !is.null(idMap)) {
#     idMap$.[leafIndex+1, depth+1] = id
#   }
  id <- paste(as.hexmode(depth), as.hexmode(leafIndex), sep="-")
  id
}

.fromMetavizNodeId=function(id, taxonomyTablePtr=NULL, 
                            depthStrSize=NULL, 
                            leafIndexStrSize=NULL) {
  tokens <- strsplit(id, "-")
  depthStr <- tokens[[1]][1]
  leafIndexStr <- tokens[[1]][2]
  depth <- as.integer(as.hexmode(depthStr))
  leafIndex <- as.integer(as.hexmode(leafIndexStr))
  return(list(depth=depth, leafIndex=leafIndex))
#   if (missing(depthStrSize)) {
#     if (missing(leafIndexStrSize)) {
#       depthStrSize <<- ceiling(log(dim(taxonomyTablePtr$.)[2], base=16))
#     } else {
#       depthStrSize = nchar(id) - leafIndexStrSize - 1
#     }
#   }
#
#   depth = as.integer(as.hexmode(substring(id, 1, depthStrSize)))
#   leafIndex = as.integer(as.hexmode(substring(id, depthStrSize + 2)))
#   return(list(depth=depth, leafIndex=leafIndex))
}
