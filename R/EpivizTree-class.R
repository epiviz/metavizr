EpivizTree <- setRefClass("EpivizTree",
  fields=list(
    rawTree="ANY",
    nodesById="ANY"
  ),
  methods=list(
    initialize=function(tree, ...) {
      rawTree <<- tree
      nodesById <<- new.env()
      .self$traverse(function(node) { nodesById[[node$id]] <<- node })
    },
    raw=function() { rawTree },
    leaves=function(node=rawTree) {
      if (node$nchildren == 0) { return(list(node)) }
      ret = list()
      for (i in 1:node$nchildren) {
        ret = c(ret, leaves(node$children[[i]]))
      }
      return(ret)
    },
    traverse=function(callback, node=rawTree) {
      callback(node)
      if (node$nchildren == 0) { return() }
      for (i in 1:node$nchildren) {
        traverse(callback, node$children[[i]])
      }
    },
    build=function(callback, node=rawTree) {
      ret = callback(node)
      if (is.null(ret) || length(ret$children) == 0) { return(ret) }
      children = c()
      n = length(ret$children)
      for (i in 1:n) {
        child = build(callback, ret$children[[i]])
        if (is.null(child)) { next }
        children = c(children, list(child))
      }
      ret$children = children
      ret
    },
    node=function(id=rawTree$id) {
      return(nodesById[[id]])
    },
    parent=function(...) {
      args = list(...)
      id = args$id
      node = args$node
      if (is.null(node)) {
        if (is.null(id)) { node = rawTree }
        else { node = .self$node(id) }
      }
      if (!is.null(node$parentId)) {
        return(.self$node(node$parentId))
      }
      return(NULL)
    }
  )
)
