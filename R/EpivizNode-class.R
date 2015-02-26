EpivizNode <- setRefClass("EpivizNode",
  fields=list(
    name="ANY",             # character
    id="character",
    parentId="ANY",         # character
    depth="numeric",
    taxonomy="ANY",         # character
    nchildren="numeric",
    selectionType="numeric",
    children="ANY",         # list of EpivizNode
    nleaves="numeric",
    order="numeric",

    childrenByName="list",
    childrenById="list"
  ),
  methods=list(
    initialize=function(id=.generatePseudoGUID(10), name=NULL, parentId=NULL, depth=0, taxonomy=NULL, nchildren=0, selectionType=1, children=NULL, nleaves=0, order=0, ...) {
      id <<- id
      name <<- name
      parentId <<- parentId
      depth <<- depth
      taxonomy <<- taxonomy
      nchildren <<- nchildren
      selectionType <<- selectionType
      children <<- children
      nleaves <<- nleaves
      order <<- order

      childrenByName <<- list()
      childrenById <<- list()

      if (length(children) > 0) {
        for (child in children) {
          childrenByName[[child$name]] <<- child
          childrenById[[child$id]] <<- child
        }
      }
    },

    childById=function(id) { childrenById[[id]] },
    childByName=function(name) { childrenByName[[name]] },

    copy=function(recursive=TRUE) {
      ret = EpivizNode$new(id=id, name=name, parentId=parentId, depth=depth, taxonomy=taxonomy, nchildren=nchildren, selectionType=selectionType, nleaves=nleaves, order=order)
      if (!recursive || length(children) == 0) { return(ret) }
      ret$children = list()
      for (child in children) {
        ret$children = c(ret$children, copy(child, recursive))
      }
      ret
    },

    raw=function(recursive=TRUE) {
      ret = list(
        name=name,
        id=id,
        parentId=parentId,
        depth=depth,
        globalDepth=depth,
        taxonomy=taxonomy,
        nchildren=nchildren,
        size=1,
        selectionType=selectionType,
        children=NULL,
        nleaves=nleaves,
        order=order
      )
      if (length(children) == 0 || !recursive) { return(ret) }
      for (child in children) {
        ret$children = c(ret$children, list(child$raw(recursive)))
      }
      return(ret)
    }
  )
)

# Used by jsonlite::toJSON
# if (!is.null(getGeneric("asJSON"))) {
#   setMethod(getGeneric("asJSON"), "EpivizNode", function(x, ...){
#     return(epivizr::toJSON(x=x$raw(), ...))
#   })
# }
