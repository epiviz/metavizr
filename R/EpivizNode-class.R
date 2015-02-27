EpivizNode <- setRefClass("EpivizNode",
  fields=list(
    name="ANY",             # character
    id="character",
    parentId="ANY",         # character
    depth="numeric",
    taxonomy="ANY",         # character
    nchildren="numeric",
    selectionType="numeric",
    childrenIds="ANY",      # list of ids
    nleaves="numeric",
    order="numeric"
  ),
  methods=list(
    initialize=function(id=.generatePseudoGUID(10), name=NULL, parentId=NULL, depth=0, taxonomy=NULL, nchildren=0, selectionType=1, childrenIds=NULL, nleaves=0, order=0, ...) {
      id <<- id
      name <<- name
      parentId <<- parentId
      depth <<- depth
      taxonomy <<- taxonomy
      nchildren <<- nchildren
      selectionType <<- selectionType
      childrenIds <<- childrenIds
      nleaves <<- nleaves
      order <<- order

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
        nleaves=nleaves,
        order=order
      )
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
