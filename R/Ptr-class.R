Ptr <- setRefClass("Ptr",
  fields=list(
    .="ANY"
  ),
  methods=list(
    initialize=function(obj) {
      if (missing(obj)) { obj = NULL }
      . <<- obj
    }
  )
)

# setMethod(getGeneric("asJSON"), "Ptr", function(x, ...){ return(toJSON(x$., ...)) })
