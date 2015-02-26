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

# setMethod("toJSON", "Ptr", function(x, container = isContainer(x, asIs, .level),
#     collapse = "\n", ..., .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
#     .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
#   toJSON(x$obj)
# })

#library("jsonlite")
# setMethod(getGeneric("asJSON"), "Ptr", function(x, ...){ return(toJSON(x$., ...)) })

# setGeneric("toJSON", signature=c("x"),
#            function(x, ...) standardGeneric("toJSON"))
#
# setMethod(toJSON, "Ptr", function (x, ...) {
#   jsonlite::toJSON(x$obj)
# })


