EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains="EpivizData",
  fields=list(taxonomy="ANY", leaves="list", samples="list"),
  methods=list(
    initialize=function(object, ...) {
      taxonomy <<- .MRexperimentToTree(object)
      leaves <<- .leaves(taxonomy)
      samples <<- list(list(name="sample1", id="s1"), list(name="sample2", id="s2")) # TODO
      callSuper(object=object, ...)
    },
    update=function(newObject, ...) {
# TODO
#       if (!is(newObject, "SummarizedExperiment"))
#         stop("'newObject' must be of class 'SummarizedExperiment'")
#
#       newObject <- reorderIfNecessary(newObject)
#
#       if(!is(rowData(newObject), "GIntervalTree"))
#         rowData(newObject) <- as(rowData(newObject), "GIntervalTree")

      callSuper(newObject, ...)
    },
    .MRexperimentToTree=function(exp) {
      filteredExp = filterData(exp,depth=6000,present=25)
      levels = colnames(fData(filteredExp))[c(3:9,1)]
      tax = fData(filteredExp)[,levels]
      tax[,1] = "Bacteria"
      table = tax[sample(1:nrow(tax),10),]

      return(.tableToTree(table))
    },
    .tableToTree=function(t) {

      tableToListTree <- function(t, colIndex=1) {
        if (colIndex > dim(t)[2]) { return(NULL) }
        groups = by(t, t[, colIndex], list, simplify=F)
        nodes = lapply(groups, function(group){
          tableToListTree(group[[1]], colIndex + 1)
        })
        return(nodes)
      }
      listTreeToJSONTree <- function(node, globalDepth=0) {
        nodeId = .generatePseudoGUID(6)
        ret = list(
          name=names(node),
          id=nodeId,
          parentId=NULL, # TODO
          depth=globalDepth, # TODO: This is the relative depth of the node (? Check in JS code)
          globalDepth=globalDepth,
          taxonomy=colnames(t)[globalDepth+1],
          nchildren=length(node[[1]]),
          size=1,
          selectionType=NULL # TODO
        )
        nleaves = 0
        if (length(node[[1]]) > 0) {
          children = c()
          for (i in 1:length(node[[1]])) {
            child = listTreeToJSONTree(node[[1]][i], globalDepth + 1)
            child$order = i - 1
            child$parentId = nodeId
            children = c(children, list(child))
            nleaves = nleaves + child$nleaves
          }
          ret$children = children
        } else {
          nleaves = 1
        }
        ret$nleaves = nleaves
        return(ret)
      }

      listTree = tableToListTree(t)
      tree = listTreeToJSONTree(listTree)
      return(tree)
    },
    .leaves=function(node) {
      if (node$nchildren == 0) {
        return(list(node))
      } else {
        ret = list()
        for (i in 1:node$nchildren) {
          ret = c(ret, .leaves(node$children[[i]]))
        }
        return(ret)
      }
    },
    .generatePseudoGUID = function(size) {
      chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
      ret = c()
      indices = sample(1:nchar(chars), size, replace=T)

      for (i in 1:size) {
        ret = c(ret, substr(chars, indices[i], indices[i]))
      }

      return(paste(ret, collapse=""))
    },
    plot=function(...) {
      # TODO
      # ms <- getMeasurements()
      # if (length(ms)<2)
      #   stop("need at least two columns to plot")

      # mgr$scatterChart(x=ms[[1]], y=ms[[2]], ...)
    }
  )
)

EpivizMetagenomicsData$methods(
  getMeasurements=function() {
    out <- lapply(samples, function(sample) {
      list(id=sample$id,
           name=sample$name,
           type="feature",
           datasourceId=id,
           datasourceGroup=id, # TODO: Something specific to this taxonomy
           defaultChartType="heatmap",
           annotation=NULL, # TODO: observationType, ageRange, etc
           minValue=0, # TODO
           maxValue=15, # TODO
           metadata=NULL) # TODO
    })
    out
  },
  getHierarchy=function() {
    taxonomy
  }
#   parseMeasurement=function(msId) {
#     column <- strsplit(msId, split="__")[[1]][2]
#     if(!.checkColumns(column)) {
#       stop("invalid parsed measurement")
#     }
#     column
#   }
)

#######################################################################

# TODO
# .valid.EpivizMetagenomicsData.object <- function(x) {
#   if(!is(x$object, "SummarizedExperiment"))
#     return("'object' must be of class 'SummarizedExperiment'")
#   if(!is(rowData(x$object), "GIntervalTree"))
#     return("'rowData(object)' must be of class 'GIntervalTree'")
#   NULL
# }
#
# .valid.EpivizMetagenomicsData.ylim <- function(x) {
#   if(!is(x$ylim, "matrix"))
#     return("'ylim' must be a matrix")
#   if(nrow(x$ylim) != 2)
#     return("'ylim' must have two rows")
#   if(ncol(x$ylim) != length(x$columns))
#     return("'ylim' must have 'length(columns)' columns")
#   NULL
# }
#
# .valid.EpivizMetagenomicsData.assay <- function(x) {
#   if (is.character(x$assay)) {
#     if(!(x$assay %in% names(assays(x$object))))
#       return("'assay' not found in 'object'")
#     return(NULL)
#   }
#
#   if (x$assay > length(assays(x$object)))
#     return("'assay' not found in 'object'")
#   NULL
# }
#
# .valid.EpivizMetagenomicsData <- function(x) {
#   c(.valid.EpivizMetagenomicsData.object(x),
#     .valid.EpivizMetagenomicsData.ylim(x),
#     .valid.EpivizMetagenomicsData.assay(x))
# }
#
# setValidity2("EpivizMetagenomicsData", .valid.EpivizMetagenomicsData)
