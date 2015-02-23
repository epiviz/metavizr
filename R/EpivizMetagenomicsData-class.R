EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains="EpivizData",
  fields=list(taxonomy="ANY", leaves="list", samples="list", maxDepth="numeric"),
  methods=list(
    initialize=function(object, ...) {
      taxonomy <<- .MRexperimentToTree(object)
      leaves <<- .leaves(taxonomy)
      samples <<- list(list(name="sample1", id="s1"), list(name="sample2", id="s2")) # TODO
      maxDepth <<- 4
      callSuper(object=object, ...)
    },
    update=function(newObject, ...) {
      # TODO
      callSuper(newObject, ...)
    },
    .MRexperimentToTree=function(exp) {
      # TODO Joe
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
  getHierarchy=function(nodeId=taxonomy$id) {
    # TODO Joe (Test)
    # grab subtree with maxDepth
    # TODO: If nodeId is NULL or not defined then return the entire tree (up to maxDepth)
    if (is.null(nodeId)) { nodeId = taxonomy$id }
    env = new.env()
    env$ret = NULL
    .traverse(taxonomy, function(node) {
      if (node$id == nodeId) { env$ret = node }
    })

    rootGlobalDepth = env$ret$globalDepth
    select = function(node) {
      if (node$globalDepth - rootGlobalDepth >= maxDepth) { return(NULL) }
      ret = node
      ret$children = NULL
      if (node$globalDepth - rootGlobalDepth < maxDepth - 1) {
        ret$children = c()
        for (i in 1:node$nchildren) {
          ret$children = c(ret$children, list(select(node$children[[i]])))
        }
      }
      return(ret)
    }

    select(env$ret)
  },
  propagateHierarchyChanges=function(selection, order) {
    # TODO Joe
    # selection is a list where the labels (names) are node ids, and values are numbers:
    #    0 means "none"
    #    1 means "leaves"
    #    2 means "node aggregation"
    # example: list(Qy2uUZ=1, abcdef=2,efghij=0, ...)
    # what you need to do here is, for each node in the taxonomy tree, change the "selectionType" field to the value given in the list

    # order has the same format as selection (labels to numbers) but numbers correpond to the order within the node parend
    # example:
    # parent: abc
    # children: x (3), y (1), z (2)
    # so in the "taxonomy" tree, we need to re-sort nodes according to these numbers
    # and also set the "order" field

    # output: taxonomy after the changes have been made, in tree format, right before it is transformed to JSON
#     .traverse(taxonomy, function(node) {
#       node$children <<- sort(node$children) # This will not work, but something with this semantic
#     })
    browser()
    taxonomy
  },
  .traverse=function(node, callback) {
    callback(node)
    if (node$nchildren == 0) { return() }
    for (i in 1:node$nchildren) {
      .traverse(node$children[[i]], callback)
    }
  }

)
