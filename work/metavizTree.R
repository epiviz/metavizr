require(msd16s)
example = filterData(msd16s,depth=6000,present=25)
tax = colnames(fData(example))[c(3:9,1)]
taxonomy = fData(example)[,tax]
taxonomy[,1] = "Bacteria"
df = taxonomy[sample(1:nrow(taxonomy),10),]

#assayData(example)[["counts"]]
counts = MRcounts(example,norm=TRUE,log=TRUE)
pd = pData(example)
pd["100489",]
l = as.list(pd["100489",])

tableToJSONTree <- function(t) {

  tableToListTree <- function(t, colIndex=1) {
    if (colIndex > dim(t)[2]) { return(NULL) }
    groups = by(t, t[, colIndex], list, simplify=F)
    nodes = lapply(groups, function(group){
      tableToListTree(group[[1]], colIndex + 1)
    })
    return(nodes)
  }

  listTreeToJSONTree <- function(node, globalDepth=0) {
    ret = list(
      name=names(node),
      id=NULL,# TODO
      parentId=NULL, # TODO
      globalDepth=globalDepth,
      taxonomy=colnames(t)[globalDepth+1],
      nchildren=length(node[[1]]),
      size=1,
      selectionType=NULL # TODO
    )
    nleaves = 0
    if (length(node[[1]]) > 0) {
      children = c()
      for (i in seq_along(node[[1]])) {
        child = listTreeToJSONTree(node[[1]][i], globalDepth + 1)
        child$order = i - 1
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
}

toTreeList <- function(t) {
  rootList = Ptr$new(list())

  for (j in 1:dim(t)[1]) {
    node = rootList
    for (i in 2:dim(t)[2]) {
      nodeName = as.character(t[j, i])
      if (is.na(nodeName)) { nodeName = "NA" }
      if (is.null(node$.[[nodeName]])) {
        node$.[[nodeName]] = Ptr$new(list())
      }
      node = node$.[[nodeName]]
    }
  }
  return(rootList)
}

leaves <- function(node) {
  if (node$nchildren == 0) {
    return(list(node))
  } else {
    ret = list()
    for (i in 1:node$nchildren) {
      ret = c(ret, leaves(node$children[[i]]))
    }
    return(ret)
  }
}

# Test method for getMeasurements() method inside EpivizMetagenomicsData
getMeasurements=function(samples, id) { # id -> datasourceId
  out <- lapply(samples, function(sample) {
    list(id=sample$id,
         name=sample$name,
         type="feature",
         datasourceId=id,
         datasourceGroup=id,
         defaultChartType="heatmap",
         annotation=NULL, # TODO: observationType, ageRange, etc
         minValue=0, # TODO
         maxValue=15, # TODO
         metadata=NULL) # TODO
  })
  out
}

generatePseudoGUID = function(size) {
  chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
  ret = c()
  indices = sample(1:nchar(chars), size, replace=T)

  for (i in 1:size) {
    ret = c(ret, substr(chars, indices[i], indices[i]))
  }

  return(paste(ret, collapse=""))
}

.tableToListTree <- function(t, colIndex=1, parent=NULL) {
  if (colIndex > dim(t)[2]) {
    if (!is.null(parent)) { parent$nleaves = 1 }
    return(NULL)
  }

  parentId = NULL
  if (!is.null(parent)) { parentId = parent$id }

  groups = by(t, t[, colIndex], list, simplify=F)
  names = names(groups)
  nodes = lapply(seq_along(names), function(i) { EpivizNode$new(name=names[i], parentId=parentId, depth=colIndex-1, taxonomy=colnames(t)[colIndex], order=i) })
  if (!is.null(parent)) {
    parent$nchildren = length(nodes)
    parent$childrenIds = lapply(nodes, function(node) { node$id })
  }

  descendants = unlist(lapply(seq_along(groups), function(i) {
    node = nodes[[i]]
    children = .tableToListTree(groups[[i]][[1]], colIndex+1, node)
    if (!is.null(parent)) {
      parent$nleaves = parent$nleaves + node$nleaves
    }
    children
  }))
  ret = c(nodes, descendants)

  return(ret)
}

tree = tableToJSONTree(df)
ls = leaves(tree)
ms = getMeasurements(ls, "my-id")
# View(RJSONIO::toJSON(tree,recursive=FALSE))

order(node$children, function(c1, c2) { c1$order - c2$order })


#####################################################################


require(msd16s)
example = filterData(msd16s,depth=6000,present=25)
tax = colnames(fData(example))[c(3:9,1)]
taxonomy = fData(example)[,tax]
taxonomy[,1] = "Bacteria"
df = taxonomy#axonomy[sample(1:nrow(taxonomy),10),]

#assayData(example)[["counts"]]
counts = MRcounts(example,norm=TRUE,log=TRUE)
pd = pData(example)
pd["100489",]
l = as.list(pd["100489",])

root = MetavizNode$new(taxonomyTablePtr=Ptr$new(df))
x = root$raw(maxDepth=0)
y = root$children()
tree = MetavizTree$new(taxonomyTablePtr=Ptr$new(df))
z = tree$selectedLeaves(0,100)


filteredExp = msd16s#filterData(msd16s,depth=3000,present=350)[,sample(1:599,10)]
fData(filteredExp) = fData(filteredExp)[,c(6,7,8,9,1)]
df = fData(filteredExp)

d = filterData(msd16s,depth=3000,present=350)[,sample(1:599,10)]
fData(d) = fData(d)[,c(6,7,8,9,1)]
t = fData(d)

