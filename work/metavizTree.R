require(msd16s)
example = filterData(msd16s,depth=6000,present=25)
tax = colnames(fData(example))[c(3:9,1)]
taxonomy = fData(example)[,tax]
taxonomy[,1] = "Bacteria"
df = taxonomy[sample(1:nrow(taxonomy),10),]

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
      for (i in 1:length(node[[1]])) {
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

tree = tableToJSONTree(df)
View(RJSONIO::toJSON(tree,recursive=FALSE))
