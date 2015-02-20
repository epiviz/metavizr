require(msd16s)
example = filterData(msd16s,depth=6000,present=25)
tax = colnames(fData(example))[c(3:9,1)]
taxonomy = fData(example)[,tax]
taxonomy[,1] = "Bacteria"
df = taxonomy[sample(1:nrow(taxonomy),10),]

metavizNode <-function(df,depth){
  list(
    list(children=metavizTreeRecursion(df[[1]][,-1,drop=FALSE],depth=depth+1)),
    nleaves = nrow(df[[1]]),
    globalDepth = depth,
    taxonomy = colnames(df)[1]
  )
}
# metavizFinalNode <-function(df,depth){
#   list(
#     df,
#     nleaves=nrow(df),
#     globalDepth=depth,
#     name="OTU",
#     taxonomy=colnames(df)[1]
#   )
# }
metavizTreeRecursion <- function(df,depth=0){
  # if(depth==1) return(list())
  if(ncol(df)==1) return(metavizNode(df,depth))
  r = by(df,df[,1],list,simplify=F)
  tr = lapply(r,function(df2){
    unlist(metavizNode(df2,depth),recursive=FALSE)
  })
  for(i in seq_along(tr)){ 
    tr[[i]]$name = names(tr)[i]; 
    names(tr[[i]])[1] = "children"
  }
  tr
}
metavizTree<-function(df){
  for(i in seq_along(df)) df[is.na(df[,i]),i] = paste(colnames(df)[i],"NA",sep=":")
  listoflists = metavizTreeRecursion(df)
  names(listoflists) = ""
  # switch to rjson
  tr = RJSONIO::toJSON(unlist(listoflists,recursive=FALSE))
  tr
}
metavizTree(df)
writeLines(metavizTree(df),"~/Desktop/test.txt")