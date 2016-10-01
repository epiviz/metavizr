#' Heatmap
#'
#' Produces a heatmap for the experiment.
#'
#' @param obj MRexperiment object or EpivizMetagenomicsData object.
#'     If EpivizMetagenomicsData object then control is ignored.
#' @param mgr manager of a metaviz session.
#' @param control List of options passed through `metavizControl`.
#' @param samples Index vector of samples to include in plot.
#' @return EpivizMetagenomicsData class object
#' @seealso \code{\link{metavizControl}} \code{\link{metaviztree}} \code{\link{metavizRank}} \code{\link{metavizOptimize}}
#' @examples
#' #todo
metavizHeatmap<-function(obj,mgr,control=metavizControl(title="heatmap"),samples=NULL) {
  # if(!class(obj)%in%c("MRexperiment","EpivizMetagenomicsData")){
  #   stop("Either a MRexperiment or EpivizMetagenomicsData")
  # }
  # if(class(obj)=="MRexperiment"){
  #   otuIndices = metavizRank(obj,control)
  #   obj = mgr$add_measurements(obj[otuIndices,],msName=control$title,control=control)
  # }
  # measurements = obj$get_measurements()
  # if(is.null(samples)){ samples = seq(measurements) }
  # heatmap = mgr$visualize("heatmap", measurements[samples])
  # invisible(obj)
  NULL
}
#' Line plot
#'
#' Produces a line plot for the experiment
#'
#' @param obj MRexperiment object or EpivizMetagenomicsData object.
#'     If EpivizMetagenomicsData object then control is ignored.
#' @param mgr manager of a metaviz session.
#' @param control List of options passed through `metavizControl`.
#' @param samples Index vector of samples to include in plot.
#' @return EpivizMetagenomicsData class object
#' @seealso \code{\link{metavizControl}} \code{\link{metaviztree}} \code{\link{metavizRank}} \code{\link{metavizOptimize}}
#' @examples
#' library(msd16s)
#' msd16s = filterData(msd16s,present=100)
#' ind = which(pData(msd16s)$Type=="Control" & pData(msd16s)$Country=="Mali")
#' obj = metavizTransformSelect(msd16s[,ind],fun=rowMeans,control=metavizControl(aggregateAtDepth="phylum",n=10))
#' # gates = metavizLine(obj,mgr,control=metavizControl(aggregateAtDepth="phylum"))
#'
#' # second example with time-series data for an individual mouse
#' data(mouseData)
#' time = order(as.Date(pData(mouseData)[,2],format="%Y-%m%-%d"))
#' mouseData = mouseData[,time]
#' status = pData(mouseData)[,"status"]
#' pm = pData(mouseData)[,"mouseID"]
#' obj = metavizTransformSelect(mouseData[,which(pm=="PM10")],fun=rowSums,control=metavizControl(aggregateAtDepth="class",n=5))
#' # metavizLine(obj,mgr=mgr,metavizControl(aggregateAtDepth="class",n=nrow(obj)))
metavizLine<-function(obj,mgr,control=metavizControl(title="line plot"),samples=NULL){
  # if(!class(obj)%in%c("MRexperiment","EpivizMetagenomicsData")){
  #   stop("Either a MRexperiment or EpivizMetagenomicsData")
  # }
  # if(class(obj)=="MRexperiment"){
  #   otuIndices = metavizRank(obj,control)
  #   obj = mgr$add_measurements(obj[otuIndices,],msName=control$title,control=control)
  # }
  # measurements = obj$get_measurements()
  # if(is.null(samples)){ samples = seq(measurements) }
  # line = mgr$visualize("lineplot", measurements[samples])
  # invisible(obj)
  NULL
}
#' Icicle
#'
#' Produces an icicle for the tree of an EpivizMetagenomicsData class object
#'
#' @param datasource An EpivizMetagenomicsData class object (i.e. output from one of the metavizPlot functions).
#' @param mgr manager of a metaviz session.
#' @return Icicle plot
#' @seealso \code{\link{startMetaviz}}
#' @examples
#' # assuming datasource, and mgr
#' # metavizTree(datasource,mgr)
#'
metavizTree<-function(datasource,mgr){
  # invisible(mgr$visualize("icicle", datasource=datasource))
  NULL
}
#' Stacked plot
#'
#' Produces the stacked plot for the experiment.
#'
#' @param obj MRexperiment object or EpivizMetagenomicsData object.
#'     If EpivizMetagenomicsData object then control is ignored.
#' @param mgr manager of a metaviz session.
#' @param control List of options passed through `metavizControl`.
#' @param samples Index vector of samples to include in plot.
#' @return EpivizMetagenomicsData class object
#' @seealso \code{\link{metavizControl}} \code{\link{metaviztree}}
#'     \code{\link{metavizRank}} \code{\link{metavizOptimize}}
#' @examples
#' # This will transform the tree to top 9 selected genera due to average abundance
#' #  of a subset of samples. Then plot the relative abundances of those 9 and 'others'.
#' # mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/",
#' #   workspace = "Ey4CYuaTjNd",useDevel=TRUE, debug=TRUE, verbose=TRUE)
#' library(msd16s)
#' msd16s = filterData(msd16s,present=100)
#' fData(msd16s) = fData(msd16s)[,c(3:9, 1)]
#' ind = which(pData(msd16s)$Type=="Control" & pData(msd16s)$Country=="Gambia")
#' obj = metavizTransformSelect(msd16s[,ind],fun=rowMeans,control=metavizControl(aggregateAtDepth="genus",n=10))
#' # gates = metavizStack(obj,mgr,control=metavizControl(aggregateFun=colSums,aggregateAtDepth=5,log=FALSE))
#'
metavizStack<-function(obj,mgr,control=metavizControl(title="stacked plot",aggregateFun=colSums),samples=NULL) {
  # if(!class(obj)%in%c("MRexperiment","EpivizMetagenomicsData")){
  #   stop("Either a MRexperiment or EpivizMetagenomicsData")
  # }
  # if(class(obj)=="MRexperiment"){
  #   otuIndices = metavizRank(obj,control)
  #   obj = mgr$add_measurements(obj[otuIndices,],control=control, datasource_name=control$title)
  # }
  # measurements = obj$get_measurements()
  # stack = mgr$visualize("stackedplot", measurements)
  # invisible(obj)
  NULL
}

#' MDS (multiple dimensional scaling) or PCA
#'
#' Produces a PCA or PCoA scatterplot of an experiment.
#'
#' @param obj MRexperiment object.
#' @param mgr manager of a metaviz session.
#' @param usePCA TRUE/FALSE whether to use PCA or MDS coordinates (TRUE isPCA).
#' @param useDist TRUE/FALSE whether to calculate distance matrix.
#' @param cord Which coordinates to compare/display.
#' @param distFun Distance function - default is stats::dist.
#' @param distMethod If useDist==TRUE, what method to calculate distances.
#' @param tran Transpose the matrix.
#' @param control List of options passed through `metavizControl`.
#' @return EpivizMetagenomicsData class object
#' @seealso \code{\link{metavizControl}} \code{\link{metavizRank}} \code{\link{cmdscale}} \code{\link{prcomp}}
#' @examples
#' data(mouseData)
#' # metavizOrd(mouseData,mgr)
#'
metavizMDS<-function(obj,mgr,usePCA=FALSE,useDist=TRUE,cord=c(1,2),
          distFun=stats::dist,distMethod="euclidian",tran=FALSE,
          control=metavizControl(title="MDS plot")){
  # if(useDist == FALSE & usePCA == FALSE)
  #   stop("Classical MDS requires distances")
  # control$n = min(nrow(obj),control$n)
  # otuIndices = metavizRank(obj,control)
  # d = MRcounts(obj[otuIndices, ],norm=TRUE,log=TRUE)
  # if(tran == FALSE) { d = t(d) }
  # if(useDist == TRUE) { d = distFun(d, method = distMethod) }
  # if(usePCA == FALSE) {
  #   ord = cmdscale(d, k = dim(as.matrix(d))[1]-1 )
  #   colnames(ord) = paste("MDS",1:ncol(ord),sep="")
  #   df = data.frame(root="MDS",components = paste("MDS",1:ncol(ord),sep=""))
  #   rownames(df) = colnames(ord)
  # } else {
  #   ord = prcomp(d)$x
  #   df = data.frame(root="PCA",PCs = paste("PC",1:ncol(ord),sep=""))
  #   rownames(df) = colnames(ord)
  # }
  # pd = pData(obj)
  # pd = cbind(Sample = "sample",pd)
  # #pd = cbind(pd,sampleNames = factor(rownames(pd)))
  # tmp = newMRexperiment(ord,featureData=AnnotatedDataFrame(pd),
  # normFactors=rep(1000,ncol(ord)),phenoData=AnnotatedDataFrame(df))
  # # override parameters
  # control$aggregateAtDepth = -1
  # control$aggregateFun = identity
  # control$minValue = min(ord)
  # control$maxValue = max(ord)
  # control$norm=FALSE
  # control$log=FALSE
  # 
  # metavizrData = mgr$add_measurements(tmp,control$title,control=control)
  # measurements = metavizrData$get_measurements()
  # scatterplot = mgr$visualize("scatterplot", list(measurements[[cord[1]]],measurements[[cord[2]]]))
  # invisible(metavizrData)
  NULL
}

#' Indices of ranked features
#'
#' Calculate top n features at the leaves. Returns indices for the rankings.
#'
#' @param obj MRexperiment object.
#' @param control List of options passed through `metavizControl`.
#' @return value
#' @seealso \code{\link{metavizControl}} \code{\link{metavizRank}}
#' @examples
#' # 10 highest mean OTUs
#' indices = metavizRank(msd16s[,1:10],control=metavizControl(n=10,rankFun=mean,log=TRUE))
#' # 5 highest sd OTUs
#' indices = metavizRank(msd16s[,1:10],control=metavizControl(n=5, rankFun=sd,log=TRUE))
#' # 20 highest max OTUs
#' indices = metavizRank(msd16s[,1:10],control=metavizControl(n=20,rankFun=max,log=TRUE))
#'
metavizRank<-function(obj,control=metavizControl(log=TRUE)){
#   norm=control$norm
#   log = control$log
#   n   = control$n
#   rankFun = control$rankFun
#   mat = MRcounts(obj,norm=norm,log=log)
#   otusToKeep = which(rowSums(mat) > 0)
#   n = min(length(otusToKeep),n)
#   otuStats = apply(mat[otusToKeep, ], 1, rankFun)
#   otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:n]]
#   otuIndices
  NULL
}
#' Tranform the tree of a MRexperiment object
#'
#' Helps in selecting the top n features at a particular level of the tree. Converts features
#' tree to a different annotation, in particular an "Other" field. This is useful,
#' when there are potentially too many leaves to display.
#'
#' @param obj MRexperiment object.
#' @param fun Row based function for ranking top OTUs.
#' @param control List of options passed through `metavizControl`.
#' @return Modified MRexperiment object.
#' @seealso \code{\link{metavizControl}} \code{\link{startMetaviz}}
#' @examples
#' library(msd16s)
#' msd16s = filterData(msd16s,present=100)
#' ind = which(pData(msd16s)$Type=="Control" & pData(msd16s)$Country=="Gambia")
#' settings = metavizControl(aggregateAtDepth="genus",n=10)
#' obj = metavizTransformSelect(msd16s[,ind],control=settings)
#'
metavizTransformSelect<-function(obj,fun=rowSums ,control=metavizControl(n=100)){
  # aggregateAtDepth = control$aggregateAtDepth
  # if(is.character(aggregateAtDepth)) aggregateAtDepth = assignValues(obj,aggregateAtDepth)
  #   n = control$n
  # tree = fData(obj)
  # 
  # treeLvl = colnames(tree)[aggregateAtDepth+1]
  # agg = aggregateByTaxonomy(obj,treeLvl,out="matrix")
  # topX = rownames(agg)[order(fun(agg),decreasing=TRUE)[1:n]]
  # ind = which(!tree[,aggregateAtDepth+1]%in%topX)
  # for(i in 1:(ncol(tree)-1)) tree[,i] = as.character(tree[,i])
  # for(i in 2:(ncol(tree)-1)) tree[ind,i] = "Other"
  # 
  # obj = newMRexperiment(MRcounts(obj),phenoData=phenoData(obj),
  #   featureData=AnnotatedDataFrame(tree),normFactors=normFactors(obj),
  #   libSize=libSize(obj))
  # return(obj)
  NULL
}

#' Check appropriate values
#'
#' Ensure options sent to addMeasurements are reasonable.
#'
#' @param obj MRexperiment object.
#' @param option Option value.
#' @return NULL
#' @seealso \code{\link{metavizControl}} \code{\link{startMetaviz}}
#'
checkValues<-function(obj,option){
  nc = ncol(fData(obj))
  if(nc<1) stop("MRexperiment must have an annotation")
  if(option>nc) stop(sprintf("value must be less than %s",nc))
  if(option<0) stop("value must be greater than zero")
}
#' Assign values of character strings
#'
#' Convert values set for aggregateAtDepth or maxDepth to integers.
#'
#' @param obj MRexperiment object.
#' @param option Option value, either character or numeric.
#' @return value
#' @seealso \code{\link{metavizControl}} \code{\link{startMetaviz}}
#'
assignValues<-function(obj,option){
  if(is.character(option)){
    option = (which(colnames(fData(obj))==option) - 1)
  }
  checkValues(obj,option)
  return(option)
}
