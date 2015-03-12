#' Indices of ranked features
#'
#' Calculate top n features at the leaves. Returns indices for the rankings.
#' 
#' @param obj MRexperiment object.
#' @param control List of options passed through `metavizControl`.
#' @return value
#' @export
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
	norm=control$norm
	log = control$log
	n   = control$n
	rankFun = control$rankFun
	mat = MRcounts(obj,norm=norm,log=log)
	otusToKeep = which(rowSums(mat) > 0)
	n = min(length(otusToKeep),n)
	otuStats = apply(mat[otusToKeep, ], 1, rankFun)
	otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:n]]
	otuIndices
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
#' @export
#' @seealso \code{\link{metavizControl}} \code{\link{startMetaviz}}
#' @examples
#' library(msd16s)
#' msd16s = filterData(msd16s,present=100)
#' ind = which(pData(msd16s)$Type=="Control" & pData(msd16s)$Country=="Gambia")
#' settings = metavizControl(aggregateAtDepth="genus",n=10)
#' obj = metavizTransformSelect(msd16s[,ind],control=settings)
#' 
metavizTransformSelect<-function(obj,fun=rowSums ,control=metavizControl(n=100)){
	aggregateAtDepth = control$aggregateAtDepth
	if(is.character(aggregateAtDepth)) aggregateAtDepth = assignValues(obj,aggregateAtDepth)
    n = control$n
	tree = fData(obj)
    
	treeLvl = colnames(tree)[aggregateAtDepth+1]
	agg = aggregateByTaxonomy(obj,treeLvl)
	topX = rownames(agg)[order(fun(agg),decreasing=TRUE)[1:n]]
	ind = which(!tree[,aggregateAtDepth+1]%in%topX)
	for(i in 1:ncol(tree)) tree[,i] = as.character(tree[,i])
	for(i in 2:ncol(tree)) tree[ind,i] = "Other"

	obj = newMRexperiment(MRcounts(obj),phenoData=phenoData(obj),
		featureData=AnnotatedDataFrame(tree),normFactors=normFactors(obj),
		libSize=libSize(obj))
	return(obj)
}
#' metavizr settings
#'
#' Default settings for the various plotting functions in metavizr.
#' 
#' @param aggregateAtDepth Level of the tree to aggregate counts at by default.
#' @param aggregateFun Function to aggregate counts by at the aggregateAtDepth level.
#' @param maxDepth Level of the tree to display by default in icicle view.
#' @param maxHistory Value for caching.
#' @param maxValue Maximum value to display.
#' @param minValue Minimum value to display.
#' @param title title.
#' @param n Number of OTUs to include in ranking.
#' @param rankFun Ranking function - single vector function.
#' @param norm Normalize MRexperiment object.
#' @param log Log tranformation of MRexperiment object.
#' @return List of setting parameters.
#' @export
#' @seealso \code{\link{metavizControl}} \code{\link{metavizRank}} \code{\link{metavizOptimize}}
#' @examples
#' settings = metavizControl()
#'
metavizControl<-function(aggregateAtDepth=3,aggregateFun=function(x) log2(1 + colSums(x)),
						maxDepth=4,maxHistory=3,maxValue=NULL,minValue=NULL,title="",
						n=10000,rankFun=sd,norm=TRUE,log=FALSE){

	list(aggregateAtDepth=aggregateAtDepth,aggregateFun=aggregateFun,
						maxDepth=maxDepth,maxHistory=maxHistory,maxValue=maxValue,minValue=minValue,title=title,
						n=n,rankFun=sd,norm=TRUE,log=TRUE)
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