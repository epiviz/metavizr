#' Line plot
#'
#' Produces a line plot for the experiment
#' 
#' @param obj MRexperiment object or EpivizMetagenomicsData object. 
#' 		If EpivizMetagenomicsData object then control is ignored.
#' @param mgr manager of a metaviz session.
#' @param control List of options passed through `metavizControl`.
#' @return EpivizMetagenomicsData class object
#' @export
#' @seealso \code{\link{metavizControl}} \code{\link{metaviztree}} \code{\link{metavizRank}} \code{\link{metavizOptimize}}
#' @examples
#' library(msd16s)
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
#' metavizLine(obj,mgr=mgr,metavizControl(aggregateAtDepth="class",n=nrow(obj)))
metavizLine<-function(obj,mgr,control=metavizControl(title="line plot")){
	if(!class(obj)%in%c("MRexperiment","EpivizMetagenomicsData")){
		stop("Either a MRexperiment or EpivizMetagenomicsData")
	}
	if(class(obj)=="MRexperiment"){
		otuIndices = metavizRank(obj,control)
		obj = mgr$addMeasurements(obj[otuIndices,],msName=control$title,control=control)
	} 
	measurements = obj$getMeasurements()
	line = mgr$visualize("lineplot", measurements)
	invisible(obj)
}
#' Icicle
#'
#' Produces an icicle for the tree of an EpivizMetagenomicsData class object
#' 
#' @param datasource An EpivizMetagenomicsData class object (i.e. output from one of the metavizPlot functions).
#' @param mgr manager of a metaviz session.
#' @return Icicle plot
#' @export
#' @seealso \code{\link{startMetaviz}}
#' @examples
#' # assuming datasource, and mgr
#' # metavizTree(datasource,mgr)
#'
metavizTree<-function(datasource,mgr){
	invisible(mgr$visualize("icicle", datasource=datasource))
}
#' Stacked plot
#'
#' Produces the stacked plot for the experiment.
#' 
#' @param obj MRexperiment object or EpivizMetagenomicsData object. 
#' 		If EpivizMetagenomicsData object then control is ignored.
#' @param mgr manager of a metaviz session.
#' @param control List of options passed through `metavizControl`.
#' @return EpivizMetagenomicsData class object
#' @export
#' @seealso \code{\link{metavizControl}} \code{\link{metaviztree}}
#' 		\code{\link{metavizRank}} \code{\link{metavizOptimize}}
#' @examples
#' # This will transform the tree to top 9 selected genera due to average abundance 
#' #	of a subset of samples. Then plot the relative abundances of those 9 and 'others'.
#' # mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", 
#' # 	workspace = "Ey4CYuaTjNd",useDevel=TRUE, debug=TRUE, verbose=TRUE)
#' library(msd16s)
#' msd16s = filterData(msd16s,present=100)
#' fData(msd16s) = fData(msd16s)[,c(3:9, 1)]
#' ind = which(pData(msd16s)$Type=="Control" & pData(msd16s)$Country=="Gambia")
#' obj = metavizTransformSelect(msd16s[,ind],fun=rowMeans,control=metavizControl(aggregateAtDepth="genus",n=10))
#' # gates = metavizStack(obj,mgr,control=metavizControl(aggregateFun=colSums,aggregateAtDepth=5,log=FALSE))
#'
metavizStack<-function(obj,mgr,control=metavizControl(title="stacked plot",aggregateFun=colSums)) {
	if(!class(obj)%in%c("MRexperiment","EpivizMetagenomicsData")){
		stop("Either a MRexperiment or EpivizMetagenomicsData")
	}
	if(class(obj)=="MRexperiment"){
		otuIndices = metavizRank(obj,control)
		obj = mgr$addMeasurements(obj[otuIndices,],msName=control$title,control=control)
	} 
	measurements = obj$getMeasurements()
	stack = mgr$visualize("stackedplot", measurements)
	invisible(obj)
}
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
#' @param valuesAnnotationFuns Function for error bars.
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
						valuesAnnotationFuns=NULL,
						maxDepth=4,maxHistory=3,maxValue=NULL,minValue=NULL,title="",
						n=10000,rankFun=sd,norm=TRUE,log=FALSE){

	list(aggregateAtDepth=aggregateAtDepth,aggregateFun=aggregateFun,
						maxDepth=maxDepth,maxHistory=maxHistory,maxValue=maxValue,minValue=minValue,title=title,
						n=n,rankFun=sd,norm=norm,log=log)
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