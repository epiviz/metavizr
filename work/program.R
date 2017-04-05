# install.packages("devtools")
# install.packages("jsonlite")
# source("http://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")
# biocLite("msd16s")
# biocLite("epivizr")
library("devtools")
#library("metagenomeSeq")
library("msd16s")
#library("stringr")

setwd("d:\\EpiViz\\github\\metavizr\\")
load_all("../epivizr")
load_all("./")

filteredExp = msd16s
fData(filteredExp) = fData(filteredExp)[,c(3:9, 1)]

#mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "6SbPOUDtKzg", useDevel=FALSE, debug=TRUE, verbose=FALSE)
#mgr = startMetaviz(localURL="http://localhost/epiviz-dev", workspace = "YGwCd2zrYOs", useDevel=FALSE, debug=TRUE, verbose=FALSE)
mgr = startMetaviz(localURL="http://localhost/epiviz-dev/index-standalone.html", workspace = "YGwCd2zrYOs", useDevel=FALSE, debug=FALSE, verbose=FALSE)
metavizrData = mgr$addMeasurements(filteredExp, "MSD 16", control=metavizControl(maxDepth=5, aggregateAtDepth=5, minValue=0, maxValue=15))

samples = metavizrData$getMeasurements()

scatter = mgr$visualize("scatterplot", list(samples[[1]], samples[[2]]))
heatmap = mgr$visualize("heatmap", list(samples[[3]]))#, samples[[2]]))
heatmap = mgr$visualize("heatmap", samples[1:400])
icicle = mgr$visualize("icicle", datasource=metavizrData)

metavizrData$changeAggregationAtDepth(mgr, depth=2, 2)

mgr$service()
#taxonomi_vis <- mgr$addDevice(msd16s, "Bacteriome Phylogenetic Tree")

mgr$stopServer()
rm(list=ls())

k = rownames(aggTax(msd16s,"genus")[order(rowSums(aggTax(msd16s,lvl="genus")),decreasing=TRUE)])[1:10]
for(i in 4:9) fData(msd16s)[which(!fData(msd16s)[,"genus"]%in%k),i] = "Other"



