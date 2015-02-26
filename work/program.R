# install.packages("devtools")
# install.packages("jsonlite")
# source("http://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")
# biocLite("msd16s")
# biocLite("epivizr")
library("devtools")
library("metagenomeSeq")
library("msd16s")

##########################################################
#library("RJSONIO")
#library("jsonlite")
#
# example = filterData(msd16s,present=25,depth=6000)
# example
#
# tax = colnames(fData(example))[c(3:9,1)]
# taxonomy = fData(example)[,tax]
# head(taxonomy)
# toJSON(taxonomy)
# getwd()

#toj = RJSONIO::toJSON
##########################################################

setwd("d:\\EpiViz\\github\\metavizr\\")
load_all("../epivizr")
load_all("./")

toj = epivizr::toJSON

mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "qSJzFdtOFPq", useDevel=FALSE, debug=TRUE, verbose=TRUE)
ms <- mgr$addMeasurements(msd16s, "Bacteriome Phylogenetic Tree")
mgr$service()
#taxonomi_vis <- mgr$addDevice(msd16s, "Bacteriome Phylogenetic Tree")

mgr$stopServer()
