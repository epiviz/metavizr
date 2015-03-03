# install.packages("devtools")
# install.packages("jsonlite")
# source("http://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")
# biocLite("msd16s")
# biocLite("epivizr")
library("devtools")
library("metagenomeSeq")
library("msd16s")
library("stringr")

load_all("../epivizr")
load_all("./")

filteredExp = msd16s
fData(filteredExp) = fData(filteredExp)[,c(3:9, 1)]

#mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "6SbPOUDtKzg", useDevel=FALSE, debug=TRUE, verbose=TRUE)
mgr = startMetaviz(localURL="http://localhost/epiviz-dev", workspace = "YGwCd2zrYOs", useDevel=FALSE, debug=TRUE, verbose=TRUE)
metavizrData = mgr$addMeasurements(filteredExp, "MSD 16", maxDepth=3, aggregateAtDepth=3)

measurements = metavizrData$getMeasurements()

scatter = mgr$visualize("scatterplot", list(measurements[[1]], measurements[[2]]))
heatmap = mgr$visualize("heatmap", list(measurements[[1]], measurements[[2]]))
icicle = mgr$visualize("icicle", datasource=metavizrData)

mgr$service()
#taxonomi_vis <- mgr$addDevice(msd16s, "Bacteriome Phylogenetic Tree")

mgr$stopServer()

