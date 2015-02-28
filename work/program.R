# install.packages("devtools")
# install.packages("jsonlite")
# source("http://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")
# biocLite("msd16s")
# biocLite("epivizr")
library("devtools")
library("metagenomeSeq")
library("msd16s")

load_all("../epivizr")
load_all("./")

# TODO Joe
filteredExp = filterData(msd16s,depth=6000,present=25)

mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "qSJzFdtOFPq", useDevel=FALSE, debug=TRUE, verbose=TRUE)
ms <- mgr$addMeasurements(filteredExp, "Bacteriome Phylogenetic Tree")
tree = buildEpivizTree(msd16s)
mgr$service()
#taxonomi_vis <- mgr$addDevice(msd16s, "Bacteriome Phylogenetic Tree")

mgr$stopServer()

str_pad(as.hexmode((1:100)*260), 0, pad="0")
