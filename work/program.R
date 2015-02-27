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

mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "qSJzFdtOFPq", useDevel=FALSE, debug=TRUE, verbose=TRUE)
ms <- mgr$addMeasurements(msd16s, "Bacteriome Phylogenetic Tree")
mgr$service()
#taxonomi_vis <- mgr$addDevice(msd16s, "Bacteriome Phylogenetic Tree")

mgr$stopServer()
