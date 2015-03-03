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

# TODO Joe
#filteredExp = filterData(msd16s,depth=6000,present=25)
# filteredExp = filterData(msd16s,depth=3000,present=350)[,sample(1:599,10)]
# fData(filteredExp) = fData(filteredExp)[,c(6,7,8,9,1)]
#counts = MRcounts(filteredExp, norm=TRUE, log=TRUE)
filteredExp = msd16s

# tax = colnames(fData(filteredExp))#[c(3:9,1)]
# taxonomy = fData(filteredExp)[,tax]
# taxonomy[,1] = "Bacteria"
# df = filteredExp
#options(error=recover)

mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "qSJzFdtOFPq", useDevel=FALSE, debug=FALSE, verbose=FALSE)
#mgr = startMetaviz(localURL="http://localhost/epiviz-dev/", workspace = "YGwCd2zrYOs", useDevel=FALSE, debug=TRUE, verbose=TRUE)
ms <- mgr$addMeasurements(filteredExp, "MSD 16", maxDepth=3, aggregateAtDepth=3)
#tree = buildEpivizTree(msd16s)
mgr$service()
#taxonomi_vis <- mgr$addDevice(msd16s, "Bacteriome Phylogenetic Tree")

mgr$stopServer()

#str_pad(as.hexmode((1:100)*260), 0, pad="0")
