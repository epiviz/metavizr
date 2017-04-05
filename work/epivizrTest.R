# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
library("devtools")

# getwd()

load_all("../epivizr")

data(tcga_colon_example)

mgr = startEpiviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "eOcpNKJ8GE", useDevel=FALSE, debug=FALSE, verbose=TRUE)
mgr$service()
blocks_dev <- mgr$addDevice(colon_blocks, "450k colon_blocks")

mgr$stopServer()
