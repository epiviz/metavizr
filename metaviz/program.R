# install.packages("devtools")
library("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite("epivizr")
load_all("../epivizr")

mgr = startEpiviz(workspace = "eOcpNKJ8GE", useDevel=FALSE, debug=FALSE, verbose=TRUE)
mgr$service()

data(tcga_colon_example)

blocks_dev <- mgr$addDevice(colon_blocks, "450k colon_blocks")

mgr$stopServer()
