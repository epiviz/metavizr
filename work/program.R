# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")
# biocLite("msd16s")
# biocLite("epivizr")
library("devtools")
library("metagenomeSeq")
library("msd16s")
library("rjson")

example = filterData(msd16s,present=25,depth=6000)
example

tax = colnames(fData(example))[c(3:9,1)]
taxonomy = fData(example)[,tax]
head(taxonomy)
toJSON(taxonomy)


getwd()
load_all("./")
#reload("../R")

mgr = startEpiviz(workspace = "eOcpNKJ8GE", useDevel=FALSE, debug=FALSE, verbose=TRUE)
mgr$service()

data(tcga_colon_example)

# blocks_dev <- mgr$addDevice(colon_blocks, "450k colon_blocks")
taxonomi_vis <- mgr$addDevice(example, "Bacteriome Phylogenetic Tree")

mgr$stopServer()


Node <- setRefClass("Node",
  fields=list(children="list")
)

node=Node$new(children=list("a", "b", "c"))
RJSONIO::toJSON(node)
