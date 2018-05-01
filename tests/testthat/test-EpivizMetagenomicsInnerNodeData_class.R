context("testing EpivizMetagenomicsInnerNodesData Class")
library(curatedMetagenomicData)

# From curatedMetagenomicData and MicrobiomeWorkshop Vignettes
zeller_file <- system.file("inst", "ZellerG_2014.RData", package = "metavizr")
load(zeller_file)
zeller.eset <- zeller[[1]]
zeller_MR <- ExpressionSet2MRexperiment(zeller.eset)
feature_order <- colnames(fData(zeller_MR))
mObj <- metavizr:::EpivizMetagenomicsDataInnerNodes$new(zeller_MR, feature_order = feature_order)

test_that("getValuesInnerNodes", {
  sampleId<- "CCIS98482370ST-3-0"
  res <- mObj$getValues(measurements = sampleId)
  
  expected <- c(0.0 , 34.1386632, 0.0, 4.9910581, 0.0, 1487.7701940, 3.2804836, 0.4248467, 0.0, 1561.4404797, 0.0, 0.0, 8.6380271, 0.0, 3.7638340, 0.0, 0.0, 0.0, 0.0, 12.3146400, 0.0, 0.0, 120.0305245, 0.3565092, 0.4894787, 0.0, 0.0, 0.0, 0.0)
  diff_result <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(diff_result), 0)
})

test_that("getHierarchyInnerNodesRoot", {
  res <- mObj$getHierarchy(nodeId = NULL)
  
  expect_equal("AllFeatures", as.character(res$tree$label))
  expect_equal(5, res$tree$nchildren)
})

test_that("getRowsInnerNodes", {
  sampleId<- "CCIS98482370ST-3-0"
  resRows <- mObj$getRows(measurement = sampleId, selectedLevels = 2)  
  expected_label <- c("Euryarchaeota", "Acidobacteria", "Actinobacteria", "Bacteroidetes", 
                      "Candidatus_Saccharibacteria", "Chlorobi", "Deferribacteres", "Deinococcus_Thermus", 
                      "Firmicutes", "Fusobacteria", "Proteobacteria", "Spirochaetes", "Synergistetes", 
                      "Verrucomicrobia", "Ascomycota", "Eukaryota_noname", "Viruses_noname")
  
  intersect_result <- intersect(resRows$metadata$label, expected_label)
  expect_equal(length(intersect_result), length(expected_label))
})
