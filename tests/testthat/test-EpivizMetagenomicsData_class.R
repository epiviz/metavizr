context("testing EpivizMetagenomicsData Class")

test_that("create EpivizMetagenomicsDataClass", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  
  expect_is(mObj, "EpivizMetagenomicsData")
  expect_is(mObj$.graph, "MetavizGraph")
})

test_that("getHierarchy", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getHierarchy(nodeId = NULL)
  
  expect_equal("Bacteria", as.character(res$tree$label))
  expect_equal(8, res$tree$nchildren)
  #Go to next level and make sure that the nchildren are correct, also feature names match
})

test_that("getValues", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  sampleId<- "PM1:20080107"
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getValues(measurement = sampleId, start=0, end=10172, selectedLevels = 3)
  
  expected <- c(2565.789474, 3125.000000, 2250.000000,6.578947,322.368421,125.000000,19.736842)
  expect_equal(nnzero(expected), nnzero(unique(res[[sampleId]])))
  
  intersect_result <- intersect(round(res[[sampleId]], digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})

test_that("getRows", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getRows(start=0, end=10172, selectedLevels = 3)
  expected_label <- c("Actinomycetales","Coriobacteriales", "Bifidobacteriales","Lactobacillales",
                      "Clostridiales","Erysipelotrichales","Rhizobiales","Campylobacterales",
                      "Enterobacteriales","Pasteurellales","Not_Annotated_order")
  
  intersect_result <- intersect(res$metadata$label, expected_label)
  expect_equal(length(intersect_result), length(expected_label))
})

#getPCA

#getAlphaDiversity

#searchTaxonomy
test_that("searchTaxonomy", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))

  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$searchTaxonomy(query = "bact", max_results = 10)

  expect_equal(10, length(res))
})
