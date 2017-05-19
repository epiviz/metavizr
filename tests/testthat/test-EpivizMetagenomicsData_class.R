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
  res <- mObj$getValues(measurement = sampleId, start=0, end=10172, selectedLevels = 4)
  
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
  res <- mObj$getRows(start=0, end=10172, selectedLevels = 4)
  expected_label <- c("Actinomycetales","Coriobacteriales", "Bifidobacteriales","Lactobacillales",
                      "Clostridiales","Erysipelotrichales","Rhizobiales","Campylobacterales",
                      "Enterobacteriales","Pasteurellales","Not_Annotated_order")
  
  intersect_result <- intersect(res$metadata$label, expected_label)
  expect_equal(length(intersect_result), length(expected_label))
})

#getPCA
test_that("getPCA", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  sampleIds<- c("PM9:20071217", "PM10:20080211")
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getPCA(measurement = sampleIds)
  
  expected_PM9_20071217 <- c(0.7798891, -0.6259177)
  expect_equal(expected_PM9_20071217, unname(c(res$data[[1]]$PC1, res$data[[1]]$PC2)), tolerance = .0001)
  
  expected_PM10_20080211 <- c(0.6259177, 0.7798891)
  expect_equal(expected_PM10_20080211, unname(c(res$data[[2]]$PC1, res$data[[2]]$PC2)), tolerance = .0001)
  
  expected_pca_variance_explained <- c(0.7798627, 0.5189621)
  expect_equal(expected_pca_variance_explained, res$pca_variance_explained, tolerance = .0001)
})

#getAlphaDiversity
test_that("getPCA", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  sampleIds<- c("PM9:20071217", "PM10:20080211")
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getAlphaDiversity(measurement = sampleIds)
  
  expected_PM9_20071217 <- 4.589864
  expect_equal(expected_PM9_20071217, res$data[[1]]$alphaDiversity, tolerance = .0001)
  
  expected_PM10_20080211 <- 3.838823
  expect_equal(expected_PM10_20080211, res$data[[2]]$alphaDiversity, tolerance = .0001)
})

#searchTaxonomy
test_that("searchTaxonomy", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))

  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$searchTaxonomy(query = "bact", max_results = 10)

  expect_equal(10, length(res))
})
