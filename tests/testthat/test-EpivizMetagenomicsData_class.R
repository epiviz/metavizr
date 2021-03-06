context("testing EpivizMetagenomicsData Class")

library(metagenomeSeq)
data(mouseData)
feature_order <- colnames(fData(mouseData))
mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)

test_that("create EpivizMetagenomicsDataClass", {
  expect_is(mObj, "EpivizMetagenomicsData")
  expect_is(mObj$.graph, "MetavizGraph")
})

test_that("getHierarchy", {
  res <- mObj$getHierarchy(nodeId = NULL)

  expect_equal("Bacteria", as.character(res$tree$label))
  expect_equal(8, res$tree$nchildren)
  #Go to next level and make sure that the nchildren are correct, also feature names match
})

test_that("getValues", {
  sampleId<- "PM1:20080107"
  res <- mObj$getValues(measurements = sampleId, start=0, end=10172, selectedLevels = 3)
  
  expected <- c(0.0, 0.0, 1812.977099, 0.0, 3.816794, 3.816794, 0.0, 1488.549618, 0.0, 3.816794, 72.519084, 187.022901, 11.450382, 76.335878, 0.0, 0.0, 0.0, 0.0, 68.702290)

  diff_result <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(diff_result), 0)
})

test_that("getRows", {
  sampleId<- "PM1:20080107"
  res <- mObj$getRows(measurement = sampleId, start=1, end=10172, selectedLevels = 3)
  expected_label <- c("Actinomycetales","Coriobacteriales", "Bifidobacteriales","Lactobacillales",
                      "Clostridiales","Erysipelotrichales","Rhizobiales","Campylobacterales",
                      "Enterobacteriales","Pasteurellales","Not_Annotated_order_Bacteria")
  
  intersect_result <- intersect(res$metadata$label, expected_label)
  expect_equal(length(intersect_result), length(expected_label))
})

#getPCA
test_that("getPCA", {
  sampleIds<- c("PM9:20071217", "PM10:20080211")
  res <- mObj$getPCA(measurement = sampleIds)
  
  expected_PM9_20071217 <- c(0.7675132, -0.6410331)
  expect_equal(expected_PM9_20071217, unname(c(res$data[[1]]$PC1, res$data[[1]]$PC2)), tolerance = .0001)
  
  expected_PM10_20080211 <- c(0.6410331, 0.7675132)
  expect_equal(expected_PM10_20080211, unname(c(res$data[[2]]$PC1, res$data[[2]]$PC2)), tolerance = .0001)
  
  expected_pca_variance_explained <- c(0.6640456, 0.4342690)
  expect_equal(expected_pca_variance_explained, res$pca_variance_explained, tolerance = .0001)
})

#getAlphaDiversity
test_that("getAlphaDiversity", {
  sampleIds<- c("PM9:20071217", "PM10:20080211")
  res <- mObj$getAlphaDiversity(measurement = sampleIds)
  
  expected_PM9_20071217 <- 4.589864
  expect_equal(expected_PM9_20071217, res$data[[1]]$alphaDiversity, tolerance = .0001)
  
  expected_PM10_20080211 <- 3.838823
  expect_equal(expected_PM10_20080211, res$data[[2]]$alphaDiversity, tolerance = .0001)
})

#searchTaxonomy
test_that("searchTaxonomy", {
  res <- mObj$searchTaxonomy(query = "bact", max_results = 10)

  expect_equal(10, length(res))
})