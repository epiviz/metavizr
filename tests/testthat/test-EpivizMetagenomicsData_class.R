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
  res <- mObj$getValues(measurements = sampleId, start=0, end=10172, selectedLevels = 3)
  
  expected <- c(0.0, 0.0, 1812.977099, 0.0, 3.816794, 3.816794, 0.0, 1488.549618, 0.0, 3.816794, 72.519084, 187.022901, 11.450382, 76.335878, 0.0, 0.0, 0.0, 0.0, 68.702290)

  diff_reuslt <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(diff_reuslt), 0)
})

test_that("getRows", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  sampleId<- "PM1:20080107"
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getRows(measurement = sampleId, start=1, end=10172, selectedLevels = 3)
  expected_label <- c("Actinomycetales","Coriobacteriales", "Bifidobacteriales","Lactobacillales",
                      "Clostridiales","Erysipelotrichales","Rhizobiales","Campylobacterales",
                      "Enterobacteriales","Pasteurellales","Not_Annotated_order_Bacteria")
  
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
  
  expected_PM9_20071217 <- c(0.7675132, -0.6410331)
  expect_equal(expected_PM9_20071217, unname(c(res$data[[1]]$PC1, res$data[[1]]$PC2)), tolerance = .0001)
  
  expected_PM10_20080211 <- c(0.6410331, 0.7675132)
  expect_equal(expected_PM10_20080211, unname(c(res$data[[2]]$PC1, res$data[[2]]$PC2)), tolerance = .0001)
  
  expected_pca_variance_explained <- c(0.6640456, 0.4342690)
  expect_equal(expected_pca_variance_explained, res$pca_variance_explained, tolerance = .0001)
})

#getAlphaDiversity
test_that("getPCA", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  sampleIds<- c("PM9:20071217", "PM10:20080211")
  mObj <- metavizr:::EpivizMetagenomicsData(cumNorm(mouseData, p =.75), feature_order = feature_order)
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

test_that("getValuesOrder", {
  library(metagenomeSeq)
  library(msd16s)
  feature_order <- colnames(fData(msd16s))[1:8]

  sampleId<- "600765"
  mObj <- metavizr:::EpivizMetagenomicsData(cumNorm(msd16s, p=.75), feature_order = feature_order)
  res <- mObj$getValues(measurements = rownames(pData(msd16s)), start=1, end=26044, selectedLevels = 4)

  expected <- c(0.9293680297, 0.0, 5.5762081784, 4382.8996282528, 0.0, 0.0, 0.0, 0.0, 0.0, 398.6988847584, 2003.717472119, 1.8587360595, 24.1635687732, 64.126394052, 0.0, 0.0, 0.0, 0.0, 0.0, 8.3643122677, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0111524164, 0.9293680297, 0.0, 0.0, 0.0, 0.0, 0.0, 1156.1338289963)

  diff_reuslt <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(diff_reuslt), 0)
})

test_that("getValuesInnerNodes", {
  library(curatedMetagenomicData)
  
  # From curatedMetagenomicData and MicrobiomeWorkshop Vignettes
  zeller.eset = ZellerG_2014.metaphlan_bugs_list.stool()
  zeller_MR <- ExpressionSet2MRexperiment(zeller.eset)
  
  feature_order <- colnames(fData(zeller_MR))
  
  sampleId<- "CCIS98482370ST-3-0"
  mObj <- metavizr:::InnerNodesEpivizMetagenomicsData$new(zeller_MR, feature_order = feature_order)
  res <- mObj$getValues(measurements = sampleId)
  
  expected <- c(0.0 , 34.1386632, 0.0, 4.9910581, 0.0, 1487.7701940, 3.2804836, 0.4248467, 0.0, 1561.4404797, 0.0, 0.0, 8.6380271, 0.0, 3.7638340, 0.0, 0.0, 0.0, 0.0, 12.3146400, 0.0, 0.0, 120.0305245, 0.3565092, 0.4894787, 0.0, 0.0, 0.0, 0.0)
  diff_result <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(diff_result), 0)
})