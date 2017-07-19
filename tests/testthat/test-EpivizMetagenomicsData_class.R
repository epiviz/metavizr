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
  res <- mObj$getValues(measurements = sampleId, start=0, end=10172, selectedLevels = 4)
  
  expected <- c(2565.789474, 3125.000000, 250.000000,6.578947,322.368421,125.000000,19.736842)
  expect_equal(nnzero(expected), nnzero(unique(res[[sampleId]])))
  
  intersect_result <- intersect(round(res[[sampleId]], digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})

test_that("getRows", {
  library(metagenomeSeq)
  data(mouseData)
  feature_order <- colnames(fData(mouseData))
  
  sampleId<- "PM1:20080107"
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData, feature_order = feature_order)
  res <- mObj$getRows(measurement = sampleId, start=1, end=10172, selectedLevels = 4)
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

test_that("getValuesOrder", {
  library(metagenomeSeq)
  library(msd16s)
  data(mouseData)
  feature_order <- colnames(fData(msd16s))[1:8]

  sampleId<- "600796"
  mObj <- metavizr:::EpivizMetagenomicsData(cumNorm(msd16s, p=.75), feature_order = feature_order)
  res <- mObj$getValues(measurements = rownames(pData(msd16s)), start=1, end=26044, selectedLevels = 4)

  expected <- c(2.5188916877,0,146.0957178841,20.1511335013,0,0,0,1319.8992443325,1201.5113350126,2544.080604534,0,513.8539042821,10.0755667506,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,496.2216624685,1365.2392947103,0,0,0,0,0,410.5793450882)

  intersect_result <- intersect(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})

test_that("getValuesFamily", {
  library(metagenomeSeq)
  library(msd16s)
  data(mouseData)
  feature_order <- colnames(fData(msd16s))[1:8]

  sampleId<- "600796"
  mObj <- metavizr:::EpivizMetagenomicsData(cumNorm(msd16s, p=.75), feature_order = feature_order)
  res <- mObj$getValues(measurements = rownames(pData(msd16s)), start=1, end=26044, selectedLevels = 5)

  expected <- c(2.5188916877,0,0,0,0,0,0,0,146.0957178841,12.5944584383,2.5188916877,5.0377833753,0,0,0,0,1319.8992443325,0,0,0,0,0,5.0377833753,2.5188916877,32.7455919395,1161.2090680101,279.59697733,0,138.5390428212,1005.0377833753,115.8690176322,5.0377833753,17.6322418136,982.3677581864,0,0,513.8539042821,10.0755667506,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,496.2216624685,1365.2392947103,0,0,0,0,0,0,410.5793450882)

  intersect_result <- intersect(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})

test_that("getValuesGenus", {
  library(metagenomeSeq)
  library(msd16s)
  data(mouseData)
  feature_order <- colnames(fData(msd16s))[1:8]

  sampleId<- "600796"
  mObj <- metavizr:::EpivizMetagenomicsData(cumNorm(msd16s, p=.75), feature_order = feature_order)
  res <- mObj$getValues(measurements = rownames(pData(msd16s)), start=1, end=26044, selectedLevels = 6)

  expected <- c(2.5188916877,0,0,0,0,0,0,0,0,0,0,0,105.7934508816,40.3022670025,0,0,12.5944584383,0,0,0,2.5188916877,5.0377833753,0,0,0,0,1319.8992443325,0,0,0,0,0,0,0,0,5.0377833753,2.5188916877,2.5188916877,30.2267002519,0,1161.2090680101,259.4458438287,20.1511335013,0,0,0,0,138.5390428212,2.5188916877,0,2.5188916877,1000,0,0,0,0,0,0,115.8690176322,0,0,5.0377833753,17.6322418136,0,62.9722921914,919.395465995,0,0,0,0,0,0,0,0,0,0,0,0,5.0377833753,508.8161209068,10.0755667506,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,2.5188916877,0,0,10.0755667506,0,95.717884131,0,367.758186398,12.5944584383,0,0,0,0,0,5.0377833753,0,0,0,0,0,0,0,1365.2392947103,0,0,0,0,0,0,0,410.5793450882)

  intersect_result <- intersect(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})

test_that("getValuesSpecies", {
  library(metagenomeSeq)
  library(msd16s)
  data(mouseData)
  feature_order <- colnames(fData(msd16s))[1:8]

  sampleId<- "600796"
  mObj <- metavizr:::EpivizMetagenomicsData(cumNorm(msd16s, p=.75), feature_order = feature_order)
  res <- mObj$getValues(measurements = rownames(pData(msd16s)), start=1, end=26044, selectedLevels = 7)

  expected <- c(0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,105.7934508816,0,0,40.3022670025,0,0,0,0,0,0,0,5.0377833753,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,2.5188916877,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1319.8992443325,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.0377833753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,2.5188916877,0,0,0,0,22.6700251889,7.556675063,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,871.5365239295,0,2.5188916877,0,0,0,0,0,0,25.1889168766,0,0,0,0,0,198.9924433249,0,0,0,0,0,2.5188916877,0,0,0,0,2.5188916877,0,0,0,0,0,2.5188916877,0,0,2.5188916877,0,0,35.2644836272,0,0,10.0755667506,0,0,0,0,5.0377833753,0,0,0,0,0,0,0,0,0,110.8312342569,0,0,0,0,0,0,0,0,73.0478589421,0,2.5188916877,0,0,0,42.8211586902,7.556675063,0,0,0,5.0377833753,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,2.5188916877,0,0,0,7.556675063,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,2.5188916877,20.1511335013,0,0,0,0,0,0,0,7.556675063,0,0,0,0,0,0,0,130.9823677582,0,0,0,0,0,0,2.5188916877,0,0,2.5188916877,0,0,0,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,115.8690176322,0,0,5.0377833753,0,0,17.6322418136,0,0,0,60.4534005038,2.5188916877,0,2.5188916877,0,0,30.2267002519,836.2720403023,0,0,2.5188916877,5.0377833753,0,0,2.5188916877,5.0377833753,0,0,0,0,0,30.2267002519,0,5.0377833753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.0377833753,0,10.0755667506,0,0,0,120.9068010076,0,55.4156171285,0,35.2644836272,0,5.0377833753,0,0,0,0,10.0755667506,0,2.5188916877,55.4156171285,45.3400503778,133.5012594458,0,2.5188916877,7.556675063,2.5188916877,17.6322418136,5.0377833753,0,7.556675063,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.5188916877,0,0,0,0,0,2.5188916877,0,0,0,10.0755667506,0,0,0,0,0,0,0,0,0,0,0,0,0,0,93.1989924433,0,2.5188916877,0,0,0,0,0,345.0881612091,15.1133501259,2.5188916877,5.0377833753,0,12.5944584383,0,0,0,0,0,0,0,0,0,5.0377833753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1322.4181360202,0,42.8211586902,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,410.5793450882)

  intersect_result <- intersect(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})