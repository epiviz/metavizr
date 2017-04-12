context("testing EpivizMetagenomicsData Class")

test_that("create EpivizMetagenomicsDataClass", {
  skip_if_not_installed("msd16s")
  require(msd16s)
  feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
  mObj <- metavizr:::EpivizMetagenomicsData(msd16s, feature_order = feature_order)
  
  expect_is(mObj, "EpivizMetagenomicsData")
  expect_is(mObj$.graph, "MetavizGraph")
})

test_that("getHierarchy", {
  skip_if_not_installed("msd16s")
  require(msd16s)

  feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
  
  mObj <- metavizr:::EpivizMetagenomicsData(msd16s, feature_order = feature_order)
  res <- mObj$getHierarchy(nodeId = NULL)
  
  expect_equal("Bacteria", as.character(res$tree$label))
  expect_equal(8, res$tree$nchildren)
  #Go to next level and make sure that the nchildren are correct, also feature names match
})

test_that("getValues", {
  skip_if_not_installed("msd16s")
  require(msd16s)

  sampleId<- "100259"
  feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
  mObj <- metavizr:::EpivizMetagenomicsData(msd16s, feature_order = feature_order)
  res <- mObj$getValues(measurement = sampleId, start=0, end=26044, selectedLevels = 4)
  
  expected <- c(22.9681978799, 28.2685512367, 2446.9964664311,560.0706713781,1010.6007067138,8.8339222615,104.2402826855,15.9010600707,1.7667844523,42.4028268551,123.6749116608,100.7067137809,3503.5335689046)
  expect_equal(nnzero(expected), nnzero(res[[sampleId]]))
  
  intersect_result <- intersect(round(res[[sampleId]], digits=3), round(expected,digits=3))
  expect_equal(length(intersect_result), length(expected))
})

test_that("getRows", {
  skip_if_not_installed("msd16s")
  require(msd16s)

  feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
  
  mObj <- metavizr:::EpivizMetagenomicsData(msd16s, feature_order = feature_order)
  res <- mObj$getRows(start=0, end=26044, selectedLevels = 4)
  expected_label <- c("Actinomycetales","Coriobacteriales","Bacteroidales","Lactobacillales","Clostridiales","Erysipelotrichales","Selenomonadales","Fusobacteriales","Rhizobiales","Campylobacterales","Enterobacteriales","Pasteurellales","Not_Annotated_order")
  
  intersect_result <- intersect(res$metadata$label, expected_label)
  expect_equal(length(intersect_result), length(expected_label))
  #Check global start index
})

#getPCA

#getAlphaDiversity

#searchTaxonomy
test_that("searchTaxonomy", {
  skip_if_not_installed("msd16s")
  require(msd16s)

  feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
  mObj <- metavizr:::EpivizMetagenomicsData(msd16s, feature_order = feature_order)
  res <- mObj$searchTaxonomy(query = "bact", max_results = 10)
  expect_equal(10, length(res))
})
