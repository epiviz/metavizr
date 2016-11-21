context("testing EpivizMetagenomicsData Class")

test_that("create EpivizMetagenomicsDataClass", {
  load(file.path(system.file("tests", package="metavizr"), "mouseData.RData"))
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData)
  
  expect_is(mObj, "EpivizMetagenomicsData")
  expect_is(mObj$.taxonomy, "MetavizTree")
})

test_that("getHierarchy", {
  load(file.path(system.file("tests", package="metavizr"), "mouseData.RData"))

  mObj <- metavizr:::EpivizMetagenomicsData(mouseData)
  res <- mObj$getHierarchy(nodeId = NULL)
  
  expect_equal("0-0", res$tree$id)
  expect_equal(8, res$tree$nchildren)

  load(file.path(system.file("tests", package="metavizr"), "expected_getHierarchy.RData")) 
  result <- epivizrServer::json_writer(res)
  expect_equal(expected, result)
})

test_that("getValues", {
  load(file.path(system.file("tests", package="metavizr"), "mouseData.RData"))
  sampleId<- "PM1:20080108"
  
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData)
  res <- mObj$getValues(measurement = "PM1:20080108", seqName = "metavizr", start=0, end=10000)
  
  expect_equal(15, length(res$values$values))
  
  expected <- "{\"globalStartIndex\":0,\"values\":{\"values\":[0,0,13.1004366812228,4004.36681222711,13.1004366812228,4.36681222707428,0,0,0,0,3349.34497816597,4.36681222707428,576.419213973805,96.0698689956341,109.170305676857]}}"
  result <- epivizrServer::json_writer(res)
  expect_equal(expected, result)
})

test_that("getRows", {
  load(file.path(system.file("tests", package="metavizr"), "mouseData.RData"))
  
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData)
  res <- mObj$getRows(seqName = "metavizr", start=0, end=10000)
  
  expect_equal(15, length(res$values$id))

  result <- epivizrServer::json_writer(res)
  expected <- "{\"globalStartIndex\":0,\"values\":{\"id\":[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],\"start\":[0,24,26,123,2449,2465,2521,2522,2524,2817,2847,8635,8686,9151,9605],\"end\":[24,26,123,2449,2465,2521,2522,2524,2817,2847,8635,8686,9151,9605,10024],\"metadata\":{\"colLabel\":[\"Actinomycetales\",\"Bifidobacteriales\",\"Coriobacteriales\",\"Bacteroidales\",\"NA\",\"NA\",\"Chloroplast\",\"Bacillales\",\"Lactobacillales\",\"NA\",\"Clostridiales\",\"NA\",\"Erysipelotrichales\",\"NA\",\"NA\"],\"ancestors\":[\"Actinomycetales,Actinobacteria,Actinobacteria,Bacteria\",\"Bifidobacteriales,Actinobacteria,Actinobacteria,Bacteria\",\"Coriobacteriales,Actinobacteria,Actinobacteria,Bacteria\",\"Bacteroidales,Bacteroidetes,Bacteroidetes,Bacteria\",\"NA,Bacteroidetes,Bacteroidetes,Bacteria\",\"NA,NA,Bacteroidetes,Bacteria\",\"Chloroplast,Cyanobacteria,Cyanobacteria,Bacteria\",\"Bacillales,Bacilli,Firmicutes,Bacteria\",\"Lactobacillales,Bacilli,Firmicutes,Bacteria\",\"NA,Bacilli,Firmicutes,Bacteria\",\"Clostridiales,Clostridia,Firmicutes,Bacteria\",\"NA,Clostridia,Firmicutes,Bacteria\",\"Erysipelotrichales,Erysipelotrichi,Firmicutes,Bacteria\",\"NA,NA,Firmicutes,Bacteria\",\"NA,NA,NA,Bacteria\"],\"lineage\":[\"3-0,2-0,1-0,0-0\",\"3-18,2-0,1-0,0-0\",\"3-1a,2-0,1-0,0-0\",\"3-7b,2-7b,1-7b,0-0\",\"3-991,2-7b,1-7b,0-0\",\"3-9a1,2-9a1,1-7b,0-0\",\"3-9d9,2-9d9,1-9d9,0-0\",\"3-9da,2-9da,1-9da,0-0\",\"3-9dc,2-9da,1-9da,0-0\",\"3-b01,2-9da,1-9da,0-0\",\"3-b1f,2-b1f,1-9da,0-0\",\"3-21bb,2-b1f,1-9da,0-0\",\"3-21ee,2-21ee,1-9da,0-0\",\"3-23bf,2-23bf,1-9da,0-0\",\"3-2585,2-2585,1-2585,0-0\"],\"OTU\":[\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\"],\"genus\":[\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\"],\"family\":[\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\",\"<NA>\"],\"order\":[\"Actinomycetales\",\"Bifidobacteriales\",\"Coriobacteriales\",\"Bacteroidales\",\"NA\",\"NA\",\"Chloroplast\",\"Bacillales\",\"Lactobacillales\",\"NA\",\"Clostridiales\",\"NA\",\"Erysipelotrichales\",\"NA\",\"NA\"],\"class\":[\"Actinobacteria\",\"Actinobacteria\",\"Actinobacteria\",\"Bacteroidetes\",\"Bacteroidetes\",\"NA\",\"Cyanobacteria\",\"Bacilli\",\"Bacilli\",\"Bacilli\",\"Clostridia\",\"Clostridia\",\"Erysipelotrichi\",\"NA\",\"NA\"],\"phylum\":[\"Actinobacteria\",\"Actinobacteria\",\"Actinobacteria\",\"Bacteroidetes\",\"Bacteroidetes\",\"Bacteroidetes\",\"Cyanobacteria\",\"Firmicutes\",\"Firmicutes\",\"Firmicutes\",\"Firmicutes\",\"Firmicutes\",\"Firmicutes\",\"Firmicutes\",\"NA\"],\"superkingdom\":[\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\",\"Bacteria\"],\"label\":[\"Actinomycetales\",\"Bifidobacteriales\",\"Coriobacteriales\",\"Bacteroidales\",\"NA\",\"NA\",\"Chloroplast\",\"Bacillales\",\"Lactobacillales\",\"NA\",\"Clostridiales\",\"NA\",\"Erysipelotrichales\",\"NA\",\"NA\"]},\"index\":[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]}}"
  expect_equal(expected, result)
})

test_that("searchTaxonomy", {
  load(file.path(system.file("tests", package="metavizr"), "mouseData.RData"))
  mObj <- metavizr:::EpivizMetagenomicsData(mouseData)
  res <- mObj$searchTaxonomy(query = "bact", max_results = 10)
  expect_equal(10, length(res))
})
