context("testing EpivizMetagenomicsInnerNodesData Class")
library(tidyr)

# From curatedMetagenomicData and MicrobiomeWorkshop Vignettes
loman_file <- system.file("extdata", "loman.RData", package = "metavizr")
load(loman_file)
loman.eset <- loman[[1]]

relative_abundance <-
  exprs(loman.eset)

number_reads <-
  loman.eset$number_reads

counts <-
  sweep(relative_abundance, 2, number_reads/100, "*")

counts <-
  round(counts)

pheno_data <-
  phenoData(loman.eset)

into_cols <-
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")

taxonomy_data <-
  data.frame(rowname = rownames(counts))

taxonomy_data <-
  separate(taxonomy_data, "rowname", into_cols, sep = "\\|", fill = "right")

rownames(counts) <-
  gsub(".+\\|", "", rownames(counts))

rownames(taxonomy_data) <-
  gsub(".+\\|", "", rownames(counts))

taxonomy_data <-
  apply(taxonomy_data, 2, function(x) {gsub("[a-z]__", "", x)})

taxonomy_data <-
  data.frame(taxonomy_data)

taxonomy_data <-
  AnnotatedDataFrame(taxonomy_data)

loman_MR <-
  metagenomeSeq::newMRexperiment(
    counts = counts,
    phenoData = pheno_data,
    featureData = taxonomy_data
  )

feature_order <- colnames(fData(loman_MR))
mObj <- metavizr:::EpivizMetagenomicsDataInnerNodes$new(loman_MR, feature_order = feature_order)

test_that("getValuesInnerNodes", {
  sampleId<- "OBK1122"
  res <- mObj$getValues(measurements = sampleId)

  expected <- c(0.000 , 0.000, 1.2948, 499.0848, 0.000, 89.1858, 0.000, 0.000, 0.000, 0.000, 7.7368, 0.000, 0.000, 0.000, 0.000)
  diff_result <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
  expect_equal(length(diff_result), 0)
})

test_that("getHierarchyInnerNodesRoot", {
  res <- mObj$getHierarchy(nodeId = NULL)

  expect_equal("AllFeatures", as.character(res$tree$label))
  expect_equal(3, res$tree$nchildren)
})

test_that("getRowsInnerNodes", {
  sampleId<- "OBK1122"
  resRows <- mObj$getRows(measurement = sampleId, selectedLevels = 2)
  expected_label <- c("Euryarchaeota", "Actinobacteria", "Bacteroidetes",
                      "Firmicutes", "Fusobacteria", "Proteobacteria",
                      "Verrucomicrobia", "Viruses_noname")

  intersect_result <- intersect(resRows$metadata$label, expected_label)
  expect_equal(length(intersect_result), length(expected_label))
})
