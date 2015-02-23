setMethod("register", "MRexperiment", function(object, ...) {
  return(EpivizMetagenomicsData$new(object=object, ...))
})
