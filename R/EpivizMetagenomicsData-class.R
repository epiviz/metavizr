#' Data container for MRexperiment objects
#' 
#' Used to serve metagenomic data (used in e.g., icicle plots and heatmaps). Wraps
#' \code{\link[metagenomeSeq]{MRexperiment-class}} objects.
#' @importClassesFrom epivizrData EpivizData
#' @importFrom methods new
#' @import RNeo4j
#' @exportClass EpivizMetagenomicsData
EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains = "EpivizData",
  fields = list(
    .taxonomy = "MetavizTree",
    .levels = "character",
    .maxDepth = "numeric",
    .aggregateAtDepth = "numeric",
    .lastRootId = "ANY",
    .feature_order = "character",

    .counts = "ANY",
    .sampleAnnotation = "ANY",

    .minValue = "ANY",
    .maxValue = "ANY",
    .aggregateFun = "ANY",
    .valuesAnnotationFuns = "ANY",

    .lastRequestRanges = "list",
    .lastLeafInfos = "list",
    .lastSelectionTypes = "list",
    .lastValues = "list",
    .maxHistory = "numeric",
    .json_query = "ANY"
  ),
  methods=list(
    # TODO use and check columns
    initialize=function(object, 
                        columns=NULL,
                        control=metavizControl(), 
                        feature_order=NULL, 
                        ...) {

      # Initialize parameters used here
      aggregateAtDepth <- control$aggregateAtDepth
      maxDepth         <- control$maxDepth
      maxHistory       <- control$maxHistory
      maxValue         <-  control$maxValue
      minValue         <-  control$minValue
      aggregateFun     <-  control$aggregateFun
      valuesAnnotationFuns <- control$valuesAnnotationFuns

      if (is.character(aggregateAtDepth)) { 
        aggregateAtDepth <- assignValues(object, aggregateAtDepth) 
      }
      
      if (is.character(maxDepth)) { 
        maxDepth <- assignValues(object, maxDepth) 
      }
      log <- control$log
      norm <- control$norm
      
      # validate MRexperiment object
      MRExpCheck <- validateObject(object)
      
      if (!MRExpCheck) {
        stop("Incompatible MRexperiment objects")
      } else {
        message("MRExperiment Object validated... PASS")
      }

      if(is.null(feature_order)) {
        .self$.feature_order = colnames(fData(object))
      }
      else {
        .self$.feature_order <- feature_order
      }

      
      # TODO: Some type checking
      .self$.taxonomy <- buildMetavizTree(object, feature_order)
      .self$.levels <- .self$.taxonomy$levels()
      .self$.maxDepth <- maxDepth
      .self$.aggregateAtDepth <- aggregateAtDepth
      .self$.lastRootId <- .self$.taxonomy$root()$id()

      taxonomy_table <- .self$.taxonomy$taxonomyTable()
      .self$.counts <- MRcounts(object[rownames(taxonomy_table),], norm=norm, log=log)

      # TODO: Make this consistent with the aggregateFun
      if (is.null(minValue)) {
        minValue <- log2(min(.self$.counts) + 1)
      }
      .self$.minValue <- minValue

      if (is.null(maxValue)) {
        maxValue <- log2(max(.self$.counts) + 1)
      }
      .self$.maxValue <- maxValue
      .self$.aggregateFun <- aggregateFun
      .self$.valuesAnnotationFuns <- valuesAnnotationFuns

      .self$.sampleAnnotation <- pData(object)

      .self$.maxHistory <- maxHistory
      .self$.lastRequestRanges <- list()
      .self$.lastLeafInfos <- list()
      .self$.lastValues <- list()
      featureSelection = control$featureSelection
      
      if(!is.null(featureSelection)){
        .self$featureSelection(featureSelection$featureNames, featureSelection$featureOrder, featureSelection$featureLevel, featureSelection$selectionType)
      }
      else if (.self$.aggregateAtDepth >= 0) {
        nodesAtDepth <- .self$.taxonomy$nodesAtDepth(.self$.aggregateAtDepth)
        node_ids <- sapply(nodesAtDepth, function(node) node$id())
        selection <- lapply(node_ids, function(node_id) SelectionType$NODE)
        names(selection) <- node_ids
        .self$.taxonomy$updateSelection(selection)
      }

      callSuper(object=object, ...)
    },

    .taxonomyLevels=function(exp) {
      .self$.levels
    },

    .getSelectedLeaves=function(start, end) {
      ret <- NULL
      if (length(.self$.lastRequestRanges) > 0) {
        for (i in rev(seq_along(.self$.lastRequestRanges))) {
          if (.self$.lastRequestRanges[[i]]$start == start || 
              .self$.lastRequestRanges[[i]]$end == end) {
            ret <- .lastLeafInfos[[i]]
            break
          }
        }
      }
      if (is.null(ret)) {
        requestRange <- list(start=start, end=end)
        ret <- Ptr$new(.self$.taxonomy$selectedLeaves(start, end))

        if (.self$.maxHistory > 0) {
          .self$.lastLeafInfos <- c(.self$.lastLeafInfos, ret)
          .self$.lastRequestRanges <- c(.self$.lastRequestRanges, list(requestRange))
          .self$.lastValues <- c(.self$.lastValues, list(NULL))

          if (length(.self$.lastRequestRanges) > .self$.maxHistory) {
            .self$.lastRequestRanges <- .self$.lastRequestRanges[2:(.self$.maxHistory+1)]
            .self$.lastLeafInfos <- .self$.lastLeafInfos[2:(.self$.maxHistory+1)]
            .self$.lastValues <- .self$.lastValues[2:(.self$.maxHistory+1)] # Not yet a value
          }
        }
      }

      ret$.
    },

    .getSelectedValues=function(measurement, start, end) {
      leafInfos <- .self$.getSelectedLeaves(start, end)
      index <- -1
      values <- NULL
      if (length(.self$.lastRequestRanges) > 0) {
        for (i in rev(seq_along(.self$.lastRequestRanges))) {
          if (.self$.lastRequestRanges[[i]]$start == start || 
              .self$.lastRequestRanges[[i]]$end == end) {
            if (!is.null(.self$.lastValues[[i]])) {
              values <- .self$.lastValues[[i]]
            } else {
              index <- i
            }
            break
          }
        }
      }

      if (is.null(values)) {
        values <- Ptr$new(list(
          values=unname(lapply(leafInfos, function(info) {
            lind <- info$node$leafIndex()
            ind = seq(lind+1, lind+info$node$nleaves())
            .self$.aggregateFun(.self$.counts[ind,, drop=FALSE])
          }))
        ))
        
        if (!is.null(.self$.valuesAnnotationFuns)) {
          for (anno in names(.self$.valuesAnnotationFuns)) {
            fun <- .self$.valuesAnnotationFuns[[anno]]
            values$.[[anno]] <- unname(lapply(leafInfos, function(info) {
              lind <- info$node$leafIndex()
              ind <- seq(lind+1, lind+info$node$nleaves())
              fun(.self$.counts[ind,, drop=FALSE])
            }))
          }
        }
      }
      
      if (index > 0) {
        .self$.lastValues[[index]] <- values
      }

      lapply(values$., function(vals) {
        lapply(vals, function(v) {
          if (length(dim(v)) > 0) { return(v[, measurement]) }
          return(v[[measurement]])
        })
      })
    },

    featureSelection=function(featureNames, featureOrder, featureLevel, selectionType){
      if(is.null(featureOrder)){
        featureOrder = .self$.feature_order
      }
      featureDepth <- (which(featureOrder == featureLevel) -1)
      
      nodes <- .self$.taxonomy$nodesAtDepth(depth=featureDepth)
      node_names <- sapply(1:length(nodes), function(i) {nodes[[i]]$.name})
      node_indices <- match(featureNames, node_names)
      node_ids <- sapply(node_indices, function(i) {nodes[[i]]$.id})
      node_all_ids <- sapply(1:length(nodes), function(i){nodes[[i]]$.id})
      
      selections <- list()
      selection_vals <- as.numeric(!is.na(match(node_names, featureNames)))
      selection_vals[which(selection_vals == 1)] <- selectionType
      selection_keys <- node_all_ids
      
      nodeSelections <- list()
      
      for(i in 1:length(nodes)) {
        # if(selection_vals[i] == 0) {
          selections[[selection_keys[i]]] <- selection_vals[i] 
        # }
          
        if(selection_vals[i] == 0) {
          nodeSelections[[selection_keys[i]]] <- selection_vals[i] 
        }
      }
      
      .self$clearSelection()
      .self$.lastSelectionTypes <- nodeSelections
      .self$.taxonomy$updateSelection(selection = selections)
      .self$.aggregateAtDepth <- featureDepth
    },

    update=function(newObject, ...) {
      # TODO
      callSuper(newObject, ...)
    },
    plot=function(...) {
      # TODO
    }
  )
)

# Data analysis features
EpivizMetagenomicsData$methods(
  nleaves=function() {
    if (is.null(dim(.self$.counts))) {
      if (length(.self$.counts) > 0) { return(1) }
      return(0)
    }

    nrow(.self$.counts)
  },
  nmeasurements=function() {
    if (is.null(dim(.self$.counts))) {
      return(length(.self$.counts))
    }

    ncol(.self$.counts)
  },

  nlevels=function() {
    length(.self$.levels)
  },
  clearSelection=function(){
    #if(!is.null(selectedLevels)) {
      # remove current selectionTypes
      nodeList <- list()
      for (node in names(.self$.taxonomy$.selectionTypes$.)) {
        nodeList[[node]] <- 1
      }
      
      .self$.lastSelectionTypes <- list()
      
      .self$.taxonomy$updateSelection(nodeList)
      
      if (!is.null(.self$.mgr)) {
        .self$.mgr$.clear_datasourceGroup_cache(.self)
      }
    #}
  },
  taxonomyTable=function() { .self$.taxonomy$taxonomyTable() },
  calcNodeId=function(rowIndex, colIndex) { .self$.taxonomy$calcNodeId(rowIndex, colIndex) },
  node=function(nodeId) { .self$.taxonomy$node(nodeId) },
  parent=function(node) { .self$.taxonomy$parent(node) },
  siblings=function(node) { .self$.taxonomy$siblings(node) },

  changeAggregation=function(nodeId, aggregationType) {
    selection <- list()
    selection[[nodeId]] <- aggregationType
    .self$.taxonomy$updateSelection(selection)
    if (!is.null(.self$.mgr)) {
      .self$.mgr$.clear_datasourceGroup_cache(.self)
    }
  },
  
  changeAggregationAtDepth=function(depth, aggregationType) {
    if (depth < 0) { return() }
    nodesAtDepth <- .self$.taxonomy$nodesAtDepth(depth)
    selection <- list()
    for (node in nodesAtDepth) {
      nodeId <- node$.id
      selection[[nodeId]] <- aggregationType
    }
    .self$.taxonomy$updateSelection(selection)
    if (!is.null(.self$.mgr)) {
      .self$.mgr$.clear_datasourceGroup_cache(.self)
    }
  }
)

# Epiviz Websockets Protocol
EpivizMetagenomicsData$methods(
  get_default_chart_type = function() { "epiviz.ui.charts.tree.Icicle" },
  get_measurements=function() {
    out <- lapply(colnames(.self$.counts), function(sample) {
      epivizrData:::EpivizMeasurement(id=sample,
           name=sample,
           type="feature",
           datasourceId=.self$.id,
           datasourceGroup=.self$.id,
           defaultChartType="heatmap",
           annotation=as.list(.self$.sampleAnnotation[sample,]),
           minValue=.minValue,
           maxValue=.maxValue,
           metadata=c(rev(.levels), "colLabel", "ancestors", "lineage", "label"))
    })
    out
  },
  getHierarchy=function(nodeId) {
    
    # clear last request ranges and values
    .self$.lastRequestRanges <- list()
    .self$.lastLeafInfos <- list()
    .self$.lastValues <- list()
    
    root <- NULL
    if (missing(nodeId) || is.null(nodeId)) { 
      root <- .self$.taxonomy$root() 
    } else {
      root <- .self$.taxonomy$parent(.self$.taxonomy$node(nodeId))
      if (is.null(root)) { root <- .self$.taxonomy$root() }
    }
    .self$.lastRootId <- nodeId
    
    resp <- list()
    resp[['tree']] <- root$raw(maxDepth=.self$.maxDepth)
    resp[['nodeSelectionTypes']] <- .self$.lastSelectionTypes
    resp[['selectionLevel']] <- .self$.aggregateAtDepth
  
    return(resp)
  },
  propagateHierarchyChanges=function(selection, order, selectedLevels) {
    
    # clear last requests and values
    .self$.lastRequestRanges <- list()
    .self$.lastLeafInfos <- list()
    .self$.lastValues <- list()
    
    if (missing(selection) && missing(order) && missing(selectedLevels)) { 
      return(getHierarchy(.self$.lastRootId)) 
    }
    
    if(!is.null(selectedLevels) && length(selectedLevels) != 0) {
      # remove current selectionTypes
      nodeList <- list()
      for (node in names(.self$.taxonomy$.selectionTypes$.)) {
        nodeList[[node]] <- 1
      }
      
      .self$.taxonomy$updateSelection(nodeList)
      
      # add new selectedLevel nodes
      for (level in names(selectedLevels)) {
        nodesAtLevel <- .self$.taxonomy$nodesAtDepth(as.numeric(level))
        nodeList <- list()
        for (node in nodesAtLevel) {
          nodeList[[node$id()]] <- selectedLevels[[level]]
        }
        .self$.taxonomy$updateSelection(nodeList)
        .self$.aggregateAtDepth <- as.numeric(level)
      }
    }

    if (!missing(selection)) {
      .self$.lastSelectionTypes <- selection
      .self$.taxonomy$updateSelection(.lastSelectionTypes)
    }

    if (!missing(order)) {
      .self$.taxonomy$updateOrder(order)
    }

    .self$.lastRequestRanges <- list()
    .self$.lastLeafInfos <- list()
    .self$.lastValues <- list()

    if (!is.null(.self$.mgr)) {
      .self$.mgr$.clear_datasourceGroup_cache(.self)
    }
    
    getHierarchy(.self$.lastRootId)
  },
  getRows=function(seqName, start, end, metadata) {
    leafInfos <- .self$.getSelectedLeaves(start, end)
    leafAncestors <- lapply(leafInfos, function(info) { 
      .self$.taxonomy$ancestors(info$node) 
    })

    leafTaxonomies <- list()
    for (i in seq_along(.self$.levels)) {
      leafTaxonomies[[.self$.levels[[i]]]] <- list()
      for (j in seq_along(leafInfos)) {
        if (i > length(leafAncestors[[j]])) { next }
        leafTaxonomies[[.self$.levels[[i]]]][[j]] = leafAncestors[[j]][[i]]
      }
    }

    ret = list(
      id = sapply(leafInfos, function(info) { info$realNodesBefore }),
      start=sapply(leafInfos, function(info) { info$start }),
      end=sapply(leafInfos, function(info) { info$start + info$node$nleaves()}),
      metadata = c(list(
        colLabel = sapply(leafInfos, function(info) { info$node$name() }),
        ancestors = sapply(leafAncestors, function(ancestors) { 
          paste(lapply(rev(ancestors), function(node) { node$name() }), collapse=",") 
        }), # TODO: Use tree .ancestryByDepth
        lineage = sapply(leafAncestors, function(ancestors) { 
          paste(lapply(rev(ancestors), function(node) { node$id() }), collapse=",") 
        })
      ), sapply(rev(.self$.levels), function(level) {
        r <- list(lapply(leafTaxonomies[[level]], function(node) { 
          if (is.null(node)) { return("<NA>") } 
          node$name() 
        }))
        
        if (length(r[[1]]) == 0) {
          r[[1]] <- lapply(seq_along(leafInfos), function(i) { "<NA>" })
        }
        return(r)
      }))
    )

    globalStartIndex = NULL
    if (length(leafInfos) > 0) { globalStartIndex = leafInfos[[1]]$realNodesBefore }
    if (length(leafInfos) == 1) {
      ret$index <- list(ret$id)
      ret$id <- list(ret$id)
      ret$start <- list(ret$start)
      ret$end <- list(ret$end)
      ret$metadata$colLabel <- list(ret$metadata$colLabel)
      ret$metadata$label <- list(ret$metadata$colLabel)
      ret$metadata$ancestors <- list(ret$metadata$ancestors)
      ret$metadata$lineage <- list(ret$metadata$lineage)
    }
    else {
      ret$index <- ret$id
      ret$metadata$label <- ret$metadata$colLabel
    }

    return(list(
      globalStartIndex = globalStartIndex,
      values = ret
    ))
  },
  getValues=function(measurement, seqName, start, end) {
    leafInfos <- .self$.getSelectedLeaves(start, end)
    globalStartIndex <- NULL
    
    if (length(leafInfos) > 0) { 
      globalStartIndex <- leafInfos[[1]]$realNodesBefore 
    }
    
    ret <- list(
      globalStartIndex = globalStartIndex,
      values = .self$.getSelectedValues(measurement, start, end)
    )
    return(ret)
  },
  getCombined=function(measurements, 
                       seqName, start, end, 
                       order, nodeSelection, selectedLevels) {
    
    # clear last request ranges and values
    .self$.lastRequestRanges <- list()
    .self$.lastLeafInfos <- list()
    .self$.lastValues <- list()
    
    # update selectedLevels on taxonomy tree
    if(!is.null(selectedLevels)) {
      
      # remove current selectionTypes
      nodeList <- list()
      for (node in names(.self$.taxonomy$.selectionTypes$.)) {
        nodeList[[node]] <- 1
      }
      
      .self$.taxonomy$updateSelection(nodeList)
      
      for (level in names(selectedLevels)) {
        
        nodesAtLevel = .self$.taxonomy$nodesAtDepth(as.numeric(level))
        
        nodeList = list()
        for (node in nodesAtLevel) {
          nodeList[[node$id()]] = selectedLevels[[level]]
        }
        
        .self$.taxonomy$updateSelection(nodeList)
        
      }
    }
    
    # update node selections types to metaviztree
    if(!is.null(nodeSelection)) {
      .self$.taxonomy$updateSelection(nodeSelection)
    }
    
    # update order
    if(!is.null(order)) {
      .self$.taxonomy$updateOrder(order)
    }
    
    # get values
    data_columns = list()
    for (m in measurements) {
      data_columns[[m]] = .self$getValues(m, '', start, end)$values$values
    }
    
    # get features for the given range
    data_rows = .self$getRows('', start, end)
    
    # convert to resp format
    result = list()
    result[['cols']] = data_columns
    result[['rows']] = data_rows$values
    result[['globalStartIndex']] = data_rows$globalStartIndex
    
    return(result)
  },
  searchTaxonomy=function(query, max_results) {
    
    results = list()
    nCount <- 1
    for (i in seq_along(.self$.levels)) {
      nodesAtLevel <- .self$.taxonomy$nodesAtDepth(i)
      for(j in seq_along(nodesAtLevel)) {
        node <- nodesAtLevel[[j]]
        if(length(grep(query, node$name(), value = TRUE)) != 0) {
          results[[nCount]] <- list("gene" = node$name(), "start" = node$leafIndex(),
                                    "end" = (node$leafIndex() + node$nleaves()), 
                                    "seqName" = "metavizr", "nodeId" = node$id(), "level" = .self$.levels[i]
          )
          nCount = nCount+1
          
          if(length(results) >= max_results) {break;}
        }  
      }
      if(length(results) >= max_results) { break;}
    }
    
    results <- unname(results, force=TRUE)
    return(results)
  }
)

 EpivizMetagenomicsData$methods(
   # TODO: do we need these MySQL functions still?
   toMySQLDb=function(con, colLabel=NULL) {
     "Save an MRexperiment object to a MySQL database."

     # TODO: Formal log
     cat("Saving column data...")
     .saveColData(con, colLabel)
     cat("Done\n")

     cat("Saving row data...")
     .saveRowData(con)
     cat("Done\n")

     cat("Saving hierarchy...")
     .saveHierarchy(con)
     cat("Done\n")

     cat("Saving counts...")
     .saveValues(con)
     cat("Done\n")

     cat("Saving levels...")
     .saveLevels(con)
     cat("Done\n")

     res <- dbSendQuery(con, "RENAME TABLE `meta_values` to `values`;")
     dbCommit(conn = con)

     cat("Saving Data Matrix...")
     .saveMatrix(con)
     cat("Done\n")

     dbDisconnect(con)
     #odbcClose(con)

   },
  .getFieldTypes = function(data) {
    colNames = colnames(data)
    colTypes = sapply(colNames, function(x) dbDataType(RMySQL::MySQL(), data[[x]]))
    colTypes = gsub("text", "VARCHAR(255)", colTypes)
    colTypes = gsub("double", "BIGINT", colTypes)
    
    tableTypes = list()
    
    for(i in seq(length(colTypes))) {
      tableTypes[[colNames[i]]] = colTypes[i]
    }
    
    tableTypes[["row_names"]] = "VARCHAR(255)"
    
    tableTypes
  },
  .saveMatrix = function(con) {
    
    locationCols = c('index', 'partition', 'start', 'end')
    
    sql_cols = paste0("SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA` = '", dbGetInfo(con)$dbname , "' AND `TABLE_NAME`='row_data';")
    
    resultSet = dbSendQuery(con, sql_cols)
    data = fetch(resultSet, n=-1)
    
    metadata_cols = setdiff(data$COLUMN_NAME, locationCols)
    metadata_cols = sapply(metadata_cols, function(x) paste0("row_data.",x))
    
    fields = c("row_data.index as row", "row_data.start", "row_data.end", "row_data.partition",
               "hierarchy.lineagelabel", "hierarchy.lineage", "hierarchy.depth",
               "`values`.val",
               "col_data.index as col", "col_data.id as measurement")
    
    fields = paste(c(fields, metadata_cols), collapse=", ")
    
    sql_data_matrix = paste0(
      "CREATE TABLE data_matrix ",
      "SELECT ",
      fields,
      " FROM (row_data JOIN hierarchy ON row_data.id = hierarchy.id) ",
      " JOIN col_data ",
      " JOIN `values` ON (row_data.index = `values`.row AND col_data.index = `values`.col) ",
      " ORDER BY row_data.index ASC; "
    )
    
    dbCommit(conn = con)
    
    dbSendQuery(con, sql_data_matrix)
    dbSendQuery(con, "ALTER TABLE `data_matrix` ENGINE=MyISAM;")
    dbSendQuery(con, "ALTER TABLE `data_matrix` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    dbSendQuery(con, "ALTER TABLE `data_matrix` ADD INDEX `row_idx` (`row` ASC);")
    dbSendQuery(con, "ALTER TABLE `data_matrix` ADD INDEX `measurement_idx` (`measurement` ASC);")
    
    dbCommit(conn = con)
  },
  .saveColData=function(con, labelCol=NULL) {
    #odbcSetAutoCommit(con, autoCommit = FALSE)
    
    cols = .sampleAnnotation
    if (!is.null(labelCol)) {
      cols$label = cols[[labelCol]]
    }
    cols$index = seq(dim(cols)[1]) - 1

    tableFieldTypes = .getFieldTypes(cols)
    
    dbSendQuery(con, "DROP TABLE IF EXISTS col_data")
    dbWriteTable(con, 
                 value = cols, 
                 name = "col_data", 
                 append = TRUE, 
                 row.names = TRUE, 
                 field.types = tableFieldTypes)
    dbCommit(conn = con)
    dbSendQuery(con, "ALTER TABLE `col_data` CHANGE COLUMN `row_names` `id` VARCHAR(255) NOT NULL;")
    #dbSendQuery(con, "ALTER TABLE `col_data` CHANGE COLUMN `index` `index` BIGINT NULL DEFAULT NULL;")
    dbSendQuery(con, "ALTER TABLE `col_data` ADD COLUMN `label` INT NOT NULL DEFAULT 1;")
    dbSendQuery(con, "UPDATE `col_data` SET `label` = `index`;")
    dbSendQuery(con, "ALTER TABLE `col_data` ADD INDEX `label_idx` USING HASH (`label` ASC);")
    dbSendQuery(con, "ALTER TABLE `col_data` ADD INDEX `index_idx` (`index` ASC) ;")
    dbCommit(conn = con)
  },
  .getValueTable=function() {
    
    #get counts
    counts = .counts
    countsIndices = expand.grid(seq(dim(counts)[1]), seq(dim(counts)[2]))
    
    df = data.frame(row=countsIndices[,1], col=countsIndices[,2], val=as.vector(counts))
    df = df[df$val != 0,]
    
    #get rows
    h = taxonomyTable()
    indexCombs = expand.grid(seq(dim(h)[1]), seq(dim(h)[2]))
    
    leafIds = lapply(seq(dim(h)[1]), function(i) {
      calcNodeId(i, dim(h)[2])
    })
    
    df$NodeId = leafIds[df$row]
    df$SampleId = rownames(.sampleAnnotation)[df$col]
    
    df
  },
  .saveRowData=function(con) {
    cat("\n  Computing taxonomy leafs...\n")
    h = taxonomyTable()
    pb = txtProgressBar(style=3, width=25)
    leafIds = lapply(seq(dim(h)[1]), function(i) {
      setTxtProgressBar(pb, i/dim(h)[1])
      #removed e$ in case
      calcNodeId(i, dim(h)[2])
    })

    cat("\n  Outputting to database...")
    leafIndices = seq(dim(h)[1]) - 1

    leafNames = h[,length(colnames(h))]


    df = data.frame(label=unlist(leafNames), index=leafIndices, start=leafIndices, end=(leafIndices+1))
    df$partition = NA
    rownames(df) = unlist(leafIds)
    
    dbSendQuery(con, "DROP TABLE IF EXISTS row_data")
    
    tableFieldTypes = .getFieldTypes(df)
    
    dbWriteTable(con, value = df, name = "row_data", append = TRUE, 
                 row.names = TRUE, 
                 field.types = tableFieldTypes)     
    dbCommit(conn = con)
    dbSendQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `row_names` `id` VARCHAR(255) NOT NULL;")

    dbSendQuery(con, "ALTER TABLE `row_data` ADD INDEX `index_idx` USING BTREE (`index` ASC);")
    dbSendQuery(con, "ALTER TABLE `row_data` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    
    #dbSendQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `label` `label` VARCHAR(255) NOT NULL;")
    dbSendQuery(con, "ALTER TABLE `row_data` ADD INDEX `label_idx` USING HASH (`label`(10)) USING HASH;")
    dbCommit(conn = con)
  },

  .saveHierarchy=function(con) {
    h = taxonomyTable()
    indexCombs = expand.grid(seq(dim(h)[1]), seq(dim(h)[2]))

    cat("\n  Extracting taxonomy nodes...\n")
    pb = txtProgressBar(style=3, width=25)
    nodeIds = lapply(seq(dim(indexCombs)[1]), function(i) {
      setTxtProgressBar(pb, i/dim(indexCombs)[1])
      calcNodeId(indexCombs[i, 1], indexCombs[i, 2])
    })
    uniqueIds = unique(nodeIds)

    cat("\n  Generating taxonomy structure...\n")
    pb = txtProgressBar(style=3, width=25)
    names = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      #node(nodeId)$name()
      pair = .fromMetavizNodeId(nodeId)
      h[pair$leafIndex+1, pair$depth+1]
    })

    cat("\n  Computing node parents...\n")
    pb = txtProgressBar(style=3, width=25)
    parentIds = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$parentId()
    })

    cat("\n  Computing lineages...\n")
    pathsList = Ptr$new(list())
    pb = txtProgressBar(style=3, width=25)
    paths = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      parentId = parentIds[[i]]
      if (is.null(parentId) || is.null(pathsList$.[[parentId]])) {
        pathsList$.[[nodeId]] = nodeId
        return(nodeId)
      }
      path = paste(pathsList$.[[parentId]], nodeId, sep=",")
      pathsList$.[[nodeId]] = path
      return(path)
    })

    cat("\n  Computing lineages labels...\n")
    pathsLabels = Ptr$new(list())
    pb = txtProgressBar(style=3, width=25)
    pathsLabels = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))

      nodeId = uniqueIds[[i]]
      pair = .fromMetavizNodeId(nodeId)
      nodeName = h[pair$leafIndex+1, pair$depth+1]

      parentId = parentIds[[i]]

      if (is.null(parentId) || is.null(pathsList$.[[parentId]])) {
        pathsList$.[[nodeId]] = nodeName
        return(nodeName)
      }

      path = paste(pathsList$.[[parentId]], nodeName, sep=",")
      pathsList$.[[nodeId]] = path
      return(path)
    })

    cat("\n  Computing index of first leaf in node subtrees...\n")
    pb = txtProgressBar(style=3, width=25)
    starts = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$leafIndex()
    })

    cat("\n  Computing leaf counts in node subtrees...\n")
    pb = txtProgressBar(style=3, width=25)
    ends = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$nleaves() + starts[[i]]
    })

    cat("\n  Computing node depths...\n")
    pb = txtProgressBar(style=3, width=25)
    depths = lapply(seq(length(uniqueIds)), function(i) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$depth()
    })

    cat("\n  Computing node children counts...")
    nchildren = list()
    for (parentId in parentIds) {
      if (is.null(parentId)) { next }
      if (is.null(nchildren[[parentId]])) {
        nchildren[[parentId]] = 0
      }
      nchildren[[parentId]] = nchildren[[parentId]] + 1
    }
    childCount = lapply(uniqueIds, function(nodeId) {
      if (is.null(nchildren[[nodeId]])) { return(0) }
      nchildren[[nodeId]]
    })
    cat("Done\n")

    cat("\n  Computing nodes order...\n")
    lastOrder = list()
    pb = txtProgressBar(style=3, width=25)
    orders = list()
    for (i in seq(length(uniqueIds))) {
      setTxtProgressBar(pb, i/length(uniqueIds))
      parentId = parentIds[[i]]
      nodeId = uniqueIds[[i]]
      if (is.null(parentId)) {
        orders[[i]] = 0
        next
      }
      o = lastOrder[[parentId]]
      if (is.null(o)) {
        lastOrder[[parentId]] = 0
        orders[[i]] = 0
        next
      }

      lastOrder[[parentId]] = o + 1
      orders[[i]] = o+1
    }

    cat("\n  Outputting to database...")
    parentIds[[1]] = NA
    df = data.frame(depth=unlist(depths), label=unlist(names), parentId=unlist(parentIds), lineage=unlist(paths), lineageLabel=unlist(pathsLabels), nchildren=unlist(childCount), start=unlist(starts), end=unlist(ends), leafIndex=unlist(starts), nleaves=unlist(ends)-unlist(starts), order=unlist(orders))
    df$partition = NA

    rownames(df) = unlist(uniqueIds)

    dbSendQuery(con, "DROP TABLE IF EXISTS hierarchy")
    
    tableFieldTypes = .getFieldTypes(df)
    
    dbWriteTable(con, value = df, name = "hierarchy", append = TRUE, row.names = TRUE, 
                 field.types = tableFieldTypes) 
    dbCommit(conn = con)
    dbSendQuery(con, "ALTER TABLE `hierarchy` CHANGE COLUMN `row_names` `id` VARCHAR(255) NOT NULL;")
    # dbSendQuery(con, "UPDATE `hierarchy` SET `id` = `row_names`;")
    dbSendQuery(con, "ALTER TABLE `hierarchy` ENGINE = MyISAM;")
    #dbSendQuery(con, "ALTER TABLE `hierarchy` CHANGE COLUMN `partition` `partition` VARCHAR(255) NULL DEFAULT NULL;")
    #Will do this later
    dbSendQuery(con, "ALTER TABLE `hierarchy` ADD INDEX `name_idx` USING HASH (`label` ASC);")
    dbSendQuery(con, "ALTER TABLE `hierarchy` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    dbSendQuery(con, "ALTER TABLE `hierarchy` ADD FULLTEXT INDEX `lineage_idx` (`lineage` ASC);")
    
    dbCommit(conn = con)
    },

  .saveValues=function(con) {
    counts = .counts
    countsIndices = expand.grid(seq(dim(counts)[1]), seq(dim(counts)[2]))

    df = data.frame(row=countsIndices[,1], col=countsIndices[,2], val=as.vector(counts))
    df = df[df$val != 0,]

    #requires user to change name of table to `values` after running this function
    dbSendQuery(con, "DROP TABLE IF EXISTS meta_values")
    dbSendQuery(con, "DROP TABLE IF EXISTS `values`")
    
    dbCommit(conn = con)
    
    tableFieldTypes = .getFieldTypes(df)
    
    dbWriteTable(con, value = df, name = "meta_values", append = TRUE, row.names = TRUE, 
                 field.types = tableFieldTypes ) 
    dbCommit(conn = con)
    
    dbSendQuery(con, "ALTER TABLE `meta_values` CHANGE COLUMN `row_names` `id` BIGINT NOT NULL  ;")
    dbSendQuery(con, "ALTER TABLE `meta_values` ADD INDEX `rowcol_idx` USING BTREE (`row` ASC, `col` ASC) ;")
    
    dbCommit(conn = con)

  },

  .saveLevels=function(con) {
    df = data.frame(depth=seq(length(.levels)) - 1, label=.levels)
    
    dbSendQuery(con, "DROP TABLE IF EXISTS levels")
    
    tableFieldTypes = .getFieldTypes(df)
    
    dbWriteTable(con, value = df, name = "levels", append = TRUE , row.names = TRUE, 
                 field.types = tableFieldTypes) 
    dbCommit(conn = con)
    
    #Do these need to be run?
    dbSendQuery(con, "ALTER TABLE `levels` ENGINE = MEMORY ;")
    dbSendQuery(con, "ALTER TABLE `levels` DROP COLUMN `row_names` ;")
  },
  
  toNEO4jDump = function(filePath = "~/", fileName="neo4j_dump.cypher") {
    "Save an MRexperiment object to a Neo4j Graph database. set the location of the output file with fileName and filePath."
    
    file = paste0(filePath, '/' ,fileName)
    
    write("begin", file=file, append = TRUE)
    .saveSampleDataNEO4J(graph=NULL, file=file)
    write(";", file=file, append = TRUE)
    write("commit", file=file, append = TRUE)
    write("begin", file=file, append = TRUE)
    .saveHierarchyNEO4J(graph=NULL, file=file)
    write(";", file=file, append = TRUE)
    write("commit", file=file, append = TRUE)
    write("begin", file=file, append = TRUE)
    .saveMatrixNEO4J(graph=NULL, file=file)
    write(";", file=file, append = TRUE)
    write("commit", file=file, append = TRUE)
    write("begin", file=file, append = TRUE)
    .neo4jUpdateProperties(graph=NULL, file=file)
    write(";", file=file, append = TRUE)
    write("commit", file=file, append = TRUE)
  },
  
  toNEO4JDb=function(graph, colLabel=NULL) {
    "Save an MRexperiment object to a Neo4j Graph database."

    cat("Saving sample data...")
    .saveSampleDataNEO4J(graph)
    cat("Done\n")
    
    cat("Saving hierarchy...")
    .saveHierarchyNEO4J(graph)
    cat("Done\n")
    
    cat("Saving Data Matrix...")
    .saveMatrixNEO4J(graph)
    cat("Done\n")
    
    cat("Saving hierarchy...")
    .neo4jUpdateProperties(graph)
    cat("Done\n")
    
  },
  
  .saveSampleDataNEO4J=function(graph, file=NULL) {
    sampleAnnotationToNeo4j = .sampleAnnotation
    sampleAnnotationToNeo4j['id'] = rownames(.sampleAnnotation)
    keys = colnames(sampleAnnotationToNeo4j)
    for (j in 1:nrow(sampleAnnotationToNeo4j)){
      row <- sampleAnnotationToNeo4j[j,]
      query = "CREATE (:Sample { "
      for (i in 1:(length(keys)-1)){
        if  (typeof(keys[i]) == "numeric")
          query = paste(query, keys[i], " : ", gsub("'", "", row[, keys[i]]), ", ", sep="")
        else
          query = paste(query, keys[i], " : '", gsub("'", "",row[, keys[i]]), "', ",sep="")
      }
      i = length(keys)
      if  (typeof(keys[i]) == "numeric")
        query = paste(query, keys[i], " : ", gsub("'", "",row[, keys[i]]), "})", sep="")
      else
        query = paste(query, keys[i], " : '", gsub("'", "",row[, keys[i]]), "'})", sep="")
      
      if(!is.null(graph)) {
        cypher(graph,query) 
      }
      else {
        write(query, file=file, append = TRUE)
      }
    }
  },
  
  .saveHierarchyNEO4J=function(graph, file=NULL) {
    h = taxonomyTable()
    indexCombs = expand.grid(seq(dim(h)[1]), seq(dim(h)[2]))
    
    # cat("\n  Extracting taxonomy nodes...\n")
    # pb = txtProgressBar(style=3, width=25)
    nodeIds = lapply(seq(dim(indexCombs)[1]), function(i) {
      # setTxtProgressBar(pb, i/dim(indexCombs)[1])
      calcNodeId(indexCombs[i, 1], indexCombs[i, 2])
    })
    uniqueIds = unique(nodeIds)
    
    # cat("\n  Generating taxonomy structure...\n")
    # pb = txtProgressBar(style=3, width=25)
    names = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      #node(nodeId)$name()
      pair = .fromMetavizNodeId(nodeId)
      h[pair$leafIndex+1, pair$depth+1]
    })
    
    # cat("\n  Computing node parents...\n")
    # pb = txtProgressBar(style=3, width=25)
    parentIds = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$parentId()
    })
    
    # cat("\n  Computing lineages...\n")
    pathsList = Ptr$new(list())
    # pb = txtProgressBar(style=3, width=25)
    paths = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      parentId = parentIds[[i]]
      if (is.null(parentId) || is.null(pathsList$.[[parentId]])) {
        pathsList$.[[nodeId]] = nodeId
        return(nodeId)
      }
      path = paste(pathsList$.[[parentId]], nodeId, sep=",")
      pathsList$.[[nodeId]] = path
      return(path)
    })
    
    # cat("\n  Computing lineages labels...\n")
    pathsLabels = Ptr$new(list())
    # pb = txtProgressBar(style=3, width=25)
    pathsLabels = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      
      nodeId = uniqueIds[[i]]
      pair = .fromMetavizNodeId(nodeId)
      nodeName = h[pair$leafIndex+1, pair$depth+1]
      
      parentId = parentIds[[i]]
      
      if (is.null(parentId) || is.null(pathsList$.[[parentId]])) {
        pathsList$.[[nodeId]] = nodeName
        return(nodeName)
      }
      
      path = paste(pathsList$.[[parentId]], nodeName, sep=",")
      pathsList$.[[nodeId]] = path
      return(path)
    })
    
    # cat("\n  Computing index of first leaf in node subtrees...\n")
    # pb = txtProgressBar(style=3, width=25)
    starts = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$leafIndex()
    })
    
    # cat("\n  Computing leaf counts in node subtrees...\n")
    # pb = txtProgressBar(style=3, width=25)
    ends = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$nleaves() + starts[[i]]
    })
    
    # cat("\n  Computing node depths...\n")
    # pb = txtProgressBar(style=3, width=25)
    depths = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$depth()
    })
    
    taxonomies = list()
    taxonomies = lapply(seq(length(uniqueIds)), function(i) {
      nodeId = uniqueIds[[i]]
      tempDepth = node(nodeId)$depth()
      .levels[tempDepth+1]
    })
    
    
    # cat("\n  Computing node children counts...")
    nchildren = list()
    for (parentId in parentIds) {
      if (is.null(parentId)) { next }
      if (is.null(nchildren[[parentId]])) {
        nchildren[[parentId]] = 0
      }
      nchildren[[parentId]] = nchildren[[parentId]] + 1
    }
    childCount = lapply(uniqueIds, function(nodeId) {
      if (is.null(nchildren[[nodeId]])) { return(0) }
      nchildren[[nodeId]]
    })
    # cat("Done\n")
    
    # cat("\n  Computing nodes order...\n")
    lastOrder = list()
    # pb = txtProgressBar(style=3, width=25)
    orders = list()
    for (i in seq(length(uniqueIds))) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      parentId = parentIds[[i]]
      nodeId = uniqueIds[[i]]
      if (is.null(parentId)) {
        orders[[i]] = 0
        next
      }
      o = lastOrder[[parentId]]
      if (is.null(o)) {
        lastOrder[[parentId]] = 0
        orders[[i]] = 0
        next
      }
      
      lastOrder[[parentId]] = o + 1
      orders[[i]] = o+1
    }
    
    # cat("\n  Outputting to database...")
    parentIds[[1]] = NA
    dfToNeo4j = data.frame(depth=unlist(depths), label=unlist(names), 
                           parentId=unlist(parentIds), lineage=unlist(paths), 
                           lineageLabel=unlist(pathsLabels), nchildren=unlist(childCount), 
                           start=unlist(starts), end=unlist(ends), leafIndex=unlist(starts), 
                           nleaves=unlist(ends)-unlist(starts), order=unlist(orders), taxonomy=unlist(taxonomies))
    dfToNeo4j$partition = NA
    rownames(dfToNeo4j) = unlist(uniqueIds)
    dfToNeo4j['id'] = unlist(uniqueIds)
    
    
    keys = colnames(dfToNeo4j)
    
    cypherCount = 0
    for (j in 1:nrow(dfToNeo4j)){
      row <- dfToNeo4j[j,]
      query = "CREATE (:Feature { "
      for (i in 1:(length(keys)-1)){
        if  (typeof(keys[i]) == "numeric")
          query = paste(query, keys[i], " : ", gsub("'", "",row[, keys[i]]), ", ",sep="")
        else
          query = paste(query, keys[i], " : '", gsub("'", "",row[, keys[i]]), "', ", sep="")
      }
      i = length(keys)
      query = paste(query, keys[i], " : '", gsub("'", "",row[, keys[i]]), "'})", sep="")
      if(!is.null(graph)) {
        cypher(graph,query) 
      }
      else {
        write(query, file=file, append = TRUE)
        
        cypherCount = cypherCount + 1
        
        # write commits if data file is too long
        if(cypherCount == 250) {
          write(";", file=file, append = TRUE)
          write("commit", file=file, append = TRUE)
          write("begin", file=file, append = TRUE)
          cypherCount = 0
        }
      }
    }
    
    
    cypherCount = 0
    for (j in 1:nrow(dfToNeo4j)){
      row <- dfToNeo4j[j,]
      query = paste("MATCH (fParent:Feature {id :'", row$parentId, "'}) MATCH (f:Feature {id:'", row$id, "'}) CREATE (fParent)-[:PARENT_OF]->(f)", sep="")
      
      if(!is.null(graph)) {
        cypher(graph,query) 
      }
      else {
        write(query, file=file, append = TRUE)
        
        cypherCount = cypherCount + 1
        
        # write commits if data file is too long
        if(cypherCount == 250) {
          write(";", file=file, append = TRUE)
          write("commit", file=file, append = TRUE)
          write("begin", file=file, append = TRUE)
          cypherCount = 0
        }
      }
    }
    
    query = paste0("MATCH (fNode:Feature)-[:PARENT_OF*]->(fLeaf:Feature {depth:'", length(.levels) - 1 , "'}) CREATE (fNode)-[:LEAF_OF]->(fLeaf)")
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
    }
    
    
    query = paste0("MATCH (fLeaf:Feature {depth:'", length(.levels) - 1,"'}) CREATE (fLeaf)-[:LEAF_OF]->(fLeaf)")
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
    }
  },
  
  .saveMatrixNEO4J = function(graph, file=NULL) {
    valuesToNeo4j = .getValueTable()
    
    cypherCount = 0

    for (j in 1:nrow(valuesToNeo4j)){
      row <- valuesToNeo4j[j,]
      # query = paste(query, paste("MATCH (f:Feature {id :'", row$NodeId, "'}) MATCH (s:Sample {id:'", row$SampleId, "'}) CREATE (s)-[:VALUE {val: ", row$val, "}]->(f)", sep=""), sep= " ")
      query = paste("MATCH (f:Feature {id :'", row$NodeId, "'}) MATCH (s:Sample {id:'", row$SampleId, "'}) CREATE (s)-[:VALUE {val: ", row$val, "}]->(f)", sep="")

      if(!is.null(graph)) {
        cypher(graph,query) 
      }
      else {
        write(query, file=file, append = TRUE)
        
        cypherCount = cypherCount + 1
        
        # write commits if data file is too long
        if(cypherCount == 250) {
          write(";", file=file, append = TRUE)
          write("commit", file=file, append = TRUE)
          write("begin", file=file, append = TRUE)
          cypherCount = 0
        }
      }
    }
  },
  
  .neo4jUpdateProperties = function(graph, file=NULL) {
    query = "MATCH (f:Feature) SET f.depth = toInt(f.depth) SET f.start = toInt(f.start) SET f.end = toInt(f.end) SET f.leafIndex = toInt(f.leafIndex) SET f.nchildren = toInt(f.nchildren) SET f.nleaves = toInt(f.nleaves) SET f.order = toInt(f.order)"
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query = "CREATE INDEX ON :Feature (depth)"
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query = "CREATE INDEX ON :Feature (start)"
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query = "CREATE INDEX ON :Feature (end)"
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query = "CREATE INDEX ON :Feature (id)"
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query = "CREATE INDEX ON :Sample (id)"
    if(!is.null(graph)) {
      cypher(graph,query) 
    }
    else {
      write(query, file=file, append = TRUE)
    }
  },
  
  toNEO4JDbHTTP =function(batch_url, neo4juser, neo4jpass, datasource, colLabel=NULL) {
    "Save an MRexperiment object to a Neo4j Graph database."
    
    cat("Saving sample data...")
    .saveSampleDataNEO4JHTTP(batch_url, neo4juser, neo4jpass)
    cat("Done\n")
    
    cat("Saving datasource data...")
    .saveDataSourceNEO4JHTTP(batch_url, neo4juser, neo4jpass, datasource)
    cat("Done\n")
    
    cat("Saving hierarchy...")
    .saveHierarchyNEO4JHTTP(batch_url, neo4juser, neo4jpass, datasource)
    cat("Done\n")
    
    cat("Saving Data Matrix...")
    .saveMatrixNEO4JHTTP(batch_url, neo4juser, neo4jpass, datasource)
    cat("Done\n")
    
    cat("Saving properties...")
    .neo4jUpdatePropertiesHTTP(batch_url, neo4juser, neo4jpass)
    cat("Done\n")
    
  },
  
  .buildBatchJSON = function(query_in, param_list, id_in=0, id_last=TRUE, full_query=TRUE,json_query_in=NULL, params_complete=FALSE){
    json_start <- "["
    method <- "{\"method\" : \"POST\","
    to <- "\"to\" : \"cypher\","
    body_start <- "\"body\" : {"
    params_start <- "\"params\" : {"
    params_end <- "}"
    body_end <- "},"
    #id <- "\"id\" : 0}"
    if(id_last){
     id <- paste("\"id\": ", as.character(id_in), "}", sep="")
    }
    else{
      id <- paste("\"id\": ", as.character(id_in), "},", sep="")
    }
    json_end <- "]"
    json_query <- ""
    if(is.null(json_query_in)){
      json_query <- ""
    }
    else{
      json_query <- json_query_in
    }
    params = ""

    if(is.null(param_list)){
      query = paste("\"query\" : \"", query_in, "\"", sep = "")
      json_query <- paste0(json_query, method, to, body_start, query, body_end, id)
      query_final <- paste0(json_start, json_query, json_end)
      print(query_final)
      .self$.json_query <- query_final
      return(query_final)
    }
    
    query = paste("\"query\" : \"", query_in, "\",", sep = "")
    for(k in seq(1, length(param_list))){
      if(k == length(param_list)){
        if(params_complete){
          params <- paste0(params, "\"", names(param_list[k]), "\" : ", unlist(unname(param_list[k])), "")
          }
        else{
          params <- paste0(params, "\"", names(param_list[k]), "\" : {", unlist(unname(param_list[k])), "}")
        }
      }
      else{
        if(params_complete){
          params <- paste0(params, "\"", names(param_list[k]), "\" : ", unlist(unname(param_list[k])), ",")
        }
        else{
          params <- paste0(params, "\"", names(param_list[k]), "\" : {\"", unlist(unname(param_list[k])), "\"}", ",")
        }
      }
    }
    json_query <- paste0(json_query, method, to, body_start, query, params_start, params, params_end, body_end, id)
    print(json_query)
    if(!full_query){
      return(json_query)
    }
    query_final <- paste0(json_start, json_query, json_end)
    print(query_final)
    .self$.json_query <- query_final
    return(query_final)
  },
  
  .saveSampleDataNEO4JHTTP =function(batch_url, neo4juser, neo4jpass, file=NULL) {
    
    json_start <- "["
    query <- "CREATE (:Sample {props})"

    json_end <- "]"
    json_query <- ""
    
    sampleAnnotationToNeo4j = .sampleAnnotation
    sampleAnnotationToNeo4j['id'] = rownames(.sampleAnnotation)
    keys = colnames(sampleAnnotationToNeo4j)
    id_counter <- 0
    for (j in 1:nrow(sampleAnnotationToNeo4j)){
      row <- sampleAnnotationToNeo4j[j,]
      props <- ""
      for (i in 1:(length(keys)-1)){
        if  (typeof(keys[i]) == "numeric")
          props <- paste(props, "\"", keys[i], "\"", " : ", "\"", gsub("'", "", row[, keys[i]]), "\"",", ", sep="")
        else
          props <- paste(props, "\"", keys[i], "\"", " : \"", gsub("'", "",row[, keys[i]]), "\", ",sep="")
      }
      i = length(keys)
      if  (typeof(keys[i]) == "numeric")
        props = paste(props, "\"", keys[i], "\"", " : ", "\"", gsub("'", "",row[, keys[i]]), "\"", sep="")
      else
        props = paste(props, "\"", keys[i], "\"", " : \"", gsub("'", "",row[, keys[i]]), "\"", sep="")
      
      params <- list("props"=props)
      if (id_counter < nrow(sampleAnnotationToNeo4j)-1){
         json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = FALSE, full_query = FALSE, json_query_in = json_query)
      }
      else {
        json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query)
      }
      id_counter <- id_counter + 1
    }
    query_final <- paste0(json_start, json_query, json_end)
    .self$.json_query <- query_final
    
    if(!is.null(batch_url)) {
      POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
    }
    else {
      write(query_final, file=file, append = TRUE)
    }
    
  },
  
  .saveDataSourceNEO4JHTTP =function(batch_url, neo4juser, neo4jpass, datasource, file=NULL) {
    query <- "CREATE (:Datasource {label: {label_param}})"
    params <- list(label_param=paste("\"", as.character(datasource), "\"", sep=""))
    
    query_final <- .buildBatchJSON(query_in = query, param_list = params, params_complete = TRUE) 
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
      write(query, file=file, append = TRUE)
    }
  },
  
  .saveHierarchyNEO4JHTTP =function(batch_url, neo4juser, neo4jpass, datasource, file=NULL) {
    h = taxonomyTable()
    indexCombs = expand.grid(seq(dim(h)[1]), seq(dim(h)[2]))
    
    # cat("\n  Extracting taxonomy nodes...\n")
    # pb = txtProgressBar(style=3, width=25)
    nodeIds = lapply(seq(dim(indexCombs)[1]), function(i) {
      # setTxtProgressBar(pb, i/dim(indexCombs)[1])
      calcNodeId(indexCombs[i, 1], indexCombs[i, 2])
    })
    uniqueIds = unique(nodeIds)
    
    # cat("\n  Generating taxonomy structure...\n")
    # pb = txtProgressBar(style=3, width=25)
    names = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      #node(nodeId)$name()
      pair = .fromMetavizNodeId(nodeId)
      h[pair$leafIndex+1, pair$depth+1]
    })
    
    # cat("\n  Computing node parents...\n")
    # pb = txtProgressBar(style=3, width=25)
    parentIds = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$parentId()
    })
    
    # cat("\n  Computing lineages...\n")
    pathsList = Ptr$new(list())
    # pb = txtProgressBar(style=3, width=25)
    paths = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      parentId = parentIds[[i]]
      if (is.null(parentId) || is.null(pathsList$.[[parentId]])) {
        pathsList$.[[nodeId]] = nodeId
        return(nodeId)
      }
      path = paste(pathsList$.[[parentId]], nodeId, sep=",")
      pathsList$.[[nodeId]] = path
      return(path)
    })
    
    # cat("\n  Computing lineages labels...\n")
    pathsLabels = Ptr$new(list())
    # pb = txtProgressBar(style=3, width=25)
    pathsLabels = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      
      nodeId = uniqueIds[[i]]
      pair = .fromMetavizNodeId(nodeId)
      nodeName = h[pair$leafIndex+1, pair$depth+1]
      
      parentId = parentIds[[i]]
      
      if (is.null(parentId) || is.null(pathsList$.[[parentId]])) {
        pathsList$.[[nodeId]] = nodeName
        return(nodeName)
      }
      
      path = paste(pathsList$.[[parentId]], nodeName, sep=",")
      pathsList$.[[nodeId]] = path
      return(path)
    })
    
    # cat("\n  Computing index of first leaf in node subtrees...\n")
    # pb = txtProgressBar(style=3, width=25)
    starts = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$leafIndex()
    })
    
    # cat("\n  Computing leaf counts in node subtrees...\n")
    # pb = txtProgressBar(style=3, width=25)
    ends = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$nleaves() + starts[[i]]
    })
    
    # cat("\n  Computing node depths...\n")
    # pb = txtProgressBar(style=3, width=25)
    depths = lapply(seq(length(uniqueIds)), function(i) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      nodeId = uniqueIds[[i]]
      node(nodeId)$depth()
    })
    
    taxonomies = list()
    taxonomies = lapply(seq(length(uniqueIds)), function(i) {
      nodeId = uniqueIds[[i]]
      tempDepth = node(nodeId)$depth()
      .levels[tempDepth+1]
    })
    
    
    # cat("\n  Computing node children counts...")
    nchildren = list()
    for (parentId in parentIds) {
      if (is.null(parentId)) { next }
      if (is.null(nchildren[[parentId]])) {
        nchildren[[parentId]] = 0
      }
      nchildren[[parentId]] = nchildren[[parentId]] + 1
    }
    childCount = lapply(uniqueIds, function(nodeId) {
      if (is.null(nchildren[[nodeId]])) { return(0) }
      nchildren[[nodeId]]
    })
    # cat("Done\n")
    
    # cat("\n  Computing nodes order...\n")
    lastOrder = list()
    # pb = txtProgressBar(style=3, width=25)
    orders = list()
    for (i in seq(length(uniqueIds))) {
      # setTxtProgressBar(pb, i/length(uniqueIds))
      parentId = parentIds[[i]]
      nodeId = uniqueIds[[i]]
      if (is.null(parentId)) {
        orders[[i]] = 0
        next
      }
      o = lastOrder[[parentId]]
      if (is.null(o)) {
        lastOrder[[parentId]] = 0
        orders[[i]] = 0
        next
      }
      
      lastOrder[[parentId]] = o + 1
      orders[[i]] = o+1
    }
    
    # cat("\n  Outputting to database...")
    parentIds[[1]] = NA
    dfToNeo4j = data.frame(depth=unlist(depths), label=unlist(names), 
                           parentId=unlist(parentIds), lineage=unlist(paths), 
                           lineageLabel=unlist(pathsLabels), nchildren=unlist(childCount), 
                           start=unlist(starts), end=unlist(ends), leafIndex=unlist(starts), 
                           nleaves=unlist(ends)-unlist(starts), order=unlist(orders), taxonomy=unlist(taxonomies))
    dfToNeo4j$partition = NA
    rownames(dfToNeo4j) = unlist(uniqueIds)
    dfToNeo4j['id'] = unlist(uniqueIds)
    
    json_start <- "["
    json_end <- "]"
    datasource_param_key <- "datasource"
    datasource_param_value <- as.character(datasource)
    
    keys = colnames(dfToNeo4j)
    query <- "CREATE (:Feature {props})"
    json_query <- ""
    cypherCount = 0
    id_counter <- 0
    
    for (j in 1:nrow(dfToNeo4j)){
      row <- dfToNeo4j[j,]
      props <- "" 
      for (i in 1:(length(keys))){
        if  (typeof(keys[i]) == "numeric")
          props <- paste(props, "\"", keys[i], "\"", " : ", "\"", gsub("'", "",row[, keys[i]]), "\"", ", ",sep="")
        else
          props <- paste(props, "\"", keys[i], "\"", " : \"", gsub("'", "",row[, keys[i]]), "\", ", sep="")
      }
      i = length(keys)
      
      props <- paste(props, "\"", datasource_param_key, "\"", " : \"", datasource_param_value, "\"", sep="")
      params <- list("props"=props)
      
      writeBatchOut <- j%%50 == 0
      if (j <= nrow(dfToNeo4j)-1){
        if (writeBatchOut){
            json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query)
        }
        else {
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = FALSE, full_query = FALSE, json_query_in = json_query)
        }
      }
      else {
        json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query)
      }
      id_counter <- id_counter + 1
      
      if (writeBatchOut){
        query_final <- paste0(json_start, json_query, json_end)
        .self$.json_query <- query_final
        
        if(!is.null(batch_url)) {
          r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
          stop_for_status(r)
          
        }
        id_counter <- 0
        json_query <- ""
      }
    }
    query_final <- paste0(json_start, json_query, json_end)
    .self$.json_query <- query_final
    
    if(!is.null(batch_url)) {
      POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
    }
    else {
      write(query, file=file, append = TRUE)
      
      cypherCount = cypherCount + 1
      
      # write commits if data file is too long
      if(cypherCount == 250) {
        write(";", file=file, append = TRUE)
        write("commit", file=file, append = TRUE)
        write("begin", file=file, append = TRUE)
        cypherCount = 0
      }
    }

    query <- "MATCH (ds:Datasource {label:{datasource_param}}) MATCH (fNode:Feature {id: {root_id}, datasource: {datasource_param}}) CREATE (ds)-[:DATASOURCE_OF]->(fNode)"
    params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), root_id=paste("\"", "0-0", "\"", sep=""))
    query_final <- .buildBatchJSON(query_in = query, param_list = params, params_complete = TRUE) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
      write(query, file=file, append = TRUE)
    }
    
    json_start <- "["
    query <- "MATCH (fParent:Feature {id : {parent_id}, datasource: {datasource_param}}) MATCH (f:Feature {id : {child_id}, datasource: {datasource_param}}) CREATE (fParent)-[:PARENT_OF]->(f)"

    json_end <- "]"
    json_query <- ""

    cypherCount = 0
    id_counter <- 0
    
    for (j in 1:nrow(dfToNeo4j)){
      row <- dfToNeo4j[j,]
      parentid <- paste("\"", row$parentId, "\"", sep="")
      childid <- paste("\"", row$id, "\"", sep="")
      writeBatchOut <- j%%50 == 0
      if (j <= nrow(dfToNeo4j)-1){
        if (writeBatchOut){
          params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), parent_id=parentid, child_id=childid)
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
        }
        else {
          params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), parent_id=parentid, child_id=childid)
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = FALSE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
        }
      }
      else {
        params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), parent_id=parentid, child_id=childid)
        json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
      }
      id_counter <- id_counter + 1
    
      if (writeBatchOut){
        query_final <- paste0(json_start, json_query, json_end)
        .self$.json_query <- query_final
        
        if(!is.null(batch_url)) {
          r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
          stop_for_status(r)
        }
        id_counter <- 0
        json_query <- ""        
      }
    }
    query_final <- paste0(json_start, json_query, json_end)
    .self$.json_query <- query_final
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
      write(query, file=file, append = TRUE)
        
      cypherCount = cypherCount + 1
      
      # write commits if data file is too long
      if(cypherCount == 250) {
        write(";", file=file, append = TRUE)
        write("commit", file=file, append = TRUE)
        write("begin", file=file, append = TRUE)
        cypherCount = 0
      }
    }

    query <- "MATCH (fNode:Feature {datasource: {datasource_param}})-[:PARENT_OF*]->(fLeaf:Feature {depth: {depth_param}, datasource: {datasource_param} }) CREATE (fNode)-[:LEAF_OF]->(fLeaf)"
    params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), depth_param=paste("\"", as.character(length(.levels) - 1), "\"", sep=""))
    query_final <- .buildBatchJSON(query_in = query, param_list = params, params_complete = TRUE) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
      write(query, file=file, append = TRUE)
    }

    query <- "MATCH (fLeaf:Feature {depth: {depth_param}, datasource: {datasource_param}}) CREATE (fLeaf)-[:LEAF_OF]->(fLeaf)"
    params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), depth_param=paste("\"", as.character(length(.levels) - 1), "\"", sep=""))
    query_final <- .buildBatchJSON(query_in = query, param_list = params, params_complete = TRUE) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
      write(query, file=file, append = TRUE)
    }

  },
  
  .saveMatrixNEO4JHTTP = function(batch_url, neo4juser, neo4jpass,  datasource, file=NULL) {
    valuesToNeo4j = .getValueTable()
    
    json_start <- "["
    query <- "MATCH (f:Feature {id : {node_id}, datasource: {datasource_param}}) MATCH (s:Sample {id: {sample_id}}) CREATE (s)-[:COUNT {val: {count_param}}]->(f)"

    json_end <- "]"
    json_query <- ""
    
    cypherCount = 0
    id_counter <- 0
    
    for (j in 1:nrow(valuesToNeo4j)){
      row <- valuesToNeo4j[j,]
      nodeid <- paste("\"", row$NodeId, "\"", sep="")
      sampleid <- paste("\"", row$SampleId, "\"", sep="")
      count <- paste("\"", as.character(row$val), "\"", sep="")
      writeBatchOut <- j%%1000 == 0
      if (j <= nrow(valuesToNeo4j)-1){
        if (writeBatchOut){
          params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), node_id=nodeid, sample_id=sampleid, count_param=count)
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
        }
        else {
          params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), node_id=nodeid, sample_id=sampleid, count_param=count)
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = FALSE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
        }
      }
      else {
        params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), node_id=nodeid, sample_id=sampleid, count_param=count)
        json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
      }
      id_counter <- id_counter + 1
      
      if (writeBatchOut){
        query_final <- paste0(json_start, json_query, json_end)
        .self$.json_query <- query_final
        
        if(!is.null(batch_url)) {
          r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
          stop_for_status(r)
        }
        id_counter <- 0
        json_query <- ""
      }

    }
    
    query_final <- paste0(json_start, json_query, json_end)
    .self$.json_query <- query_final
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
        write(query, file=file, append = TRUE)
        
        cypherCount = cypherCount + 1
        
        # write commits if data file is too long
        if(cypherCount == 250) {
          write(";", file=file, append = TRUE)
          write("commit", file=file, append = TRUE)
          write("begin", file=file, append = TRUE)
          cypherCount = 0
        }
    }
    
  },
  
  .neo4jUpdatePropertiesHTTP = function(batch_url, neo4juser, neo4jpass, file=NULL) {
    
    query <- "MATCH (f:Feature) SET f.depth = toInt(f.depth) SET f.start = toInt(f.start) SET f.end = toInt(f.end) SET f.leafIndex = toInt(f.leafIndex) SET f.nchildren = toInt(f.nchildren) SET f.nleaves = toInt(f.nleaves) SET f.order = toInt(f.order)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }

    query <- "MATCH ()-[c:COUNT]->() SET c.val = toFloat(c.val)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }

    query <- "CREATE INDEX ON :Feature (depth)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query <- "CREATE INDEX ON :Feature (start)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }

    query <- "CREATE INDEX ON :Feature (end)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query <- "CREATE INDEX ON :Feature (id)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
      write(";", file=file, append = TRUE)
      write("commit", file=file, append = TRUE)
      write("begin", file=file, append = TRUE)
    }
    
    query <- "CREATE INDEX ON :Sample (id)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    } else {
      write(query, file=file, append = TRUE)
    }
  }
)
