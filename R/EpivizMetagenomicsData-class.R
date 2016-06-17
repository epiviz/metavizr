EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains="EpivizData",
  fields=list(
    .taxonomy="MetavizTree",
    .levels="character",
    .maxDepth="numeric",
    .aggregateAtDepth="numeric",
    .lastRootId="ANY",

    .counts="ANY",
    .sampleAnnotation="ANY",

    .minValue="ANY",
    .maxValue="ANY",
    .aggregateFun="ANY",
    .valuesAnnotationFuns="ANY",

    .lastRequestRanges="list",
    .lastLeafInfos="list",
    .lastValues="list",
    .maxHistory="numeric"
  ),
  methods=list(
    initialize=function(object,control=metavizControl(), ...) {

      # Initialize parameters used here
      aggregateAtDepth = control$aggregateAtDepth
      maxDepth         = control$maxDepth
      maxHistory       = control$maxHistory
      maxValue         =  control$maxValue
      minValue         =  control$minValue
      aggregateFun     =  control$aggregateFun
      valuesAnnotationFuns=control$valuesAnnotationFuns

      if(is.character(aggregateAtDepth)){ aggregateAtDepth = assignValues(object,aggregateAtDepth) }
      if(is.character(maxDepth)){ maxDepth = assignValues(object,maxDepth) }
      log = control$log
      norm= control$norm

      # TODO: Some type checking
      .taxonomy <<- buildMetavizTree(object)
      .levels <<- .taxonomy$levels()
      .maxDepth <<- maxDepth
      .aggregateAtDepth <<- aggregateAtDepth
      .lastRootId <<- .taxonomy$root()$id()

      t = .taxonomy$taxonomyTable()
      .counts <<- MRcounts(object[rownames(t),],norm=norm,log=log)


      # TODO: Make this consistent with the aggregateFun
      if (is.null(minValue)) {
        minValue = log2(min(.counts) + 1)
      }
      .minValue <<- minValue

      if (is.null(maxValue)) {
        maxValue = log2(max(.counts) + 1)
      }
      .maxValue <<- maxValue
      .aggregateFun <<- aggregateFun
      .valuesAnnotationFuns <<- valuesAnnotationFuns

      .sampleAnnotation <<- pData(object)

      .maxHistory <<- maxHistory
      .lastRequestRanges <<- list()
      .lastLeafInfos <<- list()
      .lastValues <<- list()

      if (.aggregateAtDepth >= 0) {
        nodesAtDepth = .taxonomy$nodesAtDepth(.aggregateAtDepth)
        selection = list()
        for (node in nodesAtDepth) {
          selection[[node$id()]] = SelectionType$NODE
        }
        .taxonomy$updateSelection(selection)
      }

      callSuper(object=object, ...)
    },

    .taxonomyLevels=function(exp) {
      .levels
    },

    .getSelectedLeaves=function(start, end) {
      ret = NULL
      if (length(.lastRequestRanges) > 0) {
        for (i in rev(seq_along(.lastRequestRanges))) {
          if (.lastRequestRanges[[i]]$start == start || .lastRequestRanges[[i]]$end == end) {
            ret = .lastLeafInfos[[i]]
            break
          }
        }
      }
      if (is.null(ret)) {
        requestRange = list(start=start, end=end)
        ret = Ptr$new(.taxonomy$selectedLeaves(start, end))

        if (.maxHistory > 0) {
          .lastLeafInfos <<- c(.lastLeafInfos, ret)
          .lastRequestRanges <<- c(.lastRequestRanges, list(requestRange))
          .lastValues <<- c(.lastValues, list(NULL))

          if (length(.lastRequestRanges) > .maxHistory) {
            .lastRequestRanges <<- .lastRequestRanges[2:(.maxHistory+1)]
            .lastLeafInfos <<- .lastLeafInfos[2:(.maxHistory+1)]
            .lastValues <<- .lastValues[2:(.maxHistory+1)] # Not yet a value
          }
        }
      }

      return(ret$.)
    },

    .getSelectedValues=function(measurement, start, end) {
      leafInfos = .getSelectedLeaves(start, end)
      index = -1
      values = NULL
      if (length(.lastRequestRanges) > 0) {
        for (i in rev(seq_along(.lastRequestRanges))) {
          if (.lastRequestRanges[[i]]$start == start || .lastRequestRanges[[i]]$end == end) {
            if (!is.null(.lastValues[[i]])) {
              values = .lastValues[[i]]
            } else {
              index = i
            }
            break
          }
        }
      }

      if (is.null(values)) {
        values = Ptr$new(list(values=unname(lapply(leafInfos, function(info) {
          lind = info$node$leafIndex()
          ind = (lind+1):(lind+info$node$nleaves())
          return(.aggregateFun(.counts[ind,, drop=FALSE]))
        }))))
        if (!is.null(.valuesAnnotationFuns)) {
          for (anno in names(.valuesAnnotationFuns)) {
            fun = .valuesAnnotationFuns[[anno]]
            values$.[[anno]] = unname(lapply(leafInfos, function(info) {
              lind = info$node$leafIndex()
              ind = (lind+1):(lind+info$node$nleaves())
              return(fun(.counts[ind,, drop=FALSE]))
            }))
          }
        }
      }
      if (index > 0) {
        .lastValues[[index]] <<- values
      }

      lapply(values$., function(vals) {
        lapply(vals, function(v) {
          if (length(dim(v)) > 0) { return(v[, measurement]) }
          return(v[[measurement]])
        })
      })
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
    if (is.null(dim(.counts))) {
      if (length(.counts) > 0) { return(1) }
      return(0)
    }

    nrow(.counts)
  },
  nmeasurements=function() {
    if (is.null(dim(.counts))) {
      return(length(.counts))
    }

    ncol(.counts)
  },

  nlevels=function() {
    length(.levels)
  },

  taxonomyTable=function() { .taxonomy$taxonomyTable() },
  calcNodeId=function(rowIndex, colIndex) { .taxonomy$calcNodeId(rowIndex, colIndex) },
  node=function(nodeId) { .taxonomy$node(nodeId) },
  parent=function(node) { .taxonomy$parent(node) },
  siblings=function(node) { .taxonomy$siblings(node) },

  changeAggregation=function(mgr, nodeId, aggregationType) {
    selection = list()
    selection[[nodeId]] = aggregationType
    .taxonomy$updateSelection(selection)
    mgr$.clearDatasourceGroupCache(.self, TRUE)
  },
  changeAggregationAtDepth=function(mgr, depth, aggregationType) {
    if (depth < 0) { return() }
    nodesAtDepth = .taxonomy$nodesAtDepth(depth)
    selection = list()
    for (node in nodesAtDepth) {
      selection[[node$id()]] = aggregationType
    }
    .taxonomy$updateSelection(selection)
    mgr$.clearDatasourceGroupCache(.self, TRUE)
  }
)

# Epiviz Websockets Protocol
EpivizMetagenomicsData$methods(
  getMeasurements=function() {
    out <- lapply(colnames(.counts), function(sample) {
      list(id=sample,
           name=sample,
           type="feature",
           datasourceId=id,
           datasourceGroup=id,
           defaultChartType="heatmap",
           annotation=as.list(.sampleAnnotation[sample,]),
           minValue=.minValue,
           maxValue=.maxValue,
           metadata=c(rev(.levels), "colLabel", "ancestors", "lineage"))
    })
    out
  },
  getHierarchy=function(nodeId) {
    root = NULL
    if (missing(nodeId) || is.null(nodeId)) { root = .taxonomy$root() }
    else {
      root = .taxonomy$parent(.taxonomy$node(nodeId))
      if (is.null(root)) { root = .taxonomy$root() }
    }
    .lastRootId <<- nodeId

    return(root$raw(maxDepth=.maxDepth))
  },
  propagateHierarchyChanges=function(selection, order) {
    if (missing(selection) && missing(order)) { return(getHierarchy(.lastRootId)) }

    if (!missing(selection)) {
      .taxonomy$updateSelection(selection)
    }

    if (!missing(order)) {
      .taxonomy$updateOrder(order)
    }

    .lastRequestRanges <<- list()
    .lastLeafInfos <<- list()
    .lastValues <<- list()

    getHierarchy(.lastRootId)
  },
  getRows=function(seqName, start, end, metadata) {
    leafInfos = .self$.getSelectedLeaves(start, end)
    leafAncestors = lapply(leafInfos, function(info) { .taxonomy$ancestors(info$node) })

    leafTaxonomies = list()
    for (i in seq_along(.levels)) {
      leafTaxonomies[[.levels[[i]]]] = list()
      for (j in seq_along(leafInfos)) {
        if (i > length(leafAncestors[[j]])) { next }
        leafTaxonomies[[.levels[[i]]]][[j]] = leafAncestors[[j]][[i]]
      }
    }

    ret = list(
      id = sapply(leafInfos, function(info) { info$realNodesBefore }),
      start=sapply(leafInfos, function(info) { info$start }),
      end=sapply(leafInfos, function(info) { info$start + info$node$nleaves() - 1 }),
      metadata = c(list(
        colLabel = sapply(leafInfos, function(info) { info$node$name() }),
        ancestors = sapply(leafAncestors, function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$name() }), collapse=",") }), # TODO: Use tree .ancestryByDepth
        lineage = sapply(leafAncestors, function(ancestors) { paste(lapply(rev(ancestors), function(node) { node$id() }), collapse=",") })
      ), sapply(rev(.levels), function(level) {
        r = list(lapply(leafTaxonomies[[level]], function(node) { if(is.null(node)) { return("<NA>") }; node$name() }))
        if (length(r[[1]]) == 0) {
          r[[1]] = lapply(seq_along(leafInfos), function(i) { "<NA>" })
        }
        return(r)
      }))
    )

    globalStartIndex = NULL
    if (length(leafInfos) > 0) { globalStartIndex = leafInfos[[1]]$realNodesBefore }
    if (length(leafInfos) == 1) {
      ret$id = list(ret$id)
      ret$start = list(ret$start)
      ret$end = list(ret$end)
      ret$metadata$colLabel = list(ret$metadata$colLabel)
      ret$metadata$ancestors = list(ret$metadata$ancestors)
      ret$metadata$lineage = list(ret$metadata$lineage)
    }

    return(list(
      globalStartIndex = globalStartIndex,
      values = ret
    ))
  },
  getValues=function(measurement, seqName, start, end) {
    leafInfos = .self$.getSelectedLeaves(start, end)
    globalStartIndex = NULL
    if (length(leafInfos) > 0) { globalStartIndex = leafInfos[[1]]$realNodesBefore }
    ret = list(
      globalStartIndex = globalStartIndex,
      values = .self$.getSelectedValues(measurement, start, end)
    )
    return(ret)
  }
)

# Database serialization
EpivizMetagenomicsData$methods(
  toMySQLDb=function(con, colLabel=NULL) {
    library(RMySQL)
    # TODO: Add logging
    # con = odbcDriverConnect(connectionStr)
    # con = odbcConnect('myodbc')
    #con <- dbConnect(MySQL())
    
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
    dbDisconnect(con)
    #odbcClose(con)
    
  },
  .saveColData=function(con, labelCol=NULL) {
    #odbcSetAutoCommit(con, autoCommit = FALSE)
    
    cols = .sampleAnnotation
    if (!is.null(labelCol)) {
      cols$label = cols[[labelCol]]
    }
    cols$index = seq(dim(cols)[1]) - 1

        print(cols)
    
    #sqlDrop(con, "col_data", errors=FALSE)
    #sqlSave(con, cols, tablename="col_data", addPK=TRUE)
    #sqlQuery(con, "ALTER TABLE `col_data` CHANGE COLUMN `rownames` `id` VARCHAR(255) NOT NULL;")
    #sqlQuery(con, "ALTER TABLE `col_data` CHANGE COLUMN `index` `index` BIGINT NULL DEFAULT NULL;")
    #sqlQuery(con, "ALTER TABLE `col_data` ADD INDEX `label_idx` USING HASH (`label` ASC);")
    #sqlQuery(con, "ALTER TABLE `col_data` ADD INDEX `index_idx` (`index` ASC) ;")
    #odbcSetAutoCommit(con, autoCommit = TRUE)
    
    dbSendQuery(con, "DROP TABLE IF EXISTS col_data")
    dbWriteTable(con, value = cols, name = "col_data", append = TRUE, row.names = TRUE )     
    dbCommit(conn = con)
    dbSendQuery(con, "ALTER TABLE `col_data` CHANGE COLUMN `row_names` `id` VARCHAR(255) NOT NULL;")
    dbSendQuery(con, "ALTER TABLE `col_data` CHANGE COLUMN `index` `index` BIGINT NULL DEFAULT NULL;")
    dbSendQuery(con, "ALTER TABLE `col_data` ADD COLUMN `label` INT NOT NULL DEFAULT 1;")
    dbSendQuery(con, "UPDATE `col_data` SET `label` = `index`;")
    dbSendQuery(con, "ALTER TABLE `col_data` ADD INDEX `label_idx` USING HASH (`label` ASC);")
    dbSendQuery(con, "ALTER TABLE `col_data` ADD INDEX `index_idx` (`index` ASC) ;")
    dbCommit(conn = con)
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

    print(head(h))
    leafNames = h[,length(colnames(h))]


    df = data.frame(label=unlist(leafNames), index=leafIndices, start=leafIndices, end=(leafIndices+1))
    df$partition = NA
    rownames(df) = unlist(leafIds)

    print(head(df))
    
    #odbcSetAutoCommit(con, autoCommit = FALSE)
    #sqlDrop(con, "row_data", errors=FALSE)
    #sqlSave(con, df, tablename="row_data", addPK=TRUE)
    # sqlQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `rownames` `id` VARCHAR(255) NOT NULL  ;")
    # sqlQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `index` `index` BIGINT NULL;")
    # sqlQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `partition` `partition` VARCHAR(255) NULL DEFAULT NULL;")
    # sqlQuery(con, "ALTER TABLE `row_data` ADD INDEX `index_idx` USING BTREE (`index` ASC);")
    # sqlQuery(con, "ALTER TABLE `row_data` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    # sqlQuery(con, "ALTER TABLE `row_data` ADD INDEX `label_idx` USING HASH (`label` ASC);")
    # odbcSetAutoCommit(con, autoCommit = TRUE)
    
    dbSendQuery(con, "DROP TABLE IF EXISTS row_data")
    dbWriteTable(con, value = df, name = "row_data", append = TRUE, row.names = TRUE)     
    dbCommit(conn = con)
    dbSendQuery(con, "ALTER TABLE `row_data` ADD COLUMN `id` VARCHAR(255) NOT NULL;")
    dbSendQuery(con, "UPDATE `row_data` SET `id` = `row_names`;")
    dbSendQuery(con, "UPDATE `row_data` SET `label` = `row_names`;")
    
    # dbSendQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `label` `id` VARCHAR(255) NOT NULL  ;")
    
    dbSendQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `index` `index` BIGINT NULL;")
    dbSendQuery(con, "ALTER TABLE `row_data` CHANGE COLUMN `partition` `partition` VARCHAR(255) NULL DEFAULT NULL;")
    dbSendQuery(con, "ALTER TABLE `row_data` ADD INDEX `index_idx` USING BTREE (`index` ASC);")
    dbSendQuery(con, "ALTER TABLE `row_data` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    
    #Will need to do this at some point
    #dbSendQuery(con, "ALTER TABLE `row_data` ADD INDEX `label_idx` USING HASH (`label` ASC);")
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

    print(head(df))
    rownames(df) = unlist(uniqueIds)

    #odbcSetAutoCommit(con, autoCommit = FALSE)
    #     sqlDrop(con, "hierarchy", errors=FALSE)
    #     sqlSave(con, df, tablename="hierarchy", addPK=TRUE)
    #     sqlQuery(con, "ALTER TABLE `hierarchy` ENGINE = MyISAM;")
    #     sqlQuery(con, "ALTER TABLE `hierarchy` CHANGE COLUMN `rownames` `id` VARCHAR(255) NOT NULL;")
    #     sqlQuery(con, "ALTER TABLE `hierarchy` CHANGE COLUMN `partition` `partition` VARCHAR(255) NULL DEFAULT NULL;")
    #     sqlQuery(con, "ALTER TABLE `hierarchy` ADD INDEX `name_idx` USING HASH (`label` ASC);")
    #     sqlQuery(con, "ALTER TABLE `hierarchy` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    #     sqlQuery(con, "ALTER TABLE `hierarchy` ADD FULLTEXT INDEX `lineage_idx` (`lineage` ASC);")
    #     odbcSetAutoCommit(con, autoCommit = TRUE)

    dbSendQuery(con, "DROP TABLE IF EXISTS hierarchy")
    dbWriteTable(con, value = df, name = "hierarchy", append = TRUE, row.names = TRUE) 
    dbCommit(conn = con)
    dbSendQuery(con, "ALTER TABLE `hierarchy` ADD COLUMN `id` VARCHAR(255) NOT NULL;")
    dbSendQuery(con, "UPDATE `hierarchy` SET `id` = `row_names`;")
    dbSendQuery(con, "ALTER TABLE `hierarchy` ENGINE = MyISAM;")
    dbSendQuery(con, "ALTER TABLE `hierarchy` CHANGE COLUMN `partition` `partition` VARCHAR(255) NULL DEFAULT NULL;")
    #Will do this later
    #dbSendQuery(con, "ALTER TABLE `hierarchy` ADD INDEX `name_idx` USING HASH (`label` ASC);")
    dbSendQuery(con, "ALTER TABLE `hierarchy` ADD INDEX `location_idx` (`partition` ASC, `start` ASC, `end` ASC);")
    dbSendQuery(con, "ALTER TABLE `hierarchy` ADD FULLTEXT INDEX `lineage_idx` (`lineage` ASC);")
    
    dbCommit(conn = con)
    },

  .saveValues=function(con) {
    counts = .counts
    countsIndices = expand.grid(seq(dim(counts)[1]), seq(dim(counts)[2]))

    df = data.frame(row=countsIndices[,1], col=countsIndices[,2], val=as.vector(counts))
    df = df[df$val != 0,]

    print(head(df))

    #   odbcSetAutoCommit(con, autoCommit = FALSE)
    #   sqlDrop(con, "values", errors=FALSE)
    #   sqlSave(con, df, tablename="values", addPK=TRUE)
    #   sqlQuery(con, "ALTER TABLE `values` CHANGE COLUMN `rownames` `id` BIGINT NOT NULL  ;")
    #   sqlQuery(con, "ALTER TABLE `values` ADD INDEX `rowcol_idx` USING BTREE (`row` ASC, `col` ASC) ;")
    #   odbcSetAutoCommit(con, autoCommit = TRUE)
    
    #requires user to change name of table to `values` after running this function
    dbSendQuery(con, "DROP TABLE IF EXISTS meta_values")
    dbWriteTable(con, value = df, name = "meta_values", append = TRUE, row.names = TRUE ) 
    dbCommit(conn = con)
    
    dbSendQuery(con, "ALTER TABLE `meta_values` CHANGE COLUMN `row_names` `id` BIGINT NOT NULL  ;")
    dbSendQuery(con, "ALTER TABLE `meta_values` ADD INDEX `rowcol_idx` USING BTREE (`row` ASC, `col` ASC) ;")
    
    dbCommit(conn = con)

  },

  .saveLevels=function(con) {
    df = data.frame(depth=seq(length(.levels)) - 1, label=.levels)

    print(head(df))
    
    #   odbcSetAutoCommit(con, autoCommit = FALSE)
    #   sqlDrop(con, "levels", errors=FALSE)
    #   sqlSave(con, df, tablename="levels", addPK=FALSE)
    #   sqlQuery(con, "ALTER TABLE `levels` ENGINE = MEMORY ;")
    #   sqlQuery(con, "ALTER TABLE `levels` DROP COLUMN `rownames` ;")
    #   odbcSetAutoCommit(con, autoCommit = TRUE)
    
    dbSendQuery(con, "DROP TABLE IF EXISTS levels")
    dbWriteTable(con, value = df, name = "levels", append = TRUE , row.names = TRUE) 
    dbCommit(conn = con)
    
    #Do these need to be run?
    #dbSendQuery(con, "ALTER TABLE `levels` ENGINE = MEMORY ;")
    #dbSendQuery(con, "ALTER TABLE `levels` DROP COLUMN `rownames` ;")
    #dbCommit(conn = con)
  }
)
