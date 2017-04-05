#' Data container for MRexperiment objects
#' 
#' Used to serve metagenomic data (used in e.g., icicle plots and heatmaps). Wraps
#' \code{\link[metagenomeSeq]{MRexperiment-class}} objects.
#' @importClassesFrom epivizrData EpivizData
#' @importFrom methods new
#' @importFrom vegan diversity
#' @import data.table
#' @exportClass EpivizMetagenomicsData
EpivizMetagenomicsData <- setRefClass("EpivizMetagenomicsData",
  contains = "EpivizData",
  fields = list(
    .levels = "ANY",
    .maxDepth = "numeric",
    .feature_order = "character",
    .minValue = "numeric",
    .maxValue = "numeric",
    .sampleAnnotation = "ANY",
    .nodeSelections = "ANY",
    .levelSelected = "ANY",
    .lastRootId = "character",
    
    # Tables
    .leaf_sample_count_table = "ANY",
    .graph = "ANY"
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

      log <- control$log
      norm <- control$norm
      
      # validate MRexperiment object
      MRExpCheck <- validateObject(object)
      
      if (!MRExpCheck) {
        stop("Incompatible MRexperiment objects")
      } else {
        message("MRExperiment Object validated... PASS")
      }
      .self$.feature_order <- feature_order
      if(is.null(feature_order)) {
        .self$.feature_order = colnames(fData(object))
      }
      else {
        .self$.feature_order <- feature_order
      }

      .self$.sampleAnnotation <- pData(object)

      print("creating leaf_sample_count_table")
      .self$.leaf_sample_count_table <- .create_leaf_sample_count_table(object)

      .self$.graph <- buildMetavizGraph(object, feature_order=feature_order)
      
      .self$.minValue <- min(.self$.leaf_sample_count_table[, !c("leaf"), with=FALSE])
      .self$.maxValue <- max(.self$.leaf_sample_count_table[, !c("leaf"), with=FALSE])
      .self$.nodeSelections <- list()
      .self$.levelSelected <- 3
      .self$.lastRootId <- "1-1"
      
      featureSelection = control$featureSelection
      
      if(!is.null(featureSelection)){
        temp_selections = list()
        featureSelection <- featureSelection[which(featureSelection != "no_match")]
        for(i in seq(1,length(featureSelection))){
            temp_selections[[as.character(.self$.graph$.nodes_table[node_label==names(featureSelection)[i],child])]] <- unname(featureSelection)[i]
        }
        .self$.nodeSelections = temp_selections
      }
      
      callSuper(object=object, ...)
    },

    # Create leaf_sample_count data.table
    
    .create_leaf_sample_count_table=function(obj_in){
      normed_counts <- as.data.frame(MRcounts(cumNorm(obj_in), norm = TRUE))
      normed_counts[["leaf"]] <- rownames(normed_counts)
      ret_table <- as.data.table(normed_counts)
      return(ret_table)
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
    nrow(.self$.graph$.hierarchy_tree)
  },
  nmeasurements=function() {
    ncol(.self$.leaf_sample_count_table)-1
  },

  nlevels=function() {
    nrow(unique(.self$.graph$.nodes_table[,level]))
  }
)

# Epiviz Websockets Protocol
EpivizMetagenomicsData$methods(
  get_default_chart_type = function() { "epiviz.ui.charts.tree.Icicle" },
  get_measurements=function() {
    samplesToRet <- colnames(.self$.leaf_sample_count_table)
    samplesToRet <- samplesToRet[-(length(samplesToRet))]
    out <- lapply(samplesToRet, function(sample) {
     epivizrData:::EpivizMeasurement(id=sample,
          name=sample,
          type="feature",
          datasourceId=.self$.id,
          datasourceGroup=.self$.id,
          defaultChartType="heatmap",
          annotation="anno",
          minValue=.self$.minValue,
          maxValue=.self$.maxValue,
          metadata=c("colLabel", "ancestors", "lineage", "label"))
    })
    return(out)
  },
  
  row_to_dict=function(row){
    toRet = list()
    toRet['end'] = row['end']
    toRet['partition'] = "NA"
    toRet['leafIndex'] = row['leafIndex']
    toRet['nchildren'] = row['nchildren']
    toRet['label'] = row['label']
    toRet['name'] = row['label']
    toRet['start'] = row['start']
    toRet['depth'] = row['depth']
    toRet['globalDepth'] = row['depth']
    toRet['nleaves'] = row['nleaves']
    toRet['parentId'] = row['parentId']
    toRet['order'] = row['order']
    toRet['id'] = row['id']
    if(toRet['id'] %in% names(.self$.nodeSelections)){
      toRet['selectionType'] = .self$.nodeSelections[[as.character(toRet['id'])]]
    } else{
      toRet['selectionType'] = 1
    }
    toRet['taxonomy'] = row['taxonomy']
    toRet['size'] = 1
    toRet['children'] = NULL
    return(toRet)
  },
  
  df_to_tree=function(root, df){
    
    children = df[which(df['parentId'] == as.character(unlist(root['id']))),]
    
    if(length(children) == 0){
      root$children = NULL
      return(root)
    }
    
    otherChildren = df[which(df['parentId'] != as.character(unlist(root['id']))),]
    
    children = children[order(children['order']),]
    
    if(nrow(children) > 0){
      for(row_index in seq(1,nrow(children))){
        childDict = row_to_dict(children[row_index,])
        subDict = df_to_tree(childDict, otherChildren)
        
        if(!is.null(subDict)){
          root$children[[row_index]] = subDict
        }
        else {
          root$children = NULL
        }
      }
    }
    return(root)
  },
  
  getHierarchy=function(nodeId = NULL) {
    # getHierarchy can be called with NULL from App
    if(is.null(nodeId) || nodeId == ""){
      nodeId <- .self$.lastRootId
    }
    .self$.lastRootId <- nodeId
    root <- nodeId
  
    #Split the node id to get level and index
    split_res <- strsplit(nodeId, "-")[[1]]
    level <- as.integer(split_res[1])
    index <- which(.self$.graph$.node_ids_DT[,level, with=FALSE] == nodeId)
    
    label <- as.character(unique(.self$.graph$.hierarchy_tree[,level][index]))
    taxonomy <- colnames(.self$.graph$.hierarchy_tree)[level]
    
    level_one_down <- level+1
    level_two_down <- level+2
    level_three_down <- level+3
    
    hierarchy_slice <- unique(.self$.graph$.node_ids_DT[get(taxonomy)==nodeId, level_one_down:level_three_down])
    nodes_level_1 <- unname(unlist(unique(hierarchy_slice[,1])))
    nodes_level_2 <- unname(unlist(unique(hierarchy_slice[,2])))
    nodes_level_3 <- unname(unlist(unique(hierarchy_slice[,3])))
    
    if(!(level == "1" || nodeId == "1-1")){
      parent_of_root_taxonomy <- colnames(.self$.graph$.hierarchy_tree)[(level-1)]
      parent_of_root <- unique(.self$.graph$.node_ids_DT[get(taxonomy)==nodeId, get(parent_of_root_taxonomy)])
      nodesToRet <- c(parent_of_root,root, nodes_level_1, nodes_level_2, nodes_level_3)
    } else{
      nodesToRet <- c(root, nodes_level_1, nodes_level_2, nodes_level_3)
    }
    
    num_rows <- length(nodesToRet)
    
    starts <- rep(1, num_rows)
    labels <- rep(1, num_rows)
    leafIndexes <- rep(1, num_rows)
    parentIds <- rep(1, num_rows)
    depths <- rep(1, num_rows)
    partitions <- rep(1, num_rows)
    ends <- rep(1, num_rows)
    ids <- rep(1, num_rows)
    nchildrens <- rep(1, num_rows)
    taxonomys <- rep(1, num_rows)
    nleaves <- rep(1, num_rows)
    orders <- rep(1, num_rows)
    
    leaf_ordering_table <- as.data.table(.self$.graph$.hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    
    for(i in seq(1,num_rows)){
      nodeId <- nodesToRet[i]
      split_res <- strsplit(nodesToRet[i], "-")[[1]]
      level <- as.integer(split_res[1])
      depths[i] <- level
      
      index <- which(.self$.graph$.node_ids_DT[,level,with=FALSE] == nodeId)
      
      label <- as.character(unique(.self$.graph$.hierarchy_tree[,level][index]))
      labels[i] <- label
      taxonomy <- colnames(.self$.graph$.hierarchy_tree)[level]

      if(nodesToRet[i] != "1-1"){
        parentId_taxonomy <- colnames(.self$.graph$.hierarchy_tree)[(level-1)]
        parentId <- unique(.self$.graph$.node_ids_DT[get(taxonomy)==nodesToRet[i], get(parentId_taxonomy)])
        parentIds[i] <- parentId
      } else{
        parentIds[i] <- "NA"
      }
      
      partition <- "NA"
      partitions[i] <- partition
      
      list_of_leaves <- .self$.graph$.leaf_of_table[node_label==label,leaf]
      leaf_indexes_temp <- leaf_ordering_table[leaf %in% list_of_leaves, otu_index]
      
      if(length(leaf_indexes_temp) > 0){
        start <- min(leaf_indexes_temp)
      }    else{
        start <- nodesToRet[i]
      }
      
      if(length(leaf_indexes_temp) > 0){
        end <- max(leaf_indexes_temp)
      }    else{
        end <- nodesToRet[i]
      }
      
      starts[i] <- start
      
      leafIndex <- start
      leafIndexes[i] <- leafIndex
      ends[i] <- end
      
      id <- nodesToRet[i]
      ids[i] <- id
      
      taxonomy <- colnames(.self$.graph$.hierarchy_tree)[level]
      taxonomys[i] <- taxonomy
      
      nchildren <- length(unique(.self$.graph$.node_ids_DT[get(taxonomy)==nodesToRet[i],][[as.integer(level)+1]]))
      nchildrens[i] <- nchildren
      
      nleaves_temp <- length(unname(unlist(unique(.self$.graph$.leaf_of_table[node_label==label, leaf]))))
      nleaves[i] <- nleaves_temp
      
      if(nodesToRet[i] != "1-1"){
        orders[i] <- .self$.graph$.nodes_table[get("child")==nodesToRet[i],get("order")][[1]]
      } else {
        orders[i] <- 1
      }
    }
    ret_data_frame <- data.frame(start = starts, label = labels, leafIndex = leafIndexes, parentId = parentIds, 
                                 depth = depths, partition = partitions, end = ends, id = ids, nchildren = nchildrens, 
                                 taxonomy = taxonomys, nleaves = nleaves, order = orders)
    
    if(length(ret_data_frame) > 0){
      # convert columns to int
      ret_data_frame['start'] = as.numeric(unlist(ret_data_frame['start']))
      ret_data_frame['end'] = as.numeric(unlist(ret_data_frame['end']))
      ret_data_frame['order'] = as.numeric(unlist(ret_data_frame['order']))
      ret_data_frame['leafIndex'] = as.numeric(unlist(ret_data_frame['leafIndex']))
      ret_data_frame['nchildren'] = as.numeric(unlist(ret_data_frame['nchildren']))
      ret_data_frame['nleaves'] = as.numeric(unlist(ret_data_frame['nleaves']))
      ret_data_frame['depth'] = as.numeric(unlist(ret_data_frame['depth']))
      ret_data_frame['id'] = as.character(unlist(ret_data_frame['id']))
      

      root = ret_data_frame[1,]
      rest = ret_data_frame[-1,]
      rootDict = row_to_dict(root)
      result = df_to_tree(rootDict, rest)

      resultResp = list()
      resultResp[['nodeSelectionTypes']] <- .self$.nodeSelections
      resultResp[['selectionLevel']] <- .self$.levelSelected
      resultResp[['tree']] <- result
      return(resultResp)
    }

    return(ret_data_frame)
  },
  
  propagateHierarchyChanges=function(selection = NULL, order = NULL, selectedLevels = NULL, request_with_labels = FALSE) {
    
    if(request_with_labels && !is.null(selection)){
      temp_selections = list()
      for(i in seq(1,length(selection))){
        temp_selections[[.self$.graph$.nodes_table[node_label==names(selection)[i],child]]] <- selection[i]
      }
      selection <- temp_selections
    }
    
    # update node selections types to metaviztree
    if(!is.null(selection)) {
      for(n in names(selection)){
        .self$.nodeSelections[[n]] = selection[[n]]
      }
    }
    
    if(!is.null(selectedLevels)) {
      .self$.levelSelected <- as.integer(names(selectedLevels)[1])
    }
    .self$.mgr$.clear_datasourceGroup_cache(.self)
  },
  
  getRows=function(measurements = NULL, start = 1, end = 1000, selectedLevels = 3, selections = NULL) {
    nodes_at_level <- .self$.graph$.nodes_table[level==selectedLevels, ]
    nodes_at_level_ids <- nodes_at_level[,child]
    
    if(!is.null(selections) && !(length(selections) == 0)){
      nodes_at_level_selections <- rep(2, length(nodes_at_level_ids))
      names(nodes_at_level_selections) <- nodes_at_level_ids
      selections <- c(selections, nodes_at_level_selections)
      
      expand_selections <- which(selections == 1)
      if(!is.null(expand_selections) && length(expand_selections) > 0){
        selections <- selections[-expand_selections]
      }
      
      child_lineage <- .self$.graph$.nodes_table[child %in% names(selections),]
      remove_selections <- which(selections == 0)
      if(length(remove_selections) > 0){
        kept_nodes <- child_lineage[!grepl(paste(paste(names(remove_selections), collapse=",|"), ",",sep=""), lineage),]
        kept_nodes <- kept_nodes[!(child %in% names(remove_selections)),]
      } else {
        kept_nodes <- child_lineage
      }
      
      agg_selections <- which(selections == 2)
      if(length(agg_selections) > 0){
        kept_nodes <- as.character(kept_nodes[!grepl(paste(paste(names(agg_selections), collapse=",|"), ",",sep=""), lineage), child])
      }
      nodes_at_level <- .self$.graph$.nodes_table[child %in% kept_nodes,]
    }

    leaf_ordering_table <- as.data.table(.self$.graph$.hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    leaf_ordering_table <- leaf_ordering_table[,leaf:=as.character(leaf)]
    
    first_join <- merge(leaf_ordering_table, merge(nodes_at_level, .self$.graph$.leaf_of_table, by="node_label"), by="leaf")

    data_rows = list()
    nodes <- nodes_at_level[,node_label]

    ends <- list()
    starts <- list()
    indexes <- list()
    metadata <- list()
    metadata[['label']] <- list()

    for (i in seq(1, length(nodes))) {
      row_info <- list()
      feature_label <- nodes[i]
      res <- first_join[node_label==feature_label,otu_index]
      if(length(res) > 0) {
        ends[i] <- max(res)
        starts[i] <- min(res)
        indexes[i] <- min(res)
        metadata[['label']][i] <- feature_label
        node_level <- as.character(.self$.graph$.nodes_table[node_label==nodes[i], level])
        list_of_leaves <- .self$.graph$.leaf_of_table[node_label==nodes[i],leaf]
        leaf_indexes <- leaf_ordering_table[leaf %in% list_of_leaves, otu_index]
        id <- paste(as.hexmode(as.integer(node_level)), as.hexmode(min(leaf_indexes)), sep="-")
        metadata[['id']][i] <- id
        lineage <- .self$.graph$.nodes_table[node_label==nodes[i], lineage]
        metadata[['lineage']][i] <- lineage
      }
    }
    
    data_rows[['end']] <- ends
    data_rows[['start']] <- starts
    data_rows[['index']] <- indexes
    data_rows[['metadata']] <- metadata
    return(data_rows)
  },
  
  getValues=function(measurements = NULL, start = 1, end = 1000, selectedLevels = 3, selections = NULL) {

    nodes_at_level <- .self$.graph$.nodes_table[level==selectedLevels,]
    nodes_at_level_ids <- nodes_at_level[,child]
    
    if(!is.null(selections) && !(length(selections) == 0)){
      nodes_at_level_selections <- rep(2, length(nodes_at_level_ids))
      names(nodes_at_level_selections) <- nodes_at_level_ids
      selections <- c(selections, nodes_at_level_selections)
      
      expand_selections <- which(selections == 1)
      if(!is.null(expand_selections) && length(expand_selections) > 0){
        selections <- selections[-expand_selections]
      }
      
      child_lineage <- .self$.graph$.nodes_table[child %in% names(selections),]
      remove_selections <- which(selections == 0)
      if(length(remove_selections) > 0){
        kept_nodes <- child_lineage[!grepl(paste(paste(names(remove_selections), collapse=",|"), ",",sep=""), lineage),]
        kept_nodes <- kept_nodes[!(child %in% names(remove_selections)),]
      } else {
        kept_nodes <- child_lineage
      }
      
      agg_selections <- which(selections == 2)
      if(length(agg_selections) > 0){
        kept_nodes <- as.character(kept_nodes[!grepl(paste(paste(names(agg_selections), collapse=",|"), ",",sep=""), lineage), child])
      }
      nodes_at_level <- .self$.graph$.nodes_table[child %in% kept_nodes,]
    }
    
    data_columns = list()
    leaf_sample_count_table_sub_select <- na.omit(.self$.leaf_sample_count_table[,mget(c(measurements,"leaf"))])
    first_join <- unique(merge(na.omit(.self$.graph$.leaf_of_table),na.omit(nodes_at_level), by="node_label"))
    
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,leaf:=as.character(leaf)]
    first_join <- first_join[,leaf:=as.character(leaf)]
    
    second_join <- merge(na.omit(leaf_sample_count_table_sub_select), na.omit(first_join), by="leaf", all.x = TRUE)
    results <- second_join[,lapply(.SD, sum), .SDcols=measurements, by=node_label]
    close_results <- as.data.frame(results)
    
    data_columns = list()
    for(m in measurements){
      inner_result <- close_results[,m]
      data_columns[[m]] <- inner_result
    }
    return(data_columns)
  },
  
  getCombined=function(measurements = NULL, 
                          seqName, start = 1, end = 1000, 
                          order = NULL, nodeSelection, selectedLevels) {
    
    # update selectedLevels on taxonomy tree
    
    # update node selections types to metaviztree
    if(!is.null(nodeSelection)) {
      for(n in names(nodeSelection)){
        .self$.nodeSelections[[n]] = nodeSelection[[n]]
      }
    }
    
    selectedLevels = .self$.levelSelected
    selections = .self$.nodeSelections

    data_columns = getValues(measurements, start, end, selectedLevels, selections)
    data_rows = getRows(measurements, start, end, selectedLevels, selections)

    result <- list(
      cols = data_columns,
      rows = data_rows,
      globalStartIndex = data_rows$start[[1]]
    )

    return(result)
  },
  
  searchTaxonomy=function(query, max_results) {
    
    results = list()
    return(results)
  },
  
  getPCA=function(measurements, seqName = '', start, end) {
    init <- as.data.frame(.self$.leaf_sample_count_table[,mget(measurements)])
    if("leaf" %in% colnames(init)){
      init <- init[,-which(colnames(init)=="leaf")]
    }
    x <- t(init)
    df <- log2(x+1)
    ord <- prcomp(df)

    res <- apply(ord$x, 1, function(x) {
      return (list(PC1 = unlist(x["PC1"]), PC2 = unlist(x["PC2"])))
    })

    data <- list()
    for (row in names(res)) {
      temp <- list()
      temp$sample_id = row
      temp$PC1 = unname(res[[row]]$PC1)
      temp$PC2 = unname(res[[row]]$PC2)
      annotation = as.list(.self$.sampleAnnotation[row,])
      for (anno in names(annotation)) {
        temp[[anno]] = annotation[[anno]]
      }

      data[[row]] <- temp
    }

    result <- list()
    result$data = unname(data)
    result$pca_variance_explained = ord$sdev[1:2]
    return(result)
  },
  
  getAlphaDiversity=function(measurements, seqName = '', start, end) {
    df <- as.data.frame(.self$.leaf_sample_count_table[,mget(measurements)])
    if("leaf" %in% colnames(df)){
      df <- df[,-which(colnames(df)=="leaf")]
    }
    alpha_diversity <- vegan::diversity(t(df), index = "shannon")

    data <- list()
    for (row in 1:length(alpha_diversity)) {
      temp <- list()
      temp$sample_id = colnames(df)[row]
      temp$alphaDiversity = unname(alpha_diversity[row])
      annotation = as.list(.self$.sampleAnnotation[temp$sample_id,])
      for (anno in names(annotation)) {
        temp[[anno]] = annotation[[anno]]
      }

      data[[row]] <- temp
    }

    result <- list()
    result$data = unname(data)
    return(result)
  }
)