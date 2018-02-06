#' Data container for MRexperiment objects
#' 
#' Used to serve metagenomic data (used in e.g., icicle plots and heatmaps). Wraps
#' \code{\link[metagenomeSeq]{MRexperiment-class}} objects.
#' @importClassesFrom epivizrData EpivizData
#' @importFrom vegan diversity
#' @import data.table
#' @import digest
#' @import methods
#' @import httr
#' @exportClass InnerNodesEpivizMetagenomicsData
#' @examples
#' 
#' \dontrun{
#' library(curatedMetagenomicData)
#' zeller.eset = ZellerG_2014.metaphlan_bugs_list.stool()
#' zeller_MR <- ExpressionSet2MRexperiment(zeller.eset)
#' feature_order <- colnames(fData(zeller_MR))
#' sampleId<- "CCIS98482370ST-3-0"
#' mObj <- metavizr:::InnerNodesEpivizMetagenomicsData$new(zeller_MR, feature_order = feature_order)
#' }
#'
InnerNodesEpivizMetagenomicsData <- setRefClass("InnerNodesEpivizMetagenomicsData",
  contains = "EpivizMetagenomicsData",
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
    .json_query = "ANY",
    .graph_feature_order = "character",
                                        
    # Tables
    .leaf_sample_count_table = "ANY",
    .leaf_sample_count_table_long = "ANY",
    .graph = "ANY"
  ),
                                      
  methods=list(
    initialize=function(object, columns=NULL, control=metavizControl(), feature_order=NULL,...) {
                                          
      # Initialize parameters used here
      aggregateAtDepth <- control$aggregateAtDepth
      maxDepth <- control$maxDepth
      maxHistory <- control$maxHistory
      maxValue <-  control$maxValue
      minValue <-  control$minValue
      aggregateFun <-  control$aggregateFun
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
                                          
      if(is.null(feature_order)) {
        .self$.feature_order = colnames(fData(object))
      } else {
        .self$.feature_order <- feature_order
      }
                                          
      if(norm){
        object <- cumNorm(object, p = .75)
      }
                                          
      .self$.sampleAnnotation <- pData(object)
                                          
      .self$.graph <- buildInnerNodesMetavizGraph(object, feature_order=feature_order)
      .self$.graph_feature_order <- .self$.graph$.feature_order
                                          
      message("creating leaf_sample_count_table")
      .self$.leaf_sample_count_table <- .create_leaf_sample_count_table(object, norm=norm)
                                        
      .self$.leaf_sample_count_table_long <- .create_leaf_sample_count_table_long(object, norm=norm)
                                          
      .self$.minValue <- min(.self$.leaf_sample_count_table[, !c("leaf", "start", "end"), with=FALSE])
      .self$.maxValue <- max(.self$.leaf_sample_count_table[, !c("leaf", "start", "end"), with=FALSE])
      .self$.nodeSelections <- list()
      .self$.levelSelected <- aggregateAtDepth
      .self$.lastRootId <- "0-0"
                                          
      featureSelection = control$featureSelection
                                          
      if(!is.null(featureSelection)){
        featureSelection <- featureSelection[which(names(featureSelection) != "NA")]
        featureSelection <- featureSelection[which(names(featureSelection) != "no_match")]
                                            
        node_ids <- sapply(names(featureSelection), function(n) {
          as.character(.self$.graph$.nodes_table[node_label==n,id])
        })
                                            
        temp_selections <- unname(featureSelection)
        names(temp_selections) <- node_ids
        .self$.nodeSelections = temp_selections
      }
                                          
      callSuper(object=object, ...)
    },
                                        
    # Create leaf_sample_count data.table
    .create_leaf_sample_count_table=function(obj_in, norm = TRUE){
      normed_counts <- as.data.frame(MRcounts(obj_in, norm=norm))
      leaf_names <- rownames(normed_counts)

      level_annotated <- rep(0, nrow(fData(obj_in)))
      f_data <- fData(obj_in)
      for(i in seq(1, nrow(f_data))){
        if (length(which(is.na(f_data[i,]))) == 0){
          level_annotated[i] <- length(colnames(f_data))
        } else {
          level_annotated[i] <- min(which(is.na(f_data[i,])))-2
        }
      }
      
      normed_counts[["leaf"]] <- leaf_names
      normed_counts[["level_annotated"]] <- level_annotated
      normed_counts <- normed_counts[.self$.graph$.hierarchy_tree_order,]
      normed_counts[["start"]] <- .self$.graph$.hierarchy_tree[,c("start")]
      normed_counts[["end"]] <- .self$.graph$.hierarchy_tree[,c("end")]
      ret_table <- as.data.table(normed_counts)
      return(ret_table)
    },
                                        
    .create_leaf_sample_count_table_long=function(obj_in, norm = TRUE){
      temp_table <- .self$.leaf_sample_count_table
                                        
      temp_table_long <- melt(temp_table, id.vars = c("leaf", "start", "end", "level_annotated"), 
                                measure.vars = c(colnames(temp_table)[1:(length(colnames(temp_table))-4)]), 
                                variable.name = "sample", variable.factor = FALSE)
                                          
      ret_table_long <- temp_table_long[value != 0.0,]
      
      return(ret_table_long)
    }
                                        

  )
)

# Data analysis features
InnerNodesEpivizMetagenomicsData$methods(
  nmeasurements=function() {
    ncol(.self$.leaf_sample_count_table)-4
  }
  

)

# Epiviz Websockets Protocol
InnerNodesEpivizMetagenomicsData$methods(
    get_measurements=function() {
    "Get all annotation info for all samples
    
    \\describe{
    \\item{chart_id_or_object}{An object of class \\code{EpivizChart} or an id for
    a chart loaded to the epiviz app.}
    }
    "
    samplesToRet <- colnames(.self$.leaf_sample_count_table)
    samplesToRet <- samplesToRet[-which(samplesToRet == "leaf")]
    samplesToRet <- samplesToRet[-which(samplesToRet == "level_annotated")]
    samplesToRet <- samplesToRet[-which(samplesToRet == "start")]
    samplesToRet <- samplesToRet[-which(samplesToRet == "end")]
    out <- lapply(samplesToRet, function(sample) {
      epivizrData:::EpivizMeasurement(id=sample,
        name=sample,
        type="feature",
        datasourceId=.self$.id,
        datasourceGroup=.self$.id,
        defaultChartType="heatmap",
        annotation=as.list(.sampleAnnotation[sample,]),
        minValue=.self$.minValue,
        maxValue=.self$.maxValue,
        metadata=c("colLabel", "ancestors", "lineage", "label"))
    })
    return(out)
  },
  

  getHierarchy=function(nodeId = NULL) {
    "Retrieve feature hierarchy information for subtree with specified root
    
    \\describe{
    \\item{nodeId}{Feature identifier with level info}
    }
    "
    
    # getHierarchy can be called with NULL from App
    if(is.null(nodeId) || nodeId == ""){
      nodeId <- .self$.lastRootId
    }
    .self$.lastRootId <- nodeId
    root <- nodeId
    
    #Split the node id to get level and index
    split_res <- strsplit(nodeId, "-")[[1]]
    level <- as.integer(split_res[1])+1
    index <- which(.self$.graph$.node_ids_table[,level, with=FALSE] == nodeId)
    
    graph_tree <- .self$.graph$.hierarchy_tree[,-which(colnames(.self$.graph$.hierarchy_tree) == "start")]
    graph_tree <- graph_tree[,-which(colnames(.self$.graph$.hierarchy_tree) == "end")]
    
    label <- as.character(unique(graph_tree[,level][index]))
    taxonomy <- colnames(graph_tree)[level]
    
    if(length(.self$.graph$.feature_order) >= level+3){
      last_level_of_subtree <- level+3
    } else{
      last_level_of_subtree <- length(.self$.graph$.feature_order)
    }
    
    hierarchy_slice <- unique(.self$.graph$.node_ids_table[get(taxonomy)==nodeId, (level+1):last_level_of_subtree])
    
    nodes_of_subtree <- sapply(seq(1,length((level+1):last_level_of_subtree)), function(i) {
      unname(unlist(unique(hierarchy_slice[,i, with=FALSE])))
    })
    
    nodes_of_subtree <- unlist(nodes_of_subtree)
    nodes_of_subtree <- nodes_of_subtree[which(!is.na(nodes_of_subtree))]
    
    if(level == 0 || nodeId == "0-0"){
      nodesToRet <- c(root, unlist(nodes_of_subtree))
    } else{
      parent_of_root_taxonomy <- colnames(graph_tree)[(level-1)]
      parent_of_root <- unique(.self$.graph$.node_ids_table[get(taxonomy)==nodeId, get(parent_of_root_taxonomy)])
      nodesToRet <- c(parent_of_root,root, unlist(nodes_of_subtree))
    }
    
    num_rows <- length(nodesToRet)
    
    starts <- rep(1, num_rows)
    labels <- rep(1, num_rows)
    leafIndexes <- rep(1, num_rows)
    parentIds <- rep(1, num_rows)
    depths <- rep(0, num_rows)
    partitions <- rep(1, num_rows)
    ends <- rep(1, num_rows)
    ids <- rep(1, num_rows)
    nchildrens <- rep(1, num_rows)
    taxonomys <- rep(1, num_rows)
    nleaves <- rep(1, num_rows)
    orders <- rep(1, num_rows)
    
    for(i in seq_len(num_rows)){
      if(as.integer(strsplit(nodesToRet[i], "-")[[1]][1]) == last_level_of_subtree){
        depths[i] = length(.self$.graph$.feature_order)
        level = length(.self$.graph$.feature_order)
        index <- which(.self$.graph$.node_ids_table[,level,with=FALSE] == nodeId)
        
        label <- as.character(unique(graph_tree[,level][index]))
        labels[i] <- label
        
        partition <- "NA"
        partitions[i] <- partition
        
        starts[i] <- .self$.graph$.nodes_table[id == nodesToRet[i],start]
        leafIndexes[i] <- .self$.graph$.nodes_table[id == nodesToRet[i],start]
        ends[i] <- .self$.graph$.nodes_table[id == nodesToRet[i],end]
        
        ids[i] <- nodeId
        
        taxonomy <- colnames(graph_tree)[level]
        taxonomys[i] <- taxonomy
        
        nchildrens[i] <- 0
        nleaves[i] <- 0
        
        orders[i] <- .self$.graph$.nodes_table[get("id")==nodeId,get("order")][[1]]
      } else{
        nodeId <- nodesToRet[i]
        split_res <- strsplit(nodesToRet[i], "-")[[1]]
        depths[i] <- as.integer(split_res[1])
        level <- as.integer(split_res[1])+1
        
        index <- which(.self$.graph$.node_ids_table[,level,with=FALSE] == nodeId)
        
        label <- as.character(unique(graph_tree[,level][index]))
        labels[i] <- label
        taxonomy <- colnames(graph_tree)[level]
        
        if(nodesToRet[i] != "0-0"){
          parentId_taxonomy <- colnames(graph_tree)[(level-1)]
          parentId <- unique(.self$.graph$.node_ids_table[get(taxonomy)==nodesToRet[i], get(parentId_taxonomy)])[1]
          parentIds[i] <- parentId
        } else{
          parentIds[i] <- "NA"
        }
        
        partition <- "NA"
        partitions[i] <- partition
        
        start <- .self$.graph$.nodes_table[id == nodesToRet[i],start]
        end <- .self$.graph$.nodes_table[id == nodesToRet[i],end]
        
        starts[i] <- start
        leafIndexes[i] <- start
        ends[i] <- end
        
        id <- nodesToRet[i]
        ids[i] <- id
        
        taxonomy <- colnames(graph_tree)[level]
        taxonomys[i] <- taxonomy
        
        nchildren <- length(unique(.self$.graph$.node_ids_table[get(taxonomy)==nodesToRet[i],][[as.integer(level)+1]]))
        nchildrens[i] <- nchildren[1]
        
        nleaves[i] <- (end - start)
        
        if(nodesToRet[i] != "0-0"){
          orders[i] <- .self$.graph$.nodes_table[get("id")==nodesToRet[i],get("order")][[1]]
        } else {
          orders[i] <- 1
        }
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
      
      result[["rootTaxonomies"]] = .self$.graph$.feature_order
      lineage = .self$.graph$.nodes_table[get("id")==nodesToRet[1],get("lineage")][[1]]
      
      lineageLabel <- sapply(strsplit(lineage, ",")[[1]], function(str_id) {
        .self$.graph$.nodes_table[get("id") == str_id, get("node_label")][[1]]
      })
      
      result[["lineageLabel"]] = paste(lineageLabel, sep=", ")
      
      resultResp = list(nodeSelectionTypes = .self$.nodeSelections, 
                        selectionLevel = .self$.levelSelected, 
                        tree = result)
      
      return(resultResp)
    }
    
    return(ret_data_frame)
  },
  
  propagateHierarchyChanges=function(selection = NULL, order = NULL, selectedLevels = NULL, request_with_labels = FALSE) {
    "Update internal state for hierarchy
    
    \\describe{
    \\item{selection}{Node-id and selectionType pairs}
    \\item{order}{Ordering of features}
    \\item{selectedLevels}{Current aggregation level}
    \\item{request_with_labels}{For handling requests using fData entries from MRexperiment}
    }
    "
    
    if(request_with_labels && !is.null(selection)){
      selection_ids <- sapply(names(selection), function(i){
        .self$.graph$.nodes_table[node_label==i,id]
      })
      names(selection) <- selection_ids
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
    "Return the sample annotation and features within the specified range and level for a given sample and features
    
    \\describe{
    \\item{measurements}{Samples to retrieve for}
    \\item{start}{Start of feature range to query}
    \\item{end}{End of feature range to query}
    \\item{selections}{Node-id and selectionType pairs}
    \\item{selectedLevels}{Current aggregation level}
    
    }
    "
    
    nodes_at_level <- .self$.graph$.nodes_table[level==selectedLevels, ]
    nodes_at_level_ids <- nodes_at_level[,id]
    
    if(!is.null(selections) && !(length(selections) == 0)){
      nodes_at_level_selections <- rep(2, length(nodes_at_level_ids))
      names(nodes_at_level_selections) <- nodes_at_level_ids
      selections <- c(selections, nodes_at_level_selections)
      
      expand_selections <- which(selections == 1)
      if(!is.null(expand_selections) && length(expand_selections) > 0){
        selections <- selections[-expand_selections]
      }
      
      child_lineage <- .self$.graph$.nodes_table[id %in% names(selections),]
      remove_selections <- which(selections == 0)
      if(length(remove_selections) > 0){
        kept_nodes <- child_lineage[!grepl(paste(paste(names(remove_selections), collapse=",|"), ",",sep=""), lineage),]
        kept_nodes <- kept_nodes[!(id %in% names(remove_selections)),]
      } else {
        kept_nodes <- child_lineage
      }
      
      agg_selections <- which(selections == 2)
      if(length(agg_selections) > 0){
        kept_nodes <- as.character(kept_nodes[!grepl(paste(paste(names(agg_selections), collapse=",|"), ",",sep=""), lineage), id])
      }
      nodes_at_level <- .self$.graph$.nodes_table[id %in% kept_nodes,]
    }
    
    nodes_in_range <- nodes_at_level[start >= start,]
    nodes_in_range <- nodes_in_range[end <= end,]
    setorderv(nodes_in_range, "start")
    
    nodes_in_range_df <- as.data.frame(nodes_in_range)
    ends <- rep(0, nrow(nodes_in_range_df))
    starts <- rep(0, nrow(nodes_in_range_df))
    indexes <- rep(0, nrow(nodes_in_range_df))
    metadata <- list()
    metadata[['label']] <- list()
    metadata[['id']] <- list()
    metadata[['lineage']] <- list()
    
    for (i in seq(1,nrow(nodes_in_range_df))) {
      feature_label <- nodes_in_range_df[i, "node_label"]
      ends[i] <- as.integer(nodes_in_range_df[i, "end"])
      starts[i] <- as.integer(nodes_in_range_df[i, "start"])
      indexes[i] <- starts[i]
      metadata[['label']][i] <- feature_label
      metadata[['id']][i] <- nodes_in_range_df[i,"id"]
      metadata[['lineage']][i] <- nodes_in_range_df[i,"lineage"]
    }
    
    data_rows = list(end=ends, start=starts, index =indexes, metadata=metadata)
    return(data_rows)
  },
  
  getValues=function(measurements = NULL, start = 1, end = 1000, selectedLevels = 3, selections = NULL) {
    "Return the counts for a sample within the specified range
    
    \\describe{
    \\item{measurements}{Samples to get counts for}
    \\item{start}{Start of feature range to query}
    \\item{end}{End of feature range to query}
    \\item{selections}{Node-id and selectionType pairs}
    \\item{selectedLevels}{Current aggregation level}
    
    }
    "
    
    nodes_at_level <- .self$.graph$.nodes_table[level==selectedLevels,]
    nodes_at_level_ids <- nodes_at_level[,id]
    
    if(!is.null(selections) && !(length(selections) == 0)){
      nodes_at_level_selections <- rep(2, length(nodes_at_level_ids))
      names(nodes_at_level_selections) <- nodes_at_level_ids
      selections <- c(selections, nodes_at_level_selections)
      
      expand_selections <- which(selections == 1)
      if(!is.null(expand_selections) && length(expand_selections) > 0){
        selections <- selections[-expand_selections]
      }
      
      child_lineage <- .self$.graph$.nodes_table[id %in% names(selections),]
      remove_selections <- which(selections == 0)
      if(length(remove_selections) > 0){
        kept_nodes <- child_lineage[!grepl(paste(paste(names(remove_selections), collapse=",|"), ",",sep=""), lineage),]
        kept_nodes <- kept_nodes[!(id %in% names(remove_selections)),]
      } else {
        kept_nodes <- child_lineage
      }
      
      agg_selections <- which(selections == 2)
      if(length(agg_selections) > 0){
        kept_nodes <- as.character(kept_nodes[!grepl(paste(paste(names(agg_selections), collapse=",|"), ",",sep=""), lineage), id])
      }
      nodes_at_level <- .self$.graph$.nodes_table[id %in% kept_nodes,]
    }
    
    leaf_sample_count_table_temp <- .self$.leaf_sample_count_table_long
    
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_temp[,node_label:=as.character(leaf)]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,node_label:=sapply(strsplit(leaf_sample_count_table_sub_select[,node_label], "__"), function(i){unname(i)[2]})]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,!"leaf"]
    
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[node_label %in% nodes_at_level[,node_label],]
    
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[start >= start,]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[end <= end,]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,!"start"]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,!"end"]
    
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[sample %in% measurements,]
    
    close_results <- as.data.frame(dcast(data = leaf_sample_count_table_sub_select, formula = node_label ~ sample, value.var = "value", fill = 0.0, fun=sum))
    zero_results <- setdiff(nodes_at_level[,node_label], close_results[,"node_label"])
    if(length(zero_results) > 0){
      for(k in seq(1, length(zero_results))){
        new_index <- nrow(close_results)+1
        close_results[new_index,] <- 0.0
        close_results[new_index,]["node_label"] <- zero_results[k]
      }
    }
    
    close_results <- close_results[order(close_results[,"node_label"]),]
    close_results[is.na(close_results)] <- 0.0
    rownames(close_results) <- seq(1,nrow(close_results))
    names_to_add <- names(close_results[,1])
    data_columns = list()
    for(m in measurements){
      if(m %in% colnames(close_results)){
        inner_result <- close_results[,m]
        names(inner_result) <- names_to_add
        data_columns[[m]] <- inner_result
      } else{
        inner_result <- rep(0.0, nrow(close_results))
        names(inner_result) <- names_to_add
        data_columns[[m]] <- inner_result
      }
    }
    return(data_columns)
  },
  
 
  searchTaxonomy=function(query = NULL, max_results = 15) {
    "Return list of features matching a text-based query
    
    \\describe{
    \\item{query}{String of feature for which to search}
    \\item{max_results}{Maximum results to return}
    
    }
    "
    
    if(is.null(query)){
      return(list())
    }
    
    nodes_table_lowercase <- .self$.graph$.nodes_table
    nodes_table_lowercase <- nodes_table_lowercase[,node_label:=tolower(node_label)]
    query_lowercase <- tolower(query)
    matching_nodes <- nodes_table_lowercase[grepl(query, node_label),]
    
    if(nrow(matching_nodes) > max_results){
      num_results <- max_results
    } else{
      num_results <- nrow(matching_nodes)
    }
    matching_nodes <- matching_nodes[,head(.SD, num_results)]
    
    node_labels <- matching_nodes[,node_label]
    node_ids <- matching_nodes[,id]
    levels <- matching_nodes[,level]
    starts <- length(node_labels)
    ends <- length(node_labels)
    
    leaf_ordering_table <- as.data.table(.self$.graph$.hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    
    leaf_table_lowercase <- .self$.graph$.leaf_of_table
    leaf_table_lowercase <- leaf_table_lowercase[,node_label:=tolower(node_label)]
    
    for(i in seq_along(node_labels)){
      node <- node_labels[i]
      
      list_of_leaves <- leaf_table_lowercase[node_label==node,leaf]
      leaf_indexes_temp <- leaf_ordering_table[leaf %in% list_of_leaves, otu_index]
      
      if(length(leaf_indexes_temp) > 0){
        start <- min(leaf_indexes_temp)
      } else{
        start <- node
      }
      
      if(length(leaf_indexes_temp) > 0){
        end <- max(leaf_indexes_temp)
      } else{
        end <- node
      }
      
      starts[i] <- start
      ends[i] <- end
    }
    
    results = list()
    for(i in seq_len(num_results)){
      results[[i]] <- list("gene"=node_labels[i], "start"=starts[i], 
                           "end"=ends[i], "seqName"="metavizr", 
                           "nodeId"=node_ids[i], "level"=levels[i])
    }
    return(results)
  },
  
  getPCA=function(measurements = NULL) {
    " Compute PCA over all features for given samples
    
    \\describe{
    \\item{measurements}{Samples to compute PCA over}
    \\item{start}{Start of feature range to query }
    \\item{end}{End of feature range to query}
    }
    "
    
    if(is.null(measurements)){
      samples <- colnames(.self$.leaf_sample_count_table)
      samples <- samples[-(which(samples == "otu_index"))]
      measurements <- samples[-(which(samples == "leaf"))]
    }
    
    init <- as.data.frame(.self$.leaf_sample_count_table[,mget(measurements)])
    if("leaf" %in% colnames(init)){
      init <- init[,-which(colnames(init)=="leaf")]
    }
    x <- init
    df <- log2(x+1)
    ord <- prcomp(df)
    
    data <- list()
    for (row in rownames(ord$rotation)) {
      temp <- list(sample_id = row, PC1 = unname(ord$rotation[row,][1]), PC2 = unname(ord$rotation[row,][2]))
      annotation = as.list(.self$.sampleAnnotation[row,])
      for (anno in names(annotation)) {
        temp[[anno]] = annotation[[anno]]
      }
      
      data[[row]] <- temp
    }
    
    result <- list(data = unname(data), pca_variance_explained = ord$sdev[1:2])
    return(result)
  },
  
  getAlphaDiversity=function(measurements = NULL) {
    " Compute alpha diversity using vegan for the given samples
    
    \\describe{
    \\item{measurements}{Samples to compute alpha diversity}
    \\item{start}{Start of feature range to query }
    \\item{end}{End of feature range to query}
    }
    "
    
    if(is.null(measurements)){
      samples <- colnames(.self$.leaf_sample_count_table)
      samples <- samples[-(which(samples == "otu_index"))]
      measurements <- samples[-(which(samples == "leaf"))]
    }
    
    df <- as.data.frame(.self$.leaf_sample_count_table[,mget(measurements)])
    if("leaf" %in% colnames(df)){
      df <- df[,-which(colnames(df)=="leaf")]
    }
    alpha_diversity <- vegan::diversity(t(df), index = "shannon")
    
    data <- list()
    for (row in seq_along(alpha_diversity)) {
      temp <- list(sample_id = colnames(df)[row], alphaDiversity = unname(alpha_diversity[row]))
      annotation = as.list(.self$.sampleAnnotation[temp$sample_id,])
      for (anno in names(annotation)) {
        temp[[anno]] = annotation[[anno]]
      }
      
      data[[row]] <- temp
    }
    
    result <- list(data = unname(data))
    return(result)
  }
)
