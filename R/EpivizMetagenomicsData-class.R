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
#' @exportClass EpivizMetagenomicsData		
#' 
#' @examples
#' \dontrun{
#' library(metagenomeSeq)
#' data(mouseData)
#' obj <- metavizr:::EpivizMetagenomicsData$new(mouseData, feature_order = colnames(fData(mouseData)))
#' }
#' 
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
    .json_query = "ANY",
    .graph_feature_order = "character",
    
    # Tables
    .leaf_sample_count_table = "ANY",
    .leaf_sample_count_table_long = "ANY",
    .graph = "ANY"
  ),
  
  methods=list(
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

      if(is.null(feature_order)) {
        .self$.feature_order = colnames(fData(object))
      }
      else {
        .self$.feature_order <- feature_order
      }

      if(norm){
          object <- cumNorm(object, p = .75)
      }
      
      .self$.sampleAnnotation <- pData(object)

      .self$.graph <- buildMetavizGraph(object, feature_order=feature_order)
      .self$.graph_feature_order <- .self$.graph$.feature_order
      
      message("creating leaf_sample_count_table")
      .self$.leaf_sample_count_table <- .create_leaf_sample_count_table(object, norm=norm)
      
      .self$.leaf_sample_count_table_long <- .create_leaf_sample_count_table_long(object, norm=norm)
      
      .self$.minValue <- min(.self$.leaf_sample_count_table[, !c("leaf", "otu_index"), with=FALSE])
      .self$.maxValue <- max(.self$.leaf_sample_count_table[, !c("leaf", "otu_index"), with=FALSE])
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
      replacing_na_obj_fData <- fData(obj_in)[,.self$.feature_order]
      na_indices <- which(is.na(leaf_names))
      for(j in seq(1, length(na_indices))){
          leaf_names[na_indices[j]] <- paste("Not_Annotated", .self$.feature_order[length(.self$.feature_order)], replacing_na_obj_fData[,.self$.feature_order[1]][na_indices[j]], sep="_")
      }
      na_indices <- which(leaf_names == "NA")
      for(j in seq(1, length(na_indices))){
          leaf_names[na_indices[j]] <- paste("Not_Annotated", .self$.feature_order[length(.self$.feature_order)], replacing_na_obj_fData[,.self$.feature_order[1]][na_indices[j]], sep="_")
      }
      null_indices <- which(leaf_names == "NULL")
      for(j in seq(1, length(null_indices))){
          leaf_names[null_indices[j]] <- paste("Not_Annotated", .self$.feature_order[length(.self$.feature_order)], replacing_na_obj_fData[,.self$.feature_order[1]][null_indices[j]], sep="_")
      }
      normed_counts[["leaf"]] <- leaf_names
      normed_counts[["otu_index"]] <- .self$.graph$.hierarchy_tree[,c("otu_index")]
      ret_table <- as.data.table(normed_counts)
      return(ret_table)
    },
    
    .create_leaf_sample_count_table_long=function(obj_in, norm = TRUE){
      temp_table <- .self$.leaf_sample_count_table
      
      temp_table_long <- melt(temp_table, id.vars = c("leaf", "otu_index"), 
                             measure.vars = c(colnames(temp_table)[1:(length(colnames(temp_table))-2)]), 
                             variable.name = "sample", variable.factor = FALSE)
      
      ret_table_long <- temp_table_long[value != 0.0,]
      setorderv(ret_table_long, "otu_index")
      return(ret_table_long)
    },
    
    update=function(newObject, ...) {
      callSuper(newObject, ...)
    },
    plot=function(...) {
    }
  )
)

# Data analysis features
EpivizMetagenomicsData$methods(
  nleaves=function() {
    nrow(.self$.graph$.hierarchy_tree)
  },
  nmeasurements=function() {
    ncol(.self$.leaf_sample_count_table)-2
  },

  nlevels=function() {
    nrow(unique(.self$.graph$.nodes_table[,level]))
  }
)

# Epiviz Websockets Protocol
EpivizMetagenomicsData$methods(
  get_default_chart_type = function() { "epiviz.ui.charts.tree.Icicle" },
  get_measurements=function() {
    "Get all annotation info for all samples
    
    \\describe{
      \\item{chart_id_or_object}{An object of class \\code{EpivizChart} or an id for
        a chart loaded to the epiviz app.}
    }
    "
    samplesToRet <- colnames(.self$.leaf_sample_count_table)
    samplesToRet <- samplesToRet[-which(samplesToRet == "leaf")]
    samplesToRet <- samplesToRet[-which(samplesToRet == "otu_index")]
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
  
  row_to_dict=function(row){
    "Helper function to format each node entry for getHierarchy response
    
    \\describe{
      \\item{row}{Information for current node.}
      }
    "
    
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
    "Helper function to recursively build nested response for getHierarchy
    
    \\describe{
      \\item{root}{Root of subtree}
      \\item{df}{data.frame containing children to process}
    }
    "
    
    if(nrow(df) == 0) {
      root$children = NULL
      return(root)
    }
    
    children = df[which(df['parentId'] == as.character(unlist(root['id']))),]
    
    if(length(children) == 0){
      root$children = NULL
      return(root)
    }
    
    otherChildren = df[which(df['parentId'] != as.character(unlist(root['id']))),]
    
    children = children[order(children['order']),]
    
    if(nrow(children) > 0){
      for(row_index in seq_len(nrow(children))){
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
    
    graph_tree <- .self$.graph$.hierarchy_tree[,-which(colnames(.self$.graph$.hierarchy_tree) == "otu_index")]
    
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
    
    leaf_ordering_table <- as.data.table(.self$.graph$.hierarchy_tree[,c(.self$.graph$.feature_order[length(.self$.graph$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    
    lineage_DF <- as.data.frame(.self$.graph$.node_ids_table)
    lineage_table <- .self$.graph$.node_ids_table
    lineage_DF[,.self$.graph$.feature_order[1]] <- lineage_table[,get(.self$.graph$.feature_order[1])]
    
    for(i in seq(2,length(.self$.graph$.feature_order))){
        lineage_DF[,.self$.graph$.feature_order[i]] <- paste(lineage_DF[,.self$.graph$.feature_order[i-1]], lineage_table[,get(.self$.graph$.feature_order[i])], sep=",")
    }
    lineage_DT <- as.data.table(lineage_DF)
    
    lineage_DT_long <- melt(lineage_DT, id.vars = c("otu_index"),
                                measure.vars = c(colnames(lineage_DT)[1:(length(colnames(lineage_DT))-1)]),
                                variable.name = "taxonomy", variable.factor = FALSE)
    
    setnames(lineage_DT_long, c("otu_index", "taxonomy", "lineage"))
    
    temp_nodes_table <- merge(.self$.graph$.nodes_table, lineage_DT_long, by="lineage")
    temp_nodes_table <- temp_nodes_table[, otu_index:=as.integer(otu_index)]
      
      
    for(i in seq_len(num_rows)){
      if(as.integer(strsplit(nodesToRet[i], "-")[[1]][1]) == last_level_of_subtree){
        depths[i] = length(.self$.graph$.feature_order)
        level = length(.self$.graph$.feature_order)
        index <- which(.self$.graph$.node_ids_table[,level,with=FALSE] == nodeId)
        
        label <- as.character(unique(graph_tree[,level][index]))
        labels[i] <- label
        
        partition <- "NA"
        partitions[i] <- partition
        
        otu_index_temp <- leaf_ordering_table[leaf == nodeId, otu_index]
        starts[i] <- otu_index_temp
        leafIndexes[i] <- otu_index_temp
        ends[i] <- otu_index_temp
        
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
      
      if(nodesToRet[i] == "0-0"){
        leaf_indexes_temp <- temp_nodes_table[level == (length(.self$.graph$.feature_order)-1), otu_index,]
      } else{
        ids_match <- .self$.graph$.nodes_table[id == nodesToRet[i],]
        parents_match <- ids_match[parent == parentIds[i]]
        leaf_indexes_temp <- temp_nodes_table[lineage == parents_match[,lineage,], otu_index,]
      }
      #leaf_indexes_temp <- temp_nodes_table[lineage == .self$.graph$.nodes_table[id == nodesToRet[i],lineage,], otu_index,]
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
      
      taxonomy <- colnames(graph_tree)[level]
      taxonomys[i] <- taxonomy
      
      nchildren <- length(unique(.self$.graph$.node_ids_table[get(taxonomy)==nodesToRet[i],][[as.integer(level)+1]]))
      nchildrens[i] <- nchildren[1]
      
      nleaves_temp <- length(unname(unlist(unique(.self$.graph$.leaf_of_table[node_label==label, leaf]))))
      nleaves[i] <- nleaves_temp[1]
      
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
    
    if(!is.null(selectedLevels)){
      sanitize_selections <- list()
      for(n in names(.self$.nodeSelections)){
        n_level <- as.integer(unlist(strsplit(n, "-"))[1])
        selection_value <- as.integer(.self$.nodeSelections[[n]])
        if(!(selection_value == 2 && n_level <= .self$.levelSelected)){
          sanitize_selections[[n]] = .self$.nodeSelections[[n]]
        }
      }
      .self$.nodeSelections <- sanitize_selections
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

    leaf_ordering_table <- as.data.table(.self$.graph$.hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    leaf_ordering_table <- leaf_ordering_table[,leaf:=as.character(leaf)]
    leaf_ordering_table <- leaf_ordering_table[otu_index >= start & otu_index <= end]
    
    first_join <- merge(leaf_ordering_table, merge(nodes_at_level, .self$.graph$.leaf_of_table, by="id"), by="leaf")
    setorderv(first_join, "otu_index.x")
    
    nodes <- as.data.frame(unique(first_join[,c("node_label.x", "lineage.x")]))[,"node_label.x"]
    
    ends <- rep(0, length(nodes))
    starts <- rep(0, length(nodes))
    indexes <- rep(0, length(nodes))
    metadata <- list()
    metadata[['label']] <- list()
    metadata[['id']] <- list()
    metadata[['lineage']] <- list()

    for (i in seq_along(nodes)) {
      feature_label <- nodes[i]
      res <- first_join[node_label.x==feature_label,otu_index.x]
      if(length(res) > 0) {
        ends[i] <- max(res)
        starts[i] <- min(res)
        indexes[i] <- min(res)
        metadata[['label']][i] <- feature_label
        metadata[['id']][i] <- unique(.self$.graph$.nodes_table[node_label==nodes[i], id])[1]
        metadata[['lineage']][i] <- unique(.self$.graph$.nodes_table[node_label==nodes[i], lineage])[1]
      }
    }
    
    data_rows = list(end=ends, start=starts, index =indexes, metadata=metadata)
    return(data_rows)
  },
  
  getValues=function(measurements = NULL, start = 1, end = 1000, selectedLevels = 3, selections = NULL, row_order = NULL) { 
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
    
    leaf_sample_count_table_sub_select <- .self$.leaf_sample_count_table_long[otu_index >= start & otu_index <= end]

    first_join <- unique(merge(unique(.self$.graph$.leaf_of_table),unique(nodes_at_level), by="lineage")) 
    
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,-(which(colnames(leaf_sample_count_table_sub_select)=="otu_index")), with=FALSE]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[sample %in% measurements,]
    leaf_sample_count_table_sub_select <- leaf_sample_count_table_sub_select[,leaf:=as.character(leaf)]
    first_join <- first_join[,leaf:=as.character(leaf)]
    
    second_join <- merge(na.omit(unique(leaf_sample_count_table_sub_select)), na.omit(unique(first_join)), by="leaf", all.x = TRUE)

    results <- second_join[,sum(value), by=list(sample, node_label.x)]
    close_results <- as.data.frame(dcast(data = results, formula = node_label.x ~ sample, value.var='V1'))
    zero_results <- setdiff(nodes_at_level[,node_label], close_results[,"node_label.x"])
    if(length(zero_results) > 0){
      for(k in seq(1, length(zero_results))){
        new_index <- nrow(close_results)+1
        close_results[new_index,] <- 0.0
        close_results[new_index,]["node_label.x"] <- zero_results[k]
      }
    }
    if(is.null(row_order)){ 
      close_results <- close_results[order(close_results[,"node_label.x"]),] 
    } else { 
      rownames(close_results) <- close_results[,"node_label.x"] 
      close_results <- close_results[row_order,] 
    }
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
  
  getCombined=function(measurements = NULL, 
                          seqName, start = 1, end = 1000, 
                          order = NULL, nodeSelection = NULL, selectedLevels = NULL) {
    "Return the counts aggregated to selected nodes for the given samples
    
    \\describe{
    \\item{measurements}{Samples to get counts for}
    \\item{seqName}{name of datasource}
    \\item{start}{Start of feature range to query}
    \\item{end}{End of feature range to query}
    \\item{order}{Ordering of nodes}
    \\item{nodeSelection}{Node-id and selectionType pairs}
    \\item{selectedLevels}{Current aggregation level}
    
    }
    "
    
    # update node selections types to metaviztree
    if(!is.null(nodeSelection)) {
      for(n in names(nodeSelection)){
        .self$.nodeSelections[[n]] = nodeSelection[[n]]
      }
    }
    
    if(is.null(selectedLevels)) {
      selectedLevels = .self$.levelSelected
    }
    selections = .self$.nodeSelections
    measurements = unique(measurements)

    data_rows = getRows(measurements = measurements, start = start, end = end, selectedLevels = selectedLevels, selections = selections) 
    row_order = unlist(data_rows$metadata$label) 
    data_columns = getValues(measurements = measurements, start = start, end = end, selectedLevels = selectedLevels, selections = selections, row_order = row_order)
    
    result <- list(
      cols = data_columns,
      rows = data_rows,
      globalStartIndex = data_rows$start[[1]]
    )

    return(result)
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
  }, 
  
  updateSplineAlpha=function(alpha=NULL){ 
    " Unimplemented in EpivizMetagenomicsData 
      Implemented in EpivizMetagenomicsDataTimeSeries 
     
    \\describe{ 
    \\item{alpha}{Smoothing Spline parameter} 
    } 
    " 
  },
  
  getSpline = function(){ 
    " Unimplemented in EpivizMetagenomicsData 
      Implemented in EpivizMetagenomicsDataTimeSeries 
    " 
  } 
)


EpivizMetagenomicsData$methods(
  toNEO4JDbHTTP =function(batch_url, neo4juser, neo4jpass, datasource, description=NULL) {
    ' 
    Write an `EpivizMetagenomicsData` object to a Neo4j graph database
    
    @param batch_url (character) Neo4j database url and port for processing batch http requests
    @param neo4juser (character) Neo4j database user name
    @param neo4jpass (character) Neo4j database password
    @param datasource (character) Name of Neo4j datasource node for this `EpivizMetagenomicsData` object 
    
    @examples
    library(metagenomeSeq)
    data("mouseData")
    mobj <- metavizr:::EpivizMetagenomicsData$new(object=mouseData)
    mobj$toNEO4JDbHTTP(batch_url = "http://localhost:7474/db/data/batch", neo4juser = "neo4juser", neo4jpass = "neo4jpass", datasource = "mouse_data")
    '
    
    cat("Saving sample data...")
    .saveSampleDataNEO4JHTTP(batch_url, neo4juser, neo4jpass)
    cat("Done\n")
    
    cat("Saving datasource data...")
    .saveDataSourceNEO4JHTTP(batch_url, neo4juser, neo4jpass, datasource, description)
    cat("Done\n")
    
    cat("Saving hierarchy...")
    .saveHierarchyNEO4JHTTP(batch_url, neo4juser, neo4jpass, datasource)
    cat("Done\n")
    
    cat("Saving properties...")
    .neo4jUpdatePropertiesHTTP(batch_url, neo4juser, neo4jpass, create_id_index = TRUE)
    cat("Done\n")
    
    cat("Saving Data Matrix...")
    .saveMatrixNEO4JHTTP(batch_url, neo4juser, neo4jpass, datasource)
    cat("Done\n")
    
    cat("Saving properties...")
    .neo4jUpdatePropertiesHTTP(batch_url, neo4juser, neo4jpass, create_prop_index = TRUE)
    cat("Done\n")
  },
  
  .buildBatchJSON = function(query_in, param_list, id=0, id_last=TRUE, full_query=TRUE,json_query_in=NULL, params_complete=FALSE){
    json_start <- "["
    method <- "{\"method\" : \"POST\","
    to <- "\"to\" : \"cypher\","
    body_start <- "\"body\" : {"
    params_start <- "\"params\" : {"
    params_end <- "}"
    body_end <- "},"
    
    if(id_last){
      id <- paste("\"id\": ", as.character(id), "}", sep="")
    }
    else{
      id <- paste("\"id\": ", as.character(id), "},", sep="")
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
      .self$.json_query <- query_final
      return(query_final)
    }
    
    query = paste("\"query\" : \"", query_in, "\",", sep = "")
    for(k in seq_along(param_list)){
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
    if(!full_query){
      return(json_query)
    }
    query_final <- paste0(json_start, json_query, json_end)
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
    for (j in seq_len(nrow(sampleAnnotationToNeo4j))){
      row <- sampleAnnotationToNeo4j[j,]
      props <- ""
      for (i in seq_len(length(keys)-1)){
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
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
  },
  
  .saveDataSourceNEO4JHTTP =function(batch_url, neo4juser, neo4jpass, datasource, file=NULL, description=NULL) {
    description_temp <- as.character(description)
    if(is.null(description)){
        description_temp = as.character(datasource)
    }
    query <- "CREATE (:Datasource {label: {label_param}, description: {description_param}})"
    params <- list(label_param=paste("\"", as.character(datasource), "\"", sep=""), description_param=paste("\"", as.character(description_temp), "\"", sep=""))
      
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
    
    dfToNeo4j <- as.data.frame(.self$.graph$.nodes_table)
    dfToNeo4j <- dfToNeo4j[order(dfToNeo4j[,"id"]),]
    names(dfToNeo4j)[names(dfToNeo4j) == "node_label"] <- "label"
    names(dfToNeo4j)[names(dfToNeo4j) == "parent"] <- "parentId"
    dfToNeo4j[,"lineageLabel"] <- dfToNeo4j[,"lineage"]
    
    dfToNeo4j[,"depth"] <- (as.integer(dfToNeo4j[,"level"]))
    dfToNeo4j[,"taxonomy"] <- .self$.feature_order[(as.integer(dfToNeo4j[,"level"]) + 1)]
        
    leaf_ordering_table <- as.data.table(.self$.graph$.hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
        
    leaf_table_lowercase <- .self$.graph$.leaf_of_table
    leaf_table_lowercase <- leaf_table_lowercase[,node_label:=tolower(node_label)]
        
    starts <- rep(1,nrow(dfToNeo4j))
    ends <- rep(1,nrow(dfToNeo4j))
    nleaves <- rep(1, nrow(dfToNeo4j))
    nchildrens <- rep(1, nrow(dfToNeo4j))
    lineageLabels <- rep("", nrow(dfToNeo4j))
    
    temp_nodes_table <- merge(.self$.graph$.nodes_table, .self$.graph$.leaf_of_table, by="lineage")
    temp_nodes_table <- temp_nodes_table[, otu_index:=as.integer(otu_index)]
        
    for(i in 1:nrow(dfToNeo4j)){
      node <- dfToNeo4j[,"label"][i]
            
      leaf_indexes_temp <- temp_nodes_table[lineage == dfToNeo4j[,"lineage"][i], otu_index,]
      nleaves[i] <- length(leaf_indexes_temp)
      
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
      nchildrens[i] <- nrow(.self$.graph$.node_ids_table[get(.self$.graph_feature_order[(as.integer(dfToNeo4j[,"level"][i])+1)])==dfToNeo4j[,"id"][i],])
      lineage = .self$.graph$.nodes_table[get("id")==dfToNeo4j[,"id"][i],get("lineage")][[1]]
            
      lineageLabel <- sapply(strsplit(lineage, ",")[[1]], function(str_id) {
        .self$.graph$.nodes_table[get("id") == str_id, get("node_label")][[1]]
      })
            
      lineageLabels[i] <- paste(lineageLabel, collapse=",")
    }
        
    dfToNeo4j$start <- as.integer(starts)
    dfToNeo4j$leafIndex <- as.integer(starts)
    dfToNeo4j$end <- as.integer(ends)
    dfToNeo4j$nchildren <- nchildrens
    dfToNeo4j$nleaves <- nleaves
    dfToNeo4j$lineageLabel <- lineageLabels
    dfToNeo4j$partition = NA
    
    json_start <- "["
    json_end <- "]"
    datasource_param_key <- "datasource"
    datasource_param_value <- as.character(datasource)
    
    keys = colnames(dfToNeo4j)
    query <- "CREATE (:Feature {props})"
    json_query <- ""
    id_counter <- 0
    
    for (j in 1:nrow(dfToNeo4j)){
      row <- dfToNeo4j[j,]
      props <- "" 
      for (i in seq_along(keys)){
        if  (typeof(keys[i]) == "numeric")
          props <- paste(props, "\"", keys[i], "\"", " : ", "\"", gsub("'", "",row[, keys[i]]), "\"", ", ",sep="")
        else
          props <- paste(props, "\"", keys[i], "\"", " : \"", gsub("'", "",row[, keys[i]]), "\", ", sep="")
      }
      i = length(keys)
      
      props <- paste(props, "\"", datasource_param_key, "\"", " : \"", datasource_param_value, "\"", sep="")
      params <- list("props"=props)
      
      writeBatchOut <- j%%5000 == 0
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
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    
    query <- "CREATE INDEX ON :Feature (id)"
    query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
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
    query <- "MATCH (fParent:Feature {id : {parent_id}, datasource: {datasource_param}}) WHERE {lineage_param} CONTAINS fParent.lineage MATCH (f:Feature {id : {child_id}, datasource: {datasource_param}, lineage: {lineage_param}, order: {order_param}}) CREATE (fParent)-[:PARENT_OF]->(f)"    
    json_end <- "]"
    json_query <- ""
    
    id_counter <- 0
    
    for (j in 1:nrow(dfToNeo4j)){
      row <- dfToNeo4j[j,]
      parentid <- paste("\"", row$parentId, "\"", sep="")
      if(parentid == "None"){
        continue
      }
      childid <- paste("\"", row$id, "\"", sep="")
      lineage <- paste("\"", row$lineage, "\"", sep="")
      order <- paste("\"", row$order, "\"", sep="")
            
      writeBatchOut <- j%%5000 == 0
      if (j <= nrow(dfToNeo4j)-1){
        if (writeBatchOut){
          params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), parent_id=parentid, child_id=childid, lineage_param = lineage, order_param = order)
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = TRUE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
        }
        else {
          params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), parent_id=parentid, child_id=childid, lineage_param = lineage, order_param = order)
          json_query <- .buildBatchJSON(query_in = query, param_list = params, id= id_counter, id_last = FALSE, full_query = FALSE, json_query_in = json_query, params_complete = TRUE)
        }
      }
      else {
        params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), parent_id=parentid, child_id=childid, lineage_param = lineage, order_param = order)
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
    
    query <- "MATCH (fNode:Feature {datasource: {datasource_param}})-[:PARENT_OF*]->(fLeaf:Feature {depth: {depth_param}, datasource: {datasource_param} }) CREATE (fNode)-[:LEAF_OF]->(fLeaf)"
    params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), depth_param=paste("\"", as.character(length(.self$.feature_order)), "\"", sep=""))
    query_final <- .buildBatchJSON(query_in = query, param_list = params, params_complete = TRUE)
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
    else {
      write(query, file=file, append = TRUE)
    }
    
    query <- "MATCH (fLeaf:Feature {depth: {depth_param}, datasource: {datasource_param}}) CREATE (fLeaf)-[:LEAF_OF]->(fLeaf)"
    params <- list(datasource_param=paste("\"", as.character(datasource), "\"", sep=""), depth_param=paste("\"", as.character(length(.self$.feature_order)), "\"", sep=""))
    query_final <- .buildBatchJSON(query_in = query, param_list = params, params_complete = TRUE)
    
    if(!is.null(batch_url)) {
      r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
      stop_for_status(r)
    }
  },
  
  .saveMatrixNEO4JHTTP = function(batch_url, neo4juser, neo4jpass,  datasource, file=NULL) {
    leaf_sample_count_table_temp <- .self$.leaf_sample_count_table_long[,-(which(colnames(.self$.leaf_sample_count_table_long) == "otu_index")), with=FALSE]
    valuesToNeo4j = as.data.frame(leaf_sample_count_table_temp)
    
    json_start <- "["
    query <- "MATCH (f:Feature {label : {node_id}, datasource: {datasource_param}}) MATCH (s:Sample {id: {sample_id}}) CREATE (s)-[:COUNT {val: {count_param}}]->(f)"
    
    json_end <- "]"
    json_query <- ""
    
    cypherCount = 0
    id_counter <- 0
    
    for (j in 1:nrow(valuesToNeo4j)){
      row <- valuesToNeo4j[j,]
      nodeid <- paste("\"", row$leaf, "\"", sep="")
      sampleid <- paste("\"", row$sample, "\"", sep="")
      count <- paste("\"", as.character(row$value), "\"", sep="")
      writeBatchOut <- j%%5000 == 0
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
    
  },
  
  .neo4jUpdatePropertiesHTTP = function(batch_url, neo4juser, neo4jpass, create_id_index = FALSE, create_prop_index = FALSE) {
    
    if(create_id_index){
      query <- "CREATE INDEX ON :Sample (id)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL) 
      
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
    }
    
    if(create_prop_index){
      query <- "MATCH (f:Feature) SET f.depth = toInt(f.depth)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
      
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
      
      query <- "MATCH ()-[c:COUNT]->() SET c.val = toFloat(c.val)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
      
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
      
      query <- "MATCH (f:Feature) SET f.start = toInt(f.start)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
          
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
            
      query <- "MATCH (f:Feature) SET f.end = toInt(f.end)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
            
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
            
      query <- "MATCH (f:Feature) SET f.nleaves = toInt(f.nleaves)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
            
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
          
      query <- "MATCH (f:Feature) SET f.nchildren = toInt(f.nchildren)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
            
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
            
      query <- "MATCH (f:Feature) SET f.leafIndex = toInt(f.leafIndex)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
            
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
      
      query <- "CREATE INDEX ON :Feature (depth)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
      
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
      
      query <- "CREATE INDEX ON :Feature (start)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
      
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
      
      query <- "CREATE INDEX ON :Feature (end)"
      query_final <- .buildBatchJSON(query_in = query, param_list = NULL)
      
      if(!is.null(batch_url)) {
        r <- POST(batch_url, body = query_final, encode = "json", authenticate(user = neo4juser, password = neo4jpass))
        stop_for_status(r)
      }
    }
  }
)
