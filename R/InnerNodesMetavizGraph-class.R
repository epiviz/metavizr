#' Graph implementation to query hierarchical feature data
#'
#' Used to manage aggregation and range queries from the Metaviz app UI. 
#'  
InnerNodesMetavizGraph <- setRefClass("InnerNodesMetavizGraph",
  fields=list(
    .feature_order = "ANY",
    .hierarchy_tree = "ANY",
    .node_ids_table = "ANY",
    .nodes_table = "ANY",
    .hierarchy_tree_order = "ANY"
  ),
                            
  methods=list(
    initialize=function(object, feature_order=NULL) {
      .self$.feature_order = feature_order
                              
      message("creating hierarchy_tree")
      .self$.hierarchy_tree <- .create_hierarchy_tree(object)
                                
      message("creating node_ids_table")
      .self$.node_ids_table <- .create_node_ids()
                                
      message("creating nodes_table")
      .self$.nodes_table <- .create_nodes_table()
    },
                              
    .create_nodes_table=function(){
      "Create a data.table with information for each node in the feature hierarchy
                                
      \\describe{
      }
      \\value{ret_table}{data.table containing information for each node}
      "
                                
      feature_order <- .self$.feature_order
                                
      lineage_DF <- as.data.frame(.self$.node_ids_table)
      lineage_table <- .self$.node_ids_table
      lineage_table_df <- as.data.frame(lineage_table)
      lineage_DF[,feature_order[1]] <- lineage_table[,get(feature_order[1])]
                                
      for(i in seq(2,length(feature_order))){
        for(j in seq(1, nrow(lineage_DF))){
          if(!is.na(lineage_DF[j,feature_order[i]])){
            lineage_DF[j,feature_order[i]] <- paste(lineage_DF[j,feature_order[i-1]], lineage_table_df[j,feature_order[i]], sep=",")
          }
        }
      }
      lineage_DT <- as.data.table(lineage_DF)
                                
      root_parents <- rep("None", length(.self$.node_ids_table[,get(feature_order[1])]))
                            nodes_tab <- data.frame(id = .self$.node_ids_table[,get(feature_order[1])], parent = root_parents, 
                            lineage = .self$.node_ids_table[,get(feature_order[1])], 
                            node_label = .self$.hierarchy_tree[,1], level = rep(0, length(.self$.hierarchy_tree[,1])))
                                
      for(i in seq(2, length(feature_order))){
        temp_nodes_tab <- data.frame(id = .self$.node_ids_table[,get(feature_order[i])], 
                                       parent = .self$.node_ids_table[,get(.self$.feature_order[i-1])], 
                                       lineage = lineage_DT[,get(feature_order[i])],  node_label = .self$.hierarchy_tree[,i],  
                                       level = rep(i-1, length(.self$.hierarchy_tree[,i])))
                                  
        nodes_tab <- rbind(nodes_tab[rownames(unique(nodes_tab[,c("id","parent")])),], temp_nodes_tab[rownames(unique(temp_nodes_tab[,c("id","parent")])),])
      }
                                
      ret_table <- as.data.table(nodes_tab)
      ret_table <- ret_table[,id:=as.character(id)]
      ret_table <- ret_table[,parent:=as.character(parent)]
      ret_table <- ret_table[,lineage:=as.character(lineage)]
      ret_table <- ret_table[,node_label:=as.character(node_label)]
      ret_table <- ret_table[,level:=as.integer(level)]
      ret_table <- ret_table[node_label != "Not_Annotated",]
      ret_table <- ret_table[order(parent)]
      parent_list <- ret_table[,parent]
      orders <- rep(1, length(parent_list))
                                
      for(j in seq(2, length(parent_list))){
        if(parent_list[j] == parent_list[j-1]){
          orders[j] = orders[j-1]+1
        }
      }
      
      starts <- rep(0, nrow(ret_table))
      ends <- rep(0, nrow(ret_table))
      starts[nrow(ret_table)] <- 1
      ends[nrow(ret_table)] <- nrow(ret_table)-1
                                
      for(k in seq(2, length(feature_order)-1)){
        features_at_level <- unique(lineage_DT[,get(feature_order[k])])
        features_at_level <- features_at_level[which(!is.na(features_at_level))]
        for(f in features_at_level){
          subset_lineage_DT <- lineage_DT[get(feature_order[k]) == f,]
          starts[which(ret_table[,lineage] == f)] <- as.integer(subset_lineage_DT[is.na(get(feature_order[k+1])),start])
          ends[which(ret_table[,lineage] == f)] <- as.integer(subset_lineage_DT[is.na(get(feature_order[k+1])),end])
        }
      }
      
      k = length(feature_order)
      features_at_level <- unique(lineage_DT[,get(feature_order[k])])
      features_at_level <- features_at_level[which(!is.na(features_at_level))]
      
      for(f in features_at_level){
        starts[which(ret_table[,lineage] == f)] <- as.integer(lineage_DT[get(feature_order[k]) == f,start])
        ends[which(ret_table[,lineage] == f)] <- as.integer(lineage_DT[get(feature_order[k]) == f,end])
      }
                                
      ret_table <- ret_table[,start:=starts]
      ret_table <- ret_table[,end:=ends]
                                
      ret_table[,order:=orders]
      return(ret_table)
      },
                              
      .create_hierarchy_tree=function(obj_in){
        "Create a data.frame with the hierarchy ordered by each level of the hierarchy
                                
          \\describe{
            \\item{obj_in}{An MRexperiment object containing featureData}
            }
          \\value{ordered_fData}{data.frame with sorted hierarchy of features}
        "  
                                
        feature_order <- .self$.feature_order
                              
        level_annotated <- rep(0, nrow(fData(obj_in)))
        f_data <- fData(obj_in)
        for(i in seq(1, nrow(f_data))){
          if (length(which(is.na(f_data[i,]))) == 0){
            level_annotated[i] <- length(colnames(f_data))
          } else {
            level_annotated[i] <- min(which(is.na(f_data[i,])))-2
          }
        }
                                
        fd = fData(obj_in) 
        for( i in seq(ncol(fd))){
          fd[,i] = as.character(fd[,i]) 
        } 
        fData(obj_in) = fd
                                
        replacing_na_obj_fData <- fData(obj_in)[,feature_order]
                                
        for(i in seq(1, length(feature_order))){
          empty_indices <- which(is.na(replacing_na_obj_fData[,feature_order[i]]))
          replacing_na_obj_fData[,feature_order[i]][empty_indices] <- "Not_Annotated"
          
          na_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NA")
          replacing_na_obj_fData[,feature_order[i]][na_indices] <- "Not_Annotated"

          null_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NULL")
          replacing_na_obj_fData[,feature_order[i]][null_indices] <- "Not_Annotated"
        }
                                
        replacing_na_obj_fData[,"level_annotated"] <- level_annotated
        replacing_na_obj_fData[,"previous_order"] <- seq(1, nrow(replacing_na_obj_fData))
        
        obj_fData <- as.data.table(replacing_na_obj_fData)
        cols <- feature_order[1:length(feature_order)-1]
        order <- rep(1, length(feature_order)-1)
        ordered_fData <- setorderv(obj_fData, cols = cols, order = order)
        .self$.hierarchy_tree_order <- ordered_fData[,previous_order]
        ordered_fData <- ordered_fData[,!"previous_order"]
                                
        indexes <- seq(1:length(ordered_fData[,get(feature_order[length(feature_order)])]))
        ordered_fData <- ordered_fData[, index:=indexes]
                                
        ordered_fData <- as.data.frame(ordered_fData)
        starts <- rep(0, nrow(ordered_fData))
        ends <- rep(0, nrow(ordered_fData))
        
        for(i in seq(1, length(feature_order)-1)){
          inner_nodes_at_level <- which(ordered_fData[,feature_order[i+1]] == "Not_Annotated")
          for(n in inner_nodes_at_level){
            if(ordered_fData[,feature_order[i]][n] != "Not_Annotated"){
              insert_index <- n
              indices_for_inner_node <- ordered_fData[which(ordered_fData[,feature_order[i]] == ordered_fData[,feature_order[i]][n]),][,"index"]
              starts[insert_index] <- sort(indices_for_inner_node)[1]
              ends[insert_index] <- sort(indices_for_inner_node)[length(indices_for_inner_node)]
            }
          }
        }
                                
        inner_nodes_at_level <- which(ordered_fData[,feature_order[length(feature_order)]] != "Not_Annotated")

        for(j in seq(1, length(inner_nodes_at_level))){
          insert_index <- inner_nodes_at_level[j]
          starts[insert_index] <- inner_nodes_at_level[j]
          ends[insert_index] <- inner_nodes_at_level[j]
        }
                                
        ordered_fData <- as.data.table(ordered_fData)
        ordered_fData <- ordered_fData[,-"index"]
        ordered_fData <- ordered_fData[, start:=starts]
        ordered_fData <- ordered_fData[, end:=ends]
                                
        ordered_fData_df <- as.data.frame(ordered_fData)
                                
        if(length(unique(ordered_fData_df[,1])) > 1){
          allFeatures <- rep("AllFeatures", nrow(ordered_fData_df))
          ordered_fData_df <- cbind(allFeatures, ordered_fData_df)
          .self$.feature_order <- unlist(c("allFeatures", feature_order))
        }
                                
        return(ordered_fData_df)
      },
                              
      .create_node_ids=function(){
        "Create a data.table with unique ids used for metavizr to specify level and child for any node 
                                
          \\describe{
            \\item{feature_order}{The order of hierarchy as colnames of fData for the MRexperiment object}
            }
          \\value{table_node_ids}{data.table with node ids in sorted hierarchy}
        "  
          
        feature_order <- .self$.feature_order
                                
        table_node_ids <- .self$.hierarchy_tree
        id_list <- sapply(feature_order, function(level) {
          depth <- which(feature_order == level)
          temp_level <- data.table(table_node_ids[, c(level)])

          level_features <- as.character(table_node_ids[[level]])
          for(i in seq_len(nrow(temp_level))) {
            row <- temp_level[i,]
            if(depth==1 && i == 1){
              id <- paste(depth-1, 0, sep="-")
            } else if (row[,1] == "Not_Annotated") {
              id <- NA
            } else{
              id <- paste(depth-1, paste(digest(row[,1]), i, sep=""), sep="-")
            }
            level_features <- replace(level_features, which(level_features == row[[1]]), id)
          }
          level_features
        })
                                
        node_ids_dt <- as.data.table(id_list)
        node_ids_dt$start <- as.character(table_node_ids$start)
        node_ids_dt$end <- as.character(table_node_ids$end)
                              
        return(node_ids_dt)
      }
  )
)

#' Build a MetavizTree object from another object
#' 
#' @param object The object from which taxonomy data is extracted
#' @param ... Additional arguments
#' @return a \code{\link{InnerNodesMetavizGraph}} object
setGeneric("buildInnerNodesMetavizGraph", signature=c("object"),
           function(object, ...) standardGeneric("buildInnerNodesMetavizGraph"))

#' @describeIn buildInnerNodesMetavizGraph Build graph from a \code{\link[metagenomeSeq]{MRexperiment-class}} object
#' @importFrom Biobase fData
#' @param feature_order Ordering of leaves (features) in taxonomy tree
setMethod("buildInnerNodesMetavizGraph", "MRexperiment", function(object, feature_order, ...) {
  InnerNodesMetavizGraph$new(object, feature_order = feature_order)
})
