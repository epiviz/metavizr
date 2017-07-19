#' Graph implementation to query hierarchical feature data
#'
#' Used to manage aggregation and range queries from the Metaviz app UI. 
#'  
MetavizGraph <- setRefClass("MetavizGraph",
  fields=list(
    .feature_order = "ANY",
    .leaf_of_table = "ANY",
    .hierarchy_tree = "ANY",
    .node_ids_table = "ANY",
    .nodes_table = "ANY"
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

      message("creating leaf_of_table")
      .self$.leaf_of_table <- .create_leaf_of_table()

      .self$.leaf_of_table <- merge(unique(.self$.nodes_table[,mget(c("lineage", "id"))]),
                                    unique(.self$.leaf_of_table) , by="lineage")
      .self$.leaf_of_table <- .self$.leaf_of_table[,id:=as.character(id)]
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
      lineage_DF[,feature_order[1]] <- lineage_table[,get(feature_order[1])]
      
      for(i in seq(2,length(feature_order))){
        lineage_DF[,feature_order[i]] <- paste(lineage_DF[,feature_order[i-1]], lineage_table[,get(feature_order[i])], sep=",")
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
      
      ret_table <- ret_table[order(parent)]
      parent_list <- ret_table[,parent]
      orders <- rep(1, length(parent_list))
      
      for(j in seq(2, length(parent_list))){
        if(parent_list[j] == parent_list[j-1]){
          orders[j] = orders[j-1]+1
        }
      }
      ret_table[,order:=orders]
      return(ret_table)
    },
    
    .create_leaf_of_table=function(){
      "Create a data.table with leaf, ancestor relationship for each leaf
    
      \\describe{
          }
      \\value{ret_table}{data.table leaf of relationship for each node}
      "                               
      #            
      feature_order <- .self$.feature_order
      
      temp_hiearchy_DT <- as.data.table(.self$.hierarchy_tree)
      num_features <- length(feature_order)
      hiearchy_cols <- colnames(.self$.hierarchy_tree)
                              
      melt_res <- melt(temp_hiearchy_DT, id.vars = c(feature_order[num_features], "otu_index"), 
                      measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols)-1)]))
      label_table <- melt_res[,c(1,2,4)]
      setnames(label_table, c("leaf", "otu_index","node_label"))
                              
      label_table <- label_table[,leaf:=as.character(leaf)]
      label_table <- label_table[,otu_index:=as.character(otu_index)]
      
      lineage_DF <- as.data.frame(.self$.node_ids_table)
      lineage_table <- .self$.node_ids_table
      lineage_DF[,feature_order[1]] <- lineage_table[,get(feature_order[1])]
      
      for(i in seq(2,length(feature_order))){
        lineage_DF[,feature_order[i]] <- paste(lineage_DF[,feature_order[i-1]], lineage_table[,get(feature_order[i])], sep=",")
      }
      lineage_DT <- as.data.table(lineage_DF)
      
      melt_res_lineage <- melt(lineage_DT, id.vars = c(feature_order[num_features], "otu_index"), measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols))-1]))

      lineage_leaf_of_table <- unique(melt_res_lineage[,c(2,4)])
      setnames(lineage_leaf_of_table, c("otu_index","lineage"))

      lineage_leaf_of_table <- lineage_leaf_of_table[,otu_index:=as.character(otu_index)]
      
      lineage_df <- as.data.frame(lineage_leaf_of_table)
      leaf_node_label <- as.data.frame(label_table)[,c("leaf", "node_label")]
      
      ret_table <- as.data.table(cbind(lineage_df, leaf_node_label))
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
      
      fd = fData(obj_in) 
      for( i in seq(ncol(fd))){
        fd[,i] = as.character(fd[,i]) 
      } 
      fData(obj_in) = fd
                               
      replacing_na_obj_fData <- fData(obj_in)[,feature_order]
      for(i in seq(1, length(feature_order))){
        na_indices <- which(is.na(replacing_na_obj_fData[,feature_order[i]]))
        for(j in seq(1, length(na_indices))){
          if(i > 1) {
            replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][na_indices[j]], sep="_")
          } else {
            replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
          }
        }
        na_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NA")
        for(j in seq(1, length(na_indices))){
          if(i > 1){ 
            replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][na_indices[j]], sep="_")
          } else{
            replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
          }
        }
        null_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NULL")
        for(j in seq(1, length(null_indices))){
          if(i > 1){ 
            replacing_na_obj_fData[,feature_order[i]][null_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][null_indices[j]], sep="_")
          } else{
            replacing_na_obj_fData[,feature_order[i]][null_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
          }
        }
      }
                       
      obj_fData <- as.data.table(replacing_na_obj_fData)
      cols <- feature_order[1:length(feature_order)-1]
      order <- rep(1, length(feature_order)-1)
      ordered_fData <- setorderv(obj_fData, cols = cols, order = order)
                               
      otu_indexes <- seq(1:length(ordered_fData[,get(feature_order[length(feature_order)])]))
      ordered_fData <- ordered_fData[, otu_index:=otu_indexes]
      ordered_fData_df <- as.data.frame(ordered_fData)
      
      if(length(unique(ordered_fData_df[,1])) > 1){
        rootFeature <- rep("AllFeatures", nrow(ordered_fData_df))
        ordered_fData_df <- cbind(rootFeature, ordered_fData_df)
        .self$.feature_order <- unlist(c("rootFeature", feature_order))
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
        temp_level <- data.table(table_node_ids[, c(level, "otu_index")])
        temp_level_count <- temp_level[, .(leaf_index = .I[which.min(otu_index)], count = .N), by=eval(level)]
        
        level_features <- as.character(table_node_ids[[level]])
        for(i in seq_len(nrow(temp_level_count))) {
          row <- temp_level_count[i,]
          if(depth==1 && i == 1){
            id <- paste(depth-1, 0, sep="-")
          } else{
            id <- paste(depth-1, paste(digest(row[,1]), i, sep=""), sep="-")
          }
          level_features <- replace(level_features, which(level_features == row[[level]]), id)
        }
        level_features
      })
      
      node_ids_dt <- as.data.table(id_list)
      node_ids_dt$otu_index <- as.character(table_node_ids$otu_index)
      
      return(node_ids_dt)
    }
  )
)

#' Build a MetavizTree object from another object
#' 
#' @param object The object from which taxonomy data is extracted
#' @param ... Additional arguments
#' @return a \code{\link{MetavizGraph}} object
setGeneric("buildMetavizGraph", signature=c("object"),
           function(object, ...) standardGeneric("buildMetavizGraph"))

#' @describeIn buildMetavizGraph Build graph from a \code{\link[metagenomeSeq]{MRexperiment-class}} object
#' @importFrom Biobase fData
#' @param feature_order Ordering of leaves (features) in taxonomy tree
setMethod("buildMetavizGraph", "MRexperiment", function(object, feature_order, ...) {
  MetavizGraph$new(object, feature_order = feature_order)
})
