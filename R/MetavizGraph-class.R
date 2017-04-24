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
      .self$.hierarchy_tree <- .create_hierarchy_tree(object, feature_order = feature_order)
                               
      message("creating leaf_of_table")
      .self$.leaf_of_table <- .create_leaf_of_table(feature_order = feature_order)
                           
      message("creating node_ids_table")
      .self$.node_ids_table <- .create_node_ids(feature_order = feature_order)
                               
      message("creating nodes_table")
      .self$.nodes_table <- .create_nodes_table(feature_order = feature_order)
                               
                               
      .self$.leaf_of_table <- merge(unique(.self$.nodes_table[,mget(c("id", "node_label"))]), unique(.self$.leaf_of_table) , by="node_label")
      .self$.leaf_of_table <- .self$.leaf_of_table[,id:=as.character(id)]
    },
                             
    .create_nodes_table=function(feature_order){
      "Create a data.table with information for each node in the feature hierarchy
    
      \\describe{
      \\item{feature_order}{The order of hierarchy as colnames of fData for the MRexperiment object}
        }
      \\value{ret_table}{data.table containing information for each node}
      "
                               
      lineage_DF <- as.data.frame(.self$.node_ids_table)
      lineage_table <- .self$.node_ids_table
      lineage_DF[,feature_order[1]] <- lineage_table[,get(feature_order[1])]
                               
      for(i in seq(2,length(feature_order))){
        lineage_DF[,feature_order[i]] <- paste(lineage_DF[,feature_order[i-1]], lineage_table[,get(feature_order[i])], sep=",")
      }
      lineage_DT <- as.data.table(lineage_DF)
                               
                               
      nodes_tab_root <- .self$.node_ids_table[,get(feature_order[1])]
      root_parents <- rep("None", length(.self$.node_ids_table[,get(feature_order[1])]))
      nodes_tab_root <- data.frame(id = nodes_tab_root, parent = root_parents, lineage = .self$.node_ids_table[,get(feature_order[1])], 
                                   node_label = .self$.hierarchy_tree[,1], level = rep(1, length(.self$.hierarchy_tree[,1])))
                               
      nodes_tab <- data.frame(id = .self$.node_ids_table[,get(.self$.feature_order[2])], 
                              parent = .self$.node_ids_table[,get(.self$.feature_order[1])], 
                              lineage = lineage_DT[,get(feature_order[2])], node_label = .self$.hierarchy_tree[,2], 
                              level = rep(2, length(.self$.hierarchy_tree[,2])))
                               
      nodes_tab <- rbind(unique(nodes_tab_root), unique(nodes_tab))
                               
      for(i in seq(3, length(feature_order))){
        temp_nodes_tab <- data.frame(id = .self$.node_ids_table[,get(.self$.feature_order[i])], 
                                     parent = .self$.node_ids_table[,get(.self$.feature_order[i-1])], 
                                     lineage = lineage_DT[,get(feature_order[i])],  node_label = .self$.hierarchy_tree[,i],  
                                     level = rep(i, length(.self$.hierarchy_tree[,i])))
                                 
        nodes_tab <- rbind(unique(nodes_tab), unique(temp_nodes_tab))
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
                             
    .create_leaf_of_table=function(feature_order){
      "Create a data.table with leaf, ancestor relationship for each leaf
    
      \\describe{
        \\item{feature_order}{The order of hierarchy as colnames of fData for the MRexperiment object}
          }
      \\value{ret_table}{data.table leaf of relationship for each node}
      "                               
                               
      temp_hiearchy_DT <- as.data.table(.self$.hierarchy_tree)
      num_features <- length(feature_order)
      hiearchy_cols <- colnames(.self$.hierarchy_tree)
                               
      melt_res <- melt(temp_hiearchy_DT, id.vars = c(.self$.feature_order[num_features]), 
                       measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols)-1)]))
      ret_table <- melt_res[,c(1,3)]
      setnames(ret_table, c("leaf", "node_label"))
                               
      ret_table <- ret_table[,leaf:=as.character(leaf)]
      return(ret_table)
    },
                             
    .create_hierarchy_tree=function(obj_in, feature_order){
      "Create a data.frame with the hierarchy ordered by each level of the hierarchy
    
      \\describe{
        \\item{obj_in}{An MRexperiment object containing featureData}
        \\item{feature_order}{The order of hierarchy as colnames of fData for the MRexperiment object}
          }
      \\value{ordered_fData}{data.frame with sorted hierarchy of features}
      "  
      fd = fData(obj_in) 
      for( i in seq(ncol(fd))){
        fd[,i] = as.character(fd[,i]) 
      } 
      fData(obj_in) = fd
                               
      replacing_na_obj_fData <- fData(obj_in)[,feature_order]
      for(i in seq(1, length(feature_order))){
        na_indices <- which(is.na(replacing_na_obj_fData[,i]))
        replacing_na_obj_fData[,i][na_indices] <- paste("Not_Annotated", feature_order[i], sep="_")
        na_indices <- which(replacing_na_obj_fData[,i] == "NA")
        replacing_na_obj_fData[,i][na_indices] <- paste("Not_Annotated", feature_order[i], sep="_")
      }
                               
      obj_fData <- as.data.table(replacing_na_obj_fData)
      cols <- feature_order[1:length(feature_order)-1]
      order <- rep(1, length(feature_order)-1)
      ordered_fData <- setorderv(obj_fData, cols = cols, order = order)
                               
      otu_indexes <- seq(1:length(ordered_fData[,get(feature_order[length(feature_order)])]))
      ordered_fData <- ordered_fData[, otu_index:=otu_indexes]
                               
      return(as.data.frame(ordered_fData))
    },
                             
    .create_node_ids=function(feature_order){
      "Create a data.table with unique ids used for metavizr to specify level and child for any node 
    
       \\describe{
         \\item{feature_order}{The order of hierarchy as colnames of fData for the MRexperiment object}
           }
       \\value{table_node_ids}{data.table with node ids in sorted hierarchy}
      "  
      table_node_ids <- .self$.hierarchy_tree
      leaf_ordering_table <- as.data.table(.self$.hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
      setnames(leaf_ordering_table, c("leaf", "otu_index"))
      
      for(f in seq(1,length(feature_order))){
        nodes_at_level <- unique(.self$.hierarchy_tree[,feature_order[f]])
        for(n in nodes_at_level){
          list_of_leaves <- .self$.leaf_of_table[node_label==n,leaf]
          leaf_indexes <- leaf_ordering_table[leaf %in% list_of_leaves, otu_index]
          id <- paste(as.hexmode(f-1), as.hexmode(min(leaf_indexes)), sep="-")
          replace_indices <- which(table_node_ids[,feature_order[f]] == n)
          table_node_ids[,feature_order[f]][replace_indices] <- id
        }
      }
      
      return(as.data.table(table_node_ids))
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
