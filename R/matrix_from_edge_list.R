
# generate a binary or weighted matrix from an edge list
# needs fields "node_from","node_to", and optionally, "link_value" for weights
matrix_from_edge_list <- function(edge.list, weighted = FALSE){
  
  sp.from <- unique(edge.list$node_from)
  sp.to <- unique(edge.list$node_to)
  
  A <- matrix(0,nrow = length(sp.to),ncol = length(sp.from),dimnames = list(sp.to,sp.from))
  if(weighted){
    for(i in 1:nrow(edge.list)){
      A[edge.list[i,"node_to"],edge.list[i,"node_from"]] <- edge.list$link_value[i]
    }
  }else{
    for(i in 1:nrow(edge.list)){
      A[edge.list$node_to[i],edge.list$node_from[i]] <- 1
    }
  }
  # A[edge.list[,c("node_to","node_from")]] <- edge.list[,"link_value"]
  return(A)
  
}
