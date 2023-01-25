
# generate a interaction strength matrix from an edge list that contains the
# following fields "node_from","node_to", kingdom_from, kingdom_to. In addition,
# the following parameters (non-negative numbers) are needed: rho, delta, and 
# gamma_0

interaction_strength_matrix_from_edge_list <- function(edge.list, rho, 
                                                       delta, gamma_0){
  
  if(all(edge.list$kingdom_from == "Animalia")){
    
    animals <- unique(edge.list$node_from) %>% sort()
    plants <- unique(edge.list$node_to) %>% sort()
    
    animals_col <- 1
    plants_col <- 2
  
  }else if(all(edge.list$kingdom_from == "Plantae")){
    
    animals <- unique(edge.list$node_to) %>% sort()
    plants <- unique(edge.list$node_from) %>% sort()
    animals_col <- 2
    plants_col <- 1
    
  }
  
  species_names <- c(animals,plants)
  
  degree_nodes <- tibble(node = c(edge.list$node_from,edge.list$node_to)) %>% 
    group_by(node) %>% 
    count() %>% rename(degree = n)
  
  A <- matrix(0,nrow = length(species_names),ncol = length(species_names),
              dimnames = list(species_names,species_names))
  
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      
      if(i==j){
        
        A[i,j] <- 1
        
      }else if((species_names[i] %in% animals) == (species_names[j] %in% animals)){
        
        A[i,j] <- rho
        
      }else{
        
        if(species_names[i] %in% animals){
          
          rows_with_sp_i <- which(edge.list[,animals_col] == species_names[i])
          
          if(any(edge.list[rows_with_sp_i,plants_col]==species_names[j])){
            
            degree_sp_i <- as.numeric(
              degree_nodes$degree[degree_nodes$node == species_names[i]]
            )
            
            A[i,j] <- (-1) * gamma_0 / (degree_sp_i ^ delta)
          }
          
        }else{
          
          rows_with_sp_i <- which(edge.list[,plants_col] == species_names[i])
          
          if(any(edge.list[rows_with_sp_i,animals_col]==species_names[j])){
            
            degree_sp_i <- as.numeric(
              degree_nodes$degree[degree_nodes$node == species_names[i]]
            )
            
            A[i,j] <- (-1) * gamma_0 / (degree_sp_i ^ delta)
          }
          
        }
        
        
      }
  
    }
   
  }

  return(-A) #We add a minus sign to follow the notation in Song et al. 2018 JTB
  
}
