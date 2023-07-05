# Given an edge list for an bipartite network, this function estimates overall
# level of mutualistic strength from an edge list that contains the following
# fields: "node_from","node_to", kingdom_from, kingdom_to. In addition,
# the following parameters (non-negative numbers) are needed: rho and delta

mutualistic_strength_solver <- function(x,
                                        edge.list=matrix_edge_list_i,
                                        a=rho,b=delta){ 
  
  solution <- NULL # we store here the resulting values
  
  A_int_aux <- interaction_strength_matrix_from_edge_list(edge.list, a, b, x) 
  
  ev <- eigen(-A_int_aux) # We add a minus sign because we follow the notation in Song et al. 2018 JTB
  # extract real part of eigen values
  real_part_eigenvalues <- Re(ev$values)
  
  return(min(real_part_eigenvalues))
  
}
