
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)
library(nleqslv)

source("Scripts/aux_functions/isotropy_index.R")
source("Scripts/aux_functions/isotropy_metrics_4_int_matrix.R")

################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################



# Define an asymmetric interaction matrix: the reference------------------------
total_number_species = 6

set.seed(1234)
A_int <- -1*diag(rep(1, total_number_species))

incenter_inradius_isoprob_A_int <- anisoFun::incenter_inradius_isoprob_calculation(A_int)
incenter_A_int <- incenter_inradius_isoprob_A_int[[1]]

A_asym <- A_int
A_asym[,1] <- -(incenter_A_int)


# Extract the vertices from the asymmetric interaction matrix and compute the
# distances among them-------------------

vertices_coordinates_full_info <- vertices_unit_ball(A_asym)
vertices_coordinates <- vertices_coordinates_full_info[,(total_number_species+1):(2*total_number_species)]

arc_distance_vertices <- function(vertices_coordinates,i,j){
  
  return(acos(sum(vertices_coordinates[i,]*vertices_coordinates[j,])))
  
}

pairs_vertices <- t(combn(x=1:total_number_species, m=2)) %>% as_tibble()
pairs_vertices$arc_distance <- NA
for (pair.i in 1:nrow(pairs_vertices)) {
  pairs_vertices$arc_distance[pair.i] <- arc_distance_vertices(vertices_coordinates,
                                                               pairs_vertices$V1[pair.i],
                                                               pairs_vertices$V2[pair.i])
}

# Compute asymmetry index

set.seed(123)

number_Omega_replicates <- 100
number_boot_replicates <- number_Omega_replicates

isotropy_metrics_results <- isotropy_metrics_4_int_matrix(A_asym, number_Omega_replicates, number_boot_replicates)

J_index_aux <- isotropy_metrics_results$isotropy_index_mean
cos_variance_aux <- var(cos(pairs_vertices$arc_distance))
small_omega_mean_aux <- isotropy_metrics_results$small_omega_mean


# Estimate a new system were arc distance among nodes are proportional to the previous ones
proportions_4_tests <- seq(from = 1.1, to = .3, length.out=9)

# Variables for storaging
cos_variance <- NULL
J_index <- NULL 
J_index_lower <- NULL
J_index_upper <- NULL
small_omega_mean <- NULL
convergence_code <- NULL

for(proportion in proportions_4_tests){
  
  new_pairs_vertices <- pairs_vertices
  new_pairs_vertices$arc_distance <- proportion * pairs_vertices$arc_distance

  new_vertices_prop_triangle_V2 <- function(x){
    
    A_aux_trans <- matrix(x,ncol = total_number_species)
    A_aux <-t(A_aux_trans)
    
    A_aux_x_A_aux_trans <- A_aux %*% A_aux_trans
    
    unit_ball_conditions <- diag(A_aux_x_A_aux_trans)-1
    
    pairs_vertices_aux <- t(combn(x=1:total_number_species, m=2)) %>% as_tibble()
    
    arc_distance_conditions <- NULL
    
    for(pair.i in 1:nrow(pairs_vertices_aux)){
      
      arc_dis_aux <- A_aux_x_A_aux_trans[pairs_vertices_aux$V1[pair.i],
                                         pairs_vertices_aux$V2[pair.i]]-cos(new_pairs_vertices$arc_distance[pair.i])
      
      arc_distance_conditions <- c(arc_distance_conditions,arc_dis_aux)
    }
    
    distances_matrix <- matrix(rep(0,total_number_species*nrow(pairs_vertices_aux)),ncol = total_number_species)
    
    for(pair.i in 1:nrow(pairs_vertices_aux)){
      
      distances_matrix[pair.i,] <- A_aux[pairs_vertices_aux$V2[pair.i],] - A_aux[pairs_vertices_aux$V1[pair.i],]

    }
    
    for(i in 1:total_number_species){
      
      if(i==1){
        pair_selected <- which(pairs_vertices_aux$V1==i & pairs_vertices_aux$V2==(i+1))
        path_condition <- distances_matrix[pair_selected,]
        
      }else if((i<total_number_species)&(1<i)){
        pair_selected <- which(pairs_vertices_aux$V1==i & pairs_vertices_aux$V2==(i+1))
        path_condition <- path_condition + distances_matrix[pair_selected,]
        
      }else{
        pair_selected <- which(pairs_vertices_aux$V1==1 & pairs_vertices_aux$V2==i)
        path_condition <- path_condition - distances_matrix[pair_selected,]
        
      }
    }
    
    closed_path_condition <- path_condition[1:(total_number_species-1)]
    
    pairs_of_pairs_vertices_aux <- t(combn(x=2:total_number_species, m=2)) %>% as_tibble()
    
    angles_pairs_of_pairs_vertices <- NULL
    
    for(p_of_p.i in 1:nrow(pairs_of_pairs_vertices_aux)){
      
      row_distance_matrix_first_pair <- which(pairs_vertices_aux$V1==1 & pairs_vertices_aux$V2==pairs_of_pairs_vertices_aux$V1[p_of_p.i])
      row_distance_matrix_second_pair <- which(pairs_vertices_aux$V1==1 & pairs_vertices_aux$V2==pairs_of_pairs_vertices_aux$V2[p_of_p.i])
      
      angles_pairs_of_pairs_vertices <- c(angles_pairs_of_pairs_vertices,
                                          sum(distances_matrix[row_distance_matrix_first_pair,]*distances_matrix[row_distance_matrix_second_pair,]))
      
    }
    
    if(length(angles_pairs_of_pairs_vertices)>1){
      angles_comparisson <- t(combn(x=1:length(angles_pairs_of_pairs_vertices), m=2)) %>% as_tibble()
      
      angles_comparisson_conditions_aux <- NULL
      
      for(p_of_p.i in 1:(nrow(angles_comparisson)-1)){
        
        angles_comparisson_conditions_aux <- c(angles_comparisson_conditions_aux,
                                               angles_pairs_of_pairs_vertices[angles_comparisson$V1[p_of_p.i]]-angles_pairs_of_pairs_vertices[angles_comparisson$V2[p_of_p.i]])
        
      }
      
      number_of_needed_conditions <- total_number_species^2-1-length(arc_distance_conditions)-length(unit_ball_conditions)-length(closed_path_condition)
      angles_comparisson_conditions <- angles_comparisson_conditions_aux[1:number_of_needed_conditions]
    }else{
      angles_comparisson_conditions <- NULL
    }
    
    
    return(c(unit_ball_conditions,
             arc_distance_conditions,
             angles_comparisson_conditions,
             closed_path_condition,
             sum(A_aux[1,])-1/sqrt(total_number_species)
    ))
    
  }
  
  set.seed(1234)
  
  init_search <- rep(1,total_number_species*total_number_species) + runif(total_number_species*total_number_species)
  
  search_new_vertices_prop_triangle <- nleqslv(init_search, new_vertices_prop_triangle_V2, control = list(allowSingular=T,maxit=10000))
  search_new_vertices_prop_triangle 
  
  while(search_new_vertices_prop_triangle$termcd>3){
    search_new_vertices_prop_triangle <- nleqslv(init_search, new_vertices_prop_triangle, control = list(allowSingular=T))
    # search_new_vertices_prop_triangle 
  }
  
  new_vertices_prop_triangle_V2(search_new_vertices_prop_triangle$x)
  
  convergence_code <- c(convergence_code,search_new_vertices_prop_triangle$message)
  
  A_asym_new <- A_asym
  
  for (i in 1:nrow(A_asym_new)) {
    A_asym_new[,i] <- -search_new_vertices_prop_triangle$x[(1+(i-1)*total_number_species):(i*total_number_species)]
  }
  
  
  A_asym_new
  
  new_vertices_coordinates <- vertices_coordinates
  
  for (i in 1:nrow(new_vertices_coordinates)) {
    new_vertices_coordinates[i,] <- search_new_vertices_prop_triangle$x[(1+(i-1)*total_number_species):(i*total_number_species)]
  }
  
  
  new_pairs_vertices_aux <- t(combn(x=1:total_number_species, m=2)) %>% as_tibble()
  new_pairs_vertices_aux$arc_distance <- NA
  for (pair.i in 1:nrow(new_pairs_vertices_aux)) {
    new_pairs_vertices_aux$arc_distance[pair.i] <- arc_distance_vertices(new_vertices_coordinates,
                                                                         new_pairs_vertices_aux$V1[pair.i],
                                                                         new_pairs_vertices_aux$V2[pair.i])
  }
  
  new_pairs_vertices_aux$arc_distance - new_pairs_vertices$arc_distance
  
  # Compute new asymmetry index
  
  set.seed(1234)
  
  new_isotropy_metrics_results <- isotropy_metrics_4_int_matrix(A_asym_new, number_Omega_replicates, number_boot_replicates)
  
  J_index_new <- new_isotropy_metrics_results$isotropy_index_mean
  cos_variance_new <- var(cos(new_pairs_vertices_aux$arc_distance))
  
  
  cos_variance <- c(cos_variance,cos_variance_new)
  J_index <- c(J_index,J_index_new)
  J_index_lower <- c(J_index_lower,new_isotropy_metrics_results$isotropy_index_lowerCI)
  J_index_upper <- c(J_index_upper,new_isotropy_metrics_results$isotropy_index_upperCI)
  small_omega_mean <- c(small_omega_mean,new_isotropy_metrics_results$small_omega_mean)
  
}



all_proportions <- proportions_4_tests

data_resutls = tibble(cos_variance = cos_variance,
                      J_index = J_index,
                      proportions = all_proportions,
                      small_omega_mean = small_omega_mean,
                      convergence_solver=convergence_code)


write_csv(data_resutls,paste0("Results/results_shape_indices_comparisson_D",total_number_species,".csv"))
