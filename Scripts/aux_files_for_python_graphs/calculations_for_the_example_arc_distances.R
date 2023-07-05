library(anisoFun)
library(tidyverse)

A_int <- matrix(rep(0,9), ncol=3)
A_int[,1] <- -c(1.5, 0.53452248, 0.2672612)
A_int[,2] <- -c(0.3521804, 2, 0.6163156)
A_int[,3] <- -c(0.3698001, 0.09245003, 3)

incenter_inradius_isoprob_calculation(A_int)


vertices_unit_ball(A_int)
center <- c(0.5, 0.4, 0.6496525)
center <- center/sqrt(sum(center*center))
center
anisoFun::sp_distance_to_exclusion(center,A_int)
anisoFun::closest_to_exclusion(center, A_int)

dimensions <- ncol(A_int)
curves_data <- cone_vertices_director_vertices(A_int)

col_number <- ncol(curves_data)

dir_ver_mat <- as.matrix(curves_data[,(col_number-(dimensions-1)):col_number])

# we check the direction of the director vectors
correction_direction <- NULL

for(i in 1: dimensions){
  
  if(all((dir_ver_mat %*% as.numeric(curves_data[i,(1:dimensions)]) %>% round(12))>=0)){
    correction_direction <- c(correction_direction,1)
  }else{
    correction_direction <- c(correction_direction,-1)
  }
  
}

acr_distance_center_side <- NULL

for(i in 1:nrow(curves_data)){
  
  normal_vector_i <- correction_direction[i] * curves_data[i,c((col_number-(dimensions-1)):col_number)]
  sin_theta <- sum(normal_vector_i*center)
  acr_distance_center_side <- c(acr_distance_center_side,asin(sin_theta))
  
}

vertex_min_arc_dist_index <- c(1,2,3)

tangent_points_matrix <- matrix(rep(0,dimensions*length(vertex_min_arc_dist_index)),
                                nrow = length(vertex_min_arc_dist_index),
                                ncol = dimensions)

colnames(tangent_points_matrix) <- paste0("r_T_",1:dimensions)

for(i in 1:length(vertex_min_arc_dist_index)){
  index <- vertex_min_arc_dist_index[i]
  normal_vector_i <- correction_direction[index] * curves_data[index,c((col_number-(dimensions-1)):col_number)]
  sin_theta <- sum(normal_vector_i*center)
  tangent_aux <- center - sin_theta*normal_vector_i %>% as.numeric()
  module_tangent_aux <- sqrt(sum(tangent_aux*tangent_aux))
  tangent_points_matrix[i,] <- tangent_aux / module_tangent_aux
  
}

tangent_points_matrix
