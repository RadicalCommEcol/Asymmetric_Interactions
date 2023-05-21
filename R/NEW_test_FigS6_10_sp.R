
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)
library(nleqslv)

source("R/isotropy_index.R")
source("R/isotropy_metrics_4_int_matrix.R")

################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################



# Define an asymmetric interaction matrix: the reference------------------------
total_number_species = 10

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
# cos_variance <- c(cos_variance,cos_variance_aux)
# J_index <- c(J_index,J_index_aux)


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
  
  
  # new_vertices_prop_triangle <- function(x){
  #   
  #   V1 <- x[1:total_number_species]
  #   V2 <- x[(1+1*total_number_species):(2*total_number_species)]
  #   V3 <- x[(1+2*total_number_species):(3*total_number_species)]
  #   V4 <- x[(1+3*total_number_species):(4*total_number_species)]
  #   
  #   d12 <- V2-V1
  #   d13 <- V3-V1
  #   d14 <- V4-V1
  #   d23 <- V3-V2
  #   d24 <- V4-V2
  #   d34 <- V4-V3
  #   
  #   d21 <- -d12
  #   d31 <- -d13
  #   d41 <- -d14
  #   d32 <- -d23
  #   d42 <- -d24
  #   d43 <- -d34
  #   
  #   path <- d21+d32+d43+d41
  #   
  #   return(c(sum(V1*V1)-1, #unitball
  #            sum(V2*V2)-1, #unitball
  #            sum(V3*V3)-1, #unitball
  #            sum(V4*V4)-1, #unitball
  #            sum(V1*V2)-cos(new_pairs_vertices$arc_distance[1]), # arc distance
  #            sum(V1*V3)-cos(new_pairs_vertices$arc_distance[2]), # arc distance
  #            sum(V1*V4)-cos(new_pairs_vertices$arc_distance[3]), # arc distance
  #            sum(V2*V3)-cos(new_pairs_vertices$arc_distance[4]), # arc distance
  #            sum(V2*V4)-cos(new_pairs_vertices$arc_distance[5]), # arc distance
  #            sum(V3*V4)-cos(new_pairs_vertices$arc_distance[6]), # arc distance
  #            sum(d12*d13)-sum(d13*d14), # angle
  #            sum(d12*d13)-sum(d12*d14),# angle
  #            abs(path[1]), #closed path
  #            abs(path[2]), #closed path
  #            abs(path[3]), #closed path
  #            sum(V1)-1/sqrt(total_number_species)
  #   ))
  #   
  # }
  
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



all_proportions <- proportions_4_tests #c(1,proportions_4_tests)

data_resutls = tibble(cos_variance = cos_variance,
                      J_index = J_index,
                      proportions = all_proportions,
                      small_omega_mean = small_omega_mean,
                      convergence_solver=convergence_code)


write_csv(data_resutls,paste0("results_shape_indices_comparisson_D",total_number_species,".csv"))

coeff <- 1 # Value used to transform the data
cos_varianceeColor <- "#009E73"
J_indexColor <- "#E69F00"

ggplot(data = data_resutls,
       aes(x=proportions)) +
  
  geom_line( aes(y=cos_variance), size=2, color=cos_varianceeColor) + 
  geom_line( aes(y=J_index), size=2, color=J_indexColor) +
  geom_errorbar(aes(ymin=J_index_lower, ymax=J_index_upper), width=.1)+
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Variance of the cosine\nof the side lengths",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Asymmetry index J'(A)")
  ) + 
  scale_x_continuous(name=paste0("Ratio between the side lengths of a ",total_number_species,"-D \nspherical polytope and those of the reference polytope"))+
  theme_bw() +
  
  theme(
    axis.title.y = element_text(color = cos_varianceeColor, size=17),
    axis.title.y.right = element_text(color = J_indexColor, size=17),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold")
  )



# png("Images/cos_side_lengths_VS_J_index_2.png",
#     width = 11.69*.5, # The width of the plot in inches
#     height = 11.69*.6, units = "in", res=300*2)
ggplot(data = data_resutls,
       aes(x=proportions)) +
  
  geom_line( aes(y=cos_variance,color="Variance of the cosine of the side lengths"), size=2) + 
  geom_line( aes(y=J_index,color="Asymmetry index J'(A)"), size=2) +
  geom_errorbar(aes(ymin=J_index_lower, ymax=J_index_upper), width=.1)+
  scale_color_manual(name = NULL, values = c("Asymmetry index J'(A)" = "#E69F00", 
                                             "Variance of the cosine of the side lengths" = "#009E73"))+
  labs(x=paste0("Ratio between the side lengths of a ",total_number_species,"-D \nspherical polytope and those of the reference polytope"),
       y=NULL)+
  theme_bw() +
  theme(
    legend.position="top",
    legend.text=element_text(size=14,face="bold"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold")
  )+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

# dev.off()

