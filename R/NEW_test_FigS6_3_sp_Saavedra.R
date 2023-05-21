
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)
library(nleqslv)
library(boot)

source("R/isotropy_index.R")

################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################

# Variables
asymmetry_index_saavedra <- NULL 

# Define an asymmetric interaction matrix: the reference------------------------
total_number_species = 3

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

number_Omega_replicates <- 100000
number_boot_replicates <- number_Omega_replicates

# Estimate a new system were arc distance among nodes are proportional to the previous ones
proportions_4_tests <- seq(from = 1.0, to = .3, length.out=8)

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
  
  new_pairs_vertices_aux
  
  # Compute new asymmetry index
  
  saavedra_index_new <- var(new_pairs_vertices_aux$arc_distance)
  saavedra_index_new_lowerlimit <- var(new_pairs_vertices_aux$arc_distance)
  saavedra_index_new_upperlimit <- var(new_pairs_vertices_aux$arc_distance)
  
  saavedra_results <- tibble(saavedra_index = saavedra_index_new,
                             saavedra_index_lowerlimit = saavedra_index_new_lowerlimit,
                             saavedra_index_upperlimit = saavedra_index_new_upperlimit,
                             proportions = proportion)
  
  

  asymmetry_index_saavedra <- bind_rows(asymmetry_index_saavedra,saavedra_results)

  write_csv(asymmetry_index_saavedra,
            paste0("Data/results_saavedra_shape_indices_comparisson_D",total_number_species,".csv"))
  
}

coeff <- 1 # Value used to transform the data
saavedraColor <- "#009E73"


ggplot(data = asymmetry_index_saavedra,
       aes(x=proportions)) +
  
  geom_line( aes(y=saavedra_index), size=2, color=saavedraColor) + 
  scale_y_continuous(
    
    # Features of the first axis
    name = "saavedra",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="saaverda")
  ) + 
  scale_x_continuous(name="Ratio between the side lengths of an isosceles\nspherical triangle and \nthose of the reference triangle")+
  theme_bw() +
  
  theme(
    axis.title.y = element_text(color = saavedraColor, size=17),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold")
  )

