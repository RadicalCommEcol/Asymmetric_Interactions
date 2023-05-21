
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)


################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################



# Define an asymmetric interaction matrix: the reference------------------------
total_number_species = 4
repetitions <- 100

number_Omega_replicates <- 1000
number_boot_replicates <- number_Omega_replicates

matrix_full_info <- NULL
set.seed(1234)

for (repetition_i in 1:repetitions) {
  A_int <- matrix(-runif(total_number_species*total_number_species),
                  ncol = total_number_species)
  
  incenter_inradius_isoprob_A_int <- anisoFun::incenter_inradius_isoprob_calculation(A_int)
  incenter_A_int <- incenter_inradius_isoprob_A_int[[1]]
  
  
  # Extract the vertices from the asymmetric interaction matrix and compute the
  # distances among them-------------------
  
  vertices_coordinates_full_info <- vertices_unit_ball(A_int)
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
  
  pairs_vertices$species <- NA
  for (pair.i in 1:nrow(pairs_vertices)) {
    pairs_vertices$species[pair.i] <- paste0("sp_",
                                             which(!1:total_number_species %in% 
                                                     c(pairs_vertices$V1[pair.i],
                                                       pairs_vertices$V2[pair.i]))
    )
  }
  
  # Compute asymmetry index
  prob_exclusion_A_int <- anisoFun::prob_extinction_4_int_matrix(A_int, number_Omega_replicates, number_boot_replicates)
  
  matrix_full_info_i <- prob_exclusion_A_int %>% left_join(pairs_vertices, by = "species")
  matrix_full_info <- bind_rows(matrix_full_info,matrix_full_info_i)
}


cor.test(matrix_full_info$prob_excl_mean,matrix_full_info$arc_distance)
cor.test(matrix_full_info$prob_excl_mean,matrix_full_info$arc_distance, method = "spearman")

ggplot(matrix_full_info, aes(x=prob_excl_mean, y=arc_distance))+
  geom_point()

