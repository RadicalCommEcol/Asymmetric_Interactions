
# Estimation of asymmetry metric proposed by Grilli et al. 2017 
# (DOI: 10.1038/ncomms14389) in Caracoles, across sampling years. 

# To build the interaction matrix of the community in a given year
# from raw data we use the following auxiliary function: "Scripts/aux_functions/matrix_year_i.R"

# INPUT: 
# interaction matrices: "Data/caracoles_raw_data/alpha_heterogeneous_time.csv"

# OUTPUT:
# Asymmetry metric by Grilli et al.: "Results/caracoles_processed_data/caracoles_grilli.csv"

# -------------------------------------------------------------------------

library(matlib) # to multiply matrices
library(tidyverse)
library(pracma) # to solve n-dimensional cross products
library(anisoFun)

# Functions to run calculations on arc distances
source("Scripts/aux_functions/matrix_year_i.R")

# We create the metaweb for each year
matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


years_included <- matrix_entries_raw$year %>% unique()

caracoles_grilli_index <- NULL

for(year_i in years_included){
  
  cat(year_i,"\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  
  total_number_species <- nrow(A_int)
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
  
  grilli_i <- var(cos(pairs_vertices$arc_distance))
  
  caracoles_grilli_index <- c(caracoles_grilli_index, grilli_i)
  
  
}

caracoles_results <- tibble(year= years_included,
                            grilli_index=caracoles_grilli_index)
write_csv(caracoles_results, paste0("Results/caracoles_processed_data/caracoles_grilli.csv"))
