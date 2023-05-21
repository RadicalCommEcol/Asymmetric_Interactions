library(matlib) # to multiply matrices
library(tidyverse)
library(pracma) # to solve n-dimensional cross products
library(anisoFun)

# Functions to run calculations on arc distances
source("R/matrix_year_i.R")

# We create the metaweb for each year
matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))

growth_rates_raw_homo <- read_csv2("Data/caracoles_raw_data/lambda_heterogeneous_time.csv") %>%
  group_by(year,sp) %>% summarise(lambda = mean(lambda, na.rm =T)) %>%
  rename(focal = sp)

growth_rates_raw_hete <- read_csv2("Data/caracoles_raw_data/lambda_heterogeneous_both.csv") %>%
  group_by(year,sp) %>% summarise(mean_lambda_hete = mean(lambda, na.rm =T)) %>%
  rename(focal = sp)

abundance_raw <- read_csv("Data/caracoles_raw_data/abundance.csv") %>%
  group_by(year, species) %>% summarise(abundance = sum(individuals, na.rm =T)) %>%
  rename(focal = species)

traits_raw <- read_csv2("Data/caracoles_raw_data/01_05_plant_species_traits.csv") %>%
  rename(focal = species.code)

growth_rates_raw <- growth_rates_raw_homo %>%
  left_join(growth_rates_raw_hete, by = c("year", "focal")) %>%
  left_join(abundance_raw, by = c("year", "focal")) %>%
  left_join(traits_raw, by = c("focal")) %>%
  mutate(r_LV = log((1-(1-germination.rate)*seed.survival)/germination.rate) - lambda)

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
write_csv(caracoles_results, paste0("Data/caracoles_processed_data/caracoles_grilli.csv"))
