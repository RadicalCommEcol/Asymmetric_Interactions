

# Estimation of exclusion distances for all the species in Caracoles,
# across sampling years. To do so, we first estimate the growth vector
# of the community in a given year.

# To build the interaction matrix of the community in a given year
# from raw data we use the following auxiliary function:
# "Scripts/aux_functions/matrix_year_i.R"

# INPUT: 
# interaction matrices: "Data/caracoles_raw_data/alpha_heterogeneous_time.csv"
# growth rates: "Data/caracoles_raw_data/lambda_heterogeneous_both.csv" and 
# "Data/caracoles_raw_data/lambda_heterogeneous_both.csv"
# plant abundances: "Data/caracoles_raw_data/abundance.csv"
# plant traits: "Data/caracoles_raw_data/01_05_plant_species_traits.csv"

# OUTPUT:
# species exclusion distances: "Results/caracoles_processed_data/caracoles_arc_distance.csv"

# -------------------------------------------------------------------------


library(matlib) # to multiply matrices
library(tidyverse)
library(pracma) # to solve n-dimensional cross products
library(anisoFun)

# Functions to run calculations on exclusion distances
source("Scripts/aux_functions/matrix_year_i.R")

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

caracoles_exclusion_distances <- NULL

for(year_i in years_included){
  
  cat(year_i,"\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  
  sp_in_A_int <- colnames(A_int)
  
  data_year_i <- growth_rates_raw %>% filter(year == year_i, focal %in%  sp_in_A_int)
  
  incenter_data <- incenter_inradius_isoprob_calculation(A_int)
  incenter_arc_distance <- incenter_data[[2]]
  
  #growth_rates_i_LV <- (-1) * A_int %*% matrix(data_year_i$abundance, ncol = 1)
  
  growth_rates_i_LV <- data_year_i$r_LV
  
  r_vector_i <- growth_rates_i_LV/sqrt(sum((growth_rates_i_LV)^2))
  
  populations_growth_rates_i <- (-1)*matlib::inv(A_int) %*% r_vector_i
  
  if(!all(populations_growth_rates_i>0)){
    
    populations_growth_rates_i[populations_growth_rates_i < 0] <- 1
    
    populations_growth_rates_i <- round( populations_growth_rates_i,0)
    
    growth_rates_i_LV_adapted <- (-1)*A_int %*% populations_growth_rates_i
    
    r_vector_i <- growth_rates_i_LV_adapted/sqrt(sum((growth_rates_i_LV_adapted)^2))
    
  }
  
  arc_distance <- sp_distance_to_exclusion(r_vector_i,A_int)
  arc_distance$incenter_arc_distace = incenter_arc_distance
  arc_distance$relative_arc_distance = arc_distance$arc_distance/incenter_arc_distance
  arc_distance$year = year_i
  
  arc_distance <- arc_distance[,c("species","year","arc_distance","incenter_arc_distace","relative_arc_distance")]
  
  caracoles_exclusion_distances <- bind_rows(caracoles_exclusion_distances, arc_distance)
  
  write_csv(caracoles_exclusion_distances, paste0("Results/caracoles_processed_data/caracoles_arc_distance.csv"))

  
  
}
