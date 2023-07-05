

# Estimation of the probabilities of exclusion for all the species in Caracoles,
# across sampling years. .

# To build the interaction matrix of the community in a given year
# from raw data we use the following auxiliary function: "Scripts/aux_functions/matrix_year_i.R"

# To bootstrap the probabilities of exclusion of the community in a given year
# we use the following auxiliary function: "Scripts/aux_functions/exclusion_probabilities_par_bootstrap.R"

# INPUT: 
# interaction matrices: "Data/caracoles_raw_data/alpha_heterogeneous_time.csv"

# OUTPUT:
# species' probabilities of exclusion: "Results/caracoles_processed_data/caracoles_boot_plants_year_rep_10000.csv"



library(matlib) # to multiply matrices
library(tidyverse)
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(MultitrophicFun)
library(parallel)
library(foreach)
library(doParallel)
library(boot)
library(anisoFun)

# Functions to run calculations about the isotropic area
source("Scripts/aux_functions/matrix_year_i.R")
source("Scripts/aux_functions/exclusion_probabilities_par_bootstrap.R")

simple_mean <- function(x, indices){
  return(sum(x[indices])/length(indices))
}

# We create the metaweb for each year

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


number_Omega_replicates <- 10000
number_boot_replicates <- number_Omega_replicates

years_included <- matrix_entries_raw$year %>% unique()

caracoles_prob_exclusion_plants <- NULL

numCores <- detectCores()

registerDoParallel(numCores-4)

set.seed(123)

for(year_i in years_included){
    
    cat(year_i, "\n")
    
    A_int <- matrix_year_i(matrix_entries_raw, year_i)
    
    cat(colnames(A_int), "\n")

    dimensions <- ncol(A_int)
    
    incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
    
    I <- incenter_inradius_isoprob[[1]]
    
    center <- I
    
    # Create tibble to store year, species names, and replicate's data
    
    name_colums <- c("year","species", paste0("prob_excl_rep_",1:number_boot_replicates),
                     paste0("Omega_rep_",1:number_boot_replicates))
    
    caracoles_prob_exclusion_plants_mat <- matrix(rep(0, dimensions*(2 + 2*number_boot_replicates)),
                                                  nrow = dimensions, 
                                                  ncol = (2 + 2*number_boot_replicates))
    
    bootstrap_raw_data <-
      exclusion_probabilities_par_bootstrap(A_int, replicates = number_Omega_replicates)
    
    for (sp_i in 1:dimensions) {
      
      bootmean_area_excl_sp_i <- boot::boot(data = bootstrap_raw_data[sp_i,], 
                                      statistic = simple_mean, 
                                      R = number_boot_replicates)
      
      caracoles_prob_exclusion_plants_mat[sp_i, 3:(2+number_boot_replicates)] <- 
        bootmean_area_excl_sp_i[["t"]] 
      
    }

    Omega_values <- colSums(caracoles_prob_exclusion_plants_mat[, c(3:(2+number_boot_replicates))])
    
    caracoles_prob_exclusion_plants_mat[, c((3+number_boot_replicates):(2 + 2*number_boot_replicates))] <-
      matrix(rep(Omega_values,dimensions), nrow = number_boot_replicates) %>% t()
    
    for (sp_i in 1:nrow(caracoles_prob_exclusion_plants_mat)) {
      caracoles_prob_exclusion_plants_mat[sp_i, c(3:(2+number_boot_replicates))] <-
        caracoles_prob_exclusion_plants_mat[sp_i, c(3:(2+number_boot_replicates))]/
        caracoles_prob_exclusion_plants_mat[sp_i, c((3+number_boot_replicates):(2 + 2*number_boot_replicates))]
    }
    
    caracoles_prob_exclusion_plants_aux <- as_tibble(data.frame(caracoles_prob_exclusion_plants_mat))
    colnames(caracoles_prob_exclusion_plants_aux) <- name_colums
    
    caracoles_prob_exclusion_plants_aux$year = year_i
    caracoles_prob_exclusion_plants_aux$species = colnames(A_int)
    
    caracoles_prob_exclusion_plants <- rbind(caracoles_prob_exclusion_plants,
                                             caracoles_prob_exclusion_plants_aux)
    
    write_csv(caracoles_prob_exclusion_plants,
              paste0("Results/caracoles_processed_data/caracoles_boot_plants_year_rep_",
              number_Omega_replicates,".csv"))

  
}
