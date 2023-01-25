
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
source("R/interaction_strength_matrix_from_edge_list.R")
source("R/exclusion_probabilities_par_bootstrap.R")


simple_mean <- function(x, indices){
  return(sum(x[indices])/length(indices))
}

# Load link information
matrix_links_raw <- read_csv("Data/mangal_processed_data/NO_NAs_mangal_links_processed_kingdom.csv")
networks_included <- matrix_links_raw$network_id %>% unique()

# Sanity check
all(matrix_links_raw$kingdom_from != matrix_links_raw$kingdom_to)

# Set calculation parameters
number_Omega_replicates <- 10000
number_boot_replicates <- number_Omega_replicates

numCores <- detectCores()
registerDoParallel(numCores-4)

# Tibble to save results
mangal_prob_exclusion <- NULL

set.seed(123)

for(network_i in networks_included){
    
    cat(network_i, "\n",
        which(networks_included == network_i), "\n")
  
    matrix_edge_list_i <- matrix_links_raw %>% filter(network_id == network_i) %>%
      dplyr::select(taxonomy_name_from,taxonomy_name_to,kingdom_from,kingdom_to) %>%
      rename(node_from = taxonomy_name_from, node_to = taxonomy_name_to)
    
    rho = 0.01
    delta = 0
    
    source("R/mutualistic_strength_solver.R")
    
    xstart <- c(1)
    mutualistic_strength_results <- nleqslv(xstart,mutualistic_strength_solver)
    
    gamma_0_max <- mutualistic_strength_results$x
    
    A_int <- interaction_strength_matrix_from_edge_list(matrix_edge_list_i,
                                                        rho = rho, 
                                                        delta = delta, 
                                                        gamma_0 = 0.5*gamma_0_max)
    
    cat(colnames(A_int), "\n")

    dimensions <- ncol(A_int)
    
    incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
    
    I <- incenter_inradius_isoprob[[1]]
    
    center <- I
    
    # Create tibble to store year, species names, and replicate's data
    
    name_colums <- c("network_id","species", paste0("prob_excl_rep_",1:number_boot_replicates),
                     paste0("Omega_rep_",1:number_boot_replicates))
    
    mangal_prob_exclusion_mat <- matrix(rep(0, dimensions*(2 + 2*number_boot_replicates)),
                                                  nrow = dimensions, 
                                                  ncol = (2 + 2*number_boot_replicates))
    
    bootstrap_raw_data <-
      exclusion_probabilities_par_bootstrap(A_int, replicates = number_Omega_replicates)
    
    for (sp_i in 1:dimensions) {
      
      bootmean_area_excl_sp_i <- boot::boot(data = bootstrap_raw_data[sp_i,], 
                                      statistic = simple_mean, 
                                      R = number_boot_replicates)
      
      mangal_prob_exclusion_mat[sp_i, 3:(2+number_boot_replicates)] <- 
        bootmean_area_excl_sp_i[["t"]] 
      
    }

    Omega_values <- colSums(mangal_prob_exclusion_mat[, c(3:(2+number_boot_replicates))])
    
    mangal_prob_exclusion_mat[, c((3+number_boot_replicates):(2 + 2*number_boot_replicates))] <-
      matrix(rep(Omega_values,dimensions), nrow = number_boot_replicates) %>% t()
    
    for (sp_i in 1:nrow(mangal_prob_exclusion_mat)) {
      mangal_prob_exclusion_mat[sp_i, c(3:(2+number_boot_replicates))] <-
        mangal_prob_exclusion_mat[sp_i, c(3:(2+number_boot_replicates))]/
        mangal_prob_exclusion_mat[sp_i, c((3+number_boot_replicates):(2 + 2*number_boot_replicates))]
    }
    
    mangal_prob_exclusion_aux <- as_tibble(data.frame(mangal_prob_exclusion_mat))
    colnames(mangal_prob_exclusion_aux) <- name_colums
    
    mangal_prob_exclusion_aux$network_id = network_i
    mangal_prob_exclusion_aux$species = colnames(A_int)
    
    mangal_prob_exclusion <- rbind(mangal_prob_exclusion,
                                             mangal_prob_exclusion_aux)
    
    write_csv(mangal_prob_exclusion,
              paste0("Data/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_mangal_boot_networks_rep_",
              number_Omega_replicates,".csv"))

  
}
