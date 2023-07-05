library(tidyverse)
library(anisoFun)
source("Scripts/aux_functions/isotropy_index.R")
source("Scripts/aux_functions/isotropy_metrics_4_int_matrix.R")
source("Scripts/aux_functions/matrix_year_i.R")
################################################################################
################################################################################

# We create the metaweb for each year

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


years_included <- matrix_entries_raw$year %>% unique()

data_species <- NULL
data_kurtosis <- NULL
data_skewness <- NULL

for(year_i in years_included){
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  elements_A <- as.vector(A_int) 
  
  data_species <- c(data_species,ncol(A_int))
  data_kurtosis <- c(data_kurtosis,moments::kurtosis(elements_A))
  data_skewness <- c(data_skewness,moments::skewness(elements_A))
  
  cat("richness: ",ncol(A_int),", kurtosis: ",moments::kurtosis(elements_A),", skewness: ",moments::skewness(elements_A), "\n")
  
}
mean(data_species[data_species>8])

################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################


community_richness <- c(3,6,10)
kurtosis_index <- c(0.2,1,20)
number_repetitions <- 201

number_Omega_replicates <-  55000
number_boot_replicates <- number_Omega_replicates

set.seed(1234)

results <- NULL

for( S in community_richness){
  
  for(k in kurtosis_index ) {
    
    for(i in 1:number_repetitions){
      
      
      points_SHASH_i <- gamlss.dist::rSHASH(S*S,mu=0,sigma=1,nu=k,tau=k)
      
      # while (sd(points_SHASH_i)<0.4) {
      #   points_SHASH_i <- gamlss.dist::rSHASH(S*S,mu=0,sigma=1,nu=k,tau=k)
      # }
      
      A_int <- matrix(points_SHASH_i, ncol = S)
      
      
      repetition_i_partial_info <- tibble(richness=S,nu=k,tau=k,kurtosis=moments::kurtosis(points_SHASH_i),repetition=i)
      isotropy_metrics_i <- anisotropy_metrics_4_int_matrix(A_int, number_Omega_replicates, number_boot_replicates)
      repetition_i_full_info <- bind_cols(repetition_i_partial_info, isotropy_metrics_i)
      
      results <- bind_rows(results, repetition_i_full_info)
      
    }
    
    write_csv(results, paste0("Results/test_asymmetry_interactions_FD_",number_repetitions,".csv"))
    
  }
  
  
}
