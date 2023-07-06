
# Estimation of the size and asymmetry of random communities with 3, 6, and 10 
# species, whose interaction strengths have been randomly sampled from a 
# leptokurtic, normal, and platykurtic distribution. For each combination of 
# community size and strength distribution (9 pairs in total), we generated 
# randomly 201 interaction matrices (1,809 matrices in total). The elements of 
# the interaction matrices were sampled from the Sinh-Arcsinh (SHASH) 
# distributions in the R-package gamlss.dist v 6.0-5 
# (Stasinopoulos & Rigby, 2022) with µ = 0, σ = 1, and the following ν = τ
# values: ν = τ = 0.2 for leptokurtic distributions, ν = τ = 1 for a 
# normal distribution, and ν = τ = 20 for platykurtic distributions.

# INPUT: 
# community_richness <- c(3,6,10)
# kurtosis_index <- c(0.2,1,20)
# number_repetitions <- 201

# OUTPUT:
# Data on the size and asymmetry of the feasibility domain of 1,809 communities,
# accros a gradient of richness and interaction strengths: 
# "Results/test_asymmetry_interactions_FD_",number_201.csv"

#-----------------------------------------------------------------------



library(tidyverse)
library(anisoFun)
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
