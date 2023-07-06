

# Estimation of the average computing time (in seconds) for randomly finding 
# feasible intrinsic growth rates in symmetric communities, across a gradient of
# species richness.

# To estimate initial growth rates within the feasibility domain we use the 
# following auxiliary function:
# "Scripts/aux_functions/generate_random_feasible_growth_rates.R"


# INPUT: 
# list_total_number_species <- c(10,20,30,40,50,60,70,80,90,100)

# OUTPUT: Average computing times for the above symmetric communities:
# "Results/simulations/data_experiment_LV_times.csv"

#-----------------------------------------------------------------------

library(matlib) # to multiply matrices
library(tidyverse)
library(scales)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)

source("Scripts/aux_functions/generate_random_feasible_growth_rates.R")


################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################
# Define symmetric interaction matrix--------------------------------------------
list_total_number_species <- c(10,20,30,40,50,60,70,80,90,100)

capital_omega_list <- NULL
time_LV_list <- NULL

for (total_number_species in list_total_number_species) {
  
  A_int <- -1*diag(rep(1,total_number_species))
  
  
  species_names <- paste0("sp_",1:ncol(A_int))
  colnames(A_int) <- species_names
  rownames(A_int) <- species_names
  
  capital_omega_i <- anisoFun::Omega_bootstrap(A_int, replicates = 1e3) %>% mean()
  
  capital_omega_list <- c(capital_omega_list, capital_omega_i)
  
  # Exclusion prob from ODE
  number_random_points <- 100
  number_species <- ncol(A_int)
  
  start_time <- Sys.time()
  set.seed(1234)
  random_feasible_growth_rates <- foreach (i=1:round((number_random_points/1),0), .combine=rbind,  
                                           .packages = c("tidyverse","matlib",
                                                         "zipfR","pracma",
                                                         "CValternatives",
                                                         "anisoFun")) %dopar% {
                                                           generate_random_feasible_growth_rates(A_int, random_points=1)
                                                         }
  
  
  
  end_time <- Sys.time()
  
  time_LV_list <- c(time_LV_list, (end_time - start_time)) # Warning: this output does not inform on the time units: seconds, days,... 
  
}

time_LV_list2 <- c(0.03600192,0.07001400, 0.26200223,1.83304191,9.59698892, 60*1.09100039, 60*7.73155475,60*48.01390030,3600*4.78676046,3600*24*1.46608374)

data_experiment <- tibble(total_number_species = list_total_number_species,
                          capital_omega = capital_omega_list,
                          time_LV_unitary = time_LV_list2/number_random_points)

write_csv(data_experiment,"Results/simulations/data_experiment_LV_times.csv")

stopImplicitCluster()
