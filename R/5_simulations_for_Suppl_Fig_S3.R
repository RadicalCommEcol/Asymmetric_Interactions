
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)

source("R/matrix_year_i.R")
source("R/generate_random_feasible_growth_rates.R")
source("R/first_species_excluded_ODE.R")


# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

# Define symmetric interaction matrix--------------------------------------------
A_int <- matrix(rep(0,9), ncol=3)
A_int[,1] <- -c(1,0,0)
A_int[,2] <- -c(0.3406884710289558, 0.9328944522580684, 0.11678744219579273)
A_int[,3] <- -c(0.3406884710289558, 0.12410827817548943, 0.9319487652208257)

species_names <- paste0("sp_",1:ncol(A_int))
colnames(A_int) <- species_names
rownames(A_int) <- species_names


# Exclusion prob from Anisofun
prob_excl_AnisoFun_symmetric <- prob_extinction_4_int_matrix(A_int, 1e4, 1e4)

# Exclusion prob from ODE
number_random_points <- 5e4

number_species <- ncol(A_int)

set.seed(1234)
random_feasible_growth_rates <- generate_random_feasible_growth_rates(A_init = A_int, random_points=number_random_points)

first_sp_excluded_results_symmetric <- foreach (i=1:nrow(random_feasible_growth_rates), .combine=c,  
         .packages = c("tidyverse","matlib",
                       "zipfR","pracma",
                       "CValternatives",
                       "anisoFun","deSolve")) %dopar% {
  growth_rates_vector <- random_feasible_growth_rates[i,1:number_species] %>% as.numeric() %>% unname()
  initial_abundances <- random_feasible_growth_rates[i,(number_species+1):(2*number_species)] %>% as.numeric() %>% unname()
  species_excluded_i <- first_species_excluded_ODE(A_int, growth_rates_vector,initial_abundances, delta_t = 0.1, tmax = 1000)
  }

prob_excl_ODE_symmetric <- table(first_sp_excluded_results_symmetric)/number_random_points
prob_excl_ODE_symmetric

# prob_excl_ODE_symmetric
# sp_1    sp_2    sp_3 
# 0.33584 0.33016 0.33400

# prob_excl_AnisoFun_symmetric
# # A tibble: 3 x 4
#    species prob_excl_mean prob_excl_lowerCI prob_excl_upperCI
#   <chr>            <dbl>             <dbl>             <dbl>
# 1 sp_1             0.333             0.333             0.333
# 2 sp_2             0.333             0.333             0.333
# 3 sp_3             0.333             0.333             0.333

# Define asymmetric interaction matrix--------------------------------------------
A_int2 <- matrix(rep(0,9), ncol=3)
A_int2[,1] <- -c(1,1.5,.1)
A_int2[,2] <- -c(0.1, 1, .6)
A_int2[,3] <- -c(1.6, .5, 1)

species_names <- paste0("sp_",1:ncol(A_int2))
colnames(A_int2) <- species_names
rownames(A_int2) <- species_names


# Exclusion prob from Anisofun
prob_excl_AnisoFun_asymmetric <- prob_extinction_4_int_matrix(A_int2, 1e4, 1e4)
prob_excl_AnisoFun_asymmetric

# Exclusion prob from ODE
number_random_points_asymm <- 5e4
number_species2 <- ncol(A_int2)

set.seed(1234)
random_feasible_growth_rates_asymm <- generate_random_feasible_growth_rates(A_int2, number_random_points_asymm)
first_sp_excluded_results_asymmetric <- foreach (i=1:nrow(random_feasible_growth_rates_asymm), .combine=c,  
         .packages = c("tidyverse","matlib",
                       "zipfR","pracma",
                       "CValternatives",
                       "anisoFun","deSolve")) %dopar% {
                         growth_rates_vector_asymm <- random_feasible_growth_rates_asymm[i,1:number_species2] %>% as.numeric() %>% unname()
                         initial_abundances_asymm <- random_feasible_growth_rates_asymm[i,(number_species2+1):(2*number_species2)] %>% as.numeric() %>% unname()
                         species_excluded_i_asymm <- first_species_excluded_ODE(A_int2, growth_rates_vector_asymm,initial_abundances_asymm, delta_t = 0.1, tmax = 1000)
                         
                       }


prob_excl_ODE_asymmetric <- table(first_sp_excluded_results_asymmetric)/number_random_points_asymm
prob_excl_ODE_asymmetric

# prob_excl_ODE_asymmetric
# sp_1    sp_2    sp_3 
# 0.40510 0.32606 0.26884

# # A tibble: 3 x 4
# species    prob_excl_mean prob_excl_lowerCI prob_excl_upperCI
#   <chr>            <dbl>             <dbl>             <dbl>
# 1 sp_1             0.403             0.403             0.403
# 2 sp_2             0.327             0.327             0.327
# 3 sp_3             0.270             0.270             0.270




##########################################################################
##########################################################################
# Caracoles interaction matrix--------------------------------------------

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))

year_i <- 2018
A_int3 <- matrix_year_i(matrix_entries_raw, year_i)



# Exclusion prob from Anisofun
prob_excl_AnisoFun_caracoles <- prob_extinction_4_int_matrix(A_int3, 1e4, 1e4)
prob_excl_AnisoFun_caracoles

# Exclusion prob from ODE
number_random_points_caracoles <- 5e4
number_species3 <- ncol(A_int3)
name_species3 <- colnames(A_int3)

set.seed(1234)
random_feasible_growth_rates_caracoles <- generate_random_feasible_growth_rates(A_int3, number_random_points_caracoles)

first_sp_excluded_results_caracoles <- foreach (i=1:nrow(random_feasible_growth_rates_caracoles), .combine=c,  
                                                 .packages = c("tidyverse","matlib",
                                                               "zipfR","pracma",
                                                               "CValternatives",
                                                               "anisoFun","deSolve")) %dopar% {
                                                                 growth_rates_vector_caracoles <- random_feasible_growth_rates_caracoles[i,1:number_species3] %>% as.numeric() %>% unname()
                                                                 initial_abundances_caracoles <- random_feasible_growth_rates_caracoles[i,(number_species3+1):(2*number_species3)] %>% as.numeric() %>% unname()
                                                                 species_excluded_i_caracoles <- first_species_excluded_ODE(A_int3, growth_rates_vector_caracoles,initial_abundances_caracoles, delta_t = 0.1, tmax = 1000)
                                                                 
                                                               }


prob_excl_ODE_caracoles <- table(first_sp_excluded_results_caracoles)/number_random_points_caracoles
prob_excl_ODE_caracoles


# write_csv(random_feasible_growth_rates_caracoles,"Data/caracoles_processed_data/random_feasible_growth_rates_caracoles.csv") # Commented for security
random_feasible_growth_rates_caracoles <- read_csv("Data/caracoles_processed_data/random_feasible_growth_rates_caracoles.csv")


closest_to_exclusion_estimation_caracoles <- foreach (i=1:nrow(random_feasible_growth_rates_caracoles), .combine=c,  
                                                .packages = c("tidyverse","matlib",
                                                              "zipfR","pracma",
                                                              "CValternatives",
                                                              "anisoFun","deSolve")) %dopar% {
                                                                growth_rates_vector_caracoles <- random_feasible_growth_rates_caracoles[i,1:number_species3] %>% as.numeric() %>% unname()
                                                                closest_to_exclusion_i_caracoles <- closest_to_exclusion(growth_rates_vector_caracoles,A_int3)
                                                                number_species_excluded <- which(min(closest_to_exclusion_i_caracoles[1,(number_species3+1):(2*number_species3)])==
                                                                                                   closest_to_exclusion_i_caracoles[1,(number_species3+1):(2*number_species3)]) %>% unname()
                                                                species_excluded_i_caracoles <- name_species3[number_species_excluded]
                                                                
                                                              }

prob_excl_Aniso2_caracoles <- table(closest_to_exclusion_estimation_caracoles)/number_random_points_caracoles
prob_excl_Aniso2_caracoles

# prob_excl_AnisoFun_caracoles
# # A tibble: 8 x 4
#   species prob_excl_mean prob_excl_lowerCI prob_excl_upperCI
#   <chr>            <dbl>             <dbl>             <dbl>
# 1 BEMA            0.0952            0.0946            0.0958
# 2 CETE            0.284             0.283             0.285 
# 3 HOMA            0.110             0.109             0.110 
# 4 LEMA            0.0449            0.0446            0.0453
# 5 PAIN            0.123             0.122             0.124 
# 6 POMA            0.0810            0.0804            0.0815
# 7 PUPA            0.0532            0.0529            0.0535
# 8 SASO            0.209             0.208             0.210 

# prob_excl_ODE_caracoles
#    BEMA    CETE    HOMA    LEMA    PAIN    POMA    PUPA    SASO 
# 0.17020 0.08102 0.00762 0.07382 0.02994 0.21610 0.14054 0.28076

# closest_to_exclusion_estimation_caracoles
#   BEMA    CETE    HOMA    LEMA    PAIN    POMA    PUPA    SASO 
# 0.09934 0.26618 0.13560 0.05096 0.13112 0.08532 0.05888 0.17260

stopImplicitCluster()