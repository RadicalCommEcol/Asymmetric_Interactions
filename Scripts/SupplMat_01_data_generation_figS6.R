
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)

source("Scripts/aux_functions/generate_random_feasible_growth_rates.R")
source("Scripts/aux_functions/first_species_excluded_ODE.R")

################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################
# Define symmetric interaction matrix--------------------------------------------
total_number_species = 15
A_int <- -1*diag(rep(1,total_number_species))


species_names <- paste0("sp_",1:ncol(A_int))
colnames(A_int) <- species_names
rownames(A_int) <- species_names


# Exclusion prob from Anisofun
prob_excl_AnisoFun_symmetric <- prob_extinction_4_int_matrix(A_int, 1e4, 1e4)

# Exclusion prob from ODE
number_random_points <- 5e4

number_species <- ncol(A_int)

set.seed(1234)
random_feasible_growth_rates <- foreach (i=1:round((number_random_points/100),0), .combine=rbind,  
                                                .packages = c("tidyverse","matlib",
                                                              "zipfR","pracma",
                                                              "CValternatives",
                                                              "anisoFun")) %dopar% {
                                                               generate_random_feasible_growth_rates(A_int, random_points=100)
                                                              }

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
# sp_1   sp_10   sp_11   sp_12   sp_13   sp_14   sp_15   sp_16   sp_17   sp_18   sp_19    sp_2   sp_20   sp_21   sp_22   sp_23   sp_24   sp_25 
# 0.02062 0.02012 0.02032 0.02004 0.01946 0.02056 0.01994 0.02068 0.02130 0.01964 0.02022 0.02066 0.02034 0.02022 0.01986 0.02082 0.02044 0.01998 
# sp_26   sp_27   sp_28   sp_29    sp_3   sp_30   sp_31   sp_32   sp_33   sp_34   sp_35   sp_36   sp_37   sp_38   sp_39    sp_4   sp_40   sp_41 
# 0.02010 0.01972 0.01966 0.02084 0.01886 0.01976 0.01844 0.01960 0.01926 0.02006 0.01958 0.02020 0.02120 0.02002 0.01930 0.02114 0.02000 0.01954 
# sp_42   sp_43   sp_44   sp_45   sp_46   sp_47   sp_48   sp_49    sp_5   sp_50    sp_6    sp_7    sp_8    sp_9 
# 0.01986 0.01948 0.02056 0.01972 0.01926 0.01982 0.02006 0.01968 0.02090 0.02054 0.01962 0.01922 0.01924 0.01954 

# prob_excl_AnisoFun_symmetric
# A tibble: 50 x 4
#      species prob_excl_mean prob_excl_lowerCI prob_excl_upperCI
#     <chr>            <dbl>             <dbl>             <dbl>
#   1 sp_1            0.0200            0.0200            0.0200
#   2 sp_2            0.0200            0.0200            0.0200
#   3 sp_3            0.0200            0.0200            0.0200
#   4 sp_4            0.0200            0.0200            0.0200
#   5 sp_5            0.0200            0.0200            0.0200
#   6 sp_6            0.0200            0.0200            0.0200
#   7 sp_7            0.0200            0.0200            0.0200
#   8 sp_8            0.0200            0.0200            0.0200
#   9 sp_9            0.0200            0.0200            0.0200
#  10 sp_10           0.0200            0.0200            0.0200


write_csv(random_feasible_growth_rates,"Results/simulations/random_feasible_growth_rates_symmetric.csv")
write_csv(as.data.frame(first_sp_excluded_results_symmetric),"Results/simulations/first_sp_excluded_ODE_results_symmetric.csv")
write_csv(prob_excl_AnisoFun_symmetric,"Results/simulations/prob_excl_AnisoFun_symmetric.csv")
################################################################################
################################################################################
# Define asymmetric interaction matrix--------------------------------------------

total_number_species = 15
set.seed(1234)
A_int2 <- -1*matrix(runif(n=total_number_species^2, min=0, max=1), nrow=total_number_species)
A_int2[A_int2>-.95] <- 0
set.seed(1234)
diag(A_int2)<- runif(total_number_species, min=-1, max=-.1)
set.seed(1234)
A_int2 <- A_int2 -1*matrix(runif(n=total_number_species^2, min=0, max=0.01), nrow=total_number_species)

incenter_inradius_isoprob_A_int2 <- anisoFun::incenter_inradius_isoprob_calculation(A_int2)
incenter_A_int2 <- incenter_inradius_isoprob_A_int2[[1]]
anisoFun::closest_to_exclusion(incenter_A_int2,A_int2)
anisoFun::sp_distance_to_exclusion(incenter_A_int2,A_int2)
vertices_A_int2 <- anisoFun::vertices_unit_ball(A_int2)
vertices_A_int2 <- vertices_A_int2[,(total_number_species+1):(2*total_number_species)] %>% as.matrix()
vertices_A_int2 %*% incenter_A_int2

species_names <- paste0("sp_",1:ncol(A_int2))
colnames(A_int2) <- species_names
rownames(A_int2) <- species_names
anisoFun::Omega_bootstrap(A_int2, replicates = 1e3)


# Exclusion prob from Anisofun
prob_excl_AnisoFun_asymmetric <- prob_extinction_4_int_matrix(A_int2, 1e4, 1e4)
anisotropy_index_asymmetric <- anisoFun::anisotropy_index(prob_excl_AnisoFun_asymmetric$prob_excl_mean)
prob_excl_AnisoFun_asymmetric
anisotropy_index_asymmetric

write_csv(prob_excl_AnisoFun_asymmetric,"Results/simulations/prob_excl_AnisoFun_asymmetric.csv")

# Exclusion prob from ODE
number_random_points_asymm <- 2e4
batch_number_asymm <- 1
number_species2 <- ncol(A_int2)

set.seed(1234)
random_feasible_growth_rates_asymm <- foreach (i=1:round((number_random_points_asymm/batch_number_asymm),0),
                                               i_counter_asym=icount(), .combine=rbind,
                                               .packages = c("tidyverse","matlib",
                                                       "zipfR","pracma",
                                                       "CValternatives",
                                                       "anisoFun")) %dopar% {
                                                         write_csv(as.data.frame(i_counter_asym), "Data/simulations/iterator_random_feasible_growth_rates_asymmetric.csv")
                                                         generate_random_feasible_growth_rates(A_int2, random_points=batch_number_asymm)
                                                       }



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


write_csv(random_feasible_growth_rates_asymm,"Results/simulations/random_feasible_growth_rates_asymmetric.csv")
write_csv(as.data.frame(first_sp_excluded_results_asymmetric),"Results/simulations/first_sp_excluded_ODE_results_asymmetric.csv")


################################################################################
################################################################################
# Define block interaction matrix--------------------------------------------
A_int3 <- A_int2

i_block <- round(0.5*ncol(A_int3)+1,0)

for(i in 1:i_block){
  for(j in (i_block+1):nrow(A_int3)){
    
    A_int3[i,j] <- 0
    A_int3[j,i] <- 0
    
  }
}


species_names <- paste0("sp_",1:ncol(A_int2))
colnames(A_int3) <- species_names
rownames(A_int3) <- species_names

inv(A_int3)

anisoFun::Omega_bootstrap(A_int3, replicates = 1e3)


# Exclusion prob from Anisofun

prob_excl_AnisoFun_block <- prob_extinction_4_int_matrix(A_int3, 1e4, 1e4)
anisotropy_index_block <- anisoFun::anisotropy_index(prob_excl_AnisoFun_block$prob_excl_mean)
prob_excl_AnisoFun_block
anisotropy_index_block

write_csv(prob_excl_AnisoFun_block,"Results/simulations/prob_excl_AnisoFun_block.csv")

# Exclusion prob from ODE
number_random_points_block <- 2e4
batch_number_block <- 1
number_species3 <- ncol(A_int3)

set.seed(1234)
random_feasible_growth_rates_block <- foreach (i=1:round((number_random_points_block/batch_number_block),0), 
                                               i_counter_block=icount(), .combine=rbind,  
                                               .packages = c("tidyverse","matlib",
                                                             "zipfR","pracma",
                                                             "CValternatives",
                                                             "anisoFun")) %dopar% {
                                                               write_csv(as.data.frame(i_counter_block),
                                                                         "Data/simulations/iterator_random_feasible_growth_rates_block.csv")
                                                               generate_random_feasible_growth_rates(A_int3, random_points=batch_number_block)
                                                             }



first_sp_excluded_results_block <- foreach (i=1:nrow(random_feasible_growth_rates_block), .combine=c,  
                                                 .packages = c("tidyverse","matlib",
                                                               "zipfR","pracma",
                                                               "CValternatives",
                                                               "anisoFun","deSolve")) %dopar% {
                                                                 growth_rates_vector_block <- random_feasible_growth_rates_block[i,1:number_species3] %>% as.numeric() %>% unname()
                                                                 initial_abundances_block <- random_feasible_growth_rates_block[i,(number_species3+1):(2*number_species3)] %>% as.numeric() %>% unname()
                                                                 species_excluded_i_block <- first_species_excluded_ODE(A_int3, growth_rates_vector_block,initial_abundances_block, delta_t = 0.1, tmax = 1000)
                                                                 
                                                               }


prob_excl_ODE_block <- table(first_sp_excluded_results_block)/number_random_points_block %>% sort()
prob_excl_ODE_block


write_csv(random_feasible_growth_rates_block,"Results/simulations/random_feasible_growth_rates_block.csv")
write_csv(as.data.frame(first_sp_excluded_results_block),"Results/simulations/first_sp_excluded_results_block.csv")

################################################################################
################################################################################
# Define diagonal interaction matrix--------------------------------------------
A_int4 <- A_int2

for(i in 2:nrow(A_int4)){
  for(j in 1:(i-1)){
    
    A_int4[i,j] <- 0
    
  }
}

species_names <- paste0("sp_",1:ncol(A_int4))
colnames(A_int4) <- species_names
rownames(A_int4) <- species_names

inv(A_int4)

anisoFun::Omega_bootstrap(A_int4, replicates = 1e3)




# Exclusion prob from Anisofun
prob_excl_AnisoFun_upper <- prob_extinction_4_int_matrix(A_int4, 1e4, 1e4)
anisotropy_index_upper <- anisoFun::anisotropy_index(prob_excl_AnisoFun_upper$prob_excl_mean)
prob_excl_AnisoFun_upper
anisotropy_index_upper

write_csv(prob_excl_AnisoFun_upper,"Results/simulations/prob_excl_AnisoFun_upper.csv")

# Exclusion prob from ODE
number_random_points_upper <- 2e4
batch_number_upper <- 1
number_species4 <- ncol(A_int4)

set.seed(1234)
random_feasible_growth_rates_upper <- foreach (i=1:round((number_random_points_upper/batch_number_upper),0),
                                               i_counter_upper=icount(), .combine=rbind,  
                                               .packages = c("tidyverse","matlib",
                                                             "zipfR","pracma",
                                                             "CValternatives",
                                                             "anisoFun")) %dopar% {
                                                               write_csv(as.data.frame(i_counter_upper),
                                                                         "Data/simulations/iterator_random_feasible_growth_rates_upper.csv")
                                                               generate_random_feasible_growth_rates(A_int4, random_points=batch_number_upper)
                                                             }



first_sp_excluded_results_upper <- foreach (i=1:nrow(random_feasible_growth_rates_upper), .combine=c,  
                                            .packages = c("tidyverse","matlib",
                                                          "zipfR","pracma",
                                                          "CValternatives",
                                                          "anisoFun","deSolve")) %dopar% {
                                                            growth_rates_vector_upper <- random_feasible_growth_rates_upper[i,1:number_species4] %>% as.numeric() %>% unname()
                                                            initial_abundances_upper <- random_feasible_growth_rates_upper[i,(number_species4+1):(2*number_species4)] %>% as.numeric() %>% unname()
                                                            species_excluded_i_upper <- first_species_excluded_ODE(A_int4, growth_rates_vector_upper,initial_abundances_upper, delta_t = 0.1, tmax = 1000)
                                                            
                                                          }


prob_excl_ODE_upper <- table(first_sp_excluded_results_upper)/number_random_points_block
prob_excl_ODE_upper


write_csv(random_feasible_growth_rates_upper,"Results/simulations/random_feasible_growth_rates_upper.csv")
write_csv(as.data.frame(first_sp_excluded_results_upper),"Results/simulations/first_sp_excluded_results_upper.csv")




anisotropy_index_symmetric <- anisoFun::anisotropy_index(prob_excl_AnisoFun_symmetric$prob_excl_mean)


stopImplicitCluster()
