
# Estimation of the abundance at equilibrium of in-silico communities with
# three species.

# To estimate initial growth rates within the feasibility domain we use the 
# following auxiliary function:
# "Scripts/aux_functions/generate_random_feasible_growth_rates.R"


# INPUT: 
# interaction matrix of the in-silico community: A_int

# OUTPUT: The statistics presented in Suppl. Section 8 for species abundances
# at equilibrium.

#-----------------------------------------------------------------------


library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)

source("Scripts/aux_functions/generate_random_feasible_growth_rates.R")

# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach




# Define asymmetric interaction matrix--------------------------------------------

number_species <- 3

A_int <- matrix(rep(0,number_species*number_species), ncol=number_species)
A_int[1,] <- -c(1.308894,2.05202504,0.03489887)
A_int[2,] <- -c(1.507218,1.10978936,0.24928687)
A_int[3,] <- -c(0.883729,0.03056066,1.13730985)


# > A_int
# [,1]       [,2]         [,3]
# [1,] -1.1304281 -1.0383207 -0.009622536
# [2,] -0.7649447 -1.9958739 -1.022003557
# [3,] -0.1235011 -0.6196368 -1.467263740

#A_int <- matrix(-abs(rnorm(number_species*number_species)), ncol=number_species)
A_int

species_names <- paste0("sp_",1:ncol(A_int))
colnames(A_int) <- species_names
rownames(A_int) <- species_names

prob_excl_AnisoFun_asymmetric <- anisoFun::prob_extinction_4_int_matrix(-1*A_int,5000)
prob_excl_AnisoFun_asymmetric %>% mutate(Excl_ratio = number_species * prob_excl_mean)

# Define study parameters
number_random_points_within_FD <- 5e4
number_random_pairs_within_FD <- 1e3

# Generate points within the FD
set.seed(123)
random_feasible_growth_rates <- generate_random_feasible_growth_rates(A_int, number_random_points_within_FD)

mean(random_feasible_growth_rates$N_1)
sd(random_feasible_growth_rates$N_1)
mean(random_feasible_growth_rates$r_1)
sd(random_feasible_growth_rates$r_1)
mean(border_growth_rates$r_1[border_growth_rates$N_1 < 5e-3])
sd(border_growth_rates$r_1[border_growth_rates$N_1 < 5e-3])
mean(border_growth_rates$r_1[border_growth_rates$N_2 < 5e-3])
sd(border_growth_rates$r_1[border_growth_rates$N_2 < 5e-3])
mean(border_growth_rates$r_1[border_growth_rates$N_3 < 5e-3])
sd(border_growth_rates$r_1[border_growth_rates$N_3 < 5e-3])
# Average abundance at Equilibrium for sp1: 0.1665575 +- 0.111594
# Average growth rate within the FD fot sp1: 0.5905199 +- 0.1668483
# ER_1 = 1.46

mean(random_feasible_growth_rates$N_2)
sd(random_feasible_growth_rates$N_2)
mean(random_feasible_growth_rates$r_2)
sd(random_feasible_growth_rates$r_2)
mean(border_growth_rates$r_2)
sd(border_growth_rates$r_2)
mean(border_growth_rates$r_2[border_growth_rates$N_2 < 5e-3])
sd(border_growth_rates$r_2[border_growth_rates$N_2 < 5e-3])
mean(border_growth_rates$r_2[border_growth_rates$N_1 < 5e-3])
sd(border_growth_rates$r_2[border_growth_rates$N_1 < 5e-3])
mean(border_growth_rates$r_2[border_growth_rates$N_3 < 5e-3])
sd(border_growth_rates$r_2[border_growth_rates$N_3 < 5e-3])
# Average abundance at Equilibrium for sp2: 0.1757607 +- 0.1094629
# Average growth rate within the FD fot sp1: 0.5307306 +- 0.08103975
# ER_2 = 1.01

mean(random_feasible_growth_rates$N_3)
sd(random_feasible_growth_rates$N_3)
mean(random_feasible_growth_rates$r_3)
sd(random_feasible_growth_rates$r_3)
mean(border_growth_rates$r_3)
sd(border_growth_rates$r_3)
mean(border_growth_rates$r_3[border_growth_rates$N_3 < 5e-3])
sd(border_growth_rates$r_3[border_growth_rates$N_3 < 5e-3])
mean(border_growth_rates$r_3[border_growth_rates$N_1 < 5e-3])
sd(border_growth_rates$r_3[border_growth_rates$N_1 < 5e-3])
mean(border_growth_rates$r_3[border_growth_rates$N_2 < 5e-3])
sd(border_growth_rates$r_3[border_growth_rates$N_2 < 5e-3])
# Average abundance at Equilibrium for sp3: 0.3395075 +- 0.2166704
# Average growth rate within the FD fot sp1: 0.5386883 +- 0.2121838
# ER_3 = 0.525

