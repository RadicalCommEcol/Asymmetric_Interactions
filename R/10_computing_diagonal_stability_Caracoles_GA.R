
library(matlib) # to multiply matrices
library(tidyverse)

source("R/matrix_year_i.R")


# We create the metaweb for each year

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


years_included <- matrix_entries_raw$year %>% unique()

for(year_i in years_included){
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  
  W <- t(A_int) + A_int
  
  ev_W <- eigen(W)
  # extract components
  eig_values_W <- ev_W$values
  
  cat(all(eig_values_W<0), "\n")
  
}

# Apparently our interaction matrices are not Volterra-dissipative under the approach in 
# Grilli 2016, when using the identity matrix

# Lets find a diagonal matrix that meet d-stability criteria



fn4 <- function(x){
  
  D=diag(abs(x))
  W <- t(A_int) %*% D + D %*% A_int
  
  ev_W <- eigen(W)
  # extract components
  eig_values_W <- ev_W$values
  
  eig_values_greater_than_zero <- as.numeric(eig_values_W>0)
  
  
  #return(sum(-1*(eig_values_W*eig_values_greater_than_zero)))
  return(-1*max(eig_values_W))
  
}

eval_solution <- function(x){
  
  D=diag(abs(x))
  W <- t(A_int) %*% D + D %*% A_int
  
  ev_W <- eigen(W)
  # extract components
  eig_values_W <- ev_W$values
  
  eig_values_greater_than_zero <- as.numeric(eig_values_W>0)
  
  
  return(eig_values_W)
  #return(exp(eig_values_greater_than_zero)/sum(exp(eig_values_greater_than_zero))-1/ncol(D))
  
}

library(GA)
library(doParallel)
library(foreach)
library(iterators)

set.seed(123)

for(year_i in 2016:2021){ #years_included){
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  
  set.seed(1234)
  # GA <- GA::ga(type = "real-valued", fitness = fn4, lower = rep(0,ncol(A_int)), upper = rep(10,ncol(A_int)), 
  #              popSize = 1000, maxiter = 150000, run = 10000)
  GA <- GA::ga(type = "real-valued", fitness = fn4, lower = rep(0,ncol(A_int)), upper = rep(10,ncol(A_int)), 
               popSize = 1000, maxiter = 150000, run = 10000,seed = 12345, parallel = 4)
  
  plot(GA)
  summary(GA)
  
  eval_solution( GA@solution %>% as.numeric())
  
  write_csv(GA@solution %>% as.numeric(),paste0("Data/diagonal_matrix_year_", year_i,".csv"))
  
}
