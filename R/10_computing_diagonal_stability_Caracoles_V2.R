
library(matlib) # to multiply matrices
library(tidyverse)
library(rlist)

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


fn3 <- function(x){
  
  D=diag(abs(x))
  W <- t(A_int) %*% D + D %*% A_int
  
  ev_W <- eigen(W)
  # extract components
  eig_values_W <- ev_W$values
  
  eig_values_greater_than_zero <- as.numeric(eig_values_W>0)

  
  return(rep(max(eig_values_W),ncol(A_int))+100)
  #return(exp(eig_values_greater_than_zero)/sum(exp(eig_values_greater_than_zero))-1/ncol(D))
  
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


library("nleqslv")

set.seed(123)
diagonal_lst <- list()
for(year_i in years_included){
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  
  init_search <- abs(runif(ncol(A_int)))*1e-12
  
  x=init_search
  
  #The eigenvalues must be negative. We force the nleqslv algorithm to converge to search for eigenvalues close to
  #limit <- -abs(runif(ncol(A_int))) 
  
  
  search_diagonal <- nleqslv(init_search, fn3,control=list(maxit=10000,allowSingular=T),method = c("Newton")) #nleqslv(init_search, fn2,control=list(maxit=10000))
  search_diagonal
  fn3(search_diagonal$x)-100
  eval_solution(search_diagonal$x)
  
  while(!(all(fn3(search_diagonal$x)-100<0))){
    # cat(fn2(search_diagonal$x), "inicial\n")
    init_search <- abs(runif(ncol(A_int)))*1e-12
    search_diagonal <- nleqslv(init_search, fn3,control=list(maxit=10000,allowSingular=T),method = c("Newton"),global = "hook")
    cat(fn3(search_diagonal$x)-100, "\n")
    cat(!(all(fn3(search_diagonal$x)-100<0)), "\n")
  }
  
  cat(eval_solution(search_diagonal$x), "\n")
  cat(all(fn3(search_diagonal$x)-100<0), "\n")
  
  diagonal_lst <- list.append(diagonal_lst,search_diagonal$x) 
  
}
