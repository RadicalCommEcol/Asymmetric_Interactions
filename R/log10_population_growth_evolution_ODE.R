
# load necessary packages

library(mvtnorm)
library(mgcv)
if(!require(ggtern)) {install.packages("ggtern"); library(ggtern)}
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}



log10_population_growth_evolution_ODE <- function(A_int, initial_growth_rate,final_growth_rate,
                                                  initial_abundances, delta_t = 0.1, tmax = 1000) {
  
  number_species <- ncol(A_int)
  species_names <- colnames(A_int)
  
  alpha <- A_int
  dimnames(alpha) <- NULL
  
  time_step <- seq(0,tmax,by=delta_t)
  r <- final_growth_rate
  N0 <- initial_abundances
  parms <- list(r=r, alpha = alpha) #ODE
  model <- function(t,N,parms){ dN <- N * (parms$r + parms$alpha %*% N); list(dN)}
  sol <- ode(N0,time_step,model,parms, method = 'bdf', rtol = 1e-15, atol = 1e-15, maxsteps = 10000)# ode(N0,time_step,model,parms)
  # plot(sol)
  
  colnames(sol) <- c("time","sp1","sp2","sp3")
  
  sol_transient <- sol %>% as_tibble() %>%
    mutate(sp1=round(sp1,2),sp2=round(sp2,2),sp3=round(sp3,2)) %>% 
    filter(sp1 != lag(sp1) & sp2 != lag(sp2) & sp3 != lag(sp3)) %>%
    mutate(lag_sp1 = lag(sp1),lag_sp2 = lag(sp2),lag_sp3 = lag(sp3),
           pop_growth1 = sp1/lag_sp1, pop_growth2 = sp2/lag_sp2, pop_growth3 = sp3/lag_sp3,
           log10_pop_growth_PLUS1_1 = log10(pop_growth1+1),
           log10_pop_growth_PLUS1_2 = log10(pop_growth2+1),
           log10_pop_growth_PLUS1_3 = log10(pop_growth3+1)) 
  
  result <- slice(sol_transient, -1) %>% dplyr::select(time,log10_pop_growth_PLUS1_1,
                                                       log10_pop_growth_PLUS1_2,
                                                       log10_pop_growth_PLUS1_3)
  
  return(result)
}



full_log10_population_growth_evolution_ODE <- function(A_int, initial_growth_rate,final_growth_rate,
                                                  initial_abundances, delta_t = 0.1, tmax = 1000) {
  
  number_species <- ncol(A_int)
  species_names <- colnames(A_int)
  
  alpha <- A_int
  dimnames(alpha) <- NULL
  
  time_step <- seq(0,tmax,by=delta_t)
  r <- final_growth_rate
  N0 <- initial_abundances
  parms <- list(r=r, alpha = alpha) #ODE
  model <- function(t,N,parms){ dN <- N * (parms$r + parms$alpha %*% N); list(dN)}
  sol <- ode(N0,time_step,model,parms, method = 'bdf', rtol = 1e-15, atol = 1e-15, maxsteps = 10000)# ode(N0,time_step,model,parms)
  # plot(sol)
  
  colnames(sol) <- c("time","sp1","sp2","sp3")
  
  sol_transient <- sol %>% as_tibble() %>%
    mutate(lag_sp1 = lag(sp1),lag_sp2 = lag(sp2),lag_sp3 = lag(sp3),
           pop_growth1 = sp1/lag_sp1, pop_growth2 = sp2/lag_sp2, pop_growth3 = sp3/lag_sp3,
           log10_pop_growth_PLUS1_1 = log10(pop_growth1+1),
           log10_pop_growth_PLUS1_2 = log10(pop_growth2+1),
           log10_pop_growth_PLUS1_3 = log10(pop_growth3+1)) 
  
  result <- slice(sol_transient, -1) %>% dplyr::select(time,log10_pop_growth_PLUS1_1,
                                                       log10_pop_growth_PLUS1_2,
                                                       log10_pop_growth_PLUS1_3)
  
  return(result)
}

