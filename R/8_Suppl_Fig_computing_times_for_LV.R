
library(matlib) # to multiply matrices
library(tidyverse)
library(scales)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(iterators)

source("R/generate_random_feasible_growth_rates.R")
source("R/first_species_excluded_ODE.R")

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
  
  time_LV_list <- c(time_LV_list, (end_time - start_time)) # Cuidado esto no dice si son segundos, minutos o días 
  
}

time_LV_list2 <- c(0.03600192,0.07001400, 0.26200223,1.83304191,9.59698892, 60*1.09100039, 60*7.73155475,60*48.01390030,3600*4.78676046,3600*24*1.46608374)

data_experiment <- tibble(total_number_species = list_total_number_species,
                          capital_omega = capital_omega_list,
                          time_LV_unitary = time_LV_list2/number_random_points)

# write_csv(data_experiment,"Data/simulations/data_experiment_LV_times.csv") commented for security reasons

data_experiment <- read_csv("Data/simulations/data_experiment_LV_times.csv")

png("Images/symmetic_FD_computing_times.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)

ggplot(data_experiment, aes(x=capital_omega,y=time_LV_unitary))+
  geom_smooth(method = "lm", color = "deepskyblue", fill = "deepskyblue")+
  geom_point(size=2)+
  theme_bw()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Feasibility domain's relative volume [adimensional]",
       y="Average time to find an intrinsic\ngrowth vector in the feasibility domain [seconds]")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

dev.off()
#############

model <- lm(log10(time_LV_unitary)~log10(capital_omega),data_experiment)
summary(model)

stopImplicitCluster()
