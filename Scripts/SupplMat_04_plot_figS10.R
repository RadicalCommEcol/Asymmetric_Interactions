
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(latex2exp)
library(patchwork)

source("Scripts/aux_functions/generate_random_feasible_growth_rates.R")
source("Scripts/aux_functions/log10_population_growth_evolution_ODE.R")

# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

# Define asymmetric interaction matrix to conduct our test----------------------

number_species <- 3

A_int <- matrix(rep(0,number_species*number_species), ncol=number_species)
A_int[1,] <- -c(1.308894,2.05202504,0.03489887)
A_int[2,] <- -c(1.507218,1.10978936,0.24928687)
A_int[3,] <- -c(0.883729,0.03056066,1.13730985)


species_names <- paste0("sp_",1:ncol(A_int))
colnames(A_int) <- species_names
rownames(A_int) <- species_names

prob_excl_AnisoFun_asymmetric <- anisoFun::prob_extinction_4_int_matrix(-1*A_int,5000)
prob_excl_AnisoFun_asymmetric %>% mutate(Excl_ratio = number_species * prob_excl_mean)

# Define study parameters
number_random_points_within_FD <- 5e4

# Generate points within the FD
set.seed(123)
random_feasible_growth_rates <- generate_random_feasible_growth_rates(A_int, number_random_points_within_FD)

#############################################################################
# TEST ONE STEP SEQUENCE FROM THE INCENTER TO THE BORDER: SHORT-TERM CHANGES
#############################################################################

border_growth_rates <- random_feasible_growth_rates %>% filter(N_1 < 5e-3 | N_2 < 5e-3 | N_3 < 5e-3)

# Extract growth rates for each step
incenter_inradius_isoprob_data <- anisoFun::incenter_inradius_isoprob_calculation(A_int)
incenter <- incenter_inradius_isoprob_data[[1]]
inradius <- incenter_inradius_isoprob_data[[2]]

population_incenter <- (-1)*Inverse(A_int) %*% matrix(incenter,ncol=1) %>% as.numeric()

population_growth_evolution_one_step <- foreach (i=1:nrow(border_growth_rates), .combine=rbind,  
                                                 .packages = c("tidyverse","matlib",
                                                               "zipfR","pracma",
                                                               "CValternatives",
                                                               "anisoFun","deSolve")) %dopar% {
                                                                 
                                                                 initial_growth_rate <- incenter
                                                                 final_growth_rate <- border_growth_rates[i,1:number_species] %>% 
                                                                   as.numeric() %>% unname()
                                                                 initial_abundances <- population_incenter
                                                                 species_excluded_i_asymm <- log10_population_growth_evolution_ODE(A_int, initial_growth_rate,final_growth_rate,
                                                                                                                                   initial_abundances, delta_t = .05, tmax = 500)
                                                                 
                                                               }


results_one_step <- tibble(time = c(rep(population_growth_evolution_one_step$time,3)) %>%
                             as.numeric(),
                           
                           log10_pop_growth = c(population_growth_evolution_one_step$log10_pop_growth_PLUS1_1,
                                                population_growth_evolution_one_step$log10_pop_growth_PLUS1_2,
                                                population_growth_evolution_one_step$log10_pop_growth_PLUS1_3) %>%
                             as.numeric(),
                           species = c(rep("sp1",nrow(population_growth_evolution_one_step)),
                                       rep("sp2",nrow(population_growth_evolution_one_step)),
                                       rep("sp3",nrow(population_growth_evolution_one_step))),
                           exclusion_ratio = c(rep(number_species*prob_excl_AnisoFun_asymmetric$prob_excl_lowerCI[1],
                                                   nrow(population_growth_evolution_one_step)),
                                               rep(number_species*prob_excl_AnisoFun_asymmetric$prob_excl_lowerCI[2],
                                                   nrow(population_growth_evolution_one_step)),
                                               rep(number_species*prob_excl_AnisoFun_asymmetric$prob_excl_lowerCI[3],
                                                   nrow(population_growth_evolution_one_step)))
                           
)

plot_from_incenter_to_border <- ggplot(data = results_one_step, aes(x=round(exclusion_ratio,3),y=log10_pop_growth,group=round(exclusion_ratio,3)))+
  geom_boxplot(outlier.alpha = 0.1)+
  #geom_jitter(alpha=0.2)+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(TeX("Exclusion ratio:  $ER_i$"))+
  theme_bw()+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

png("Images/figS10.png",
    width = 11.69*.6, # The width of the plot in inches
    height = 11.69*.55, units = "in", res=300*2)

plot_from_incenter_to_border

dev.off()