
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(deSolve)
library(foreach)
library(doParallel)
library(latex2exp)
library(patchwork)

source("R/generate_random_feasible_growth_rates.R")
source("R/log10_population_growth_evolution_ODE.R")

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

border_growth_rates <- random_feasible_growth_rates %>% filter(N_1 < 5e-3 | N_2 < 5e-3 | N_3 < 5e-3)


# Plot FD 
ggplot(random_feasible_growth_rates, aes(x=r_1,y=r_2))+
  geom_point(alpha=0.3)+
  geom_point(data = border_growth_rates, aes(x=r_1,y=r_2), color="red")

ggplot(random_feasible_growth_rates, aes(x=r_1,y=r_3))+
  geom_point(alpha=0.3)+
  geom_point(data = border_growth_rates, aes(x=r_1,y=r_3), color="red")

ggplot(random_feasible_growth_rates, aes(x=r_2,y=r_3))+
  geom_point(alpha=0.3)+
  geom_point(data = border_growth_rates, aes(x=r_2,y=r_3), color="red")


plot_N1 <- ggplot(random_feasible_growth_rates, aes(x=r_2,y=r_3, color = N_1))+
  geom_point()+
  ylab(TeX("$r_2$"))+
  xlab(TeX("$r_3$"))+
  theme_bw()+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  scale_color_continuous(TeX("$N_1^*$"))

plot_N2 <- ggplot(random_feasible_growth_rates, aes(x=r_2,y=r_3, color = N_2))+
  geom_point()+
  ylab(TeX("$r_2$"))+
  xlab(TeX("$r_3$"))+
  theme_bw()+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  scale_color_continuous(TeX("$N_2^*$"))

plot_N3 <- ggplot(random_feasible_growth_rates, aes(x=r_2,y=r_3, color = N_3))+
  geom_point()+
  ylab(TeX("$r_2$"))+
  xlab(TeX("$r_3$"))+
  theme_bw()+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  scale_color_continuous(TeX("$N_3^*$"))



png("Images/population_equilibrium_unitball.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*350/1300, units = "in", res=300*2)

(plot_N1 | plot_N2 | plot_N3)
dev.off()

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

# Extract growth rates for each step
incenter_inradius_isoprob_data <- anisoFun::incenter_inradius_isoprob_calculation(A_int)
incenter <- incenter_inradius_isoprob_data[[1]]
inradius <- incenter_inradius_isoprob_data[[2]]


population_incenter <- (-1)*Inverse(A_int) %*% matrix(incenter,ncol=1) %>% as.numeric()



#########################################################
# TEST ONE STEP SEQUENCE FROM THE INCENTER TO THE BORDER
#########################################################

border_growth_rates <- random_feasible_growth_rates %>% filter(N_1 < 5e-3 | N_2 < 5e-3 | N_3 < 5e-3)

ggplot(random_feasible_growth_rates, aes(x=r_1,y=r_2))+
  geom_point(alpha=0.3)+
  geom_point(data = border_growth_rates, aes(x=r_1,y=r_2), color="red")

ggplot(random_feasible_growth_rates, aes(x=r_1,y=r_3))+
  geom_point(alpha=0.3)+
  geom_point(data = border_growth_rates, aes(x=r_1,y=r_3), color="red")


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

png("Images/results_from_incenter_to_border.png",
    width = 11.69*.6, # The width of the plot in inches
    height = 11.69*.55, units = "in", res=300*2)

plot_from_incenter_to_border

dev.off()

# Interquartilic ranges
IQR(results_one_step$log10_pop_growth[results_one_step$species=="sp1"])
IQR(results_one_step$log10_pop_growth[results_one_step$species=="sp2"])
IQR(results_one_step$log10_pop_growth[results_one_step$species=="sp3"])
prob_excl_AnisoFun_asymmetric %>% mutate(Excl_ratio = number_species * prob_excl_mean)
