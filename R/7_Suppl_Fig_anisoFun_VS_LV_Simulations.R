

library(matlib) # to multiply matrices
library(tidyverse)
library(scales)
library(anisoFun)
library(patchwork)
library(cowplot)

MAE <- function(x,y){
  
  return(sum((x-y)^2)/len(x))
}

# Symmetric results

first_sp_excluded_results_symmetric <- read_csv("Data/simulations/first_sp_excluded_ODE_results_symmetric.csv")
samples_symmetric <-  nrow(first_sp_excluded_results_symmetric)
prob_excl_ODE_symmetric <- first_sp_excluded_results_symmetric %>% 
  rename(species=first_sp_excluded_results_symmetric) %>%
  group_by(species) %>% count() %>%
  mutate(prob_excl_ODE = n / samples_symmetric) %>% select(-n)

prob_excl_AnisoFun_symmetric <- 
  read_csv("Data/simulations/prob_excl_AnisoFun_symmetric.csv") %>%
  left_join(prob_excl_ODE_symmetric, by = "species")

number_species_symmetric <- nrow(prob_excl_AnisoFun_symmetric)

plot_symm <- ggplot(prob_excl_AnisoFun_symmetric, aes(x = prob_excl_ODE,y = prob_excl_mean))+
  geom_abline(slope = 1, intercept = 0,
              color = "#009E73", linetype = "dashed", 
              size = 1.2, alpha =.5)+
  geom_abline(slope = 0, intercept = 1/number_species_symmetric,
              color = "#E69F00", linetype = "dotted", 
              size = 1.2, alpha =.5)+
  geom_vline(xintercept = 1/number_species_symmetric,
              color = "#E69F00", linetype = "dotted", 
              size = 1.2, alpha =.5)+
  geom_point(size=2,alpha=0.5)+
  geom_errorbar(aes(ymin=prob_excl_lowerCI, ymax=prob_excl_upperCI), width=.0005)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlim(0.033,.1)+ylim(0.033,.1)+
  labs(x = NULL, y =NULL,
       title = paste0("Asymmetry index J'(A) = ",round(1,2),
                      "\nMean square error = ", 
                      round(MAE(x = prob_excl_AnisoFun_symmetric$prob_excl_ODE,
                          y = prob_excl_AnisoFun_symmetric$prob_excl_mean),5)))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 14, face = "bold"))
  

# Block matrix results

first_sp_excluded_results_block <- read_csv("Data/simulations/first_sp_excluded_results_block.csv")
samples_block <-  nrow(first_sp_excluded_results_block)
prob_excl_ODE_block <- first_sp_excluded_results_block %>% 
  rename(species=first_sp_excluded_results_block) %>%
  group_by(species) %>% count() %>%
  mutate(prob_excl_ODE = n / samples_block) %>% select(-n)

prob_excl_AnisoFun_block <- 
  read_csv("Data/simulations/prob_excl_AnisoFun_block.csv") %>%
  left_join(prob_excl_ODE_block, by = "species")

number_species_block <- nrow(prob_excl_AnisoFun_block)

plot_block <- ggplot(prob_excl_AnisoFun_block, aes(x = prob_excl_ODE,y = prob_excl_mean))+
  geom_abline(slope = 1, intercept = 0,
              color = "#009E73", linetype = "dashed", 
              size = 1.2, alpha =.5)+
  geom_point(size=2,alpha=0.5)+
  geom_errorbar(aes(ymin=prob_excl_lowerCI, ymax=prob_excl_upperCI), width=.0005)+
  labs(x = NULL, y ="Prob. of exclusion based on\nfeasibility domain's shape",
       title = paste0("Asymmetry index J'(A) = ",
                      round(anisoFun::anisotropy_index(prob_excl_AnisoFun_block$prob_excl_mean),2),
                      "\nMean square error = ", 
                      round(MAE(x = prob_excl_AnisoFun_block$prob_excl_ODE,
                                y = prob_excl_AnisoFun_block$prob_excl_mean),5)))+
  theme_bw()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 14, face = "bold"))


# Upper matrix results

first_sp_excluded_results_upper <- read_csv("Data/simulations/first_sp_excluded_results_upper.csv")
samples_upper <-  nrow(first_sp_excluded_results_upper)
prob_excl_ODE_upper <- first_sp_excluded_results_upper %>% 
  rename(species=first_sp_excluded_results_upper) %>%
  group_by(species) %>% count() %>%
  mutate(prob_excl_ODE = n / samples_upper) %>% select(-n)

prob_excl_AnisoFun_upper <- 
  read_csv("Data/simulations/prob_excl_AnisoFun_upper.csv") %>%
  left_join(prob_excl_ODE_upper, by = "species")

number_species_upper <- nrow(prob_excl_AnisoFun_upper)

plot_upper <- ggplot(prob_excl_AnisoFun_upper, aes(x = prob_excl_ODE,y = prob_excl_mean))+
  geom_abline(slope = 1, intercept = 0,
              color = "#009E73", linetype = "dashed", 
              size = 1.2, alpha =.5)+
  geom_point(size=2,alpha=0.5)+
  geom_errorbar(aes(ymin=prob_excl_lowerCI, ymax=prob_excl_upperCI), width=.0005)+
  labs(x = "Prob. of exclusion from LV estimation", y =NULL,
       title = paste0("Asymmetry index J'(A) = ",
                      round(anisoFun::anisotropy_index(prob_excl_AnisoFun_upper$prob_excl_mean),2),
                      "\nMean square error = ", 
                      round(MAE(x = prob_excl_AnisoFun_upper$prob_excl_ODE,
                                y = prob_excl_AnisoFun_upper$prob_excl_mean),5)))+
  theme_bw()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 14, face = "bold"))


png("Images/Anisofun_ODE_asymmetry.png",
    width = 11.69*1.1, # The width of the plot in inches
    height = 11.69*0.4, units = "in", res=300*2)
plot_block | plot_upper | plot_symm
dev.off()
