
library(matlib) # to multiply matrices
library(tidyverse)
library(anisoFun)
library(latex2exp)
library(patchwork)

source("Scripts/aux_functions/generate_random_feasible_growth_rates.R")




# Define asymmetric interaction matrix--------------------------------------------

number_species <- 3

A_int <- matrix(rep(0,number_species*number_species), ncol=number_species)
A_int[1,] <- -c(1.308894,2.05202504,0.03489887)
A_int[2,] <- -c(1.507218,1.10978936,0.24928687)
A_int[3,] <- -c(0.883729,0.03056066,1.13730985)


# Define study parameters
number_random_points_within_FD <- 5e4
number_random_pairs_within_FD <- 1e3

# Generate points within the FD
set.seed(123)
random_feasible_growth_rates <- generate_random_feasible_growth_rates(A_int, number_random_points_within_FD)


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



png("Images/figS9.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*350/1300, units = "in", res=300*2)

(plot_N1 | plot_N2 | plot_N3)
dev.off()
