library(matlib) # to multiply matrices
library(tidyverse)
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(MultitrophicFun)
library(foreach)
library(doParallel)
library(boot) # to bootstrap
library(CValternatives) # to estimate PV index
library(ggpmisc)

# Functions to run the isotropic-area calculations
# source("R/vertices_unit_ball.R")
# source("R/cone_vertices_director_vertices.R")
# source("R/incenter_inradius_isoprob_calculation.R")
# source("R/incenter_tangent_points_solver.R")
# source("R/tangent_points_spherical_triangle.R")
# source("R/furthest_from_extinction.R")
# source("R/closest_to_extinction.R")
# source("R/Omega_song_et_al_par_bootstrap.R")
# source("R/exclusion_probabilities_par_bootstrap.R")
# source("R/boot_prob_excl_Omega_raw.R")
source("R/isotropy_index.R")
source("R/isotropy_metrics_4_int_matrix.R")
# source("R/prob_extinction_4_int_matrix.R")
# source("R/small_omega_iso.R")

library(anisoFun)

numCores <- detectCores()

registerDoParallel(numCores-4)

# Create an interaction matrix
number_species <- 3

isosceles_proportions <- read.csv("Images/isosceles_proportions.csv")

A_int <- matrix(rep(0,9), ncol = 3)

cos_variance <- NULL
L2_norm_sd <- NULL
J_index <- NULL
J_index_lower <- NULL
J_index_upper <- NULL
prob_excl_var <- NULL

for (i in 1:nrow(isosceles_proportions)) {
  
  A_int[,1]=-c(isosceles_proportions[i,1], isosceles_proportions[i,2], isosceles_proportions[i,3])
  A_int[,2]=-c(isosceles_proportions[i,4], isosceles_proportions[i,5], isosceles_proportions[i,6])
  A_int[,3]=-c(isosceles_proportions[i,7], isosceles_proportions[i,8], isosceles_proportions[i,9])
  
  colnames(A_int) <- paste0("sp_",1:number_species)
  rownames(A_int) <- paste0("sp_",1:number_species)
  
  cos_variance <- c(cos_variance, var(c(sum(A_int[,1]*A_int[,2]),sum(A_int[,1]*A_int[,3]),sum(A_int[,2]*A_int[,3]))))
  L2_norm_sd <- c(L2_norm_sd, sd(c(sqrt(sum(A_int[,1]*A_int[,1])),sqrt(sum(A_int[,2]*A_int[,2])),sqrt(sum(A_int[,3]*A_int[,3])))))
  
  # Estimation of the incenter position and of the inradius for all the three vertices (of the 
  # simplicial cone) combinations
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
  
  I <- incenter_inradius_isoprob[[1]]
  
  set.seed(1234)
  
  number_Omega_replicates <- 1000
  number_boot_replicates <- number_Omega_replicates
  
  prob_excl_var_result <- anisoFun::prob_extinction_4_int_matrix(A_int) %>% select(prob_excl_mean) %>% pull() %>% var()
  
  isotropy_metrics_results <- isotropy_metrics_4_int_matrix(A_int, number_Omega_replicates, number_boot_replicates)
  
  J_index <- c(J_index,isotropy_metrics_results$isotropy_index_mean)
  J_index_lower <- c(J_index_lower,isotropy_metrics_results$isotropy_index_lowerCI)
  J_index_upper <- c(J_index_upper,isotropy_metrics_results$isotropy_index_upperCI)
  prob_excl_var <- c(prob_excl_var, prob_excl_var_result)
}


png("Images/cos_side_lengths_VS_J_index.png",
    width = 11.69*.4, # The width of the plot in inches
    height = 11.69*.3, units = "in", res=300*2)

ggplot(data = tibble(cos_variance = cos_variance,
                     J_index = J_index,
                     J_index_lower = J_index_lower,
                     J_index_upper = J_index_upper,
                     prob_excl_var = prob_excl_var,
                     proportions = isosceles_proportions$prop),
       aes(x = cos_variance, y = J_index,color=proportions))+
  geom_point()+
  geom_errorbar(aes(ymin=J_index_lower, ymax=J_index_upper), width=.005)+
  labs(x="Variance of the cosine of the side lengths", y= "J'(A)", color=NULL)+
  theme_bw()


dev.off()

# Value used to transform the data
coeff <- 1

# A few constants
cos_varianceeColor <- "#009E73"
J_indexColor <- "#E69F00"

ggplot(data = tibble(cos_variance = cos_variance,
                     J_index = J_index,
                     J_index_lower = J_index_lower,
                     J_index_upper = J_index_upper,
                     prob_excl_var = prob_excl_var,
                     proportions = isosceles_proportions$prop),
       aes(x=proportions)) +
  
  geom_line( aes(y=cos_variance), size=2, color=cos_varianceeColor) + 
  geom_line( aes(y=J_index), size=2, color=J_indexColor) +
  #geom_errorbar(aes(ymin=J_index_lower, ymax=J_index_upper), width=.1)+
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Variance of the cosine\nof the side lengths",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Asymmetry index J'(A)")
  ) + 
  scale_x_continuous(name="Ratio between the side lengths of an isosceles\nspherical triangle and \nthose of the reference triangle")+
  theme_bw() +
  
  theme(
    axis.title.y = element_text(color = cos_varianceeColor, size=17),
    axis.title.y.right = element_text(color = J_indexColor, size=17),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold")
  )



png("Images/cos_side_lengths_VS_J_index_2.png",
    width = 11.69*.5, # The width of the plot in inches
    height = 11.69*.6, units = "in", res=300*2)
ggplot(data = tibble(cos_variance = cos_variance,
                     J_index = J_index,
                     J_index_lower = J_index_lower,
                     J_index_upper = J_index_upper,
                     prob_excl_var = prob_excl_var,
                     proportions = isosceles_proportions$prop),
       aes(x=proportions)) +
  
  geom_line( aes(y=cos_variance,color="Variance of the cosine of the side lengths"), size=2) + 
  geom_line( aes(y=J_index,color="Asymmetry index J'(A)"), size=2) +
  #geom_errorbar(aes(ymin=J_index_lower, ymax=J_index_upper), width=.1)+
  scale_color_manual(name = NULL, values = c("Asymmetry index J'(A)" = "#E69F00", 
                                             "Variance of the cosine of the side lengths" = "#009E73"))+
labs(x="Ratio between the side lengths of an isosceles\nspherical triangle and those of the reference triangle",
     y=NULL)+
  theme_bw() +
  theme(
    legend.position="top",
    legend.text=element_text(size=14,face="bold"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold")
  )+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

dev.off()




ggplot(data = tibble(cos_variance = cos_variance,
                     J_index = J_index,
                     J_index_lower = J_index_lower,
                     J_index_upper = J_index_upper,
                     prob_excl_var = prob_excl_var,
                     proportions = isosceles_proportions$prop),
       aes(y = prob_excl_var, x = J_index,color=proportions))+
  geom_point()+
  geom_errorbar(aes(xmin=J_index_lower, xmax=J_index_upper), width=.0001)+
  labs(y="Variance of of exclusion prob.", x= "J'(A)", color=NULL)+
  theme_bw()



png("Images/variance_cos_side_lengths_VS_variance_excl_prob.png",
    width = 11.69*.4, # The width of the plot in inches
    height = 11.69*.3, units = "in", res=300*2)

ggplot(data = tibble(cos_variance = cos_variance,
                     J_index = J_index,
                     J_index_lower = J_index_lower,
                     J_index_upper = J_index_upper,
                     prob_excl_var = prob_excl_var,
                     proportions = isosceles_proportions$prop),
       aes(x = cos_variance, y = prob_excl_var, color=proportions))+
  stat_poly_line(method="lm", formula = y~x+I(x^.5)) +
  # stat_poly_eq(method="lm",formula = y~I(x^.5)+x,aes(label = paste(after_stat(eq.label),
  #                                after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point()+
  # geom_smooth(method="lm")+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(x="Variance of the cosine of the side lengths", y= "Variance of exclusion prob.", color=NULL)+
  theme_bw()

dev.off()

model <- lm(prob_excl_var~(cos_variance)+I(cos_variance^0.5), data = tibble(cos_variance = cos_variance,
                                             J_index = J_index,
                                             J_index_lower = J_index_lower,
                                             J_index_upper = J_index_upper,
                                             prob_excl_var = prob_excl_var,
                                             proportions = isosceles_proportions$prop))

summary(model)
