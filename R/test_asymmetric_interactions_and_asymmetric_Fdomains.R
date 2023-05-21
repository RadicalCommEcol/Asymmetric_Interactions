library(tidyverse)
# library(gamlss.dist)
# library(moments)
library(anisoFun)
source("R/isotropy_index.R")
source("R/isotropy_metrics_4_int_matrix.R")
source("R/matrix_year_i.R")
################################################################################
################################################################################

# We create the metaweb for each year

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


years_included <- matrix_entries_raw$year %>% unique()

data_species <- NULL
data_kurtosis <- NULL
data_skewness <- NULL

for(year_i in years_included){
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  elements_A <- as.vector(A_int) 
  
  data_species <- c(data_species,ncol(A_int))
  data_kurtosis <- c(data_kurtosis,moments::kurtosis(elements_A))
  data_skewness <- c(data_skewness,moments::skewness(elements_A))
  
  cat("richness: ",ncol(A_int),", kurtosis: ",moments::kurtosis(elements_A),", skewness: ",moments::skewness(elements_A), "\n")
  
}
mean(data_species[data_species>8])

################################################################################
################################################################################
# Init parallelization
workers <- 8
cl <- makeCluster(workers)
registerDoParallel(cl) # register the cluster for using foreach

################################################################################
################################################################################


community_richness <- c(3,6,10)
kurtosis_index <- c(0.2,1,20)
number_repetitions <- 201

number_Omega_replicates <-  55000
number_boot_replicates <- number_Omega_replicates

set.seed(1234)

results <- NULL

for( S in community_richness){
  
  for(k in kurtosis_index ) {
    
    for(i in 1:number_repetitions){
      
      
      points_SHASH_i <- gamlss.dist::rSHASH(S*S,mu=0,sigma=1,nu=k,tau=k)
      
      # while (sd(points_SHASH_i)<0.4) {
      #   points_SHASH_i <- gamlss.dist::rSHASH(S*S,mu=0,sigma=1,nu=k,tau=k)
      # }
      
      A_int <- matrix(points_SHASH_i, ncol = S)
      
      
      repetition_i_partial_info <- tibble(richness=S,nu=k,tau=k,kurtosis=moments::kurtosis(points_SHASH_i),repetition=i)
      isotropy_metrics_i <- anisotropy_metrics_4_int_matrix(A_int, number_Omega_replicates, number_boot_replicates)
      repetition_i_full_info <- bind_cols(repetition_i_partial_info, isotropy_metrics_i)
      
      results <- bind_rows(results, repetition_i_full_info)
      
    }
    
    # write_csv(results, paste0("Data/test_asymmetry_interactions_FD_",number_repetitions,".csv"))
    
  }
  
  
}


results <- read_csv("Data/test_asymmetry_interactions_FD_201_omega_rep_55000.csv")

facet_names <- c(
  `3` = "3 species",
  `6` = "6 species",
  `10` = "10 species"
)

results_plot <- results %>% mutate(big_omega=small_omega_mean^richness)
results_plot$type <- "Platykurtic"
results_plot$type[results_plot$nu==1] <- "Mesokurtic"
results_plot$type[results_plot$nu<1] <- "Leptokurtic"

df <- group_by(results_plot, richness, type)
df_isotropy <- summarise(df, avg = median(anisotropy_index_mean))
df_isotropy

df_capital_omega <- summarise(df, avg = median(big_omega))
df_capital_omega

ggplot(results_plot, aes(x=big_omega, fill=as.factor(type))) +
  geom_histogram(binwidth=.025, alpha=.5, position="identity")+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  geom_vline(data=df_capital_omega, aes(xintercept=avg,  colour=as.factor(type)),
             linetype="dashed", size=1)+
  theme_bw()+
  labs(x="FD size", y="Frequency")+
  theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )

ggplot(results_plot, aes(x=anisotropy_index_mean, fill=as.factor(type))) +
  geom_histogram(binwidth=.025, alpha=.5, position="identity")+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  geom_vline(data=df_isotropy, aes(xintercept=avg,  colour=as.factor(type)),
             linetype="dashed", size=1)+
  labs(x="Asymmetry index, J'", y="Frequency")+
  theme_bw()+
  theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )



number_points <- 10000
points_SHASH_normal <- gamlss.dist::rSHASH(number_points,mu=0,sigma=1,nu=1,tau=1)#rSHASH(10, mu = 0, sigma = 1, nu = 0.0, tau = 0.5)
hist(points_SHASH_normal, 50)

points_SHASH_platy <- gamlss.dist::rSHASH(number_points,mu=0,sigma=1,nu=20,tau=20)#rSHASH(10, mu = 0, sigma = 1, nu = 0.0, tau = 0.5)
hist(points_SHASH_platy, 50)

points_SHASH_lepto <- gamlss.dist::rSHASH(number_points,mu=0,sigma=1,nu=.2,tau=.2)#rSHASH(10, mu = 0, sigma = 1, nu = 0.0, tau = 0.5)
hist(points_SHASH_lepto, 50)

distr_SHASH_normal <- tibble(value=points_SHASH_normal,type="Mesokurtic")
distr_SHASH_platy <- tibble(value=points_SHASH_platy,type="Platykurtic")
distr_SHASH_lepto <- tibble(value=points_SHASH_lepto,type="Leptokurtic")

distr_results <- bind_rows(distr_SHASH_platy,distr_SHASH_normal,distr_SHASH_lepto)



ggplot(distr_results, aes(x=value, fill=as.factor(type))) +
  geom_density(alpha=0.5)+
  #facet_wrap(type ~ .,scales = "free_x")+
  xlim(-5,5)+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()+
  theme(
    axis.text.x= element_text(color="black",size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )



ggplot(results_plot, aes(x=anisotropy_index_mean, fill=as.factor(type))) +
  geom_density(alpha=.3)+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  geom_vline(data=df_isotropy, aes(xintercept=avg,  colour=as.factor(type)),
             linetype="dashed", size=1)+
  theme_bw()+
  theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )


ggplot(results_plot, aes(y=anisotropy_index_mean, x = kurtosis, color=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  #geom_smooth()+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  theme_bw()

png("Images/asymmetry_interactions_domain.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.45, units = "in", res=300*2)
ggplot(results_plot, aes(x=as.factor(type), y = anisotropy_index_mean,color=as.factor(type))) +
  geom_boxplot(alpha=.3,size=1)+
  geom_hline(aes(yintercept=1),  color="deepskyblue1",
             linetype="dashed", size=1)+
  labs(y="Asymmetry index, J'(A)",x=NULL)+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+ 
  theme(legend.position="none") +
  theme(legend.title=element_blank())
dev.off()

library(scales)
png("Images/asymmetry_test_NACHO.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.45, units = "in", res=300*2)
ggplot(results_plot, aes(x=kurtosis, y = anisotropy_index_mean,color=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  geom_smooth(method = "lm",formula = y~x)+
  geom_hline(aes(yintercept=1),  color="deepskyblue1",
             linetype="dashed", size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(y="Asymmetry index, J'(A)",x=NULL)+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+ 
  theme(legend.title=element_blank())
dev.off()


library(latex2exp)
library(scales)

png("Images/asymmetry_interactions_domain_FD_sizes.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.45, units = "in", res=300*2)
ggplot(results_plot, aes(x=as.factor(type), y = big_omega,color=as.factor(type))) +
  geom_boxplot(alpha=.3,size=1)+
  # geom_hline(aes(yintercept=0.5),  color="deepskyblue1",
             # linetype="dashed", size=1)+
  ylab(TeX("Feasibility domain's size:  $\\Omega$"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  # labs(y="Feasibility domain's size",x=NULL)+
  xlab(NULL)+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+ 
  theme(legend.position="none") +
  theme(legend.title=element_blank())
dev.off()


library(scales)

ggplot(results_plot, aes(y=anisotropy_index_mean, x = big_omega, color=as.factor(type), fill=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth()+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(x="FD Size", y = "Asymmetry index, J'")+
  theme(
    axis.text.x= element_text(color="black",size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )

library("ggExtra")


asym_size_plot_16_sp <- ggplot(results_plot %>% filter(richness==16), 
                         aes(y=anisotropy_index_mean, x = small_omega_mean, color=as.factor(type), fill=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth()+
  # geom_smooth(method = "lm")+
  theme_bw()+
  labs(x="FD Size", y = "Asymmetry index, J'", title="16 species")+
  theme(
    axis.text.x= element_text(color="black",size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )

asym_size_plot_4_sp <- ggplot(results_plot %>% filter(richness==4), 
                               aes(y=anisotropy_index_mean, x = small_omega_mean, color=as.factor(type), fill=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth()+
  # geom_smooth(method = "lm")+
  theme_bw()+
  labs(x="FD Size", y = "Asymmetry index, J'", title="4 species")+
  theme(
    axis.text.x= element_text(color="black",size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )


asym_size_plot_8_sp <- ggplot(results_plot %>% filter(richness==8), 
                               aes(y=anisotropy_index_mean, x = small_omega_mean, color=as.factor(type), fill=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth()+
  # geom_smooth(method = "lm")+
  theme_bw()+
  labs(x="FD Size", y = "Asymmetry index, J'", title="8 species")+
  theme(
    axis.text.x= element_text(color="black",size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )



ggMarginal_sp4 <- ggMarginal(asym_size_plot_4_sp, groupColour = TRUE, groupFill = TRUE, alpha= 0.5)
ggMarginal_sp8 <- ggMarginal(asym_size_plot_8_sp, groupColour = TRUE, groupFill = TRUE, alpha= 0.5)
ggMarginal_sp16 <- ggMarginal(asym_size_plot_16_sp, groupColour = TRUE, groupFill = TRUE, alpha= 0.5)

library(patchwork)

png("Images/asymmetry_interactions_all.png",
    width = 11.69*0.8, # The width of the plot in inches
    height = 11.69*0.45, units = "in", res=300*2)
patchwork::wrap_elements(ggMarginal_sp4) +  patchwork::wrap_elements(ggMarginal_sp8) + patchwork::wrap_elements(ggMarginal_sp16 )
dev.off()

library(ggpubr)
ggarrange(ggMarginal_sp4, ggMarginal_sp8, ggMarginal_sp16, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

ggplot(results_plot, aes(y=anisotropy_index_mean, x = big_omega, color=as.factor(type), fill=as.factor(type))) +
  geom_point(alpha=.3,size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #geom_smooth()+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  # geom_smooth(method = "lm",formula = y~x+I(x^2))+
  theme_bw()+
  labs(x="FD Size", y = "Asymmetry index, J'")+
  theme(
    axis.text.x= element_text(color="black",size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1),
    panel.background = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
    legend.position="bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA)
  )

