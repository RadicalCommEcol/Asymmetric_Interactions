
library(tidyverse)
library(ggpubr)

total_number_species <- 3
data_results_3 <- read_csv(paste0("Results/results_shape_indices_comparisson_D",total_number_species,".csv"))
cos_ref_3 <- data_results_3$cos_variance[data_results_3$proportions==1]
J_ref_3 <- data_results_3$J_index[data_results_3$proportions==1]

data_saavedra_3 <- read_csv(paste0("Results/results_saavedra_shape_indices_comparisson_D",total_number_species,".csv"))
saavedra_ref_3 <- data_saavedra_3$saavedra_index[data_saavedra_3$proportions==1]

data_results_3_adi <- data_results_3 %>% mutate(cos_adi=cos_variance/cos_ref_3,
                                                  J_adi = J_index / J_ref_3,
                                                saavedra_adi = data_saavedra_3$saavedra_index/saavedra_ref_3,
                                                  total_number_species=total_number_species)


total_number_species <- 6
data_results_6 <- read_csv(paste0("Results/results_shape_indices_comparisson_D",total_number_species,".csv"))
cos_ref_6 <- data_results_6$cos_variance[data_results_6$proportions==1]
J_ref_6 <- data_results_6$J_index[data_results_6$proportions==1]

data_saavedra_6 <- read_csv(paste0("Results/results_saavedra_shape_indices_comparisson_D",total_number_species,".csv"))
saavedra_ref_6 <- data_saavedra_6$saavedra_index[data_saavedra_6$proportions==1]


data_results_6_adi <- data_results_6 %>% mutate(cos_adi=cos_variance/cos_ref_6,
                                                  J_adi = J_index / J_ref_6,
                                                saavedra_adi = data_saavedra_6$saavedra_index/saavedra_ref_6,
                                                  total_number_species=total_number_species)

total_number_species <- 10
data_results_10 <- read_csv(paste0("Results/results_shape_indices_comparisson_D",total_number_species,".csv")) %>%
  filter(proportions<=1)
cos_ref_10 <- data_results_10$cos_variance[data_results_10$proportions==1]
J_ref_10 <- data_results_10$J_index[data_results_10$proportions==1]

data_saavedra_10 <- read_csv(paste0("Results/results_saavedra_shape_indices_comparisson_D",total_number_species,".csv"))
saavedra_ref_10 <- data_saavedra_10$saavedra_index[data_saavedra_10$proportions==1]


data_results_10_adi <- data_results_10 %>% mutate(cos_adi=cos_variance/cos_ref_10,
                                                  J_adi = J_index / J_ref_10,
                                                  saavedra_adi = data_saavedra_10$saavedra_index/saavedra_ref_10,
                                                  total_number_species=total_number_species)




coeff <- 1 # Value used to transform the data


data_plot_3sp <- data_results_3_adi %>% filter(proportions<=1,convergence_solver=="Function criterion near zero") %>%
  select("proportions","cos_adi","J_adi","saavedra_adi") %>%
  gather(key = "metric",value = "value",-"proportions")

data_plot_3sp$metric[data_plot_3sp$metric=="cos_adi"] <- "Variance of the cosine of the side lengths"
data_plot_3sp$metric[data_plot_3sp$metric=="J_adi"] <- "Asymmetry index J'(A)"
data_plot_3sp$metric[data_plot_3sp$metric=="saavedra_adi"] <- "Variance of the niche differences"

results_3 <- ggplot(data = data_plot_3sp,
       aes(x=1-proportions, y = value,col=metric,linetype=metric)) +
  geom_line(lwd = 1.5)+
  scale_linetype_manual(values = c(rep("solid",1),
                                   rep("dotted",1),
                                   rep("dashed",1)))+
  scale_color_manual(values = c(rep("#E69F00",1),
                                rep("#56B4E9",1),
                                rep("#009E73",1)))+
  labs(title="3 species",
       x="Relative decrease in\nfeasibility domain size",
       y="Index relative change")+
  ylim(c(0,1.25))+
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.text=element_text(size=14,face="bold"),
    legend.title=element_blank(),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.key.width = unit(1.7,"cm")
  )+
  guides(linetype=guide_legend(nrow=3,byrow=TRUE))


data_plot_6sp <- data_results_6_adi %>% filter(proportions<=1,convergence_solver=="Function criterion near zero") %>%
  select("proportions","cos_adi","J_adi","saavedra_adi") %>%
  gather(key = "metric",value = "value",-"proportions")

data_plot_6sp$metric[data_plot_6sp$metric=="cos_adi"] <- "Variance of the cosine of the side lengths"
data_plot_6sp$metric[data_plot_6sp$metric=="J_adi"] <- "Asymmetry index J'(A)"
data_plot_6sp$metric[data_plot_6sp$metric=="saavedra_adi"] <- "Variance of the niche differences"

results_6 <- ggplot(data = data_plot_6sp,
                    aes(x=1-proportions, y = value,col=metric,linetype=metric)) +
  geom_line(lwd = 1.5)+
  scale_linetype_manual(values = c(rep("solid",1),
                                   rep("dotted",1),
                                   rep("dashed",1)))+
  scale_color_manual(values = c(rep("#E69F00",1),
                                rep("#56B4E9",1),
                                rep("#009E73",1)))+
  labs(title="6 species",
       x="Relative decrease in\nfeasibility domain size",
       y="Index relative change")+
  ylim(c(0,1.25))+
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.text=element_text(size=14,face="bold"),
    legend.title=element_blank(),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.key.width = unit(1.7,"cm")
  )+
  guides(linetype=guide_legend(nrow=3,byrow=TRUE))



data_plot_10sp <- data_results_10_adi %>% filter(proportions<=1,convergence_solver=="Function criterion near zero") %>%
  select("proportions","cos_adi","J_adi","saavedra_adi") %>%
  gather(key = "metric",value = "value",-"proportions")

data_plot_10sp$metric[data_plot_10sp$metric=="cos_adi"] <- "Variance of the cosine of the side lengths"
data_plot_10sp$metric[data_plot_10sp$metric=="J_adi"] <- "Asymmetry index J'(A)"
data_plot_10sp$metric[data_plot_10sp$metric=="saavedra_adi"] <- "Variance of the niche differences"

results_10 <- ggplot(data = data_plot_10sp,
                    aes(x=1-proportions, y = value,col=metric,linetype=metric)) +
  geom_line(lwd = 1.5)+
  scale_linetype_manual(values = c(rep("solid",1),
                                   rep("dotted",1),
                                   rep("dashed",1)))+
  scale_color_manual(values = c(rep("#E69F00",1),
                                rep("#56B4E9",1),
                                rep("#009E73",1)))+
  labs(title="10 species",
       x="Relative decrease in\nfeasibility domain size",
       y="Index relative change")+
  ylim(c(0,1.25))+
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.text=element_text(size=14,face="bold"),
    legend.title=element_blank(),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.key.width = unit(1.7,"cm")
  )+
  guides(linetype=guide_legend(nrow=3,byrow=TRUE))



png("Images/fig3.png",
    width = 11.69*1, # The width of the plot in inches
    height = 11.69*.5, units = "in", res=300*2)
ggarrange(results_3, results_6, results_10, ncol=3, nrow=1, common.legend = TRUE, legend="top")
dev.off()


