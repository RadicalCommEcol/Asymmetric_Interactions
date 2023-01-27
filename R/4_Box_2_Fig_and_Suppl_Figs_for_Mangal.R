

################################################################################
# Before running this script, you should run the following scripts:
# - 3a_Process_Prob_exclusion_Mangal_NO_mutualistic_tradeoff.R
# - 3b_Process_Prob_exclusion_Mangal_NO_interspec_competition.R
################################################################################


library(tidyverse)
library(anisoFun)
# Functions to run calculations about the isotropic area
source("R/interaction_strength_matrix_from_edge_list.R")

number_Omega_replicates <- 10000
p_excl_data_aux_NO_mutualistic_tradeoff <- read_csv(
          paste0("Data/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_p_excl_data_aux_rep_",
                 number_Omega_replicates,".csv")) %>% mutate(type = "no_mutua_tradeoff")

p_excl_data_aux_NO_interspec_competition <- read_csv(
          paste0("Data/mangal_processed_data/NO_NAs_NO_interspec_competition_p_excl_data_aux_rep_",
                 number_Omega_replicates,".csv")) %>% mutate(type = "no_intrasp_competition")

isoprob_data <-p_excl_data_aux_NO_mutualistic_tradeoff %>% select(network_id, isoprob) %>% unique()

p_excl_data_aux <- rbind(p_excl_data_aux_NO_mutualistic_tradeoff,
                         p_excl_data_aux_NO_interspec_competition)

number_animals_plants <- p_excl_data_aux %>% filter(repetition == 1) %>%
  group_by(network_id,kingdom) %>% count()  %>% spread(kingdom, n)


p_excl_data_aux <- p_excl_data_aux %>% 
  left_join(number_animals_plants, by = c("network_id")) 

p_excl_data_aux <-  mutate(p_excl_data_aux, 
                           specialization = if_else(kingdom == "Animalia",
                                                    (Plantae-degree)/(Plantae-1),
                                                    (Animalia-degree)/(Animalia-1)))


mean_p_excl_data <- p_excl_data_aux %>% 
  group_by(network_id,species,type,kingdom, metric) %>% summarise_all(mean)

# Fig. 5 (box 2)

mean_p_excl_data_5019 <- mean_p_excl_data  %>% filter(network_id=="mangal_5019",
                                                      type != "no_intrasp_competition")


mean_p_excl_data_5019_animal <- mean_p_excl_data_5019  %>% filter(kingdom!="Plantae")
mean_p_excl_data_5019_plantae <- mean_p_excl_data_5019  %>% filter(kingdom=="Plantae")

cor.test(mean_p_excl_data_5019$value,mean_p_excl_data_5019$degree)

mean(mean_p_excl_data_5019_animal$value)
sd(mean_p_excl_data_5019_animal$value)

mean(mean_p_excl_data_5019_plantae$value)
sd(mean_p_excl_data_5019_plantae$value)

cor.test(mean_p_excl_data$value,mean_p_excl_data$degree, method = "spearman")

png("Images/mangal_example_no_intrasp_competition.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)

ggplot(data = mean_p_excl_data  %>% filter(network_id=="mangal_5019",
                                           type != "no_intrasp_competition"),
       aes(x = reorder(as.factor(species), -value, na.rm = T), y = value))+
  #geom_boxplot()+
  # stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.3)+
  stat_summary(fun = mean, geom="point", shape=23, 
               aes(size=degree, fill = as.factor(kingdom)))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  #geom_hline(aes(yintercept = isoprob), 
  #           color = "deepskyblue", linetype = "dashed", size = 1.2, alpha =.5)+
  labs(  y="Probability of exclusion",x="Species", fill = "Taxa", size = "Degree")+
  theme_bw()+theme(axis.text.x=element_blank())+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  theme(legend.text = element_text(size=18),legend.title=element_text(size=20))+ 
  guides(fill = guide_legend(override.aes = list(size=5)))

dev.off()

png("Images/mangal_no_intrasp_competition.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)

ggplot(data = mean_p_excl_data  %>% filter(type == "no_intrasp_competition"),
       aes(x = reorder(as.factor(species), -value, na.rm = T), y = value))+
  #geom_boxplot()+
  #stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.3)+
  stat_summary(fun = mean, geom="point", shape=23, 
               aes(size=degree, fill = as.factor(kingdom)))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  facet_wrap(~network_id, ncol = 8, scales = "free")+
  # geom_hline(data = isoprob_data, aes(yintercept = isoprob), 
  #            color = "deepskyblue", linetype = "dashed", size = 1.2, alpha =.5)+
  labs(y="Probability of exclusion",x="Species", fill = "Taxa", size = "Degree")+
  theme_bw()+theme(axis.text.x=element_blank(),legend.position="bottom")+
  guides(fill = guide_legend(override.aes = list(size=5)))

dev.off()

png("Images/mangal_no_mutualistic_tradeoff.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)

ggplot(data = mean_p_excl_data  %>% filter(type != "no_intrasp_competition"),
       aes(x = reorder(as.factor(species), -value, na.rm = T), y = value))+
  #geom_boxplot()+
  # stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.3)+
  stat_summary(fun = mean, geom="point", shape=23, 
               aes(size=degree, fill = as.factor(kingdom)))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  facet_wrap(~network_id, ncol = 8, scales = "free")+
  # geom_hline(data = isoprob_data, aes(yintercept = isoprob), 
  #            color = "deepskyblue", linetype = "dashed", size = 1.2, alpha =.5)+
  labs(  y="Probability of exclusion",x="Species", fill = "Taxa", size = "Degree")+
  theme_bw()+theme(axis.text.x=element_blank(),legend.position="bottom")+
  guides(fill = guide_legend(override.aes = list(size=5)))

  

dev.off()

