
# Statistical information for the 88 networks from mangal, included in
# the text Box 2.

# INPUT: 
# species' probabilities of exclusion when there is NO interspecific competition:
# "Results/mangal_processed_data/NO_NAs_NO_interspec_competition_p_excl_data_aux_rep_100000.csv"
# species' probabilities of exclusion when there is NO mutualistic tradeoff:
# "Results/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_p_excl_data_aux_rep_100000.csv"

################################################################################
# Before running this script, you should run the Box2_02_process_mangal scripts
# Note: Some of the inputs are a 3GB files and are not available in our public 
# repository due to storage constraints.
################################################################################

library(tidyverse)
library(anisoFun)


number_Omega_replicates <- 10000
p_excl_data_aux_NO_mutualistic_tradeoff <- read_csv(
          paste0("Results/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_p_excl_data_aux_rep_",
                 number_Omega_replicates,".csv")) %>% mutate(type = "no_mutua_tradeoff")

p_excl_data_aux_NO_interspec_competition <- read_csv(
          paste0("Results/mangal_processed_data/NO_NAs_NO_interspec_competition_p_excl_data_aux_rep_",
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


# Correlation between prob. of excl. and degree for 88 networks

cor.test(mean_p_excl_data$value,mean_p_excl_data$degree)
mean(mean_p_excl_data$value[mean_p_excl_data$kingdom=="Animalia"])
mean(mean_p_excl_data$value[mean_p_excl_data$kingdom=="Plantae"])
sd(mean_p_excl_data$value[mean_p_excl_data$kingdom=="Animalia"])
sd(mean_p_excl_data$value[mean_p_excl_data$kingdom=="Plantae"])

