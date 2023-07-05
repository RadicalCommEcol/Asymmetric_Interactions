
# Extraction of the probabilities of exclusion for the 88 networks in mangal,
# by using a parameterization without mutualistic tradeoff (see Saavedra et al.
# 2016 for further details [doi: 10.1002/ece3.1930])

# To build the interaction matrix from raw data we use the following auxiliary 
# function: "Scripts/aux_functions/interaction_strength_matrix_from_edge_list.R"

# To estimate mutualistic strength of each network and parametererization,
# we use the function: "Scripts/aux_functions/mutualistic_strength_solver.R"

# INPUT: 
# Pre-processed species' probabilities of exclusion: 
# "Results/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_mangal_boot_networks_rep_100000.csv"
# nodes data: "Data/mangal_raw_data/mangal_pollination_quant_nodes.csv"
# links data: "Data/mangal_raw_data/mangal_pollination_quant_links.csv"
# other information: "Data/mangal_raw_data/mangal_pollination_quant_nets.csv"

# OUTPUT:
# species' probabilities of exclusion:
# "Results/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_p_excl_data_aux_rep_100000.csv"

# Note: Since the output is a 3GB file, it is not available in our public repository
# due to storage constraints.

#-----------------------------------------------------------------------


library(tidyverse)
library(anisoFun)
library(nleqslv)

# Functions to run calculations about the isotropic area
source("Scripts/aux_functions/interaction_strength_matrix_from_edge_list.R")

# Load link information
matrix_links_raw <- read_csv("Results/mangal_processed_data/NO_NAs_mangal_links_processed_kingdom.csv")
lat_long_data <- read_csv2("Data/mangal_raw_data/mangal_pollination_quant_nets.csv") %>%
  select(network_id,network_lat,network_lon)
networks_included <- matrix_links_raw$network_id %>% unique()

kingdom_data_aux1 <- read_csv("Results/mangal_processed_data/NO_NAs_mangal_links_processed_kingdom.csv") %>%
  select(network_id,taxonomy_name_from,kingdom_from) %>%  
  group_by(network_id,taxonomy_name_from,kingdom_from) %>% count() %>%
  rename(species = taxonomy_name_from, kingdom = kingdom_from, degree = n)

kingdom_data_aux2 <- read_csv("Results/mangal_processed_data/NO_NAs_mangal_links_processed_kingdom.csv") %>%
  select(network_id,taxonomy_name_to,kingdom_to) %>%  
  group_by(network_id,taxonomy_name_to,kingdom_to) %>% count() %>%
  rename(species = taxonomy_name_to, kingdom = kingdom_to, degree = n)

kingdom_data <- rbind(kingdom_data_aux1,kingdom_data_aux2) %>% 
  arrange(network_id,kingdom)

network_aux1 <- matrix_links_raw %>% select(network_id,taxonomy_name_from) %>% 
  unique() %>% arrange(network_id,taxonomy_name_from) %>% rename(node=taxonomy_name_from)
network_aux2 <- matrix_links_raw %>% select(network_id,taxonomy_name_to) %>% 
  unique() %>% arrange(network_id,taxonomy_name_to) %>% rename(node=taxonomy_name_to)

networks_sizes <- rbind(network_aux1,network_aux2) %>% group_by(network_id) %>% count()

#####################################
# Extract maximum mean mutualistic strength

# Load link information
max_mutualistic_strength <- tibble(network_id = networks_included,
                                   gamma_0_max = as.numeric(NA))


for(network_i in networks_included){
  
  cat(network_i, "\n",
      which(networks_included == network_i), "\n")
  
  matrix_edge_list_i <- matrix_links_raw %>% filter(network_id == network_i) %>%
    dplyr::select(taxonomy_name_from,taxonomy_name_to,kingdom_from,kingdom_to) %>%
    rename(node_from = taxonomy_name_from, node_to = taxonomy_name_to)
  
  rho = 0.01
  delta = 0
  
  source("Scripts/aux_functions/mutualistic_strength_solver.R")
  
  xstart <- c(1)
  mutualistic_strength_results <- nleqslv(xstart,mutualistic_strength_solver)
  
  gamma_0_max <- mutualistic_strength_results$x
  
  max_mutualistic_strength[which(networks_included==network_i),2] <- gamma_0_max
  
}


#####################################

number_Omega_replicates <-  10000
number_boot_replicates <- number_Omega_replicates
mangal_plot_network_raw <- read_csv(paste0("Results/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_mangal_boot_networks_rep_",
                                           number_Omega_replicates,".csv"))

number_plant_sp_network <- mangal_plot_network_raw %>% group_by(network_id) %>% count() %>%
  rename(number_plant_sp = n)

mangal_plot_network <- mangal_plot_network_raw %>% 
  left_join(number_plant_sp_network, by = c("network_id")) %>%
  mutate(iso_prob_exclusion = 1/number_plant_sp)

name_colums <- c("network_id","number_plant_sp","iso_prob_exclusion","species", 
                 paste0("prob_excl_rep_",1:number_boot_replicates),
                 paste0("Omega_rep_",1:number_boot_replicates))

mangal_plot_network <- mangal_plot_network[,name_colums]

diff_prob_mat <- as.matrix(mangal_plot_network[,(5:(4+number_boot_replicates))])
colnames(diff_prob_mat) <- paste0("diff_isoprob_excl_",1:number_boot_replicates)

for(i in 1:nrow(diff_prob_mat)){
  diff_prob_mat[i,] <- diff_prob_mat[i,] - mangal_plot_network$iso_prob_exclusion[i]
}

mangal_plot_network_final_aux <- bind_cols(mangal_plot_network, as_tibble(diff_prob_mat))
mangal_plot_network_final <- mangal_plot_network_final_aux %>% 
  gather("metric", "value", -network_id,-number_plant_sp,-iso_prob_exclusion,-species) %>%
  mutate(repetition = as.numeric(gsub("[^0-9.-]", "", metric)),
         metric = gsub('[[:digit:]]+', '', metric),
         metric = substring(metric,1, nchar(metric)-1))

mangal_plot_network_final$metric %>% unique()

diff_isoprob_excl_data <- mangal_plot_network_final %>% filter(metric == "diff_isoprob_excl")%>% 
  arrange(network_id, repetition, value) %>% 
  left_join(kingdom_data,by = c("network_id","species")) %>% 
  left_join(max_mutualistic_strength, by = "network_id") %>%
  mutate(species = paste0(network_id,species), isoprob = 1/ number_plant_sp,
         mutualistic_strength = 0.5*gamma_0_max/(degree^0.0),
         total_mutualistic_strength = degree*0.5*gamma_0_max/(degree^0.0))
  


p_excl_data_aux <- mangal_plot_network_final %>% filter(metric == "prob_excl_rep") %>% 
  arrange(network_id, repetition, value) %>% 
  left_join(kingdom_data,by = c("network_id","species")) %>% 
  left_join(max_mutualistic_strength, by = "network_id") %>%
  mutate(species = paste0(network_id,species), isoprob = 1/ number_plant_sp,
       mutualistic_strength = 0.5*gamma_0_max/(degree^0.0),
       total_mutualistic_strength = degree*0.5*gamma_0_max/(degree^0.0))

isoprob_data <- p_excl_data_aux %>% select(network_id, isoprob) %>% unique()

write_csv(p_excl_data_aux,
          paste0("Results/mangal_processed_data/NO_NAs_NO_mutualistic_tradeoff_p_excl_data_aux_rep_",
                 number_Omega_replicates,".csv"))