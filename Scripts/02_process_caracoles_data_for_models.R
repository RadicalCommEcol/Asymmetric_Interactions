library(tidyverse)
############################################################
# Load species' arc distances

arc_distance <- read_csv("Results/caracoles_processed_data/caracoles_arc_distance.csv")
grilli_index <- read_csv("Results/caracoles_processed_data/caracoles_grilli.csv")
############################################################
# Load probabilities

number_Omega_replicates <-  10000
number_boot_replicates <- number_Omega_replicates
caracoles_plot_year_raw <- read_csv(paste0("Results/caracoles_processed_data/caracoles_boot_plants_year_rep_",
                                           number_Omega_replicates,".csv"))

number_plant_sp_year <- caracoles_plot_year_raw %>% group_by(year) %>% count() %>%
  rename(number_plant_sp = n)

caracoles_plot_year <- caracoles_plot_year_raw %>% 
  left_join(number_plant_sp_year, by = c("year")) %>%
  mutate(iso_prob_exclusion = 1/number_plant_sp)

name_colums <- c("year","number_plant_sp","iso_prob_exclusion","species", 
                 paste0("prob_excl_rep_",1:number_boot_replicates),
                 paste0("Omega_rep_",1:number_boot_replicates))

caracoles_plot_year <- caracoles_plot_year[,name_colums]

diff_prob_mat <- as.matrix(caracoles_plot_year[,(5:(4+number_boot_replicates))])
colnames(diff_prob_mat) <- paste0("diff_isoprob_excl_",1:number_boot_replicates)

for(i in 1:nrow(diff_prob_mat)){
  diff_prob_mat[i,] <- diff_prob_mat[i,] - caracoles_plot_year$iso_prob_exclusion[i]
}

caracoles_plot_year_final_aux <- bind_cols(caracoles_plot_year, as_tibble(diff_prob_mat))
caracoles_plot_year_final <- caracoles_plot_year_final_aux %>% 
  gather("metric", "value", -year,-number_plant_sp,-iso_prob_exclusion,-species) %>%
  mutate(repetition = as.numeric(gsub("[^0-9.-]", "", metric)),
         metric = gsub('[[:digit:]]+', '', metric),
         metric = substring(metric,1, nchar(metric)-1))

caracoles_plot_year_final$metric %>% unique()

mean_prob_excl <- caracoles_plot_year_final %>% filter(metric == "prob_excl_rep") %>% # "prob_excl_rep") %>% 
  select(year, species, value) %>%
  group_by(year, species) %>% summarise_all(mean) %>% rename(prob_excl = value)

mean_omega <- caracoles_plot_year_final %>% filter(metric == "Omega_rep") %>% # "prob_excl_rep") %>% 
  select(year, species, value) %>%
  group_by(year, species) %>% summarise_all(mean) %>% rename(Omega = value) %>%
  left_join(number_plant_sp_year, by = "year") %>%
  mutate(iso_prob = Omega^(1/number_plant_sp))


############################################################
# Load abundances

abundances_raw <- read_csv("Data/caracoles_raw_data/abundance.csv")
abundances_year <- abundances_raw %>% group_by(year,species) %>%
  count(wt = individuals) %>% filter(n>0) %>%
  rename(total_individuals = n)

full_pool_species <- abundances_raw %>% dplyr::select(species) %>% unique() %>% pull() %>% sort()

abundance_years <- abundances_raw %>% group_by(species, year) %>% count(wt =individuals) %>% spread(year, n)

pool_species <- abundances_year %>% dplyr::select(species,year) %>% distinct() %>% group_by(species) %>%
  count() %>% filter(n>4) %>% dplyr::select(species) %>% pull()

abundances_year_pool <- abundances_year %>% filter(species %in% c(pool_species))

abundances_year_pool %>% group_by(species) %>%
  count(wt = total_individuals)

# In 2018 there are some species that are virtually absent from Caracoles
abundances_year_pool %>% filter(total_individuals < 10)

############################################################
# Data for models

abundances_prev_year_data <- abundances_year_pool %>%
  mutate(year = year + 1, 
         total_individuals_prev_year = total_individuals) %>%
  dplyr::select(species,year,total_individuals_prev_year)

arc_distance_4_model <- arc_distance %>%
  mutate(year = year + 1, 
         arc_distance_prev_year = arc_distance,
         relative_arc_distance_prev_year = relative_arc_distance) %>%
  dplyr::select(species,year,arc_distance_prev_year,relative_arc_distance_prev_year)

mean_prob_excl_4_model <- mean_prob_excl %>%
  left_join(mean_omega, by = c("year","species")) %>%
  mutate(year = year + 1,
         exclusion_ratio_prev_year = prob_excl/iso_prob,
         prob_excl_prev_year = prob_excl)

abundance_model_data_aux <- abundances_year_pool %>% 
  left_join(abundances_prev_year_data, by = c("year","species")) %>%
  left_join(mean_prob_excl_4_model, by = c("year","species")) %>%
  left_join(arc_distance_4_model, by = c("year","species")) %>%
  filter(!year %in% c(2015,2018),species != "PUPA")


species_absent_previous_year <- abundance_model_data_aux$species[is.na(abundance_model_data_aux$prob_excl)] %>%
  unique()

abundance_model_data <- abundance_model_data_aux %>%
  filter(!species %in% species_absent_previous_year) %>%
  mutate(year=year,species = as.factor(species))  %>%
  mutate(dif_total_individuals = total_individuals - total_individuals_prev_year)


abundance_model_data$species %>% unique()

# Save data for the pool of 7 species
write_csv(abundance_model_data,"Results/abundance_model_data.csv")
