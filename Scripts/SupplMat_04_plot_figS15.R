
library(tidyverse)
library(corrplot)
############################################################
# Load species' arc distances

arc_distance <- read_csv("Data/caracoles_processed_data/caracoles_arc_distance.csv")
############################################################
# Load exclusion probabilities

number_Omega_replicates <-  10000
number_boot_replicates <- number_Omega_replicates
caracoles_plot_year_raw <- read_csv(paste0("Data/caracoles_processed_data/caracoles_boot_plants_year_rep_",
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

ggplot(abundances_year_pool,
       aes(x=as.factor(year), y=log10(total_individuals), fill=as.factor(species))) +
  geom_bar(stat="identity", position=position_dodge())

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
  mutate(year=as.factor(year),species = as.factor(species))  %>%
  mutate(dif_total_individuals = total_individuals - total_individuals_prev_year)


abundance_model_data$species %>% unique()


# Correlograms for plant species' exclusion metrics

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

final_pool_species <- abundance_model_data$species %>% unique() %>% levels()
data_exclusion_ratio <- spread(data = mean_prob_excl %>% filter(species %in% final_pool_species) %>% 
                           left_join(number_plant_sp_year, by = "year") %>%
                           mutate(exclusion_ratio = number_plant_sp*prob_excl), key=species, value=exclusion_ratio, fill = 0)

M_excl_ratio<- cor(data_exclusion_ratio %>% ungroup() %>% dplyr::select(-year,-prob_excl,-number_plant_sp))
# matrix of the p-value of the correlation
p.mat <- cor.mtest(data_exclusion_ratio %>% ungroup() %>% dplyr::select(-year,-prob_excl,-number_plant_sp))

png("Images/corr_plot_excl_ratio.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)
corrplot(M_excl_ratio, method="number",order="hclust",col=brewer.pal(n=10, name="RdYlBu"),p.mat = p.mat, sig.level = 0.05)

dev.off()


data_prob_excl <- spread(data = mean_prob_excl %>% filter(species %in% final_pool_species) %>% 
                                 left_join(number_plant_sp_year, by = "year"), key=species, value=prob_excl, fill = 0)

M_prob_excl<- cor(data_prob_excl %>% ungroup() %>% dplyr::select(-year,-number_plant_sp))
# matrix of the p-value of the correlation
p.mat <- cor.mtest(data_prob_excl %>% ungroup() %>% dplyr::select(-year,-number_plant_sp))

png("Images/corr_plot_prob_excl.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)
corrplot(M_prob_excl, method="number",order="hclust",col=brewer.pal(n=10, name="RdYlBu"),p.mat = p.mat, sig.level = 0.05)

dev.off()


data_arc_distance <- spread(data = arc_distance %>% filter(species %in% final_pool_species), key=species, value=arc_distance, fill = 0)

M_arc_distance<- cor(data_arc_distance %>% ungroup() %>% dplyr::select(-year,-relative_arc_distance,-incenter_arc_distace))
# matrix of the p-value of the correlation
p.mat.arc <- cor.mtest(data_arc_distance %>% ungroup() %>% dplyr::select(-year,-relative_arc_distance,-incenter_arc_distace))

png("Images/corr_plot_arc_distances.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)
corrplot(M_arc_distance, method="number",order="hclust",col=brewer.pal(n=10, name="RdYlBu"),p.mat = p.mat.arc, sig.level = 0.05)


dev.off()