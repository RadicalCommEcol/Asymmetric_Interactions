
library(tidyverse)
library(scales)
############################################################
# Load species' arc distances

arc_distance <- read_csv("Data/caracoles_processed_data/caracoles_arc_distance.csv")
############################################################
# Load probabilities

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

#######################################################
library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(DHARMa)
library(performance)
library(visreg)
library(RColorBrewer)
library(nlme)
library(usdm)

# Correlation among explanatory variables and response variables

cor.test(abundance_model_data$total_individuals,
         abundance_model_data$exclusion_ratio_prev_year)

cor.test(abundance_model_data$total_individuals,
         abundance_model_data$exclusion_ratio_prev_year, method = "spearman")

cor.test(abundance_model_data$total_individuals,
         abundance_model_data$arc_distance_prev_year)
cor.test(abundance_model_data$total_individuals,
         abundance_model_data$arc_distance_prev_year, method = "spearman")

cor.test(abundance_model_data$exclusion_ratio_prev_year,
         abundance_model_data$arc_distance_prev_year)
cor.test(abundance_model_data$exclusion_ratio_prev_year,
         abundance_model_data$arc_distance_prev_year, method = "spearman")

cor.test(abundance_model_data$total_individuals,
         abundance_model_data$total_individuals_prev_year)
cor.test(abundance_model_data$total_individuals,
         abundance_model_data$total_individuals_prev_year, method = "spearman")


vif(as.data.frame(dplyr::select(abundance_model_data %>% ungroup(),
                                total_individuals_prev_year,prob_excl)))
cor.test(abundance_model_data$total_individuals_prev_year,
         abundance_model_data$prob_excl_prev_year)

vif(as.data.frame(dplyr::select(abundance_model_data %>% ungroup(),
                                total_individuals_prev_year,exclusion_ratio_prev_year)))
cor.test(abundance_model_data$total_individuals_prev_year,
         abundance_model_data$exclusion_ratio_prev_year)

cor.test(abundance_model_data$total_individuals,
         abundance_model_data$exclusion_ratio_prev_year, method = "spearman")

vif(as.data.frame(dplyr::select(abundance_model_data %>% ungroup(),
                                prob_excl_prev_year,exclusion_ratio_prev_year)))
cor.test(abundance_model_data$prob_excl_prev_year,
         abundance_model_data$exclusion_ratio_prev_year)

vif(as.data.frame(dplyr::select(abundance_model_data %>% ungroup(),
                                total_individuals_prev_year,arc_distance_prev_year)))
cor.test(abundance_model_data$total_individuals_prev_year,abundance_model_data$arc_distance_prev_year)


vif(as.data.frame(dplyr::select(abundance_model_data %>% ungroup(),
                                prob_excl_prev_year,arc_distance_prev_year)))
cor.test(abundance_model_data$prob_excl_prev_year,abundance_model_data$arc_distance_prev_year)

vif(as.data.frame(dplyr::select(abundance_model_data %>% ungroup(),
                                total_individuals_prev_year,exclusion_ratio_prev_year,arc_distance_prev_year)))

# Caracoles models

abundance_model <- glmmTMB((total_individuals) ~ scale(total_individuals_prev_year)+
                             species,#(1|species),
                           #ziformula = ~1,
                           family = nbinom1(),
                           data = abundance_model_data)

abundance_excl_ratio_model <- glmmTMB((total_individuals) ~ scale(total_individuals_prev_year)+scale(exclusion_ratio_prev_year)+
                             species,#(1|species),
                           #ziformula = ~1,
                           family = nbinom1(),
                           data = abundance_model_data)

abundance_arc_distance_model <- glmmTMB((total_individuals) ~ scale(total_individuals_prev_year)+scale(arc_distance_prev_year)+
                                        species,#(1|species),
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = abundance_model_data)

abundance_prob_excl_model <- glmmTMB((total_individuals) ~ scale(total_individuals_prev_year)+scale(prob_excl_prev_year)+
                                          species,#(1|species),
                                        #ziformula = ~1,
                                        family = nbinom1(),
                                        data = abundance_model_data)

abundance_excl_ratio_arc_distance_model <- glmmTMB((total_individuals) ~ scale(total_individuals_prev_year)+scale(arc_distance_prev_year)+
                                                     scale(exclusion_ratio_prev_year)+species,#(1|species),
                                        #ziformula = ~1,
                                        family = nbinom1(),
                                        data = abundance_model_data)

abundance_prob_excl_arc_distance_model <- glmmTMB((total_individuals) ~ scale(total_individuals_prev_year)+scale(arc_distance_prev_year)+
                                                     scale(prob_excl_prev_year)+species,#(1|species),
                                                   #ziformula = ~1,
                                                   family = nbinom1(),
                                                   data = abundance_model_data)


##########

summary(abundance_model)
summary(abundance_excl_ratio_model)
summary(abundance_prob_excl_model)
summary(abundance_arc_distance_model)
summary(abundance_excl_ratio_arc_distance_model)
summary(abundance_prob_excl_arc_distance_model)


# AIC comparison

AIC(abundance_model,
    abundance_excl_ratio_model,
    abundance_prob_excl_model,
    abundance_arc_distance_model,
    abundance_excl_ratio_arc_distance_model,
    abundance_prob_excl_arc_distance_model)


# Best model: abundance_arc_distance_model

excl_distance_predictions <- tibble(real_data = abundance_model_data$total_individuals,
                                predicted_data = exp(predict(abundance_arc_distance_model, abundance_model_data)),
                                species=abundance_model_data$species,
                                year=abundance_model_data$year)


# Figure 4------------------------------------

species_names <- c(
  `BEMA` = "Beta\nmacrocarpa",
  `CETE` = "Centaurium\ntenuiflorum",
  `HOMA` = "Hordeum\nmarinum",
  `LEMA` = "Leontodon\nmaroccanus",
  `PAIN` = "Parapholis\nincurva",
  `POMA` = "Polypogon\nmaritimus",
  `SASO` = "Salsola\nsoda",
  `SCLA` = "Scorzonera\nlaciniata"
)

excl_distance_predictions_sp_names <- excl_distance_predictions

excl_distance_predictions_sp_names$species_name <- NA

excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "BEMA"] <- "Beta\nmacrocarpa"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "CETE"] <- "Centaurium\ntenuiflorum"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "HOMA"] <- "Hordeum\nmarinum"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "LEMA"] <- "Leontodon\nmaroccanus"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "PAIN"] <- "Parapholis\nincurva"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "POMA"] <- "Polypogon\nmaritimus"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "SASO"] <- "Salsola\nsoda"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "SCLA"] <- "Scorzonera\nlaciniata"
excl_distance_predictions_sp_names$species_name[excl_distance_predictions_sp_names$species == "PUPA"] <- "Pulicaria\npaludosa"

png("Images/caracoles_predictions_log.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.55, units = "in", res=300*2)

ggplot(excl_distance_predictions_sp_names %>% filter(species!="PUPA"), aes(x=real_data, y=predicted_data, color = species_name))+
  geom_abline(slope = 1, intercept = 0,color = "deepskyblue", linetype = "dashed", size = 1.2, alpha =.5)+
  geom_point(size=5,alpha=0.5)+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+
  theme_bw()+
  labs(x="Observed abundance",y="Predicted abundance",color=NULL)+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

dev.off()


# Check residuals with Dharma

# Simulating residuals with Dahrma
library(DHARMa)
res_abundance_model <- simulateResiduals(fittedModel = abundance_model, n = 1500)
res_abundance_excl_ratio_model <- simulateResiduals(fittedModel = abundance_excl_ratio_model, n = 1500)
res_abundance_arc_distance_model <- simulateResiduals(fittedModel = abundance_arc_distance_model, n = 1500)
res_abundance_excl_ratio_arc_distance_model <- simulateResiduals(fittedModel = abundance_excl_ratio_arc_distance_model, n = 1500)

res_dif_abundance_model <- simulateResiduals(fittedModel = dif_abundance_model, n = 1500)
res_dif_abundance_excl_ratio_model <- simulateResiduals(fittedModel = dif_abundance_excl_ratio_model, n = 1500)
res_dif_abundance_arc_distance_model <- simulateResiduals(fittedModel = dif_abundance_arc_distance_model, n = 1500)
res_dif_abundance_excl_ratio_arc_distance_model <- simulateResiduals(fittedModel = dif_abundance_excl_ratio_arc_distance_model, n = 1500)

# Checking Residuals 
testZeroInflation(res_abundance_model)
testDispersion(res_abundance_model) # overdispersion: OK
plot(res_abundance_model) # dispersion: deviation significant

testZeroInflation(res_abundance_excl_ratio_model)
testDispersion(res_abundance_excl_ratio_model) # overdispersion: OK
plot(res_abundance_excl_ratio_model) 

testZeroInflation(res_abundance_arc_distance_model)
testDispersion(res_abundance_arc_distance_model) # overdispersion: OK
plot(res_abundance_arc_distance_model)

testZeroInflation(res_abundance_excl_ratio_arc_distance_model)
testDispersion(res_abundance_excl_ratio_arc_distance_model) # overdispersion: OK
plot(res_abundance_excl_ratio_arc_distance_model) # quantile deviation significant


testZeroInflation(res_dif_abundance_model)
testDispersion(res_dif_abundance_model) # overdispersion: OK
plot(res_dif_abundance_model) # dispersion: deviation significant

testZeroInflation(res_dif_abundance_excl_ratio_model)
testDispersion(res_dif_abundance_excl_ratio_model) # overdispersion: OK
plot(res_dif_abundance_excl_ratio_model) 

testZeroInflation(res_dif_abundance_arc_distance_model)
testDispersion(res_dif_abundance_arc_distance_model) # overdispersion: OK
plot(res_dif_abundance_arc_distance_model)

testZeroInflation(res_dif_abundance_excl_ratio_arc_distance_model)
testDispersion(res_dif_abundance_excl_ratio_arc_distance_model) # overdispersion: OK
plot(res_dif_abundance_excl_ratio_arc_distance_model) # quantile deviation significant


png("Images/model_exc_ratio_residuals.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*0.6, units = "in", res=300*2)

plot(res_abundance_excl_ratio_model) 

dev.off()

png("Images/model_arc_distance_residuals.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*0.6, units = "in", res=300*2)

plot(res_abundance_arc_distance_model) 

dev.off()