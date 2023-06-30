
library(tidyverse)
library(scales)
############################################################
# Load species' arc distances

arc_distance <- read_csv("Data/caracoles_processed_data/caracoles_arc_distance.csv")
grilli_index <- read_csv("Data/caracoles_processed_data/caracoles_grilli.csv")
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
  mutate(year=year,species = as.factor(species))  %>%
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


# MODELS
# Create data frame with input data
data = abundance_model_data %>% mutate(ratio_tplus1_over_t=total_individuals/total_individuals_prev_year,
                                       difference_tplus1_minus_t=total_individuals-total_individuals_prev_year,
                                       number_species_prev_year = round(exclusion_ratio_prev_year/prob_excl_prev_year,0))

library(GGally)
# Pairwise visualization
ggpairs(data %>% ungroup() %>% dplyr::select("total_individuals","total_individuals_prev_year",
                               "prob_excl_prev_year","exclusion_ratio_prev_year",
                               "arc_distance_prev_year"), title="correlogram population size") 

ggpairs(data  %>% ungroup() %>% dplyr::select("ratio_tplus1_over_t", "total_individuals_prev_year",
                               "prob_excl_prev_year","exclusion_ratio_prev_year",
                               "arc_distance_prev_year"), title="correlogram population growth") 


# png("Images/correlogram log10(population size).png",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27, units = "in", res=300*2)
# ggpairs(data %>% ungroup()  %>% mutate(log10_population_size=log10(total_individuals+1)) %>% 
#           rename(exclusion_distance_prev_year = arc_distance_prev_year) %>% 
#           dplyr::select("log10_population_size","total_individuals_prev_year",
#                                              "prob_excl_prev_year","exclusion_ratio_prev_year",
#                                              "exclusion_distance_prev_year","species"), title="correlogram log10(population size)",
#         mapping=ggplot2::aes(colour = species),
#         lower=list(combo=wrap("facethist",binwidth=1))
#         ) 
# 
# dev.off()

# png("Images/correlogram log10(population growth+1).png",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27, units = "in", res=300*2)
# 
# ggpairs(data  %>% ungroup() %>% mutate(log10_ratio_tplus1_over_t = log10(ratio_tplus1_over_t+1))%>% 
#           rename(exclusion_distance_prev_year = arc_distance_prev_year) %>% 
#           dplyr::select("log10_ratio_tplus1_over_t", "total_individuals_prev_year",
#                                               "prob_excl_prev_year","exclusion_ratio_prev_year",
#                                               "exclusion_distance_prev_year","species"), title="correlogram log10(population growth+1)",
#         mapping=ggplot2::aes(colour = species),
#         lower=list(combo=wrap("facethist",binwidth=1)))
# 
# 
# dev.off()


##################################
##################################
# GAMMA FIT ######################
##################################
##################################

ggplot(data,aes(x=total_individuals_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

ggplot(data,aes(x=prob_excl_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

ggplot(data,aes(x=exclusion_ratio_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

ggplot(data,aes(x=arc_distance_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

gamma_individuals_prev_year_model <- glm(ratio_tplus1_over_t ~ total_individuals_prev_year,data = data,
                                         family = Gamma(link = log),
                                         control = list(maxit = 50)) # Error
gamma_exclusion_ratio_model <- glm(ratio_tplus1_over_t ~ exclusion_ratio_prev_year,data = data,
                                   family = Gamma(link = log),
                                   control = list(maxit = 50)) # Converged
gamma_arc_distance_model <- glm(ratio_tplus1_over_t ~ arc_distance_prev_year,data = data,
                                family = Gamma(link = log),
                                control = list(maxit = 50))# Converged
gamma_prob_excl_model <- glm(ratio_tplus1_over_t ~ prob_excl_prev_year,data = data,
                             family = Gamma(link = log),
                             control = list(maxit = 100))# Converged

summary(gamma_exclusion_ratio_model)
summary(gamma_prob_excl_model)
summary(gamma_arc_distance_model)

##################################
##################################
# GAUSSIAN FIT ###################
##################################
##################################
gaussian_individuals_prev_year_model <- glm(ratio_tplus1_over_t ~ total_individuals_prev_year,data = data,
                                         family = gaussian(link = log),
                                         control = list(maxit = 50))# Converged
gaussian_exclusion_ratio_model <- glm(ratio_tplus1_over_t ~ exclusion_ratio_prev_year,data = data,
                                   family = gaussian(link = log),
                                   control = list(maxit = 50)) # Did not converged
gaussian_arc_distance_model <- glm(ratio_tplus1_over_t ~ arc_distance_prev_year,data = data,
                                family = gaussian(link = log),
                                control = list(maxit = 50))# Converged
gaussian_prob_excl_model <- glm(ratio_tplus1_over_t ~ prob_excl_prev_year,data = data,
                             family = gaussian(link = log),
                             control = list(maxit = 100))# Converged

summary(gaussian_individuals_prev_year_model)
summary(gaussian_exclusion_ratio_model) # Algorithm did not converge
summary(gaussian_prob_excl_model)
summary(gaussian_arc_distance_model)


AIC(
  gaussian_individuals_prev_year_model,
  gaussian_exclusion_ratio_model,
  gaussian_arc_distance_model,
  gaussian_prob_excl_model
)


##################################
##################################
# EXPO. FIT ######################
##################################
##################################
# Update data for exponential fit
data_growth <- data  %>% ungroup() %>% mutate(log10_ratio_tplus1_over_t = log10(ratio_tplus1_over_t+1))%>% 
  rename(exclusion_distance_prev_year = arc_distance_prev_year) %>% 
  dplyr::select("year","log10_ratio_tplus1_over_t", "total_individuals_prev_year",
                "prob_excl_prev_year","exclusion_ratio_prev_year",
                "exclusion_distance_prev_year","species") %>% mutate(number_species_prev_year= exclusion_ratio_prev_year/prob_excl_prev_year)


# Data visualization
ggplot(data_growth,aes(x=total_individuals_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()

ggplot(data_growth,aes(x=prob_excl_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()

ggplot(data_growth,aes(x=exclusion_ratio_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()

ggplot(data_growth,aes(x=exclusion_distance_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()


data_growth$type_exclusion_ratio <- NA

for (i in 1:nrow(data_growth)) {
  if(data_growth$exclusion_ratio_prev_year[i] < .75){
    data_growth$type_exclusion_ratio[i] <- "0 < Excl. ratio < 0.75"
  }else if(data_growth$exclusion_ratio_prev_year[i] < 1.25){
    data_growth$type_exclusion_ratio[i] <- "0.75 < Excl. ratio < 1.25"
  }else{
    data_growth$type_exclusion_ratio[i] <- "1.25 < Excl. ratio"
  }
}



library(latex2exp)

real_data <- ggplot(data_growth,aes(x=type_exclusion_ratio,y=log10_ratio_tplus1_over_t))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("$\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(title="Observed data")+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=0, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )

real_data

plot_prob_exp_neg <- ggplot(data_growth, aes(x = prob_excl_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ c*exp(x*d), 
              method.args = list(start = list(c = 0.5,d=1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ c*exp(x*d)",y="r0")


plot_excl_ratio_exp_neg <- ggplot(data_growth, aes(x = exclusion_ratio_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ c*exp(x*d), 
              method.args = list(start = list(c = 0.5,d=-1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ c*exp(x*d)",y="r0")

plot_excl_ratio_exp_neg

plot_excl_ratio_hyp <- ggplot(data_growth, aes(x = exclusion_ratio_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ a/(1+b*x), 
              method.args = list(start = list(a = 0.5,b=1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ a/(1+b*x)",y="r0")
plot_excl_ratio_hyp 

plot_excl_ratio_exp_neg <- ggplot(data_growth, aes(x = exclusion_ratio_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ c*exp(x*d), 
              method.args = list(start = list(c = 0.1,d=-1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ c*exp(x*d)",y="r0")

plot_excl_ratio_exp_neg

plot_exclusion_distance_exp_neg <- ggplot(data_growth, aes(x = exclusion_distance_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ c*exp(x*d), 
              method.args = list(start = list(c = 0.1,d=-1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ c*exp(x*d)",y="r0")

plot_exclusion_distance_exp_neg

plot_abundance_exp_neg <- ggplot(data_growth, aes(x = total_individuals_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ c*exp(x*d), 
              method.args = list(start = list(c = 0.1,d=-1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ c*exp(x*d)")

plot_abundance_exp_neg # ERROR


library(patchwork)
png("Images/non_linear_fits.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.55, units = "in", res=300)

(plot_prob_exp_neg)/(plot_excl_ratio_exp_neg|plot_excl_ratio_hyp )/plot_prob_exp_neg/plot_exclusion_distance_exp_neg#/plot_abundance_exp_neg
dev.off()



#############################
#############################
###### FIT MODELS EXP NEG ###
#############################
#############################


fit_excl_ratio_hyp <- nls(log10_ratio_tplus1_over_t ~ a/(1+b*exclusion_ratio_prev_year), 
                              data = data_growth,
                              start = list(a = 0.5, b=1))
summary(fit_excl_ratio_hyp)


fit_excl_ratio_exp_neg <- nls(log10_ratio_tplus1_over_t ~ c*exp(exclusion_ratio_prev_year*d), 
                              data = data_growth,
                              start = list(c = 0.5, d=-1))
summary(fit_excl_ratio_exp_neg)


fit_abundance_exp_neg <- nls(log10_ratio_tplus1_over_t ~ c*exp(total_individuals_prev_year*d), 
                                   data = data_growth,
                                   start = list(c = 0.5, d=-1)) # ERROR

fit_prob_excl_exp_neg <- nls(log10_ratio_tplus1_over_t ~ c*exp(prob_excl_prev_year*d), 
                              data = data_growth,
                              start = list(c = 0.5, d=-1))
summary(fit_prob_excl_exp_neg)

fit_excl_dist_exp_neg <- nls(log10_ratio_tplus1_over_t ~ c*exp(exclusion_distance_prev_year*d), 
                             data = data_growth,
                             start = list(c = 0.5, d=-1))
summary(fit_excl_dist_exp_neg)

# AIC comparisson

AIC(fit_excl_ratio_exp_neg,fit_excl_ratio_hyp,fit_prob_excl_exp_neg,fit_excl_dist_exp_neg,gamma_exclusion_ratio_model,
    gamma_arc_distance_model,
    gamma_prob_excl_model,
    gaussian_individuals_prev_year_model,
    gaussian_exclusion_ratio_model,
    gaussian_arc_distance_model,
    gaussian_prob_excl_model)


# Test the ACF and PACF of the best model's residuals
nlm_excl_ratio_exp_neg_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                    predicted_data = (predict(fit_excl_ratio_exp_neg, data_growth)),
                                    residuals = residuals(fit_excl_ratio_exp_neg),
                                    species=data$species,
                                    year=data$year,
                                    exclusion_ratio_prev_year=data$exclusion_ratio_prev_year)

BEMA_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="BEMA",c(5,3)]
BEMA_series <- bind_rows(BEMA_aux_series[1:2,],tibble(year=2018,residuals=NA),BEMA_aux_series[3:5,])
residuals_BEMA <- ts(BEMA_series[,2],start = 2016)

CETE_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="CETE",c(5,3)]
CETE_series <- bind_rows(CETE_aux_series[1:2,],tibble(year=2018,residuals=NA),CETE_aux_series[3:5,])
residuals_CETE <- ts(CETE_series[,2],start = 2016)

HOMA_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="HOMA",c(5,3)]
HOMA_series <- bind_rows(HOMA_aux_series[1:2,],tibble(year=2018,residuals=NA),HOMA_aux_series[3:5,])
residuals_HOMA <- ts(HOMA_series[,2],start = 2016)


LEMA_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="LEMA",c(5,3)]
LEMA_series <- bind_rows(LEMA_aux_series[1:2,],tibble(year=2018,residuals=NA),LEMA_aux_series[3:5,])
residuals_LEMA <- ts(LEMA_series[,2],start = 2016)


PAIN_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="PAIN",c(5,3)]
PAIN_series <- bind_rows(PAIN_aux_series[1:2,],tibble(year=2018,residuals=NA),PAIN_aux_series[3:5,])
residuals_PAIN <- ts(PAIN_series[,2],start = 2016)

POMA_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="POMA",c(5,3)]
POMA_series <- bind_rows(POMA_aux_series[1:2,],tibble(year=2018,residuals=NA),POMA_aux_series[3:5,])
residuals_POMA <- ts(POMA_series[,2],start = 2016)


SASO_aux_series <- nlm_excl_ratio_exp_neg_predictions[nlm_excl_ratio_exp_neg_predictions$species=="SASO",c(5,3)]
SASO_series <- bind_rows(SASO_aux_series[1:2,],tibble(year=2018,residuals=NA),SASO_aux_series[3:5,])
residuals_SASO <- ts(SASO_series[,2],start = 2016)



library(TSstudio)
acf(residuals_BEMA,na.action = na.pass) #OK
pacf(residuals_BEMA,na.action = na.pass) #No OK

acf(residuals_CETE,na.action = na.pass) #OK
pacf(residuals_CETE,na.action = na.pass) #NO OK

acf(residuals_HOMA,na.action = na.pass) #OK
pacf(residuals_HOMA,na.action = na.pass) #OK

acf(residuals_LEMA,na.action = na.pass) #OK
pacf(residuals_LEMA, na.action = na.pass) #OK

acf(residuals_PAIN,na.action = na.pass) #OK
pacf(residuals_PAIN,na.action = na.pass) #OK

acf(residuals_POMA,na.action = na.pass) #OK
pacf(residuals_POMA,na.action = na.pass) #OK

acf(residuals_SASO,na.action = na.pass) #OK
pacf(residuals_SASO,na.action = na.pass) # No OK


# Plot results for the model with the smallest AIC
nlm_excl_ratio_exp_neg_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                             predicted_data = (predict(fit_excl_ratio_exp_neg, data_growth)),
                                             residuals = residuals(fit_excl_ratio_exp_neg),
                                             species=data$species,
                                             year=data$year,
                                             exclusion_ratio_prev_year=data$exclusion_ratio_prev_year)

nlm_excl_ratio_exp_neg_predictions$type_exclusion_ratio <- NA

for (i in 1:nrow(nlm_excl_ratio_exp_neg_predictions)) {
  if(nlm_excl_ratio_exp_neg_predictions$exclusion_ratio_prev_year[i] < .75){
    nlm_excl_ratio_exp_neg_predictions$type_exclusion_ratio[i] <- "0 < Excl. ratio < 0.75"
  }else if(nlm_excl_ratio_exp_neg_predictions$exclusion_ratio_prev_year[i] < 1.25){
    nlm_excl_ratio_exp_neg_predictions$type_exclusion_ratio[i] <- "0.75 < Excl. ratio < 1.25"
  }else{
    nlm_excl_ratio_exp_neg_predictions$type_exclusion_ratio[i] <- "1.25 < Excl. ratio"
  }
}



library(latex2exp)

modeled_data <- ggplot(nlm_excl_ratio_exp_neg_predictions,aes(x=type_exclusion_ratio,y=predicted_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("$\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(title="Modeled data")+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=0, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+
  ylim(c(-1,2.5))

modeled_data
real_data


data_boxplot <- bind_rows(data_growth %>% mutate(type_data = "Observed") %>% rename(y = log10_ratio_tplus1_over_t),
                          nlm_excl_ratio_exp_neg_predictions %>% mutate(type_data = "Modeled") %>% rename(y = predicted_data)) 


population_growth_boxplots <- ggplot(data = data_boxplot, aes(x=type_exclusion_ratio,y=y, fill = type_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  scale_fill_manual(values=c("grey89", "white"))+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(fill=NULL)+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_boxplots


data_growth$species_name[data_growth$species == "BEMA"] <- "Beta\nmacrocarpa"
data_growth$species_name[data_growth$species == "CETE"] <- "Centaurium\ntenuiflorum"
data_growth$species_name[data_growth$species == "HOMA"] <- "Hordeum\nmarinum"
data_growth$species_name[data_growth$species == "LEMA"] <- "Leontodon\nmaroccanus"
data_growth$species_name[data_growth$species == "PAIN"] <- "Parapholis\nincurva"
data_growth$species_name[data_growth$species == "POMA"] <- "Polypogon\nmaritimus"
data_growth$species_name[data_growth$species == "SASO"] <- "Salsola\nsoda"
data_growth$species_name[data_growth$species == "SCLA"] <- "Scorzonera\nlaciniata"
data_growth$species_name[data_growth$species == "PUPA"] <- "Pulicaria\npaludosa"


# Extract nls mean curve
mean_curve_fit <- tibble(exclusion_ratio_prev_year = seq(min(data_growth$exclusion_ratio_prev_year), 2.5, by = 0.01),
                         log10_ratio_tplus1_over_t = 0.6462*exp(-0.3181*exclusion_ratio_prev_year))

# Extract mean curve fir and confidence intervals from nls
library(investr)

new.data <- data.frame(exclusion_ratio_prev_year = seq(min(data_growth$exclusion_ratio_prev_year),2.5,by = 0.01))
interval <- as_tibble(investr::predFit(fit_excl_ratio_exp_neg, newdata = new.data, interval = "confidence", level= 0.95)) %>% 
  mutate(exclusion_ratio_prev_year = new.data$exclusion_ratio_prev_year)


population_growth_plants <- ggplot(mean_curve_fit, aes(x=exclusion_ratio_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_line(lwd=2,color="grey")+
  geom_ribbon(data=interval, aes(x=exclusion_ratio_prev_year, ymin=lwr, ymax=upr), alpha=0.3, inherit.aes=F, fill="gray")+
  geom_point(data=data_growth,aes(x=exclusion_ratio_prev_year,y=log10_ratio_tplus1_over_t,color=species_name),
             size=5,alpha=0.5)+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  xlim(c(0,2.5))+
  labs(x="Exclusion ratio",color=NULL)+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_plants

library(patchwork)
population_growth_plants + population_growth_boxplots

png("Images/neg_exp_population_growth_model_ER.png",
    width = 11.69*1.5, # The width of the plot in inches
    height = 11.69*0.8, units = "in", res=300*2)
(population_growth_plants+
    ylim(c(-.2,2.5))) + (population_growth_boxplots+
                           ylim(c(-.2,2.5)))
dev.off()


##############################
##############################
##############################
##############################
#  Model that corrects temporal correlation and random groups for the previous best model

lmer_simple_fit <- lme(log(log10_ratio_tplus1_over_t) ~ exclusion_ratio_prev_year, 
                       random= ~ 1|species,
                       data=data_growth%>% mutate(species=as.factor(species)),
                       correlation=corCAR1(form= ~ year|species)
)

summary(lmer_simple_fit)

ranef(lmer_simple_fit)

AIC(fit_excl_ratio_exp_neg,lmer_simple_fit) # This new model has twice AIC


###########################
# lmer model visualization
lmer_simple_fit_exp_neg_predictions_2 <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                                predicted_data = exp(predict(lmer_simple_fit, data_growth)),
                                                residuals = residuals(lmer_simple_fit),
                                                species=data$species,
                                                year=data$year,
                                                exclusion_ratio_prev_year=data$exclusion_ratio_prev_year)
lmer_simple_fit_exp_neg_predictions_2$type_exclusion_ratio <- NA

for (i in 1:nrow(lmer_simple_fit_exp_neg_predictions_2)) {
  if(lmer_simple_fit_exp_neg_predictions_2$exclusion_ratio_prev_year[i] < .75){
    lmer_simple_fit_exp_neg_predictions_2$type_exclusion_ratio[i] <- "0 < Excl. ratio < 0.75"
  }else if(lmer_simple_fit_exp_neg_predictions_2$exclusion_ratio_prev_year[i] < 1.25){
    lmer_simple_fit_exp_neg_predictions_2$type_exclusion_ratio[i] <- "0.75 < Excl. ratio < 1.25"
  }else{
    lmer_simple_fit_exp_neg_predictions_2$type_exclusion_ratio[i] <- "1.25 < Excl. ratio"
  }
}



library(latex2exp)

modeled_data <- ggplot(lmer_simple_fit_exp_neg_predictions_2,aes(x=type_exclusion_ratio,y=predicted_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("$\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(title="Modeled data")+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=0, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+
  ylim(c(0,2.5))

modeled_data
real_data

############################
# check residualsexp(predict(lmer_simple_fit,mean_HOMA_fit))

lmer_exp_neg_predictions_2 <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                     predicted_data = exp(predict(lmer_simple_fit,data_growth)),
                                     residuals = residuals(lmer_simple_fit),
                                     species=data$species,
                                     year=data$year,
                                     exclusion_ratio_prev_year=data$exclusion_ratio_prev_year)


BEMA_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="BEMA",c(5,3)]
BEMA_series <- bind_rows(BEMA_aux_series[1:2,],tibble(year=2018,residuals=NA),BEMA_aux_series[3:5,])
residuals_BEMA <- ts(BEMA_series[,2],start = 2016)

CETE_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="CETE",c(5,3)]
CETE_series <- bind_rows(CETE_aux_series[1:2,],tibble(year=2018,residuals=NA),CETE_aux_series[3:5,])
residuals_CETE <- ts(CETE_series[,2],start = 2016)

HOMA_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="HOMA",c(5,3)]
HOMA_series <- bind_rows(HOMA_aux_series[1:2,],tibble(year=2018,residuals=NA),HOMA_aux_series[3:5,])
residuals_HOMA <- ts(HOMA_series[,2],start = 2016)


LEMA_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="LEMA",c(5,3)]
LEMA_series <- bind_rows(LEMA_aux_series[1:2,],tibble(year=2018,residuals=NA),LEMA_aux_series[3:5,])
residuals_LEMA <- ts(LEMA_series[,2],start = 2016)


PAIN_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="PAIN",c(5,3)]
PAIN_series <- bind_rows(PAIN_aux_series[1:2,],tibble(year=2018,residuals=NA),PAIN_aux_series[3:5,])
residuals_PAIN <- ts(PAIN_series[,2],start = 2016)

POMA_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="POMA",c(5,3)]
POMA_series <- bind_rows(POMA_aux_series[1:2,],tibble(year=2018,residuals=NA),POMA_aux_series[3:5,])
residuals_POMA <- ts(POMA_series[,2],start = 2016)


SASO_aux_series <- lmer_exp_neg_predictions_2[lmer_exp_neg_predictions_2$species=="SASO",c(5,3)]
SASO_series <- bind_rows(SASO_aux_series[1:2,],tibble(year=2018,residuals=NA),SASO_aux_series[3:5,])
residuals_SASO <- ts(SASO_series[,2],start = 2016)


library(TSstudio)
acf(residuals_BEMA,na.action = na.pass) #OK 
pacf(residuals_BEMA,na.action = na.pass) # NO OK

acf(residuals_CETE,na.action = na.pass) #OK
pacf(residuals_CETE,na.action = na.pass) #NO OK

acf(residuals_HOMA,na.action = na.pass) #OK
pacf(residuals_HOMA,na.action = na.pass) #OK

acf(residuals_LEMA,na.action = na.pass) #OK
pacf(residuals_LEMA, na.action = na.pass) #OK

acf(residuals_PAIN,na.action = na.pass) #OK
pacf(residuals_PAIN,na.action = na.pass) # NO OK

acf(residuals_POMA,na.action = na.pass) #OK
pacf(residuals_POMA,na.action = na.pass) #OK

acf(residuals_SASO,na.action = na.pass) #OK
pacf(residuals_SASO,na.action = na.pass) # NO OK

###################################
###################################

# Plot results for exclusion distance exponential NLS fit

nlm_excl_distance_exp_neg_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                             predicted_data = (predict(fit_excl_dist_exp_neg, data_growth)),
                                             residuals = residuals(fit_excl_ratio_exp_neg),
                                             species=data$species,
                                             year=data$year,
                                             exclusion_distance_prev_year=data$arc_distance_prev_year)


ggplot(nlm_excl_distance_exp_neg_predictions, aes(x=exclusion_distance_prev_year,
                                                   y=real_data)) +
  geom_point()

nlm_excl_distance_exp_neg_predictions$type_exclusion_ratio <- NA

for (i in 1:nrow(nlm_excl_distance_exp_neg_predictions)) {
  if(nlm_excl_distance_exp_neg_predictions$exclusion_distance_prev_year[i] < .025){
    nlm_excl_distance_exp_neg_predictions$type_exclusion_ratio[i] <- "0 < Excl. distance < 0.025"
  }else{
    nlm_excl_distance_exp_neg_predictions$type_exclusion_ratio[i] <- "0.025 < Excl. distance"
  }
}



library(latex2exp)

modeled_data <- ggplot(nlm_excl_distance_exp_neg_predictions,aes(x=type_exclusion_ratio,y=predicted_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("$\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(title="Modeled data")+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=0, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+
  ylim(c(-1,2.5))

modeled_data

data_growth_excl_distance <- data  %>% ungroup() %>% mutate(log10_ratio_tplus1_over_t = log10(ratio_tplus1_over_t+1))%>% 
  rename(exclusion_distance_prev_year = arc_distance_prev_year) %>% 
  dplyr::select("year","log10_ratio_tplus1_over_t", 
                "prob_excl_prev_year","exclusion_ratio_prev_year",
                "exclusion_distance_prev_year","species") %>% mutate(number_species_prev_year= exclusion_ratio_prev_year/prob_excl_prev_year)

ggplot(data_growth_excl_distance,aes(x=exclusion_distance_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()


data_growth_excl_distance$type_exclusion_ratio <- NA

for (i in 1:nrow(data_growth_excl_distance)) {
  if(data_growth_excl_distance$exclusion_distance_prev_year[i] < .025){
    data_growth_excl_distance$type_exclusion_ratio[i] <- "0 < Excl. distance < 0.025"
  }else{
    data_growth_excl_distance$type_exclusion_ratio[i] <- "0.025 < Excl. distance"
  }
}



data_boxplot_excl_dist <- bind_rows(data_growth_excl_distance %>% mutate(type_data = "Observed") %>% rename(y = log10_ratio_tplus1_over_t),
                                    nlm_excl_distance_exp_neg_predictions %>% mutate(type_data = "Modeled") %>% rename(y = predicted_data)) 


population_growth_boxplots_exclusion_distance <- ggplot(data = data_boxplot_excl_dist, aes(x=type_exclusion_ratio,y=y, fill = type_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  scale_fill_manual(values=c("grey89", "white"))+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(fill=NULL)+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_boxplots_exclusion_distance


data_growth_excl_distance$species_name[data_growth_excl_distance$species == "BEMA"] <- "Beta\nmacrocarpa"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "CETE"] <- "Centaurium\ntenuiflorum"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "HOMA"] <- "Hordeum\nmarinum"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "LEMA"] <- "Leontodon\nmaroccanus"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "PAIN"] <- "Parapholis\nincurva"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "POMA"] <- "Polypogon\nmaritimus"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "SASO"] <- "Salsola\nsoda"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "SCLA"] <- "Scorzonera\nlaciniata"
data_growth_excl_distance$species_name[data_growth_excl_distance$species == "PUPA"] <- "Pulicaria\npaludosa"


# Extract nls mean curve
mean_curve_fit_exclusion_distance <- tibble(exclusion_distance_prev_year = seq(min(data_growth_excl_distance$exclusion_distance_prev_year), 
                                                            max(data_growth_excl_distance$exclusion_distance_prev_year), by = 0.01),
                         log10_ratio_tplus1_over_t = 0.5312*exp(3.3297*exclusion_distance_prev_year))

# Extract mean curve fir and confidence intervals from nls
library(investr)

new.data_exclusion_distance <- data.frame(exclusion_distance_prev_year = seq(min(data_growth_excl_distance$exclusion_distance_prev_year), 
                                                                             max(data_growth_excl_distance$exclusion_distance_prev_year),by = 0.001))
interval_exclusion_distance <- as_tibble(investr::predFit(fit_excl_dist_exp_neg, 
                                                          newdata = new.data_exclusion_distance, 
                                                          interval = "confidence", level= 0.95)) %>% 
  mutate(exclusion_distance_prev_year = new.data_exclusion_distance$exclusion_distance_prev_year)


population_growth_plants_exclusion_distance <- ggplot(mean_curve_fit_exclusion_distance, aes(x=exclusion_distance_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_line(lwd=2,color="grey")+
  geom_ribbon(data=interval_exclusion_distance, aes(x=exclusion_distance_prev_year, ymin=lwr, ymax=upr), alpha=0.3, inherit.aes=F, fill="gray")+
  geom_point(data=data_growth_excl_distance,aes(x=exclusion_distance_prev_year,y=log10_ratio_tplus1_over_t,color=species_name),
             size=5,alpha=0.5)+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  #xlim(c(0,2.5))+
  labs(x="Exclusion distance",color=NULL)+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_plants_exclusion_distance

library(patchwork)
population_growth_plants_exclusion_distance + population_growth_boxplots_exclusion_distance

png("Images/caracoles_excl_distance_negative_exp_simple_COMBO.png",
    width = 11.69*1.5, # The width of the plot in inches
    height = 11.69*0.8, units = "in", res=300*2)
(population_growth_plants_exclusion_distance+
    ylim(c(-.1,2.5))) + (population_growth_boxplots_exclusion_distance+
                           ylim(c(-.1,2.5)))
dev.off()

###################################
###################################

# Plot results for probability of exclusion exponential NLS fit

nlm_prob_excl_exp_neg_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                                predicted_data = (predict(fit_prob_excl_exp_neg, data_growth)),
                                                residuals = residuals(fit_excl_ratio_exp_neg),
                                                species=data$species,
                                                year=data$year,
                                            prob_excl_prev_year=data$prob_excl_prev_year)


ggplot(nlm_prob_excl_exp_neg_predictions, aes(x=prob_excl_prev_year,
                                                  y=real_data)) +
  geom_point()

nlm_prob_excl_exp_neg_predictions$type_exclusion_ratio <- NA

for (i in 1:nrow(nlm_excl_distance_exp_neg_predictions)) {
  if(nlm_prob_excl_exp_neg_predictions$prob_excl_prev_year[i] < .15){
    nlm_prob_excl_exp_neg_predictions$type_exclusion_ratio[i] <- "0 < Prob. of excl. < 0.15"
  }else{
    nlm_prob_excl_exp_neg_predictions$type_exclusion_ratio[i] <- "0.15 < Prob. of excl."
  }
}



library(latex2exp)

modeled_data_prob_excl <- ggplot(nlm_prob_excl_exp_neg_predictions,aes(x=type_exclusion_ratio,y=predicted_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("$\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(title="Modeled data")+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=0, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+
  ylim(c(-1,2.5))

modeled_data_prob_excl

data_growth_prob_excl <- data  %>% ungroup() %>% mutate(log10_ratio_tplus1_over_t = log10(ratio_tplus1_over_t+1))%>% 
  rename(prob_excl_prev_year = prob_excl_prev_year) %>% 
  dplyr::select("year","log10_ratio_tplus1_over_t", 
                "prob_excl_prev_year","exclusion_ratio_prev_year",
                "species") %>% mutate(number_species_prev_year= exclusion_ratio_prev_year/prob_excl_prev_year)

ggplot(data_growth_prob_excl,aes(x=prob_excl_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()


data_growth_prob_excl$type_exclusion_ratio <- NA

for (i in 1:nrow(data_growth_excl_distance)) {
  if(data_growth_prob_excl$prob_excl_prev_year[i] < .15){
    data_growth_prob_excl$type_exclusion_ratio[i] <- "0 < Prob. of excl. < 0.15"
  }else{
    data_growth_prob_excl$type_exclusion_ratio[i] <- "0.15 < Prob. of excl."
  }
}



data_boxplot_prob_excl <- bind_rows(data_growth_prob_excl %>% mutate(type_data = "Observed") %>% rename(y = log10_ratio_tplus1_over_t),
                                    nlm_prob_excl_exp_neg_predictions %>% mutate(type_data = "Modeled") %>% rename(y = predicted_data)) 


population_growth_boxplots_prob_excl <- ggplot(data = data_boxplot_prob_excl, aes(x=type_exclusion_ratio,y=y, fill = type_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  scale_fill_manual(values=c("grey89", "white"))+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(fill=NULL)+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_boxplots_exclusion_distance


data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "BEMA"] <- "Beta\nmacrocarpa"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "CETE"] <- "Centaurium\ntenuiflorum"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "HOMA"] <- "Hordeum\nmarinum"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "LEMA"] <- "Leontodon\nmaroccanus"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "PAIN"] <- "Parapholis\nincurva"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "POMA"] <- "Polypogon\nmaritimus"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "SASO"] <- "Salsola\nsoda"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "SCLA"] <- "Scorzonera\nlaciniata"
data_boxplot_prob_excl$species_name[data_boxplot_prob_excl$species == "PUPA"] <- "Pulicaria\npaludosa"


# Extract nls mean curve
mean_curve_fit_prob_excl <- tibble(prob_excl_prev_year = seq(min(data_boxplot_prob_excl$prob_excl_prev_year), 
                                                                               max(data_boxplot_prob_excl$prob_excl_prev_year), by = 0.01),
                                            log10_ratio_tplus1_over_t = 0.5926*exp(-0.8000*prob_excl_prev_year))

# Extract mean curve fir and confidence intervals from nls
library(investr)

new.data_prob_excl <- data.frame(prob_excl_prev_year = seq(min(data_boxplot_prob_excl$prob_excl_prev_year), 
                                                                             max(data_boxplot_prob_excl$prob_excl_prev_year),by = 0.001))
interval_prob_excl <- as_tibble(investr::predFit(fit_prob_excl_exp_neg, 
                                                          newdata = new.data_prob_excl, 
                                                          interval = "confidence", level= 0.95)) %>% 
  mutate(prob_excl_prev_year = new.data_prob_excl$prob_excl_prev_year)


population_growth_plants_prob_excl <- ggplot(mean_curve_fit_prob_excl, aes(x=prob_excl_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_line(lwd=2,color="grey")+
  geom_ribbon(data=interval_prob_excl, aes(x=prob_excl_prev_year, ymin=lwr, ymax=upr), alpha=0.3, inherit.aes=F, fill="gray")+
  geom_point(data=data_boxplot_prob_excl,aes(x=prob_excl_prev_year,y=y,color=species_name),
             size=5,alpha=0.5)+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  #xlim(c(0,2.5))+
  labs(x="Prob. of excl.",color=NULL)+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_plants_prob_excl

library(patchwork)
population_growth_plants_prob_excl + population_growth_boxplots_prob_excl

png("Images/caracoles_prob_excl_negative_exp_simple_COMBO.png",
    width = 11.69*1.5, # The width of the plot in inches
    height = 11.69*0.8, units = "in", res=300*2)
(population_growth_plants_prob_excl+
    ylim(c(-.1,2.5))) + (population_growth_boxplots_prob_excl+
                           ylim(c(-.1,2.5)))
dev.off()

#################################################
#################################################
#################################################

# Plot (fixed effect) results for GLMM

lmer_simple_fit_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                                predicted_data = exp(predict(lmer_simple_fit, data_growth)),
                                                residuals = residuals(fit_excl_ratio_exp_neg),
                                                species=data$species,
                                                year=data$year,
                                      exclusion_ratio_prev_year=data$exclusion_ratio_prev_year)


ggplot(lmer_simple_fit_predictions, aes(x=exclusion_ratio_prev_year,
                                                  y=real_data)) +
  geom_point()


lmer_simple_fit_predictions$type_exclusion_ratio <- NA

for (i in 1:nrow(lmer_simple_fit_predictions)) {
  if(lmer_simple_fit_predictions$exclusion_ratio_prev_year[i] < .75){
    lmer_simple_fit_predictions$type_exclusion_ratio[i] <- "0 < Excl. ratio < 0.75"
  }else if(lmer_simple_fit_predictions$exclusion_ratio_prev_year[i] < 1.25){
    lmer_simple_fit_predictions$type_exclusion_ratio[i] <- "0.75 < Excl. ratio < 1.25"
  }else{
    lmer_simple_fit_predictions$type_exclusion_ratio[i] <- "1.25 < Excl. ratio"
  }
}



library(latex2exp)

modeled_data_lmer <- ggplot(lmer_simple_fit_predictions,aes(x=type_exclusion_ratio,y=predicted_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  ylab(TeX("$\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(title="Modeled data")+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=0, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+
  ylim(c(0,2.5))

modeled_data_lmer


data_growth_lmer <- data  %>% ungroup() %>% mutate(log10_ratio_tplus1_over_t = log10(ratio_tplus1_over_t+1))%>% 
  rename(exclusion_distance_prev_year = arc_distance_prev_year) %>% 
  dplyr::select("year","log10_ratio_tplus1_over_t", 
                "prob_excl_prev_year","exclusion_ratio_prev_year",
                "exclusion_distance_prev_year","species") %>% mutate(number_species_prev_year= exclusion_ratio_prev_year/prob_excl_prev_year)

ggplot(data_growth_lmer,aes(x=exclusion_ratio_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_point()


data_growth_lmer$type_exclusion_ratio <- NA

for (i in 1:nrow(data_growth_lmer)) {
  if(data_growth_lmer$exclusion_ratio_prev_year[i] < .75){
    data_growth_lmer$type_exclusion_ratio[i] <- "0 < Excl. ratio < 0.75"
  }else if(data_growth_lmer$exclusion_ratio_prev_year[i] < 1.25){
    data_growth_lmer$type_exclusion_ratio[i] <- "0.75 < Excl. ratio < 1.25"
  }else{
    data_growth_lmer$type_exclusion_ratio[i] <- "1.25 < Excl. ratio"
  }
}



data_boxplot_lmer <- bind_rows(data_growth_lmer %>% mutate(type_data = "Observed") %>% rename(y = log10_ratio_tplus1_over_t),
                               lmer_simple_fit_predictions %>% mutate(type_data = "Modeled") %>% rename(y = predicted_data)) 


population_growth_boxplots_lmer <- ggplot(data = data_boxplot_lmer, aes(x=type_exclusion_ratio,y=y, fill = type_data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  scale_fill_manual(values=c("grey89", "white"))+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  xlab(NULL)+
  labs(fill=NULL)+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_boxplots_lmer


data_boxplot_lmer$species_name[data_boxplot_lmer$species == "BEMA"] <- "Beta\nmacrocarpa"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "CETE"] <- "Centaurium\ntenuiflorum"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "HOMA"] <- "Hordeum\nmarinum"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "LEMA"] <- "Leontodon\nmaroccanus"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "PAIN"] <- "Parapholis\nincurva"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "POMA"] <- "Polypogon\nmaritimus"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "SASO"] <- "Salsola\nsoda"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "SCLA"] <- "Scorzonera\nlaciniata"
data_boxplot_lmer$species_name[data_boxplot_lmer$species == "PUPA"] <- "Pulicaria\npaludosa"


# Extract nls mean curve
exp(predict(lmer_simple_fit, data_growth))
mean_curve_fit_lmer <- tibble(exclusion_ratio_prev_year = seq(min(data_boxplot_lmer$exclusion_ratio_prev_year), 
                                                                               max(data_boxplot_lmer$exclusion_ratio_prev_year), by = 0.01),
                                            log10_ratio_tplus1_over_t = exp(-1.2622388-0.13181617*exclusion_ratio_prev_year))

# Extract mean curve fir and confidence intervals from nls
library(investr)

new.data_lmer <- data.frame(exclusion_ratio_prev_year = seq(min(data_boxplot_lmer$exclusion_ratio_prev_year), 
                                                                             max(data_boxplot_lmer$exclusion_ratio_prev_year),by = 0.001))
interval_lmer <- as_tibble(investr::predFit(lmer_simple_fit, 
                                                          newdata = new.data_lmer, 
                                                          interval = "confidence", level= 0.95)) %>% 
  mutate(exclusion_ratio_prev_year = new.data_lmer$exclusion_ratio_prev_year)


population_growth_plants_lmer <- ggplot(mean_curve_fit_lmer, aes(x=exclusion_ratio_prev_year,y=log10_ratio_tplus1_over_t))+
  geom_line(lwd=2,color="grey")+
  #geom_ribbon(data=interval_lmer, aes(x=exclusion_ratio_prev_year, ymin=lwr, ymax=upr), alpha=0.3, inherit.aes=F, fill="gray")+
  geom_point(data=data_boxplot_lmer,aes(x=exclusion_ratio_prev_year,y=y,color=species_name),
             size=5,alpha=0.5)+
  geom_hline(aes(yintercept=log10(2)),  color="deepskyblue1",
             linetype="dashed", size=1)+
  #xlim(c(0,2.5))+
  labs(x="Exclusion ratio",color=NULL)+
  ylab(TeX("Population growth:  $\\log_{10}\\left ( \\frac{N_{\\,t+1}}{N_{\\,t}} +1 \\right )$"))+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

population_growth_plants_lmer

library(patchwork)
population_growth_plants_lmer + population_growth_boxplots_lmer

png("Images/caracoles_excl_ratio_lmer_COMBO.png",
    width = 11.69*1.5, # The width of the plot in inches
    height = 11.69*0.8, units = "in", res=300*2)
(population_growth_plants_lmer+
    ylim(c(-.1,2.5))) + (population_growth_boxplots_lmer+
                           ylim(c(-.1,2.5)))
dev.off()
