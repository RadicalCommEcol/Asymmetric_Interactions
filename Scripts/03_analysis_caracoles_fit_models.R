
library(tidyverse)
library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(performance)
library(visreg)
library(RColorBrewer)
library(nlme)
library(usdm)
library(GGally)
library(TSstudio)

abundance_model_data <- read_csv("Results/abundance_model_data.csv") %>%
  mutate(year=year,species = as.factor(species))

# Preliminary analyses-----------------------------------------------------------
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


# MODELS------------------------------------------------------------------------
# Create data frame with input data
data = abundance_model_data %>% mutate(ratio_tplus1_over_t=total_individuals/total_individuals_prev_year,
                                       difference_tplus1_minus_t=total_individuals-total_individuals_prev_year,
                                       number_species_prev_year = round(exclusion_ratio_prev_year/prob_excl_prev_year,0))


# Pairwise visualization
ggpairs(data %>% ungroup() %>% dplyr::select("total_individuals","total_individuals_prev_year",
                                             "prob_excl_prev_year","exclusion_ratio_prev_year",
                                             "arc_distance_prev_year"), title="correlogram population size") 

ggpairs(data  %>% ungroup() %>% dplyr::select("ratio_tplus1_over_t", "total_individuals_prev_year",
                                              "prob_excl_prev_year","exclusion_ratio_prev_year",
                                              "arc_distance_prev_year"), title="correlogram population growth") 

# Aditional previsualizations
ggplot(data,aes(x=total_individuals_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

ggplot(data,aes(x=prob_excl_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

ggplot(data,aes(x=exclusion_ratio_prev_year,y=ratio_tplus1_over_t))+
  geom_point()

ggplot(data,aes(x=arc_distance_prev_year,y=ratio_tplus1_over_t))+
  geom_point()


##################################
##################################
# GAMMA FIT ######################
##################################
##################################

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


plot_prob_exp_neg <- ggplot(data_growth, aes(x = prob_excl_prev_year, y = log10_ratio_tplus1_over_t)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ c*exp(x*d), 
              method.args = list(start = list(c = 0.5,d=1)), se = F, #starting values obtained from fit above
              color = "dark red")+
  labs(title = "y ~ c*exp(x*d)",y="r0")

plot_prob_exp_neg

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

AIC(fit_excl_ratio_exp_neg,
    fit_excl_ratio_hyp,
    fit_prob_excl_exp_neg,
    fit_excl_dist_exp_neg)

##############################################################################
##############################################################################

# AIC comparisson for every fitted model

AIC(fit_excl_ratio_exp_neg,
    fit_excl_ratio_hyp,
    fit_prob_excl_exp_neg,
    fit_excl_dist_exp_neg,
    gamma_exclusion_ratio_model,
    gamma_arc_distance_model,
    gamma_prob_excl_model,
    gaussian_individuals_prev_year_model,
    gaussian_exclusion_ratio_model,
    gaussian_arc_distance_model,
    gaussian_prob_excl_model)


#                                     df        AIC
# fit_excl_ratio_exp_neg                3   63.86213
# fit_excl_ratio_hyp                    3   64.08529
# fit_prob_excl_exp_neg                 3   65.02887
# fit_excl_dist_exp_neg                 3   64.86826
# gamma_exclusion_ratio_model           3  197.71752
# gamma_arc_distance_model              3  198.27493
# gamma_prob_excl_model                 3  197.69274
# gaussian_individuals_prev_year_model  3  345.62198
# gaussian_exclusion_ratio_model        3 2314.43864
# gaussian_arc_distance_model           3  369.90875
# gaussian_prob_excl_model              3  369.92612


##############################################################################
##############################################################################
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


# Check temporal (partial)-autocorrelation
acf(residuals_BEMA,na.action = na.pass) #OK
pacf(residuals_BEMA,na.action = na.pass) #No OK

acf(residuals_CETE,na.action = na.pass) #OK
pacf(residuals_CETE,na.action = na.pass) #OK

acf(residuals_HOMA,na.action = na.pass) #OK
pacf(residuals_HOMA,na.action = na.pass) #OK

acf(residuals_LEMA,na.action = na.pass) #OK
pacf(residuals_LEMA, na.action = na.pass) #OK

acf(residuals_PAIN,na.action = na.pass) #OK
pacf(residuals_PAIN,na.action = na.pass) #OK

acf(residuals_POMA,na.action = na.pass) #OK
pacf(residuals_POMA,na.action = na.pass) #OK

acf(residuals_SASO,na.action = na.pass) #OK
pacf(residuals_SASO,na.action = na.pass) # OK


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

# check residuals

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


# Check temporal (partial)-autocorrelation
acf(residuals_BEMA,na.action = na.pass) #OK 
pacf(residuals_BEMA,na.action = na.pass) # NO OK

acf(residuals_CETE,na.action = na.pass) #OK
pacf(residuals_CETE,na.action = na.pass) #OK

acf(residuals_HOMA,na.action = na.pass) #OK
pacf(residuals_HOMA,na.action = na.pass) #OK

acf(residuals_LEMA,na.action = na.pass) #OK
pacf(residuals_LEMA, na.action = na.pass) #OK

acf(residuals_PAIN,na.action = na.pass) #OK
pacf(residuals_PAIN,na.action = na.pass) # NO OK

acf(residuals_POMA,na.action = na.pass) #OK
pacf(residuals_POMA,na.action = na.pass) #OK

acf(residuals_SASO,na.action = na.pass) #OK
pacf(residuals_SASO,na.action = na.pass) # OK
