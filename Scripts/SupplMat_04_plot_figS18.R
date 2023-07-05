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
library(latex2exp)
library(investr)
library(patchwork)

abundance_model_data <- read_csv("Results/abundance_model_data.csv") %>%
  mutate(year=year,species = as.factor(species))

# Create data frame with input data
data = abundance_model_data %>% mutate(ratio_tplus1_over_t=total_individuals/total_individuals_prev_year,
                                       difference_tplus1_minus_t=total_individuals-total_individuals_prev_year,
                                       number_species_prev_year = round(exclusion_ratio_prev_year/prob_excl_prev_year,0))
# Update data for exponential fit
data_growth <- data  %>% ungroup() %>% mutate(log10_ratio_tplus1_over_t = log10(ratio_tplus1_over_t+1))%>% 
  rename(exclusion_distance_prev_year = arc_distance_prev_year) %>% 
  dplyr::select("year","log10_ratio_tplus1_over_t", "total_individuals_prev_year",
                "prob_excl_prev_year","exclusion_ratio_prev_year",
                "exclusion_distance_prev_year","species") %>% mutate(number_species_prev_year= exclusion_ratio_prev_year/prob_excl_prev_year)

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


#############################
#############################
###### FIT MODELS EXP NEG ###
#############################
#############################

fit_excl_dist_exp_neg <- nls(log10_ratio_tplus1_over_t ~ c*exp(exclusion_distance_prev_year*d), 
                             data = data_growth,
                             start = list(c = 0.5, d=-1))


# Plot results for exclusion distance exponential NLS fit

nlm_excl_distance_exp_neg_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                                predicted_data = (predict(fit_excl_dist_exp_neg, data_growth)),
                                                residuals = residuals(fit_excl_dist_exp_neg),
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


png("Images/figS18.png",
    width = 11.69*1.5, # The width of the plot in inches
    height = 11.69*0.8, units = "in", res=300*2)
(population_growth_plants_exclusion_distance+
    ylim(c(-.1,2.5))) + (population_growth_boxplots_exclusion_distance+
                           ylim(c(-.1,2.5)))
dev.off()
