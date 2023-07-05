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
###### FIT GLMM #############
#############################
#############################

#  Model that corrects temporal correlation and random groups for the previous best model

lmer_simple_fit <- lme(log(log10_ratio_tplus1_over_t) ~ exclusion_ratio_prev_year, 
                       random= ~ 1|species,
                       data=data_growth%>% mutate(species=as.factor(species)),
                       correlation=corCAR1(form= ~ year|species)
)


# Plot (fixed effect) results for GLMM

lmer_simple_fit_predictions <- tibble(real_data = data_growth$log10_ratio_tplus1_over_t,
                                      predicted_data = exp(predict(lmer_simple_fit, data_growth)),
                                      residuals = residuals(lmer_simple_fit),
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


png("Images/figS17.png",
    width = 11.69*1.5, # The width of the plot in inches
    height = 11.69*0.8, units = "in", res=300*2)
(population_growth_plants_lmer+
    ylim(c(-.1,2.5))) + (population_growth_boxplots_lmer+
                           ylim(c(-.1,2.5)))
dev.off()