
library(tidyverse)
library(CValternatives)
library(patchwork)
library(scales)
library(latex2exp)
library(anisoFun)

# Load data

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


# SPECIES METRICS -----------

png("Images/evolution_isoprobability_of_exclusion.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*0.6, units = "in", res=300*2)

iso_plot <- ggplot(caracoles_plot_year_final %>% select(year, iso_prob_exclusion) %>%
         distinct(),
       aes(x = as.factor(year), y = iso_prob_exclusion, group=1))+
  geom_line()+
  labs(y="Isoprobability of exclusion",x="Year")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sp_plot <- ggplot(caracoles_plot_year_final %>% select(year, number_plant_sp) %>%
                     distinct(),
                   aes(x = as.factor(year), y = number_plant_sp, group=1))+
  geom_line()+scale_y_continuous(breaks= pretty_breaks())+
  labs(y="Number plant species in the metacommunity",x="Year")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

iso_plot + sp_plot

dev.off()

png("Images/evolution_exclsuion_rates_paper.png",
    width = 11.69*0.8, # The width of the plot in inches
    height = 11.69*0.8*0.5, units = "in", res=300*2)

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


ggplot(caracoles_plot_year_final %>% filter(metric == "prob_excl_rep",
                                            species %in% c("BEMA","CETE", "HOMA",
                                                           "LEMA","PAIN","POMA",
                                                           "SASO")) %>%
         mutate(exclusion_rate = value/iso_prob_exclusion),
       aes(x = as.factor(year),y = exclusion_rate))+
  geom_boxplot()+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.3)+
  stat_summary(fun = mean, geom="point", shape=23, size=3,fill="red")+
  geom_hline(yintercept = 1, color = "deepskyblue", linetype = "dashed", size = 1.3)+
  facet_wrap(~species, ncol = 4, labeller = as_labeller(species_names))+ #, scales = "free")+
  # scale_y_continuous(labels = scientific)+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  ylab("Exclusion ratio")+
  labs(x=NULL)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text = element_text(face = "italic"))

dev.off()

png("Images/evolution_arc_distances_paper.png",
    width = 11.69*0.8, # The width of the plot in inches
    height = 11.69*0.8*0.5, units = "in", res=300*2)

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


arc_distance <- read_csv("Data/caracoles_processed_data/caracoles_arc_distance.csv")

ggplot(caracoles_plot_year_final %>% filter(metric == "prob_excl_rep",
                                            species %in% c("BEMA","CETE", "HOMA",
                                                           "LEMA","PAIN","POMA",
                                                           "SASO")) %>%
         left_join(arc_distance, by = c("year","species")),
       aes(x = as.factor(year),y = arc_distance))+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.3)+
  stat_summary(fun = mean, geom="point", shape=23, size=3,fill="red")+
  facet_wrap(~species, ncol = 4, labeller = as_labeller(species_names))+ #, scales = "free")+
  # scale_y_continuous(labels = scientific)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Exclusion distance")+
  labs(x=NULL)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text = element_text(face = "italic"))

dev.off()


# COMMUNITY METRICS -----------

list_years <- caracoles_plot_year_final_aux$year %>% unique()

anisotropy_agg_info <- tibble(year = list_years)
small_omega_mean <- NA
small_omega_lowerCI <- NA
small_omega_upperCI <- NA
anisotropy_index_mean <- NA
anisotropy_index_lowerCI <- NA
anisotropy_index_upperCI <- NA
PV_index_mean <- NA
PV_index_lowerCI <- NA
PV_index_upperCI <- NA


for (year_i in 1:length(list_years)) {
  
  data_year_i <- caracoles_plot_year_final_aux %>% filter(year == list_years[year_i])
  
  Omega_rep <- data_year_i[1,(5+number_boot_replicates):(4+2*number_boot_replicates)] %>%
    as.numeric()
  small_omega_rep <- Omega_rep^(1/data_year_i$number_plant_sp[1])
  anisotropy_rep <- NULL
  PV_rep <- NULL
  
  for(i in 1:number_boot_replicates){
    
    probability_vector <- data_year_i[,(4+i)] %>% pull()
    
    anisotropy_rep <-  c(anisotropy_rep, anisotropy_index(probability_vector))
    PV_rep <- c(PV_rep, PV(probability_vector))
    
  }
  
  small_omega_rep_CI <- quantile(small_omega_rep, prob=c(.025,.975)) %>% as.numeric()
  anisotropy_rep_CI <- quantile(anisotropy_rep, prob=c(.025,.975)) %>% as.numeric()
  PV_rep_CI <- quantile(PV_rep, prob=c(.025,.975)) %>% as.numeric()
  
  anisotropy_agg_info$small_omega_mean[year_i] <- mean(small_omega_rep, na.rm = T)
  anisotropy_agg_info$small_omega_lowerCI[year_i] <- small_omega_rep_CI[1]
  anisotropy_agg_info$small_omega_upperCI[year_i] <- small_omega_rep_CI[2]
  anisotropy_agg_info$anisotropy_index_mean[year_i] <- mean(anisotropy_rep, na.rm = T)
  anisotropy_agg_info$anisotropy_index_lowerCI[year_i] <- anisotropy_rep_CI[1]
  anisotropy_agg_info$anisotropy_index_upperCI[year_i] <- anisotropy_rep_CI[2]
  anisotropy_agg_info$PV_index_mean[year_i] <- mean(PV_rep, na.rm = T)
  anisotropy_agg_info$PV_index_lowerCI[year_i] <- PV_rep_CI[1]
  anisotropy_agg_info$PV_index_upperCI[year_i] <- PV_rep_CI[2]
  
}

all_anisotropy_agg_info <- anisotropy_agg_info %>% left_join(number_plant_sp_year, by="year") %>%
  mutate(capital_omega_mean=small_omega_mean^number_plant_sp,
         capital_omega_lowerCI=small_omega_lowerCI^number_plant_sp,
         capital_omega_upperCI=small_omega_upperCI^number_plant_sp)

cor.test(all_anisotropy_agg_info$small_omega_mean,all_anisotropy_agg_info$number_plant_sp)
cor.test(all_anisotropy_agg_info$small_omega_mean,all_anisotropy_agg_info$anisotropy_index_mean, method = "spearman")
cor.test(all_anisotropy_agg_info$small_omega_mean,all_anisotropy_agg_info$anisotropy_index_mean)

png("Images/asimmetry index.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)

ggplot(data = anisotropy_agg_info,
       aes(x = as.factor(year),
           y = anisotropy_index_mean,
           color = as.factor(year))) + 
  scale_color_manual(values=c("#E69F00", "#009E73","#56B4E9","#0072B2","#D55E00","#CC79A7","#000000"))+
  geom_point(size=1.5) + 
  geom_errorbar(aes(ymin = anisotropy_index_lowerCI,ymax = anisotropy_index_upperCI)) + 
  # scale_color_brewer(palette = "Dark3")+
  labs(y="Feasibility asymmetry index, J'(A)", x=NULL, color = NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+theme(legend.text=element_text(size=15))+ 
  guides(colour = guide_legend(override.aes = list(size=1.5)))
dev.off()

png("Images/Small omega.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)

ggplot(data = anisotropy_agg_info,aes(x = as.factor(year),
                                      y = small_omega_mean, 
                                      color = as.factor(year))) + 
  scale_color_manual(values=c("#E69F00", "#009E73","#56B4E9","#0072B2","#D55E00","#CC79A7","#000000"))+
  geom_point(size=1.5) + 
  geom_errorbar(aes(ymin = small_omega_lowerCI,ymax = small_omega_upperCI)) + 
  # scale_color_brewer(palette = "Dark3")+
  labs(y="Feasibility domain's relative volume per species", x=NULL, color = NULL)+
  theme_bw()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+theme(legend.text=element_text(size=15))+ 
  guides(colour = guide_legend(override.aes = list(size=1.5)))

dev.off()