library(matlib) # to multiply matrices
library(tidyverse)
library(moments)
library(RColorBrewer)

source("R/matrix_year_i.R")


# We create the metaweb for each year

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


selected_species <- c("BEMA","CETE", "HOMA", "LEMA","PAIN","POMA", "SASO","SCLA")

years_included <- matrix_entries_raw$year %>% unique()

A_int_elements <- NULL

for(i in 1:length(years_included)){
  
  year_i <- years_included[i]
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  A_int_elements_aux <- tibble(elements = as.numeric(A_int), Year = year_i, type ="All\nspecies")
  
  columns_with_elements_selected_species <- which(colnames(A_int) %in% selected_species)
  
  A_int_elements_aux_selected <- tibble(elements = as.numeric(A_int[,columns_with_elements_selected_species]), 
                                        Year = year_i, type ="Selected\nspecies")
  
  A_int_elements_aux <- bind_rows(A_int_elements_aux, A_int_elements_aux_selected)
  
  A_int_elements <- bind_rows(A_int_elements,A_int_elements_aux)
  
  
  cat(skewness(as.numeric(A_int)), "\n")
  cat(skewness(A_int_elements_aux %>% select(elements) %>% pull %>% abs()), "\n")
  
}

png("Images/skewness.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*0.6, units = "in", res=300*2)

ggplot(A_int_elements, aes(x= abs(elements)))+
  geom_density(aes(fill= as.factor(type)))+ scale_fill_manual(values=c("#E69F00", "#009E73"))+
  facet_wrap(~Year+type, nrow = 3, ncol=6, scale="free")+
  theme_bw()+
  labs(x="| Interaction strength |",y="Density")+
  theme(legend.position="bottom")+theme(legend.position = "none")#+
  #theme(strip.background =element_rect(fill="grey"))


dev.off()