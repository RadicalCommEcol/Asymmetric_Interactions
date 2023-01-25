library("rmangal")
library(tidyverse) 

search_datasets(query = "908")
search_datasets(list(id = ))


all_datasets <- search_datasets("", verbose = FALSE)
all_datasets$id_network <- NA

id_and_id_network <- NULL

for(i in 1:nrow(all_datasets)){
  
  if(any(!is.na( all_datasets[[10]][[i]]))){
    
    id_and_id_network_aux <- tibble(id=all_datasets$id[i],
                                    id_network = all_datasets[[10]][[i]]$id)
    
    id_and_id_network <- bind_rows(id_and_id_network,id_and_id_network_aux)
  }
  
  
}

mangal_data <- read_csv("Data/mangal_processed_data/NO_NAs_mangal_links_processed_kingdom.csv")

id_network_mangal <- gsub("[^0-9.-]", "", mangal_data$network_id) %>% as.numeric()
id_and_id_network_mangal <- id_and_id_network %>% filter(id_network %in% id_network_mangal)

number_networks_in_study <- id_and_id_network_mangal %>% group_by(id) %>% count() %>%
  rename(number_networks = n)

id_mangal <- unique(id_and_id_network_mangal$id)

mangal_datasets <- all_datasets[all_datasets$id %in% id_mangal,]

mangal_table <- NULL

for(i.mangal in 1:nrow(mangal_datasets)){
  
  mangal_table_aux <- mangal_datasets[[9]][[i.mangal]]
  mangal_table_aux$
  mangal_table <- bind_rows(mangal_table, mangal_table_aux)
  
}

mangal_table$number_networks <- number_networks_in_study$number_networks
mangal_table_final <- mangal_table[,c("first_author","number_networks","year","doi")]
mangal_table_final
