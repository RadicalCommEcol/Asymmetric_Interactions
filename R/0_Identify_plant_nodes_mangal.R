
library("taxize")
library(tidyverse)

matrix_nodes_raw <- read_csv2("Data/mangal_raw_data/mangal_pollination_quant_nodes.csv")
matrix_links_raw <- read_csv2("Data/mangal_raw_data/mangal_pollination_quant_links.csv")
network_raw <- read_csv2("Data/mangal_raw_data/mangal_pollination_quant_nets.csv")


# Number of networks
matrix_links_raw$network_id %>% unique()

# Remove networks with nodes without a toxonomy name

matrix_nodes_raw$node_id[is.na(matrix_nodes_raw$taxonomy_name)]

matrix_nodes_without_taxonomy <- matrix_links_raw %>%
  filter(node_from == "mangal_node_0090" | node_to == "mangal_node_0090") %>% 
  select(network_id) %>% pull() %>% unique()

# We study networks where all the nodes have their taxonomy name

matrix_links_processed <- matrix_links_raw %>% filter(! network_id %in% matrix_nodes_without_taxonomy)

matrix_links_processed$network_id %>% unique() %>% length() # 88 networks

matrix_nodes_selected <- unique(c(matrix_links_processed$node_from,matrix_links_processed$node_to))
matrix_nodes_processed <- matrix_nodes_raw %>% 
  filter(node_id %in% matrix_nodes_selected) # 339 nodes selected

all(matrix_nodes_selected %in% matrix_nodes_processed$node_id )

write_csv(matrix_nodes_processed, "Data/mangal_processed_data/mangal_pollination_quant_nodes.csv")
write_csv(matrix_links_processed, "Data/mangal_processed_data/mangal_pollination_quant_links.csv")


# nodes classification in plantae and animalia

matrix_nodes_processed$kingdom <- NA

for(i in 1:nrow(matrix_nodes_processed)){
  
  species_i <- matrix_nodes_processed$taxonomy_name[i]
  
  if(!is.na(species_i)){
    
    classification_results <- classification(species_i, db = 'itis')
    
    if(!is.na(classification_results[[species_i]])){
      matrix_nodes_processed$kingdom[i] <- classification_results[[species_i]][["name"]][1]
      cat(matrix_nodes_processed$kingdom[i],"\n","\n","\n")
    }
  }
  
}

names_2_review <- matrix_nodes_processed$taxonomy_name[!matrix_nodes_processed$kingdom %in% c("Animalia", "Plantae")]

for(i in 276:1){
  
  species_i <- matrix_nodes_processed$taxonomy_name[i]
  
  if((is.na(matrix_nodes_processed$kingdom[i]) | (!matrix_nodes_processed$kingdom[i] %in% c("Animalia", "Plantae")))& !is.na(species_i)){
    
    classification_results <- classification(species_i, db = 'eol')
    
    if(!is.na(classification_results[[species_i]])){
      matrix_nodes_processed$kingdom[i] <- classification_results[[species_i]][["name"]][1]
      cat(matrix_nodes_processed$kingdom[i],"\n","\n","\n")
    }
  }
  
}

matrix_nodes_processed$kingdom[matrix_nodes_processed$kingdom %in% c("Insecta","Arthropoda",
                                                         "Coleoptera",)] <- "Animalia"

matrix_nodes_processed$kingdom[matrix_nodes_processed$kingdom %in% c("Dipsacaceae")] <- "Plantae"

write_csv(matrix_nodes_processed, "Data/mangal_processed_data/mangal_pollination_quant_nodes.csv")


names_2_review <- matrix_nodes_processed$taxonomy_name[!matrix_nodes_processed$kingdom %in% c("Animalia", "Plantae")]







matrix_nodes_processed$kingdom[matrix_nodes_processed$taxonomy_name %in% c("Cheilosia aenea","Melissodes lupina",
                                                                           "Polistes olivaceus","Eristalinus vicarians",
                                                                           "Cratopus aurostriatus","Oedemera lurida",
                                                                           "Oedemera flavipes","Anthidium sticticum",
                                                                           "Psilothrix viridicoerulea","Apis melifera",
                                                                           "Polistes gallicus","Chalicodoma pyrenaica",
                                                                           "Milichidae","Mordella bipunctata",
                                                                           "Helorus coruscus","Oxythyrea funesta",
                                                                           "Callophrys rubi","Mordellistena pumila",
                                                                           "Dasytidae","Oedemera nobilis",
                                                                           "Scolia hirta","Ceratina curcurbitina",
                                                                           "Oedemera barbara","Camponotus aetiops"
                                                                           )] <- "Animalia"





matrix_nodes_processed$kingdom[matrix_nodes_processed$taxonomy_name %in% c("Galactites tomentosa","Polygonum arenastrum",
                                                                           "Grindelia camporum","Aster chilensis",
                                                                           "Croton setigerus","Dillenia ferruginea",
                                                                           "Sonchus tenerrinus","Aetheorrhiza bulbosa",
                                                                           "Alyssum maritimum","Thymelaea hirsuta",
                                                                           "Carpobrotus affinis acinaciformis",
                                                                           "Echium sabulicola","Daucus carota carota",
                                                                           "Sedum sediforme","Helianthemum guttatum",
                                                                           "Psoralea bituminosa","Artemesia")] <- "Plantae"


# Sanity check
matrix_nodes_processed$taxonomy_name[!matrix_nodes_processed$kingdom %in% c("Animalia", "Plantae")]


write_csv(matrix_nodes_processed, "Data/mangal_processed_data/mangal_pollination_quant_nodes.csv")

#####################################
# BIPARTITE NETWORKS??

matrix_links_processed_kingdom <- matrix_links_processed %>%
  left_join(matrix_nodes_processed %>% rename(node_from = node_id,
                                              taxonomy_name_from =taxonomy_name,
                                              taxonomy_rank_from = taxonomy_rank,
                                              kingdom_from = kingdom),
            by = "node_from") %>%
  
  left_join(matrix_nodes_processed %>% rename(node_to = node_id,
                                              taxonomy_name_to =taxonomy_name,
                                              taxonomy_rank_to = taxonomy_rank,
                                              kingdom_to = kingdom),
            by = "node_to")

matrix_links_processed_kingdom  %>% filter(kingdom_to == kingdom_from)

# It seems that networks are bipartite

networks_with_NA_NA_int <- matrix_links_raw %>%
  filter(node_from == "mangal_node_0090" & node_to == "mangal_node_0090") %>% 
  select(network_id) %>% pull() %>% unique()

# We study networks where all the nodes have their taxonomy name

matrix_links_processed_2 <- matrix_links_raw %>% filter(!network_id %in% unique(c(networks_with_NA_NA_int,
                                                                                  unique(matrix_links_processed$network_id))))

matrix_links_processed_2$network_id %>% unique() %>% length() # 209 networks

matrix_nodes_selected_2_aux <- unique(c(matrix_links_processed_2$node_from,matrix_links_processed_2$node_to))
matrix_nodes_selected_2 <- matrix_nodes_selected_2_aux[!matrix_nodes_selected_2_aux %in% matrix_nodes_selected]
matrix_nodes_processed_2 <- matrix_nodes_raw %>% 
  filter(node_id %in% matrix_nodes_selected_2) # 339 nodes selected

all(matrix_nodes_selected_2 %in% matrix_nodes_processed_2$node_id )
all(matrix_nodes_processed_2$node_id %in% matrix_nodes_selected_2)

write_csv(matrix_nodes_processed_2, "Data/mangal_processed_data/mangal_pollination_quant_nodes_2.csv")
write_csv(matrix_links_processed_2, "Data/mangal_processed_data/mangal_pollination_quant_links_2.csv")


# nodes classification in plantae and animalia

matrix_nodes_processed_2$kingdom <- NA

for(i in 1:nrow(matrix_nodes_processed_2)){
  
  species_i <- matrix_nodes_processed_2$taxonomy_name[i]
  
  if(!is.na(species_i)){
    
    classification_results <- classification(species_i, db = 'itis')
    
    if(!is.na(classification_results[[species_i]])){
      matrix_nodes_processed_2$kingdom[i] <- classification_results[[species_i]][["name"]][1]
      cat(matrix_nodes_processed_2$kingdom[i],"\n","\n","\n")
    }
  }
  
}

names_2_review_2 <- matrix_nodes_processed_2$taxonomy_name[!matrix_nodes_processed_2$kingdom %in% c("Animalia", "Plantae")]

# USE GBIF, Brazil, Japan, Findland, UK classification to retrieve the kingdom

for(i in nrow(matrix_nodes_processed_2):1){
  
  species_i <- matrix_nodes_processed_2$taxonomy_name[i]
  
  if((is.na(matrix_nodes_processed_2$kingdom[i]) | (!matrix_nodes_processed_2$kingdom[i] %in% c("Animalia", "Plantae")))& !is.na(species_i)){
    
    classification_results <- classification(species_i, db = 'eol')
    
    if(!is.na(classification_results[[species_i]])){
      matrix_nodes_processed_2$kingdom[i] <- classification_results[[species_i]][["name"]][1]
      cat(matrix_nodes_processed_2$kingdom[i],"\n","\n","\n")
    }
  }
  
}

matrix_nodes_processed_2$kingdom[matrix_nodes_processed_2$kingdom %in% c("Insecta","Arthropoda",
                                                                     "Coleoptera")] <- "Animalia"

matrix_nodes_processed_2$kingdom[matrix_nodes_processed_2$kingdom %in% c("Dipsacaceae")] <- "Plantae"

write_csv(matrix_nodes_processed_2, "Data/mangal_processed_data/mangal_pollination_quant_nodes_2.csv")


names_2_review_2 <- matrix_nodes_processed_2$taxonomy_name[!matrix_nodes_processed_2$kingdom %in% c("Animalia", "Plantae")]







matrix_nodes_processed_2$kingdom[matrix_nodes_processed_2$taxonomy_name %in% c("Camptocladius brachytomus","Orthocentrinae",
                                                                               "Bradysia giraudi","Rhamphodon naevius",
                                                                               "Micropeza thoracicum","Eristalix tenax",
                                                                               "Anothomyia pluvialis","Anthrax tinctus",
                                                                               "Copestylum pusillum","Homoneura quadrivitta",
                                                                               "Polistes hebraeus","Campopleginae",
                                                                               "Spodoptera littoralis","Sarrothripinae",
                                                                               "Thysanoplusia orichalcea","Parasarcophaga hirtipes",
                                                                               "Simosyrphus aegyptius","Plagiostenopterina cyanosoma",
                                                                               "Physiphora azurea","Pachycerina seychellensis",
                                                                               "Tricimba seychellensis","Oedemeronia lucidicollis",
                                                                               "Oxycetonia jucunda","Helophilus virgatus",
                                                                               "Dalopius tamui","Megacampsomeris grossa",
                                                                               "Coptotermes gambrinus","Macrolagria robusticeps",
                                                                               "Trypherus niponicus", "Fischeria")] <- "Animalia"





matrix_nodes_processed_2$kingdom[matrix_nodes_processed_2$taxonomy_name %in% c("Acanthostachys strobilaceae","Erigeron tenuifolius",
                                                                               "Lotus purshianus","Eremocarpus setigerus",
                                                                               "Nephrosperma vanhoutteanum","Ochna kirkii",
                                                                               "Patrinia scabiosaefolia","Artemesia")] <- "Plantae"


# Sanity check
matrix_nodes_processed_2$taxonomy_name[!matrix_nodes_processed_2$kingdom %in% c("Animalia", "Plantae")]


write_csv(matrix_nodes_processed_2, "Data/mangal_processed_data/mangal_pollination_quant_nodes_2.csv")

matrix_nodes_processed_12 <- rbind(matrix_nodes_processed,matrix_nodes_processed_2)

matrix_links_processed_kingdom_2 <- matrix_links_processed_2 %>%
  left_join(matrix_nodes_processed_12 %>% rename(node_from = node_id,
                                              taxonomy_name_from =taxonomy_name,
                                              taxonomy_rank_from = taxonomy_rank,
                                              kingdom_from = kingdom),
            by = "node_from") %>%
  
  left_join(matrix_nodes_processed_12 %>% rename(node_to = node_id,
                                              taxonomy_name_to =taxonomy_name,
                                              taxonomy_rank_to = taxonomy_rank,
                                              kingdom_to = kingdom),
            by = "node_to")

# Sanity check
bad_results_check1 <- matrix_links_processed_kingdom_2  %>% filter(kingdom_to == kingdom_from)
bad_results_check2 <- matrix_links_processed_kingdom_2  %>% filter(is.na(kingdom_from),is.na(kingdom_to))

NA_results <- matrix_links_processed_kingdom_2 %>% filter(is.na(kingdom_from) | is.na(kingdom_to))



# Save links data for 88 + 209 networks

matrix_links_processed_kingdom_2_review <- matrix_links_processed_kingdom_2

for (i in 1:nrow(matrix_links_processed_kingdom_2_review)){
  
  if(is.na(matrix_links_processed_kingdom_2_review$kingdom_to[i]) & 
     (matrix_links_processed_kingdom_2_review$kingdom_from[i] == "Plantae")){
    
    matrix_links_processed_kingdom_2_review$kingdom_to[i] <- "Animalia"
    matrix_links_processed_kingdom_2_review$taxonomy_name_to[i] <- "Unknown Animalia"
    
  }else if(is.na(matrix_links_processed_kingdom_2_review$kingdom_to[i]) & 
           (matrix_links_processed_kingdom_2_review$kingdom_from[i] == "Animalia")){
    
    matrix_links_processed_kingdom_2_review$kingdom_to[i] <- "Plantae"
    matrix_links_processed_kingdom_2_review$taxonomy_name_to[i] <- "Unknown Plantae"
    
  }else if(is.na(matrix_links_processed_kingdom_2_review$kingdom_from[i]) & 
           (matrix_links_processed_kingdom_2_review$kingdom_to[i] == "Plantae")){
    
    matrix_links_processed_kingdom_2_review$kingdom_from[i] <- "Animalia"
    matrix_links_processed_kingdom_2_review$taxonomy_name_from[i] <- "Unknown Animalia"
    
  }else if(is.na(matrix_links_processed_kingdom_2_review$kingdom_from[i]) & 
           (matrix_links_processed_kingdom_2_review$kingdom_to[i] == "Animalia")){
    
    matrix_links_processed_kingdom_2_review$kingdom_from[i] <- "Plantae"
    matrix_links_processed_kingdom_2_review$taxonomy_name_from[i] <- "Unknown Plantae"
  }
  
}

resulting_matrix_links_processed_kingdom <- rbind(matrix_links_processed_kingdom,
                                                  matrix_links_processed_kingdom_2_review)

write_csv(matrix_links_processed_kingdom, 
          "Data/mangal_processed_data/NO_NAs_mangal_links_processed_kingdom.csv")

write_csv(resulting_matrix_links_processed_kingdom, 
          "Data/mangal_processed_data/ALL_mangal_links_processed_kingdom.csv")
