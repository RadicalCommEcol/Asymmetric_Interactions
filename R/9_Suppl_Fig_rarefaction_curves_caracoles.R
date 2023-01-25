
library(tidyverse)
library(vegan)
library(iNEXT)

# Load caracoles data

interactions_year_raw <- read_csv("Data/caracoles_raw_data/competition.csv") %>% 
  select(-comments, -observers, -fruit, -seed,-FRPU,-SUSP,-COSQ,-CRCR,-DA,
         -MEPO,-RAPE,-ANAR,-ACHI,-ARTE,-NEWGRASS)

focal <- c("PAIN","PUPA","POMO","SCLA","HOMA","CETE","POMA","MESU","SPRU","CHMI","MEEL","LEMA","LYTR","PLCO","BEMA","SUSP","SASO","SOAS","CHFU")#interactions_year_raw$focal %>% unique()
counts_colnames <- colnames(interactions_year_raw)
common_colnames <- counts_colnames[1:5]
focal_colnames <- counts_colnames[7:length(counts_colnames)]

interactions_caracoles <- NULL

for(focal.i in focal){
  
  interactions_focal_i <- interactions_year_raw %>% filter(focal==focal.i) %>%
    select(-focal)
  colnames(interactions_focal_i) <- c(common_colnames,paste0(focal.i,"-",focal_colnames))
  interactions_caracoles <- bind_rows(interactions_caracoles,interactions_focal_i)
  
}

#################
# RAREFACTION CURVES

x <- interactions_caracoles %>% select(-day, -month, -plot, -subplot) %>% 
  group_by(year) %>% summarise_all(sum,na.rm = TRUE)

col_names_x <- paste0("year ",x$year)

x <- x %>% ungroup() %>% select(-year)

rownames(x) <- col_names_x

x.list <- setNames(split(x, seq(nrow(x))), rownames(x))
str(x.list)

year_i = 1
out1 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction1 <- ggiNEXT(out1, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2015")
rarefaction1


year_i = 2
out2 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction2 <- ggiNEXT(out2, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2016")
rarefaction2

year_i = 3
out3 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction3 <- ggiNEXT(out3, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2017")
rarefaction3


year_i = 4
out4 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction4 <- ggiNEXT(out4, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2018")
rarefaction4


year_i = 5
out5 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction5 <- ggiNEXT(out5, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2019")
rarefaction5


year_i = 6
out6 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction6 <- ggiNEXT(out6, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2020")
rarefaction6


year_i = 7
out7 <- iNEXT(unlist(x.list[[year_i]]), q=0, datatype="abundance", endpoint=100000, nboot=1000)
rarefaction7 <- ggiNEXT(out7, type=2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 18) +
  theme(legend.position="none")+
  labs(x="Number of interactions",
       title="2021")
rarefaction7

library(patchwork)

png("Images/plant_plant_rarefaction_curves.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69*0.842, units = "in", res=300*2)

(rarefaction1|rarefaction2|rarefaction3)/(rarefaction4|rarefaction5|rarefaction6)/(plot_spacer()|rarefaction7|plot_spacer())

dev.off()
