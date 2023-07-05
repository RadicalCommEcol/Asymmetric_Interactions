library(tidyverse)
library(latex2exp)
library(scales)
#################################################################################
################################################################################


results <- read_csv("Results/test_asymmetry_interactions_FD_201_omega_rep_55000.csv")

facet_names <- c(
  `3` = "3 species",
  `6` = "6 species",
  `10` = "10 species"
)

results_plot <- results %>% mutate(big_omega=small_omega_mean^richness)
results_plot$type <- "Platykurtic"
results_plot$type[results_plot$nu==1] <- "Mesokurtic"
results_plot$type[results_plot$nu<1] <- "Leptokurtic"

df <- group_by(results_plot, richness, type)
df_isotropy <- summarise(df, avg = median(anisotropy_index_mean))
df_isotropy

df_capital_omega <- summarise(df, avg = median(big_omega))
df_capital_omega


png("Images/figS1.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.45, units = "in", res=300*2)
ggplot(results_plot, aes(x=as.factor(type), y = anisotropy_index_mean,color=as.factor(type))) +
  geom_boxplot(alpha=.3,size=1)+
  geom_hline(aes(yintercept=1),  color="deepskyblue1",
             linetype="dashed", size=1)+
  labs(y="Asymmetry index, J'(A)",x=NULL)+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+ 
  theme(legend.position="none") +
  theme(legend.title=element_blank())
dev.off()



png("Images/figS2.png",
    width = 11.69*0.72, # The width of the plot in inches
    height = 11.69*0.45, units = "in", res=300*2)
ggplot(results_plot, aes(x=as.factor(type), y = big_omega,color=as.factor(type))) +
  geom_boxplot(alpha=.3,size=1)+
  # geom_hline(aes(yintercept=0.5),  color="deepskyblue1",
             # linetype="dashed", size=1)+
  ylab(TeX("Feasibility domain's size:  $\\Omega$"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  # labs(y="Feasibility domain's size",x=NULL)+
  xlab(NULL)+
  facet_wrap(richness ~ ., labeller = as_labeller(facet_names))+
  theme_bw()+theme(
    axis.text.x= element_text(color="black",angle=45, size=10, vjust=0.5),
    axis.text.y= element_text(color="black", size=12, vjust=0.5),
    axis.title.y = element_text(color="black",size=12, vjust=0.5),
    plot.title = element_text(color="black",face="bold",size=14, hjust=0.5,vjust=1)
  )+ 
  theme(legend.position="none") +
  theme(legend.title=element_blank())
dev.off()


