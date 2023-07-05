
library(tidyverse)
library(scales)

data_experiment <- read_csv("Results/simulations/data_experiment_LV_times.csv")

png("Images/figS8.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)

ggplot(data_experiment, aes(x=capital_omega,y=time_LV_unitary))+
  geom_smooth(method = "lm", color = "deepskyblue", fill = "deepskyblue")+
  geom_point(size=2)+
  theme_bw()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Feasibility domain's relative volume [adimensional]",
       y="Average time to find an intrinsic\ngrowth vector in the feasibility domain [seconds]")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

dev.off()
#############