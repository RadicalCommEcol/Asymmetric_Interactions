library(tidyverse)
library(scales)
library(latex2exp)

png("Images/fig4b.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.5, units = "in", res=300*2)


expectations <- tibble(x = 5+c(5,10,5,10),
                       y = 5+c(7.5,10,7.5,5),
                       name=c("up","up","down","down"))

ggplot(expectations, aes(x=x, y=y, color = name))+
  geom_point(size = 8, alpha =.5, linetype ="dashed")+
  geom_line(size = 2, alpha =.5, linetype ="dashed")+
  scale_color_discrete(labels = unname(TeX(c("$\\;\\;\\;\\;\\theta_A\\uparrow\\uparrow,\\;P_A^E\\downarrow\\downarrow,\\;ER_A\\downarrow\\downarrow$",
                                             "$\\;\\;\\;\\;\\theta_A\\downarrow\\downarrow,\\;P_A^E\\uparrow\\uparrow,\\;ER_A\\uparrow\\uparrow$"))))+
  theme_bw()+
  labs(x="Time",y="Expected population growth",color=NULL)+
  theme(legend.position="none")+
  annotate("text", x=12.5, y=15.5, size= 8, colour = "cyan3", label=TeX(c("$\\;\\;\\;\\;\\theta_A\\uparrow\\uparrow,\\;P_A^E\\downarrow\\downarrow,\\;ER_A\\downarrow\\downarrow$")))+
  annotate("text", x=12.5, y=9.5, size= 8, colour = "salmon", label=TeX(c("$\\;\\;\\;\\;\\theta_A\\downarrow\\downarrow,\\;P_A^E\\uparrow\\uparrow,\\;ER_A\\uparrow\\uparrow$")))+
  ylim(9,15.6)+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+ 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_x_continuous(breaks=c(10,15),
                     labels=c(TeX("$t$"), TeX("$t+\\Delta t$")))+
  theme(axis.text.x = element_text(size = 22))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


dev.off()
