library(cowplot)
#library(tidyverse)
library(dplyr)
library(ggplot2)
library(geomtextpath)

# Read in data for bar plot
data = read.csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/5_plot_B_inputs/evidence_B.csv", sep="\t" , header=TRUE)

# Separate into the four datasets for the figure
CS_S_main <- data %>%
  filter(platforms=="Clingen_vs_SEQ" & orientation == "up" & sub=="main")

CS_S <- data %>%
  filter(platforms=="Clingen_vs_SEQ" & orientation == "up")

CS_C_main <- data %>%
  filter(platforms=="Clingen_vs_SEQ" & orientation == "down" & sub=="main")

CS_C <- data %>%
  filter(platforms=="Clingen_vs_SEQ" & orientation == "down")

# Order terms
## CS_C
CS_C_main_order_levels <- CS_C_main %>%
  arrange(codes) %>%
  pull(codes)

CS_C$codes <- factor(CS_C$codes, levels = (CS_C_main_order_levels))

## CS_S
CS_S_main_order_levels <- CS_C_main %>%
  arrange(codes) %>%
  pull(codes)

CS_S$codes <- factor(CS_S$codes, levels = (CS_S_main_order_levels))



##CS_C graph (reverse oriented)
p1 <- ggplot(CS_C, aes(fill=factor(sub,levels=c("total","main")) , y=variant_count, x=codes, color=sub)) +
  geom_bar(position="stack", stat="identity")+
  
  scale_colour_manual(values = c("black",adjustcolor( "grey95", alpha.f = 0.00001))) +
  
  scale_fill_manual(values=c('grey95', '#baffc9')) +
  
  geom_text(aes(label = main, y = main+225), hjust = -0.3,size = 3, color="darkgreen") +
  geom_text(aes(label = total_label, y = total+200), hjust = 1.3,size = 3, color="grey")+
  geom_text(aes(label = freq, y = 3200),size = 3.5, color="red")+
  
  coord_flip(xlim= NULL, ylim = rev(c(0, 3500)), expand=FALSE, clip = "on")+
  ggtitle("Clingen Specific") +
  xlab("Evidence Code") + ylab("Variant Counts")+
  scale_x_discrete(name = "",position = "top",expand=c(0,0)) +
  scale_y_reverse(expand=c(0, 0),breaks = c(0, 500, 1000, 1500,2000,2500,3000,3500))


p1 <- p1 + theme(legend.box="horizontal",
                 legend.key=element_blank(),
                 legend.title=element_blank(),
                 legend.position="none",
                 panel.border = element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.length.y = unit(0.5, "cm"),
                 axis.ticks.y = element_line(linetype="dotted"),
                 axis.text.y.right=element_text(color="black", face="bold", hjust=1),
                 panel.background = element_blank(),
                 plot.margin = unit(c(2, 1, 0, 0), "cm"),
                 axis.line = element_line(colour = "black"),
                 plot.title = element_text(color="black", size=10, face="italic"),
                 axis.title.x = element_text(color="black", size=14),
                 #axis.title.y = element_text(color="#993333", size=14, face="bold")
)

##CS_S graph
p2 <- ggplot(CS_S, aes(fill=factor(sub,levels=c("total","main")) , y=variant_count, x=codes, color=sub)) +
  geom_bar(position="stack", stat="identity")+
  
  scale_colour_manual(values = c("black",adjustcolor( "grey95", alpha.f = 0.00001))) +
  
  scale_fill_manual(values=c('grey95', '#bae1ff')) +
  
  geom_text(aes(label = main, y = main+225), hjust = 1.3,size = 3, color="darkblue") +
  geom_text(aes(label = total_label, y = total+200), hjust = -0.3,size = 3, color="grey")+
  geom_text(aes(label = freq, y = 3200),size = 3.5, color="red")+
  
  coord_flip(xlim= NULL, ylim = (c(0, 3500)), expand=FALSE, clip = "on")+
  
  ggtitle("Genomize-SEQ Specific") +
  xlab("Evidence Code") + ylab("Variant Counts") +
  scale_x_discrete(name = "",position = "bottom",expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),breaks = c(0, 500, 1000, 1500, 2000,2500,3000,3500))


p2 <- p2 + theme(legend.box="horizontal",
                 legend.key=element_blank(),
                 legend.title=element_blank(),
                 legend.position="none",
                 panel.border = element_blank(),
                 axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.length.y = unit(0.5, "cm"),
                 axis.ticks.y = element_line(linetype="dotted"),
                 #axis.ticks.y=element_blank(),
                 panel.background = element_blank(),
                 plot.margin = unit(c(2, 0, 1, -0.2 ), "cm"),
                 #plot.title="abc",
                 axis.line = element_line(colour = "black"),
                 plot.title = element_text(color="black", size=10, face="italic", hjust=1),
                 axis.title.x = element_text(color="black", size=14),
                 #axis.title.y = element_text(color="#993333", size=14, face="bold")
)


p1_2 <- plot_grid(
  p1, NULL, p2,
  nrow = 1, ncol=3,
  rel_widths = c(1, -0.08, 1),
  #labels = c("A","",""),
  align="h") +
  draw_label(label = "Evidence\n Codes", x = 0.465, y =0.05, size = 8)+
  draw_label(label = "Specific Evidence Codes\n(Clingen vs. Genomize-SEQ)", 
             x = 0.465,
             y = 1,
             hjust = 0.5,
             vjust = 1,
             size = 15,
             lineheight = 1,
             fontface = "bold") 

p1_2 <- p1_2 + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2, linetype="dotted"),
                     plot.margin = unit(c(1, 1, 1, 1 ), "cm")
                     )



# Separate into the four datasets for the figure
CF_F_main <- data %>%
  filter(platforms=="Clingen_vs_Franklin" & orientation == "up" & sub=="main")

CF_F <- data %>%
  filter(platforms=="Clingen_vs_Franklin" & orientation == "up")

CF_C_main <- data %>%
  filter(platforms=="Clingen_vs_Franklin" & orientation == "down" & sub=="main")

CF_C <- data %>%
  filter(platforms=="Clingen_vs_Franklin" & orientation == "down")

# Order terms
## CF_C
CF_C_main_order_levels <- CF_C_main %>%
  arrange(codes) %>%
  pull(codes)

CF_C$codes <- factor(CF_C$codes, levels = (CF_C_main_order_levels))

## CF_F
CF_F_main_order_levels <- CF_C_main %>%
  arrange(codes) %>%
  pull(codes)

CF_F$codes <- factor(CF_F$codes, levels = (CF_F_main_order_levels))


##CF_C graph (reverse oriented)
p3 <- ggplot(CF_C, aes(fill=factor(sub,levels=c("total","main")) , y=variant_count, x=codes, color=sub)) +
  geom_bar(position="stack", stat="identity")+
  
  scale_colour_manual(values = c("black",adjustcolor( "grey95", alpha.f = 0.00001))) +
  
  scale_fill_manual(values=c('grey95', '#baffc9')) +
  
  geom_text(aes(label = main, y = main+225), hjust = -0.3,size = 3, color="darkgreen") +
  geom_text(aes(label = total_label, y = total+200), hjust = 1.3,size = 3, color="grey")+
  geom_text(aes(label = freq, y = 3200),size = 3.5, color="red")+
  
  coord_flip(xlim= NULL, ylim = rev(c(0, 3500)), expand=FALSE, clip = "on")+
  ggtitle("Clingen Specific") +
  ylab("Variant Counts")+
  scale_x_discrete(name = "",position = "top",expand=c(0,0)) +
  scale_y_reverse(expand=c(0, 0),breaks = c(0, 500, 1000, 1500,2000,2500,3000,3500))


p3 <- p3 + theme(legend.box="horizontal",
                 legend.key=element_blank(),
                 legend.title=element_blank(),
                 legend.position="none",
                 panel.border = element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.length.y = unit(0.5, "cm"),
                 axis.ticks.y = element_line(linetype="dotted"),
                 axis.text.y.right=element_text(color="black", face="bold", hjust=1),
                 panel.background = element_blank(),
                 plot.margin = unit(c(2, 0, 1, 0), "cm"),
                 axis.line = element_line(colour = "black"),
                 plot.title = element_text(color="black", size=10, face="italic"),
                 axis.title.x = element_text(color="black", size=14),
                 #axis.title.y = element_text(color="#993333", size=14, face="bold")
)

##CF_F graph
p4 <- ggplot(CF_F, aes(fill=factor(sub,levels=c("total","main")) , y=variant_count, x=codes, color=sub)) +
  geom_bar(position="stack", stat="identity")+
  
  scale_colour_manual(values = c("black",adjustcolor( "grey95", alpha.f = 0.00001))) +
  
  scale_fill_manual(values=c('grey95', '#ffb3ba')) +
  
  geom_text(aes(label = main, y = main+225), hjust = 1.3,size = 3, color="deeppink") +
  geom_text(aes(label = total_label, y = total+200), hjust = -0.3,size = 3, color="grey")+
  geom_text(aes(label = freq, y = 3200),size = 3.5, color="red")+
  
  coord_flip(xlim= NULL, ylim = (c(0, 3500)), expand=FALSE, clip = "on")+
  
  ggtitle("Franklin Specific") +
  ylab("Variant Counts") +
  scale_x_discrete(name = "",position = "bottom",expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),breaks = c(0, 500, 1000, 1500, 2000,2500,3000,3500))


p4 <- p4 + theme(legend.box="horizontal",
                 legend.key=element_blank(),
                 legend.title=element_blank(),
                 legend.position="none",
                 panel.border = element_blank(),
                 axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.length.y = unit(0.5, "cm"),
                 axis.ticks.y = element_line(linetype="dotted"),
                 #axis.ticks.y=element_blank(),
                 panel.background = element_blank(),
                 plot.margin = unit(c(2, 0, 1, -0.2 ), "cm"),
                 #plot.title="abc",
                 axis.line = element_line(colour = "black"),
                 plot.title = element_text(color="black", size=10, face="italic", hjust=1),
                 axis.title.x = element_text(color="black", size=14),
                 #axis.title.y = element_text(color="#993333", size=14, face="bold")
)

p3_4 <- plot_grid(
  p3, NULL, p4,
  nrow = 1, ncol=3,
  rel_widths = c(1, -0.027, 1),
  rel_heights = c(1,1),
  #labels = c("B","",""),
  align="h") +
  draw_label(label = "Evidence\n Codes", x = 0.48, y =0.05, size = 8)+
  draw_label(label = "Specific Evidence Codes\n(Clingen vs. Franklin)", 
           x = 0.48,
           y = 1,
           hjust = 0.5,
           vjust = 1,
           size = 15,
           lineheight = 1,
           fontface = "bold") 

p3_4 <- p3_4 + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2, linetype="dotted"),
                     plot.margin = unit(c(1, 1, 1, 1 ), "cm")
)






library(ggrepel)
library(tidyverse)
library(ggplot2)
library(geomtextpath)
data = read.csv("/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/5_plot_B_inputs/common_B.csv", sep="\t" , header=TRUE)
CS_data <- data %>%
  filter(platform=="CS_B")

CF_data <- data %>%
  filter(platform=="CF_B")


CS_order_levels <- CS_data %>%
  arrange(count) %>%
  pull(codes)
CS_data$codes <- factor(CS_data$codes, levels = CS_order_levels)


CF_order_levels <- CF_data %>%
  arrange(count) %>%
  pull(codes)

CF_data$codes <- factor(CF_data$codes, levels = CF_order_levels)


donut1 <- ggplot(CS_data, aes(x = codes, y = count, fill = codes, label=label)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values = c(
    "#8AD6E4",
    "#E5A8A0",
    "#B6EEE2",
    "#CAB8E5",
    "#C6EDC5",
    "#A1BDE6",
    "#EBD4AB",
    "#A0CDE2",
    "#CFB192",
    "#87C4B8",
    "#E6BBC7",
    "#B2C69A",
    "#E2C5E0"
  )) +
  geom_labelpath(aes(label=label,y=count), angle=-90, hjust=-0.2, size = 2.5,  text_only = TRUE, fill="#F6F6FF") +
  coord_polar(theta = "y", direction = 1, clip = "on")+
  ggtitle("Common Evidence Codes \n(Clingen and Genomize-SEQ)")


donut1 <- donut1 + theme(
  #panel.background = element_rect(fill = "white"),
  legend.position="none",
  panel.grid = element_line(color="black",size=0.2, linetype = "dotted"),
  #panel.grid.major = element_line(color="black",size=0.2, linetype = "dotted"),
  #panel.grid.minor = element_line(color="black",size=0.2, linetype = "dotted"),
  panel.background = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.ticks.length.y = unit(0.5, "cm"),
  #panel.border = element_rect(colour = "black", fill=NA, size=1)
  plot.title = element_text(color="black", size=15, face="bold", hjust=0.5),
)


donut2 <- ggplot(CF_data, aes(x = codes, y = count, fill = codes, label=label)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values = c(
    "#8AD6E4",
    "#E5A8A0",
    "#B6EEE2",
    "#CAB8E5",
    "#C6EDC5",
    "#A1BDE6",
    "#EBD4AB",
    "#A0CDE2",
    "#CFB192",
    "#87C4B8",
    "#E6BBC7",
    "#B2C69A",
    "#E2C5E0",
    "#A3C0A6",
    "#B5B4C8",
    "#DBE7C5",
    "#A4BEB8",
    "#F0DCD0",
    "#CDE4E9"
  )) +
  geom_labelpath(aes(label=label,y=count), angle=-90, hjust=-0.2, size = 2.5,  text_only = TRUE, fill="#F6F6FF") +
  coord_polar(theta = "y", direction = 1, clip = "on")+
  ggtitle("Common Evidence Codes \n(Clingen and Franklin)")


donut2 <- donut2 + theme(
  #panel.background = element_rect(fill = "white"),
  legend.position="none",
  panel.grid = element_line(color="black",size=0.2, linetype = "dotted"),
  #panel.grid.major = element_line(color="black",size=0.2, linetype = "dotted"),
  #panel.grid.minor = element_line(color="black",size=0.2, linetype = "dotted"),
  panel.background = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.ticks.length.y = unit(0.5, "cm"),
  #panel.border = element_rect(colour = "black", fill=NA, size=1)
  plot.title = element_text(color="black", size=14, face="bold", hjust=0.5),
)


first <- plot_grid(
  p1_2, donut1,
  ncol=2,
  nrow=1,
  rel_heights = c(1, 1),
  rel_widths = c(1, 0.5),
  labels = c("A","B"),
  align="h")


second <- plot_grid(
  p3_4, donut2,
  ncol=2,
  nrow=1,
  rel_heights = c(1, 1),
  rel_widths = c(1, 0.5),
  labels = c("C","D"),
  align="h")


#save image to file

ggsave(filename = "/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/5_plot_B_outputs/B_CS.pdf", 
       plot = first,
       units = "in", 
       width = 20,
       height= 9,
       dpi = 600)


ggsave(filename = "/Users/islekbro/Desktop/genomize_makale/2_evidence_codes/evidence_outputs/5_plot_B_outputs/B_CF.pdf", 
       plot = second,
       units = "in", 
       width = 20,
       height= 9,
       dpi = 600)
