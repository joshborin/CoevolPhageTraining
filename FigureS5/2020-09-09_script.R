
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(vegan)
library(stringr)

# This script is going to generate a 3 panel Figure S5 from data on cI26 and cI26+ phages
setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/FigureS5")
rawdata <- read.csv(file = "DataFormatted.csv") 
platemap <- read.csv(file = "2020-07-29_platemap.csv")
# Note thats this has scrap data from other genotypes that I'm not considering from this data
# e.g. future phages and 2811 phages which were then conducted in separate experiments

reshaped <- melt(rawdata, id=c("Time"), variable.name="Well", 
                 value.name="OD600")
annotated <- inner_join(reshaped, platemap, by="Well")
write.csv(annotated, "2020-09-09_annotated.csv")

# First, check Blank Controls
ggplot(data= subset(annotated, Bacteria == "Blank"),
       mapping = aes(x=Time/3600, y=OD600,
                     group=Well,
                     color = Well)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 22)) +
  scale_y_continuous(limits = c(0, 0.7))

# Lets get rid of that scrap data and focus those to be analyzed here...
annotated <- subset(annotated, subset = Treatment == 'cI26' &
                      Bacteria != 'REL606' &
                      Status != 'Future')
# All wells with phage
plot_cI26 <- ggplot(data=annotated,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=Well,
                                  color=interaction(Phage, Bacteria))) +
  geom_line(size = 1) +
  facet_wrap(Bacteria ~ ., ncol = 4) +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_cI26)
is.numeric(annotated$Time)
# NOW I WANT TO SPLIT UP THE DATA INTO DIFFERENT POPULATIONS AND I WILL MAKE THE PLOTS SEPARATELY
# And add in Blank controls
P1 <- subset(annotated, Bacteria == 'P1-T3-2-R')
#####################################################################################
#### Plots for Figures ####
A <- ggplot(data=P1,
                 mapping = aes(x=Time/3600, y=OD600,
                               group=Well,
                               color=Phage)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = c(0.2,0.6),
        legend.key=element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        #axis.line = element_line(colour = "black"),
        # Change axes title sizes
        axis.title.x=element_text(size=12, color = 'black'),
        axis.title.y=element_text(size=12, color = 'black'),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = 'Time (Hours)', y = 'Bacterial Density (OD600)') +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="A", size = 6) +
  scale_color_manual(values = c('#00BFC4', '#440154FF'))
print(A)
ggsave(plot=A,
       filename="FigureS5A.pdf",
       width = 6, height = 4)

####################################################################################
#############################   POPULATION 2 #######################################
####################################################################################
P2 <- subset(annotated, Bacteria == 'P2-T3-1-R')
B <- ggplot(data=P2,
            mapping = aes(x=Time/3600, y=OD600,
                          group=Well,
                          color=Phage)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = c(0.2,0.6),
        legend.key=element_blank(),
        #axis.line = element_line(colour = "black"),
        # Change axes title sizes
        axis.title.x=element_text(size=12, color = 'black'),
        axis.title.y=element_text(size=12, color = 'black'),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = 'Time (Hours)', y = 'Bacterial Density (OD600)') +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="B", size = 6) +
  scale_color_manual(values = c('#00BFC4', '#440154FF'))
print(B)
ggsave(plot=B,
       filename="FigureS5B.pdf",
       width = 6, height = 4)


####################################################################################
#############################   POPULATION 3 #######################################
####################################################################################
P3 <- subset(annotated, Bacteria == 'P3-T3-1-R')
C <- ggplot(data=P3,
            mapping = aes(x=Time/3600, y=OD600,
                          group=Well,
                          color=Phage)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = c(0.2,0.6),
        legend.key=element_blank(),
        #axis.line = element_line(colour = "black"),
        # Change axes title sizes
        axis.title.x=element_text(size=12, color = 'black'),
        axis.title.y=element_text(size=12, color = 'black'),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=10, color = 'black'),
        axis.text.y = element_text(size=10, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = 'Time (Hours)', y = 'Bacterial Density (OD600)') +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="C", size = 6) +
  scale_color_manual(values = c('#00BFC4', '#440154FF'))
print(C)
ggsave(plot=C,
       filename="FigureS5C.pdf",
       width = 6, height = 4)

###############################################################################
###############################################################################
###############################################################################
plot_matrix <- ggarrange(A, B, C,
                         common.legend = FALSE,
                         ncol = 3, nrow = 1,
                         widths = 6, heights = 4)
plot_matrix
ggsave(filename = 'FigureS5A-C.pdf', plot_matrix,
       width = 12, height = 4)










