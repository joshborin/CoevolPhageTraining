
# Figure 2
# 
# This script uses a master metadata set
# To plot the resistance status of different populations across various time-points
# This script will generate individual plots as well as assemble them into Figure 2
#  into stacked bar charts
#
# The metadata set was created from looking at growth rate trajectories for all isolates
# These trajectories are included in the folder
# Raw data for these trajectories will be made available upon request

library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
# Set working directory to point to github repo folder and script
# Example:
setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/Figure2")
master <- read.csv(file = "MasterData_forR.csv") 

# Clean up the dataset and remove missing data
dat_nb <- subset(master, LiquidStatus == 'R' | LiquidStatus == 'S' | LiquidStatus == 'P')
dat_nb$Timepoint <- as.factor(dat_nb$Timepoint)

TEST <- dat_nb %>% 
  group_by(Population,Timepoint,LiquidStatus) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))
TEST$Timepoint <- as.numeric(as.character(TEST$Timepoint))
length(levels(TEST$Population))

######## cI26 POPULAION 1
cp1 <- subset(TEST, subset = Population == 'cI26 Hi 1')
cp1_plot <- ggplot(cp1, aes(x = as.numeric(Timepoint), y = perc,
                fill = factor(LiquidStatus, levels=c('R','P','S'))))+
  geom_bar(stat='identity', width = 0.7) +
  labs(y = 'Proportion of Population', fill = "Resistance") +
  scale_x_continuous(limits = c(0,31),
                     breaks = c(0,3,5,10,15,20,25,30),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.03)) +
  scale_fill_manual(values = c(
    # complete
    '#FF0000',
    # partial
    '#ffaa00',
    # sensitive
    '#00aaff'),
    labels = c('Completely Resistant',
               'Partially Resistant',
               'Sensitive')) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            # Legend alterations
            legend.position = 'top',
            legend.title = element_blank(),
            legend.text = element_text(margin = margin(r = 0.5, unit = 'inch'),
                                       size = 12),
        # Change axes title sizes
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12, color = 'black'),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0,0.5,0.25,1), "lines")) +
  annotate('text', x = 1.4, y = .95, label = 'A', size = 6)
print(cp1_plot)
ggsave(filename = 'Figure2A.pdf', cp1_plot,
       width = 4.22, height = 2.5)

########################################################################
########################################################################
######## cI26 POPULAION 2
cp2 <- subset(TEST, subset = Population == 'cI26 Hi 2')
cp2_plot <- ggplot(cp2, aes(x = as.numeric(Timepoint), y = perc,
                            fill = factor(LiquidStatus, levels=c('R','P','S'))))+
  geom_bar(stat='identity', width = 0.7) +
  labs(y = 'Proportion of Population', fill = "Resistance") +
  scale_x_continuous(limits = c(0,31),
                     breaks = c(0,3,5,10,15,20,25,30),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.03)) +
  scale_fill_manual(values = c(
    # complete
    '#FF0000',
    # partial
    '#ffaa00',
    # sensitive
    '#00aaff')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'none',
        # Change axes title sizes
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12, color = 'black'),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0.25,0.5,0.25,1), "lines")) +
  annotate('text', x = 1.4, y = .95, label = 'B', size = 6)
print(cp2_plot)
ggsave(filename = 'Figure2B.pdf', cp2_plot,
       width = 4.22, height = 2.5)

########################################################################
########################################################################
######## cI26 POPULAION 3
cp3 <- subset(TEST, subset = Population == 'cI26 Hi 3')
cp3_plot <- ggplot(cp3, aes(x = as.numeric(Timepoint), y = perc,
                            fill = factor(LiquidStatus, levels=c('R','P','S'))))+
  geom_bar(stat='identity', width = 0.7) +
  labs(x = 'Day', y = 'Proportion of Population', fill = "Resistance") +
  scale_x_continuous(limits = c(0,31),
                     breaks = c(0,3,5,10,15,20,25,30),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.03)) +
  scale_fill_manual(values = c(
    # complete
    '#FF0000',
    # partial
    '#ffaa00',
    # sensitive
    '#00aaff')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'none',
        # Change axes title sizes
        axis.title.x=element_text(size=12, color = 'black'),
        axis.title.y=element_text(size=12, color = 'black'),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0.25,0.5,0,1), "lines")) +
  annotate('text', x = 1.4, y = .95, label = 'C', size = 6)
print(cp3_plot)
ggsave(filename = 'Figure2C.pdf', cp3_plot,
       width = 4.22, height = 2.6)

########################################################################
########################################################################
######## 2811 POPULAION 1
twop1 <- subset(TEST, subset = Population == '2811 Hi 1')
twop1_plot <- ggplot(twop1, aes(x = as.numeric(Timepoint), y = perc,
                            fill = factor(LiquidStatus, levels=c('R','P','S'))))+
  geom_bar(stat='identity', width = 0.7) +
  labs(y = 'Proportion of Population', fill = "Resistance") +
  scale_x_continuous(limits = c(0,31),
                     breaks = c(0,3,5,10,15,20,25,30),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.03)) +
  scale_fill_manual(values = c(
    # complete
    #'#FF0000',
    # partial
    '#ffaa00',
    # sensitive
    '#00aaff')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'none',
        # Change axes title sizes
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0,1,0.25,0.5), "lines")) +
  annotate('text', x = 1.4, y = .95, label = 'D', size = 6)
print(twop1_plot)
ggsave(filename = 'Figure2D.pdf', twop1_plot,
       width = 4, height = 2.5)

########################################################################
########################################################################
######## 2811 POPULAION 2
twop2 <- subset(TEST, subset = Population == '2811 Hi 2')
twop2_plot <- ggplot(twop2, aes(x = as.numeric(Timepoint), y = perc,
                                fill = factor(LiquidStatus, levels=c('R','P','S'))))+
  geom_bar(stat='identity', width = 0.7) +
  labs(y = 'Proportion of Population', fill = "Resistance") +
  scale_x_continuous(limits = c(0,31),
                     breaks = c(0,3,5,10,15,20,25,30),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.03)) +
  scale_fill_manual(values = c(
    # complete
    '#FF0000',
    # partial
    '#ffaa00',
    # sensitive
    '#00aaff'),
    labels = c('Complete Resistance','Partial Resistance','Sensitive')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'none',
        # Change axes title sizes
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0.25,1,0.25,0.5), "lines")) +
  annotate('text', x = 1.4, y = .95, label = 'E', size = 6)
print(twop2_plot)
ggsave(filename = 'Figure2E.pdf', twop2_plot,
       width = 4, height = 2.5)

########################################################################
########################################################################
######## 2811 POPULAION 3
twop3 <- subset(TEST, subset = Population == '2811 Hi 3')
twop3_plot <- ggplot(twop3, aes(x = as.numeric(Timepoint), y = perc,
                                fill = factor(LiquidStatus, levels=c('R','P','S'))))+
  geom_bar(stat='identity', width = 0.7) +
  labs(x = 'Day', y = 'Proportion of Population', fill = "Resistance") +
  scale_x_continuous(limits = c(0,31),
                     breaks = c(0,3,5,10,15,20,25,30),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.03)) +
  scale_fill_manual(values = c(
    # complete
    '#FF0000',
    # partial
    '#ffaa00',
    # sensitive
    '#00aaff')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'none',
        # Change axes title sizes
        axis.title.x=element_text(size=12, color = 'black'),
        axis.title.y=element_blank(),
        # Change axes ticks
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color= 'black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.margin = unit(c(0.25,1,0,0.5), "lines")) +
  annotate('text', x = 1.4, y = .95, label = 'F', size = 6)
print(twop3_plot)
ggsave(filename = 'Figure2F.pdf', twop3_plot,
       width = 4, height = 2.6)

plot_matrix <- ggarrange(cp1_plot, twop1_plot, cp2_plot, twop2_plot, cp3_plot, twop3_plot,
          legend = 'top', common.legend = TRUE,
          ncol = 2, nrow = 3,
          widths = 3, heights = 2.5)
plot_matrix
ggsave(filename = 'Figure2A-F.pdf', plot_matrix,
       width = 8.5, height = 8)












