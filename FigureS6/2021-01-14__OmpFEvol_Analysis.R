
# 2021-01-14
# This script is appropriated from that used to make Figure 1
#  However, instead of overall population dynamics
#  I will plot the titer of phage in untrained phage treatment
#  that have evolved to use the new, OmpF receptor.

library(ggplot2)
library(ggpubr)
setwd("CoevolReplay_FastResistance/Manuscript/FigureS6")
raw.data <- read.csv(file = "2020-04-02_CoevolutionCounts.csv")

colnames(raw.data)
# Will just focus on Hi inoculum data from unstrained phage lysates
dater <- subset(raw.data, Inoculum == 'Hi'
                  & Phage_Strain == 'cI26'
                  & Plated == 'Lysate')
  
################################################################################
  # Plot for lysate phage titers
  lysate.plot <- ggplot(data = subset(dater, Plated == 'Lysate'),
                           aes(x = Timepoint,
                               y = log10(as.numeric(TiterOn19A)),
                               color = as.factor(Rep))) +
    geom_path(aes(linetype = as.factor(Rep)), size = 0.5) +
    geom_point(size = 1.5) +
    annotate('rect', xmin = -0.3, xmax = 30.3, ymin = 0, ymax = 2.7,
             color = NA, fill = 'grey70', alpha = 0.5) +
    annotate('text', label = 'Limit of Detection',
             x = 15, y = 1.35, size = 5, color = 'grey20') +
    labs(y = bquote('OmpF+ Phage Titer ('~log[10]~'PFU/mL)'), x = 'Day') +
    scale_y_continuous(breaks = c(1,2,3,4,5),
                       limits = c(0,6),
                       expand = c(0,0)) +
    scale_x_continuous(breaks = c(0,17,19,22,30),
                       expand = c(0.01,0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Legend alterations
          legend.position = 'none',
          #axis.line = element_line(colour = "black"),
          # Change axes title sizes
          axis.title.x=element_text(size = 12, color = 'black'),
          axis.title.y=element_text(size = 12, color = 'black'),
          # Change axes ticks
          axis.line = element_line(color='black'),
          axis.ticks = element_line(color= 'black'),
          # Change axes text sizes
          axis.text.x = element_text(size=12, color = 'black'),
          axis.text.y = element_text(size=12, color = 'black'),
          panel.border = element_rect(fill = NA, colour = "black")) +
    geom_hline(yintercept = 2.7, linetype = 'dashed') +
    guides(linetype = FALSE)
  lysate.plot
  
  ggsave(filename = "OmpF Plot.pdf", 
         plot = lysate.plot, width = 6, height = 4, bg = 'transparent')

