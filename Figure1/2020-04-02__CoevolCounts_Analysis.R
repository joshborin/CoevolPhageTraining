
# 2020-04-02
# This script will analyze the CFU and PFU counts from plating populations
#  during the Coevolution Experiment
# Then it will generate each panel of Figure 1 independently, as well as assembled.

library(ggplot2)
library(ggpubr)
# Set working directory to point to github repo folder and script
# Example:
setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/Figure1")
raw.data <- read.csv(file = "2020-04-02_CoevolutionCounts.csv")

colnames(raw.data)
Hi.data <- subset(raw.data, Inoculum == 'Hi')

# Plot coculture bacteria titers and remove NAs
coculture.plot <- ggplot(data = subset(Hi.data, Plated == 'Coculture' & Titer !='NA'),
                      aes(x = Timepoint,
                          y = log10(as.numeric(Titer)),
                          color = interaction(Phage_Strain, Inoculum))) +
  geom_path(aes(linetype = as.factor(Rep)), size = 0.5) +
  geom_point(size = 1.5) +
  annotate('rect', xmin = -0.3, xmax = 30.3, ymin = 0, ymax = 2,
            color = NA, fill = 'grey70', alpha = 0.5) +
  annotate('text', label = 'Limit of Detection',
           x = 15, y = 1, size = 5, color = 'grey20') +
  annotate('text', label = "A",
           x = 0.5, y = 10.3, size = 6, color = 'black') +
  labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),
                     limits = c(0,11),
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_linetype_discrete(name='Population') +
  scale_color_discrete(name='Evolved with:',
                       labels = c('Pre-Coevolved Phage', 'Ancestral Phage')) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'none',
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
  # hline is the limit of detection
  geom_hline(yintercept = 2, linetype = 'dashed') +
  # Removes the Population number from the legends
  guides(linetype = FALSE)
coculture.plot
  ggsave(filename = "Fig1A.pdf", 
       plot = coculture.plot, width = 6, height = 4, bg = 'transparent')
  
################################################################################
lysate.dat <- subset(Hi.data, Plated=='Lysate')
  
  # Plot for lysate phage titers
  lysate.plot <- ggplot(data = subset(Hi.data, Plated == 'Lysate'),
                           aes(x = Timepoint,
                               y = log10(as.numeric(Titer)),
                               color = interaction(Phage_Strain, Inoculum))) +
    geom_path(aes(linetype = as.factor(Rep)), size = 0.5) +
    geom_point(size = 1.5) +
    annotate('rect', xmin = -0.3, xmax = 30.3, ymin = 0, ymax = 2.7,
             color = NA, fill = 'grey70', alpha = 0.5) +
    annotate('text', label = 'Limit of Detection',
             x = 15, y = 1.35, size = 5, color = 'grey20') +
    annotate('text', label = "B",
             x = 0.5, y = 11.3, size = 6, color = 'black') +
    labs(y = bquote('Phage Titer ('~log[10]~'PFU/mL)'), x = 'Day') +
    scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11),
                       limits = c(0,12),
                       expand = c(0,0)) +
    scale_x_continuous(expand = c(0.01,0)) +
    scale_linetype_discrete(name='Population') +
    scale_color_discrete(name='Evolved with:',
                         labels = c('Pre-Coevolved Phage', 'Ancestral Phage')) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Legend alterations
          legend.position = 'none',
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
    # hline is the limit of detection
    geom_hline(yintercept = 2.7, linetype = 'dashed')
    # Removes the Population number from the legends
  lysate.plot
  ggsave(filename = "Fig1B.pdf", 
         plot = lysate.plot, width = 6, height = 4, bg = 'transparent')
  plot_matrix <- ggarrange(coculture.plot, lysate.plot,
                           ncol = 1, nrow = 2)
  plot_matrix
  ggsave(filename = 'Fig1AB.pdf', plot_matrix,
         width = 6, height = 8)
  
