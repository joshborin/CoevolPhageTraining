

library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)

# This script is to prepare Figure S1 for the Coevolutionary Phage Training Manuscript
# I read in a data file with raw plate reader trajectories
# This file was made by manually combining the raw data of previous, independent runs of representative strains of interest
# I've made a plate map so that these runs can be matched with respective strain info
# Note that the Well does not correspond to the actual well position the strains were run in the original plate

setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/FigureS1")
rawdata <- read.csv(file = "DataFormatted.csv") 
platemap <- read.csv(file = "2020-08-28_platemap.csv")
platemap$Interaction <- paste(platemap$Treatment, platemap$Strain, sep = '-')
platemap$Label <- paste(platemap$Strain, platemap$Status, sep = '-')

# Reshape the data. Instead of rows containing the Time and readings for each Well,
# rows will contain the Time and Well ID, and the reading at that Well.
reshaped <- melt(rawdata, id=c("Time"), variable.name="Well", 
                 value.name="OD600")
annotated <- inner_join(reshaped, platemap, by="Well")
########################################################################
# Now plotting all of the strains in plots overlapping w/ and w/o phage
# First ALL TOGETHER
overlap_all <- ggplot(data= annotated,
                      mapping = aes(x=Time/3600, y=OD600,
                                    group=interaction(Well, Strain, Phage),
                                    color=Phage)) +
  geom_line(size = 1) +
  facet_wrap(interaction(Treatment, Strain) ~ ., ncol = 3) +
  scale_x_continuous(limits = c(0, 22)) +
  scale_y_continuous(limits = c(0,0.9))
print(overlap_all)
########################################################################
# THEN to separately generate plots for cI26 and 2811
dat_c <- subset(platemap, subset = Treatment == 'cI26')
length(unique(dat_c$Interaction))
unique(dat_c$Interaction)
plot_list_c = list()

for (i in 1:length(unique(dat_c$Interaction))) {
  print(unique(dat_c$Interaction)[i])
  sub_lab <- unique(dat_c$Interaction)[i]
  sub_title <- unique(dat_c$Label)[i]
  sub_dat <- subset(annotated, subset = Interaction == sub_lab)
  plot_list_c[[i]] <- ggplot(data=sub_dat,
                           mapping = aes(x=Time/3600, y=OD600,
                                         group=Well,
                                         color=Phage,
                                         linetype = Phage)) +
    geom_line(size = 0.7) +
    scale_x_continuous(limits = c(0, 20), breaks = c(0,5,10,15,20)) +
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Legend alterations
          legend.position = 'none',
          #axis.line = element_line(colour = "black"),
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
          plot.margin = unit(c(0.25,1,0,0.5), "lines")) +
    annotate('text', x = 10, y = 0.9, label = sub_title, size = 4) +
    labs(x = 'Time (Hours)', y = 'OD600') +
    scale_color_manual(breaks = c('No','Yes'), values = c('black','#00BFC4')) +
    scale_linetype_manual(breaks = c('No','Yes'), values = c(2,1))
}
plot_list_c[[1]]
################################
dat_2 <- subset(platemap, subset = Treatment == '2811')
length(unique(dat_2$Interaction))
unique(dat_2$Interaction)
plot_list_2 = list()

for (i in 1:length(unique(dat_2$Interaction))) {
  print(unique(dat_2$Interaction)[i])
  sub_lab <- unique(dat_2$Interaction)[i]
  sub_title <- unique(dat_2$Label)[i]
  sub_dat <- subset(annotated, subset = Interaction == sub_lab)
  plot_list_2[[i]] <- ggplot(data=sub_dat,
                             mapping = aes(x=Time/3600, y=OD600,
                                           group=Well,
                                           color=Phage,
                                           linetype = Phage)) +
    geom_line(size = 0.7) +
    scale_x_continuous(limits = c(0, 20), breaks = c(0,5,10,15,20)) +
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Legend alterations
          legend.position = 'none',
          #axis.line = element_line(colour = "black"),
          # Change axes title sizes
          axis.title.x=element_text(size=12, color = 'black'),
          axis.title.y=element_text(size=12, color = 'black', vjust = 0.5, hjust = 0.5),
          # Change axes ticks
          axis.line = element_line(color='black'),
          axis.ticks = element_line(color= 'black'),
          # Change axes text sizes
          axis.text.x = element_text(size=12, color = 'black'),
          axis.text.y = element_text(size=12, color = 'black'),
          panel.border = element_rect(fill = NA, colour = "black"),
          plot.margin = unit(c(0.25,1,0,0.5), "lines")) +
    annotate('text', x = 10, y = 0.9, label = sub_title, size = 4) +
    labs(x = 'Time (Hours)', y = expression(paste('OD600'))) +
    scale_color_manual(breaks = c('No','Yes'), values = c('black','#F8766D')) +
    scale_linetype_manual(breaks = c('No','Yes'), values = c(2,1))
}
plot_list_2[[2]]
plot_matrix <- ggarrange(plot_list_c[[1]], plot_list_2[[1]], plot_list_2[[5]],
                         plot_list_c[[4]], plot_list_2[[2]], plot_list_2[[6]],
                         plot_list_c[[3]], plot_list_2[[3]], plot_list_2[[7]],
                         plot_list_c[[2]], plot_list_2[[4]], plot_list_2[[8]],
                         ncol = 3, nrow = 4,
                         widths = 3, heights = 3)
plot_matrix
ggsave(filename = 'FigureS1.pdf', plot_matrix,
       width = 8, height = 8)






