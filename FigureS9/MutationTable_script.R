
# 2021-03-07
# Script adopted from that for NS Binding manuscript

library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(viridis)

# Import the data file and map which includes information about Environment
setwd("CoevolReplay_FastResistance/Manuscript/FigureS9")
rawdata <- read.csv(file = "MutationTable.csv") 
colnames(rawdata) <- c('Mutation','9+','11+','14+','21+','23+','24+',
                    '1-','3-','6-','10-','13-','15-')
map <- read.csv(file = "map.csv")

# Melt data and map together into long format. Then make sure all
# data types are factors, clean up names, and order them for plotting.
reshaped <- melt(rawdata, id=c("Mutation"), variable.name="Phage", 
                 value.name="mut")
reshaped$mut <- as.factor(as.numeric(reshaped$mut))
reshaped$Phage <- factor(reshaped$Phage,
                           levels = c('15-','13-','10-','6-','3-','1-',
                                      '24+','23+','21+','14+','11+','9+'))
reshaped$Phage

reshaped$Mutation <- factor(reshaped$Mutation,
                           levels = c('C2879T','C2938A','C2969T','C2988A','A2989G',
                                      'T2991G','C2999A','C2999T','C3033A','A3034G',
                                      'C3119T','C3147A','C3227A','C3227T','G3229A',
                                      'T3230C','A3233G','T3248C','C3310T','G3319A',
                                      'A3320G','T3321A'))
reshaped$Mutation

ggplot(reshaped, aes(x=Mutation, y=Phage)) +
  geom_tile(aes(fill = mut, color = 'poop'), show.legend = FALSE) +
  scale_fill_viridis(option = 'cividis', discrete=TRUE, begin = 0.5, end = 0.1) +
  scale_color_manual(values = 'white') +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(labels=c('OmpF- 6','OmpF- 5','OmpF- 4','OmpF- 3',
                            'OmpF- 2','OmpF- 1','OmpF+ 6','OmpF+ 5',
                            'OmpF+ 4','OmpF+ 3','OmpF+ 2','OmpF+ 1')) +
  guides(fill = FALSE) + ylab("") + xlab("") +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.top = element_text(color = 'black', angle=45, vjust = 0, hjust = 0.1),
        axis.text.y = element_text(color = 'black'))

ggsave('FigureS9_MutTable.pdf', height = 4, width = 6)


val <- c(1,2,3,4,5,6)
t <- c(1,1,1,1,1,1)
fac <- c('1','2','3','4','5','6')
symdat <- data.frame(cbind(val,t,fac))

symplot <- ggplot(symdat, aes(x=t,y=val,shape=fac)) +
  geom_point() +
  theme(legend.key = element_rect(fill='NA',color=NA))
ggsave('symbols.pdf', height = 4, width = 6)
