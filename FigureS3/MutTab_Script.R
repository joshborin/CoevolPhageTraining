
# This script uses consolidated mutation data in MutTab
#  to generate the matrix for figure S3
#  The data that was used to create MutTab.csv is available in the 'Sequencing' folder

# Please note that some formatting of the figure was done in post using Inkscape
#  The matrix figure generated here was spliced together with the table (also located in the repo)
#  Rows and strains were colored in post
#  Lastly, strain JB47 was added in post

library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(viridis)
library(ggpubr)

# Import the data file and map which includes information about Environment
setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/FigureS3")
rawdata <- read.csv(file = "MutTab.csv", check.names = FALSE) 
# Melt data and map together into long format. Then make sure all
# data types are factors, clean up names, and order them for plotting.
reshaped <- melt(rawdata, id=c("Position"), variable.name="Isolate", 
                 value.name="mut")
reshaped$mut <- as.factor(as.numeric(reshaped$mut))
reshaped$Isolate <- as.factor(reshaped$Isolate)
reshaped$Isolate
levels(as.factor(reshaped$Position))
order(reshaped$Position)

MutTable <- ggplot(reshaped, aes(x=Isolate, y=as.factor(Position))) +
  geom_tile(aes(fill = mut, color = 'poop'),
            show.legend = FALSE) +
  scale_fill_viridis(option = 'cividis', discrete=TRUE, begin = 0.5, end = 0.1) +
  scale_color_manual(values = 'white') +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev(levels(as.factor(reshaped$Position)))) +
  guides(fill = FALSE) + ylab("") + xlab("") +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.top = element_text(color = 'black', angle=90, vjust = 0.5, hjust = -0.1),
        axis.text.y = element_text(color = 'black'))
print(MutTable)
ggsave(MutTable, filename = 'MutTable.pdf', height = 7, width = 4)
  



