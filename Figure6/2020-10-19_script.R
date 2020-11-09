
# This script uses the computed growth rates of phage pre- and post-recombination
#  to generate Figure 6
# For this script, I read in a cleaned up data set but the original data is
#  also in the same Figure 6 folder in the repository

library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)

setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/Figure6")
dat <- read.csv(file = "DataForFig.csv") 

levels(dat$Host)
dat <- dat %>%
  arrange(PercentGrowth) %>%
  # Note: I am arranging these strains by treatment
  # and then PC1 value as determined in Fig. S3
  mutate(Host = factor(Host, levels=c('wt',
                                      'malT-',
                                      'malT-ompF-')))
levels(dat$Host)
bplot2 <- ggplot(dat, aes(x = Host, y = RelativeGrowth)) +
  geom_boxplot() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(y = 'Relative Growth Rate of Recombinant', x = 'Host Genotype') +
  scale_y_continuous(limits = c(0.9,2.8), breaks = c(1, 1.5, 2, 2.5)) +
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
        axis.text.x = element_text(size=12, color = 'black', face = 'italic'),
        axis.text.y = element_text(size=12, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x = 1, y = 2.75, label = 'A', size = 6) +
  annotate('text', x = 2, y = 2.75, label = 'A', size = 6) +
  annotate('text', x = 3, y = 2.75, label = 'B', size = 6)
print(bplot2)
ggsave(filename = 'Figure6.pdf', bplot2,
       width = 4, height = 4)









