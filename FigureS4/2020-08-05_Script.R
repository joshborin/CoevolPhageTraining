
# This script is similar to that for Figure 3
# It reads in data of OD differences between cultures without and with phages
# Then it conducts a PCA
# And generates Figure S4 of this PCA
# In post (Inkscape), I moved the strain labels on the plot to be adjacent to points

library(reshape2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(vegan)
library(stringr)
library(jcolors)

setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/FigureS4")
dat <- read.csv(file = "CompiledData_ODdiff.csv", header = T) 

pca_test <- rda(dat[,4:244])
biplot(pca_test)
summary(pca_test)
# PC1 explains 92.25% of the variation
# PC2 explains 5.90% of the variation
# save the PCA summary as well as the 'scores'
a <- summary(pca_test)
site.scores <- a$sites
# save these scores (=Principal components) and corresponding strains to data frame
pca <- as.data.frame(site.scores[,c(1,2,3,4,5,6)])
Strains <- as.data.frame(dat$Strain)
Labels <- as.data.frame(dat$Label)
Treatment <- as.data.frame(dat$Treatment)
colnames(Strains) <- c('Strain')
colnames(Labels) <- c('Label')
colnames(Treatment) <- c('Treatment')
dater <- bind_cols(Strains, Labels, Treatment, pca)

## Plot
pca.plot1 <-  ggplot(data = subset(dater, subset = Treatment != 'Ftest'),
                     aes(x = PC1, y = PC2, color = Strain, label = Strain)) +
  geom_point(size = 4) +
  geom_text(show.legend = FALSE) +
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
        plot.margin = unit(c(0.25,1,0,0.5), "lines"),
        legend.key=element_blank()) +
  labs(x = 'PC1 (92.25%)', y = 'PC2 (5.90%)') +
  scale_color_jcolors(palette = 'pal8')
pca.plot1
ggsave(pca.plot1, file = 'FigureS4.pdf', width = 6, height = 4)
