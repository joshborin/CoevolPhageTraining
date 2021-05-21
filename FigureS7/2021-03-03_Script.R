
# This code builds on that for Figure 3
# However, instead of regressing PC1 (resistance) with # mutations
# Here, I regress PC1 with EOP, another indicator of resistance

# This will be used to create supplemental figure S7

library(reshape2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(vegan)
library(stringr)
library(pals)

setwd("CoevolReplay_FastResistance/Manuscript/FigureS7")
dat <- read.csv(file = "CompiledData_ODdiff.csv", header = T) 

##################################################
# This chunk of code is used to generate Figure S3
# But is a prerequisite for generating the data that is used for Figure 3
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
# This is just to verify that PCA is working correctly, and matches Fig S4
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
##################################################

# read in data file for EOP
mut_dat <- read.csv(file = "2021-03-03_EOPdat.csv", header = T)
# join the data and clean up (only care about PC1)
joined_dat <- merge(dater, mut_dat, by = 'Strain')



corr_plot <- ggplot(data = joined_dat,
                        aes(y = EOP, x = PC1, color = Strain, label = Strain, alpha(0.6))) +
  geom_jitter(size = 3, width = 0.001, height = 0.001) +
  labs(x = 'Principal Component 1 (PC1, 92.25%)', y = 'Efficiency of Plating (EOP)', color = 'Bacterial Isolate') +
  scale_color_manual(values = as.vector(tol(n=11))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'right',
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
        legend.key = element_blank())
print(corr_plot)
ggsave(corr_plot, file = 'PC1_x_EOP_corrplot.pdf', width = 6, height = 4)

cor.test(joined_dat$PC1, joined_dat$EOP, method = c('spearman'))
# S = 413.35
# p = 0.00037
# rho = -0.879

















