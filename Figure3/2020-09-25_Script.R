
# This script reads in two metadata files
##  One contains the number of total and putative resistance mutations in bacterial isolates
##  The other contains computed differences between OD for representative bacterial isolates
##   growing without and with phage. The OD trajectories will be made available upon request.
# Finally, the script will generate each panel of Figure 3 independently, as well as assembled.

library(reshape2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(vegan)
library(stringr)
library(jcolors)
library(ggpubr)

setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/Figure3")
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
ggsave(pca.plot1, file = 'PCA_plot.pdf', width = 8, height = 4)
##################################################

# read in data file for numbers of mutations in each strain
mut_dat <- read.csv(file = "2020-07-20_MutMetadata.csv", header = T)
# join the data and clean up (only care about PC1)
joined_dat <- merge(dater, mut_dat, by = 'Strain')
# now want to separate 2811 and cI26 data

############################################################
# LINEAR REGRESSIONS WITH PUTATIVE MUTATIONS
############################################################

#cI26
lm_put_cI26 <- lm(PC1 ~ NumPutativeMut, data = subset(joined_dat, Treatment == 'cI26' | Treatment == 'Ancestor'))
summary(lm_put_cI26)
# R-sq = 0.996
# positive slope SIGNIFICANT (p=0.002)
put_cI26_plot <- ggplot(data = subset(joined_dat, Treatment == 'cI26' | Treatment == 'Ancestor'),
                        aes(y = PC1, x = NumPutativeMut, color = Strain, label = Strain)) +
  geom_jitter(size = 3, width = 0.01) +

  geom_abline(slope = lm_put_cI26$coef[2],
              intercept = lm_put_cI26$coef[1]) +
  scale_x_continuous(breaks = c(0,1), labels = c('0','1')) +
  labs(x = 'Putative Resistance Mutations', y = 'Phage Resistance (PC1, 92.25%)', color = 'Bacterial Isolate') +
  scale_color_jcolors(palette = 'pal2') +
  annotate('text', x = min(subset(joined_dat, Treatment == 'cI26' | Treatment == 'Ancestor')$NumPutativeMut),
           y = max(subset(joined_dat, Treatment == 'cI26' | Treatment == 'Ancestor')$PC1), label = 'A', size = 6) +
  annotate('text', x = 0.6, y = -0.5,
           label = expression(paste('R'^2, ' = 0.9961')), size = 4) +
  annotate('text', x = 0.578, y = -0.7,
           label = expression(paste('p = 0.002')), size = 4) +  
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
        legend.key = element_blank())
print(put_cI26_plot)
ggsave(put_cI26_plot, file = 'Figure3A.pdf', width = 6, height = 4)

# 2811 without Ftest strains
lm_put_2811 <- lm(PC1 ~ NumPutativeMut, data = subset(joined_dat, Treatment != 'cI26' & Treatment != 'Ftest'))
summary(lm_put_2811)
# R-sq = 0.746
# positive slope SIGNIFICANT (p=0.003)
put_2811_plot <- ggplot(data = subset(joined_dat, Treatment != 'cI26' & Treatment != 'Ftest'),
                        aes(y = PC1, x = NumPutativeMut, color = Strain, label = Strain)) +
  geom_point(size = 3) +
  geom_abline(slope = lm_put_2811$coef[2],
              intercept = lm_put_2811$coef[1]) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  labs(x = 'Putative Resistance Mutations', y = 'Phage Resistance (PC1, 92.25%)', color = 'Bacterial Isolate') +
  scale_color_jcolors(palette = 'rainbow') +
  annotate('text', x = min(subset(joined_dat, Treatment == 'cI26' | Treatment == 'Ancestor')$NumPutativeMut),
           y = max(subset(joined_dat, Treatment == 'cI26' | Treatment == 'Ancestor')$PC1), label = 'B', size = 6) +
  annotate('text', x = 1.6, y = 0.5,
           label = expression(paste('R'^2, ' = 0.7462')), size = 4) +
  annotate('text', x = 1.47, y = 0.3,
           label = expression(paste('p = 0.003')), size = 4) +  
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
        legend.key = element_blank())
print(put_2811_plot)
ggsave(put_2811_plot, file = 'Figure3B.pdf', width = 6, height = 4)

plot_matrix <- ggarrange(put_cI26_plot, put_2811_plot,
                         legend = 'right', common.legend = FALSE,
                         ncol = 1, nrow = 2,
                         align = 'v')

plot_matrix
ggsave(filename = 'Figure3AB.pdf', plot_matrix,
       width = 6, height = 6)

