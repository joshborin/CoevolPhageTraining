
# This script will analyze fitness data from competitions
#  and generate Figure 4
#  it also includes the statistical tests that were conducted

library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/Figure4")
dat <- read.csv(file = "2020-07-23_Results_forR.csv") 

levels(dat$Label)
dat <- dat %>%
  arrange(W) %>%
  # Note: I am arranging these strains by treatment
  # and then PC1 value as determined in Fig. S3
  mutate(Label = factor(Label, levels=c('P2-T3-1-R',
                                      'P3-T3-1-R',
                                      'P1-T3-2-R',
                                      'P1-T30-1-S',
                                      'P3-T5-1-P',
                                      'P1-T25-1-P',
                                      'P2-T15-1-P',
                                      'P2-T20-1-P',
                                      'P2-T25-10-P',
                                      'P2-T25-9-R',
                                      'P3-T15-2-R',
                                      'JB39-P',
                                      'JB42-P',
                                      'JB43-P',
                                      'JB47-R',
                                      'JB48-P')))
  all_plot <- ggplot(dat, aes(x=Label, y=W, color = Treatment)) +
    geom_boxplot() +
    labs(x= 'Strain',y= 'Relative Fitness') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.4))
  all_plot

sub_dat <- subset(dat, Treatment != 'Ftest')
sub_plot <- ggplot(sub_dat, aes(x=Label, y=W, color = Treatment)) +
  geom_boxplot() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x= 'Bacterial Strain',y= 'Relative Fitness (W)') +
  scale_y_continuous(limits = c(0,1.25),
                     expand = c(0,0)) +
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
        axis.text.x = element_text(size=10, color = 'black',
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10, color = 'black'),
        panel.border = element_rect(fill = NA, colour = "black"))
sub_plot
ggsave(sub_plot, filename = 'Figure4.pdf', width = 6, height = 4)

t.test(W~Treatment, data = sub_dat)
# Two sample t-test comparing mean fitness of cI26 to 2811 resistance
# p= 1.998 * 10^-8
# t = -8.1924
# Mean W in cI26 = 0.972425
# Mean W in 2811 = 0.595780

t.test(x = subset(sub_dat, Treatment == 'cI26')$W, mu = 1)
# One sample t-test comparing mean fitness of cI26 to 1
# Significantly different
# p = 0.005
# thats annoying

# Lastly, a test to compare fitness of partially and completely resistant strains against 2811
R <- subset(sub_dat, Treatment == '2811' & Isolate == '9' | Strain == '2811 P3 T15 2')
R$Resistance <- 'R'
P <- subset(sub_dat, Treatment == '2811'
            & Isolate != '9'
            & Strain != '2811 P3 T15 2'
            & Strain != '2811 P1 T30 1')
P$Resistance <- 'P'
test <- rbind(R,P)
t.test(W~Resistance, data = test)
# Two sample t-test comparing fitness cost of partial and complete resistance in 2811
# NOT different
# p = 0.67





