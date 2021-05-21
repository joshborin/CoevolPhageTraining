
# 2021-02-17
# This script is adapted from the code used for Figure 1A
#  of the Coevol Phage Training manuscript

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidystats)
setwd('CoevolReplay_FastResistance/Manuscript/Figure7')
raw.data <- read.csv(file = "2021-02-14_CFU_titers.csv")

colnames(raw.data)
as.factor(raw.data$Category)

raw.data$Category <- factor(raw.data$Category, levels = c('2811','OmpF+','cI26','OmpF-'))
raw.data$Category  # notice the changed order of factor levels

# Plot bacteria titer, ALL TREATMENTS
# NOTE THAT I HAVE REMOVED SOME NAs WHICH FILL IN THE GAP OF MISSING VALUES IN THE PLOT
coculture.plot <- ggplot(data = subset(raw.data, Titer != 'NA'),
                      aes(x = Day,
                          y = log10(as.numeric(Titer)),
                          group = interaction(Phage_Strain, Rep),
                          color = Category)) +
  geom_path(aes(linetype = as.factor(Rep)), size = 1) +
  geom_point(size = 1.5) +
  labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10),
                     limits = c(0,10),
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0.01,0), breaks = c(0,1,2,3,4,5,6,7)) +
  scale_linetype_discrete(name='Population') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'right',
        #axis.line = element_line(colour = "black"),
        # Change axes ticks
        axis.line = element_line(color='black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black')) +
  # hline is the limit of detection
  geom_hline(yintercept = 2, linetype = 'dashed') +
  # Removes the Population number from the legends
  guides(linetype = FALSE)
coculture.plot

ggsave(filename = "2021-02-17_bacteriaplot.png", 
       plot = coculture.plot, width = 6, height = 4, bg = 'transparent')

########################################################################
# cI26 vs 2811
###############
dat1 <- subset(raw.data, 
               Category == 'cI26' |
                 Category == '2811' &
                 !is.na(Titer))
dat1$Category <- factor(dat1$Category, levels=c("cI26", "2811"))
plot_cv2 <- ggplot(dat1,
                         aes(x = Day,
                             y = log10(as.numeric(Titer)),
                             group = interaction(Category, Rep),
                             color = Category)) +
  geom_path(size = 1) +
  geom_point(aes(shape = factor(Rep)), size = 2) +
  labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),
                     limits = c(0,10),
                     expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05), breaks = c(0,1,2,3,4,5,6,7)) +
  scale_linetype_discrete(name='Population') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = 'white'),
        legend.key.width = unit(1, 'cm'),
        #axis.line = element_line(colour = "black"),
        # Change axes ticks
        axis.line = element_line(color='black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black')) +
  # hline is the limit of detection
  geom_hline(yintercept = 2, linetype = 'dashed') +
  # Removes the Population number from the legends
  guides(shape = FALSE) +
  annotate(geom = 'text', x = 0, y=9.5, label = 'A', size = 5) +
  scale_color_manual(values = c('#00BFC4','#F8766D')) +
  annotate(geom = 'text', x=1, y=9.8, label = '**') +
  annotate(geom = 'text', x=2, y=9.8, label = '**') +
  annotate(geom = 'text', x=3, y=9.8, label = '**') +
  annotate(geom = 'text', x=4, y=9.8, label = '**') +
  annotate(geom = 'text', x=5, y=9.8, label = '**') +
  annotate(geom = 'text', x=6, y=10, label = 'ns') +
  annotate(geom = 'text', x=7, y=10, label = 'ns')
plot_cv2

vec1 <- vector(length = 7)
for (i in 1:7) {
  x = subset(dat1, Day == i)
  print(i)
  y = wilcox.test(data=x, Titer~Category)$p.value
  print(y)
  vec1[i] <- y
  }
vec1
# p-values
# day 1 -- 0.002 **
# day 2 -- 0.002 **
# day 3 -- 0.004 **
# day 4 -- 0.005 **
# day 5 -- 0.002 **
# day 6 -- 0.093 ns
# day 7 -- 0.065 ns


########################################################################
# cI26 vs trained OmpF+
###############
dat2 <- subset(raw.data,
       Category == 'cI26' |
         Category == 'OmpF+' &
         !is.na(Titer))
dat2$Repli <- rep(1:6, times = length(dat2)/6)
dat2$Category <- factor(dat2$Category, levels=c("cI26", "OmpF+"))
plot_cvOplus <- ggplot(data = dat2,
                   aes(x = Day,
                       y = log10(as.numeric(Titer)),
                       group = interaction(Phage_Strain, Rep),
                       color = Category,
                       shape = factor(Repli))) +
  geom_path(size = 1) +
  geom_point(size = 2) +
  labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),
                     limits = c(0,10),
                     expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05), breaks = c(0,1,2,3,4,5,6,7)) +
  scale_linetype_discrete(name='Population') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = 'white'),
        legend.key.width = unit(1, 'cm'),
        #axis.line = element_line(colour = "black"),
        # Change axes ticks
        axis.line = element_line(color='black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black')) +
  # hline is the limit of detection
  geom_hline(yintercept = 2, linetype = 'dashed') +
  # Removes the Population number from the legends
  guides(shape = FALSE) +
  scale_color_manual(values = c('#00BFC4','#00DF0F')) +
  annotate(geom = 'text', x = 0, y=9.5, label = 'B', size = 5) +
  annotate(geom = 'text', x=1, y=9.8, label = '**') +
  annotate(geom = 'text', x=2, y=9.8, label = '**') +
  annotate(geom = 'text', x=3, y=9.8, label = '**') +
  annotate(geom = 'text', x=4, y=9.8, label = '**') +
  annotate(geom = 'text', x=5, y=9.8, label = '**') +
  annotate(geom = 'text', x=6, y=9.8, label = '**') +
  annotate(geom = 'text', x=7, y=9.8, label = '**')
plot_cvOplus

vec2 <- vector(length = 7)
for (i in 1:7) {
  x = subset(dat2, Day == i)
  print(i)
  y = wilcox.test(data=x, Titer~Category)$p.value
  print(y)
  vec2[i] <- y
}
vec2
# p-values
# day 1 -- 0.005 **
# day 2 -- 0.004 **
# day 3 -- 0.002 **
# day 4 -- 0.002 **
# day 5 -- 0.002 **
# day 6 -- 0.002 **
# day 7 -- 0.002 **


########################################################################
# cI26 vs trained OmpF-
###############
dat3 <- subset(raw.data,
               Category == 'cI26' |
                 Category == 'OmpF-')
dat3$Repli <- rep(1:6, times = length(dat3)/6)
dat3$Category <- factor(dat3$Category, levels=c("cI26", "OmpF-"))
dat3 <- subset(dat3, !is.na(Titer))
plot_cvOminus <- ggplot(data = dat3,
                       aes(x = Day,
                           y = log10(as.numeric(Titer)),
                           group = interaction(Phage_Strain, Rep),
                           color = Category,
                           shape = factor(Repli))) +
  geom_path(size = 1) +
  geom_point(size = 2) +
  labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),
                     limits = c(0,10),
                     expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05), breaks = c(0,1,2,3,4,5,6,7)) +
  scale_linetype_discrete(name='Population') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = 'white'),
        legend.key.width = unit(1, 'cm'),
        #axis.line = element_line(colour = "black"),
        # Change axes ticks
        axis.line = element_line(color='black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black')) +
  # hline is the limit of detection
  geom_hline(yintercept = 2, linetype = 'dashed') +
  # Removes the Population number from the legends
  guides(shape = FALSE) +
  scale_color_manual(values = c('#00BFC4','#003BCE')) +
  annotate(geom = 'text', x = 0, y=9.5, label = 'C', size = 5) +
  annotate(geom = 'text', x=1, y=10, label = 'ns') +
  annotate(geom = 'text', x=2, y=9.8, label = '**') +
  annotate(geom = 'text', x=3, y=10, label = 'ns') +
  annotate(geom = 'text', x=4, y=10, label = 'ns') +
  annotate(geom = 'text', x=5, y=9.8, label = '*') +
  annotate(geom = 'text', x=6, y=10, label = 'ns') +
  annotate(geom = 'text', x=7, y=9.8, label = '**')
plot_cvOminus

vec3 <- vector(length = 7)
for (i in 1:7) {
  x = subset(dat3, Day == i)
  print(i)
  y = wilcox.test(data=x, Titer~Category)$p.value
  print(y)
  vec3[i] <- y
}
vec3
# p-values
# day 1 -- 0.589 ns
# day 2 -- 0.004 **
# day 3 -- 0.200 ns
# day 4 -- 0.937 ns
# day 5 -- 0.041 *
# day 6 -- 0.099 ns
# day 7 -- 0.004 **


########################################################################
# 2811 vs trained OmpF+
###############
dat4 <- subset(raw.data,
               Category == '2811' |
                 Category == 'OmpF+')
dat4$Repli <- rep(1:6, times = length(dat4)/6)
dat4 <- subset(dat4, !is.na(Titer))
dat4$Category <- factor(dat4$Category, levels=c("2811", "OmpF+"))
plot_2vOplus <- ggplot(data = dat4,
                        aes(x = Day,
                            y = log10(as.numeric(Titer)),
                            group = interaction(Phage_Strain, Rep),
                            color = Category,
                            shape = factor(Repli))) +
  geom_path(size = 1) +
  geom_point(size = 2) +
  labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),
                     limits = c(0,10),
                     expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0.05,0.05), breaks = c(0,1,2,3,4,5,6,7)) +
  scale_linetype_discrete(name='Population') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = 'white'),
        legend.key.width = unit(1, 'cm'),
        #axis.line = element_line(colour = "black"),
        # Change axes ticks
        axis.line = element_line(color='black'),
        # Change axes text sizes
        axis.text.x = element_text(size=12, color = 'black'),
        axis.text.y = element_text(size=12, color = 'black')) +
  # hline is the limit of detection
  geom_hline(yintercept = 2, linetype = 'dashed') +
  # Removes the Population number from the legends
  guides(shape = FALSE) +
  scale_color_manual(values = c('#F8766D','#00DF0F')) +
  annotate(geom = 'text', x = 0, y=9.5, label = 'D', size = 5) +
  annotate(geom = 'text', x=1, y=10, label = 'ns') +
  annotate(geom = 'text', x=2, y=10, label = 'ns') +
  annotate(geom = 'text', x=3, y=10, label = 'ns') +
  annotate(geom = 'text', x=4, y=10, label = 'ns') +
  annotate(geom = 'text', x=5, y=10, label = 'ns') +
  annotate(geom = 'text', x=6, y=10, label = 'ns') +
  annotate(geom = 'text', x=7, y=10, label = 'ns')
plot_2vOplus

vec4 <- vector(length = 7)
for (i in 1:7) {
  x = subset(dat4, Day == i)
  print(i)
  y = wilcox.test(data=x, Titer~Category)$p.value
  print(y)
  vec4[i] <- y
}
vec4
# p-values
# day 1 -- 0.297 ns
# day 2 -- 1.000 ns
# day 3 -- 0.873 ns
# day 4 -- 0.572 ns
# day 5 -- 0.818 ns
# day 6 -- 0.393 ns
# day 7 -- 0.937 ns


############################################
############################################
plot_matrix <- ggarrange(plot_cv2, plot_cvOplus, plot_cvOminus, plot_2vOplus,
                         legend = 'top',
                         ncol = 2, nrow = 2,
                         widths = 3, heights = 3)
plot_matrix

ggsave(filename = 'all_plots.pdf', plot_matrix,
       width = 8, height = 6)

