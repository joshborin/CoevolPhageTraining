
# This rather beastly script is going to generate a 9 panel Figure 5 from several different data sets
# I have split up the script to process data for each of the different phage populations
# NOTE, the script needs to be run sequentially to work properly

library(reshape2)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(vegan)
library(stringr)

setwd("CoevolReplay_FastResistance/Manuscript/Github_repo/Figure5")
####################################################################################
#############################   POPULATION 1 #######################################
####################################################################################
rawdata_P1 <- read.csv(file = "DataFormatted_P1.csv") 
platemap_P1 <- read.csv(file = "2020-08-14_P1_platemap.csv") 

reshaped <- melt(rawdata_P1, id=c("Time"), variable.name="Well", 
                 value.name="OD600")
annotated1 <- inner_join(reshaped, platemap_P1, by="Well")
write.csv(annotated1, "2020-09-03_annotated_P1.csv")

# First, check Blank Controls
ggplot(data= subset(annotated1, Bacteria == "Blank"),
       mapping = aes(x=Time/3600, y=OD600,
                     group=Well,
                     color = Well)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 22)) +
  scale_y_continuous(limits = c(0, 0.7))

annotated1_nb <- subset(annotated1, Bacteria !="Blank")

# Plot bacteria without phage as positive control for growth
bact_data <- subset(annotated1_nb, Phage == 'None')
plot_bact <- ggplot(data=bact_data,
       mapping = aes(x=Time/3600, y=OD600,
                     group=(Well))) +
  geom_line(size = 1) +
  facet_wrap(Bacteria ~ ., ncol = 2) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_bact)
# All wells with phage
data_2 <- subset(annotated1_nb, Phage != 'None' & Phage != ' Blank')
plot_2811 <- ggplot(data=data_2,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=Well,
                                  color=Phage)) +
  geom_line(size = 1) +
  facet_wrap(Dil ~ ., ncol = 2) +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_2811)

# Now equalize initial phage densities by choosing appropriate wells/dilutions
# This is based on spot titering of phages on REL606
dat1 <- subset(annotated1, Phage == '2811' & Dil == '0.0125')
dat2 <- subset(annotated1, Phage == '2811-P1-T25' & Dil == 'undil')
dat3 <- subset(annotated1, Phage == 'None')
dat4 <- subset(annotated1, Phage == 'Blank')

dat_sub <- rbind(dat1, dat2, dat3, dat4)
plot_sub <- ggplot(data=dat_sub,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=Well,
                                  color=interaction(Phage,Dil))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_sub)

# Now subset the raw data from my annotated dataset
# And transform the data so that I can conduct a PCA to look diff in suppression

vec <- unique(dat_sub$Status)
vec
dat_anc <- subset(dat_sub, Status == vec[1])
dat_contemp <- subset(dat_sub, Status == vec[2])
dat_none <- subset(dat_sub, Status == vec[3])
dat_blank <- subset(dat_sub, Status == vec[4])

unique(dat_sub$Well)
length(unique(dat_sub$Well))

########
mat_sub <- matrix(nrow = 24, ncol = 241)
for (i in 1:length(unique(dat_sub$Well))) {
  well <- unique(dat_sub$Well)[i]
  print(well)
  dat <- subset(dat_sub, Well == well)
  mat_sub[i,] <- dat[,3]
}
pca_sub <- rda(mat_sub[,1:241])
biplot(pca_sub)
summary(pca_sub)

a <- summary(pca_sub)
site.scores <- a$sites
# save these scores (=Principal components) and corresponding strains to data frame
pc <- as.data.frame(site.scores[,c(1,2,3,4,5,6)])
length(unique(dat_sub$Phage))
Strains <- data.frame('Strain' = c(rep('2811', times = 8),
                                   rep('2811-P1-T25', times = 8),
                                   rep('None', times = 4),
                                   rep('Blank', times = 4)))
colnames(Strains) <- c('Strain')
dater <- bind_cols(Strains, pc)
## Plot
pc.plot <-  ggplot(data = dater, aes(x = PC1, y = PC2, color = Strain)) +
  geom_point() +
  stat_ellipse(geom='polygon', aes(fill = Strain),
               alpha = 0.2,
               show.legend = FALSE,
               level = 0.95) +
  theme_classic()
pc.plot
box.plot <- ggplot(data = dater, aes(x = Strain, y = PC1, color = Strain)) +
  geom_boxplot() +
  theme_classic()
print(box.plot)

t.test(PC1~Strain, data = subset(dater, Strain == '2811' | Strain == '2811-P1-T25'))
# 2811 and 2811+ P1 are sig different in PC1 // t = 3.7531 // p=0.004)
#####################################################################################
#### Plots for Figures ####
P1_dat <- subset(dat_sub, Phage == '2811' | Phage == '2811-P1-T25')
A <- ggplot(data=P1_dat,
                 mapping = aes(x=Time/3600, y=OD600,
                               group=Well,
                               color=Phage)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = c(0.3,0.6),
        legend.key=element_blank(),
        legend.background=element_blank(),
        # legend.title = element_text(size = 12),
        # legend.text = element_text(size = 12),
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = 'Time (Hours)', y = 'Bacterial Density (OD600)') +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="A", size = 6) +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A')
print(A)
ggsave(plot=A,
       filename="Figure5A.pdf",
       width = 6, height = 4)

B <- ggplot(data = subset(dater, Strain == '2811' | Strain == '2811-P1-T25'),
                 aes(x = Strain, y = PC1, color = Strain)) +
  geom_boxplot() +
  geom_point() +
  geom_bracket(inherit.aes = FALSE, xmin = '2811', xmax = '2811-P1-T25', y.position = 0.18,
    label = "**") +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A') +
  scale_y_continuous(limits = c(-0.4, 0.24)) +
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="B", size = 6) +
  labs(x = 'Phage', y = 'Resistance (PC1 89.94%)')
print(B)
ggsave(B, file = 'Figure5B.pdf', width = 3, height = 6)

# EOP Figure
EOP <- read.csv(file = 'Consolidated_EOP.csv')

P1_EOP <- subset(EOP, Host == 'P1-T25-1-P')
t.test(data=P1_EOP, EOP~Phage)
# Significantly different EOP, p = 0.022, t = -3.799
# mean 2811 is 0.14
# mean 2811-P1-T25 is 0.45
C <- ggplot(P1_EOP, aes(x = Phage, y = as.numeric(EOP), color = Phage)) +
  geom_boxplot() +
  geom_point() + 
  geom_bracket(inherit.aes = FALSE, xmin = '2811', xmax = '2811-P1-T25', y.position = 0.68,
               label = "*") +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A') +
  scale_y_continuous(limits = c(0, 0.75)) +
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="C", size = 6) +
  labs(x = 'Phage', y = 'Efficiency of Plating (EOP)')
print(C)
ggsave(C, file = 'Figure5C.pdf', width = 3, height = 6)

####################################################################################
#############################   POPULATION 2 #######################################
####################################################################################
rawdata_P2 <- read.csv(file = "DataFormatted_P2.csv") 
platemap_P2 <- read.csv(file = "2020-08-21_P2_platemap.csv") 

reshaped <- melt(rawdata_P2, id=c("Time"), variable.name="Well", 
                 value.name="OD600")
annotated2 <- inner_join(reshaped, platemap_P2, by="Well")
write.csv(annotated2, "2020-09-04_annotated_P2.csv")

# First, check Blank Controls
ggplot(data= subset(annotated2, Bacteria == "Blank"),
       mapping = aes(x=Time/3600, y=OD600,
                     group=Well,
                     color = Well)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 22)) +
  scale_y_continuous(limits = c(0, 0.7))
annotated2_nb <- subset(annotated2, Bacteria !="Blank")

# Plot bacteria without phage as positive control for growth
bact_data <- subset(annotated2_nb, Phage == 'None')
plot_bact <- ggplot(data=bact_data,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=(Well))) +
  geom_line(size = 1) +
  facet_wrap(Bacteria ~ ., ncol = 2) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_bact)
# All wells with phage
data_2 <- subset(annotated2_nb, Phage != 'None' & Phage != ' Blank')
plot_2811 <- ggplot(data=data_2,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=Well,
                                  color=Phage)) +
  geom_line(size = 1) +
  facet_wrap(Dil ~ ., ncol = 2) +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_2811)

# Now equalize initial phage densities by choosing appropriate wells/dilutions
# This is based on spot titering of phages on REL606
dat1 <- subset(annotated2, Phage == '2811' & Dil == '0.003125')
dat2 <- subset(annotated2, Phage == '2811-P2-T15' & Dil == 'undil')
dat3 <- subset(annotated2, Phage == 'None')
dat4 <- subset(annotated2, Phage == 'Blank')

dat_sub <- rbind(dat1, dat2, dat3, dat4)
plot_sub <- ggplot(data=dat_sub,
                   mapping = aes(x=Time/3600, y=OD600,
                                 group=Well,
                                 color=interaction(Phage,Dil))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_sub)
# Now subset the raw data from my annotated dataset
# And transform the data so that I can conduct a PCA to look diff in suppression
vec <- unique(dat_sub$Status)
vec
dat_anc <- subset(dat_sub, Status == vec[1])
dat_contemp <- subset(dat_sub, Status == vec[2])
dat_none <- subset(dat_sub, Status == vec[3])
dat_blank <- subset(dat_sub, Status == vec[4])

unique(dat_sub$Well)
length(unique(dat_sub$Well))

########
mat_sub <- matrix(nrow = 24, ncol = 241)
for (i in 1:length(unique(dat_sub$Well))) {
  well <- unique(dat_sub$Well)[i]
  print(well)
  dat <- subset(dat_sub, Well == well)
  mat_sub[i,] <- dat[,3]
}
pca_sub <- rda(mat_sub[,1:241])
biplot(pca_sub)
summary(pca_sub)

a <- summary(pca_sub)
site.scores <- a$sites
# save these scores (=Principal components) and corresponding strains to data frame
pc <- as.data.frame(site.scores[,c(1,2,3,4,5,6)])
length(unique(dat_sub$Phage))
Strains <- data.frame('Strain' = c(rep('2811', times = 8),
                                   rep('2811-P2-T15', times = 8),
                                   rep('None', times = 4),
                                   rep('Blank', times = 4)))
colnames(Strains) <- c('Strain')
dater <- bind_cols(Strains, pc)
## Plot
pc.plot <-  ggplot(data = dater, aes(x = PC1, y = PC2, color = Strain)) +
  geom_point() +
  stat_ellipse(geom='polygon', aes(fill = Strain),
               alpha = 0.2,
               show.legend = FALSE,
               level = 0.95) +
  theme_classic()
pc.plot
box.plot <- ggplot(data = dater, aes(x = Strain, y = PC1, color = Strain)) +
  geom_boxplot() +
  theme_classic()
print(box.plot)
t.test(PC1~Strain, data = subset(dater, Strain == '2811' | Strain == '2811-P2-T15'))
# 2811 and 2811+ P2 are sig different in PC1 / t = 3.0404 / p=0.017)
# Note that for P2 phage are more dense (not sig) and sig more suppression
#  however when 2811 more dense (not sig) there is not sig more suppression
#####################################################################################
#### Plots for Figures ####
P2_dat <- subset(dat_sub, Phage == '2811' | Phage == '2811-P2-T15')
D <- ggplot(data=P2_dat,
            mapping = aes(x=Time/3600, y=OD600,
                          group=Well,
                          color=Phage)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = c(0.3,0.6),
        legend.key=element_blank(),
        legend.background=element_blank(),
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = 'Time (Hours)', y = 'Bacterial Density (OD600)') +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="D", size = 6) +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A')
print(D)
ggsave(plot=D,
       filename="Figure5D.pdf",
       width = 6, height = 4)

E <- ggplot(data = subset(dater, Strain == '2811' | Strain == '2811-P2-T15'),
            aes(x = Strain, y = PC1, color = Strain)) +
  geom_boxplot() +
  geom_point() +
  geom_bracket(inherit.aes = FALSE, xmin = '2811', xmax = '2811-P2-T15', y.position = 0.08,
               label = "*") +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A') +
  scale_y_continuous(limits = c(-0.4, 0.2)) +
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="E", size = 6) +
  labs(x = 'Phage', y = 'Resistance (PC1 90.38%)')
print(E)
ggsave(E, file = 'Figure5E.pdf', width = 3, height = 6)

# EOP Figure
P2_EOP <- subset(EOP, Host == 'P2-T15-1-P')
t.test(data=P2_EOP, EOP~Phage)
# Significantly different EOP, p = 0.045, t = -3.0109
# mean 2811 is 0.30
# mean 2811-P1-T25 is 0.87

pF <- ggplot(P2_EOP, aes(x = Phage, y = as.numeric(EOP), color = Phage)) +
  geom_boxplot() +
  geom_point() + 
  geom_bracket(inherit.aes = FALSE, xmin = '2811', xmax = '2811-P2-T15', y.position = 1.45,
               label = "*") +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A') +
  scale_y_continuous(limits = c(0, 1.8)) +
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="F", size = 6) +
  labs(x = 'Phage', y = 'Efficiency of Plating (EOP)')
print(pF)
ggsave(pF, file = 'Figure5F.pdf', width = 3, height = 6)

####################################################################################
#############################   POPULATION 3 #######################################
####################################################################################
rawdata_P3 <- read.csv(file = "DataFormatted_P3.csv") 
platemap_P3 <- read.csv(file = "2020-08-20_P3_platemap.csv") 

reshaped <- melt(rawdata_P3, id=c("Time"), variable.name="Well", 
                 value.name="OD600")
annotated3 <- inner_join(reshaped, platemap_P3, by="Well")
write.csv(annotated3, "2020-09-04_annotated_P3.csv")

# First, check Blank Controls
ggplot(data= subset(annotated3, Bacteria == "Blank"),
       mapping = aes(x=Time/3600, y=OD600,
                     group=Well,
                     color = Well)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 22)) +
  scale_y_continuous(limits = c(0, 0.7))

annotated3_nb <- subset(annotated3, Bacteria !="Blank")

# Plot bacteria without phage as positive control for growth
bact_data <- subset(annotated3_nb, Phage == 'None')
plot_bact <- ggplot(data=bact_data,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=(Well))) +
  geom_line(size = 1) +
  facet_wrap(Bacteria ~ ., ncol = 2) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_bact)
# All wells with phage
data_2 <- subset(annotated3_nb, Phage != 'None' & Phage != ' Blank')
plot_2811 <- ggplot(data=data_2,
                    mapping = aes(x=Time/3600, y=OD600,
                                  group=Well,
                                  color=Phage)) +
  geom_line(size = 1) +
  facet_wrap(Dil ~ ., ncol = 2) +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_2811)

# Now equalize initial phage densities by choosing appropriate wells/dilutions
# This is based on spot titering of phages on REL606
dat1 <- subset(annotated3, Phage == '2811' & Dil == '0.025')
dat2 <- subset(annotated3, Phage == '2811-P3-T5' & Dil == '-1')
dat3 <- subset(annotated3, Phage == 'None')
dat4 <- subset(annotated3, Phage == 'Blank')

dat_sub <- rbind(dat1, dat2, dat3, dat4)
plot_sub <- ggplot(data=dat_sub,
                   mapping = aes(x=Time/3600, y=OD600,
                                 group=Well,
                                 color=interaction(Phage,Dil))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.7)) +
  scale_color_brewer(type = 'qual')
print(plot_sub)

# Now subset the raw data from my annotated dataset
# And transform the data so that I can conduct a PCA to look diff in suppression

vec <- unique(dat_sub$Status)
vec
dat_anc <- subset(dat_sub, Status == vec[1])
dat_contemp <- subset(dat_sub, Status == vec[2])
dat_none <- subset(dat_sub, Status == vec[3])
dat_blank <- subset(dat_sub, Status == vec[4])

unique(dat_sub$Well)
length(unique(dat_sub$Well))

########
mat_sub <- matrix(nrow = 24, ncol = 241)
for (i in 1:length(unique(dat_sub$Well))) {
  well <- unique(dat_sub$Well)[i]
  print(well)
  dat <- subset(dat_sub, Well == well)
  mat_sub[i,] <- dat[,3]
}
pca_sub <- rda(mat_sub[,1:241])
biplot(pca_sub)
summary(pca_sub)

a <- summary(pca_sub)
site.scores <- a$sites
# save these scores (=Principal components) and corresponding strains to data frame
pc <- as.data.frame(site.scores[,c(1,2,3,4,5,6)])
length(unique(dat_sub$Phage))
Strains <- data.frame('Strain' = c(rep('2811', times = 8),
                                   rep('2811-P3-T5', times = 8),
                                   rep('None', times = 4),
                                   rep('Blank', times = 4)))
colnames(Strains) <- c('Strain')
dater <- bind_cols(Strains, pc)
## Plot
pc.plot <-  ggplot(data = dater, aes(x = PC1, y = PC2, color = Strain)) +
  geom_point() +
  stat_ellipse(geom='polygon', aes(fill = Strain),
               alpha = 0.2,
               show.legend = FALSE,
               level = 0.95) +
  theme_classic()
pc.plot
box.plot <- ggplot(data = dater, aes(x = Strain, y = PC1, color = Strain)) +
  geom_boxplot() +
  theme_classic()
print(box.plot)
t.test(PC1~Strain, data = subset(dater, Strain == '2811' | Strain == '2811-P3-T5'))
# 2811 and 2811+ P3 are sig different in PC1 / t = 2.7621 / p=0.017
# Note that for P3 phage is sig less dense and has sig more suppression!
#####################################################################################
#### Plots for Figures ####
P3_dat <- subset(dat_sub, Phage == '2811' | Phage == '2811-P3-T5')
G <- ggplot(data=P3_dat,
            mapping = aes(x=Time/3600, y=OD600,
                          group=Well,
                          color=Phage)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Legend alterations
        legend.position = c(0.3,0.6),
        legend.key=element_blank(),
        legend.background=element_blank(),
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = 'Time (Hours)', y = 'Bacterial Density (OD600)') +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="G", size = 6) +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A')
print(G)
ggsave(plot=G,
       filename="Figure5G.pdf",
       width = 6, height = 4)

H <- ggplot(data = subset(dater, Strain == '2811' | Strain == '2811-P3-T5'),
            aes(x = Strain, y = PC1, color = Strain)) +
  geom_boxplot() +
  geom_point() +
  geom_bracket(inherit.aes = FALSE, xmin = '2811', xmax = '2811-P3-T5', y.position = -0.125,
               label = "*") +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A') +
  scale_y_continuous(limits = c(-0.35, -0.08)) +
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x=-Inf, y=Inf, hjust = -1, vjust = 2, label="H", size = 6) +
  labs(x = 'Phage', y = 'Resistance (PC1 98.12%)')
print(H)
ggsave(H, file = 'Figure5H.pdf', width = 3, height = 6)

# EOP Figure
P3_EOP <- subset(EOP, Host == 'P3-T5-1-P')
t.test(data=P3_EOP, EOP~Phage)
# Significantly different EOP, p = 0.037, t = -3.5377
# mean 2811 is 0.18
# mean 2811-P1-T25 is 0.68

I <- ggplot(P3_EOP, aes(x = Phage, y = as.numeric(EOP), color = Phage)) +
  geom_boxplot() +
  geom_point() + 
  geom_bracket(inherit.aes = FALSE, xmin = '2811', xmax = '2811-P3-T5', y.position = 1.2,
               label = "*") +
  scale_color_viridis_d(begin = 0.68, end = 0.3, option = 'A') +
  scale_y_continuous(limits = c(0, 1.5)) +
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
        panel.border = element_rect(fill = NA, colour = "black")) +
  annotate('text', x=-Inf, y=Inf, hjust = -3, vjust = 2, label="I", size = 6) +
  labs(x = 'Phage', y = 'Efficiency of Plating (EOP)')
print(I)
ggsave(I, file = 'Figure5I.pdf', width = 3, height = 6)


###############################################################################
###############################################################################
###############################################################################

plot_matrix <- ggarrange(A, B, C,
                         D, E, pF,
                         G, H, I,
                         common.legend = FALSE,
                         ncol = 3, nrow = 3,
                         widths = c(5,4,4), heights = c(4,4,4))
plot_matrix

ggsave(filename = 'Figure5A-I.pdf', plot_matrix,
       width = 8, height = 9)










