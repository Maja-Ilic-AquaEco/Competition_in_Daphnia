############################################################
##                                                        ##
##  Essential dietary fatty acids affect intraspecific    ## 
##  competition in herbivorous zooplankton                ## 
##                                                        ##
##  Authors of the study:                                 ##
##  Maja Ilic, Sina Brehm, Maria Stockenreiter,           ##
##  Eric von Elert, Patrick Fink                          ##
##                                                        ##
##  Corresponding author & author of the scripts:         ##
##  Maja Ilic                                             ##
##  maja.ilic.bio@gmail.com                               ##
##                                                        ##
##  Data and scripts can be found here:                   ##
##  github.com/Maja-Ilic-AquaEco/Competition_in_Daphnia   ##
##                                                        ##
##  Script 1:                                             ##
##  Chlorophyll a vs Daphnia                              ##
##                                                        ##
##  This script contains the code used to produce Fig. 1  ##
##                                                        ##
############################################################

#============================#
## Packages ----

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(RColorBrewer)
library(cowplot)

#============================#
## Set working directory ----

setwd("~/My Documents/Paper Mesocosm experiment with D. longi/Third submission/Data GitHub")

path.figures <- "~/My Documents/Paper Mesocosm experiment with D. longi/Third submission/Figures"

#============================#
## Import data ----

df <- read_excel("Data_1_Chlorophyll_a_Daphnia_entire_duration.xlsx", sheet = "Data",
                 col_types = c(rep("numeric",2),rep("text",2),rep("numeric",4)))

head(df)
str(df)

## Column names, units and explanation:
# Day	              day	      Day of the experiment (Start = Day 0)
# Sampling		                Sampling event (numbered)
# Mesocosm		                ID of the mesocosm (1-15)
# Inoculum	        µl	      Volume of the natural phytoplankton community from lake Klostersee used for the initial inoculation of the experimental phytoplankton communities
# TChla	            μg L-1	  Total chlorophyll a
# Adult_per_L	      ind L-1	  Number of adult Daphnia per L
# Juvenile_per_L	  ind L-1	  Number of juvenile Daphnia per L
# Daphnia_per_L	    ind L-1	  Total number of Daphnia per L (sum of adult and juvenile Daphnia)

#============================#
## Summary statistics ----

# Calculate mean and standard deviation for Chlorophyll a and Daphnia abundancies
# See also Supplementary Information, Tab. S4

df.summary <- df %>% group_by(Day) %>%
  summarize(Mean_TChla = mean(TChla, na.rm = T),
            SD_TChla = sd(TChla, na.rm = T),
            Mean_Daphnia = mean(Daphnia_per_L, na.rm = T),
            SD_Daphnia = sd(Daphnia_per_L, na.rm = T),
            Mean_Adult = mean(Adult_per_L, na.rm = T),
            SD_Adult = sd(Adult_per_L, na.rm = T),
            Mean_Juvenile = mean(Juvenile_per_L, na.rm = T),
            SD_Juvenile = sd(Juvenile_per_L, na.rm = T))

#============================#
## Plots ----

# Define colors and labels for the plot 

cols <- c(rgb(51,160,44,max = 255),
          rgb(31,120,180,max = 255))

ylab1 <- expression(bold("Chlorophyll ")*bolditalic("a")*bold(" (µg "*L^"-1"*")"))
ylab2 <- expression(bolditalic("D. longispina")*bold(" (ind "*L^"-1"*")"))

## Chlorophyll a ----

p1 <- ggplot(df.summary, aes(x = Day, y = Mean_TChla)) +
  geom_rect(aes(xmin = -Inf, xmax = 28, ymin = -Inf, ymax = Inf), fill = "grey70") +
  geom_line(color = cols[1], linewidth = 1.2) +
  geom_errorbar(aes(ymin = Mean_TChla - SD_TChla, ymax = Mean_TChla + SD_TChla),
                color = cols[1], width = 1) +
  geom_point(shape = 21, color = cols[1], fill = cols[1], size = 3) +
  theme_bw() +
  labs(x = "Day", y = ylab1) +
  scale_x_continuous(limits = c(-1.5,99.5), breaks = c(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98)) +
  scale_y_continuous(limits = c(-2,20), breaks = c(0,5,10,15,20)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

p1

## Daphnia abundances ----

# First, reorganize the data

# All means 

daphnia.mean <- df.summary %>% 
  select(Day, Mean_Daphnia, Mean_Adult, Mean_Juvenile) %>% 
  gather(key = "Group", value = "Mean", -Day)

daphnia.mean$Group <- gsub("Mean_", "", 
                           as.character(daphnia.mean$Group))

# All standard deviations

daphnia.sd <- df.summary %>% 
  select(Day, SD_Daphnia, SD_Adult, SD_Juvenile) %>% 
  gather(key = "Group", value = "SD", -Day)

daphnia.sd$Group <- gsub("SD_", "", 
                         as.character(daphnia.sd$Group))

# Combine

df.daphnia <- left_join(daphnia.mean,
                        daphnia.sd,
                        by = c("Day","Group"))

df.daphnia$Group <- as.factor(df.daphnia$Group)
df.daphnia$Group <- factor(df.daphnia$Group, levels = c("Daphnia","Adult","Juvenile"))

# Define colors

cols.daphnia <- c("black",colorRampPalette(brewer.pal(8, "Paired"))(8)[(c(2,1))])

adult.juvenile <- df.daphnia %>% 
  filter(Group != "Daphnia")

total.daphnia <- df.daphnia %>% 
  filter(Group == "Daphnia")

p2 <- ggplot(total.daphnia, aes(x = Day, y = Mean)) +
  geom_rect(aes(xmin = -Inf, xmax = 28, ymin = -Inf, ymax = Inf), fill = "grey70", color = NA) +
  geom_area(data = adult.juvenile, aes(x = Day, y = Mean, color = Group, fill = Group)) +
  geom_vline(xintercept = c(56,63,70,77,84), linetype = "dashed", color = "grey70") +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 2) +
  geom_point(shape = 21, size = 3, color = "black", fill = "white") +
  theme_bw() +
  labs(x = "Day", y = ylab2) +
  scale_color_manual(values = cols.daphnia[2:3], 
                     name = "") +
  scale_fill_manual(values = cols.daphnia[2:3], 
                    name = "") +
  scale_x_continuous(limits = c(-1.5,99.5), breaks = c(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98)) +
  scale_y_continuous(limits = c(-2,80), breaks = c(0,20,40,60,80)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1, "cm"),
        legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

p2

## Join plots ----

p1.new <- p1 + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())

plot.final <- plot_grid(p1.new,p2,
                        align = "v",
                        nrow = 2,
                        labels = c("A","B"),
                        label_x = c(0.945,0.945),
                        label_y = c(0.965,0.965),
                        rel_heights = c(0.75,1))

ggsave(paste0(path.figures,"/For manuscript - Figure 1.png"), 
       plot.final, width = 7, height = 8)
