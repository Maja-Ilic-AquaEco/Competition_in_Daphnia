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
##  Script 3:                                             ##
##  Daphnia genotype frequencies                          ##
##                                                        ##
##  This script contains the code used to produce Fig. 5  ##
##                                                        ##
############################################################

#============================#
## Packages ----

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(scales)

#============================#
## Set working directory ----

setwd("~/My Documents/Paper Mesocosm experiment with D. longi/Third submission/Data GitHub")

path.figures <- "~/My Documents/Paper Mesocosm experiment with D. longi/Third submission/Figures"

#============================#
## Import data ----

df <- read_excel("Data_4_Daphnia_genotype_frequences_entire_duration.xlsx", sheet = "Data", 
                 col_types = c(rep("numeric",2),"text", rep("numeric",5)))

head(df)
str(df)

## Column names, units and explanation:
# Day	                day	  Day of the experiment
# Mesocosm		              ID of the mesocosm (1-15)
# Inoculum	          µl	  Volume of the natural phytoplankton community from lake Klostersee used for the initial inoculation of the experimental phytoplankton communities
# KL14		                  Frequency of the genotype KL14 (0-1)
# KL83		                  Frequency of the genotype KL83 (0-1)
# KL93		                  Frequency of the genotype KL93 (0-1)
# H_Daphnia		              Shannon-Diversity based on Daphnia genotype frequencies
# spec_rich_Daphnia		      Species richness (number of Daphnia genotypes present in the community)

#============================#
## Summary statistics ----

# First, transform the data set from wide to long format, 
# so that genotype frequencies are in one single columns

df.long <- df %>% gather(key = "Clone", value = "Frequency",
                         c(KL14,KL83,KL93))

# Calculate mean and standard deviation for the genotype frequencies per day
# See also Supplementary Information, Tab. S8

df.summary <- df.long %>% 
  group_by(Day,Clone) %>% 
  summarize(Mean = mean(Frequency,na.rm = T),
            SD = sd(Frequency,na.rm = T))

df.summary$Clone <- factor(df.summary$Clone, levels = c("KL93", "KL83", "KL14"))

#============================#
## Plot ----

# Define color-vision-deficiency friendly colors

palette("okabe-ito")
cvd.cols <- palette("okabe-ito")
show_col(cvd.cols)

gen.cols <- cvd.cols[c(4,2,3)]

# Note: The plot is showing mean frequency (0-1) of each genotype per day
# (average of all 15 mesocosms). White dashed lines represent days on which 
# Daphnia individuals were removed to avoid crowding. Note that the genotype 
# frequencies were not accessed on day 35 (one week upon Daphnia introduction). 

ggplot(df.summary ,aes(x = Day, y = Mean, fill = Clone)) +
  geom_area() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold",
                                    margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold",
                                    margin = margin(0,10,0,0)),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none") +
  scale_fill_manual(values = gen.cols) +
  labs(x = "Day of experiment", y = "Frequency") +
  scale_x_continuous(limits = c(28,98),breaks = c(28,42,56,70,84,98), exp = c(0,0)) +
  scale_y_continuous(limits = c(0,1.0001), breaks = c(0,0.2,0.4,0.6,0.8,1), exp = c(0,0)) +
  annotate("text", x = 50, y = 0.5, label = "KL93", hjust = 0,
           color = "white", size = 4) +
  annotate("text", x = 50, y = 0.25, label = "KL83", hjust = 0,
           color = "white", size = 4) +
  annotate("text", x = 50, y = 0.08, label = "KL14", hjust = 0,
           color = "white", size = 4) +
  geom_vline(aes(xintercept = 56), color = "white", linetype = "dashed") +
  geom_vline(aes(xintercept = 63), color = "white", linetype = "dashed") +
  geom_vline(aes(xintercept = 70), color = "white", linetype = "dashed") +
  geom_vline(aes(xintercept = 77), color = "white", linetype = "dashed") +
  geom_vline(aes(xintercept = 84), color = "white", linetype = "dashed")

ggsave(paste0(path.figures,"/For manuscript - Figure 5.png"), 
              width = 15, height = 10, units = "cm")
