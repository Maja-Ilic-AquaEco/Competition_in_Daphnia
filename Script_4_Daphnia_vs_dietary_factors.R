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
##  Script 4:                                             ##
##  Daphnia vs dietary factors                            ##
##                                                        ##
##  This script contains the code used to produce Fig.    ##
##  6-7 and run linear models and dbRDA based on Daphnia  ##
##  longispina genotype frequencies and dietary factors   ##
##  (diversity or abundance of fatty acids).              ##
##                                                        ##
############################################################

#============================#
## Packages ----

library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(vegan)
library(readxl)
library(scales)

#============================#
## Set working directory ----

setwd("~/My Documents/Paper Mesocosm experiment with D. longi/Third submission/Data GitHub")

path.figures <- "~/My Documents/Paper Mesocosm experiment with D. longi/Third submission/Figures"

#============================#
## Import data ----

df <- read_excel("Data_5_Daphnia_vs_dietary_factors.xlsx",
                 sheet = "Data")


head(df)
str(df)

## Column names, units and explanation:
# Mesocosm		              ID of the mesocosm (1-15)
# Inoculum    µl	          Volume of the natural phytoplankton community from lake Klostersee used for the initial inoculation of the experimental phytoplankton communities
# KL14		                  Frequency of the genotype KL14 (0 - 1) on day 56 of the experiment
# KL83		                  Frequency of the genotype KL83 (0 - 1) on day 56 of the experiment
# KL93    		              Frequency of the genotype KL93 (0 - 1) on day 56 of the experiment
# H_PUFA		                Shannon-Diversity of polyunsaturated fatty acids (PUFAs) measured in the phytoplankton communities on day 28
# PUFA	      % total FA	  Abundance of PUFAs relative to total FAs measured in the phytoplankton communities on day 28
# Omega.3	    % total FA	  Abundance of omega-3 PUFAs relative to total FAs measured in the phytoplankton communities on day 28
# Omega.6	    % total FA	  Abundance of omega-6 PUFAs relative to total FAs measured in the phytoplankton communities on day 28
# ALA	        % total FA	  Abundance of α-linolenic acid ALA relative to total FAs measured in the phytoplankton communities on day 28
# ARA	        % total FA	  Abundance of arachidonic acid ARA relative to total FAs measured in the phytoplankton communities on day 28
# DHA	        % total FA	  Abundance of docosahexaenoic acid DHA relative to total FAs measured in the phytoplankton communities on day 28
# EPA	        % total FA	  Abundance of eicosapentaenoic acid EPA relative to total FAs measured in the phytoplankton communities on day 28

# Note: Shannon-Diversity Indices H_FA and H_PUFA were calculated previously 
# using the function diversity() from the package vegan. Complete data set 
# containing the fatty acid profiles of phytoplankton community in each 
# mesocosm is given in Data_3_Fatty_acid_composition_day_28.xlsx.

#============================#
## Prepare the dataset for further analysis ----

# Gather the data so that all three genotypes are in only one column

# Frequency

df.rel <- df %>% 
  gather(key = "Genotype", value = "Frequency", 
         KL14,KL83,KL93)

# Gather the data so that all FA-derived parameters are in only one column

df.FA <- df.rel %>% 
  gather(key = "FA.parameter", value = "Value",
         H_PUFA, PUFA, Omega.3, Omega.6,
         ALA, ARA, DHA, EPA)

# Claculate summary statistics (Mean and Standard deviation)
# See also Supplementary Information, Tab. S6 and S7

df.summary <- df.FA %>% 
  group_by(FA.parameter) %>% 
  summarize(Mean = mean(Value, na.rm = T),
            SD = sd(Value, na.rm = T)) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

#============================#
## Nest data and run multiple linear regressions ----

lm.FA <- df.FA %>% 
  filter(!is.na(Frequency)) %>% 
  nest(-Genotype,-FA.parameter) %>% 
  mutate(
    fit = map(data, ~lm(Frequency ~ Value, data = .x)),
    tidied = map(fit,tidy),
    glanced = map(fit,glance)
  ) %>% 
  unnest(tidied) 

lm.FA.glanced <- lm.FA %>%
  filter(term == "Value") %>% 
  select(-statistic,-p.value) %>% 
  unnest(glanced)

lm.FA <- left_join(lm.FA,
                   lm.FA.glanced[,c("Genotype","FA.parameter","r.squared","adj.r.squared")],
                   by = c("Genotype","FA.parameter"))

lm.FA[which(lm.FA$p.value < 0.05), "Significant"] <- "Yes"
lm.FA[which(lm.FA$p.value >= 0.05), "Significant"] <- "No"

View(lm.FA %>% 
       filter(term == "Value" & Significant == "Yes") %>% 
       select(Genotype,FA.parameter,Significant))

# Add significance to the original data set

df.new <- left_join(df.FA,
                    lm.FA[lm.FA$term == "Value", c("Genotype","FA.parameter","Significant")],
                    by = c("Genotype","FA.parameter"))

df.new$Significant <- factor(df.new$Significant, levels = c("Yes","No"))

## Create final output table (see Tab. S9)

# Extract intercept (a)

inter.FA <- lm.FA %>% 
  filter(term == "(Intercept)") %>% 
  select(Genotype,FA.parameter,estimate)

names(inter.FA)[which(names(inter.FA) == "estimate")] <- "Intercept"

# Extract slope, p and R2 

slope.FA <- lm.FA %>% 
  filter(term == "Value") %>% 
  select(Genotype,FA.parameter,estimate,p.value,r.squared)

names(slope.FA)[which(names(slope.FA) == "estimate")] <- "Slope"
names(slope.FA)[which(names(slope.FA) == "p.value")] <- "p"
names(slope.FA)[which(names(slope.FA) == "r.squared")] <- "R2"

# Combine

output.lm <- left_join(inter.FA,
                       slope.FA)

# Round all numeric variables

output.lm <- output.lm %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  arrange(Genotype)

write.table(output.lm, "Output_Daphnia_vs_dietary_factors_linear_regressions.csv",
            sep = ";", row.names = F)

#============================#
## Plots ----

# Define color-vision-deficiency friendly colors

palette("okabe-ito")
cvd.cols <- palette("okabe-ito")
show_col(cvd.cols)

gen.cols <- cvd.cols[c(3,2,4)]

#============================#
## H_PUFA ----

g1 <- df.new %>% 
  filter(FA.parameter == "H_PUFA") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T,
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0.3,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(1,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(paste(bolditalic("H'")[bold("PUFA")]))) +
  scale_y_continuous(limits = c(-10,100)) +
  coord_cartesian(xlim = c(min(df.new$Value[df.new$FA.parameter == "H_PUFA"]),
                           max(df.new$Value[df.new$FA.parameter == "H_PUFA"])),
                  ylim = c(0,100)) 

g1

#============================#
## PUFA ----

g2 <- df.new %>% 
  filter(FA.parameter == "PUFA") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0.3,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(1,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = "PUFA (% total FA)") +
  coord_cartesian(xlim = c(20,45), ylim = c(0,100))

g2

#============================#
## omega3-PUFA ----

g3 <- df.new %>% 
  filter(FA.parameter == "Omega.3") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0.3,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(1,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(bold(omega*"3-PUFA (% total FA)"))) +
  coord_cartesian(xlim = c(12,32), ylim = c(0,100))

g3

#============================#
## omega6-PUFA ----

g4 <- df.new %>% 
  filter(FA.parameter == "Omega.6") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(2,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(paste(bold(omega*"6-PUFA (% total FA)")))) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(6,14))

g4

#============================#
## ALA ----

g5 <- df.new %>% 
  filter(FA.parameter == "ALA") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0.3,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(1,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(paste(bold("ALA, C18:3 "*omega*"3 (% total FA)")))) +
  coord_cartesian(xlim = c(7,27), ylim = c(0,100))

g5

#============================#
## ARA ----

g6 <- df.new %>% 
  filter(FA.parameter == "ARA") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0.3,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(1,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(paste(bold("ARA, C20:4 "*omega*"6 (% total FA)")))) +
  coord_cartesian(xlim = c(0,1.5), ylim = c(0,100))

g6

#============================#
## EPA ----

g7 <- df.new %>% 
  filter(FA.parameter == "EPA") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(2,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(paste(bold("EPA, C20:5 "*omega*"3 (% total FA)")))) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,5))

g7

#============================#
## DHA ----

g8 <- df.new %>% 
  filter(FA.parameter == "DHA") %>% 
  ggplot(., aes(x = Value, y = Frequency*100, fill = Genotype, color = Genotype)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, 
              aes(alpha = Significant, linetype = Significant)) +
  scale_fill_manual(values = gen.cols) +
  scale_color_manual(values = gen.cols) +
  scale_alpha_manual(values = c(0,0), name = "p < 0.05") +
  scale_linetype_manual(values = c(2,2), name = "p < 0.05") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(10,0,0,0)),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(8,1,4,6)) +
  labs(x = expression(paste(bold("DHA, C22:6 "*omega*"3 (% total FA)")))) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,6))

g8

#============================#
## Plot 1 ----

plot1 <- plot_grid(g1,g2,g3,g4,
                   ncol = 2,
                   align = "hv",
                   labels = c("A","B","C","D"),
                   label_fontface = "bold",
                   label_x = 0.915, label_y = 0.975)

y.label <- textGrob(expression(bolditalic("          D. longispina")*bold("genotype (%)")), 
                    gp = gpar(fontsize = 14), rot = 90)

plot1 <- arrangeGrob(plot1, left = y.label)

grid.arrange(plot1)

ggsave(paste0(path.figures,"/For manuscript - Figure 6.png"),
       grid.arrange(plot1),
       width = 17.5, height = 17, units = "cm")

#============================#
## Plot 2 ----

plot2 <- plot_grid(g5,g6,g7,g8,
                   ncol = 2,
                   align = "hv",
                   labels = c("A","B","C","D"),
                   label_fontface = "bold",
                   label_x = 0.915, label_y = 0.975)

y.label <- textGrob(expression(bolditalic("          D. longispina")*bold("genotype (%)")), 
                    gp = gpar(fontsize = 14), rot = 90)

plot2 <- arrangeGrob(plot2, left = y.label)

grid.arrange(plot2)

ggsave(paste0(path.figures,"/For manuscript - Figure 7.png"),
       grid.arrange(plot2),
       width = 17.5, height = 17, units = "cm")

