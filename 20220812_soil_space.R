#--------------------------------------------------------------------------#
#     .______.
#   __| _/\_ |__    ____
#  / __ |  | __ \ _/ ___\
# / /_/ |  | \_\ \\  \___
# \____ |  |___  / \___  >
#      \/      \/      \/
#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2022-03-16 10:13:57
# @Last Modified by:
# @Last Modified time: 2022-08-12 14:54:52
# @Description: Anova Tables
#--------------------------------------------------------------------------#
# ==========================================================================
#  preparation
# ==========================================================================
# cleaning memory
cat("\014")
rm(list=ls())
gc()

# Loading packages
# You have loaded plyr after dplyr - this is likely to cause problems.
# If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
# library(plyr); library(dplyr)

library(broom)
library(car)
library(ciTools)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpubr)
library(nlme)
library(patchwork)
library(purrr)
library(readr)
library(readxl)
library(skimr)
library(tidyr)
library(dplyr)

# set the current working directory
setwd("G:/我的坚果云/高硕硕/NEW_RESULT")
getwd()

# Loading data of soil space
data4gao <- read.csv("./data4gao.csv")
skim(data4gao)
View(data4gao)

# Convert input to a factor
data4gao$Function_group <- as.factor(data4gao$Function_group)
data4gao$Species_code <- as.factor(data4gao$Species_code)
data4gao$Space_code <- as.factor(data4gao$Space_code)
str(data4gao)

# ==========================================================================
# Data analyses - 01
# ==========================================================================
# ref:
# https://fukamilab.github.io/BIO202/04-B-binary-data.html#glm_for_proportional_data

# data range
range(data4gao$Total_mass, na.rm = T)
range(data4gao$Shoot_mass, na.rm = T)
range(data4gao$Root_mass, na.rm = T)
range(data4gao$Rmass_per_vol, na.rm = T)
range(data4gao$RM_ratio, na.rm = T)

# normality test
ggdensity.Total_mass <- ggdensity(log(data4gao$Total_mass), main ="density plot", xlab ="Total_mass")
ggdensity.Shoot_mass <- ggdensity(log(data4gao$Shoot_mass), main ="density plot", xlab ="Shoot_mass")
ggdensity.Root_mass <- ggdensity(log(data4gao$Root_mass), main ="density plot", xlab ="Root_mass")
ggdensity.Rmass_per_vol <- ggdensity(log(data4gao$Rmass_per_vol), main ="density plot", xlab ="Rmass_per_vol")
ggdensity.RM_ratio <- ggdensity(sqrt(data4gao$RM_ratio), main ="density plot", xlab ="RM_ratio")

ggdensity.Total_mass + ggdensity.Shoot_mass + ggdensity.Root_mass + ggdensity.Rmass_per_vol + ggdensity.RM_ratio

# qq plot
ggqqplot.Total_mass <- ggqqplot(log(data4gao$Total_mass), main ="density plot", xlab ="Total_mass")
ggqqplot.Shoot_mass <- ggqqplot(log(data4gao$Shoot_mass), main ="density plot", xlab ="Shoot_mass")
ggqqplot.Root_mass <- ggqqplot(log(data4gao$Root_mass), main ="density plot", xlab ="Root_mass")
ggqqplot.Rmass_per_vol <- ggqqplot(log(data4gao$Rmass_per_vol), main ="density plot", xlab ="Rmass_per_vol")
ggqqplot.RM_ratio <- ggqqplot(sqrt(data4gao$RM_ratio), main ="density plot", xlab ="RM_ratio")

ggqqplot.Total_mass + ggqqplot.Shoot_mass + ggqqplot.Root_mass + ggqqplot.Rmass_per_vol + ggqqplot.RM_ratio

# ==========================================================================
# Data summarization
# ==========================================================================
summary.data4gao <- data4gao %>% group_by(Function_group, Species_code) %>% dplyr::summarise(Total_mass.mean = mean(Total_mass, na.rm = TRUE),
                                                                             Shoot_mass.mean = mean(Shoot_mass, na.rm = TRUE),
                                                                             Root_mass.mean = mean(Root_mass, na.rm = TRUE),
                                                                             Rmass_per_vol.mean = mean(Rmass_per_vol, na.rm = TRUE),
                                                                             RM_ratio.mean = mean(RM_ratio, na.rm = TRUE))

# save data summarization
write.csv(summary.data4gao, "summary.data4gao.csv", row.names = FALSE)

# ==========================================================================
# Data transformation
# ==========================================================================
data4gao$Total_mass.tf    <- log(data4gao$Total_mass)
data4gao$Shoot_mass.tf    <- log(data4gao$Shoot_mass)
data4gao$Root_mass.tf     <- log(data4gao$Root_mass)
data4gao$Rmass_per_vol.tf <- log(data4gao$Rmass_per_vol)
data4gao$RM_ratio.tf      <- sqrt(data4gao$RM_ratio)

# library(car)
# leveneTest_Total_mass.tf    <- leveneTest(Total_mass.tf ~ Function_group * Space_code, data = data4gao)
# leveneTest_Shoot_mass.tf    <- leveneTest(Shoot_mass.tf ~ Function_group * Space_code, data = data4gao)
# leveneTest_Root_mass.tf     <- leveneTest(Root_mass.tf ~ Function_group * Space_code, data = data4gao)
# leveneTest_Rmass_per_vol.tf <- leveneTest(Rmass_per_vol.tf ~ Function_group * Space_code, data = data4gao)
# leveneTest_RM_ratio.tf      <- leveneTest(RM_ratio.tf ~ Function_group * Space_code, data = data4gao)

# leveneTest_Total_mass.tf
# leveneTest_Shoot_mass.tf
# leveneTest_Root_mass.tf
# leveneTest_Rmass_per_vol.tf
# leveneTest_RM_ratio.tf

# ==========================================================================
#
# Functional groups scale
# nested ANOVAs via lme function
#
# ==========================================================================
# nested anova
fit_Total_mass.tf    <- lme(Total_mass.tf ~ Function_group + Space_code + Function_group:Space_code, random = ~1|Species_code, data = data4gao, na.action = na.omit, method = "REML")
fit_Shoot_mass.tf    <- lme(Shoot_mass.tf ~ Function_group + Space_code + Function_group:Space_code, random = ~1|Species_code, data = data4gao, na.action = na.omit, method = "REML")
fit_Root_mass.tf     <- lme(Root_mass.tf ~ Function_group + Space_code + Function_group:Space_code, random = ~1|Species_code, data = data4gao, na.action = na.omit, method = "REML")
fit_Rmass_per_vol.tf <- lme(Rmass_per_vol.tf ~ Function_group + Space_code + Function_group:Space_code, random = ~1|Species_code, data = data4gao, na.action = na.omit, method = "REML")
fit_RM_ratio.tf      <- lme(RM_ratio.tf ~ Function_group + Space_code + Function_group:Space_code, random = ~1|Species_code, data = data4gao, na.action = na.omit, method = "REML")

fit_Total_mass.tf
fit_Shoot_mass.tf
fit_Root_mass.tf
fit_Rmass_per_vol.tf
fit_RM_ratio.tf

Anova(fit_Total_mass.tf, type = "III")
Anova(fit_Shoot_mass.tf, type = "III")
Anova(fit_Root_mass.tf, type = "III")
Anova(fit_Rmass_per_vol.tf, type = "III")
Anova(fit_RM_ratio.tf, type = "III")

# Output the cleaned results of ANOVAs
lme_aov.Total_mass    <- tidy(Anova(fit_Total_mass.tf, type = "III"))
lme_aov.Shoot_mass    <- tidy(Anova(fit_Shoot_mass.tf, type = "III"))
lme_aov.Root_mass     <- tidy(Anova(fit_Root_mass.tf, type = "III"))
lme_aov.Rmass_per_vol <- tidy(Anova(fit_Rmass_per_vol.tf, type = "III"))
lme_aov.RM_ratio      <- tidy(Anova(fit_RM_ratio.tf, type = "III"))
library(rlist)

#
list.lme_aovs <- list(
                     lme_aov.Total_mass = lme_aov.Total_mass,
                     lme_aov.Shoot_mass = lme_aov.Shoot_mass,
                     lme_aov.Root_mass = lme_aov.Root_mass,
                     lme_aov.Rmass_per_vol = lme_aov.Rmass_per_vol,
                     lme_aov.RM_ratio = lme_aov.RM_ratio
                     ) %>% list.rbind()

# save the cleaned results of ANOVAs
write.csv(list.lme_aovs, "list.lme_aovs.csv", row.names = TRUE)

# ==========================================================================
# Post-Hoc pairwise comparisons
# ==========================================================================
# loading the package
library(emmeans)

# post-hoc pairwise comparisons
emmeans(fit_Total_mass.tf, pairwise ~ Space_code)
emmeans(fit_Shoot_mass.tf, pairwise ~ Space_code)
emmeans(fit_Root_mass.tf, pairwise ~ Space_code)
emmeans(fit_Rmass_per_vol.tf, pairwise ~ Space_code)
emmeans(fit_RM_ratio.tf, pairwise ~ Space_code)

# ==========================================================================
#
# Species scale
# ANOVAs via lm function
#
# ==========================================================================
#
type3 <- list(Species_code = contr.sum, Space_code = contr.sum)

#nested anova
# fit.lm_Total_mass.tf    <- lm(Total_mass.tf ~ Species_code + Space_code + Species_code:Space_code, data = data4gao, na.action = na.omit)
fit.lm_Total_mass.tf    <- lm(Total_mass.tf ~ Species_code + Space_code + Species_code:Space_code, data = data4gao, na.action = na.omit, contrasts = type3)
fit.lm_Shoot_mass.tf    <- lm(Shoot_mass.tf ~ Species_code + Space_code + Species_code:Space_code, data = data4gao, na.action = na.omit, contrasts = type3)
fit.lm_Root_mass.tf     <- lm(Root_mass.tf ~ Species_code + Space_code + Species_code:Space_code, data = data4gao, na.action = na.omit, contrasts = type3)
fit.lm_Rmass_per_vol.tf <- lm(Rmass_per_vol.tf ~ Species_code + Space_code + Species_code:Space_code, data = data4gao, na.action = na.omit, contrasts = type3)
fit.lm_RM_ratio.tf      <- lm(RM_ratio.tf ~ Species_code + Space_code + Species_code:Space_code, data = data4gao, na.action = na.omit, contrasts = type3)

fit.lm_Total_mass.tf
fit.lm_Shoot_mass.tf
fit.lm_Root_mass.tf
fit.lm_Rmass_per_vol.tf
fit.lm_RM_ratio.tf

Anova(fit.lm_Total_mass.tf, type = "III")
Anova(fit.lm_Shoot_mass.tf, type = "III")
Anova(fit.lm_Root_mass.tf, type = "III")
Anova(fit.lm_Rmass_per_vol.tf, type = "III")
Anova(fit.lm_RM_ratio.tf, type = "III")

lm_aov.Total_mass    <- tidy(Anova(fit.lm_Total_mass.tf, type = "III"))
lm_aov.Shoot_mass    <- tidy(Anova(fit.lm_Shoot_mass.tf, type = "III"))
lm_aov.Root_mass     <- tidy(Anova(fit.lm_Root_mass.tf, type = "III"))
lm_aov.Rmass_per_vol <- tidy(Anova(fit.lm_Rmass_per_vol.tf, type = "III"))
lm_aov.RM_ratio      <- tidy(Anova(fit.lm_RM_ratio.tf, type = "III"))

library(rlist)

list.lm_aovs <- list(
                     lm_aov.Total_mass = lm_aov.Total_mass,
                     lm_aov.Shoot_mass = lm_aov.Shoot_mass,
                     lm_aov.Root_mass = lm_aov.Root_mass,
                     lm_aov.Rmass_per_vol = lm_aov.Rmass_per_vol,
                     lm_aov.RM_ratio = lm_aov.RM_ratio
                     ) %>% list.rbind()

write.csv(list.lm_aovs, "list.lm_aovs.csv", row.names = TRUE)

# ==========================================================================
# Post-Hoc pairwise comparisons
# ==========================================================================
library(emmeans)

# Species_code
emmeans(fit.lm_Total_mass.tf, pairwise ~ Species_code)
emmeans(fit.lm_Shoot_mass.tf, pairwise ~ Species_code)
emmeans(fit.lm_Root_mass.tf, pairwise ~ Species_code)
emmeans(fit.lm_Rmass_per_vol.tf, pairwise ~ Species_code)
emmeans(fit.lm_RM_ratio.tf, pairwise ~ Species_code)

# Space_code
emmeans(fit.lm_Total_mass.tf, pairwise ~ Space_code)
emmeans(fit.lm_Shoot_mass.tf, pairwise ~ Space_code)
emmeans(fit.lm_Root_mass.tf, pairwise ~ Space_code)
emmeans(fit.lm_Rmass_per_vol.tf, pairwise ~ Space_code)
emmeans(fit.lm_RM_ratio.tf, pairwise ~ Space_code)

# interaction
tukey_Total_mass.tf    <- emmeans(fit.lm_Total_mass.tf, pairwise ~ Space_code|Species_code)
tukey_Shoot_mass.tf    <- emmeans(fit.lm_Shoot_mass.tf, pairwise ~ Space_code|Species_code)
tukey_Root_mass.tf     <- emmeans(fit.lm_Root_mass.tf, pairwise ~ Space_code|Species_code)
tukey_Rmass_per_vol.tf <- emmeans(fit.lm_Rmass_per_vol.tf, pairwise ~ Space_code|Species_code)
tukey_RM_ratio.tf      <- emmeans(fit.lm_RM_ratio.tf, pairwise ~ Space_code|Species_code)

tukey_Total_mass.tf$contrasts
tukey_Shoot_mass.tf$contrasts
tukey_Root_mass.tf$contrasts
tukey_Rmass_per_vol.tf$contrasts
tukey_RM_ratio.tf$contrasts

#--------------------------------------------------------------------------#
#     .______.
#   __| _/\_ |__    ____
#  / __ |  | __ \ _/ ___\
# / /_/ |  | \_\ \\  \___
# \____ |  |___  / \___  >
#      \/      \/      \/
#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2022-03-16 10:13:57
# @Last Modified by:
# @Last Modified time: 2022-08-12 14:54
# @Description: Plots
#--------------------------------------------------------------------------#
# ==========================================================================
# Functional groups scale
# ==========================================================================
library(Rmisc)

# total mass
# --------------------------------------------------------------------------
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
dt_fun.Total_mass.tf <- summarySE(data4gao, measurevar = "Total_mass.tf", groupvars = c("Function_group", "Space_code"), na.rm = TRUE)
head(dt_fun.Total_mass.tf)

fg_plot.Total_mass.tf <- ggplot(dt_fun.Total_mass.tf, aes(x = Function_group, y = Total_mass.tf, fill = as.factor(Space_code))) +
                          geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Total_mass.tf - se, ymax = Total_mass.tf + se)) +
                          geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                                               scale_fill_manual(name = "Soil Space:",
                                                                  breaks = c("H15D15","H15D30","H60D15"),
                                                                  labels = c("D15H15","D30H15","D15H60"),
                                                                  values = c("white", "grey", "black")) +
                                               scale_x_discrete(breaks = c("grass", "forb", "legume"), labels = c("Grasses", "Forbs", "Legumes")) +
                                               theme(panel.background = element_rect(fill = NA),
                                                    legend.position=c(0.5, 0.85),
                                                    legend.text = element_text(family = "serif", colour = "black", size = 12),
                                                    legend.title = element_text(family = "serif", colour = "black", size = 12),
                                                    axis.line = element_line(size = 1, linetype = "solid"),
                                                    axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                                                    axis.text = element_text(family = "serif", colour = "black", size = 12),
                                                    axis.title = element_text(family = "serif", colour = "black", size = 12),
                                                    # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                                                    ) +
                                                labs(x = "Functional groups", y = "Log(total mass, g)") +
                                                scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1))

fg_plot.Total_mass.tf

# aboveground mass
# --------------------------------------------------------------------------
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

dt_fun.Shoot_mass.tf <- summarySE(data4gao, measurevar = "Shoot_mass.tf", groupvars = c("Function_group", "Space_code"), na.rm = TRUE)
head(dt_fun.Shoot_mass.tf)

fg_plot.Shoot_mass.tf <- ggplot(dt_fun.Shoot_mass.tf, aes(x = Function_group, y = Shoot_mass.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Shoot_mass.tf - se, ymax = Shoot_mass.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         scale_x_discrete(breaks = c("grass", "forb", "legume"), labels = c("Grasses", "Forbs", "Legumes")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 12),
                              axis.title = element_text(family = "serif", colour = "black", size = 12),
                              # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "Functional groups", y = "Log(shoot mass, g)") +
                          scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1))

fg_plot.Shoot_mass.tf


# root mass
# --------------------------------------------------------------------------
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

dt_fun.Root_mass.tf <- summarySE(data4gao, measurevar = "Root_mass.tf", groupvars = c("Function_group", "Space_code"), na.rm = TRUE)
head(dt_fun.Root_mass.tf)

fg_plot.Root_mass.tf <- ggplot(dt_fun.Root_mass.tf, aes(x = Function_group, y = Root_mass.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Root_mass.tf - se, ymax = Root_mass.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         scale_x_discrete(breaks = c("grass", "forb", "legume"), labels = c("Grasses", "Forbs", "Legumes")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 12),
                              axis.title = element_text(family = "serif", colour = "black", size = 12),
                              # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "Functional groups", y = "Log(root mass, g)") +
                          scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 1))

fg_plot.Root_mass.tf


# root density
# --------------------------------------------------------------------------
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

dt_fun.Rmass_per_vol.tf <- summarySE(data4gao, measurevar = "Rmass_per_vol.tf", groupvars = c("Function_group", "Space_code"), na.rm = TRUE)
head(dt_fun.Rmass_per_vol.tf)

fg_plot.Rmass_per_vol.tf <- ggplot(dt_fun.Rmass_per_vol.tf, aes(x = Function_group, y = Rmass_per_vol.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Rmass_per_vol.tf - se, ymax = Rmass_per_vol.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         scale_x_discrete(breaks = c("grass", "forb", "legume"), labels = c("Grasses", "Forbs", "Legumes")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 12),
                              axis.title = element_text(family = "serif", colour = "black", size = 12),
                              # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "Functional groups", y = "Log(root density, g/L)") +
                          scale_y_continuous(limits = c(-4, 1), breaks = seq(-4, 1, by = 1))

fg_plot.Rmass_per_vol.tf


# root mass fraction
# --------------------------------------------------------------------------
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

dt_fun.RM_ratio.tf <- summarySE(data4gao, measurevar = "RM_ratio.tf", groupvars = c("Function_group", "Space_code"), na.rm = TRUE)
head(dt_fun.RM_ratio.tf)

fg_plot.RM_ratio.tf <- ggplot(dt_fun.RM_ratio.tf, aes(x = Function_group, y = RM_ratio.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = RM_ratio.tf - se, ymax = RM_ratio.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         scale_x_discrete(breaks = c("grass", "forb", "legume"), labels = c("Grasses", "Forbs", "Legumes")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 12),
                              axis.title = element_text(family = "serif", colour = "black", size = 12),
                              # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "Functional groups", y = "Sqrt(root mass fraction)") +
                          scale_y_continuous(limits = c(0.3, 0.7), breaks = seq(0.3, 0.7, by = 0.1))

fg_plot.RM_ratio.tf

# --------------------------------------------------------------------------
# combine separate plots
# --------------------------------------------------------------------------
# ref:
# https://blog.csdn.net/weixin_42350411/article/details/113316500
library(patchwork)

# figure 1
fg_plots_p1 <- fg_plot.Total_mass.tf +
               fg_plot.Shoot_mass.tf +
               fg_plot.Root_mass.tf +
               fg_plot.Rmass_per_vol.tf +
               fg_plot.RM_ratio.tf +
               plot_layout(ncol = 2, byrow = FALSE) +
               plot_annotation(tag_levels = 'A')

(fg_plots_p1)

# --------------------------------------------------------------------------
# export the combined figure
# --------------------------------------------------------------------------
library(ggpubr)
ggexport(fg_plots_p1, filename = "./20220508fg_plots.png",
         width = 1800,
         height = 2700,
         pointsize = 12,
         res = 300)

# ==========================================================================
# species scale
# ==========================================================================
# loading packages
library(Rmisc)
library(dplyr)
library(grid)

# data summarization
data4gao.summary <- data4gao %>% group_by(Function_group, Species_code) %>% dplyr::summarise(n = sum(!is.na(Total_mass)))
data4gao.summary

levels.order <- data4gao.summary$Species_code

# total mass
# --------------------------------------------------------------------------
# data summarization
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
dt_sp.Total_mass.tf <- summarySE(data4gao, measurevar = "Total_mass.tf", groupvars = c("Species_code", "Space_code"), na.rm = TRUE)
head(dt_sp.Total_mass.tf)

# adjust the order of Species_code levels
dt_sp.Total_mass.tf$Species_code <- factor(dt_sp.Total_mass.tf$Species_code, order = TRUE, levels = levels.order)
levels(dt_sp.Total_mass.tf$Species_code)

# plot
# add text below the x-axis, ref:
# https://stackoverflow.com/questions/61549243/secondary-x-axis-labels-for-sample-size-with-ggplot2-on-r
sp_plot.Total_mass.tf <- ggplot(dt_sp.Total_mass.tf, aes(x = Species_code, y = Total_mass.tf, fill = as.factor(Space_code))) +
                          geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Total_mass.tf - se, ymax = Total_mass.tf + se)) +
                          geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                                               scale_fill_manual(name = "Soil Space",
                                                                  breaks = c("H15D15","H15D30","H60D15"),
                                                                  labels = c("D15H15","D30H15","D15H60"),
                                                                  values = c("white", "grey", "black")) +
                                               theme(panel.background = element_rect(fill = NA),
                                                    legend.position=c(0.2, 0.85),
                                                    legend.text = element_text(family = "serif", colour = "black", size = 14),
                                                    legend.title = element_text(family = "serif", colour = "black", size = 14),
                                                    axis.line = element_line(size = 1, linetype = "solid"),
                                                    axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                                                    axis.text = element_text(family = "serif", colour = "black", size = 14),
                                                    axis.title = element_text(family = "serif", colour = "black", size = 14),
                                                    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                                                    ) +
                                                labs(x = "", y = "Log(total mass, g)") +
                                                scale_y_continuous(limits = c(-2, 8), breaks = seq(-2, 8, by = 2))

# aboveground mass
# --------------------------------------------------------------------------
# data summarization
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
dt_sp.Shoot_mass.tf <- summarySE(data4gao, measurevar = "Shoot_mass.tf", groupvars = c("Species_code", "Space_code"), na.rm = TRUE)
head(dt_sp.Shoot_mass.tf)

# adjust the order of Species_code levels
dt_sp.Shoot_mass.tf$Species_code <- factor(dt_sp.Shoot_mass.tf$Species_code, order = TRUE, levels = levels.order)
levels(dt_sp.Shoot_mass.tf$Species_code)

# plot
sp_plot.Shoot_mass.tf <- ggplot(dt_sp.Shoot_mass.tf, aes(x = Species_code, y = Shoot_mass.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Shoot_mass.tf - se, ymax = Shoot_mass.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 14),
                              axis.title = element_text(family = "serif", colour = "black", size = 14),
                              axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "", y = "Log(shoot mass, g)") +
                          scale_y_continuous(limits = c(-2, 8), breaks = seq(-2, 8, by = 2))

sp_plot.Shoot_mass.tf

# root mass
# --------------------------------------------------------------------------
# data summarization
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
dt_sp.Root_mass.tf <- summarySE(data4gao, measurevar = "Root_mass.tf", groupvars = c("Species_code", "Space_code"), na.rm = TRUE)
head(dt_sp.Root_mass.tf)

# adjust the order of Species_code levels
dt_sp.Root_mass.tf$Species_code <- factor(dt_sp.Root_mass.tf$Species_code, order = TRUE, levels = levels.order)
levels(dt_sp.Root_mass.tf$Species_code)

# plot
sp_plot.Root_mass.tf <- ggplot(dt_sp.Root_mass.tf, aes(x = Species_code, y = Root_mass.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Root_mass.tf - se, ymax = Root_mass.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 14),
                              axis.title = element_text(family = "serif", colour = "black", size = 14),
                              axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "", y = "Log(root mass, g)") +
                          scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2))


# root density
# --------------------------------------------------------------------------
# data summarization
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

dt_sp.Rmass_per_vol.tf <- summarySE(data4gao, measurevar = "Rmass_per_vol.tf", groupvars = c("Species_code", "Space_code"), na.rm = TRUE)
head(dt_sp.Rmass_per_vol.tf)

# adjust the order of Species_code levels
dt_sp.Rmass_per_vol.tf$Species_code <- factor(dt_sp.Rmass_per_vol.tf$Species_code, order = TRUE, levels = levels.order)
levels(dt_sp.Rmass_per_vol.tf$Species_code)

# plot
sp_plot.Rmass_per_vol.tf <- ggplot(dt_sp.Rmass_per_vol.tf, aes(x = Species_code, y = Rmass_per_vol.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = Rmass_per_vol.tf - se, ymax = Rmass_per_vol.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(name = "Soil Space",
                                            breaks = c("H15D15","H15D30","H60D15"),
                                            labels = c("D15H15","D30H15","D15H60"),
                                            values = c("white", "grey", "black")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position=c(0.2, 0.85),
                              legend.text = element_text(family = "serif", colour = "black", size = 14),
                              legend.title = element_text(family = "serif", colour = "black", size = 14),
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 14),
                              axis.title = element_text(family = "serif", colour = "black", size = 14),
                              axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "", y = "Log(root density, g/L)") +
                          scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2))

sp_plot.Rmass_per_vol.tf

# root mass fraction
# --------------------------------------------------------------------------
# data summarization
# summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

dt_sp.RM_ratio.tf <- summarySE(data4gao, measurevar = "RM_ratio.tf", groupvars = c("Species_code", "Space_code"), na.rm = TRUE)
head(dt_sp.RM_ratio.tf)

# adjust the order of Species_code levels
dt_sp.RM_ratio.tf$Species_code <- factor(dt_sp.RM_ratio.tf$Species_code, order = TRUE, levels = levels.order)
levels(dt_sp.RM_ratio.tf$Species_code)

sp_plot.RM_ratio.tf <- ggplot(dt_sp.RM_ratio.tf, aes(x = Species_code, y = RM_ratio.tf, fill = as.factor(Space_code))) +
    geom_errorbar(position = position_dodge(.8), width = .2, size = 0.5, aes(ymin = RM_ratio.tf - se, ymax = RM_ratio.tf + se)) +
    geom_point(position = position_dodge(.8), shape = 21, size = 4) +
                         scale_fill_manual(values = c("white", "grey", "black")) +
                         theme(panel.background = element_rect(fill = NA),
                              legend.position = "none",
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 14),
                              axis.title = element_text(family = "serif", colour = "black", size = 14),
                              axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
                              ) +
                          labs(x = "", y = "Sqrt(root mass fraction)") +
                          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

sp_plot.RM_ratio.tf

# ==========================================================================
# combine separate plots
# ==========================================================================
# add reference line
# --------------------------------------------------------------------------
sp_plot.Root_mass.tf01 <- sp_plot.Root_mass.tf + coord_cartesian(clip = "off") +
annotation_custom(linesGrob(x = unit(c(0.03, 0.29), 'npc'), y = unit(c(- 0.345, - 0.345), 'npc'), gp = gpar(lwd = 2))) +
annotation_custom(linesGrob(x = unit(c(0.34, 0.66), 'npc'), y = unit(c(- 0.345, - 0.345), 'npc'), gp = gpar(lwd = 2))) +
annotation_custom(linesGrob(x = unit(c(0.71, 0.97), 'npc'), y = unit(c(- 0.345, - 0.345), 'npc'), gp = gpar(lwd = 2))) +
annotation_custom(textGrob(label = "Forbs", gp = gpar(fontfamily = "serif", fontsize = 14)), xmin = 3.5, xmax = 3.5, ymin = -8, ymax = -8) +
annotation_custom(textGrob(label = "Grasses", gp = gpar(fontfamily = "serif", fontsize = 14)), xmin = 10, xmax = 10, ymin = -8, ymax = -8) +
annotation_custom(textGrob(label = "Legumes", gp = gpar(fontfamily = "serif", fontsize = 14)), xmin = 16.5, xmax = 16.5, ymin = -8, ymax = -8)

sp_plot.Root_mass.tf01


sp_plot.RM_ratio.tf01 <- sp_plot.RM_ratio.tf + coord_cartesian(clip = "off") +
annotation_custom(linesGrob(x = unit(c(0.03, 0.29), 'npc'), y = unit(c(- 0.345, - 0.345), 'npc'), gp = gpar(lwd = 2))) +
annotation_custom(linesGrob(x = unit(c(0.34, 0.66), 'npc'), y = unit(c(- 0.345, - 0.345), 'npc'), gp = gpar(lwd = 2))) +
annotation_custom(linesGrob(x = unit(c(0.71, 0.97), 'npc'), y = unit(c(- 0.345, - 0.345), 'npc'), gp = gpar(lwd = 2))) +
annotation_custom(textGrob(label = "Forbs", gp = gpar(fontfamily = "serif", fontsize = 14)), xmin = 3.5, xmax = 3.5, ymin = -0.5, ymax = -0.5) +
annotation_custom(textGrob(label = "Grasses", gp = gpar(fontfamily = "serif", fontsize = 14)), xmin = 10, xmax = 10, ymin = -0.5, ymax = -0.5) +
annotation_custom(textGrob(label = "Legumes", gp = gpar(fontfamily = "serif", fontsize = 14)), xmin = 16.5, xmax = 16.5, ymin = -0.5, ymax = -0.5)

sp_plot.RM_ratio.tf01

# add p values
# --------------------------------------------------------------------------
# using levels.order, determine the position of levels in the x-axis
df.levels.order <- data.frame(Species_code = as.character(levels.order), levels.order = seq(levels.order))
df.levels.order

# using emmeans pacakge to check the p-values
df.tukey_Total_mass.tf    <- as.data.frame(tukey_Total_mass.tf$contrasts)
df.tukey_Shoot_mass.tf    <- as.data.frame(tukey_Shoot_mass.tf$contrasts)
df.tukey_Root_mass.tf     <- as.data.frame(tukey_Root_mass.tf$contrasts)
df.tukey_Rmass_per_vol.tf <- as.data.frame(tukey_Rmass_per_vol.tf$contrasts)
df.tukey_RM_ratio.tf      <- as.data.frame(tukey_RM_ratio.tf$contrasts)

df.tukey_Total_mass.tf
df.tukey_Shoot_mass.tf
df.tukey_Root_mass.tf
df.tukey_Rmass_per_vol.tf
df.tukey_RM_ratio.tf

#
# total mass
# --------------------------------------------------------------------------
# remain the data with p < 0.05, and the corresponding species name
tukey_Total_mass.species_code <- df.tukey_Total_mass.tf %>% filter(p.value <= 0.05) %>% dplyr::select(Species_code) %>% unique()
tukey_Total_mass.species_code$Species_code <- as.character(tukey_Total_mass.species_code$Species_code)

# combine df.levels.order and tukey_Total_mass.species_code
tukey_Total_mass.species_code01 <- tukey_Total_mass.species_code %>% dplyr::left_join(df.levels.order, by = "Species_code")
tukey_Total_mass.species_code01

# combine df.levels.order01 and dt_sp.Total_mass.tf
dt_sp.Total_mass.tf$Species_code <- as.character(dt_sp.Total_mass.tf$Species_code)
tukey_Total_mass.species_code02 <- tukey_Total_mass.species_code01 %>% dplyr::left_join(dt_sp.Total_mass.tf, by = "Species_code")
tukey_Total_mass.species_code02

# A new column is added, marked by letters, to show the significance
# df.tukey_Total_mass.tf %>% filter(p.value <=0.05)
tukey_Total_mass.species_code03 <- tukey_Total_mass.species_code02 %>%
                        dplyr::select(Species_code, Space_code, levels.order, Total_mass.tf) %>%
                        mutate(p_lab = c("a", "a", "b",
                                         "a", "b", "b",
                                         "a", "b", "c",
                                         "a", "ab", "b",
                                         "a", "ab", "b",
                                         "a", "ab", "b",
                                         "a", "ab", "b",
                                         "a", "a", "b"))

tukey_Total_mass.species_code03

# add significance letter in sp_plot.Total_mass.tf
sp_plot.Total_mass.tf01 <- sp_plot.Total_mass.tf +
                           geom_text(data = tukey_Total_mass.species_code03, aes(x = levels.order, y = Total_mass.tf + 1.4, label = p_lab), size = 5, position = position_dodge(.8))
sp_plot.Total_mass.tf01

# aboveground mass
# --------------------------------------------------------------------------
# remain the data with p < 0.05, and the corresponding species name
tukey_Shoot_mass.species_code <- df.tukey_Shoot_mass.tf %>% filter(p.value <= 0.05) %>% dplyr::select(Species_code) %>% unique()
tukey_Shoot_mass.species_code$Species_code <- as.character(tukey_Shoot_mass.species_code$Species_code)

# combine df.levels.order and tukey_Shoot_mass.species_code
tukey_Shoot_mass.species_code01 <- tukey_Shoot_mass.species_code %>% dplyr::left_join(df.levels.order, by = "Species_code")
tukey_Shoot_mass.species_code01

# combine df.levels.order01 and dt_sp.Shoot_mass.tf
dt_sp.Shoot_mass.tf$Species_code <- as.character(dt_sp.Shoot_mass.tf$Species_code)
tukey_Shoot_mass.species_code02 <- tukey_Shoot_mass.species_code01 %>% dplyr::left_join(dt_sp.Shoot_mass.tf, by = "Species_code")
tukey_Shoot_mass.species_code02

# A new column is added, marked by letters, to show the significance
# df.tukey_Total_mass.tf %>% filter(p.value <=0.05)
tukey_Shoot_mass.species_code03 <- tukey_Shoot_mass.species_code02 %>%
                        dplyr::select(Species_code, Space_code, levels.order, Shoot_mass.tf) %>%
                        mutate(p_lab = c("a", "a", "b",
                                         "a", "b", "b",
                                         "a", "b", "c",
                                         "a", "ab", "b",
                                         "a", "b", "b",
                                         "a", "ab", "b",
                                         "a", "ab", "b",
                                         "a", "b", "a",
                                         "a", "ab", "b",
                                         "a", "a", "b"))

tukey_Shoot_mass.species_code03

# add significance letter in sp_plot.Shoot_mass.tf
sp_plot.Shoot_mass.tf01 <- sp_plot.Shoot_mass.tf +
                           geom_text(data = tukey_Shoot_mass.species_code03, aes(x = levels.order, y = Shoot_mass.tf + 1.4, label = p_lab), size = 5, position = position_dodge(.8))
sp_plot.Shoot_mass.tf01

# --------------------------------------------------------------------------
# combine separate plots
# --------------------------------------------------------------------------
# loading packages
library(patchwork)

# figure 2
sp_plots_p1 <- sp_plot.Total_mass.tf01 +
               sp_plot.Shoot_mass.tf01 +
               sp_plot.Root_mass.tf01 +
               plot_layout(ncol = 1, byrow = FALSE) +
               plot_annotation(tag_levels = 'A')

# figure 3
sp_plots_p2 <- sp_plot.Rmass_per_vol.tf +
               sp_plot.RM_ratio.tf01 +
               plot_layout(ncol = 1, byrow = FALSE) +
               plot_annotation(tag_levels = 'A')

(sp_plots_p1)
(sp_plots_p2)

# --------------------------------------------------------------------------
# export the combined figure
# --------------------------------------------------------------------------
library(ggpubr)
ggexport(sp_plots_p1, filename = "./20220508_sp_plots01.png",
         width = 3500,
         height = 3000,
         pointsize = 12,
         res = 300)

ggexport(sp_plots_p2, filename = "./20220508_sp_plots02.png",
         width = 3500,
         height = 2000,
         pointsize = 12,
         res = 300)


