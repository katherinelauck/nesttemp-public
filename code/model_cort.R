### Create all objects needed for cort analysis

### make prez figures
require(tidyverse)
require(ggplot2)
require(lubridate)
require(scales)
require(viridis)
require(ggeffects)
require(ggimage)
require(magick)
require(cowplot)
require(ggggeffects)
require(emmeans)
require(multcomp)
require(gt)
library(glmmTMB)
library(dagitty)
library(ggdag)
library(ggplot2)
library(igraph)
library(tidyverse)
library(performance)
library(modelsummary)
library(car)
library(lme4)
library(hms)
library(weathermetrics)
library(future)
# modelsummary::get_gof(my_model)
# performance::model_performance(my_model)
# wmean <- read_rds("data/wmean.rds")  # unused

g <- read_rds("data/growth_cort_provis_manytempmeasures.rds") %>%
  mutate(year_fct = as.factor(year))

## Objectives

# 1. New responses:
#   - s2 - s1: absolute change
# - (s2-s1)/s1: proportional change
# 1. Do responses satisfy normality?
#   1. Produce one model with max and one model with min, look at graphically to see what looks nicest
# 1. Evaluate model:
#   - fit
# - qqplots
# - outliers
# - residuals vs temp
#
# ## Questions for lab meeting
#
# 1. If we think that age of nestling mediates their ability to respond to stress (we do), should I test for a 3-way interaction between age, temp, and habitat? Does this make more sense as a formal mediation analysis?
#   1. I like the idea of parsing out the effect of max temp from min temp, could be a useful comparison to make with previous work looking at the cort effects of cold, and matches our pre-existing expectation that cold and hot operate differently. Based on your understanding of bird physiology, what are the squared temp vs the cold and hot temp models telling us?
#   1. Reflect on finding that abs change has more of an interactive signal than baseline cort

## Do responses satisfy normality?

## WEBL

### Baseline cort


# hist(g %>% filter(Species == "WEBL") %>% pull(cort_s1))
# qqnorm(g %>% filter(Species == "WEBL") %>% pull(cort_s1))
# qqline(g %>% filter(Species == "WEBL") %>% pull(cort_s1))
#
#
# May be basically okay but a transformation may help.
#
# ### Square root transformation
#
#
# hist(sqrt(g %>% filter(Species == "WEBL") %>% pull(cort_s1)))
# qqnorm(sqrt(g %>% filter(Species == "WEBL") %>% pull(cort_s1)))
# qqline(sqrt(g %>% filter(Species == "WEBL") %>% pull(cort_s1)))
#
#
# This looks like enough for cort_s1.
#
# ### Absolute difference cort_s2-cort_s1
#
#
# hist(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort))
# qqnorm(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort))
# qqline(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort))
#
#
# Not great! Needs to be transformed!
#
#   ### Square root transformation
#
#
# hist(sqrt(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort)))
# qqnorm(sqrt(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort)))
# qqline(sqrt(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort)))
#
#
# Closer... Maybe try natural log?
#
#   ### Natural log
#
#
# hist(log(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort)))
# qqnorm(log(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort)))
# qqline(log(g %>% filter(Species == "WEBL") %>% pull(abs_change_cort)))
#
# Too far! Looks like sqrt is the best transformation.
#
# ### Proportional difference (cort_s2-cort_s1)/cort_s1
#
#
# hist(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort))
# qqnorm(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort))
# qqline(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort))
#
#
# Needs to be transformed.
#
# ### Square root transformation
#
#
# hist(sqrt(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort)))
# qqnorm(sqrt(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort)))
# qqline(sqrt(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort)))
#
#
# Not quite enough... maybe the natural log? But this does look better.
#
# ### Natural log
#
#
# hist(log(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort)))
# qqnorm(log(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort)))
# qqline(log(g %>% filter(Species == "WEBL") %>% pull(prop_change_cort)))
#
#
# Too far, looks like the square root transformation is the one for all three responses for WEBL.
#
# ### Conclusions
#
# Square root transformation for all three cort response variables helps with normality of response variable.
#
# ## TRES
#
# ### Baseline cort
#
#
# hist(g %>% filter(Species == "TRES") %>% pull(cort_s1))
# qqnorm(g %>% filter(Species == "TRES") %>% pull(cort_s1))
# qqline(g %>% filter(Species == "TRES") %>% pull(cort_s1))
#
#
# May be basically okay but a transformation may help.
#
# ### Square root transformation
#
#
# hist(sqrt(g %>% filter(Species == "TRES") %>% pull(cort_s1)))
# qqnorm(sqrt(g %>% filter(Species == "TRES") %>% pull(cort_s1)))
# qqline(sqrt(g %>% filter(Species == "TRES") %>% pull(cort_s1)))
#
#
# This looks like enough for cort_s1.
#
# ### Absolute difference cort_s2-cort_s1
#
#
# hist(g %>% filter(Species == "TRES") %>% pull(abs_change_cort))
# qqnorm(g %>% filter(Species == "TRES") %>% pull(abs_change_cort))
# qqline(g %>% filter(Species == "TRES") %>% pull(abs_change_cort))
#
#
# Honestly looks mostly okay to me.
#
# ### Square root transformation
#
#
# hist(sqrt(g %>% filter(Species == "TRES") %>% pull(abs_change_cort)))
# qqnorm(sqrt(g %>% filter(Species == "TRES") %>% pull(abs_change_cort)))
# qqline(sqrt(g %>% filter(Species == "TRES") %>% pull(abs_change_cort)))
#
#
# That's enough for abs change!
#
# ### Proportional difference (cort_s2-cort_s1)/cort_s1
#
#
# hist(g %>% filter(Species == "TRES") %>% pull(prop_change_cort))
# qqnorm(g %>% filter(Species == "TRES") %>% pull(prop_change_cort))
# qqline(g %>% filter(Species == "TRES") %>% pull(prop_change_cort))
#
#
# Needs to be transformed.
#
# ### Square root transformation
#
#
# hist(sqrt(g %>% filter(Species == "TRES") %>% pull(prop_change_cort)))
# qqnorm(sqrt(g %>% filter(Species == "TRES") %>% pull(prop_change_cort)))
# qqline(sqrt(g %>% filter(Species == "TRES") %>% pull(prop_change_cort)))
#
#
# Honestly looks pretty much good. Will try the natural log just in case.
#
# ### Natural log
#
#
# hist(log(g %>% filter(Species == "TRES") %>% pull(prop_change_cort)))
# qqnorm(log(g %>% filter(Species == "TRES") %>% pull(prop_change_cort)))
# qqline(log(g %>% filter(Species == "TRES") %>% pull(prop_change_cort)))
#
#
# Too far, looks like the square root transformation is the one for all three responses for TRES.
#
# ### Conclusions
#
# Square root transformation for all three cort response variables helps with normality of response variable.
#
# ## For each response, produce 2 models:
#
# 1. linear max and min
# 1. Squared max
#
# ## WEBL
#
# ### sqrt(cort_s1)
#
#
# s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# s1_sqmax <- lmerTest::lmer(sqrt(cort_s1) ~ poly(meanmaxtempI_scaled,2) * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                   ~ scale(.x)[,1],
#                                                   .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# #### Plots
#
#
# (pl <- ggpredict(s1_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#    plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#    facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Blood cort concentration") +
#    scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max and min") +
#    # ylim(0,100) +
#    annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
#    )
#
#
#
# summary(s1_lintemp)
#
#
#
#
# (pl <- ggpredict(s1_sqmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#    plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#    facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Blood cort concentration") +
#    scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; sq max") +
#    # ylim(0,100) +
#    annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
#    )
#
#
#
# summary(s1_sqmax)
#
#
#
#
#
# performance::check_predictions(s1_lintemp)
# performance::check_predictions(s1_sqmax)
#
#
#
# performance::check_collinearity(s1_lintemp)
#
#
#
# performance::check_collinearity(s1_sqmax)
#
#
# High VIFs all seem to come from interactions so probably not a huge issue. But maybe the second model is the best one because there is a moderate VIF for max and min temp.
#
# #### What happens to VIFs for linear model without interactions?
#
# This should tell us whether max and min are worryingly collinear with each other.
#
#
# s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + (1|year) + (1|site/attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
# check_collinearity(s1_lintemp)
#
#
# It looks like they aren't, which is what I thought. So we should be good there.
#
# ### sqrt(abs_change_cort)
#
#
# abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                 mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                               ~ scale(.x)[,1],
#                                               .names = "{.col}_scaled"),
#                                        meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# abs_sqmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ poly(meanmaxtempI_scaled,2) * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                             ~ scale(.x)[,1],
#                                             .names = "{.col}_scaled"),
#                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# #### Plots
#
#
# (pl <- ggpredict(abs_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Absolute increase S1 -> S2") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max and min") +
#     # ylim(0,100) +
#     annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
# )
#
#
# <!--   -->
#   <!-- (pl <- ggpredict(abs_sqmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>% -->
#           <!--    plot(line_size = 1.5,alpha = .2,show_data = TRUE) + -->
#           <!--     theme_classic() + -->
#           <!--    facet_wrap(~ group, ncol = 2) + -->
#           <!--     xlab("Mean daily max temp over preceding week") + -->
#           <!--     ylab("Absolute increase S1 -> S2") + -->
#           <!--    scale_fill_viridis(discrete = TRUE) + -->
#           <!--     scale_color_viridis(discrete = TRUE) + -->
#           <!--     theme(text = element_text(size = 16)) + -->
#           <!--     labs(title = "WEBL; sq max") + -->
#           <!--    # ylim(0,100) + -->
#           <!--    annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp))) -->
#           <!--    ) -->
#   <!--   -->
#
#
#   <!--   -->
#   <!-- performance::check_predictions(abs_lintemp) -->
#   <!--   -->
#
#
#   <!--   -->
#   <!-- performance::check_predictions(abs_sqmax) -->
#   <!--   -->
#
#
#
# performance::check_collinearity(abs_lintemp)
#
#
#
# performance::check_collinearity(abs_sqmax)
#
#
# This response seems a little bit rough in the model dept. Tough to tell what exactly is going on - got several errors when graphing and checking performance.
#
# #### sqrt(prop_change_cort)
#
#
# prop_lintemp <- lmerTest::lmer(sqrt(prop_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                ~ scale(.x)[,1],
#                                                .names = "{.col}_scaled"),
#                                         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# prop_sqmax <- lmerTest::lmer(sqrt(prop_change_cort) ~ poly(meanmaxtempI_scaled,2) * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                              ~ scale(.x)[,1],
#                                              .names = "{.col}_scaled"),
#                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# #### Plots
#
#
# (pl <- ggpredict(prop_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Absolute increase S1 -> S2") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max and min") +
#     # ylim(0,100) +
#     annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
# )
#
#
# <!--   -->
#   <!-- (pl <- ggpredict(prop_sqmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>% -->
#           <!--    plot(line_size = 1.5,alpha = .2,show_data = TRUE) + -->
#           <!--     theme_classic() + -->
#           <!--    facet_wrap(~ group, ncol = 2) + -->
#           <!--     xlab("Mean daily max temp over preceding week") + -->
#           <!--     ylab("Absolute increase S1 -> S2") + -->
#           <!--    scale_fill_viridis(discrete = TRUE) + -->
#           <!--     scale_color_viridis(discrete = TRUE) + -->
#           <!--     theme(text = element_text(size = 16)) + -->
#           <!--     labs(title = "WEBL; sq max") + -->
#           <!--    # ylim(0,100) + -->
#           <!--    annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp))) -->
#           <!--    ) -->
#   <!--   -->
#
#
#   <!--   -->
#   <!-- performance::check_predictions(abs_lintemp) -->
#   <!--   -->
#
#
#   <!--   -->
#   <!-- performance::check_predictions(abs_sqmax) -->
#   <!--   -->
#
#
#
# performance::check_collinearity(abs_lintemp)
#
#
#
# performance::check_collinearity(abs_sqmax)
#
#
# This response seems a little bit rough in the model dept. Tough to tell what exactly is going on because I got so many errors.
#
# ### Conclusions:
#
# I like sqrt(cort_s1) as a response better than the others. There seems to be more signal there, and it's the most normal response.
#
# Revisiting those two models, I think an argument could be made for either. The VIFs are inflated by the interaction terms for both models, and when I run the linear model without the max * habitat interaction, its VIF is under the threshold.
#
# I like the idea of parsing out the effect of max temp from min temp, could be a useful comparison to make with previous work looking at the cort effects of cold. This is something I'd like to bring up in lab meeting - maybe crowdsource some thoughts about what those two models tell us.
#
# Looks like row crop is significantly different from forest with linear temp, but not with squared temp. Another less rigorous argument for linear temp. It is also so much easier to interpret the coefficients with linear temp!!!
#
#
#
#   ## TRES
#
#   ### sqrt(cort_s1)
#
#
# s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                              ~ scale(.x)[,1],
#                                              .names = "{.col}_scaled"),
#                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# s1_sqmax <- lmerTest::lmer(sqrt(cort_s1) ~ poly(meanmaxtempI_scaled,2) * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                              mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                            ~ scale(.x)[,1],
#                                            .names = "{.col}_scaled"),
#                                     meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# #### Plots
#
#
# (pl <- ggpredict(s1_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Blood cort concentration") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max and min") +
#     # ylim(0,100) +
#     annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
# )
#
#
#
# summary(s1_lintemp)
#
#
#
#
# (pl <- ggpredict(s1_sqmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Blood cort concentration") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; sq max") +
#     ylim(0,100) +
#     annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
# )
#
#
#
# summary(s1_sqmax)
#
#
#
#
#
# performance::check_predictions(s1_lintemp)
# performance::check_predictions(s1_sqmax)
#
#
#
# performance::check_collinearity(s1_lintemp)
#
#
#
# performance::check_collinearity(s1_sqmax)
#
#
# High VIFs all seem to come from interactions so probably not a huge issue. But maybe the second model is the best one because there is a moderate VIF for max and min temp.
#
# #### What happens to VIFs for without interactions?
#
# This should tell us whether max and min are worryingly collinear with each other.
#
#
# s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + (1|year) + (1|site/attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                              ~ scale(.x)[,1],
#                                              .names = "{.col}_scaled"),
#                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
# check_collinearity(s1_lintemp)
#
#
# It looks like they aren't, which is what I thought. So we should be good there.
#
#
# s1_sqmax <- lmerTest::lmer(sqrt(cort_s1) ~ poly(meanmaxtempI_scaled,2) + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                   ~ scale(.x)[,1],
#                                                   .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
# check_collinearity(s1_sqmax)
#
#
# ### sqrt(abs_change_cort)
#
#
# abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# abs_sqmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ poly(meanmaxtempI_scaled,2) * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                   ~ scale(.x)[,1],
#                                                   .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# #### Plots
#
#
# (pl <- ggpredict(abs_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#    plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#    facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Absolute increase S1 -> S2") +
#    scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max and min") +
#    ylim(0,100) +
#    annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
#    )
#
#
#
# (pl <- ggpredict(abs_sqmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#    plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#    facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Absolute increase S1 -> S2") +
#    scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; sq max") +
#    ylim(0,120) +
#    annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
#    )
#
#
# I feel like there's just not enough data in forest and orchard to model those. I wonder if filtering to just grassland and row crop would help here.
#
#
#
# performance::check_predictions(abs_lintemp)
#
#
#
#
# performance::check_predictions(abs_sqmax)
#
#
#
#
# performance::check_collinearity(abs_lintemp)
#
#
#
# abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + habitat + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                 mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                               ~ scale(.x)[,1],
#                                               .names = "{.col}_scaled"),
#                                        meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# performance::check_collinearity(abs_lintemp)
#
#
#
# performance::check_collinearity(abs_sqmax)
#
#
#
#
# abs_sqmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ poly(meanmaxtempI_scaled,2) + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                             ~ scale(.x)[,1],
#                                             .names = "{.col}_scaled"),
#                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# check_collinearity(abs_sqmax)
#
#
#
# #### sqrt(prop_change_cort)
#
#
# prop_lintemp <- lmerTest::lmer(sqrt(prop_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                                ~ scale(.x)[,1],
#                                                .names = "{.col}_scaled"),
#                                         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# prop_sqmax <- lmerTest::lmer(sqrt(prop_change_cort) ~ poly(meanmaxtempI_scaled,2) * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
#                                              ~ scale(.x)[,1],
#                                              .names = "{.col}_scaled"),
#                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
# #### Plots
#
#
# (pl <- ggpredict(prop_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Proportional increase S1 -> S2") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max and min") +
#     # ylim(0,100) +
#     annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
# )
#
#
#
# (pl <- ggpredict(prop_sqmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week") +
#     ylab("Proportional increase S1 -> S2") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; sq max") +
#     # ylim(0,100) +
#     annotate("text", x = -1.5, y = 15, size = 5, label = paste("N =",nobs(s1_lintemp)))
# )
#
#
#
#
# performance::check_predictions(abs_lintemp)
#
#
#
#
# performance::check_predictions(abs_sqmax)
#
#
#
#
# performance::check_collinearity(abs_lintemp)
#
#
#
# performance::check_collinearity(abs_sqmax)
#
#
# This response seems a little bit rough in the model dept. Tough to tell what exactly is going on, but again I do like the simplicity of the linear model.
#
# ### Conclusions:
#
# For TRES, I like both cort_s1 and absolute increase. prop increase has a singular fit and looks sketchier when graphed. Could consider comparing just grassland and row crop for TRES cort as forest and orchard have very few data points.
#
# With how little data TRES cort has, I think I'd argue for the simplicity of max and min temp. The VIFs indicate including both is not a problem.
#
# I like the idea of parsing out the effect of max temp from min temp, could be a useful comparison to make with previous work looking at the cort effects of cold, and matches our pre-existing expectation that cold and hot operate differently. This is something I'd like to bring up in lab meeting - maybe crowdsource some thoughts about what those two models tell us.
#
# There might not be much going on (or enough data to identify trends) for TRES cort and temp.
#
# ## Now that we have determined the temp structure for cort, does the interaction with habitat help the model?

# Will check interaction with max and min.

### WEBL meanmaxtempI

#### s1_cort
s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_webl,"figures/int_tab_s1_webl.html")

s1_webl <- s1_lintemp_noint

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_webl,"figures/s1byhabitat_webl.html")


## trends
##
(webls1trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1trendmax,"figures/webls1trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


# ss_year_webl_s1 <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
# ss_year_webl_s1 %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
#   gtsave("figures/ss_year_webl_s1.html")


ss_year_webl_s1 <- s1_lintemp_noint@frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total = `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_webl_s1 <- ss_year_webl_s1 %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
    })
)


t_ss_year_webl_s1


samp_s1_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_webl %>% gt()


dat_text_s1_webl <- data.frame(
  label = paste("N =",samp_s1_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_webl = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_s1_webl <- mean(data_s1_webl %>% pull(meanmaxtempI))
sd_temp_s1_webl <- sd(data_s1_webl %>% pull(meanmaxtempI))


temp_trans_s1_webl <- trans_new("temp_trans_s1_webl",
                                transform = function(x){(x * sd_temp_s1_webl) + mean_temp_s1_webl},
                                inverse = function(x){x})

(figs2_webl <- ggpredict(s1_lintemp_noint,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_webl,
                       breaks = c((20-mean_temp_s1_webl)/sd_temp_s1_webl,
                                  (25-mean_temp_s1_webl)/sd_temp_s1_webl,
                                  (30-mean_temp_s1_webl)/sd_temp_s1_webl,
                                  (35-mean_temp_s1_webl)/sd_temp_s1_webl,
                                  (40-mean_temp_s1_webl)/sd_temp_s1_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bymaxtempxhab_WEBL.png",plot =  figs2_webl, width = 10, height = 6.6)




data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp <- mean(data %>% pull(meanmintempI))
sd_temp <- sd(data %>% pull(meanmintempI))


temp_trans <- trans_new("temp_trans",
                        transform = function(x){(x * sd_temp) + mean_temp},
                        inverse = function(x){x})

(pl <- ggpredict(s1_lintemp_noint,terms = c("meanmintempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily min temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans,
                       breaks = c((8-mean_temp)/sd_temp,
                                  (10-mean_temp)/sd_temp,
                                  (12-mean_temp)/sd_temp,
                                  (14-mean_temp)/sd_temp,
                                  (16-mean_temp)/sd_temp,
                                  (18-mean_temp)/sd_temp,
                                  (20-mean_temp)/sd_temp),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bymintempxhab_WEBL.png",plot =  pl, width = 10, height = 6.6)






#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_webl,"figures/int_tab_abs_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl <- abs_lintemp

# For abs_change in WEBL, there is an interaction with both min and max.

### sample size


# ss_year_webl_abs <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
# ss_year_webl_abs %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
#   gtsave("figures/ss_year_webl_abs.html")


ss_year_webl_abs <- abs_lintemp@frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total = `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_webl_abs <- ss_year_webl_abs %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
    })
)


t_ss_year_webl_abs


samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(absbyhabitat_webl <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(absbyhabitat_webl,"figures/absbyhabitat_webl.html")



summary(abs_lintemp_noint)





data_webl = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_webl <- mean(data_webl %>% pull(meanmaxtempI))
sd_temp_webl <- sd(data_webl %>% pull(meanmaxtempI))


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

(fig4_webl <- ggpredict(abs_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl,
                       breaks = c((20-mean_temp_webl)/sd_temp_webl,
                                  (25-mean_temp_webl)/sd_temp_webl,
                                  (30-mean_temp_webl)/sd_temp_webl,
                                  (35-mean_temp_webl)/sd_temp_webl,
                                  (40-mean_temp_webl)/sd_temp_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    ylim(0,60) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/absbymaxtempxhab_WEBL.png",plot =  fig4_webl, width = 10, height = 6.6)




data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp <- mean(data %>% pull(meanmintempI))
sd_temp <- sd(data %>% pull(meanmintempI))


temp_trans <- trans_new("temp_trans",
                        transform = function(x){(x * sd_temp) + mean_temp},
                        inverse = function(x){x})

(pl <- ggpredict(abs_lintemp,terms = c("meanmintempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily min temp over preceding week (\u00b0C)") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL; min temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans,
                       breaks = c((8-mean_temp)/sd_temp,
                                  (10-mean_temp)/sd_temp,
                                  (12-mean_temp)/sd_temp,
                                  (14-mean_temp)/sd_temp,
                                  (16-mean_temp)/sd_temp,
                                  (18-mean_temp)/sd_temp,
                                  (20-mean_temp)/sd_temp),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/absbymintempxhab_WEBL.png",plot =  pl, width = 10, height = 6.6)



summary(abs_lintemp)



(weblabstrendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabstrendmax,"figures/weblabstrendmax.html")

data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

# (t <- emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(meanmaxtempI_scaled = (meanmaxtempI_scaled * sd(data$meanmaxtempI)) + mean(data$meanmaxtempI),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            meanmaxtempI_scaled = round(meanmaxtempI_scaled),
#            meanmaxtempI_scaled = paste0(meanmaxtempI_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`Max temperature` = "meanmaxtempI_scaled",`Predicted delta cort` = "response",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     group_by(`Max temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`Max temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`Max temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`Maximum TA_P` = if_else(`Maximum TA_P` == 0.000,"<0.001",as.character(`Maximum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"../figures/weblabsdeltamax.html")

# ((emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(abs_lintemp,specs = pairwise ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(abs_lintemp,formula = habitat ~ meanmaxtempI_scaled, at = list(meanmaxtempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()


### Minimum temperature


(t <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(t,"figures/weblabstrendmin.html")

data = dplyr::filter(g,Species == "WEBL",!is.na(meanmintempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmintempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmintempI_scaled_sq = meanmintempI_scaled * meanmintempI_scaled)

# (t <- emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(meanmintempI_scaled = (meanmintempI_scaled * sd(data$meanmintempI)) + mean(data$meanmintempI),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            meanmintempI_scaled = round(meanmintempI_scaled),
#            meanmintempI_scaled = paste0(meanmintempI_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`min temperature` = "meanmintempI_scaled",`Predicted delta cort` = "response",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     group_by(`min temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`min temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`min temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"../figures/weblabsdeltamin.html")

# ((emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(abs_lintemp,specs = pairwise ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(abs_lintemp,formula = habitat ~ meanmintempI_scaled, at = list(meanmintempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()

#### s2_cort
s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_webl,"figures/int_tab_s2_webl.html")

s2_webl <- s2_lintemp

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_webl,"figures/s2byhabitat_webl.html")


## trends
##
(webls2trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2trendmax,"figures/webls2trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


# ss_year_webl_s2 <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
# ss_year_webl_s2 %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
#   gtsave("figures/ss_year_webl_s2.html")


ss_year_webl_s2 <- s2_lintemp@frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total = `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_webl_s2 <- ss_year_webl_s2 %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
    })
)


t_ss_year_webl_s2

samp_s2_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_webl %>% gt()


dat_text_s2_webl <- data.frame(
  label = paste("N =",samp_s2_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_webl = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_s2_webl <- mean(data_s2_webl %>% pull(meanmaxtempI))
sd_temp_s2_webl <- sd(data_s2_webl %>% pull(meanmaxtempI))


temp_trans_s2_webl <- trans_new("temp_trans_s2_webl",
                                transform = function(x){(x * sd_temp_s2_webl) + mean_temp_s2_webl},
                                inverse = function(x){x})

(fig_webl_s2 <- ggpredict(s2_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_webl,
                       breaks = c((20-mean_temp_s2_webl)/sd_temp_s2_webl,
                                  (25-mean_temp_s2_webl)/sd_temp_s2_webl,
                                  (30-mean_temp_s2_webl)/sd_temp_s2_webl,
                                  (35-mean_temp_s2_webl)/sd_temp_s2_webl,
                                  (40-mean_temp_s2_webl)/sd_temp_s2_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bymaxtempxhab_WEBL.png",plot =  fig_webl_s2, width = 10, height = 6.6)


summary(s2_lintemp)

#### use prior day temp to predict cort instead


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_priordayt_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_priordayt_webl,"figures/int_tab_s1_priordayt_webl.html")

s1_prior_day_webl <- s1_lintemp_addmax
check_collinearity(s1_lintemp_addmax)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_priordayt_webl <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_priordayt_webl,"figures/s1byhabitat_priordayt_webl.html")


## trends
##
(webls1_priordayt_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_priordayt_trendmax,"figures/webls1_priordayt_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_priordayt <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_priordayt.html")

samp_s1_priordayt_webl <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordayt_webl %>% gt()


dat_text_s1_priordayt_webl <- data.frame(
  label = paste("N =",samp_s1_priordayt_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_priordayt_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled)


mean_temp_s1_priordayt_webl <- mean(data_s1_priordayt_webl %>% pull(maxt_prior))
sd_temp_s1_priordayt_webl <- sd(data_s1_priordayt_webl %>% pull(maxt_prior))


temp_trans_s1_priordayt_webl <- trans_new("temp_trans_s1_priordayt_webl",
                                transform = function(x){(x * sd_temp_s1_priordayt_webl) + mean_temp_s1_priordayt_webl},
                                inverse = function(x){x})

(figs2_priordayt_webl <- ggpredict(s1_lintemp_addmax,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_priordayt_webl,
                       breaks = c((20-mean_temp_s1_priordayt_webl)/sd_temp_s1_priordayt_webl,
                                  (25-mean_temp_s1_priordayt_webl)/sd_temp_s1_priordayt_webl,
                                  (30-mean_temp_s1_priordayt_webl)/sd_temp_s1_priordayt_webl,
                                  (35-mean_temp_s1_priordayt_webl)/sd_temp_s1_priordayt_webl,
                                  (40-mean_temp_s1_priordayt_webl)/sd_temp_s1_priordayt_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bypriordaytxhab_WEBL.png",plot =  figs2_priordayt_webl, width = 10, height = 6.6)

#### use prior day temp to predict cort instead


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_priordayt_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_priordayt_webl,"figures/int_tab_s2_priordayt_webl.html")

s2_prior_day_webl <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_priordayt_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_priordayt_webl,"figures/s2byhabitat_priordayt_webl.html")


## trends
##
(webls2_priordayt_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_priordayt_trendmax,"figures/webls2_priordayt_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_priordayt <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_priordayt.html")

samp_s2_priordayt_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordayt_webl %>% gt()


dat_text_s2_priordayt_webl <- data.frame(
  label = paste("N =",samp_s2_priordayt_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_priordayt_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled)


mean_temp_s2_priordayt_webl <- mean(data_s2_priordayt_webl %>% pull(maxt_prior))
sd_temp_s2_priordayt_webl <- sd(data_s2_priordayt_webl %>% pull(maxt_prior))


temp_trans_s2_priordayt_webl <- trans_new("temp_trans_s2_priordayt_webl",
                                          transform = function(x){(x * sd_temp_s2_priordayt_webl) + mean_temp_s2_priordayt_webl},
                                          inverse = function(x){x})

(fig_priordayt_webl_s2 <- ggpredict(s2_lintemp_addmax,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_priordayt_webl,
                       breaks = c((20-mean_temp_s2_priordayt_webl)/sd_temp_s2_priordayt_webl,
                                  (25-mean_temp_s2_priordayt_webl)/sd_temp_s2_priordayt_webl,
                                  (30-mean_temp_s2_priordayt_webl)/sd_temp_s2_priordayt_webl,
                                  (35-mean_temp_s2_priordayt_webl)/sd_temp_s2_priordayt_webl,
                                  (40-mean_temp_s2_priordayt_webl)/sd_temp_s2_priordayt_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bypriordaytxhab_WEBL.png",plot =  fig_priordayt_webl_s2, width = 10, height = 6.6)


#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_priordayt_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_priordayt_webl,"figures/int_tab_abs_priordayt_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_priordayt <- abs_lintemp_noint

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_priordayt <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_priordayt.html")

samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abspriordayt_byhabitat_webl <- emmeans(abs_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abspriordayt_byhabitat_webl,"figures/abspriordayt_byhabitat_webl.html")



summary(abs_lintemp_noint)





data_webl_priordayt = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled)


mean_temp_webl_priordayt <- mean(data_webl_priordayt %>% pull(maxt_prior))
sd_temp_webl_priordayt <- sd(data_webl_priordayt %>% pull(maxt_prior))


temp_trans_webl_priordayt <- trans_new("temp_trans_webl_priordayt",
                             transform = function(x){(x * sd_temp_webl_priordayt) + mean_temp_webl_priordayt},
                             inverse = function(x){x})

(fig4_webl_priordayt <- ggpredict(abs_lintemp_noint,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl_priordayt,
                       breaks = c((20-mean_temp_webl_priordayt)/sd_temp_webl_priordayt,
                                  (25-mean_temp_webl_priordayt)/sd_temp_webl_priordayt,
                                  (30-mean_temp_webl_priordayt)/sd_temp_webl_priordayt,
                                  (35-mean_temp_webl_priordayt)/sd_temp_webl_priordayt,
                                  (40-mean_temp_webl_priordayt)/sd_temp_webl_priordayt),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_priordayt_xhab_WEBL.png",plot =  fig4_webl_priordayt, width = 10, height = 6.6)




(weblabs_priordayt_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_priordayt_trendmax,"figures/weblabs_priordayt_trendmax.html")


#### use prior day heat index to predict cort instead


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                               mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled + minhi_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                     mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_priordaymaxhhi_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_priordaymaxhhi_webl,"figures/int_tab_s1_priordaymaxhhi_webl.html")

s1_priordaymaxhhi_webl <- s1_lintemp_addmax
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_priordaymaxhhi_webl <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_priordaymaxhhi_webl,"figures/s1byhabitat_priordaymaxhhi_webl.html")


## trends
##
(webls1_priordaymaxhhi_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_priordaymaxhhi_trendmax,"figures/webls1_priordaymaxhhi_trendmax.html")


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_priordaymaxhhi <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_priordaymaxhhi.html")

samp_s1_priordaymaxhhi_webl <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordaymaxhhi_webl %>% gt()


dat_text_s1_priordaymaxhhi_webl <- data.frame(
  label = paste("N =",samp_s1_priordaymaxhhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_priordaymaxhhi_webl =  dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_s1_priordaymaxhhi_webl <- mean(data_s1_priordaymaxhhi_webl %>% pull(maxhi_prior))
sd_temp_s1_priordaymaxhhi_webl <- sd(data_s1_priordaymaxhhi_webl %>% pull(maxhi_prior))


temp_trans_s1_priordaymaxhhi_webl <- trans_new("temp_trans_s1_priordaymaxhhi_webl",
                                          transform = function(x){(x * sd_temp_s1_priordaymaxhhi_webl) + mean_temp_s1_priordaymaxhhi_webl},
                                          inverse = function(x){x})

(figs2_priordaymaxhhi_webl <- ggpredict(s1_lintemp_addmax,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior day (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_priordaymaxhhi_webl,
                       breaks = c((20-mean_temp_s1_priordaymaxhhi_webl)/sd_temp_s1_priordaymaxhhi_webl,
                                  (40-mean_temp_s1_priordaymaxhhi_webl)/sd_temp_s1_priordaymaxhhi_webl,
                                  (60-mean_temp_s1_priordaymaxhhi_webl)/sd_temp_s1_priordaymaxhhi_webl,
                                  (80-mean_temp_s1_priordaymaxhhi_webl)/sd_temp_s1_priordaymaxhhi_webl,
                                  (100-mean_temp_s1_priordaymaxhhi_webl)/sd_temp_s1_priordaymaxhhi_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_priordaymaxhhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bypriordaymaxhhixhab_WEBL.png",plot =  figs2_priordaymaxhhi_webl, width = 10, height = 6.6)

#### s2 cort

s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                               mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled + minhi_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                     mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_priordaymaxhhi_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_priordaymaxhhi_webl,"figures/int_tab_s2_priordaymaxhhi_webl.html")

s2_priordaymaxhhi_webl <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_priordaymaxhhi_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_priordaymaxhhi_webl,"figures/s2byhabitat_priordaymaxhhi_webl.html")


## trends
##
(webls2_priordaymaxhhi_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_priordaymaxhhi_trendmax,"figures/webls2_priordaymaxhhi_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_priordaymaxhhi <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_priordaymaxhhi.html")

samp_s2_priordaymaxhhi_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordaymaxhhi_webl %>% gt()


dat_text_s2_priordaymaxhhi_webl <- data.frame(
  label = paste("N =",samp_s2_priordaymaxhhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_priordaymaxhhi_webl =  dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_s2_priordaymaxhhi_webl <- mean(data_s2_priordaymaxhhi_webl %>% pull(maxhi_prior))
sd_temp_s2_priordaymaxhhi_webl <- sd(data_s2_priordaymaxhhi_webl %>% pull(maxhi_prior))


temp_trans_s2_priordaymaxhhi_webl <- trans_new("temp_trans_s2_priordaymaxhhi_webl",
                                               transform = function(x){(x * sd_temp_s2_priordaymaxhhi_webl) + mean_temp_s2_priordaymaxhhi_webl},
                                               inverse = function(x){x})

(fig_priordaymaxhhi_webl_s2 <- ggpredict(s2_lintemp_noint,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior day (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_priordaymaxhhi_webl,
                       breaks = c((20-mean_temp_s2_priordaymaxhhi_webl)/sd_temp_s2_priordaymaxhhi_webl,
                                  (40-mean_temp_s2_priordaymaxhhi_webl)/sd_temp_s2_priordaymaxhhi_webl,
                                  (60-mean_temp_s2_priordaymaxhhi_webl)/sd_temp_s2_priordaymaxhhi_webl,
                                  (80-mean_temp_s2_priordaymaxhhi_webl)/sd_temp_s2_priordaymaxhhi_webl,
                                  (100-mean_temp_s2_priordaymaxhhi_webl)/sd_temp_s2_priordaymaxhhi_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_priordaymaxhhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bypriordaymaxhhixhab_WEBL.png",plot =  fig_priordaymaxhhi_webl_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                       mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                       mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled + minhi_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_priordaymaxhhi_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_priordaymaxhhi_webl,"figures/int_tab_abs_priordaymaxhhi_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_priordaymaxhhi <- abs_lintemp_noint

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_priordaymaxhhi <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_priordaymaxhhi.html")

samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_priordaymaxhhi_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abspriordaymaxhhi_byhabitat_webl <- emmeans(abs_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abspriordaymaxhhi_byhabitat_webl,"figures/abspriordaymaxhhi_byhabitat_webl.html")



summary(abs_lintemp_noint)





data_webl_priordaymaxhhi = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(minhi_prior)) %>%
  mutate(across(c(gweight,maxhi_prior,minhi_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_webl_priordaymaxhhi <- mean(data_webl_priordaymaxhhi %>% pull(maxhi_prior))
sd_temp_webl_priordaymaxhhi <- sd(data_webl_priordaymaxhhi %>% pull(maxhi_prior))


temp_trans_webl_priordaymaxhhi <- trans_new("temp_trans_webl_priordaymaxhhi",
                                       transform = function(x){(x * sd_temp_webl_priordaymaxhhi) + mean_temp_webl_priordaymaxhhi},
                                       inverse = function(x){x})

(fig4_webl_priordaymaxhhi <- ggpredict(abs_lintemp_noint,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior day (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl_priordaymaxhhi,
                       breaks = c((20-mean_temp_webl_priordaymaxhhi)/sd_temp_webl_priordaymaxhhi,
                                  (40-mean_temp_webl_priordaymaxhhi)/sd_temp_webl_priordaymaxhhi,
                                  (60-mean_temp_webl_priordaymaxhhi)/sd_temp_webl_priordaymaxhhi,
                                  (80-mean_temp_webl_priordaymaxhhi)/sd_temp_webl_priordaymaxhhi,
                                  (100-mean_temp_webl_priordaymaxhhi)/sd_temp_webl_priordaymaxhhi),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_priordaymaxhhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_priordaymaxhhi_xhab_WEBL.png",plot =  fig4_webl_priordaymaxhhi, width = 10, height = 6.6)




(weblabs_priordaymaxhhi_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_priordaymaxhhi_trendmax,"figures/weblabs_priordaymaxhhi_trendmax.html")

#### use prior week heat index to predict cort instead


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_weekhi_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_weekhi_webl,"figures/int_tab_s1_weekhi_webl.html")

s1_weekhi_webl <- s1_lintemp_noint
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_weekhi_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_weekhi_webl,"figures/s1byhabitat_weekhi_webl.html")


## trends
##
(webls1_weekhi_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_weekhi_trendmax,"figures/webls1_weekhi_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_weekhi <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_weekhi.html")

samp_s1_weekhi_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_weekhi_webl %>% gt()


dat_text_s1_weekhi_webl <- data.frame(
  label = paste("N =",samp_s1_weekhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_weekhi_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_s1_weekhi_webl <- mean(data_s1_weekhi_webl %>% pull(maxhi_week))
sd_temp_s1_weekhi_webl <- sd(data_s1_weekhi_webl %>% pull(maxhi_week))


temp_trans_s1_weekhi_webl <- trans_new("temp_trans_s1_weekhi_webl",
                                          transform = function(x){(x * sd_temp_s1_weekhi_webl) + mean_temp_s1_weekhi_webl},
                                          inverse = function(x){x})

(figs2_weekhi_webl <- ggpredict(s1_lintemp_noint,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_weekhi_webl,
                       breaks = c((20-mean_temp_s1_weekhi_webl)/sd_temp_s1_weekhi_webl,
                                  (30-mean_temp_s1_weekhi_webl)/sd_temp_s1_weekhi_webl,
                                  (40-mean_temp_s1_weekhi_webl)/sd_temp_s1_weekhi_webl,
                                  (50-mean_temp_s1_weekhi_webl)/sd_temp_s1_weekhi_webl,
                                  (60-mean_temp_s1_weekhi_webl)/sd_temp_s1_weekhi_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_weekhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1byweekhixhab_WEBL.png",plot =  figs2_weekhi_webl, width = 10, height = 6.6)


#### s2 cort


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_weekhi_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_weekhi_webl,"figures/int_tab_s2_weekhi_webl.html")

s2_weekhi_webl <- s2_lintemp
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_weekhi_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_weekhi_webl,"figures/s2byhabitat_weekhi_webl.html")


## trends
##
(webls2_weekhi_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_weekhi_trendmax,"figures/webls2_weekhi_trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_weekhi <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_weekhi.html")

samp_s2_weekhi_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_weekhi_webl %>% gt()


dat_text_s2_weekhi_webl <- data.frame(
  label = paste("N =",samp_s2_weekhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_weekhi_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_s2_weekhi_webl <- mean(data_s2_weekhi_webl %>% pull(maxhi_week))
sd_temp_s2_weekhi_webl <- sd(data_s2_weekhi_webl %>% pull(maxhi_week))


temp_trans_s2_weekhi_webl <- trans_new("temp_trans_s2_weekhi_webl",
                                       transform = function(x){(x * sd_temp_s2_weekhi_webl) + mean_temp_s2_weekhi_webl},
                                       inverse = function(x){x})

(fig_weekhi_webl_s2 <- ggpredict(s2_lintemp,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_weekhi_webl,
                       breaks = c((20-mean_temp_s2_weekhi_webl)/sd_temp_s2_weekhi_webl,
                                  (30-mean_temp_s2_weekhi_webl)/sd_temp_s2_weekhi_webl,
                                  (40-mean_temp_s2_weekhi_webl)/sd_temp_s2_weekhi_webl,
                                  (50-mean_temp_s2_weekhi_webl)/sd_temp_s2_weekhi_webl,
                                  (60-mean_temp_s2_weekhi_webl)/sd_temp_s2_weekhi_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_weekhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2byweekhixhab_WEBL.png",plot =  fig_weekhi_webl_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_weekhi_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_weekhi_webl,"figures/int_tab_abs_weekhi_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_weekhi <- abs_lintemp

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_weekhi <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_weekhi.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl_weekhi <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(absweekhi_byhabitat_webl <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(absweekhi_byhabitat_webl,"figures/absweekhi_byhabitat_webl.html")



summary(abs_lintemp)





data_webl_weekhi = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_webl_weekhi <- mean(data_webl_weekhi %>% pull(maxhi_week))
sd_temp_webl_weekhi <- sd(data_webl_weekhi %>% pull(maxhi_week))


temp_trans_webl_weekhi <- trans_new("temp_trans_webl_weekhi",
                                       transform = function(x){(x * sd_temp_webl_weekhi) + mean_temp_webl_weekhi},
                                       inverse = function(x){x})

(fig4_webl_weekhi <- ggpredict(abs_lintemp,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl_weekhi,
                       breaks = c((20-mean_temp_webl_weekhi)/sd_temp_webl_weekhi,
                                  (30-mean_temp_webl_weekhi)/sd_temp_webl_weekhi,
                                  (40-mean_temp_webl_weekhi)/sd_temp_webl_weekhi,
                                  (50-mean_temp_webl_weekhi)/sd_temp_webl_weekhi,
                                  (60-mean_temp_webl_weekhi)/sd_temp_webl_weekhi),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_webl_weekhi, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_weekhi_xhab_WEBL.png",plot =  fig4_webl_weekhi, width = 10, height = 6.6)




(weblabs_weekhi_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_weekhi_trendmax,"figures/weblabs_weekhi_trendmax.html")

#### use cumulative prior day hi to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumhiday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumhiday_webl,"figures/int_tab_s1_cumhiday_webl.html")

s1_cumhiday_webl <- s1_lintemp_noint
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_cumhiday_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumhiday_webl,"figures/s1byhabitat_cumhiday_webl.html")


## trends
##
(webls1_cumhiday_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_cumhiday_trendmax,"figures/webls1_cumhiday_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_cumhiday <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_cumhiday.html")

samp_s1_cumhiday_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiday_webl %>% gt()


dat_text_s1_cumhiday_webl <- data.frame(
  label = paste("N =",samp_s1_cumhiday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumhiday_webl = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled)


mean_temp_s1_cumhiday_webl <- mean(data_s1_cumhiday_webl %>% pull(hihours_over_30hi_priorday))
sd_temp_s1_cumhiday_webl <- sd(data_s1_cumhiday_webl %>% pull(hihours_over_30hi_priorday))


temp_trans_s1_cumhiday_webl <- trans_new("temp_trans_s1_cumhiday_webl",
                                       transform = function(x){(x * sd_temp_s1_cumhiday_webl) + mean_temp_s1_cumhiday_webl},
                                       inverse = function(x){x})

(figs2_cumhiday_webl <- ggpredict(s1_lintemp_noint,terms = c("hihours_over_30hi_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumhiday_webl,
                       breaks = c((0-mean_temp_s1_cumhiday_webl)/sd_temp_s1_cumhiday_webl,
                                  (200-mean_temp_s1_cumhiday_webl)/sd_temp_s1_cumhiday_webl,
                                  (400-mean_temp_s1_cumhiday_webl)/sd_temp_s1_cumhiday_webl,
                                  (600-mean_temp_s1_cumhiday_webl)/sd_temp_s1_cumhiday_webl,
                                  (800-mean_temp_s1_cumhiday_webl)/sd_temp_s1_cumhiday_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumhiday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumhidayxhab_WEBL.png",plot =  figs2_cumhiday_webl, width = 10, height = 6.6)


#### use cumulative prior day hi to predict cort


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumhiday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumhiday_webl,"figures/int_tab_s2_cumhiday_webl.html")

s2_cumhiday_webl <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_cumhiday_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumhiday_webl,"figures/s2byhabitat_cumhiday_webl.html")


## trends
##
(webls2_cumhiday_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_cumhiday_trendmax,"figures/webls2_cumhiday_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_cumhiday <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_cumhiday.html")

samp_s2_cumhiday_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiday_webl %>% gt()


dat_text_s2_cumhiday_webl <- data.frame(
  label = paste("N =",samp_s2_cumhiday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumhiday_webl = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled)


mean_temp_s2_cumhiday_webl <- mean(data_s2_cumhiday_webl %>% pull(hihours_over_30hi_priorday))
sd_temp_s2_cumhiday_webl <- sd(data_s2_cumhiday_webl %>% pull(hihours_over_30hi_priorday))


temp_trans_s2_cumhiday_webl <- trans_new("temp_trans_s2_cumhiday_webl",
                                         transform = function(x){(x * sd_temp_s2_cumhiday_webl) + mean_temp_s2_cumhiday_webl},
                                         inverse = function(x){x})

(fig_cumhiday_webl_s2 <- ggpredict(s2_lintemp_noint,terms = c("hihours_over_30hi_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumhiday_webl,
                       breaks = c((0-mean_temp_s2_cumhiday_webl)/sd_temp_s2_cumhiday_webl,
                                  (200-mean_temp_s2_cumhiday_webl)/sd_temp_s2_cumhiday_webl,
                                  (400-mean_temp_s2_cumhiday_webl)/sd_temp_s2_cumhiday_webl,
                                  (600-mean_temp_s2_cumhiday_webl)/sd_temp_s2_cumhiday_webl,
                                  (800-mean_temp_s2_cumhiday_webl)/sd_temp_s2_cumhiday_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumhiday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumhidayxhab_WEBL.png",plot =  fig_cumhiday_webl_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumhiday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumhiday_webl,"figures/int_tab_abs_cumhiday_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_cumhiday <- abs_lintemp_noint

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_cumhiday <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_cumhiday.html")

samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl_cumhiday <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumhiday_byhabitat_webl <- emmeans(abs_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumhiday_byhabitat_webl,"figures/abscumhiday_byhabitat_webl.html")



summary(abs_lintemp_noint)





data_webl_cumhiday = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled)


mean_temp_webl_cumhiday <- mean(data_webl_cumhiday %>% pull(hihours_over_30hi_priorday))
sd_temp_webl_cumhiday <- sd(data_webl_cumhiday %>% pull(hihours_over_30hi_priorday))


temp_trans_webl_cumhiday <- trans_new("temp_trans_webl_cumhiday",
                                    transform = function(x){(x * sd_temp_webl_cumhiday) + mean_temp_webl_cumhiday},
                                    inverse = function(x){x})

(fig4_webl_cumhiday <- ggpredict(abs_lintemp_noint,terms = c("hihours_over_30hi_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index-hours >25\u00b0C over prior day") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl_cumhiday,
                       breaks = c((0-mean_temp_webl_cumhiday)/sd_temp_webl_cumhiday,
                                  (200-mean_temp_webl_cumhiday)/sd_temp_webl_cumhiday,
                                  (400-mean_temp_webl_cumhiday)/sd_temp_webl_cumhiday,
                                  (600-mean_temp_webl_cumhiday)/sd_temp_webl_cumhiday,
                                  (800-mean_temp_webl_cumhiday)/sd_temp_webl_cumhiday),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "solid")) +
    # ylim(0,60) +
    geom_text(data = dat_text_webl_cumhiday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumhiday_xhab_WEBL.png",plot =  fig4_webl_cumhiday, width = 10, height = 6.6)




(weblabs_cumhiday_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_cumhiday_trendmax,"figures/weblabs_cumhiday_trendmax.html")


#### use cumulative prior week hi to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumhiweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumhiweek_webl,"figures/int_tab_s1_cumhiweek_webl.html")

s1_cumhiweek_webl <- s1_lintemp_noint
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_cumhiweek_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumhiweek_webl,"figures/s1byhabitat_cumhiweek_webl.html")


## trends
##
(webls1_cumhiweek_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_cumhiweek_trendmax,"figures/webls1_cumhiweek_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_cumhiweek <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_cumhiweek.html")

samp_s1_cumhiweek_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiweek_webl %>% gt()


dat_text_s1_cumhiweek_webl <- data.frame(
  label = paste("N =",samp_s1_cumhiweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumhiweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_temp_s1_cumhiweek_webl <- mean(data_s1_cumhiweek_webl %>% pull(hihours_over_30hi_priorweek))
sd_temp_s1_cumhiweek_webl <- sd(data_s1_cumhiweek_webl %>% pull(hihours_over_30hi_priorweek))


temp_trans_s1_cumhiweek_webl <- trans_new("temp_trans_s1_cumhiweek_webl",
                                         transform = function(x){(x * sd_temp_s1_cumhiweek_webl) + mean_temp_s1_cumhiweek_webl},
                                         inverse = function(x){x})

(figs2_cumhiweek_webl <- ggpredict(s1_lintemp_noint,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumhiweek_webl,
                       breaks = c((0-mean_temp_s1_cumhiweek_webl)/sd_temp_s1_cumhiweek_webl,
                                  (1000-mean_temp_s1_cumhiweek_webl)/sd_temp_s1_cumhiweek_webl,
                                  (2000-mean_temp_s1_cumhiweek_webl)/sd_temp_s1_cumhiweek_webl,
                                  (3000-mean_temp_s1_cumhiweek_webl)/sd_temp_s1_cumhiweek_webl,
                                  (4000-mean_temp_s1_cumhiweek_webl)/sd_temp_s1_cumhiweek_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumhiweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumhiweekxhab_WEBL.png",plot =  figs2_cumhiweek_webl, width = 10, height = 6.6)


#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumhiweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumhiweek_webl,"figures/int_tab_s2_cumhiweek_webl.html")

s2_cumhiweek_webl <- s2_lintemp
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_cumhiweek_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumhiweek_webl,"figures/s2byhabitat_cumhiweek_webl.html")


## trends
##
(webls2_cumhiweek_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_cumhiweek_trendmax,"figures/webls2_cumhiweek_trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_cumhiweek <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_cumhiweek.html")

samp_s2_cumhiweek_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiweek_webl %>% gt()


dat_text_s2_cumhiweek_webl <- data.frame(
  label = paste("N =",samp_s2_cumhiweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumhiweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_temp_s2_cumhiweek_webl <- mean(data_s2_cumhiweek_webl %>% pull(hihours_over_30hi_priorweek))
sd_temp_s2_cumhiweek_webl <- sd(data_s2_cumhiweek_webl %>% pull(hihours_over_30hi_priorweek))


temp_trans_s2_cumhiweek_webl <- trans_new("temp_trans_s2_cumhiweek_webl",
                                          transform = function(x){(x * sd_temp_s2_cumhiweek_webl) + mean_temp_s2_cumhiweek_webl},
                                          inverse = function(x){x})

(fig_cumhiweek_webl_s2 <- ggpredict(s2_lintemp,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumhiweek_webl,
                       breaks = c((0-mean_temp_s2_cumhiweek_webl)/sd_temp_s2_cumhiweek_webl,
                                  (1000-mean_temp_s2_cumhiweek_webl)/sd_temp_s2_cumhiweek_webl,
                                  (2000-mean_temp_s2_cumhiweek_webl)/sd_temp_s2_cumhiweek_webl,
                                  (3000-mean_temp_s2_cumhiweek_webl)/sd_temp_s2_cumhiweek_webl,
                                  (4000-mean_temp_s2_cumhiweek_webl)/sd_temp_s2_cumhiweek_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dashed","solid","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumhiweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumhiweekxhab_WEBL.png",plot =  fig_cumhiweek_webl_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumhiweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumhiweek_webl,"figures/int_tab_abs_cumhiweek_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_cumhiweek <- abs_lintemp

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_cumhiweek <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_cumhiweek.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_abs_text_webl_cumhiweek <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumhiweek_byhabitat_webl <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumhiweek_byhabitat_webl,"figures/abscumhiweek_byhabitat_webl.html")



summary(abs_lintemp)





data_abs_webl_cumhiweek = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_abs_temp_webl_cumhiweek <- mean(data_abs_webl_cumhiweek %>% pull(hihours_over_30hi_priorweek))
sd_abs_temp_webl_cumhiweek <- sd(data_abs_webl_cumhiweek %>% pull(hihours_over_30hi_priorweek))


temp_abs_trans_webl_cumhiweek <- trans_new("temp_trans_webl_cumhiweek",
                                      transform = function(x){(x * sd_abs_temp_webl_cumhiweek) + mean_abs_temp_webl_cumhiweek},
                                      inverse = function(x){x})

(fig4_webl_cumhiweek <- ggpredict(abs_lintemp,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index > 25\u00b0C over prior week") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_abs_trans_webl_cumhiweek,
                       breaks = c((0-mean_abs_temp_webl_cumhiweek)/sd_abs_temp_webl_cumhiweek,
                                  (1000-mean_abs_temp_webl_cumhiweek)/sd_abs_temp_webl_cumhiweek,
                                  (2000-mean_abs_temp_webl_cumhiweek)/sd_abs_temp_webl_cumhiweek,
                                  (3000-mean_abs_temp_webl_cumhiweek)/sd_abs_temp_webl_cumhiweek,
                                  (4000-mean_abs_temp_webl_cumhiweek)/sd_abs_temp_webl_cumhiweek),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_abs_text_webl_cumhiweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumhiweek_xhab_WEBL.png",plot =  fig4_webl_cumhiweek, width = 10, height = 6.6)




(weblabs_cumhiweek_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_cumhiweek_trendmax,"figures/weblabs_cumhiweek_trendmax.html")

#### use cumulative prior day temp to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumdegreeday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumdegreeday_webl,"figures/int_tab_s1_cumdegreeday_webl.html")

s1_cumdegreeday_webl <- s1_lintemp_noint
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeday_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumdegreeday_webl,"figures/s1byhabitat_cumdegreeday_webl.html")


## trends
##
(webls1_cumdegreeday_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_cumdegreeday_trendmax,"figures/webls1_cumdegreeday_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_cumdegreeday <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_cumdegreeday.html")

samp_s1_cumdegreeday_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeday_webl %>% gt()


dat_text_s1_cumdegreeday_webl <- data.frame(
  label = paste("N =",samp_s1_cumdegreeday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumdegreeday_webl = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled)


mean_temp_s1_cumdegreeday_webl <- mean(data_s1_cumdegreeday_webl %>% pull(degreehours_over_30C_priorday))
sd_temp_s1_cumdegreeday_webl <- sd(data_s1_cumdegreeday_webl %>% pull(degreehours_over_30C_priorday))


temp_trans_s1_cumdegreeday_webl <- trans_new("temp_trans_s1_cumdegreeday_webl",
                                         transform = function(x){(x * sd_temp_s1_cumdegreeday_webl) + mean_temp_s1_cumdegreeday_webl},
                                         inverse = function(x){x})

(figs2_cumdegreeday_webl <- ggpredict(s1_lintemp_noint,terms = c("degreehours_over_30C_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumdegreeday_webl,
                       breaks = c((0-mean_temp_s1_cumdegreeday_webl)/sd_temp_s1_cumdegreeday_webl,
                                  (100-mean_temp_s1_cumdegreeday_webl)/sd_temp_s1_cumdegreeday_webl,
                                  (200-mean_temp_s1_cumdegreeday_webl)/sd_temp_s1_cumdegreeday_webl,
                                  (300-mean_temp_s1_cumdegreeday_webl)/sd_temp_s1_cumdegreeday_webl,
                                  (400-mean_temp_s1_cumdegreeday_webl)/sd_temp_s1_cumdegreeday_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumdegreeday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumdegreedayxhab_WEBL.png",plot =  figs2_cumdegreeday_webl, width = 10, height = 6.6)

#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumdegreeday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumdegreeday_webl,"figures/int_tab_s2_cumdegreeday_webl.html")

s2_cumdegreeday_webl <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeday_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumdegreeday_webl,"figures/s2byhabitat_cumdegreeday_webl.html")


## trends
##
(webls2_cumdegreeday_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_cumdegreeday_trendmax,"figures/webls2_cumdegreeday_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_cumdegreeday <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_cumdegreeday.html")

samp_s2_cumdegreeday_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeday_webl %>% gt()


dat_text_s2_cumdegreeday_webl <- data.frame(
  label = paste("N =",samp_s2_cumdegreeday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumdegreeday_webl = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled)


mean_temp_s2_cumdegreeday_webl <- mean(data_s2_cumdegreeday_webl %>% pull(degreehours_over_30C_priorday))
sd_temp_s2_cumdegreeday_webl <- sd(data_s2_cumdegreeday_webl %>% pull(degreehours_over_30C_priorday))


temp_trans_s2_cumdegreeday_webl <- trans_new("temp_trans_s2_cumdegreeday_webl",
                                             transform = function(x){(x * sd_temp_s2_cumdegreeday_webl) + mean_temp_s2_cumdegreeday_webl},
                                             inverse = function(x){x})

(fig_cumdegreeday_webl_s2 <- ggpredict(s2_lintemp_noint,terms = c("degreehours_over_30C_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumdegreeday_webl,
                       breaks = c((0-mean_temp_s2_cumdegreeday_webl)/sd_temp_s2_cumdegreeday_webl,
                                  (100-mean_temp_s2_cumdegreeday_webl)/sd_temp_s2_cumdegreeday_webl,
                                  (200-mean_temp_s2_cumdegreeday_webl)/sd_temp_s2_cumdegreeday_webl,
                                  (300-mean_temp_s2_cumdegreeday_webl)/sd_temp_s2_cumdegreeday_webl,
                                  (400-mean_temp_s2_cumdegreeday_webl)/sd_temp_s2_cumdegreeday_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumdegreeday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumdegreedayxhab_WEBL.png",plot =  fig_cumdegreeday_webl_s2, width = 10, height = 6.6)


#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumdegreeday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumdegreeday_webl,"figures/int_tab_abs_cumdegreeday_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_cumdegreeday <- abs_lintemp

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_cumdegreeday <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_cumdegreeday.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl_cumdegreeday <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumdegreeday_byhabitat_webl <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumdegreeday_byhabitat_webl,"figures/abscumdegreeday_byhabitat_webl.html")



summary(abs_lintemp)





data_webl_cumdegreeday = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled)


mean_temp_webl_cumdegreeday <- mean(data_webl_cumdegreeday %>% pull(degreehours_over_30C_priorday))
sd_temp_webl_cumdegreeday <- sd(data_webl_cumdegreeday %>% pull(degreehours_over_30C_priorday))


temp_trans_webl_cumdegreeday <- trans_new("temp_trans_webl_cumdegreeday",
                                      transform = function(x){(x * sd_temp_webl_cumdegreeday) + mean_temp_webl_cumdegreeday},
                                      inverse = function(x){x})

(fig4_webl_cumdegreeday <- ggpredict(abs_lintemp,terms = c("degreehours_over_30C_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl_cumdegreeday,
                       breaks = c((0-mean_temp_webl_cumdegreeday)/sd_temp_webl_cumdegreeday,
                                  (100-mean_temp_webl_cumdegreeday)/sd_temp_webl_cumdegreeday,
                                  (200-mean_temp_webl_cumdegreeday)/sd_temp_webl_cumdegreeday,
                                  (300-mean_temp_webl_cumdegreeday)/sd_temp_webl_cumdegreeday,
                                  (400-mean_temp_webl_cumdegreeday)/sd_temp_webl_cumdegreeday),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_webl_cumdegreeday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumdegreeday_xhab_WEBL.png",plot =  fig4_webl_cumdegreeday, width = 10, height = 6.6)




(weblabs_cumdegreeday_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_cumdegreeday_trendmax,"figures/weblabs_cumdegreeday_trendmax.html")

#### use cumulative prior week temp to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumdegreeweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumdegreeweek_webl,"figures/int_tab_s1_cumdegreeweek_webl.html")

s1_cumdegreeweek_webl <- s1_lintemp_addmin
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeweek_webl <- emmeans(s1_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumdegreeweek_webl,"figures/s1byhabitat_cumdegreeweek_webl.html")


## trends
##
(webls1_cumdegreeweek_trendmax <- emtrends(s1_lintemp_addmin,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls1_cumdegreeweek_trendmax,"figures/webls1_cumdegreeweek_trendmax.html")


summary(s1_lintemp_addmin)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_webl_s1_cumdegreeweek <- s1_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s1_cumdegreeweek.html")

samp_s1_cumdegreeweek_webl <- s1_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeweek_webl %>% gt()


dat_text_s1_cumdegreeweek_webl <- data.frame(
  label = paste("N =",samp_s1_cumdegreeweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumdegreeweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_s1_cumdegreeweek_webl <- mean(data_s1_cumdegreeweek_webl %>% pull(degreehours_over_30C_priorweek))
sd_temp_s1_cumdegreeweek_webl <- sd(data_s1_cumdegreeweek_webl %>% pull(degreehours_over_30C_priorweek))


temp_trans_s1_cumdegreeweek_webl <- trans_new("temp_trans_s1_cumdegreeweek_webl",
                                         transform = function(x){(x * sd_temp_s1_cumdegreeweek_webl) + mean_temp_s1_cumdegreeweek_webl},
                                         inverse = function(x){x})

(figs2_cumdegreeweek_webl <- ggpredict(s1_lintemp_addmin,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumdegreeweek_webl,
                       breaks = c((0-mean_temp_s1_cumdegreeweek_webl)/sd_temp_s1_cumdegreeweek_webl,
                                  (1000-mean_temp_s1_cumdegreeweek_webl)/sd_temp_s1_cumdegreeweek_webl,
                                  (2000-mean_temp_s1_cumdegreeweek_webl)/sd_temp_s1_cumdegreeweek_webl,
                                  (3000-mean_temp_s1_cumdegreeweek_webl)/sd_temp_s1_cumdegreeweek_webl,
                                  (4000-mean_temp_s1_cumdegreeweek_webl)/sd_temp_s1_cumdegreeweek_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("solid","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumdegreeweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumdegreeweekxhab_WEBL.png",plot =  figs2_cumdegreeweek_webl, width = 10, height = 6.6)


#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumdegreeweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumdegreeweek_webl,"figures/int_tab_s2_cumdegreeweek_webl.html")

s2_cumdegreeweek_webl <- s2_lintemp
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling WEBL.

## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeweek_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumdegreeweek_webl,"figures/s2byhabitat_cumdegreeweek_webl.html")


## trends
##
(webls2_cumdegreeweek_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(webls2_cumdegreeweek_trendmax,"figures/webls2_cumdegreeweek_trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_webl_s2_cumdegreeweek <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_s2_cumdegreeweek.html")

samp_s2_cumdegreeweek_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeweek_webl %>% gt()


dat_text_s2_cumdegreeweek_webl <- data.frame(
  label = paste("N =",samp_s2_cumdegreeweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumdegreeweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_s2_cumdegreeweek_webl <- mean(data_s2_cumdegreeweek_webl %>% pull(degreehours_over_30C_priorweek))
sd_temp_s2_cumdegreeweek_webl <- sd(data_s2_cumdegreeweek_webl %>% pull(degreehours_over_30C_priorweek))


temp_trans_s2_cumdegreeweek_webl <- trans_new("temp_trans_s2_cumdegreeweek_webl",
                                              transform = function(x){(x * sd_temp_s2_cumdegreeweek_webl) + mean_temp_s2_cumdegreeweek_webl},
                                              inverse = function(x){x})

(fig_cumdegreeweek_webl_s2 <- ggpredict(s2_lintemp,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumdegreeweek_webl,
                       breaks = c((0-mean_temp_s2_cumdegreeweek_webl)/sd_temp_s2_cumdegreeweek_webl,
                                  (1000-mean_temp_s2_cumdegreeweek_webl)/sd_temp_s2_cumdegreeweek_webl,
                                  (2000-mean_temp_s2_cumdegreeweek_webl)/sd_temp_s2_cumdegreeweek_webl,
                                  (3000-mean_temp_s2_cumdegreeweek_webl)/sd_temp_s2_cumdegreeweek_webl,
                                  (4000-mean_temp_s2_cumdegreeweek_webl)/sd_temp_s2_cumdegreeweek_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumdegreeweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumdegreeweekxhab_WEBL.png",plot =  fig_cumdegreeweek_webl_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumdegreeweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumdegreeweek_webl,"figures/int_tab_abs_cumdegreeweek_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl_cumdegreeweek <- abs_lintemp

# For abs_change in WEBL, there are no interactions

### sample size


ss_year_webl_abs_cumdegreeweek <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_webl_abs_cumdegreeweek.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl_cumdegreeweek <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumdegreeweek_byhabitat_webl <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumdegreeweek_byhabitat_webl,"figures/abscumdegreeweek_byhabitat_webl.html")



summary(abs_lintemp)





data_webl_cumdegreeweek = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_webl_cumdegreeweek <- mean(data_webl_cumdegreeweek %>% pull(degreehours_over_30C_priorweek))
sd_temp_webl_cumdegreeweek <- sd(data_webl_cumdegreeweek %>% pull(degreehours_over_30C_priorweek))


temp_trans_webl_cumdegreeweek <- trans_new("temp_trans_webl_cumdegreeweek",
                                      transform = function(x){(x * sd_temp_webl_cumdegreeweek) + mean_temp_webl_cumdegreeweek},
                                      inverse = function(x){x})

(fig4_webl_cumdegreeweek <- ggpredict(abs_lintemp,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "WEBL cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_webl_cumdegreeweek,
                       breaks = c((0-mean_temp_webl_cumdegreeweek)/sd_temp_webl_cumdegreeweek,
                                  (1000-mean_temp_webl_cumdegreeweek)/sd_temp_webl_cumdegreeweek,
                                  (2000-mean_temp_webl_cumdegreeweek)/sd_temp_webl_cumdegreeweek,
                                  (3000-mean_temp_webl_cumdegreeweek)/sd_temp_webl_cumdegreeweek,
                                  (4000-mean_temp_webl_cumdegreeweek)/sd_temp_webl_cumdegreeweek),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_webl_cumdegreeweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumdegreeweek_xhab_WEBL.png",plot =  fig4_webl_cumdegreeweek, width = 10, height = 6.6)




(weblabs_cumdegreeweek_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblabs_cumdegreeweek_trendmax,"figures/weblabs_cumdegreeweek_trendmax.html")

### TRES

#### cort_s1


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_tres,"figures/int_tab_s1_tres.html")

s1_tres <- s1_lintemp_addmax

### Sample sizes:


# ss_year_tres_s1 <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
# ss_year_tres_s1 %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
#   gtsave("figures/ss_year_tres_s1.html")


ss_year_tres_s1 <- s1_lintemp_addmax@frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total = `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_tres_s1 <- ss_year_tres_s1 %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
    })
)


t_ss_year_tres_s1

samp_s1_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_s1_tres <- data.frame(
  label = paste("N =",samp_s1_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


# For cort_s1 in TRES, there is an interaction with min but not with max.
#
# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)

## trends
##
(tress1trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1trendmax,"figures/tress1trendmax.html")



data_s1_tres = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_s1_tres <- mean(data_s1_tres %>% pull(meanmaxtempI))
sd_temp_s1_tres <- sd(data_s1_tres %>% pull(meanmaxtempI))


temp_trans_s1_tres <- trans_new("temp_trans_s1_tres",
                                transform = function(x){(x * sd_temp_s1_tres) + mean_temp_s1_tres},
                                inverse = function(x){x})

(figs2_tres <- ggpredict(s1_lintemp_addmax,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_tres,
                       breaks = c((20-mean_temp_s1_tres)/sd_temp_s1_tres,
                                  (25-mean_temp_s1_tres)/sd_temp_s1_tres,
                                  (30-mean_temp_s1_tres)/sd_temp_s1_tres,
                                  (35-mean_temp_s1_tres)/sd_temp_s1_tres,
                                  (40-mean_temp_s1_tres)/sd_temp_s1_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bymaxtempxhab_TRES.png",plot =  figs2_tres, width = 10, height = 6.6)








data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp <- mean(data %>% pull(meanmintempI))
sd_temp <- sd(data %>% pull(meanmintempI))


temp_trans <- trans_new("temp_trans",
                        transform = function(x){(x * sd_temp) + mean_temp},
                        inverse = function(x){x})

(pl <- ggpredict(s1_lintemp_addmax,terms = c("meanmintempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily min temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans,
                       breaks = c((8-mean_temp)/sd_temp,
                                  (10-mean_temp)/sd_temp,
                                  (12-mean_temp)/sd_temp,
                                  (14-mean_temp)/sd_temp,
                                  (16-mean_temp)/sd_temp,
                                  (18-mean_temp)/sd_temp,
                                  (20-mean_temp)/sd_temp),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bymintempxhab_TRES.png",plot =  pl, width = 10, height = 6.6)


## Emmeans to check for effect of habitat


(s1byhabitat_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_tres,"figures/s1byhabitat_tres.html")


### Minimum temperature


(t <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(t,"figures/tresabstrendmin.html")

# data = dplyr::filter(g,Species == "TRES",!is.na(meanmintempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmintempI,meanmintempI,juliandate,brood_size,age),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          meanmintempI_scaled_sq = meanmintempI_scaled * meanmintempI_scaled)
#
# (t <- emmeans(s1_lintemp_addmax,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(meanmintempI_scaled = (meanmintempI_scaled * sd(data$meanmintempI)) + mean(data$meanmintempI),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            meanmintempI_scaled = round(meanmintempI_scaled),
#            meanmintempI_scaled = paste0(meanmintempI_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`min temperature` = "meanmintempI_scaled",`Predicted delta cort` = "response",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     group_by(`min temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`min temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`min temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"figures/tresabsdeltamin.html")

# ((emmeans(s1_lintemp_addmax,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(s1_lintemp_addmax,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(s1_lintemp_addmax,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(s1_lintemp_addmax,specs = pairwise ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(s1_lintemp_addmax,formula = habitat ~ meanmintempI_scaled, at = list(meanmintempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()

#### cort_s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_tres,"figures/int_tab_s2_tres.html")

s2_tres <- s2_lintemp_addmin

### Sample sizes:


# ss_year_tres_s2 <- s2_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
# ss_year_tres_s2 %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
#   gtsave("figures/ss_year_tres_s2.html")


ss_year_tres_s2 <- s2_lintemp_addmin@frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total = `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_tres_s2 <- ss_year_tres_s2 %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
    })
)


t_ss_year_tres_s2

samp_s2_tres <- s2_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_s2_tres <- data.frame(
  label = paste("N =",samp_s2_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

(s2byhabitat_tres <- emmeans(s2_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_tres,"figures/s2byhabitat_tres.html")

# For cort_s2 in TRES, there is an interaction with min but not with max.
#
# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.


summary(s2_lintemp_addmin)
check_collinearity(s2_lintemp_noint)

## trends
##
(tress2trendmax <- emtrends(s2_lintemp_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2trendmax,"figures/tress2trendmax.html")



data_s2_tres = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_s2_tres <- mean(data_s2_tres %>% pull(meanmaxtempI))
sd_temp_s2_tres <- sd(data_s2_tres %>% pull(meanmaxtempI))


temp_trans_s2_tres <- trans_new("temp_trans_s2_tres",
                                transform = function(x){(x * sd_temp_s2_tres) + mean_temp_s2_tres},
                                inverse = function(x){x})

(fig_tres_s2 <- ggpredict(s2_lintemp_addmin,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_tres,
                       breaks = c((20-mean_temp_s2_tres)/sd_temp_s2_tres,
                                  (25-mean_temp_s2_tres)/sd_temp_s2_tres,
                                  (30-mean_temp_s2_tres)/sd_temp_s2_tres,
                                  (35-mean_temp_s2_tres)/sd_temp_s2_tres,
                                  (40-mean_temp_s2_tres)/sd_temp_s2_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bymaxtempxhab_TRES.png",plot =  fig_tres_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_tres,"figures/int_tab_abs_tres.html")

abs_tres <- abs_lintemp

# For abs_change in TRES, there is an interaction with max and marginally min.


### Sample sizes:


# ss_year_tres_abs <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
# ss_year_tres_abs %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
#   gtsave("figures/ss_year_tres_abs.html")


ss_year_tres_abs <- abs_lintemp@frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total = `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_tres_abs <- ss_year_tres_abs %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
    })
)


t_ss_year_tres_abs

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



summary(abs_lintemp)

check_collinearity(abs_lintemp_noint)




data_tres = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_tres <- mean(data_tres %>% pull(meanmaxtempI))
sd_temp_tres <- sd(data_tres %>% pull(meanmaxtempI))


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

(fig4_tres <- ggpredict(abs_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres,
                       breaks = c((20-mean_temp_tres)/sd_temp_tres,
                                  (25-mean_temp_tres)/sd_temp_tres,
                                  (30-mean_temp_tres)/sd_temp_tres,
                                  (35-mean_temp_tres)/sd_temp_tres,
                                  (40-mean_temp_tres)/sd_temp_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dashed","Row crop" = "dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/absbymaxtempxhab_TRES.png",plot = fig4_tres, width = 10, height = 6.6)









data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp <- mean(data %>% pull(meanmintempI))
sd_temp <- sd(data %>% pull(meanmintempI))


temp_trans <- trans_new("temp_trans",
                        transform = function(x){(x * sd_temp) + mean_temp},
                        inverse = function(x){x})

(pl <- ggpredict(abs_lintemp,terms = c("meanmintempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily min temp over preceding week (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans,
                       breaks = c((8-mean_temp)/sd_temp,
                                  (10-mean_temp)/sd_temp,
                                  (12-mean_temp)/sd_temp,
                                  (14-mean_temp)/sd_temp,
                                  (16-mean_temp)/sd_temp,
                                  (18-mean_temp)/sd_temp,
                                  (20-mean_temp)/sd_temp),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/absbymintempxhab_TRES.png",plot =  pl, width = 10, height = 6.6)


## Emmeans to check for effect of habitat


(absbyhabitat_tres <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(absbyhabitat_tres,"figures/absbyhabitat_tres.html")


#I think given this finding it might be helpful to subset the TRES data down to just grassland and row crop to compare them. Sample size for forest and orchard are so small.

### Maximum temperature


(tresabstrendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabstrendmax,"figures/tresabstrendmax.html")

data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmintempI_scaled_sq = meanmintempI_scaled * meanmintempI_scaled)

# (t <- emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(meanmaxtempI_scaled = (meanmaxtempI_scaled * sd(data$meanmaxtempI)) + mean(data$meanmaxtempI),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            meanmaxtempI_scaled = round(meanmaxtempI_scaled),
#            meanmaxtempI_scaled = paste0(meanmaxtempI_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`Max temperature` = "meanmaxtempI_scaled",`Predicted delta cort` = "response",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     group_by(`Max temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`Max temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`Max temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"../figures/tresabsdeltamax.html")

# ((emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(abs_lintemp,specs = pairwise ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(abs_lintemp,formula = habitat ~ meanmaxtempI_scaled, at = list(meanmaxtempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



### Minimum temperature


(t <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(t,"figures/tresabstrendmin.html")

data = dplyr::filter(g,Species == "TRES",!is.na(meanmintempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmintempI,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmintempI_scaled_sq = meanmintempI_scaled * meanmintempI_scaled)

# (t <- emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(meanmintempI_scaled = (meanmintempI_scaled * sd(data$meanmintempI)) + mean(data$meanmintempI),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            meanmintempI_scaled = round(meanmintempI_scaled),
#            meanmintempI_scaled = paste0(meanmintempI_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`Min temperature` = "meanmintempI_scaled",`Predicted delta cort` = "response",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     group_by(`Min temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`Min temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`Min temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"figures/tresabsdeltamin.html")

# ((emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(abs_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(abs_lintemp,specs = pairwise ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(abs_lintemp,formula = habitat ~ meanmintempI_scaled, at = list(meanmintempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()


# ## Sample sizes
#
#
# c %>% filter(Age == "L") %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "WEBL")
#
# c %>% filter(Age == "L") %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "WEBL") %>% dplyr::select(attempt_id) %>% unique() %>% nrow()
#
# c %>% filter(Age == "L") %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "TRES")
#
# c %>% filter(Age == "L") %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "TRES") %>% dplyr::select(attempt_id) %>% unique() %>% nrow()
#
# c %>% filter(Age == "L") %>% dplyr::select(attempt_id,blood_num,Species,blood_num_mother) %>% distinct(.keep_all = TRUE) %>% filter(is.na(blood_num_mother))
#
#
# ## Sample sizes
#
#
# c %>% filter(Age %in% c("SY","ASY","AHY")) %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "WEBL")
#
# c %>% filter(Age %in% c("SY","ASY","AHY")) %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "WEBL") %>% dplyr::select(attempt_id) %>% unique() %>% nrow()
#
# c %>% filter(Age %in% c("SY","ASY","AHY")) %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "TRES")
#
# c %>% filter(Age %in% c("SY","ASY","AHY")) %>% dplyr::select(attempt_id,blood_num,Species) %>% distinct(.keep_all = TRUE) %>% filter(Species == "TRES") %>% dplyr::select(attempt_id) %>% unique() %>% nrow()
#
# ### hist of abs change
#
#
# c %>% mutate(abs_change = mean_s2-mean_s1) %>% filter(abs_change < 10) %>% pull(abs_change) %>% hist()

#### use prior day temp to predict cort instead


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_priordayt_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_priordayt_tres,"figures/int_tab_s1_priordayt_tres.html")

check_collinearity(s1_lintemp_noint)
s1_prior_day_tres <- s1_lintemp_addmin

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_priordayt_tres <- emmeans(s1_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_priordayt_tres,"figures/s1byhabitat_priordayt_tres.html")


## trends
##
(tress1_priordayt_trendmax <- emtrends(s1_lintemp_addmin,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_priordayt_trendmax,"figures/tress1_priordayt_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_priordayt <- s1_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_priordayt.html")

samp_s1_priordayt_tres <- s1_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordayt_tres %>% gt()


dat_text_s1_priordayt_tres <- data.frame(
  label = paste("N =",samp_s1_priordayt_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_priordayt_tres = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled)


mean_temp_s1_priordayt_tres <- mean(data_s1_priordayt_tres %>% pull(maxt_prior))
sd_temp_s1_priordayt_tres <- sd(data_s1_priordayt_tres %>% pull(maxt_prior))


temp_trans_s1_priordayt_tres <- trans_new("temp_trans_s1_priordayt_tres",
                                          transform = function(x){(x * sd_temp_s1_priordayt_tres) + mean_temp_s1_priordayt_tres},
                                          inverse = function(x){x})

(figs2_priordayt_tres <- ggpredict(s1_lintemp_addmin,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_priordayt_tres,
                       breaks = c((25-mean_temp_s1_priordayt_tres)/sd_temp_s1_priordayt_tres,
                                  (35-mean_temp_s1_priordayt_tres)/sd_temp_s1_priordayt_tres,
                                  (45-mean_temp_s1_priordayt_tres)/sd_temp_s1_priordayt_tres,
                                  (55-mean_temp_s1_priordayt_tres)/sd_temp_s1_priordayt_tres,
                                  (65-mean_temp_s1_priordayt_tres)/sd_temp_s1_priordayt_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bypriordaytxhab_TRES.png",plot =  figs2_priordayt_tres, width = 10, height = 6.6)

#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_priordayt_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_priordayt_tres,"figures/int_tab_s2_priordayt_tres.html")

check_collinearity(s2_lintemp_noint)
s2_prior_day_tres <- s2_lintemp

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_priordayt_tres <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_priordayt_tres,"figures/s2byhabitat_priordayt_tres.html")


## trends
##
(tress2_priordayt_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_priordayt_trendmax,"figures/tress2_priordayt_trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_priordayt <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_priordayt.html")

samp_s2_priordayt_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordayt_tres %>% gt()


dat_text_s2_priordayt_tres <- data.frame(
  label = paste("N =",samp_s2_priordayt_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_priordayt_tres = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled)


mean_temp_s2_priordayt_tres <- mean(data_s2_priordayt_tres %>% pull(maxt_prior))
sd_temp_s2_priordayt_tres <- sd(data_s2_priordayt_tres %>% pull(maxt_prior))


temp_trans_s2_priordayt_tres <- trans_new("temp_trans_s2_priordayt_tres",
                                          transform = function(x){(x * sd_temp_s2_priordayt_tres) + mean_temp_s2_priordayt_tres},
                                          inverse = function(x){x})

(fig_priordayt_tres_s2 <- ggpredict(s2_lintemp,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_priordayt_tres,
                       breaks = c((25-mean_temp_s2_priordayt_tres)/sd_temp_s2_priordayt_tres,
                                  (35-mean_temp_s2_priordayt_tres)/sd_temp_s2_priordayt_tres,
                                  (45-mean_temp_s2_priordayt_tres)/sd_temp_s2_priordayt_tres,
                                  (55-mean_temp_s2_priordayt_tres)/sd_temp_s2_priordayt_tres,
                                  (65-mean_temp_s2_priordayt_tres)/sd_temp_s2_priordayt_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dashed")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bypriordaytxhab_TRES.png",plot =  fig_priordayt_tres_s2, width = 10, height = 6.6)


#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_priordayt_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_priordayt_tres,"figures/int_tab_abs_priordayt_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_priordayt <- abs_lintemp

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_priordayt <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_priordayt.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abspriordayt_byhabitat_tres <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abspriordayt_byhabitat_tres,"figures/abspriordayt_byhabitat_tres.html")



summary(abs_lintemp)





data_tres_priordayt = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled)


mean_temp_tres_priordayt <- mean(data_tres_priordayt %>% pull(maxt_prior))
sd_temp_tres_priordayt <- sd(data_tres_priordayt %>% pull(maxt_prior))


temp_trans_tres_priordayt <- trans_new("temp_trans_tres_priordayt",
                                       transform = function(x){(x * sd_temp_tres_priordayt) + mean_temp_tres_priordayt},
                                       inverse = function(x){x})

(fig4_tres_priordayt <- ggpredict(abs_lintemp,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres_priordayt,
                       breaks = c((20-mean_temp_tres_priordayt)/sd_temp_tres_priordayt,
                                  (25-mean_temp_tres_priordayt)/sd_temp_tres_priordayt,
                                  (30-mean_temp_tres_priordayt)/sd_temp_tres_priordayt,
                                  (35-mean_temp_tres_priordayt)/sd_temp_tres_priordayt,
                                  (40-mean_temp_tres_priordayt)/sd_temp_tres_priordayt),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_priordayt_xhab_TRES.png",plot =  fig4_tres_priordayt, width = 10, height = 6.6)




(tresabs_priordayt_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_priordayt_trendmax,"figures/tresabs_priordayt_trendmax.html")

#### use prior day heat index to predict cort instead


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_priordaymaxhhi_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_priordaymaxhhi_tres,"figures/int_tab_s1_priordaymaxhhi_tres.html")

s1_priordaymaxhhi_tres <- s1_lintemp_addmax
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_priordaymaxhhi_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_priordaymaxhhi_tres,"figures/s1byhabitat_priordaymaxhhi_tres.html")


## trends
##
(tress1_priordaymaxhhi_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_priordaymaxhhi_trendmax,"figures/tress1_priordaymaxhhi_trendmax.html")


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_priordaymaxhhi <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_priordaymaxhhi.html")

samp_s1_priordaymaxhhi_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordaymaxhhi_tres %>% gt()


dat_text_s1_priordaymaxhhi_tres <- data.frame(
  label = paste("N =",samp_s1_priordaymaxhhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_priordaymaxhhi_tres =  dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_s1_priordaymaxhhi_tres <- mean(data_s1_priordaymaxhhi_tres %>% pull(maxhi_prior))
sd_temp_s1_priordaymaxhhi_tres <- sd(data_s1_priordaymaxhhi_tres %>% pull(maxhi_prior))


temp_trans_s1_priordaymaxhhi_tres <- trans_new("temp_trans_s1_priordaymaxhhi_tres",
                                               transform = function(x){(x * sd_temp_s1_priordaymaxhhi_tres) + mean_temp_s1_priordaymaxhhi_tres},
                                               inverse = function(x){x})

(figs2_priordaymaxhhi_tres <- ggpredict(s1_lintemp_addmax,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior day (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_priordaymaxhhi_tres,
                       breaks = c((20-mean_temp_s1_priordaymaxhhi_tres)/sd_temp_s1_priordaymaxhhi_tres,
                                  (40-mean_temp_s1_priordaymaxhhi_tres)/sd_temp_s1_priordaymaxhhi_tres,
                                  (60-mean_temp_s1_priordaymaxhhi_tres)/sd_temp_s1_priordaymaxhhi_tres,
                                  (80-mean_temp_s1_priordaymaxhhi_tres)/sd_temp_s1_priordaymaxhhi_tres,
                                  (100-mean_temp_s1_priordaymaxhhi_tres)/sd_temp_s1_priordaymaxhhi_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("solid","solid","solid","solid")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_priordaymaxhhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bypriordaymaxhhixhab_TRES.png",plot =  figs2_priordaymaxhhi_tres, width = 10, height = 6.6)


#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_priordaymaxhhi_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_priordaymaxhhi_tres,"figures/int_tab_s2_priordaymaxhhi_tres.html")

s2_priordaymaxhhi_tres <- s2_lintemp
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_priordaymaxhhi_tres <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_priordaymaxhhi_tres,"figures/s2byhabitat_priordaymaxhhi_tres.html")


## trends
##
(tress2_priordaymaxhhi_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_priordaymaxhhi_trendmax,"figures/tress2_priordaymaxhhi_trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_priordaymaxhhi <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_priordaymaxhhi.html")

samp_s2_priordaymaxhhi_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordaymaxhhi_tres %>% gt()


dat_text_s2_priordaymaxhhi_tres <- data.frame(
  label = paste("N =",samp_s2_priordaymaxhhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_priordaymaxhhi_tres =  dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_s2_priordaymaxhhi_tres <- mean(data_s2_priordaymaxhhi_tres %>% pull(maxhi_prior))
sd_temp_s2_priordaymaxhhi_tres <- sd(data_s2_priordaymaxhhi_tres %>% pull(maxhi_prior))


temp_trans_s2_priordaymaxhhi_tres <- trans_new("temp_trans_s2_priordaymaxhhi_tres",
                                               transform = function(x){(x * sd_temp_s2_priordaymaxhhi_tres) + mean_temp_s2_priordaymaxhhi_tres},
                                               inverse = function(x){x})

(fig_priordaymaxhhi_tres_s2 <- ggpredict(s2_lintemp,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior day (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_priordaymaxhhi_tres,
                       breaks = c((20-mean_temp_s2_priordaymaxhhi_tres)/sd_temp_s2_priordaymaxhhi_tres,
                                  (40-mean_temp_s2_priordaymaxhhi_tres)/sd_temp_s2_priordaymaxhhi_tres,
                                  (60-mean_temp_s2_priordaymaxhhi_tres)/sd_temp_s2_priordaymaxhhi_tres,
                                  (80-mean_temp_s2_priordaymaxhhi_tres)/sd_temp_s2_priordaymaxhhi_tres,
                                  (100-mean_temp_s2_priordaymaxhhi_tres)/sd_temp_s2_priordaymaxhhi_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","solid","dotted","dashed")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_priordaymaxhhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bypriordaymaxhhixhab_TRES.png",plot =  fig_priordaymaxhhi_tres_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_priordaymaxhhi_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_priordaymaxhhi_tres,"figures/int_tab_abs_priordaymaxhhi_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_priordaymaxhhi <- abs_lintemp

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_priordaymaxhhi <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_priordaymaxhhi.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_priordaymaxhhi_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abspriordaymaxhhi_byhabitat_tres <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abspriordaymaxhhi_byhabitat_tres,"figures/abspriordaymaxhhi_byhabitat_tres.html")



summary(abs_lintemp)




data_tres_priordaymaxhhi = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(mint_prior)) %>%
  mutate(across(c(gweight,maxhi_prior,mint_prior,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_tres_priordaymaxhhi <- mean(data_tres_priordaymaxhhi %>% pull(maxhi_prior))
sd_temp_tres_priordaymaxhhi <- sd(data_tres_priordaymaxhhi %>% pull(maxhi_prior))


temp_trans_tres_priordaymaxhhi <- trans_new("temp_trans_tres_priordaymaxhhi",
                                            transform = function(x){(x * sd_temp_tres_priordaymaxhhi) + mean_temp_tres_priordaymaxhhi},
                                            inverse = function(x){x})

(fig4_tres_priordaymaxhhi <- ggpredict(abs_lintemp,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior day (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres_priordaymaxhhi,
                       breaks = c((20-mean_temp_tres_priordaymaxhhi)/sd_temp_tres_priordaymaxhhi,
                                  (40-mean_temp_tres_priordaymaxhhi)/sd_temp_tres_priordaymaxhhi,
                                  (60-mean_temp_tres_priordaymaxhhi)/sd_temp_tres_priordaymaxhhi,
                                  (80-mean_temp_tres_priordaymaxhhi)/sd_temp_tres_priordaymaxhhi,
                                  (100-mean_temp_tres_priordaymaxhhi)/sd_temp_tres_priordaymaxhhi),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_priordaymaxhhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_priordaymaxhhi_xhab_TRES.png",plot =  fig4_tres_priordayt, width = 10, height = 6.6)




(tresabs_priordaymaxhhi_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_priordaymaxhhi_trendmax,"figures/tresabs_priordaymaxhhi_trendmax.html")


## Use prior week heat index to test against cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_weekhi_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_weekhi_tres,"figures/int_tab_s1_weekhi_tres.html")

s1_weekhi_tres <- s1_lintemp_noint
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_weekhi_tres <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_weekhi_tres,"figures/s1byhabitat_weekhi_tres.html")


## trends
##
(tress1_weekhi_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_weekhi_trendmax,"figures/tress1_weekhi_trendmax.html")


summary(s1_lintemp_noint)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_weekhi <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_weekhi.html")

samp_s1_weekhi_tres <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_weekhi_tres %>% gt()


dat_text_s1_weekhi_tres <- data.frame(
  label = paste("N =",samp_s1_weekhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_weekhi_tres = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_s1_weekhi_tres <- mean(data_s1_weekhi_tres %>% pull(maxhi_week))
sd_temp_s1_weekhi_tres <- sd(data_s1_weekhi_tres %>% pull(maxhi_week))


temp_trans_s1_weekhi_tres <- trans_new("temp_trans_s1_weekhi_tres",
                                       transform = function(x){(x * sd_temp_s1_weekhi_tres) + mean_temp_s1_weekhi_tres},
                                       inverse = function(x){x})

(figs2_weekhi_tres <- ggpredict(s1_lintemp_noint,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_weekhi_tres,
                       breaks = c((20-mean_temp_s1_weekhi_tres)/sd_temp_s1_weekhi_tres,
                                  (30-mean_temp_s1_weekhi_tres)/sd_temp_s1_weekhi_tres,
                                  (40-mean_temp_s1_weekhi_tres)/sd_temp_s1_weekhi_tres,
                                  (50-mean_temp_s1_weekhi_tres)/sd_temp_s1_weekhi_tres,
                                  (60-mean_temp_s1_weekhi_tres)/sd_temp_s1_weekhi_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_weekhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1byweekhixhab_TRES.png",plot =  figs2_weekhi_tres, width = 10, height = 6.6)


## s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_weekhi_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_weekhi_tres,"figures/int_tab_s2_weekhi_tres.html")

s2_weekhi_tres <- s2_lintemp_addmin
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_weekhi_tres <- emmeans(s2_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_weekhi_tres,"figures/s2byhabitat_weekhi_tres.html")


## trends
##
(tress2_weekhi_trendmax <- emtrends(s2_lintemp_addmin,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_weekhi_trendmax,"figures/tress2_weekhi_trendmax.html")


summary(s2_lintemp_addmin)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_weekhi <- s2_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_weekhi.html")

samp_s2_weekhi_tres <- s2_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_weekhi_tres %>% gt()


dat_text_s2_weekhi_tres <- data.frame(
  label = paste("N =",samp_s2_weekhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_weekhi_tres = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_s2_weekhi_tres <- mean(data_s2_weekhi_tres %>% pull(maxhi_week))
sd_temp_s2_weekhi_tres <- sd(data_s2_weekhi_tres %>% pull(maxhi_week))


temp_trans_s2_weekhi_tres <- trans_new("temp_trans_s2_weekhi_tres",
                                       transform = function(x){(x * sd_temp_s2_weekhi_tres) + mean_temp_s2_weekhi_tres},
                                       inverse = function(x){x})

(fig_weekhi_tres_s2 <- ggpredict(s2_lintemp_addmin,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_weekhi_tres,
                       breaks = c((20-mean_temp_s2_weekhi_tres)/sd_temp_s2_weekhi_tres,
                                  (30-mean_temp_s2_weekhi_tres)/sd_temp_s2_weekhi_tres,
                                  (40-mean_temp_s2_weekhi_tres)/sd_temp_s2_weekhi_tres,
                                  (50-mean_temp_s2_weekhi_tres)/sd_temp_s2_weekhi_tres,
                                  (60-mean_temp_s2_weekhi_tres)/sd_temp_s2_weekhi_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_weekhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2byweekhixhab_TRES.png",plot =  fig_weekhi_tres_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_weekhi_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_weekhi_tres,"figures/int_tab_abs_weekhi_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_weekhi <- abs_lintemp

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_weekhi <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_weekhi.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres_weekhi <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(absweekhi_byhabitat_tres <- emmeans(abs_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(absweekhi_byhabitat_tres,"figures/absweekhi_byhabitat_tres.html")



summary(abs_lintemp_noint)





data_tres_weekhi = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_tres_weekhi <- mean(data_tres_weekhi %>% pull(maxhi_week))
sd_temp_tres_weekhi <- sd(data_tres_weekhi %>% pull(maxhi_week))


temp_trans_tres_weekhi <- trans_new("temp_trans_tres_weekhi",
                                    transform = function(x){(x * sd_temp_tres_weekhi) + mean_temp_tres_weekhi},
                                    inverse = function(x){x})

(fig4_tres_weekhi <- ggpredict(abs_lintemp,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Max heat index over prior week (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres_weekhi,
                       breaks = c((30-mean_temp_tres_weekhi)/sd_temp_tres_weekhi,
                                  (40-mean_temp_tres_weekhi)/sd_temp_tres_weekhi,
                                  (50-mean_temp_tres_weekhi)/sd_temp_tres_weekhi,
                                  (60-mean_temp_tres_weekhi)/sd_temp_tres_weekhi,
                                  (70-mean_temp_tres_weekhi)/sd_temp_tres_weekhi),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres_weekhi, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_weekhi_xhab_TRES.png",plot =  fig4_tres_weekhi, width = 10, height = 6.6)




(tresabs_weekhi_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_weekhi_trendmax,"figures/tresabs_weekhi_trendmax.html")



#### use cumulative prior day hi to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumhiday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumhiday_tres,"figures/int_tab_s1_cumhiday_tres.html")

s1_cumhiday_tres <- s1_lintemp_addmax
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_cumhiday_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumhiday_tres,"figures/s1byhabitat_cumhiday_tres.html")


## trends
##
(tress1_cumhiday_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_cumhiday_trendmax,"figures/tress1_cumhiday_trendmax.html")


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_cumhiday <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_cumhiday.html")

samp_s1_cumhiday_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiday_tres %>% gt()


dat_text_s1_cumhiday_tres <- data.frame(
  label = paste("N =",samp_s1_cumhiday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumhiday_tres = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled)


mean_temp_s1_cumhiday_tres <- mean(data_s1_cumhiday_tres %>% pull(hihours_over_30hi_priorday))
sd_temp_s1_cumhiday_tres <- sd(data_s1_cumhiday_tres %>% pull(hihours_over_30hi_priorday))


temp_trans_s1_cumhiday_tres <- trans_new("temp_trans_s1_cumhiday_tres",
                                         transform = function(x){(x * sd_temp_s1_cumhiday_tres) + mean_temp_s1_cumhiday_tres},
                                         inverse = function(x){x})

(figs2_cumhiday_tres <- ggpredict(s1_lintemp_addmax,terms = c("hihours_over_30hi_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index hours >25\u00b0C over prior day") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumhiday_tres,
                       breaks = c((0-mean_temp_s1_cumhiday_tres)/sd_temp_s1_cumhiday_tres,
                                  (200-mean_temp_s1_cumhiday_tres)/sd_temp_s1_cumhiday_tres,
                                  (400-mean_temp_s1_cumhiday_tres)/sd_temp_s1_cumhiday_tres,
                                  (600-mean_temp_s1_cumhiday_tres)/sd_temp_s1_cumhiday_tres,
                                  (800-mean_temp_s1_cumhiday_tres)/sd_temp_s1_cumhiday_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumhiday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumhidayxhab_TRES.png",plot =  figs2_cumhiday_tres, width = 10, height = 6.6)

#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumhiday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumhiday_tres,"figures/int_tab_s2_cumhiday_tres.html")

s2_cumhiday_tres <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_cumhiday_tres <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumhiday_tres,"figures/s2byhabitat_cumhiday_tres.html")


## trends
##
(tress2_cumhiday_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_cumhiday_trendmax,"figures/tress2_cumhiday_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_cumhiday <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_cumhiday.html")

samp_s2_cumhiday_tres <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiday_tres %>% gt()


dat_text_s2_cumhiday_tres <- data.frame(
  label = paste("N =",samp_s2_cumhiday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumhiday_tres = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled)


mean_temp_s2_cumhiday_tres <- mean(data_s2_cumhiday_tres %>% pull(hihours_over_30hi_priorday))
sd_temp_s2_cumhiday_tres <- sd(data_s2_cumhiday_tres %>% pull(hihours_over_30hi_priorday))


temp_trans_s2_cumhiday_tres <- trans_new("temp_trans_s2_cumhiday_tres",
                                         transform = function(x){(x * sd_temp_s2_cumhiday_tres) + mean_temp_s2_cumhiday_tres},
                                         inverse = function(x){x})

(fig_cumhiday_tres_s2 <- ggpredict(s2_lintemp_noint,terms = c("hihours_over_30hi_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index hours >25\u00b0C over prior day") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumhiday_tres,
                       breaks = c((0-mean_temp_s2_cumhiday_tres)/sd_temp_s2_cumhiday_tres,
                                  (200-mean_temp_s2_cumhiday_tres)/sd_temp_s2_cumhiday_tres,
                                  (400-mean_temp_s2_cumhiday_tres)/sd_temp_s2_cumhiday_tres,
                                  (600-mean_temp_s2_cumhiday_tres)/sd_temp_s2_cumhiday_tres,
                                  (800-mean_temp_s2_cumhiday_tres)/sd_temp_s2_cumhiday_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumhiday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumhidayxhab_TRES.png",plot =  fig_cumhiday_tres_s2, width = 10, height = 6.6)


#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumhiday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumhiday_tres,"figures/int_tab_abs_cumhiday_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_cumhiday <- abs_lintemp_noint

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_cumhiday <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_cumhiday.html")

samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres_cumhiday <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumhiday_byhabitat_tres <- emmeans(abs_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumhiday_byhabitat_tres,"figures/abscumhiday_byhabitat_tres.html")



summary(abs_lintemp_noint)





data_tres_cumhiday = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorday_scaled_sq = hihours_over_30hi_priorday_scaled * hihours_over_30hi_priorday_scaled)


mean_temp_tres_cumhiday <- mean(data_tres_cumhiday %>% pull(hihours_over_30hi_priorday))
sd_temp_tres_cumhiday <- sd(data_tres_cumhiday %>% pull(hihours_over_30hi_priorday))


temp_trans_tres_cumhiday <- trans_new("temp_trans_tres_cumhiday",
                                      transform = function(x){(x * sd_temp_tres_cumhiday) + mean_temp_tres_cumhiday},
                                      inverse = function(x){x})

(fig4_tres_cumhiday <- ggpredict(abs_lintemp_noint,terms = c("hihours_over_30hi_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index hours >25\u00b0C over prior day") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres_cumhiday,
                       breaks = c((0-mean_temp_tres_cumhiday)/sd_temp_tres_cumhiday,
                                  (200-mean_temp_tres_cumhiday)/sd_temp_tres_cumhiday,
                                  (400-mean_temp_tres_cumhiday)/sd_temp_tres_cumhiday,
                                  (600-mean_temp_tres_cumhiday)/sd_temp_tres_cumhiday,
                                  (800-mean_temp_tres_cumhiday)/sd_temp_tres_cumhiday),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dashed","Grassland" = "dashed","Row crop" = "dashed")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres_cumhiday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumhiday_xhab_TRES.png",plot =  fig4_tres_cumhiday, width = 10, height = 6.6)




(tresabs_cumhiday_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_cumhiday_trendmax,"figures/tresabs_cumhiday_trendmax.html")


#### use cumulative prior week hi to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumhiweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumhiweek_tres,"figures/int_tab_s1_cumhiweek_tres.html")

s1_cumhiweek_tres <- s1_lintemp_addmax
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_cumhiweek_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumhiweek_tres,"figures/s1byhabitat_cumhiweek_tres.html")


## trends
##
(tress1_cumhiweek_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_cumhiweek_trendmax,"figures/tress1_cumhiweek_trendmax.html")


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_cumhiweek <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_cumhiweek.html")

samp_s1_cumhiweek_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiweek_tres %>% gt()


dat_text_s1_cumhiweek_tres <- data.frame(
  label = paste("N =",samp_s1_cumhiweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumhiweek_tres = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_temp_s1_cumhiweek_tres <- mean(data_s1_cumhiweek_tres %>% pull(hihours_over_30hi_priorweek))
sd_temp_s1_cumhiweek_tres <- sd(data_s1_cumhiweek_tres %>% pull(hihours_over_30hi_priorweek))


temp_trans_s1_cumhiweek_tres <- trans_new("temp_trans_s1_cumhiweek_tres",
                                          transform = function(x){(x * sd_temp_s1_cumhiweek_tres) + mean_temp_s1_cumhiweek_tres},
                                          inverse = function(x){x})

(figs2_cumhiweek_tres <- ggpredict(s1_lintemp_addmax,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index hours >25\u00b0C over prior week") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumhiweek_tres,
                       breaks = c((0-mean_temp_s1_cumhiweek_tres)/sd_temp_s1_cumhiweek_tres,
                                  (1000-mean_temp_s1_cumhiweek_tres)/sd_temp_s1_cumhiweek_tres,
                                  (2000-mean_temp_s1_cumhiweek_tres)/sd_temp_s1_cumhiweek_tres,
                                  (3000-mean_temp_s1_cumhiweek_tres)/sd_temp_s1_cumhiweek_tres,
                                  (4000-mean_temp_s1_cumhiweek_tres)/sd_temp_s1_cumhiweek_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumhiweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumhiweekxhab_TRES.png",plot =  figs2_cumhiweek_tres, width = 10, height = 6.6)


#### s2


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumhiweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumhiweek_tres,"figures/int_tab_s2_cumhiweek_tres.html")

s2_cumhiweek_tres <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_cumhiweek_tres <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumhiweek_tres,"figures/s2byhabitat_cumhiweek_tres.html")


## trends
##
(tress2_cumhiweek_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_cumhiweek_trendmax,"figures/tress2_cumhiweek_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_cumhiweek <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_cumhiweek.html")

samp_s2_cumhiweek_tres <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiweek_tres %>% gt()


dat_text_s2_cumhiweek_tres <- data.frame(
  label = paste("N =",samp_s2_cumhiweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumhiweek_tres = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_temp_s2_cumhiweek_tres <- mean(data_s2_cumhiweek_tres %>% pull(hihours_over_30hi_priorweek))
sd_temp_s2_cumhiweek_tres <- sd(data_s2_cumhiweek_tres %>% pull(hihours_over_30hi_priorweek))


temp_trans_s2_cumhiweek_tres <- trans_new("temp_trans_s2_cumhiweek_tres",
                                          transform = function(x){(x * sd_temp_s2_cumhiweek_tres) + mean_temp_s2_cumhiweek_tres},
                                          inverse = function(x){x})

(fig_cumhiweek_tres_s2 <- ggpredict(s2_lintemp_noint,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index hours >25\u00b0C over prior week") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumhiweek_tres,
                       breaks = c((0-mean_temp_s2_cumhiweek_tres)/sd_temp_s2_cumhiweek_tres,
                                  (1000-mean_temp_s2_cumhiweek_tres)/sd_temp_s2_cumhiweek_tres,
                                  (2000-mean_temp_s2_cumhiweek_tres)/sd_temp_s2_cumhiweek_tres,
                                  (3000-mean_temp_s2_cumhiweek_tres)/sd_temp_s2_cumhiweek_tres,
                                  (4000-mean_temp_s2_cumhiweek_tres)/sd_temp_s2_cumhiweek_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumhiweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumhiweekxhab_TRES.png",plot =  fig_cumhiweek_tres_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumhiweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumhiweek_tres,"figures/int_tab_abs_cumhiweek_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_cumhiweek <- abs_lintemp

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_cumhiweek <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_cumhiweek.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_abs_text_tres_cumhiweek <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumhiweek_byhabitat_tres <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumhiweek_byhabitat_tres,"figures/abscumhiweek_byhabitat_tres.html")



summary(abs_lintemp)





data_abs_tres_cumhiweek = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_abs_temp_tres_cumhiweek <- mean(data_abs_tres_cumhiweek %>% pull(hihours_over_30hi_priorweek))
sd_abs_temp_tres_cumhiweek <- sd(data_abs_tres_cumhiweek %>% pull(hihours_over_30hi_priorweek))


temp_abs_trans_tres_cumhiweek <- trans_new("temp_trans_tres_cumhiweek",
                                           transform = function(x){(x * sd_abs_temp_tres_cumhiweek) + mean_abs_temp_tres_cumhiweek},
                                           inverse = function(x){x})

(fig4_tres_cumhiweek <- ggpredict(abs_lintemp,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index hours >25\u00b0C over prior week") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_abs_trans_tres_cumhiweek,
                       breaks = c((0-mean_abs_temp_tres_cumhiweek)/sd_abs_temp_tres_cumhiweek,
                                  (1000-mean_abs_temp_tres_cumhiweek)/sd_abs_temp_tres_cumhiweek,
                                  (2000-mean_abs_temp_tres_cumhiweek)/sd_abs_temp_tres_cumhiweek,
                                  (3000-mean_abs_temp_tres_cumhiweek)/sd_abs_temp_tres_cumhiweek,
                                  (4000-mean_abs_temp_tres_cumhiweek)/sd_abs_temp_tres_cumhiweek),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "solid","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_abs_text_tres_cumhiweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumhiweek_xhab_TRES.png",plot =  fig4_tres_cumhiweek, width = 10, height = 6.6)




(tresabs_cumhiweek_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_cumhiweek_trendmax,"figures/tresabs_cumhiweek_trendmax.html")

#### use cumulative prior day temp to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumdegreeday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumdegreeday_tres,"figures/int_tab_s1_cumdegreeday_tres.html")

s1_cumdegreeday_tres <- s1_lintemp_addmax
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeday_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumdegreeday_tres,"figures/s1byhabitat_cumdegreeday_tres.html")


## trends
##
(tress1_cumdegreeday_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_cumdegreeday_trendmax,"figures/tress1_cumdegreeday_trendmax.html")


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_cumdegreeday <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_cumdegreeday.html")

samp_s1_cumdegreeday_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeday_tres %>% gt()


dat_text_s1_cumdegreeday_tres <- data.frame(
  label = paste("N =",samp_s1_cumdegreeday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumdegreeday_tres = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled)


mean_temp_s1_cumdegreeday_tres <- mean(data_s1_cumdegreeday_tres %>% pull(degreehours_over_30C_priorday))
sd_temp_s1_cumdegreeday_tres <- sd(data_s1_cumdegreeday_tres %>% pull(degreehours_over_30C_priorday))


temp_trans_s1_cumdegreeday_tres <- trans_new("temp_trans_s1_cumdegreeday_tres",
                                             transform = function(x){(x * sd_temp_s1_cumdegreeday_tres) + mean_temp_s1_cumdegreeday_tres},
                                             inverse = function(x){x})

(figs2_cumdegreeday_tres <- ggpredict(s1_lintemp_addmax,terms = c("degreehours_over_30C_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior day") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumdegreeday_tres,
                       breaks = c((0-mean_temp_s1_cumdegreeday_tres)/sd_temp_s1_cumdegreeday_tres,
                                  (100-mean_temp_s1_cumdegreeday_tres)/sd_temp_s1_cumdegreeday_tres,
                                  (200-mean_temp_s1_cumdegreeday_tres)/sd_temp_s1_cumdegreeday_tres,
                                  (300-mean_temp_s1_cumdegreeday_tres)/sd_temp_s1_cumdegreeday_tres,
                                  (400-mean_temp_s1_cumdegreeday_tres)/sd_temp_s1_cumdegreeday_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumdegreeday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumdegreedayxhab_TRES.png",plot =  figs2_cumdegreeday_tres, width = 10, height = 6.6)


#### use cumulative prior day temp to predict cort


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumdegreeday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumdegreeday_tres,"figures/int_tab_s2_cumdegreeday_tres.html")

s2_cumdegreeday_tres <- s2_lintemp
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeday_tres <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumdegreeday_tres,"figures/s2byhabitat_cumdegreeday_tres.html")


## trends
##
(tress2_cumdegreeday_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_cumdegreeday_trendmax,"figures/tress2_cumdegreeday_trendmax.html")


summary(s2_lintemp)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_cumdegreeday <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_cumdegreeday.html")

samp_s2_cumdegreeday_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeday_tres %>% gt()


dat_text_s2_cumdegreeday_tres <- data.frame(
  label = paste("N =",samp_s2_cumdegreeday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumdegreeday_tres = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled)


mean_temp_s2_cumdegreeday_tres <- mean(data_s2_cumdegreeday_tres %>% pull(degreehours_over_30C_priorday))
sd_temp_s2_cumdegreeday_tres <- sd(data_s2_cumdegreeday_tres %>% pull(degreehours_over_30C_priorday))


temp_trans_s2_cumdegreeday_tres <- trans_new("temp_trans_s2_cumdegreeday_tres",
                                             transform = function(x){(x * sd_temp_s2_cumdegreeday_tres) + mean_temp_s2_cumdegreeday_tres},
                                             inverse = function(x){x})

(fig_cumdegreeday_tres_s2 <- ggpredict(s2_lintemp,terms = c("degreehours_over_30C_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior day") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumdegreeday_tres,
                       breaks = c((0-mean_temp_s2_cumdegreeday_tres)/sd_temp_s2_cumdegreeday_tres,
                                  (100-mean_temp_s2_cumdegreeday_tres)/sd_temp_s2_cumdegreeday_tres,
                                  (200-mean_temp_s2_cumdegreeday_tres)/sd_temp_s2_cumdegreeday_tres,
                                  (300-mean_temp_s2_cumdegreeday_tres)/sd_temp_s2_cumdegreeday_tres,
                                  (400-mean_temp_s2_cumdegreeday_tres)/sd_temp_s2_cumdegreeday_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumdegreeday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumdegreedayxhab_TRES.png",plot =  fig_cumdegreeday_tres_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumdegreeday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumdegreeday_tres,"figures/int_tab_abs_cumdegreeday_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_cumdegreeday <- abs_lintemp

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_cumdegreeday <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_cumdegreeday.html")

samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres_cumdegreeday <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumdegreeday_byhabitat_tres <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumdegreeday_byhabitat_tres,"figures/abscumdegreeday_byhabitat_tres.html")



summary(abs_lintemp)





data_tres_cumdegreeday = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorday),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorday,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorday_scaled_sq = degreehours_over_30C_priorday_scaled * degreehours_over_30C_priorday_scaled)


mean_temp_tres_cumdegreeday <- mean(data_tres_cumdegreeday %>% pull(degreehours_over_30C_priorday))
sd_temp_tres_cumdegreeday <- sd(data_tres_cumdegreeday %>% pull(degreehours_over_30C_priorday))


temp_trans_tres_cumdegreeday <- trans_new("temp_trans_tres_cumdegreeday",
                                          transform = function(x){(x * sd_temp_tres_cumdegreeday) + mean_temp_tres_cumdegreeday},
                                          inverse = function(x){x})

(fig4_tres_cumdegreeday <- ggpredict(abs_lintemp,terms = c("degreehours_over_30C_priorday_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior day") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres_cumdegreeday,
                       breaks = c((0-mean_temp_tres_cumdegreeday)/sd_temp_tres_cumdegreeday,
                                  (100-mean_temp_tres_cumdegreeday)/sd_temp_tres_cumdegreeday,
                                  (200-mean_temp_tres_cumdegreeday)/sd_temp_tres_cumdegreeday,
                                  (300-mean_temp_tres_cumdegreeday)/sd_temp_tres_cumdegreeday,
                                  (400-mean_temp_tres_cumdegreeday)/sd_temp_tres_cumdegreeday),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres_cumdegreeday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumdegreeday_xhab_TRES.png",plot =  fig4_tres_cumdegreeday, width = 10, height = 6.6)




(tresabs_cumdegreeday_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_cumdegreeday_trendmax,"figures/tresabs_cumdegreeday_trendmax.html")

#### use cumulative prior week temp to predict cort


s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(s1_lintemp,s1_lintemp_addmin,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s1_lintemp,s1_lintemp_addmax,s1_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s1_cumdegreeweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s1_cumdegreeweek_tres,"figures/int_tab_s1_cumdegreeweek_tres.html")

s1_cumdegreeweek_tres <- s1_lintemp_addmax
check_collinearity(s1_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeweek_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s1byhabitat_cumdegreeweek_tres,"figures/s1byhabitat_cumdegreeweek_tres.html")


## trends
##
(tress1_cumdegreeweek_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress1_cumdegreeweek_trendmax,"figures/tress1_cumdegreeweek_trendmax.html")


summary(s1_lintemp_addmax)
check_collinearity(s1_lintemp_noint)


### Sample sizes:


ss_year_tres_s1_cumdegreeweek <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s1_cumdegreeweek.html")

samp_s1_cumdegreeweek_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeweek_tres %>% gt()


dat_text_s1_cumdegreeweek_tres <- data.frame(
  label = paste("N =",samp_s1_cumdegreeweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s1_cumdegreeweek_tres = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_s1_cumdegreeweek_tres <- mean(data_s1_cumdegreeweek_tres %>% pull(degreehours_over_30C_priorweek))
sd_temp_s1_cumdegreeweek_tres <- sd(data_s1_cumdegreeweek_tres %>% pull(degreehours_over_30C_priorweek))


temp_trans_s1_cumdegreeweek_tres <- trans_new("temp_trans_s1_cumdegreeweek_tres",
                                              transform = function(x){(x * sd_temp_s1_cumdegreeweek_tres) + mean_temp_s1_cumdegreeweek_tres},
                                              inverse = function(x){x})

(figs2_cumdegreeweek_tres <- ggpredict(s1_lintemp_addmax,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s1_cumdegreeweek_tres,
                       breaks = c((0-mean_temp_s1_cumdegreeweek_tres)/sd_temp_s1_cumdegreeweek_tres,
                                  (1000-mean_temp_s1_cumdegreeweek_tres)/sd_temp_s1_cumdegreeweek_tres,
                                  (2000-mean_temp_s1_cumdegreeweek_tres)/sd_temp_s1_cumdegreeweek_tres,
                                  (3000-mean_temp_s1_cumdegreeweek_tres)/sd_temp_s1_cumdegreeweek_tres,
                                  (4000-mean_temp_s1_cumdegreeweek_tres)/sd_temp_s1_cumdegreeweek_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_cumdegreeweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s1bycumdegreeweekxhab_TRES.png",plot =  figs2_cumdegreeweek_tres, width = 10, height = 6.6)


#### use cumulative prior week temp to predict cort


s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(s2_lintemp,s2_lintemp_addmin,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(s2_lintemp,s2_lintemp_addmax,s2_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_s2_cumdegreeweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_s2_cumdegreeweek_tres,"figures/int_tab_s2_cumdegreeweek_tres.html")

s2_cumdegreeweek_tres <- s2_lintemp_noint
check_collinearity(s2_lintemp_noint)

# Looks like there is not an interactive effect of temp and habitat for cort in nestling TRES.

## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeweek_tres <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(s2byhabitat_cumdegreeweek_tres,"figures/s2byhabitat_cumdegreeweek_tres.html")


## trends
##
(tress2_cumdegreeweek_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tress2_cumdegreeweek_trendmax,"figures/tress2_cumdegreeweek_trendmax.html")


summary(s2_lintemp_noint)
check_collinearity(s2_lintemp_noint)


### Sample sizes:


ss_year_tres_s2_cumdegreeweek <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_s2_cumdegreeweek.html")

samp_s2_cumdegreeweek_tres <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeweek_tres %>% gt()


dat_text_s2_cumdegreeweek_tres <- data.frame(
  label = paste("N =",samp_s2_cumdegreeweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_s2_cumdegreeweek_tres = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_s2_cumdegreeweek_tres <- mean(data_s2_cumdegreeweek_tres %>% pull(degreehours_over_30C_priorweek))
sd_temp_s2_cumdegreeweek_tres <- sd(data_s2_cumdegreeweek_tres %>% pull(degreehours_over_30C_priorweek))


temp_trans_s2_cumdegreeweek_tres <- trans_new("temp_trans_s2_cumdegreeweek_tres",
                                              transform = function(x){(x * sd_temp_s2_cumdegreeweek_tres) + mean_temp_s2_cumdegreeweek_tres},
                                              inverse = function(x){x})

(fig_cumdegreeweek_tres_s2 <- ggpredict(s2_lintemp_noint,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_s2_cumdegreeweek_tres,
                       breaks = c((0-mean_temp_s2_cumdegreeweek_tres)/sd_temp_s2_cumdegreeweek_tres,
                                  (1000-mean_temp_s2_cumdegreeweek_tres)/sd_temp_s2_cumdegreeweek_tres,
                                  (2000-mean_temp_s2_cumdegreeweek_tres)/sd_temp_s2_cumdegreeweek_tres,
                                  (3000-mean_temp_s2_cumdegreeweek_tres)/sd_temp_s2_cumdegreeweek_tres,
                                  (4000-mean_temp_s2_cumdegreeweek_tres)/sd_temp_s2_cumdegreeweek_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_cumdegreeweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/s2bycumdegreeweekxhab_TRES.png",plot =  fig_cumdegreeweek_tres_s2, width = 10, height = 6.6)

#### abs_change_cort


abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(abs_lintemp,abs_lintemp_addmin,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(abs_lintemp,abs_lintemp_addmax,abs_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_abs_cumdegreeweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 3)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_abs_cumdegreeweek_tres,"figures/int_tab_abs_cumdegreeweek_tres.html")

check_collinearity(abs_lintemp_noint)

abs_tres_cumdegreeweek <- abs_lintemp_noint

# For abs_change in TRES, there are no interactions

### sample size


ss_year_tres_abs_cumdegreeweek <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_tres_abs_cumdegreeweek.html")

samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres_cumdegreeweek <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


## Emmeans to check for effect of habitat


(abscumdegreeweek_byhabitat_tres <- emmeans(abs_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(abscumdegreeweek_byhabitat_tres,"figures/abscumdegreeweek_byhabitat_tres.html")



summary(abs_lintemp_noint)





data_tres_cumdegreeweek = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_tres_cumdegreeweek <- mean(data_tres_cumdegreeweek %>% pull(degreehours_over_30C_priorweek))
sd_temp_tres_cumdegreeweek <- sd(data_tres_cumdegreeweek %>% pull(degreehours_over_30C_priorweek))


temp_trans_tres_cumdegreeweek <- trans_new("temp_trans_tres_cumdegreeweek",
                                           transform = function(x){(x * sd_temp_tres_cumdegreeweek) + mean_temp_tres_cumdegreeweek},
                                           inverse = function(x){x})

(fig4_tres_cumdegreeweek <- ggpredict(abs_lintemp_noint,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES cort; max temp * habitat interaction") +
    scale_x_continuous(trans = temp_trans_tres_cumdegreeweek,
                       breaks = c((0-mean_temp_tres_cumdegreeweek)/sd_temp_tres_cumdegreeweek,
                                  (1000-mean_temp_tres_cumdegreeweek)/sd_temp_tres_cumdegreeweek,
                                  (2000-mean_temp_tres_cumdegreeweek)/sd_temp_tres_cumdegreeweek,
                                  (3000-mean_temp_tres_cumdegreeweek)/sd_temp_tres_cumdegreeweek,
                                  (4000-mean_temp_tres_cumdegreeweek)/sd_temp_tres_cumdegreeweek),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres_cumdegreeweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/abs_cumdegreeweek_xhab_TRES.png",plot =  fig4_tres_cumdegreeweek, width = 10, height = 6.6)




(tresabs_cumdegreeweek_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresabs_cumdegreeweek_trendmax,"figures/tresabs_cumdegreeweek_trendmax.html")




save(list = ls(), file = "data/models_cort.RData")
