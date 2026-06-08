### Create all objects needed for survival analysis

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
library(emmeans)
# modelsummary::get_gof(my_model)
# performance::model_performance(my_model)
# wmean <- read_rds("../data/wmean.rds")

s <- read_rds("data/survival_attempt.rds") %>%
  mutate(juliandate_inc = yday(inc_date),
         juliandate_hatch = yday(hatch_date),
         year_fct = factor(year),
         site = case_match(site,
                           c("MBNC","MBNG","MBNR") ~ "MB",
                           c("RRC","RRRG","RRRR") ~ "RR",
                           c("PCC","PCE","PICG","PICO","PICR") ~ "PG",
                           ("WI" ~ "WI")),
         site = factor(site),
         eggs_hatched = `Eggs Hatched_num`,
         clutch_size = `Clutch Size_num`,
         nest_fledged = `Nestlings Fledging_num`,
         brood_size = `Brood Size_num`)




# ## Objectives
#
# 1. Exploration
# 1. Determine model to use for binomial outcome variable (number of successes in n trials)
# 1. Evaluate model:
#   - fit
# - qqplots
# - outliers
# - residuals vs temp
# 1. Test for interaction with land cover and temp
# 1. Produce graphs of relationship

## Exploration first


# ggplot(filter(s, Species == "WEBL"),aes(x = habitat, y = p_eggs_hatched)) + geom_boxplot()
#
# ggplot(filter(s, Species == "WEBL"),aes(x = meanmaxt_inc, y = p_eggs_hatched)) + geom_point() + facet_wrap(~habitat)
# # ggplot(filter(s, Species == "WEBL"),aes(x = meanmint_inc, y = p_eggs_hatched)) + geom_point() + facet_wrap(~habitat)
#
# ggplot(filter(s, Species == "WEBL"),aes(x = meanmaxt_nest, y = p_nest_fledged)) + geom_point() + facet_wrap(~habitat)
#
# ggplot(filter(s, Species == "WEBL"),aes(x = meanmaxt_nestpd, y = p_eggs_fledged)) + geom_point() + facet_wrap(~habitat)
#
#
#
#
# ggplot(filter(s, Species == "TRES"),aes(x = habitat, y = p_eggs_hatched)) + geom_boxplot()
# ggplot(filter(s, Species == "TRES"),aes(x = habitat, y = p_nest_fledged)) + geom_boxplot()
# ggplot(filter(s, Species == "TRES"),aes(x = habitat, y = p_eggs_fledged)) + geom_boxplot()
#
#
# ggplot(filter(s, Species == "TRES"),aes(x = meanmaxt_inc, y = p_eggs_hatched)) + geom_point() + facet_wrap(~habitat)
# # ggplot(filter(s, Species == "TRES"),aes(x = meanmint_inc, y = p_eggs_hatched)) + geom_point() + facet_wrap(~habitat)
#
# ggplot(filter(s, Species == "TRES"),aes(x = meanmaxt_nest, y = p_nest_fledged)) + geom_point() + facet_wrap(~habitat)
#
# ggplot(filter(s, Species == "TRES"),aes(x = meanmaxt_nestpd, y = p_eggs_fledged)) + geom_point() + facet_wrap(~habitat)
#
#
#
#
# ## Using glmer with family = "binomial"
#
# ### WEBL
#
# #### Incubation model
#
#
# s_inc_WEBL <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled * habitat + meanmint_inc_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                   family = binomial(link = "logit"),
#                   data = dplyr::filter(s,
#                                        Species == "WEBL",
#                                        !is.na(eggs_hatched),
#                                        !is.na(clutch_size)) %>%
#                     mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                   ~ scale(.x)[,1],
#                                   .names = "{.col}_scaled")
#                     )
# )
#
#
#
# Using site as a random effect results in a singular fit, so I left it as a fixed effect.
#
#
# summary(s_inc_WEBL)
#
#
#
#
# data = dplyr::filter(s,
#                      Species == "WEBL",
#                      !is.na(eggs_hatched),
#                      !is.na(clutch_size)) %>%
#   mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled")
#   )
#
# mean_temp <- mean(data %>% pull(meanmaxt_inc),na.rm = TRUE)
# sd_temp <- sd(data %>% pull(meanmaxt_inc),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_inc_WEBL$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
#
#
# (pl <- predict_response(s_inc_WEBL,terms = c("meanmaxt_inc_scaled [all]","habitat"),bias_correction = TRUE,margin = "empirical") %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over incubation period (\u00b0C)") +
#     ylab("Predicted proportion of eggs hatching") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = .3,label = label),hjust = -.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
#
#
#
#
#
# #### Nestling model
#
#
# s_nest_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled * habitat + meanmint_nest_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                    family = binomial(link = "logit"),
#                    data = dplyr::filter(s,
#                                         Species == "WEBL",
#                                         !is.na(nest_fledged),
#                                         !is.na(brood_size)) %>%
#                      mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_inc),
#                                    ~ scale(.x)[,1],
#                                    .names = "{.col}_scaled")
#                      )
# )
#
#
#
# Using site as a random effect results in a singular fit, so I left it as a fixed effect.
#
#
# summary(s_inc_WEBL)
#
#
#
#
# data = dplyr::filter(s,
#                      Species == "WEBL",
#                      !is.na(nest_fledged),
#                      !is.na(brood_size)) %>%
#   mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_inc),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled")
#   )
#
# mean_temp <- mean(data %>% pull(meanmaxt_nest),na.rm = TRUE)
# sd_temp <- sd(data %>% pull(meanmaxt_nest),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_nest_WEBL$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
#
#
# (pl <- predict_response(s_nest_WEBL,terms = c("meanmaxt_nest_scaled [all]","habitat"),bias_correction = TRUE,margin = "empirical") %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over nestling period (\u00b0C)") +
#     ylab("Predicted proportion of nestlings fledging") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = .3,label = label),hjust = -.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
#
#
#
#
#
# #### Nest period model
#
#
# s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                      family = binomial(link = "logit"),
#                      data = dplyr::filter(s,
#                                           Species == "WEBL",
#                                           !is.na(nest_fledged),
#                                           !is.na(clutch_size)) %>%
#                        mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_inc),
#                                      ~ scale(.x)[,1],
#                                      .names = "{.col}_scaled")
#                        )
# )
#
#
#
# Using site as a random effect results in a singular fit, so I left it as a fixed effect.
#
#
# summary(s_nestpd_WEBL)
#
#
#
#
#
# mean_temp <- mean(s_nestpd_WEBL$data %>% pull(meanmaxt_nestpd),na.rm = TRUE)
# sd_temp <- sd(s_nestpd_WEBL$data %>% pull(meanmaxt_nestpd),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_nestpd_WEBL$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
#
#
# (pl <- predict_response(s_nestpd_WEBL,terms = c("meanmaxt_nestpd_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp incubation to fledging (\u00b0C)") +
#     ylab("Predicted proportion of eggs fledging") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = .3,label = label),hjust = -.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
#
#
#
#
#
#
# ### TRES
#
# #### Incubation model
#
#
# s_inc_TRES <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled * habitat + meanmint_inc_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                   family = binomial(link = "logit"),
#                   data = dplyr::filter(s,
#                                        Species == "TRES",
#                                        !is.na(eggs_hatched),
#                                        !is.na(clutch_size),
#                                        meanmaxt_inc < 50) %>%
#                     mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                   ~ scale(.x)[,1],
#                                   .names = "{.col}_scaled")
#                     )
# )
#
#
#
# Using site as a random effect results in a singular fit, so I left it as a fixed effect.
#
#
# summary(s_inc_TRES)
#
#
#
#
# data = s_inc_TRES$data
#
# mean_temp <- mean(data %>% pull(meanmaxt_inc),na.rm = TRUE)
# sd_temp <- sd(data %>% pull(meanmaxt_inc),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_inc_TRES$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
#
#
# (pl <- predict_response(s_inc_TRES,terms = c("meanmaxt_inc_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over incubation period (\u00b0C)") +
#     ylab("Predicted proportion of eggs hatching") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = .3,label = label),hjust = -.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
#
#
#
#
#
# #### Nestling model
#
#
# s_nest_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled * habitat + meanmint_nest_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                    family = binomial(link = "logit"),
#                    data = dplyr::filter(s,
#                                         Species == "TRES",
#                                         !is.na(nest_fledged),
#                                         !is.na(brood_size),
#                                         meanmaxt_nest < 50) %>%
#                      mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_inc),
#                                    ~ scale(.x)[,1],
#                                    .names = "{.col}_scaled")
#                      )
# )
#
#
#
# Using site as a random effect results in a singular fit, so I left it as a fixed effect.
#
#
# summary(s_inc_TRES)
#
#
#
#
# data = s_inc_TRES$data
#
# mean_temp <- mean(data %>% pull(meanmaxt_nest),na.rm = TRUE)
# sd_temp <- sd(data %>% pull(meanmaxt_nest),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_nest_TRES$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
#
#
# (pl <- predict_response(s_nest_TRES,terms = c("meanmaxt_nest_scaled [all]","habitat"),bias_correction = TRUE) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over nestling period (\u00b0C)") +
#     ylab("Predicted proportion of nestlings fledging") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = .3,label = label),hjust = -.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
#
#
#
#
#
# #### Nest period model
#
#
# s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                      family = binomial(link = "logit"),
#                      data = dplyr::filter(s,
#                                           Species == "TRES",
#                                           !is.na(nest_fledged),
#                                           !is.na(clutch_size),
#                                           meanmaxt_nestpd < 50) %>%
#                        mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_inc),
#                                      ~ scale(.x)[,1],
#                                      .names = "{.col}_scaled")
#                        )
# )
#
#
#
# Using site as a random effect results in a singular fit, so I left it as a fixed effect.
#
#
# summary(s_nestpd_TRES)
#
#
#
#
#
# mean_temp <- mean(s_nestpd_TRES$data %>% pull(meanmaxt_nestpd),na.rm = TRUE)
# sd_temp <- sd(s_nestpd_TRES$data %>% pull(meanmaxt_nestpd),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_nestpd_TRES$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
#
#
# (pl <- predict_response(s_nestpd_TRES,terms = c("meanmaxt_nestpd_scaled [all]","habitat")) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp incubation to fledging (\u00b0C)") +
#     ylab("Predicted proportion of eggs fledging") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = .3,label = label),hjust = -.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
#
#



## Test for interaction with land cover and temp

### WEBL

#### Incubation model


# s_inc_WEBL <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled * habitat + meanmint_inc_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                   family = binomial(link = "logit"),
#                   data = dplyr::filter(s,
#                                        Species == "WEBL",
#                                        !is.na(eggs_hatched),
#                                        !is.na(clutch_size)) %>%
#                     mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                   ~ scale(.x)[,1],
#                                   .names = "{.col}_scaled")
#                     )
# )
#
# s_inc_WEBL_addmax <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled + meanmint_inc_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                          family = binomial(link = "logit"),
#                          data = dplyr::filter(s,
#                                               Species == "WEBL",
#                                               !is.na(eggs_hatched),
#                                               !is.na(clutch_size)) %>%
#                            mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                          ~ scale(.x)[,1],
#                                          .names = "{.col}_scaled")
#                            )
# )
#
# s_inc_WEBL_addmin <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled * habitat + meanmint_inc_scaled + juliandate_inc_scaled + year_fct + site,
#                          family = binomial(link = "logit"),
#                          data = dplyr::filter(s,
#                                               Species == "WEBL",
#                                               !is.na(eggs_hatched),
#                                               !is.na(clutch_size)) %>%
#                            mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                          ~ scale(.x)[,1],
#                                          .names = "{.col}_scaled")
#                            )
# )
# s_inc_WEBL_noint <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled + meanmint_inc_scaled + habitat + juliandate_inc_scaled + year_fct + site,
#                         family = binomial(link = "logit"),
#                         data = dplyr::filter(s,
#                                              Species == "WEBL",
#                                              !is.na(eggs_hatched),
#                                              !is.na(clutch_size)) %>%
#                           mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                         ~ scale(.x)[,1],
#                                         .names = "{.col}_scaled")
#                           )
# )
#
# c1 <- anova(s_inc_WEBL,s_inc_WEBL_addmin,s_inc_WEBL_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# c2 <- anova(s_inc_WEBL,s_inc_WEBL_addmax,s_inc_WEBL_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# (int_tab_survival_webl <- bind_rows(c1,c2) %>%
#     as_tibble() %>%
#     mutate(across(where(is.numeric),~round(.x,digits = 4)),
#            P = `Pr(>Chi)`) %>%
#     dplyr::select(Model,Deviance,P) %>%
#     mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
#            across(c(Deviance), ~ round(.x, digits = 2)),
#            P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
#     group_by(max_or_min) %>% gt())
# gtsave(int_tab_survival_webl,"figures/int_tab_survival_webl.html")
#
#
#
# #Conclusion: For WEBL incubation, temp does not interact with habitat.
#
#
# summary(s_inc_WEBL_noint)
#
#
# #Also no direct effect of max or min temp on hatching success; and no habitat effects.
#
#
# mean_temp <- mean(s_inc_WEBL_noint$data %>% pull(meanmaxt_inc),na.rm = TRUE)
# sd_temp <- sd(s_inc_WEBL_noint$data %>% pull(meanmaxt_inc),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- s_inc_WEBL_noint$data %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count)
# )
#
#
#
# (pl <- predict_response(s_inc_WEBL_noint,terms = c("meanmaxt_inc_scaled [all]")) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     #facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp incubation to fledging (\u00b0C)") +
#     ylab("Predicted proportion of eggs fledging") +
#     # scale_fill_viridis(discrete = TRUE) +
#     #  scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )
#
#
#
# #### Nestling model
#
#
# s_nest_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled * habitat + meanmint_nest_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
#                    family = binomial(link = "logit"),
#                    data = dplyr::filter(s,
#                                         Species == "WEBL",
#                                         !is.na(nest_fledged),
#                                         !is.na(brood_size)) %>%
#                      mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                    ~ scale(.x)[,1],
#                                    .names = "{.col}_scaled")
#                      )
# )
#
# s_nest_WEBL_addmax <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled + meanmint_nest_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
#                           family = binomial(link = "logit"),
#                           data = dplyr::filter(s,
#                                                Species == "WEBL",
#                                                !is.na(nest_fledged),
#                                                !is.na(brood_size)) %>%
#                             mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                           ~ scale(.x)[,1],
#                                           .names = "{.col}_scaled")
#                             )
# )
#
# s_nest_WEBL_addmin <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled * habitat + meanmint_nest_scaled + juliandate_hatch_scaled + year_fct + site,
#                           family = binomial(link = "logit"),
#                           data = dplyr::filter(s,
#                                                Species == "WEBL",
#                                                !is.na(nest_fledged),
#                                                !is.na(brood_size)) %>%
#                             mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                           ~ scale(.x)[,1],
#                                           .names = "{.col}_scaled")
#                             )
# )
# s_nest_WEBL_noint <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled + meanmint_nest_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
#                          family = binomial(link = "logit"),
#                          data = dplyr::filter(s,
#                                               Species == "WEBL",
#                                               !is.na(nest_fledged),
#                                               !is.na(brood_size)) %>%
#                            mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                          ~ scale(.x)[,1],
#                                          .names = "{.col}_scaled")
#                            )
# )
#
# c1 <- anova(s_nest_WEBL,s_nest_WEBL_addmin,s_nest_WEBL_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# c2 <- anova(s_nest_WEBL,s_nest_WEBL_addmax,s_nest_WEBL_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# (tab <- bind_rows(c1,c2) %>%
#     as_tibble() %>%
#     mutate(across(where(is.numeric),~round(.x,digits = 4)),
#            P = `Pr(>Chi)`) %>%
#     dplyr::select(Model,Deviance,P) %>%
#     mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
#            across(c(Deviance), ~ round(.x, digits = 2)),
#            P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
#     group_by(max_or_min) %>% gt())
# gtsave(tab,"../figures/int_tab_survival_nest_webl.html")
#
#
#
# Conclusion: For nestling period, max temp interacts with land cover.
#
#
# summary(s_nest_WEBL_addmin)
#
#
# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).
#
# #### Emtrends to calculate effect of max temp in each habitat
#
#
# (t <- emtrends(s_nest_WEBL_addmin,specs = pairwise ~ habitat, var = c("meanmaxt_nest_scaled")) %>% test() %>% pluck("emtrends") %>%
#     mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
#            df = round(df),
#            p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
#     rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nest_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
#     gt())
#
#
# gtsave(t,"../figures/weblsurvival_nest_trendmax.html")
#
# data = s_nest_WEBL_addmin$data
#
# (t <- emmeans(s_nest_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nest_scaled"), at = list(meanmaxt_nest_scaled = c(-2,0,2)),type = "response") %>%
#     as.tibble() %>% #gt() %>%
#     mutate(meanmaxt_nest_scaled = (meanmaxt_nest_scaled * sd(data$meanmaxt_nest,na.rm = TRUE)) + mean(data$meanmaxt_nest,na.rm = TRUE),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            meanmaxt_nest_scaled = round(meanmaxt_nest_scaled),
#            meanmaxt_nest_scaled = paste0(meanmaxt_nest_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`Max temperature` = "meanmaxt_nest_scaled",`Predicted P(survival)` = "prob",`2.5%` = "asymp.LCL",`97.5%` = "asymp.UCL" ) %>%
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
# gtsave(t,"../figures/weblsurvival_nest_deltamax.html")
#
# ((emmeans(s_nest_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nest_scaled"), at = list(meanmaxt_nest_scaled = c(2))) %>% as.tibble() %>% pull(emmean))-(emmeans(s_nest_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nest_scaled"), at = list(meanmaxt_nest_scaled = c(-2))) %>% as.tibble() %>% pull(emmean)))/(emmeans(s_nest_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nest_scaled"), at = list(meanmaxt_nest_scaled = c(-2))) %>% as.tibble() %>% pull(emmean))
#
# emmeans(s_nest_WEBL_addmin,specs = pairwise ~ habitat,by = c("meanmaxt_nest_scaled"), at = list(meanmaxt_nest_scaled = c(-2,0,2))) %>% plot()
# emmip(s_nest_WEBL_addmin,formula = habitat ~ meanmaxt_nest_scaled, at = list(meanmaxt_nest_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()
#
#
#
# ## Emmeans to check for effect of habitat
#
#
# (t <- emmeans(s_nest_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
#     mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
#            across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nest_byhabitat_webl.html")
#
#
# ### Check for effect of temperature
#
#
# (t <- summary(s_nest_WEBL_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#     # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#     mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
#            across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nest_byhabitat_summary_webl.html")
#
#
#
#
# data = s_nest_WEBL_addmin$data
#
#
# mean_temp <- mean(data %>% pull(meanmaxt_nest),na.rm = TRUE)
# sd_temp <- sd(data %>% pull(meanmaxt_nest),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                         transform = function(x){(x * sd_temp) + mean_temp},
#                         inverse = function(x){x})
#
# samp <- data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
#   label = paste("N =",samp$count),
#   group   = factor(samp$habitat)
# )
#
# (pl <- predict_response(s_nest_WEBL_addmin,terms = c("meanmaxt_nest_scaled [all]","habitat")) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp over preceding week (\u00b0C)") +
#     ylab("Predicted proportion of nestlings fledging") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "WEBL; linear max interacts with land cover") +
#     scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                        ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#     ) +
#     # ylim(0,100) +
#     geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
#     theme(legend.position = "none")
# )





#### Nest period model


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
                                          Species == "WEBL",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_WEBL_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(s,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_WEBL_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(s,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_WEBL_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(s,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

c1 <- anova(s_nestpd_WEBL,s_nestpd_WEBL_addmin,s_nestpd_WEBL_noint,test="Chisq") %>% tibble() %>%
  rename(Chisq = Deviance) %>%
  mutate(AIC = c(AIC(s_nestpd_WEBL),AIC(s_nestpd_WEBL_addmin),AIC(s_nestpd_WEBL_noint)),.before = Df) %>%
  mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)

c2 <- anova(s_nestpd_WEBL,s_nestpd_WEBL_addmax,s_nestpd_WEBL_noint,test="Chisq") %>% tibble() %>%
  rename(Chisq = Deviance) %>%
  mutate(AIC = c(AIC(s_nestpd_WEBL),
                 AIC(s_nestpd_WEBL_addmax),
                 AIC(s_nestpd_WEBL_noint)),.before = Df) %>%
  mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)

(int_tab_survival_nestpd_webl <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_survival_nestpd_webl,"figures/int_tab_survival_nestpd_webl.html")



#Conclusion: For nest attempt overall, land cover interacts with either max or min temp but not both together. It looks like the mean temp interaction model is slightly more explanatory so we'll go with that.


summary(s_nestpd_WEBL_addmin)
(survivalbyhabitat_summary_webl <- summary(s_nestpd_WEBL_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
   # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
   mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
          across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(survivalbyhabitat_summary_webl,"figures/survivalbyhabitat_summary_webl.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax <- emtrends(s_nestpd_WEBL_addmin,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
   gt())


# gtsave(weblsurvival_nestpd_trendmax,"figures/weblsurvival_nestpd_trendmax.html")

data = s_nestpd_WEBL_addmin$data

# (t <- emmeans(s_nestpd_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nestpd_scaled"), at = list(meanmaxt_nestpd_scaled = c(-2,0,2)),type = "response") %>%
#     as.tibble() %>% #gt() %>%
#   mutate(meanmaxt_nestpd_scaled = (meanmaxt_nestpd_scaled * sd(data$meanmaxt_nestpd,na.rm = TRUE)) + mean(data$meanmaxt_nestpd,na.rm = TRUE),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     meanmaxt_nestpd_scaled = round(meanmaxt_nestpd_scaled),
#     meanmaxt_nestpd_scaled = paste0(meanmaxt_nestpd_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "meanmaxt_nestpd_scaled",`Predicted P(survival)` = "prob",`2.5%` = "asymp.LCL",`97.5%` = "asymp.UCL" ) %>%
#   group_by(`Max temperature`) %>%
#   mutate(row=row_number()) %>%
#   pivot_longer(-c(`Max temperature`, row,Habitat)) %>%
#   pivot_wider(names_from=c(`Max temperature`, name), values_from=value) %>%
#   dplyr::select(-row) %>%
#   #mutate(`Maximum TA_P` = if_else(`Maximum TA_P` == 0.000,"<0.001",as.character(`Maximum TA_P`))) %>%
#   #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#   gt() %>% tab_options(data_row.padding = px(1)) %>%
#   tab_spanner_delim(
#     delim="_"
#   ))
#
# gtsave(t,"figures/weblsurvival_nestpd_deltamax.html")
#
# ((emmeans(s_nestpd_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nestpd_scaled"), at = list(meanmaxt_nestpd_scaled = c(2))) %>% as.tibble() %>% pull(emmean))-(emmeans(s_nestpd_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nestpd_scaled"), at = list(meanmaxt_nestpd_scaled = c(-2))) %>% as.tibble() %>% pull(emmean)))/(emmeans(s_nestpd_WEBL_addmin,specs = ~ habitat,by = c("meanmaxt_nestpd_scaled"), at = list(meanmaxt_nestpd_scaled = c(-2))) %>% as.tibble() %>% pull(emmean))
#
# emmeans(s_nestpd_WEBL_addmin,specs = pairwise ~ habitat,by = c("meanmaxt_nestpd_scaled"), at = list(meanmaxt_nestpd_scaled = c(-2,0,2))) %>% plot()
# emmip(s_nestpd_WEBL_addmin,formula = habitat ~ meanmaxt_nestpd_scaled, at = list(meanmaxt_nestpd_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl <- emmeans(s_nestpd_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
   mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
          across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(survival_nestpd_byhabitat_webl,"figures/survival_nestpd_byhabitat_webl.html")


# Forest survival is lower than in the other land covers.

### Check for effect of temperature


# (t <- summary(s_nestpd_WEBL_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
#           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nestpd_byhabitat_summary_webl.html")
#
#
# Hot temps and low temps reduce survival overall.


data_webl = s_nestpd_WEBL_addmin$data


mean_temp_webl <- mean(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                          transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                          inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt()

ss_year_survival_webl <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

  ss_year_survival_webl %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_survival_webl.html")


dat_text_webl <- data.frame(
  label = paste("N =",samp_webl$count),
  group   = factor(samp_webl$habitat)
)

(fig3_webl <- predict_response(s_nestpd_WEBL_addmin,terms = c("meanmaxt_nestpd_scaled [all]","habitat")) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Predicted percent of eggs fledging") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
   scale_x_continuous(trans = temp_trans_webl,
                       breaks = c((20-mean_temp_webl)/sd_temp_webl,
                                  (25-mean_temp_webl)/sd_temp_webl,
                                  (30-mean_temp_webl)/sd_temp_webl,
                                  (35-mean_temp_webl)/sd_temp_webl,
                                  (40-mean_temp_webl)/sd_temp_webl#,
                                  #(45-mean_temp)/sd_temp
                                  ),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
                       ) +
   # ylim(0,100) +
   scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dashed","Grassland" = "dotted","Row crop" = "dashed")) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )

# ggsave("figures/survbytempxhab_webl.png",fig3_webl,width = 10, height = 6.6)


### TRES

#### Incubation model


# s_inc_TRES <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled * habitat + meanmint_inc_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(eggs_hatched),
#                                     !is.na(clutch_size)) %>%
#                  mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
#
# s_inc_TRES_addmax <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled + meanmint_inc_scaled * habitat + juliandate_inc_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(eggs_hatched),
#                                     !is.na(clutch_size)) %>%
#                  mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
#
# s_inc_TRES_addmin <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled * habitat + meanmint_inc_scaled + juliandate_inc_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(eggs_hatched),
#                                     !is.na(clutch_size)) %>%
#                  mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
# s_inc_TRES_noint <- glm(cbind(eggs_hatched,clutch_size - eggs_hatched) ~ meanmaxt_inc_scaled + meanmint_inc_scaled + habitat + juliandate_inc_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(eggs_hatched),
#                                     !is.na(clutch_size)) %>%
#                  mutate(across(c(meanmaxt_inc,meanmint_inc,juliandate_inc),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
#
# c1 <- anova(s_inc_TRES,s_inc_TRES_addmin,s_inc_TRES_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# c2 <- anova(s_inc_TRES,s_inc_TRES_addmax,s_inc_TRES_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# (tab <- bind_rows(c1,c2) %>%
#   as_tibble() %>%
#   mutate(across(where(is.numeric),~round(.x,digits = 4)),
#          P = `Pr(>Chi)`) %>%
#   dplyr::select(Model,Deviance,P) %>%
#   mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
#          across(c(Deviance), ~ round(.x, digits = 2)),
#          P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
#   group_by(max_or_min) %>% gt())
# gtsave(tab,"../figures/int_tab_survival_tres.html")
#
#
#
# Land cover does not interact with max or min temp during incubation.
#
#
# summary(s_inc_TRES_noint)
#
#
# However, there is a direct negative effect of high temps. Also, survival is generally higher in Grassland and Row crop than Forest.
#
#
# (t <- emmeans(s_inc_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
#    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
#           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_inc_byhabitat_tres.html")
#
#
# Posthoc test does not support the difference between row crop and forest?
#
#
# mean_temp <- mean(s_inc_TRES_noint$data %>% pull(meanmaxt_inc),na.rm = TRUE)
# sd_temp <- sd(s_inc_TRES_noint$data %>% pull(meanmaxt_inc),na.rm = TRUE)
#
#
# temp_trans <- trans_new("temp_trans",
#                           transform = function(x){(x * sd_temp) + mean_temp},
#                            inverse = function(x){x})
#
# samp <- s_inc_TRES_noint$data %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
# label = paste("N =",samp$count)
# )
#
#
#
# (pl <- predict_response(s_inc_TRES_noint,terms = c("meanmaxt_inc_scaled [all]")) %>%
#    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#    theme_classic() +
#    #facet_wrap(~ group, ncol = 2) +
#     xlab("Mean daily max temp incubation to fledging (\u00b0C)") +
#     ylab("Predicted proportion of eggs fledging") +
#    # scale_fill_viridis(discrete = TRUE) +
#    #  scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max") +
#    scale_x_continuous(trans = temp_trans,
#                        breaks = c((20-mean_temp)/sd_temp,
#                                   (25-mean_temp)/sd_temp,
#                                   (30-mean_temp)/sd_temp,
#                                   (35-mean_temp)/sd_temp,
#                                   (40-mean_temp)/sd_temp#,
#                                   #(45-mean_temp)/sd_temp
#                                   ),
#                        # breaks = c(20,30,40,50),
#                        # labels = c("20","30","40","50"),
#                        # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#                        #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#                        # limits = c(18,5)
#                        ) +
#    # ylim(0,100) +
#    geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
#     theme(legend.position = "none")
#  )
#
#
#
#
# samp <- s_inc_TRES_noint$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
# label = paste("N =",samp$count),
# group = factor(samp$habitat)
# )
#
#
#
# dat <- predict_response(s_inc_TRES_noint,terms = c("habitat"))
#
# (pl <- ggplot(data = dat, mapping = aes(x = x, color = x, y = predicted,ymin = conf.low,ymax = conf.high)) +
#    # plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE,colors = viridis(4)) +
#    geom_point(size = 4) +
#     geom_linerange(linewidth = 2) +
#    theme_classic() +
#    #facet_wrap(~ group, ncol = 2) +
#     xlab("Land cover") +
#     ylab("Predicted proportion of eggs fledging") +
#    # scale_fill_viridis(discrete = TRUE) +
#    scale_color_manual(values = viridis(4)) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max") +
#    # scale_x_continuous(trans = temp_trans,
#    #                     breaks = c((20-mean_temp)/sd_temp,
#    #                                (25-mean_temp)/sd_temp,
#    #                                (30-mean_temp)/sd_temp,
#    #                                (35-mean_temp)/sd_temp,
#    #                                (40-mean_temp)/sd_temp#,
#    #                                #(45-mean_temp)/sd_temp
#    #                                ),
#    #                     # breaks = c(20,30,40,50),
#    #                     # labels = c("20","30","40","50"),
#    #                     # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#    #                     #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#    #                     # limits = c(18,5)
#    #                     ) +
#    # ylim(0,100) +
#    geom_text(data = dat_text, mapping = aes(x = c(1,2,3,4), y = -Inf,label = label),hjust = .5,vjust = -.7,inherit.aes = FALSE) +
#     theme(legend.position = "none") +
#     annotate(geom = "text", x = c(1,2,3,4),y = Inf,hjust = .5,vjust = 1,label = c("a","ab","b","ab"))
#  )
#
#
#
# #### Nestling model
#
#
# s_nest_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled * habitat + meanmint_nest_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(nest_fledged),
#                                     !is.na(brood_size)) %>%
#                  mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
#
# s_nest_TRES_addmax <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled + meanmint_nest_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(nest_fledged),
#                                     !is.na(brood_size)) %>%
#                  mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
#
# s_nest_TRES_addmin <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled * habitat + meanmint_nest_scaled + juliandate_hatch_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(nest_fledged),
#                                     !is.na(brood_size)) %>%
#                  mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
# s_nest_TRES_noint <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meanmaxt_nest_scaled + meanmint_nest_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
#                family = binomial(link = "logit"),
#                data = dplyr::filter(s,
#                                     Species == "TRES",
#                                     !is.na(nest_fledged),
#                                     !is.na(brood_size)) %>%
#                  mutate(across(c(meanmaxt_nest,meanmint_nest,juliandate_hatch),
#                                ~ scale(.x)[,1],
#                                .names = "{.col}_scaled")
#                       )
#                )
#
# c1 <- anova(s_nest_TRES,s_nest_TRES_addmin,s_nest_TRES_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# c2 <- anova(s_nest_TRES,s_nest_TRES_addmax,s_nest_TRES_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)
#
# (tab <- bind_rows(c1,c2) %>%
#   as_tibble() %>%
#   mutate(across(where(is.numeric),~round(.x,digits = 4)),
#          P = `Pr(>Chi)`) %>%
#   dplyr::select(Model,Deviance,P) %>%
#   mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
#          across(c(Deviance), ~ round(.x, digits = 2)),
#          P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
#   group_by(max_or_min) %>% gt())
# gtsave(tab,"../figures/int_tab_survival_nest_tres.html")
#
#
#
# For TRES nestlings, neither max nor mean temp interacts with land cover.
#
#
# summary(s_nest_TRES_noint)
#
#
# No main effects of temp but potentially main effects of habitat- yes, grassland and row crop have higher survival than forest.
#
#
# (t <- emmeans(s_nest_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
#    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
#           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nest_byhabitat_tres.html")
#
#
#
#
#
#
# samp <- s_nest_TRES_noint$data %>% group_by(habitat) %>% summarize(count = n())
# # samp %>% gt()
#
#
# dat_text <- data.frame(
# label = paste("N =",samp$count),
# group = factor(samp$habitat)
# )
#
#
#
# dat <- predict_response(s_nest_TRES_noint,terms = c("habitat"))
#
# (pl <- ggplot(data = dat, mapping = aes(x = x, color = x, y = predicted,ymin = conf.low,ymax = conf.high)) +
#    # plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE,colors = viridis(4)) +
#    geom_point(size = 4) +
#     geom_linerange(linewidth = 2) +
#    theme_classic() +
#    #facet_wrap(~ group, ncol = 2) +
#     xlab("Land cover") +
#     ylab("Predicted proportion of nestlings fledging") +
#    # scale_fill_viridis(discrete = TRUE) +
#    scale_color_manual(values = viridis(4)) +
#     theme(text = element_text(size = 16)) +
#     labs(title = "TRES; linear max") +
#    # scale_x_continuous(trans = temp_trans,
#    #                     breaks = c((20-mean_temp)/sd_temp,
#    #                                (25-mean_temp)/sd_temp,
#    #                                (30-mean_temp)/sd_temp,
#    #                                (35-mean_temp)/sd_temp,
#    #                                (40-mean_temp)/sd_temp#,
#    #                                #(45-mean_temp)/sd_temp
#    #                                ),
#    #                     # breaks = c(20,30,40,50),
#    #                     # labels = c("20","30","40","50"),
#    #                     # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
#    #                     #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
#    #                     # limits = c(18,5)
#    #                     ) +
#    ylim(0,1.1) +
#    geom_text(data = dat_text, mapping = aes(x = c(1,2,3,4), y = -Inf,label = label),hjust = .5,vjust = -1,inherit.aes = FALSE) +
#     theme(legend.position = "none") +
#     annotate(geom = "text", x = c(1,2,3,4),y = Inf,hjust = .5,vjust = 1,label = c("a","ab","b","b"))
#  )


#### Nest period model


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
               family = binomial(link = "logit"),
               data = dplyr::filter(s,
                                    Species == "TRES",
                                    !is.na(nest_fledged),
                                    !is.na(clutch_size)) %>%
                 mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                               ~ scale(.x)[,1],
                               .names = "{.col}_scaled")
                      )
               )

s_nestpd_TRES_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
               family = binomial(link = "logit"),
               data = dplyr::filter(s,
                                    Species == "TRES",
                                    !is.na(nest_fledged),
                                    !is.na(clutch_size)) %>%
                 mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                               ~ scale(.x)[,1],
                               .names = "{.col}_scaled")
                      )
               )

s_nestpd_TRES_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
               family = binomial(link = "logit"),
               data = dplyr::filter(s,
                                    Species == "TRES",
                                    !is.na(nest_fledged),
                                    !is.na(clutch_size)) %>%
                 mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                               ~ scale(.x)[,1],
                               .names = "{.col}_scaled")
                      )
               )
s_nestpd_TRES_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
               family = binomial(link = "logit"),
               data = dplyr::filter(s,
                                    Species == "TRES",
                                    !is.na(nest_fledged),
                                    !is.na(clutch_size)) %>%
                 mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                               ~ scale(.x)[,1],
                               .names = "{.col}_scaled")
                      )
               )


c1 <- anova(s_nestpd_TRES,s_nestpd_TRES_addmin,s_nestpd_TRES_noint,test="Chisq") %>% tibble() %>%
  rename(Chisq = Deviance) %>%
  mutate(AIC = c(AIC(s_nestpd_TRES),AIC(s_nestpd_TRES_addmin),AIC(s_nestpd_TRES_noint)),.before = Df) %>%
  mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)

c2 <- anova(s_nestpd_TRES,s_nestpd_TRES_addmax,s_nestpd_TRES_noint,test="Chisq") %>% tibble() %>%
  rename(Chisq = Deviance) %>%
  mutate(AIC = c(AIC(s_nestpd_TRES),
                 AIC(s_nestpd_TRES_addmax),
                 AIC(s_nestpd_TRES_noint)),.before = Df) %>%
  mutate(Model = c("no interaction","single interaction","both interacting"),.before = Df)

(int_tab_survival_nestpd_tres <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_survival_nestpd_tres,"figures/int_tab_survival_nestpd_tres.html")




#No interaction of land cover with max or min temp for TRES nest period.

(tressurvival_nestpd_trendmax <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


# gtsave(tressurvival_nestpd_trendmax,"figures/tressurvival_nestpd_trendmax.html")


summary(s_nestpd_TRES_noint)


#No main effects of max or min temp. Effects of habitat?


(survival_nestpd_byhabitat_tres <- emmeans(s_nestpd_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
   mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
          across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(survival_nestpd_byhabitat_tres,"figures/survival_nestpd_byhabitat_tres.html")



# Yes - survival is higher in grassland and row crop than in forest.



samp_tres <- s_nestpd_TRES_noint$data %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt()
ss_year_survival_tres <- s_nestpd_TRES_noint$data %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
ss_year_survival_tres %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))
  # gtsave("figures/ss_year_survival_tres.html")


dat_text_tres <- data.frame(
label = paste("N =",samp_tres$count),
group = factor(samp_tres$habitat)
)

prop_trans_tres <- trans_new("prop_trans_tres",
                             transform = function(x){(x * 100)},
                             inverse = function(x){x})

dat_tres <- predict_response(s_nestpd_TRES_noint,terms = c("habitat"))

(fig3_tres <- ggplot(data = dat_tres, mapping = aes(x = x, color = x, y = predicted,ymin = conf.low,ymax = conf.high)) +
   # plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE,colors = viridis(4)) +
   geom_point(size = 4) +
    geom_linerange(linewidth = 2) +
   theme_classic() +
   #facet_wrap(~ group, ncol = 2) +
    xlab("Land cover") +
    ylab("Predicted percent of eggs fledging") +
   # scale_fill_viridis(discrete = TRUE) +
   scale_color_manual(values = viridis(4)) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_y_continuous(trans = prop_trans_tres,
                       breaks = c(.2,.4,.6,.8),
                       labels = c("20%","40%","60%","80%")) +
   # scale_x_continuous(trans = temp_trans,
   #                     breaks = c((20-mean_temp)/sd_temp,
   #                                (25-mean_temp)/sd_temp,
   #                                (30-mean_temp)/sd_temp,
   #                                (35-mean_temp)/sd_temp,
   #                                (40-mean_temp)/sd_temp#,
   #                                #(45-mean_temp)/sd_temp
   #                                ),
   #                     # breaks = c(20,30,40,50),
   #                     # labels = c("20","30","40","50"),
   #                     # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
   #                     #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
   #                     # limits = c(18,5)
   #                     ) +
   # ylim(0,1.1) +
   geom_text(data = dat_text_tres, mapping = aes(x = c(1,2,3,4), y = -Inf,label = label),hjust = .5,vjust = -.2,inherit.aes = FALSE) +
    theme(legend.position = "none") +
    annotate(geom = "text", x = c(1,2,3,4),y = Inf,hjust = .5,vjust = 1,label = c("a","ab","b","b"))
 )

# ggsave("figures/survbytempxhab_tres.png",fig3_tres,width = 10, height = 6.6)




## combined survival results for manuscript


library(egg)
ggplot_build(fig3_webl)$layout$panel_scales_y
(p_full <- ggarrange(fig3_webl + theme(text = element_text(size = 12)),fig3_tres + theme(text = element_text(size = 12),axis.title.y = element_blank()),ncol = 2,
          labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig3_survival_by_temp_hab.png",p_full,width = 6.25,height = 4)




# Effect of provis and cort on survival - sample size is 26 WEBL nestlings and 3 TRES nestlings??? yep, mostly because cort coverage is bad.


# s <- read_rds("../data/growth_with_survival_indiv.rds")
#
#
#
# s_webl <- glm(survive_to_next_week ~ cort_s1_scaled + occ3_mean_provis_scaled + occ3_mean_meanmaxt_scaled + occ3_mean_meanmint_scaled + age_scaled + condition_scaled,
#          family = binomial(link="logit"),
#          data = dplyr::filter(s,
#                                                      Species == "WEBL",
#                                                      !is.na(cort_s1),
#                                                      !is.na(occ3_mean_provis),
#                                                      !is.na(occ3_mean_meanmaxt),
#                                                      !is.na(occ3_mean_meanmint),
#                                                      !is.na(age),
#                                                      !is.na(condition)) %>%
#   mutate(across(c(cort_s1,occ3_mean_provis,occ3_mean_meanmaxt,occ3_mean_meanmint,age,condition),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled")))
#
# summary(s_webl)
#
#
#
#
# s_webl <- glm(survive_to_next_week ~ occ3_mean_provis_scaled,
#          family = binomial(link="logit"),
#          data = dplyr::filter(s,
#                                                      Species == "WEBL",
#                                                      !is.na(occ3_mean_provis)) %>%
#   mutate(across(c(cort_s1,occ3_mean_provis,occ3_mean_meanmaxt,occ3_mean_meanmint,age,condition),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled")))
#
# summary(s_webl)


#There's just not enough data to draw any meaningful conclusions here.


### What about trying to estimate the effect of provis on survival? Will need to return to DAG to know what to put in the model.
#DAG says that to estimate the direct effect of provis on survival, I'll need to include cort. So we're back to very low sample size. Shoot.



# s_webl <- glm(survive_to_next_week ~ mean_provis_scaled + _scaled + occ3_mean_meanmint_scaled + age_scaled + condition_scaled,
#               family = binomial(link="logit"),
#               data = dplyr::filter(s,
#                                    Species == "WEBL",
#                                    !is.na(cort_s1),
#                                    !is.na(occ3_mean_provis),
#                                    !is.na(occ3_mean_meanmaxt),
#                                    !is.na(occ3_mean_meanmint),
#                                    !is.na(age),
#                                    !is.na(condition)) %>%
#                 mutate(across(c(cort_s1,occ3_mean_provis,occ3_mean_meanmaxt,occ3_mean_meanmint,age,condition),
#                               ~ scale(.x)[,1],
#                               .names = "{.col}_scaled")))
#
# summary(s_webl)



# The following is attempt-level, but Danny and I discussed trying an individual model.

# I think it only really makes sense to estimate the relationship between nestling survival and cort/provis.
#
#
# p <- read_rds("../data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
#   mutate(year = year(date),
#          year_fct = as.factor(year))
#
# g <- read_rds("../data/growth_and_provis_mobilenetv3-original_dataset.h5.rds") %>%
#   mutate(abs_change_cort = cort_s2-cort_s1,
#          prop_change_cort = (cort_s2-cort_s1)/cort_s1,
#          year_fct = as.factor(year))
#
#
# add_provis_and_cort <- function(row){
#   # print(row)
#   nestbox <- row$Nestbox
#   egg_start <- row$`First Egg Date_lubridate`
#   inc_start <- row$`Incubation Begins Date_lubridate`
#   egg_end <- row$`Hatch Date_lubridate`
#   if(is.na(egg_end)){
#     egg_end <- egg_start + ddays(14)
#   }
#   nest_end <- row$`Date Banded_lubridate`
#   if(is.na(nest_end)){
#     nest_end <- egg_start + ddays(30)
#   }
#   temp_provis <- filter(p,nestbox == nestbox,date >= egg_end & date <= nest_end)
#   temp_cort <- filter(g,nestbox == Nestbox,date >= egg_end & date <= nest_end)
#   if(nrow(temp_provis) == 0){
#     meanprovis <- NA
#     meanaftprovis <- NA
#     meanmornprovis <- NA
#   } else {
#     meanprovis <- mean(temp_provis$valid_detections,na.rm = TRUE)
#     meanaftprovis <- filter(temp_provis,tod_h >= 12) %>% pull(valid_detections) %>% mean(na.rm = TRUE)
#     meanmornprovis <- filter(temp_provis,tod_h < 12) %>% pull(valid_detections) %>% mean(na.rm = TRUE)
#   }
#   if(nrow(temp_cort) == 0){
#     meancort_s1 <- NA
#     meancort_absdiff <- NA
#     meancort_propdiff <- NA
#     meanage <- NA
#     meancond <- NA
#     meanh <- NA
#   } else {
#     temp <- filter(temp_cort,occasion == max(occasion,na.rm = TRUE))
#     meancort_s1 <- pull(temp_cort,cort_s1) %>% mean(na.rm = TRUE)
#     meancort_absdiff <- mutate(temp_cort,absdiff = cort_s2 - cort_s1) %>% pull(absdiff) %>% mean(na.rm = TRUE)
#     meancort_propdiff <- mutate(temp_cort,propdiff = (cort_s2 - cort_s1)/cort_s1) %>% pull(propdiff) %>% mean(na.rm = TRUE)
#     meanage <- pull(temp_cort,age) %>% mean(na.rm = TRUE)
#     meancond <- pull(temp_cort,condition) %>% mean(na.rm = TRUE)
#     meanh <- pull(temp_cort,meanh) %>% mean(na.rm = TRUE)
#   }
#   row <- bind_cols(row,
#                    meanprovis = meanprovis,
#                    meanaftprovis = meanaftprovis,
#                    meanmornprovis = meanmornprovis,
#                    meancort_s1 = meancort_s1,
#                    meancort_absdiff = meancort_absdiff,
#                    meancort_propdiff = meancort_propdiff,
#                    meanage = meanage,
#                    meancond = meancond,
#                    meanh = meanh)
#   return(row)
# }
#
# future::plan(multisession,workers = availableCores()-2)
#
# s_with_provis_and_cort <- s %>% rowwise() %>% group_split() %>% future_map(~ add_provis_and_cort(.x),.progress = TRUE) %>% list_rbind()
#
#
#
# ## WEBL
#
# ### Overall provis
#
#
#
# s_provis_cort_s1_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_s1_scaled + meanprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanage_scaled + meancond_scaled,
#                              family = binomial(link="logit"),
#                              data = dplyr::filter(s_with_provis_and_cort,
#                                                   Species == "WEBL",
#                                                   !is.na(meanmaxt_nest),
#                                                   !is.na(meanmint_nest),
#                                                   !is.na(meancort_s1),
#                                                   !is.na(meanprovis),
#                                                   !is.na(meanage),
#                                                   !is.na(meancond)) %>%
#                                mutate(across(c(meancort_s1,meanprovis,meanmaxt_nest,meanmint_nest,meanage,meancond),
#                                              ~ scale(.x)[,1],
#                                              .names = "{.col}_scaled")))
# summary(s_provis_cort_s1_WEBL)
#
#
#
#
# s_provis_cort_absdiff_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_absdiff_scaled + meanprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanage_scaled + meancond_scaled,
#                                   family = binomial(link="logit"),
#                                   data = dplyr::filter(s_with_provis_and_cort,
#                                                        Species == "WEBL",
#                                                        !is.na(meanmaxt_nest),
#                                                        !is.na(meanmint_nest),
#                                                        !is.na(meancort_absdiff),
#                                                        !is.na(meanprovis),
#                                                        !is.na(meanage),
#                                                        !is.na(meancond)) %>%
#                                     mutate(across(c(meancort_absdiff,meanprovis,meanmaxt_nest,meanmint_nest,meanage,meancond),
#                                                   ~ scale(.x)[,1],
#                                                   .names = "{.col}_scaled")))
# summary(s_provis_cort_absdiff_WEBL)
#
#
#
#
#
# s_provis_cort_propdiff_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_propdiff_scaled + meanprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                    family = binomial(link="logit"),
#                                    data = dplyr::filter(s_with_provis_and_cort,
#                                                         Species == "WEBL",
#                                                         !is.na(meanmaxt_nest),
#                                                         !is.na(meanmint_nest),
#                                                         !is.na(meancort_propdiff),
#                                                         !is.na(meanprovis),
#                                                         !is.na(meanh),
#                                                         !is.na(meanage),
#                                                         !is.na(meancond)) %>%
#                                      mutate(across(c(meancort_propdiff,meanprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                    ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled")))
# summary(s_provis_cort_propdiff_WEBL)
#
#
#
# ### Morning cort
#
#
#
# s_mornprovis_cort_s1_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_s1_scaled + meanmornprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                  family = binomial(link="logit"),
#                                  data = dplyr::filter(s_with_provis_and_cort,
#                                                       Species == "WEBL",
#                                                       !is.na(meanmaxt_nest),
#                                                       !is.na(meanmint_nest),
#                                                       !is.na(meancort_s1),
#                                                       !is.na(meanmornprovis),
#                                                       !is.na(meanh),
#                                                       !is.na(meanage),
#                                                       !is.na(meancond)) %>%
#                                    mutate(across(c(meancort_s1,meanmornprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                  ~ scale(.x)[,1],
#                                                  .names = "{.col}_scaled")))
# summary(s_mornprovis_cort_s1_WEBL)
#
#
#
#
# s_mornprovis_cort_absdiff_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_absdiff_scaled + meanmornprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                       family = binomial(link="logit"),
#                                       data = dplyr::filter(s_with_provis_and_cort,
#                                                            Species == "WEBL",
#                                                            !is.na(meanmaxt_nest),
#                                                            !is.na(meanmint_nest),
#                                                            !is.na(meancort_absdiff),
#                                                            !is.na(meanmornprovis),
#                                                            !is.na(meanh),
#                                                            !is.na(meanage),
#                                                            !is.na(meancond)) %>%
#                                         mutate(across(c(meancort_absdiff,meanmornprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                       ~ scale(.x)[,1],
#                                                       .names = "{.col}_scaled")))
# summary(s_mornprovis_cort_absdiff_WEBL)
#
#
#
#
#
# s_mornprovis_cort_propdiff_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_propdiff_scaled + meanmornprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                        family = binomial(link="logit"),
#                                        data = dplyr::filter(s_with_provis_and_cort,
#                                                             Species == "WEBL",
#                                                             !is.na(meanmaxt_nest),
#                                                             !is.na(meanmint_nest),
#                                                             !is.na(meancort_propdiff),
#                                                             !is.na(meanmornprovis),
#                                                             !is.na(meanh),
#                                                             !is.na(meanage),
#                                                             !is.na(meancond)) %>%
#                                          mutate(across(c(meancort_propdiff,meanmornprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                        ~ scale(.x)[,1],
#                                                        .names = "{.col}_scaled")))
# summary(s_mornprovis_cort_propdiff_WEBL)
#
#
#
# ### Afternoon cort
#
#
#
# s_aftprovis_cort_s1_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_s1_scaled + meanaftprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                 family = binomial(link="logit"),
#                                 data = dplyr::filter(s_with_provis_and_cort,
#                                                      Species == "WEBL",
#                                                      !is.na(meanmaxt_nest),
#                                                      !is.na(meanmint_nest),
#                                                      !is.na(meancort_s1),
#                                                      !is.na(meanaftprovis),
#                                                      !is.na(meanh),
#                                                      !is.na(meanage),
#                                                      !is.na(meancond)) %>%
#                                   mutate(across(c(meancort_s1,meanaftprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                 ~ scale(.x)[,1],
#                                                 .names = "{.col}_scaled")))
# summary(s_aftprovis_cort_s1_WEBL)
#
#
#
#
# s_aftprovis_cort_absdiff_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_absdiff_scaled + meanaftprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                      family = binomial(link="logit"),
#                                      data = dplyr::filter(s_with_provis_and_cort,
#                                                           Species == "WEBL",
#                                                           !is.na(meanmaxt_nest),
#                                                           !is.na(meanmint_nest),
#                                                           !is.na(meancort_absdiff),
#                                                           !is.na(meanaftprovis),
#                                                           !is.na(meanh),
#                                                           !is.na(meanage),
#                                                           !is.na(meancond)) %>%
#                                        mutate(across(c(meancort_absdiff,meanaftprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                      ~ scale(.x)[,1],
#                                                      .names = "{.col}_scaled")))
# summary(s_aftprovis_cort_absdiff_WEBL)
#
#
#
#
#
# s_aftprovis_cort_propdiff_WEBL <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_propdiff_scaled + meanaftprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                       family = binomial(link="logit"),
#                                       data = dplyr::filter(s_with_provis_and_cort,
#                                                            Species == "WEBL",
#                                                            !is.na(meanmaxt_nest),
#                                                            !is.na(meanmint_nest),
#                                                            !is.na(meancort_propdiff),
#                                                            !is.na(meanaftprovis),
#                                                            !is.na(meanh),
#                                                            !is.na(meanage),
#                                                            !is.na(meancond)) %>%
#                                         mutate(across(c(meancort_propdiff,meanaftprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                       ~ scale(.x)[,1],
#                                                       .names = "{.col}_scaled")))
# summary(s_aftprovis_cort_propdiff_WEBL)
#
#
#
# ## TRES
#
# ### Overall provis
#
#
#
# s_provis_cort_s1_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_s1_scaled + meanprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                              family = binomial(link="logit"),
#                              data = dplyr::filter(s_with_provis_and_cort,
#                                                   Species == "TRES",
#                                                   !is.na(meanmaxt_nest),
#                                                   !is.na(meanmint_nest),
#                                                   !is.na(meancort_s1),
#                                                   !is.na(meanprovis),
#                                                   !is.na(meanh),
#                                                   !is.na(meanage),
#                                                   !is.na(meancond)) %>%
#                                mutate(across(c(meancort_s1,meanprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                              ~ scale(.x)[,1],
#                                              .names = "{.col}_scaled")))
# summary(s_provis_cort_s1_TRES)
#
#
#
#
# s_provis_cort_absdiff_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_absdiff_scaled + meanprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                   family = binomial(link="logit"),
#                                   data = dplyr::filter(s_with_provis_and_cort,
#                                                        Species == "TRES",
#                                                        !is.na(meanmaxt_nest),
#                                                        !is.na(meanmint_nest),
#                                                        !is.na(meancort_absdiff),
#                                                        !is.na(meanprovis),
#                                                        !is.na(meanh),
#                                                        !is.na(meanage),
#                                                        !is.na(meancond)) %>%
#                                     mutate(across(c(meancort_absdiff,meanprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                   ~ scale(.x)[,1],
#                                                   .names = "{.col}_scaled")))
# summary(s_provis_cort_absdiff_TRES)
#
#
#
#
#
# s_provis_cort_propdiff_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_propdiff_scaled + meanprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                    family = binomial(link="logit"),
#                                    data = dplyr::filter(s_with_provis_and_cort,
#                                                         Species == "TRES",
#                                                         !is.na(meanmaxt_nest),
#                                                         !is.na(meanmint_nest),
#                                                         !is.na(meancort_propdiff),
#                                                         !is.na(meanprovis),
#                                                         !is.na(meanh),
#                                                         !is.na(meanage),
#                                                         !is.na(meancond)) %>%
#                                      mutate(across(c(meancort_propdiff,meanprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                    ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled")))
# summary(s_provis_cort_propdiff_TRES)
#
#
#
# ### Morning cort
#
#
#
# s_mornprovis_cort_s1_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_s1_scaled + meanmornprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                  family = binomial(link="logit"),
#                                  data = dplyr::filter(s_with_provis_and_cort,
#                                                       Species == "TRES",
#                                                       !is.na(meanmaxt_nest),
#                                                       !is.na(meanmint_nest),
#                                                       !is.na(meancort_s1),
#                                                       !is.na(meanmornprovis),
#                                                       !is.na(meanh),
#                                                       !is.na(meanage),
#                                                       !is.na(meancond)) %>%
#                                    mutate(across(c(meancort_s1,meanmornprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                  ~ scale(.x)[,1],
#                                                  .names = "{.col}_scaled")))
# summary(s_mornprovis_cort_s1_TRES)
#
#
#
#
# s_mornprovis_cort_absdiff_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_absdiff_scaled + meanmornprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                       family = binomial(link="logit"),
#                                       data = dplyr::filter(s_with_provis_and_cort,
#                                                            Species == "TRES",
#                                                            !is.na(meanmaxt_nest),
#                                                            !is.na(meanmint_nest),
#                                                            !is.na(meancort_absdiff),
#                                                            !is.na(meanmornprovis),
#                                                            !is.na(meanh),
#                                                            !is.na(meanage),
#                                                            !is.na(meancond)) %>%
#                                         mutate(across(c(meancort_absdiff,meanmornprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                       ~ scale(.x)[,1],
#                                                       .names = "{.col}_scaled")))
# summary(s_mornprovis_cort_absdiff_TRES)
#
#
#
#
#
# s_mornprovis_cort_propdiff_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_propdiff_scaled + meanmornprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                        family = binomial(link="logit"),
#                                        data = dplyr::filter(s_with_provis_and_cort,
#                                                             Species == "TRES",
#                                                             !is.na(meanmaxt_nest),
#                                                             !is.na(meanmint_nest),
#                                                             !is.na(meancort_propdiff),
#                                                             !is.na(meanmornprovis),
#                                                             !is.na(meanh),
#                                                             !is.na(meanage),
#                                                             !is.na(meancond)) %>%
#                                          mutate(across(c(meancort_propdiff,meanmornprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                        ~ scale(.x)[,1],
#                                                        .names = "{.col}_scaled")))
# summary(s_mornprovis_cort_propdiff_TRES)
#
#
#
# ### Afternoon cort
#
#
#
# s_aftprovis_cort_s1_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_s1_scaled + meanaftprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                 family = binomial(link="logit"),
#                                 data = dplyr::filter(s_with_provis_and_cort,
#                                                      Species == "TRES",
#                                                      !is.na(meanmaxt_nest),
#                                                      !is.na(meanmint_nest),
#                                                      !is.na(meancort_s1),
#                                                      !is.na(meanaftprovis),
#                                                      !is.na(meanh),
#                                                      !is.na(meanage),
#                                                      !is.na(meancond)) %>%
#                                   mutate(across(c(meancort_s1,meanaftprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                 ~ scale(.x)[,1],
#                                                 .names = "{.col}_scaled")))
# summary(s_aftprovis_cort_s1_TRES)
#
#
#
#
# s_aftprovis_cort_absdiff_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_absdiff_scaled + meanaftprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                      family = binomial(link="logit"),
#                                      data = dplyr::filter(s_with_provis_and_cort,
#                                                           Species == "TRES",
#                                                           !is.na(meanmaxt_nest),
#                                                           !is.na(meanmint_nest),
#                                                           !is.na(meancort_absdiff),
#                                                           !is.na(meanaftprovis),
#                                                           !is.na(meanh),
#                                                           !is.na(meanage),
#                                                           !is.na(meancond)) %>%
#                                        mutate(across(c(meancort_absdiff,meanaftprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                      ~ scale(.x)[,1],
#                                                      .names = "{.col}_scaled")))
# summary(s_aftprovis_cort_absdiff_TRES)
#
#
#
#
#
# s_aftprovis_cort_propdiff_TRES <- glm(cbind(nest_fledged,brood_size - nest_fledged) ~ meancort_propdiff_scaled + meanaftprovis_scaled + meanmaxt_nest_scaled + meanmint_nest_scaled + meanh_scaled + meanage_scaled + meancond_scaled,
#                                       family = binomial(link="logit"),
#                                       data = dplyr::filter(s_with_provis_and_cort,
#                                                            Species == "TRES",
#                                                            !is.na(meanmaxt_nest),
#                                                            !is.na(meanmint_nest),
#                                                            !is.na(meancort_propdiff),
#                                                            !is.na(meanaftprovis),
#                                                            !is.na(meanh),
#                                                            !is.na(meanage),
#                                                            !is.na(meancond)) %>%
#                                         mutate(across(c(meancort_propdiff,meanaftprovis,meanmaxt_nest,meanmint_nest,meanh,meanage,meancond),
#                                                       ~ scale(.x)[,1],
#                                                       .names = "{.col}_scaled")))
# summary(s_aftprovis_cort_propdiff_TRES)




# Summary chart of findings?
