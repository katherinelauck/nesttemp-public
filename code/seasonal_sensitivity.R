## ----setup, include=FALSE-----------------------------------------------------
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

# g <- read_rds("data/growth_and_provis_combined_mobilenetv3-original_dataset.h5.rds") %>%
#   mutate(abs_change_cort = cort_s2-cort_s1,
#          prop_change_cort = (cort_s2-cort_s1)/cort_s1,
#          year_fct = as.factor(year))
p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))


g <- read_rds("data/growth_cort_provis_manytempmeasures.rds")

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

p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))



## -----------------------------------------------------------------------------

prior_model <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_webl <- bind_rows(c1,c2) %>%
  as.tibble() %>%
  mutate(across(where(is.numeric),~round(.x,digits = 4)),
         P = `Pr(>Chisq)`) %>%
  dplyr::select(Model,AIC,Chisq,P) %>%
  mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
         across(c(AIC,Chisq), ~ round(.x, digits = 2)),
         P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
  group_by(max_or_min) %>% gt())

summary(g_lintemp)
check_collinearity(g_lintemp_noint)

anova(g_lintemp,prior_model)

g_lintemp_webl <- g_lintemp



## -----------------------------------------------------------------------------
weblgrowthtrendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
   mutate(Model = "No temp * juliandate int") %>%
  gt()

weblgrowthtrendmax <- emtrends(g_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    mutate(Model = "interaction") %>%
  gt()




## -----------------------------------------------------------------------------

samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())



dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





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

(fig2_webl <- predict_response(prior_model,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
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
   #ylim(-5,5) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )




## -----------------------------------------------------------------------------
(fig2_webl <- predict_response(g_lintemp,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
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
   #ylim(-5,5) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )




## -----------------------------------------------------------------------------

(plott <- predict_response(g_lintemp,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    # labs(title = element_blank()) +
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
    ) #+
    #ylim(-5,5) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------

prior_model <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_tres <- bind_rows(c1,c2) %>%
  as.tibble() %>%
  mutate(across(where(is.numeric),~round(.x,digits = 4)),
         P = `Pr(>Chisq)`) %>%
  dplyr::select(Model,AIC,Chisq,P) %>%
  mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
         across(c(AIC,Chisq), ~ round(.x, digits = 2)),
         P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
  group_by(max_or_min) %>% gt())

check_collinearity(g_lintemp_noint)

summary(g_lintemp_noint)

g_lintemp_addmin_tres <- g_lintemp_addmin

(anova(g_lintemp_addmin,prior_model))



## -----------------------------------------------------------------------------
tresgrowthtrendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
   mutate(Model = "No temp * juliandate int") %>%
  gt()

tresgrowthtrendmax <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    mutate(Model = "interaction") %>%
  gt()

(trend_comp <- rbind(tresgrowthtrendmax_prior$`_data`, tresgrowthtrendmax$`_data`) %>% gt())



## -----------------------------------------------------------------------------
samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

data_tres = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_tres <- mean(data_tres %>% pull(meanmaxtempI))
sd_temp_tres <- sd(data_tres %>% pull(meanmaxtempI))


temp_trans_tres <- trans_new("temp_trans",
                          transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                          inverse = function(x){x})


(fig2_tres <- predict_response(prior_model,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
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
   #ylim(-5,5) +
   #geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )




## -----------------------------------------------------------------------------
(fig2_tres <- predict_response(g_lintemp_addmin,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
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
   #ylim(-5,5) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )




## -----------------------------------------------------------------------------

(plott <- predict_response(g_lintemp_addmin,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    # labs(title = element_blank()) +
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
    ) #+
    #ylim(-5,5) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------

prior_model <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
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

s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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

s_nestpd_WEBL_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled * habitat + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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

s_nestpd_WEBL_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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
s_nestpd_WEBL_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled + habitat + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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



#Conclusion: For nest attempt overall, land cover interacts with either max or min temp but not both together. It looks like the max temp interaction model is slightly more explanatory so we'll go with that.


summary(s_nestpd_WEBL_addmin)

anova(s_nestpd_WEBL_addmin, prior_model)


## -----------------------------------------------------------------------------
(weblsurvival_nestpd_trendmax_prior <- emtrends(prior_model,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
   gt())

(weblsurvival_nestpd_trendmax <- emtrends(s_nestpd_WEBL_addmin,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
   gt())


## -----------------------------------------------------------------------------
data_webl = s_nestpd_WEBL_addmin$data


mean_temp_webl <- mean(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                          transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                          inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_webl.html")


dat_text_webl <- data.frame(
  label = paste("N =",samp_webl$count),
  group   = factor(samp_webl$habitat)
)

(fig3_webl <- predict_response(prior_model,terms = c("meanmaxt_nestpd_scaled [all]","habitat")) %>%
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


## -----------------------------------------------------------------------------
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
   scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


## -----------------------------------------------------------------------------
(plott <- predict_response(s_nestpd_WEBL_addmin,terms = c("meanmaxt_nestpd_scaled [all]","juliandate_hatch_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Survival") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    # labs(title = element_blank()) +
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
    ) #+
    #ylim(-5,5) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------

prior_model <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
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

s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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

s_nestpd_TRES_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled * habitat + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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

s_nestpd_TRES_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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
s_nestpd_TRES_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled + habitat + meanmaxt_nestpd_scaled * juliandate_hatch_scaled + year_fct + site,
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



#Conclusion: For nest attempt overall, land cover interacts with either max or min temp but not both together. It looks like the max temp interaction model is slightly more explanatory so we'll go with that.


summary(s_nestpd_TRES_noint)

anova(s_nestpd_TRES_noint, prior_model)


## -----------------------------------------------------------------------------
(tressurvival_nestpd_trendmax_prior <- emtrends(prior_model,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
   gt())

(tressurvival_nestpd_trendmax <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
   gt())


## -----------------------------------------------------------------------------
data_tres = s_nestpd_TRES_noint$data


mean_temp_tres <- mean(data_tres %>% pull(meanmaxt_nestpd),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(meanmaxt_nestpd),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                          transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                          inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_tres.html")


dat_text_tres <- data.frame(
  label = paste("N =",samp_tres$count),
  group   = factor(samp_tres$habitat)
)

(fig3_tres <- predict_response(prior_model,terms = c("meanmaxt_nestpd_scaled [all]","habitat")) %>%
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
   scale_x_continuous(trans = temp_trans_tres,
                       breaks = c((20-mean_temp_tres)/sd_temp_tres,
                                  (25-mean_temp_tres)/sd_temp_tres,
                                  (30-mean_temp_tres)/sd_temp_tres,
                                  (35-mean_temp_tres)/sd_temp_tres,
                                  (40-mean_temp_tres)/sd_temp_tres#,
                                  #(45-mean_temp)/sd_temp
                                  ),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
                       ) +
   # ylim(0,100) +
   scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


## -----------------------------------------------------------------------------
(fig3_tres <- predict_response(s_nestpd_TRES_noint,terms = c("meanmaxt_nestpd_scaled [all]","habitat")) %>%
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
   scale_x_continuous(trans = temp_trans_tres,
                       breaks = c((20-mean_temp_tres)/sd_temp_tres,
                                  (25-mean_temp_tres)/sd_temp_tres,
                                  (30-mean_temp_tres)/sd_temp_tres,
                                  (35-mean_temp_tres)/sd_temp_tres,
                                  (40-mean_temp_tres)/sd_temp_tres#,
                                  #(45-mean_temp)/sd_temp
                                  ),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
                       ) +
   # ylim(0,100) +
   scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


## -----------------------------------------------------------------------------
(plott <- predict_response(s_nestpd_TRES_noint,terms = c("meanmaxt_nestpd_scaled [all]","juliandate_hatch_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Survival") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    # labs(title = element_blank()) +
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
    ) #+
    #ylim(-5,5) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
g <- read_rds("data/growth_cort_provis_manytempmeasures.rds") %>%
  mutate(year_fct = as.factor(year))

prior_model <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + habitat + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
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

anova(s1_lintemp_noint,prior_model)


## -----------------------------------------------------------------------------
(webls1trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())

(webls1trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s1_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
# samp_s1_webl %>% gt() %>% gtsave("figures/ss_webl_s1.html")


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

(figs2_webl <- ggpredict(s1_lintemp_noint,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    #scale_fill_viridis(discrete = TRUE) +
    #scale_color_viridis(discrete = TRUE) +
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
    ) #+
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    #geom_text(data = dat_text_s1_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    #theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
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
#gtsave(int_tab_s1_tres,"figures/int_tab_s1_tres.html")

s1_tres <- s1_lintemp_addmax

anova(s1_lintemp_addmax,prior_model)


## -----------------------------------------------------------------------------
(tress1trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(tress1trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s1_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_tres_s1.html")


dat_text_s1_tres <- data.frame(
  label = paste("N =",samp_s1_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

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

(figs2_tres <- ggpredict(prior_model,terms = c("meanmaxtempI_scaled [all]"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    # geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
samp_s1_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_tres_s1.html")


dat_text_s1_tres <- data.frame(
  label = paste("N =",samp_s1_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

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

(figs2_tres <- ggpredict(s1_lintemp_addmax,terms = c("meanmaxtempI_scaled [all]"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    # geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(figs2_tres <- ggpredict(s1_lintemp_addmax,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) #+
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    # geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    #theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
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
#gtsave(int_tab_abs_webl,"figures/int_tab_abs_webl.html")

check_collinearity(abs_lintemp_noint)

abs_webl <- abs_lintemp
summary(abs_lintemp)

anova(abs_lintemp,prior_model)


## -----------------------------------------------------------------------------
(weblabstrendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(weblabstrendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_webl_abs.html")


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

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

(fig4_webl <- ggpredict(prior_model,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------

(fig4_webl <- ggpredict(abs_lintemp,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    # ylim(0,60) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

prior_model_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                       mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
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

summary(abs_lintemp)


anova(abs_lintemp,prior_model)
summary(abs_lintemp_addmin)

anova(abs_lintemp_addmin,prior_model_addmin)




## -----------------------------------------------------------------------------
(tresabstrendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())

(tresabstrendmax <- emtrends(prior_model_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())

(tresabstrendmax <- emtrends(abs_lintemp_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_tres_abs.html")


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

(fig4_tres <- ggpredict(prior_model,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
(fig4_tres <- ggpredict(abs_lintemp_addmin,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
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
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig4_tres <- ggpredict(abs_lintemp_addmin,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    #ylim(-5,5) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
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

summary(s2_lintemp)

anova(s2_lintemp,prior_model)


## -----------------------------------------------------------------------------
(webls2trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())

(webls2trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s2_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
# samp_s2_webl %>% gt() %>% gtsave("figures/ss_webl_s2.html")


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

(fig_webl_s2 <- ggpredict(prior_model,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
(fig_webl_s2 <- ggpredict(s2_lintemp,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s2_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
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

summary(s2_lintemp_addmin)

anova(s2_lintemp_addmin,prior_model)


## -----------------------------------------------------------------------------
(tress2trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(tress2trendmax <- emtrends(s2_lintemp_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s2_tres <- s2_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_tres_s2.html")


dat_text_s2_tres <- data.frame(
  label = paste("N =",samp_s2_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


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

(fig_tres_s2 <- ggpredict(prior_model,terms = c("meanmaxtempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
(fig_tres_s2 <- ggpredict(s2_lintemp_addmin,terms = c("meanmaxtempI_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    #ylim(-5,5) +
    # geom_text(data = dat_text_s2_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------

prior_model <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
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

summary(s1_lintemp_addmax)

anova(s1_lintemp_addmax,prior_model)


## -----------------------------------------------------------------------------
(webls1_priordayt_trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(webls1_priordayt_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s1_priordayt_webl <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
# samp_s1_priordayt_webl %>% gt() %>% gtsave("figures/ss_webl_s1_priordayt_.html")


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

(figs2_priordayt_webl <- ggpredict(prior_model,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
(figs2_priordayt_webl <- ggpredict(s1_lintemp_addmax,terms = c("maxt_prior_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s1_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s1_lintemp_noint <- lmerTest::lmer(sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
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

summary(s1_lintemp_addmin)

anova(s1_lintemp_addmax,prior_model)


## -----------------------------------------------------------------------------
(tress1_priordayt_trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(tress1_priordayt_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s1_priordayt_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
# samp_s1_priordayt_tres %>% gt() %>% gtsave("figures/ss_tres_s1_priordayt_.html")


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

(figs2_priordayt_tres <- ggpredict(prior_model,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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
    scale_linetype_manual(values = c("solid","solid","solid","solid")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(figs2_priordayt_tres <- ggpredict(s1_lintemp_addmax,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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
    scale_linetype_manual(values = c("solid","solid","solid","solid")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(figs2_priordayt_tres <- ggpredict(s1_lintemp_addmax,terms = c("maxt_prior_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
#     scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
#     #ylim(-5,5) +
#     geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
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

summary(abs_lintemp_noint)

abs_webl_priordayt <- abs_lintemp_noint
anova(abs_lintemp_noint,prior_model)


## -----------------------------------------------------------------------------
(weblabs_priordayt_trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(weblabs_priordayt_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_webl_abs_priordayt.html")


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

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

(fig4_webl_priordayt <- ggpredict(prior_model,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
(fig4_webl_priordayt <- ggpredict(abs_lintemp_noint,terms = c("maxt_prior_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced - Baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # # ylim(0,60) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +  maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                              ~ scale(.x)[,1],
                                              .names = "{.col}_scaled"),
                                       maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmax <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled +  maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_addmin <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled +  maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                       mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                     ~ scale(.x)[,1],
                                                     .names = "{.col}_scaled"),
                                              maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

abs_lintemp_noint <- lmerTest::lmer(sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled +  maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
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

summary(abs_lintemp)

abs_tres_priordayt <- abs_lintemp

anova(abs_lintemp,prior_model)


## -----------------------------------------------------------------------------
(tresabs_priordayt_trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(tresabs_priordayt_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_tres_abs_priordayt.html")


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

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

(fig4_tres_priordayt <- ggpredict(prior_model,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
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
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,60) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig4_tres_priordayt <- ggpredict(abs_lintemp,terms = c("maxt_prior_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced - baseline corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    # # ylim(0,60) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                     mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxt_prior),!is.na(mint_prior)) %>%
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

summary(s2_lintemp_noint)

anova(s2_lintemp_noint,prior_model)


## -----------------------------------------------------------------------------
(webls2_priordayt_trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(webls2_priordayt_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s2_priordayt_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
# samp_s2_priordayt_webl %>% gt() %>% gtsave("figures/ss_webl_s2_priordayt_.html")


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

(fig_priordayt_webl_s2 <- ggpredict(prior_model,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
(fig_priordayt_webl_s2 <- ggpredict(s2_lintemp_addmax,terms = c("maxt_prior_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s2_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                               mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmax <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_addmin <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
                                      mutate(across(c(gweight,maxt_prior,mint_prior,juliandate,brood_size,age),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             maxt_prior_scaled_sq = maxt_prior_scaled * maxt_prior_scaled))

s2_lintemp_noint <- lmerTest::lmer(sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxt_prior),!is.na(mint_prior)) %>%
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

summary(s2_lintemp)

anova(s2_lintemp,prior_model)


## -----------------------------------------------------------------------------
(tress2_priordayt_trendmax_prior <- emtrends(prior_model,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())
(tress2_priordayt_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## -----------------------------------------------------------------------------
samp_s2_priordayt_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
# samp_s2_priordayt_tres %>% gt() %>% gtsave("figures/ss_tres_s2_priordayt_.html")


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

(fig_priordayt_tres_s2 <- ggpredict(prior_model,terms = c("maxt_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
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


## -----------------------------------------------------------------------------
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
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig_priordayt_tres_s2 <- ggpredict(s2_lintemp,terms = c("maxt_prior_scaled [all]","juliandate_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Max temp over prior day (\u00b0C)") +
    ylab("Stress-induced corticosterone (ng/mL)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dashed")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------

prior_model <- glmmTMB(valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(p,!is.na(mean_temp),
                                         !is.na(julian_date),
                                         !is.na(mean_nestling_age),
                                         !is.na(tod_h),
                                         !is.na(site),
                                         !is.na(attempt_id),
                                         !is.na(year),
                                         species == "WEBL"# ,
                                         # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                         # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                         # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                    ) %>%
                      mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                                    ~ scale(.x)[,1],
                                    .names = "{.col}_scaled")))

m <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) * habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(p,!is.na(mean_temp),
                                  !is.na(julian_date),
                                  !is.na(mean_nestling_age),
                                  !is.na(tod_h),
                                  !is.na(site),
                                  !is.na(attempt_id),
                                  !is.na(year),
                                  species == "WEBL"# ,
                                  # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                  # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                  # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
             ) %>%
               mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                             ~ scale(.x)[,1],
                             .names = "{.col}_scaled")))

emmeans(m,specs = pairwise ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
emmip(m,formula = habitat ~ mean_temp_scaled, at = list(mean_temp_scaled = seq(from = -2.5, to = 2.5, by = .1)))

m_linint <- glmmTMB(valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(p,!is.na(mean_temp),
                                         !is.na(julian_date),
                                         !is.na(mean_nestling_age),
                                         !is.na(tod_h),
                                         !is.na(site),
                                         !is.na(attempt_id),
                                         !is.na(year),
                                         species == "WEBL"# ,
                                         # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                         # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                         # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                    ) %>%
                      mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                                    ~ scale(.x)[,1],
                                    .names = "{.col}_scaled")))

m_noint <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) + habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(p,!is.na(mean_temp),
                                        !is.na(julian_date),
                                        !is.na(mean_nestling_age),
                                        !is.na(tod_h),
                                        !is.na(site),
                                        !is.na(attempt_id),
                                        !is.na(year),
                                        species == "WEBL"# ,
                                        # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                        # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                        # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                   ) %>%
                     mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                                   ~ scale(.x)[,1],
                                   .names = "{.col}_scaled")))
c <- anova(m,m_linint,m_noint) %>% tibble() %>% mutate(Model = c("No temp * LC interaction","linear temp * LC interaction", "quadratic temp * LC  interaction"),.before = Df)

(int_tab_provis_webl <- c %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>% gt())
# gtsave(int_tab_provis_webl,"figures/int_tab_provis_webl.html")

summary(m_linint)
check_collinearity(m_noint)

provis_webl <- m_linint

summary(m_linint)

anova(m_linint,prior_model)


## -----------------------------------------------------------------------------
(weblprovistrend_prior <- emtrends(prior_model,specs = ~ habitat, var = c("mean_temp_scaled"),max.degree = 1) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    #group_by(degree) %>%
    gt())
(weblprovistrend <- emtrends(m_linint,specs = ~ habitat, var = c("mean_temp_scaled"),max.degree = 1) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    #group_by(degree) %>%
    gt())


## -----------------------------------------------------------------------------
samp <- m_linint$frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_webl_provis.html")


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

data_webl = dplyr::filter(p,!is.na(mean_temp),
                          !is.na(julian_date),
                          !is.na(mean_nestling_age),
                          !is.na(tod_h),
                          !is.na(site),
                          !is.na(attempt_id),
                          !is.na(year),
                          species == "WEBL"# ,
                          # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                          # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                          # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
) %>%
  mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


mean_temp_webl <- mean(data_webl %>% pull(mean_temp))
sd_temp_webl <- sd(data_webl %>% pull(mean_temp))


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

(fig5_webl <- predict_response(prior_model,terms = c("mean_temp_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean temperature (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title =  element_blank()) +
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
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig5_webl <- predict_response(m_linint,terms = c("mean_temp_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean temperature (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title =  element_blank()) +
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
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig5_webl <- predict_response(m_linint,terms = c("mean_temp_scaled [all]","julian_date_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean temperature (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title =  element_blank()) +
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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
prior_model <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(p,!is.na(mean_temp),
                                  !is.na(julian_date),
                                  !is.na(mean_nestling_age),
                                  !is.na(tod_h),
                                  !is.na(site),
                                  !is.na(attempt_id),
                                  !is.na(year),
                                  valid_detections < 40,
                                  species == "TRES"# ,
                                  # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                  # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                  # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
             ) %>%
               mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                             ~ scale(.x)[,1],
                             .names = "{.col}_scaled")))

m <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) * habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(p,!is.na(mean_temp),
                                  !is.na(julian_date),
                                  !is.na(mean_nestling_age),
                                  !is.na(tod_h),
                                  !is.na(site),
                                  !is.na(attempt_id),
                                  !is.na(year),
                                  valid_detections < 40,
                                  species == "TRES"# ,
                                  # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                  # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                  # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
             ) %>%
               mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                             ~ scale(.x)[,1],
                             .names = "{.col}_scaled")))

m_linint <- glmmTMB(valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(p,!is.na(mean_temp),
                                         !is.na(julian_date),
                                         !is.na(mean_nestling_age),
                                         !is.na(tod_h),
                                         !is.na(site),
                                         !is.na(attempt_id),
                                         !is.na(year),
                                         valid_detections < 40,
                                         species == "TRES"# ,
                                         # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                         # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                         # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                    ) %>%
                      mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                                    ~ scale(.x)[,1],
                                    .names = "{.col}_scaled")))

m_noint <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) + habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(p,!is.na(mean_temp),
                                        !is.na(julian_date),
                                        !is.na(mean_nestling_age),
                                        !is.na(tod_h),
                                        !is.na(site),
                                        !is.na(attempt_id),
                                        !is.na(year),
                                        valid_detections < 40,
                                        species == "TRES"# ,
                                        # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                        # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                        # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                   ) %>%
                     mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                                   ~ scale(.x)[,1],
                                   .names = "{.col}_scaled")))

c <- anova(m,m_linint,m_noint) %>% tibble() %>% mutate(Model = c("No temp * LC interaction","Linear temp * LC interaction","Squared temp * LC interaction"),.before = Df)

(int_tab_provis_tres <- c %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>% gt())
# gtsave(int_tab_provis_tres,"figures/int_tab_provis_tres.html")

summary(m)
check_collinearity(m)


provis_tres <- m

summary(m)
anova(m,prior_model)


## -----------------------------------------------------------------------------
(tresprovistrend_prior <- emtrends(prior_model,specs = ~ degree | habitat, var = "mean_temp_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())
(tresprovistrend <- emtrends(m,specs = ~ degree | habitat, var = "mean_temp_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


## -----------------------------------------------------------------------------
samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("figures/ss_tres_provis.html")


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_tres = dplyr::filter(p,!is.na(mean_temp),
                          !is.na(julian_date),
                          !is.na(mean_nestling_age),
                          !is.na(tod_h),
                          !is.na(site),
                          !is.na(attempt_id),
                          !is.na(year),
                          valid_detections < 40,
                          species == "TRES"# ,
                          # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                          # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                          # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
) %>%
  mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


mean_temp_tres <- mean(data_tres %>% pull(mean_temp))
sd_temp_tres <- sd(data_tres %>% pull(mean_temp))


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

(fig5_tres <- ggpredict(prior_model,terms = c("mean_temp_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean temperature (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES provis: maxsq interaction with LC, but not with linear max") +
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
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    ylim(0,50) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig5_tres <- ggpredict(m,terms = c("mean_temp_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean temperature (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES provis: maxsq interaction with LC, but not with linear max") +
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
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    ylim(0,50) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
(fig5_tres <- ggpredict(prior_model,terms = c("mean_temp_scaled [all]","julian_date_scaled"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    # aes(linetype = .data[["group"]]) +
    theme_classic() +
    # facet_wrap(~ group, ncol = 2) +
    xlab("Mean temperature (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    # scale_fill_viridis(discrete = TRUE) +
    # scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = "TRES provis: maxsq interaction with LC, but not with linear max") +
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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # ylim(0,50) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


## -----------------------------------------------------------------------------
growth <- rbind(weblgrowthtrendmax_prior$`_data`,tresgrowthtrendmax_prior$`_data`,
                weblgrowthtrendmax$`_data`,tresgrowthtrendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Growth ~ Mean max temperature over prior week",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

survival <- rbind(weblsurvival_nestpd_trendmax_prior$`_data`,tressurvival_nestpd_trendmax_prior$`_data`,
                weblsurvival_nestpd_trendmax$`_data`,tressurvival_nestpd_trendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Survival ~ Mean max temperature over prior week",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

s1_priorday <- rbind(webls1_priordayt_trendmax_prior$`_data`,tress1_priordayt_trendmax_prior$`_data`,
                webls1_priordayt_trendmax$`_data`,tress1_priordayt_trendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Baseline corticosterone ~ Maximum temperature of the prior day",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

s1_priorweek <- rbind(webls1trendmax_prior$`_data`,tress1trendmax_prior$`_data`,
                webls1trendmax$`_data`,tress1trendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Baseline corticosterone ~ Mean max temperature of the prior week",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

abs_priorday <- rbind(weblabs_priordayt_trendmax_prior$`_data`,tresabs_priordayt_trendmax_prior$`_data`,
                weblabs_priordayt_trendmax$`_data`,tresabs_priordayt_trendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Difference between stress-induced and baseline corticosterone ~ Maximum temperature of the prior day",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

abs_priorweek <- rbind(weblabstrendmax_prior$`_data`,tresabstrendmax_prior$`_data`,
                weblabstrendmax$`_data`,tresabstrendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Difference between stress-induced and baseline corticosterone ~ Mean max temperature of the prior week",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

s2_priorday <- rbind(webls2_priordayt_trendmax_prior$`_data`,tress2_priordayt_trendmax_prior$`_data`,
                webls2_priordayt_trendmax$`_data`,tress2_priordayt_trendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Stres-induced corticosterone ~ Maximum temperature of the prior day",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

s2_priorweek <- rbind(webls2trendmax_prior$`_data`,tress2trendmax_prior$`_data`,
                webls2trendmax$`_data`,tress2trendmax$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "T-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",4)),2),
         Response = "Stress-induced corticosterone ~ Mean max temperature of the prior week",
         Model = c(rep("No maximum temperature * day of year interaction",8),
                       rep("Maximum temperature * day of year interaction",8)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)

provis <- rbind(weblprovistrend_prior$`_data` %>% mutate(degree = "linear"),tresprovistrend_prior$`_data`,
                weblprovistrend$`_data` %>% mutate(degree = "linear"),tresprovistrend$`_data`) %>%
  rename(Trend = "Max temp trend",`Statistic` = "Z-ratio",`P-value` = "P") %>%
  mutate(Species = rep(c(rep("Western Bluebird",4),rep("Tree Swallow",8)),2),
         Response = "Provisioning ~ Average hourly temperature",
         Model = c(rep("No maximum temperature * day of year interaction",12),
                       rep("Maximum temperature * day of year interaction",12)),
         sig = if_else(`P-value` == "<0.001","***",
                       if_else(as.numeric(`P-value`) < .01, "**",
                               if_else(as.numeric(`P-value`) < .05,"*",""))),
         Estimate = paste0(as.character(Trend),sig),
         Predictor = "Mean max temperature over prior week") %>%
  dplyr::select(Species, Response, Model, Habitat, Estimate, degree) %>%
  pivot_wider(names_from = degree,values_from = Estimate) %>%
  mutate(quadratic = replace_na(quadratic, "NA"),
         Estimate = paste0(linear,", ",quadratic)) %>%
  dplyr::select(-c(linear,quadratic)) %>%
  group_by(Response) %>%
  pivot_longer(-c(Species, Response, Model,Habitat)) %>%
  pivot_wider(names_from=c(Species, Habitat), values_from=value) %>%
  dplyr::select(-name)



## -----------------------------------------------------------------------------
(trend <- rbind(growth,survival,s1_priorday,s1_priorweek,abs_priorday,abs_priorweek,s2_priorday,s2_priorweek,provis) %>% gt(rowname_col = "Model",groupname_col = "Response") %>% tab_options(data_row.padding = px(1)) %>%

  tab_spanner_delim(
    delim="_"))

# gtsave(trend,"../figures/val_tmaxxjday_tbl.html")


save(list = ls(), file = "data/models_seasonal.RData")
