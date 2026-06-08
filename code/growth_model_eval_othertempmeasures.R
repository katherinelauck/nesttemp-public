## Create all objects needed for growth analysis using maxhi_week & cumulative week measures to predict growth

#knitr::opts_chunk$set(echo = TRUE,warning = FALSE,message = FALSE)
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
library(egg)
library(igraph)
library(tidyverse)
library(performance)
library(modelsummary)
library(car)
library(lme4)
library(hms)
library(emmeans)
library(patchwork)
# modelsummary::get_gof(my_model)
# performance::model_performance(my_model)
# wmean <- read_rds("data/wmean.rds")  # unused; data pre-computed in growth_cort_provis_manytempmeasures.rds

# g <- read_rds("data/growth_and_provis_combined_mobilenetv3-original_dataset.h5.rds") %>%
#   mutate(abs_change_cort = cort_s2-cort_s1,
#          prop_change_cort = (cort_s2-cort_s1)/cort_s1,
#          year_fct = as.factor(year))
p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

# b <- read_csv("data/banding-and-morphometrics_proofed.csv")  # unused


g <- read_rds("data/growth_cort_provis_manytempmeasures.rds")


### WEBL
#### maxhi_week


g_lintemp <- lmerTest::lmer(gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_maxhiweek_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_maxhiweek_webl,"figures/int_tab_growth_maxhiweek_webl.html")
summary(g_lintemp)
check_collinearity(g_lintemp_noint)

g_lintemp_maxhiweek_webl <- g_lintemp

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_maxhiweek_webl <- g_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_maxhiweek_webl <- samp_year_maxhiweek_webl %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_maxhiweek_webl

samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_maxhiweek_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_maxhiweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_maxhiweek_webl <- mean(data_maxhiweek_webl %>% pull(maxhi_week))
sd_temp_maxhiweek_webl <- sd(data_maxhiweek_webl %>% pull(maxhi_week))


temp_trans_maxhiweek_webl <- trans_new("temp_trans_maxhiweek_webl",
                             transform = function(x){(x * sd_temp_maxhiweek_webl) + mean_temp_maxhiweek_webl},
                             inverse = function(x){x})

(fig2_maxhiweek_webl <- predict_response(g_lintemp,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max heat index (\u00b0C) over preceding week") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_maxhiweek_webl,
                       breaks = c((20-mean_temp_maxhiweek_webl)/sd_temp_maxhiweek_webl,
                                  (25-mean_temp_maxhiweek_webl)/sd_temp_maxhiweek_webl,
                                  (30-mean_temp_maxhiweek_webl)/sd_temp_maxhiweek_webl,
                                  (35-mean_temp_maxhiweek_webl)/sd_temp_maxhiweek_webl,
                                  (40-mean_temp_maxhiweek_webl)/sd_temp_maxhiweek_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$maxhi_week,na.rm = TRUE))/sd(g$maxhi_week,na.rm = TRUE),
                       #            (55-mean(g$maxhi_week,na.rm = TRUE))/sd(g$maxhi_week,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_maxhiweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbymaxhiweekxhab_WEBL.png",plot =  fig2_maxhiweek_webl, width = 10, height = 6.6)



summary(g_lintemp)



(weblgrowthtrendmaxhiweek <- emtrends(g_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblgrowthtrendmaxhiweek,"figures/weblgrowthtrendmaxhiweek.html")

data_maxhiweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(maxhi_week_scaled = (maxhi_week_scaled * sd(data$maxhi_week)) + mean(data$maxhi_week),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     maxhi_week_scaled = round(maxhi_week_scaled),
#     maxhi_week_scaled = paste0(maxhi_week_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "maxhi_week_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/weblgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ maxhi_week_scaled, at = list(maxhi_week_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiweek_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_maxhiweek_webl,"figures/growthbyhabitat_maxhiweek_webl.html")


### Check for effect of temperature


(growthbyhabitat_summary_maxhiweek_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_maxhiweek_webl,"figures/growthbyhabitat_summary_maxhiweek_webl.html")



#### Check whether controlling for brood size explains the habitat * temp differences.

# Because higher temps reduce survival, especially in forest, we are thinking that brood size is the outcome of the habitat * temp interaction. So, if we include brood size and the relationship growth ~ habitat * temp disappears, we think that's what's going on.

## create brood size by week measure

#
# brood_size <- b %>% group_by(`Banding Date`,Nestbox) %>% summarize(brood_size_weekly = n()) %>% ungroup() %>% mutate(date = mdy(`Banding Date`)) %>% dplyr::select(!`Banding Date`)
#
# g <- left_join(g %>% ungroup(),brood_size, by = c("date","Nestbox"))
#
#
# # If this is so, then the main driver of growth and survival may be predation.
#
#
# g_lintemp <- lmerTest::lmer(gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + brood_size_weekly_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
# g_lintemp_addmax <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled * habitat + brood_size_weekly_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
# g_lintemp_addmin <- lmerTest::lmer(gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled + brood_size_weekly_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
# g_lintemp_noint <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled + habitat + brood_size_weekly_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
# c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)
#
# c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)
#
# (tab <- bind_rows(c1,c2) %>%
#   as.tibble() %>%
#   mutate(across(where(is.numeric),~round(.x,digits = 4)),
#          P = `Pr(>Chisq)`) %>%
#   dplyr::select(Model,AIC,Chisq,P) %>%
#   mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
#          across(c(AIC,Chisq), ~ round(.x, digits = 2)),
#          P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
#   group_by(max_or_min) %>% gt())
# gtsave(tab,"../figures/int_tab_growth_with_broodsize_webl.html")
#
# summary(g_lintemp)
#
#
# # Is brood_size significant if I remove temp and habitat?
#
#
# g_lintemp <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + brood_size_weekly_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
#
#
# summary(g_lintemp)
#
# g_lintemp <- lmerTest::lmer(gweight ~ age_scaled + brood_size_weekly_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
#
#
# summary(g_lintemp)

#### maxhi_prior


g_lintemp <- lmerTest::lmer(gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ maxhi_prior_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ maxhi_prior_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_maxhiday_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_maxhiday_webl,"figures/int_tab_growth_maxhiday_webl.html")
summary(g_lintemp)
check_collinearity(g_lintemp_noint)

g_lintemp_maxhiday_webl <- g_lintemp

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_maxhiday_webl <- g_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_maxhiday_webl <- samp_year_maxhiday_webl %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_maxhiday_webl

samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_maxhiday_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_maxhiday_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_maxhiday_webl <- mean(data_maxhiday_webl %>% pull(maxhi_prior))
sd_temp_maxhiday_webl <- sd(data_maxhiday_webl %>% pull(maxhi_prior))


temp_trans_maxhiday_webl <- trans_new("temp_trans_maxhiday_webl",
                                       transform = function(x){(x * sd_temp_maxhiday_webl) + mean_temp_maxhiday_webl},
                                       inverse = function(x){x})

(fig2_maxhiday_webl <- predict_response(g_lintemp,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max heat index (\u00b0C) over preceding day") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "solid","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_maxhiday_webl,
                       breaks = c((20-mean_temp_maxhiday_webl)/sd_temp_maxhiday_webl,
                                  (25-mean_temp_maxhiday_webl)/sd_temp_maxhiday_webl,
                                  (30-mean_temp_maxhiday_webl)/sd_temp_maxhiday_webl,
                                  (35-mean_temp_maxhiday_webl)/sd_temp_maxhiday_webl,
                                  (40-mean_temp_maxhiday_webl)/sd_temp_maxhiday_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$maxhi_prior,na.rm = TRUE))/sd(g$maxhi_prior,na.rm = TRUE),
                       #            (55-mean(g$maxhi_prior,na.rm = TRUE))/sd(g$maxhi_prior,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_maxhiday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbymaxhidayxhab_WEBL.png",plot =  fig2_maxhiday_webl, width = 10, height = 6.6)



summary(g_lintemp)



(weblgrowthtrendmaxhiday <- emtrends(g_lintemp,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblgrowthtrendmaxhiday,"figures/weblgrowthtrendmaxhiday.html")

data_maxhiday_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(maxhi_prior_scaled = (maxhi_prior_scaled * sd(data$maxhi_prior)) + mean(data$maxhi_prior),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     maxhi_prior_scaled = round(maxhi_prior_scaled),
#     maxhi_prior_scaled = paste0(maxhi_prior_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "maxhi_prior_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/weblgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ maxhi_prior_scaled, at = list(maxhi_prior_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiday_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_maxhiday_webl,"figures/growthbyhabitat_maxhiday_webl.html")


### Check for effect of temperature


(growthbyhabitat_summary_maxhiday_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_maxhiday_webl,"figures/growthbyhabitat_summary_maxhiday_webl.html")


#### degreehours_over_30C_priorweek


g_lintemp <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(`degreehours_over_30C_priorweek`),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_deghr30week_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_deghr30week_webl,"figures/int_tab_growth_deghr30week_webl.html")
summary(g_lintemp)
check_collinearity(g_lintemp_noint)

g_lintemp_deghr30week_webl <- g_lintemp

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_deghr30week_webl <- g_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_deghr30week_webl <- samp_year_deghr30week_webl %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_deghr30week_webl

samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_deghr30week_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_deghr30week_webl = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_deghr30week_webl <- mean(data_deghr30week_webl %>% pull(degreehours_over_30C_priorweek))
sd_temp_deghr30week_webl <- sd(data_deghr30week_webl %>% pull(degreehours_over_30C_priorweek))


temp_trans_deghr30week_webl <- trans_new("temp_trans_deghr30week_webl",
                                      transform = function(x){(x * sd_temp_deghr30week_webl) + mean_temp_deghr30week_webl},
                                      inverse = function(x){x})

(fig2_deghr30week_webl <- predict_response(g_lintemp,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week day") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dashed")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_deghr30week_webl,
                       breaks = c((0-mean_temp_deghr30week_webl)/sd_temp_deghr30week_webl,
                                  (1000-mean_temp_deghr30week_webl)/sd_temp_deghr30week_webl,
                                  (2000-mean_temp_deghr30week_webl)/sd_temp_deghr30week_webl,
                                  (3000-mean_temp_deghr30week_webl)/sd_temp_deghr30week_webl,
                                  (4000-mean_temp_deghr30week_webl)/sd_temp_deghr30week_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$degreehours_over_30C_priorweek,na.rm = TRUE))/sd(g$degreehours_over_30C_priorweek,na.rm = TRUE),
                       #            (55-mean(g$degreehours_over_30C_priorweek,na.rm = TRUE))/sd(g$degreehours_over_30C_priorweek,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_deghr30week_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbydeghr30weekxhab_WEBL.png",plot =  fig2_deghr30week_webl, width = 10, height = 6.6)



summary(g_lintemp)



(weblgrowthtrenddeghr30week <- emtrends(g_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblgrowthtrenddeghr30week,"figures/weblgrowthtrenddeghr30week.html")

data_deghr30week_webl = dplyr::filter(g,Species == "WEBL",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(degreehours_over_30C_priorweek_scaled = (degreehours_over_30C_priorweek_scaled * sd(data$degreehours_over_30C_priorweek)) + mean(data$degreehours_over_30C_priorweek),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     degreehours_over_30C_priorweek_scaled = round(degreehours_over_30C_priorweek_scaled),
#     degreehours_over_30C_priorweek_scaled = paste0(degreehours_over_30C_priorweek_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "degreehours_over_30C_priorweek_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/weblgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ degreehours_over_30C_priorweek_scaled, at = list(degreehours_over_30C_priorweek_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_deghr30week_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_deghr30week_webl,"figures/growthbyhabitat_deghr30week_webl.html")


### Check for effect of temperature


(growthbyhabitat_summary_deghr30week_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_deghr30week_webl,"figures/growthbyhabitat_summary_deghr30week_webl.html")

#### hihours_over_30hi_priorweek


g_lintemp <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(`hihours_over_30hi_priorweek`),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_hihr25week_webl <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_hihr25week_webl,"figures/int_tab_growth_hihr25week_webl.html")
summary(g_lintemp)
check_collinearity(g_lintemp_noint)

g_lintemp_hihr25week_webl <- g_lintemp

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_hihr25week_webl <- g_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_hihr25week_webl <- samp_year_hihr25week_webl %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_hihr25week_webl

samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_hihr25week_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_hihr25week_webl = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_temp_hihr25week_webl <- mean(data_hihr25week_webl %>% pull(hihours_over_30hi_priorweek))
sd_temp_hihr25week_webl <- sd(data_hihr25week_webl %>% pull(hihours_over_30hi_priorweek))


temp_trans_hihr25week_webl <- trans_new("temp_trans_hihr25week_webl",
                                         transform = function(x){(x * sd_temp_hihr25week_webl) + mean_temp_hihr25week_webl},
                                         inverse = function(x){x})

(fig2_hihr25week_webl <- predict_response(g_lintemp,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index-hours >25\u00b0C over prior day") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dashed")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_hihr25week_webl,
                       breaks = c((0-mean_temp_hihr25week_webl)/sd_temp_hihr25week_webl,
                                  (1000-mean_temp_hihr25week_webl)/sd_temp_hihr25week_webl,
                                  (2000-mean_temp_hihr25week_webl)/sd_temp_hihr25week_webl,
                                  (3000-mean_temp_hihr25week_webl)/sd_temp_hihr25week_webl,
                                  (4000-mean_temp_hihr25week_webl)/sd_temp_hihr25week_webl),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$hihours_over_30hi_priorweek,na.rm = TRUE))/sd(g$hihours_over_30hi_priorweek,na.rm = TRUE),
                       #            (55-mean(g$hihours_over_30hi_priorweek,na.rm = TRUE))/sd(g$hihours_over_30hi_priorweek,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_hihr25week_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbyhihr25weekxhab_WEBL.png",plot =  fig2_hihr25week_webl, width = 10, height = 6.6)



summary(g_lintemp)



(weblgrowthtrendhihr25week <- emtrends(g_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(weblgrowthtrendhihr25week,"figures/weblgrowthtrendhihr25week.html")

data_hihr25week_webl = dplyr::filter(g,Species == "WEBL",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(hihours_over_30hi_priorweek_scaled = (hihours_over_30hi_priorweek_scaled * sd(data$hihours_over_30hi_priorweek)) + mean(data$hihours_over_30hi_priorweek),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     hihours_over_30hi_priorweek_scaled = round(hihours_over_30hi_priorweek_scaled),
#     hihours_over_30hi_priorweek_scaled = paste0(hihours_over_30hi_priorweek_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "hihours_over_30hi_priorweek_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/weblgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ hihours_over_30hi_priorweek_scaled, at = list(hihours_over_30hi_priorweek_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_hihr25week_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_hihr25week_webl,"figures/growthbyhabitat_hihr25week_webl.html")


### Check for effect of temperature


(growthbyhabitat_summary_hihr25week_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_hihr25week_webl,"figures/growthbyhabitat_summary_hihr25week_webl.html")



### TRES
#### maxhi_week

g_lintemp <- lmerTest::lmer(gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_maxhiweek_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_maxhiweek_tres,"figures/int_tab_growth_maxhiweek_tres.html")
summary(g_lintemp)
check_collinearity(g_lintemp_noint)

g_lintemp_maxhiweek_tres <- g_lintemp_addmin

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_maxhiweek_tres <- g_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_maxhiweek_tres <- samp_year_maxhiweek_tres %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_maxhiweek_tres

samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_maxhiweek_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_maxhiweek_tres = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)


mean_temp_maxhiweek_tres <- mean(data_maxhiweek_tres %>% pull(maxhi_week))
sd_temp_maxhiweek_tres <- sd(data_maxhiweek_tres %>% pull(maxhi_week))


temp_trans_maxhiweek_tres <- trans_new("temp_trans_maxhiweek_tres",
                                       transform = function(x){(x * sd_temp_maxhiweek_tres) + mean_temp_maxhiweek_tres},
                                       inverse = function(x){x})

(fig2_maxhiweek_tres <- predict_response(g_lintemp_addmin,terms = c("maxhi_week_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max heat index (\u00b0C) over preceding week") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_maxhiweek_tres,
                       breaks = c((20-mean_temp_maxhiweek_tres)/sd_temp_maxhiweek_tres,
                                  (25-mean_temp_maxhiweek_tres)/sd_temp_maxhiweek_tres,
                                  (30-mean_temp_maxhiweek_tres)/sd_temp_maxhiweek_tres,
                                  (35-mean_temp_maxhiweek_tres)/sd_temp_maxhiweek_tres,
                                  (40-mean_temp_maxhiweek_tres)/sd_temp_maxhiweek_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$maxhi_week,na.rm = TRUE))/sd(g$maxhi_week,na.rm = TRUE),
                       #            (55-mean(g$maxhi_week,na.rm = TRUE))/sd(g$maxhi_week,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_maxhiweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbymaxhiweekxhab_TRES.png",plot =  fig2_maxhiweek_tres, width = 10, height = 6.6)



summary(g_lintemp_addmin)



(tresgrowthtrendmaxhiweek <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresgrowthtrendmaxhiweek,"figures/tresgrowthtrendmaxhiweek.html")

data_maxhiweek_tres = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(maxhi_week_scaled = (maxhi_week_scaled * sd(data$maxhi_week)) + mean(data$maxhi_week),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     maxhi_week_scaled = round(maxhi_week_scaled),
#     maxhi_week_scaled = paste0(maxhi_week_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "maxhi_week_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/tresgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("maxhi_week_scaled"), at = list(maxhi_week_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ maxhi_week_scaled, at = list(maxhi_week_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiweek_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_maxhiweek_tres,"figures/growthbyhabitat_maxhiweek_tres.html")


### Check for effect of temperature


(growthbyhabitat_summary_maxhiweek_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_maxhiweek_tres,"figures/growthbyhabitat_summary_maxhiweek_tres.html")


#### maxhi_prior


g_lintemp <- lmerTest::lmer(gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ maxhi_prior_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ maxhi_prior_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_maxhiday_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_maxhiday_tres,"figures/int_tab_growth_maxhiday_tres.html")
summary(g_lintemp_addmin)
check_collinearity(g_lintemp_noint)

g_lintemp_maxhiday_tres <- g_lintemp_addmin

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_maxhiday_tres <- g_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_maxhiday_tres <- samp_year_maxhiday_tres %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_maxhiday_tres

samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_maxhiday_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_maxhiday_tres = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)


mean_temp_maxhiday_tres <- mean(data_maxhiday_tres %>% pull(maxhi_prior))
sd_temp_maxhiday_tres <- sd(data_maxhiday_tres %>% pull(maxhi_prior))


temp_trans_maxhiday_tres <- trans_new("temp_trans_maxhiday_tres",
                                      transform = function(x){(x * sd_temp_maxhiday_tres) + mean_temp_maxhiday_tres},
                                      inverse = function(x){x})

(fig2_maxhiday_tres <- predict_response(g_lintemp_addmin,terms = c("maxhi_prior_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max heat index (\u00b0C) over preceding day") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_maxhiday_tres,
                       breaks = c((20-mean_temp_maxhiday_tres)/sd_temp_maxhiday_tres,
                                  (25-mean_temp_maxhiday_tres)/sd_temp_maxhiday_tres,
                                  (30-mean_temp_maxhiday_tres)/sd_temp_maxhiday_tres,
                                  (35-mean_temp_maxhiday_tres)/sd_temp_maxhiday_tres,
                                  (40-mean_temp_maxhiday_tres)/sd_temp_maxhiday_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$maxhi_prior,na.rm = TRUE))/sd(g$maxhi_prior,na.rm = TRUE),
                       #            (55-mean(g$maxhi_prior,na.rm = TRUE))/sd(g$maxhi_prior,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_maxhiday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbymaxhidayxhab_TRES.png",plot =  fig2_maxhiday_tres, width = 10, height = 6.6)



summary(g_lintemp_addmin)



(tresgrowthtrendmaxhiday <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresgrowthtrendmaxhiday,"figures/tresgrowthtrendmaxhiday.html")

data_maxhiday_tres = dplyr::filter(g,Species == "TRES",!is.na(maxhi_prior),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,maxhi_prior,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         maxhi_prior_scaled_sq = maxhi_prior_scaled * maxhi_prior_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(maxhi_prior_scaled = (maxhi_prior_scaled * sd(data$maxhi_prior)) + mean(data$maxhi_prior),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     maxhi_prior_scaled = round(maxhi_prior_scaled),
#     maxhi_prior_scaled = paste0(maxhi_prior_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "maxhi_prior_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/tresgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("maxhi_prior_scaled"), at = list(maxhi_prior_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ maxhi_prior_scaled, at = list(maxhi_prior_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiday_tres <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_maxhiday_tres,"figures/growthbyhabitat_maxhiday_tres.html")


### Check for effect of temperature


(growthbyhabitat_summary_maxhiday_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_maxhiday_tres,"figures/growthbyhabitat_summary_maxhiday_tres.html")


#### degreehours_over_30C_priorweek


g_lintemp <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(`degreehours_over_30C_priorweek`),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_deghr30week_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_deghr30week_tres,"figures/int_tab_growth_deghr30week_tres.html")
summary(g_lintemp)
check_collinearity(g_lintemp_noint)

g_lintemp_deghr30week_tres <- g_lintemp_addmin

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_deghr30week_tres <- g_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_deghr30week_tres <- samp_year_deghr30week_tres %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_deghr30week_tres

samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_deghr30week_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_deghr30week_tres = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)


mean_temp_deghr30week_tres <- mean(data_deghr30week_tres %>% pull(degreehours_over_30C_priorweek))
sd_temp_deghr30week_tres <- sd(data_deghr30week_tres %>% pull(degreehours_over_30C_priorweek))


temp_trans_deghr30week_tres <- trans_new("temp_trans_deghr30week_tres",
                                         transform = function(x){(x * sd_temp_deghr30week_tres) + mean_temp_deghr30week_tres},
                                         inverse = function(x){x})

(fig2_deghr30week_tres <- predict_response(g_lintemp_addmin,terms = c("degreehours_over_30C_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C over prior week day") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dotted","Grassland" = "solid","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_deghr30week_tres,
                       breaks = c((0-mean_temp_deghr30week_tres)/sd_temp_deghr30week_tres,
                                  (1000-mean_temp_deghr30week_tres)/sd_temp_deghr30week_tres,
                                  (2000-mean_temp_deghr30week_tres)/sd_temp_deghr30week_tres,
                                  (3000-mean_temp_deghr30week_tres)/sd_temp_deghr30week_tres,
                                  (4000-mean_temp_deghr30week_tres)/sd_temp_deghr30week_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$degreehours_over_30C_priorweek,na.rm = TRUE))/sd(g$degreehours_over_30C_priorweek,na.rm = TRUE),
                       #            (55-mean(g$degreehours_over_30C_priorweek,na.rm = TRUE))/sd(g$degreehours_over_30C_priorweek,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_deghr30week_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbydeghr30weekxhab_TRES.png",plot =  fig2_deghr30week_tres, width = 10, height = 6.6)



summary(g_lintemp_addmin)



(tresgrowthtrenddeghr30week <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresgrowthtrenddeghr30week,"figures/tresgrowthtrenddeghr30week.html")

data_deghr30week_tres = dplyr::filter(g,Species == "TRES",!is.na(degreehours_over_30C_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,degreehours_over_30C_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         degreehours_over_30C_priorweek_scaled_sq = degreehours_over_30C_priorweek_scaled * degreehours_over_30C_priorweek_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(degreehours_over_30C_priorweek_scaled = (degreehours_over_30C_priorweek_scaled * sd(data$degreehours_over_30C_priorweek)) + mean(data$degreehours_over_30C_priorweek),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     degreehours_over_30C_priorweek_scaled = round(degreehours_over_30C_priorweek_scaled),
#     degreehours_over_30C_priorweek_scaled = paste0(degreehours_over_30C_priorweek_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "degreehours_over_30C_priorweek_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/tresgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("degreehours_over_30C_priorweek_scaled"), at = list(degreehours_over_30C_priorweek_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ degreehours_over_30C_priorweek_scaled, at = list(degreehours_over_30C_priorweek_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_deghr30week_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_deghr30week_tres,"figures/growthbyhabitat_deghr30week_tres.html")


### Check for effect of temperature


(growthbyhabitat_summary_deghr30week_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_deghr30week_tres,"figures/growthbyhabitat_summary_deghr30week_tres.html")

#### hihours_over_30hi_priorweek


g_lintemp <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(`hihours_over_30hi_priorweek`),!is.na(meanmintempI)) %>%
                              mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled"),
                                     hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
                                    mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled))

c1 <- anova(g_lintemp,g_lintemp_addmin,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

c2 <- anova(g_lintemp,g_lintemp_addmax,g_lintemp_noint) %>% tibble() %>% mutate(Model = c("no interaction","single interaction","both interacting"),.before = npar)

(int_tab_growth_hihr25week_tres <- bind_rows(c1,c2) %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
# gtsave(int_tab_growth_hihr25week_tres,"figures/int_tab_growth_hihr25week_tres.html")
summary(g_lintemp_addmin)
check_collinearity(g_lintemp_noint)

g_lintemp_hihr25week_tres <- g_lintemp_addmin

### Conclusion: growth interacts with max and min.

### Sample sizes:


samp_year_hihr25week_tres <- g_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
(t_samp_year_growth_hihr25week_tres <- samp_year_hihr25week_tres %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))


t_samp_year_growth_hihr25week_tres

samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_hihr25week_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)





data_hihr25week_tres = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)


mean_temp_hihr25week_tres <- mean(data_hihr25week_tres %>% pull(hihours_over_30hi_priorweek))
sd_temp_hihr25week_tres <- sd(data_hihr25week_tres %>% pull(hihours_over_30hi_priorweek))


temp_trans_hihr25week_tres <- trans_new("temp_trans_hihr25week_tres",
                                        transform = function(x){(x * sd_temp_hihr25week_tres) + mean_temp_hihr25week_tres},
                                        inverse = function(x){x})

(fig2_hihr25week_tres <- predict_response(g_lintemp_addmin,terms = c("hihours_over_30hi_priorweek_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index-hours >25\u00b0C over prior day") +
    ylab("Growth (g/day)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "solid")) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_hihr25week_tres,
                       breaks = c((0-mean_temp_hihr25week_tres)/sd_temp_hihr25week_tres,
                                  (1000-mean_temp_hihr25week_tres)/sd_temp_hihr25week_tres,
                                  (2000-mean_temp_hihr25week_tres)/sd_temp_hihr25week_tres,
                                  (3000-mean_temp_hihr25week_tres)/sd_temp_hihr25week_tres,
                                  (4000-mean_temp_hihr25week_tres)/sd_temp_hihr25week_tres),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$hihours_over_30hi_priorweek,na.rm = TRUE))/sd(g$hihours_over_30hi_priorweek,na.rm = TRUE),
                       #            (55-mean(g$hihours_over_30hi_priorweek,na.rm = TRUE))/sd(g$hihours_over_30hi_priorweek,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    #ylim(-5,5) +
    geom_text(data = dat_text_hihr25week_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

# ggsave("figures/growthbyhihr25weekxhab_TRES.png",plot =  fig2_hihr25week_tres, width = 10, height = 6.6)



summary(g_lintemp_addmin)



(tresgrowthtrendhihr25week <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


# gtsave(tresgrowthtrendhihr25week,"figures/tresgrowthtrendhihr25week.html")

data_hihr25week_tres = dplyr::filter(g,Species == "TRES",!is.na(hihours_over_30hi_priorweek),!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,hihours_over_30hi_priorweek,meanmintempI,juliandate,brood_size,age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         hihours_over_30hi_priorweek_scaled_sq = hihours_over_30hi_priorweek_scaled * hihours_over_30hi_priorweek_scaled)

# (t <- emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(hihours_over_30hi_priorweek_scaled = (hihours_over_30hi_priorweek_scaled * sd(data$hihours_over_30hi_priorweek)) + mean(data$hihours_over_30hi_priorweek),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     hihours_over_30hi_priorweek_scaled = round(hihours_over_30hi_priorweek_scaled),
#     hihours_over_30hi_priorweek_scaled = paste0(hihours_over_30hi_priorweek_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "hihours_over_30hi_priorweek_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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
# gtsave(t,"figures/tresgrowthdeltamax.html")

# ((emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("hihours_over_30hi_priorweek_scaled"), at = list(hihours_over_30hi_priorweek_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp,formula = habitat ~ hihours_over_30hi_priorweek_scaled, at = list(hihours_over_30hi_priorweek_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_hihr25week_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_hihr25week_tres,"figures/growthbyhabitat_hihr25week_tres.html")


### Check for effect of temperature


(growthbyhabitat_summary_hihr25week_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_summary_hihr25week_tres,"figures/growthbyhabitat_summary_hihr25week_tres.html")


## Combined WEBL and TRES growth plots


# ggplot_build(fig2_webl)$layout$panel_scales_y
# (p_full <- ggarrange(fig2_webl + theme(text = element_text(size = 12),axis.title.x = element_blank()),fig2_tres + ylim(-1.09,3.42) + theme(axis.text.y = element_blank(),
#                                                                                                                                            axis.ticks.y = element_blank(),
#                                                                                                                                            text = element_text(size = 12),
#                                                                                                                                            axis.title.y = element_blank(),
#                                                                                                                                            axis.title.x = element_text(hjust = 2.8)),ncol = 2,
#                      labels = c("(a): Western Bluebird","(b): Tree Swallow")))
#
# ggsave("figures/fig2_growth_by_temp_hab.png",p_full,width = 6.25,height = 4)




# Models to estimate the unbiased causal effect of cort and provisioning on growth

#Model structure pulled from DAG: will use a minimal adjustment set to estimate the total effect of cort and provisioning on growth.
#{ attempt, humidity, nest_age, nestcond, surftemp }

#I could subtract the direct from the total effect to get the indirect effect. I'm still vacillating between thinking that I should be measuring the total effect or the direct effect. I think because nestcond is really the only mediator in the DAG that I actually want total effect?

#  The question I'd like to ask: Does the direct effect of these two factors on growth change depending on habitat? I'm not sure that this causal structure is the correct one for that question. Is that a mediation analysis? After reading some papers, I think that actually the question is, do cort and provisioning mediate the effect of habitat and temperature on growth? I still think there may be other questions in this dataset. Maybe the answer is to do this paper just linear modeling with causal inference support, and then get Sage in on a paper that delves more deeply into a structural causal modeling framework.

#The following focuses on the more simple question: how do cort and provisioning affect growth (total effect)? I will do the same model structure for survival in that notebook.

#adjustment set: { attempt, humidity, nest_age, nestcond, surftemp }



# g_provis_cort <- lmerTest::lmer(gweight ~ cort_s1_scaled + provis_mean_scaled + poly(maxhi_week_scaled,2) + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45) %>%
#                                   mutate(across(c(gweight,maxhi_week,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                                                 ~ scale(.x)[,1],
#                                                 .names = "{.col}_scaled"),
#                                          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
#
# g_provis_cort_maxhiweek_webl <- g_provis_cort
# check_collinearity(g_provis_cort)
# summary(g_provis_cort)
#
# (growth_by_cort_provis_summary_maxhiweek_webl <- summary(g_provis_cort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#     # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#     mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
#            across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_cort_provis_summary_maxhiweek_webl,"figures/growth_by_cort_provis_summary_maxhiweek_webl.html")
#
# data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45) %>%
#   mutate(across(c(gweight,maxhi_week,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
#
# # (t <- emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
# #     mutate(provis_mean_scaled = (provis_mean_scaled * sd(data$provis_mean, na.rm = TRUE)) + mean(data$provis_mean,na.rm = TRUE),
# #            across(where(is.numeric), ~ round(.x, digits = 2)),
# #            provis_mean_scaled = round(provis_mean_scaled)) %>%
# #     #tibble() %>%
# #     dplyr::select(-df) %>%
# #     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
# #     rename(`Provisioning` = "provis_mean_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
# #     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
# #     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
# #     gt() %>% tab_options(data_row.padding = px(1)))
# # gtsave(t,"figures/growth_by_basecort_provis_delta_webl.html")
# #
# # ((emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
#
#
#
#
# g_provis_abscort <- lmerTest::lmer(gweight ~ abs_change_cort_scaled + provis_mean_scaled + maxhi_week_scaled + meanmintempI_scaled + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45,!is.na(meanmintempI)) %>%
#                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
#                                                    ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
# g_provis_abscort_maxhiweek_webl <- g_provis_abscort
#
# summary(g_provis_abscort)
#
# (growth_by_abscort_provis_summary_maxhiweek_webl <- summary(g_provis_abscort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#     # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#     mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
#            across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_abscort_provis_summary_maxhiweek_webl,"figures/growth_by_abscort_provis_summary_maxhiweek_webl.html")
#
# data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45,!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
#
# (t <- emmeans(g_provis_abscort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(provis_mean_scaled = (provis_mean_scaled * sd(data$provis_mean,na.rm = TRUE)) + mean(data$provis_mean,na.rm = TRUE),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            provis_mean_scaled = round(provis_mean_scaled)) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(`max temperature` = "provis_mean_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)))
# gtsave(t,"figures/growth_by_abscort_provis_delta_maxhiweek_webl.html")
#
#
#
# #### Sample size
#
#
# data_s1_maxhiweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45) %>%
#   mutate(across(c(gweight,maxhi_week,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
# ss_year_webl_growthbybasecortprovis_maxhiweek <- data_s1_maxhiweek_webl %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_webl_growthbybasecortprovis_maxhiweek %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_webl_growthbybasecortprovis_maxhiweek.html")
# data_abs_maxhiweek_webl = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45,!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
# ss_year_webl_growthbyabscortprovis_maxhiweek <- data_abs_maxhiweek_webl %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_webl_growthbyabscortprovis_maxhiweek %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_webl_growthbyabscortprovis_maxhiweek.html")
#
#
#
#
#
# mean_provis_maxhiweek_webl <- mean(data_s1_maxhiweek_webl %>% pull(provis_mean))
# sd_provis_maxhiweek_webl <- sd(data_s1_maxhiweek_webl %>% pull(provis_mean))
#
#
# provis_trans_maxhiweek_webl <- trans_new("provis_trans_maxhiweek_webl",
#                                transform = function(x){(x * sd_provis_maxhiweek_webl) + mean_provis_maxhiweek_webl},
#                                inverse = function(x){x})
#
# mean_s1_maxhiweek_webl <- mean(data_s1_maxhiweek_webl %>% pull(cort_s1))
# sd_s1_maxhiweek_webl <- sd(data_s1_maxhiweek_webl %>% pull(cort_s1))
#
#
# s1_trans_maxhiweek_webl <- trans_new("s1_trans_maxhiweek_webl",
#                            transform = function(x){(x * sd_s1_maxhiweek_webl) + mean_s1_maxhiweek_webl},
#                            inverse = function(x){x})
#
# mean_abs_maxhiweek_webl <- mean(data_abs_maxhiweek_webl %>% pull(abs_change_cort))
# sd_abs_maxhiweek_webl <- sd(data_abs_maxhiweek_webl %>% pull(abs_change_cort))
#
#
# abs_trans_maxhiweek_webl <- trans_new("abs_trans_maxhiweek_webl",
#                             transform = function(x){(x * sd_abs_maxhiweek_webl) + mean_abs_maxhiweek_webl},
#                             inverse = function(x){x})
#
# # (fig6_provis_maxhiweek_webl <- predict_response(g_provis_cort,terms = c("provis_mean_scaled"),bias_correction = TRUE,margin = "empirical") %>%
# #     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
# #     theme_classic() +
# #     #facet_wrap(~ group, ncol = 2) +
# #     xlab("Provisioning") +
# #     ylab("Growth (g/day)") +
# #     #scale_fill_viridis(discrete = TRUE) +
# #     #scale_color_viridis(discrete = TRUE) +
# #     theme(text = element_text(size = 16)) +
# #     labs(title = element_blank())  +
# #     scale_x_continuous(trans = provis_trans_maxhiweek_webl,
# #                        breaks = c(#(0-mean_provis_maxhiweek_webl)/sd_provis_maxhiweek_webl,
# #                                   (5-mean_provis_maxhiweek_webl)/sd_provis_maxhiweek_webl,
# #                                   (10-mean_provis_maxhiweek_webl)/sd_provis_maxhiweek_webl,
# #                                   (15-mean_provis_maxhiweek_webl)/sd_provis_maxhiweek_webl,
# #                                   (20-mean_provis_maxhiweek_webl)/sd_provis_maxhiweek_webl)
# #     ) +
# #     #ylim(-5,5) +
# #     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
# #     theme(legend.position = "none") +
# #     annotate(geom = "text",label = "N = 40",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5)
# # )
# #
# # (fig6_corts1_maxhiweek_webl <- predict_response(g_provis_cort,terms = c("cort_s1_scaled [all]"),bias_correction = TRUE,margin = "empirical") %>%
# #     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
# #     aes(linetype = .data[["group"]]) +
# #     theme_classic() +
# #     #facet_wrap(~ group, ncol = 2) +
# #     xlab("Baseline corticosterone (ng/\U00B5L)") +
# #     ylab("Growth (g/day)") +
# #     #scale_fill_viridis(discrete = TRUE) +
# #     #scale_color_viridis(discrete = TRUE) +
# #     theme(text = element_text(size = 16)) +
# #     labs(title = element_blank())  +
# #     scale_x_continuous(trans = s1_trans_webl,
# #                        breaks = c((2-mean_s1_webl)/sd_s1_webl,
# #                                   (4-mean_s1_webl)/sd_s1_webl,
# #                                   (6-mean_s1_webl)/sd_s1_webl,
# #                                   (8-mean_s1_webl)/sd_s1_webl,
# #                                   (10-mean_s1_webl)/sd_s1_webl,
# #                                   (12-mean_s1_webl)/sd_s1_webl,
# #                                   (14-mean_s1_webl)/sd_s1_webl)
# #     ) +
# #     scale_linetype_manual(values = "dotted") +
# #     #ylim(-5,5) +
# #     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
# #     theme(legend.position = "none") +
# #     annotate(geom = "text",label = "N = 40",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5)
# # )
# #
# # (fig6_abscort_webl <- predict_response(g_provis_abscort,terms = c("abs_change_cort_scaled"),bias_correction = TRUE,margin = "empirical") %>%
# #     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
# #     aes(linetype = .data[["group"]]) +
# #     theme_classic() +
# #     #facet_wrap(~ group, ncol = 2) +
# #     xlab("Stress-induced corticosterone (ng/\U00B5L)") +
# #     ylab("Growth (g/day)") +
# #     #scale_fill_viridis(discrete = TRUE) +
# #     #scale_color_viridis(discrete = TRUE) +
# #     theme(text = element_text(size = 16)) +
# #     labs(title = element_blank())  +
# #     scale_x_continuous(trans = abs_trans_webl,
# #                        breaks = c((0-mean_abs_webl)/sd_abs_webl,
# #                                   (10-mean_abs_webl)/sd_abs_webl,
# #                                   (20-mean_abs_webl)/sd_abs_webl,
# #                                   (30-mean_abs_webl)/sd_abs_webl,
# #                                   (40-mean_abs_webl)/sd_abs_webl,
# #                                   (50-mean_abs_webl)/sd_abs_webl)
# #     ) +
# #     scale_linetype_manual(values = "dotted") +
# #     #ylim(-5,5) +
# #     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
# #     theme(legend.position = "none") +
# #     annotate(geom = "text",label = "N = 35",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5)
# # )
# #
# # ggsave("figures/growthbyprovis_WEBL.png",plot =  fig6_provis_webl, width = 10, height = 6.6)
#
#
#
# ## Canonical adjustment set overloads the data:
#
# #{ attempt, broodsize, habitat, humidity, jul_date, mother_s1, nest_age, site, surftemp, year }
#
#
#
# # g_provis_cort_canon <- lmerTest::lmer(gweight ~ cort_s1_scaled + provis_mean_scaled + (1|attempt_id) + brood_size_scaled + habitat + meanh_scaled + juliandate_scaled + cort_s1_mother_scaled + age_scaled + maxhi_week_scaled + meanmintempI_scaled + year_fct,data = dplyr::filter(g,Species == "WEBL",!is.na(maxhi_week),maxhi_week < 45,!is.na(meanmintempI)) %>%
# #                                         mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,cort_s1_mother,condition,provis_mean),
# #                                                       ~ scale(.x)[,1],
# #                                                       .names = "{.col}_scaled"),
# #                                                maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
# # summary(g_provis_cort_canon)
# #
# #
#
# ## TRES
#
#
#
# g_provis_cort <- lmerTest::lmer(gweight ~ cort_s1_scaled + provis_mean_scaled + poly(maxhi_week_scaled,2) + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#                                   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                                                 ~ scale(.x)[,1],
#                                                 .names = "{.col}_scaled"),
#                                          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
# check_collinearity(g_provis_cort)
# summary(g_provis_cort)
#
# g_provis_cort_tres <- g_provis_cort
#
# (growth_by_cort_provis_summary_tres <- summary(g_provis_cort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#     # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#     mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
#            across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_cort_provis_summary_tres,"figures/growth_by_cort_provis_summary_tres.html")
#
# data = dplyr::filter(g,Species == "TRES",!is.na(gweight),!is.na(cort_s1),!is.na(provis_mean),!is.na(maxhi_week),
#                      !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
#
# # (t <- emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
# #     mutate(cort_s1_scaled = (cort_s1_scaled * sd(data$cort_s1, na.rm = TRUE)) + mean(data$cort_s1,na.rm = TRUE),
# #            across(where(is.numeric), ~ round(.x, digits = 2)),
# #            cort_s1_scaled = round(cort_s1_scaled)) %>%
# #     #tibble() %>%
# #     dplyr::select(-df) %>%
# #     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
# #     rename(`Baseline cort` = "cort_s1_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
# #     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
# #     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
# #     gt() %>% tab_options(data_row.padding = px(1)))
# # gtsave(t,"figures/growth_by_basecort_provis_delta_tres.html")
# #
# # ((emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
#
#
# #### Abs diff cort
#
#
#
# g_provis_abscort <- lmerTest::lmer(gweight ~ abs_change_cort_scaled + provis_mean_scaled + maxhi_week_scaled + meanmintempI_scaled + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(g,Species == "TRES",!is.na(maxhi_week),!is.na(meanmintempI)) %>%
#                                      mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
#                                                    ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                             maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled))
# summary(g_provis_cort)
# g_provis_abscort_tres <- g_provis_abscort
#
# (growth_by_abscort_provis_summary_tres <- summary(g_provis_abscort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#     # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#     mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
#            across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_abscort_provis_summary_tres,"figures/growth_by_abscort_provis_summary_tres.html")
#
# data = dplyr::filter(g,Species == "TRES",!is.na(gweight),!is.na(abs_change_cort),!is.na(provis_mean),!is.na(maxhi_week),
#                      !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
#
# (t <- emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(abs_change_cort_scaled = (abs_change_cort_scaled * sd(data$abs_change_cort, na.rm = TRUE)) + mean(data$abs_change_cort,na.rm = TRUE),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            abs_change_cort_scaled = round(abs_change_cort_scaled)) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(`Baseline cort` = "abs_change_cort_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)))
# gtsave(t,"figures/growth_by_abscort_provis_delta_tres.html")
#
# ((emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
#
#
# #### Sample size
#
#
# data_s1_tres = dplyr::filter(g,Species == "TRES",!is.na(gweight),!is.na(cort_s1),!is.na(provis_mean),!is.na(maxhi_week),
#                              !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
# ss_year_TRES_growthbybasecortprovis <- data_s1_tres %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_TRES_growthbybasecortprovis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_TRES_growthbybasecortprovis.html")
# data_abs_tres = dplyr::filter(g,Species == "TRES",!is.na(gweight),!is.na(abs_change_cort),!is.na(provis_mean),!is.na(maxhi_week),
#                               !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
#   mutate(across(c(gweight,maxhi_week,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
#                 ~ scale(.x)[,1],
#                 .names = "{.col}_scaled"),
#          maxhi_week_scaled_sq = maxhi_week_scaled * maxhi_week_scaled)
# ss_year_TRES_growthbyabscortprovis <- data_abs_tres %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_TRES_growthbyabscortprovis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_TRES_growthbyabscortprovis.html")
#
#
#
#
#
# mean_provis_tres <- mean(data_s1_tres %>% pull(provis_mean))
# sd_provis_tres <- sd(data_s1_tres %>% pull(provis_mean))
#
#
# provis_trans_tres <- trans_new("provis_trans_tres",
#                                transform = function(x){(x * sd_provis_tres) + mean_provis_tres},
#                                inverse = function(x){x})
#
# mean_s1_tres <- mean(data_s1_tres %>% pull(cort_s1))
# sd_s1_tres <- sd(data_s1_tres %>% pull(cort_s1))
#
#
# s1_trans_tres <- trans_new("s1_trans_tres",
#                            transform = function(x){(x * sd_s1_tres) + mean_s1_tres},
#                            inverse = function(x){x})
#
# mean_abs_tres <- mean(data_abs_tres %>% pull(abs_change_cort))
# sd_abs_tres <- sd(data_abs_tres %>% pull(abs_change_cort))
#
#
# abs_trans_tres <- trans_new("abs_trans_tres",
#                             transform = function(x){(x * sd_abs_tres) + mean_abs_tres},
#                             inverse = function(x){x})
#
# # (fig6_tres <- predict_response(g_provis_cort,terms = c("provis_mean_scaled","cort_s1_scaled [-2,2]"),bias_correction = TRUE) %>%
# #     #filter(group != 0.35) %>%
# #     # mutate(group = as_factor(group) %>% fct_collapse(Low = "-0.93",High = "1.63")) %>%
# #     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE,colors = "viridis") +
# #     #scale_fill_discrete(labels = c("Low","High")) +
# #     #scale_color_discrete(labels = c("Low","High")) +
# #     # geom_point(show.legend = FALSE) +
# #     theme_classic() +
# #     #facet_wrap(~ group, ncol = 2) +
# #     xlab("Provisioning (mean visits/hr)") +
# #     ylab("Growth (g/day)") +
# #     #scale_fill_viridis(discrete = TRUE) +
# #     #scale_color_viridis(discrete = TRUE) +
# #     theme(text = element_text(size = 16)) +
# #     labs(title = element_blank(),color = "Baseline cort") +
# #     scale_x_continuous(trans = provis_trans_tres,
# #                        breaks = c((0-mean_provis_tres)/sd_provis_tres,
# #                                   (10-mean_provis_tres)/sd_provis_tres,
# #                                   (20-mean_provis_tres)/sd_provis_tres,
# #                                   (30-mean_provis_tres)/sd_provis_tres,
# #                                   (40-mean_provis_tres)/sd_provis_tres,
# #                                   (50-mean_provis_tres)/sd_provis_tres)
# #     ) +
# #     #ylim(-5,5) #+
# #     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
# #     theme(legend.position = "none") +
# #     annotate("text",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5,label = "N = 16")
# # )
#
# (fig6_provis_tres <- predict_response(g_provis_cort,terms = c("provis_mean_scaled"),bias_correction = TRUE,margin = "empirical") %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     aes(linetype = .data[["group"]]) +
#     theme_classic() +
#     #facet_wrap(~ group, ncol = 2) +
#     xlab("Provisioning (visits/hr)") +
#     ylab("Growth (g/day)") +
#     #scale_fill_viridis(discrete = TRUE) +
#     #scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = element_blank())  +
#     scale_x_continuous(trans = provis_trans_tres,
#                        breaks = c((0-mean_provis_tres)/sd_provis_tres,
#                                   (10-mean_provis_tres)/sd_provis_tres,
#                                   (20-mean_provis_tres)/sd_provis_tres,
#                                   (30-mean_provis_tres)/sd_provis_tres,
#                                   (40-mean_provis_tres)/sd_provis_tres,
#                                   (50-mean_provis_tres)/sd_provis_tres)
#     ) +
#     scale_linetype_manual(values = "dotted") +
#     #ylim(-5,5) +
#     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
#     theme(legend.position = "none") +
#     annotate(geom = "text",label = "N = 16",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5)
# )
#
# (fig6_corts1_tres <- predict_response(g_provis_cort,terms = c("cort_s1_scaled [all]"),bias_correction = TRUE,margin = "empirical") %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     theme_classic() +
#     #facet_wrap(~ group, ncol = 2) +
#     xlab("Baseline corticosterone (ng/\U00B5L)") +
#     ylab("Growth (g/day)") +
#     #scale_fill_viridis(discrete = TRUE) +
#     #scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = element_blank())  +
#     scale_x_continuous(trans = s1_trans_tres,
#                        breaks = c((5-mean_s1_tres)/sd_s1_tres,
#                                   (10-mean_s1_tres)/sd_s1_tres,
#                                   (15-mean_s1_tres)/sd_s1_tres,
#                                   (20-mean_s1_tres)/sd_s1_tres,
#                                   (25-mean_s1_tres)/sd_s1_tres)
#     ) +
#     #ylim(-5,5) +
#     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
#     theme(legend.position = "none") +
#     annotate(geom = "text",label = "N = 16",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5)
# )
#
# (fig6_abscort_tres <- predict_response(g_provis_abscort,terms = c("abs_change_cort_scaled"),bias_correction = TRUE,margin = "empirical") %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
#     aes(linetype = .data[["group"]]) +
#     theme_classic() +
#     #facet_wrap(~ group, ncol = 2) +
#     xlab("Stress-induced corticosterone (ng/\U00B5L)") +
#     ylab("Growth (g/day)") +
#     #scale_fill_viridis(discrete = TRUE) +
#     #scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = element_blank())  +
#     scale_x_continuous(trans = abs_trans_tres,
#                        breaks = c((20-mean_abs_tres)/sd_abs_tres,
#                                   (40-mean_abs_tres)/sd_abs_tres,
#                                   (60-mean_abs_tres)/sd_abs_tres,
#                                   (80-mean_abs_tres)/sd_abs_tres,
#                                   (100-mean_abs_tres)/sd_abs_tres)
#     ) +
#     scale_linetype_manual(values = "dotted") +
#     #ylim(-5,5) +
#     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
#     theme(legend.position = "none") +
#     annotate(geom = "text",label = "N = 14",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5)
# )
#
# ggsave("figures/growthbycort+provis_TRES.png",plot =  fig6_abscort_tres, width = 10, height = 6.6)
#
#
# ## Combined WEBL and TRES growth plots
#
#
# # ggplot_build(fig6_provis_webl)$layout$panel_scales_y
# # ggplot_build(fig6_provis_tres)$layout$panel_scales_y
# #
# # ggplot_build(fig6_corts1_webl)$layout$panel_scales_y
# # ggplot_build(fig6_corts1_tres)$layout$panel_scales_y
# #
# # ggplot_build(fig6_abscort_webl)$layout$panel_scales_y
# # ggplot_build(fig6_abscort_tres)$layout$panel_scales_y
# # (p_full <- ggarrange(fig6_provis_webl + theme(text = element_text(size = 12),
# #                                               legend.position = "none",
# #                                               axis.title.x = element_blank(),
# #                                               axis.title.y = element_text(hjust = -3)),
# #                      fig6_provis_tres + theme(axis.text.y = element_blank(),
# #                                               axis.ticks.y = element_blank(),
# #                                               text = element_text(size = 12),
# #                                               axis.title.y = element_blank(),
# #                                               legend.position = 'none',
# #                                               axis.title.x = element_text(hjust = -.7)) +
# #                        ylim(-.183,2.73),
# #                      fig6_corts1_webl + theme(text = element_text(size = 12),
# #                                               legend.position = "none",
# #                                               axis.title.x = element_blank(),
# #                                               axis.title.y = element_blank()),
# #                      fig6_corts1_tres + theme(axis.text.y = element_blank(),
# #                                               axis.ticks.y = element_blank(),
# #                                               text = element_text(size = 12),
# #                                               axis.title.y = element_blank(),
# #                                               legend.position = 'none',
# #                                               axis.title.x = element_text(hjust = -2.7)) +
# #                        ylim(-.183,2.73),
# #                      fig6_abscort_webl + theme(text = element_text(size = 12),
# #                                                legend.position = "none",
# #                                                axis.title.x = element_blank(),
# #                                                axis.title.y = element_blank()),
# #                      fig6_abscort_tres + xlab("Stress-induced - Baseline corticosterone (ng/\U00B5L)") + theme(axis.text.y = element_blank(),
# #                                                                                                                axis.ticks.y = element_blank(),
# #                                                                                                                text = element_text(size = 12),
# #                                                                                                                axis.title.y = element_blank(),
# #                                                                                                                legend.position = 'none',
# #                                                                                                                axis.title.x = element_text(hjust = 2.5)) +
# #                        ylim(-.183,2.73),
# #                      ncol = 2 ,
# #                      labels = c("(a): Western Bluebird","(b): Tree Swallow","","","","")
# # ))
# #
# # ggsave("figures/fig6_growth_by_temp_hab.png",p_full,width = 6.25,height = 8)

