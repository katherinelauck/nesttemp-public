## ================================================================
## SETUP: Libraries and data loading
## ================================================================

library(tidyverse)
library(ggplot2)
library(lubridate)
library(scales)
library(viridis)
library(ggeffects)
library(ggimage)
library(magick)
library(cowplot)
library(emmeans)
library(multcomp)
library(gt)
library(glmmTMB)
library(dagitty)
library(ggdag)
library(igraph)
library(performance)
library(modelsummary)
library(car)
library(lme4)
library(hms)
library(weathermetrics)
library(future)
library(egg)
library(patchwork)

# --- Abbreviations -----------------------------------------------------------
# Species:  WEBL = Western Bluebird, TRES = Tree Swallow
# Cort:     s1 = baseline (pre-stress) corticosterone sample
#           s2 = stress-induced corticosterone sample
#           abs = absolute change (s2 - s1)
# Temp vars:
#   I suffix (e.g. meanmaxtempI) = Fahrenheit / Imperial units
#   hi                           = heat index
#   cum                          = cumulative over the nest/sampling period
#   priorday / priorweek         = rolling window ending the day/week before sampling
#   degreehours_over_30C         = accumulated degree-hours exceeding 30 C
# Habitats: Forest, Orchard, Grassland, Row crop
# Data frames: dat_growth = growth & cort records; dat_provis = provisioning;
#              dat_surv = nest survival
# Scaled vars: original name + _scaled suffix; squared term + _scaled_sq suffix
# -----------------------------------------------------------------------------

dat_provis <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

dat_growth <- read_rds("data/growth_cort_provis_manytempmeasures.rds") %>%
  mutate(year_fct = as.factor(year))

dat_surv <- read_rds("data/survival_attempt.rds") %>%
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


## ================================================================
## HELPER FUNCTIONS
## ================================================================

scale_cols <- function(data, cols) {
  data %>% mutate(across(all_of(cols), ~ scale(.x)[,1], .names = "{.col}_scaled"))
}
# Filter dat to species, drop NA rows for temp_var and min_temp_var, scale
# cols_to_scale, and append a squared column for temp_var (add_sq = TRUE).
prep_model_data <- function(dat, species, temp_var, min_temp_var, cols_to_scale,
                             add_sq = TRUE) {
  d <- dat %>%
    dplyr::filter(
      Species == species,
      !is.na(.data[[temp_var]]),
      !is.na(.data[[min_temp_var]])
    ) %>%
    scale_cols(cols_to_scale)
  if (add_sq) {
    sq_col <- paste0(temp_var, "_scaled_sq")
    d[[sq_col]] <- d[[paste0(temp_var, "_scaled")]]^2
  }
  d
}

make_temp_trans <- function(data, col) {
  m <- mean(data[[col]], na.rm = TRUE)
  s <- sd(data[[col]], na.rm = TRUE)
  list(
    mean = m, sd = s,
    trans = scales::trans_new(col,
      transform = function(x) (x - m) / s,
      inverse   = function(x) x * s + m)
  )
}

anova_int_tab <- function(m_full, m_addmin, m_addmax, m_noint, digits = 4) {
  if (inherits(m_full, "glm") && !inherits(m_full, "lmerMod")) {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint, test = "Chisq") %>% tibble() %>%
        rename(Chisq = Deviance) %>%
        mutate(AIC = c(AIC(m_full), AIC(m_single), AIC(m_noint)), .before = Df) %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = Df)
    }
    p_col    <- "Pr(>Chi)"
    rnd_cols <- "Chisq"
  } else {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint) %>% tibble() %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = npar)
    }
    p_col    <- "Pr(>Chisq)"
    rnd_cols <- c("AIC","Chisq")
  }
  bind_rows(mk_row(m_addmin), mk_row(m_addmax)) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric), ~round(.x, digits = digits)),
           P = .data[[p_col]]) %>%
    dplyr::select(Model, AIC, Chisq, P) %>%
    mutate(max_or_min = c(rep("Max temp", times = 3), rep("Min temp", times = 3)),
           across(all_of(rnd_cols), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001, "<0.001", as.character(P))) %>%
    group_by(max_or_min) %>% gt()
}


## ================================================================
## SECTION 1: GROWTH MODELS
## ================================================================


data_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

(int_tab_growth_webl <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

g_lintemp_webl <- g_lintemp

### Conclusion: growth interacts with max and min.

### Sample sizes:


  samp_year_webl <- g_lintemp@frame %>%
    group_by(habitat,year_fct) %>%
    summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
    rowwise() %>%
    mutate(count = tibble(count_ind,count_attempt)) %>%
    dplyr::select(-c(count_ind, count_attempt)) %>%
    pivot_wider(values_from = count,names_from = year_fct) %>%
    as.tibble() %>% rename(Habitat = 'habitat') %>%
    ungroup() %>%
    mutate(Total = `2021` + `2022` + `2023`) %>%
    rowwise() %>%
    mutate(`2021` = paste0(`2021`$count_ind,", ",`2021`$count_attempt),
           `2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
           `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
           Total = paste0(Total$count_ind,", ",Total$count_attempt)
           )

  (t_samp_year_growth_webl <- samp_year_webl %>% gt() %>%
      grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
        t <- str_split(.,", ",simplify = TRUE)
        t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
        })
    )


  t_samp_year_growth_webl

samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



mean_temp_webl <- mean(data_webl %>% pull(meanmaxtempI))
sd_temp_webl <- sd(data_webl %>% pull(meanmaxtempI))


temp_trans_webl <- trans_new("temp_trans_webl",
                          transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                          inverse = function(x){x})

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
                       ) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )

## visualize temp * julian date interaction

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
    ) #+
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


(weblgrowthtrendmax <- emtrends(g_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
   gt())


data <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
   mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
          across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
   # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
   mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
          across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


## min temp


data <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


mean_temp <- mean(data %>% pull(meanmintempI))
sd_temp <- sd(data %>% pull(meanmintempI))


temp_trans <- trans_new("temp_trans",
                          transform = function(x){(x * sd_temp) + mean_temp},
                          inverse = function(x){x})

(pl <- ggpredict(g_lintemp,terms = c("meanmintempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily min temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans,
                       breaks = c((8-mean_temp)/sd_temp,
                                  (12-mean_temp)/sd_temp,
                                  (16-mean_temp)/sd_temp,
                                  (20-mean_temp)/sd_temp),
                       ) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


(t <- emtrends(g_lintemp,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
   gt())


data <- prep_model_data(
  dat_growth, "WEBL", "meanmintempI", "meanmintempI",
  c("gweight","meanmintempI","meanmintempI","juliandate","brood_size","age"))

(t <- emmeans(g_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
  mutate(meanmintempI_scaled = (meanmintempI_scaled * sd(data$meanmintempI)) + mean(data$meanmintempI),
    across(where(is.numeric), ~ round(.x, digits = 2)),
    meanmintempI_scaled = round(meanmintempI_scaled),
    meanmintempI_scaled = paste0(meanmintempI_scaled,"\u00b0C")) %>%
  #tibble() %>%
  dplyr::select(-df) %>%
  #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
  rename(Habitat = "habitat",`min temperature` = "meanmintempI_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
  group_by(`min temperature`) %>%
  mutate(row=row_number()) %>%
  pivot_longer(-c(`min temperature`, row,Habitat)) %>%
  pivot_wider(names_from=c(`min temperature`, name), values_from=value) %>%
  dplyr::select(-row) %>%
  #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
  #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
  gt() %>% tab_options(data_row.padding = px(1)) %>%
  tab_spanner_delim(
    delim="_"
  ))


((emmeans(g_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))

emmeans(g_lintemp,specs = pairwise ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
emmip(g_lintemp,formula = habitat ~ meanmintempI_scaled, at = list(meanmintempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()


#### Check whether controlling for brood size explains the habitat * temp differences.


## create brood size by week measure


g_lintemp_addmin_tres <- g_lintemp_addmin


### Conclusion: growth interacts with max but not min.


#### Check whether controlling for brood size explains the habitat * temp differences.


## create brood size by week measure

data_maxhiweek_webl <- prep_model_data(
  dat_growth, "WEBL", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiweek_webl)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiweek_webl)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiweek_webl)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiweek_webl)

(int_tab_growth_maxhiweek_webl <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_maxhiweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblgrowthtrendmaxhiweek <- emtrends(g_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiweek_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_maxhiweek_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


#### Check whether controlling for brood size explains the habitat * temp differences.


## create brood size by week measure

data_maxhiday_webl <- prep_model_data(
  dat_growth, "WEBL", "maxhi_prior", "meanmintempI",
  c("gweight","maxhi_prior","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiday_webl)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiday_webl)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiday_webl)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_maxhiday_webl)

(int_tab_growth_maxhiday_webl <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_maxhiday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblgrowthtrendmaxhiday <- emtrends(g_lintemp,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiday_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_maxhiday_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


#### degreehours_over_30C_priorweek


data_deghr30week_webl <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_deghr30week_webl)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_deghr30week_webl)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_deghr30week_webl)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_deghr30week_webl)

(int_tab_growth_deghr30week_webl <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_deghr30week_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblgrowthtrenddeghr30week <- emtrends(g_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data_deghr30week_webl <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_deghr30week_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_deghr30week_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())

#### hihours_over_30hi_priorweek


data_hihr25week_webl <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_hihr25week_webl)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_hihr25week_webl)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_hihr25week_webl)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_hihr25week_webl)

(int_tab_growth_hihr25week_webl <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_hihr25week_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblgrowthtrendhihr25week <- emtrends(g_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data_hihr25week_webl <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_hihr25week_webl <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_hihr25week_webl <- summary(g_lintemp) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


### TRES
#### maxhi_week

data_maxhiweek_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiweek_tres)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiweek_tres)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiweek_tres)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiweek_tres)

(int_tab_growth_maxhiweek_tres <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_maxhiweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresgrowthtrendmaxhiweek <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data_maxhiweek_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiweek_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_maxhiweek_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


#### maxhi_prior


data_maxhiday_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_prior", "meanmintempI",
  c("gweight","maxhi_prior","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiday_tres)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiday_tres)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiday_tres)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ maxhi_prior_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_maxhiday_tres)

(int_tab_growth_maxhiday_tres <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_maxhiday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresgrowthtrendmaxhiday <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data_maxhiday_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_prior", "meanmintempI",
  c("gweight","maxhi_prior","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_maxhiday_tres <- emmeans(g_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_maxhiday_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


#### degreehours_over_30C_priorweek


data_deghr30week_tres <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_deghr30week_tres)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_deghr30week_tres)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_deghr30week_tres)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_deghr30week_tres)

(int_tab_growth_deghr30week_tres <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_deghr30week_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresgrowthtrenddeghr30week <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data_deghr30week_tres <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_deghr30week_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_deghr30week_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())

#### hihours_over_30hi_priorweek


data_hihr25week_tres <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

g_lintemp <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_hihr25week_tres)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_hihr25week_tres)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_hihr25week_tres)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_hihr25week_tres)

(int_tab_growth_hihr25week_tres <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))

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
    ) +
    geom_text(data = dat_text_hihr25week_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresgrowthtrendhihr25week <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data_hihr25week_tres <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))


## Emmeans to check for effect of habitat


(growthbyhabitat_hihr25week_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(growthbyhabitat_summary_hihr25week_tres <- summary(g_lintemp_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())


### TRES


g_lintemp <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
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
# gtsave(int_tab_growth_tres,"figures/int_tab_growth_tres.html")


check_collinearity(g_lintemp_noint)

summary(g_lintemp_noint)

g_lintemp_addmin_tres <- g_lintemp_addmin


### Conclusion: growth interacts with max but not min.



#### Check whether controlling for brood size explains the habitat * temp differences.

#Because higher temps reduce survival, especially in forest, we are thinking that brood size is the outcome of the habitat * temp interaction. So, if we include brood size and the relationship growth ~ habitat * temp disappears, we think that's what's going on.

## create brood size by week measure


# brood_size <- b %>% group_by(`Banding Date`,Nestbox) %>% summarize(brood_size_weekly = n()) %>% ungroup() %>% mutate(date = mdy(`Banding Date`)) %>% dplyr::select(!`Banding Date`)
#
# g <- left_join(g %>% ungroup(),brood_size, by = c("date","Nestbox"))
#
#
# If this is so, then the main driver of growth and survival may be predation.
#
#
# g_lintemp <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + brood_size_weekly_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# g_lintemp_addmax <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + brood_size_weekly_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# g_lintemp_addmin <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + brood_size_weekly_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
# g_lintemp_noint <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + brood_size_weekly_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
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
# gtsave(tab,"../figures/int_tab_growth_with_broodsize_TRES.html")
#
# summary(g_lintemp)
#
#
# Is brood_size significant if I remove temp and habitat?
#
#
# g_lintemp <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + brood_size_weekly_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
#
# summary(g_lintemp)
#
# g_lintemp <- lmerTest::lmer(gweight ~ age_scaled + brood_size_weekly_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
#   mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size_weekly,age),
#                                                   ~ scale(.x)[,1],
#                                                    .names = "{.col}_scaled"),
#                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
#
#
#
# summary(g_lintemp)


### Sample sizes:


# samp_year_tres <- g_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
# (t_samp_year_growth_tres <- samp_year_tres %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)))
#   t_samp_year_growth_tres



  samp_year_tres <- g_lintemp_addmin@frame %>%
    group_by(habitat,year_fct) %>%
    summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
    rowwise() %>%
    mutate(count = tibble(count_ind,count_attempt)) %>%
    dplyr::select(-c(count_ind, count_attempt)) %>%
    pivot_wider(values_from = count,names_from = year_fct) %>%
    as.tibble() %>% rename(Habitat = 'habitat') %>%
    ungroup() %>%
    mutate(Total = `2021` + `2022` + `2023`) %>%
    rowwise() %>%
    mutate(`2021` = paste0(`2021`$count_ind,", ",`2021`$count_attempt),
           `2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
           `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
           Total = paste0(Total$count_ind,", ",Total$count_attempt)
    )

  (t_samp_year_growth_tres <- samp_year_tres %>% gt() %>%
      grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
        t <- str_split(.,", ",simplify = TRUE)
        t[,1] %>% as.numeric() %>% sum() %>% paste(t[,2] %>% as.numeric() %>% sum(),sep = ", ")
      })
  )

  t_samp_year_growth_tres

samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


#### Combined sample sizes for growth



# rbind(t_samp_year_growth_webl$`_data`,t_samp_year_growth_tres$`_data`) %>% dplyr::select(Habitat, `2021`, `2022`, `2023`) %>%
#   mutate(Species = c(rep("WEBL",4),rep("TRES",4)),
#          Response = "Growth") %>%
#   # mutate(row=row_number()) %>%
#   group_by(Response) %>%
#   pivot_longer(-c(Species, Response, Habitat)) %>%
#   pivot_wider(names_from=c(Species, name), values_from=value) %>%
#   # dplyr::select(-row) %>%
#   gt() %>% tab_options(data_row.padding = px(1)) %>%
#   tab_spanner_delim(
#     delim="_")






data_tres = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp_tres <- mean(data_tres %>% pull(meanmaxtempI))
sd_temp_tres <- sd(data_tres %>% pull(meanmaxtempI))


temp_trans_tres <- trans_new("temp_trans",
                          transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                          inverse = function(x){x})




summary(g_lintemp_addmin)



(tresgrowthtrendmax <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
   gt())


# gtsave(tresgrowthtrendmax,"figures/tresgrowthtrendmax.html")

data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

# (t <- emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#   mutate(meanmaxtempI_scaled = (meanmaxtempI_scaled * sd(data$meanmaxtempI)) + mean(data$meanmaxtempI),
#     across(where(is.numeric), ~ round(.x, digits = 2)),
#     meanmaxtempI_scaled = round(meanmaxtempI_scaled),
#     meanmaxtempI_scaled = paste0(meanmaxtempI_scaled,"\u00b0C")) %>%
#   #tibble() %>%
#   dplyr::select(-df) %>%
#   #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#   rename(Habitat = "habitat",`Max temperature` = "meanmaxtempI_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
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

# ((emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))
#
# emmeans(g_lintemp_addmin,specs = pairwise ~ habitat,by = c("meanmaxtempI_scaled"), at = list(meanmaxtempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(g_lintemp_addmin,formula = habitat ~ meanmaxtempI_scaled, at = list(meanmaxtempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(growthbyhabitat_tres <- emmeans(g_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
   mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
          across(p.value,~round(.x,digits = 3))) %>% gt())
# gtsave(growthbyhabitat_tres,"figures/growthbyhabitat_tres.html")


## min temp



data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


mean_temp <- mean(data %>% pull(meanmintempI))
sd_temp <- sd(data %>% pull(meanmintempI))


temp_trans <- trans_new("temp_trans",
                          transform = function(x){(x * sd_temp) + mean_temp},
                          inverse = function(x){x})

(pl <- ggpredict(g_lintemp_addmin,terms = c("meanmintempI_scaled [all]","habitat"),bias_correction = TRUE) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily min temp over preceding week (\u00b0C)") +
    ylab("Growth (g/day)") +
   scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans,
                       breaks = c((8-mean_temp)/sd_temp,
                                  (12-mean_temp)/sd_temp,
                                  (16-mean_temp)/sd_temp,
                                  (20-mean_temp)/sd_temp),
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmintempI,na.rm = TRUE))/sd(g$meanmintempI,na.rm = TRUE),
                       #            (55-mean(g$meanmintempI,na.rm = TRUE))/sd(g$meanmintempI,na.rm = TRUE))
                       # limits = c(18,5)
                       ) +
   #ylim(-5,5) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )

# ggsave("figures/growthbymintempxhab_TRES.png",plot =  pl, width = 10, height = 6.6)



summary(g_lintemp_addmin)



(t <- emtrends(g_lintemp_addmin,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
   gt())


# gtsave(t,"figures/tresgrowthtrendmin.html")

data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmintempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmintempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmintempI_scaled_sq = meanmintempI_scaled * meanmintempI_scaled)

(t <- emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
  mutate(meanmintempI_scaled = (meanmintempI_scaled * sd(data$meanmintempI)) + mean(data$meanmintempI),
    across(where(is.numeric), ~ round(.x, digits = 2)),
    meanmintempI_scaled = round(meanmintempI_scaled),
    meanmintempI_scaled = paste0(meanmintempI_scaled,"\u00b0C")) %>%
  #tibble() %>%
  dplyr::select(-df) %>%
  #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
  rename(Habitat = "habitat",`min temperature` = "meanmintempI_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
  group_by(`min temperature`) %>%
  mutate(row=row_number()) %>%
  pivot_longer(-c(`min temperature`, row,Habitat)) %>%
  pivot_wider(names_from=c(`min temperature`, name), values_from=value) %>%
  dplyr::select(-row) %>%
  #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
  #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
  gt() %>% tab_options(data_row.padding = px(1)) %>%
  tab_spanner_delim(
    delim="_"
  ))

# gtsave(t,"figures/tresgrowthdeltamin.html")

((emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_lintemp_addmin,specs = ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))

emmeans(g_lintemp_addmin,specs = pairwise ~ habitat,by = c("meanmintempI_scaled"), at = list(meanmintempI_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
emmip(g_lintemp_addmin,formula = habitat ~ meanmintempI_scaled, at = list(meanmintempI_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Combined WEBL and TRES growth plots





# Models to estimate the unbiased causal effect of cort and provisioning on growth

#Model structure pulled from DAG: will use a minimal adjustment set to estimate the total effect of cort and provisioning on growth.
#{ attempt, humidity, nest_age, nestcond, surftemp }

#I could subtract the direct from the total effect to get the indirect effect. I'm still vacillating between thinking that I should be measuring the total effect or the direct effect. I think because nestcond is really the only mediator in the DAG that I actually want total effect?

#  The question I'd like to ask: Does the direct effect of these two factors on growth change depending on habitat? I'm not sure that this causal structure is the correct one for that question. Is that a mediation analysis? After reading some papers, I think that actually the question is, do cort and provisioning mediate the effect of habitat and temperature on growth? I still think there may be other questions in this dataset. Maybe the answer is to do this paper just linear modeling with causal inference support, and then get Sage in on a paper that delves more deeply into a structural causal modeling framework.

#The following focuses on the more simple question: how do cort and provisioning affect growth (total effect)? I will do the same model structure for survival in that notebook.

#adjustment set: { attempt, humidity, nest_age, nestcond, surftemp }



g_provis_cort <- lmerTest::lmer(gweight ~ cort_s1_scaled + provis_mean_scaled + poly(meanmaxtempI_scaled,2) + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "WEBL",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
                                  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                                                ~ scale(.x)[,1],
                                                .names = "{.col}_scaled"),
                                         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_provis_cort_webl <- g_provis_cort
check_collinearity(g_provis_cort)
summary(g_provis_cort)

(growth_by_cort_provis_summary_webl <- summary(g_provis_cort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_cort_provis_summary_webl,"figures/growth_by_cort_provis_summary_webl.html")

data = dplyr::filter(dat_growth,Species == "WEBL",!is.na(gweight),!is.na(cort_s1),!is.na(provis_mean),!is.na(meanmaxtempI),
                     !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

# (t <- emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(provis_mean_scaled = (provis_mean_scaled * sd(data$provis_mean, na.rm = TRUE)) + mean(data$provis_mean,na.rm = TRUE),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            provis_mean_scaled = round(provis_mean_scaled)) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(`Provisioning` = "provis_mean_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)))
# gtsave(t,"figures/growth_by_basecort_provis_delta_webl.html")
#
# ((emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_provis_cort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))





g_provis_abscort <- lmerTest::lmer(gweight ~ abs_change_cort_scaled + provis_mean_scaled + meanmaxtempI_scaled + meanmintempI_scaled + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "WEBL",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
g_provis_abscort_webl <- g_provis_abscort

summary(g_provis_abscort)

(growth_by_abscort_provis_summary_webl <- summary(g_provis_abscort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_abscort_provis_summary_webl,"figures/growth_by_abscort_provis_summary_webl.html")

data = dplyr::filter(dat_growth,Species == "WEBL",!is.na(gweight),!is.na(abs_change_cort),!is.na(provis_mean),!is.na(meanmaxtempI),
                     !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

(t <- emmeans(g_provis_abscort,specs = ~ provis_mean_scaled,by = c("provis_mean_scaled"), at = list(provis_mean_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
    mutate(provis_mean_scaled = (provis_mean_scaled * sd(data$provis_mean,na.rm = TRUE)) + mean(data$provis_mean,na.rm = TRUE),
           across(where(is.numeric), ~ round(.x, digits = 2)),
           provis_mean_scaled = round(provis_mean_scaled)) %>%
    #tibble() %>%
    dplyr::select(-df) %>%
    #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
    rename(`max temperature` = "provis_mean_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
    #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
    #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
    gt() %>% tab_options(data_row.padding = px(1)))
# gtsave(t,"figures/growth_by_abscort_provis_delta_webl.html")



#### Sample size


data_s1_webl = dplyr::filter(dat_growth,Species == "WEBL",!is.na(gweight),!is.na(cort_s1),!is.na(provis_mean),!is.na(meanmaxtempI),
                          !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)
# ss_year_webl_growthbybasecortprovis <- data_s1_webl %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_webl_growthbybasecortprovis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_webl_growthbybasecortprovis.html")

ss_year_webl_growthbybasecortprovis <- data_s1_webl %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total =`2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_webl_growthbybasecortprovis <- ss_year_webl_growthbybasecortprovis %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    })
)


t_ss_year_webl_growthbybasecortprovis





data_abs_webl = dplyr::filter(dat_growth,Species == "WEBL",!is.na(gweight),!is.na(abs_change_cort),!is.na(provis_mean),!is.na(meanmaxtempI),
                     !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)
# ss_year_webl_growthbyabscortprovis <- data_abs_webl %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_webl_growthbyabscortprovis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_webl_growthbyabscortprovis.html")


ss_year_webl_growthbyabscortprovis <- data_abs_webl %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total =`2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_webl_growthbyabscortprovis <- ss_year_webl_growthbyabscortprovis %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    })
)


t_ss_year_webl_growthbyabscortprovis





mean_provis_webl <- mean(data_s1_webl %>% pull(provis_mean))
sd_provis_webl <- sd(data_s1_webl %>% pull(provis_mean))


provis_trans_webl <- trans_new("provis_trans_webl",
                               transform = function(x){(x * sd_provis_webl) + mean_provis_webl},
                               inverse = function(x){x})

mean_s1_webl <- mean(data_s1_webl %>% pull(cort_s1))
sd_s1_webl <- sd(data_s1_webl %>% pull(cort_s1))


s1_trans_webl <- trans_new("s1_trans_webl",
                               transform = function(x){(x * sd_s1_webl) + mean_s1_webl},
                               inverse = function(x){x})

mean_abs_webl <- mean(data_abs_webl %>% pull(abs_change_cort))
sd_abs_webl <- sd(data_abs_webl %>% pull(abs_change_cort))


abs_trans_webl <- trans_new("abs_trans_webl",
                               transform = function(x){(x * sd_abs_webl) + mean_abs_webl},
                               inverse = function(x){x})




## Canonical adjustment set overloads the data:

#{ attempt, broodsize, habitat, humidity, jul_date, mother_s1, nest_age, site, surftemp, year }



# g_provis_cort_canon <- lmerTest::lmer(gweight ~ cort_s1_scaled + provis_mean_scaled + (1|attempt_id) + brood_size_scaled + habitat + meanh_scaled + juliandate_scaled + cort_s1_mother_scaled + age_scaled + meanmaxtempI_scaled + meanmintempI_scaled + year_fct,data = dplyr::filter(dat_growth,Species == "WEBL",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
#                                         mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,cort_s1_mother,condition,provis_mean),
#                                                       ~ scale(.x)[,1],
#                                                       .names = "{.col}_scaled"),
#                                                meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
# summary(g_provis_cort_canon)
#
#

## TRES



g_provis_cort <- lmerTest::lmer(gweight ~ cort_s1_scaled + provis_mean_scaled + poly(meanmaxtempI_scaled,2) + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                                                ~ scale(.x)[,1],
                                                .names = "{.col}_scaled"),
                                         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
check_collinearity(g_provis_cort)
summary(g_provis_cort)

g_provis_cort_tres <- g_provis_cort

(growth_by_cort_provis_summary_tres <- summary(g_provis_cort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_cort_provis_summary_tres,"figures/growth_by_cort_provis_summary_tres.html")

data = dplyr::filter(dat_growth,Species == "TRES",!is.na(gweight),!is.na(cort_s1),!is.na(provis_mean),!is.na(meanmaxtempI),
                     !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

# (t <- emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(cort_s1_scaled = (cort_s1_scaled * sd(data$cort_s1, na.rm = TRUE)) + mean(data$cort_s1,na.rm = TRUE),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            cort_s1_scaled = round(cort_s1_scaled)) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(`Baseline cort` = "cort_s1_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
#     #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)))
# gtsave(t,"figures/growth_by_basecort_provis_delta_tres.html")
#
# ((emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_provis_cort,specs = ~ cort_s1_scaled,by = c("cort_s1_scaled"), at = list(cort_s1_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))



#### Abs diff cort



g_provis_abscort <- lmerTest::lmer(gweight ~ abs_change_cort_scaled + provis_mean_scaled + meanmaxtempI_scaled + meanmintempI_scaled + meanh_scaled + age_scaled + condition_scaled + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),!is.na(meanmintempI)) %>%
                                     mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
                                                   ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                            meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
summary(g_provis_cort)
g_provis_abscort_tres <- g_provis_abscort

(growth_by_abscort_provis_summary_tres <- summary(g_provis_abscort) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
           across(`Pr(>|t|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(growth_by_abscort_provis_summary_tres,"figures/growth_by_abscort_provis_summary_tres.html")

data = dplyr::filter(dat_growth,Species == "TRES",!is.na(gweight),!is.na(abs_change_cort),!is.na(provis_mean),!is.na(meanmaxtempI),
                     !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,abs_change_cort,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

(t <- emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
    mutate(abs_change_cort_scaled = (abs_change_cort_scaled * sd(data$abs_change_cort, na.rm = TRUE)) + mean(data$abs_change_cort,na.rm = TRUE),
           across(where(is.numeric), ~ round(.x, digits = 2)),
           abs_change_cort_scaled = round(abs_change_cort_scaled)) %>%
    #tibble() %>%
    dplyr::select(-df) %>%
    #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("minimum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
    rename(`Baseline cort` = "abs_change_cort_scaled",`Predicted growth` = "emmean",`2.5%` = "lower.CL",`97.5%` = "upper.CL" ) %>%
    #mutate(`minimum TA_P` = if_else(`minimum TA_P` == 0.000,"<0.001",as.character(`minimum TA_P`))) %>%
    #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
    gt() %>% tab_options(data_row.padding = px(1)))
# gtsave(t,"figures/growth_by_abscort_provis_delta_tres.html")

((emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(emmean))-(emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean)))/(emmeans(g_provis_abscort,specs = ~ abs_change_cort_scaled,by = c("abs_change_cort_scaled"), at = list(abs_change_cort_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(emmean))



#### Sample size


data_s1_tres = dplyr::filter(dat_growth,Species == "TRES",!is.na(gweight),!is.na(cort_s1),!is.na(provis_mean),!is.na(meanmaxtempI),
                          !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


# ss_year_TRES_growthbybasecortprovis <- data_s1_tres %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_TRES_growthbybasecortprovis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_TRES_growthbybasecortprovis.html")
#


ss_year_tres_growthbybasecortprovis <- data_s1_tres %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total =`2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_tres_growthbybasecortprovis <- ss_year_tres_growthbybasecortprovis %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    })
)


t_ss_year_tres_growthbybasecortprovis




data_abs_tres = dplyr::filter(dat_growth,Species == "TRES",!is.na(gweight),!is.na(abs_change_cort),!is.na(provis_mean),!is.na(meanmaxtempI),
                     !is.na(meanmintempI),!is.na(meanh),!is.na(age),!is.na(condition)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age,meanh,cort_s1,condition,provis_mean),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)


# ss_year_TRES_growthbyabscortprovis <- data_abs_tres %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   # relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2022`+`2023`)
# ss_year_TRES_growthbyabscortprovis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_TRES_growthbyabscortprovis.html")


ss_year_tres_growthbyabscortprovis <- data_abs_tres %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(Total =`2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_tres_growthbyabscortprovis <- ss_year_tres_growthbyabscortprovis %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    })
)


t_ss_year_tres_growthbyabscortprovis


mean_provis_tres <- mean(data_s1_tres %>% pull(provis_mean))
sd_provis_tres <- sd(data_s1_tres %>% pull(provis_mean))


provis_trans_tres <- trans_new("provis_trans_tres",
                               transform = function(x){(x * sd_provis_tres) + mean_provis_tres},
                               inverse = function(x){x})

mean_s1_tres <- mean(data_s1_tres %>% pull(cort_s1))
sd_s1_tres <- sd(data_s1_tres %>% pull(cort_s1))


s1_trans_tres <- trans_new("s1_trans_tres",
                               transform = function(x){(x * sd_s1_tres) + mean_s1_tres},
                               inverse = function(x){x})

mean_abs_tres <- mean(data_abs_tres %>% pull(abs_change_cort))
sd_abs_tres <- sd(data_abs_tres %>% pull(abs_change_cort))


abs_trans_tres <- trans_new("abs_trans_tres",
                           transform = function(x){(x * sd_abs_tres) + mean_abs_tres},
                           inverse = function(x){x})

# (fig6_tres <- predict_response(g_provis_cort,terms = c("provis_mean_scaled","cort_s1_scaled [-2,2]"),bias_correction = TRUE) %>%
#     #filter(group != 0.35) %>%
#     # mutate(group = as_factor(group) %>% fct_collapse(Low = "-0.93",High = "1.63")) %>%
#     plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE,colors = "viridis") +
#     #scale_fill_discrete(labels = c("Low","High")) +
#     #scale_color_discrete(labels = c("Low","High")) +
#     # geom_point(show.legend = FALSE) +
#     theme_classic() +
#     #facet_wrap(~ group, ncol = 2) +
#     xlab("Provisioning (mean visits/hr)") +
#     ylab("Growth (g/day)") +
#     #scale_fill_viridis(discrete = TRUE) +
#     #scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     labs(title = element_blank(),color = "Baseline cort") +
#     scale_x_continuous(trans = provis_trans_tres,
#                        breaks = c((0-mean_provis_tres)/sd_provis_tres,
#                                   (10-mean_provis_tres)/sd_provis_tres,
#                                   (20-mean_provis_tres)/sd_provis_tres,
#                                   (30-mean_provis_tres)/sd_provis_tres,
#                                   (40-mean_provis_tres)/sd_provis_tres,
#                                   (50-mean_provis_tres)/sd_provis_tres)
#     ) +
#     #ylim(-5,5) #+
#     # geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) # +
#     theme(legend.position = "none") +
#     annotate("text",x = -Inf,y = -Inf,hjust = -.2,vjust = -.5,label = "N = 16")
# )





save(list = ls(), file = "data/models_growth.RData")
rm(list = ls()); gc()
dat_growth <- read_rds("data/growth_cort_provis_manytempmeasures.rds") %>%
  mutate(year_fct = as.factor(year))
scale_cols <- function(data, cols) {
  data %>% mutate(across(all_of(cols), ~ scale(.x)[,1], .names = "{.col}_scaled"))
}
prep_model_data <- function(dat, species, temp_var, min_temp_var, cols_to_scale,
                             add_sq = TRUE) {
  d <- dat %>%
    dplyr::filter(Species == species, !is.na(.data[[temp_var]]), !is.na(.data[[min_temp_var]])) %>%
    scale_cols(cols_to_scale)
  if (add_sq) { sq_col <- paste0(temp_var, "_scaled_sq"); d[[sq_col]] <- d[[paste0(temp_var, "_scaled")]]^2 }
  d
}
anova_int_tab <- function(m_full, m_addmin, m_addmax, m_noint, digits = 4) {
  if (inherits(m_full, "glm") && !inherits(m_full, "lmerMod")) {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint, test = "Chisq") %>% tibble() %>%
        rename(Chisq = Deviance) %>%
        mutate(AIC = c(AIC(m_full), AIC(m_single), AIC(m_noint)), .before = Df) %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = Df)
    }
    p_col <- "Pr(>Chi)"; rnd_cols <- "Chisq"
  } else {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint) %>% tibble() %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = npar)
    }
    p_col <- "Pr(>Chisq)"; rnd_cols <- c("AIC","Chisq")
  }
  bind_rows(mk_row(m_addmin), mk_row(m_addmax)) %>% as_tibble() %>%
    mutate(across(where(is.numeric), ~round(.x, digits = digits)), P = .data[[p_col]]) %>%
    dplyr::select(Model, AIC, Chisq, P) %>%
    mutate(max_or_min = c(rep("Max temp", times = 3), rep("Min temp", times = 3)),
           across(all_of(rnd_cols), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001, "<0.001", as.character(P))) %>%
    group_by(max_or_min) %>% gt()
}


## ================================================================
## SECTION 2: CORTICOSTERONE MODELS
## ================================================================

data_s1_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_webl)

(int_tab_s1_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_webl <- s1_lintemp_noint


## Emmeans to check for effect of habitat


(s1byhabitat_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


data <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


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
    ) +
    geom_text(data = dat_text_s1_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

(int_tab_abs_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl <- abs_lintemp


### sample size


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    ylim(0,60) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


data <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


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
    ) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabstrendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


### Minimum temperature


(t <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data <- prep_model_data(
  dat_growth, "WEBL", "meanmintempI", "meanmintempI",
  c("gweight","meanmintempI","meanmintempI","juliandate","brood_size","age"))


#### s2_cort
data_s2_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_webl)

(int_tab_s2_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_webl <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### use prior day temp to predict cort instead


data_s1_priordayt_webl <- prep_model_data(
  dat_growth, "WEBL", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_webl)

(int_tab_s1_priordayt_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_prior_day_webl <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_priordayt_webl <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_priordayt_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_priordayt <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_priordayt_webl <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordayt_webl %>% gt()


dat_text_s1_priordayt_webl <- data.frame(
  label = paste("N =",samp_s1_priordayt_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### use prior day temp to predict cort instead


data_s2_priordayt_webl <- prep_model_data(
  dat_growth, "WEBL", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_webl)

(int_tab_s2_priordayt_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_prior_day_webl <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_priordayt_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_priordayt_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_priordayt <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_priordayt_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordayt_webl %>% gt()


dat_text_s2_priordayt_webl <- data.frame(
  label = paste("N =",samp_s2_priordayt_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl_priordayt <- prep_model_data(
  dat_growth, "WEBL", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordayt)

(int_tab_abs_priordayt_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_priordayt <- abs_lintemp_noint


### sample size


ss_year_webl_abs_priordayt <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_priordayt_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use prior day heat index to predict cort instead


data_webl_priordaymaxhhi <- prep_model_data(
  dat_growth, "WEBL", "maxhi_prior", "minhi_prior",
  c("gweight","maxhi_prior","minhi_prior","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordaymaxhhi)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordaymaxhhi)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordaymaxhhi)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled + minhi_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordaymaxhhi)

(int_tab_s1_priordaymaxhhi_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_priordaymaxhhi_webl <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_priordaymaxhhi_webl <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_priordaymaxhhi_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_priordaymaxhhi <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_priordaymaxhhi_webl <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordaymaxhhi_webl %>% gt()


dat_text_s1_priordaymaxhhi_webl <- data.frame(
  label = paste("N =",samp_s1_priordaymaxhhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


data_s1_priordaymaxhhi_webl <- prep_model_data(
  dat_growth, "WEBL", "maxhi_prior", "mint_prior",
  c("gweight","maxhi_prior","mint_prior","juliandate","brood_size","age"))


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_priordaymaxhhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2 cort

data_webl_priordaymaxhhi <- prep_model_data(
  dat_growth, "WEBL", "maxhi_prior", "minhi_prior",
  c("gweight","maxhi_prior","minhi_prior","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordaymaxhhi)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordaymaxhhi)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordaymaxhhi)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled + minhi_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordaymaxhhi)

(int_tab_s2_priordaymaxhhi_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_priordaymaxhhi_webl <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_priordaymaxhhi_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_priordaymaxhhi_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_priordaymaxhhi <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_priordaymaxhhi_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordaymaxhhi_webl %>% gt()


dat_text_s2_priordaymaxhhi_webl <- data.frame(
  label = paste("N =",samp_s2_priordaymaxhhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


data_s2_priordaymaxhhi_webl <- prep_model_data(
  dat_growth, "WEBL", "maxhi_prior", "mint_prior",
  c("gweight","maxhi_prior","mint_prior","juliandate","brood_size","age"))


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_priordaymaxhhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl_priordaymaxhhi <- prep_model_data(
  dat_growth, "WEBL", "maxhi_prior", "minhi_prior",
  c("gweight","maxhi_prior","minhi_prior","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordaymaxhhi)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled + minhi_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordaymaxhhi)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + minhi_prior_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordaymaxhhi)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled + minhi_prior_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_priordaymaxhhi)

(int_tab_abs_priordaymaxhhi_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_priordaymaxhhi <- abs_lintemp_noint


### sample size


ss_year_webl_abs_priordaymaxhhi <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_priordaymaxhhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_priordaymaxhhi_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use prior week heat index to predict cort instead


data_s1_weekhi_webl <- prep_model_data(
  dat_growth, "WEBL", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_weekhi_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_weekhi_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_weekhi_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_weekhi_webl)

(int_tab_s1_weekhi_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_weekhi_webl <- s1_lintemp_noint


## Emmeans to check for effect of habitat


(s1byhabitat_weekhi_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_weekhi_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_weekhi <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_weekhi_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_weekhi_webl %>% gt()


dat_text_s1_weekhi_webl <- data.frame(
  label = paste("N =",samp_s1_weekhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_weekhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2 cort


data_s2_weekhi_webl <- prep_model_data(
  dat_growth, "WEBL", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_weekhi_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_weekhi_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_weekhi_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_weekhi_webl)

(int_tab_s2_weekhi_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_weekhi_webl <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_weekhi_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_weekhi_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_weekhi <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_weekhi_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_weekhi_webl %>% gt()


dat_text_s2_weekhi_webl <- data.frame(
  label = paste("N =",samp_s2_weekhi_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_weekhi_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl_weekhi <- prep_model_data(
  dat_growth, "WEBL", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_weekhi)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_weekhi)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_weekhi)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl_weekhi)

(int_tab_abs_weekhi_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_weekhi <- abs_lintemp


### sample size


ss_year_webl_abs_weekhi <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl_weekhi, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_weekhi_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior day hi to predict cort


data_s1_cumhiday_webl <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorday", "meanmintempI",
  c("gweight","hihours_over_30hi_priorday","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumhiday_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiday_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiday_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiday_webl)

(int_tab_s1_cumhiday_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumhiday_webl <- s1_lintemp_noint


## Emmeans to check for effect of habitat


(s1byhabitat_cumhiday_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_cumhiday_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_cumhiday <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumhiday_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiday_webl %>% gt()


dat_text_s1_cumhiday_webl <- data.frame(
  label = paste("N =",samp_s1_cumhiday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumhiday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### use cumulative prior day hi to predict cort


data_s2_cumhiday_webl <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorday", "meanmintempI",
  c("gweight","hihours_over_30hi_priorday","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumhiday_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiday_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiday_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiday_webl)

(int_tab_s2_cumhiday_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumhiday_webl <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_cumhiday_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_cumhiday_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_cumhiday <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumhiday_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiday_webl %>% gt()


dat_text_s2_cumhiday_webl <- data.frame(
  label = paste("N =",samp_s2_cumhiday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_cumhiday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl_cumhiday <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorday", "meanmintempI",
  c("gweight","hihours_over_30hi_priorday","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumhiday)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumhiday)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumhiday)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumhiday)

(int_tab_abs_cumhiday_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_cumhiday <- abs_lintemp_noint


### sample size


ss_year_webl_abs_cumhiday <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "solid")) +
    geom_text(data = dat_text_webl_cumhiday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_cumhiday_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior week hi to predict cort


data_s1_cumhiweek_webl <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_webl)

(int_tab_s1_cumhiweek_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumhiweek_webl <- s1_lintemp_noint


## Emmeans to check for effect of habitat


(s1byhabitat_cumhiweek_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_cumhiweek_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_cumhiweek <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumhiweek_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiweek_webl %>% gt()


dat_text_s1_cumhiweek_webl <- data.frame(
  label = paste("N =",samp_s1_cumhiweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumhiweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_cumhiweek_webl <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_webl)

(int_tab_s2_cumhiweek_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumhiweek_webl <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_cumhiweek_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_cumhiweek_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_cumhiweek <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumhiweek_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiweek_webl %>% gt()


dat_text_s2_cumhiweek_webl <- data.frame(
  label = paste("N =",samp_s2_cumhiweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dashed","solid","dotted")) +
    geom_text(data = dat_text_s2_cumhiweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_abs_webl_cumhiweek <- prep_model_data(
  dat_growth, "WEBL", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_webl_cumhiweek)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_webl_cumhiweek)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_webl_cumhiweek)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_webl_cumhiweek)

(int_tab_abs_cumhiweek_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_cumhiweek <- abs_lintemp


### sample size


ss_year_webl_abs_cumhiweek <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_abs_text_webl_cumhiweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_cumhiweek_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior day temp to predict cort


data_s1_cumdegreeday_webl <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorday", "meanmintempI",
  c("gweight","degreehours_over_30C_priorday","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_webl)

(int_tab_s1_cumdegreeday_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumdegreeday_webl <- s1_lintemp_noint


## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeday_webl <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_cumdegreeday_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_cumdegreeday <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumdegreeday_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeday_webl %>% gt()


dat_text_s1_cumdegreeday_webl <- data.frame(
  label = paste("N =",samp_s1_cumdegreeday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumdegreeday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_cumdegreeday_webl <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorday", "meanmintempI",
  c("gweight","degreehours_over_30C_priorday","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_webl)

(int_tab_s2_cumdegreeday_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumdegreeday_webl <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeday_webl <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_cumdegreeday_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_cumdegreeday <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumdegreeday_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeday_webl %>% gt()


dat_text_s2_cumdegreeday_webl <- data.frame(
  label = paste("N =",samp_s2_cumdegreeday_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_cumdegreeday_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl_cumdegreeday <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorday", "meanmintempI",
  c("gweight","degreehours_over_30C_priorday","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeday)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeday)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeday)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeday)

(int_tab_abs_cumdegreeday_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_cumdegreeday <- abs_lintemp


### sample size


ss_year_webl_abs_cumdegreeday <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl_cumdegreeday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_cumdegreeday_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior week temp to predict cort


data_s1_cumdegreeweek_webl <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_webl)

(int_tab_s1_cumdegreeweek_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumdegreeweek_webl <- s1_lintemp_addmin


## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeweek_webl <- emmeans(s1_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls1_cumdegreeweek_trendmax <- emtrends(s1_lintemp_addmin,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s1_cumdegreeweek <- s1_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s1_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumdegreeweek_webl <- s1_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeweek_webl %>% gt()


dat_text_s1_cumdegreeweek_webl <- data.frame(
  label = paste("N =",samp_s1_cumdegreeweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("solid","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumdegreeweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_cumdegreeweek_webl <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_webl)

(int_tab_s2_cumdegreeweek_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumdegreeweek_webl <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeweek_webl <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(webls2_cumdegreeweek_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_webl_s2_cumdegreeweek <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_s2_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumdegreeweek_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeweek_webl %>% gt()


dat_text_s2_cumdegreeweek_webl <- data.frame(
  label = paste("N =",samp_s2_cumdegreeweek_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_cumdegreeweek_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_webl_cumdegreeweek <- prep_model_data(
  dat_growth, "WEBL", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeweek)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeweek)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeweek)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_cumdegreeweek)

(int_tab_abs_cumdegreeweek_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_cumdegreeweek <- abs_lintemp


### sample size


ss_year_webl_abs_cumdegreeweek <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_webl_abs_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl_cumdegreeweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblabs_cumdegreeweek_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### TRES

#### cort_s1


data_s1_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_tres)

(int_tab_s1_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_tres <- s1_lintemp_addmax

### Sample sizes:


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


## trends
##
(tress1trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


data <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


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
    ) +
    geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## Emmeans to check for effect of habitat


(s1byhabitat_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Minimum temperature


(t <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### cort_s2


data_s2_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_tres)

(int_tab_s2_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_tres <- s2_lintemp_addmin

### Sample sizes:


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


## trends
##
(tress2trendmax <- emtrends(s2_lintemp_addmin,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres)

(int_tab_abs_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))

abs_tres <- abs_lintemp


### Sample sizes:


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dashed","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


data <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


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
    ) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## Emmeans to check for effect of habitat


(absbyhabitat_tres <- emmeans(abs_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Maximum temperature


(tresabstrendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmaxtempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxtempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


### Minimum temperature


(t <- emtrends(abs_lintemp,specs = ~ habitat, var = c("meanmintempI_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `min temp trend` = "meanmintempI_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


data <- prep_model_data(
  dat_growth, "TRES", "meanmintempI", "meanmintempI",
  c("gweight","meanmintempI","meanmintempI","juliandate","brood_size","age"))


#### use prior day temp to predict cort instead


data_s1_priordayt_tres <- prep_model_data(
  dat_growth, "TRES", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_tres)

(int_tab_s1_priordayt_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_prior_day_tres <- s1_lintemp_addmin


## Emmeans to check for effect of habitat


(s1byhabitat_priordayt_tres <- emmeans(s1_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_priordayt_trendmax <- emtrends(s1_lintemp_addmin,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_priordayt <- s1_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_priordayt_tres <- s1_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordayt_tres %>% gt()


dat_text_s1_priordayt_tres <- data.frame(
  label = paste("N =",samp_s1_priordayt_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
    geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_priordayt_tres <- prep_model_data(
  dat_growth, "TRES", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_tres)

(int_tab_s2_priordayt_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_prior_day_tres <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_priordayt_tres <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_priordayt_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_priordayt <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_priordayt_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordayt_tres %>% gt()


dat_text_s2_priordayt_tres <- data.frame(
  label = paste("N =",samp_s2_priordayt_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dashed")) +
    geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres_priordayt <- prep_model_data(
  dat_growth, "TRES", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_tres_priordayt)

(int_tab_abs_priordayt_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_priordayt <- abs_lintemp


### sample size


ss_year_tres_abs_priordayt <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_priordayt %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_priordayt_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxt_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxt_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use prior day heat index to predict cort instead


data_s1_priordaymaxhhi_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_prior", "mint_prior",
  c("gweight","maxhi_prior","mint_prior","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordaymaxhhi_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordaymaxhhi_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordaymaxhhi_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordaymaxhhi_tres)

(int_tab_s1_priordaymaxhhi_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_priordaymaxhhi_tres <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_priordaymaxhhi_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_priordaymaxhhi_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_priordaymaxhhi <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_priordaymaxhhi_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_priordaymaxhhi_tres %>% gt()


dat_text_s1_priordaymaxhhi_tres <- data.frame(
  label = paste("N =",samp_s1_priordaymaxhhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("solid","solid","solid","solid")) +
    geom_text(data = dat_text_s1_priordaymaxhhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_priordaymaxhhi_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_prior", "mint_prior",
  c("gweight","maxhi_prior","mint_prior","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordaymaxhhi_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordaymaxhhi_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordaymaxhhi_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordaymaxhhi_tres)

(int_tab_s2_priordaymaxhhi_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_priordaymaxhhi_tres <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_priordaymaxhhi_tres <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_priordaymaxhhi_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_priordaymaxhhi <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_priordaymaxhhi_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_priordaymaxhhi_tres %>% gt()


dat_text_s2_priordaymaxhhi_tres <- data.frame(
  label = paste("N =",samp_s2_priordaymaxhhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","solid","dotted","dashed")) +
    geom_text(data = dat_text_s2_priordaymaxhhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres_priordaymaxhhi <- prep_model_data(
  dat_growth, "TRES", "maxhi_prior", "mint_prior",
  c("gweight","maxhi_prior","mint_prior","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_priordaymaxhhi)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_tres_priordaymaxhhi)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled * habitat + mint_prior_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_tres_priordaymaxhhi)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_tres_priordaymaxhhi)

(int_tab_abs_priordaymaxhhi_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_priordaymaxhhi <- abs_lintemp


### sample size


ss_year_tres_abs_priordaymaxhhi <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_priordaymaxhhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_priordaymaxhhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_priordaymaxhhi_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxhi_prior_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_prior_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


## Use prior week heat index to test against cort


data_s1_weekhi_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_weekhi_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_weekhi_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_weekhi_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_weekhi_tres)

(int_tab_s1_weekhi_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_weekhi_tres <- s1_lintemp_noint


## Emmeans to check for effect of habitat


(s1byhabitat_weekhi_tres <- emmeans(s1_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_weekhi_trendmax <- emtrends(s1_lintemp_noint,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_weekhi <- s1_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_weekhi_tres <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_weekhi_tres %>% gt()


dat_text_s1_weekhi_tres <- data.frame(
  label = paste("N =",samp_s1_weekhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_weekhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## s2


data_s2_weekhi_tres <- prep_model_data(
  dat_growth, "TRES", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_weekhi_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_weekhi_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_weekhi_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_weekhi_tres)

(int_tab_s2_weekhi_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_weekhi_tres <- s2_lintemp_addmin


## Emmeans to check for effect of habitat


(s2byhabitat_weekhi_tres <- emmeans(s2_lintemp_addmin,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_weekhi_trendmax <- emtrends(s2_lintemp_addmin,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_weekhi <- s2_lintemp_addmin@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_weekhi_tres <- s2_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_weekhi_tres %>% gt()


dat_text_s2_weekhi_tres <- data.frame(
  label = paste("N =",samp_s2_weekhi_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
    geom_text(data = dat_text_s2_weekhi_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres_weekhi <- prep_model_data(
  dat_growth, "TRES", "maxhi_week", "meanmintempI",
  c("gweight","maxhi_week","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_weekhi)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_weekhi)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_weekhi)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxhi_week_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_weekhi)

(int_tab_abs_weekhi_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_weekhi <- abs_lintemp


### sample size


ss_year_tres_abs_weekhi <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_weekhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres_weekhi, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_weekhi_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("maxhi_week_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "maxhi_week_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior day hi to predict cort


data_s1_cumhiday_tres <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorday", "meanmintempI",
  c("gweight","hihours_over_30hi_priorday","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumhiday_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiday_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiday_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiday_tres)

(int_tab_s1_cumhiday_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumhiday_tres <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_cumhiday_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_cumhiday_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_cumhiday <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumhiday_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiday_tres %>% gt()


dat_text_s1_cumhiday_tres <- data.frame(
  label = paste("N =",samp_s1_cumhiday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumhiday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_cumhiday_tres <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorday", "meanmintempI",
  c("gweight","hihours_over_30hi_priorday","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumhiday_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiday_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiday_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiday_tres)

(int_tab_s2_cumhiday_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumhiday_tres <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_cumhiday_tres <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_cumhiday_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_cumhiday <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumhiday_tres <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiday_tres %>% gt()


dat_text_s2_cumhiday_tres <- data.frame(
  label = paste("N =",samp_s2_cumhiday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_cumhiday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres_cumhiday <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorday", "meanmintempI",
  c("gweight","hihours_over_30hi_priorday","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumhiday)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumhiday)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumhiday)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumhiday)

(int_tab_abs_cumhiday_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_cumhiday <- abs_lintemp_noint


### sample size


ss_year_tres_abs_cumhiday <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumhiday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dashed","Grassland" = "dashed","Row crop" = "dashed")) +
    geom_text(data = dat_text_tres_cumhiday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_cumhiday_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior week hi to predict cort


data_s1_cumhiweek_tres <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_cumhiweek_tres)

(int_tab_s1_cumhiweek_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumhiweek_tres <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_cumhiweek_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_cumhiweek_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_cumhiweek <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumhiweek_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumhiweek_tres %>% gt()


dat_text_s1_cumhiweek_tres <- data.frame(
  label = paste("N =",samp_s1_cumhiweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumhiweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### s2


data_s2_cumhiweek_tres <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_cumhiweek_tres)

(int_tab_s2_cumhiweek_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumhiweek_tres <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_cumhiweek_tres <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_cumhiweek_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_cumhiweek <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumhiweek_tres <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumhiweek_tres %>% gt()


dat_text_s2_cumhiweek_tres <- data.frame(
  label = paste("N =",samp_s2_cumhiweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_cumhiweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_abs_tres_cumhiweek <- prep_model_data(
  dat_growth, "TRES", "hihours_over_30hi_priorweek", "meanmintempI",
  c("gweight","hihours_over_30hi_priorweek","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_tres_cumhiweek)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_tres_cumhiweek)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_tres_cumhiweek)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ hihours_over_30hi_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_abs_tres_cumhiweek)

(int_tab_abs_cumhiweek_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_cumhiweek <- abs_lintemp


### sample size


ss_year_tres_abs_cumhiweek <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumhiweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_abs_text_tres_cumhiweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_cumhiweek_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("hihours_over_30hi_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihours_over_30hi_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior day temp to predict cort


data_s1_cumdegreeday_tres <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorday", "meanmintempI",
  c("gweight","degreehours_over_30C_priorday","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeday_tres)

(int_tab_s1_cumdegreeday_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumdegreeday_tres <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeday_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_cumdegreeday_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_cumdegreeday <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumdegreeday_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeday_tres %>% gt()


dat_text_s1_cumdegreeday_tres <- data.frame(
  label = paste("N =",samp_s1_cumdegreeday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumdegreeday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### use cumulative prior day temp to predict cort


data_s2_cumdegreeday_tres <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorday", "meanmintempI",
  c("gweight","degreehours_over_30C_priorday","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeday_tres)

(int_tab_s2_cumdegreeday_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumdegreeday_tres <- s2_lintemp


## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeday_tres <- emmeans(s2_lintemp,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_cumdegreeday_trendmax <- emtrends(s2_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_cumdegreeday <- s2_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumdegreeday_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeday_tres %>% gt()


dat_text_s2_cumdegreeday_tres <- data.frame(
  label = paste("N =",samp_s2_cumdegreeday_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
    geom_text(data = dat_text_s2_cumdegreeday_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres_cumdegreeday <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorday", "meanmintempI",
  c("gweight","degreehours_over_30C_priorday","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeday)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeday)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeday)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorday_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeday)

(int_tab_abs_cumdegreeday_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_cumdegreeday <- abs_lintemp


### sample size


ss_year_tres_abs_cumdegreeday <- abs_lintemp@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumdegreeday %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres_cumdegreeday, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_cumdegreeday_trendmax <- emtrends(abs_lintemp,specs = ~ habitat, var = c("degreehours_over_30C_priorday_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorday_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


#### use cumulative prior week temp to predict cort


data_s1_cumdegreeweek_tres <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_cumdegreeweek_tres)

(int_tab_s1_cumdegreeweek_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_cumdegreeweek_tres <- s1_lintemp_addmax


## Emmeans to check for effect of habitat


(s1byhabitat_cumdegreeweek_tres <- emmeans(s1_lintemp_addmax,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress1_cumdegreeweek_trendmax <- emtrends(s1_lintemp_addmax,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s1_cumdegreeweek <- s1_lintemp_addmax@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s1_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s1_cumdegreeweek_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s1_cumdegreeweek_tres %>% gt()


dat_text_s1_cumdegreeweek_tres <- data.frame(
  label = paste("N =",samp_s1_cumdegreeweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_cumdegreeweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### use cumulative prior week temp to predict cort


data_s2_cumdegreeweek_tres <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_cumdegreeweek_tres)

(int_tab_s2_cumdegreeweek_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_cumdegreeweek_tres <- s2_lintemp_noint


## Emmeans to check for effect of habitat


(s2byhabitat_cumdegreeweek_tres <- emmeans(s2_lintemp_noint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:t.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


## trends
##
(tress2_cumdegreeweek_trendmax <- emtrends(s2_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


### Sample sizes:


ss_year_tres_s2_cumdegreeweek <- s2_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_s2_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

samp_s2_cumdegreeweek_tres <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())
samp_s2_cumdegreeweek_tres %>% gt()


dat_text_s2_cumdegreeweek_tres <- data.frame(
  label = paste("N =",samp_s2_cumdegreeweek_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_cumdegreeweek_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


#### abs_change_cort


data_tres_cumdegreeweek <- prep_model_data(
  dat_growth, "TRES", "degreehours_over_30C_priorweek", "meanmintempI",
  c("gweight","degreehours_over_30C_priorweek","meanmintempI","juliandate","brood_size","age"))

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeweek)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeweek)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled * habitat + meanmintempI_scaled + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeweek)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ degreehours_over_30C_priorweek_scaled + meanmintempI_scaled + habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_cumdegreeweek)

(int_tab_abs_cumdegreeweek_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_cumdegreeweek <- abs_lintemp_noint


### sample size


ss_year_tres_abs_cumdegreeweek <- abs_lintemp_noint@frame %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2022` + `2023`)
ss_year_tres_abs_cumdegreeweek %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres_cumdegreeweek, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresabs_cumdegreeweek_trendmax <- emtrends(abs_lintemp_noint,specs = ~ habitat, var = c("degreehours_over_30C_priorweek_scaled")) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "degreehours_over_30C_priorweek_scaled.trend",Df = "df", `T-ratio` = "t.ratio", P = "p.value") %>%
    gt())


save(list = ls(), file = "data/models_cort.RData")
rm(list = ls()); gc()
dat_surv <- read_rds("data/survival_attempt.rds") %>%
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
anova_int_tab <- function(m_full, m_addmin, m_addmax, m_noint, digits = 4) {
  if (inherits(m_full, "glm") && !inherits(m_full, "lmerMod")) {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint, test = "Chisq") %>% tibble() %>%
        rename(Chisq = Deviance) %>%
        mutate(AIC = c(AIC(m_full), AIC(m_single), AIC(m_noint)), .before = Df) %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = Df)
    }
    p_col <- "Pr(>Chi)"; rnd_cols <- "Chisq"
  } else {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint) %>% tibble() %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = npar)
    }
    p_col <- "Pr(>Chisq)"; rnd_cols <- c("AIC","Chisq")
  }
  bind_rows(mk_row(m_addmin), mk_row(m_addmax)) %>% as_tibble() %>%
    mutate(across(where(is.numeric), ~round(.x, digits = digits)), P = .data[[p_col]]) %>%
    dplyr::select(Model, AIC, Chisq, P) %>%
    mutate(max_or_min = c(rep("Max temp", times = 3), rep("Min temp", times = 3)),
           across(all_of(rnd_cols), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001, "<0.001", as.character(P))) %>%
    group_by(max_or_min) %>% gt()
}


## ================================================================
## SECTION 3: SURVIVAL MODELS
## ================================================================


#### Nest period model


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
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
                            data = dplyr::filter(dat_surv,
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
                            data = dplyr::filter(dat_surv,
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
                           data = dplyr::filter(dat_surv,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_webl <- anova_int_tab(s_nestpd_WEBL, s_nestpd_WEBL_addmin, s_nestpd_WEBL_addmax, s_nestpd_WEBL_noint))


(survivalbyhabitat_summary_webl <- summary(s_nestpd_WEBL_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
   mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
          across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax <- emtrends(s_nestpd_WEBL_addmin,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
   mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
          df = round(df),
          p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
   rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
   gt())


data = s_nestpd_WEBL_addmin$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl <- emmeans(s_nestpd_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
   mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
          across(p.value,~round(.x,digits = 3))) %>% gt())


data_webl = s_nestpd_WEBL_addmin$data

mean_temp_webl <- mean(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)

temp_trans_webl <- trans_new("temp_trans_webl",
                          transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                          inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())

ss_year_survival_webl <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))

dat_text_webl <- data.frame(
  label = paste("N =",samp_webl$count),
  group   = factor(samp_webl$habitat)
)

(fig3_webl <- predict_response(s_nestpd_WEBL_addmin,terms = c("meanmaxt_nestpd_scaled [all]","habitat")) %>%
   plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
   facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max temp over preceding week (°C)") +
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
                                  (40-mean_temp_webl)/sd_temp_webl
                                  )) +
   scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dashed","Grassland" = "dotted","Row crop" = "dashed")) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


### TRES

#### Nest period model


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
               family = binomial(link = "logit"),
               data = dplyr::filter(dat_surv,
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
               data = dplyr::filter(dat_surv,
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
               data = dplyr::filter(dat_surv,
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
               data = dplyr::filter(dat_surv,
                                    Species == "TRES",
                                    !is.na(nest_fledged),
                                    !is.na(clutch_size)) %>%
                 mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                               ~ scale(.x)[,1],
                               .names = "{.col}_scaled")
                      )
               )


(int_tab_survival_nestpd_tres <- anova_int_tab(s_nestpd_TRES, s_nestpd_TRES_addmin, s_nestpd_TRES_addmax, s_nestpd_TRES_noint))


(tressurvival_nestpd_trendmax <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("meanmaxt_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxt_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


(survival_nestpd_byhabitat_tres <- emmeans(s_nestpd_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
   mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
          across(p.value,~round(.x,digits = 3))) %>% gt())


samp_tres <- s_nestpd_TRES_noint$data %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_tres <- s_nestpd_TRES_noint$data %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)
ss_year_survival_tres %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_tres <- data.frame(
label = paste("N =",samp_tres$count),
group = factor(samp_tres$habitat)
)

prop_trans_tres <- trans_new("prop_trans_tres",
                             transform = function(x){(x * 100)},
                             inverse = function(x){x})

dat_tres <- predict_response(s_nestpd_TRES_noint,terms = c("habitat"))

(fig3_tres <- ggplot(data = dat_tres, mapping = aes(x = x, color = x, y = predicted,ymin = conf.low,ymax = conf.high)) +
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
   #                     breaks = c((20-mean_temp)/sd_temp,
   #                                ),
   #                     # breaks = c(20,30,40,50),
   #                     # labels = c("20","30","40","50"),
   #                     # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
   #                     # limits = c(18,5)
   #                     ) +
   geom_text(data = dat_text_tres, mapping = aes(x = c(1,2,3,4), y = -Inf,label = label),hjust = .5,vjust = -.2,inherit.aes = FALSE) +
    theme(legend.position = "none") +
    annotate(geom = "text", x = c(1,2,3,4),y = Inf,hjust = .5,vjust = 1,label = c("a","ab","b","b"))
 )


## combined survival results for manuscript


library(egg)
ggplot_build(fig3_webl)$layout$panel_scales_y
(p_full <- ggarrange(fig3_webl + theme(text = element_text(size = 12)),fig3_tres + theme(text = element_text(size = 12),axis.title.y = element_blank()),ncol = 2,
          labels = c("(a): Western Bluebird","(b): Tree Swallow")))


### What about trying to estimate the effect of provis on survival? Will need to return to DAG to know what to put in the model.


## Other temperature measures
## WEBL

### meanmaxhi_nestpd


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled * habitat + meanminhi_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
                                          Species == "WEBL",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_WEBL_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled + meanminhi_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_WEBL_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled * habitat + meanminhi_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_WEBL_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled + meanminhi_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_webl_meanmaxhi <- anova_int_tab(s_nestpd_WEBL, s_nestpd_WEBL_addmin, s_nestpd_WEBL_addmax, s_nestpd_WEBL_noint))


(survivalbyhabitat_summary_webl_meanmaxhi <- summary(s_nestpd_WEBL_addmax) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax_meanmaxhi <- emtrends(s_nestpd_WEBL_addmax,specs = pairwise ~ habitat, var = c("meanmaxhi_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxhi_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


data = s_nestpd_WEBL_addmax$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl_meanmaxhi <- emmeans(s_nestpd_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


data_webl = s_nestpd_WEBL_addmax$data


mean_temp_webl <- mean(data_webl %>% pull(meanmaxhi_nestpd),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(meanmaxhi_nestpd),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_webl_meanmaxhi <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl_meanmaxhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_webl <- data.frame(
  label = paste("N =",samp_webl$count),
  group   = factor(samp_webl$habitat)
)

(fig3_webl_meanmaxhi <- predict_response(s_nestpd_WEBL_addmax,terms = c("meanmaxhi_nestpd_scaled [all]","habitat")) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max heat index over preceding week (\u00b0C)") +
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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


### deghr_30


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
                                          Species == "WEBL",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_WEBL_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_WEBL_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_WEBL_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_webl_deghr_30 <- anova_int_tab(s_nestpd_WEBL, s_nestpd_WEBL_addmin, s_nestpd_WEBL_addmax, s_nestpd_WEBL_noint))


(survivalbyhabitat_summary_webl_deghr_30 <- summary(s_nestpd_WEBL_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax_deghr_30 <- emtrends(s_nestpd_WEBL_addmin,specs = pairwise ~ habitat, var = c("deghr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "deghr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


data = s_nestpd_WEBL_addmin$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl_deghr_30 <- emmeans(s_nestpd_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


data_webl = s_nestpd_WEBL_addmin$data


mean_temp_webl <- mean(data_webl %>% pull(deghr_30),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(deghr_30),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_webl_deghr_30 <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl_deghr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_webl <- data.frame(
  label = paste("N =",samp_webl$count),
  group   = factor(samp_webl$habitat)
)

(fig3_webl_deghr_30 <- predict_response(s_nestpd_WEBL_addmin,terms = c("deghr_30_scaled [all]","habitat")) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C") +
    ylab("Predicted percent of eggs fledging") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_webl,
                       breaks = c((0-mean_temp_webl)/sd_temp_webl,
                                  (2000-mean_temp_webl)/sd_temp_webl,
                                  (4000-mean_temp_webl)/sd_temp_webl,
                                  (6000-mean_temp_webl)/sd_temp_webl,
                                  (8000-mean_temp_webl)/sd_temp_webl#,
                                  #(10000-mean_temp_webl)/sd_temp_webl#,
                       ),
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


### hihr_30


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
                                          Species == "WEBL",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_WEBL_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_WEBL_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "WEBL",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_WEBL_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_webl_hihr_30 <- anova_int_tab(s_nestpd_WEBL, s_nestpd_WEBL_addmin, s_nestpd_WEBL_addmax, s_nestpd_WEBL_noint))


(survivalbyhabitat_summary_webl_hihr_30 <- summary(s_nestpd_WEBL_addmax) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax_hihr_30 <- emtrends(s_nestpd_WEBL_addmax,specs = pairwise ~ habitat, var = c("hihr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


data = s_nestpd_WEBL_addmax$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl_hihr_30 <- emmeans(s_nestpd_WEBL_addmax,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


data_webl = s_nestpd_WEBL_addmax$data


mean_temp_webl <- mean(data_webl %>% pull(hihr_30),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(hihr_30),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_webl_hihr_30 <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl_hihr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_webl <- data.frame(
  label = paste("N =",samp_webl$count),
  group   = factor(samp_webl$habitat)
)

(fig3_webl_hihr_30 <- predict_response(s_nestpd_WEBL_addmax,terms = c("hihr_30_scaled [all]","habitat")) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index-hours >25\u00b0C") +
    ylab("Predicted percent of eggs fledging") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_webl,
                       breaks = c((0-mean_temp_webl)/sd_temp_webl,
                                  (3000-mean_temp_webl)/sd_temp_webl,
                                  (6000-mean_temp_webl)/sd_temp_webl,
                                  (9000-mean_temp_webl)/sd_temp_webl,
                                  (12000-mean_temp_webl)/sd_temp_webl,
                                  (15000-mean_temp_webl)/sd_temp_webl#,
                       ),
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## TRES

### meanmaxhi_nestpd


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled * habitat + meanminhi_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
                                          Species == "TRES",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_TRES_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled + meanminhi_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "TRES",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_TRES_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled * habitat + meanminhi_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "TRES",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_TRES_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled + meanminhi_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_tres_meanmaxhi <- anova_int_tab(s_nestpd_TRES, s_nestpd_TRES_addmin, s_nestpd_TRES_addmax, s_nestpd_TRES_noint))


(survivalbyhabitat_summary_tres_meanmaxhi <- summary(s_nestpd_TRES_addmax) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(tressurvival_nestpd_trendmax_meanmaxhi <- emtrends(s_nestpd_TRES_addmax,specs = pairwise ~ habitat, var = c("meanmaxhi_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxhi_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


data = s_nestpd_TRES_addmax$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_tres_meanmaxhi <- emmeans(s_nestpd_TRES_addmax,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


data_tres = s_nestpd_TRES_addmax$data


mean_temp_tres <- mean(data_tres %>% pull(meanmaxhi_nestpd),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(meanmaxhi_nestpd),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_tres_meanmaxhi <- data_tres %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_tres_meanmaxhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_tres <- data.frame(
  label = paste("N =",samp_tres$count),
  group   = factor(samp_tres$habitat)
)

(fig3_tres_meanmaxhi <- predict_response(s_nestpd_TRES_addmax,terms = c("meanmaxhi_nestpd_scaled [all]","habitat")) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean daily max heat index over preceding week (\u00b0C)") +
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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


### deghr_30


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
                                          Species == "TRES",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_TRES_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "TRES",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_TRES_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "TRES",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_TRES_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_tres_deghr_30 <- anova_int_tab(s_nestpd_TRES, s_nestpd_TRES_addmin, s_nestpd_TRES_addmax, s_nestpd_TRES_noint))


(survivalbyhabitat_summary_tres_deghr_30 <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(tressurvival_nestpd_trendmax_deghr_30 <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("deghr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "deghr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


data = s_nestpd_TRES_noint$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_tres_deghr_30 <- emmeans(s_nestpd_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


data_tres = s_nestpd_TRES_noint$data


mean_temp_tres <- mean(data_tres %>% pull(deghr_30),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(deghr_30),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_tres_deghr_30 <- data_tres %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_tres_deghr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_tres <- data.frame(
  label = paste("N =",samp_tres$count),
  group   = factor(samp_tres$habitat)
)

(fig3_tres_deghr_30 <- predict_response(s_nestpd_TRES_noint,terms = c("deghr_30_scaled [all]","habitat")) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative degree-hours >30\u00b0C") +
    ylab("Predicted percent of eggs fledging") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_tres,
                       breaks = c((0-mean_temp_tres)/sd_temp_tres,
                                  (2000-mean_temp_tres)/sd_temp_tres,
                                  (4000-mean_temp_tres)/sd_temp_tres,
                                  (6000-mean_temp_tres)/sd_temp_tres,
                                  (8000-mean_temp_tres)/sd_temp_tres#,
                                  #(10000-mean_temp_tres)/sd_temp_tres#,
                       ),
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "solid")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


### hihr_30


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(dat_surv,
                                          Species == "TRES",
                                          !is.na(nest_fledged),
                                          !is.na(clutch_size)) %>%
                       mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                     ~ scale(.x)[,1],
                                     .names = "{.col}_scaled")
                       )
)

s_nestpd_TRES_addmax <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "TRES",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)

s_nestpd_TRES_addmin <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
                                                 Species == "TRES",
                                                 !is.na(nest_fledged),
                                                 !is.na(clutch_size)) %>%
                              mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                            ~ scale(.x)[,1],
                                            .names = "{.col}_scaled")
                              )
)
s_nestpd_TRES_noint <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_tres_hihr_30 <- anova_int_tab(s_nestpd_TRES, s_nestpd_TRES_addmin, s_nestpd_TRES_addmax, s_nestpd_TRES_noint))


(survivalbyhabitat_summary_tres_hihr_30 <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


#### Emtrends to calculate effect of max temp in each habitat


(tressurvival_nestpd_trendmax_hihr_30 <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("hihr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


data = s_nestpd_TRES_noint$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_tres_hihr_30 <- emmeans(s_nestpd_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


data_tres = s_nestpd_TRES_noint$data


mean_temp_tres <- mean(data_tres %>% pull(hihr_30),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(hihr_30),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
ss_year_survival_tres_hihr_30 <- data_tres %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_tres_hihr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.))


dat_text_tres <- data.frame(
  label = paste("N =",samp_tres$count),
  group   = factor(samp_tres$habitat)
)

(fig3_tres_hihr_30 <- predict_response(s_nestpd_TRES_noint,terms = c("hihr_30_scaled [all]","habitat")) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Cumulative heat index-hours >25\u00b0C") +
    ylab("Predicted percent of eggs fledging") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title = element_blank()) +
    scale_x_continuous(trans = temp_trans_tres,
                       breaks = c((0-mean_temp_tres)/sd_temp_tres,
                                  (3000-mean_temp_tres)/sd_temp_tres,
                                  (6000-mean_temp_tres)/sd_temp_tres,
                                  (9000-mean_temp_tres)/sd_temp_tres,
                                  (12000-mean_temp_tres)/sd_temp_tres,
                                  (15000-mean_temp_tres)/sd_temp_tres#,
                       ),
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


save(list = ls(), file = "data/models_survival.RData")
rm(list = ls()); gc()
dat_provis <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

## ================================================================
## SECTION 4: PROVISIONING MODELS
## ================================================================


dat_provis %>% dplyr::select(c(species, attempt_id)) %>% distinct() %>% summarize(.by = species,n = n())


## Sample size:


## Is interaction significant?

### WEBL


m <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_linint <- glmmTMB(
  valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_noint <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(dat_provis,!is.na(mean_temp),
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


provis_webl <- m_linint

### Sample size


ss_year_webl_provis <- m_linint$frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(`2021` = replace_na(`2021`,list(count_ind = 0,count_attempt = 0))) %>%
  mutate(Total = `2021` + `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2021` = paste0(`2021`$count_ind,", ",`2021`$count_attempt),
         `2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_webl_provis <- ss_year_webl_provis %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    })
)


t_ss_year_webl_provis

samp <- m_linint$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

## num nestboxes
m_linint$frame %>% pull(attempt_id) %>% unique() %>% length()


data_webl = dplyr::filter(dat_provis,!is.na(mean_temp),
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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblprovistrend <- emtrends(m_linint,specs = ~ habitat, var = c("mean_temp_scaled"),max.degree = 1) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    #group_by(degree) %>%
    gt())


data = data = dplyr::filter(dat_provis,!is.na(mean_temp),
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


## Emmeans to check for effect of habitat


(provisbyhabitat_webl <- emmeans(m_linint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(provisbyhabitat_summary_webl <- summary(m_noint) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


### TRES


m <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_linint <- glmmTMB(
  valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_noint <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(dat_provis,!is.na(mean_temp),
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


provis_tres <- m
### Sample size


ss_year_tres_provis <- m$frame %>%
  group_by(habitat,year_fct) %>%
  summarize(count_ind = n(),count_attempt = n_distinct(attempt_id)) %>%
  rowwise() %>%
  mutate(count = tibble(count_ind,count_attempt)) %>%
  dplyr::select(-c(count_ind, count_attempt)) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>%
  as.tibble() %>% rename(Habitat = 'habitat') %>%
  ungroup() %>%
  mutate(`2021` = replace_na(`2021`,list(count_ind = 0,count_attempt = 0)),
         `2022` = replace_na(`2022`,list(count_ind = 0,count_attempt = 0))) %>%
  mutate(Total = `2021` + `2022` + `2023`) %>%
  rowwise() %>%
  mutate(`2021` = paste0(`2021`$count_ind,", ",`2021`$count_attempt),
         `2022` = paste0(`2022`$count_ind,", ",`2022`$count_attempt),
         `2023` = paste0(`2023`$count_ind,", ",`2023`$count_attempt),
         Total = paste0(Total$count_ind,", ",Total$count_attempt)
  )

(t_ss_year_tres_provis <- ss_year_tres_provis %>% gt() %>%
    grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ {
      t <- str_split(.,", ",simplify = TRUE)
      t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    })
)


t_ss_year_tres_provis

samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

## num nestboxes
m$frame %>% pull(attempt_id) %>% unique() %>% length()


data_tres = dplyr::filter(dat_provis,!is.na(mean_temp),
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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    ylim(0,50) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


## combined figure for fig5


library(egg)
ggplot_build(fig5_webl)$layout$panel_scales_y
ggplot_build(fig5_tres)$layout$panel_scales_y
(p_full <- ggarrange(fig5_webl + ylim(0,40) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                    axis.title.x = element_blank()),fig5_tres + labs(title = element_blank()) + ylim(0,40) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                                                                                                   axis.text.y = element_blank(),
                                                                                                                                                                                   axis.title.y = element_blank(),
                                                                                                                                                                                   axis.title.x = element_text(hjust = -.8)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))


## Emmeans to check for effect of habitat


(provisbyhabitat_tres <- emmeans(m,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(provisbyhabitat_summary_tres <- summary(m_noint) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


(tresprovistrend <- emtrends(m,specs = ~ degree | habitat, var = "mean_temp_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


data = dplyr::filter(dat_provis,!is.na(mean_temp),
                     !is.na(julian_date),
                     !is.na(mean_nestling_age),
                     !is.na(tod_h),
                     !is.na(site),
                     !is.na(attempt_id),
                     !is.na(year),
                     valid_detections < 30,
                     species == "TRES"# ,
                     # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                     # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                     # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
) %>%
  mutate(across(c(tod_h,mean_temp,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


## Other temperature measures

library(future)
instant_temp <- read_rds("data/provis_manytempmeasures.rds")

## Is interaction significant?

### WEBL
#### maxhi of closest 30 min period to start time of hour


m <- glmmTMB(
  valid_detections ~ poly(hi_30min_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(instant_temp,!is.na(hi_30min),
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
               mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                             ~ scale(.x)[,1],
                             .names = "{.col}_scaled")))

emmeans(m,specs = pairwise ~ habitat,by = c("hi_30min_scaled"), at = list(hi_30min_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
emmip(m,formula = habitat ~ hi_30min_scaled, at = list(hi_30min_scaled = seq(from = -2.5, to = 2.5, by = .1)))

m_linint <- glmmTMB(
  valid_detections ~ hi_30min_scaled * habitat + poly(hi_30min_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(instant_temp,!is.na(hi_30min),
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
                      mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                                    ~ scale(.x)[,1],
                                    .names = "{.col}_scaled")))

m_noint <- glmmTMB(
  valid_detections ~ poly(hi_30min_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(instant_temp,!is.na(hi_30min),
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
                     mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                                   ~ scale(.x)[,1],
                                   .names = "{.col}_scaled")))
c <- anova(m,m_linint,m_noint) %>% tibble() %>% mutate(Model = c("No temp * LC interaction","linear temp * LC interaction", "quadratic temp * LC  interaction"),.before = Df)

(int_tab_provis_webl_hi <- c %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           instant_temp = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,instant_temp) %>%
    mutate(across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           instant_temp = if_else(instant_temp < 0.001,"<0.001",as.character(instant_temp))) %>% gt())


provis_webl_hi <- m

### Sample size


ss_year_webl_provis_hi <- m$frame %>%
  group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
  relocate(`2021`,.before = `2022`) %>%
  ungroup() %>%
  mutate(across(c(`2021`,`2022`,`2023`),~ replace_na(.x,0)),
         Total = `2021`+`2022`+`2023`)
ss_year_webl_provis_hi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE))

samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_webl_hi <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

## num nestboxes
m$frame %>% pull(attempt_id) %>% unique() %>% length()


data_webl = dplyr::filter(instant_temp,!is.na(hi_30min),
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
  mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


mean_hi_30min_webl <- mean(data_webl %>% pull(hi_30min))
sd_hi_30min_webl <- sd(data_webl %>% pull(hi_30min))


hi_30min_trans_webl <- trans_new("hi_30min_trans_webl",
                             transform = function(x){(x * sd_hi_30min_webl) + mean_hi_30min_webl},
                             inverse = function(x){x})

(fig5_webl_hi <- predict_response(m,terms = c("hi_30min_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean heat index (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title =  "WEBL provis: max * LC interaction") +
    scale_x_continuous(trans = hi_30min_trans_webl,
                       breaks = c((20-mean_hi_30min_webl)/sd_hi_30min_webl,
                                  (25-mean_hi_30min_webl)/sd_hi_30min_webl,
                                  (30-mean_hi_30min_webl)/sd_hi_30min_webl,
                                  (35-mean_hi_30min_webl)/sd_hi_30min_webl,
                                  (40-mean_hi_30min_webl)/sd_hi_30min_webl),
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl_hi, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(weblprovistrend_hi <- emtrends(m,specs = ~ degree | habitat, var = "hi_30min_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hi_30min_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


data = dplyr::filter(instant_temp,!is.na(hi_30min),
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
  mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


## Emmeans to check for effect of habitat


(provisbyhabitat_webl_hi <- emmeans(m,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(provisbyhabitat_summary_webl_hi <- summary(m) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(hi_30min_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())

### TRES
#### maxhi of closest 30 min period to start time of hour


m <- glmmTMB(
  valid_detections ~ poly(hi_30min_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(instant_temp,!is.na(hi_30min),
                                  !is.na(julian_date),
                                  !is.na(mean_nestling_age),
                                  !is.na(tod_h),
                                  !is.na(site),
                                  !is.na(attempt_id),
                                  !is.na(year),
                                  species == "TRES"# ,
                                  # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                  # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                  # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
             ) %>%
               mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                             ~ scale(.x)[,1],
                             .names = "{.col}_scaled")))


m_linint <- glmmTMB(
  valid_detections ~ hi_30min_scaled * habitat + poly(hi_30min_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(instant_temp,!is.na(hi_30min),
                                         !is.na(julian_date),
                                         !is.na(mean_nestling_age),
                                         !is.na(tod_h),
                                         !is.na(site),
                                         !is.na(attempt_id),
                                         !is.na(year),
                                         species == "TRES"# ,
                                         # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                         # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                         # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                    ) %>%
                      mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                                    ~ scale(.x)[,1],
                                    .names = "{.col}_scaled")))

m_noint <- glmmTMB(
  valid_detections ~ poly(hi_30min_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(instant_temp,!is.na(hi_30min),
                                        !is.na(julian_date),
                                        !is.na(mean_nestling_age),
                                        !is.na(tod_h),
                                        !is.na(site),
                                        !is.na(attempt_id),
                                        !is.na(year),
                                        species == "TRES"# ,
                                        # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                                        # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                                        # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
                   ) %>%
                     mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                                   ~ scale(.x)[,1],
                                   .names = "{.col}_scaled")))
c <- anova(m,m_linint,m_noint) %>% tibble() %>% mutate(Model = c("No temp * LC interaction","linear temp * LC interaction", "quadratic temp * LC  interaction"),.before = Df)

(int_tab_provis_tres_hi <- c %>%
    as.tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           instant_temp = `Pr(>Chisq)`) %>%
    dplyr::select(Model,AIC,Chisq,instant_temp) %>%
    mutate(across(c(AIC,Chisq), ~ round(.x, digits = 2)),
           instant_temp = if_else(instant_temp < 0.001,"<0.001",as.character(instant_temp))) %>% gt())


provis_tres_hi <- m

### Sample size


ss_year_tres_provis_hi <- m$frame %>%
  group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
  relocate(`2021`,.before = `2022`) %>%
  ungroup() %>%
  mutate(across(c(`2021`,`2022`,`2023`),~ replace_na(.x,0)),
         Total = `2021`+`2022`+`2023`)
ss_year_tres_provis_hi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE))

samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt()


dat_text_tres_hi <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

## num nestboxes
m$frame %>% pull(attempt_id) %>% unique() %>% length()


data_tres = dplyr::filter(instant_temp,!is.na(hi_30min),
                          !is.na(julian_date),
                          !is.na(mean_nestling_age),
                          !is.na(tod_h),
                          !is.na(site),
                          !is.na(attempt_id),
                          !is.na(year),
                          species == "TRES"# ,
                          # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                          # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                          # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
) %>%
  mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


mean_hi_30min_tres <- mean(data_tres %>% pull(hi_30min))
sd_hi_30min_tres <- sd(data_tres %>% pull(hi_30min))


hi_30min_trans_tres <- trans_new("hi_30min_trans_tres",
                             transform = function(x){(x * sd_hi_30min_tres) + mean_hi_30min_tres},
                             inverse = function(x){x})

(fig5_tres_hi <- predict_response(m,terms = c("hi_30min_scaled [all]","habitat"),bias_correction = TRUE) %>%
    plot(line_size = 1.5,alpha = .2,show_data = TRUE,limit_range = TRUE) +
    aes(linetype = .data[["group"]]) +
    theme_classic() +
    facet_wrap(~ group, ncol = 2) +
    xlab("Mean heat index (\u00b0C)") +
    ylab("Provisioning (count/hour)") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16)) +
    labs(title =  "TRES provis: max * LC interaction") +
    scale_x_continuous(trans = temp_trans_tres,
                       breaks = c((20-mean_hi_30min_tres)/sd_hi_30min_tres,
                                  (25-mean_hi_30min_tres)/sd_hi_30min_tres,
                                  (30-mean_hi_30min_tres)/sd_hi_30min_tres,
                                  (35-mean_hi_30min_tres)/sd_hi_30min_tres,
                                  (40-mean_hi_30min_tres)/sd_hi_30min_tres),
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    ylim(0,60) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


(tresprovistrend_hi <- emtrends(m,specs = ~ degree | habitat, var = "hi_30min_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hi_30min_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


data = data = dplyr::filter(instant_temp,!is.na(hi_30min),
                            !is.na(julian_date),
                            !is.na(mean_nestling_age),
                            !is.na(tod_h),
                            !is.na(site),
                            !is.na(attempt_id),
                            !is.na(year),
                            species == "TRES"# ,
                            # tod >= hms::as.hms('11:00:00', tz = "America/Los_Angeles"),
                            # tod_end <= hms::as.hms('17:00:00',tz = "America/Los_Angeles"),
                            # tod_end >= hms::as.hms('11:00:00', tz = "America/Los_Angeles")
) %>%
  mutate(across(c(tod_h,hi_30min,julian_date,mean_nestling_age),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))


## Emmeans to check for effect of habitat


(provisbyhabitat_tres_hi <- emmeans(m,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())


### Check for effect of temperature


(provisbyhabitat_summary_tres_hi <- summary(m) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(hi_30min_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())


## combined figure for fig5


library(egg)
ggplot_build(fig5_webl_hi)$layout$panel_scales_y
ggplot_build(fig5_tres_hi)$layout$panel_scales_y
(p_full <- ggarrange(fig5_webl_hi + ylim(0,55.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                    axis.title.x = element_blank()),fig5_tres_hi + labs(title = element_blank()) + ylim(0,55.9) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                                                                                                   axis.text.y = element_blank(),
                                                                                                                                                                                   axis.title.y = element_blank(),
                                                                                                                                                                                   axis.title.x = element_text(hjust = -.8)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))


save(list = ls(), file = "data/models_provis.RData")
rm(list = ls()); gc()
dat_provis <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))
dat_growth <- read_rds("data/growth_cort_provis_manytempmeasures.rds") %>%
  mutate(year_fct = as.factor(year))
dat_surv <- read_rds("data/survival_attempt.rds") %>%
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
scale_cols <- function(data, cols) {
  data %>% mutate(across(all_of(cols), ~ scale(.x)[,1], .names = "{.col}_scaled"))
}
prep_model_data <- function(dat, species, temp_var, min_temp_var, cols_to_scale,
                             add_sq = TRUE) {
  d <- dat %>%
    dplyr::filter(Species == species, !is.na(.data[[temp_var]]), !is.na(.data[[min_temp_var]])) %>%
    scale_cols(cols_to_scale)
  if (add_sq) { sq_col <- paste0(temp_var, "_scaled_sq"); d[[sq_col]] <- d[[paste0(temp_var, "_scaled")]]^2 }
  d
}
anova_int_tab <- function(m_full, m_addmin, m_addmax, m_noint, digits = 4) {
  if (inherits(m_full, "glm") && !inherits(m_full, "lmerMod")) {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint, test = "Chisq") %>% tibble() %>%
        rename(Chisq = Deviance) %>%
        mutate(AIC = c(AIC(m_full), AIC(m_single), AIC(m_noint)), .before = Df) %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = Df)
    }
    p_col <- "Pr(>Chi)"; rnd_cols <- "Chisq"
  } else {
    mk_row <- function(m_single) {
      anova(m_full, m_single, m_noint) %>% tibble() %>%
        mutate(Model = c("no interaction","single interaction","both interacting"), .before = npar)
    }
    p_col <- "Pr(>Chisq)"; rnd_cols <- c("AIC","Chisq")
  }
  bind_rows(mk_row(m_addmin), mk_row(m_addmax)) %>% as_tibble() %>%
    mutate(across(where(is.numeric), ~round(.x, digits = digits)), P = .data[[p_col]]) %>%
    dplyr::select(Model, AIC, Chisq, P) %>%
    mutate(max_or_min = c(rep("Max temp", times = 3), rep("Min temp", times = 3)),
           across(all_of(rnd_cols), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001, "<0.001", as.character(P))) %>%
    group_by(max_or_min) %>% gt()
}


## ================================================================
## SECTION 5: SEASONAL SENSITIVITY
## Does adding temp * julian date interaction improve model fit?
## ================================================================


data_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl)

g_lintemp <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

g_lintemp_addmax <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

g_lintemp_addmin <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

g_lintemp_noint <- lmerTest::lmer(
  gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_webl)

(int_tab_growth_webl <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))


anova(g_lintemp,prior_model)

g_lintemp_webl <- g_lintemp


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


samp <- g_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
                       ) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
                       ) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
    ) #+
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


prior_model <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmax <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_addmin <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

g_lintemp_noint <- lmerTest::lmer(gweight ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled + meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),data = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size,age),
                                                  ~ scale(.x)[,1],
                                                   .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

(int_tab_growth_tres <- anova_int_tab(g_lintemp, g_lintemp_addmin, g_lintemp_addmax, g_lintemp_noint))


g_lintemp_addmin_tres <- g_lintemp_addmin

(anova(g_lintemp_addmin,prior_model))


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


samp <- g_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

data_tres = dplyr::filter(dat_growth,Species == "TRES",!is.na(meanmaxtempI),meanmaxtempI < 45,!is.na(meanmintempI)) %>%
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
                       ) +
   #geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
                       ) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
    ) #+
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


prior_model <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled * habitat + meanmint_nestpd_scaled + juliandate_hatch_scaled + year_fct + site,
                            family = binomial(link = "logit"),
                            data = dplyr::filter(dat_surv,
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
                     data = dplyr::filter(dat_surv,
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
                            data = dplyr::filter(dat_surv,
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
                            data = dplyr::filter(dat_surv,
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
                           data = dplyr::filter(dat_surv,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_webl <- anova_int_tab(s_nestpd_WEBL, s_nestpd_WEBL_addmin, s_nestpd_WEBL_addmax, s_nestpd_WEBL_noint))


anova(s_nestpd_WEBL_addmin, prior_model)


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


data_webl = s_nestpd_WEBL_addmin$data


mean_temp_webl <- mean(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(meanmaxt_nestpd),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                          transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                          inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())


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
                       ) +
   scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dashed","Grassland" = "dotted","Row crop" = "dashed")) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
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
                       ) +
   scale_linetype_manual(values = c("Forest" = "dashed","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
   geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
    ) #+
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


prior_model <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxt_nestpd_scaled + meanmint_nestpd_scaled + habitat + juliandate_hatch_scaled + year_fct + site,
                           family = binomial(link = "logit"),
                           data = dplyr::filter(dat_surv,
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
                     data = dplyr::filter(dat_surv,
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
                            data = dplyr::filter(dat_surv,
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
                            data = dplyr::filter(dat_surv,
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
                           data = dplyr::filter(dat_surv,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxt_nestpd,meanmint_nestpd,juliandate_hatch),
                                           ~ scale(.x)[,1],
                                           .names = "{.col}_scaled")
                             )
)

(int_tab_survival_nestpd_tres <- anova_int_tab(s_nestpd_TRES, s_nestpd_TRES_addmin, s_nestpd_TRES_addmax, s_nestpd_TRES_noint))


anova(s_nestpd_TRES_noint, prior_model)


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


data_tres = s_nestpd_TRES_noint$data


mean_temp_tres <- mean(data_tres %>% pull(meanmaxt_nestpd),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(meanmaxt_nestpd),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                          transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                          inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())


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
                       ) +
   scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
                       ) +
   scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
   geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
   )


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
    ) #+
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


dat_growth <- read_rds("data/growth_cort_provis_manytempmeasures.rds") %>%
  mutate(year_fct = as.factor(year))

data_s1_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + habitat + meanmintempI_scaled + habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_webl)

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_webl)

(int_tab_s1_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_webl <- s1_lintemp_noint

anova(s1_lintemp_noint,prior_model)


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


samp_s1_webl <- s1_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s1_webl <- data.frame(
  label = paste("N =",samp_s1_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) #+
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    #geom_text(data = dat_text_s1_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    #theme(legend.position = "none")
)


data_s1_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_tres)

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_tres)

(int_tab_s1_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_tres <- s1_lintemp_addmax

anova(s1_lintemp_addmax,prior_model)


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


samp_s1_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s1_tres <- data.frame(
  label = paste("N =",samp_s1_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


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
    ) +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


samp_s1_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s1_tres <- data.frame(
  label = paste("N =",samp_s1_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

data_s1_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))


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
    ) +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) #+
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # geom_text(data = dat_text_s1_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    #theme(legend.position = "none")
)


data_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl)

(int_tab_abs_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl <- abs_lintemp

anova(abs_lintemp,prior_model)


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


samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    ylim(0,60) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    ylim(0,60) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres)

prior_model_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres)

(int_tab_abs_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))

abs_tres <- abs_lintemp


anova(abs_lintemp,prior_model)

anova(abs_lintemp_addmin,prior_model_addmin)


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


samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dashed","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_s2_webl <- prep_model_data(
  dat_growth, "WEBL", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_webl)

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_webl)

(int_tab_s2_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_webl <- s2_lintemp


anova(s2_lintemp,prior_model)


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


samp_s2_webl <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s2_webl <- data.frame(
  label = paste("N =",samp_s2_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s2_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_s2_tres <- prep_model_data(
  dat_growth, "TRES", "meanmaxtempI", "meanmintempI",
  c("gweight","meanmaxtempI","meanmintempI","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_tres)

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled * habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ meanmaxtempI_scaled + meanmintempI_scaled + habitat + age_scaled +
    meanmaxtempI_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_tres)

(int_tab_s2_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_tres <- s2_lintemp_addmin


anova(s2_lintemp_addmin,prior_model)


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


samp_s2_tres <- s2_lintemp_addmin@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s2_tres <- data.frame(
  label = paste("N =",samp_s2_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    geom_text(data = dat_text_s2_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","solid","dotted")) +
    # geom_text(data = dat_text_s2_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_s1_priordayt_webl <- prep_model_data(
  dat_growth, "WEBL", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_webl)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_webl)

(int_tab_s1_priordayt_webl <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_prior_day_webl <- s1_lintemp_addmax


anova(s1_lintemp_addmax,prior_model)


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


samp_s1_priordayt_webl <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s1_priordayt_webl <- data.frame(
  label = paste("N =",samp_s1_priordayt_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s1_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s1_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_s1_priordayt_tres <- prep_model_data(
  dat_growth, "TRES", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_tres)

s1_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s1) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s1_priordayt_tres)

(int_tab_s1_priordayt_tres <- anova_int_tab(s1_lintemp, s1_lintemp_addmin, s1_lintemp_addmax, s1_lintemp_noint, digits = 3))

s1_prior_day_tres <- s1_lintemp_addmin


anova(s1_lintemp_addmax,prior_model)


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


samp_s1_priordayt_tres <- s1_lintemp_addmax@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s1_priordayt_tres <- data.frame(
  label = paste("N =",samp_s1_priordayt_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("solid","solid","solid","solid")) +
    geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("solid","solid","solid","solid")) +
    geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
#     scale_linetype_manual(values = c("dotted","dashed","dotted","dotted")) +
#     #ylim(-5,5) +
#     geom_text(data = dat_text_s1_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
#     theme(legend.position = "none")
)


data_webl_priordayt <- prep_model_data(
  dat_growth, "WEBL", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_priordayt)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_webl_priordayt)

(int_tab_abs_priordayt_webl <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_webl_priordayt <- abs_lintemp_noint
anova(abs_lintemp_noint,prior_model)


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


samp <- abs_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # # ylim(0,60) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_tres_priordayt <- prep_model_data(
  dat_growth, "TRES", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +
     maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp_addmax <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled +
     maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp_addmin <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled +
     maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_priordayt)

abs_lintemp_noint <- lmerTest::lmer(
  sqrt(abs_change_cort) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled +
     maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_tres_priordayt)

(int_tab_abs_priordayt_tres <- anova_int_tab(abs_lintemp, abs_lintemp_addmin, abs_lintemp_addmax, abs_lintemp_noint, digits = 3))


abs_tres_priordayt <- abs_lintemp

anova(abs_lintemp,prior_model)


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


samp <- abs_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "solid","Grassland" = "dotted","Row crop" = "dotted")) +
    # # ylim(0,60) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_s2_priordayt_webl <- prep_model_data(
  dat_growth, "WEBL", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + juliandate_scaled + year_fct +
    (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_webl)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_webl)

(int_tab_s2_priordayt_webl <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_prior_day_webl <- s2_lintemp_noint


anova(s2_lintemp_noint,prior_model)


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


samp_s2_priordayt_webl <- s2_lintemp_noint@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s2_priordayt_webl <- data.frame(
  label = paste("N =",samp_s2_priordayt_webl$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s2_priordayt_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


data_s2_priordayt_tres <- prep_model_data(
  dat_growth, "TRES", "maxt_prior", "mint_prior",
  c("gweight","maxt_prior","mint_prior","juliandate","brood_size","age"))

prior_model <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled + juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled * habitat + age_scaled +
    maxt_prior_scaled * juliandate_scaled + year_fct + (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp_addmax <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled * habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp_addmin <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled * habitat + mint_prior_scaled + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_tres)

s2_lintemp_noint <- lmerTest::lmer(
  sqrt(cort_s2) ~ maxt_prior_scaled + mint_prior_scaled + habitat + age_scaled + maxt_prior_scaled * juliandate_scaled +
    year_fct + (1|attempt_id),
  data = data_s2_priordayt_tres)

(int_tab_s2_priordayt_tres <- anova_int_tab(s2_lintemp, s2_lintemp_addmin, s2_lintemp_addmax, s2_lintemp_noint, digits = 3))

s2_prior_day_tres <- s2_lintemp


anova(s2_lintemp,prior_model)


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


samp_s2_priordayt_tres <- s2_lintemp@frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_s2_priordayt_tres <- data.frame(
  label = paste("N =",samp_s2_priordayt_tres$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)



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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dashed")) +
    geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("dotted","dotted","dotted","dotted")) +
    geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("dotted","dotted","dotted","dashed")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_s2_priordayt_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


prior_model <- glmmTMB(
  valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) * habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_linint <- glmmTMB(
  valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + mean_temp_scaled * julian_date_scaled +
    poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_noint <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) + habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(dat_provis,!is.na(mean_temp),
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


provis_webl <- m_linint


anova(m_linint,prior_model)


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


samp <- m_linint$frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

data_webl = dplyr::filter(dat_provis,!is.na(mean_temp),
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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    # #ylim(-5,5) +
    # geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


prior_model <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) * habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
             ziformula = ~ 1,
             family = nbinom2(),
             data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_linint <- glmmTMB(
  valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + mean_temp_scaled * julian_date_scaled +
    poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                    ziformula = ~ 1,
                    family = nbinom2(),
                    data = dplyr::filter(dat_provis,!is.na(mean_temp),
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

m_noint <- glmmTMB(
  valid_detections ~ poly(mean_temp_scaled,2) + habitat + mean_temp_scaled * julian_date_scaled + poly(tod_h_scaled,2) +
    mean_nestling_age_scaled + (1|attempt_id) + year_fct,
                   ziformula = ~ 1,
                   family = nbinom2(),
                   data = dplyr::filter(dat_provis,!is.na(mean_temp),
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


provis_tres <- m

anova(m,prior_model)


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


samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)


data_tres = dplyr::filter(dat_provis,!is.na(mean_temp),
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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    ylim(0,50) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    ylim(0,50) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)


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
    ) # +
    # scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    # geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    # theme(legend.position = "none")
)


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


(trend <- rbind(growth,survival,s1_priorday,s1_priorweek,abs_priorday,abs_priorweek,s2_priorday,s2_priorweek,provis) %>% gt(rowname_col = "Model",groupname_col = "Response") %>% tab_options(data_row.padding = px(1)) %>%

  tab_spanner_delim(
    delim="_"))


save(list = ls(), file = "data/models_seasonal.RData")
