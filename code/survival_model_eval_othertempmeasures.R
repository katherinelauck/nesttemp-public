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
wmean <- read_rds("data/wmean.rds")

# g <- read_rds("data/growth_and_provis_combined_mobilenetv3-original_dataset.h5.rds") %>%
#   mutate(abs_change_cort = cort_s2-cort_s1,
#          prop_change_cort = (cort_s2-cort_s1)/cort_s1,
#          year_fct = as.factor(year))
p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

b <- read_csv("data/banding-and-morphometrics_proofed.csv")


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


## WEBL

### meanmaxhi_nestpd


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled * habitat + meanminhi_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                           data = dplyr::filter(s,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
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

(int_tab_survival_nestpd_webl_meanmaxhi <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
gtsave(int_tab_survival_nestpd_webl_meanmaxhi,"figures/int_tab_survival_nestpd_webl_meanmaxhi.html")



#Conclusion: For nest attempt overall, land cover interacts with min temp.


summary(s_nestpd_WEBL_addmax)
(survivalbyhabitat_summary_webl_meanmaxhi <- summary(s_nestpd_WEBL_addmax) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(survivalbyhabitat_summary_webl_meanmaxhi,"figures/survivalbyhabitat_summary_webl_meanmaxhi.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax_meanmaxhi <- emtrends(s_nestpd_WEBL_addmax,specs = pairwise ~ habitat, var = c("meanmaxhi_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxhi_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


gtsave(weblsurvival_nestpd_trendmax_meanmaxhi,"figures/weblsurvival_nestpd_trendmax_meanmaxhi.html")

data = s_nestpd_WEBL_addmax$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl_meanmaxhi <- emmeans(s_nestpd_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(survival_nestpd_byhabitat_webl_meanmaxhi,"figures/survival_nestpd_byhabitat_webl_meanmaxhi.html")


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


data_webl = s_nestpd_WEBL_addmax$data


mean_temp_webl <- mean(data_webl %>% pull(meanmaxhi_nestpd),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(meanmaxhi_nestpd),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_webl.html")
#
ss_year_survival_webl_meanmaxhi <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl_meanmaxhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
  gtsave("figures/ss_year_survival_webl_meanmaxhi.html")


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
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    # ylim(0,100) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/survbytempxhab_webl_meanmaxhi.png",fig3_webl_meanmaxhi,width = 10, height = 6.6)


### deghr_30


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                           data = dplyr::filter(s,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
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

(int_tab_survival_nestpd_webl_deghr_30 <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
gtsave(int_tab_survival_nestpd_webl_deghr_30,"figures/int_tab_survival_nestpd_webl_deghr_30.html")



#Conclusion: For nest attempt overall, land cover interacts with cumulative degree hours.


summary(s_nestpd_WEBL_addmin)
(survivalbyhabitat_summary_webl_deghr_30 <- summary(s_nestpd_WEBL_addmin) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(survivalbyhabitat_summary_webl_deghr_30,"figures/survivalbyhabitat_summary_webl_deghr_30.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax_deghr_30 <- emtrends(s_nestpd_WEBL_addmin,specs = pairwise ~ habitat, var = c("deghr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "deghr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


gtsave(weblsurvival_nestpd_trendmax_deghr_30,"figures/weblsurvival_nestpd_trendmax_deghr_30.html")

data = s_nestpd_WEBL_addmin$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl_deghr_30 <- emmeans(s_nestpd_WEBL_addmin,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(survival_nestpd_byhabitat_webl_deghr_30,"figures/survival_nestpd_byhabitat_webl_deghr_30.html")


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


mean_temp_webl <- mean(data_webl %>% pull(deghr_30),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(deghr_30),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_webl.html")
#
ss_year_survival_webl_deghr_30 <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl_deghr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
  gtsave("figures/ss_year_survival_webl_deghr_30.html")


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
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    # ylim(0,100) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/survbytempxhab_webl_deghr_30.png",fig3_webl_deghr_30,width = 10, height = 6.6)

### hihr_30


s_nestpd_WEBL <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                           data = dplyr::filter(s,
                                                Species == "WEBL",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
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

(int_tab_survival_nestpd_webl_hihr_30 <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
gtsave(int_tab_survival_nestpd_webl_hihr_30,"figures/int_tab_survival_nestpd_webl_hihr_30.html")



#Conclusion: For nest attempt overall, land cover interacts with cumulative degree hours.


summary(s_nestpd_WEBL_addmax)
(survivalbyhabitat_summary_webl_hihr_30 <- summary(s_nestpd_WEBL_addmax) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(survivalbyhabitat_summary_webl_hihr_30,"figures/survivalbyhabitat_summary_webl_hihr_30.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(weblsurvival_nestpd_trendmax_hihr_30 <- emtrends(s_nestpd_WEBL_addmax,specs = pairwise ~ habitat, var = c("hihr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


gtsave(weblsurvival_nestpd_trendmax_hihr_30,"figures/weblsurvival_nestpd_trendmax_hihr_30.html")

data = s_nestpd_WEBL_addmax$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_webl_hihr_30 <- emmeans(s_nestpd_WEBL_addmax,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(survival_nestpd_byhabitat_webl_hihr_30,"figures/survival_nestpd_byhabitat_webl_hihr_30.html")


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


data_webl = s_nestpd_WEBL_addmax$data


mean_temp_webl <- mean(data_webl %>% pull(hihr_30),na.rm = TRUE)
sd_temp_webl <- sd(data_webl %>% pull(hihr_30),na.rm = TRUE)


temp_trans_webl <- trans_new("temp_trans_webl",
                             transform = function(x){(x * sd_temp_webl) + mean_temp_webl},
                             inverse = function(x){x})

samp_webl <- data_webl %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_webl.html")
#
ss_year_survival_webl_hihr_30 <- data_webl %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_webl_hihr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
  gtsave("figures/ss_year_survival_webl_hihr_30.html")


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
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    # ylim(0,100) +
    scale_linetype_manual(values = c("Forest" = "dotted","Orchard" = "dotted","Grassland" = "dotted","Row crop" = "dotted")) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/survbytempxhab_webl_hihr_30.png",fig3_webl_hihr_30,width = 10, height = 6.6)

## TRES

### meanmaxhi_nestpd


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ meanmaxhi_nestpd_scaled * habitat + meanminhi_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                           data = dplyr::filter(s,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(meanmaxhi_nestpd,meanminhi_nestpd,juliandate_hatch),
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

(int_tab_survival_nestpd_tres_meanmaxhi <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
gtsave(int_tab_survival_nestpd_tres_meanmaxhi,"figures/int_tab_survival_nestpd_tres_meanmaxhi.html")



#Conclusion: For nest attempt overall, land cover interacts with min temp.


summary(s_nestpd_TRES_addmax)
(survivalbyhabitat_summary_tres_meanmaxhi <- summary(s_nestpd_TRES_addmax) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(survivalbyhabitat_summary_tres_meanmaxhi,"figures/survivalbyhabitat_summary_tres_meanmaxhi.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(tressurvival_nestpd_trendmax_meanmaxhi <- emtrends(s_nestpd_TRES_addmax,specs = pairwise ~ habitat, var = c("meanmaxhi_nestpd_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "meanmaxhi_nestpd_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


gtsave(tressurvival_nestpd_trendmax_meanmaxhi,"figures/tressurvival_nestpd_trendmax_meanmaxhi.html")

data = s_nestpd_TRES_addmax$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_tres_meanmaxhi <- emmeans(s_nestpd_TRES_addmax,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(survival_nestpd_byhabitat_tres_meanmaxhi,"figures/survival_nestpd_byhabitat_tres_meanmaxhi.html")


# Forest survival is lower than in the other land covers.

### Check for effect of temperature


# (t <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
#           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nestpd_byhabitat_summary_tres.html")
#
#
# Hot temps and low temps reduce survival overall.


data_tres = s_nestpd_TRES_addmax$data


mean_temp_tres <- mean(data_tres %>% pull(meanmaxhi_nestpd),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(meanmaxhi_nestpd),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_tres.html")
#
ss_year_survival_tres_meanmaxhi <- data_tres %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_tres_meanmaxhi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
  gtsave("figures/ss_year_survival_tres_meanmaxhi.html")


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

ggsave("figures/survbytempxhab_tres_meanmaxhi.png",fig3_tres_meanmaxhi,width = 10, height = 6.6)


### deghr_30


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ deghr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                           data = dplyr::filter(s,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(deghr_30,meanmint_nestpd,juliandate_hatch),
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

(int_tab_survival_nestpd_tres_deghr_30 <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
gtsave(int_tab_survival_nestpd_tres_deghr_30,"figures/int_tab_survival_nestpd_tres_deghr_30.html")



#Conclusion: For nest attempt overall, land cover interacts with cumulative degree hours.


summary(s_nestpd_TRES_noint)
(survivalbyhabitat_summary_tres_deghr_30 <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(survivalbyhabitat_summary_tres_deghr_30,"figures/survivalbyhabitat_summary_tres_deghr_30.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(tressurvival_nestpd_trendmax_deghr_30 <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("deghr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "deghr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


gtsave(tressurvival_nestpd_trendmax_deghr_30,"figures/tressurvival_nestpd_trendmax_deghr_30.html")

data = s_nestpd_TRES_noint$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_tres_deghr_30 <- emmeans(s_nestpd_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(survival_nestpd_byhabitat_tres_deghr_30,"figures/survival_nestpd_byhabitat_tres_deghr_30.html")


# Forest survival is lower than in the other land covers.

### Check for effect of temperature


# (t <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
#           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nestpd_byhabitat_summary_tres.html")
#
#
# Hot temps and low temps reduce survival overall.


data_tres = s_nestpd_TRES_noint$data


mean_temp_tres <- mean(data_tres %>% pull(deghr_30),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(deghr_30),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_tres.html")
#
ss_year_survival_tres_deghr_30 <- data_tres %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_tres_deghr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
  gtsave("figures/ss_year_survival_tres_deghr_30.html")


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
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    # ylim(0,100) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "solid")) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = -Inf,label = label),hjust = -.2,vjust = -.7,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/survbytempxhab_tres_deghr_30.png",fig3_tres_deghr_30,width = 10, height = 6.6)

### hihr_30


s_nestpd_TRES <- glm(cbind(nest_fledged,clutch_size - nest_fledged) ~ hihr_30_scaled * habitat + meanmint_nestpd_scaled * habitat + juliandate_hatch_scaled + year_fct + site,
                     family = binomial(link = "logit"),
                     data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                            data = dplyr::filter(s,
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
                           data = dplyr::filter(s,
                                                Species == "TRES",
                                                !is.na(nest_fledged),
                                                !is.na(clutch_size)) %>%
                             mutate(across(c(hihr_30,meanmint_nestpd,juliandate_hatch),
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

(int_tab_survival_nestpd_tres_hihr_30 <- bind_rows(c1,c2) %>%
    as_tibble() %>%
    mutate(across(where(is.numeric),~round(.x,digits = 4)),
           P = `Pr(>Chi)`) %>%
    dplyr::select(Model,AIC,Chisq,P) %>%
    mutate(max_or_min = c(rep("Max temp",times = 3),rep("Min temp",times = 3)),
           across(c(Chisq), ~ round(.x, digits = 2)),
           P = if_else(P < 0.001,"<0.001",as.character(P))) %>%
    group_by(max_or_min) %>% gt())
gtsave(int_tab_survival_nestpd_tres_hihr_30,"figures/int_tab_survival_nestpd_tres_hihr_30.html")



#Conclusion: For nest attempt overall, land cover interacts with cumulative degree hours.


summary(s_nestpd_TRES_noint)
(survivalbyhabitat_summary_tres_hihr_30 <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(survivalbyhabitat_summary_tres_hihr_30,"figures/survivalbyhabitat_summary_tres_hihr_30.html")



# There is a direct effect of min temp in the direction we expect (higher min temp = higher survival).

#### Emtrends to calculate effect of max temp in each habitat


(tressurvival_nestpd_trendmax_hihr_30 <- emtrends(s_nestpd_TRES_noint,specs = pairwise ~ habitat, var = c("hihr_30_scaled")) %>% test() %>% pluck("emtrends") %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hihr_30_scaled.trend",Df = "df", `T-ratio` = "z.ratio", P = "p.value") %>%
    gt())


gtsave(tressurvival_nestpd_trendmax_hihr_30,"figures/tressurvival_nestpd_trendmax_hihr_30.html")

data = s_nestpd_TRES_noint$data


## Emmeans to check for effect of habitat


(survival_nestpd_byhabitat_tres_hihr_30 <- emmeans(s_nestpd_TRES_noint,"habitat") %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(survival_nestpd_byhabitat_tres_hihr_30,"figures/survival_nestpd_byhabitat_tres_hihr_30.html")


# Forest survival is lower than in the other land covers.

### Check for effect of temperature


# (t <- summary(s_nestpd_TRES_noint) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
#    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
#    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
#           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
# gtsave(t,"../figures/survival_nestpd_byhabitat_summary_tres.html")
#
#
# Hot temps and low temps reduce survival overall.


data_tres = s_nestpd_TRES_noint$data


mean_temp_tres <- mean(data_tres %>% pull(hihr_30),na.rm = TRUE)
sd_temp_tres <- sd(data_tres %>% pull(hihr_30),na.rm = TRUE)


temp_trans_tres <- trans_new("temp_trans_tres",
                             transform = function(x){(x * sd_temp_tres) + mean_temp_tres},
                             inverse = function(x){x})

samp_tres <- data_tres %>% group_by(habitat) %>% summarize(count = n())
# samp %>% gt() %>% gtsave("../figures/ss_tres.html")
#
ss_year_survival_tres_hihr_30 <- data_tres %>% group_by(habitat,year_fct) %>% summarize(count = n()) %>%
  pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>% ungroup() %>% mutate(Total = `2021` + `2022` + `2023`)

ss_year_survival_tres_hihr_30 %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.)) %>%
  gtsave("figures/ss_year_survival_tres_hihr_30.html")


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

ggsave("figures/survbytempxhab_tres_hihr_30.png",fig3_tres_hihr_30,width = 10, height = 6.6)

