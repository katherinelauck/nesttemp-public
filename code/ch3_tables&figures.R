#### Make tables and figures for ch3 manuscript

setwd("code")
source(knitr::purl("ch3_analysis_reorg_after_meeting_with_Danny.Rmd"))
setwd("..")

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

### Tables

#### Sample size

provis <- rbind(t_ss_year_webl_provis$`_data`,t_ss_year_tres_provis$`_data`) %>% dplyr::select(`2021`,`2022`, `2023`) %>%
  replace_na(replace = list(`2021` = 0,`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Provisioning ~ Mean hourly temp\U00B2",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_priorday <- rbind(t_ss_year_webl_provis_q2h1_maxt_prior$`_data`,t_ss_year_tres_provis_q2h1_maxt_prior$`_data`) %>% dplyr::select(`2021`,`2022`, `2023`) %>%
  replace_na(replace = list(`2021` = 0,`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Provisioning ~ Mean hourly temp\U00B2 * Prior day max temp",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_res <- rbind(t_ss_year_webl_provis_q2h2_res$`_data`,t_ss_year_tres_provis_q2h2_res$`_data`) %>% dplyr::select(`2022`, `2023`) %>%
  replace_na(replace = list(`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Provisioning ~ Mean hourly temp\U00B2 * Female condition",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_s1 <- rbind(t_ss_year_webl_provis_q2h3$`_data`,t_ss_year_tres_provis_q2h3$`_data`) %>% dplyr::select(`2022`, `2023`) %>%
  replace_na(replace = list(`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Provisioning ~ Mean hourly temp\U00B2 * Female baseline corticosterone",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_s2 <- rbind(t_ss_year_webl_provis_q4h2$`_data`,t_ss_year_tres_provis_q4h2$`_data`) %>% dplyr::select(`2022`, `2023`) %>%
  replace_na(replace = list(`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Provisioning ~ Mean hourly temp\U00B2 * Female stress-induced corticosterone",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

condition_by_temp <- rbind(t_ss_year_webl_q3h1_priorday$`_data`,t_ss_year_tres_q3h1_priorday$`_data`) %>% dplyr::select(`2022`, `2023`) %>%
  replace_na(replace = list(`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Female condition ~ Temp\U00B2",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

s1_by_temp <- rbind(t_ss_year_webl_q3h2$`_data`,t_ss_year_tres_q3h2$`_data`) %>% dplyr::select(`2022`, `2023`) %>%
  replace_na(replace = list(`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Baseline corticosterone ~ Temp\U00B2",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

s2_by_temp <- rbind(t_ss_year_webl_q3h3$`_data`,t_ss_year_tres_q3h3$`_data`) %>% dplyr::select(`2022`, `2023`) %>%
  replace_na(replace = list(`2022` = 0, `2023` = 0)) %>%
  mutate(Species = c(rep("Western Bluebird",1),rep("Tree Swallow",1)),
         Model = "Stress-induced corticosterone ~ Temp\U00B2",
         across(everything(),~ as.character(.))) %>%
  group_by(Model) %>%
  pivot_longer(-c(Species, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

(samp_size_onecol <- rbind(provis,provis_by_priorday,provis_by_res,provis_by_s1,provis_by_s2,condition_by_temp,s1_by_temp,s2_by_temp) %>%
    mutate(across(everything(),~replace_na(.x,"0, 0")),
           across(everything(),~str_replace_all(.x,"NA","0"))) %>%
    ungroup() %>%
    gt() %>% tab_options(data_row.padding = px(1)) %>%
    tab_spanner_delim(
      delim="_") # %>%
    # summary_rows(fns = list(id = "Total") ~ {
    #   t <- str_split(.,", ",simplify = TRUE)
    #   if(ncol(t) == 2){
    #     t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% paste(t[,2] %>% as.numeric() %>% sum(na.rm = TRUE),sep = ", ")
    #   } else {
    #     t[,1] %>% as.numeric() %>% sum(na.rm = TRUE) %>% as.character()
    #   }
    #
    # })
)

gtsave(samp_size_onecol,"figures/samp_size_tbl_ch3.html")

#### Effect size table

provis_by_priorday <- rbind(weblprovistrend_q2h1_maxt_prior$`_data`,tresprovistrend_q2h1_maxt_prior$`_data`) %>%
  rename(`Mediator value` = `Prior day maximum temperatures`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         Mediator = "Prior day maximum temperature",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, Mediator,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_cond <- rbind(weblprovistrend_q2h2_res$`_data`,tresprovistrend_q2h2_res$`_data`) %>%
  rename(`Mediator value` = `Female condition`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         Mediator = "Female condition",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, Mediator,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_s1 <- rbind(weblprovistrend_q2h3$`_data`,tresprovistrend_q2h3$`_data`) %>%
  rename(`Mediator value` = `Female baseline corticosterone`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         Mediator = "Female baseline corticosterone",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, Mediator,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

provis_by_s2 <- rbind(weblprovistrend_q4h2$`_data`,tresprovistrend_q4h2$`_data`) %>%
  rename(`Mediator value` = `Female stress-induced corticosterone`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         Mediator = "Female stress-induced corticosterone",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, Mediator,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

(effect_size_q2 <- rbind(provis_by_priorday,provis_by_cond,provis_by_s1,provis_by_s2) %>% gt(rowname_col = "Degree",groupname_col = "Mediator") %>% tab_options(data_row.padding = px(1)) %>%

    tab_spanner_delim(
      delim="_"))

gtsave(effect_size_q2,"figures/effect_size_q2_ch3.html")

#### Effect size table for sensitivity analyses

priorday <- rbind(weblprovistrend_q2h1_maxt_prior$`_data`,tresprovistrend_q2h1_maxt_prior$`_data`) %>%
  rename(`Mediator value` = `Prior day maximum temperatures`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         `Temp measure` = "Prior day maximum temperature",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, `Temp measure`,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

deghr_above30 <- rbind(weblprovistrend_q2h1_deghrover30$`_data`,TRESprovistrend_q2h1_deghrover30$`_data`) %>%
  rename(`Mediator value` = `Prior day degree-hours difference from 30C`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         `Temp measure` = "Prior day degree-hours above 30\u00b0C",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, `Temp measure`,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

deghr_above35 <- rbind(weblprovistrend_q2h1_deghrover35$`_data`,TRESprovistrend_q2h1_deghrover35$`_data`) %>%
  rename(`Mediator value` = `Prior day degree-hours difference from 30C`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         `Temp measure` = "Prior day degree-hours above 35\u00b0C",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, `Temp measure`,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

deghr_diff30 <- rbind(weblprovistrend_q2h1_deghrdiff30$`_data`,TRESprovistrend_q2h1_deghrdiff30$`_data`) %>%
  rename(`Mediator value` = `Prior day degree-hours difference from 30C`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         `Temp measure` = "Prior day degree-hours difference from 30\u00b0C",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, `Temp measure`,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

deghr_diff35 <- rbind(weblprovistrend_q2h1_deghrdiff35$`_data`,tresprovistrend_q2h1_deghrdiff35$`_data`) %>%
  rename(`Mediator value` = `Prior day degree-hours difference from 35C`,Degree = "degree",`Statistic` = `Z-ratio`,`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",4),rep("Tree Swallow",4)),
         `Temp measure` = "Prior day degree-hours difference from 35\u00b0C",
         across(everything(),~ as.character(.)),
         `Mediator level` = rep(c("-2 SD","+2 SD"),4)) %>%
  arrange(desc(`Mediator level`)) %>%
  # mutate(row=row_number()) %>%
  pivot_longer(-c(Species, `Temp measure`,`Mediator level`,Degree)) %>%
  #print(n = Inf)
  pivot_wider(names_from=c(Species, name), values_from=value)

(effect_size_sensitivity <- rbind(priorday,deghr_above30,deghr_above35) %>% gt(rowname_col = "Degree",groupname_col = "Temp measure") %>% tab_options(data_row.padding = px(1)) %>%

    tab_spanner_delim(
      delim="_"))

gtsave(effect_size_sensitivity,"figures/effect_size_sensitivity_ch3.html")

#### Coefficient table q1

webl <- summary(m_webl_q1h1) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird")

tres <- summary(m_tres_q1h1) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow")

name_key_q1h1 <- tibble(old = c(webl$Covariate,tres$Covariate) %>% unique(),
                          new = c("Intercept","Mean hourly temperature (scaled)",
                                  "Mean hourly temperature\U00B2 (scaled)",
                                  "Nestling age",
                                  "Brood size",
                                  "Year: 2022",
                                  "Year: 2023",
                                  "Time of day",
                                  "Time of day\U00B2"))

(q1h1_coef <- rbind(webl,tres) %>%
    rename(SE = `Std. Error`,`Z-value` = `z value`,P = `Pr(>|z|)`) %>%
    filter(!Covariate=="(Intercept)") %>%
    mutate(Covariate = name_key_q1h1$new[match(Covariate, name_key_q1h1$old)]) %>%
    pivot_longer(-c(Species,Covariate)) %>%
    mutate(value = as.character(value),
           value = if_else(name == "P" & value == "0","<0.001",value)) %>%
    pivot_wider(names_from = c(Species,name),values_from = value) %>% gt() %>%
    tab_options(data_row.padding = px(1)) %>%
    tab_spanner_delim(
      delim="_"))

gtsave(q1h1_coef,"figures/q1h1_coef_ch3.html")

#### Coefficient table q2

webl_h1 <- summary(m_webl_q2h1_maxt_prior) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Mediator = "Prior day maximum temperature")

tres_h1 <- summary(m_tres_q2h1_maxt_prior) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Mediator = "Prior day maximum temperature")

webl_h2 <- summary(m_webl_q2h2_res) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Mediator = "Female condition") %>%
  filter(Covariate != "poly(mean_temp_scaled, degree = 2)1")

tres_h2 <- summary(m_tres_q2h2_res) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Mediator = "Female condition") %>%
  filter(Covariate != "poly(mean_temp_scaled, degree = 2)1")

webl_h3 <- summary(m_webl_q2h3) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Mediator = "Baseline corticosterone")

tres_h3 <- summary(m_tres_q2h3) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Mediator = "Baseline corticosterone")

webl_h4 <- summary(m_webl_q4h2) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Mediator = "Stress-induced corticosterone")

tres_h4 <- summary(m_tres_q4h2) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
         across(`Pr(>|z|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Mediator = "Stress-induced corticosterone")

name_key_q2 <- tibble(old = c(webl_h1$Covariate,tres_h1$Covariate,webl_h2$Covariate,tres_h2$Covariate,webl_h3$Covariate,tres_h3$Covariate,webl_h4$Covariate,tres_h4$Covariate) %>% unique(),
                        new = c("Intercept","Mean hourly temperature (scaled)",
                                "Mean hourly temperature\U00B2 (scaled)",
                                "Prior day maximum temperature",
                                "Time of day",
                                "Time of day\U00B2",
                                "Nestling age",
                                "Brood size",
                                "Year: 2022",
                                "Year: 2023",
                                "Mean hourly temperature * prior day maximum temperature",
                                "Mean hourly temperature\U00B2 * prior day maximum temperature",
                                "Mean hourly temperature (scaled)",
                                "Female condition",
                                "Mean hourly temperature\U00B2 (scaled)",
                                "Mean hourly temperature * female condition",
                                "Baseline corticosterone",
                                "Mean hourly temperature * baseline corticosterone",
                                "Mean hourly temperature\U00B2 * baseline corticosterone",
                                "Mean hourly temperature (scaled)",
                                "Stress-induced corticosterone"))

(q2h1_coef <- rbind(webl_h1,tres_h1,webl_h2,tres_h2,webl_h3,tres_h3,webl_h4,tres_h4) %>%
    rename(SD = `Std. Error`,`Z-value` = `z value`,P = `Pr(>|z|)`) %>%
    filter(!Covariate=="(Intercept)") %>%
    mutate(Covariate = name_key_q2$new[match(Covariate, name_key_q2$old)]) %>%
    pivot_longer(-c(Species,Covariate,Mediator)) %>%
    mutate(value = as.character(value),
           value = if_else(name == "P" & value == "0","<0.001",value)) %>%
    pivot_wider(names_from = c(Species,name),values_from = value) %>% gt(groupname_col = "Mediator") %>%
    tab_options(data_row.padding = px(1)) %>%
    tab_spanner_delim(
      delim="_"))

gtsave(q2h1_coef,"../figures/q2h1_coef_ch3.html")


#### Coefficient table q3

webl_h1 <- summary(m_webl_q3h1_priorday) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Model = "Female condition ~ prior day maximum temperature")

tres_h1 <- summary(m_tres_q3h1_priorday) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Model = "Female condition ~ prior day maximum temperature")

webl_h1_priorweek <- summary(m_webl_q3h1) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Model = "Female condition ~ prior week average daily maximum temperature")

tres_h1_priorweek <- summary(m_tres_q3h1) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Model = "Female condition ~ prior week average daily maximum temperature")

webl_h2 <- summary(m_webl_q3h2) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Model = "Baseline corticosterone ~ prior day maximum temperature")

tres_h2 <- summary(m_tres_q3h2) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Model = "Baseline corticosterone ~ prior day maximum temperature")

webl_h2_priorweek <- summary(m_webl_q3h2_priorweek) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Model = "Baseline corticosterone ~ prior week average daily maximum temperature")

tres_h2_priorweek <- summary(m_tres_q3h2_priorweek) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Model = "Baseline corticosterone ~ prior week average daily maximum temperature")

webl_h3 <- summary(m_webl_q3h3_priorday) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Model = "Stress-induced corticosterone ~ prior day maximum temperature")

tres_h3 <- summary(m_tres_q3h3_priorday) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Model = "Stress-induced corticosterone ~ prior day maximum temperature")

webl_h3_priorweek <- summary(m_webl_q3h3) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Western Bluebird",
         Model = "Stress-induced corticosterone ~ prior week average daily maximum temperature")

tres_h3_priorweek <- summary(m_tres_q3h3) %>% coef() %>% as_tibble(rownames = "Covariate") %>%
  # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
  mutate(across(Estimate:`t value`,~round(.x,digits = 2)),
         across(`Pr(>|t|)`,~round(.x,digits = 3))) %>%
  mutate(Species = "Tree Swallow",
         Model = "Stress-induced corticosterone ~ prior week average daily maximum temperature")

name_key_q3 <- tibble(old = c(webl_h1$Covariate,tres_h1$Covariate,webl_h1_priorweek$Covariate,tres_h1_priorweek$Covariate,webl_h2$Covariate,tres_h2$Covariate,webl_h2_priorweek$Covariate,tres_h2_priorweek$Covariate,webl_h3$Covariate,tres_h3$Covariate,webl_h3_priorweek$Covariate,tres_h3_priorweek$Covariate) %>% unique(),
                      new = c("Intercept","Prior day maximum temperature",
                              "Prior day maximum temperature\U00B2",
                              "Brood size",
                              "Day of year",
                              "Year: 2023",
                              "Site: Picnic Grounds",
                              "Site: Russell Ranch",
                              "Prior week average daily maximum temperature",
                              "Prior week average daily maximum temperature\U00B2"))

(q3_coef <- rbind(webl_h1,tres_h1,webl_h1_priorweek,tres_h1_priorweek,webl_h2,tres_h2,webl_h2_priorweek,tres_h2_priorweek,webl_h3,tres_h3,webl_h3_priorweek,tres_h3_priorweek) %>%
    rename(SD = `Std. Error`,`T-value` = `t value`,P = `Pr(>|t|)`) %>%
    filter(!Covariate=="(Intercept)") %>%
    mutate(Covariate = name_key_q3$new[match(Covariate, name_key_q3$old)]) %>%
    pivot_longer(-c(Species,Covariate,Model)) %>%
    mutate(value = as.character(value),
           value = if_else(name == "P" & value == "0","<0.001",value)) %>%
    pivot_wider(names_from = c(Species,name),values_from = value) %>% gt(groupname_col = "Model") %>%
    tab_options(data_row.padding = px(1)) %>%
    tab_spanner_delim(
      delim="_"))

gtsave(q3_coef,"../figures/q3_coef_ch3.html")

#### Maximum likelihood ratio tests

priorday <- rbind(int_tab_ch3_q2h1_maxt_prior_webl$`_data`,int_tab_ch3_q2h1_maxt_prior_tres$`_data`) %>%
  rename(`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",3),rep("Tree Swallow",3)),
         Mediator = "Prior day maximum temperature",
         AIC = as.character(AIC),
         Chisq = as.character(Chisq),
         `P-value` = as.character(`P-value`),
         Model = rep(c("No temp * mediator interaction","Linear temp * mediator interaction","Temp\U00B2 * mediator interaction"),2)) %>%
  group_by(Mediator) %>%
  pivot_longer(-c(Species, Mediator, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

cond <- rbind(int_tab_webl_q2h2$`_data`,int_tab_tres_q2h2$`_data`) %>%
  rename(`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",3),rep("Tree Swallow",3)),
         Mediator = "Female condition",
         AIC = as.character(AIC),
         Chisq = as.character(Chisq),
         `P-value` = as.character(`P-value`),
         Model = rep(c("No temp * mediator interaction","Linear temp * mediator interaction","Temp\U00B2 * mediator interaction"),2)) %>%
  group_by(Mediator) %>%
  pivot_longer(-c(Species, Mediator, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

s1 <- rbind(int_tab_webl_q2h3$`_data`,int_tab_tres_q2h3$`_data`) %>%
  rename(`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",3),rep("Tree Swallow",3)),
         Mediator = "Baseline corticosterone",
         AIC = as.character(AIC),
         Chisq = as.character(Chisq),
         `P-value` = as.character(`P-value`),
         Model = rep(c("No temp * mediator interaction","Linear temp * mediator interaction","Temp\U00B2 * mediator interaction"),2)) %>%
  group_by(Mediator) %>%
  pivot_longer(-c(Species, Mediator, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

s2 <- rbind(int_tab_webl_q2h4$`_data`,int_tab_tres_q2h4$`_data`) %>%
  rename(`P-value` = "P") %>%
  mutate(Species = c(rep("Western Bluebird",3),rep("Tree Swallow",3)),
         Mediator = "Stress-induced corticosterone",
         AIC = as.character(AIC),
         Chisq = as.character(Chisq),
         `P-value` = as.character(`P-value`),
         Model = rep(c("No temp * mediator interaction","Linear temp * mediator interaction","Temp\U00B2 * mediator interaction"),2)) %>%
  group_by(Mediator) %>%
  pivot_longer(-c(Species, Mediator, Model)) %>%
  pivot_wider(names_from=c(Species, name), values_from=value)

(int_tab <- rbind(priorday,cond,s1,s2) %>% gt(rowname_col = "Model",groupname_col = "Mediator") %>% tab_options(data_row.padding = px(1)) %>%

    tab_spanner_delim(
      delim="_"))

gtsave(int_tab,"../figures/int_tbl_ch3.html")



#### Figure 2: provisioning by hourly temp

ggplot_build(figq1h1_webl)$layout$panel_scales_y
ggplot_build(figq1h1_tres)$layout$panel_scales_y
(p_full <- ggarrange(figq1h1_webl+ labs(title = NULL) + theme(text = element_text(size = 12),axis.title.x = element_blank()) + ylim(0,65),figq1h1_tres + ylim(0,65) + labs(title = NULL) + theme(axis.text.y = element_blank(),
                                                                                                                                           axis.ticks.y = element_blank(),
                                                                                                                                           text = element_text(size = 12),
                                                                                                                                           axis.title.y = element_blank(),
                                                                                                                                           axis.title.x = element_text(hjust = -.5)
                                                                                                                                           ),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/ch3_fig2_provisbytemp.png",p_full,width = 6.25,height = 4)

#### Figure 3: WEBL Provis by prior day max temp, condition, baseline cort

ggplot_build(figq2h1_maxt_prior_webl)$layout$panel_scales_y
ggplot_build(figq2h2_res_webl)$layout$panel_scales_y
ggplot_build(figq2h3_webl)$layout$panel_scales_y
(p_full <- ggarrange(figq2h1_maxt_prior_webl +
                       theme(legend.position = "inside",legend.position.inside = c(.65,.83),legend.direction = "horizontal",
                             legend.background = element_rect(fill = "transparent"),text = element_text(size = 12),
                             axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
                       labs(title = NULL,color = "",fill = "") +
                       # ylab(NULL) +
                       xlab(NULL),
                     figq2h2_res_webl +
                       theme(legend.position = "inside",legend.position.inside = c(.65,.83),legend.direction = "horizontal",
                             legend.background = element_rect(fill = "transparent"),text = element_text(size = 12),
                             axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
                       labs(title = NULL,color = "",fill = "") +
                       xlab(NULL),
                     figq2h3_webl +
                       theme(legend.position = "inside",legend.position.inside = c(.65,.83),legend.direction = "horizontal",
                             legend.background = element_rect(fill = "transparent"),text = element_text(size = 12)) +
                       labs(title = NULL,color = "",fill = "") #+
                       # ylab(NULL)
                       ,
                     ncol = 1,
labels = c("(a): Prior day maximum temperature","(b): Female body condition","(c): Female baseline corticosterone")))


ggsave("figures/ch3_fig3_provisbytempxmod_webl.png",p_full,width = 4,height = 8)

#### Figure 4: TRES Provis by prior day max temp, condition, baseline cort (no interaction)

ggplot_build(figq2h1_maxt_prior_tres)$layout$panel_scales_y
ggplot_build(figq2h2_res_tres)$layout$panel_scales_y
ggplot_build(figq2h3_tres)$layout$panel_scales_y

(p_full <- ggarrange(figq2h1_maxt_prior_tres +
                       theme(legend.position = "none",legend.position.inside = c(.65,.83),legend.direction = "horizontal",
                             legend.background = element_rect(fill = "transparent"),text = element_text(size = 12)) +
                       labs(title = NULL,color = "",fill = "") #+
                       #ylab(NULL)
                       ,
                     figq2h2_res_tres +
                       theme(legend.position = "none",legend.position.inside = c(.65,.83),legend.direction = "horizontal",
                             legend.background = element_rect(fill = "transparent"),text = element_text(size = 12)) +
                       labs(title = NULL,color = "",fill = ""),
                     figq2h3_tres +
                       theme(legend.position = "inside",legend.position.inside = c(.65,.83),legend.direction = "horizontal",
                             legend.background = element_rect(fill = "transparent"),text = element_text(size = 12)) +
                       labs(title = NULL,color = "",fill = "")# +
                       #ylab(NULL)
                       ,
                     ncol = 1,
                     labels = c("(a): Prior day maximum temperature","(b): Female body condition","(c): Female baseline corticosterone")))

ggsave("figures/ch3_fig4_provisbytempxmod_tres.png",p_full,width = 4,height = 8)

# #### OG figures 3, 4, 5 (when WEBL and TRES were in the same graph)
#
# ggplot_build(figq2h1_maxt_prior_webl)$layout$panel_scales_y
# ggplot_build(figq2h1_maxt_prior_tres)$layout$panel_scales_y
# (p_full <- ggarrange(figq2h1_maxt_prior_webl + labs(title = NULL) + theme(legend.position = "none",text = element_text(size = 12),axis.title.x = element_blank()) + ylim(0,50),figq2h1_maxt_prior_tres + labs(title = NULL) + ylim(0,50) + theme(legend.position = "none",axis.text.y = element_blank(),
#                                                                                                                                                                                                                                                  axis.ticks.y = element_blank(),
#                                                                                                                                                                                                                                                  text = element_text(size = 12),
#                                                                                                                                                                                                                                                  axis.title.y = element_blank(),
#                                                                                                                                                                                                                                                  axis.title.x = element_text(hjust = -.5)
# ),ncol = 2,
# labels = c("(a): Western Bluebird","(b): Tree Swallow")))
#
# ggsave("figures/ch3_fig3_provisbytemppriorday.png",p_full,width = 6.25,height = 4)
#
# ggplot_build(figq2h2_res_webl)$layout$panel_scales_y
# ggplot_build(figq2h2_res_tres)$layout$panel_scales_y
# (p_full <- ggarrange(figq2h2_res_webl + labs(title = NULL) + theme(legend.position = "none",text = element_text(size = 12),axis.title.x = element_blank()) + ylim(0,100),figq2h2_res_tres + ylim(0,100) + labs(title = NULL) + theme(legend.position = "none",axis.text.y = element_blank(),
#                                                                                                                                                                                                                                      axis.ticks.y = element_blank(),
#                                                                                                                                                                                                                                      text = element_text(size = 12),
#                                                                                                                                                                                                                                      axis.title.y = element_blank(),
#                                                                                                                                                                                                                                      axis.title.x = element_text(hjust = -.5)
# ),ncol = 2,
# labels = c("(a): Western Bluebird","(b): Tree Swallow")))
#
# ggsave("figures/ch3_fig4_provisbytempcond.png",p_full,width = 6.25,height = 4)

# ggplot_build(figq2h3_webl)$layout$panel_scales_y
# ggplot_build(figq2h3_tres)$layout$panel_scales_y
# (p_full <- ggarrange(figq2h3_webl + labs(title = NULL) + theme(legend.position = "none",text = element_text(size = 12),axis.title.x = element_blank()) + ylim(0,100),figq2h3_tres + ylim(0,100) + labs(title = NULL) + theme(legend.position = "none",axis.text.y = element_blank(),
#                                                                                                                                                                                                                    axis.ticks.y = element_blank(),
#                                                                                                                                                                                                                    text = element_text(size = 12),
#                                                                                                                                                                                                                    axis.title.y = element_blank(),
#                                                                                                                                                                                                                    axis.title.x = element_text(hjust = -.5)
# ),ncol = 2,
# labels = c("(a): Western Bluebird","(b): Tree Swallow")))
#
# ggsave("figures/ch3_fig5_provisbytemps1.png",p_full,width = 6.25,height = 4)

#### Figure S1: Provis by stress-induced cort

ggplot_build(figq4h2_webl)$layout$panel_scales_y
ggplot_build(figq4h2_tres)$layout$panel_scales_y
(p_full <- ggarrange(figq4h2_webl + theme(legend.position = "none",text = element_text(size = 12),axis.title.x = element_blank()) + labs(title = NULL) + ylim(-.5,110),figq4h2_tres + ylim(-.5,110) + labs(title = NULL) + theme(legend.position = "none",axis.text.y = element_blank(),
                                                                                                                                                                                                           axis.ticks.y = element_blank(),
                                                                                                                                                                                                           text = element_text(size = 12),
                                                                                                                                                                                                           axis.title.y = element_blank(),
                                                                                                                                                                                                           axis.title.x = element_text(hjust = 2.5)
),ncol = 2,
labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/ch3_figs1_provisbytemps2.png",p_full,width = 6.25,height = 4)

#### Figure S2: condition by prior week, s1 by prior day, s2 by prior week

ggplot_build(figq3h1_webl)$layout$panel_scales_y
ggplot_build(figq3h1_tres)$layout$panel_scales_y
ggplot_build(figq3h2_webl)$layout$panel_scales_y
ggplot_build(figq3h2_tres)$layout$panel_scales_y
ggplot_build(figq3h3_webl)$layout$panel_scales_y
ggplot_build(figq3h3_tres)$layout$panel_scales_y


(p_full <- ggarrange(figq3h1_webl + labs(title = NULL) + ylim(-16,12) + theme(axis.title.x = element_blank(),
                                                                                text = element_text(size = 12)),
                     figq3h1_tres + labs(title = NULL) + ylim(-16,12) + theme(axis.text.y = element_blank(),
                                                                                axis.ticks.y = element_blank(),
                                                                                text = element_text(size = 12),
                                                                                axis.title.y = element_blank(),
                                                                                axis.title.x = element_text(hjust = 10)),
                     figq3h2_webl + labs(title = NULL) + ylim(0,16.3) + theme(axis.title.x = element_blank(),
                                                                                  text = element_text(size = 12)),
                     figq3h2_tres + labs(title = NULL) + ylim(0,16.3) + theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  text = element_text(size = 12),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.title.x = element_text(hjust = -1)),
                     figq3h3_webl + labs(title = NULL) + ylim(0,93.8) + theme(axis.title.x = element_blank(),
                                                                                          text = element_text(size = 12)),
                     figq3h3_tres + labs(title = NULL) + ylim(0,93.8) + theme(axis.text.y = element_blank(),
                                                                                          axis.ticks.y = element_blank(),
                                                                                          text = element_text(size = 12),
                                                                                          axis.title.y = element_blank(),
                                                                                          axis.title.x = element_text(hjust = 10)),
                     ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow",
                                "(c)","(d)",
                                "(e)","(f)")
                     ))


ggsave("figures/ch3_figs2_condbytemp.png",p_full,width = 8,height = 10)

#### Figure 6: condition by prior day, s1 by prior week, s2 by prior day

ggplot_build(figq3h1_priorday_webl)$layout$panel_scales_y
ggplot_build(figq3h1_priorday_tres)$layout$panel_scales_y
ggplot_build(figq3h2_priorweek_webl)$layout$panel_scales_y
ggplot_build(figq3h2_priorweek_tres)$layout$panel_scales_y
ggplot_build(figq3h3_priorday_webl)$layout$panel_scales_y
ggplot_build(figq3h3_priorday_tres)$layout$panel_scales_y


(p_full <- ggarrange(figq3h1_priorday_webl + labs(title = NULL) + ylim(-10.6,12) + theme(axis.title.x = element_blank(),
                                                                                text = element_text(size = 12)),
                     figq3h1_priorday_tres + labs(title = NULL) + ylim(-10.6,12) + theme(axis.text.y = element_blank(),
                                                                                axis.ticks.y = element_blank(),
                                                                                text = element_text(size = 12),
                                                                                axis.title.y = element_blank(),
                                                                                axis.title.x = element_text(hjust = -1)),
                     figq3h2_priorweek_webl + labs(title = NULL) + ylim(0,30.5) + theme(axis.title.x = element_blank(),
                                                                                  text = element_text(size = 12)),
                     figq3h2_priorweek_tres + labs(title = NULL) + ylim(0,30.5) + theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  text = element_text(size = 12),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.title.x = element_text(hjust = 10)),
                     figq3h3_priorday_webl + labs(title = NULL) + ylim(2.59,78.9) + theme(axis.title.x = element_blank(),
                                                                                          text = element_text(size = 12)),
                     figq3h3_priorday_tres + labs(title = NULL) + ylim(2.59,78.9) + theme(axis.text.y = element_blank(),
                                                                                          axis.ticks.y = element_blank(),
                                                                                          text = element_text(size = 12),
                                                                                          axis.title.y = element_blank(),
                                                                                          axis.title.x = element_text(hjust = -1)),
                     ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow",
                                "(c)","(d)",
                                "(e)","(f)")
))


ggsave("figures/ch3_fig6_condbytemp.png",p_full,width = 8,height = 10)


#### AOS figures

## provis by mean temp

f1_aos <- figq1h1_webl + theme(text = element_text(size = 24)) + labs(title = NULL)

ggsave("figures/provisbymaxtempxhab_WEBL_q1h1_aosprez.png",plot =  f1_aos, width = 10, height = 6.6)

summary(m_webl_q1h1)

## Provis by temp * prior day t

ggplot_build(figq2h1_maxt_prior_webl)$layout$panel_scales_y
ggplot_build(figq2h2_res_webl)$layout$panel_scales_y
ggplot_build(figq2h3_webl)$layout$panel_scales_y

(full <- figq2h1_maxt_prior_webl +
  labs(title = NULL) +
  theme(legend.position = "inside",legend.position.inside = c(.27,.89),legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size = 24)))

ggsave("figures/ch3_figq2h1_maxt_prior_webl_full.png",plot = full, width = 10, height = 6.6)

(low <- figq2h1_maxt_prior_webl +
  labs(title = NULL) +
  scale_fill_manual(values = c("#0d0887","transparent"),guide = guide_legend(override.aes = list(alpha = 0) )) +
  scale_color_manual(values = c("#0d0887","transparent"),guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.position = "inside",legend.position.inside = c(.27,.89),legend.direction = "horizontal",
        text = element_text(size = 24,color = "transparent"),legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent")))

ggsave("figures/ch3_figq2h1_maxt_prior_webl_low.png",plot = low, width = 10, height = 6.6)

summary(m_webl_q2h1_maxt_prior)

## condition

(full <- figq2h2_res_webl +
  labs(title = NULL) +
  theme(legend.position = "inside",legend.position.inside = c(.7,.9),legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size = 24)))

ggsave("figures/ch3_figq2h2_res_webl_full.png",plot = full, width = 10, height = 6.6)

(low <- figq2h2_res_webl +
  labs(title = NULL) +
  scale_fill_manual(values = c("black","transparent"),guide = guide_legend(override.aes = list(alpha = 0) )) +
  scale_color_manual(values = c("black","transparent"),guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.position = "inside",legend.position.inside = c(.7,.83),legend.direction = "horizontal",
        text = element_text(size = 24,color = "transparent"),legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent")))

ggsave("figures/ch3_figq2h2_res_webl_low.png",plot = low, width = 10, height = 6.6)

summary(m_webl_q2h2_res)


## Baseline cort

(full <- figq2h3_webl +
  labs(title = NULL) +
  theme(legend.position = "inside",legend.position.inside = c(.55,.83),legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size = 24)))

ggsave("figures/ch3_figq2h3_webl_full.png",plot = full, width = 10, height = 6.6)

(high <- figq2h3_webl +
  labs(title = NULL) +
  scale_fill_manual(values = c("transparent","black"),guide = guide_legend(override.aes = list(alpha = 0) )) +
  scale_color_manual(values = c("transparent","black"),guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.position = "inside",legend.position.inside = c(.55,.83),legend.direction = "horizontal",
        text = element_text(size = 24,color = "transparent"),legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent")))

ggsave("figures/ch3_figq2h3_webl_high.png",plot = high, width = 10, height = 6.6)

summary(m_webl_q2h3)

