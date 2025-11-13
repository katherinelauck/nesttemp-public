## Create all objects needed for provis analysis
##

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
library(future)

instant_temp <- read_rds("../data/provis_manytempmeasures.rds")

b <- read_csv("data/banding-and-morphometrics_proofed.csv")


g <- read_rds("data/growth_cort_provis_manytempmeasures.rds")


## Is interaction significant?

### WEBL
#### maxhi of closest 30 min period to start time of hour


m <- glmmTMB(valid_detections ~ poly(hi_30min_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_linint <- glmmTMB(valid_detections ~ hi_30min_scaled * habitat + poly(hi_30min_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_noint <- glmmTMB(valid_detections ~ poly(hi_30min_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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
gtsave(int_tab_provis_webl_hi,"figures/int_tab_provis_webl_hi.html")

summary(m)
check_collinearity(m_noint)

provis_webl_hi <- m

### Sample size


ss_year_webl_provis_hi <- m$frame %>%
  group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
  relocate(`2021`,.before = `2022`) %>%
  ungroup() %>%
  mutate(across(c(`2021`,`2022`,`2023`),~ replace_na(.x,0)),
         Total = `2021`+`2022`+`2023`)
ss_year_webl_provis_hi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
  gtsave("figures/ss_year_webl_provis_hi.html")

samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt() %>% gtsave("figures/ss_webl_provis_hi.html")


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
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "solid","Grassland" = "solid","Row crop" = "dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_webl_hi, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/provisbymaxtempxhab_WEBL_hi.png",plot =  fig5_webl_hi, width = 10, height = 6.6)



(weblprovistrend_hi <- emtrends(m,specs = ~ degree | habitat, var = "hi_30min_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hi_30min_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


gtsave(weblprovistrend_hi,"figures/weblprovistrend_hi.html")

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
gtsave(provisbyhabitat_webl_hi,"figures/provisbyhabitat_webl_hi.html")


### Check for effect of temperature


(provisbyhabitat_summary_webl_hi <- summary(m) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(hi_30min_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(provisbyhabitat_summary_webl_hi,"figures/provisbyhabitat_summary_webl_hi.html")

### TRES
#### maxhi of closest 30 min period to start time of hour


m <- glmmTMB(valid_detections ~ poly(hi_30min_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

#emmeans(m,specs = pairwise ~ habitat,by = c("hi_30min_scaled"), at = list(hi_30min_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
#emmip(m,formula = habitat ~ hi_30min_scaled, at = list(hi_30min_scaled = seq(from = -2.5, to = 2.5, by = .1)))

m_linint <- glmmTMB(valid_detections ~ hi_30min_scaled * habitat + poly(hi_30min_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_noint <- glmmTMB(valid_detections ~ poly(hi_30min_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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
gtsave(int_tab_provis_tres_hi,"figures/int_tab_provis_tres_hi.html")

summary(m)
check_collinearity(m_noint)

provis_tres_hi <- m

### Sample size


ss_year_tres_provis_hi <- m$frame %>%
  group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
  relocate(`2021`,.before = `2022`) %>%
  ungroup() %>%
  mutate(across(c(`2021`,`2022`,`2023`),~ replace_na(.x,0)),
         Total = `2021`+`2022`+`2023`)
ss_year_tres_provis_hi %>% gt() %>%
  grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
  gtsave("figures/ss_year_tres_provis_hi.html")

samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt() %>% gtsave("figures/ss_tres_provis_hi.html")


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
                       # breaks = c(20,30,40,50),
                       # labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       # limits = c(18,5)
    ) +
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dotted","Grassland" = "solid","Row crop" = "dotted")) +
    ylim(0,60) +
    geom_text(data = dat_text_tres, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/provisbymaxtempxhab_TRES_hi.png",plot =  fig5_tres_hi, width = 10, height = 6.6)



(tresprovistrend_hi <- emtrends(m,specs = ~ degree | habitat, var = "hi_30min_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "hi_30min_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


gtsave(tresprovistrend_hi,"figures/tresprovistrend_hi.html")

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
gtsave(provisbyhabitat_tres_hi,"figures/provisbyhabitat_tres_hi.html")


### Check for effect of temperature


(provisbyhabitat_summary_tres_hi <- summary(m) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(hi_30min_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(provisbyhabitat_summary_tres_hi,"figures/provisbyhabitat_summary_tres_hi.html")


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

ggsave("figures/fig5_provis_by_hi_hab.png",p_full,width = 6.25,height = 4)
