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
# modelsummary::get_gof(my_model)
# performance::model_performance(my_model)
wmean <- read_rds("data/wmean.rds")

g <- read_rds("data/growth_and_provis_mobilenetv3-original_dataset.h5.rds")
p2 <- read_rds("data/provis_with_attempt_1h_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year)) %>%
  mutate(model = "mobilenetv3",.before = start)
p1 <- read_rds("data/provis_with_attempt_1h.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year)) %>%
  mutate(model = "squeezenet",.before = start)

p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

p %>% dplyr::select(c(species, attempt_id)) %>% distinct() %>% summarize(.by = species,n = n())


## Sample size:


# p2 %>% dplyr::select(folder) %>% unique() %>% nrow()
#
#
#
# ## Objectives:
#
# # * Check videos associated with high provisioning rates
# # * Prepare for D-RUG meeting- simplify Danny meeting notes
#
# ## Check vids w/ high provis rates
#
#
# p1 %>% filter(valid_detections > 30) #%>% write_csv("../data/high_provis.csv")
# p2 %>% filter(valid_detections > 30)
#
#
#
#
# filter(p1,is.na(valid_detections))
# rbind(p1,p2) %>% arrange(nestbox,start) %>% filter(folder == "MBNG7_card9")
#
#
# What happens if I take the lower of the two valid_detections score?
#
#
# rbind(p1,p2) %>% filter(valid_detections == min(valid_detections,na.rm = TRUE) & row_number() == 1,.by = c(folder,start)) %>% filter(valid_detections > 40) %>% arrange(folder) %>% pull(species)
#
#
# I think I've solved it?? All of the high detections are in TRES boxes, where the nestlings do actually just like hang out at the edge of the box cuz they're hot.
#
# So the strategy is now, take the minimum valid_detections of these two models. They are generally very similar in terms of the number of detections, except in the cases where they have different tendencies to hallucinate. I think by combining the two models we've eliminated the very high detections where a huge series of videos is all tagged as bird. The ones left are where a TRES chick is hanging out in the box hole the whole time. I think that does mean that the models I built do not work for TRES so to do that we'd need to change our approach I think. Which is probably not worth it now.


# p <- rbind(p1,p2) %>% dplyr::filter(!is.na(valid_detections)) %>% group_by(folder,start) %>% dplyr::filter(valid_detections == min(valid_detections)) %>% dplyr::filter(row_number() == 1) %>% ungroup()




## Is interaction significant?

### WEBL


m <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_linint <- glmmTMB(valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_noint <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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
gtsave(int_tab_provis_webl,"figures/int_tab_provis_webl.html")

summary(m_linint)
check_collinearity(m_noint)

provis_webl <- m_linint

### Sample size


# ss_year_webl_provis <- m_linint$frame %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2021`,`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2021`+`2022`+`2023`)
# ss_year_webl_provis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_webl_provis.html")


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


t_ss_year_webl_provis %>% gtsave("figures/ss_year_webl_provis.html")

samp <- m_linint$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt() %>% gtsave("figures/ss_webl_provis.html")


dat_text_webl <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

## num nestboxes
m_linint$frame %>% pull(attempt_id) %>% unique() %>% length()


#
# (pl <- predict_response(m,terms = c("mean_temp_scaled [all]","habitat"),bias_correction = TRUE) %>%
#    plot(line_size = 1.5,alpha = .2,show_data = TRUE) +
#    facet_wrap(~ group,ncol = 2) +
#    # theme_classic() +
#     xlab("Average hourly temp") +
#     ylab("Nest visits per hour") +
#     labs(title = NULL,color = "Cover type") +
#    #ggtitle("Predicted growth across habitat and mean max daily temp") +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
#     theme(text = element_text(size = 16)) +
#     theme_classic() +
#     xlab("Mean temp over 60 mins") +
#     ylab("Visits/60 mins") +
#     labs(title = "WEBL") +
#    #ylim(0,80) +
#    annotate("text", x = -1.5, y = 50, size = 5, label = paste("N =",nobs(m)))
#    )
# ggsave("../figures/provisbytempxhab_WEBL.png",plot =  pl, width = 10, height = 6.6)
#




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
    scale_linetype_manual(values = c("Forest" = "solid","Orchard" = "dashed","Grassland" = "solid","Row crop" = "dotted")) +
    #ylim(-5,5) +
    geom_text(data = dat_text_webl, mapping = aes(x = -Inf, y = Inf,label = label),hjust = -.2, vjust = 1.2,inherit.aes = FALSE) +
    theme(legend.position = "none")
)

ggsave("figures/provisbymaxtempxhab_WEBL.png",plot =  fig5_webl, width = 10, height = 6.6)



(weblprovistrend <- emtrends(m_linint,specs = ~ habitat, var = c("mean_temp_scaled"),max.degree = 1) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    #group_by(degree) %>%
    gt())


gtsave(weblprovistrend,"figures/weblprovistrend.html")

data = data = dplyr::filter(p,!is.na(mean_temp),
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

# (t <- emmeans(m_linint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(mean_temp_scaled = (mean_temp_scaled * sd(data$mean_temp)) + mean(data$mean_temp),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            mean_temp_scaled = round(mean_temp_scaled),
#            mean_temp_scaled = paste0(mean_temp_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`Temperature` = "mean_temp_scaled",`Predicted provisioning` = "response",`2.5%` = "asymp.LCL",`97.5%` = "asymp.UCL" ) %>%
#     group_by(`Temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`Temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`Temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`Maximum TA_P` = if_else(`Maximum TA_P` == 0.000,"<0.001",as.character(`Maximum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"figures/weblprovisdelta.html")

# ((emmeans(m_linint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(m_linint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(m_linint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(m_linint,specs = pairwise ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(m_linint,formula = habitat ~ mean_temp_scaled, at = list(mean_temp_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()



## Emmeans to check for effect of habitat


(provisbyhabitat_webl <- emmeans(m_linint,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(provisbyhabitat_webl,"figures/provisbyhabitat_webl.html")


### Check for effect of temperature


(provisbyhabitat_summary_webl <- summary(m_noint) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(provisbyhabitat_summary_webl,"figures/provisbyhabitat_summary_webl.html")



### TRES


m <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) * habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_linint <- glmmTMB(valid_detections ~ mean_temp_scaled * habitat + poly(mean_temp_scaled,2) + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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

m_noint <- glmmTMB(valid_detections ~ poly(mean_temp_scaled,2) + habitat + julian_date_scaled + poly(tod_h_scaled,2) + mean_nestling_age_scaled + (1|attempt_id) + year_fct,
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
gtsave(int_tab_provis_tres,"figures/int_tab_provis_tres.html")

summary(m)
check_collinearity(m)


provis_tres <- m
### Sample size


# ss_year_tres_provis <- m$frame %>%
#   group_by(habitat,year_fct) %>% summarize(count = n()) %>% pivot_wider(values_from = count,names_from = year_fct) %>% as.tibble() %>% rename(Habitat = 'habitat') %>%
#   relocate(`2021`,.before = `2022`) %>%
#   ungroup() %>%
#   mutate(across(c(`2021`,`2022`,`2023`),~ replace_na(.x,0)),
#          Total = `2021`+`2022`+`2023`)
# ss_year_tres_provis %>% gt() %>%
#   grand_summary_rows(columns = -c(Habitat),fns = list(id = "Total") ~ sum(.,na.rm = TRUE)) %>%
#   gtsave("figures/ss_year_tres_provis.html")



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


t_ss_year_tres_provis %>% gtsave("figures/ss_year_tres_provis.html")

samp <- m$frame %>% group_by(habitat) %>% summarize(count = n())
samp %>% gt() %>% gtsave("figures/ss_tres_provis.html")


dat_text_tres <- data.frame(
  label = paste("N =",samp$count),
  group   = factor(c("Forest","Orchard","Grassland","Row crop"))
)

## num nestboxes
m$frame %>% pull(attempt_id) %>% unique() %>% length()





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

ggsave("figures/provisbymeantempxhab_TRES.png",plot =  fig5_tres, width = 10, height = 6.6)



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

ggsave("figures/fig5_provis_by_temp_hab.png",p_full,width = 6.25,height = 4)


## Emmeans to check for effect of habitat


(provisbyhabitat_tres <- emmeans(m,"habitat") %>% regrid() %>% pairs() %>% as_tibble() %>%
    mutate(across(estimate:z.ratio,~round(.x,digits = 2)),
           across(p.value,~round(.x,digits = 3))) %>% gt())
gtsave(provisbyhabitat_tres,"figures/provisbyhabitat_tres.html")


### Check for effect of temperature


(provisbyhabitat_summary_tres <- summary(m_noint) %>% coef() %>% pluck("cond") %>% as_tibble(rownames = "Covariate") %>%
    # dplyr::filter(Covariate != "poly(mean_temp_scaled, 2)1") %>%
    mutate(across(Estimate:`z value`,~round(.x,digits = 2)),
           across(`Pr(>|z|)`,~round(.x,digits = 3))) %>% gt())
gtsave(provisbyhabitat_summary_tres,"figures/provisbyhabitat_summary_tres.html")




(tresprovistrend <- emtrends(m,specs = ~ degree | habitat, var = "mean_temp_scaled",max.degree = 2) %>% test() %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)),
           df = round(df),
           p.value = if_else(p.value == 0.000,"<0.001",as.character(p.value))) %>%
    rename(Habitat = "habitat", `Max temp trend` = "mean_temp_scaled.trend",Df = "df", `Z-ratio` = "z.ratio", P = "p.value") %>%
    group_by(Habitat) %>%
    gt())


gtsave(tresprovistrend,"figures/tresprovistrend.html")

data = dplyr::filter(p,!is.na(mean_temp),
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

# (t <- emmeans(m_noint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2,0,2)),type = "response") %>% as.tibble() %>% #gt() %>%
#     mutate(mean_temp_scaled = (mean_temp_scaled * sd(data$mean_temp)) + mean(data$mean_temp),
#            across(where(is.numeric), ~ round(.x, digits = 2)),
#            mean_temp_scaled = round(mean_temp_scaled),
#            mean_temp_scaled = paste0(mean_temp_scaled,"\u00b0C")) %>%
#     #tibble() %>%
#     dplyr::select(-df) %>%
#     #mutate(Model = rep(c("TA2 * LU + TA * LU", "TA2 + TA * LU", "TA2 + TA + LU"),2), type = c(rep("Maximum TA", 3),rep("Minimum TA",3)), .before = AIC) %>%
#     rename(Habitat = "habitat",`Temperature` = "mean_temp_scaled",`Predicted provisioning` = "response",`2.5%` = "asymp.LCL",`97.5%` = "asymp.UCL" ) %>%
#     group_by(`Temperature`) %>%
#     mutate(row=row_number()) %>%
#     pivot_longer(-c(`Temperature`, row,Habitat)) %>%
#     pivot_wider(names_from=c(`Temperature`, name), values_from=value) %>%
#     dplyr::select(-row) %>%
#     #mutate(`Maximum TA_P` = if_else(`Maximum TA_P` == 0.000,"<0.001",as.character(`Maximum TA_P`))) %>%
#     #mutate(`Minimum TA_P` = if_else(`Minimum TA_P` == 0.000,"<0.001",as.character(`Minimum TA_P`))) %>%
#     gt() %>% tab_options(data_row.padding = px(1)) %>%
#     tab_spanner_delim(
#       delim="_"
#     ))
#
# gtsave(t,"figures/tresprovisdelta.html")

# ((emmeans(m_noint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(2)),type = "response") %>% as.tibble() %>% pull(response))-(emmeans(m_noint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response)))/(emmeans(m_noint,specs = ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2)),type = "response") %>% as.tibble() %>% pull(response))
#
# emmeans(m_noint,specs = pairwise ~ habitat,by = c("mean_temp_scaled"), at = list(mean_temp_scaled = c(-2,0,2))) %>% plot(comparisons = TRUE)
# emmip(m_noint,formula = habitat ~ mean_temp_scaled, at = list(mean_temp_scaled = seq(from = -2.5, to = 2.5, by = .1)),CIs = TRUE, plotit = FALSE) %>% emmip_ggplot() + theme_classic()


