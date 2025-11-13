# Cort figures for publication
#
library(tidyverse)
library(ggplot2)
library(egg)
library(ggeffects)

source("code/cort_model_eval.R")


## combined figs2 (baseline cort) prior week max temp


library(egg)
ggplot_build(figs2_webl)$layout$panel_scales_y
ggplot_build(figs2_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_webl + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                            axis.title.x = element_blank()),
                     figs2_tres + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                            axis.text.y = element_blank(),
                                                                                            axis.title.y = element_blank(),
                                                                                            axis.title.x = element_text(hjust = 2.6)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_temp_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior week max temp


library(egg)

(p_full <- ggarrange(fig_webl_s2 + ylim(1.3,139) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                            axis.title.x = element_blank()),
                     fig_tres_s2 + ylim(1.3,139) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                            axis.text.y = element_blank(),
                                                                                            axis.title.y = element_blank(),
                                                                                            axis.title.x = element_text(hjust = 2.6)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab.png",p_full,width = 6.25,height = 4)

## combined fig4 prior week max temp

library(egg)
ggplot_build(fig4_tres)$layout$panel_scales_y
ggplot_build(fig4_webl)$layout$panel_scales_y
ggplot_build(fig_webl_s2)$layout$panel_scales_y
ggplot_build(fig_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl + ylim(0.0011,106) + ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                                                                          axis.title.x = element_blank()),fig4_tres + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                                                                                                                                                            axis.text.y = element_blank(),
                                                                                                                                                                                                                                            axis.title.y = element_blank(),
                                                                                                                                                                                                                                            axis.title.x = element_text(hjust = 2.)),
                     fig_webl_s2 + ylim(1.3,139) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                         axis.title.x = element_blank()),
                     fig_tres_s2 + ylim(1.3,139) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                         axis.text.y = element_blank(),
                                                                                         axis.title.y = element_blank(),
                                                                                         axis.title.x = element_text(hjust = 2.6)),
                     ncol = 2,
                     labels = c("(a)   Western Bluebird","(b)   Tree Swallow","(c)","(d)")))

ggsave("figures/fig4_abscort_by_temp_hab.png",p_full,width = 6.25,height = 8)

## combined figure for figs2 prior day temp


library(egg)
ggplot_build(figs2_priordayt_webl)$layout$panel_scales_y
ggplot_build(figs2_priordayt_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_priordayt_webl + ylim(0.0812,17.9) + ylab("Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                                                                       axis.title.x = element_blank()),
                     figs2_priordayt_tres +
                       ylim(0.0812,17.9) +
                       labs(title = element_blank()) +
                       theme(text = element_text(size = 12),
                             axis.ticks.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.title.x = element_text(hjust = -2)),
                     ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_priordayt_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior week max temp


library(egg)
ggplot_build(fig_priordayt_webl_s2)$layout$panel_scales_y
ggplot_build(fig_priordayt_tres_S2)$layout$panel_scales_y
(p_full <- ggarrange(fig_priordayt_webl_s2 + ylim(1.3,133) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                         axis.title.x = element_blank()),
                     fig_priordayt_tres_S2 + ylim(1.3,133) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                         axis.text.y = element_blank(),
                                                                                         axis.title.y = element_blank(),
                                                                                         axis.title.x = element_text(hjust = -1)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_priordayt.png",p_full,width = 6.25,height = 4)

## combined figure for fig4 prior day temp


library(egg)
ggplot_build(fig4_webl_priordayt)$layout$panel_scales_y
ggplot_build(fig4_tres_priordayt)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_priordayt + ylim(0.0065,125) + ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                                                                       axis.title.x = element_blank()),
                     fig4_tres_priordayt +
                       ylim(0.0065,125) +
                       labs(title = element_blank()) +
                       theme(text = element_text(size = 12),
                             axis.ticks.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.title.x = element_text(hjust = -2)),
                     ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_priordayt_hab.png",p_full,width = 6.25,height = 4)

## combined figure for figs2 prior week max heat index

library(egg)
ggplot_build(figs2_weekhi_webl)$layout$panel_scales_y
ggplot_build(figs2_weekhi_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_weekhi_webl + ylim(0.0812,17.9) + ylab("Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                                                                    axis.title.x = element_blank()),
                     figs2_weekhi_tres +
                       ylim(0.0812,17.9) +
                       labs(title = element_blank()) +
                       theme(text = element_text(size = 12),
                             axis.ticks.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.title.x = element_text(hjust = -20)),
                     ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_priorweekmaxhhi_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior week max heat index


library(egg)
ggplot_build(fig_weekhi_webl_s2)$layout$panel_scales_y
ggplot_build(fig_weekhi_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig_weekhi_webl_s2 + ylim(1.3,112) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                   axis.title.x = element_blank()),
                     fig_weekhi_tres_s2 + ylim(1.3,112) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                   axis.text.y = element_blank(),
                                                                                                   axis.title.y = element_blank(),
                                                                                                   axis.title.x = element_text(hjust = -1)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_weekhi.png",p_full,width = 6.25,height = 4)

## combined figure for fig4 prior week max heat index

library(egg)
ggplot_build(fig4_webl_weekhi)$layout$panel_scales_y
ggplot_build(fig4_tres_weekhi)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_weekhi + ylim(0.0065,300) + ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                                                                    axis.title.x = element_blank()),
                     fig4_tres_weekhi +
                       ylim(0.0065,300) +
                       labs(title = element_blank()) +
                       theme(text = element_text(size = 12),
                             axis.ticks.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.title.x = element_text(hjust = -25)),
                     ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_priorweekmaxhhi_hab.png",p_full,width = 6.25,height = 4)

#combined figs2 for prior day heat index

library(egg)
ggplot_build(figs2_priordaymaxhhi_webl)$layout$panel_scales_y
ggplot_build(figs2_priordaymaxhhi_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_priordaymaxhhi_webl + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                           axis.title.x = element_blank()),
                     figs2_priordaymaxhhi_tres + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                           axis.text.y = element_blank(),
                                                                                                           axis.title.y = element_blank(),
                                                                                                           axis.title.x = element_text(hjust = -5)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_priordaymaxhhi_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior week max heat index


library(egg)
ggplot_build(fig_priordaymaxhhi_webl_s2)$layout$panel_scales_y
ggplot_build(fig_priordaymaxhhi_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig_priordaymaxhhi_webl_s2 + ylim(.0173,181) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                axis.title.x = element_blank()),
                     fig_priordaymaxhhi_tres_s2 + ylim(.0173,181) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                axis.text.y = element_blank(),
                                                                                                axis.title.y = element_blank(),
                                                                                                axis.title.x = element_text(hjust = -1)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_priordaymaxhhi.png",p_full,width = 6.25,height = 4)


#combined fig4 for prior day heat index

library(egg)
ggplot_build(fig4_webl_priordaymaxhhi)$layout$panel_scales_y
ggplot_build(fig4_tres_priordaymaxhhi)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_priordaymaxhhi + ylim(0.0065,168) + ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                           axis.title.x = element_blank()),
                     fig4_tres_priordaymaxhhi + ylim(0.0065,168) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                           axis.text.y = element_blank(),
                                                                                                           axis.title.y = element_blank(),
                                                                                                           axis.title.x = element_text(hjust = -5)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_priordaymaxhhi_hab.png",p_full,width = 6.25,height = 4)

#combined figs2 for prior day cumulative heat index-hours over 25C

library(egg)
ggplot_build(figs2_cumhiday_webl)$layout$panel_scales_y
ggplot_build(figs2_cumhiday_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_cumhiday_webl + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                           axis.title.x = element_blank()),
                     figs2_cumhiday_tres + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                           axis.text.y = element_blank(),
                                                                                                           axis.title.y = element_blank(),
                                                                                                           axis.title.x = element_text(hjust = 2.2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_priordaycumhi_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior day cumulative heat index-hours over 25C


library(egg)
ggplot_build(fig_cumhiday_webl_s2)$layout$panel_scales_y
ggplot_build(fig_cumhiday_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig_cumhiday_webl_s2 + ylim(1.3,118) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                          axis.title.x = element_blank()),
                     fig_cumhiday_tres_s2 + ylim(1.3,118) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                          axis.text.y = element_blank(),
                                                                                                          axis.title.y = element_blank(),
                                                                                                          axis.title.x = element_text(hjust = 2.6)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_cumhiday.png",p_full,width = 6.25,height = 4)

#combined fig4 for prior day cumulative heat index-hours over 25C

library(egg)
ggplot_build(fig4_webl_cumhiday)$layout$panel_scales_y
ggplot_build(fig4_tres_cumhiday)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_cumhiday + ylim(0.0065,106)+ ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                     axis.title.x = element_blank()),
                     fig4_tres_cumhiday + ylim(0.0065,106) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                     axis.text.y = element_blank(),
                                                                                                     axis.title.y = element_blank(),
                                                                                                     axis.title.x = element_text(hjust = 2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_priordaycumhi_hab.png",p_full,width = 6.25,height = 4)

#combined figs2 for prior week cumulative heat index-hours over 25C

library(egg)
ggplot_build(figs2_cumhiweek_webl)$layout$panel_scales_y
ggplot_build(figs2_cumhiweek_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_cumhiweek_webl + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                     axis.title.x = element_blank()),
                     figs2_cumhiweek_tres + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                     axis.text.y = element_blank(),
                                                                                                     axis.title.y = element_blank(),
                                                                                                     axis.title.x = element_text(hjust = -10)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_priorweekcumhi_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior week cumulative heat index-hours over 25C


library(egg)
ggplot_build(fig_cumhiweek_webl_s2)$layout$panel_scales_y
ggplot_build(fig_cumhiweek_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig_cumhiweek_webl_s2 + ylim(.84,123) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                  axis.title.x = element_blank()),
                     fig_cumhiweek_tres_s2 + ylim(.84,123) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                  axis.text.y = element_blank(),
                                                                                                  axis.title.y = element_blank(),
                                                                                                  axis.title.x = element_text(hjust = 2.6)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_cumhiweek.png",p_full,width = 6.25,height = 4)

#combined fig4 for prior week cumulative heat index-hours over 25C

library(egg)
ggplot_build(fig4_webl_cumhiweek)$layout$panel_scales_y
ggplot_build(fig4_tres_cumhiweek)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_cumhiweek + ylim(0.0065,141)+ ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                   axis.title.x = element_blank()),
                     fig4_tres_cumhiweek + ylim(0.0065,141) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                   axis.text.y = element_blank(),
                                                                                                   axis.title.y = element_blank(),
                                                                                                   axis.title.x = element_text(hjust = 2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_priorweekcumhi_hab.png",p_full,width = 6.25,height = 4)

#combined figs2 for prior day cumulative degree-hours over 30C

library(egg)
ggplot_build(figs2_cumdegreeday_webl)$layout$panel_scales_y
ggplot_build(figs2_cumdegreeday_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_cumdegreeday_webl + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                      axis.title.x = element_blank()),
                     figs2_cumdegreeday_tres + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                      axis.text.y = element_blank(),
                                                                                                      axis.title.y = element_blank(),
                                                                                                      axis.title.x = element_text(hjust = 2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_cumdegreeday_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior day cumulative degree-hours over 30C


library(egg)
ggplot_build(fig_cumdegreeday_webl_s2)$layout$panel_scales_y
ggplot_build(fig_cumdegreeday_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig_cumdegreeday_webl_s2 + ylim(.0812,128) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                   axis.title.x = element_blank()),
                     fig_cumdegreeday_tres_s2 + ylim(.0812,128) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                   axis.text.y = element_blank(),
                                                                                                   axis.title.y = element_blank(),
                                                                                                   axis.title.x = element_text(hjust = 2.6)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_cumdegreeday.png",p_full,width = 6.25,height = 4)

#combined fig4 for prior day cumulative degree-hours over 30C

library(egg)
ggplot_build(fig4_webl_cumdegreeday)$layout$panel_scales_y
ggplot_build(fig4_tres_cumdegreeday)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_cumdegreeday + ylim(0.0065,130)+ ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                    axis.title.x = element_blank()),
                     fig4_tres_cumdegreeday + ylim(0.0065,130) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                    axis.text.y = element_blank(),
                                                                                                    axis.title.y = element_blank(),
                                                                                                    axis.title.x = element_text(hjust = 2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_cumdegreeday_hab.png",p_full,width = 6.25,height = 4)

#combined figs2 for prior week cumulative degree-hours over 30C

library(egg)
ggplot_build(figs2_cumdegreeweek_webl)$layout$panel_scales_y
ggplot_build(figs2_cumdegreeweek_tres)$layout$panel_scales_y
(p_full <- ggarrange(figs2_cumdegreeweek_webl + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                         axis.title.x = element_blank()),
                     figs2_cumdegreeweek_tres + ylim(0.0812,17.9) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                         axis.text.y = element_blank(),
                                                                                                         axis.title.y = element_blank(),
                                                                                                         axis.title.x = element_text(hjust = 2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/figs2_s1cort_by_cumdegreeweek_hab.png",p_full,width = 6.25,height = 4)

## combined cort s2 prior week cumulative degree-hours over 30C


library(egg)
ggplot_build(fig_cumdegreeweek_webl_s2)$layout$panel_scales_y
ggplot_build(fig_cumdegreeweek_tres_s2)$layout$panel_scales_y
(p_full <- ggarrange(fig_cumdegreeweek_webl_s2 + ylim(1.3,118) + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                        axis.title.x = element_blank()),
                     fig_cumdegreeweek_tres_s2 + ylim(1.3,118) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                        axis.text.y = element_blank(),
                                                                                                        axis.title.y = element_blank(),
                                                                                                        axis.title.x = element_text(hjust = 2.6)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig_s2cort_by_temp_hab_cumdegreeweek.png",p_full,width = 6.25,height = 4)

#combined fig4 for prior week cumulative degree-hours over 30C

library(egg)
ggplot_build(fig4_webl_cumdegreeweek)$layout$panel_scales_y
ggplot_build(fig4_tres_cumdegreeweek)$layout$panel_scales_y
(p_full <- ggarrange(fig4_webl_cumdegreeweek + ylim(0,106)+ ylab("Stress-induced - Baseline corticosterone (ng/\u00B5L)") + labs(title = element_blank()) + theme(text = element_text(size = 12),
                                                                                                                                                                      axis.title.x = element_blank()),
                     fig4_tres_cumdegreeweek + ylim(0,106) + labs(title = element_blank()) + theme(text = element_text(size = 12),axis.ticks.y = element_blank(),
                                                                                                       axis.text.y = element_blank(),
                                                                                                       axis.title.y = element_blank(),
                                                                                                       axis.title.x = element_text(hjust = 2)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig4_abscort_by_cumdegreeweek_hab.png",p_full,width = 6.25,height = 4)
