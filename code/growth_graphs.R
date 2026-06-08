# growth figures for publication
#
library(tidyverse)
library(ggplot2)
library(egg)
library(ggeffects)

source("code/growth_model_eval.R")

## Combined WEBL and TRES growth plots meanmaxtempI


ggplot_build(fig2_webl)$layout$panel_scales_y
ggplot_build(fig2_tres)$layout$panel_scales_y
(p_full <- ggarrange(fig2_webl + theme(text = element_text(size = 12),axis.title.x = element_blank()),fig2_tres + ylim(-1.09,3.42) + theme(axis.text.y = element_blank(),
                                                                                                                                           axis.ticks.y = element_blank(),
                                                                                                                                           text = element_text(size = 12),
                                                                                                                                           axis.title.y = element_blank(),
                                                                                                                                           axis.title.x = element_text(hjust = 2.8)),ncol = 2,
                     labels = c("(a): Western Bluebird","(b): Tree Swallow")))

ggsave("figures/fig2_growth_by_temp_hab.png",p_full,width = 6.25,height = 4)

## Combined WEBL and TRES growth plots maxhi_week

# ggplot_build(fig2_maxhiweek_webl)$layout$panel_scales_y
# ggplot_build(fig2_maxhiweek_tres)$layout$panel_scales_y
# (p_full <- ggarrange(fig2_maxhiweek_webl + theme(text = element_text(size = 12),axis.title.x = element_blank()),fig2_maxhiweek_tres + ylim(-.767,3.65) + theme(axis.text.y = element_blank(),
#                                                                                                                                            axis.ticks.y = element_blank(),
#                                                                                                                                            text = element_text(size = 12),
#                                                                                                                                            axis.title.y = element_blank(),
#                                                                                                                                            axis.title.x = element_text(hjust = 2)),ncol = 2,
#                      labels = c("(a): Western Bluebird","(b): Tree Swallow")))
# ggsave("figures/fig2_growth_by_maxhiweek_hab.png",p_full,width = 6.25,height = 4)

## Combined WEBL and TRES growth plots maxhi_day

# ggplot_build(fig2_maxhiday_webl)$layout$panel_scales_y
# ggplot_build(fig2_maxhiday_tres)$layout$panel_scales_y
# (p_full <- ggarrange(fig2_maxhiday_webl + theme(text = element_text(size = 12),axis.title.x = element_blank()),fig2_maxhiday_tres + ylim(-.767,3.76) + theme(axis.text.y = element_blank(),
#                                                                                                                                                                axis.ticks.y = element_blank(),
#                                                                                                                                                                text = element_text(size = 12),
#                                                                                                                                                                axis.title.y = element_blank(),
#                                                                                                                                                                axis.title.x = element_text(hjust = 2)),ncol = 2,
#                      labels = c("(a): Western Bluebird","(b): Tree Swallow")))
# ggsave("figures/fig2_growth_by_maxhiday_hab.png",p_full,width = 6.25,height = 4)

## Combined WEBL and TRES growth plots deghr30week

# ggplot_build(fig2_deghr30week_webl)$layout$panel_scales_y
# ggplot_build(fig2_deghr30week_tres)$layout$panel_scales_y
# (p_full <- ggarrange(fig2_deghr30week_webl + theme(text = element_text(size = 12),axis.title.x = element_blank()),fig2_deghr30week_tres + ylim(-.767,3.46) + theme(axis.text.y = element_blank(),
#                                                                                                                                                              axis.ticks.y = element_blank(),
#                                                                                                                                                              text = element_text(size = 12),
#                                                                                                                                                              axis.title.y = element_blank(),
#                                                                                                                                                              axis.title.x = element_text(hjust = 2)),ncol = 2,
#                      labels = c("(a): Western Bluebird","(b): Tree Swallow")))
# ggsave("figures/fig2_growth_by_deghr30week_hab.png",p_full,width = 6.25,height = 4)

## Combined WEBL and TRES growth plots hihr25week

# ggplot_build(fig2_hihr25week_webl)$layout$panel_scales_y
# ggplot_build(fig2_hihr25week_tres)$layout$panel_scales_y
# (p_full <- ggarrange(fig2_hihr25week_webl + theme(text = element_text(size = 12),axis.title.x = element_blank()),fig2_hihr25week_tres + ylim(-.768,3.27) + theme(axis.text.y = element_blank(),
#                                                                                                                                                                    axis.ticks.y = element_blank(),
#                                                                                                                                                                    text = element_text(size = 12),
#                                                                                                                                                                    axis.title.y = element_blank(),
#                                                                                                                                                                    axis.title.x = element_text(hjust = 2)),ncol = 2,
#                      labels = c("(a): Western Bluebird","(b): Tree Swallow")))
# ggsave("figures/fig2_growth_by_hihr25week_hab.png",p_full,width = 6.25,height = 4)

