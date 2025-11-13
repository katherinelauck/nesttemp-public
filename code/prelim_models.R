##### Preliminary modeling

library(tidyverse)
require(ggeffects)
require(lme4)
require(viridis)
library(car)

g <- read_rds("data/growth.rds")

summary(lm(meanh ~ habitat,data = g))
summary(lm(meanmaxtempI ~ habitat,data = g))

# Does habitat affect growth interactively with temp?

#### WEBL

model <- lmerTest::lmer(gweight_scaled ~ meanmaxtempI_scaled * habitat + meanmaxtempI_scaled_sq * habitat + meanmintempI_scaled + juliandate_scaled + site + year + (1|attempt_id),data = filter(g,Species == "WEBL",
                !is.na(meanmaxtempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

data <- filter(g,Species == "WEBL",
       !is.na(meanmaxtempI)) %>%
  mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"),
         meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled)

temp <- range(data$meanmaxtempI_scaled) # c(-2, 2) # define range of predictor in standard scale
temp <- seq(from = temp[1],to = temp[2], length.out = 200) # build dummy temp evenly spread across range

grid=expand.grid(meanmaxtempI_scaled = temp, habitat = c("Forest","Orchard","Grassland","Row Crop"),meanmintempI_scaled = 0, juliandate_scaled=0,gweight_scaled=0, year = factor("2022",levels = c("2021","2022","2023"),ordered = FALSE),site = factor("WI",levels = "WI","MB","PG","RR","D"),ordered = FALSE) # not sure how to make this go when model includes year and site

grid=cbind(grid[,1],grid[,1]^2,grid[,2:5])
colnames(grid)[1]<-"meanmaxtempI_scaled"

colnames(grid)[2]="meanmaxtempI_scaled_sq"

temp_predict <- function(model,grid) {
  mm=model.matrix(terms(model), grid)
  predicted = mm %*% fixef(model)
  pvar1 <- diag(mm %*% vcov(model) %*% t(mm))
  return(bind_cols(habitat = grid$habitat,temp = grid$meanmaxtempI_scaled, predicted = predicted, lower = predicted-2*sqrt(pvar1), upper = predicted+2*sqrt(pvar1))) # bind into dataframe & return
}

fig_data <- temp_predict(meanI_WEBL,grid) %>% mutate(facet = factor(habitat,c("Forest","Orchard","Grassland","Row Crop"),ordered = TRUE))

grid.quant_sd=expand.grid(Tmax_std_gridmet = c(-2,0,2),
                          NewLU1= c("Forest","Ag","Natural_open","Human"),
                          pcpbefore_raw_gridmet = 0,NLCD_p_forest= 0, NLCD_p_human=0, NLCD_p_ag=0,elevation=0,
                          substrate_binary=1,laydate_scaled=0,at_least_one_success=1)
grid.quant_sd=cbind(grid.quant_sd[,1],grid.quant_sd[,1]^2,grid.quant_sd[,2:10])
colnames(grid.quant_sd)[1]<-"Tmax_std_gridmet"

colnames(grid.quant_sd)[2]="Tmax_std_gridmet_sq"

quant_sd <- temp_predict(m,grid.quant_sd)

### which interactions are significant?

# ag <- read_rds("results/revisions/mainv1_withregion_ag.rds")
# summary(ag)
# forest <- read_rds("results/revisions/mainv1_withregion_forest.rds")
# summary(forest)
# human <- read_rds("results/revisions/mainv1_withregion_human.rds")
# summary(human)
# open <- read_rds("results/revisions/mainv1_withregion_open.rds")
# summary(open)

# f_labels <- data.frame(lu = factor(c("Forest", "Ag", "Natural_open", "Human")), label = c("**", "*", "**","*")) ## significance labels for each facet

p <- ggplot(data = fig_data, aes(temp, predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, color = habitat, fill = habitat), linetype = 2, alpha = .2) + # plot confidence intervals
  geom_line(aes(color = habitat), linewidth = 1) + # overlay line
  # facet_wrap(~habitat, nrow = 4, ncol = 1) + # facet by land use
  ylab("Daily growth (g/day)") +
  # xlab("Mean maximum temperature over nesting period (z-scaled)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values = c("#7CAE00","#F8766D","#00BFC4","#C77CFF")) +
  scale_fill_manual(values = c("#7CAE00","#F8766D","#00BFC4","#C77CFF"))
  # geom_text(x = -2, y = .25, aes(label = label), data = f_labels) +
  # ylim(0,.95)

hist <- ggplot(data = data, aes(x = Tmax_std_gridmet)) +
  geom_histogram() +
  xlab("Maximum temperature anomaly") +
  ylab("Frequency") +
  theme_classic()

plot <- ggarrange(p, hist, heights = c(4,1), ncol = 1, nrow = 2, align = "v")

meanI_WEBL_poly <- lmerTest::lmer(gweight_scaled ~ poly(meanmaxtempI_scaled,2,raw = TRUE) * habitat + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "WEBL",
                                                                                                                                                                                                      !is.na(meanmaxtempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

(p <- ggpredict(meanI_WEBL_poly,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                      alpha = .2,show_data = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13))
)
ggsave("figures/growthbymeanI_WEBL.png",p,width = 6,height = 4)



meanI_WEBL_linear <- lmerTest::lmer(gweight_scaled ~ poly(meanmaxtempI_scaled,1,raw = TRUE) * habitat + poly(meanmaxtempI_scaled,2,raw = TRUE) + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "WEBL",
                                                                                                                                                                                                 !is.na(meanmaxtempI)) %>%
                                    mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

(p <- ggpredict(meanI_WEBL_linear,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                           alpha = .2,show_data = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13))
)
ggsave("figures/growthbymeanI_linearT_WEBL.png",p,width = 6,height = 4)

meanI_WEBL_noint <- lmerTest::lmer(gweight_scaled ~ poly(meanmaxtempI_scaled,1,raw = TRUE) + habitat + poly(meanmaxtempI_scaled,2,raw = TRUE) + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "WEBL",
                                                                                                                                                                                                   !is.na(meanmaxtempI)) %>%
                                      mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                                    ~ scale(.x)[,1],
                                                    .names = "{.col}_scaled"),
                                             meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))


(p <- ggpredict(meanI_WEBL_noint,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                             alpha = .2,show_data = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13))
)
ggsave("figures/growthbymeanI_noint_WEBL.png",p,width = 6,height = 4)

# interaction significant?
anova(meanI_WEBL_poly,meanI_WEBL_linear,meanI_WEBL_noint)

# sq term necessary?

meanI_WEBL_nosq <- lmerTest::lmer(gweight_scaled ~ meanmaxtempI_scaled * habitat + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "WEBL",
                                                                                                                                                                                           !is.na(meanmaxtempI)) %>%
                                    mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

meanI_WEBL_notemp <- lmerTest::lmer(gweight_scaled ~ habitat + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "WEBL",
                                                                                                                                                                        !is.na(meanmaxtempI)) %>%
                                    mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))
anova(meanI_WEBL_poly,meanI_WEBL_nosq,meanI_WEBL_notemp)

anova(meanI_WEBL_nosq,meanI_WEBL_linear)

(p <- ggpredict(meanI_WEBL_nosq,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                            alpha = .2,show_data = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13))
)
ggsave("figures/growthbymeanI_nosq_WEBL.png",p,width = 6,height = 4)


#### TRES

meanI_TRES <- lmerTest::lmer(gweight_scaled ~ poly(meanmaxtempI_scaled,2,raw = TRUE) * habitat + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "TRES",
                                                                                                                                                                                            !is.na(meanmaxtempI)) %>%
                               filter(meanmaxtempI < 45) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

(p <- ggpredict(meanI_TRES,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                           alpha = .2,show_data = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13))
)
ggsave("figures/growthbymeanI_WEBL.png",p,width = 6,height = 4)

meanI_TRES_linear <- lmerTest::lmer(gweight_scaled ~ poly(meanmaxtempI_scaled,1,raw = TRUE) * habitat + poly(meanmaxtempI_scaled,2,raw = TRUE) + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "TRES",
                                                                                                                                                                                            !is.na(meanmaxtempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

meanI_TRES_noint <- lmerTest::lmer(gweight_scaled ~ poly(meanmaxtempI_scaled,2,raw = TRUE) + habitat + meanmintempI_scaled + juliandate_scaled + year + site + (1|attempt_id),data = filter(g,Species == "TRES",
                                                                                                                                                                                            !is.na(meanmaxtempI)) %>%
                               mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                             ~ scale(.x)[,1],
                                             .names = "{.col}_scaled"),
                                      meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))





meanI_WEBL <- lmerTest::lmer(gweight_scaled ~ meanmaxtempI_scaled * habitat + juliandate_scaled + year_factor + site + (1|attempt_id),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_scaled),gweight_scaled > -1, gweight_scaled < 3))

hist(g$gweight_scaled)

filter(g,Species == "WEBL",!is.nan(meanmaxtempI_scaled))$gweight %>% hist()

filter(g,Species == "WEBL",!is.nan(meanmaxtempI_scaled),gweight < 0)

meanI_WEBL <- lmerTest::lmer(pweight_scaled~meanmaxtempI_scaled * habitat + meanmaxtempI_scaled_sq + juliandate_scaled + site + (1|Nestbox),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_scaled)))


meanI_WEBL_nosq <- lmerTest::lmer(gweight_1~meanmaxtempI_1 * habitat + juliandate_1 + (1|Nestbox),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_1)))
anova(meanI_WEBL,meanI_WEBL_nosq)
meanI_WEBL_nosqint <- lmerTest::lmer(gweight_1~meanmaxtempI_1 * habitat + poly(meanmaxtempI_1,2) + juliandate_1 + (1|Nestbox),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_1)))
anova(meanI_WEBL,meanI_WEBL_nosqint)
meanI_WEBL_noint <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) + habitat + juliandate_1 + (1|Nestbox),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_1)))
anova(meanI_WEBL,meanI_WEBL_noint)



(p <- ggpredict(meanI_WEBL,terms = c("meanmaxtempI_1 [all]","habitat")) %>% plot(line.size = 1.5,
                                                                    alpha = .2) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_manual(values = c("transparent","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    scale_color_manual(values = c("transparent","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    theme(text = element_text(size = 13),legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent"))
)
ggsave("figures/growthbymeanI_WEBL_rowcrop.png",p,bg = "transparent",width = 6,height = 4)
(p <- ggpredict(meanI_WEBL,terms = c("meanmaxtempI_1 [all]","habitat")) %>% plot(line.size = 1.5,
                                                                    alpha = .2) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled)") +
    ylab("Nestling growth (z-scaled)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_manual(values = c("#440154FF","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    scale_color_manual(values = c("#440154FF","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    theme(text = element_text(size = 16),legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent"))
)
ggsave("figures/growthbymeanI_WEBL_rowcropforest.png",p,bg = "transparent",width = 6,height = 4)
(p <- ggpredict(meanI_WEBL,terms = c("meanmaxtempI_1 [all]","habitat")) %>% plot(line.size = 1.5,
                                                                    alpha = .2) +
    theme_classic() +
    xlab("Weekly mean maximum temperature") +
    ylab("Weekly nestling growth") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_manual(values = c("#440154FF","transparent","#35B779FF","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    scale_color_manual(values = c("#440154FF","transparent","#35B779FF","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    theme(text = element_text(size = 16),legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent"))
)
ggsave("figures/growthbymeanI_WEBL_rowcropforestgrassland.png",p,bg = "transparent",width = 6,height = 4)

# canopy?

meanI_WEBL_canopy <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) * canopy_cover_1 + juliandate_1 + (1|Nestbox),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_1)))
meanI_WEBL_canopy_noint <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) + canopy_cover_1 + juliandate_1 + (1|Nestbox),data = filter(g,Species == "WEBL",!is.nan(meanmaxtempI_1)))
anova(meanI_WEBL_canopy,meanI_WEBL_canopy_noint)
(p <- ggpredict(meanI_WEBL_canopy_noint,terms = c("meanmaxtempI_1","canopy_cover_1")) %>%
    filter(group != "0.34") %>%
    plot(line.size = 1.5,alpha = .2) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "canopy_cover_1") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13)) +
    labs(title = element_blank(),color = "Canopy cover") +
    scale_fill_manual(values = c("#E14209FF","#30123BFF"),labels = c("Low","High")) +
    scale_color_manual(values = c("#E14209FF","#30123BFF"),labels = c("Low","High"))
)

summary(meanI_WEBL_canopy_noint)

ggsave("figures/growthbymeanI_WEBL_canopy.png",p,width = 6,height = 4)

# back to original q with tres

meanI_TRES <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) * habitat + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1)))
meanI_TRES_nosqint <- lmerTest::lmer(gweight_1~meanmaxtempI_1 * habitat + poly(meanmaxtempI_1,2) + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1)))
meanI_TRES_nosq <- lmerTest::lmer(gweight_1~meanmaxtempI_1 * habitat + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1)))
anova(meanI_TRES,meanI_TRES_nosq)
meanI_TRES_noint <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) + habitat + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1)))
anova(meanI_TRES,meanI_TRES_noint)
(p <- ggpredict(meanI_TRES_nosq,terms = c("meanmaxtempI_1","habitat")) %>% plot(line.size = 1.5,
                                                                         alpha = .2) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled)") +
    ylab("Nestling growth (z-scaled)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 16))
)
ggsave("figures/growthbymeanI_TRES.png",p,width = 6,height = 4)


(p <- ggpredict(meanI_TRES_nosq,terms = c("meanmaxtempI_1","habitat")) %>% plot(line.size = 1.5,
                                                                         alpha = .2) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled)") +
    ylab("Nestling growth (z-scaled)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_manual(values = c("transparent","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    scale_color_manual(values = c("transparent","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    theme(text = element_text(size = 16),legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent"))
)
ggsave("figures/growthbymeanI_TRES_rowcrop.png",p,bg = "transparent",width = 6,height = 4)
(p <- ggpredict(meanI_TRES_nosq,terms = c("meanmaxtempI_1","habitat")) %>% plot(line.size = 1.5,
                                                                         alpha = .2) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled)") +
    ylab("Nestling growth (z-scaled)") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_manual(values = c("#440154FF","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    scale_color_manual(values = c("#440154FF","transparent","transparent","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    theme(text = element_text(size = 16),legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent"))
)
ggsave("figures/growthbymeanI_TRES_rowcropforest.png",p,bg = "transparent",width = 6,height = 4)
(p <- ggpredict(meanI_TRES_nosq,terms = c("meanmaxtempI_1","habitat")) %>% plot(line.size = 1.5,
                                                                         alpha = .2) +
    theme_classic() +
    xlab("Weekly mean maximum temperature") +
    ylab("Weekly nestling growth") +
    labs(title = element_blank(),color = "Cover type") +
    scale_fill_manual(values = c("#440154FF","transparent","#35B779FF","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    scale_color_manual(values = c("#440154FF","transparent","#35B779FF","#FDE725FF"),guide = guide_legend(override.aes = list(alpha = 0) )) +
    theme(text = element_text(size = 16),legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),axis.text = element_text(color = "transparent"),axis.ticks = element_line(color = "transparent"),axis.line = element_line(color = "transparent"),axis.title = element_text(color = "transparent"),panel.background = element_rect(fill = "transparent"),plot.background = element_rect(fill = "transparent"),panel.border = element_rect(color = "transparent",fill = "transparent"),legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(color = "transparent",fill = "transparent"))
)
ggsave("figures/growthbymeanI_TRES_rowcropforestgrassland.png",p,bg = "transparent",width = 6,height = 4)

# canopy?
meanI_TRES_canopy <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) * canopy_cover_1 + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1),gweight_1 < 2))
meanI_TRES_canopy_nosq <- lmerTest::lmer(gweight_1~meanmaxtempI_1 * canopy_cover_1 + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1),gweight_1 < 2))
anova(meanI_TRES_canopy,meanI_TRES_canopy_nosq)
meanI_TRES_canopy_nosqint <- lmerTest::lmer(gweight_1~meanmaxtempI_1 * canopy_cover_1 + poly(meanmaxtempI_1,2) + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1),gweight_1 < 2))
anova(meanI_TRES_canopy,meanI_TRES_canopy_nosqint)
meanI_TRES_canopy_noint <- lmerTest::lmer(gweight_1~poly(meanmaxtempI_1,2) + canopy_cover_1 + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1),gweight_1 < 2))
anova(meanI_TRES_canopy,meanI_TRES_canopy_noint)
meanI_TRES_canopy_noint_nosq <- lmerTest::lmer(gweight_1~meanmaxtempI_1 + canopy_cover_1 + juliandate_1 + (1|Nestbox),data = filter(g,Species == "TRES",!is.nan(meanmaxtempI_1),gweight_1 < 2))
anova(meanI_TRES_canopy,meanI_TRES_canopy_noint_nosq)
(p <- ggpredict(meanI_TRES_canopy,terms = c("meanmaxtempI_1","canopy_cover_1")) %>%
    plot(line.size = 1.5,alpha = .2,rawdata = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled)") +
    ylab("Nestling growth (z-scaled)") +
    labs(title = element_blank(),color = "Canopy cover") +
    scale_fill_viridis() +
    scale_color_viridis() +
    theme(text = element_text(size = 16))
)
ggsave("figures/growthbymeanI_TRES_canopy_sq.png",p,width = 6,height = 4)


(p <- ggpredict(meanI_canopy,terms = c("meanmaxtempI","canopy_cover")) %>%
    filter(group != "0.03") %>%
    plot(line.size = 1.5,alpha = .2) +
    theme_classic() +
    xlab("Weekly mean maximum temperature") +
    ylab("Weekly nestling growth") +
    labs(title = element_blank(),color = "Canopy cover") +
    scale_fill_manual(values = c("#E14209FF","#30123BFF"),labels = c("Low","High")) +
    scale_color_manual(values = c("#E14209FF","#30123BFF"),labels = c("Low","High")))
ggsave("../figures/growthbymeanI_canopy.png",p,width = 6,height = 4)
summary(meanI_canopy)
