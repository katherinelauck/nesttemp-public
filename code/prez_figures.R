### make prez figures
library(tidyverse)
library(ggplot2)
library(lubridate)
library(scales)
library(viridis)
wmean <- read_rds("data/wmean.rds")

g <- read_rds("data/growth.rds") %>%
  filter(!is.na(gweight),!is.infinite(gweight),!is.infinite(meanmaxtempI),
         !is.infinite(meanmaxtempO))

canopy_trans <- trans_new("canopy_trans",
                          transform = function(x){(x * sd(g$canopy_cover,na.rm = TRUE)) +
                              mean(g$canopy_cover,na.rm = TRUE)},
                          inverse = function(x){x})
temp_trans <- trans_new("temp_trans",
                          transform = function(x){(x * sd(g$meanmaxtempI,na.rm = TRUE)) +
                              mean(g$meanmaxtempI,na.rm = TRUE)},
                          inverse = function(x){x})

(out <- filter(wmean,logger_position == "O") %>%
    ggplot(mapping = aes(x = habitat, y = resid)) +
    geom_boxplot(aes(fill = habitat)) +
    xlab("Cover type") +
    ylab("Daily max temp above daily mean max (\u00B0C)") +
    #stat_summary(fun.y = median, fun.ymax = length,
    #geom = "text", aes(label = ..ymax..), vjust = -1) +
    labs(title = element_blank(),fill = "Cover type") +
    annotate("text",x = c(1:4),y = 16,label = c("a","b","c","d")) +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic() +
    theme(text = element_text(size = 13)))
ggsave("figures/max-weightedmean_outside.png",out,width = 6,height = 4)
TukeyHSD(aov(resid~habitat, filter(wmean,logger_position == "O")),conf.level = .95)


TukeyHSD(aov(canopy_cover~habitat,g,conf.level = .95))

(box <- g %>% ggplot(mapping = aes(x = habitat,y = canopy_cover)) +
    geom_boxplot(aes(fill = habitat)) +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic() +
    annotate("label",x = c(1:4),y = 1.4,label = c("a","b","c","c"),fill = "white",label.size = NA) +
    labs(x = "Cover type",
         y = "Canopy cover",
         title = element_blank(),
         fill = "Cover type") +
    theme(text = element_text(size = 16)) +
    scale_y_continuous(trans = canopy_trans,
                       #breaks = c((0-mean(g$canopy_cover,na.rm = TRUE))/sd(g$canopy_cover,na.rm = TRUE),
                        #          (50-mean(g$canopy_cover,na.rm = TRUE))/sd(g$canopy_cover,na.rm = TRUE),
                         #         (100-mean(g$canopy_cover,na.rm = TRUE))/sd(g$canopy_cover,na.rm = TRUE)),
                       breaks = c(0,50,100),
                       labels = c("0%","50%","100%")))


ggsave("figures/canopybyhabitat.png",box,width = 6,height = 4)

(box <- g %>% ggplot(aes(x = canopy_cover,y = meanmaxtempO)) +
    geom_smooth(method = "lm", se = TRUE,color = "black") +
    geom_point() +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic()
  +
    labs(x = "Canopy cover",
         y = "Max temperature (\u00B0C)",
         title = element_blank()) +
    theme(text = element_text(size = 13)) +
    scale_x_continuous(trans = canopy_trans,
                       breaks = c(0,50,100),
                       labels = c("0%","50%","100%")) +
    scale_y_continuous(trans = temp_trans,
                       # breaks = c((20-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (30-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (40-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (50-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE)),
                       breaks = c(20,30,40,50),
                       labels = c("20","30","40","50"),
                       # limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
                       #            (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE))
                       limits = c(18,52)
                       )
)
lm(meanmaxtempO ~ canopy_cover,data = g) %>% summary()

ggsave("figures/tempbycanopy.png",box,width = 6,height = 4)

(box <- g %>% filter(Species == "WEBL",gweight>-5,gweight<5,meanmaxtempI<50  ) %>% ggplot(aes(x = meanmaxtempI, y = gweight)) +
    geom_smooth(method = "lm",formula = y ~ x + I(x^2), se = TRUE,color = "black") +
    geom_point() +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic()
  +
    labs(x = "Maximum temperature (C)",
         y = "Nestling growth (g)",
         title = "WEBL") +
    theme(text = element_text(size = 13))
  # scale_y_continuous(trans = temp_trans,
  #                    breaks = c((20-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (30-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (40-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (50-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE)),
  #                    labels = c("20","30","40","50"),
  #                    limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE)))
)
lm(gweight ~ meanmaxtempI + I(meanmaxtempI^2),data = g%>% filter(Species=="WEBL")) %>% summary()

ggsave("figures/weblgrowthbytemp.png",box,width = 6,height = 4)

(box <- g %>% filter(Species == "TRES",gweight<3.7,meanmaxtempI<50) %>% ggplot(aes(x = meanmaxtempI, y = gweight)) +
    geom_smooth(method = "lm",formula = y ~ x + I(x^2), se = TRUE,color = "black") +
    geom_point() +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic()
  +
    labs(x = "Maximum temperature (C)",
         y = "Nestling growth (g)",
         title = "TRES") +
    theme(text = element_text(size = 13))
  # scale_y_continuous(trans = temp_trans,
  #                    breaks = c((20-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (30-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (40-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (50-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE)),
  #                    labels = c("20","30","40","50"),
  #                    limits = c((15-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE),
  #                               (55-mean(g$meanmaxtempI,na.rm = TRUE))/sd(g$meanmaxtempI,na.rm = TRUE)))
)
lm(gweight ~ meanmaxtempI + I( meanmaxtempI^2),data = g%>% filter(Species=="TRES")) %>% summary()

ggsave("figures/tresgrowthbytemp.png",box,width = 6,height = 4)
