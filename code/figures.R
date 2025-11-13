### important figures collected

library(tidyverse)
require(ggeffects)
require(lme4)
require(viridis)
library(car)
library(gt)
g <- read_rds("data/growth.rds")


## plot showing interactive effect in WEBL and TRES

meanI_WEBL_nosq <- lmerTest::lmer(gweight_scaled ~ meanmaxtempI_scaled * habitat + poly(meanmintempI_scaled,2) * habitat + age + juliandate_scaled + (1|year) + (1|site/attempt_id),
                                  data = filter(g,Species == "WEBL",!is.na(meanmaxtempI)) %>%
                                    mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

model <- meanI_WEBL_nosq

p <- ggpredict(model,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                alpha = .2,show_data = TRUE) +
  theme_classic() +
  xlab("Maximum temperature (\u00B0C)") +
  ylab("Nestling growth (g/day)") +
  labs(title = element_blank(),color = "Cover type") +
  # ggtitle("Predicted growth across habitat and mean max daily temp") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(text = element_text(size = 13))

unscale_x <- function(axis_break) {
  if(is.na(axis_break)) {
    return(NA)
  } else {
    sd_x = getData(model) %>% pull(meanmaxtempI) %>% sd(na.rm = TRUE)
    mean_x = getData(model) %>% pull(meanmaxtempI) %>% mean(na.rm = TRUE)
    out <- (axis_break * sd_x + mean_x) %>% round(digits = 0)
  }
}

unscale_y <- function(axis_break) {
  if(is.na(axis_break)) {
    return(NA)
  } else {
    sd_y = getData(model) %>% pull(gweight) %>% sd(na.rm = TRUE)
    mean_y = getData(model) %>% pull(gweight) %>% mean(na.rm = TRUE)
    out <- (axis_break * sd_y + mean_y) %>% round(digits = 1)
  }
}



(p_rescaled <- p +
    scale_x_continuous(labels = map(layer_scales(p)$x$break_positions(),unscale_x) %>% unlist()) +
    scale_y_continuous(labels = map(layer_scales(p)$y$break_positions(),unscale_y) %>% unlist()))

ggsave("figures/rescaled_meanI_WEBL_nosq.png",p_rescaled,height = 4,width = 6)

### TRES

meanI_TRES_nosq <- lmerTest::lmer(gweight_scaled ~ meanmaxtempI_scaled * habitat + poly(meanmintempI_scaled,2) * habitat + age + juliandate_scaled + (1|year) + (1|site/attempt_id),data = filter(g,Species == "TRES",
                                                                                                                                                                                                  !is.na(meanmaxtempI)) %>%
                                    filter(meanmaxtempI < 45) %>%
                                    mutate(across(c(gweight,meanmaxtempI,meanmintempI,juliandate,brood_size),
                                                  ~ scale(.x)[,1],
                                                  .names = "{.col}_scaled"),
                                           meanmaxtempI_scaled_sq = meanmaxtempI_scaled * meanmaxtempI_scaled))

(p <- ggpredict(meanI_TRES_nosq,terms = c("meanmaxtempI_scaled [all]","habitat")) %>% plot(line.size = 1.5,
                                                                                           alpha = .2,show_data = TRUE) +
    theme_classic() +
    xlab("Maximum temperature (z-scaled; ~20-55C)") +
    ylab("Daily nestling growth (z-scaled; ~-2-8g)") +
    labs(title = element_blank(),color = "Cover type") +
    ggtitle("Predicted growth across habitat and mean max daily temp") +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    theme(text = element_text(size = 13))
)
