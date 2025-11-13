#### Exploring 2021 nest checks
#### Author: Katherine Lauck
#### Last updated: 4/18/2022

nest <- get_nest_checks()

nest %>% filter(`First Egg Date` < "2023-01-01", Species %in% c("WEBL","TRES")) %>%
  count(Species)

t <- get_temp_data()

# Ideas for temp visualization
# 1) calculate time over 40 C - implemented in get_temp_data
# 2) standardize each logger's max temp over the average max across all loggers per day (weighted to account for sample size differences)
# 3) regress max temp on day of year and plot residuals
# 4) regress time over 40 C on day of year and plot residuals
#
# 2)

wmean <- t %>%
  group_by(date,habitat,logger_position) %>%
  summarize(meanmax = mean(max),n = n(),.groups = "keep") %>%
  ungroup(habitat) %>%
  summarize(wmean = weighted.mean(meanmax,n)) %>%
  right_join(t,by = c("date","logger_position")) %>%
  mutate(resid = max - wmean)

filter(wmean,logger_position == "O") %>% ggplot(mapping = aes(x = habitat, y = resid)) +
  geom_boxplot() +
  xlab("Cover type") +
  ylab("Max temp standardized by daily max") +
  labs(title = "Outside loggers") +
  theme_classic()

ggsave("figures/max-weightedmean_outside.png",width = 6,height = 4)

TukeyHSD(aov(resid~habitat, filter(wmean,logger_position == "O")),conf.level = .95)

filter(wmean,logger_position == "I") %>% ggplot(mapping = aes(x = habitat, y = resid)) +
  geom_boxplot() +
  xlab("Cover type") +
  ylab("Max temp standardized by daily max") +
  labs(title = "Inside loggers") +
  theme_classic()

ggsave("figures/max-weightedmean_inside.png",width = 6,height = 4)

TukeyHSD(aov(resid~habitat, filter(wmean,logger_position == "I")),conf.level = .95)

# 3)

lmresid <- tibble(lmresid = resid(aov(max~yday(date),wmean)))
lmresid <- wmean %>% ungroup() %>% mutate(lmresid = pull(lmresid,lmresid))

filter(lmresid,logger_position == "O") %>% ggplot(mapping = aes(x = habitat, y = lmresid)) +
  geom_boxplot() +
  xlab("Cover type") +
  ylab("Residual of max temp ~ day of year") +
  labs(title = "Outside loggers") +
  theme_classic()

ggsave("figures/resid_outside.png",width = 6,height = 4)

TukeyHSD(aov(lmresid~habitat, filter(lmresid,logger_position == "O")),conf.level = .95)

filter(lmresid,logger_position == "I") %>% ggplot(mapping = aes(x = habitat, y = lmresid)) +
  geom_boxplot() +
  xlab("Cover type") +
  ylab("Residual of max temp ~ day of year") +
  labs(title = "Inside loggers") +
  theme_classic()

ggsave("figures/resid_inside.png",width = 6,height = 4)

TukeyHSD(aov(lmresid~habitat, filter(lmresid,logger_position == "I")),conf.level = .95)

# 4)

resid40 <- tibble(resid40 = resid(aov(time_above_40~yday(date),lmresid)))
lmresid <- lmresid %>% mutate(resid40 = pull(resid40,resid40))

filter(lmresid,logger_position == "O") %>% ggplot(mapping = aes(x = habitat, y = time_above_40)) +
  geom_boxplot() +
  xlab("Cover type") +
  ylab("Residual of time above 40C ~ day of year") +
  labs(title = "Outside loggers") +
  theme_classic()

ggsave("figures/resid_timeabove40_outside.png",width = 6,height = 4)

TukeyHSD(aov(resid40~habitat, filter(lmresid,logger_position == "O")),conf.level = .95)

filter(lmresid,logger_position == "I") %>% ggplot(mapping = aes(x = habitat, y = time_above_40)) +
  geom_boxplot() +
  xlab("Cover type") +
  ylab("Residual of time above 40C ~ day of year") +
  labs(title = "Inside loggers") +
  theme_classic()

ggsave("figures/resid_timeabove40_inside.png",width = 6,height = 4)

TukeyHSD(aov(resid40~habitat, filter(lmresid,logger_position == "I")),conf.level = .95)
