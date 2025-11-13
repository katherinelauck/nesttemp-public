#### Preliminary cort modeling

library(tidyverse)
require(ggeffects)
require(lme4)
require(viridis)
require(insight)
source("code/helper_functions.R")

c <- read_rds("data/cort_with_id.rds")
cort <- read_rds("data/growth_cort.rds")

only_adultf <- c %>% filter(Sex == "F" & Age %in% c("AHY","ASY","SY")) %>%
  mutate()

habitat_s1 <- lm(s1_scaled ~ habitat + brood_size_scaled,data = filter(only_adultf,Species == "TRES")) #swap julian date for temp
summary(habitat_s1)
habitat_s2 <- lm(s2_scaled ~ habitat + juliandate_scaled + brood_size_scaled,data = filter(only_adultf,Species == "WEBL"))
summary(habitat_s2)

habitat_s2_s1 <- lm(s2_scaled ~ s1_scaled + habitat + juliandate_scaled + brood_size_scaled, data = only_adultf)
summary(habitat_s2_s1)

nestling <- c %>% filter(Age == "L")

ms1_s1 <- lmerTest::lmer(s1_scaled ~ s1_mother_scaled + juliandate_scaled + brood_size_scaled + (1|Nestbox/color), data = nestling)
summary(ms1_s1)

ms2_s2 <- lmerTest::lmer(s2_scaled ~ s2_mother_scaled + juliandate_scaled + brood_size_scaled + (1|Nestbox/color), data = nestling)
summary(ms2_s2)

ms2_s1 <- lmerTest::lmer(s2_scaled ~ s1_mother_scaled + juliandate_scaled + brood_size_scaled + (1|Nestbox/color), data = nestling)
summary(ms2_s1)

get_variance_intercept(ms2_s1)

ms1 <- lm(s1_scaled ~ habitat * juliandate_scaled * condition_scaled, data = filter(n, Species == "WEBL"))
summary(ms1) # replace julian date with temp here too

ms2 <- lm(s2_scaled ~ habitat * juliandate_scaled * condition_scaled, data = filter(n, Species == "WEBL"))
summary(ms2) # replace julian date with temp here too

ms1 <- lm(s1_scaled ~ habitat * juliandate_scaled * condition_scaled, data = filter(n, Species == "TRES"))
summary(ms1) # replace julian date with temp here too

ms2 <- lm(s2_scaled ~ habitat * juliandate_scaled * condition_scaled, data = filter(n, Species == "TRES"))
summary(ms2) # replace julian date with temp here too

n <- mutate(nestling,
            condition = weight/tarsus,
            condition_scaled = scale(condition))
n <- n %>% group_by(year, Nestbox, attempt) %>%
  summarize(mean_s1 = mean(s1), mean_s2 = mean(s2))
