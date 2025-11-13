##### Temp differences among habitats, across canopy cover, inside vs outside
##### Author: Katherine Lauck
##### Last updated: 10 January 2022

source("code/helper_functions.R")
install.packages("lme4")
install.packages("tidyverse")
install.packages("lmerTest")
install.packages("lubridate")
install.packages("map_dfr")
library(lme4)
library(tidyverse)
library(lmerTest)

t <- get_temp_data()
minus <- function(x) sum(x[1],na.rm=T) - sum(x[2],na.rm=T)
hab_key <- tibble(site = unique(t$site),habitat = c("orchard","forest","grassland","row crop","orchard","orchard","forest","forest","forest","row crop","orchard","grassland","orchard","row crop","forest","forest","grassland","row crop","row crop"))
t_grouped <- t  %>% left_join(hab_key,by = "site") %>% group_by(date,box,site,habitat)

diffs <- t_grouped %>%
  arrange(logger_position,.by_group = TRUE) %>%
  summarize(diff = -minus(max))


  ungroup(diffs) %>%
  summarize(mean(diff))

  ungroup(diffs)
summarize(mean(diff(row crop))

mod.dif<-lmer(formula = diff ~ habitat + (1 | site), data = diffs)
summary(mod.dif)

ggplot(diffs,aes(x = habitat,y = Temperature)) +
  geom_boxplot()

mod.max.I <- lmer(formula = max ~ habitat + as.numeric(date) + (1|site) , data = filter(t_grouped,logger_position == "I"))
mod.max.O <- lmer(formula = max ~ habitat + as.numeric(date) + (1|site) , data = filter(t_grouped,logger_position == "O"))
mod.max <- lmer(formula = max ~ habitat + logger_position + as.numeric(date) + (1|site) , data = t_grouped)

ggplot(filter(t_grouped,logger_position == "I"),aes(x = date,y = max,color = habitat)) +
  geom_point()
ggplot(filter(t_grouped,logger_position == "O"),aes(x = date,y = max,color = habitat)) +
  geom_point()

