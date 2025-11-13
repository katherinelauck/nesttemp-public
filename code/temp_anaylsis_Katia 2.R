##### Temp differences among habitats, across canopy cover, inside vs outside
##### Author: Katherine Lauck
##### Last updated: 10 January 2022

source("code/helper_functions.R")

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

mod.dif<-lmer(formula = diff ~ habitat + (1 | site), data = diffs)
summary(mod.dif)

ggplot(diffs,aes(x = habitat,y = diff)) +
  geom_boxplot()

mod.max.I <- lmer(formula = max ~ habitat + as.numeric(date) + (1|site) , data = filter(t_grouped,logger_position == "I"))
mod.max.O <- lmer(formula = max ~ habitat + as.numeric(date) + (1|site) , data = filter(t_grouped,logger_position == "O"))
mod.max <- lmer(formula = max ~ habitat + logger_position + as.numeric(date) + (1|site) , data = t_grouped)

ggplot(filter(t_grouped,logger_position == "I"),aes(x = date,y = max,color = habitat)) +
  geom_point()
ggplot(filter(t_grouped,logger_position == "O"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(diffs,aes(x = habitat,y = diff)) data == filter(diffs,habitat == "forest") +
  geom_boxplot()


data == filter(diffs,habitat == "forest")

ggplot(filter(t_grouped,logger_position == "O"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(diffs,habitat == "orchard"),aes(x = date,y = diff,color = habitat)) +
  geom_point()

ggplot(filter(diffs,habitat == "forest"),aes(x = date,y = diff,color = habitat)) +
  geom_point()

ggplot(filter(diffs,habitat == "row crop"),aes(x = date,y = diff,color = habitat)) +
  geom_point()

ggplot(filter(diffs,habitat == "grassland"),aes(x = date,y = diff,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "PICG"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "PICR"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "PCC"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "PICO"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "RRRR"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "RRRG"),aes(x = date,y = max,color = habitat)) +
  geom_point()

ggplot(filter(t_grouped,site == "RRC"),aes(x = date,y = max,color = habitat)) +
  geom_point()
