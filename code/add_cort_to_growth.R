### cort modeling

library(tidyverse)
require(ggeffects)
require(lme4)
require(viridis)
source("code/helper_functions.R")

g <- read_rds("data/growth.rds")
c <- get_cort()

c <- c %>% pivot_wider(names_from = sample, values_from = mean) %>% arrange(date)

gc <- left_join(g,
                select(c,c(blood_num,s1,s2)),
                by = c("blood_num"),
                multiple = "first")
write_rds(gc,"data/growth_cort.rds")


