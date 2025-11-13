require(tidyverse)
require(lubridate)
source("code/helper_functions.R")

t <- get_temp_data("data/temp-loggers") %>% list_rbind() %>%
  mutate(site = factor(str_extract(box,"([:alpha:]|[:punct:])+")),
         habitat = fct_collapse(site,
                                Orchard = c("BRO","MCE","PICO","MCO","FBFO","PG-O-"),
                                Grassland = c("MBNG","RRRG","PICG","FBFG"),
                                `Row crop` = c("MBNR","PICR","RRRR","FBFR","PG-C-","RU-C-"),
                                Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","PCE")),
         habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")))

write_rds(t,"data/temp.rds")

h <- get_humidity_data("data/humidity sensors") %>% list_rbind() %>%
  mutate(site = factor(str_extract(box,"([:alpha:]|[:punct:])+")),
         habitat = fct_collapse(site,
                                Orchard = c("BRO","MCE","PICO","MCO","FBFO","PG-O-"),
                                Grassland = c("MBNG","RRRG","PICG","FBFG"),
                                `Row crop` = c("MBNR","PICR","RRRR","FBFR","PG-C-","RU-C-"),
                                Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","PCE")),
         habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")))

write_rds(h,"data/humidity.rds")


wmean <- t %>%
  group_by(date,habitat,logger_position) %>%
  summarize(meanmax = mean(maxt),n = n(),.groups = "keep") %>%
  ungroup(habitat) %>%
  summarize(wmean = weighted.mean(meanmax,n)) %>%
  right_join(t,by = c("date","logger_position")) %>%
  mutate(resid = maxt - wmean)

write_rds(wmean,"data/wmean.rds")
