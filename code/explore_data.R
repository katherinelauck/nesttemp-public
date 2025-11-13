# data cleaning and exploration
require(tidyverse)


# source("code/cat_clean_data.R")
g <- read_rds("data/growth.rds")
c <- read_rds("data/cort_with_id.rds")

pull(dg,"Nestbox") %>% unique()
pull(dg,"color") %>% unique()
hist(dg$juliandate)
pull(dg,"Species") %>% unique()
pull(dg,"Age") %>% unique()
pull(dg,"Sex") %>% unique()
pull(dg,"blood_num") %>% unique()
pull(dg,"life_stage") %>% unique()
pull(dg,"attempt_id") %>% unique()
pull(dg,"ind_id") %>% unique()
pull(dg,"blood_num_mother") %>% unique()
pull(dg,"site") %>% unique()
pull(dg,"habitat") %>% unique()

hist(dg$weight)
hist(dg$gweight)
hist(dg$pweight)
plot(dg$weight ~ dg$age)
plot(dg$gweight ~ dg$age)
plot(dg$gweight ~ dg$occasion)
hist(dg$wing)
hist(dg$gwing)
hist(dg$gtail)
hist(dg$gskull)
hist(dg$brood_size)
hist(dg$age)
hist(dg$canopy_cover)
hist(dg$meanmaxtempI)
hist(dg$meanh)
plot(dg$meanh ~ dg$habitat)
hist(filter(dg,cort_s1_mother < 100)$cort_s2_mother)
hist(filter(dg,cort_s2_mother < 100)$cort_s2_mother)
hist(filter(dg,cort_s1 < 100)$cort_s1)
hist(filter(dg,cort_s2 < 100)$cort_s2)

plot(gweight ~ brood_size, data = dg)
plot(meanmaxtempI ~ meanmintempI, data = dg)
plot(gweight ~ juliandate, data = dg)

  filter(dg,ind_id %in% (filter(dg,gweight < 0) %>% pull("ind_id"))) %>% write_csv("data/belowzerogweights.csv")
filter(dg,ind_id %in% (filter(dg,gwing < 0) %>% pull("ind_id"))) %>% write_csv("data/belowzerogwing.csv")
filter(dg,ind_id %in% (filter(dg,gskull < 0) %>% pull("ind_id"))) %>% write_csv("data/belowzerogskull.csv")

filter(dg,ind_id %in% (filter(dg,gweight < 0) %>% pull("ind_id"))) %>% write_csv("data/belowzerogweights.csv")

hist(filter(c,s1 < 15,Species == "WEBL")$s1)
hist(filter(c,s1 < 15,Species == "TRES")$s1)
filter(c,s1 > 15) %>% write_csv("data/extremes1.csv")

hist(filter(c,s2 < 70,Species == "WEBL")$s2)
filter(c,s2 > 70,Species == "WEBL") %>% write_csv("data/extremeWEBLs2.csv")

hist(filter(c,s2 < 90,Species == "TRES")$s2)
filter(c,s2 > 90,Species == "TRES") %>% write_csv("data/extremeTRESs2.csv")

hist(c$mean_s2)

filter(c,mean_s1 > 20)


c <- c %>% mutate(diff_s2_s1 = mean_s2-mean_s1)

hist(c$diff_s2_s1)

filter(c,diff_s2_s1 < 0) %>% write_csv("data/mislabelled_cort.csv")
