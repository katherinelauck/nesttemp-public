# concatenate canopy cover measurements

require(tidyverse)
require(lubridate)

pix <- list.files("data/canopy-cover-pictures/",pattern = ".*rds$",full.names = TRUE) %>%
  map(read_rds) %>%
  list_rbind()

cc <- read_csv("data/canopy_cover_picture_log.csv") %>%
  select(date,nestbox,picture_file_name,canopy_cover) %>%
  mutate(date = mdy(date,tz = "America/Los_Angeles"),
         picture_file_name = ifelse(str_detect(picture_file_name,"100"),str_c(picture_file_name,".JPG"),NA),
         picture_file_name = str_replace(picture_file_name,"100-","IMG_"),
         canopy_cover = as.numeric(canopy_cover)) %>%
  left_join(pix,by = join_by(picture_file_name == id)) %>%
  mutate(canopy_cover = ifelse(is.na(canopy_cover),100-DIFN,canopy_cover),
         year = year(date)) %>%
  select(year,nestbox,canopy_cover)

write_rds(cc,"data/cc.rds")

