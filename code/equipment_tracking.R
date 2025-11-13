# How many cams were put on nests?
#
#
library(tidyverse)

g <- read_csv("data/equipment_tracking.csv") %>%
  filter(equipment_type == "camera", action == "deploy") %>%
  dplyr::select(date,box) %>%
  distinct() %>%
  mutate(date = mdy(date),year = year(date)) %>%
  dplyr::select(year,box)

b <- read_csv("data/banding-and-morphometrics_proofed.csv") %>%
  dplyr::select(`Banding Date`,Species,Nestbox) %>%
  mutate(date = mdy(`Banding Date`),year = year(date)) %>%
  dplyr:: select(year, Species , Nestbox) %>%
  distinct() %>%
  rename(box = Nestbox)

g_sp <- left_join(g,b) %>%
  dplyr::select(Species) %>%
  group_by(Species) %>%
  summarize(count = n())
