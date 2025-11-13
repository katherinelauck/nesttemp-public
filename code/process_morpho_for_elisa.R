##### Process morpho for blood number accounting
#### Author: Katherine Lauck
#### Last updated: 23 August 2023

require(tidyverse)
require(lubridate)
require(ggeffects)
require(lme4)
library(hemispheR)
library(viridis)
library(ggpubr)
library(patchwork)
source("helper_functions.R")

g <- nestbox_drive_get("banding-and-morphometrics") %>%
  mutate(`Blood number` = ifelse(str_detect(`Blood number`,"(A|B)([:digit:]){3}"),`Blood number`,str_extract(`Comments`,"(A|B)([:digit:]){3}"))) %>%
  arrange(`Blood number`) %>%
  select(c(`Blood number`,Species)) %>%
  filter(!is.na(`Blood number`)) %>%
  mutate(plate = "",
         wells1.1 = "",
         wells1.2 = "",
         wells2.1 = "",
         wells2.2 = "")

write_csv(g,"data/blood_num_elisa.csv")
