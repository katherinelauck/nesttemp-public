require(tidyverse)
require(lubridate)
library(hemispheR)
source("helper_functions.R")

list <- list.files("data/canopy-cover-pictures",pattern = "JPG",full.names = TRUE)

write_rds(list,"data/cc_list.csv")
