##### analyze canopy cover files

require(tidyverse)
library(hemispheR)

args <- commandArgs(trailingOnly = TRUE)

pix <- paste0("data/canopy-cover-pictures/",args[1]) %>%
  import_fisheye() %>%
  binarize_fisheye() %>%
  gapfrac_fisheye() %>%
  canopy_fisheye() %>%
  select(id,DIFN) %>%
  write_rds(file = paste0("data/canopy-cover-pictures/",args[1],".rds"))
