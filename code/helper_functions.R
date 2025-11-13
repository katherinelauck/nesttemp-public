#### Functions to load data from nestbox-temp-monitoring shared drive
#### Author: Katherine Lauck
#### Last updated: 4/18/2022


nestbox_drive_get <- function(sheet_name) {
  ### Return tibble of google sheet from shared drive nestbox-temp-monitoring given the quoted name of the google sheet without the file extension.
  ### Also checks for the correct package googledrive version.
  require(tidyverse)
  require(googledrive)
  require(googlesheets4)
  require(rlang)

  if(compareVersion(as.character(packageVersion("googledrive")),"1.9.0.9000") == -1){
    abort("Oopsies, your version of package googledrive is out of date! Update with devtools::install_github(\"tidyverse/googledrive\")")
  }
  drive_get(sheet_name, shared_drive = "nestbox-temp-monitoring") %>%
    read_sheet(col_types = "????????cccc????")
}


gear_in_field <- function() {
  ### Return tibble of gear remaining in field.

  gear <- nestbox_drive_get("equipment_tracking")
  gear_out <- filter(gear, action == "deploy" | action == "collect") %>%
    count(equipment_name) %>%
    filter(n < 2)

  return(filter(gear, equipment_name %in% gear_out$equipment_name))

}

get_temp_data <- function(dir){

  require(lubridate)
  require(tidyverse)

  ### Return tibble of daily high, daily low, and daily range temperature from HOBO loggers. Other columns included are date, habitat, canopy cover, box, site, and whether the logger is on the inside or outside of the box.

  f <- list.files(dir,pattern = "csv",full.names = TRUE)
  collapse <- function(file_name) {
    # Given relative file name of temperature log from one logger, collapse into tibble with one row for each complete day of temperature data.
    if(read_csv(file_name) %>% problems() %>% nrow() > 100){
      t <- (readLines(file_name) %>% tibble() %>% slice(2:n()) %>% pull() %>%
              I() %>% read_csv(col_names = c("x1","x2","x3")))[,2:3]
    } else {
      t <- read_csv(file_name)[,2:3]
    }
    logger_name <- str_extract(file_name,paste0("(?<=^",dir,"/)[:graph:]+(?=\\s)"))
    names(t) <- c("datetime","temp")
    t$datetime <- mdy_hms(t$datetime, tz = "America/Los_Angeles")
    t$date <- t$datetime %>% date()
    above_40 <- function(d,...){
      ifelse(nrow(filter(d,temp > 40,.preserve = TRUE)) == 0,
             0,
             interval(filter(d,temp > 40,.preserve = TRUE) %>% pull(datetime) %>% min(na.rm = TRUE),
                      filter(d,temp > 40,.preserve = TRUE) %>% pull(datetime) %>% max(na.rm = TRUE)) / dhours(1))
    }
    time_above_40 <- group_map(t %>% group_by(date),above_40)
    sum <- t %>%
      group_by(date) %>%
      summarize(meant = mean(temp,na.rm = TRUE),
                maxt = quantile(temp,probs = 0.95,na.rm = TRUE),
                mint = quantile(temp,probs = 0.05,na.rm = TRUE),
                ranget = maxt-mint) %>%
      mutate(time_above_40 = unlist(time_above_40),
             box = str_extract(logger_name,"[:graph:]*(?=-[:alpha:])"),
             logger_position = str_extract(logger_name,"(?<=-)[:alpha:]+$"),
             identity = paste(box,logger_position,as.character(date(date[1])),sep = "_"))
  }

  return(map(f,collapse))

}


get_humidity_data <- function(dir){

  require(lubridate)
  require(tidyverse)

  ### Return tibble of daily high, daily low, and daily range temperature from HOBO loggers. Other columns included are date, habitat, canopy cover, box, site, and whether the logger is on the inside or outside of the box.

  f <- list.files(dir,pattern = "csv",full.names = TRUE)
  collapse <- function(file_name) {
    # Given relative file name of temperature log from one logger, collapse into tibble with one row for each complete day of temperature data.
    if(read_csv(file_name) %>% problems() %>% nrow() > 100){
      t <- read_csv(file_name)[,c(2,4)]
      # t <- (readLines(file_name) %>% tibble() %>% slice(2:n()) %>% pull() %>%
      #         I() %>% read_csv(col_names = c("x1","x2","x3")))[,2:3]
    } else {
      t <- read_csv(file_name)[,c(2,4)]
    }
    logger_name <- str_extract(file_name,paste0("(?<=^",dir,"/)[:graph:]+(?=\\s)"))
    names(t) <- c("datetime","humidity")
    t$datetime <- mdy_hms(t$datetime, tz = "America/Los_Angeles")
    t$date <- t$datetime %>% date()
    # above_40 <- function(d,...){
    #   ifelse(nrow(filter(d,temp > 40,.preserve = TRUE)) == 0,
    #          0,
    #          interval(filter(d,temp > 40,.preserve = TRUE) %>% pull(datetime) %>% min(na.rm = TRUE),
    #                   filter(d,temp > 40,.preserve = TRUE) %>% pull(datetime) %>% max(na.rm = TRUE)) / dhours(1))
    # }
    # time_above_40 <- group_map(t %>% group_by(date),above_40)
    sum <- t %>%
      group_by(date) %>%
      mutate(humidity = as.numeric(humidity)) %>%
      summarize(meanh = mean(humidity,na.rm = TRUE),
                maxh = quantile(humidity,probs = 0.95,na.rm = TRUE),
                minh = quantile(humidity,probs = 0.05,na.rm = TRUE),
                range = maxh-minh) %>%
      mutate(box = str_extract(logger_name,"[:graph:]*(?=-[:alpha:])"))
  }

  return(map(f,collapse))

}


get_cort <- function() {
  # Construct cort data from raw assay results. Steps:
  # 1. Load sample identity data
  # 2. Load cort data as a series of matrices, the names of which are the assay number
  # 3. Use date and assay number to match sample identity to results

  library(tidyverse)
  library(lubridate)
  id <- read_csv("data/elisa_assay_log.csv") %>%
    pivot_longer(cols = c("s1.1","s1.2","s2.1","s2.2"),names_to = "sample", values_to = "cell")

  cort_files <- list.files("data/cort-assay-results/",pattern = "^A\\d{2}_clean.csv",full.names = TRUE)

  c <- cort_files %>% map(read_csv,col_types = "cc") %>%
    set_names(nm = str_extract(cort_files,pattern = "A\\d{2}")) %>%
    list_rbind(names_to = "plate")

  cort <- left_join(id,c,by = c("plate" = "plate","cell" = "well"),na_matches = "never",keep = FALSE) %>%
    mutate(conc = as.numeric(conc)) %>%
    filter(!is.na(blood_num),!is.na(species),!is.na(conc)) %>%
    dplyr::select(!c("cell")) %>%
    mutate(conc = conc/1000) %>%
    mutate(replicate = as.factor(sample) %>% fct_collapse(r1 = c("s1.1","s2.1"),r2 = c("s1.2","s2.2")),
           sample = as.factor(sample) %>% fct_collapse(s1 = c("s1.1","s1.2"),s2 = c("s2.1","s2.2"))) %>%
    pivot_wider(names_from = replicate,values_from = conc) %>%
    rowwise() %>%
    mutate(mean = mean(c(r1,r2),na.rm = TRUE),d = abs(r1-r2)) %>%
    ungroup() %>%
    pivot_wider(names_from = sample,values_from = c(r1,r2,mean,d)) %>%
    mutate(date = mdy(date))


  return(cort)

}
