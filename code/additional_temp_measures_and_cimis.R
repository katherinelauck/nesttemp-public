
library(tidyverse)
library(lubridate)
library(furrr)
library(fuzzyjoin)
library(weathermetrics)

### Describe historical temp regime in Davis with CIMIS data

#### 1. Read in CIMIS data

# dir <- "data/cimis"
# f <- list.files(dir,pattern = "csv",full.names = TRUE)
#
# hist_temp <- map(f,read_csv) %>% list_rbind() %>% dplyr::select(-c(`qc...8`,`qc...10`,`qc...12`)) %>%
#   mutate(across(c(`Rel Hum (%)`,`Wind Speed (m/s)`),~ if_else(.x == -9997,NA,.x)))
#
# names(hist_temp) <- c("station_id","station_name","cimis_region","date","hour","julian_date","air_temp","rh","wind_speed")
#
# hist_temp_date_corrected <- hist_temp %>%
#   mutate(datetime = str_c(date,hour,sep = " "),
#          datetime = mdy_hm(datetime,tz = "America/Los_Angeles"),
#          julian_date = as.numeric(julian_date),
#          date = mdy(date,tz = "America/Los_Angeles")) %>%
#   filter(month(date) > 2,month(date) < 9)
#
# write_rds(hist_temp_date_corrected,"data/cimis_1982-2012_hourly.rds")
#
# hist_temp_date_corrected_hsum <- hist_temp_date_corrected %>%
#   mutate(month = month(datetime),
#          day = mday(datetime),
#          hour = hour(datetime)) %>%
#   summarize(hmeant = mean(air_temp,na.rm = TRUE),
#             .by = c(month, day, hour))
#
# write_rds(hist_temp_date_corrected_hsum,"data/cimis_1982-2012_hourlymeant.rds")
#
# hist_temp_sum <- hist_temp_date_corrected %>%
#   summarize(maxt = max(air_temp,na.rm = TRUE),
#             .by = c(date)) %>%
#   mutate(maxt = if_else(is.infinite(maxt),NA_integer_,maxt))
#
# write_rds(hist_temp_sum,"data/cimis_1982-2012_dailymax.rds")
#
# cimis_1883 <- read_csv("data/Davis_CIMIS_1883.csv") %>%
#   mutate(date = parse_date_time(DATE, orders = c("ymd","mdy"),tz = "America/Los_Angeles"),
#          year = year(date),
#          month = month(date),
#          day = mday(date),
#          maxt = fahrenheit.to.celsius(TMAX)) %>%
#   dplyr::select(c(date,year,month,day,maxt)) %>%
#   filter(month %in% c(3:8),
#          year <= 1913)
#
# write_rds(cimis_1883,"data/cimis_1883_dailymax.rds")
#
# str(cimis_1883)

### compare results: average across a given week per year, then across years, or across years and then a given week.

# start <- mdy("07/07/2022",tz = "America/Los_Angeles")
# end <- mdy("07/14/2022",tz = "America/Los_Angeles")
#
# wk <- tibble(date = seq(start,end,by = "days")) %>%
#   mutate(m = month(date),
#          d = mday(date))
#
# hist_temp_sum %>%
#   mutate(year = year(date),
#          month = month(date),
#          day = mday(date)) %>%
#   filter(month %in% wk$m, day %in% wk$d) %>%
#   summarize(wkmeanmax = mean(maxt,na.rm = TRUE),
#             .by = year) %>%
#   summarize(mean = mean(wkmeanmax,na.rm = TRUE))
#
# hist_temp_sum %>%
#   mutate(year = year(date),
#          month = month(date),
#          day = mday(date)) %>%
#   filter(month %in% wk$m, day %in% wk$d) %>%
#   summarize(daymeanmax = mean(maxt,na.rm = TRUE),
#             .by = c(month,day)) %>%
#   summarize(mean = mean(daymeanmax,na.rm = TRUE))

# What this means is that I can create a mean max for each month/day combo across years and then average across different intervals and get the same answer as if I averaged across the interval for each year and then across years. I suspected this would be the case but just wanted to confirm it and be able to prove it to others before implementing it.

## summarize over each day


# hist_temp_sum <- hist_temp_sum %>%
#   mutate(month = month(date),
#          day = mday(date)) %>%
#   summarize(meanmax = mean(maxt,na.rm = TRUE),
#             sdmax = sd(maxt,na.rm = TRUE),
#             .by = c(month,day))
#
# write_rds(hist_temp_sum,"data/cimis_1982-2012_dailymeanmax.rds")
#
# hist_temp_sum_1883 <- cimis_1883 %>%
#   summarize(meanmax = mean(maxt,na.rm = TRUE),
#             sdmax = sd(maxt,na.rm = TRUE),
#             .by = c(month,day))
#
# write_rds(hist_temp_sum_1883,"data/cimis_1883_dailymeanmax.rds")


### Calculate instantaneous heat index by matching each temp measure with a relevant humidity measure (datetime and site)
###


# dir <- "data/temp-loggers"
# f <- list.files(dir,pattern = "csv",full.names = TRUE)
#
# collapse <- function(file_name) {
#   # Given relative file name of temperature log from one logger, return all rows.
#   if(read_csv(file_name) %>% problems() %>% nrow() > 100){
#     t <- (readLines(file_name) %>% tibble() %>% slice(2:n()) %>% pull() %>%
#             I() %>% read_csv(col_names = c("x1","x2","x3")))[,2:3]
#   } else {
#     t <- read_csv(file_name)[,2:3]
#   }
#   logger_name <- str_extract(file_name,paste0("(?<=^",dir,"/)[:graph:]+(?=\\s)"))
#   names(t) <- c("datetime","temp")
#   t$datetime <- mdy_hms(t$datetime, tz = "America/Los_Angeles")
#   t$date <- t$datetime %>% date()
#   t <- t %>%
#     mutate(box = str_extract(logger_name,"[:graph:]*(?=-[:alpha:])"),
#            logger_position = str_extract(logger_name,"(?<=-)[:alpha:]+$"),
#            identity = paste(box,logger_position,as.character(date(date[1])),sep = "_"),
#            site = factor(str_extract(box,"([:alpha:]|[:punct:])+")), # extract site from nestbox
#            habitat = fct_collapse(site, # collapse sites into habitat classes
#                                   Orchard = c("BRO","MCE","PICO","MCO","FBFO","PG-O-"),
#                                   Grassland = c("MBNG","RRRG","PICG","FBFG"),
#                                   `Row crop` = c("MBNR","PICR","RRRR","FBFR","PG-C-","RU-C-"),
#                                   Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","RRRC")),
#            habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")),
#            site = str_replace_all(site,"MCE","MCO"), # Collapse site vectors when there are multiple site abbreviations per actual site
#            site = str_replace_all(site,"PCT","PCC"),
#            site = str_replace_all(site,"PCE","PCC"),
#            site = str_replace_all(site,"RRE","RRC"),
#            site = str_replace_all(site,"SFP","MBNC"),
#            site = str_replace_all(site,"PG-O-","PICO"),
#            site = str_replace_all(site,"PG-C-","PICR"),
#            site = str_replace_all(site,"RU-C-","RRRR"),
#            site = as_factor(site)
#            )
# }
#
# future::plan(multisession,workers = availableCores()-2)
#
# t <- future_map(f,collapse,.progress = TRUE) %>% list_rbind()
#
# dir <- "data/humidity sensors"
# f <- list.files(dir,pattern = "csv",full.names = TRUE)
#
# collapse <- function(file_name) {
#   # Given relative file name of temperature log from one logger, return all rows.
#   if(read_csv(file_name) %>% problems() %>% nrow() > 100){
#     t <- (readLines(file_name) %>% tibble() %>% slice(2:n()) %>% pull() %>%
#             I() %>% read_csv(col_names = c("x1","x2","x3")))[,2:4]
#   } else {
#     t <- read_csv(file_name)[,2:4]
#   }
#   logger_name <- str_extract(file_name,paste0("(?<=^",dir,"/)[:graph:]+(?=\\s)"))
#   names(t) <- c("datetime","temp","humidity")
#   t$datetime <- mdy_hms(t$datetime, tz = "America/Los_Angeles")
#   t$date <- t$datetime %>% date()
#   t <- t %>%
#     mutate(box = str_extract(logger_name,"[:graph:]*(?=-[:alpha:])"),
#            humidity = as.numeric(humidity),
#            site = factor(str_extract(box,"([:alpha:]|[:punct:])+")), # extract site from nestbox
#            habitat = fct_collapse(site, # collapse sites into habitat classes
#                                   Orchard = c("BRO","MCE","PICO","MCO","FBFO"),
#                                   Grassland = c("MBNG","RRRG","PICG","FBFG"),
#                                   `Row crop` = c("MBNR","PICR","RRRR","FBFR"),
#                                   Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","RRRC")),
#            habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")),
#            site = str_replace_all(site,"MCE","MCO"), # Collapse site vectors when there are multiple site abbreviations per actual site
#            site = str_replace_all(site,"PCT","PCC"),
#            site = str_replace_all(site,"PCE","PCC"),
#            site = str_replace_all(site,"RRE","RRC"),
#            site = str_replace_all(site,"SFP","MBNC"),
#            site = str_replace_all(site,"PG-O-","PICO"),
#            site = str_replace_all(site,"PG-C-","PICR"),
#            site = str_replace_all(site,"RU-C-","RRRR"),
#            site = as_factor(site))
# }
#
# future::plan(multisession,workers = availableCores()-2)
#
# h <- future_map(f,collapse,.progress = TRUE) %>% list_rbind()
#
#
# h_grouped <- h %>%
#   mutate(site = as.character(site),
#          datetime_hour = hour(datetime),
#          datetime_min = minute(datetime),
#          datetime_min = if_else(datetime_min < 30,0,30)) %>%
#   group_by(date,box,site,habitat,datetime_hour,datetime_min) %>%
#   summarize(meanh_30min = mean(humidity,na.rm = TRUE))
# t_grouped <- t %>%
#   mutate(site = as.character(site)) %>%
#   arrange(box,datetime) %>%
#   mutate(datetime_hour = hour(datetime),
#          datetime_min = minute(datetime),
#          datetime_min = if_else(datetime_min < 30,0,30)) %>%
#   group_by(date,box,site,habitat,logger_position,datetime_hour,datetime_min) %>%
#   summarize(meant_30min = mean(temp,na.rm = TRUE))
#
# ## Identify replacements for bad data
#
# ### date, box, site, position combinations with anomalously high data
# match <- dplyr::filter(t_grouped %>% ungroup(),meant_30min>55) %>% dplyr::select(date,box,site,datetime_hour,datetime_min,logger_position) %>% distinct()
#
# ### candidates to replace
# matched <- rowwise(match) %>% group_map(~ dplyr::filter(t_grouped %>% ungroup(),date == .x$date,site == .x$site,datetime_hour == .x$datetime_hour, datetime_min == .x$datetime_min, logger_position == "I"))
#
# box_match <- tibble(old = pull(match,box) %>% unique(),new = c("RRRG1","RRC12","MBNG6"))
#
# match <- match %>%
#   left_join(box_match,by = join_by(box == old)) %>%
#   left_join(t_grouped %>% ungroup() %>% dplyr::select(meant_30min,datetime_hour,datetime_min,logger_position,date,box),
#             by = join_by(new == box,datetime_hour == datetime_hour,datetime_min == datetime_min, date == date,logger_position == logger_position))
#
# t_grouped <- t_grouped %>% ungroup() %>%
#   left_join(match %>% dplyr::select(-c(site,new)),by = join_by(box == box,datetime_hour == datetime_hour,datetime_min == datetime_min, date == date,logger_position == logger_position)) %>%
#   mutate(meant_30min = if_else(is.na(meant_30min.y),meant_30min.x,meant_30min.y)) %>%
#   dplyr::select(-c(meant_30min.x,meant_30min.y))
#
# t_h <- left_join(t_grouped %>% ungroup(),
#                  h_grouped %>% ungroup(),
#                  by = join_by(date,site,datetime_hour,datetime_min),
#                  na_matches = "never",
#                  multiple = "first") %>%
#   dplyr::select(-c(box.y,habitat.y)) %>%
#   rename(box = box.x,
#          habitat = habitat.x)
#
#
# t_h_hi <- t_h %>%
#   mutate(hi_30min = heat.index(t = meant_30min,rh = meanh_30min,temperature.metric = "celsius",output.metric = "celsius",round = 2),
#          over_40C = if_else(meant_30min >= 40,1,0),
#          over_30C = if_else(meant_30min >= 30,1,0),
#          over_35C = if_else(meant_30min >= 35,1,0),
#          over_75hi = if_else(hi_30min >= 75,1,0),
#          over_50hi = if_else(hi_30min >= 50,1,0),
#          over_30hi = if_else(hi_30min >= 30,1,0),
#          over_35hi = if_else(hi_30min >= 35,1,0),
#          deg_diff_from_35C = abs(35-meant_30min),
#          deg_sq_diff_from_35C = (abs(35-meant_30min))^2,
#          hi_diff_from_35C = abs(35-hi_30min),
#          hi_sq_diff_from_35C = (abs(35-hi_30min))^2,
#          deg_diff_from_30C = abs(30-meant_30min),
#          deg_sq_diff_from_30C = (abs(30-meant_30min))^2,
#          hi_diff_from_30C = abs(30-hi_30min),
#          hi_sq_diff_from_30C = (abs(30-hi_30min))^2,
#          box = str_replace(box,"PG-C-0","PICR") %>%
#                   str_replace("PG-O-0","PICO") %>%
#                   str_replace("RU-C-0","RRRR"))
#
# # Output unsummarized heat index measures for provisioning
# nrow(t_h_hi)
#
# write_rds(t_h_hi[1:400000,],"data/heatindex_30min_1.rds")
# write_rds(t_h_hi[400001:nrow(t_h_hi),],"data/heatindex_30min_2.rds")
#
#
# t_h_hi_sum <- t_h_hi %>%
#   dplyr::summarize(maxt = max(meant_30min,na.rm = TRUE),
#             mint = min(meant_30min,na.rm = TRUE),
#             meant = mean(meant_30min,na.rm = TRUE),
#             maxh = max(meanh_30min,na.rm = TRUE),
#             minh = min(meanh_30min,na.rm = TRUE),
#             meanh = mean(meanh_30min,na.rm = TRUE),
#             maxhi = max(hi_30min,na.rm = TRUE),
#             minhi = min(hi_30min,na.rm = TRUE),
#             meanhi = mean(hi_30min,na.rm = TRUE),
#             degreehours_over_40 = sum((meant_30min * over_40C) * .5,na.rm = TRUE),
#             degreehours_over_30 = sum((meant_30min * over_30C) * .5,na.rm = TRUE),
#             degreehours_over_35 = sum((meant_30min * over_35C) * .5,na.rm = TRUE),
#             hihours_over_75 =  sum((hi_30min * over_75hi) * .5,na.rm = TRUE),
#             hihours_over_50 =  sum((hi_30min * over_50hi) * .5,na.rm = TRUE),
#             hihours_over_30 =  sum((hi_30min * over_30hi) * .5,na.rm = TRUE),
#             hihours_over_35 =  sum((hi_30min * over_35hi) * .5,na.rm = TRUE),
#             deghours_diff_35C = sum(deg_diff_from_35C * .5,na.rm = TRUE),
#             deghours_sqdiff_35C = sum(deg_sq_diff_from_35C * .5,na.rm = TRUE),
#             hihours_diff_35C = sum(hi_diff_from_35C * .5,na.rm = TRUE),
#             hihours_sqdiff_35C = sum(hi_sq_diff_from_35C * .5,na.rm = TRUE),
#             deghours_diff_30C = sum(deg_diff_from_30C * .5,na.rm = TRUE),
#             deghours_sqdiff_30C = sum(deg_sq_diff_from_30C * .5,na.rm = TRUE),
#             hihours_diff_30C = sum(hi_diff_from_30C * .5,na.rm = TRUE),
#             hihours_sqdiff_30C = sum(hi_sq_diff_from_30C * .5,na.rm = TRUE),
#             .by = c(date,box,site,habitat,logger_position)) %>%
#   distinct() %>%
#   mutate(box = str_replace(box,"PG-C-0","PICR") %>%
#            str_replace("PG-O-0","PICO") %>%
#            str_replace("RU-C-0","RRRR"),
#          maxh = if_else(is.infinite(maxh),NA_integer_,maxh),
#          minh = if_else(is.infinite(minh),NA_integer_,minh),
#          meanh = if_else(is.nan(meanh),NA_integer_,meanh),
#          maxhi = if_else(is.infinite(maxhi),NA_integer_,maxhi),
#          minhi = if_else(is.infinite(minhi),NA_integer_,minhi),
#          meanhi = if_else(is.nan(meanhi),NA_integer_,meanhi),
#          degreehours_over_40 = if_else(is.nan(degreehours_over_40) | is.na(maxt),NA_integer_,degreehours_over_40),
#          degreehours_over_30 = if_else(is.nan(degreehours_over_30) | is.na(maxt),NA_integer_,degreehours_over_30),
#          degreehours_over_35 = if_else(is.nan(degreehours_over_35) | is.na(maxt),NA_integer_,degreehours_over_35),
#          hihours_over_75 = if_else(is.nan(hihours_over_75) | is.na(maxh),NA_integer_,hihours_over_75),
#          hihours_over_50 = if_else(is.nan(hihours_over_50) | is.na(maxh),NA_integer_,hihours_over_50),
#          hihours_over_30 = if_else(is.nan(hihours_over_30) | is.na(maxh),NA_integer_,hihours_over_30),
#          hihours_over_35 = if_else(is.nan(hihours_over_35) | is.na(maxh),NA_integer_,hihours_over_35),
#          deghours_diff_35C = if_else(is.nan(deghours_diff_35C) | is.na(maxt),NA_integer_,deghours_diff_35C),
#          deghours_sqdiff_35C = if_else(is.nan(deghours_sqdiff_35C) | is.na(maxt),NA_integer_,deghours_sqdiff_35C),
#          hihours_diff_35C = if_else(is.nan(hihours_diff_35C) | is.na(maxh),NA_integer_,hihours_diff_35C),
#          hihours_sqdiff_35C = if_else(is.nan(hihours_sqdiff_35C) | is.na(maxh),NA_integer_,hihours_sqdiff_35C),
#          deghours_diff_30C = if_else(is.nan(deghours_diff_30C) | is.na(maxt),NA_integer_,deghours_diff_30C),
#          deghours_sqdiff_30C = if_else(is.nan(deghours_sqdiff_30C) | is.na(maxt),NA_integer_,deghours_sqdiff_30C),
#          hihours_diff_30C = if_else(is.nan(hihours_diff_30C) | is.na(maxh),NA_integer_,hihours_diff_30C),
#          hihours_sqdiff_30C = if_else(is.nan(hihours_sqdiff_30C) | is.na(maxh),NA_integer_,hihours_sqdiff_30C)
#          )
#
# write_rds(t_h_hi_sum,"data/heatindex.rds")

# future::plan(multisession,workers = availableCores()-2)
#
# join_and_write <- function(x,idx){
#   h_sub <- dplyr::filter(h,
#                          site == site,
#                          date(datetime) == date(x$datetime))
#   gc()
#   distance_left_join(x, h_sub,
#                   by = c("site" = "site",
#                          'datetime' = "h_datetime_minus2min",
#                          'datetime' = "h_datetime_plus2min"),
#                   match_fun = list(`==`, `>`, `<`)) %>%
#     write_rds(paste0("data/temp_and_humidity/slice",idx,".rds"))
#   print(paste0("rownum: ",idx))
#   gc()
# }
#
# t %>% future_iwalk(\(x,idx) join_and_write(x,idx),
#                 .progress = TRUE
#                 )

g <- read_rds("data/growth_and_provis_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(abs_change_cort = cort_s2-cort_s1,
         prop_change_cort = (cort_s2-cort_s1)/cort_s1,
         year_fct = as.factor(year),
         site_box = factor(str_extract(Nestbox,"([:alpha:]|[:punct:])+")))
p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

c <- read_rds("data/cort_with_id.rds")

## add day prior temp to g

### Access temp object; filter to make the mapping go faster
hi <- read_rds("data/heatindex.rds")
histtmax <- read_rds("data/cimis_1883_dailymax.rds")
# histtmeanmax <- read_rds("data/cimis_1883_dailymeanmax.rds")
future::plan(multisession,workers = availableCores()-2)

filterhi_dayprior <- future_map2(pull(g,date),pull(g,Nestbox),function(x,y){filter(hi,
                                                                     date == x - days(1),
                                                                     box == y,
                                                                     logger_position == "I")})
filterhi_week <- future_map2(pull(g,interval),pull(g,Nestbox),function(x,y){filter(hi,
                                                                               date %within% x,
                                                                               box == y,
                                                                               logger_position == "I")})

# With the stipulation that the average should be over a 2 week period centered on the date, the prior day and prior week historical max should be the same for each record.

# filterhisttday_dayprior <- future_map2(pull(g,date),pull(g,Nestbox),function(x,y){filter(histtmeanmax,
#                                                                                          month == month(x),
#                                                                                          day == mday(x)-1)})

filterhisttday_week <- future_map2(pull(g,interval),pull(g,Nestbox),function(x,y){
  s <- tibble(dates = seq(int_end(x)-days(7),int_end(x)+days(7),by = "days")) %>%
    mutate(month = month(dates),
           day = mday(dates))
  filter(histtmax,
         month(date) %in% (pull(s,month) %>% unique()),
         mday(date) %in% (pull(s,day) %>% unique())) %>%
    mutate(year = year(date)) %>%
    summarize(meanmax_year = mean(maxt,na.rm = TRUE),.by = year) %>%
    summarize(meanmax = mean(meanmax_year,na.rm = TRUE),
              sdmax = sd(meanmax_year,na.rm = TRUE))
  })

### For each date - 1 day and nestbox in g, pull mean, max, min, and timeover40.

g <- g %>% mutate(
  maxt_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(maxt) %>% mean(na.rm = TRUE)}),
  mint_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(mint) %>% mean(na.rm = TRUE)}),
  meant_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(meant) %>% mean(na.rm = TRUE)}),
  maxhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(maxhi) %>% mean(na.rm = TRUE)}),
  minhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(minhi) %>% mean(na.rm = TRUE)}),
  meanhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(meanhi) %>% mean(na.rm = TRUE)}),
  maxhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(maxhi) %>% mean(na.rm = TRUE)}),
  minhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(minhi) %>% mean(na.rm = TRUE)}),
  meanhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(meanhi) %>% mean(na.rm = TRUE)}),
  hihours_over_75hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_75) %>% mean(na.rm = TRUE)}),
  hihours_over_75hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_75) %>% sum()}),
  degreehours_over_40C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_40) %>% mean(na.rm = TRUE)}),
  degreehours_over_40C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_40) %>% sum()}),
  hihours_over_50hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_50) %>% mean(na.rm = TRUE)}),
  hihours_over_50hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_50) %>% sum()}),
  hihours_over_30hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_30) %>% mean(na.rm = TRUE)}),
  hihours_over_30hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_30) %>% sum()}),
  degreehours_over_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_30) %>% mean(na.rm = TRUE)}),
  degreehours_over_30C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_30) %>% sum()}),
  degreehours_over_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_35) %>% mean(na.rm = TRUE)}),
  degreehours_over_35C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_35) %>% sum()}),
  deghours_diff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_diff_35C) %>% mean(na.rm = TRUE)}),
  deghours_sqdiff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_sqdiff_35C) %>% mean(na.rm = TRUE)}),
  hihours_diff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_diff_35C) %>% mean(na.rm = TRUE)}),
  hihours_sqdiff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_sqdiff_35C) %>% mean(na.rm = TRUE)}),
  deghours_diff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_diff_30C) %>% mean(na.rm = TRUE)}),
  deghours_sqdiff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_sqdiff_30C) %>% mean(na.rm = TRUE)}),
  hihours_diff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_diff_30C) %>% mean(na.rm = TRUE)}),
  hihours_sqdiff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_sqdiff_30C) %>% mean(na.rm = TRUE)}),
  histmaxt_priorweek_mean = filterhisttday_week %>%
    map_dbl(function(x){x %>% pull(meanmax) %>% mean(na.rm = TRUE)}),
  histmaxt_priorweek_sd = filterhisttday_week %>%
    map_dbl(function(x){x %>% pull(sdmax) %>% mean(na.rm = TRUE)}),
  maxt_priorday_minus_histmeanmax = maxt_prior - histmaxt_priorweek_mean,
  maxt_priorday_histzscale = maxt_priorday_minus_histmeanmax/histmaxt_priorweek_sd,
  maxt_priorweek_minus_histmeanmax = meanmaxtempI - histmaxt_priorweek_mean,
  maxt_priorweek_histzscale = maxt_priorweek_minus_histmeanmax/histmaxt_priorweek_sd,
  across(maxt_prior:maxt_priorweek_histzscale,~if_else(is.nan(.x),NA_integer_,.x))
)

# plot(provis_mean ~ hihours_diff_35C_priorday, data = g)


write_rds(g,"data/growth_cort_provis_manytempmeasures.rds")


p <- read_rds("data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds") %>%
  mutate(year = year(date),
         year_fct = as.factor(year))

## pull in heat index, harmonize with needs of match_temp function
hi <- read_rds("data/heatindex_30min_1.rds") %>% bind_rows(read_rds("data/heatindex_30min_2.rds")) %>%
  mutate(year = year(date),
         datetime_conv = paste0(as.character(date)," ",as.character(datetime_hour),":",
                                if_else(as.character(datetime_min) == "0","00",as.character(datetime_min)),":00"),
         datetime_conv = ymd_hms(datetime_conv),
         datetime_conv = force_tz(datetime_conv,tz = "America/Los_Angeles")
  )

match_temp <- function(vid_data,temp_data){

  # testing
  # vid_data <- instant_temp[[1]] #%>%
  #mutate(realtime_num = as.POSIXct(realtime) %>% as.numeric())
  # temp_data <- filter(temp_data,
  # 			 year == unique(vid_data$year)[1],
  # 			 box == unique(vid_data$nestbox)[1],
  # 			 logger_position == "O") %>%
  # 	select(datetime,temp) #%>%
  #mutate(datetime_num = as.POSIXct(datetime) %>% as.numeric())

  # str(vid_data)
  #str(temp_data)
  #names(temp_data)
  #names(vid_data)

  # function

  temp_data <- temp_data %>%
    dplyr::select(datetime_conv,hi_30min)

  # out <- vid_data %>%
  # 	left_join(vid_data,
  # 						temp_data,
  # 						by = join_by(closest(realtime <= datetime_conv)))
  #
  out <- difference_left_join(
    vid_data, temp_data,
    by = c("start" = "datetime_conv"),
    max_dist = as.difftime(31, units = "mins"),
    distance_col = "time_difference") %>%
    group_by(start) %>%
    slice_min(order_by = time_difference,n = 1) %>% ungroup()

}

future::plan(multisession,workers = availableCores()-2)
instant_temp <- p %>%
  group_by(year,nestbox) %>%
  group_split() %>%
  future_map(~ match_temp(.x,filter(hi,
                                    year == unique(.x$year)[1],
                                    box == unique(.x$nestbox)[1],
                                    logger_position == "O")),
             .progress = TRUE) %>%
  list_rbind()

## Add prior day deg & hi difference from 35C, also sum of squared differences

temp <- read_rds("data/heatindex.rds")

filterhi_dayprior <- future_map2(pull(instant_temp,date),pull(instant_temp,nestbox),function(x,y){filter(temp,
                                                                                   date == x - days(1),
                                                                                   box == y,
                                                                                   logger_position == "O")})

filterhi_week <- future_map2(pull(instant_temp,date),pull(instant_temp,nestbox),function(x,y){filter(temp,
                                                                                   date %within% interval(x - days(1),x - days(7)),
                                                                                   box == y,
                                                                                   logger_position == "O")})

p <- instant_temp %>% mutate(
  maxt_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(maxt) %>% mean(na.rm = TRUE)}),
  mint_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(mint) %>% mean(na.rm = TRUE)}),
  meant_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(meant) %>% mean(na.rm = TRUE)}),
  maxt_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(maxt) %>% mean(na.rm = TRUE)}),
  mint_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(mint) %>% mean(na.rm = TRUE)}),
  meant_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(meant) %>% mean(na.rm = TRUE)}),
  maxhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(maxhi) %>% mean(na.rm = TRUE)}),
  minhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(minhi) %>% mean(na.rm = TRUE)}),
  meanhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(meanhi) %>% mean(na.rm = TRUE)}),
  maxhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(maxhi) %>% mean(na.rm = TRUE)}),
  minhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(minhi) %>% mean(na.rm = TRUE)}),
  meanhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(meanhi) %>% mean(na.rm = TRUE)}),
  hihours_over_75hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_75) %>% mean(na.rm = TRUE)}),
  hihours_over_75hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_75) %>% sum()}),
  degreehours_over_40C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_40) %>% mean(na.rm = TRUE)}),
  degreehours_over_40C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_40) %>% sum()}),
  hihours_over_50hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_50) %>% mean(na.rm = TRUE)}),
  hihours_over_50hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_50) %>% sum()}),
  hihours_over_30hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_30) %>% mean(na.rm = TRUE)}),
  hihours_over_30hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_30) %>% sum()}),
  degreehours_over_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_30) %>% mean(na.rm = TRUE)}),
  degreehours_over_30C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_30) %>% sum()}),
  degreehours_over_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_35) %>% mean(na.rm = TRUE)}),
  degreehours_over_35C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_35) %>% sum()}),
  deghours_diff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_diff_35C) %>% mean(na.rm = TRUE)}),
  deghours_sqdiff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_sqdiff_35C) %>% mean(na.rm = TRUE)}),
  hihours_diff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_diff_35C) %>% mean(na.rm = TRUE)}),
  hihours_sqdiff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_sqdiff_35C) %>% mean(na.rm = TRUE)}),
  deghours_diff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_diff_30C) %>% mean(na.rm = TRUE)}),
  deghours_sqdiff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_sqdiff_30C) %>% mean(na.rm = TRUE)}),
  hihours_diff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_diff_30C) %>% mean(na.rm = TRUE)}),
  hihours_sqdiff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_sqdiff_30C) %>% mean(na.rm = TRUE)}),
  across(maxt_prior:hihours_sqdiff_30C_priorday,~if_else(is.nan(.x),NA_integer_,.x))
)

write_rds(p,"data/provis_manytempmeasures.rds")

c <- read_rds("data/cort_with_id.rds") %>%
  dplyr::filter(Age %in% c("AHY","ASY","SY"),Sex == "F") %>%
  rename(assay_date = "date.x",
         sample_date = "date.y") %>%
  mutate(year_fct = as.factor(year))

str(c)


hi <- read_rds("data/heatindex.rds")
library(future)
library(furrr)
histtmax <- read_rds("data/cimis_1883_dailymax.rds")

future::plan(multisession,workers = availableCores()-2)

filterhi_dayprior <- future_map2(pull(c,sample_date),pull(c,Nestbox),function(x,y){filter(hi,
                                                                                          date == x - days(1),
                                                                                          box == y,
                                                                                          logger_position == "O")})
filterhi_week <- future_map2(pull(c,sample_date),pull(c,Nestbox),function(x,y){filter(hi,
                                                                                      date %within% interval(x-days(7),x),
                                                                                      box == y,
                                                                                      logger_position == "O")})

filterhisttday_week <- future_map2(pull(c,sample_date),pull(c,Nestbox),function(x,y){
  s <- tibble(dates = seq(x-days(7),x+days(7),by = "days")) %>%
    mutate(month = month(dates),
           day = mday(dates))
  filter(histtmax,
         month(date) %in% (pull(s,month) %>% unique()),
         mday(date) %in% (pull(s,day) %>% unique())) %>%
    mutate(year = year(date)) %>%
    summarize(meanmax_year = mean(maxt,na.rm = TRUE),.by = year) %>%
    summarize(meanmax = mean(meanmax_year,na.rm = TRUE),
              sdmax = sd(meanmax_year,na.rm = TRUE))
})

### For each date - 1 day and nestbox in c, pull mean, max, min, and timeover40.

c <- c %>% mutate(
  maxt_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(maxt) %>% mean(na.rm = TRUE)}),
  mint_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(mint) %>% mean(na.rm = TRUE)}),
  meant_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(meant) %>% mean(na.rm = TRUE)}),
  maxt_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(maxt) %>% mean(na.rm = TRUE)}),
  mint_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(mint) %>% mean(na.rm = TRUE)}),
  meant_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(meant) %>% mean(na.rm = TRUE)}),
  maxhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(maxhi) %>% mean(na.rm = TRUE)}),
  minhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(minhi) %>% mean(na.rm = TRUE)}),
  meanhi_prior = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(meanhi) %>% mean(na.rm = TRUE)}),
  maxhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(maxhi) %>% mean(na.rm = TRUE)}),
  minhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(minhi) %>% mean(na.rm = TRUE)}),
  meanhi_week = filterhi_week %>%
    map_dbl(function(x){x %>% pull(meanhi) %>% mean(na.rm = TRUE)}),
  hihours_over_75hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_75) %>% mean(na.rm = TRUE)}),
  hihours_over_75hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_75) %>% sum()}),
  degreehours_over_40C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_40) %>% mean(na.rm = TRUE)}),
  degreehours_over_40C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_40) %>% sum()}),
  hihours_over_50hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_50) %>% mean(na.rm = TRUE)}),
  hihours_over_50hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_50) %>% sum()}),
  hihours_over_30hi_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_over_30) %>% mean(na.rm = TRUE)}),
  hihours_over_30hi_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(hihours_over_30) %>% sum()}),
  degreehours_over_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_30) %>% mean(na.rm = TRUE)}),
  degreehours_over_30C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_30) %>% sum()}),
  degreehours_over_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(degreehours_over_35) %>% mean(na.rm = TRUE)}),
  degreehours_over_35C_priorweek = filterhi_week %>%
    map_dbl(function(x){x %>% pull(degreehours_over_35) %>% sum()}),
  deghours_diff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_diff_35C) %>% mean(na.rm = TRUE)}),
  deghours_sqdiff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_sqdiff_35C) %>% mean(na.rm = TRUE)}),
  hihours_diff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_diff_35C) %>% mean(na.rm = TRUE)}),
  hihours_sqdiff_35C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_sqdiff_35C) %>% mean(na.rm = TRUE)}),
  deghours_diff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_diff_30C) %>% mean(na.rm = TRUE)}),
  deghours_sqdiff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(deghours_sqdiff_30C) %>% mean(na.rm = TRUE)}),
  hihours_diff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_diff_30C) %>% mean(na.rm = TRUE)}),
  hihours_sqdiff_30C_priorday = filterhi_dayprior %>%
    map_dbl(function(x){x %>% pull(hihours_sqdiff_30C) %>% mean(na.rm = TRUE)}),
  histmaxt_priorweek_mean = filterhisttday_week %>%
    map_dbl(function(x){x %>% pull(meanmax) %>% mean(na.rm = TRUE)}),
  histmaxt_priorweek_sd = filterhisttday_week %>%
    map_dbl(function(x){x %>% pull(sdmax) %>% mean(na.rm = TRUE)}),
  maxt_priorday_minus_histmeanmax = maxt_prior - histmaxt_priorweek_mean,
  maxt_priorday_histzscale = maxt_priorday_minus_histmeanmax/histmaxt_priorweek_sd,
  maxt_priorweek_minus_histmeanmax = meanmaxtempI - histmaxt_priorweek_mean,
  maxt_priorweek_histzscale = maxt_priorweek_minus_histmeanmax/histmaxt_priorweek_sd,
  across(maxt_prior:maxt_priorweek_histzscale,~if_else(is.nan(.x),NA_integer_,.x))
)

# plot(provis_mean ~ hihours_over_30hi_priorweek, data = g)


write_rds(c,"data/cort_manytempmeasures.rds")
