##### Generate comp log from nest checks
##### Author: Katherine Lauck

require(tidyverse)
require(lubridate)
require(DescTools)
source('code/helper_functions.R')

identify_end <- function(checks,nest_start) {
  ### Identify the end of an attempt, given a series of checks and the index of the beginning of the attempt.
  # First end is one of:
  # Skip (NA) after first eggs or nestlings
  # 0 nestlings and 0 eggs after first eggs
  # 0 nestlings after first nestlings
  # Decrease and then increase in eggs
  # Decrease in eggs with no nestlings, then no nestlings after 3 checks
  # No nestlings at fifth check of eggs
  # There are more than 21 days between two checks

  for(i in seq(nest_start,nrow(checks))){
    if(i == nrow(checks)){
      nest_end <- i
      break
    }
    if((checks$Eggs[i] > 0 | checks$Nestlings[i] > 0) &
       is.na(checks$Nestlings[i+1])) {
      nest_end <- min(i + max(which(is.na(checks$Eggs[(i+1):(i+3)]))),nrow(checks))
      break
    } else if((is.na(checks$Nestlings[i]) | checks$Nestlings[i] > 0) &
              identical(checks$Nestlings[i+1],0)){
      nest_end <- i + 1
      break
    } else if(checks$Eggs[i] > 0 &
              identical(checks$Eggs[i+1],0) &
              identical(checks$Nestlings[i+1],0)){
      nest_end <- i + 1
      break
    } else if(checks$Eggs[i] > 0 &
              checks$Eggs[i+1] < checks$Eggs[i] &
              any(checks$Eggs[(i+2):(i+3)] > checks$Eggs[i+1],na.rm = TRUE)){
      nest_end <- i + 1
      break
    } else if(checks$Eggs[i] > 0 &
              checks$Eggs[i+1] < checks$Eggs[i] &
              checks$Nestlings[i+1] < 1 &
              all(identical(checks$Nestlings[(i+2):(i+3)],0),na.rm = TRUE)){
      nest_end <- i + 1
      break
    } else if(checks$Eggs[i] > 0 &
              all(identical(checks$Nestlings[(i+1):(i+3)],0),na.rm = TRUE)){
      nest_end <- i + 1
      break
    } else if((as_date(checks$Date[i+1]) - as_date(checks$Date[i])) %>% as.numeric() > 21 &
              !is.na(as_date(checks$Date[i+1]))){
      nest_end <- i
      break
    } else {nest_end <- NA;next}
  }

  return(nest_end)

}

identify_start <- function(checks, end){
  ### Identify the start of an attempt, given a series of checks and the index of the end of the last attempt.
  start <- end + min(which(checks$Eggs[(end + 1):(nrow(checks))] > 0),na.rm = TRUE)
  return(start)
}

collapse <- function(checks,...){
  ### Collapse a set of checks from one nest collected as part of the Nestbox Highway Project into information about each attempted nest with eggs constructed, one attempt per line.

  checks <- arrange(checks,Date) # order checks by date to allow parsing of attempts

  start <- 0
  end <- 0

  while(!is.na(end[length(end)]) & end[length(end)] < nrow(checks)){
    first_egg <- identify_start(checks,end[length(end)])
    start <- c(start,first_egg)
    if(identical(start[length(start)],Inf)) {end <- c(end,NA);break}
    nest_end <- identify_end(checks,first_egg)
    end <- c(end,nest_end)
  }

  attempt_map <- tibble(attempt = seq(0,length(start)-1),start = start, end = end) %>%
    filter(start != Inf,attempt != 0)

  attempt <- rep(NA,times = nrow(checks))

  if(nrow(attempt_map) == 0) {} else {
    if(is.na(attempt_map$end[nrow(attempt_map)])) {attempt_map$end[nrow(attempt_map)] <- nrow(checks)}
    for(i in attempt_map$attempt){
      attempt[attempt_map$start[i]:attempt_map$end[i]] <- i
    }
  }

  checks <- mutate(checks, attempt = attempt) %>%
    filter(!is.na(attempt)) %>%
    group_by(attempt)

  if(nrow(checks) == 0){
    hobo_logger = NA
    site = NA
    box = NA
    attempt = NA
    species = NA
    male = NA
    female = NA
    first_egg_date = NA
    first_egg_certainty = NA
    clutch_size = NA
    incubation_begins = NA
    eggs_hatched = NA
    eggs_failed = NA
    eggs_abandoned = NA
    eggs_eaten = NA
    eggs_destroyed = NA
    egg_unk = NA
    brood_size = NA
    hatch_date = NA
    hd_certainty = NA
    early_banding = NA
    late_banding = NA
    band_date = NA
    chicks_fledged = NA
    nest_perish = NA
    nest_abandoned = NA
    nest_pred = NA
    band_nums = NA
  } else {

    eggs_unhatched <- function(checks,...){
      if(max(checks$Nestlings,na.rm = TRUE) == 0){return(0)
        } else if(any(!is.na(checks$Eggs[which(checks$Nestlings == max(checks$Nestlings,
                                                                   na.rm = TRUE))]))) {
        return(checks$Eggs[which(checks$Nestlings == max(checks$Nestlings,na.rm = TRUE))] %>%
                 min(na.rm = TRUE)
        )
      } else {return(NA)}
    }

    eggs_aband <- function(checks,...){
      if(max(checks$Nestlings,na.rm = TRUE) == 0){return(checks$Eggs[nrow(checks)])} else {return(0)}
    }

    wing_tbl <- tibble(
      age = c(rep(seq(0,18,1),3),0,0),
      species = c(rep(c("WEBL","TRES","ATFL"),each = 19),"skipped","Skipped"),
      wing = c(6,7,8,9,11,13,16,19,23,27,31,35,39,43,47,51,55,59,63,6,6,7,8,10,13,16,20,25,30,35,41,47,52,57,62,66,72,78,8,8,9,12,14,17,21,25,32,38,44,50,55,60,65,70,75,77,79,0,0)
    ) %>% group_by(species)

    nest_age <- function(...,wing=wing_tbl) {
      x <- tibble(...)
      if(!is.na(x$`Wing chords`) & !is.null(x$`Wing chords` %>% unlist()) &
         x$Species %in% wing_tbl$species) {
        mean_wing <- str_extract_all(x$`Wing chords`,"[[:digit:]]+") %>%
          unlist() %>%
          # str_split(",") %>%
          # map(function(x){if(identical(x,character(0))){NA} else {x}}) %>%
          # map(function(x){if(!is.na(x)){as.numeric(x)} else {NA}}) %>%
          as.numeric() %>%
          mean(na.rm = TRUE)

        age <- wing %>% filter(abs(wing-mean_wing)==min(abs(wing-mean_wing))) %>%
          # map2(.y = checks$Species,.f = function(x,y){filter(x,species == y)}) %>%
          filter(species == x$Species,age == min(age)) %>%
          # ungroup() %>%
          # map(function(x){select(x,age)}) %>%
          # map(function(x){return(x[1,])}) %>%
          pull(age) %>%
          # dplyr::bind_rows() %>%
          # pull(age) %>%
          as.integer()
      } else {age <- NA_integer_}

      return(age)
    }

    count_fledged <- function(checks,...) {
      active_checks <- filter(checks, Nestlings > 1)
      if(nrow(active_checks) < 2){
        fledged = NA
      } else if(pmap(active_checks,nest_age) %>% unlist() %>% is.na() %>% all()) {
        fledged = NA
      } else {
        ages <- pmap(active_checks,nest_age) %>% unlist()
        hatch_date <- active_checks$Date[which.min(ages)] %>% as_date() - min(ages,na.rm = TRUE)
        if(as_date(max(active_checks$Date,na.rm = TRUE)) - hatch_date > 7) {
          fledged <- active_checks$Nestlings[which.max(active_checks$Date)]
          dead_chicks <- sum(checks$dead_nestlings,na.rm = TRUE)
          fledged <- fledged[1] - dead_chicks[1]
        } else {fledged <- NA}
      }
      return(fledged)
    }

    get_hatch_date <- function(checks,...){
      active_checks <- filter(checks, Nestlings > 1)
      if(pmap(active_checks,nest_age) %>% unlist() %>% is.na() %>% all()){
        hatch_date = NA
      } else {
        ages <- pmap(active_checks,nest_age) %>% unlist()
        hatch_date <- active_checks$Date[which.min(ages)] %>% as_date() - min(ages,na.rm = TRUE)
      }
      return(hatch_date)
    }

    get_hatch_date_uncertainty <- function(checks,...){
      active_checks <- filter(checks, Nestlings > 1)
      if(pmap(active_checks,nest_age) %>% unlist() %>% is.na() %>% all()){
        hatch_date = NA
        hd_certainty = NA
      } else {
        ages <- pmap(active_checks,nest_age) %>% unlist()
        hatch_date <- active_checks$Date[which.min(ages)] %>% as_date() - min(ages,na.rm = TRUE)
        hd_certainty = ifelse(min(ages,na.rm = TRUE) < 1,
                              'Confirmed',
                              'Estimated')
      }
      return(hd_certainty)
    }

    get_early_banding_date <- function(checks,...){
      hd <- get_hatch_date(checks)
      if(nrow(checks) < 1){
        early_banding <- NA
      } else {
        if((checks$Species %>% unique())[1] == 'WEBL'){
          early_banding <- as_date(hd) + 13
        } else if((checks$Species %>% unique())[1] == 'TRES'){
          early_banding <- as_date(hd) + 11
        } else if((checks$Species %>% unique())[1] == 'ATFL'){
          early_banding <- as_date(hd) + 9
        } else if ((checks$Species %>% unique())[1] == 'HOWR'){
          early_banding <- as_date(hd) + 8
        }
      }
      return(early_banding)
    }

    get_species <- function(checks,...){
      select(checks,Species) %>% slice(1) %>% as.character()
    }

    get_logger_and_id <- function(int,nestbox){
      g <- read_rds('data/morph.rds') %>%
        filter(as_date(`Banding Date`) %within% int,
               Nestbox == nestbox)
      t <- read_rds('data/temp.rds') %>%
        filter(date %within% int,
               box == nestbox)
      hobo_logger <- if_else(nrow(t) > 0,
                             'HOBO loggers',
                             NA)
      fe <- g %>%
        filter(Sex == 'F',Age %in% c('AHY','SY','ASY'))
      female <- if_else(nrow(fe) == 0,NA,fe$`Band Number`[1])
      ma <- g %>%
        filter(Sex == 'M',Age %in% c('AHY','SY','ASY'))
      male <- if_else(nrow(ma) == 0,NA,ma$`Band Number`[1])
      band_date <- g %>% filter(Age == 'L',
                                `Band Number` != 'NA') %>%
        slice(1) %>% pull(`Banding Date`) %>% as_date()
      band_nums <- g %>% filter(Age == 'L',
                                `Band Number` != 'NA') %>%
        pull(`Band Number`) %>% unique() %>%
        paste(collapse = ', ',sep = ', ')

      d <- tibble(hobo_logger = hobo_logger,
                  female = female,
                  male = male,
                  band_date = band_date,
                  band_nums = band_nums)
      if(nrow(d) == 0){
        d <- tibble(hobo_logger = NA,
                    female = NA,
                    male = NA,
                    band_date = NA,
                    band_nums = NA)
      }
      return(d)
    }

    site <- rep(checks$site[1],nrow(attempt_map))
    box <- rep(checks$Nestbox[1],nrow(attempt_map))
    attempt <- attempt_map$attempt %>% as.numeric()
    species <- group_map(checks,get_species) %>% unlist()
    first_egg_date <- group_map(checks, ~ if_else(.x$Eggs[1] == 1,
                                                 as_date(.x$Date[1]),
                                                 as_date(.x$Date[1]) - .x$Eggs[1])) %>%
      unlist() %>%
      as_date()
    first_egg_certainty <- group_map(checks, ~ if_else(.x$Eggs[1] == 1,
                                                      'Confirmed',
                                                      'Estimated')) %>%
      unlist()
    eggs_hatched <- group_map(checks, ~if_else(nrow(.x) == 0,
                                               NA,
                                               max(.x$Nestlings,na.rm = TRUE))) %>%
      unlist()
    clutch_size <- group_map(checks, ~ if_else(any(!is.na(.x$Eggs)),
                                               max(.x$Eggs[which.max(.x$Eggs)] +
                                                 .x$Nestlings[which.max(.x$Eggs)],
                                                 max(.x$Nestlings,na.rm = TRUE),
                                                 na.rm = TRUE),
                                               NA_integer_)) %>% unlist()
    incubation_begins <- first_egg_date + clutch_size
    eggs_failed <- group_map(checks,eggs_unhatched) %>% unlist()
    eggs_failed <- if_else(eggs_hatched == clutch_size,0,eggs_failed)
    eggs_abandoned <- group_map(checks, eggs_aband) %>% unlist()
    eggs_eaten <- rep(0,nrow(attempt_map))
    eggs_destroyed <- rep(0,nrow(attempt_map)) # placeholder; we may have HOSP this season so may need to use this column more correctly
    egg_unk <- clutch_size - eggs_hatched - eggs_failed - eggs_abandoned - eggs_eaten - eggs_destroyed
    brood_size <- if_else(sapply(eggs_hatched,is.numeric),eggs_hatched,NA)
    hatch_date <- group_map(checks,get_hatch_date) %>%
      unlist() %>% as_date()
    hd_certainty <- group_map(checks,get_hatch_date_uncertainty) %>%
      unlist()
    early_banding <- group_map(checks,get_early_banding_date) %>% unlist() %>% as_date()
    late_banding <- as_date(early_banding + 2)
    chicks_fledged <- group_map(checks, count_fledged) %>% unlist()
    nest_perish <- if_else(chicks_fledged > 0,
                           brood_size - chicks_fledged,
                           0)
    nest_abandoned <- if_else(chicks_fledged == 0,
                           brood_size,
                           0)
    nest_pred <- if_else(sapply(chicks_fledged,is.numeric),
                         0,
                         NA)
    nest_interval <- interval(first_egg_date,late_banding,tz = 'US/Pacific')
    logger_and_id <- map2(.x = nest_interval,
                                .y = box,
                                get_logger_and_id) %>%
      list_rbind()
    hobo_logger <- pull(logger_and_id,hobo_logger)
    male <- pull(logger_and_id,male)
    female <- pull(logger_and_id,female)
    band_date <- pull(logger_and_id,band_date)
    band_nums <- pull(logger_and_id,band_nums)

  }
  # return tibble of columns

  return(tibble(Comments = hobo_logger,
                `OK to get fecal sample?` = NA,
                `Date of fecal sample` = NA,
                `Comments/info from fecal sample visit` = NA,
                Site = site,
                Nestbox = box,
                `Attempt No.` = attempt,
                Species = species,
                `Male ID` = male,
                `Female ID` = female,
                `First Egg Date` = first_egg_date,
                `First Egg Date Certainty` = first_egg_certainty,
                `Clutch Size` = clutch_size,
                `Incubation Begins Date` = incubation_begins,
                `Eggs Hatched` = eggs_hatched,
                `Eggs Failed-to-Hatch` = eggs_failed,
                `Eggs Abandoned` = eggs_abandoned,
                `Eggs Depredated` = eggs_eaten,
                `Eggs Destroyed` = eggs_destroyed,
                `Egg Unknown Fate` = egg_unk,
                `Brood Size` = brood_size,
                `Hatch Date` = hatch_date,
                `HD Observed or Estimated?` = hd_certainty,
                `Earliest Target Banding Date` = early_banding,
                `Latest Target Banding Date` = late_banding,
                `Date Banded` = band_date,
                `Nestlings Fledging` = chicks_fledged,
                `Nestlings Perishing` = nest_perish,
                `Nestlings Abandoned` = nest_abandoned,
                `Nestlings Depredated` = nest_pred,
                `Nestlings Destroyed` = 0,
                `BAND numbers of Nestlings` = band_nums,
                Notes = NA,
  ))

}

eggs_nestlings_to_numeric <- function(string) {
  strip <- str_match(string,"([0-9]+)|(Skipped)|(skipped)") %>% replace_na("0")
  strip[which(strip %in% c("Skipped","skipped"))] <- NA
  return(as.numeric(strip[,1]))
}

d <- nestbox_drive_get("nest_checks") %>%
  mutate(Eggs = eggs_nestlings_to_numeric(Eggs),Nestlings = eggs_nestlings_to_numeric(Nestlings),
         year = year(Date)) %>%
  group_by(Nestbox,year) %>%
  mutate(Species = ifelse(Species == "ATFL!","ATFL",Species),
         Species = ifelse(Species == "HOFI (?)","HOFI",Species),
         `Wing chords` = as.character(`Wing chords`),
         `Wing chords` = ifelse(`Wing chords` == "NULL",NA,`Wing chords`)) %>%
  mutate(site = str_extract(Nestbox,"([:alpha:]|[:punct:])+")) %>%
  filter(site %in% c('MBNG','MBNR','PICG','PICO','PICR',"RRRG","RRRR"))


comp <- group_split(d) %>% map(collapse) %>% map_dfr(bind_rows) %>%
  filter(!is.na(Site)) %>%
  mutate(Site = fct_collapse(Site,
                             `Russell Ranch` = c("RRRG","RRRR"),
                             `Mace Blvd` = c('MBNG','MBNR'),
                             `Picnic Grounds` = c('PICG','PICO','PICR'))) %>%
  arrange(year(`First Egg Date`),Site,Nestbox,`Attempt No.`)

write_csv(comp,"data/comprehensive_logs_Lauck.csv")

## Prepare nest checks for delivery to Melanie

require(hms)

n <- nestbox_drive_get("nest_checks")

nest <- n %>%
  mutate(site = str_extract(Nestbox,"([:alpha:]|[:punct:])+"),
         Site = fct_collapse(site,
                              `Russell Ranch` = c("RRRG","RRRR"),
                              `Mace Blvd` = c('MBNG','MBNR'),
                              `Picnic Grounds` = c('PICG','PICO','PICR')),
         `Start Time` = strftime(`Start Time`, format="%H:%M"),
         `End Time` = strftime(`End Time`, format="%H:%M"),
         Comments = NA,
         Nest = paste0(Nest,"\t")) %>%
  select(Comments,Date,`Start Time`,`End Time`,`Observer(s)`,Site,Nestbox,Species,`Adult(s)`,
         Nest,Eggs,Nestlings,Notes) %>%
  arrange(Date,Nestbox)

write_csv(nest,"data/nest_checks_Lauck.csv",na = '')

## Prepare morph data for delivery to Melanie

g <- read_csv("data/banding-and-morphometrics_proofed.csv") %>%
  mutate(site = str_extract(Nestbox,"([:alpha:]|[:punct:])+"),
         `Band Number` = band_number) %>%
  filter(!((site %in% c("BRO","MBNC","MCE","MCO","D","PCC","PCT","RRC","RRE","SFP")) &
             !is.na(band_number))) %>%
  select(`Band Number`,Disposition:Comments) %>%
  arrange(year(`Banding Date` %>% parse_date_time(orders = "mdy")),Site, Nestbox)

write_csv(g,"data/banding_data_Lauck.csv")
