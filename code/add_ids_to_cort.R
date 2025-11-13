#### add ids to cort

library(tidyverse)
require(ggeffects)
require(lme4)
require(viridis)
source("code/helper_functions.R")

c <- get_cort() %>% pivot_wider(names_from = sample, values_from = mean) %>% arrange(date)

attempt <- function(d,...) {
  if(nrow(d) > 1){
    attempt <- rep(NA,nrow(d))
    for(i in 2:nrow(d)){
      attempt[1] <- 1
      attempt[i] <- if_else(d$date[i]-d$date[1] < 21, 1,2)
    }
    d2 <- d %>% mutate(attempt = attempt) %>%
      filter(attempt == 2)
    if(nrow(d2) > 1) {
      attempt2 <- pull(d2,attempt)
      for(i in 2:nrow(d2)){
        attempt2[1] <- 2
        attempt2[i] <- if_else(d2$date[i]-d2$date[1] < 21, 2,3)
      }
      d <- mutate(d,attempt = attempt[which(attempt == 1)] %>% c(attempt2))
    } else {
      d <- mutate(d,attempt = attempt)
    }
  } else{
    d <- mutate(d,attempt = rep(1,nrow(d)))
  }
  return(d)
}

mom_blood_num <- function(d,...) {
  if(nrow(d) > 0){
    if("L" %in% d$Age){
      dg <- d %>% filter(Age == "L") %>%
        group_by(date) %>%
        summarize(brood_size = n())
      d <- mutate(d,brood_size = max(dg$brood_size))

      if(nrow(d) > 1 && "F" %in% filter(d,life_stage == "adult") %>% pull(Sex)){
        d <- mutate(d,blood_num_mother = (filter(d,life_stage == "adult",Sex == "F") %>% pull(blood_num))[1])
      } else {
      d <- mutate(d,blood_num_mother = NA)
      }

    } else {
      d <- mutate(d,blood_num_mother = NA,brood_size = NA)
    }
  } else {
    d <- mutate(d,blood_num_mother = NA, brood_size = NA)
  }
  return(d)
}

m <- read_rds("data/g.rds") %>%
  mutate(blood_num = ifelse(str_detect(Comments,"B([:digit:]){3}"), # is there a blood number? if so,
                            str_extract(Comments,"B( )?([:digit:]){3}"), # extract it
                            `Blood number`), # if not, use the listed blood number, which may be NA
         blood_num = if_else(blood_num == "-", NA, blood_num),
         color = Comments %>% str_to_lower() %>% str_extract("(blue|green|red|yellow|orange|pink|purple)2?"),
         year = year(`Banding Date`), # Extract year
         date = as_date(`Banding Date`), # Extract date
         weight = as.character(`Weight (g)`) %>% as.numeric(), # make all growth columns numeric
         wing = as.character(`Right Wing Chord (mm)`) %>% as.numeric(),
         tail = as.character(`Tail Length (mm)`) %>% as.numeric(),
         skull = as.character(`Skull Length (mm)`) %>% as.numeric(),
         bill = as.character(`Nares-Tip (mm)`) %>% as.numeric(),
         tarsus = as.character(`Tarsus (mm)`) %>% as.numeric(),
         site = factor(str_extract(Nestbox,"([:alpha:]|[:punct:])+")), # extract site from nestbox
         habitat = fct_collapse(site, # collapse sites into habitat classes
                                Orchard = c("BRO","MCE","PICO","MCO","FBFO"),
                                Grassland = c("MBNG","RRRG","PICG","FBFG"),
                                `Row crop` = c("MBNR","PICR","RRRR","FBFR"),
                                Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","FBF","RRRC")),
         habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")), # order factor correctly
         juliandate = yday(date)) %>%
  arrange(year,Nestbox,date,color) %>% # Order by year, date,nestbox, color
  group_by(Nestbox, year) %>% # group by nestbox, year.
  group_map(attempt,.keep = TRUE) %>% # label first and second (and potentially third) attempts
  list_rbind() %>%
  group_by(Nestbox,year,attempt) %>%
  group_map(mom_blood_num, .keep = TRUE) %>%
  list_rbind() %>%
  left_join(select(c,c(blood_num,s1,s2)),
            by = join_by(blood_num_mother == blood_num),
            multiple = "first") %>%
  rename(s1_mother = s1,s2_mother = s2)

cort <- c %>% left_join(select(m,c(blood_num,blood_num_mother,brood_size,s1_mother,s2_mother,`Band Number`,Species, Nestbox, year,attempt, date, juliandate, habitat,site, Age,Sex,color,
                                   weight:tarsus)),by = c("blood_num")) %>%
  mutate(across(c(s1, # scale numeric responses and predictors
                  s2,
                  brood_size,
                  juliandate,
                  s1_mother,
                  s2_mother,
                  weight:tarsus),
                ~ scale(.x)[,1],
                .names = "{.col}_scaled"))

write_rds(cort,"data/cort_with_id.rds")

