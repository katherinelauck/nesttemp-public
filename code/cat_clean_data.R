#### Build growth data for SEM modelling, including cort

# Author: Katherine Lauck
# Last updated: 15 Feb 2024

require(tidyverse)
require(lubridate)
source("code/helper_functions.R")

### todo
# 1. Look for extreme high cort samples: take the reasonable replicate if there is one
# 2. Are any samples mislabeled? If any S1's are higher than S2's, they're likely swapped
# 3. Missing data WEBL s2?
# 4. Need to upload final bits of logger data

# pull data from google drive (not yet fully proofed)
g <- read_csv("data/banding-and-morphometrics_proofed.csv")
t <- read_rds("data/temp.rds")
cc <- read_rds("data/cc_full.rds")
h <- read_rds("data/humidity.rds")
raw_cort <- get_cort()

wmean <- t %>%
  mutate(site = factor(str_extract(box,"([:alpha:]|[:punct:])+")),
         habitat = fct_collapse(site,
                                Orchard = c("BRO","MCE","PICO","MCO","FBFO","PG-O-"),
                                Grassland = c("MBNG","RRRG","PICG","FBFG"),
                                `Row crop` = c("MBNR","PICR","RRRR","FBFR","PG-C-","RU-C-"),
                                Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","PCE")),
         habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop"))) %>%
  group_by(date,habitat,logger_position) %>%
  summarize(meanmax = mean(maxt,na.rm = TRUE),n = n(),.groups = "keep") %>%
  ungroup(habitat) %>%
  summarize(wmean = weighted.mean(meanmax,n)) %>%
  right_join(t,by = c("date","logger_position")) %>%
  mutate(resid = maxt - wmean, year = year(date),
        box = if_else(box == "PG-C-03", "PICR3", box),
        box = if_else(box == "RU-C-02", "RRRR2", box),
        box = if_else(box == "PG-O-01", "PICO1", box)) %>%
  filter(resid != 0) %>%
  left_join(cc, by = join_by(year,box == nestbox), multiple = "first")

# check boxes with no canopy cover measurement

# filter(wmean,is.na(canopy_cover)) %>%
#   ungroup(date) %>%
#   select(c(year,box)) %>%
#   distinct() %>%
#   write_csv("data/no_cc.csv")
#
# # manually entered in data for many boxes, now want to check for cc measurements in other years to fill in remainder
#
# no_cc <- read_csv("data/no_cc.csv")
#
# add_cc <- filter(wmean,box %in% (filter(no_cc,is.na(canopy_cover)) %>%
#                                    pull(box))) %>%
#   ungroup(date) %>%
#   select(c(year,box,canopy_cover)) %>%
#   distinct() %>% arrange(box) %>%
#   filter(!is.na(canopy_cover))
#
# # add these canopy cover replacements to the canopy cover database
#
# cc_full <- left_join(no_cc,select(add_cc,c(box,canopy_cover)),
#                         by = join_by(box == box)) %>%
#   arrange(box) %>%
#   mutate(canopy_cover = if_else(is.na(canopy_cover.x),
#                                 canopy_cover.y,
#                                 canopy_cover.x)) %>%
#   select(-c(canopy_cover.x,canopy_cover.y)) %>%
#   rename(nestbox = box) %>%
#   add_row(cc) %>%
#   write_rds("data/cc_full.rds")
#
# # replace canopy cover in wmean with the correct canopy cover
#
# wmean <- wmean %>% select(-canopy_cover) %>%
#   left_join(cc_full, by = join_by(year,box == nestbox), multiple = "first")
#
write_rds(wmean, "data/wmean.rds")
#
# # make updated canopy_cover the default canopy cover
#
# cc <- cc_full

c <- raw_cort %>%
  arrange(date) %>%
  filter(plate != "A02")

prob_c <- c %>% mutate(diff_s2_s1 = mean_s2-mean_s1)
hist(prob_c$diff_s2_s1)
filter(prob_c,diff_s2_s1 < 0) %>% write_csv("data/mislabelled_cort.csv")

c$r2_s1[which(c$blood_num == "B486")] <- NA
c$r2_s1[which(c$blood_num == "B582")] <- NA
c$r1_s2[which(c$blood_num == "B582")] <- NA
c$r1_s1[which(c$blood_num == "B584")] <- NA
c$r2_s2[which(c$blood_num == "B584")] <- NA
c$r2_s2[which(c$blood_num == "B584")] <- NA
c$r2_s2[which(c$blood_num == "B088")] <- NA
b088_r1_s1 <- c$r1_s1[which(c$blood_num == "B088")]
b088_r2_s1 <- c$r2_s1[which(c$blood_num == "B088")]
b088_r1_s2 <- c$r1_s2[which(c$blood_num == "B088")]
c$r1_s1[which(c$blood_num == "B088")] <- b088_r1_s2
c$r1_s2[which(c$blood_num == "B088")] <- b088_r1_s1
c$r2_s2[which(c$blood_num == "B088")] <- b088_r2_s1

B107_r1_s1 <- c$r1_s1[which(c$blood_num == "B107")]
B107_r2_s1 <- c$r2_s1[which(c$blood_num == "B107")]
B107_r1_s2 <- c$r1_s2[which(c$blood_num == "B107")]
B107_r2_s2 <- c$r2_s2[which(c$blood_num == "B107")]
c$r1_s1[which(c$blood_num == "B107")] <- B107_r1_s2
c$r2_s1[which(c$blood_num == "B107")] <- B107_r2_s2
c$r1_s2[which(c$blood_num == "B107")] <- B107_r1_s1
c$r2_s2[which(c$blood_num == "B107")] <- B107_r2_s1

B409_r1_s1 <- c$r1_s1[which(c$blood_num == "B409")]
B409_r2_s1 <- c$r2_s1[which(c$blood_num == "B409")]
B409_r1_s2 <- c$r1_s2[which(c$blood_num == "B409")]
B409_r2_s2 <- c$r2_s2[which(c$blood_num == "B409")]
c$r1_s1[which(c$blood_num == "B409")] <- B409_r1_s2
c$r2_s1[which(c$blood_num == "B409")] <- B409_r2_s2
c$r1_s2[which(c$blood_num == "B409")] <- B409_r1_s1
c$r2_s2[which(c$blood_num == "B409")] <- B409_r2_s1

B338_r1_s1 <- c$r1_s1[which(c$blood_num == "B338")]
B338_r2_s1 <- c$r2_s1[which(c$blood_num == "B338")]
B338_r1_s2 <- c$r1_s2[which(c$blood_num == "B338")]
B338_r2_s2 <- c$r2_s2[which(c$blood_num == "B338")]
c$r1_s1[which(c$blood_num == "B338")] <- B338_r1_s2
c$r2_s1[which(c$blood_num == "B338")] <- B338_r2_s2
c$r1_s2[which(c$blood_num == "B338")] <- B338_r1_s1
c$r2_s2[which(c$blood_num == "B338")] <- B338_r2_s1

B342_r1_s1 <- c$r1_s1[which(c$blood_num == "B342")]
B342_r2_s1 <- c$r2_s1[which(c$blood_num == "B342")]
B342_r1_s2 <- c$r1_s2[which(c$blood_num == "B342")]
B342_r2_s2 <- c$r2_s2[which(c$blood_num == "B342")]
c$r1_s1[which(c$blood_num == "B342")] <- B342_r1_s2
c$r2_s1[which(c$blood_num == "B342")] <- B342_r2_s2
c$r1_s2[which(c$blood_num == "B342")] <- B342_r1_s1
c$r2_s2[which(c$blood_num == "B342")] <- B342_r2_s1

B129_r1_s1 <- c$r1_s1[which(c$blood_num == "B129")]
B129_r2_s1 <- c$r2_s1[which(c$blood_num == "B129")]
B129_r1_s2 <- c$r1_s2[which(c$blood_num == "B129")]
B129_r2_s2 <- c$r2_s2[which(c$blood_num == "B129")]
c$r1_s1[which(c$blood_num == "B129")] <- B129_r1_s2
c$r2_s1[which(c$blood_num == "B129")] <- B129_r2_s2
c$r1_s2[which(c$blood_num == "B129")] <- B129_r1_s1
c$r2_s2[which(c$blood_num == "B129")] <- B129_r2_s1

B233_r1_s1 <- c$r1_s1[which(c$blood_num == "B233")]
B233_r2_s1 <- c$r2_s1[which(c$blood_num == "B233")]
B233_r1_s2 <- c$r1_s2[which(c$blood_num == "B233")]
B233_r2_s2 <- c$r2_s2[which(c$blood_num == "B233")]
c$r1_s1[which(c$blood_num == "B233")] <- B233_r1_s2
c$r2_s1[which(c$blood_num == "B233")] <- B233_r2_s2
c$r1_s2[which(c$blood_num == "B233")] <- B233_r1_s1
c$r2_s2[which(c$blood_num == "B233")] <- B233_r2_s1

B302_r1_s1 <- c$r1_s1[which(c$blood_num == "B302")]
B302_r2_s1 <- c$r2_s1[which(c$blood_num == "B302")]
B302_r1_s2 <- c$r1_s2[which(c$blood_num == "B302")]
B302_r2_s2 <- c$r2_s2[which(c$blood_num == "B302")]
c$r1_s1[which(c$blood_num == "B302")] <- B302_r1_s2
c$r2_s1[which(c$blood_num == "B302")] <- B302_r2_s2
c$r1_s2[which(c$blood_num == "B302")] <- B302_r1_s1
c$r2_s2[which(c$blood_num == "B302")] <- B302_r2_s1



problem_cort <- filter(c,d_s1 > 5 | d_s2 > 5)
problem_cort <- filter(c,(mean_s1 > 15) | (mean_s2 > 80))
problem_cort <- filter(c,d_s1 > 5 | d_s2 > 5 | (mean_s1 > 15) | (mean_s2 > 80))

write_csv(problem_cort,"data/problem_cort.csv")

c$r2_s1[which(c$blood_num == "B571")] <- NA
c$r1_s1[which(c$blood_num == "B546")] <- NA
c$r1_s2[which(c$blood_num == "B579")] <- NA
c$r1_s2[which(c$blood_num == "B580")] <- NA
c$r2_s1[which(c$blood_num == "B581")] <- NA
c$r1_s2[which(c$blood_num == "B583")] <- NA
c$r1_s2[which(c$blood_num == "B583")] <- NA
c[which(c$blood_num == "B036"),c("r1_s2","r2_s2")] <- c[which(c$blood_num == "B036"),c("r1_s1","r2_s1")]
c[which(c$blood_num == "B036"),c("r1_s1","r2_s1")] <- NA
c[which(c$blood_num == "B087"),c("r1_s2","r2_s2")] <- c[which(c$blood_num == "B087"),c("r1_s1","r2_s1")]
c[which(c$blood_num == "B087"),c("r1_s1","r2_s1")] <- NA
c$r2_s1[which(c$blood_num == "B088")] <- NA
c$r2_s1[which(c$blood_num == "B517")] <- NA
c$r1_s1[which(c$blood_num == "B107")] <- NA
c$r2_s1[which(c$blood_num == "B107")] <- NA
c[which(c$blood_num == "B412"),c("r1_s1","r2_s1")] <- c[which(c$blood_num == "B412"),c("r1_s2","r2_s2")]
c[which(c$blood_num == "B412"),c("r1_s2","r2_s2")] <- NA
c$r2_s2[which(c$blood_num == "B463")] <- NA
c$r2_s2[which(c$blood_num == "B329")] <- NA
c$r1_s1[which(c$blood_num == "B384")] <- NA
c$r1_s1[which(c$blood_num == "B130")] <- NA
c$r2_s1[which(c$blood_num == "B130")] <- NA
c[which(c$blood_num == "B202"),c("r1_s2","r2_s2")] <- c[which(c$blood_num == "B202"),c("r1_s1","r2_s1")]
c[which(c$blood_num == "B202"),c("r1_s1","r2_s1")] <- NA
c$r1_s1[which(c$blood_num == "B301")] <- NA
c[which(c$blood_num == "B308"),c("r1_s2","r2_s2")] <- c[which(c$blood_num == "B308"),c("r1_s1","r2_s1")]
c[which(c$blood_num == "B308"),c("r1_s1","r2_s1")] <- NA
c$r1_s2[which(c$blood_num == "B273")] <- NA
c[which(c$blood_num == "B516"),c("r1_s2","r2_s2")] <- c[which(c$blood_num == "B516"),c("r1_s1","r2_s1")]
c[which(c$blood_num == "B516"),c("r1_s1","r2_s1")] <- NA
c[which(c$blood_num == "B581"),c("r1_s1","r2_s1")] <- NA
c[which(c$blood_num == "B228"),c("r1_s1","r2_s1")] <- NA

c <- c %>% rowwise() %>%
  mutate(mean_s1 = mean(c(r1_s1,r2_s1),na.rm = TRUE),
         mean_s2 = mean(c(r1_s2,r2_s2),na.rm = TRUE),
         d_s1 = abs(r1_s1-r2_s1),
         d_s2 = abs(r1_s2-r2_s2)) %>%
  ungroup()


wing_tbl <- tibble( # recreate wing to age chart from Museum of Fish and Wildlife Biology
  age = c(rep(seq(0,18,1),3),0,0),
  species = c(rep(c("WEBL","TRES","ATFL"),each = 19),"skipped","Skipped"),
  wing = c(6,7,8,9,11,13,16,19,23,27,31,35,39,43,47,51,55,59,63,6,6,7,8,10,13,16,20,25,30,35,41,47,52,57,62,66,72,78,8,8,9,12,14,17,21,25,32,38,44,50,55,60,65,70,75,77,79,0,0)
) %>% group_by(species)

# function to calculate growth by matching color and nestbox

growth <- function(d,...) {
  out <- d %>% select(band_number,year,Nestbox,color,date,Species,Sex,blood_num,Age,life_stage,attempt,attempt_id, # select relevant columns
                      weight,
                      wing,
                      tail,
                      skull,
                      bill,
                      tarsus) %>%
    arrange(date) # order by date
  if(nrow(out)==2){ # if we have more than one measurement of chick per nest, meaning we are able to calculate growth, then:
    out <- add_column(out, # calculate growth in grams or millimeters per day (so a standardized rate)
                      gweight = c(NA,(out$weight[2]-out$weight[1])/as.numeric(out$date[2]-out$date[1])),
                      gwing = c(NA,(out$wing[2]-out$wing[1])/as.numeric(out$date[2]-out$date[1])),
                      gtail = c(NA,(out$tail[2]-out$tail[1])/as.numeric(out$date[2]-out$date[1])),
                      gskull = c(NA,(out$skull[2]-out$skull[1])/as.numeric(out$date[2]-out$date[1])),
                      gbill = c(NA,(out$bill[2]-out$bill[1])/as.numeric(out$date[2]-out$date[1])),
                      gtarsus = c(NA,(out$tarsus[2]-out$tarsus[1])/as.numeric(out$date[2]-out$date[1])))
    out <- add_column(out,
                      pweight = c(NA,out$gweight[2]/out$weight[1]/as.numeric(out$date[2]-out$date[1])),
                      pwing = c(NA,out$gwing[2]/out$wing[1]/as.numeric(out$date[2]-out$date[1])),
                      ptail = c(NA,out$gtail[2]/out$tail[1]/as.numeric(out$date[2]-out$date[1])),
                      pskull = c(NA,out$gskull[2]/out$skull[1]/as.numeric(out$date[2]-out$date[1])),
                      pbill = c(NA,out$gbill[2]/out$bill[1]/as.numeric(out$date[2]-out$date[1])),
                      ptarsus = c(NA,out$gtarsus[2]/out$tarsus[1]/as.numeric(out$date[2]-out$date[1])),
                      occasion = c(1,2),
                      ind_id = if_else(any(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}")),
                                       out$band_number[which(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}"))][1],
                                       paste0(out$attempt_id,"_",out$color)[1],
                                       missing = paste0(out$attempt_id,"_",out$color)[1]))
  } else if(nrow(out)==3){ # if we have three measurements, calculate two growth intervals per day
    out <- add_column(out,
                      gweight = c(NA,(out$weight[2]-out$weight[1])/as.numeric(out$date[2]-out$date[1]),
                                  (out$weight[3]-out$weight[2])/as.numeric(out$date[3]-out$date[2])),
                      gwing = c(NA,(out$wing[2]-out$wing[1])/as.numeric(out$date[2]-out$date[1]),
                                (out$wing[3]-out$wing[2])/as.numeric(out$date[3]-out$date[2])),
                      gtail = c(NA,(out$tail[2]-out$tail[1])/as.numeric(out$date[2]-out$date[1]),
                                (out$tail[3]-out$tail[2])/as.numeric(out$date[3]-out$date[2])),
                      gskull = c(NA,(out$skull[2]-out$skull[1])/as.numeric(out$date[2]-out$date[1]),
                                 (out$skull[3]-out$skull[2])/as.numeric(out$date[3]-out$date[2])),
                      gbill = c(NA,(out$bill[2]-out$bill[1])/as.numeric(out$date[2]-out$date[1]),
                                (out$bill[3]-out$bill[2])/as.numeric(out$date[3]-out$date[2])),
                      gtarsus = c(NA,(out$tarsus[2]-out$tarsus[1])/as.numeric(out$date[2]-out$date[1]),
                                  (out$tarsus[3]-out$tarsus[2])/as.numeric(out$date[3]-out$date[2])))
    out <- add_column(out,
                      pweight = c(NA,(out$gweight[2]/out$weight[1])/as.numeric(out$date[2]-out$date[1]),
                                  (out$gweight[3]/out$weight[2])/as.numeric(out$date[3]-out$date[2])),
                      pwing = c(NA,(out$gwing[2]/out$wing[1])/as.numeric(out$date[2]-out$date[1]),
                                (out$gwing[3]/out$wing[2])/as.numeric(out$date[3]-out$date[2])),
                      ptail = c(NA,(out$gtail[2]/out$tail[1])/as.numeric(out$date[2]-out$date[1]),
                                (out$gtail[3]/out$tail[2])/as.numeric(out$date[3]-out$date[2])),
                      pskull = c(NA,(out$gskull[2]/out$skull[1])/as.numeric(out$date[2]-out$date[1]),
                                 (out$gskull[3]/out$skull[2])/as.numeric(out$date[3]-out$date[2])),
                      pbill = c(NA,(out$gbill[2]/out$bill[1])/as.numeric(out$date[2]-out$date[1]),
                                (out$gbill[3]/out$bill[2])/as.numeric(out$date[3]-out$date[2])),
                      ptarsus = c(NA,(out$gtarsus[2]/out$tarsus[1])/as.numeric(out$date[2]-out$date[1]),
                                  (out$gtarsus[3]/out$tarsus[2])/as.numeric(out$date[3]-out$date[2])),
                      occasion = c(1,2,3),
                      ind_id = if_else(any(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}")),
                                       out$band_number[which(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}"))][1],
                                       paste0(out$attempt_id,"_",out$color)[1],
                                       missing = paste0(out$attempt_id,"_",out$color)[1]))
  } else if(nrow(out) > 3) {
    out <- add_column(out,
                      gweight = NA,
                      gwing = NA,
                      gtail = NA,
                      gskull = NA,
                      gbill = NA,
                      gtarsus = NA,
                      pweight = NA,
                      pwing = NA,
                      ptail = NA,
                      pskull = NA,
                      pbill = NA,
                      ptarsus = NA,
                      occasion = 99,
                      ind_id = if_else(any(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}")),
                                       out$band_number[which(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}"))][1],
                                       paste0(out$attempt_id,"_",out$color)[1],
                                       missing = paste0(out$attempt_id,"_",out$color)[1]))
  } else { # otherwise, growth can't be calculated.
    out <- add_column(out,
                      gweight = NA,
                      gwing = NA,
                      gtail = NA,
                      gskull = NA,
                      gbill = NA,
                      gtarsus = NA,
                      pweight = NA,
                      pwing = NA,
                      ptail = NA,
                      pskull = NA,
                      pbill = NA,
                      ptarsus = NA,
                      occasion = 1,
                      ind_id = if_else(any(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}")),
                                       out$band_number[which(str_detect(out$band_number,"[:digit:]{4}-[:digit:]{5}"))][1],
                                       paste0(out$attempt_id,"_",out$color)[1],
                                       missing = paste0(out$attempt_id,"_",out$color)[1]))
  }
}

attempt <- function(d,...) {
  if(nrow(d) > 1){
    attempt <- rep(NA,nrow(d))
    for(i in 2:nrow(d)){
      attempt[1] <- 1
      attempt[i] <- if_else(d$date[i]-d$date[1] < 28, 1,2)
    }
    d2 <- d %>% mutate(attempt = attempt) %>%
      filter(attempt == 2)
    if(nrow(d2) > 1) {
      attempt2 <- pull(d2,attempt)
      for(i in 2:nrow(d2)){
        attempt2[1] <- 2
        attempt2[i] <- if_else(d2$date[i]-d2$date[1] < 28, 2,3)
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

brood_size <- function(d,...) {
  if(nrow(d) > 0){
    if("L" %in% d$Age){
      dg <- d %>% filter(Age == "L") %>%
        group_by(date) %>%
        summarize(brood_size = n(),
                  age = wing_tbl$age[which.min(abs(wing_tbl$wing-mean(wing,na.rm = TRUE)))])
      d <- mutate(d,brood_size = max(dg$brood_size,na.rm = TRUE)) %>%
        left_join(select(dg,c("date","age")),by = c("date"))

      if(nrow(d) > 1 && "F" %in% (filter(d,life_stage == "adult") %>% pull("Sex"))){
        d <- mutate(d,blood_num_mother = (filter(d,life_stage == "adult",Sex == "F") %>% pull("blood_num"))[1])
      } else {
        d <- mutate(d,blood_num_mother = NA)
      }

    } else {
      d <- mutate(d,brood_size = NA,age = NA,blood_num_mother = NA)
    }
  } else {
    d <- mutate(d, brood_size = NA,age = NA,blood_num_mother = NA)
  }
  return(d)
}

gr <- g %>%
  mutate(
    # blood_num = ifelse(str_detect(Comments,"B([:digit:]){3}"), # is there a blood number? if so,
    #                         str_extract(Comments,"B( )?([:digit:]){3}"), # extract it
    #                         `Blood number`), # if not, use the listed blood number, which may be NA
         blood_num = if_else(`Blood number` == "-", NA, `Blood number`), # if - in blood num, insert NA
    #      color = Comments %>% str_to_lower() %>%
    #        str_extract("(blue|green|red|yellow|orange|pink|purple)2?"), # Extract color from comments
         date = mdy(`Banding Date`), # Extract date
         year = year(date), # Extract year
         weight = as.character(`Weight (g)`) %>% as.numeric(), # make all growth columns numeric
         wing = as.character(`Right Wing Chord (mm)`) %>% as.numeric(),
         tail = as.character(`Tail Length (mm)`) %>% as.numeric(),
         skull = as.character(`Skull Length (mm)`) %>% as.numeric(),
         bill = as.character(`Nares-Tip (mm)`) %>% as.numeric(),
         tarsus = as.character(`Tarsus (mm)`) %>% as.numeric()) %>% # make all growth columns numeric
  arrange(year,Nestbox,date,color) %>% # Order by year, date,nestbox, color
  group_by(Nestbox, year) %>% # group by nestbox, year.
  group_map(attempt,.keep = TRUE) %>% # label first and second (and potentially third) attempts
  list_rbind() %>%
  mutate(attempt_id = paste0(year,"_",Nestbox,"_",attempt)) %>%
  group_by(attempt_id,color) %>% # add final grouping variables
  group_map(growth,.keep = TRUE) %>% # map growth over each group. See todo.
  list_rbind() %>% # bind into tibble
  group_by(attempt_id) %>%
  group_map(brood_size,.keep = TRUE) %>%
  list_rbind() # %>%
  # filter(!is.na(gweight),!is.infinite(gweight)) # drop NAs and infinites

dg <- gr %>%
  mutate(site = factor(str_extract(Nestbox,"([:alpha:]|[:punct:])+")), # extract site from nestbox
         habitat = fct_collapse(site, # collapse sites into habitat classes
                                Orchard = c("BRO","MCE","PICO","MCO","FBFO"),
                                Grassland = c("MBNG","RRRG","PICG","FBFG"),
                                `Row crop` = c("MBNR","PICR","RRRR","FBFR"),
                                Forest = c("MBNC","PCC","PCT","RRC","RRE","SFP","RRRC")),
         habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")), # order factor correctly
         interval = interval(date-ddays(7),date,tz = "America/Los_Angeles"), # specify interval to extract temp measurements
         juliandate = yday(date)) %>% # extract julian date
  left_join(cc,by = join_by(year == year,Nestbox == nestbox),multiple = "first") # join canopy cover measurements

filtertemp <- map2(pull(dg,interval),pull(dg,Nestbox),function(x,y){filter(t,date %within% x,box == y)}) # filter temp data to the applicable dates and boxes
filterh <- map2(pull(dg,interval),pull(dg,site),function(x,y){filter(h,date %within% x,as.character(site) == as.character(y))})

dg <- dg %>%
  mutate(meanmaxtempI = filtertemp %>% # calculate mean max internal temp
           map_dbl(function(x){filter(x,logger_position == "I") %>% pull(maxt) %>% mean()}),
         maxmaxtempI = filtertemp %>% # max max internal temp
           map_dbl(function(x){filter(x,logger_position == "I") %>% pull(maxt) %>% max(na.rm = TRUE)}),
         cumulativeover40I = filtertemp %>% # time over 40 degrees C internal
           map_dbl(function(x){filter(x,logger_position == "I") %>% pull(time_above_40) %>% sum()}),
         meanmaxtempO = filtertemp %>% # mean max external temp
           map_dbl(function(x){filter(x,logger_position == "O") %>% pull(maxt) %>% mean()}),
         maxmaxtempO = filtertemp %>% # max max external temp
           map_dbl(function(x){filter(x,logger_position == "O") %>% pull(maxt) %>% max(na.rm = TRUE)}),
         cumulativeover40O = filtertemp %>% # time over 40 degrees C external
           map_dbl(function(x){filter(x,logger_position == "O") %>% pull(time_above_40) %>% sum}),
         meanh = filterh %>% # mean humidity
           map_dbl(function(x){pull(x,meanh) %>% mean()}),
         meanmaxh = filterh %>% # mean max humidity
           map_dbl(function(x){pull(x,maxh) %>% mean()}),
         meanminh = filterh %>% # mean min humidity
           map_dbl(function(x){pull(x,minh) %>% mean()}),
         meanmintempI = filtertemp %>% # mean min internal temp
           map_dbl(function(x){filter(x,logger_position == "I") %>% pull(mint) %>% mean()}),
         minmintempI = filtertemp %>% # min min internal temp
           map_dbl(function(x){filter(x,logger_position == "I") %>% pull(mint) %>% min(na.rm = TRUE)}),
         meanmintempO = filtertemp %>% # mean min external temp
           map_dbl(function(x){filter(x,logger_position == "O") %>% pull(mint) %>% mean()}),
         minmintempO = filtertemp %>% # min min external temp
           map_dbl(function(x){filter(x,logger_position == "O") %>% pull(mint) %>% min(na.rm = TRUE)})) %>%
  # filter(!is.infinite(gweight)) %>% # filter out infinite gweights - need to figure out why this is happening
  filter(!(site %in% c("FBFG","FBFR","FBFO")), # remove Full Belly Farm nests
         # life_stage == "nestling"
         # ,
         # !(year == 2023) # remove 2023 nests, but this line can be dropped soon
         ) %>%
  mutate(site = str_replace_all(site,"MCE","MCO"), # Collapse site vectors when there are multiple site abbreviations per actual site
         site = str_replace_all(site,"PCT","PCC"),
         site = str_replace_all(site,"RRE","RRC"),
         site = str_replace_all(site,"SFP","MBNC"),
         site = as_factor(site),
         site = fct_collapse(site,
                             PG = c("PCC","PICG","PICO","PICR"),
                             RR = c("RRC","RRRG","RRRR"),
                             MB = c("MBNC","MBNG","MBNR"),
                             WI = c("BRO","MCO")), # factorize site
         year = as_factor(year), # factorize year
         Species = as_factor(Species),
         condition = weight/tarsus) %>% # factorize species
  left_join(select(c,c(blood_num,mean_s1,mean_s2)),
            by = join_by(blood_num_mother == blood_num),
            multiple = "last") %>%
  rename(cort_s1_mother = mean_s1,cort_s2_mother = mean_s2) %>%
  left_join(select(c,blood_num:d_s2),
            by = "blood_num",
            multiple = "last") %>%
  rename(cort_s1 = mean_s1, cort_s2 = mean_s2)

cort <- c %>% left_join(select(dg,c(blood_num,blood_num_mother,brood_size,cort_s1_mother,cort_s2_mother,Species, Nestbox, year, habitat,canopy_cover,site,attempt_id, date, juliandate,ind_id,band_number,occasion,Age,age,Sex,meanh,meanmaxtempI,
                                   condition,weight:tarsus,gweight:gtarsus)),
                        by = c("blood_num"),multiple = "last", na_matches = "never") %>%
  filter(!is.na(Nestbox))

dg <- dg %>%
  left_join(select(cort,c(blood_num,condition)) %>% rename(condition_mother = condition),
            by = join_by(blood_num_mother == blood_num),
            multiple = "last")

write_rds(cort,"data/cort_with_id.rds")

write_rds(filter(dg,life_stage == "nestling"),"data/growth.rds")

trigly <- read_csv("data/trigly_assay/plate_b1_triglycerides.csv") %>% left_join(cort,by = join_by(sample_num==blood_num),multiple = "last")

write_rds(trigly, "data/trigly_with_id.rds")

##### survival data

## To-do:
## - filter to attempts that I monitored
## - make table with data sheet

# s <- get_nest_checks() %>% write_rds("data/survival_katieboxes.rds")
#
# g <- read_rds("data/growth.rds")
#
# read_survival_data_csvs <- function(file) {
#   col_names <- read_csv("data/pcnh_comp_logs/mbn_2022.csv")[,c(5:22,27)] %>% names()
#   read_csv(file) %>% select(col_names) %>%
#     mutate(across(c(`First Egg Date`,`Incubation Begins Date`,`Hatch Date`),mdy),
#            across(c(`Clutch Size`,`Eggs Hatched`:`Brood Size`,`Nestlings Fledging`),
#                   ~ as.character(.x) %>% as.numeric()))
# }
#
# files <- c("data/pcnh_comp_logs/mbn_2022.csv","data/pcnh_comp_logs/sfp_2023.csv","data/pcnh_comp_logs/pic_2022.csv","data/pcnh_comp_logs/pic_2023.csv","data/pcnh_comp_logs/rrr_2022.csv","data/pcnh_comp_logs/rrr_2023.csv","data/pcnh_comp_logs/bro_2022.csv","data/pcnh_comp_logs/bro_2023.csv","data/pcnh_comp_logs/mcn_2022.csv","data/pcnh_comp_logs/all_2021.csv")
#
# s <- files %>% map(read_survival_data_csvs) %>% bind_rows() %>% bind_rows(read_rds("data/survival_katieboxes.rds")) %>%
#   mutate(year = year(`Hatch Date`),
#          boxyear = str_c(Nestbox,year %>% as.character())) %>%
#   filter(boxyear %in% str_c(g$Nestbox,g$year)) %>%
#   mutate(site = str_extract(s$Nestbox,"([:alpha:]|[:punct:])+") %>% factor() %>%
#            fct_collapse(PCC = c("PCC","PCT"),
#                         RRC = c("RRE","RRC"),
#                         MBNC = c("MBNC","SFP"),
#                         WI = c("BRO","MCE","MCO")), # extract site from nestbox
#          habitat = fct_collapse(site,
#                                 Orchard = c("WI","PICO"),
#                                 Grassland = c("MBNG","RRRG","PICG"),
#                                 `Row crop` = c("MBNR","PICR","RRRR"),
#                                 Forest = c("MBNC","PCC","RRC")),
#          habitat = factor(habitat,levels = c("Forest","Orchard","Grassland","Row crop")),
#          )
# levels(s$site)
#
# write_rds(s,"data/survival.rds")
