### number of chicks fledged per attempt by site, habitat, year

require(tidyverse)
require(lubridate)
require(gt)

s <- read_rds("data/survival.rds") # the data I used, attached to this email. obv pls delete after porting this code to your own data

s %>% # data in the table will just be Nestlings Fledging
  mutate(`Nestlings Fledging` = replace_na(`Nestlings Fledging`,0)) %>% # get rid of NAs. They're not relevant here, but they would be relevant if I were averaging per nest.
  select(`Nestlings Fledging`,site,habitat,year) %>% group_by(habitat,site,year) %>% # Select data + grouping columns, then group
  summarize(fledgelings = sum(`Nestlings Fledging`)) %>% # Summarize over habitat, site, year
  pivot_wider(names_from = c(habitat,site), values_from = fledgelings) %>% # pivot into shape of table - years on rows, habitat and site on columns
  mutate(across(-year, ~ replace_na(.x,0))) %>% # There's one site with no attempts in 2023 - replace that NA with a 0 so it will be properly colored. I found out in the process of coloring the table that you can change how values display after using gt() but the underlying value doesn't change - and the color function uses the underlying values.
  ungroup() %>% # remove grouping
  group_by(year) %>% # add the grouping I want for summary
  mutate(Total = sum(Forest_MBNC:`Row crop_RRRR`)) %>% # Create summary column adding together all the nestlings from each year
  ungroup() %>% # drop grouping again because gt formats with respect to it and I don't want that. I think gt can actually do this computation but I prefer to do it manually because I can see and adjust each step
  gt() %>% tab_spanner_delim(delim = "_") %>% # split column names to group sites under habitat. First part of name will end up as a column spanner, second will be the column name
  data_color(columns = -c(year,Total),fn = scales::col_numeric( # manually color the entire dataset. In this line, indicate which columns to avoid and set the binning algorithm - in this case it's just linear.
    palette = viridis::viridis(seq(select(.$`_data`,-c(year,Total)) %>% min(),
                                   select(.$`_data`,-c(year,Total)) %>% max()) %>% length()), # set the palette, and tell it how many colors to generate
    domain = c(select(.$`_data`,-c(year,Total)) %>% min(),select(.$`_data`,-c(year,Total)) %>% max()) # what table values to color
  )) %>% # I bet there's a way to generate the number of colors and the domain more dynamically, but I didn't try to figure out how. As is I figured these out by just reading the values in the uncolored table.
  grand_summary_rows(columns = -year, # create a grand summary row, don't use year
                     fns = list(Total = ~sum(., na.rm = TRUE))) %>% # Sum over each column selected.
  # tab_options(data_row.padding = px(1)) %>% # reduce unnecessary spacing if you're saving the html to copy into a word doc.
  gtsave("figures/chicksfledgedbysite.png") # save as a png, result attached as well. Can also do pdf, jpg, and html. For word docs I save as html, open in browser, and then copy/paste the rendered table to word (although note that you can also just make word docs directly using rmarkdown so the document is fully dynamic. #goals for upcoming manuscripts). It's pretty faithful, just needs some slight adjustments. Apparently word developers put a lot of effort towards keeping html to word fairly painless so might as well take advantage.
