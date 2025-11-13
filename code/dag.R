##### DAG of relationships among measured variables
##### Author: Katherine Lauck
##### Last updated: 24 January 2022

#  set theme of all DAGs to `theme_dag()`
library(dagitty)
library(ggdag)
library(ggplot2)
library(tidyverse)
theme_set(theme_dag())

growth <- dagify(fitness ~ nestgrowth,
                 nestgrowth ~ foodprovis + neststress,
                 foodprovis ~ habitat + surftemp,
                 neststress ~ surftemp + foodprovis,
                 surftemp ~ climate + habitat,
                 labels = c("fitness" = "fitness",
                            "nestgrowth" = "nestling growth\n& survival",
                            "foodprovis" = "food provisioning",
                            "neststress" = "nestling stress",
                            "adultstress" = "adult stress",
                            "foodsupply" = "food supply",
                            "surftemp" = "surface temperature",
                            "habitat" = "land cover",
                            "climate" = "climate"),
                 exposure = c("surftemp","habitat"),
                 outcome = "nestgrowth")

ggdag_status(growth, text = FALSE, use_labels = 'label', layout = "nicely") + geom_dag_edges()
ggsave("figures/dag.png")

g <- dagify(fitness ~ nestgrowth,
            nestgrowth ~ s1_cort + nestcond + foodprovis + surftemp_humint_growth,
            s1_cort ~ site + surftemp + mother_s1 + habitat + foodprovis,
            nestcond ~ s1_cort + foodprovis,
            foodprovis ~ habitat + year + mother_s1 + site + surftemp_humint_foodprovis,
            mother_s1 ~ site + surftemp + attempt + habitat  + jul_date,
            surftemp ~ year + habitat + jul_date,
            attempt ~ jul_date,
            surftemp_humint_growth ~ surftemp + humidity,
            surftemp_humint_foodprovis ~ surftemp + humidity,
            exposure = c("s1_cort","foodprovis"),
            outcome = "nestgrowth")
ggdag_status(g,layout = "nicely") + geom_dag_edges()
adjSet_suff <- adjustmentSets(g,type = 'all')
adjSet_minimal <- adjustmentSets(g)


