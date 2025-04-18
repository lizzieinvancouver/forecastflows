## Started 17 April 2025 ##
## By Lizzie so far ... ## 

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# set your directory if you're me
setwd("~/Documents/git/projects/misc/miscmisc/workflowsPhilTrans/analyses")

# packages 
library(data.table)
# What?! These two are because I was cheap and copied from Johnson et al. 2024
library(countrycode)
library(tidyverse)


lpi <- fread("input/LPD2022_public.csv", na = c("", "NA", "NULL"))

# Below is copied from :
# Johnson_etal_2024_code_GitTFJ-correlated_effect_model-21c0d95/code/data_compile.Rmd
custom_match <- c(`Cura??ao` = "CUW", `?Óland Islands` = "ALA", `International Waters` = "ABNJ")

lpi <- subset(lpi, Replicate == 0)
lpi <- lpi %>% 
  mutate(
    citation_year = str_replace_all(Citation, "\\\"|\\.|-|\\\'", ""),
    citation_year = parse_number(citation_year),
    site = paste(word(Citation, 1, sep = " |,"),
                 citation_year,
                 Location,
                 sep = "_")
  ) %>%
  select(site, Country, Latitude, Longitude, Binomial, Units, `1950`:`2020`) %>% 
  pivot_longer(cols = `1950`:`2020`,
               names_to = "date",
               values_to = "abundance")
lpi <- as.data.frame(lpi)


## Thinking we need to subset for sanity ... 
length(unique(lpi$site))
sort(unique(lpi$Country))
table(lpi$Country, lpi$site)

sa <- subset(lpi, Country=="South Africa")

ggplot(sa, aes(x=date, y=abundance)) +
	geom_point() + 
	facet_wrap(.~site)