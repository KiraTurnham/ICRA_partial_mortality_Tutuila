#merging NCRMP site-level density with ESA site-level density using updated ICRA colony totals


#load libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggridges)

rm(list=ls())

####PREP NCMRP DATA---------------

#load data 
ICRA.dat <- read.csv("data/NCRMP_COlony_level_TUT_filtered.csv")%>% mutate_if(is.character,as.factor)%>%
  select(YEAR, DATE_, SITE, SPCODE, COLONYID)


site <- read.csv("data/CoralBelt_Adults_raw_CLEANED_2023.csv")%>% mutate_if(is.character,as.factor)%>%
  filter(ISLANDCODE == "TUT", REEF_ZONE == "Forereef", OBS_YEAR != "2020", DEPTH_BIN == "Mid")%>%
  rename(YEAR = OBS_YEAR)%>%
  select(YEAR, DATE_, SITE, LATITUDE, LONGITUDE, MAX_DEPTH_M, TRANSECTAREA)%>%
  distinct()%>%
  droplevels()

#sum transect area per site
site <- site %>% 
  group_by(YEAR, SITE, LATITUDE, LONGITUDE, MAX_DEPTH_M)%>%
  mutate(SURVEYAREA = sum(TRANSECTAREA))%>%
  select(-TRANSECTAREA)%>%
  distinct()
#115 total sites


#colonies per site
col <-  ICRA.dat%>%
  group_by(YEAR, SITE) %>%
  summarise(COL_COUNT = n_distinct(COLONYID))

#join colony and site data
NCRMP <- left_join(site, col)%>%
  replace(is.na(.), 0) %>%
mutate(
  DATE_ = as.character(DATE_),                     # factor → character
  DATETIME = ymd_hms(DATE_),                       # parse full datetime properly
  DATE_only = as.Date(DATETIME),                   # just the date
  TIME_only = format(DATETIME, "%H:%M:%S"),        # just the time
  DATE_formatted = format(DATE_only, "%m/%d/%Y") # custom format
)

NCRMP$YEAR <- as.factor(NCRMP$YEAR)

####read in and merge 2025 ESA data to NCRMP data. add year column
ESA <- read_csv("data/ESA_Corals_Site_Density_2025.csv") %>% mutate_if(is.character,as.factor)
as.factor(ESA$YEAR <- "2025")

#rename columns to match
ESA <- ESA %>% 
  rename(DATE_formatted=Date, SITE=Site, LATITUDE =Lat, LONGITUDE = Long, SURVEYAREA = Survey_area, MAX_DEPTH_M = Max_depth_m)%>%
select(YEAR, DATE_formatted, SITE, LATITUDE, LONGITUDE, MAX_DEPTH_M, SURVEYAREA, COL_COUNT)


NCRMP <- NCRMP %>% select(YEAR, DATE_formatted, SITE, LATITUDE, LONGITUDE, MAX_DEPTH_M, SURVEYAREA, COL_COUNT)

#double check column names are the same before merging
colnames(NCRMP)
colnames(ESA)

SITE_COUNT <- rbind(NCRMP, ESA)

#calculate density by dividing colonies by survey area to get colonies per m2
ALL_COLONY_DENSITY<- SITE_COUNT %>%
  mutate(DENSITY=(COL_COUNT/SURVEYAREA))%>%
  ungroup()

save(ALL_COLONY_DENSITY, file="data/ALL_COLONY_DENSITY.RData")
save(ALL_COLONY_DENSITY, file="data/ALL_COLONY_DENSITY.csv")

#subset these data to only include data points from South side of Tutuila to remove surveying bias (ICRA is rare on north side)
#this was done by importing this data into ArcGIS and selectign points manually. See map script for visualization.
south_den<-read.csv("data/South_ICRA_Site_Counts.csv")


#calculate density
SOUTH_COLONY_DENSITY<- south_den %>%
  mutate(DENSITY=(COL_COUNT/SURVEYAREA))%>%
    ungroup()

#load in the colonies that were too big to be counted (from colony-level script)
load("data/colonies_removed_due_to_size.RData")

SOUTH_COLONY_DENSITY_filtered<- SOUTH_COLONY_DENSITY %>%
  left_join(removed_summary, by = "SITE") %>%
  mutate(
    adjusted_colony_count = COL_COUNT - coalesce(removed_count, 0),
    adjusted_density = adjusted_colony_count / SURVEYAREA
  )

#Use adjusted_density in future analyses:
save(SOUTH_COLONY_DENSITY_filtered, file="data/SOUTH_COLONY_DENSITY_filtered.RData")
write.csv(SOUTH_COLONY_DENSITY_filtered, file="data/SOUTH_COLONY_DENSITY_filtered.csv")
