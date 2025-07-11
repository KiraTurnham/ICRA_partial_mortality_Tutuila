#Cleaning ICRA data colony level
rm(list=ls())

#load libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(readxl)
library(dplyr)


####PREP DATA---------------

#load data and filter data 

ncrmp <- read.csv("data/CoralBelt_Adults_raw_CLEANED_2023.csv")%>% mutate_if(is.character,as.factor) %>%  
  filter(ISLANDCODE == "TUT", REEF_ZONE == "Forereef", OBS_YEAR != "2020", DEPTH_BIN == "Mid", SPCODE %in% c("ISSP", "ICRA")) %>%
  filter(!MORPHOLOGY %in% c("Branching", "Columnar")) %>%
  mutate(PER_DEAD = OLDDEAD + RDEXTENT1 + RDEXTENT2,
         Area_surveyed_m2 = 10) %>%
  # filter(as.numeric(COLONYLENGTH) > 4.9) %>% #only use colonies >5
  rename(YEAR = OBS_YEAR)%>%
  droplevels()

#confirm no corals below 5cm are included. And note largest corals measured by NCRMP, which will be used as cutoff for 2025 (to correct for survey differences biasing larger colonies) 
ncrmp %>%
  group_by(YEAR) %>%
  summarise(
    min_size = min(COLONYLENGTH, na.rm = TRUE),
    max_size = max(COLONYLENGTH, na.rm = TRUE)
  )

#Calculate quintiles of pre-bleaching year data to assign size classes 
#check normality
shapiro.test(ncrmp$COLONYLENGTH) #not normal
#log transform colony length, check normality
ncrmp$logCOLONYLENGTH = log(ncrmp$COLONYLENGTH)
shapiro.test(ncrmp$logCOLONYLENGTH) #not normal even after transformation. Use untransformed data
q<- data.frame(quantile(ncrmp$COLONYLENGTH, probs =  c(20,80)/100, na.rm = FALSE, names = TRUE, type = 9, digits = 4))
#(5-12 cm = "small", >40 = "brood stock")

#add quintiles to dataframe in new column TAIL_BINS
#ncrmp$logTAIL_BINS=cut(ncrmp$logCOLONYLENGTH,c(-Inf,q1[1,1],q1[2,1],Inf),labels=c('Q20','QMED','Q80'))
ncrmp$TAIL_BINS=cut(ncrmp$COLONYLENGTH,c(-Inf,q[1,1],q[2,1],Inf),labels=c('Q20','QMED','Q80'))

#note: tail bins were the same log tansforming or not. 

#write.csv(ncrmp, "NCRMP_COlony_level_TUT_filtered.csv")

#make sure columns are read as correct format
ncrmp$YEAR <- as.factor(ncrmp$YEAR)
ncrmp$COLONYLENGTH <- as.numeric(ncrmp$COLONYLENGTH)
ncrmp$PER_DEAD <- as.numeric(ncrmp$PER_DEAD)
ncrmp$DATE_ <- as.Date(ncrmp$DATE_)
ncrmp2 <- dplyr::select(ncrmp, DATE_, COLONYLENGTH, Area_surveyed_m2, MAX_DEPTH_M, SITE, PER_DEAD, LATITUDE, LONGITUDE, YEAR, TAIL_BINS)


#read in 2025 survey data
ICRA_2025 <- read.csv("data/2025_ICRA_colony_level_TUT.csv") %>%
  mutate(ICRA_size_cm = as.numeric(ICRA_size_cm),
         ICRA_partial_mortality = as.numeric(ICRA_partial_mortality),
         Year=as.factor(Year),
         Date = as.Date(Date),
         MAX_DEPTH_M = MAX_depth_ft * 0.3048) %>%
  rename(PER_DEAD = ICRA_partial_mortality,
         COLONYLENGTH = ICRA_size_cm,
         YEAR = Year,
         LATITUDE = Lat,
         LONGITUDE = Long,
         SITE = Site,
         DATE_=Date)

#apply quintiles from 2015 calculations
ICRA_2025$TAIL_BINS <- cut(
  ICRA_2025$COLONYLENGTH, 
  breaks = c(-Inf, 12, 40, Inf), 
  labels = c('Q20', 'QMED', 'Q80'))

#subset the data in order to merge w NCRMP
esa <- dplyr::select(ICRA_2025, DATE_, COLONYLENGTH, MAX_DEPTH_M, Area_surveyed_m2, SITE, PER_DEAD, LATITUDE, LONGITUDE, YEAR, TAIL_BINS)

#merge ncmrp and esa data: double check colnames are same   
colnames(esa)
colnames(ncrmp2)
ALL_ICRA_SIZE_PM <- rbind(
  mutate(esa, Data_Source = "esa"),      
  mutate(ncrmp2, Data_Source = "ncrmp2") 
) %>%
  mutate(
    Bleaching_Period = ifelse(YEAR %in% c(2015, 2018, 2023), "Pre-Bleaching", "Post-Bleaching")
  )
ALL_ICRA_SIZE_PM$YEAR <- ordered(ALL_ICRA_SIZE_PM$YEAR, levels = c("2015", "2018", "2023", "2025"))

#save this raw data
save(ALL_ICRA_SIZE_PM, file ="data/ICRA_PM_NS.csv")
save(ALL_ICRA_SIZE_PM, file ="data/All_ICRA_SIZE_PM.RData")

#manually removed all ICRA data from north side of island in arcGIS. both are plotted in map script for transparency and visualization.
SOUTH_COLONY_SIZE_PM<-read.csv("data/south_only_ICRA_Colony_level_data.csv")

save(south_ICRA_survey_data, file ="data/ICRA_PM_SIZE_USE.Rdata") #not for temporal comparisons

#identify corals outside the size range for temporal analyses (>116cm, to account for differences in survey methods in 2025 vs. ncrmp)
removed_summary <- SOUTH_COLONY_SIZE_PM %>% 
  filter(COLONYLENGTH < 4.9 | COLONYLENGTH > 116.1) %>%
  group_by(SITE) %>%
  summarise(
    removed_count = n(),
    min_size = min(COLONYLENGTH, na.rm = TRUE),
    max_size = max(COLONYLENGTH, na.rm = TRUE)
  ) %>%
  arrange(desc(removed_count))

#view how many were removed and from what site
print(removed_summary)
#save(removed_summary, file = "data/colonies_removed_due_to_size.RData")

#remove those corals from dataset
ICRA_SIZE_PM_SOUTH_filtered <- SOUTH_COLONY_SIZE_PM %>%
  filter(COLONYLENGTH >= 4.9 & COLONYLENGTH <= 116.1)

#save as csv and rdata
save(ICRA_SIZE_PM_SOUTH_filtered, file ="data/ICRA_SIZE_PM_SOUTH_sizefiltered_colony_level.RData") #remember this has subset of 2025 data that were sized
write.csv(ICRA_SIZE_PM_SOUTH_filtered, "data/south_ICRA_Colony_level_data_sizefiltered.csv", row.names = FALSE)