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
t<-read.csv("data/CoralBelt_Adults_raw_CLEANED_2023.csv")

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
#   YEAR min_size max_size
#  2015        5       73
#  2018        5       70
#  2023        5      116

#Calculate quintiles of pre-bleaching year data to assign size classes 
#log transorm colony length, check normality
ncrmp$logCOLONYLENGTH = log(ncrmp$COLONYLENGTH)

shapiro.test(ncrmp$logCOLONYLENGTH)
#W = 0.98347, p-value = 0.002419

#q1<- data.frame(quantile(ncrmp$logCOLONYLENGTH, probs =  c(20,80)/100, na.rm = FALSE, names = TRUE, type = 9, digits = 4))

q<- data.frame(quantile(ncrmp$COLONYLENGTH, probs =  c(20,80)/100, na.rm = FALSE, names = TRUE, type = 9, digits = 4))
#(5-12 cm = "small", >40 = "brood stock")

#add quintiles to dataframe in new column TAIL_BINS
#ncrmp$logTAIL_BINS=cut(ncrmp$logCOLONYLENGTH,c(-Inf,q1[1,1],q1[2,1],Inf),labels=c('Q20','QMED','Q80'))
ncrmp$TAIL_BINS=cut(ncrmp$COLONYLENGTH,c(-Inf,q[1,1],q[2,1],Inf),labels=c('Q20','QMED','Q80'))

#tail bins were the same log tansforming or not. 

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

esa <- dplyr::select(ICRA_2025, DATE_, COLONYLENGTH, MAX_DEPTH_M, Area_surveyed_m2, SITE, PER_DEAD, LATITUDE, LONGITUDE, YEAR, TAIL_BINS)

#merge ncmrp and esa data    
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

summary_by_year_site <- ALL_ICRA_SIZE_PM %>%
  group_by(YEAR) %>%
  summarise(num_sites = n_distinct(SITE))

save(ALL_ICRA_SIZE_PM, file ="data/ICRA_PM_NS.RData")

ICRA_2025_unfiltered<-ALL_ICRA_SIZE_PM%>%
  filter(YEAR=="2025")

all_ICRA_PM_site<-ICRA_2025_unfiltered%>%
  group_by(SITE)%>%
  summarise(date = min(DATE_), #keep date column
            mean_PM = mean(PER_DEAD, na.rm = TRUE),
            sd_PM = sd(PER_DEAD, na.rm = TRUE),
            n = sum(!is.na(PER_DEAD)),
            se = sd_PM / sqrt(n),
            .groups = "drop")


save(all_ICRA_PM_site, file ="data/ICRA_2025_PM.RData")

save(ALL_ICRA_SIZE_PM, file ="data/All_ICRA_SIZE_PM.RData")
#write.csv(ALL_ICRA_SIZE_PM, "data/All_ICRA_Colony_level_data_filtered.csv", row.names = FALSE)



#manually removed all ICRA data from north side of island in arcGIS. both are plotted in map script for transparency and visualization.
SOUTH_COLONY_SIZE_PM<-read.csv("~/GitHub/ICRA/data/south_only_ICRA_Colony_level_data.csv")

south_ICRA_survey_data <- SOUTH_COLONY_SIZE_PM %>%
  filter(COLONYLENGTH >= 4.9)
save(south_ICRA_survey_data, file ="data/ICRA_PM_SIZE_USE.Rdata")

#identify corals outside the size range (<5cm or >116cm, to account for differences in survey methods in 2025 vs. ncrmp) (to update density dataset)
#this automaticall excludes any feb data (sicne we didn't measure size)
removed_data <- SOUTH_COLONY_SIZE_PM %>% 
  filter(!(as.numeric(COLONYLENGTH) > 4.9 )) %>%
  dplyr::select(SITE, COLONYLENGTH)

removed_summary <- SOUTH_COLONY_SIZE_PM %>% 
  filter(COLONYLENGTH < 4.9 | COLONYLENGTH > 116.1) %>%
  group_by(SITE) %>%
  summarise(
    removed_count = n(),
    min_size = min(COLONYLENGTH, na.rm = TRUE),
    max_size = max(COLONYLENGTH, na.rm = TRUE)
  ) %>%
  arrange(desc(removed_count))

# View the summary
print(removed_summary)
#save(removed_summary, file = "data/colonies_removed_due_to_size.RData")


#remove those corals from dataset
ICRA_SIZE_PM_SOUTH_filtered <- SOUTH_COLONY_SIZE_PM %>%
  filter(COLONYLENGTH >= 4.9 & COLONYLENGTH <= 116.1)

save(ICRA_SIZE_PM_SOUTH_filtered, file ="data/ICRA_SIZE_PM_SOUTH_sizefiltered.RData") #remember this has subset of 2025 data that were sized
write.csv(ICRA_SIZE_PM_SOUTH_filtered, "data/south_ICRA_Colony_level_data_sizefiltered.csv", row.names = FALSE)

#calculate average PM per site
avg_PM_south2025_nofeb<-ICRA_SIZE_PM_SOUTH_filtered%>%
  filter(YEAR==2025)%>%
  group_by(SITE) %>%
  summarise(
    n = sum(!is.na(PER_DEAD)),
    mean_PM = mean(PER_DEAD, na.rm = TRUE),
    sd_PM = sd(PER_DEAD, na.rm = TRUE),
    max_PM = max(PER_DEAD, na.rm = TRUE)
  )

save(avg_PM_south2025_nofeb, file ="data/ICRA_2025_SIZE_PM_nofeb.RData")
write.csv(avg_PM_south2025_nofeb, "data/south_ICRA_Colony_level_2025_data_filtered_nofeb.csv", row.names = FALSE)