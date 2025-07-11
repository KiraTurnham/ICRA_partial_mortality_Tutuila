rm(list = ls())
library(tidyverse)
library(corrplot)   
library(GGally)    
library(eds)
library(betareg)
library(statmod)
library(lmtest)
library(ggplot2)
library(dplyr)

select = dplyr::select
rename  = dplyr::rename

#load("data/eds_output.Rdata")
load("eds_output_Sep12025.Rdata")
load("eds_output_Sep12025_UPDATEDJPL.Rdata")
load(file ="data/ICRA_PM_SIZE_USE.Rdata")
merged_PM_site_all_YR01.csv

#subset 2025 (can later add this to clean-colony level script)
s<-south_ICRA_survey_data %>%
  filter( YEAR == '2025')




#subset variables you want to use:
sub_eds_1yr <- eds %>%
  select(SITE,
         SST_range = mean_annual_range_Sea_Surface_Temperature_jplMUR_Daily_YR01,
         SST_Mean = mean_Sea_Surface_Temperature_jplMUR_Daily_YR01
         DHW_meandur = DHW.MeanDur_Major_Degree_Heating_Weeks_jplMUR_Daily_YR01
         DHW_meanmax = DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR01
         DHW_mean
         
  )

range(sub_eds$SST_Mean)
range(sub_eds$SST_range)

# year_eds <- eds %>%
#   select(SITE,
#          SST_range = mean_annual_range_Sea_Surface_Temperature_jplMUR_Daily_,
#          SST_Mean = mean_Sea_Surface_Temperature_jplMUR_Daily_
#   )
# 
# range(year_eds$SST_Mean)
# range(year_eds$SST_range)

#next is merging variables of interest back with PM colony data
merged2025_eds_PM_S_colony6m<-sub_eds%>%
  left_join(south_ICRA_survey_data, by = "SITE")%>%
  filter(!is.na(PER_DEAD))
write.csv(merged2025_eds_PM_S_colony6m, file ="merged2025_eds_PM_S_colony6mUPDATEDJPL.csv") 
save(merged2025_eds_PM_S_colony6m, file="paper/merged2025_eds_PM_S_colony6mUPDATEDJPL.Rdata")



#check if SST mean and range are not correlated
sub_numeric <- sub_eds %>%
  select(where(is.numeric))

#make a matrix
sub_numeric_matrix <- as.matrix(sub_numeric)
sub_numeric_matrix[!is.finite(sub_numeric_matrix)] <- NA
sub_numeric_matrix <- na.omit(sub_numeric_matrix)

#correlation plot 
M <- cor(sub_numeric_matrix, use = "pairwise.complete.obs") #pearsons
res1 <- cor.mtest(sub_numeric_matrix, conf.level = 0.95)
corrplot(M, p.mat = res1$p, 
         sig.level = 0.05, 
         insig = "p-value", # if you want to print nonsig pvalues
         order = 'hclust', 
         addrect = 2, 
         tl.srt = 45, 
         tl.cex = 0.6,  
         pch.cex = 0.8, 
         type = 'upper')

