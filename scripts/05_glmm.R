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
library(glmmTMB)
library(DHARMa)
library(performance)
library(emmeans)
library(ggeffects)

select = dplyr::select
rename  = dplyr::rename

load("data/eds_output.Rdata")
load(file ="data/ICRA_PM_SIZE_USE.Rdata")

#subset 2025 (can later add this to clean-colony level script)
s<-SOUTH_COLONY_SIZE_PM %>%
  filter( YEAR == '2025')

select = dplyr::select
rename  = dplyr::rename

#subset 2025 from eds
#save the colnames as a file for ease in viewing variable names
column_names <- colnames(eds)
column_names_df <- data.frame(column_names)

#view just the 6month 1km data 
YR01_columns <- column_names[grepl("_YR01$", column_names)]
YR01_columns_df<- data.frame(YR01_columns)
non_zero_YR01_columns <- YR01_columns[colSums(eds[, YR01_columns] != 0) > 0]
print(non_zero_YR01_columns)

colnames(eds)
#subset jpm 2025, MO06 data
# eds_YR01 <- eds %>%
#   select(SITE, ends_with("jplMUR_Daily_YR01"))


eds_YR01 <- eds %>%
  select(SITE,
         DHW_MeanMax = DHW.MeanMax_Degree_Heating_Weeks_jplMUR_Daily_YR01,
         # DHW_MeanMax_Major = DHW.MeanMax_Major_Degree_Heating_Weeks_jplMUR_Daily_YR01,
         DHW_Dur = DHW.MeanDur_Degree_Heating_Weeks_jplMUR_Daily_YR01,
         # DHW_Dur_Major = DHW.MeanDur_Major_Degree_Heating_Weeks_jplMUR_Daily_YR01,
         # DHW_Max_Major = DHW.MaxMax_Major_Degree_Heating_Weeks_jplMUR_Daily_YR01,
         # DHW_Max_Major5 = DHW.MaxMax_Major_Degree_Heating_Weeks_jplMUR_Daily_YR01,
         SST_biweekly = mean_biweekly_range_Sea_Surface_Temperature_jplMUR_Daily_YR01,
         #  SST_AnnRange = mean_annual_range_Sea_Surface_Temperature_jplMUR_Daily_YR01,
         #SST_MonthRange = mean_monthly_range_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         SST_Mean = mean_Sea_Surface_Temperature_jplMUR_Daily_YR01,
         #SST_Q05 = q05_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         #SST_Q95 = q95_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         #SST_SD = sd_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         # SST_BiweekRange = mean_biweekly_range_Sea_Surface_Temperature_jplMUR_Daily_MO06)
         CRW_DHW_MeanMax = DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR01,
         CRW_DHW_Dur = DHW.MeanDur_Major_Degree_Heating_Weeks_CRW_Daily_YR01
  )
#   select(-any_of(c("DHW.YearsToLast_Degree_Heating_Weeks_jplMUR_Daily_MO06", 
#                    "DHW.YearsToLast_Major_Degree_Heating_Weeks_jplMUR_Daily_MO06", 
#                    "DHW.CI95Max_Major_Degree_Heating_Weeks_jplMUR_Daily_MO06", 
#                    "DHW.CI95Max_Degree_Heating_Weeks_jplMUR_Daily_MO06",
#                    "DHW.Np10y_Degree_Heating_Weeks_jplMUR_Daily_MO06",             
#                    "DHW.Np10y_Major_Degree_Heating_Weeks_jplMUR_Daily_MO06")
# ))
#save(eds_MO06, file = "eds_MO06_JPL_updated_DHW.Rdata")

sub_numeric <- eds_YR01 %>%
  select(where(is.numeric))

#make a matrix
sub_numeric_matrix <- as.matrix(sub_numeric)
sub_numeric_matrix[!is.finite(sub_numeric_matrix)] <- NA
sub_numeric_matrix <- na.omit(sub_numeric_matrix)

M <- cor(sub_numeric_matrix, method = "pearson", use = "pairwise.complete.obs") #pearsons
res1 <- cor.mtest(sub_numeric_matrix, conf.level = 0.95)

#correlation plot showing the correlation coefficient
png("plots/corrplot_output_w_CRW.png", width = 800, height = 600, res = 150)
corrplot(M, p.mat = res1$p, 
         sig.level = 0.05, 
         insig = "p-value",
         order = 'hclust', 
         addrect = 2, 
         tl.srt = 45, 
         tl.cex = 0.6,  
         pch.cex = 0.8, 
         type = 'upper')
dev.off()
###############################################################################################
#subset variables you want to use:
sub_eds <- eds_YR01 %>%
  select(SITE,
         SST_biweekly,
         SST_Mean,
         DHW_MeanMax,
         DHW_Dur,
         CRW_DHW_MeanMax,
         CRW_DHW_Dur
  )

#look at range
range(sub_eds$SST_Mean)
range(sub_eds$SST_biweekly)
range(sub_eds$DHW_MeanMax_Major)
range(sub_eds$DHW_Dur)
range(sub_eds$CRW_DHW_Dur)
range(sub_eds$CRW_DHW_MeanMax)

#next is merging variables of interest back with PM colony data
######
merged2025_eds_PM_S_colony1yr_sep<-sub_eds%>%
  left_join(SOUTH_COLONY_SIZE_PM, by = "SITE")%>%
  filter(!is.na(PER_DEAD))

#################################################################################################################
#adjust 0s fo beta regression before averaging, then take mean PM, convert to proportion.

colnames(merged2025_eds_PM_S_colony1yr_sep)
n <- nrow(merged2025_eds_PM_S_colony1yr_sep)


site<-merged2025_eds_PM_S_colony1yr_sep%>%
  mutate(
    PER_DEAD_TRUE = (PER_DEAD * (n - 1) + 0.5) / n)%>%                       
  group_by(SITE, TAIL_BINS)%>%
  summarise(
    n = n(),
    mean_PM = mean(PER_DEAD_TRUE, na.rm = TRUE),
    prop_mean_PM = mean_PM / 100,
    SST_Mean = first(SST_Mean),
    SST_biweekly = first(SST_biweekly),
    DHW_MeanMax_Major = first(DHW_MeanMax),
    DHW_Dur = first(DHW_Dur),
    CRW_DHW_Dur = first(CRW_DHW_Dur),
    CRW_DHW_MeanMax = first(CRW_DHW_MeanMax),
    Latitude= first(LATITUDE),
    LONGITUDE = first(LONGITUDE),
    .groups = "drop")

#histogram
hist(site$prop_mean_PM, breaks = 30)

## check which model is best
glm.sstmean_biweekly <- glmmTMB(prop_mean_PM ~ SST_Mean * TAIL_BINS + SST_biweekly * TAIL_BINS,
                 data = site,
                 family = beta_family()
)

glm.sstmean <- glmmTMB(prop_mean_PM ~ SST_Mean * TAIL_BINS,
                  data = site,
                  family = beta_family()
)

glm.sstbiweekly<-glmmTMB(prop_mean_PM ~ SST_biweekly * TAIL_BINS,
               data = site,
               family = beta_family()
)

glm.dhwmeanmax<-glmmTMB(prop_mean_PM ~ DHW_MeanMax_Major * TAIL_BINS,
                        data = site,
                        family = beta_family()
)

glm.dur<-glmmTMB(prop_mean_PM ~ DHW_Dur * TAIL_BINS,
                 data = site,
                 family = beta_family()
)
glm.crw.dur<-glmmTMB(prop_mean_PM ~ CRW_DHW_Dur * TAIL_BINS,
                            data = site,
                            family = beta_family()
)
glm.crw.meanmax<-glmmTMB(prop_mean_PM ~ CRW_DHW_MeanMax * TAIL_BINS,
                            data = site,
                            family = beta_family()
)


AIC(glm.sstmean_biweekly, glm.sstmean, glm.sstbiweekly, glm.dhwmeanmax, glm.dur, glm.crw.dur, glm.crw.meanmax)
compare_performance(glm.sstmean_biweekly, glm.sstmean, glm.sstbiweekly, glm.dhwmeanmax, glm.dur, glm.crw.dur, glm.crw.meanmax, rank = TRUE)

# 
# AIC(glm.1, glm.1a, glm.1b)
# compare_performance(glm.1, glm.1a, glm.1b, rank = TRUE)

# best is just JPL SST mean
model<- glmmTMB(
  prop_mean_PM ~ (SST_Mean * TAIL_BINS), #+ (SST_biweekly * TAIL_BINS),
  data = site,
  family = beta_family()
)
summary(model)

#Ferrari's R2
performance::r2(model)



#emmeans estimated mean PM across SST for each size class
emmeans(model, ~ SST_Mean | TAIL_BINS, at = list(SST_Mean = seq(min(site$SST_Mean), max(site$SST_Mean), length.out = 5)))

#which sizes differ
mean_sst <- mean(site$SST_Mean, na.rm = TRUE)
# # Estimated means at mean SST
# emms_fixedSST <- emmeans(model, ~ TAIL_BINS)
# pairs(emms_fixedSST, adjust = "tukey")

#test if slope differs by size. estiamte slope of SST mean for each size class
emtrends(model, ~ TAIL_BINS, var = "SST_Mean")

pairs(emtrends(model, ~ TAIL_BINS, var = "SST_Mean"))

trends<- emtrends(model, specs = "TAIL_BINS", var = "SST_Mean")
summary(trends, infer = c(TRUE, TRUE))
pairs(trends)

# Pairwise tests of slope differences
pairs(trends)

############################################
#look at residuals. first simulate w DHARMa#
############################################

sim_res <- simulateResiduals(fittedModel = model, plot = TRUE)
#no significant deviations or problems detected
hist(sim_res$scaledResiduals, main = "Histogram of Scaled Residuals", xlab = "Residuals")

plotResiduals(sim_res, site$SST_Mean)

testUniformity(sim_res)     # Are residuals uniformly distributed?
#D = ns, KS test p = 0.53
testDispersion(sim_res)     # Is there overdispersion?
#p=0.528
testOutliers(sim_res)       # Are there extreme values?
#ns

#plot resudials by site
plotResiduals(sim_res, site$SITE)
#ns


################
#plot raw data##
################


site$Size_cat <- recode(site$TAIL_BINS,
                        "Q20" = "Small",
                        "QMED" = "Medium",
                        "Q80" = "Large")

colors<-c("goldenrod", "gold", "darkgrey" )

site$Size_cat <- factor(site$Size_cat, levels = c("Small", "Medium", "Large"))

#plot raw data
ggplot(site, aes(x = SST_Mean, y = prop_mean_PM, color = Size_cat)) +
  geom_point(color = "black", size = 2, alpha = 0.6) +
  geom_smooth(method = "lm",  se = TRUE) +
  #facet_wrap(~ Size_cat) +
  theme_bw() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    x = expression(bold("SST Mean ("*~degree*C*")")),
    y = expression(bold("Observed Partial Mortality")),
    # title = "Observed PM by DHW Mean Across Size Classes"
  ) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

# predict PM using ggpredict
site$TAIL_BINS <- factor(site$TAIL_BINS)
model <- glmmTMB(prop_mean_PM ~ SST_Mean * TAIL_BINS, family = beta_family, data = site)

preds <- ggpredict(model, terms = c("SST_Mean [all]", "TAIL_BINS"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, color = NA) +
  geom_point(data = site, aes(x = SST_Mean, y = prop_mean_PM, color = TAIL_BINS), 
             size = 2, alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    x = expression("SST mean ("*~degree*C*")"),
    y = "Predicted partial mortality (proportion)",
    color = "Size class",
    fill = "Size class"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 14),
    text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
      ylab("Predicted partial mortality \n(proportion)") 
  

#to plot scaled x axis

site$SST_scaled <- scale(site$SST_Mean)
model <- glmmTMB(prop_mean_PM ~ SST_scaled * TAIL_BINS, data = site, family = beta_family)

preds <- ggpredict(model, terms = c("SST_scaled [all]", "TAIL_BINS"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, color = NA) +
  geom_point(data = site, aes(x = SST_scaled, y = prop_mean_PM, color = TAIL_BINS), 
             size = 2, alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    x = expression("SST mean ("*~degree*C*")"),
    y = "Predicted partial mortality (proportion)",
    color = "Size class",
    fill = "Size class"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 14),
    text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
  ylab("Predicted partial mortality \n(proportion)")+
  xlab("Scaled SST Mean")
)

#run 2 category model:
#= smalls vs all (>12 cm), pooling med/large colonies
#adjust 0s fo beta regression before averaging, then take mean PM, convert to proportion.

colnames(merged2025_eds_PM_S_colony1yr_sep)
n <- nrow(merged2025_eds_PM_S_colony1yr_sep)

merged2025_eds_PM_S_colony1yr_sep<-merged2025_eds_PM_S_colony1yr_sep%>%
  mutate(size_cat = case_when(
    TAIL_BINS == "Q20" ~ "Q20",
    TAIL_BINS %in% c("Q80", "QMED") ~ "big"
  ))

site<-merged2025_eds_PM_S_colony1yr_sep%>%
  mutate(
    PER_DEAD_TRUE = (PER_DEAD * (n - 1) + 0.5) / n)%>%                       
  group_by(SITE, size_cat)%>%
  summarise(
    n = n(),
    mean_PM = mean(PER_DEAD_TRUE, na.rm = TRUE),
    prop_mean_PM = mean_PM / 100,
    SST_Mean = first(SST_Mean),
    SST_biweekly = first(SST_biweekly),
    DHW_MeanMax_Major = first(DHW_MeanMax),
    DHW_Dur = first(DHW_Dur),
    CRW_DHW_Dur = first(CRW_DHW_Dur),
    CRW_DHW_MeanMax = first(CRW_DHW_MeanMax),
    Latitude= first(LATITUDE),
    LONGITUDE = first(LONGITUDE),
    .groups = "drop")

#histogram
hist(site$prop_mean_PM, breaks = 30)

## check which model is best
glm.sstmean_biweekly <- glmmTMB(prop_mean_PM ~ SST_Mean * size_cat + SST_biweekly * size_cat,
                                data = site,
                                family = beta_family()
)

glm.sstmean <- glmmTMB(prop_mean_PM ~ SST_Mean * size_cat,
                       data = site,
                       family = beta_family()
)

glm.sstbiweekly<-glmmTMB(prop_mean_PM ~ SST_biweekly * size_cat,
                         data = site,
                         family = beta_family()
)

glm.dhwmeanmax<-glmmTMB(prop_mean_PM ~ DHW_MeanMax_Major * size_cat,
                        data = site,
                        family = beta_family()
)

glm.dur<-glmmTMB(prop_mean_PM ~ DHW_Dur * size_cat,
                 data = site,
                 family = beta_family()
)
glm.crw.dur<-glmmTMB(prop_mean_PM ~ CRW_DHW_Dur * size_cat,
                     data = site,
                     family = beta_family()
)
glm.crw.meanmax<-glmmTMB(prop_mean_PM ~ CRW_DHW_MeanMax * size_cat,
                         data = site,
                         family = beta_family()
)


AIC(glm.sstmean_biweekly, glm.sstmean, glm.sstbiweekly, glm.dhwmeanmax, glm.dur, glm.crw.dur, glm.crw.meanmax)
compare_performance(glm.sstmean_biweekly, glm.sstmean, glm.sstbiweekly, glm.dhwmeanmax, glm.dur, glm.crw.dur, glm.crw.meanmax, rank = TRUE)

# best is just JPL SST mean
model<- glmmTMB(
  prop_mean_PM ~ (SST_Mean * size_cat), #+ (SST_biweekly * TAIL_BINS),
  data = site,
  family = beta_family()
)
summary(model)

#Ferrari's R2
performance::r2(model)



#emmeans estimated mean PM across SST for each size class
emmeans(model, ~ SST_Mean | size_cat, at = list(SST_Mean = seq(min(site$SST_Mean), max(site$SST_Mean), length.out = 5)))

#which sizes differ
mean_sst <- mean(site$SST_Mean, na.rm = TRUE)
# # Estimated means at mean SST
# emms_fixedSST <- emmeans(model, ~ TAIL_BINS)
# pairs(emms_fixedSST, adjust = "tukey")

#test if slope differs by size. estiamte slope of SST mean for each size class
emtrends(model, ~ size_cat, var = "SST_Mean")

pairs(emtrends(model, ~ size_cat, var = "SST_Mean"))

trends<- emtrends(model, specs = "size_cat", var = "SST_Mean")
summary(trends, infer = c(TRUE, TRUE))
pairs(trends)

# predict PM using ggpredict
site$size_cat <- factor(site$size_cat)
model <- glmmTMB(prop_mean_PM ~ SST_Mean * size_cat, family = beta_family, data = site)

preds <- ggpredict(model, terms = c("SST_Mean [all]", "size_cat"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, color = NA) +
  geom_point(data = site, aes(x = SST_Mean, y = prop_mean_PM, color = size_cat), 
             size = 2, alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    x = expression("SST mean ("*~degree*C*")"),
    y = "Predicted partial mortality (proportion)",
    color = "Size class",
    fill = "Size class"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 14),
    text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
  ylab("Predicted partial mortality \n(proportion)") 
)

