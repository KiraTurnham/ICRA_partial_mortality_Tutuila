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
eds_YR01 <- eds %>%
  select(SITE, ends_with("jplMUR_Daily_YR01"))


eds_YR01 <- eds_YR01 %>%
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
         SST_Mean = mean_Sea_Surface_Temperature_jplMUR_Daily_YR01
         #SST_Q05 = q05_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         #SST_Q95 = q95_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         #SST_SD = sd_Sea_Surface_Temperature_jplMUR_Daily_MO06,
         # SST_BiweekRange = mean_biweekly_range_Sea_Surface_Temperature_jplMUR_Daily_MO06)
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
png("plots/corrplot_output.png", width = 800, height = 600, res = 150)
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
         # SST_AnnRange,
         SST_biweekly,
         SST_Mean,
         DHW_Dur,
         DHW_MeanMax
  )

#look at range
range(sub_eds$SST_Mean)
range(sub_eds$SST_biweekly)
range(sub_eds$DHW_Dur)
range(sub_eds$DHW_MeanMax)

#next is merging variables of interest back with PM colony data
######
merged2025_eds_PM_S_colony1yr_sep<-sub_eds%>%
  left_join(SOUTH_COLONY_SIZE_PM, by = "SITE")%>%
  filter(!is.na(PER_DEAD))
# write.csv(merged2025_eds_PM_S_colony1yr_DHW_sep, file ="merged2025_eds_PM_S_colony1yr_JPLDHW_Sep.csv") 
# save(merged2025_eds_PM_S_colony1yr_DHW_sep, file="paper/merged2025_eds_PM_S_colony1yr_JPLDHW_Sep.Rdata")

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
    #SST_AnnRange = first(SST_AnnRange),
    SST_biweekly = first(SST_biweekly),
    DHW_Dur = first(DHW_Dur),
    DHW_MeanMax = first(DHW_MeanMax),
    Latitude= first(LATITUDE),
    LONGITUDE = first(LONGITUDE),
    .groups = "drop")

#histogram
hist(site$prop_mean_PM, breaks = 30)

## check which model is best
glm.1 <- glmmTMB(prop_mean_PM ~ SST_Mean * TAIL_BINS + SST_biweekly * TAIL_BINS,
                 data = site,
                 family = beta_family()
)
glm.1a <- glmmTMB(prop_mean_PM ~ DHW_MeanMax * TAIL_BINS + DHW_Dur * TAIL_BINS,
                  data = site,
                  family = beta_family()
)
glm.1b <- glmmTMB(prop_mean_PM ~ SST_Mean * TAIL_BINS + DHW_Dur * TAIL_BINS + SST_biweekly * TAIL_BINS,
                  data = site,
                  family = beta_family() 
)        
glm.1c <- glmmTMB(prop_mean_PM ~ SST_Mean * TAIL_BINS,
                  data = site,
                  family = beta_family()
)
glm.1d <- glmmTMB(prop_mean_PM ~ DHW_MeanMax * TAIL_BINS,
                  data = site,
                  family = beta_family()
)
glm1e<-glmmTMB(prop_mean_PM ~ SST_biweekly * TAIL_BINS,
               data = site,
               family = beta_family()
)

glm.1f <- glmmTMB(prop_mean_PM ~ SST_biweekly * TAIL_BINS + DHW_MeanMax * TAIL_BINS,
                  data = site,
                  family = beta_family())


AIC(glm.1, glm.1a, glm.1b, glm.1c, glm.1d, glm1e, glm.1f)
compare_performance(glm.1, glm.1a, glm.1b, glm.1c, glm.1d, glm1e, glm.1f, rank = TRUE)

# best is just SST mean
model<- glmmTMB(
  prop_mean_PM ~ (SST_Mean * TAIL_BINS), #+ (SST_biweekly * TAIL_BINS),
  data = site,
  family = beta_family()
)
summary(model)

#Ferrari's R2
performance::r2(model)

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

# predict PM
# Create prediction sequences per size bin
make_newdata <- function(bin_label, site_data) {
  subset_data <- subset(site_data, TAIL_BINS == bin_label)
  SST_seq <- seq(min(subset_data$SST_Mean), max(subset_data$SST_Mean), length.out = 100)
  data.frame(
    SST_Mean = SST_seq,
    TAIL_BINS = bin_label
  )
}

newdata_small <- make_newdata("Q20", site)
newdata_med   <- make_newdata("QMED", site)
newdata_large <- make_newdata("Q80", site)

predict_and_bind <- function(newdata, label) {
  p <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
  p_df <- as.data.frame(p)
  colnames(p_df) <- c("Predicted_PM", "SE_PM")
  newdata <- cbind(newdata, p_df)
  newdata$Predict.lwr <- newdata$Predicted_PM - 1.96 * newdata$SE_PM
  newdata$Predict.upr <- newdata$Predicted_PM + 1.96 * newdata$SE_PM
  newdata$Size_cat <- label
  return(newdata)
}

newdata_small <- predict_and_bind(newdata_small, "Small")
newdata_med   <- predict_and_bind(newdata_med, "Medium")
newdata_large <- predict_and_bind(newdata_large, "Large")

all.newdata <- rbind(newdata_small, newdata_med, newdata_large)

# Reorder size factor
all.newdata$Size_cat <- factor(all.newdata$Size_cat, levels = c("Large", "Medium", "Small"))
site$Size_cat <- recode(site$TAIL_BINS,
                        "Q20" = "Small",
                        "QMED" = "Medium",
                        "Q80" = "Large")
site$Size_cat <- factor(site$Size_cat, levels = c("Large", "Medium", "Small"))


p<-ggplot() +
  geom_line(data = all.newdata, aes(x = SST_Mean, y = Predicted_PM, color = Size_cat), size = 1) +
  geom_ribbon(data = all.newdata, aes(x = SST_Mean, ymin = Predict.lwr, ymax = Predict.upr, fill = Size_cat), alpha = 0.1) +
  geom_point(data = site, aes(x = SST_Mean, y = prop_mean_PM, color = Size_cat), size = 2, alpha = 0.7) +
  theme_bw() +
  theme(
    plot.margin = unit(c(2, 1, 1, 1), "cm"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 14),
    text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  ylab("Predicted partial mortality \n(proportion)") +
  xlab(expression("SST mean ("*~degree*C*")"))+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)
ggsave("plots/SST_mean_regression.png", plot = p, width = 8, height = 6, dpi = 300)
ggsave("plots/SST_mean_regression.pdf", plot = p, width = 8, height = 6, dpi = 300)
ggsave("plots/SST_mean_regression.png")
