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

select = dplyr::select
rename  = dplyr::rename

load("data/eds_MO06_JPL_updated_DHW.Rdata")
load(file ="data/ICRA_PM_SIZE_USE.Rdata")

#subset 2025 (can later add this to clean-colony level script)
s<-SOUTH_COLONY_SIZE_PM %>%
  filter( YEAR == '2025')

select = dplyr::select
rename  = dplyr::rename


eds_MO06 <- eds_MO06 %>%
  select(SITE,
         SST_BiweekRange,
         SST_Mean,
         DHW_MeanMax_Major,
         DHW_Dur,
         CRW_DHW_MeanMax,
         CRW_DHW_Dur
  )
#CRW_DHW_MeanMax = DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR01,
#CRW_DHW_Dur = DHW.MeanDur_Major_Degree_Heating_Weeks_CRW_Daily_YR01

sub_numeric <- eds_MO06 %>%
  select(where(is.numeric))

#make a matrix
sub_numeric_matrix <- as.matrix(sub_numeric)
sub_numeric_matrix[!is.finite(sub_numeric_matrix)] <- NA
sub_numeric_matrix <- na.omit(sub_numeric_matrix)

M <- cor(sub_numeric_matrix, method = "pearson", use = "pairwise.complete.obs") #pearsons
res1 <- cor.mtest(sub_numeric_matrix, conf.level = 0.95)

#correlation plot showing the correlation coefficient
#png("plots/corrplot_output.png", width = 800, height = 600, res = 150)
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
sub_eds <- eds_MO06 

#look at range
range(sub_eds$SST_Mean)
range(sub_eds$SST_BiweekRange)
range(sub_eds$DHW_MeanMax_Major)
range(sub_eds$DHW_Dur)
range(sub_eds$CRW_DHW_Dur)
range(sub_eds$CRW_DHW_MeanMax)

#next is merging variables of interest back with PM colony data
######
merged2025_eds_PM_S_colony_6mjuly<-sub_eds%>%
  left_join(SOUTH_COLONY_SIZE_PM, by = "SITE")%>%
  filter(!is.na(PER_DEAD))

#################################################################################################################
#adjust 0s fo beta regression before averaging, then take mean PM, convert to proportion.

colnames(merged2025_eds_PM_S_colony_6mjuly)
n <- nrow(merged2025_eds_PM_S_colony_6mjuly)


site<-merged2025_eds_PM_S_colony_6mjuly%>%
  mutate(
    PER_DEAD_TRUE = (PER_DEAD * (n - 1) + 0.5) / n)%>%                       
  group_by(SITE, TAIL_BINS)%>%
  summarise(
    n = n(),
    mean_PM = mean(PER_DEAD_TRUE, na.rm = TRUE),
    prop_mean_PM = mean_PM / 100,
    SST_Mean = first(SST_Mean),
    SST_biweekly = first(SST_BiweekRange),
    #DHW_MeanMax_Major = first(DHW_MeanMax_Major),
    #DHW_Dur = first(DHW_Dur),
   # CRW_DHW_Dur = first(CRW_DHW_Dur),
  #  CRW_DHW_MeanMax = first(CRW_DHW_MeanMax),
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


# best is just SST mean
model<- glmmTMB(
  prop_mean_PM ~ (SST_Mean * TAIL_BINS), #+ (SST_biweekly * TAIL_BINS),
  data = site,
  family = beta_family()
)
summary(model)

#Ferrari's R2
performance::r2(model)



#emmeans across SST_Mean by size class
emmeans(model, ~ SST_Mean | TAIL_BINS, at = list(SST_Mean = seq(min(site$SST_Mean), max(site$SST_Mean), length.out = 5)))

#which sizes differ
mean_sst <- mean(site$SST_Mean, na.rm = TRUE)
# Estimated means at mean SST
emms_fixedSST <- emmeans(model, ~ TAIL_BINS)

# Pairwise comparisons
pairs(emms_fixedSST, adjust = "tukey")

#test if slope differs by size
emtrends(model, ~ TAIL_BINS, var = "SST_Mean")

pairs(emtrends(model, ~ TAIL_BINS, var = "SST_Mean"))


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

# # calculate percent change
# calculate_percent_change <- function(data, size_label) {
#   min_val <- min(data$Predicted_PM, na.rm = TRUE)
#   max_val <- max(data$Predicted_PM, na.rm = TRUE)
#   percent_change <- ((max_val - min_val) / min_val) * 100
#   cat(sprintf("Percent change in predicted PM for %s colonies: %.2f%%\n", size_label, percent_change))
# }
# 
# calculate_percent_change(newdata_med, "Medium")
# calculate_percent_change(newdata_large, "Large")

# Reorder size factor
all.newdata$Size_cat <- factor(all.newdata$Size_cat, levels = c("Large", "Medium", "Small"))
site$Size_cat <- recode(site$TAIL_BINS,
                        "Q20" = "Small",
                        "QMED" = "Medium",
                        "Q80" = "Large")
site$Size_cat <- factor(site$Size_cat, levels = c("Large", "Medium", "Small"))


p<-ggplot() +
  geom_line(data = all.newdata, aes(x = SST_Mean, y = Predicted_PM, color = Size_cat), linewidth = 1) +
  geom_ribbon(data = all.newdata, aes(x = SST_Mean, ymin = Predict.lwr, ymax = Predict.upr, fill = Size_cat), alpha = 0.1) +
  geom_point(data = site, aes(x = SST_Mean, y = prop_mean_PM, color = Size_cat), size = 2, alpha = 0.7) +
  theme_bw() +
  theme(
    plot.margin = unit(c(2, 1, 1, 1), "cm"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 18),
    text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  ylab("Predicted partial mortality \n(proportion)") +
  xlab(expression("SST mean ("*~degree*C*")"))+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)+
  ggtitle ( "6 month prior to July 1, 2024")
ggsave("plots/supp_SST_mean_regression_6m.png", plot = p, width = 8, height = 6, dpi = 300)
ggsave("plots/supp_SST_mean_regression_6m.pdf", plot = p, width = 8, height = 6, dpi = 300)
#ggsave("plots/supp_SST_mean_regression_6m.png")
