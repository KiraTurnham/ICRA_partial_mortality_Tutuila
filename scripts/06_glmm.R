rm(list = ls())
library(tidyverse)
library(corrplot)   
library(GGally)    
library(betareg)
library(statmod)
library(lmtest)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(broom)
library(tidyverse)

#exploring if spatial variation in SST mean and duration correlates w PM
#adding adjusitng based on PALS (y * (n - 1) + 0.5) / n to the COLONIES 
load("paper/merged2025_eds_PM_S_colony.Rdata")
load("paper/merged2025_eds_PM_S_colony6m.Rdata")

sep<-merged2025_eds_PM_S_colony
july<-merged2025_eds_PM_S_colony6m


colnames(merged2025_eds_PM_S_colony6m)
n <- nrow(merged2025_eds_PM_S_colony6m)
#adjust 0s fo beta regression before averaging, then take mean PM, convert to proportion.
site<-merged2025_eds_PM_S_colony6m%>%
  mutate(
    PER_DEAD_TRUE = (PER_DEAD * (n - 1) + 0.5) / n)%>%                       
  group_by(SITE, TAIL_BINS)%>%
  summarise(
    n = n(),
    mean_PM = mean(PER_DEAD_TRUE, na.rm = TRUE),
    prop_mean_PM = mean_PM / 100,
    SST_Mean = first(SST_Mean),
    SST_range = first(SST_range),
    Latitude= first(LATITUDE),
    LONGITUDE = first(LONGITUDE),
    .groups = "drop")

range(site$prop_mean_PM)
range(site$SST_Mean)
range(site$SST_range)

hist(site$prop_mean_PM, breaks = 30)

model<- glmmTMB(
  prop_mean_PM ~ (SST_Mean * TAIL_BINS) + (SST_range * TAIL_BINS),
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
#plotResiduals(sim_res, site2$SST_range)

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
#colors<-c("cyan4","purple3","gray67")
#colors<-c( "darkgrey", "#E1AD01",  "gold")
#colors<-c("grey", "orangered", "darkred" )
colors<-c("goldenrod", "gold", "darkgrey" )

site$Size_cat <- factor(site$Size_cat, levels = c("Small", "Medium", "Large"))

#plot
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
    # title = "Observed PM by SST Mean Across Size Classes"
  ) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

############################
## Visualize interactions ##
############################


####Step 1: create 3 new data frames for time since heat stress interval (or size bin our case) and hold everything but variable of interest constant####
#set SST range to mean (holding it constant)
small_data <- subset(site, TAIL_BINS == "Q20")
newdata_small<-small_data
newdata_small$SST_range <- mean(small_data$SST_range)  # Set SST_range to mean
newdata_small$SST_mean <- seq(min(small_data$SST_Mean), 
                              max(small_data$SST_Mean), 
                              length.out = nrow(newdata_small)) # making the sequence matches the number of rows. courtney used by = round(diff(range(small_data$SST_Mean)), 3) / nrow(small_data))

med_data <- subset(site, TAIL_BINS == "QMED")
newdata_med<-med_data
newdata_med$SST_range <- mean(med_data$SST_range)  # Set SST_range to mean
newdata_med$SST_mean <- seq(min(med_data$SST_Mean),
                            max(med_data$SST_Mean),
                            length.out = nrow(newdata_med))


large_data <- subset(site, TAIL_BINS == "Q80")
newdata_large<-large_data
newdata_large$SST_range <- mean(large_data$SST_range)  # Set SST_range to mean
newdata_large$SST_mean <- seq(min(large_data$SST_Mean), 
                              max(large_data$SST_Mean), 
                              length.out = nrow(newdata_large))

####Step 2: predict PM as a function of SST_mean using the model output (e.g. model in this case)####
#this would be predicted PM , create lower and upper CI
p <- predict(model, newdata = newdata_small, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_PM","SE_PM")
newdata_small<-cbind(newdata_small,p)
newdata_small$Predict.lwr <- newdata_small$Predicted_PM - 1.96 * newdata_small$SE_PM # confidence interval upper bound
newdata_small$Predict.upr <- newdata_small$Predicted_PM + 1.96 * newdata_small$SE_PM # confidence interval lower bound
newdata_small$Size_cat<-"Small"

p <- predict(model, newdata = newdata_med, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_PM","SE_PM")
newdata_med<-cbind(newdata_med,p)
newdata_med$Predict.lwr <- newdata_med$Predicted_PM - 1.96 * newdata_med$SE_PM # confidence interval upper bound
newdata_med$Predict.upr <- newdata_med$Predicted_PM + 1.96 * newdata_med$SE_PM # confidence interval lower bound
newdata_med$Size_cat<-"Medium"

p <- predict(model, newdata = newdata_large, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_PM","SE_PM")
newdata_large<-cbind(newdata_large,p)
newdata_large$Predict.lwr <- newdata_large$Predicted_PM - 1.96 * newdata_large$SE_PM # confidence interval upper bound
newdata_large$Predict.upr <- newdata_large$Predicted_PM + 1.96 * newdata_large$SE_PM # confidence interval lower bound
newdata_large$Size_cat<-"Large"

names(newdata_small)
names(newdata_med)
names(newdata_large)
#Merge into 1 dataframe and add column for each size category.  THIS didnt work bc diff rows : 
all.newdata<-rbind(newdata_small,newdata_med,newdata_large)


#Reorder size variables
all.newdata$Size_cat <- factor(all.newdata$Size_cat, levels = c("Large","Medium","Small"))
all.newdata<- all.newdata[order(all.newdata$Size_cat),];head(all.newdata)


#ver1
ggplot() +
  geom_line(data = all.newdata, aes(x = SST_Mean, y = Predicted_PM, color = Size_cat), size = 1) +
  #geom_ribbon(data = all.newdata, aes(x = SST_Mean, ymin = Predict.lwr, ymax = Predict.upr), alpha = 0.1) +
  geom_ribbon(data = all.newdata, aes(x = SST_Mean, ymin = Predict.lwr, ymax = Predict.upr, fill = Size_cat), alpha = 0.1) +
  #facet_wrap(~Size_cat)+
  geom_point(data = site, aes(x = SST_Mean, y = prop_mean_PM, color = Size_cat), size = 2, alpha = 0.7) +
  #geom_rug(data = all.newdata, mapping = aes(x = SST_Mean, y = 0, color = Size_cat), sides = "b", size = 1.5)+
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.7, 'cm'),
    legend.text = element_text(size = 14),
    text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  ) +
  ylab(expression(bold("Predicted PM"))) +
  xlab(expression(bold("SST Mean ("*~degree*C*")"))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)


ggsave("paper/regression.png")
