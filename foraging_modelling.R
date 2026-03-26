# Load packages
library(lme4)
library(lmerTest)
library(here)
library(mgcv)
library(tidyverse)


master <- read.csv(here("csv/processed/foraging_master.csv"))

master$unique_trial_ID <- as.factor(master$unique_trial_ID)
master$trial <- as.factor(master$trial)
master$season <- as.factor(master$season)


model1 <- glmer(first_island_time ~ trial + season + (1 | season_ID), data = master, family = Gamma)
summary(model1)
# Residual plot
plot(model1)
# Normality of residuals
qqnorm(resid(model1))
qqline(resid(model1))
hist(master$first_island_time)
summary(master$first_island_time)
min(master$first_island_time, na.rm=TRUE)
ggplot(master, aes(x = season, y = first_island_time, color = trial)) +
  geom_boxplot()
  
# Total distance
model_dist <- lmer(distance_rate ~ trial + season + (1 | season_ID), data = master)
summary(model_dist)
plot(model_dist)
qqnorm(resid(model_dist))
qqline(resid(model_dist))

# Time moving
model_move <- lmer(moving_time ~ trial + season + (1 | season_ID), data = master)
summary(model_move)
plot(model_move)
qqnorm(resid(model_move))
qqline(resid(model_move))

#tot visits
model_visit <- lmer(tot_visit ~ trial + season + (1 | season_ID), data = master)
summary(model_visit)
plot(model_visit)
qqnorm(resid(model_visit))
qqline(resid(model_visit))
