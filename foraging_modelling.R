# Load packages
library(lme4)
library(lmerTest)
library(here)
library(mgcv)
library(tidyverse)

master <- read.csv(here("csv/processed/foraging_master.csv"))
path <- read.csv(here("csv/processed/foraging_similarities.csv"))
efficiency <- read.csv(here("csv/processed/interactions_counts.csv"))
edges <- read.csv(here("csv/foraging_edges_new.csv"))

master$unique_trial_ID <- as.factor(master$unique_trial_ID)
master$trial <- as.factor(master$trial)
master$season <- as.factor(master$season)

#Ordering the factors
master$season <- factor(master$season, 
                        levels = c("summer", "winter", "spring"), 
                        ordered = TRUE)  
#master$trial <- factor(master$trial,
                       #levels = c("T1S1", "T1S2", "T2S1", "T2S2"),
                       #ordered = TRUE)

#Models for exploration tendency using "master" and "edges"

#Model first island time (master). I want to test if the time to reach the first island changes among trials and season
hist(master$first_island_time)

model_istime <- lmer(first_island_time ~ trial + season + (1 | season_ID),data = master)
summary(model_istime)
plot(model_istime)
qqline(resid(model))

ggplot(master, aes(x = season, y = first_island_time, color = trial)) +
  geom_boxplot()
  
#Model on the distance rate (master), to test if there is a change in distance covered among trials and seasons.
#Distance rate= distance covered in cm /time moving
hist(master$distance_rate)

model_dist <- lmer(distance_rate ~ trial + season + (1 | season_ID), data = master)
summary(model_dist)
plot(model_dist)
qqnorm(resid(model_dist))
qqline(resid(model_dist))

ggplot(master, aes(x = season, y = distance_rate, color = trial)) +
  geom_boxplot()

#Model on time moving (master) among trial and season
hist(master$moving_time)

model_move <- lmer(moving_time ~ trial + season + (1 | season_ID), data = master)
summary(model_move)
plot(model_move)
qqline(resid(model_move))

ggplot(master, aes(x = season, y = moving_time, color = trial)) +
  geom_boxplot()

#Model on the total number of visits at the islands (master) per trial and season
hist(master$tot_visit)

model_visit <- lmer(tot_visit ~ trial + season + (1 | season_ID), data = master)
summary(model_visit)
plot(model_visit)
qqnorm(resid(model_visit))
qqline(resid(model_visit))

ggplot(master, aes(x = season, y = tot_visit, color = trial)) +
  geom_boxplot()

#Models for path consistency using "path"

hist(path$mean_sim_noncons)
hist(path$sd_noncons)
