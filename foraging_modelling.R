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

edges<- edges %>% 
  mutate(season_ID = paste(season, ID, sep = "_")) %>% 
  mutate(prop_time_center = time_center/time_total)
  
master$unique_trial_ID <- as.factor(master$unique_trial_ID)
master$trial <- as.factor(master$trial)
master$season <- as.factor(master$season)

#Ordering the factors in master
master$season <- factor(master$season, 
                        levels = c("summer", "winter", "spring"), 
                        ordered = TRUE)  
master$trial <- factor(master$trial,
                       levels = c("T1S1", "T1S2", "T2S1", "T2S2"),
                       ordered = TRUE)

#Models for exploration tendency using "master" and "edges"

#Model first island time (master). I want to test if the time to reach the first island changes among trials and season
hist(master$first_island_time)
hist(log(master$first_island_time + 1))

model_istime <- lmer(log(first_island_time) ~ season + trial + (1 | season_ID),data = master)
summary(model_istime)
plot(model_istime)
qqnorm(resid(model_istime))

ggplot(master, aes(x = season, y = log(first_island_time), color = trial)) +
  geom_boxplot()
  
#Model on the distance rate (master), to test if there is a change in distance covered among trials and seasons.
#Distance rate= distance covered in cm /time moving
hist(master$distance_rate)
hist(log(master$distance_rate+1))

model_dist <- lmer(log(distance_rate) ~ season + trial + (1 | season_ID), data = master)
summary(model_dist)
plot(model_dist)
qqnorm(resid(model_dist))
qqline(resid(model_dist))

ggplot(master, aes(x = season, y = distance_rate, color = trial)) +
  geom_boxplot()

#Model on time moving (master) among trial and season
hist(master$moving_time)
master$moving_time_sqrt <- sqrt(master$moving_time)
hist(master$moving_time_sqrt)

model_move <- lmer(moving_time_sqrt ~ trial + season + (1 | season_ID), data = master)
summary(model_move)
plot(model_move)
qqline(resid(model_move))

ggplot(master, aes(x = season, y = moving_time, color = trial)) +
  geom_boxplot()

#Model on the total number of visits at the islands (master) per trial and season
hist(master$tot_visit)
hist(log(master$tot_visit+1))

master$tot_visit_log <- log(master$tot_visit + 1)
model_visit <- lmer(tot_visit_log ~ trial + season + (1 | season_ID), data = master)
summary(model_visit)
plot(model_visit)
qqnorm(resid(model_visit))
qqline(resid(model_visit))

ggplot(master, aes(x = season, y = tot_visit_log, color = trial)) +
  geom_boxplot()

#Ordering the factors in edges
edges$season <- factor(master$season, 
                        levels = c("summer", "winter", "spring"), 
                        ordered = TRUE)  
edges$trial <- factor(master$trial,
                       levels = c("T1S1", "T1S2", "T2S1", "T2S2"),
                       ordered = TRUE)
#Models on the proportional time spent on the edge per trial and season
model_edges<- lmer(prop_time_edge ~ trial + season + (1 | season_ID), data = edges)
summary(model_edges)
plot(model_edges)
qqline(resid(model_edges))


hist(edges$prop_time_edges)
ggplot(edges, aes(x = trial, y = prop_time_edge))+
  geom_boxplot()+
  facet_wrap (~ season)


#Models for path consistency using "path"
