# Load packages
library(lme4)
library(lmerTest)
library(here)
library(mgcv)
library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(glmmTMB)
library(ggeffects)
library(emmeans)

#Frequentist modelling ----
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

#Models using interaction counts to investigate if shrews have a strategy on door visiting
#Rearrange the data
# Door columns
door_cols <- c(
  "A_1","A_2","A_3","A_4","A_5","A_6",
  "B_1","B_2","B_3","B_4","B_5","B_6",
  "C_1","C_2","C_3","C_4","C_5","C_6",
  "D_1","D_2","D_3","D_4","D_5","D_6")

# Convert to long format
df_long <- pivot_longer(data = efficiency, cols = all_of(door_cols), names_to = "door", values_to = "visits")

# Remove missing values
df_long <- na.omit(df_long)

# Factors
df_long$trial <- factor(df_long$trial)
df_long$season <- factor(df_long$season)
df_long$ID <- factor(df_long$ID)
df_long$door <- factor(df_long$door)

#Model 
#Whether visit counts differ across trials, and whether this pattern is different in different seasons. 
#It does NOT tell about learning or strategy, it just describes the counts.
df_long$season <- factor(df_long$season, levels = c("summer", "winter", "spring"))

model_nb <- glmmTMB(visits ~ trial * season + (1 | ID), data = df_long, family = nbinom2) 
#Negative binomial: visit counts can be zero, and the variance is much larger than the mean (overdispersion — some doors get visited many times, most get zero)
summary(model_nb)

#This should compare the number of visits in each trial and check for significant differences
emmeans(model_nb, pairwise ~ trial)
emmeans(model_nb,
        pairwise ~ trial,
        at = list(trial = c("T2S1","T2S2")))

#Model with door as a covariate
model_nb2 <- glmmTMB(visits ~ trial * season + door + (1 | ID), data = df_long, family = nbinom2)
summary (model_nb2)


pred <- ggpredict(model_nb, terms = c("trial", "season"))

ggplot(pred, aes(x = x, y = predicted, color = group, group = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  theme_classic() +
  labs(
    x = "Trial",
    y = "Predicted number of visits",
    color = "Season",
    title = "Predicted visitation across trials and seasons")
