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
library(dplyr)


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
df_long$season <- factor(df_long$season, levels = c("summer", "winter", "spring"))

model_nb <- glmmTMB(visits ~ trial * season + (1 | ID), data = df_long, family = nbinom2)
summary(model_nb)

#check from here on
library(car)
Anova(model_nb)
library(emmeans)
emmeans(model_nb, pairwise ~ trial)
emmeans(model_nb,
        pairwise ~ trial,
        at = list(trial = c("T2S1","T2S2")))

model_nb2 <- glmmTMB(
  visits ~ trial * season + door + (1 | ID),
  data = df_long,
  family = nbinom2
)
library(DHARMa)

sim <- simulateResiduals(model_nb2)

plot(sim)
---
pred <- ggpredict(model_nb2,
                  terms = c("trial", "season"))
ggplot(pred,
       aes(x = x,
           y = predicted,
           color = group,
           group = group)) +
  
  geom_line(size = 1.2) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.1) +
  
  theme_bw() +
  
  labs(
    x = "Trial",
    y = "Predicted visits",
    color = "Season"
  )
---
pred <- ggpredict(
  model_nb,
  terms = c("trial", "season")
)

ggplot(pred,
       aes(x = x,
           y = predicted,
           color = group,
           group = group)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = conf.low,
        ymax = conf.high),
    width = 0.1
  ) +
  
  theme_bw() +
  
  labs(
    x = "Trial",
    y = "Predicted number of visits",
    color = "Season",
    title = "Predicted visitation across trials and seasons"
  )


# For each individual, compare T1S1 → T1S2 (food depletes within day 1)
# and T2S1 → T2S2 (food depletes within day 2)

# First create a binary: was door visited (1) or not (0)?
df_binary <- df_long %>%
  mutate(visited = as.integer(visits > 0))

# Pivot so each row = ID × door, with columns for each trial
df_wide_doors <- df_binary %>%
  select(ID, season, door, trial, visited) %>%
  pivot_wider(names_from = trial,
              values_from = visited)

# Create the key predictor: was this door visited in the PREVIOUS trial?
# Pair 1: T1S1 (full) → T1S2 (empty same day)
pair1 <- df_wide_doors %>%
  select(ID, season, door, T1S1, T1S2) %>%
  rename(visited_prev = T1S1,
         visited_next = T1S2) %>%
  mutate(pair = "day1_S1toS2",
         food_available = 1)   # food was present in previous trial

# Pair 2: T2S1 (full) → T2S2 (empty same day)  
pair2 <- df_wide_doors %>%
  select(ID, season, door, T2S1, T2S2) %>%
  rename(visited_prev = T2S1,
         visited_next = T2S2) %>%
  mutate(pair = "day2_S1toS2",
         food_available = 1)

# Pair 3: T1S1 → T2S1 (same food state, different days = learning/memory)
pair3 <- df_wide_doors %>%
  select(ID, season, door, T1S1, T2S1) %>%
  rename(visited_prev = T1S1,
         visited_next = T2S1) %>%
  mutate(pair = "S1_day1today2",
         food_available = 1)

df_pairs <- bind_rows(pair1, pair2, pair3) %>%
  drop_na(visited_prev, visited_next) %>%
  mutate(
    ID     = factor(ID),
    season = factor(season),
    door   = factor(door),
    pair   = factor(pair),
    island = str_sub(door, 1, 1)
  )


df_pairs$season <- factor(df_pairs$season, levels = c("summer", "winter", "spring"))
m1 <- brm(
  visited_next ~ visited_prev * season + (1 | ID) + (1 | door),
  data   = df_pairs,
  family = bernoulli(link = "logit"),
  prior  = c(
    prior(normal(0, 1.5), class = b),
    prior(normal(0, 1.5), class = Intercept),
    prior(exponential(1), class = sd)
  ),
  chains = 4,
  iter   = 4000,
  warmup = 1000,
  cores  = 4,
  seed   = 42
)

summary(m1)
pp_check(m1, ndraws = 4000)

m2 <- brm(
  visited_next ~ visited_prev * season * pair + (1 | ID) + (1 | door),
  data   = df_pairs,
  family = bernoulli(link = "logit"),
  prior  = c(
    prior(normal(0, 1.5), class = b),
    prior(normal(0, 1.5), class = Intercept),
    prior(exponential(1), class = sd)
  ),
  chains = 4,
  iter   = 4000,
  warmup = 1000,
  cores  = 4,
  seed   = 42
)
summary(m2)
pp_check(m2, ndraws = 4000)

# Collapse to island level
df_island <- df_pairs %>%
  group_by(ID, season, pair, island) %>%
  summarise(
    visited_prev_isl = as.integer(sum(visited_prev, na.rm=T) > 0),
    visited_next_isl = as.integer(sum(visited_next, na.rm=T) > 0),
    .groups = "drop"
  )

m3 <- brm(
  visited_next_isl ~ visited_prev_isl * season + (1 | ID) + (1 | island),
  data   = df_island,
  family = bernoulli(link = "logit"),
  prior  = c(
    prior(normal(0, 1.5), class = b),
    prior(normal(0, 1.5), class = Intercept),
    prior(exponential(1), class = sd)
  ),
  chains = 4, iter = 4000, warmup = 1000, cores = 4, seed = 42
)

summary(m3)
pp_check(m3, ndraws = 4000)

--- #Plots 

# Get predicted probability for visited_prev=1 and visited_prev=0
# then compute the difference (the learning signal)

newdata <- expand.grid(
  visited_prev = c(0, 1),
  season       = c("spring", "summer", "winter")
)

draws <- m1 %>%
  epred_draws(newdata   = newdata,
              re_formula = NA)

# Compute contrast: P(revisit | visited before) - P(revisit | not visited before)
contrast <- draws %>%
  ungroup() %>%
  select(season, visited_prev, .draw, .epred) %>%
  pivot_wider(names_from  = visited_prev,
              values_from = .epred,
              names_prefix = "prev_") %>%
  mutate(
    contrast   = prev_1 - prev_0,
    season     = factor(season,
                        levels = c("winter","summer","spring"))
  )

ggplot(contrast,
       aes(x = contrast, y = season, fill = season)) +
  
  stat_halfeye(
    .width       = c(0.66, 0.95),
    point_interval = median_qi,
    alpha        = 0.8,
    normalize    = "xy"
  ) +
  
  geom_vline(xintercept = 0,
             linetype   = "dashed",
             color      = "grey40") +
  
  scale_fill_manual(
    values = c(winter = "#5B8DB8",
               summer = "#E8A838",
               spring = "#6BAF6B")
  ) +
  
  labs(
    x     = "P(revisit | visited before) − P(revisit | not visited before)",
    y     = NULL,
    title = "Does prior visitation predict revisitation?",
    subtitle = "Positive = animal returns to previously visited doors\nVertical line = no effect"
  ) +
  
  theme_bw() +
  theme(legend.position = "none",
        plot.subtitle    = element_text(color = "grey40"))

draws_summary <- draws %>%
  group_by(season, visited_prev) %>%
  median_qi(.epred, .width = c(0.66, 0.95)) %>%
  mutate(
    visited_prev = factor(visited_prev,
                          levels = c(0, 1),
                          labels = c("not visited\nbefore",
                                     "visited\nbefore")),
    season = factor(season,
                    levels = c("winter","summer","spring"))
  )

ggplot(draws_summary,
       aes(x = visited_prev,
           y = .epred,
           color = season,
           group = season)) +
  
  geom_line(linewidth = 1.2, alpha = 0.8) +
  
  geom_point(size = 4) +
  
  geom_errorbar(
    aes(ymin = .lower, ymax = .upper),
    width     = 0.08,
    linewidth = 0.8
  ) +
  
  facet_wrap(~ season) +
  
  scale_color_manual(
    values = c(winter = "#5B8DB8",
               summer = "#E8A438",
               spring = "#6BAF6B")
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  
  labs(
    x     = NULL,
    y     = "Predicted probability of visiting door",
    title = "Revisitation probability by season",
    subtitle = "Thick bars = 66% CI, thin bars = 95% CI"
  ) +
  
  theme_bw() +
  theme(legend.position = "none")

# Get individual-level predictions (include random effects)
newdata_ind <- df_pairs %>%
  distinct(ID, season) %>%
  mutate(visited_prev = 1)   # probability of revisiting a door visited before

draws_ind <- m1 %>%
  epred_draws(newdata    = newdata_ind,
              re_formula = ~ (1 | ID),
              allow_new_levels = FALSE)

draws_ind_sum <- draws_ind %>%
  group_by(ID, season) %>%
  median_qi(.epred, .width = 0.89) %>%
  mutate(season = factor(season,
                         levels = c("winter","summer","spring")))

ggplot(draws_ind_sum,
       aes(x     = .epred,
           y     = reorder(ID, .epred),
           color = season)) +
  
  geom_pointrange(
    aes(xmin = .lower, xmax = .upper),
    linewidth = 0.6
  ) +
  
  geom_vline(xintercept = 0.5,
             linetype   = "dashed",
             color      = "grey50") +
  
  facet_wrap(~ season,
             scales = "free_y",
             ncol   = 3) +
  
  scale_color_manual(
    values = c(winter = "#5B8DB8",
               summer = "#E8A438",
               spring = "#6BAF6B")
  ) +
  
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::percent) +
  
  labs(
    x        = "P(revisit a previously visited door)",
    y        = NULL,
    title    = "Individual revisitation probabilities",
    subtitle = "89% credible intervals · dashed line = chance"
  ) +
  
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y      = element_text(size = 8))
