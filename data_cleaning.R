library(tidyverse)
library(here)
library(ggdist)
library(trajr)

result <- read.csv(here("csv/foraging_results.csv"))

foraging_master <- result %>%
  distinct(unique_trial_ID) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  separate(unique_trial_ID, into = c("season", "ID", "trial"),sep = "_", remove = FALSE)

# ----- First frame and time to baited island interaction and join-----
first_success <- result %>%
  filter(island_debug %in% c("A", "D")) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>% 
  mutate(frame = as.numeric(frame)) %>%
  group_by(unique_trial_ID) %>%
  slice_min(frame, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(unique_trial_ID, first_AD_frame = frame, first_AD_island = island_debug, first_AD_door = place)

foraging_master <- foraging_master %>%
  left_join(first_success, by = "unique_trial_ID") %>%
  mutate(first_AD_time = round(first_AD_frame / 30, 2))

#Islands visited
letter_food_counts <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  filter(!is.na(island_debug), !is.na(food)) %>%
  filter(food %in% c(0, 1)) %>%
  group_by(unique_trial_ID, island_debug, food) %>%
  summarise(count = n(), .groups = "drop") %>%
  unite(letter_food, island_debug, food, sep = "_") %>%
  pivot_wider(names_from = letter_food, values_from = count, values_fill = 0)

foraging_master <- foraging_master %>%
  left_join(letter_food_counts, by = "unique_trial_ID")

#Time management
#Travel time
travelling_counts <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>% 
  group_by(unique_trial_ID) %>%
  summarise(travelling_frames = sum(journey == "travelling"), .groups = "drop")

foraging_master <- foraging_master %>%
  left_join(travelling_counts, by = "unique_trial_ID") %>%
  mutate(travelling_time = round(travelling_frames / 30, 2))
#Island time
result <- result %>%
  mutate(frame = as.numeric(frame))

island_frame <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%  
  group_by(unique_trial_ID) %>%
  summarise(last_frame = max(frame, na.rm = TRUE),.groups = "drop")

foraging_master <- foraging_master %>%
  left_join(island_frame, by = "unique_trial_ID") 

foraging_master <- foraging_master %>%  
  mutate(island_frame = last_frame - travelling_frames)

foraging_master <- foraging_master %>%
  mutate(islands_time = round(island_frame * (1/30), 2))
#Time not moving
foraging_master <- foraging_master %>% 
  mutate(nonmoving_frames = 43200 - last_frame)

foraging_master <- foraging_master %>% 
  mutate (nonmoving_time = round(nonmoving_frames * (1/30), 2))
#Time moving
foraging_master <- foraging_master %>% 
  mutate(moving_time = round(last_frame * 30, 2))

#Ordering the columns
foraging_master <- foraging_master %>%
  select(unique_trial_ID, season, ID, trial, first_AD_time, first_AD_island, first_AD_door, A_0, A_1, B_0, C_0, D_0, D_1,
         travelling_time, islands_time, nonmoving_time, moving_time, last_frame)

# ----- Island and door interactions (new df) -----
position_counts <- result %>%
  distinct(unique_trial_ID) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  separate(unique_trial_ID, into = c("season", "ID", "trial"),sep = "_", remove = FALSE)

place_counts <- result %>%
  filter(!is.na(island_debug), !is.na(place)) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  count(unique_trial_ID, island_debug, place, name = "n") %>%
  mutate(col_name = paste0(island_debug, "_", place)) %>%
  select(unique_trial_ID, col_name, n) %>%
  pivot_wider(names_from = col_name, values_from = n, values_fill = 0)

position_counts <- position_counts %>%
  left_join(place_counts, by = "unique_trial_ID")

#Save csv
write.csv( position_counts, here("csv", "position_counts.csv"),row.names = FALSE)

# ----- Loop trial -----
df_trajectory <- result %>%
  select(unique_trial_ID, season, ID, trial, time, x, y, journey) %>%
  mutate(time = as.numeric(time)) %>%
  filter(!is.na(x), !is.na(y), !is.na(time)) %>%
  arrange(unique_trial_ID, time)

#List by unique_trial_ID
df_list <- df_trajectory %>%
  group_split(unique_trial_ID) %>%
  setNames(unique(df_trajectory$unique_trial_ID))

#Trajectory function
safe_traj <- function(df) {
  df <- df %>% filter(!is.na(x), !is.na(y), !is.na(time))
  if (nrow(df) < 2) return(NULL)
  TrajFromCoords(df, xCol = "x", yCol = "y", fps = 30, spatialUnits = "cm")}

#Metrics calculation
compute_metrics <- function(df_trial) {
  trial_id <- df_trial$unique_trial_ID[1]
  season   <- df_trial$season[1]
  ID       <- df_trial$ID[1]
  trial    <- df_trial$trial[1]
#Total trajectory
  traj_all <- safe_traj(df_trial)
  if (is.null(traj_all)) return(NULL)
  distance_tot      <- TrajLength(traj_all)
  straightness_tot  <- TrajStraightness(traj_all)
#Travelling trajectory
  df_travel <- df_trial %>%
    filter(journey == "travelling")
  traj_travel <- safe_traj(df_travel)
  straightness_travel <- ifelse(
    is.null(traj_travel), NA, TrajStraightness(traj_travel))
#Island trajectory
  df_island <- df_trial %>%
    filter(grepl("^at_", journey))
  traj_island <- safe_traj(df_island)
  straightness_island <- ifelse(
    is.null(traj_island), NA, TrajStraightness(traj_island))
#Output
  data.frame(
    unique_trial_ID       = trial_id,
    season                = season,
    ID                    = ID,
    trial                 = trial,
    distance_total_cm     = distance_tot,
    straightness_total    = straightness_tot,
    straightness_travel   = straightness_travel,
    straightness_island   = straightness_island)}

#Run for all the trials
results_df <- df_list %>%
  lapply(compute_metrics) %>%
  bind_rows()

results_df <- results_df %>%
  mutate(distance_total_cm   = round(distance_total_cm, 2),
         straightness_total  = round(straightness_total, 5),
         straightness_travel = round(straightness_travel, 5),
         straightness_island = round(straightness_island, 5))

#Join to the master
metrics_to_join <- results_df %>%
  select(unique_trial_ID, distance_total_cm, straightness_total, straightness_travel, straightness_island)
foraging_master <- foraging_master %>%
  left_join(metrics_to_join, by = "unique_trial_ID")

#Save csv
write.csv( foraging_master, here("csv", "foraging_master.csv"),row.names = FALSE)

# ----- Plots -------
#Plot the island door interactions  
all_letters_long <- position_counts %>%
  pivot_longer(cols = matches("^[A-D]_[0-6]$"),
    names_to = "letter_position",
    values_to = "count")
ggplot(all_letters_long, aes(x = letter_position, y = count, fill = season)) +
  geom_col(position = "dodge") +
  labs(x = "Island door interaction", y = "Count", fill = "Season") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot of the moving time
ggplot(foraging_master, aes(x = trial, y = moving_time, color = season))+
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1, alpha = 0.4) +
  facet_wrap (~ season)+
  theme_classic()

#Plot the time to reach the first successful island per season
counts <- foraging_master %>%
  group_by(season) %>%
  summarise(n = sum(!is.na(first_AD_time))) 

ggplot(foraging_master, aes(x = season, y = first_AD_time)) +
  geom_violin(fill = "skyblue", color = "black", alpha = 0.3, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, color = "darkblue") +
  labs(x = "Season", y = "Time to first baited island (seconds)", title = "Raincloud plot of First AD Time by Season") +
  theme_minimal(base_size = 14)


#Plot the time to reach the first successful island per trial per season
foraging_master <- foraging_master %>%
  mutate(trial = factor(trial, levels = c("T1S1", "T1S2", "T2S1", "T2S2")),
    season_trial = factor(paste(season, trial, sep = "_"),
    levels = c("spring_T1S1", "spring_T1S2", "spring_T2S1", "spring_T2S2",
               "summer_T1S1", "summer_T1S2", "summer_T2S1", "summer_T2S2",
               "winter_T1S1", "winter_T1S2", "winter_T2S1", "winter_T2S2")))
ggplot(foraging_master, aes(x = season_trial, y = first_AD_time)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, color = "darkblue") +
  labs(x = "Season and Trial", y = "Time to first baited island (seconds)", title = "First AD Time per Season and Trial") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Plot the first island A/D and which door is opened
ggplot(foraging_master %>% filter(!is.na(first_AD_island)),
  aes(x = first_AD_island, y = first_AD_time, fill = factor(first_AD_door))) +
  # half violin
  #geom_violin(position = position_dodge(width = 0.8), alpha = 0.5, trim = FALSE) +
  # boxplot
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  # jittered points ("rain")
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1, alpha = 0.4) +
  labs(x = "First A/D Island", y = "Time to First A/D (s)", fill = "Door") +
  facet_wrap (~ season)+
  theme_classic()

#Plot for the distance covered
ggplot(foraging_master, aes(x = season_trial, y = distance_total_cm)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, color = "darkblue") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
