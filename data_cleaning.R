library(tidyverse)
library(here)
library(ggdist)

result <- read.csv(here("csv/foraging_results.csv"))

foraging_master <- result %>%
  distinct(unique_trial_ID) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  separate(unique_trial_ID, into = c("season", "ID", "trial"),sep = "_", remove = FALSE)

#First frame and time to baited island interaction and join
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
  mutate (nonmoving_time = round(nonmoving_frames * (1/30),2))


#Ordering the columns
foraging_master <- foraging_master %>%
  select(unique_trial_ID, season, ID, trial, first_AD_frame, first_AD_time, first_AD_island, first_AD_door, A_0, A_1, B_0, C_0, D_0, D_1,travelling_frames,
         travelling_time,island_frame, islands_time, nonmoving_frames, nonmoving_time, last_frame)

#Save csv
write.csv( foraging_master, here("csv", "foraging_master.csv"),row.names = FALSE)

#Island and door interactions
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

#Plot the position interaction with the trial
# pivot A_1, A_2, etc. into long format
A_long <- position_counts %>%
  select(trial, season, starts_with("A_")) %>%
  pivot_longer(
    cols = starts_with("A_"),
    names_to = "A_position",
    values_to = "count")

ggplot(A_long, aes(x = A_position, y = trial, color = season, size = count)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("spring" = "green", "summer" = "orange", "winter" = "purple")) +
  labs(
    x = "Position",
    y = "Trial",
    color = "Season",
    size = "Count"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Plots
counts <- foraging_master %>%
  group_by(season) %>%
  summarise(n = sum(!is.na(first_AD_time))) 
#Plot the time to reach the first successful island per season
ggplot(foraging_master, aes(x = season, y = first_AD_time)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  geom_text(data = counts, aes(x = season, y = max(foraging_master$first_AD_time, na.rm = TRUE) + 1, 
                              label = paste0("n=", n)), inherit.aes = FALSE, size = 5) + 
  labs(x = "Season", y = "Time to first baited island (seconds)",
  title = "Distribution of First AD Time by Season") +
  theme_minimal(base_size = 14)


foraging_master <- foraging_master %>%
  mutate(trial = factor(trial, levels = c("T1S1", "T1S2", "T2S1", "T2S2")),
    season_trial = factor(paste(season, trial, sep = "_"),
    levels = c("spring_T1S1", "spring_T1S2", "spring_T2S1", "spring_T2S2",
               "summer_T1S1", "summer_T1S2", "summer_T2S1", "summer_T2S2",
               "winter_T1S1", "winter_T1S2", "winter_T2S1", "winter_T2S2")))
#Plot the time to reach the first successful island per trial per season
ggplot(foraging_master, aes(x = season_trial, y = first_AD_time)) +
  geom_boxplot(fill = "skyblue", color = "black") +
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
  theme_classic()