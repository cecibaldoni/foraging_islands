library(tidyverse)
library(here)

result <- read.csv(here("csv/foraging_results.csv"))

foraging_master <- result %>%
  distinct(unique_trial_ID) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  separate(unique_trial_ID, into = c("season", "ID", "trial"),sep = "_", remove = FALSE)

#First frame and time to baited island interaction and join
first_success <- result %>%
  filter(island_debug %in% c("A", "D")) %>%
  group_by(unique_trial_ID) %>%
  summarise(first_AD_frame = as.numeric (min(frame, na.rm = TRUE)),.groups = "drop")

foraging_master <- foraging_master %>%
  left_join(first_success, by = "unique_trial_ID")

foraging_master <- foraging_master %>%
  left_join(first_success, by = "unique_trial_ID") %>%
  mutate(first_AD_time  = round(first_AD_frame * (1/24),2))

# #This does not work for me?
# Error in `mutate()`:
#   ℹ In argument: `first_AD_time = round(first_AD_frame * (1/24), 2)`.
# Caused by error:
#   ! object 'first_AD_frame' not found
# Run `rlang::last_trace()` to see where the error occurred.

#Islands visited
letter_food_counts <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  filter(!is.na(island_debug), !is.na(food)) %>%
  filter(food %in% c(0, 1)) %>%
  group_by(unique_trial_ID, island_debug, food) %>%
  summarise(count = n(), .groups = "drop") %>%
  unite(letter_food, island_debug, food, sep = "_") %>%
  pivot_wider(names_from = letter_food, values_from = count, values_fill = 0)

letter_food_counts <- letter_food_counts %>%
  select(-B_1)  # Removes the column named B_1 because there should not be any (typo)

foraging_master <- foraging_master %>%
  left_join(letter_food_counts, by = "unique_trial_ID")

#Time management
#Travel time
travelling_counts <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>% 
  group_by(unique_trial_ID) %>%
  summarise(travelling_frames = sum(journey == "travelling"), .groups = "drop")

foraging_master <- foraging_master %>%
  left_join(travelling_counts, by = "unique_trial_ID")

foraging_master <- foraging_master %>%
  mutate(travelling_time = round(travelling_frames * (1/24), 2))

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
  mutate(islands_time = round(island_frame * (1/24), 2))

#Time not moving
foraging_master <- foraging_master %>% 
  mutate(nonmoving_frames = 43200 - last_frame)

foraging_master <- foraging_master %>% 
  mutate (nonmoving_time = round(nonmoving_frames * (1/24),2))


#Ordering the columns
foraging_master <- foraging_master %>%
  select(unique_trial_ID, season, ID, trial, first_AD_frame, first_AD_time, A_0, A_1, B_0, C_0, D_0, D_1,travelling_frames,
         travelling_time,island_frame, islands_time, nonmoving_frames, nonmoving_time, last_frame)

#Save csv
write.csv( foraging_master, here("csv", "foraging_master.csv"),row.names = FALSE)


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

