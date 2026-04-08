library(tidyverse)
library(here)
library(ggdist)
library(trajr)
library(purrr)
library(SimilarityMeasures)

result <- read.csv(here("csv/processed/foraging_results.csv"))

# Island and door interactions df -----
# This dataframe will produce a table with IDs, number of interaction for each island door, and success rate.
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

position_counts <- position_counts %>%
  mutate(unique_ID = sub("_[^_]+$", "", unique_trial_ID))

#Loop for trials 

expected_trials <- c("T1S1", "T1S2", "T2S1", "T2S2")

unique_ids <- unique(position_counts$unique_ID)

trial_ls <- lapply(unique_ids, function(uid) {
df_id <- position_counts %>% 
  filter (unique_ID == uid) %>% 
  arrange (trial)

message('Processing: ', uid)

found_trials <- df_id$trial
n_trials <- length(found_trials)

if (n_trials != 4) {missing <- setdiff(expected_trials, found_trials)
  message("   -> Warning: ", n_trials, " trial(s) found. Missing trial(s): ", paste(missing, collapse = ", "))}

return(df_id)})

names(trial_ls) <- unique_ids

#Success rate
#Find AD
doors_AD <- function(df) {
  grep("^(A|D)_", names(df), value = TRUE)
}
#Process each id
trial_ls_processed <- map(trial_ls, function(df) {
  #df = trial_ls[[17]]
   ad_cols <- doors_AD(df)
#~ .x & !lag(.x) This means: count it if it is used now and if it was NOT used in the previous session 
  df <- df %>%
    arrange(trial) %>%
    mutate(across(all_of(ad_cols), ~replace_na(.x, 0))) %>% 
    mutate(across(all_of(ad_cols), ~ .x != 0)) %>%
    mutate(across(all_of(ad_cols), ~ case_when(
          trial %in% c("T1S1", "T2S1") ~ .x,
          trial %in% c("T1S2", "T2S2") ~ .x & !lag(.x),
          TRUE ~ FALSE),.names = "count_{.col}")) %>%
    rowwise() %>%
    mutate(door_count = sum(c_across(starts_with("count_")), na.rm = TRUE)) %>%
    ungroup()
   
  df <- df  %>%
    arrange(trial) %>%
    mutate(success_rate = case_when(
        trial %in% c("T1S1", "T2S1") ~ door_count / 12,
        trial %in% c("T1S2", "T2S2") &
          (12 - lag(door_count)) > 0 ~
          door_count / (12 - lag(door_count)),TRUE ~ NA_real_))
})

#Join the success rate to position count
doors_summary_df <- bind_rows(trial_ls_processed)
doors_summary_df <- doors_summary_df %>%
  select(unique_trial_ID, door_count, success_rate)
position_counts <- position_counts %>%
  select(-any_of(c("door_count", "success_rate"))) %>%
  left_join(doors_summary_df, by = "unique_trial_ID") %>% 
  mutate (door_count = replace_na(door_count, 0), success_rate = replace_na(success_rate, 0))

#Save csv
write.csv( position_counts, here("csv/processed", "interactions_counts.csv"),row.names = FALSE)

#Plot success rate per trial
ggplot(position_counts, aes(x =trial, y = success_rate)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, color = "darkblue") +
  #geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", linewidth = 1) +
  facet_wrap (~ season)+
  labs(x = "Season and Trial", y = "Success rate") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

a <- ggplot(position_counts, aes(x = success_rate)) +
  geom_density(fill = "blue", alpha = 0.4)
b <- ggplot(position_counts, aes(x = success_rate)) +
  geom_density(fill = "blue", alpha = 0.4) +
  facet_wrap(~season)
a
b

# Master df -----
foraging_master <- result %>%
  distinct(unique_trial_ID) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  separate(unique_trial_ID, into = c("season", "ID", "trial"),sep = "_", remove = FALSE) %>% 
  mutate(season_ID = paste(season, ID, sep = "_"))

# --- First frame and time of island interaction and join ---
first_success <- result %>%
  filter(island_debug %in% c("A", 'B', 'C', "D")) %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  mutate(frame = as.numeric(frame)) %>%
  group_by(unique_trial_ID) %>%
  slice_min(frame, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(unique_trial_ID, first_island_frame = frame, first_island = island_debug, first_door = place)

foraging_master <- foraging_master %>%
  left_join(first_success, by = "unique_trial_ID") %>%
  mutate(first_island_time = round(first_island_frame / 30, 2))

## ---- Number of visits per island ----
letter_food_counts <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>%
  filter(!is.na(island_debug)) %>%
  group_by(unique_trial_ID, island_debug) %>%
  summarise(count = n(), .groups = "drop") %>%
  unite(letter_food, island_debug, sep = "_") %>%
  pivot_wider(names_from = letter_food, values_from = count, values_fill = 0)

foraging_master <- foraging_master %>%
  left_join(letter_food_counts, by = "unique_trial_ID") %>% 
  mutate (tot_visit = rowSums(across(c(A, B, C, D)), na.rm = TRUE))

## ----- Time management -----
#Travel time
travelling_counts <- result %>%
  filter(unique_trial_ID != "unique_trial_ID") %>% 
  group_by(unique_trial_ID) %>%
  summarise(travelling_frames = sum(journey == "travelling"), .groups = "drop")

foraging_master <- foraging_master %>%
  left_join(travelling_counts, by = "unique_trial_ID") %>%
  mutate(travelling_time = round(travelling_frames* (1/ 30), 2))
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
  mutate(nonmoving_frames = (54000 - last_frame))

foraging_master <- foraging_master %>% 
  mutate (nonmoving_time = pmax(round(nonmoving_frames * (1/30), 2), 0))
#Time moving
foraging_master <- foraging_master %>% 
  mutate(moving_time = round(last_frame * (1/30), 2))


##--- Distance covered in cm ----
# Clean trajectory data (time, x and y)
df_trajectory <- result %>%
  select(unique_trial_ID, season, ID, trial, time, x, y) %>%
  filter(!is.na(x), !is.na(y), !is.na(time)) %>%
  arrange(unique_trial_ID, time)

# Function to compute total distance
compute_distance <- function(df) {
  if(nrow(df) < 2) return(NULL)
  traj <- TrajFromCoords(df[,c("x","y")],
                         fps = 30,
                         spatialUnits = "cm")
  data.frame(
    unique_trial_ID = df$unique_trial_ID[1],
    season          = df$season[1],
    ID              = df$ID[1],
    trial           = df$trial[1],
    distance_total_cm = TrajLength(traj)
  )
}

# Apply to each trial
results_df <- df_trajectory %>%
  group_split(unique_trial_ID) %>%
  map_dfr(compute_distance) %>%
  mutate(distance_total_cm = round(distance_total_cm, 2))

#Join to the master
metrics_to_join <- results_df %>%
  select(unique_trial_ID, distance_total_cm)
foraging_master <- foraging_master %>%
  left_join(metrics_to_join, by = "unique_trial_ID")
#normalisation total_distance_cm and plot
foraging_master <- foraging_master %>%
  mutate(distance_rate = round(distance_total_cm / moving_time, 2))

#Ordering the columns
foraging_master <- foraging_master %>%
  select(unique_trial_ID, season, ID, trial, season_ID, first_island_time, first_island, A, B, C, D, tot_visit,
         travelling_time, islands_time, nonmoving_time, moving_time, distance_total_cm, distance_rate)

#Save csv
write.csv( foraging_master, here("csv/processed", "foraging_master.csv"),row.names = FALSE)

# Similarities in the areas covered ----
grid_size <- 30
# Mh. A 30×30 grid = 900 cells
#This means that  most will be zero for any given trial
# Correlating two very sparse vectors inflates similarity
# if they avoid 850 cells, they look very similar just by sharing zeros
# 
# result_grid <- result %>%
#   mutate(
#     x_bin = cut(x, breaks = grid_size, labels = FALSE),
#     y_bin = cut(y, breaks = grid_size, labels = FALSE))
# 
# occupancy <- result_grid %>% 
#   mutate(season_ID = paste(season, ID, sep = "_")) %>% 
#   count(season_ID, trial, x_bin, y_bin) 
#   
# maps <- occupancy %>%
#   pivot_wider(names_from = c(x_bin, y_bin), values_from = n, values_fill = 0)
# 
# maps_norm <- maps %>%
#   rowwise() %>%
#   # mutate(across(-c(season_ID, trial), ~ . / sum(c_across(-c(season_ID, trial)))))
#   # This is incorrect
#   mutate(row_total = sum(c_across(-c(season_ID, trial)))) %>%
#   mutate(across(-c(season_ID, trial, row_total), ~ . / row_total)) %>%
#   select(-row_total)
# 
# # Comparison among T1-T2 T2-T3 T3-T4
# # Function for the correlation of two vectors
# similarity_fun <- function(a,b){
#   cor(a,b)
# }

grid_size <- 10  # reduced from 30 (so 10x10 = 100 cells, much less sparse)

result_grid <- result %>%
  filter(!is.na(x), !is.na(y)) %>%          # guard against NAs in coordinates
  mutate(
    x_bin = cut(x, breaks = grid_size, labels = FALSE),
    y_bin = cut(y, breaks = grid_size, labels = FALSE)
  )

occupancy <- result_grid %>%
  mutate(season_ID = paste(season, ID, sep = "_")) %>%
  count(season_ID, trial, x_bin, y_bin)

maps <- occupancy %>%
  pivot_wider(names_from = c(x_bin, y_bin), values_from = n, values_fill = 0)

# compute row total first, then divide
maps_norm <- maps %>%
  rowwise() %>%
  mutate(row_total = sum(c_across(-c(season_ID, trial)))) %>%
  mutate(across(-c(season_ID, trial, row_total), ~ . / row_total)) %>%
  select(-row_total) %>%
  ungroup()

# Bhattacharyya coefficient instead of Pearson correlation
similarity_fun <- function(a, b) {
  sum(sqrt(a * b))  # BC: ranges 0 (no overlap) to 1 (identical)
}

# The same three comparisons 
similarities_1to4 <- maps_norm %>%
  group_by(season_ID) %>%
  arrange(trial, .by_group = TRUE) %>% #need to turn maps_norm into matrix
  group_modify(~{
    mat <- as.matrix(.x[,-c(1,2)]) #remove first 2 columns
    sims <- sapply(1:(nrow(mat)-1), function(i){ 
      similarity_fun(mat[i, ], mat[i+1, ]) #compare the 3 pairs of trials
    })
    tibble(mean_similarity = mean(sims)) #calculate mean similarity
  })

#se ho capito bene, questo code compara gni trial con il suo vicino
#quidi T1S1 con T2S1 con T1S2 con T2S2. T1S1 and T2S2 are never directly compared to each other.
#not necessarily wrong, but just need to be sure that that is what you want

# Comparison of first trial (T1S1) vs last trial (T2S2) only
similarities_firstlast <- maps_norm %>%
  group_by(season_ID) %>%
  arrange(trial, .by_group = TRUE) %>%
  group_modify(~{
    mat <- as.matrix(.x[,-c(1,2)])
    if(nrow(mat) >= 4){
      sim <- similarity_fun(mat[1,], mat[4,])
    } else {
      sim <- NA_real_
    }
    tibble(mean_similarity = sim)
  })

#Comparison of 1-2 and 3-4
similarities_2days <- maps_norm %>%
  group_by(season_ID) %>%
  arrange(trial, .by_group = TRUE) %>%
  group_modify(~{
    mat <- as.matrix(.x[,-c(1,2)])
    idx <- seq(1, nrow(mat)-1, by = 2) #only 2 pair of trials 
    sims <- sapply(idx, function(i){
      similarity_fun(mat[i, ], mat[i+1, ])
    })
    tibble(mean_similarity = mean(sims, na.rm = TRUE))
  })

#Comparison of 1-3 and 2-4
similarities_noncons <- maps_norm %>%
  group_by(season_ID) %>%
  arrange(trial, .by_group = TRUE) %>%
  group_modify(~{
    mat <- as.matrix(.x[,-c(1,2)])
    sims <- c()
    if(nrow(mat) >= 3){
      sims <- c(sims, similarity_fun(mat[1,], mat[3,]))
    }
    if(nrow(mat) >= 4){
      sims <- c(sims, similarity_fun(mat[2,], mat[4,]))
    }
    tibble(mean_similarity = mean(sims, na.rm = TRUE))
  })

#Rename the columns 
similarities_1to4<- similarities_1to4 %>% rename(sim_1to4 = mean_similarity)
similarities_2days <- similarities_2days %>% rename(sim_2days = mean_similarity)
similarities_noncons <- similarities_noncons %>% rename(sim_noncons = mean_similarity)
#Merge
foraging_similarities <- similarities_1to4 %>%
  left_join(similarities_2days, by = "season_ID") %>%
  left_join(similarities_noncons, by = "season_ID") %>% 
  separate(season_ID, into = c("season", "ID"),sep = "_", remove = FALSE)
#Save csv
write.csv(foraging_similarities, here("csv/processed", "foraging_similarities.csv"),row.names = FALSE)

# ----- Plots -------

#Plot the island interactions count
#Here we can switch trial and season 
all_letters_long <- foraging_master %>%
  pivot_longer(cols = c(A, B, C, D),
    names_to = "letter",
    values_to = "count")
ggplot(all_letters_long, aes(x = letter, y = count, fill = season)) +
  geom_col(position = "dodge") +
  facet_wrap (~ trial) +
  labs(x = "Island door interaction", y = "Count", fill = "Season") +
  theme_classic() 

#Plot the distance rate
ggplot(foraging_master, aes(x = season, y = distance_rate, color = season)) +
  #geom_violin(fill = "skyblue", color = "black") +
  geom_boxplot(fill = "skyblue", color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap (~ trial)

#Plot of the moving time
ggplot(foraging_master, aes(x = trial, y = moving_time, color = season))+
  geom_violin(width = 0.15, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1, alpha = 0.4) +
  facet_wrap (~ season)+
  theme_classic()
                                                                                                                                                                                   
#Plot for the distance covered
ggplot(foraging_master, aes(x = trial, y = distance_total_cm, color = season)) +
  geom_violin(fill = "skyblue", color = "black") +
  #geom_boxplot(fill = "skyblue", color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap (~ season) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot the time to the first island
ggplot(foraging_master, aes(x = season, y = first_island_time, color = trial)) +
  geom_violin(fill = "skyblue", color = "black", alpha = 0.3, trim = FALSE) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(x = "Door", y = "Time to first island (s)") +
  facet_wrap (~ trial) +
  theme_minimal(base_size = 14)

# Summarize occupancy per square per season_ID
heatmap_data <- occupancy %>%
  mutate(season_ID = as.character(season_ID), season = sub("_.*$", "", season_ID)) %>%
  group_by(season_ID, season, x_bin, y_bin) %>%
  summarise(count = sum(n), .groups = "drop") %>%
  mutate(x_bin = as.numeric(x_bin), y_bin = as.numeric(y_bin))
# Split data by season_ID
heatmap_list <- split(heatmap_data, heatmap_data$season_ID)
#plot for each season_ID summed trials
plots <- lapply(names(heatmap_list), function(season) {
  ggplot(heatmap_list[[season]], aes(x = x_bin, y = y_bin, fill = count)) +
    geom_tile(color = "grey80") +
    scale_fill_gradientn(colours = c("#ffffcc", "#ffeda0", "#feb24c", "#fd8d3c", "#f03b20"))+
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = paste("Animal Movement Heatmap -", season),
      x = "X Grid",
      y = "Y Grid",
      fill = "Pass Count"
    ) +
    theme_minimal()
})
plots[[15]]
#Total plot per season
ggplot(heatmap_data, aes(x = x_bin, y = y_bin, fill = count)) +
  geom_tile(color = "grey80") +
  scale_fill_gradientn(colours = c("#ffffcc", "#ffeda0", "#feb24c", "#fd8d3c", "#f03b20"))+
  coord_fixed() +
  scale_y_reverse() +
  labs(
    title = "Animal movement heatmap by season",
    x = "X grid",
    y = "Y grid",
    fill = "Pass count"
  ) +
  facet_wrap(~season)+
  theme_minimal()
  
#Plot for similarities and anova
ggplot(foraging_similarities, aes(x = season, y = sim_1to4, color = season))+ #y=sim_noncons or y=sim_2days
  geom_boxplot(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) 

anova_model_sim <- aov(sim_1to4 ~ season, data = foraging_similarities)
summary(anova_model_sim)
plot(anova_model_sim)
TukeyHSD(anova_model_sim)
