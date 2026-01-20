# Load libraries ----

library(tidyverse)
library(sf)
# library(mapview)
library(parallel)
library(gganimate)
library(here)
library(trajr)

# Load and tidy data ----
tracking <- read.csv(here("csv/merged.csv")) %>%
 # rename(season = ID, ID = trial, trial = season) %>%
  mutate(season = tolower(season))  
  
cmperpixel = 0.187192

coords <- read.csv(here("csv/islands.csv")) %>%
  separate(A, into = c("A_x", "A_y"), sep = ";", convert = TRUE) %>% 
  separate(B, into = c("B_x", "B_y"), sep = ";", convert = TRUE) %>% 
  separate(C, into = c("C_x", "C_y"), sep = ";", convert = TRUE) %>% 
  separate(D, into = c("D_x", "D_y"), sep = ";", convert = TRUE) %>% 
  separate(door, into = c("door_x", "door_y"), sep = ";", convert = TRUE) %>% 
  mutate(unique_trial_ID = as.factor(paste(season, ID, trial, sep = "_")),
         ID = as.factor(ID),
         season = as.factor(season),
         trial = as.factor(trial),
         across(contains("_x"), ~ . * cmperpixel),
         across(contains("_y"), ~ . * cmperpixel))

tracking[c('ANGLE')][sapply(tracking[c('ANGLE')], is.infinite)] <- NA #transform inf values in NA, then drop them
tracking <- tracking %>% 
  drop_na(ANGLE) %>% 
  mutate(unique_trial_ID = paste(season, ID, trial, sep = "_"))

foraging_ls <- split(tracking, tracking$unique_trial_ID)

island_visit <- read.csv(here("csv/island_visit.csv"))

#foraging append is the tracking with manual observation. To do before the cleaning to avoid frame mismatch.
foraging_append <- map(foraging_ls, function(df) {
  df %>%
    mutate(frame = as.integer(frame)) %>%
    left_join(island_visit, by = c("unique_trial_ID", "frame"))
})

# Smooth and clean tracking data ----
clean_trajectory <- function(data, door_x, door_y, max_jump = 20) {
  data <- data %>%
    mutate(dist_from_door = sqrt((x - door_x)^2 + (y - door_y)^2),
           before_close = cumsum(dist_from_door <= 15) == 0) %>%
    filter(!before_close) %>%
    select(-dist_from_door, -before_close) %>%
    mutate(frame = row_number(), # re-number frames after filtering
           time = (frame - 1) / 30) #re-number time after filtering 
  
          
  
  # detect and fix clusters of jumps
  data <- data %>%
    mutate(x_lag = lag(x, default = first(x)),
           y_lag = lag(y, default = first(y)),
           dist = sqrt((x - x_lag)^2 + (y - y_lag)^2),
           is_jump = dist > max_jump) %>%
    mutate(is_jump = replace_na(is_jump, FALSE))
  
  anchor_x <- data$x[1]
  anchor_y <- data$y[1]
  
  for (i in seq_len(nrow(data))) {
    if (data$is_jump[i]) {
      for (j in i:nrow(data)) {
        dist_to_anchor <- sqrt((data$x[j] - anchor_x)^2 + (data$y[j] - anchor_y)^2)
        if (dist_to_anchor <= max_jump) {
          break
        } else {
          data$x[j] <- anchor_x
          data$y[j] <- anchor_y
        }
      }
    } else {
      anchor_x <- data$x[i]
      anchor_y <- data$y[i]
    }
  }
  data <- data %>% select(-x_lag, -y_lag, -dist, -is_jump)
  return(data)
}

#trajectory calculations function

safe_traj <- function(df) {
  df <- df %>% filter(!is.na(x), !is.na(y), !is.na(time))
  if (nrow(df) < 2) return(NULL)
  
  TrajFromCoords(df, xCol = "x", yCol = "y", fps = 30)
}

#define all possible combinations of islands and locations (from A to D, from 1 to 6: 24 combinations)
islands_ref <- c("A", "B", "C", "D")
food_ref    <- 1:6
all_possible_islands <- expand.grid(food = food_ref, island = islands_ref) %>%
  mutate(col_name = paste0(island, "_", food)) %>%
  pull(col_name)

# Main Function ----

#define paths
results_file <- here("csv/processed", "foraging_results.csv")
master_file <- here("csv/processed", "foraging_master.csv")

all_ls <- lapply(foraging_append, function(x){
  ## to try the loop with only one list element, run the next line
  # x = foraging_append[[1]]
  message("Processing: ", unique(x$unique_trial_ID))
  
  door <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("door_x", "door_y")) #%>%
    #st_as_sf(coords = c("door_x", "door_y"))
  
  x <- clean_trajectory(x, door$door_x[[1]], door$door_y[[1]])
  
  if (nrow(x) == 0) {
    message("   ->  No data after cleaning")
    return(NULL)
  }
  
  door <- door %>% 
    st_as_sf(coords = c("door_x", "door_y"))
  
  door_buffer <- door %>%
    st_buffer(dist = 4)
  
  hex_ls <- list()
  islands_buffer <- list()
  
  for (island in islands_ref) {
    filtered_coords <- coords %>%
      filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
      dplyr::select(paste0(island, "_x"), paste0(island, "_y")) %>%
      mutate(!!paste0(island, "_x") := as.numeric(!!sym(paste0(island, "_x"))),
        !!paste0(island, "_y") := as.numeric(!!sym(paste0(island, "_y"))))
    center <- as.numeric(filtered_coords[1, ])
    distance <- 9
    hex_coords <- matrix(nrow = 6, ncol = 2)
    for (i in 1:6) {
      angle <- (i - 1) * 2 * pi / 6
      hex_coords[i, ] <- center + distance * c(cos(angle), sin(angle))
    }
    hex_closed <- rbind(hex_coords, hex_coords[1, ])
    hex_sf <- st_polygon(list(hex_closed))
    hex_ls[[island]] <- hex_sf
    island_buffer <- hex_sf %>%
      st_buffer(dist = 4)
    islands_buffer[[island]] <- island_buffer
  }
  
  x <- x %>%
    arrange(frame) %>%
    mutate(x_lag = lag(x), 
           y_lag = lag(y),
           dist = sqrt((x - x_lag)^2 + (y - y_lag)^2)) %>%
    filter(is.na(dist) | dist < 10) %>%
    select(-x_lag, -y_lag, -dist)
  
  track_sf <- x %>%
    st_as_sf(coords = c("x", "y"))
  
  intersection_df <- data.frame()
  for (island in islands_ref) {
    at_island <- track_sf %>%
      #st_intersection(hex_ls[[island]]) %>%
      st_intersection(islands_buffer[[island]]) %>% 
      as.data.frame() %>%
      arrange(frame) %>%
      mutate(island = island) %>%
      mutate(timediff = frame - lag(frame)) %>%
      mutate(new_timediff = ifelse(is.na(timediff) | timediff != 1, 1, 0)) %>%
      mutate(visit_seq = cumsum(new_timediff))
    intersection_df <- rbind(intersection_df, at_island)
  }
  
  track_sf_2 <- track_sf %>%
    full_join(intersection_df[c("frame", "island", "visit_seq")], by = "frame") %>% 
    arrange(frame) %>% 
    mutate(journey = ifelse(!is.na(island), paste0("at_", island, "_", visit_seq), "travelling"))
  
  
  track_save <- track_sf_2 %>% 
    extract(geometry, c('x', 'y'), '\\((.*), (.*)\\)', convert = TRUE) %>% 
    relocate(x, .after = frame) %>% 
    relocate(y, .after = x) %>% 
    relocate(unique_trial_ID, .before = ID)
  
  
  write.table(track_save, file = results_file, 
              append = TRUE, sep = ",", col.names = !file.exists(results_file),
              row.names = FALSE)
  
  ## TRAJECTORY/TIME/... CALCULATIONS ----
  
  # first time at A or D
  first_ad <- track_save %>%
    filter(island %in% c("A", "D")) %>%
    slice_min(frame, n = 1, with_ties = FALSE)
  
  # Island/Food Counts (A_0, A_1, etc.)
  island_counts <- track_save %>%
    filter(!is.na(island), !is.na(food)) %>%
    count(island, food) %>%
    mutate(col_name = paste0(island, "_", food)) %>%
    select(col_name, n) %>%
    tidyr::pivot_wider(names_from = col_name, values_from = n, values_fill = 0)
  
  #add all the 24 possible combinations into 24 columns, even if they were not visited
  missing_cols <- setdiff(all_possible_islands, names(island_counts))
  if (length(missing_cols) > 0) {
    island_counts[missing_cols] <- 0
  }
  
  #fix the column order, to avoid errors later
  island_counts <- island_counts %>% select(all_of(all_possible_islands))
  
  # trajectory metrics (with the function defined before the loop)
  traj_all <- safe_traj(track_save)
  
  #compile the Master df
  master_row <- data.frame(
    #add info columns
    unique_trial_ID    = x$unique_trial_ID[1],
    season             = x$season[1],
    ID                 = x$ID[1],
    trial              = x$trial[1],
    #add metrics
    first_AD_time      = if(nrow(first_ad) > 0) round(first_ad$frame[1] / 30, 2) else NA,
    distance_total_cm  = if(!is.null(traj_all)) round(TrajLength(traj_all), 2) else NA,
    straightness_total = if(!is.null(traj_all)) round(TrajStraightness(traj_all), 5) else NA,
    travelling_time    = round(sum(track_save$journey == "travelling") / 30, 2),
    moving_time        = round(max(track_save$frame) / 30, 2),
    nonmoving_time     = round((43200 - max(track_save$frame)) / 30, 2))
  
   
  # attache the island_counts column, with all 24 combinations
  if (nrow(island_counts) > 0) {
    master_row <- bind_cols(master_row, island_counts %>% select(-any_of("unique_trial_ID")))
  }
  
  #TODO: Think if you need more information. 
  # EXAMPLE 1: Instead of First AD, why not check the first island visited regardless? If you have that information in a column, it's easy to filter A and D inside that. While if you have it already filtered, you don't know if they went to the empty islands before going in the baited one
  # EXAMPLE 2: would it be possible to investigate the instances at each position (A1, vs C4) and how many times it happened and when? this maybe is already possible, I didn't think about it yet
  
  
  # save to Master CSV
  write.table(master_row, file = master_file, append = TRUE, sep = ",", 
              col.names = !file.exists(master_file), row.names = FALSE)
  
  return(master_row)
})

#Processing: summer_20210802-3_T2S1
# ->  No data after cleaning
# Processing: summer_20210803-1_T1S1
# ->  No data after cleaning


#Write result csv
result <- read.csv(here("csv/processed/foraging_results.csv"))
result_ls <- split(result, result$unique_trial_ID)

df <- result_ls[[1]]


# Plots (example code) ----

#Run this instead of lapply for a specific ID
target_id <- "summer_20210803-4_T1S2"
i <- result_ls[[target_id]]
i <- result_ls[[which(sapply(result_ls, function(x)
  unique(x$unique_trial_ID) == target_id))]]

lapply(head(result_ls), function(i) { 
  islands <- coords %>%
    filter(unique_trial_ID == unique(i$unique_trial_ID)) %>%
    select(4:11, 14) %>% 
    mutate(across(everything(), ~ as.numeric(as.character(.))))
  
  i <- i %>%
    mutate(across(c(x, y), ~ as.numeric(as.character(.))))
  
  plot <- i %>% 
    ggplot(aes(x, y, group = 1, colour = frame)) +
    ggtitle(i$unique_trial_ID) +
    geom_point(data = islands, aes(x = A_x, y = A_y), size = 10, colour = "red") +
    geom_point(data = islands, aes(x = B_x, y = B_y), size = 10, colour = "blue") +
    geom_point(data = islands, aes(x = C_x, y = C_y), size = 10, colour = "green") +
    geom_point(data = islands, aes(x = D_x, y = D_y), size = 10, colour = "gold") +
    geom_path() +
    
    geom_point(aes(x = x, y = y), shape = 4, size = 4, 
               color = ifelse(i$food == 1, "green", "red")) +
    
    scale_y_reverse() +
    
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  print(plot)
})

#Animated plot
library(gganimate)

df <- result_final[[which(sapply(result_final, function(df) df$unique_trial_ID[1] == "summer_20210803-4_T1S2"))]]

islands <- coords %>%
  filter(unique_trial_ID == df$unique_trial_ID[1]) %>%
  select(A_x, A_y, B_x, B_y, C_x, C_y, D_x, D_y) %>%
  mutate(across(everything(), as.numeric))

df_food <- df %>% filter(!is.na(food) & food %in% c(0, 1))

df <- df %>%
  mutate(across(c(x, y), ~ as.numeric(as.character(.))))

plot <- ggplot() +
  geom_point(aes(x = islands$A_x, y = islands$A_y), color = "red", size = 10) +
  geom_point(aes(x = islands$B_x, y = islands$B_y), color = "blue", size = 10) +
  geom_point(aes(x = islands$C_x, y = islands$C_y), color = "green", size = 10) +
  geom_point(aes(x = islands$D_x, y = islands$D_y), color = "gold", size = 10) +
  geom_path(data = df, aes(x = x, y = y, group = 1, color = frame)) +
  geom_point(data = df_food, aes(x = x, y = y, color = factor(food)),
             shape = 4, size = 4) +
  scale_color_manual(values = c("red", "green")) +
  scale_y_reverse() +
  labs(title = 'Frame: {frame}') +
  theme_void() +
  transition_reveal(along = frame) 
 #+ ease_aes('linear')

animate(plot, fps = 50)


#Various checks
str(island_visit)
island_visit_ls <- split(island_visit, island_visit$unique_trial_ID)
df_island <- island_visit_ls$`summer_20210803-3_T2S1`
df_2 <- island_visit_ls[[12]]
str(df_island)
str(df)
island_visit <- read.csv(here("csv/island_visit_FR.csv"))
island_separate <- island_visit %>%
  separate(unique_trial_ID, into = c("season", "ID", "session"), sep = "_")

#check for the observation per session in each season
obs_count <- island_visit %>% 
  separate(unique_trial_ID, into = c("season", "ID", "session"), sep = "_") %>%
  count(season, session, name = "n_observations")

#check if I analized 4 session for each individual
sessions_count <- island_visit %>%
  separate(unique_trial_ID, into = c("season", "id", "session"), sep = "_") %>%
  distinct(season, id, session) %>%        
  count(season, id, name = "n_sessions")   

##check if the csv is merged with the tracking
check <- view(result_final[["any unique trial ID"]])

#check for the mismatches among my observations and the pvs
mismatch_rows <- lapply(result_ls, function(df) {
  df[df$journey == "travelling" & !is.na(df$island_debug) & df$island_debug != "", ]
})

mismatch_table <- data.frame(
  df = seq_along(result_ls),
  mismatches = sapply(result_ls, function(df) {
    sum(df$island_debug != df$island, na.rm = TRUE)
  })
)

#List of df with more than 10 mismatch
mismatch_rows_filtered <- mismatch_rows[sapply(mismatch_rows, nrow) > 10]

#check if the observations from the pv have been all manually analysed
result %>% summarise(unique_trial_ID = n_distinct(unique_trial_ID))
island_visit %>% summarise(unique_trial_ID = n_distinct(unique_trial_ID))
missing_ids <- setdiff(result$unique_trial_ID, island_visit$unique_trial_ID)
missing_table <- data.frame(unique_trial_ID = missing_ids)
write.csv(missing_table, "missing_ids.csv", row.names = FALSE)
