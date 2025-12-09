# Load libraries ----
library(tidyverse)
library(sf)
# library(mapview)
library(parallel)
library(ggplot2)
library(gganimate)
library(here)
library(dplyr)

# Load data ----
tracking <- read.csv(here("csv/merged.csv")) %>%
  rename(season = ID, ID = trial, trial = season) %>%
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

# Main Function ----

all_ls <- lapply(foraging_ls, function(x){
  door <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("door_x", "door_y")) #%>%
  #st_as_sf(coords = c("door_x", "door_y"))
  
  #x <- clean_trajectory(x, door$door_x[[1]], door$door_y[[1]])
  
  door <- door %>% 
    st_as_sf(coords = c("door_x", "door_y"))
  
  door_buffer <- door %>%
    st_buffer(dist = 4)
  
  islands <- c("A", "B", "C", "D")
  hex_ls <- list()
  islands_buffer <- list()
  for (island in islands) {
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
  for (island in islands) {
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
    full_join(intersection_df[c("frame", "island", "visit_seq")]) %>% 
    arrange(frame) %>% 
    mutate(journey = ifelse(!is.na(island), paste0("at_", island, "_", visit_seq), "travelling"))
  
  track_save <- track_sf_2 %>% 
    extract(geometry, c('x', 'y'), '\\((.*), (.*)\\)', convert = TRUE) %>% 
    relocate(x, .after = frame) %>% 
    relocate(y, .after = unique_trial_ID) %>% 
    relocate(unique_trial_ID, .before = ID)
  
  file_path <- here("csv", "foraging_results.csv")
  #write.csv(track_save, file = paste0("~/data/foraging/results/", unique(x$unique_trial_ID),".csv"))
  write.table(track_save, file = file_path, 
              append = TRUE, sep = ",", col.names = !file.exists("file_path"),
              row.names = FALSE)
  
})

#Write result csv
result <- read.csv(here("csv/foraging_results.csv"))
result_ls <- split(result, result$unique_trial_ID)


island_visit <- read.csv(here("csv/island_visit_FR.csv"))

result_final <- map(result_ls, function(df) {
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
    mutate(frame = row_number())  # Re-number frames after filtering
  
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

# check for mismatches
mismatch_rows <- lapply(result_final, function(df) {
  df[df$journey == "travelling" & !is.na(df$island_debug) & df$island_debug != "", ]
})


# trial to adjust the frameshifts but it's not working
move_island_debug <- function(df, travel_letters = c("A", "B", "C", "D")) {
  
  # Copy the debug column so we modify a new one
  new_debug <- result_final$island_debug
  
  for (i in seq_len(nrow(df))) {
    
    letter <- result_final$island_debug[i]
    
    if (!is.na(letter) && letter %in% travel_letters) {
      
      # search window ±10 frames
      start <- max(1, i - 20)
      end   <- min(nrow(df), i + 20)
      
      # find matching letter in island
      match_idx <- which(df$island[start:end] == letter)
      
      if (length(match_idx) > 0) {
        # convert window indices to absolute frame indices
        target_frame <- start + match_idx[1] - 1
        
        # move the debug letter:
        new_debug[target_frame] <- letter
        new_debug[i] <- NA
      }
    }
  }
  
  df$island_debug_aligned <- new_debug
  return(df)
}

# plot
df <- result_final[[which(sapply(result_final, function(df) df$unique_trial_ID[1] == "spring_20201103-1_T1S1"))]]

islands <- coords %>%
  filter(unique_trial_ID == unique(df$unique_trial_ID)) %>%
  select(4:11, 14) %>% 
  mutate(across(everything(), ~ as.numeric(as.character(.))))

df <- df %>%
  mutate(across(c(x, y), ~ as.numeric(as.character(.))))

plot <- df %>% 
  ggplot(aes(x, y, colour = frame)) +
  ggtitle(unique(df$unique_trial_ID)) +
  geom_point(data = islands, aes(x = A_x, y = A_y), size = 10, colour = "red") +
  geom_point(data = islands, aes(x = B_x, y = B_y), size = 10, colour = "blue") +
  geom_point(data = islands, aes(x = C_x, y = C_y), size = 10, colour = "green") +
  geom_point(data = islands, aes(x = D_x, y = D_y), size = 10, colour = "gold") +
  geom_path() +
  geom_point(aes(x = x, y = y), shape = 4, size = 4, 
             color = ifelse(df$food == 1, "green", "red")) +
  scale_y_reverse() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

print(plot)

