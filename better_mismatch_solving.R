# Load libraries ----

library(tidyverse)
library(sf)
# library(mapview)
library(parallel)
library(ggplot2)
library(gganimate)
library(here)
library(dplyr)

# Load and tidy data ----
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

#At this point I have tracking whit the positions at each frame and the coordinates of the islands
#Next step is to create the islands and the buffer, then merge and add the column for travelling

x <- foraging_ls[["summer_20210803-3_T2S1"]] #Just one object

# Get door coordinates
door <- coords %>%
  filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
  dplyr::select(c("door_x", "door_y"))

# set the door
x <- clean_trajectory(x, door$door_x[[1]], door$door_y[[1]])

# convert door to sf + buffer
door_sf <- door %>%
  st_as_sf(coords = c("door_x", "door_y"))

door_buffer <- door_sf %>%
  st_buffer(dist = 4)

#Build hexagons and buffers for each island
islands <- c("A", "B", "C", "D")
hex_ls <- list()
islands_buffer <- list()

for (island in islands) {
  
  filtered_coords <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(paste0(island, "_x"), paste0(island, "_y")) %>%
    mutate(
      !!paste0(island, "_x") := as.numeric(!!sym(paste0(island, "_x"))),
      !!paste0(island, "_y") := as.numeric(!!sym(paste0(island, "_y")))
    )
  
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
  islands_buffer[[island]] <- st_buffer(hex_sf, dist = 4)
}

# --- Clean jumps in the track ---
# x <- x %>%
#   arrange(frame) %>%
#   mutate(
#     x_lag = lag(x),
#     y_lag = lag(y),
#     dist = sqrt((x - x_lag)^2 + (y - y_lag)^2)
#   ) %>%
#   filter(is.na(dist) | dist < 10) %>%
#   select(-x_lag, -y_lag, -dist)

# convert track to sf
track_sf <- x %>%
  st_as_sf(coords = c("x", "y"))

# --- Intersection with island buffers ---
intersection_df <- data.frame()

for (island in islands) {
  at_island <- track_sf %>%
    st_intersection(islands_buffer[[island]]) %>%
    as.data.frame() %>%
    arrange(frame) %>%
    mutate(
      island = island,
      timediff = frame - lag(frame),
      new_timediff = ifelse(is.na(timediff) | timediff != 1, 1, 0),
      visit_seq = cumsum(new_timediff)
    )
  
  intersection_df <- rbind(intersection_df, at_island)
}

# --- Join back to track and classify journey ---
track_sf_2 <- track_sf %>%
  full_join(intersection_df[c("frame", "island", "visit_seq")]) %>%
  arrange(frame) %>%
  mutate(journey = ifelse(!is.na(island),
                          paste0("at_", island, "_", visit_seq),
                          "travelling"))

# --- Extract coordinates back out of sf ---
track_save <- track_sf_2 %>%
  extract(geometry, c('x', 'y'), '\\((.*), (.*)\\)', convert = TRUE) %>%
  relocate(x, .after = frame) %>%
  relocate(y, .after = unique_trial_ID) %>%
  relocate(unique_trial_ID, .before = ID)

# --- Save output ---
file_path <- here("csv", "foraging_results.csv")

# Here I join my manual observations
island_visit <- read.csv(here("csv/island_visit_FR.csv"))

final_result <- track_save %>%
  left_join(island_visit, by = c("unique_trial_ID", "frame"))

# Plot

# i = a single dataframe
i <- final_result 

# Prepare islands coordinates
islands <- coords %>%
  filter(unique_trial_ID == unique(i$unique_trial_ID)) %>%
  select(4:11, 14) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

# Ensure x and y are numeric
i <- i %>%
  mutate(across(c(x, y), ~ as.numeric(as.character(.))))

# Make the plot
ggplot(i, aes(x = x, y = y, colour = frame)) +
  ggtitle(unique(i$unique_trial_ID)) +
  
  # Add island points
  geom_point(data = islands, aes(x = A_x, y = A_y), size = 10, colour = "red") +
  geom_point(data = islands, aes(x = B_x, y = B_y), size = 10, colour = "blue") +
  geom_point(data = islands, aes(x = C_x, y = C_y), size = 10, colour = "green") +
  geom_point(data = islands, aes(x = D_x, y = D_y), size = 10, colour = "gold") +
  
  # Trajectory path
  geom_path() +
  
  # Points with food info
  geom_point(aes(x = x, y = y),
             shape = 4, size = 4,
             color = ifelse(i$food == 1, "green", "red")) +
  
  # Flip y axis if needed
  scale_y_reverse() +
  
  # Remove axis labels/ticks
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
