library(tidyverse)
library(sf)
library(here)

# -----------------------------
# 1. Load tracking data
# -----------------------------
foraging_tracking <- read_csv(here("csv/processed/foraging_results.csv")) %>% 
  mutate(
    ID = as.factor(ID),
    season = as.factor(season),
    trial = as.factor(trial),
    unique_trial_ID = as.factor(paste(season, ID, trial, sep = "_"))
  ) %>% 
  relocate(unique_trial_ID, .before = ID)

# cm per pixel conversion
fr_cm_per_pixel <- 0.187192

# -----------------------------
# 2. Load arena coordinates
# -----------------------------
fr_coords <- read_csv(here("csv/islands.csv")) %>%
  separate(A, into = c("A_x", "A_y"), sep = ";", convert = TRUE) %>% 
  separate(B, into = c("B_x", "B_y"), sep = ";", convert = TRUE) %>% 
  separate(C, into = c("C_x", "C_y"), sep = ";", convert = TRUE) %>% 
  separate(D, into = c("D_x", "D_y"), sep = ";", convert = TRUE) %>% 
  separate(door, into = c("door_x", "door_y"), sep = ";", convert = TRUE) %>% 
  mutate(
    unique_trial_ID = as.factor(paste(season, ID, trial, sep = "_")),
    across(contains("_x"), ~ . * fr_cm_per_pixel),
    across(contains("_y"), ~ . * fr_cm_per_pixel)
  )

# -----------------------------
# 3. Load corners data
# -----------------------------
fr_corners <- read_csv(here("csv/foraging_corners.csv")) %>% 
  mutate(
    unique_trial_ID = as.factor(paste(season, ID, trial, sep = "_")),
    across(contains("_x"), ~ . * fr_cm_per_pixel),
    across(contains("_y"), ~ . * fr_cm_per_pixel)
  )

# Keep only trials with tracking data
valid_ids <- unique(fr_coords$unique_trial_ID)
fr_corners <- fr_corners %>% filter(unique_trial_ID %in% valid_ids) %>% droplevels()

# Split tracking data by trial
foraging_ls <- split(foraging_tracking, foraging_tracking$unique_trial_ID)

# -----------------------------
# 4. Edge calculation per trial
# -----------------------------
foraging_edges_file <- here("csv/foraging_edges_new.csv")

calculate_distance <- function(data) {
  data %>%
    st_drop_geometry() %>%
    arrange(frame) %>%
    mutate(
      x_lag = lag(x),
      y_lag = lag(y),
      dist = sqrt((x - x_lag)^2 + (y - y_lag)^2)
    ) %>%
    summarise(total = sum(dist, na.rm = TRUE)) %>%
    pull(total)
}

lapply(foraging_ls, function(x) {
  # Convert tracking points to sf
  track_sf <- x %>% st_as_sf(coords = c("x", "y")) %>%
    mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2])
  
  # Get corners for the trial
  corners <- fr_corners %>%
    filter(as.character(unique_trial_ID) == as.character(unique(x$unique_trial_ID))) %>% 
    select(topleft_x, topleft_y, topright_x, topright_y, downright_x, downright_y, downleft_x, downleft_y)
  
  # Create polygon lines from corners
  corner_coords <- rbind(
    c(corners$topleft_x, corners$topleft_y),
    c(corners$topright_x, corners$topright_y),
    c(corners$downright_x, corners$downright_y),
    c(corners$downleft_x, corners$downleft_y),
    c(corners$topleft_x, corners$topleft_y) # close polygon
  )
  side_lines_sf <- st_sf(geometry = st_sfc(st_linestring(corner_coords)))
  
  # Buffer around edges
  buffer_distance <- 4
  edges_buffer <- st_buffer(side_lines_sf, dist = buffer_distance)
  
  # Detect points at edge
  track_sf <- track_sf %>%
    mutate(at_edge = as.vector(st_intersects(geometry, edges_buffer, sparse = FALSE)))
  
  # Time calculations
  time_edge <- track_sf %>% filter(at_edge == TRUE) %>% summarise(total_time = sum(time, na.rm = TRUE)) %>% pull(total_time)
  time_center <- track_sf %>% filter(at_edge == FALSE) %>% summarise(total_time = sum(time, na.rm = TRUE)) %>% pull(total_time)
  time_total <- time_edge + time_center
  prop_time_edge <- ifelse(time_total > 0, time_edge / time_total, NA)
  
  # Distance calculations
  path_length_edge <- calculate_distance(track_sf %>% filter(at_edge == TRUE))
  path_length_center <- calculate_distance(track_sf %>% filter(at_edge == FALSE))
  total_path_length <- path_length_edge + path_length_center
  prop_path_length_edge <- ifelse(total_path_length > 0, path_length_edge / total_path_length, NA)
  
  # Combine results
  df <- data.frame(
    task = "foraging",
    x[1,1], x[1,2], x[1,3], x[1,4],
    time_total = time_total, time_edge = time_edge, time_center = time_center, prop_time_edge = prop_time_edge,
    path_length_edge = path_length_edge, path_length_center = path_length_center, prop_path_length_edge = prop_path_length_edge
  ) %>% droplevels()
  
  # Save to CSV
  write.table(df, file = foraging_edges_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(foraging_edges_file))
})

