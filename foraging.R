library(tidyverse)
library(sf)
library(mapview)
library(parallel)
library(ggplot2)
library(gganimate)

# tracking <- read.csv("~/data/foraging/csv/merged.csv") %>%
  mutate(season = tolower(season))
cmperpixel = 0.187192
coords <- read.csv("~/data/foraging/islands.csv") %>%
  separate(A, into = c("A_x", "A_y"), sep = ";", convert = TRUE) %>% 
  separate(B, into = c("B_x", "B_y"), sep = ";", convert = TRUE) %>% 
  separate(C, into = c("C_x", "C_y"), sep = ";", convert = TRUE) %>% 
  separate(D, into = c("D_x", "D_y"), sep = ";", convert = TRUE) %>% 
  separate(door, into = c("door_x", "door_y"), sep = ";", convert = TRUE) %>% 
  mutate(unique_trial_ID = as.factor(paste(season, trial, ID, sep = "_")),
         ID = as.factor(ID),
         season = as.factor(season),
         trial = as.factor(trial),
         across(contains("_x"), ~ . * cmperpixel),
         across(contains("_y"), ~ . * cmperpixel))

tracking[c('ANGLE')][sapply(tracking[c('ANGLE')], is.infinite)] <- NA #transform inf values in NA, then drop them
tracking <- tracking %>% 
  drop_na(ANGLE) %>% 
  mutate(unique_trial_ID = paste(season, trial, ID, sep = "_"))

foraging_ls <- split(tracking, tracking$unique_trial_ID)

x = foraging_ls[['spring_T1S1_20201103-5']]
door <- coords %>%
  filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
  dplyr::select(c("door_x", "door_y"))
str(door)

x1 <- x %>%
  # Step 1: Remove frames at the beginnin, far from the door
  mutate(dist_from_door = sqrt((x - door$door_x[[1]])^2 + (y - door$door_y[[1]])^2),
         before_close = cumsum(dist_from_door <= 15) == 0) %>%
  filter(!before_close) %>%
  select(-dist_from_door, -before_close) %>% 
  mutate(frame = row_number())

ggplot(data = x1[100:300,], aes(x, y, colour = frame)) +
  geom_path()



ggplot(data = try, aes(x, y, colour = frame)) +
  geom_path()

x %>% 
  ggplot(aes(x, y, colour = frame)) +
  ggtitle(x$unique_trial_ID) +
  geom_point(data = coords %>% filter(unique_trial_ID == unique(x$unique_trial_ID)), aes(x = A_x, y = A_y), size = 10, colour = "red") +
  geom_point(data = coords %>% filter(unique_trial_ID == unique(x$unique_trial_ID)), aes(x = B_x, y = B_y), size = 10, colour = "blue") +
  geom_point(data = coords %>% filter(unique_trial_ID == unique(x$unique_trial_ID)), aes(x = C_x, y = C_y), size = 10, colour = "green") +
  geom_point(data = coords %>% filter(unique_trial_ID == unique(x$unique_trial_ID)), aes(x = D_x, y = D_y), size = 10, colour = "gold") +
  geom_path()
#Need to smooth movement

all_ls <- lapply(foraging_ls, function(x){
  x = foraging_ls[['spring_T1S1_20201103-5']]
  door <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("door_x", "door_y")) %>%
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
    relocate(y, .after = x) %>% 
    relocate(unique_trial_ID, .before = ID)
  
  #write.csv(track_save, file = paste0("~/data/foraging/results/", unique(x$unique_trial_ID),".csv"))
  write.table(track_save, file = "~/data/foraging/results/all_trials.csv", 
              append = TRUE, sep = ",", col.names = !file.exists("~/data/foraging/results/all_trials.csv"),  # Write headers only once
              row.names = FALSE)
  
})

result <- read.csv("~/data/foraging/results/all_trials.csv")
result_ls <- split(result, result$unique_trial_ID)

#Plots

lapply(head(result_ls), function(i) { 
  islands <- coords %>%
    filter(unique_trial_ID == unique(i$unique_trial_ID)) %>%
    select(4:11, 14)
  
  plot <- i %>% 
    ggplot(aes(x, y, colour = frame)) +
    ggtitle(i$unique_trial_ID) +
    geom_point(data = islands, aes(x = A_x, y = A_y), size = 10, colour = "red") +
    geom_point(data = islands, aes(x = B_x, y = B_y), size = 3, colour = "blue") +
    geom_point(data = islands, aes(x = C_x, y = C_y), size = 3, colour = "green") +
    geom_point(data = islands, aes(x = D_x, y = D_y), size = 3, colour = "gold") +
    geom_path() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  print(plot)
})

