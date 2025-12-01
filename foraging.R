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
  #rename(season = ID, ID = trial, trial = season) %>%
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

# TO DO ----

#Add merging logic here, once the island_visit.csv is ready
#Once the merging is done, then the smoothing and cleaning of the dataset can work out


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

# Main Function ----

all_ls <- lapply(foraging_ls, function(x){
  door <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("door_x", "door_y")) #%>%
    #st_as_sf(coords = c("door_x", "door_y"))
  
  x <- clean_trajectory(x, door$door_x[[1]], door$door_y[[1]])
  
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

# TO DO ---- 
#move all the code down here to before the cleaning. The logic should still work, but keep an eyeif everything is merged correctly

## Left join information from the visits of the island to the tracking results

island_visit <- read.csv(here("csv/island_visit_FR.csv"))

## Check if the types are the same and in case we have to edit IDs in all character columns
island_visit <- island_visit %>%
   mutate(across(where(is.character), ~ sub("_(\\d)_", "-\\1_", .)))

## Write back to CSV (overwrite original)
write_csv(island_visit, here("csv/island_visit_FR.csv"))

result_final <- map(result_ls, function(df) {
  df %>%
    mutate(frame = as.integer(frame)) %>%
    left_join(island_visit, by = c("unique_trial_ID", "frame"))
})

##check if the csv is merged with the tracking
#check <- view(result_final[["any unique trial ID"]])

#check for the mismatches among my observations and the pvs
mismatch_rows <- lapply(result_final, function(df) {
  df[df$journey == "travelling" & !is.na(df$island_debug) & df$island_debug != "", ]
})

#List of df with more than 10 mismatch
mismatch_rows_filtered <- mismatch_rows[sapply(mismatch_rows, nrow) > 10]

#check if the observations from the pv have been all manually analysed
result %>% summarise(unique_trial_ID = n_distinct(unique_trial_ID))
island_visit %>% summarise(unique_trial_ID = n_distinct(unique_trial_ID))
missing_ids <- setdiff(result$unique_trial_ID, island_visit$unique_trial_ID)
missing_table <- data.frame(unique_trial_ID = missing_ids)

write.csv(missing_table, "missing_ids.csv", row.names = FALSE)


# Plots (example code) ----

lapply(head(result_final), function(i) { 
  islands <- coords %>%
    filter(unique_trial_ID == unique(i$unique_trial_ID)) %>%
    select(4:11, 14) %>% 
    mutate(across(everything(), ~ as.numeric(as.character(.))))
  
  i <- i %>%
    mutate(across(c(x, y), ~ as.numeric(as.character(.))))
  
  plot <- i %>% 
    ggplot(aes(x, y, colour = frame)) +
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

library(gganimate)

df <- result_final[[which(sapply(result_final, function(df) df$unique_trial_ID[1] == "summer_20210803-3_T2S1"))]]

#df <- result_final[[which(sapply(result_final, function(df) df$unique_trial_ID[1] == "spring_20201103-5_T2S1"))]]

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
  # geom_point(data = df_food, aes(x = x, y = y, color = factor(food)), 
  #            shape = 4, size = 4) +
  # scale_color_manual(values = c("red", "green")) +
  scale_y_reverse() +
  labs(title = 'Frame: {frame}') +
  theme_void() +
  transition_reveal(along = frame) 
 #+ ease_aes('linear')

animate(plot, fps = 1)

str(island_visit)
island_visit_ls <- split(island_visit, island_visit$unique_trial_ID)


df_island <- island_visit_ls$`summer_20210803-3_T2S1`
df_2 <- island_visit_ls[[12]]
str(df_island)
str(df)
