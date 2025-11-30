# Load libraries ----

library(tidyverse)
library(sf)
library(mapview)
library(parallel)
library(ggplot2)
library(gganimate)
library(here)
library(dplyr)

# 1. Merge island visits with tracking 

 ##tracking
tracking <- read.csv(here("csv/merged.csv")) %>%
  #rename(season = ID, ID = trial, trial = season) %>%
  mutate(season = tolower(season)) %>% 
  unite("unique_trial_ID", season, ID, trial, sep = "_")

summer <- tracking %>% 
  filter(unique_trial_ID == "summer_20210803-3_T2S1")
 
 ##island visit
island_visit <- read.csv(here("csv/island_visit_FR.csv"))%>% 
  filter(unique_trial_ID == "summer_20210803-3_T2S1")
 
 ##merging island visits and tracking
merged <- summer %>%
  left_join(island_visit, by = c("unique_trial_ID", "frame"))

# 2. include the shape of the arena, the feeders and merge it to the tracking

cmperpixel = 0.187192

arena_coords <- read.csv(here("csv/islands.csv")) %>%
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

arena_summer <- arena_coords %>% 
  filter(unique_trial_ID =="summer_20210803-3_T2S1") %>% 
  unite("unique_trial_ID", season, ID, trial, sep = "_")

## Convert position points (frames & x/y) to sf points
merged_sf <- st_as_sf(merged, coords = c("x", "y"), crs = NA)  

## Convert arena points to sf points
hex_df <- arena_summer %>%
  select(unique_trial_ID, A_x, A_y, B_x, B_y, C_x, C_y, D_x, D_y) %>%
  pivot_longer(cols = -unique_trial_ID, names_to = c("island", ".value"), names_pattern = "([A-D])_(x|y)")

hex_sf <- hex_df %>%
  st_as_sf(coords = c("x", "y"), crs = NA)

## Define buffer radius (adjust as needed) and hexagonal buffers
buffer_radius <- 13  
hex_buffers <- st_buffer(hex_sf, dist = buffer_radius, nQuadSegs = 6)

## Join merged points with hex buffers
merged_with_islands <- st_join(merged_sf, hex_buffers, join = st_intersects)

# 3. Plot
ggplot() +
  ## plot hexagonal buffers
  geom_sf(data = hex_buffers, fill = "lightblue", color = "blue", alpha = 0.4) +
  ## plot island points
  geom_sf(data = hex_sf, color = "darkblue", size = 3) +
  ## plot animal positioning
  geom_point(data = merged %>% 
               filter (!is.na(food)),
             aes(x = x, y = y, color = factor(food)), size = 1, shape = 4) +
  scale_color_manual(values = c("0" = "red", "1" = "green")) +

  theme_minimal() +
  labs(title = "Animal visits and Island Hex Buffers",
       x = "X", y = "Y")




door <- coords %>%
  filter(unique_trial_ID == "summer_20210803-3_T2S1") %>%
  dplyr::select(c("door_x", "door_y")) #%>%

x <- clean_trajectory(merged, door$door_x[[1]], door$door_y[[1]])
str(x)
str(merged)

ggplot() +
  geom_path(data = merged, aes(x=x, y=y, colour= frame)) +
  geom_point(aes(x = islands$A_x, y = islands$A_y), color = "red", size = 10) +
  geom_point(aes(x = islands$B_x, y = islands$B_y), color = "blue", size = 10) +
  geom_point(aes(x = islands$C_x, y = islands$C_y), color = "green", size = 10) +
  geom_point(aes(x = islands$D_x, y = islands$D_y), color = "gold", size = 10)

ggplot() +
  geom_path(data = x, aes(x=x, y=y, colour= frame)) +
  geom_point(aes(x = islands$A_x, y = islands$A_y), color = "red", size = 10) +
  geom_point(aes(x = islands$B_x, y = islands$B_y), color = "blue", size = 10) +
  geom_point(aes(x = islands$C_x, y = islands$C_y), color = "green", size = 10) +
  geom_point(aes(x = islands$D_x, y = islands$D_y), color = "gold", size = 10)

ggplot() +
  ## plot hexagonal buffers
  geom_sf(data = hex_buffers, fill = "lightblue", color = "blue", alpha = 0.4) +
  ## plot island points
  geom_sf(data = hex_sf, color = "darkblue", size = 3) +
  ## plot animal positioning
  geom_point(data = x %>% 
               filter (!is.na(food)),
             aes(x = x, y = y, color = factor(food)), size = 1, shape = 4) +
  scale_color_manual(values = c("0" = "red", "1" = "green")) +
  labs(title = "Animal visits and Island Hex Buffers",
       x = "X", y = "Y")
