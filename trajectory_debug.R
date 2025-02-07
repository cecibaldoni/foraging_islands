library(dplyr)
library(ggplot2)

foraging_ls <- split(tracking, tracking$unique_trial_ID)

x =foraging_ls[['spring_T1S1_20201103-5']]
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


clean_trajectory <- function(data, max_jump = 20) {
  data <- data %>%
    mutate(x_lag = lag(x, default = first(x)),  #check if the first row has a valid lag
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
          # End of the cluster
          break
        } else {
          #overwrite
          data$x[j] <- anchor_x
          data$y[j] <- anchor_y
        }
      }
    } else {
      # update anchr
      anchor_x <- data$x[i]
      anchor_y <- data$y[i]
    }
  }
  data <- data %>% select(-x_lag, -y_lag, -dist, -is_jump)
  return(data)
}


try <- clean_trajectory(x1)

ggplot(data = x, aes(x, y, colour = frame)) +
   geom_path()
ggplot(data = try, aes(x, y, colour = frame)) +
   geom_path()

