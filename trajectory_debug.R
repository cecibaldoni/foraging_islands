library(dplyr)
library(ggplot2)

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

