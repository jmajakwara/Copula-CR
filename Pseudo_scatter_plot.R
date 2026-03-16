library(survival)
library(copula)
library(ggplot2)
library(dplyr)
library(gridExtra)

my_data <- read.csv("Breast Cancer Data.csv")
data <- my_data |>
   dplyr::select(time = follow_up_duration, status = event_death) |>
   dplyr::mutate(status = as.integer(status), time = as.numeric(time))


# Quick check
cat("n =", nrow(data), "\n")
cat("Event rate =", mean(data$status), "\n")
cat("Censoring rate =", mean(1 - data$status), "\n")


#  Compute censoring-adjusted pseudo-observations (ranks)


# Kaplan-Meier survival estimator (for censoring adjustment)
km_fit <- survfit(Surv(time, status) ~ 1, data = data)

# Inverse probability of censoring weights (IPCW) style ranks
# Rank for event times: use KM survival at time just before event
# Rank for censored times: use KM survival at censoring time


times_ordered <- sort(unique(data$time))

# KM survival at each unique time (left-continuous)
S_km <- summary(km_fit, times = times_ordered - 1e-6)$surv

# Map back to each observation
data <- data %>%
  dplyr::mutate(
    # Approximate KM survival just before this time
    S_before = approx(times_ordered, S_km, xout = time, method = "constant", rule = 2)$y,
    
    # Pseudo-observation (rank adjustment for censoring)
    u = ifelse(
      status == 1,
      # For events: rank among all times up to this point
      (rank(time, ties.method = "average") - 0.5) / nrow(data),
      # For censored: use survival probability at censoring time
      S_before
    )
  )



# Split into pseudo-observations for events and censored
u_event <- data$u[data$status == 1]
u_cens  <- data$u[data$status == 0]



# Scatter plot of pseudo-observations

ggplot(data, aes(x = rank(time)/nrow(data), y = u)) +
  geom_point(aes(color = factor(status), shape = factor(status)), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("0" = "grey60", "1" = "darkred")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(
    title = "Pseudo observations (censoring adjusted ranks)",
    subtitle = "Events (red) vs Censored (grey)",
    x = "Rank based pseudo-value (marginal)",
    y = "Censoring adjusted pseudo-value",
    color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  coord_equal(xlim = c(0,1), ylim = c(0,1))





