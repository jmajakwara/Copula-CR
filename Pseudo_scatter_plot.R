#--------------------------------------------------
# Competing risks pseudo-observations
#--------------------------------------------------

library(survival)
library(dplyr)
library(ggplot2)

# Data
data <- my_data %>%
  dplyr::select(time = follow_up_duration,
                status = event_death) %>%
  dplyr::mutate(time = as.numeric(time),
         status = as.integer(status))

n <- nrow(data)

#--------------------------------------------------
# KM for each "cause"
#--------------------------------------------------

# Cause 1: Event
km_event <- survfit(Surv(time, status == 1) ~ 1, data = data)

# Cause 2: Censoring (treated as competing event)
km_cens  <- survfit(Surv(time, status == 0) ~ 1, data = data)

#--------------------------------------------------
# Extract survival just before time
#--------------------------------------------------

S_event <- summary(km_event, times = data$time, extend = TRUE)$surv
S_cens  <- summary(km_cens,  times = data$time, extend = TRUE)$surv

#--------------------------------------------------
# Construct pseudo-observations
#--------------------------------------------------

data <- data %>%
  dplyr::mutate(
    U1 = 1-S_event,   # event margin
    U2 = 1-S_cens     # censoring margin
  )

#--------------------------------------------------
# Plot (true copula space)
#--------------------------------------------------
data$event <- factor(data$status)
ggplot(data, aes(x = U1, y = U2)) +
  geom_point(aes(color = event, shape = event),
             alpha = 0.6,
             size = 2,
             position = position_jitter(width = 0.01, height = 0.01)) +
  scale_color_manual(values = c("0" = "grey60", "1" = "darkred")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(
    title = "Pseudo-observations under competing risks representation",
    subtitle = "Event (red) vs Censoring treated as competing event (grey)",
    x = expression(hat(F)[1](t^{`-`})),
    y = expression(hat(F)[2](t^{`-`}))) +
  theme_minimal() +
  coord_equal(xlim = c(0,1), ylim = c(0,1))

# Create matrix with U1 and U2 (left-continuous version)
U <- cbind(U1, U2)


#--------------------------------------------------
# Plot using splom2 from copula package
#--------------------------------------------------

splom2(U, 
       varnames = c("Event", "Censoring"),
       cex = 0.4,
       col.mat = "red",
       main = "Pseudo-observations under competing risks representation")


