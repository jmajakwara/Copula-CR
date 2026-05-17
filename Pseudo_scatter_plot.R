# --------------------------------------------------
# Enhanced Wide Pseudo-Observations Plot
# --------------------------------------------------
library(survival)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(viridis)


# Data preparation
data <- my_data %>%
  dplyr::select(time = follow_up_duration, status = event_death) %>%
  dplyr::mutate(time = as.numeric(time),
                status = as.integer(status))

# KM estimates
km_event <- survfit(Surv(time, status == 1) ~ 1, data = data)
km_cens  <- survfit(Surv(time, status == 0) ~ 1, data = data)

S_event <- summary(km_event, times = data$time, extend = TRUE)$surv
S_cens  <- summary(km_cens,  times = data$time, extend = TRUE)$surv

data <- data %>%
  mutate(
    U1 = 1 - S_event,
    U2 = 1 - S_cens,
    Event = factor(status, levels = c(0,1), 
                   labels = c("Censored", "Event (Death)"))
  )

# ----------------------
# Main Plot - Wider & Bigger
# ----------------------
p_main <- ggplot(data, aes(x = U1, y = U2, color = Event)) +
  geom_point(aes(shape = Event), 
             alpha = 0.6, #0.75
             size = 2, #2.08
             position = position_jitter(width = 0.01, height = 0.01)) +
  scale_color_manual(values = c("Censored" = "#1f77b4", 
                                "Event (Death)" = "#d62728")) +
  scale_shape_manual(values = c(1, 16)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "grey50", alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
  labs(title = "Pseudo-observations under competing risks representation",
       subtitle = "Treating right-censoring as a competing risk",
       x = expression(hat(F)[1](t) ~ "(Pseudo-observations for event)"),
       y = expression(hat(F)[2](t) ~ "(Pseudo-observations for censoring)"),
       color = "Status", 
       shape = "Status") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.margin = margin(t = 25, r = 25, b = 25, l = 25)
  ) +
  coord_equal(xlim = c(0,1), ylim = c(0,1))

# ----------------------------------------------------------
# Add marginal densities
# ----------------------------------------------------------
p_final <- ggMarginal(
  p_main,
  type = "density",
  groupColour = TRUE,
  groupFill = TRUE
)

print(p_final)

ggsave(
  filename = "Pseudo_scatter_plot.pdf",
  plot = p_final,
  width = 12,
  height = 10,
  units = "in"
)


