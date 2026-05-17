library(survival)
library(copula)
library(ggplot2)
library(gridExtra)


# ----------------------
# Clayton / Gumbel RW Estimator
# ----------------------

rw_estimator <- function(Y, delta, theta, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  # Clayton generator: phi(u) = (u^{-theta} - 1)/theta
  # phi^{-1}(s) = (1 + theta*s)^(-1/theta)
  if (family == "clayton") {
    copula_gen <- function(u) (u^(-theta) - 1) / theta
    copula_gen_inv <- function(s) (1 + theta * s)^(-1/theta)
  } else {
    # Gumbel: generator phi(u) = (-log u)^theta, phi^{-1}(s) = exp(- s^(1/theta))
    copula_gen <- function(u) (-log(u))^theta
    copula_gen_inv <- function(s) exp(- s^(1/theta))
  }
  
  n <- length(Y)
  # Use unique event times only
  event_times <- sort(unique(Y[delta == 1]))
  if (length(event_times) == 0) {
    stop("No events observed (delta==1). RW estimator cannot be computed.")
  }
  
  # Empirical survive function for pi_hat(t) = P(Y > t)
  pi_hat_vals <- sapply(event_times, function(t) sum(Y > t) / n)
  
  sum_term <- 0

rw_surv <- numeric(length(event_times))

  for (i in seq_along(event_times)) {
    pi0 <- pi_hat_vals[i]
    pi1 <- pi0 + 1 / n
    # keep arguments in (0,1)
    pi0 <- min(max(pi0, 1e-12), 1 - 1e-12)
    pi1 <- min(max(pi1, 1e-12), 1 - 1e-12)
    # increment = - [phi(pi1) - phi(pi0)]  (negate the difference)
    increment <- copula_gen(pi0) - copula_gen(pi1)
    # accumulate (should be nonnegative in theory)
    sum_term <- sum_term + increment

    # guard numerical negativity
    s_val <- pmax(sum_term, 0)
    rw_surv[i] <- copula_gen_inv(s_val)  # apply inverse  
  }
    
  return(data.frame(time = event_times, surv = rw_surv))
}


#### Simulation ####

simulate_data <- function(n, weibull_T, weibull_C, theta, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (family == "clayton") {
    cop <- archmCopula(family = "clayton", param = theta, dim = 2)
  } else {
    # Gumbel in copula package uses 'gumbel' name; parameter >= 1
    cop <- archmCopula(family = "gumbel", param = theta, dim = 2)
  }
  U <- rCopula(n, cop) # n x 2 matrix of uniforms
  U1 <- U[, 1]
  U2 <- U[, 2]
  
  # Weibull inverse CDF (parametrisation: F(t)=1-exp(-(t/scale)^shape))
  T <- weibull_T$scale * (-log(U1))^(1 / weibull_T$shape)
  C <- weibull_C$scale * (-log(U2))^(1 / weibull_C$shape)
  
  Y <- pmin(T, C)
  delta <- as.numeric(T <= C)
  return(data.frame(T = T, C = C, Y = Y, delta = delta))
}


# ----------------------
# Single sensitivity (MSE) using estimator
# ----------------------

single_sensitivity_analysis <- function(n, true_theta, weibull_T, weibull_C,
                                        theta_range = NULL,
                                        family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (is.null(theta_range)) {
    if (family == "clayton") theta_range <- seq(0.1, 10, length.out = 20)
    if (family == "gumbel") theta_range <- seq(1.1, 10, length.out = 20)
  }
  
  sim_df <- simulate_data(n, weibull_T, weibull_C, true_theta, family)
  
# sorted unique times
time_grid <- sort(unique(sim_df$Y))
  true_surv <- pweibull(time_grid, shape = weibull_T$shape,
                        scale = weibull_T$scale, lower.tail = FALSE)
  
  mse_results <- data.frame(theta = theta_range, mse = NA_real_)
  
  for (i in seq_along(theta_range)) {
    theta_test <- theta_range[i]
    # compute RW with test theta and the same family
    rw_fit <- tryCatch(
      rw_estimator(sim_df$Y, sim_df$delta, theta = theta_test, family = family),
      error = function(e) return(NULL)
    )
    if (is.null(rw_fit) || nrow(rw_fit) == 0) {
      mse_results$mse[i] <- NA_real_
      next
    }
    # linear interpolation (survival defined for times not in event grid)
    rw_surv_interp <- approx(rw_fit$time, rw_fit$surv, xout = time_grid,
                             yleft = 1, yright = tail(rw_fit$surv, 1))$y
    mse_results$mse[i] <- mean((rw_surv_interp - true_surv)^2, na.rm = TRUE)
  }
  
  return(list(mse_results = mse_results,
              censoring_rate = 1 - mean(sim_df$delta)))
}


theta_to_tau <- function(theta, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (family == "clayton") {
    return(theta / (theta + 2))  # Clayton
  } else {
    return(1 - 1/theta)          # Gumbel
  }
}

tau_to_theta <- function(tau, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (family == "clayton") {
    return(2 * tau / (1 - tau))  # Inverse: 
  } else {
    return(1 / (1 - tau))        # Inverse:
  }
}


# ----------------------
# Kendall's Tau-based Sensitivity Analysis
# ----------------------

single_sensitivity_analysis_tau <- function(n, true_tau, weibull_T, weibull_C,
                                           tau_range = NULL,
                                           family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  
  # Convert true tau to theta for data generation
  true_theta <- tau_to_theta(true_tau, family)
  
  if (is.null(tau_range)) {
    # Use meaningful tau range: 0.1 to 0.9 (weak to strong dependence)
    tau_range <- seq(0.01, 0.99, by = 0.01)
  }
  
  sim_df <- simulate_data(n, weibull_T, weibull_C, true_theta, family)
  
  # sorted unique times
  time_grid <- sort(unique(sim_df$Y))

  true_surv <- pweibull(time_grid, shape = weibull_T$shape,
                        scale = weibull_T$scale, lower.tail = FALSE)
  
  mse_results <- data.frame(tau = tau_range, mse = NA_real_)


  for (i in seq_along(tau_range)) {
    tau_test <- tau_range[i]
    # Convert test tau to theta for RW estimation
    theta_test <- tau_to_theta(tau_test, family)
    
    # compute RW with test theta
    rw_fit <- tryCatch(
      rw_estimator(sim_df$Y, sim_df$delta, theta = theta_test, family = family),
      error = function(e) return(NULL)
    )
    if (is.null(rw_fit) || nrow(rw_fit) == 0) {
      mse_results$mse[i] <- NA_real_
      next
    }
    
    rw_surv_interp <- approx(rw_fit$time, rw_fit$surv, xout = time_grid,
                             yleft = 1, yright = tail(rw_fit$surv, 1))$y
    mse_results$mse[i] <- mean((rw_surv_interp - true_surv)^2, na.rm = TRUE)
  }

 # Calculate KM MSE for baseline comparison
  km_fit <- survfit(Surv(Y, delta) ~ 1, data = sim_df)
  km_surv_interp <- approx(km_fit$time, km_fit$surv, xout = time_grid,
                          yleft = 1, yright = tail(km_fit$surv, 1))$y
  km_mse <- mean((km_surv_interp - true_surv)^2, na.rm = TRUE)
  
  return(list(mse_results = mse_results,
	      km_mse = km_mse,
              censoring_rate = 1 - mean(sim_df$delta),
              true_tau = true_tau,
              true_theta = true_theta))
}



# ----------------------
# Comparison curves 
# ----------------------

plot_comparison_curves <- function(n, tau, weibull_T, weibull_C, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  
  # Convert tau to theta for data generation and RW estimation
theta <- tau_to_theta(tau, family)
  
# Simulate data
 
  sim_df <- simulate_data(n, weibull_T, weibull_C, theta, family)
  
  # Kaplan-Meier estimator
  km_fit <- survfit(Surv(Y, delta) ~ 1, data = sim_df)
  
  # Rivest-Wells estimator with converted theta
  rw_fit_df <- rw_estimator(sim_df$Y, sim_df$delta, theta = theta, family = family)
  
  # True survival curve
  time_grid <- sort(unique(sim_df$Y))
  true_surv <- pweibull(time_grid, shape = weibull_T$shape,
                        scale = weibull_T$scale, lower.tail = FALSE)
  
  # Create comparison data frame
  km_df <- data.frame(time = km_fit$time, surv = km_fit$surv, method = "Kaplan-Meier")
  rw_df <- data.frame(time = rw_fit_df$time, surv = rw_fit_df$surv, method = "Rivest-Wells")
  true_df <- data.frame(time = time_grid, surv = true_surv, method = "True Survival")
  
  plot_df <- rbind(km_df, rw_df, true_df)
  
  censoring_rate <- 1 - mean(sim_df$delta)
  
  p <- ggplot(plot_df, aes(x = time, y = surv, color = method, linetype = method)) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    labs(title = paste("Survival Curve Comparison -", toupper(family), "Copula"),
         subtitle = paste("Tau =", tau, "(Theta =", round(theta, 2), ")", 
                          "| Censoring =", round(censoring_rate * 100, 1), "%"),
         x = "Time",
         y = "Survival Probability",
         color = "Method",
         linetype = "Method") +
    scale_color_manual(values = c("True Survival" = "black",
                                  "Kaplan-Meier" = "blue",
                                  "Rivest-Wells" = "darkgreen")) +
    scale_linetype_manual(values = c("True Survival" = "solid",
                                     "Kaplan-Meier" = "dashed",
                                     "Rivest-Wells" = "dotted")) +
    theme(legend.position = "bottom")
  
  return(p)
}

# ----------------------
# Run Kendall's Tau Sensitivity Analysis
# ----------------------

# =========================================================
# PARAMETERS
# =========================================================
set.seed(2025)
N <- c(200, 500)

tau_values <- c(0.1, 0.8)

weibull_params_T <- list(shape = 1.5, scale = 1)

weibull_params_C <- list(shape = 1.2, scale = 1)



# =========================================================
# 2x2 PANEL FUNCTION - SINGLE SHARED LEGEND
# =========================================================
make_panel_plots <- function(family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  plot_list <- list()
  k <- 1
  
  for (n in N) {
    for (tau in tau_values) {
      theta <- tau_to_theta(tau, family)
      
      # Simulate data
      sim_df <- simulate_data(n = n, 
                              weibull_T = weibull_params_T,
                              weibull_C = weibull_params_C, 
                              theta = theta, 
                              family = family)
      
      # KM estimator
      km_fit <- survfit(Surv(Y, delta) ~ 1, data = sim_df)
      km_df <- data.frame(time = km_fit$time, surv = km_fit$surv, 
                          method = "KM")
      
      # Copula-Graphic (RW) estimator
      rw_fit_df <- rw_estimator(sim_df$Y, sim_df$delta, 
                                theta = theta, family = family)
      rw_df <- data.frame(time = rw_fit_df$time, surv = rw_fit_df$surv, 
                          method = "CG")
      
      # True survival
      time_grid <- seq(0, max(sim_df$Y), length.out = 400)
      true_surv <- pweibull(time_grid, 
                            shape = weibull_params_T$shape,
                            scale = weibull_params_T$scale, 
                            lower.tail = FALSE)
      true_df <- data.frame(time = time_grid, surv = true_surv, 
                            method = "True Survival")
      
      plot_df <- rbind(km_df, rw_df, true_df)
      censoring_rate <- 1 - mean(sim_df$delta)
      
      # ================== PLOT WITHOUT LEGEND ==================
      p <- ggplot(plot_df, aes(x = time, y = surv, color = method, linetype = method)) +
        geom_step(data = subset(plot_df, method != "True Survival"), linewidth = 1.1) +
        geom_line(data = subset(plot_df, method == "True Survival"), linewidth = 1.1) +
        theme_minimal(base_size = 13) +
        labs(title = paste("n =", n, "| tau =", tau, 
                          "| censoring =", round(100 * censoring_rate, 1), "%"),
             x = "Time", y = "Survival Probability") +
        scale_color_manual(values = c("True Survival" = "black",
                                      "KM" = "blue",
                                      "CG" = "darkgreen")) +
        scale_linetype_manual(values = c("True Survival" = "solid",
                                         "KM" = "dashed",
                                         "CG" = "dotted")) +
        theme(legend.position = "none")   # ← Remove legend from individual plots
     
      plot_list[[k]] <- p
      k <- k + 1
    }
  }
  

  
  # ---------------------- Create 2x2 Grid ----------------------
  grid_plots <- grid.arrange(grobs = plot_list, ncol = 2, nrow = 2,
                             top = textGrob(paste("Survival curves under", 
                                                  tools::toTitleCase(family), "copula"),
                                            gp = gpar(fontsize = 16, fontface = "bold")))
  
  # ---------------------- Extract Shared Legend ----------------------
  # Take legend from the first plot
  legend_plot <- plot_list[[1]] + 
    theme(legend.position = "bottom") +
    guides(
      color = guide_legend(
        nrow = 1
      ),
      linetype = guide_legend(
        nrow = 1,
        override.aes = list(linewidth = 1.5)
      )
    )
  
  shared_legend <- get_legend(legend_plot)
  
  # Final panel
  final_panel <- grid.arrange(grid_plots, shared_legend, 
                              ncol = 1, 
                              heights = c(1, 0.10))   # slightly more space for legend
  
  #print(final_panel)
  return(invisible(final_panel))
}

# =========================================================
# Run for both families
# =========================================================
make_panel_plots("clayton")
make_panel_plots("gumbel")
