# Required libraries
library(survival)
library(copula)
library(ggplot2)

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

# ----------------------
# Reliable copula-based sampler using 'copula' package
# ----------------------
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
# Single sensitivity (MSE) using corrected estimator
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
# Kendall's Tau MC Sensitivity Analysis
# ----------------------

mc_sensitivity_analysis_tau <- function(num_mc = 1000, n, true_tau,
                                       weibull_T, weibull_C,
                                       tau_range = NULL,
                                       family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  
  if (is.null(tau_range)) {
    tau_range <- seq(0.01, 0.99, by = 0.01)
  }
  
  avg_mse_df <- data.frame(tau = tau_range, avg_mse = 0, num_valid = 0)
  censoring_rates <- numeric(num_mc)
  km_mse_values <- numeric(num_mc)
  
  for (mc in seq_len(num_mc)) {
    set.seed(123 + mc)
    sens_res <- single_sensitivity_analysis_tau(n, true_tau, weibull_T, weibull_C,
                                               tau_range = tau_range, family = family)
    valid_idx <- !is.na(sens_res$mse_results$mse)
    avg_mse_df$avg_mse[valid_idx] <- avg_mse_df$avg_mse[valid_idx] + sens_res$mse_results$mse[valid_idx]
    avg_mse_df$num_valid[valid_idx] <- avg_mse_df$num_valid[valid_idx] + 1
    km_mse_values[mc] <- sens_res$km_mse
    censoring_rates[mc] <- sens_res$censoring_rate
    if (mc %% 50 == 0) message("Completed MC iter: ", mc)
  }
  
  avg_mse_df$avg_mse <- ifelse(avg_mse_df$num_valid > 0, 
                               avg_mse_df$avg_mse / avg_mse_df$num_valid, NA)
  
  best_tau <- avg_mse_df$tau[which.min(avg_mse_df$avg_mse)]
  best_theta <- tau_to_theta(best_tau, family)
  
  return(list(avg_mse_results = avg_mse_df[, c("tau", "avg_mse")],
              best_tau = best_tau,
              best_theta = best_theta,
              min_avg_mse = min(avg_mse_df$avg_mse, na.rm = TRUE),
	      avg_km_mse = mean(km_mse_values),
              true_tau = true_tau,
              true_theta = sens_res$true_theta,
              copula_family = family,
              avg_censoring_rate = mean(censoring_rates)))
}

# ----------------------
# Kendall's Tau Plotting
# ----------------------

plot_sensitivity_results_tau <- function(mc_results) {
  df <- mc_results$avg_mse_results
  p <- ggplot(df, aes(x = tau, y = avg_mse)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_vline(xintercept = mc_results$true_tau,
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = mc_results$best_tau,
               color = "darkgreen", linetype = "dashed", linewidth = 1) +
    theme_minimal() +
    labs(title = paste("MC Sensitivity Analysis -", toupper(mc_results$copula_family), "Copula"),
         subtitle = paste("True tau =", round(mc_results$true_tau, 3), 
                          "(theta =", round(mc_results$true_theta, 3), ")",
                          "| Best tau =", round(mc_results$best_tau, 3),
                          "| Avg Censoring =", round(mc_results$avg_censoring_rate * 100, 1), "%"),
         x = "Kendall's Tau",
         y = "Average MSE") +
    annotate("text", x = mc_results$true_tau, y = max(df$avg_mse, na.rm = TRUE) * 0.9,
             label = "True Tau", color = "red", hjust = -0.1, size = 3) +
    annotate("text", x = mc_results$best_tau, y = max(df$avg_mse, na.rm = TRUE) * 0.8,
             label = "Best Tau", color = "darkgreen", hjust = -0.1, size = 3)
  return(p)
}

# ----------------------
# Comparison curves 
# ----------------------

plot_comparison_curves <- function(n, tau, weibull_T, weibull_C, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  
  # Convert tau to theta for data generation and RW estimation
  theta <- tau_to_theta(tau, family)
  
  # Simulate data
  set.seed(2025)
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

# Define parameters
N <- 200
weibull_params_T <- list(shape = 1.5, scale = 1)
weibull_params_C <- list(shape = 1.2, scale = 1)
runs <- 1000

# Use comparable tau values across families
tau_values <- c(0.1, 0.5, 0.8)  # Weak, moderate, strong dependence

# Clayton Copula MC Sensitivity Analysis
cat("Running Clayton and Gumbel Copulas MC Sensitivity Analysis (Kendall's Tau)...\n")
cat("Running Clayton Copula MC Sensitivity Analysis (Low Tau)...\n")
clayton_tau_low <- mc_sensitivity_analysis_tau(runs, N, tau_values[1], 
                                              weibull_params_T, weibull_params_C, 
                                              family = "clayton")
cat("Running Clayton Copula MC Sensitivity Analysis (Med Tau)...\n")
clayton_tau_med <- mc_sensitivity_analysis_tau(runs, N, tau_values[2], 
                                              weibull_params_T, weibull_params_C, 
                                              family = "clayton")
cat("Running Clayton Copula MC Sensitivity Analysis (high Tau)...\n")
clayton_tau_high <- mc_sensitivity_analysis_tau(runs, N, tau_values[3], 
                                               weibull_params_T, weibull_params_C, 
                                               family = "clayton")

# Gumbel Copula MC Sensitivity Analysis
cat("Running Gumbel Copula MC Sensitivity Analysis (Low Tau)...\n")
gumbel_tau_low <- mc_sensitivity_analysis_tau(runs, N, tau_values[1], 
                                             weibull_params_T, weibull_params_C, 
                                             family = "gumbel")
cat("Running Gumbel Copula MC Sensitivity Analysis (Med Tau)...\n")
gumbel_tau_med <- mc_sensitivity_analysis_tau(runs, N, tau_values[2], 
                                             weibull_params_T, weibull_params_C, 
                                             family = "gumbel")
cat("Running Gumbel Copula MC Sensitivity Analysis (High Tau)...\n")
gumbel_tau_high <- mc_sensitivity_analysis_tau(runs, N, tau_values[3], 
                                              weibull_params_T, weibull_params_C, 
                                              family = "gumbel")

# Plot results
pdf(file="clayton_low_200.pdf")
print(plot_sensitivity_results_tau(clayton_tau_low))
dev.off()
pdf(file="clayton_med_200.pdf")
print(plot_sensitivity_results_tau(clayton_tau_med))
dev.off()
pdf(file="clayton_high_200.pdf")
print(plot_sensitivity_results_tau(clayton_tau_high))
dev.off()
pdf(file="gumbel_low_200.pdf")
print(plot_sensitivity_results_tau(gumbel_tau_low))
dev.off()
pdf(file="gumbel_med_200.pdf")
print(plot_sensitivity_results_tau(gumbel_tau_med))
dev.off()
pdf(file="gumbel_high_200.pdf")
print(plot_sensitivity_results_tau(gumbel_tau_high))
dev.off()

# Plot comparison curves for selected scenarios (single run examples)
pdf(file="clayton_200.pdf")
print(plot_comparison_curves(N, 0.8, weibull_params_T, weibull_params_C, "clayton"))
dev.off()
pdf(file="gumbel_200.pdf")
print(plot_comparison_curves(N, 0.8, weibull_params_T, weibull_params_C, "gumbel"))
dev.off()

cat("\n=== CLAYTON COPULA MC SENSITIVITY ANALYSIS SUMMARY ===\n")
cat("Number of MC runs per scenario:", runs, "\n")
cat("True tau = 0.1: KM MSE =", round(clayton_tau_low$avg_km_mse, 4), "\n")
cat("True tau = 0.5: KM MSE =", round(clayton_tau_med$avg_km_mse, 4), "\n")
cat("True tau = 0.8: KM MSE =", round(clayton_tau_high$avg_km_mse, 4), "\n")
cat("True tau = 0.1: Best tau =", clayton_tau_low$best_tau, "Avg Min MSE =", round(clayton_tau_low$min_avg_mse, 4), "\n")
cat("True tau = 0.5: Best tau =", clayton_tau_med$best_tau, "Avg Min MSE =", round(clayton_tau_med$min_avg_mse, 4), "\n")
cat("True tau = 0.8: Best tau =", clayton_tau_high$best_tau, "Avg Min MSE =", round(clayton_tau_high$min_avg_mse, 4), "\n")

cat("\n=== GUMBEL COPULA MC SENSITIVITY ANALYSIS SUMMARY ===\n")
cat("Number of MC runs per scenario:", runs, "\n")
cat("True tau = 0.1: KM MSE =", round(gumbel_tau_low$avg_km_mse, 4), "\n")
cat("True tau = 0.5: KM MSE =", round(gumbel_tau_med$avg_km_mse, 4), "\n")
cat("True tau = 0.8: KM MSE =", round(gumbel_tau_high$avg_km_mse, 4), "\n")
cat("True tau = 0.1: Best tau =", gumbel_tau_low$best_tau, "Avg Min MSE =", round(gumbel_tau_low$min_avg_mse, 4), "\n")
cat("True tau = 0.5: Best tau =", gumbel_tau_med$best_tau, "Avg Min MSE =", round(gumbel_tau_med$min_avg_mse, 4), "\n")
cat("True tau = 0.8: Best tau =", gumbel_tau_high$best_tau, "Avg Min MSE =", round(gumbel_tau_high$min_avg_mse, 4), "\n")