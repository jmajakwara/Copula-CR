# Required libraries
library(survival)
library(copula)
library(ggplot2)

# ----------------------
# Clayton / Gumbel RW Estimator (YOUR ORIGINAL)
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
# Data sampler (YOUR ORIGINAL)
# ----------------------
simulate_data <- function(n, weibull_T, weibull_C, theta, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (family == "clayton") {
    cop <- archmCopula(family = "clayton", param = theta, dim = 2)
  } else {
    cop <- archmCopula(family = "gumbel", param = theta, dim = 2)
  }
  U <- rCopula(n, cop)
  U1 <- U[, 1]
  U2 <- U[, 2]
  
  T <- weibull_T$scale * (-log(U1))^(1 / weibull_T$shape)
  C <- weibull_C$scale * (-log(U2))^(1 / weibull_C$shape)
  
  Y <- pmin(T, C)
  delta <- as.numeric(T <= C)
  return(data.frame(T = T, C = C, Y = Y, delta = delta))
}

# ----------------------
# Helper functions
# ----------------------
tau_to_theta <- function(tau, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (family == "clayton") {
    return(2 * tau / (1 - tau))
  } else {
    return(1 / (1 - tau))
  }
}

theta_to_tau <- function(theta, family = c("clayton", "gumbel")) {
  family <- match.arg(family)
  if (family == "clayton") {
    return(theta / (theta + 2))
  } else {
    return(1 - 1/theta)
  }
}

# =====================================================
# THE KEY FUNCTION WITH MISSPECIFICATION
# =====================================================
# NOTE: family_sim and family_fit are DIFFERENT for misspecification!
# =====================================================

single_misspec_analysis <- function(n, true_tau, weibull_T, weibull_C,
                                    tau_range = NULL,
                                    family_sim,   # <-- Family for DATA GENERATION
                                    family_fit) { # <-- Family for FITTING (DIFFERENT!)
  
  # Convert true tau to theta using SIMULATION family
  true_theta <- tau_to_theta(true_tau, family_sim)
  
  if (is.null(tau_range)) {
    tau_range <- seq(0.01, 0.99, by = 0.01)
  }
  
  # STEP 1: Simulate data using family_sim
  sim_df <- simulate_data(n, weibull_T, weibull_C, true_theta, family_sim)
  
  # True survival curve
  time_grid <- sort(unique(sim_df$Y))
  true_surv <- pweibull(time_grid, shape = weibull_T$shape,
                        scale = weibull_T$scale, lower.tail = FALSE)
  
  mse_results <- data.frame(tau = tau_range, mse = NA_real_)
  
  for (i in seq_along(tau_range)) {
    tau_test <- tau_range[i]
    # Convert test tau to theta using FITTING family
    theta_test <- tau_to_theta(tau_test, family_fit)
    
    # STEP 2: Fit using family_fit (DIFFERENT from simulation!)
    rw_fit <- tryCatch(
      rw_estimator(sim_df$Y, sim_df$delta, theta = theta_test, family = family_fit),
      error = function(e) return(NULL)
    )
    
    if (is.null(rw_fit) || nrow(rw_fit) == 0) next
    
    rw_surv_interp <- approx(rw_fit$time, rw_fit$surv, xout = time_grid,
                             yleft = 1, yright = tail(rw_fit$surv, 1))$y
    mse_results$mse[i] <- mean((rw_surv_interp - true_surv)^2, na.rm = TRUE)
  }
  
  # KM baseline
  km_fit <- survfit(Surv(Y, delta) ~ 1, data = sim_df)
  km_surv_interp <- approx(km_fit$time, km_fit$surv, xout = time_grid,
                           yleft = 1, yright = tail(km_fit$surv, 1))$y
  km_mse <- mean((km_surv_interp - true_surv)^2, na.rm = TRUE)
  
  return(list(mse_results = mse_results,
              km_mse = km_mse,
              censoring_rate = 1 - mean(sim_df$delta),
              true_tau = true_tau,
              true_theta_sim = true_theta))
}

# =====================================================
# Monte Carlo wrapper for misspecification
# =====================================================
mc_misspec_analysis <- function(num_mc = 1000, n, true_tau,
                                weibull_T, weibull_C,
                                tau_range = NULL,
                                family_sim,
                                family_fit) {
  
  if (is.null(tau_range)) {
    tau_range <- seq(0.01, 0.99, by = 0.01)
  }
  
  avg_mse_df <- data.frame(tau = tau_range, avg_mse = 0, num_valid = 0)
  censoring_rates <- numeric(num_mc)
  km_mse_values <- numeric(num_mc)
  
  for (mc in seq_len(num_mc)) {
    set.seed(123 + mc)
    sens_res <- single_misspec_analysis(n, true_tau, weibull_T, weibull_C,
                                        tau_range = tau_range,
                                        family_sim = family_sim,
                                        family_fit = family_fit)
    
    valid_idx <- !is.na(sens_res$mse_results$mse)
    avg_mse_df$avg_mse[valid_idx] <- avg_mse_df$avg_mse[valid_idx] + sens_res$mse_results$mse[valid_idx]
    avg_mse_df$num_valid[valid_idx] <- avg_mse_df$num_valid[valid_idx] + 1
    km_mse_values[mc] <- sens_res$km_mse
    censoring_rates[mc] <- sens_res$censoring_rate
    
    if (mc %% 100 == 0) message("MC iter: ", mc, " | Sim: ", family_sim, " | Fit: ", family_fit)
  }
  
  avg_mse_df$avg_mse <- ifelse(avg_mse_df$num_valid > 0,
                               avg_mse_df$avg_mse / avg_mse_df$num_valid, NA)
  
  best_tau <- avg_mse_df$tau[which.min(avg_mse_df$avg_mse)]
  best_theta <- tau_to_theta(best_tau, family_fit)
  
  return(list(avg_mse_results = avg_mse_df[, c("tau", "avg_mse")],
              best_tau = best_tau,
              best_theta = best_theta,
              min_avg_mse = min(avg_mse_df$avg_mse, na.rm = TRUE),
              avg_km_mse = mean(km_mse_values),
              true_tau = true_tau,
              true_theta_sim = tau_to_theta(true_tau, family_sim),
              family_sim = family_sim,
              family_fit = family_fit,
              avg_censoring_rate = mean(censoring_rates)))
}

# ----------------------
# Plotting function
# ----------------------
plot_misspec_results <- function(mc_results) {
  df <- mc_results$avg_mse_results
  
  scenario_label <- paste0("Simulate: ", toupper(mc_results$family_sim),
                           " | Fit: ", toupper(mc_results$family_fit))
  
  p <- ggplot(df, aes(x = tau, y = avg_mse)) +
    geom_line(linewidth = 1.2, color = "darkgreen") +
    geom_point(size = 2, color = "darkgreen") +
    geom_hline(yintercept = mc_results$avg_km_mse,
               color = "blue", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = mc_results$true_tau,
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = mc_results$best_tau,
               color = "darkgreen", linetype = "dashed", linewidth = 1) +
    theme_minimal() +
    labs(title = paste("MC Sensitivity Analysis - MISSPECIFICATION"),
         subtitle = paste(scenario_label,
                          "\nTrue tau =", round(mc_results$true_tau, 3),
                          "| Best tau =", round(mc_results$best_tau, 3),
                          "| KM MSE =", round(mc_results$avg_km_mse, 5),
                          "| Min MSE =", round(mc_results$min_avg_mse, 5),
                          "| Censoring =", round(mc_results$avg_censoring_rate * 100, 1), "%"),
         x = "Assumed Kendall's Tau",
         y = "Average MSE") +
    annotate("text", x = mc_results$true_tau, y = max(df$avg_mse, na.rm = TRUE) * 0.9,
             label = "True Tau", color = "red", hjust = -0.1, size = 3) +
    annotate("text", x = mc_results$best_tau, y = max(df$avg_mse, na.rm = TRUE) * 0.8,
             label = "Best Tau", color = "darkgreen", hjust = -0.1, size = 3)
  
  return(p)
}

# =====================================================
# RUN MISSPECIFICATION ANALYSES
# =====================================================

N <- 200
weibull_params_T <- list(shape = 1.5, scale = 1)
weibull_params_C <- list(shape = 1.2, scale = 1)
runs <- 1000
tau_values <- c(0.1, 0.5, 0.8)

# =====================================================
# SCENARIO A: Simulate Clayton, Fit Gumbel (MISSPECIFICATION)
# =====================================================
cat("\n========================================\n")
cat("MISSPECIFICATION: Simulate Clayton, Fit Gumbel\n")
cat("========================================\n")

# THE MISMATCH IS HERE: family_sim = "clayton", family_fit = "gumbel"
clayton_sim_gumbel_fit_low <- mc_misspec_analysis(runs, N, tau_values[1],
                                                  weibull_params_T, weibull_params_C,
                                                  family_sim = "clayton",
                                                  family_fit = "gumbel")

clayton_sim_gumbel_fit_med <- mc_misspec_analysis(runs, N, tau_values[2],
                                                  weibull_params_T, weibull_params_C,
                                                  family_sim = "clayton",
                                                  family_fit = "gumbel")

clayton_sim_gumbel_fit_high <- mc_misspec_analysis(runs, N, tau_values[3],
                                                   weibull_params_T, weibull_params_C,
                                                   family_sim = "clayton",
                                                   family_fit = "gumbel")

# =====================================================
# SCENARIO B: Simulate Gumbel, Fit Clayton (MISSPECIFICATION)
# =====================================================
cat("\n========================================\n")
cat("MISSPECIFICATION: Simulate Gumbel, Fit Clayton\n")
cat("========================================\n")

# THE MISMATCH IS HERE: family_sim = "gumbel", family_fit = "clayton"
gumbel_sim_clayton_fit_low <- mc_misspec_analysis(runs, N, tau_values[1],
                                                  weibull_params_T, weibull_params_C,
                                                  family_sim = "gumbel",
                                                  family_fit = "clayton")

gumbel_sim_clayton_fit_med <- mc_misspec_analysis(runs, N, tau_values[2],
                                                  weibull_params_T, weibull_params_C,
                                                  family_sim = "gumbel",
                                                  family_fit = "clayton")

gumbel_sim_clayton_fit_high <- mc_misspec_analysis(runs, N, tau_values[3],
                                                   weibull_params_T, weibull_params_C,
                                                   family_sim = "gumbel",
                                                   family_fit = "clayton")

# =====================================================
# SAVE PLOTS
# =====================================================

pdf(file="misspec_clayton_sim_gumbel_fit_low_200.pdf")
print(plot_misspec_results(clayton_sim_gumbel_fit_low))
dev.off()
pdf(file="misspec_clayton_sim_gumbel_fit_med_200.pdf")
print(plot_misspec_results(clayton_sim_gumbel_fit_med))
dev.off()
pdf(file="misspec_clayton_sim_gumbel_fit_high_200.pdf")
print(plot_misspec_results(clayton_sim_gumbel_fit_high))
dev.off()

pdf(file="misspec_gumbel_sim_clayton_fit_low_200.pdf")
print(plot_misspec_results(gumbel_sim_clayton_fit_low))
dev.off()
pdf(file="misspec_gumbel_sim_clayton_fit_med_200.pdf")
print(plot_misspec_results(gumbel_sim_clayton_fit_med))
dev.off()
pdf(file="misspec_gumbel_sim_clayton_fit_high_200.pdf")
print(plot_misspec_results(gumbel_sim_clayton_fit_high))
dev.off()

# =====================================================
# SUMMARY OUTPUT
# =====================================================

cat("\n\n========================================\n")
cat("MISSPECIFICATION ANALYSIS SUMMARY\n")
cat("========================================\n")

cat("\n--- Scenario A: Simulate Clayton, Fit Gumbel ---\n")
cat("True tau = 0.1: Best tau =", round(clayton_sim_gumbel_fit_low$best_tau, 3),
    "| Min MSE =", round(clayton_sim_gumbel_fit_low$min_avg_mse, 6),
    "| KM MSE =", round(clayton_sim_gumbel_fit_low$avg_km_mse, 6), "\n")
cat("True tau = 0.5: Best tau =", round(clayton_sim_gumbel_fit_med$best_tau, 3),
    "| Min MSE =", round(clayton_sim_gumbel_fit_med$min_avg_mse, 6),
    "| KM MSE =", round(clayton_sim_gumbel_fit_med$avg_km_mse, 6), "\n")
cat("True tau = 0.8: Best tau =", round(clayton_sim_gumbel_fit_high$best_tau, 3),
    "| Min MSE =", round(clayton_sim_gumbel_fit_high$min_avg_mse, 6),
    "| KM MSE =", round(clayton_sim_gumbel_fit_high$avg_km_mse, 6), "\n")

cat("\n--- Scenario B: Simulate Gumbel, Fit Clayton ---\n")
cat("True tau = 0.1: Best tau =", round(gumbel_sim_clayton_fit_low$best_tau, 3),
    "| Min MSE =", round(gumbel_sim_clayton_fit_low$min_avg_mse, 6),
    "| KM MSE =", round(gumbel_sim_clayton_fit_low$avg_km_mse, 6), "\n")
cat("True tau = 0.5: Best tau =", round(gumbel_sim_clayton_fit_med$best_tau, 3),
    "| Min MSE =", round(gumbel_sim_clayton_fit_med$min_avg_mse, 6),
    "| KM MSE =", round(gumbel_sim_clayton_fit_med$avg_km_mse, 6), "\n")
cat("True tau = 0.8: Best tau =", round(gumbel_sim_clayton_fit_high$best_tau, 3),
    "| Min MSE =", round(gumbel_sim_clayton_fit_high$min_avg_mse, 6),
    "| KM MSE =", round(gumbel_sim_clayton_fit_high$avg_km_mse, 6), "\n")

cat("\nAll misspecification results saved to PDF files.\n")