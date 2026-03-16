#==================================================
# Breast Cancer Survival Analysis: Full Workflow + CG & Copula-Cox by Tau
#==================================================

#setwd("C:/Users/User/Downloads/")
#-------------------------
# 1) Working directory & data
#-------------------------
my_data <- read.csv("Breast Cancer Data.csv")

#-------------------------
# 2) Packages
#-------------------------
library(survival)
library(survminer)
library(dplyr)
library(knitr)
library(psych)
library(compound.Cox)
library(ggplot2)
library(gridExtra)
library(joint.Cox)
library(GGally)


#-------------------------
# 3) Quick data check
#-------------------------
cols_needed <- c("event_death","follow_up_duration","HBB","ARID3A","PKM2", "CCNA2", "G6PD")
print(colSums(is.na(my_data[cols_needed])))
print(head(my_data))
print(summary(my_data[cols_needed]))
print(psych::describe(my_data[cols_needed]))
#pairs(my_data[, c("HBB","ARID3A","PKM2", "CCNA2", "G6PD")],col=7,pch=1)
ggpairs(my_data[, c("HBB","ARID3A","PKM2", "CCNA2", "G6PD")], diag=list(continuous="barDiag"),
        upper = list(continuous = wrap(ggally_cor, method = "kendall")))+
  labs(x="log_2(gene expression)",y="log_2(gene expression)")

#-------------------------
# 4) Survival object
#-------------------------
surv_object <- Surv(time = my_data$follow_up_duration, event = my_data$event_death)

#-------------------------
# 5) Kaplan-Meier (overall)
#-------------------------
fit_overall <- survfit(surv_object ~ 1)
ggsurvplot(fit_overall, data = my_data, risk.table = TRUE,
           title = "Kaplan-Meier: Overall Survival")

#-------------------------
# 6) Median-based High/Low grouping
#-------------------------
genes <- c("HBB","ARID3A","PKM2", "CCNA2", "G6PD")
grouping <- function(gene_name){
  v <- my_data[[gene_name]]
  ifelse(v >= median(v, na.rm = TRUE), "High", "Low")
}
for (g in genes) my_data[[paste0(g,"_group")]] <- grouping(g)

#-------------------------
# 7) KM by gene group
#-------------------------
for (g in genes) {
  fit_g <- survfit(surv_object ~ my_data[[paste0(g,"_group")]])
  print(fit_g)
  print(
    ggsurvplot(fit_g, data = my_data, pval = TRUE, pval.method = TRUE,
               conf.int = TRUE, risk.table = TRUE,
               title = paste("Kaplan-Meier:", g, "High vs Low"),
               legend.labs = c("High","Low"))
  )
}

#==================================================
# 8) Copula-Graphic (CG) curves per copula × gene × τ
#==================================================
cg_results <- data.frame()

# Convert Kendall's tau to copula parameter (theta)
tau_to_theta <- function(copula_name, tau_value){
  if (copula_name %in% c("Frank", "Frank_pos", "Frank_neg")) copula_name <- "Frank"
  switch(copula_name,
         "Clayton" = {
           if (tau_value < 0) NA else 2 * tau_value / (1 - tau_value)
         },
         "Gumbel" = {
           if (tau_value < 0) NA else tau_value / (1 - tau_value)
         },
         "Frank" = {
           if (tau_value == 0) 0 else {
             tau_from_theta <- function(theta){
               integral_value <- integrate(function(t) t / (exp(t) - 1), 0, theta)$value
               1 - (4 / theta) * (1 - (integral_value / theta))
             }
             uniroot(function(theta) tau_from_theta(theta) - tau_value, c(-50, 50))$root
           }
         }
  )
}

# --- Function to plot Copula-Graphic curve ---
plot_CG <- function(copula_fun, copula_name, tau_value, time_vec, event_vec, group_vec, gene_name) {
  theta <- tau_to_theta(copula_name, tau_value)
  if (is.na(theta) || is.infinite(theta)) stop("Invalid theta for ", copula_name, " tau=", tau_value)
  
  t1 <- time_vec[group_vec == "Low"];  d1 <- event_vec[group_vec == "Low"]
  t2 <- time_vec[group_vec == "High"]; d2 <- event_vec[group_vec == "High"]
  tau_cut <- min(tapply(time_vec, group_vec, max), na.rm = TRUE)
  
  CG1 <- copula_fun(t1, d1, alpha = theta, S.plot = FALSE)
  CG2 <- copula_fun(t2, d2, alpha = theta, S.plot = FALSE)
  CG_test <- CG.test(time_vec, event_vec, group_vec,
                     copula = copula_fun, alpha = theta, S.plot = FALSE)
  
  res <- data.frame(
    Gene = gene_name,
    Copula = copula_name,
    Tau = tau_value,
    Theta = round(theta, 2),
    D = round(CG_test$test["Survival.diff"], 4),
    Pvalue = round(CG_test$test["P.value"], 4)
  )
  assign("cg_results", rbind(get("cg_results", envir = .GlobalEnv), res), envir = .GlobalEnv)
  
  p_title <- paste0(copula_name, " | θ=", round(theta,3),
                    " | τ=", round(tau_value,2),
                    " | p=", round(CG_test$test["P.value"],4))
  
  df <- data.frame(
    time = c(sort(t1), sort(t2)),
    surv = c(CG1$surv, CG2$surv),
    group = rep(c("Low", "High"), c(length(CG1$surv), length(CG2$surv)))
  )
  df$group_label <- paste(gene_name, df$group)   # “PKM2 Low” / “PKM2 High”
  
  ggplot(df, aes(x = time, y = surv, color = group_label)) +
    geom_step(size = 0.8) +
    geom_vline(xintercept = tau_cut, linetype = "dotted", linewidth = 0.8) +
    coord_cartesian(xlim = c(0, max(c(t1, t2), na.rm = TRUE)), ylim = c(0, 1)) +
    labs(title = p_title, x = "Time (years)", y = "Survival Probability", color = NULL) +
    scale_color_manual(values = c("blue", "red")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = c(0.2, 0.2)
    )
}

# --- Define copulas ---
copulas <- list(
  "Clayton"   = CG.Clayton,
  "Gumbel"    = CG.Gumbel,
  "Frank_pos" = CG.Frank,
  "Frank_neg" = CG.Frank
)

# --- Tau grids ---
cg_tau_grid <- list(
  Clayton   = c(0, 0.2, 0.5, 0.8),
  Gumbel    = c(0, 0.2, 0.5, 0.8),
  Frank_pos = c(0, 0.2, 0.5, 0.8),
  Frank_neg = c(0, -0.2, -0.5, -0.8)
)


# Copula-Graphic (CG) curves per copula × gene × τ

plots <- list()

for (copula_name in names(copulas)) {
  copula_fun <- copulas[[copula_name]]
  
  for (gene in genes) {
    group_var <- my_data[[paste0(gene, "_group")]]
    subset_plots <- list()  # temporary list per gene × copula
    
    # --- generate plots for all τ ---
    for (tau in cg_tau_grid[[copula_name]]) {
      plot_id <- paste(copula_name, gene, tau, sep = "_")
      message("Generating plot: ", plot_id)
      
      p <- tryCatch({
        plot_CG(
          copula_fun, copula_name, tau,
          my_data$follow_up_duration, my_data$event_death,
          group_var, gene
        )
      }, error = function(e) {
        message("Skipping ", plot_id, " (", e$message, ")")
        NULL
      })
      
      if (!is.null(p)) {
        subset_plots[[plot_id]] <- p
        message(plot_id, " done")  # show the plot is done
      }
    }
    
    # Arrange plots per gene × copula
    subset_plots <- Filter(Negate(is.null), subset_plots)
    
    if (length(subset_plots) > 0) {
      # dynamically choose layout based on how many plots were created
      ncol_val <- ifelse(length(subset_plots) >= 4, 2, length(subset_plots))
      nrow_val <- ceiling(length(subset_plots) / ncol_val)
      
      grid.arrange(
        grobs = subset_plots,
        ncol = ncol_val,
        nrow = nrow_val,
        top = paste(gene, "-", copula_name, "(Copula-Graphic by τ)")
      )
    } else {
      message("No plots generated for ", gene, " under ", copula_name)
    }
    
    # store subset in master list
    plots[[paste0(gene, "_", copula_name)]] <- subset_plots
    message("Plot for ", gene, " under ", copula_name, " done")
  }
}


#==================================================
# 9) Standard Cox regression (uni & multi)
#==================================================
cfit_HBB    <- coxph(surv_object ~ HBB, data = my_data)
cfit_ARID3A <- coxph(surv_object ~ ARID3A, data = my_data)
cfit_PKM2   <- coxph(surv_object ~ PKM2, data = my_data)
cfit_CCNA2  <- coxph(surv_object ~ CCNA2, data = my_data)
cfit_G6PD   <- coxph(surv_object ~ G6PD, data = my_data)
cfit_multi  <- coxph(surv_object ~ HBB + ARID3A + PKM2 + CCNA2 + G6PD, data = my_data)

cat("===== Standard Cox Regression =====\n")
print(summary(cfit_HBB))
print(summary(cfit_ARID3A))
print(summary(cfit_PKM2))
print(summary(cfit_CCNA2))
print(summary(cfit_G6PD))
print(summary(cfit_multi))

#==================================================
# 10) Copula-Cox regression (per-gene HR + joint p-value)
#==================================================

# Define gene sets
genes_all <- genes                      # all genes in your data
genes_sel <- c("ARID3A", "CCNA2", "G6PD")  # selected 3 genes for PI_select

# Per-gene Cox regression for all genes
copula_summary <- data.frame()
for (gene in genes_all){
  fit <- coxph(as.formula(paste0("Surv(follow_up_duration, event_death) ~ ", gene)), 
               data = my_data)
  temp <- data.frame(
    Gene = gene,
    HR = summary(fit)$coef[,2],
    CI_lower = summary(fit)$conf.int[,3],
    CI_upper = summary(fit)$conf.int[,4],
    P_value = summary(fit)$coef[,5]
  )
  copula_summary <- rbind(copula_summary, temp)
}
cat("===== Per-gene Cox Regression (all genes) =====\n")
print(copula_summary)

# Per-gene Cox regression for selection genes
copula_summary <- data.frame()
for (gene in genes_sel){
  fit <- coxph(as.formula(paste0("Surv(follow_up_duration, event_death) ~ ", gene)), 
               data = my_data)
  temp <- data.frame(
    Gene = gene,
    HR = summary(fit)$coef[,2],
    CI_lower = summary(fit)$conf.int[,3],
    CI_upper = summary(fit)$conf.int[,4],
    P_value = summary(fit)$coef[,5]
  )
  copula_summary <- rbind(copula_summary, temp)
}
cat("===== Per-gene Cox Regression (selected genes) =====\n")
print(copula_summary)


#==================================================
# 11) Prognostic Index 1 (PI_genes_all) via all genes
#==================================================

X1.mat <- as.matrix(my_data[, genes_all])
PI_genes_all_value <- uni.score(my_data$follow_up_duration, my_data$event_death, X1.mat)
beta1.vec <- as.numeric(PI_genes_all_value$beta)

my_data$PI_genes_all <- as.vector(X1.mat %*% beta1.vec)
my_data$PI_genes_all_group <- ifelse(my_data$PI_genes_all >= median(my_data$PI_genes_all, na.rm = TRUE), "High", "Low")

# KM Plot for PI_genes_all
fit_PI_genes_all <- survfit(surv_object ~ PI_genes_all_group, data = my_data)
print(ggsurvplot(fit_PI_genes_all, data = my_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                 title = "Kaplan-Meier: PI (All Genes)"))

cfit_PI_genes_all <- coxph(surv_object ~ PI_genes_all, data = my_data)
print(summary(cfit_PI_genes_all))


#==================================================
# 12) Prognostic Index 2 (PI_select) via selected 3 genes (ARID3A, CCNA2, G6PD)
#==================================================

X2.mat <- as.matrix(my_data[, genes_sel])
PI_select_value <- uni.score(my_data$follow_up_duration, my_data$event_death, X2.mat)
beta2.vec <- as.numeric(PI_select_value$beta)

my_data$PI_select <- as.vector(X2.mat %*% beta2.vec)
my_data$PI_select_group <- ifelse(my_data$PI_select >= median(my_data$PI_select, na.rm = TRUE), "High", "Low")

# KM Plot for PI_select
fit_PI_select <- survfit(surv_object ~ PI_select_group, data = my_data)
print(ggsurvplot(fit_PI_select, data = my_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                 title = "Kaplan-Meier: PI (ARID3A + CCNA2 + G6PD)"))

cfit_PI_select <- coxph(surv_object ~ PI_select, data = my_data)
print(summary(cfit_PI_select))


#==================================================
# 13) Copula-Graphic (CG) Curves for PI_genes_all and PI_select
#==================================================

plots_PI <- list()
for (copula_name in names(copulas)) {
  copula_fun <- copulas[[copula_name]]
  
  for (PI_name in c("PI_genes_all", "PI_select")) {
    group_var <- my_data[[paste0(PI_name, "_group")]]
    copula_plots <- list()
    
    for (tau in cg_tau_grid[[copula_name]]) {
      plot_id <- paste(copula_name, PI_name, tau, sep = "_")
      message("Generating plot: ", plot_id)
      
      p <- tryCatch({
        plot_CG(
          copula_fun, copula_name, tau,
          my_data$follow_up_duration, my_data$event_death,
          group_var, PI_name
        )
      }, error = function(e) {
        message("Skipping ", plot_id, " (", e$message, ")")
        NULL
      })
      
      if (!is.null(p)) {
        copula_plots[[plot_id]] <- p
        plots_PI[[plot_id]] <- p
      }
    }
    
    valid_plots <- Filter(Negate(is.null), copula_plots)
    if (length(valid_plots) > 0) {
      grid.arrange(
        grobs = valid_plots,
        ncol = 2,
        top = paste(PI_name, "-", copula_name, "(Copula-Graphic by τ)")
      )
    }
  }
}

#==================================================
# 14) Optional – Combine PI_genes_all and PI_select Summary Table
#==================================================

PI_summary <- data.frame(
  PI = c("PI_genes_all (All Genes)", "PI_select (ARID3A + CCNA2 + G6PD)"),
  HR = c(summary(cfit_PI_genes_all)$coef[2], summary(cfit_PI_select)$coef[2]),
  CI_lower = c(summary(cfit_PI_genes_all)$conf.int[3], summary(cfit_PI_select)$conf.int[3]),
  CI_upper = c(summary(cfit_PI_genes_all)$conf.int[4], summary(cfit_PI_select)$conf.int[4]),
  P_value = c(summary(cfit_PI_genes_all)$coef[5], summary(cfit_PI_select)$coef[5])
)

cat("===== Prognostic Index Summary =====\n")
print(PI_summary)

#-------------------------------------------------------
# Show the PI formulas for PI_genes_all and PI_select
#-------------------------------------------------------

# Formula for PI_genes_all
cat("\n===== Formula: PI_genes_all =====\n")
PI_genes_all_formula <- paste0(
  "PI_genes_all = ",
  paste(
    sprintf("%.4f * %s", beta1.vec, genes_all),
    collapse = " + "
  )
)
cat(PI_genes_all_formula, "\n")

# Formula for PI_select (3 genes: ARID3A, CCNA2, G6PD)
cat("\n===== Formula: PI_select =====\n")
PI_select_formula <- paste0(
  "PI_select = ",
  paste(
    sprintf("%.4f * %s", beta2.vec, genes_sel),
    collapse = " + "
  )
)
cat(PI_select_formula, "\n")

#==================================================================
# 15) Copula-Dependent (Clayton) HR Estimates for All Genes + PI Models
#==================================================================
Cdepend_list <- list()
tau_1 <- 0.1 / (2 + 0.1)

for (gene in c(genes, "PI_genes_all", "PI_select")) {
  Cdepend_list[[gene]] <- list()
  
  if (length(unique(na.omit(my_data[[gene]]))) < 2) next  # skip constant vars
  
  for (tau in c(cg_tau_grid[["Clayton"]], 0.6, tau_1)) {
    res <- tryCatch(
      dependCox.reg(my_data$follow_up_duration, my_data$event_death, my_data[[gene]], tau,
                    var = TRUE, censor.reg = TRUE, baseline = TRUE),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      Cdepend_list[[gene]][[as.character(tau)]] <- list(
        surv.reg   = res$surv.reg,
        censor.reg = res$censor.reg
      )
    }
  }
}

print(Cdepend_list)


#==================================================
# 16) Copula-Graphic Test and HR Summary Tables
#==================================================
cg_results <- cg_results %>% mutate(Pvalue = as.numeric(Pvalue), D = as.numeric(D))

cat("===== Copula-Graphic Test Results (All) =====\n")
print(kable(cg_results, row.names = FALSE, digits = 3, caption = "CG Test Results: Difference Metric D & P-value"))

Cdepend_tbl <- purrr::map_dfr(names(Cdepend_list), function(gene) {
  purrr::map_dfr(names(Cdepend_list[[gene]]), function(tau) {
    res <- Cdepend_list[[gene]][[tau]]$surv.reg
    beta <- as.numeric(res["beta"])
    
    tibble::tibble(
      Gene   = gene,
      Tau    = as.numeric(tau),
      HR     = exp(beta),
      Pvalue = as.numeric(res["P"])
    )
  })
})

cat("===== Hazard Ratio (HR) Estimates Results (Clayton copula) =====\n")
print(kable(Cdepend_tbl, digits = 4,
            caption = "Hazard (HR) Estimates Results: HR & P-value under Clayton copula"))