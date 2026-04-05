#==================================================
# Breast Cancer Survival Analysis (ARID3A ONLY)
#==================================================


library(survival)
library(survminer)
library(dplyr)
library(knitr)
library(psych)
library(compound.Cox)
library(ggplot2)
library(gridExtra)
library(joint.Cox)

my_data <- read.csv("Breast Cancer Data.csv")

#-------------------------
#  Quick data check
#-------------------------
cols_needed <- c("event_death","follow_up_duration","ARID3A")

print(colSums(is.na(my_data[cols_needed])))
print(psych::describe(my_data[cols_needed]))

surv_object <- Surv(time = my_data$follow_up_duration,
                    event = my_data$event_death)

#-------------------------
# Kaplan-Meier (overall)
#-------------------------
fit_overall <- survfit(surv_object ~ 1)

ggsurvplot(fit_overall, data = my_data,
           risk.table = TRUE,
           title = "Kaplan-Meier: Overall Survival")

#-------------------------
# Median-based grouping (ARID3A ONLY)
#-------------------------
my_data$ARID3A_group <- ifelse(
  my_data$ARID3A >= median(my_data$ARID3A, na.rm = TRUE),
  "High", "Low"
)

#-------------------------
# KM by ARID3A group
#-------------------------
fit_ARID3A <- survfit(surv_object ~ ARID3A_group, data = my_data)

print(fit_ARID3A)

ggsurvplot(fit_ARID3A,
           data = my_data,
           pval = TRUE,
           pval.method = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           title = "Kaplan-Meier: ARID3A High vs Low",
           legend.labs = c("High","Low"))

#==================================================
# Copula-Graphic (CG) Analysis
#==================================================

cg_results <- data.frame()

#-------------------------
# Tau to Theta
#-------------------------
tau_to_theta <- function(copula_name, tau_value){

  switch(copula_name,

    "Clayton" = ifelse(tau_value < 0, NA,
                       2*tau_value/(1 - tau_value)),

    "Gumbel"  = ifelse(tau_value < 0, NA,
                       tau_value/(1 - tau_value)) # due to reparameterisation of tau=alpha/(alpha+1) in the compound.Cox package
  )
}

#-------------------------
# Copulas
#-------------------------
copulas <- list(
  Clayton = CG.Clayton,
  Gumbel  = CG.Gumbel
)

tau_grid <- list(
  Clayton = c(0,0.2,0.5,0.8), #c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), 
  Gumbel  = c(0,0.2,0.5,0.8)  #c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) 
)

#-------------------------
# CG  Plot Function 
#-------------------------
plot_CG <- function(t, d, PI, copula_fun, copula_name, tau){

  theta <- tau_to_theta(copula_name, tau)
  if (is.na(theta) || is.infinite(theta)) return(NULL)

  idx <- complete.cases(t, d, PI)
  t <- t[idx]; d <- d[idx]; PI <- PI[idx]

  cutoff <- median(PI)

  #-------------------------
  # CG Test
  #-------------------------
  CG_res <- CG.test(
    t.vec = t,
    d.vec = d,
    PI    = PI,
    cutoff = cutoff,
    alpha = theta,
    copula = copula_fun,
    S.plot = FALSE
  )

  D_val <- as.numeric(CG_res$test["Survival.diff"])
  P_val <- as.numeric(CG_res$test["P.value"])

  result_row <- data.frame(
    Copula = copula_name,
    Tau = tau,
    Theta = round(theta,3),
    D = round(D_val,4),
    Pvalue = round(P_val,4)
  )

  #==================================================
  # CG STRATIFIED CURVES
  #==================================================
  idx_low  <- PI < cutoff
  idx_high <- PI >= cutoff

  CG_low <- copula_fun(t[idx_low], d[idx_low],
                       alpha = theta, S.plot = FALSE)

  CG_high <- copula_fun(t[idx_high], d[idx_high],
                        alpha = theta, S.plot = FALSE)

  df_cg <- rbind(
    data.frame(time = CG_low$time,
               surv = CG_low$surv,
               group = "Low"),
    data.frame(time = CG_high$time,
               surv = CG_high$surv,
               group = "High")
  )

  group <- factor(ifelse(idx_low, "Low (KM)", "High (KM)"))



  #==================================================
  # PLOT
  #==================================================
  p <- ggplot() +
    geom_step(data = df_cg,
              aes(time, surv, color = group),
              linewidth = 0.9) +
    labs(
      title = paste0(
        copula_name,
        " (tau=", tau,
        ", theta=", round(theta,2),
        ", p=", round(P_val,3), ")"
      ),
      x = "Time",
      y = "Survival",
      color = "Group"
    ) +
    theme_minimal()

  return(list(plot = p, result = result_row))
}

plots <- list()

for (copula_name in names(copulas)) {

  subset_plots <- list()

  for (tau in tau_grid[[copula_name]]) {

    res <- plot_CG(
      my_data$follow_up_duration,
      my_data$event_death,
      my_data$ARID3A,   
      copulas[[copula_name]],
      copula_name,
      tau
    )

    if (!is.null(res)) {
      subset_plots[[paste0("tau_",tau)]] <- res$plot
      cg_results <- rbind(cg_results, res$result)
    }
  }

  grid.arrange(grobs = subset_plots,
               ncol = 2,
               top = paste("CG", copula_name))
}

cg_results <- cg_results %>%
  distinct(Copula, Tau, .keep_all = TRUE)

cat("===== CG Results ===== \n")
print(kable(cg_results, digits=3))

