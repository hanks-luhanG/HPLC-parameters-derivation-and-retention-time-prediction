################################################################################
# Project: HPLC parameters derivation and retention time prediction
# Script:  Quadratic_Model_Fitting.R
# Author:  Luhan Guan (Shenzhen Chipscreen Biosciences / Chongqing University)
# Date:    2026-05-13
# Version: 1.0.0
# Description: This script performs non-linear least squares optimization 
#              based on the Quadratic model for RP-HPLC data.
# Contact: [glhmail2000@foxmail.com]
################################################################################

library(numDeriv)
library(dplyr)
library(parallel)
library(deSolve)
library(readxl)
library(writexl)

# =======================================================================
# ---- 1. Core function definition ----
# =======================================================================

# [1.1 training ODE]
ode_func <- function(t, L, parms) {
  t_prime <- t - parms$td
  if (t_prime <= 0) phi <- parms$phi0
  else if (t_prime <= parms$t_grad_end) phi <- parms$phi0 + parms$B * t_prime
  else phi <- parms$phi1
  
  exponent <- -parms$S1 * phi + parms$S2 * phi^2
  exp_term <- exp(pmin(exponent, 700))
  dLdt <- parms$u_linear / (1 + parms$kw * exp_term)
  return(list(c(L = dLdt)))
}

# [1.2 prediction ODE]
ode_func_prediction <- function(t, L, parms) {
  t_prime <- pmax(t - parms$td, 0)
  stage_idx <- findInterval(t_prime, parms$stage_boundaries, left.open = TRUE) + 1
  stage_idx <- pmin(stage_idx, length(parms$B_values))
  phi_eff <- parms$phi_start_of_stage[stage_idx] + parms$B_values[stage_idx] * (t_prime - parms$stage_start_times[stage_idx])
  
  exponent <- -parms$S1 * phi_eff + parms$S2 * phi_eff^2
  exp_term <- exp(pmin(exponent, 700))
  dLdt <- parms$u_linear / (1 + parms$kw * exp_term)
  return(list(c(L = dLdt)))
}

# [1.3 Integral calculation]
stable_integrate <- function(par, model_type, single_row_data) {
  lnkw <- par[1]; S1 <- par[2]; S2 <- if (model_type == "three_par") par[3] else 0
  
  parms <- list(kw = exp(lnkw), S1 = S1, S2 = S2, 
                phi0 = single_row_data$phi0, B = single_row_data$B, 
                td = single_row_data$td, u_linear = single_row_data$u_linear,
                phi1 = single_row_data$phi1, t_grad_end = single_row_data$t_grad_end)
  
  t_end <- as.numeric(single_row_data$t_R)
  ode_res <- tryCatch({
    deSolve::ode(y = c(L = 0), times = c(0, t_end), func = ode_func, 
                 parms = parms, method = "lsoda", atol = 1e-8, rtol = 1e-8)
  }, error = function(e) NA)
  
  if (is.matrix(ode_res)) return(ode_res[nrow(ode_res), "L"]) else return(NA)
}

# [1.4 2-param fitting]
objective <- function(par, model_type, fit_data) {
  lnkw <- par[1]; S1 <- par[2]; S2 <- if (model_type == "three_par") par[3] else 0
  if (lnkw < -15 || lnkw > 15 || S1 < 0.01 || S1 > 30 ||
      (model_type == "three_par" && (S2 < -30 || S2 > 30))) return(1e9)
  
  err <- 0
  for(j in 1:nrow(fit_data)) {
    current_row <- as.list(fit_data[j, , drop = FALSE])
    integral <- stable_integrate(par, model_type, current_row)
    if(is.na(integral)) return(1e6)
    err <- err + (integral - fit_data$L[j])^2
  }
  return(err)
}

# [1.5 3-param optimization]
objective_stabilized <- function(par, lambda_reg, fit_data, lss_lnkw, lss_S1) {
  err <- objective(par, "three_par", fit_data)
  if (err >= 1e6) return(err)
  penalty <- lambda_reg * ((par[1] - lss_lnkw)^2 + (par[2] - lss_S1)^2 + 100 * par[3]^2)
  return(err + penalty)
}

# [1.6 Prediction function]
predict_t_R_robust_new <- function(lnkw, S1, S2, target_L, Q_flow, td_delay, phi0_start, phi_values_end, stage_boundaries) {
  A_cm2 <- pi * (4.6 / (2 * 10))^2
  u_linear <- Q_flow / A_cm2
  num_stages <- length(phi_values_end)
  stage_start_times <- c(0, stage_boundaries[1:(num_stages - 1)])
  phi_start_of_stage <- c(phi0_start, if(num_stages > 1) phi_values_end[1:(num_stages - 1)])
  B_values <- (phi_values_end - phi_start_of_stage) / (stage_boundaries - stage_start_times)
  
  parms_pred <- list(kw = exp(lnkw), S1 = S1, S2 = S2, u_linear = u_linear, td = td_delay,
                     B_values = B_values, phi_start_of_stage = phi_start_of_stage,
                     stage_boundaries = stage_boundaries, stage_start_times = stage_start_times)
  
  obj_func <- function(t_R) {
    if (t_R <= 0) return(target_L)
    res <- deSolve::ode(y = c(L = 0), times = c(0, t_R), func = ode_func_prediction, parms = parms_pred)
    return(res[nrow(res), "L"] - target_L)
  }
  return(tryCatch(uniroot(obj_func, lower = 0.1, upper = 100, extendInt = "upX")$root, error = function(e) NA))
}

# =======================================================================
# ---- 2. batch Processing ----
# =======================================================================

path_to_file <- "./data/example_data.xlsx"
raw_data_input <- read_excel(path_to_file, sheet = "acetonitrile-ACE AR")

# function starting
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, library(deSolve))
clusterExport(cl, c("ode_func", "stable_integrate", "objective", "objective_stabilized", "ode_func_prediction"))

final_output_list <- list()

for (i in 1:nrow(raw_data_input)) {
  row_val <- raw_data_input[i, ]
  cat(sprintf("\n[%d/%d] Processing: %s\n", i, nrow(raw_data_input), as.character(row_val[[1]])))
  
  # Construct a standalone dataframe for compound
  current_fit_data <- data.frame(
    t_R = as.numeric(c(row_val[[3]], row_val[[4]], row_val[[5]])),
    phi0 = 0.1, phi1 = 1.0, t_grad_end = c(25, 30, 40),
    L = 15, Q = 1, td = 0.822, A = pi * (4.6 / 2)^2
  )
  current_fit_data$B <- (current_fit_data$phi1 - current_fit_data$phi0) / current_fit_data$t_grad_end
  current_fit_data$u_linear <- current_fit_data$Q / (current_fit_data$A / 100)

  # [step 1: LSS model fitting（2-param）]
  results_lss <- parLapply(cl, 1:10, function(idx, f_data) {
    set.seed(123 + idx)
    optim(par = c(runif(1, 5, 12), runif(1, 1, 20)), 
          fn = objective, model_type = "two_par", fit_data = f_data,
          method = "L-BFGS-B", lower = c(-15, 0.01), upper = c(15, 30))
  }, f_data = current_fit_data)
  
  best_lss <- results_lss[[which.min(sapply(results_lss, `[[`, "value"))]]
  l_opt <- best_lss$par[1]; s_opt <- best_lss$par[2]

  # [step 2: Quadratic model optimization（3-param）]
  results_quad <- parLapply(cl, 1:10, function(idx, f_data, l_ref, s_ref) {
    set.seed(123 + idx)
    optim(par = c(l_ref, s_ref, runif(1, -2, 2)), 
          fn = objective_stabilized, lambda_reg = 1e-4, 
          fit_data = f_data, lss_lnkw = l_ref, lss_S1 = s_ref,
          method = "L-BFGS-B", lower = c(-15, 0.01, -30), upper = c(15, 30, 30))
  }, f_data = current_fit_data, l_ref = l_opt, s_ref = s_opt)
  
  best_quad <- results_quad[[which.min(sapply(results_quad, `[[`, "value"))]]

  # [step 3: prediction]
  pred_tR <- predict_t_R_robust_new(
    lnkw = best_quad$par[1], S1 = best_quad$par[2], S2 = best_quad$par[3],
    target_L = 15, Q_flow = 1, td_delay = 0.822, phi0_start = 0.1,
    phi_values_end = c(1.0, 1.0), stage_boundaries = c(35, 45)  
  )

  # [step 4:storage]
  final_output_list[[i]] <- data.frame(
    Name = as.character(row_val[[1]]),
    lnkw = best_quad$par[1], S1 = best_quad$par[2], S2 = best_quad$par[3],
    SSE = best_quad$value, Predicted_tR = pred_tR
  )
}

stopCluster(cl)
write_xlsx(bind_rows(final_output_list), "./Quadratic_Model_results.xlsx")
cat("\n✅ Fully calculated！")