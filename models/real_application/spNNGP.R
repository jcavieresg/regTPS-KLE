library(spNNGP)
library(loo)
library(abind)

# ============================================================
# Reading the ata
# ============================================================
y_sqrt <- sqrt(sp_data$y_obs)
hist(y_sqrt)
coords <- as.matrix(sp_data[, c("s1", "s2")])

# Same prediction grid as SPDE and regTPS-KLE
grid_total  <- regTPS_KLE_tmb$grid_total
grid_matrix <- as.matrix(grid_total)


# ============================================================
# Starting values and priors for spNNGP
# ============================================================
sigma_sq_start <- var(y_sqrt) * 0.8
tau_sq_start   <- var(y_sqrt) * 0.2

phi_start <- 3 / (0.3 * diff(range(coords[, 1])))

starting <- list(phi = phi_start, sigma.sq = sigma_sq_start, tau.sq   = tau_sq_start, nu = 1.0)
tuning <- list(phi = 0.5, sigma.sq = 0.5, tau.sq = 0.5, nu = 1.0)

priors <- list(phi.Unif = c(3 / (2 * diff(range(coords[, 1]))), 3 / (0.01 * diff(range(coords[, 1])))),
               sigma.sq.IG = c(2, sigma_sq_start),
               tau.sq.IG   = c(2, tau_sq_start),
               nu.Unif     = c(0.1, 2.0))


# ============================================================
# 3. Fit 3 independent chains
# ============================================================
n_chains  <- 3
n_samples <- 3000

spnngp_fits <- vector("list", n_chains)
spnngp_times <- numeric(n_chains)

for (chain in seq_len(n_chains)) {
  
  cat("\n========================================\n")
  cat("Running spNNGP chain", chain, "of", n_chains, "\n")
  cat("========================================\n")
  
  set.seed(1000 + chain)
  t0 <- Sys.time()
  spnngp_fits[[chain]] <- spNNGP(
    y_sqrt ~ 1,
    coords        = coords,
    method        = "latent",
    family        = "gaussian",
    n.neighbors   = 15,
    starting      = starting,
    tuning        = tuning,
    priors        = priors,
    cov.model     = "matern",
    n.samples     = n_samples,
    n.omp.threads = no_cores,
    n.report      = 500)
  t1 <- Sys.time()
  spnngp_times[chain] <- as.numeric(difftime(t1, t0, units = "min"))
  cat("Chain", chain, "fit time (min):", round(spnngp_times[chain], 2), "\n")
}


cat("\nTotal spNNGP fitting time (min):", round(sum(spnngp_times), 2), "\n")
cat("Mean chain fitting time (min):", round(mean(spnngp_times), 2), "\n")


saveRDS(spnngp_fits, file='spnngp_fits.RDS')

# ============================================================
# Apply the SAME burn-in trimming used in your prediction
# script (burn_in = 1000 per chain) before computing pointwise log-lik,
# ============================================================

compute_spnngp_loglik <- function(spnngp_fit, X, y, burn_in = 0) {
  
  w_samples     <- spnngp_fit$p.w.samples
  beta_samples  <- spnngp_fit$p.beta.samples
  theta_samples <- spnngp_fit$p.theta.samples
  
  n_draws_total <- ncol(w_samples)
  keep_idx <- seq.int(from = min(burn_in + 1, n_draws_total), to = n_draws_total)
  
  w_samples     <- w_samples[, keep_idx, drop = FALSE]
  beta_samples  <- beta_samples[keep_idx, , drop = FALSE]
  theta_samples <- theta_samples[keep_idx, , drop = FALSE]
  
  n_draws <- ncol(w_samples)
  n_obs   <- length(y)
  
  tau_sq_col <- grep("^tau\\.sq$|tau\\.sq", colnames(theta_samples))
  if (length(tau_sq_col) != 1) stop("Could not uniquely identify tau.sq in theta samples.")
  tau_sq_samples <- theta_samples[, tau_sq_col]
  
  log_lik <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (d in seq_len(n_draws)) {
    mu_d <- as.numeric(X %*% beta_samples[d, ]) + w_samples[, d]
    sigma_d <- sqrt(tau_sq_samples[d])
    log_lik[d, ] <- dnorm(y, mean = mu_d, sd = sigma_d, log = TRUE)
  }
  log_lik
}

# ============================================================
# Recompute with burn-in removed, matching your prediction script
# ============================================================

burn_in <- 700   # same value used in the prediction/mapping script

n_obs <- length(y_sqrt)
X_full <- matrix(1, nrow = n_obs, ncol = 1)
log_lik_chains <- vector("list", n_chains)

for (chain in seq_len(n_chains)) {
  cat("Computing log-likelihood for chain", chain, "(post burn-in)\n")
  log_lik_chains[[chain]] <- compute_spnngp_loglik(spnngp_fits[[chain]], X = X_full, y = y_sqrt,
                                                   burn_in = burn_in)
  cat("  Dimensions:", paste(dim(log_lik_chains[[chain]]), collapse = " x "), "\n")
}

log_lik_3_matrix <- rbind(log_lik_chains[[1]], log_lik_chains[[2]], log_lik_chains[[3]])
cat("Total post-burn-in draws used for LOO:", nrow(log_lik_3_matrix), "\n")

loo_3 <- loo(log_lik_3_matrix, cores = no_cores)
print(loo_3)

# ============================================================
# Also check Pareto-k diagnostics for ALL THREE models -- strong
# spatial correlation can stress PSIS-LOO's importance-sampling
# assumptions; worth confirming none of the three show a
# disproportionate number of problematic k values.
# ============================================================
cat("\nSPDE Pareto-k summary:\n");       print(pareto_k_table(loo_1))
cat("\nregTPS-KLE Pareto-k summary:\n"); print(pareto_k_table(loo_2))
cat("\nspNNGP Pareto-k summary:\n");     print(pareto_k_table(loo_3))

cmp <- loo_compare(loo_1, loo_2, loo_3)
rownames(cmp) <- c("SPDE","regTPS-KLE","spNNGP")[match(rownames(cmp), c("model1","model2","model3"))]
print(cmp)








# ============================================================
# Bulk ESS for spNNGP posterior samples
# ============================================================
library(posterior)

# ============================================================
# 1. Extract post-burn-in theta samples from the 3 chains
# ============================================================

burn_in <- 700

theta_list <- lapply(spnngp_fits, function(fit) {
  theta <- fit$p.theta.samples
  keep <- (burn_in + 1):nrow(theta)
  theta[keep, , drop = FALSE]
})


# Check dimensions
lapply(theta_list, dim)

# Check parameter names
colnames(theta_list[[1]])


# ============================================================
# 2. Bulk ESS for each hyperparameter
# ============================================================

theta_names <- colnames(theta_list[[1]])

ess_theta <- numeric(length(theta_names))
names(ess_theta) <- theta_names

for (p in seq_along(theta_names)) {
  # Combine the same parameter across chains
  x <- lapply(theta_list, function(z) z[, p])
  # Convert to draws_array:
  # iterations x chains x variables
  x_array <- array(unlist(x),dim = c(nrow(theta_list[[1]]), n_chains, 1))
  dimnames(x_array)[[3]] <- theta_names[p]
  x_draws <- posterior::as_draws_array(x_array)
  ess_theta[p] <- posterior::ess_bulk(x_draws)[1]
}

ess_theta




# ============================================================
# 3. Median Bulk ESS
# ============================================================
median_bulk_ess <- median(ess_theta, na.rm = TRUE)
cat("Median Bulk ESS =", median_bulk_ess, "\n")




# ============================================================
# 4. Sampling efficiency
# ============================================================
sampling_time_min <- sum(spnngp_times)
ess_per_min <- median_bulk_ess / sampling_time_min
log_ess_per_min <- log(ess_per_min)
cat("Total sampling time (min):", sampling_time_min, "\n")
cat("Median Bulk ESS:", median_bulk_ess, "\n")
cat("Median Bulk ESS / min:", ess_per_min, "\n")
cat("log(Median Bulk ESS / min):", log_ess_per_min, "\n")






# ============================================================
# spNNGP: posterior prediction + prediction time
# ============================================================
library(spNNGP)

# ------------------------------------------------------------
# 1. Prediction locations
# ------------------------------------------------------------
grid_total  <- regTPS_KLE_tmb$grid_total
grid_matrix <- as.matrix(grid_total)
n_grid <- nrow(grid_matrix)

X_pred <- matrix(1, nrow = n_grid, ncol = 1)
cat("Number of prediction locations:", n_grid, "\n")


# ------------------------------------------------------------
# 2. Burn-in
# ------------------------------------------------------------
burn_in <- 700
cat("Burn-in per chain:", burn_in, "\n")


# ============================================================
# 3. Prediction for each chain
# ============================================================

pred_spnngp_list <- vector("list", n_chains)
prediction_times <- numeric(n_chains)

for (chain in seq_len(n_chains)) {
  cat("\n========================================\n")
  cat("Predicting spNNGP chain", chain, "of", n_chains, "\n")
  cat("========================================\n")
  # Start timing ONLY the prediction step
  t0 <- Sys.time()
  pred_spnngp_list[[chain]] <- predict(spnngp_fits[[chain]], 
                                       X.0 = X_pred, coords.0 = grid_matrix,
                                       n.omp.threads = no_cores)
  t1 <- Sys.time()
  prediction_times[chain] <- as.numeric(difftime(t1, t0, units = "mins"))
  cat("Chain", chain, "prediction time (min):", round(prediction_times[chain], 4), "\n")
}


# ============================================================
# 4. Total and mean prediction time
# ============================================================
total_prediction_time <- sum(prediction_times)
mean_prediction_time <- mean(prediction_times)
sd_prediction_time <- sd(prediction_times)


cat("\n========================================\n")
cat("spNNGP prediction time summary\n")
cat("========================================\n")

cat("Chain prediction times (min):", paste(round(prediction_times, 4), collapse = ", "), "\n")
cat("Total prediction time (min):", round(total_prediction_time, 4), "\n")
cat("Mean prediction time per chain (min):", round(mean_prediction_time, 4), "\n")
cat("SD prediction time per chain (min):", round(sd_prediction_time, 4), "\n")


# ============================================================
# 5. Extract posterior predictive samples
# ============================================================
pred_draws_list <- lapply(pred_spnngp_list,
  function(pred) {
    keep_idx <- seq(burn_in + 1, ncol(pred$p.y.0))
    pred$p.y.0[, keep_idx, drop = FALSE]
  }
)


# ============================================================
# 6. Combine chains
# ============================================================
combined_p_y_0 <- do.call(cbind, pred_draws_list)

cat("\n========================================\n")
cat("Posterior predictive samples\n")
cat("========================================\n")

cat("Prediction locations:", nrow(combined_p_y_0), "\n")
cat("Posterior predictive draws:", ncol(combined_p_y_0), "\n")


# ============================================================
# 7. Posterior predictive quantiles
# ============================================================
y_summary_spnngp <- t(apply(combined_p_y_0, 1, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE))
colnames(y_summary_spnngp) <- c("q025", "median", "q975")


# ============================================================
# 8. Posterior predictive mean
# ============================================================
y_mean_spnngp <- rowMeans(combined_p_y_0, na.rm = TRUE)


# ============================================================
# 9. Final prediction data frame
# ============================================================
df_spnngp_prediction <- data.frame(s1 = grid_total[, 1],
                                   s2 = grid_total[, 2],
                                   mean = y_mean_spnngp,
                                   q025 = y_summary_spnngp[, "q025"],
                                   median = y_summary_spnngp[, "median"],
                                   q975 = y_summary_spnngp[, "q975"])


# ============================================================
# 10. Prediction-time summary table
# ============================================================
prediction_time_spnngp <- data.frame(
  Model = "spNNGP",
  Prediction_time_total_min = total_prediction_time,
  Prediction_time_mean_min = mean_prediction_time,
  Prediction_time_sd_min = sd_prediction_time,
  n_prediction_locations = n_grid,
  n_chains = n_chains,
  n_posterior_draws = ncol(combined_p_y_0))

print(prediction_time_spnngp)


# ============================================================
# 11. Save results
# ============================================================

saveRDS(list(predictions = df_spnngp_prediction,
             posterior_draws = combined_p_y_0,
             prediction_times = prediction_times,
             prediction_time_summary = prediction_time_spnngp),
  file = "outputs/spnngp_predictions_timing.RDS")
