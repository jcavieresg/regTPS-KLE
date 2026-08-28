rm(list = ls())
options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, reshape2,
               patchwork, purrr, kableExtra)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading TMB models
tmb_spde <- readRDS('fits_TMB_spde.RDS')
tmb_tps <- readRDS('fits_TMB_tps.RDS')

#=================
# TPS/KLE
y_obs_tps   <- tmb_tps[[1]]$tmb_data$y_obs
field_tps_sp   <- tmb_tps[[1]]$u_true_sp
field_tps_grid   <- tmb_tps[[1]]$u_true_grid


# Extract
y_obs_spde  <- tmb_spde[[1]]$tmb_data$y
field_spde_sp  <- tmb_spde[[1]]$u_true_sp
field_spde_grid  <- tmb_spde[[1]]$u_true_grid

#----------------------------------------------------------
# 5) Compare results
cat("TPS vs SPDE  (y_obs):        ", all.equal(y_obs_tps, y_obs_spde), "\n")
cat("TPS vs SPDE  (field_sp):        ", all.equal(field_tps_sp, field_spde_sp), "\n")
cat("TPS vs SPDE  (field_grid):        ", all.equal(field_tps_grid, field_spde_grid), "\n")


# MCMC models (Stan)
mcmc_spde1 <- readRDS('stan_spde_1.RDS')
mcmc_spde2 <- readRDS('stan_spde_2.RDS')
mcmc_spde3 <- readRDS('stan_spde_3.RDS')
mcmc_spde4 <- readRDS('stan_spde_4.RDS')

mcmc_tps1 <- readRDS('stan_tps_1.RDS')
mcmc_tps2 <- readRDS('stan_tps_2.RDS')
mcmc_tps3 <- readRDS('stan_tps_3.RDS')
mcmc_tps4 <- readRDS('stan_tps_4.RDS')



# #===================================================
# # Extracting the posteriors for the latent field
# # SPDE
# #===================================================
extract_spde_field <- function(stan_fit, A_grid) {
  post <- rstan::extract(stan_fit)

  u_draws       <- post$u_tilde         # iter x n_spde
  rho_draws     <- post$rho       # iter
  sigma_u_draws <- post$sigma_u   # iter

  n_iter <- nrow(u_draws)
  n_grid <- nrow(A_grid)

  # ---- Storage ----
  field_mean <- numeric(n_grid)
  for (iter in seq_len(n_iter)) {
   # Transform hyperparameters
    kappa_iter <- sqrt(8) / rho_draws[iter]
    tau_iter   <- 1.0 / (kappa_iter * sigma_u_draws[iter])
   # Reconstruct centered field
    u_iter <- u_draws[iter, ] / tau_iter
    field_mean <- field_mean + as.vector(A_grid %*% u_iter)
  }
  field_mean <- field_mean / n_iter
  return(list(field_grid_mean = field_mean))
}


spde_field_grid <- vector("list", 4)

for (i in 1:4) {
  spde_field_grid[[i]] <- extract_spde_field( stan_fit = get(paste0("mcmc_spde", i)), A_grid   = tmb_spde[[i]][[4]]$A_grid)$field_grid_mean
}

#===================================================
# Extracting the posteriors for the latent field
# TPS
#===================================================
extract_tps_field <- function(stan_fit, Phi_kle_grid, S_diag_truncated, M_P_null_space) {

  z_tilde_draws <- rstan::extract(stan_fit)$z_tilde   # iter x M
  alpha_draws  <- exp(rstan::extract(stan_fit)$logalpha)
  n_iter <- nrow(z_tilde_draws)
  M      <- ncol(z_tilde_draws)
  n_grid <- nrow(Phi_kle_grid)

  field_mean <- numeric(n_grid)
  for (iter in seq_len(n_iter)) {

    alpha_iter <- alpha_draws[iter]

    # Scaling vector
    scale_iter <- numeric(M)
    scale_iter[1:M_P_null_space] <- 1.0
    scale_iter[(M_P_null_space + 1):M] <- 1.0 / sqrt(1.0 + alpha_iter * S_diag_truncated[(M_P_null_space + 1):M])
    # Reconstruct Z
    Z_iter <- scale_iter * z_tilde_draws[iter, ]
    # Accumulate posterior mean of field
    field_mean <- field_mean + as.vector(Phi_kle_grid %*% Z_iter)
}

  field_mean <- field_mean / n_iter
  return(list(field_grid_mean = field_mean))
}



tps_field_grid <- vector("list", 4)

for (i in 1:4) {
  tps_field_grid[[i]] <- extract_tps_field(stan_fit = get(paste0("mcmc_tps", i)),
                                           Phi_kle_grid    = tmb_tps[[i]][[4]]$Phi_kle_grid,
                                           S_diag_truncated= tmb_tps[[i]][[4]]$S_diag_truncated,
                                           M_P_null_space  = tmb_tps[[i]][[4]]$M_P_null_space)$field_grid_mean
}



#=====================================================================================
# Plots of the paper
#=====================================================================================
alpha_post1 <- exp(rstan::extract(mcmc_tps1)$logalpha)  # Posterior samples of KL coefficients
alpha_post2 <- exp(rstan::extract(mcmc_tps2)$logalpha)  # Posterior samples of KL coefficients
alpha_post3 <- exp(rstan::extract(mcmc_tps3)$logalpha)  # Posterior samples of KL coefficients
alpha_post4 <- exp(rstan::extract(mcmc_tps4)$logalpha)  # Posterior samples of KL coefficients

mcmc_tps_list <- list("Sce.1" = mcmc_tps1,
                      "Sce.2" = mcmc_tps2,
                      "Sce.3" = mcmc_tps3,
                      "Sce.4" = mcmc_tps4)

# Extract logalpha and then exponentiate to get alpha.
alpha_posteriors <- lapply(mcmc_tps_list, function(x) {
  exp(rstan::extract(x)$logalpha)
})

# Combine all posterior samples into a single data frame for ggplot.
# The `bind_rows` function will create a column named '.' when it combines unnamed vectors.
# We can rename that column directly inside the `bind_rows` call.
plot_df <- bind_rows(lapply(alpha_posteriors, function(x) {
  data.frame(alpha = x)
}), .id = "Scenarios")

#==================================================
#               Plots
#==================================================
# Prior on alpha
# Log-Cauchy density function
dlogcauchy <- function(x, location, scale) {
  1 / (pi * scale * (1 + ((log(x) - location) / scale)^2)) * (1 / x)
}

# Prior parameters for log(alpha)
mu0 <- 0     # location in log-space
sigma0 <- 1  # scale in log-space

# Create prior density data for each scenario
alpha_range <- seq(min(plot_df$alpha), max(plot_df$alpha), length.out = 500)

prior_df <- do.call(rbind, lapply(names(mcmc_tps_list), function(scen) {
  data.frame(alpha = alpha_range,
             density = dlogcauchy(alpha_range, location = mu0, scale = sigma0), Scenarios = scen, type = "Prior")}))

posterior_df <- plot_df %>%
  group_by(Scenarios) %>%
  do({dens <- density(.$alpha)
  data.frame(alpha = dens$x, density = dens$y)}) %>%
  ungroup() %>%
  mutate(type = "Posterior")


combined_df <- bind_rows(posterior_df, prior_df)

#==================================================
# Add posterior quantiles (20% and 80%)
#==================================================
quantiles_df <- plot_df %>%
  group_by(Scenarios) %>%
  summarise(q20 = quantile(alpha, 0.20), q80 = quantile(alpha, 0.80), .groups = "drop")

#==================================================
# Plot with quantile markers
#==================================================
plot1 <- ggplot(combined_df, aes(x = alpha, y = density, color = type, linetype = type)) +
  geom_line(size = 1) +
  facet_wrap(~Scenarios, scales = "free") +
  # Add vertical lines for quantiles
  geom_vline(data = quantiles_df, aes(xintercept = q20),
             linetype = "dashed", color = "black", size = 0.4, inherit.aes = FALSE) +
  geom_vline(data = quantiles_df, aes(xintercept = q80),
             linetype = "dashed", color = "black", size = 0.4, inherit.aes = FALSE) +
  labs(title = expression(bold("Posterior and Prior Distributions of") ~ alpha),
                          x = expression(alpha ~ "values"),y = "Density", color = "", linetype = "") +
  scale_color_manual(values = c("Prior" = "red", "Posterior" = "blue")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "top", 
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))

plot1



#===================================================
# Extract Posterior Means for Z Coefficients
# (Non-centered parameterization)
#===================================================
extract_z_posterior <- function(mcmc_fit, S_diag_truncated, M_P_null_space) {
  # Extract MCMC samples
  post <- rstan::extract(mcmc_fit)
  
  z_tilde_draws <- post$z_tilde  # iter x M (whitened)
  alpha_draws   <- exp(post$logalpha)
  
  n_iter <- nrow(z_tilde_draws)
  M      <- ncol(z_tilde_draws)
  
  # Storage for transformed Z
  Z_samples <- matrix(0, n_iter, M)
  
  # Transform each MCMC sample
  for (iter in 1:n_iter) {
    alpha_iter <- alpha_draws[iter]
    
    # NON-CENTERED TRANSFORMATION
    z_iter <- numeric(M)
    
    # Null space (unpenalized)
    if (M_P_null_space > 0) {
      z_iter[1:M_P_null_space] <- z_tilde_draws[iter, 1:M_P_null_space]
    }
    
    # Penalized modes
    if (M > M_P_null_space) {
      idx_penalized <- (M_P_null_space + 1):M
      S_k <- S_diag_truncated[idx_penalized]
      
      scale_factor <- 1.0 / sqrt(1.0 + alpha_iter * S_k + 1e-10)
      z_iter[idx_penalized] <- z_tilde_draws[iter, idx_penalized] * scale_factor
    }
    
    Z_samples[iter, ] <- z_iter
  }
  
  # Return posterior mean and samples
  return(list(mean = colMeans(Z_samples), samples = Z_samples, alpha_mean = mean(alpha_draws)))
}

#===================================================
# Extract for all scenarios
#===================================================
z_results <- list()

for (i in 1:4) {
  z_results[[i]] <- extract_z_posterior(mcmc_fit = get(paste0("mcmc_tps", i)),
    S_diag_truncated = tmb_tps[[i]][[4]]$S_diag_truncated,
    M_P_null_space   = tmb_tps[[i]][[4]]$M_P_null_space)
}

#===================================================
# Compute Prior SDs using posterior mean of alpha
#===================================================
compute_prior_sd <- function(alpha_mean, S_diag_truncated, M_P_null_space) {
  M <- length(S_diag_truncated)
  prior_sd <- sqrt(1 / (1 + alpha_mean * S_diag_truncated + 1e-10))
  # Null space modes have prior SD = 1
  if (M_P_null_space > 0) {
    prior_sd[1:M_P_null_space] <- 1
  }
  return(prior_sd)
}

prior_sds <- list()
for (i in 1:4) {
  prior_sds[[i]] <- compute_prior_sd(alpha_mean = z_results[[i]]$alpha_mean,
    S_diag_truncated = tmb_tps[[i]][[4]]$S_diag_truncated,
    M_P_null_space   = tmb_tps[[i]][[4]]$M_P_null_space)
}

#===================================================
# Create plotting dataframe
#===================================================
plot_data <- bind_rows(lapply(1:4, function(i) {
    M_P_null <- tmb_tps[[i]][[4]]$M_P_null_space
    M <- length(prior_sds[[i]])
    data.frame(PriorSD = prior_sds[[i]],
      z_post_abs = abs(z_results[[i]]$mean),
      z_post = z_results[[i]]$mean,  # Keep sign for interpretation
      ModeIndex = 1:M,
      ModeType = ifelse(1:M <= M_P_null, "Null Space", "Penalized"),
      Scenario = paste0("Sce. ", i))
  })
)

#===================================================
# Enhanced Diagnostic Plot
#===================================================
plot2 <- ggplot(plot_data, aes(x = PriorSD, y = z_post_abs)) +
  geom_ribbon(aes(ymin = 0, ymax = 2 * PriorSD),  fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 2, color = "red", linetype = "dotted", linewidth = 0.6) +
  geom_point(aes(color = ModeType, shape = ModeType), size = 2.5, alpha = 0.7) +
  facet_wrap(~Scenario, scales = "free", ncol = 2) +
  labs(title = "Posterior Mean Z vs Prior Standard Deviation",
    x = "Prior Standard Deviation (SD)", y = "Posterior Mean of Z", color = "Mode Type", shape = "Mode Type") +
  scale_color_manual(values = c("Null Space" = "#E41A1C", "Penalized" = "#377EB8")) +
  scale_shape_manual(values = c("Null Space" = 16, "Penalized" = 17)) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

print(plot2)



#===================================================
# Enhanced Diagnostic Plot - LOG-LOG SCALE
#===================================================
plot2 <- ggplot(plot_data, aes(x = PriorSD, y = z_post_abs)) +
  geom_ribbon(aes(ymin = 0, ymax = 2 * PriorSD),  fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1,  color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 2,  color = "red", linetype = "dotted", linewidth = 0.6) +
  geom_point(aes(color = ModeType, shape = ModeType), size = 2.5, alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  # Facet by scenario
  facet_wrap(~Scenario, scales = "free", ncol = 2) +
  labs(title = bquote(bold("Posterior Mean Z vs Prior SD (" * log[10] * " scale)")), 
       x = "Prior Standard Deviation (SD)", y = "Posterior Mean of |Z|", color = "Mode Type", shape = "Mode Type") +
  scale_color_manual(values = c("Null Space" = "#E41A1C", "Penalized" = "#377EB8")) +
  scale_shape_manual(values = c("Null Space" = 16, "Penalized" = 17)) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank())

#===============================================
# PLOT3
#===============================================
compute_kle_eigenvalues <- function(scenario_idx, tmb_tps, z_results) {
  v_k <- tmb_tps[[scenario_idx]][[4]]$S_diag_truncated
  alpha_mean <- z_results[[scenario_idx]]$alpha_mean
  lambda_k <- 1 / (1 + alpha_mean * v_k + 1e-10)
  M_P_null <- tmb_tps[[scenario_idx]][[4]]$M_P_null_space
  if (M_P_null > 0) {
    lambda_k[1:M_P_null] <- 1  # Null space has λ = 1
  }
  M_truncation <- tmb_tps[[scenario_idx]][[6]]
  return(list(v_k = v_k, lambda_k = lambda_k, alpha = alpha_mean, M_truncation = M_truncation, M_P_null = M_P_null))
}
# Compute for all scenarios
eigenvalue_data <- lapply(1:4, function(i) {
  compute_kle_eigenvalues(i, tmb_tps, z_results)
})

#===================================================
# Create combined dataframe
#===================================================
combined_data <- bind_rows(
  lapply(1:4, function(i) {
    data.frame(k = 1:length(eigenvalue_data[[i]]$lambda_k),
      penalty_eigenvalue = eigenvalue_data[[i]]$v_k,  # v_k from penalty matrix
      kle_eigenvalue = eigenvalue_data[[i]]$lambda_k,  # \lambda_k for KLE
      alpha = eigenvalue_data[[i]]$alpha,
      M_null = eigenvalue_data[[i]]$M_P_null,
      mode_type = ifelse(1:length(eigenvalue_data[[i]]$lambda_k) <= eigenvalue_data[[i]]$M_P_null, "Null Space", "Penalized"),
      scenario = paste0("Sce. ", i))
  })
)

#===================================================
# Plot: KLE Eigenvalue Decay ($\lambda_k$)
#===================================================
plot_data <- bind_rows(lapply(seq_along(tmb_tps), function(i) {
    fit <- tmb_tps[[i]]
    # Safely compute cumulative variance
    lam <- fit$lambda_k
    cum_var <- cumsum(lam) / sum(lam)
    data.frame(scenario  = paste0("Sce. ", i),
      k  = seq_along(lam), lambda_k = lam,
      cumulative_variance = cum_var,
      M_truncation  = fit$M_truncation,
      k_basis  = fit$k_basis)
  })
)

# Extract precise truncation point per scenario
truncation_points <- plot_data |>
  group_by(scenario) |>
  filter(k == first(M_truncation)) |>
  ungroup()

#=========================================================
# KLE Eigenvalue Decay
#=========================================================
plot_lambda <- ggplot(plot_data, aes(x = k, y = lambda_k)) +
  geom_line(colour = "steelblue", linewidth = 1.1) +
  geom_vline(data = truncation_points, aes(xintercept = M_truncation), colour = "red", linetype = 2, linewidth = 0.7) +
  geom_point(data = truncation_points, aes(x = k, y = lambda_k), colour = "red", size = 2.5) +
  # Multiplicative nudge for log scale to avoid clipping/errors
  geom_label(data = truncation_points, aes(x = k, y = lambda_k * 2.2, label = paste0("M = ", M_truncation)),
    size = 3.5, fill = "white", colour = "red", label.size = 0.25) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), expand = expansion(mult = c(0.05, 0.25))) +
  facet_wrap(~ scenario, scales = "free_x", ncol = 2) +
  labs(title = "Decay of KLE Eigenvalues", x = "Basis function index (k)", y = expression(lambda[k])) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))

#=========================================================
# Cumulative Variance Explained
#=========================================================
plot_cumvar <- ggplot(plot_data, aes(x = k, y = cumulative_variance)) +
  geom_line(colour = "steelblue", linewidth = 1.1) +
  geom_hline(yintercept = 0.99, colour = "darkgrey", linetype = 3, linewidth = 0.8) +
  geom_vline(data = truncation_points, aes(xintercept = M_truncation), colour = "red", linetype = 2, linewidth = 0.7) +
  geom_point(data = truncation_points, aes(x = k, y = cumulative_variance), colour = "red", size = 2.5) +
  geom_label(data = truncation_points, aes(x = k, y = cumulative_variance,  label = paste0("M = ", M_truncation)),
    size = 3.3, fill = "white", colour = "red", label.size = 0.25, vjust = 1.35 ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ scenario, scales = "free_x", ncol = 2) +
  labs(title = "Cumulative Variance Explained", x = "Number of Basis Functions (k)", y = "Cumulative Proportion of Variance") +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "black"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank())

#===================================================
# Plot 3: Penalty Eigenvalue Spectrum (v_k)
#===================================================
plot4 <- ggplot(combined_data, aes(x = k, y = penalty_eigenvalue)) +
  geom_point(aes(color = mode_type), size = 3, alpha = 1.2, shape = 1) +
  geom_line(aes(group = 1), linewidth = 0.8, color = "blue") +
  facet_wrap(~ scenario, scales = "free", ncol = 2) +
  scale_y_log10() +
  scale_color_manual(values = c("Null Space" = "#E41A1C", 
                                "Penalized" = "#377EB8")) +
  labs(title = "Eigenvalues of the TPS roughness penalty matrix S",
    subtitle = expression("v"[k] * " from S = ψΛψ'"), x = "Mode Index (increasing roughness)", y = expression("Penalty Eigenvalue (v"[k]*")"), color = "Mode Type") +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "black"),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank())


#=====================================================================================
# EXTRACT SPECTRAL RESULTS FUNCTION 
#=====================================================================================
extract_spectral_results <- function(scenario_obj, scenario_name = "Scenario") {
  # Extract components using $ notation
  opt <- scenario_obj$opt
  tmb_data <- scenario_obj$tmb_data
  M_truncation <- scenario_obj$M_truncation
  # Get estimated parameters
  alpha_est <- exp(opt$par["logalpha"])
  sigma_est <- exp(opt$par["logsigma"])
  # Get S_diag from tmb_data
  S_diag <- tmb_data$S_diag_truncated
  M_P_null_space <- tmb_data$M_P_null_space
  
  # Compute spectral quantities (Hankel transform perspective for d=2)
  # For biharmonic operator: μ_k ≈ ω_k^4
  omega_k <- (pmax(S_diag, 1e-10))^(1/4)
  wavelength_k <- 2*pi / omega_k
  # Regularized kernel eigenvalues with estimated alpha
  lambda_k_est <- 1 / (1 + alpha_est * pmax(S_diag, 0))
  lambda_k_est[S_diag <= 1e-10] <- 1.0
  
  # Explained variance
  total_var <- sum(lambda_k_est)
  cumulative_var <- cumsum(lambda_k_est)
  explained_var_ratio <- cumulative_var / total_var
  # Spectral density approximation (from Hankel transform theory)
  spectral_density <- lambda_k_est
  # Power spectrum (variance per mode)
  power_spectrum <- lambda_k_est
  # Create results dataframe
  spectral_df <- data.frame(mode = 1:M_truncation, mu_k = S_diag, omega_k = omega_k, wavelength = wavelength_k, lambda_k = lambda_k_est,
                            spectral_density = spectral_density, power_spectrum = power_spectrum, explained_var_cumulative = explained_var_ratio,
    explained_var_incremental = c(explained_var_ratio[1], diff(explained_var_ratio)), is_null_space = 1:M_truncation <= M_P_null_space,
    scenario = scenario_name, alpha = alpha_est, sigma = sigma_est, M_truncation = M_truncation, M_null_space = M_P_null_space)
 return(spectral_df)
}


#====================================================
#             PLOT DEL PAPER (Figure 5)
#====================================================
fits_TMB_spde <- readRDS('fits_TMB_spde.RDS')
fits_TMB_tps <- readRDS('fits_TMB_tps.RDS')


# Function to perform all calculations for a given scenario
process_scenario <- function(scenario_number) {
  # Reconstruct the true covariance and SPDE covariance
  cov_true <- fits_TMB_spde[[scenario_number]]$Cov_true
  
  rho <- fits_TMB_spde[[scenario_number]]$opt$par[2]
  sigma_u <- fits_TMB_spde[[scenario_number]]$opt$par[3]
  mesh <- fits_TMB_spde[[scenario_number]]$mesh
  spde <- fits_TMB_spde[[scenario_number]]$spde
  Q_spde <- inla.spde2.precision(spde, theta = c(log(rho), log(sigma_u)))
  cov_spde <- as.matrix(solve(Q_spde))
  
  # Reconstruct and project the regTPS-KLE covariance
  S_diag <- fits_TMB_tps[[scenario_number]]$S_diag_full
  evectors <- fits_TMB_tps[[scenario_number]]$evectors
  alpha_est <- exp(fits_TMB_tps[[scenario_number]]$rep$par.fixed["logalpha"])
  cov_regTPS <- evectors %*% diag(1 / (1 + alpha_est * S_diag)) %*% t(evectors)
  
  sm_obj <- fits_TMB_tps[[scenario_number]]$sm
  mesh_loc <- mesh$loc
  X_mesh <- mgcv::PredictMat(sm_obj, data.frame(s1 = mesh_loc[, 1], s2 = mesh_loc[, 2]))
  cov_regTPS_projected <- X_mesh %*% cov_regTPS %*% t(X_mesh)
  
  # Calculate the difference matrices
  diff_tps_true <- cov_regTPS_projected - cov_true
  diff_spde_true <- cov_spde - cov_true
  
  # Get the number of nodes for indexing
  n_nodes <- nrow(cov_true)
  
  # Function to convert a matrix to a long data frame
  prep_data <- function(mat, model_name, scenario_num) {
    as.data.frame(mat) %>%
      mutate(s1_idx = 1:n_nodes) %>%
      pivot_longer(-s1_idx, names_to = "s2_idx", values_to = "value") %>%
      mutate(model = model_name, scenario = paste0("Sce.", scenario_num), s2_idx = as.numeric(s2_idx))
  }
  
  # Prepare data frames for the two difference matrices
  plot_data_diff_tps <- prep_data(diff_tps_true, "regTPS-KLE - Cov True", scenario_number)
  plot_data_diff_spde <- prep_data(diff_spde_true, "SPDE - Cov True", scenario_number)
  
  # Combine and normalize the data for this scenario
  combined_data_scenario <- bind_rows(plot_data_diff_tps, plot_data_diff_spde)
  return(combined_data_scenario)
}

# Process both scenarios and combine the results
combined_all_data <- bind_rows(process_scenario(1), 
                               process_scenario(2),
                               process_scenario(3),
                               process_scenario(4)) %>%
  group_by(model, scenario) %>%  # Scale within each model-scenario combination
  mutate(scaled_value = value / max(abs(value))) %>%
  ungroup()

filter(combined_all_data, model == "regTPS-KLE - Cov True")
filter(combined_all_data, model == "SPDE - Cov True")

# Plot the combined heatmaps with facet_grid and free scales
combined_all_data$model <- factor(combined_all_data$model, levels = c("SPDE - Cov True", "regTPS-KLE - Cov True"))
plot5 <- ggplot(combined_all_data, aes(x = s1_idx, y = s2_idx, fill = scaled_value)) +
  geom_tile() + facet_grid(scenario ~ model, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name = "Scaled\nDifference") +
  labs(title = "Comparison of Covariance Model Differences", x = "Mesh Node Index", y = "Mesh Node Index") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black")) +
        ggh4x::facet_grid2(scenario ~ model, scales = "free", independent = "all")

plot5

#=========================================
# Field Comparison Plot
#=========================================
library(reshape2)

#===================================================
# Reshape function
#===================================================
reshape_data <- function(true_data,  spde_model, tps_model,
                         scenario, A_grid, Phi_kle_grid,
                         S_diag_truncated, M_P_null_space,
                         dim_grid = 30) {
  
  #===================================================
  # SPDE: Apply non-centered transformation
  #===================================================
  spde_post <- rstan::extract(spde_model)
  u_tilde_draws <- as.matrix(spde_post$u_tilde)
  rho_draws     <- as.vector(spde_post$rho)
  sigma_u_draws <- as.vector(spde_post$sigma_u)
  n_iter <- nrow(u_tilde_draws)
  n_grid <- nrow(A_grid)
  
  # Transform and project each iteration
  spde_field_sum <- numeric(n_grid)
  for (iter in 1:n_iter) {
    kappa_iter <- sqrt(8) / rho_draws[iter]
    tau_iter   <- 1.0 / (kappa_iter * sigma_u_draws[iter])
    u_iter <- u_tilde_draws[iter, ] / tau_iter
    spde_field_sum <- spde_field_sum + as.vector(A_grid %*% u_iter)
  }
  spde_field <- spde_field_sum / n_iter
  
  #===================================================
  # regTPS-KLE: Apply non-centered transformation
  #===================================================
  tps_post <- rstan::extract(tps_model)
  z_tilde_draws <- as.matrix(tps_post$z_tilde)
  alpha_draws   <- exp(as.vector(tps_post$logalpha))
  
  n_iter_tps <- nrow(z_tilde_draws)
  M <- ncol(z_tilde_draws)
  n_grid_tps <- nrow(Phi_kle_grid)
  
  # Transform and project each iteration
  tps_field_sum <- numeric(n_grid_tps)
  for (iter in 1:n_iter_tps) {
    alpha_iter <- alpha_draws[iter]
    
    z_iter <- numeric(M)
    
    # Null space
    if (M_P_null_space > 0) {
      z_iter[1:M_P_null_space] <- z_tilde_draws[iter, 1:M_P_null_space]
    }
    
    # Penalized modes
    if (M > M_P_null_space) {
      idx_penalized <- (M_P_null_space + 1):M
      S_k <- S_diag_truncated[idx_penalized]
      scale_factor <- 1.0 / sqrt(1.0 + alpha_iter * S_k + 1e-10)
      z_iter[idx_penalized] <- z_tilde_draws[iter, idx_penalized] * scale_factor
    }
    
    tps_field_sum <- tps_field_sum + as.vector(Phi_kle_grid %*% z_iter)
  }
  tps_field <- tps_field_sum / n_iter_tps
  

  true_df <- melt(matrix(true_data, dim_grid, dim_grid)) %>% mutate(type = "True GRF")
  spde_df <- melt(matrix(spde_field, dim_grid, dim_grid)) %>% mutate(type = "SPDE")
  tps_df <- melt(matrix(tps_field, dim_grid, dim_grid)) %>% mutate(type = "regTPS-KLE")
  bind_rows(true_df, spde_df, tps_df) %>% mutate(scenario = paste0("Sce. ", scenario))
}

#===================================================
# Apply to all scenarios
#===================================================
all_data <- bind_rows(reshape_data(true_data = tmb_tps[[1]][[12]], spde_model = mcmc_spde1, tps_model = mcmc_tps1, 
    scenario = 1,A_grid = tmb_spde[[1]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[1]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[1]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[1]]$tmb_data$M_P_null_space),
  reshape_data(true_data = tmb_tps[[2]][[12]], spde_model = mcmc_spde2,  tps_model = mcmc_tps2, 
    scenario = 2, A_grid = tmb_spde[[2]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[2]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[2]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[2]]$tmb_data$M_P_null_space),
  reshape_data(true_data = tmb_tps[[3]][[12]], spde_model = mcmc_spde3, tps_model = mcmc_tps3, 
    scenario = 3, A_grid = tmb_spde[[3]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[3]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[3]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[3]]$tmb_data$M_P_null_space),
  reshape_data(true_data = tmb_tps[[4]][[12]], spde_model = mcmc_spde4, tps_model = mcmc_tps4, 
    scenario = 4, A_grid = tmb_spde[[4]]$tmb_data$A_grid,
    Phi_kle_grid = tmb_tps[[4]]$tmb_data$Phi_kle_grid,
    S_diag_truncated = tmb_tps[[4]]$tmb_data$S_diag_truncated,
    M_P_null_space = tmb_tps[[4]]$tmb_data$M_P_null_space)
)

#===================================================
# Order factor levels
#===================================================
all_data$type <- factor(all_data$type, levels = c("True GRF", "SPDE", "regTPS-KLE"))

#===================================================
# Per-scenario normalization
#===================================================
all_data <- all_data %>%
  group_by(scenario, type) %>%
  mutate(value_norm = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

#===================================================
# Main plot
#===================================================
plot6 <- ggplot(all_data, aes(x = Var1, y = Var2, fill = value_norm)) +
  geom_tile() +
  facet_grid(scenario ~ type) +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Comparison of True and Approximated GRFs",
    x = "X-coordinate", y = "Y-coordinate",  fill = "Field Values") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 14, colour = "black"),
    strip.text.y = element_text(size = 14, colour = "black"),
    strip.background = element_rect(fill = "gray90"),
    axis.ticks = element_line(color = "black"),
    panel.spacing = unit(0.5, "lines")) +
  coord_fixed(ratio = 1)

#===============================
#            Plot 7
#===============================
# Simulating from the SIMULATE function
# #==================================
# # Compile the model and load it
dyn.load(dynlib("C:/Users/jcavi/OneDrive/Escritorio/KLE/spde"))

#==================================
# Compile the model and load it
dyn.load(dynlib("C:/Users/jcavi/OneDrive/Escritorio/KLE/tps_kle"))



# Creating a NA matrix values
set.seed(1234)
mat_sim = matrix(data=NA, nrow=length(tmb_spde[[4]][[1]]$simulate()$y_sim), ncol=100)
mat_sim
for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = tmb_spde[[4]][[1]]$simulate()$y_sim
  }
}

# mat_sim
df1_spde <- data.frame(tmb_spde[[1]][[1]]$simulate()$y_sim, mat_sim)
names(df1_spde)[names(df1_spde) == 'tmb_spde..1....1...simulate...y_sim'] <- 'y_sim'

df2_spde <- data.frame(tmb_spde[[2]][[1]]$simulate()$y_sim, mat_sim)
names(df2_spde)[names(df2_spde) == 'tmb_spde..2....1...simulate...y_sim'] <- 'y_sim'

df3_spde <- data.frame(tmb_spde[[3]][[1]]$simulate()$y_sim, mat_sim)
names(df3_spde)[names(df3_spde) == 'tmb_spde..3....1...simulate...y_sim'] <- 'y_sim'

df4_spde <- data.frame(tmb_spde[[4]][[1]]$simulate()$y_sim, mat_sim, mat_sim)
names(df4_spde)[names(df4_spde) == 'tmb_spde..4....1...simulate...y_sim'] <- 'y_sim'




# Histogram with kernel density
p1_spde <- ggplot(df1_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 1", colour = "blue", size = 10)

for (i in 2:ncol(df1_spde)) {
  p1_spde <- p1_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1_spde[, i]), 
                                       sd = sd(df1_spde[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p2_spde <- ggplot(df2_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 2", colour = "blue", size = 10)

for (i in 2:ncol(df2_spde)) {
  p2_spde <- p2_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2_spde[, i]), 
                                       sd = sd(df2_spde[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p3_spde <- ggplot(df3_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 3", colour = "blue", size = 10)

for (i in 2:ncol(df3_spde)) {
  p3_spde <- p3_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df3_spde[, i]), 
                                       sd = sd(df3_spde[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p4_spde <- ggplot(df4_spde, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 4", colour = "blue", size = 10)

for (i in 2:ncol(df4_spde)) {
  p4_spde <- p4_spde + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4_spde[, i]), 
                                       sd = sd(df4_spde[, i])), lwd = 1, col = 'orange')}


library(gridExtra)
grid.arrange(p1_spde, p2_spde, p3_spde, p4_spde, ncol = 2)



#=======================================================
# Using face_wrap for plot SPDE
#=======================================================

# 1) Hist data: keep only y_sim + scenario (different lengths are fine)
df_hist <- bind_rows(df1_spde %>% transmute(y_sim, scenario = "Sce. 1"),
  df2_spde %>% transmute(y_sim, scenario = "Sce. 2"),
  df3_spde %>% transmute(y_sim, scenario = "Sce. 3"),
  df4_spde %>% transmute(y_sim, scenario = "Sce. 4"))

# Helper: compute mu/sigma for all simulated columns (except y_sim) in one df
compute_params <- function(df, scenario_label) {
  df <- dplyr::as_tibble(df)
  # Ensure unique column names (handles duplicated mat_sim in df4)
  names(df) <- make.unique(names(df))
  sim_cols <- setdiff(names(df), "y_sim")
  if (length(sim_cols) == 0) return(tibble())  # no extra sim cols
  
  tibble(
    scenario = scenario_label,
    sim_id   = sim_cols,
    mu       = sapply(df[sim_cols], function(x) mean(x, na.rm = TRUE)),
    sigma    = sapply(df[sim_cols], function(x) sd(x,   na.rm = TRUE))
  ) %>%
    filter(is.finite(mu), is.finite(sigma), !is.na(mu), !is.na(sigma), sigma > 0)
}

# 2) Normal params per scenario (no row alignment issues)
normal_params <- bind_rows(
  compute_params(df1_spde, "Sce. 1"),
  compute_params(df2_spde, "Sce. 2"),
  compute_params(df3_spde, "Sce. 3"),
  compute_params(df4_spde, "Sce. 4")
)

# 3) Build precomputed normal curves over each scenario's x-range
make_curves <- function(params_df, hist_df) {
  if (nrow(params_df) == 0) return(tibble())
  curves_list <- lapply(unique(params_df$scenario), function(sc) {
    psc <- dplyr::filter(params_df, scenario == sc)
    rng <- range(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
    # if degenerate range, widen a bit:
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      m <- mean(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
      rng <- c(m - 1, m + 1)
    }
    x <- seq(rng[1], rng[2], length.out = 400)
    # Cartesian product to compute densities for every sim_id at every x
    merge(psc, data.frame(x = x), by = NULL) |>
      mutate(density = dnorm(x, mean = mu, sd = sigma))
  })
  bind_rows(curves_list)
}

curves_df <- make_curves(normal_params, df_hist)

# 4) Plot: faceted histogram + KDE + overlaid normal fits
ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  # geom_density(linewidth = 1) +
  geom_line(
    data = curves_df,
    aes(x = x, y = density, group = sim_id),
    linewidth = 0.7, alpha = 0.9, color = "orange"
  ) +
  facet_wrap(~ scenario, ncol = 2) +
  theme_bw(base_size = 14) +
  labs(title = "Simulated values using SPDE", x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))





## Posterior predictive check
# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(tmb_tps[[4]][[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = tmb_tps[[4]][[1]]$simulate()$y_sim
  }
}
# mat_sim


df1_tps <- data.frame(tmb_tps[[1]][[1]]$simulate()$y_sim, mat_sim)
names(df1_tps)[names(df1_tps) == 'tmb_tps..1....1...simulate...y_sim'] <- 'y_sim'

df2_tps <- data.frame(tmb_tps[[2]][[1]]$simulate()$y_sim, mat_sim)
names(df2_tps)[names(df2_tps) == 'tmb_tps..2....1...simulate...y_sim'] <- 'y_sim'

df3_tps <- data.frame(tmb_tps[[3]][[1]]$simulate()$y_sim, mat_sim)
names(df3_tps)[names(df3_tps) == 'tmb_tps..3....1...simulate...y_sim'] <- 'y_sim'

df4_tps <- data.frame(tmb_tps[[4]][[1]]$simulate()$y_sim, mat_sim, mat_sim)
names(df4_tps)[names(df4_tps) == 'tmb_tps..4....1...simulate...y_sim'] <- 'y_sim'




# Histogram with kernel density
p1_tps <- ggplot(df1_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 1", colour = "blue", size = 10)

for (i in 2:ncol(df1_tps)) {
  p1_tps <- p1_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1_tps[, i]), 
                                       sd = sd(df1_tps[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p2_tps <- ggplot(df2_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 2", colour = "blue", size = 10)

for (i in 2:ncol(df2_tps)) {
  p2_tps <- p2_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2_tps[, i]), 
                                       sd = sd(df2_tps[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p3_tps <- ggplot(df3_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 3", colour = "blue", size = 10)

for (i in 2:ncol(df3_tps)) {
  p3_tps <- p3_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df3_tps[, i]), 
                                       sd = sd(df3_tps[, i])), lwd = 1, col = 'orange')}





# Histogram with kernel density
p4_tps <- ggplot(df4_tps, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Sce. 4", colour = "blue", size = 10)

for (i in 2:ncol(df4_tps)) {
  p4_tps <- p4_tps + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4_tps[, i]), 
                                       sd = sd(df4_tps[, i])), lwd = 1, col = 'orange')}


library(gridExtra)
grid.arrange(p1_tps, p2_tps, p3_tps, p4_tps, ncol = 2)


#=======================================================
# Using face_wrap

# 1) Hist data: keep only y_sim + scenario (different lengths are fine)
df_hist <- bind_rows(
  df1_tps %>% transmute(y_sim, scenario = "Sce. 1"),
  df2_tps %>% transmute(y_sim, scenario = "Sce. 2"),
  df3_tps %>% transmute(y_sim, scenario = "Sce. 3"),
  df4_tps %>% transmute(y_sim, scenario = "Sce. 4")
)

# Helper: compute mu/sigma for all simulated columns (except y_sim) in one df
compute_params <- function(df, scenario_label) {
  df <- dplyr::as_tibble(df)
  # Ensure unique column names (handles duplicated mat_sim in df4)
  names(df) <- make.unique(names(df))
  sim_cols <- setdiff(names(df), "y_sim")
  if (length(sim_cols) == 0) return(tibble())  # no extra sim cols
  
  tibble(
    scenario = scenario_label,
    sim_id   = sim_cols,
    mu       = sapply(df[sim_cols], function(x) mean(x, na.rm = TRUE)),
    sigma    = sapply(df[sim_cols], function(x) sd(x,   na.rm = TRUE))
  ) %>%
    filter(is.finite(mu), is.finite(sigma), !is.na(mu), !is.na(sigma), sigma > 0)
}

# 2) Normal params per scenario (no row alignment issues)
normal_params <- bind_rows(
  compute_params(df1_tps, "Sce. 1"),
  compute_params(df2_tps, "Sce. 2"),
  compute_params(df3_tps, "Sce. 3"),
  compute_params(df4_tps, "Sce. 4")
)

# 3) Build precomputed normal curves over each scenario's x-range
make_curves <- function(params_df, hist_df) {
  if (nrow(params_df) == 0) return(tibble())
  curves_list <- lapply(unique(params_df$scenario), function(sc) {
    psc <- dplyr::filter(params_df, scenario == sc)
    rng <- range(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
    # if degenerate range, widen a bit:
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      m <- mean(hist_df$y_sim[hist_df$scenario == sc], na.rm = TRUE)
      rng <- c(m - 1, m + 1)
    }
    x <- seq(rng[1], rng[2], length.out = 400)
    # Cartesian product to compute densities for every sim_id at every x
    merge(psc, data.frame(x = x), by = NULL) |>
      mutate(density = dnorm(x, mean = mu, sd = sigma))
  })
  bind_rows(curves_list)
}

curves_df <- make_curves(normal_params, df_hist)

# 4) Plot: faceted histogram + KDE + overlaid normal fits
ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  # geom_density(linewidth = 1) +
  geom_line(
    data = curves_df,
    aes(x = x, y = density, group = sim_id),
    linewidth = 0.7, alpha = 0.9, color = "orange"
  ) +
  facet_wrap(~ scenario, ncol = 2) +
  theme_bw(base_size = 14) +
  labs(title = "Simulated values using regTPS-KLE", x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))






#=====================================
#            PLOT TOTAL 
#=====================================
prep_df <- function(df, scenario_label, model_label) {
  df <- dplyr::as_tibble(df)
  names(df) <- make.unique(names(df))   # fix duplicates
  df %>%
    mutate(scenario = scenario_label, model = model_label)
}
compute_params <- function(df, scenario_label, model_label) {
  df <- dplyr::as_tibble(df)
  names(df) <- make.unique(names(df))
  sim_cols <- setdiff(names(df), "y_sim")
  if (length(sim_cols) == 0) return(tibble())
  tibble(
    scenario = scenario_label,
    model    = model_label,
    sim_id   = sim_cols,
    mu       = sapply(df[sim_cols], function(x) mean(x, na.rm = TRUE)),
    sigma    = sapply(df[sim_cols], function(x) sd(x,   na.rm = TRUE))
  ) %>%
    filter(is.finite(mu), is.finite(sigma), !is.na(mu), !is.na(sigma), sigma > 0)
}
make_curves <- function(params_df, hist_df) {
  if (nrow(params_df) == 0) return(tibble())
  curves_list <- lapply(unique(params_df$scenario), function(sc) {
    for_model <- params_df %>% filter(scenario == sc)
    rng <- range(hist_df$y_sim[hist_df$scenario == sc & 
                                 hist_df$model == unique(for_model$model)], na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      m <- mean(hist_df$y_sim[hist_df$scenario == sc &
                                hist_df$model == unique(for_model$model)], na.rm = TRUE)
      rng <- c(m - 1, m + 1)
    }
    x <- seq(rng[1], rng[2], length.out = 400)
    merge(for_model, data.frame(x = x), by = NULL) |>
      mutate(density = dnorm(x, mean = mu, sd = sigma))
  })
  bind_rows(curves_list)
}
# --- Build data frames (SPDE + TPS) ---
df_hist <- bind_rows(
  prep_df(df1_spde, "Sce. 1", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df2_spde, "Sce. 2", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df3_spde, "Sce. 3", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df4_spde, "Sce. 4", "SPDE") %>% dplyr::select(y_sim, scenario, model),
  prep_df(df1_tps,  "Sce. 1", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model),
  prep_df(df2_tps,  "Sce. 2", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model),
  prep_df(df3_tps,  "Sce. 3", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model),
  prep_df(df4_tps,  "Sce. 4", "regTPS-KLE")  %>% dplyr::select(y_sim, scenario, model)
)
normal_params <- bind_rows(
  compute_params(df1_spde, "Sce. 1", "SPDE"),
  compute_params(df2_spde, "Sce. 2", "SPDE"),
  compute_params(df3_spde, "Sce. 3", "SPDE"),
  compute_params(df4_spde, "Sce. 4", "SPDE"),
  compute_params(df1_tps,  "Sce. 1", "regTPS-KLE"),
  compute_params(df2_tps,  "Sce. 2", "regTPS-KLE"),
  compute_params(df3_tps,  "Sce. 3", "regTPS-KLE"),
  compute_params(df4_tps,  "Sce. 4", "regTPS-KLE")
)
curves_df <- make_curves(normal_params, df_hist)
# --- Plot ---
# Changed factor levels order to put regTPS-KLE first (left column)
df_hist$model <- factor(df_hist$model, levels = c("SPDE", "regTPS-KLE"))
df_hist$scenario <- factor(df_hist$scenario, levels = c("Sce. 1", "Sce. 2", "Sce. 3", "Sce. 4"))

# Also update curves_df factor levels if needed
curves_df$model <- factor(curves_df$model, levels = c("SPDE", "regTPS-KLE"))

plot7 <- ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  # geom_density(linewidth = 1, color = "blue") +
  geom_line(data = curves_df,
            aes(x = x, y = density, group = sim_id),
            color = "orange", linewidth = 0.6, alpha = 0.7) +
  facet_grid(scenario ~ model) +   
  theme_bw(base_size = 14) +
  labs(title = "Simulated Values From Posteriors", x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        # plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))



# Compute mean per scenario/model
mean_df <- df_hist %>%
  group_by(scenario, model) %>%
  summarise(mean_val = mean(y_sim, na.rm = TRUE), .groups = "drop")

# Add vertical lines for means
plot7 <- ggplot(df_hist, aes(x = y_sim)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.35, fill = "grey", color = "black", alpha = 0.6) +
  geom_line(data = curves_df,
            aes(x = x, y = density, group = sim_id),
            color = "orange", linewidth = 0.6, alpha = 0.7) +
  geom_vline(data = mean_df, aes(xintercept = mean_val),
             color = "red", linetype = "dashed", linewidth = 1) +
  facet_grid(scenario ~ model) +   
  theme_bw(base_size = 14) +
  labs(title = "Simulated Values From Posteriors", 
       x = "Simulated Values", 
       y = "Density") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        # plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(color = "black"))

plot7

#====================================================
#                      Plot
#====================================================

# Tidy simulations function 
get_sim_df <- function(tmb_obj, model_name, n_scenarios = 4, n_reps = 100) {
  df_list <- list()
  
  for (k in 1:n_scenarios) {
    mat_sim <- replicate(n_reps, tmb_obj[[k]][[1]]$simulate()$y_sim)
    df_tmp <- as.data.frame(mat_sim)
    df_tmp$y_obs <- tmb_obj[[k]][[1]]$simulate()$y_sim
    
    df_tmp <- df_tmp %>%
      pivot_longer(-y_obs, names_to = "rep", values_to = "sim") %>%
      mutate(scenario = paste0("Sce. ", k),
             model = model_name)
    
    df_list[[k]] <- df_tmp
  }
  
  bind_rows(df_list)
}

# --- Build SPDE + TPS ---
df_spde <- get_sim_df(tmb_spde, "SPDE")
df_tps  <- get_sim_df(tmb_tps,  "regTPS-KLE")

df_all <- bind_rows(df_spde, df_tps)

df_all$model     <- factor(df_all$model, levels = c("SPDE", "regTPS-KLE"))
df_all$scenario  <- factor(df_all$scenario, levels = paste0("Sce. ", 1:4))

# --- Compute densities per replicate ---
dens_list <- df_all %>%
  group_split(model, scenario, rep) %>%
  map(~{
    d <- density(.x$sim, from = min(df_all$sim), to = max(df_all$sim))
    data.frame(
      x = d$x,
      y = d$y,
      model = .x$model[1],
      scenario = .x$scenario[1],
      rep = .x$rep[1]
    )
  })

dens_df <- bind_rows(dens_list)

# --- Average density across replicates ---
dens_mean <- dens_df %>%
  group_by(model, scenario, x) %>%
  summarise(y = mean(y), .groups = "drop")

# --- Plot ---
p_all <- ggplot(df_all, aes(x = sim, group = rep)) +
  # observed histogram
  geom_histogram(aes(x = y_obs, y = after_stat(density)), 
                 fill = "grey80", color = "black", binwidth = 0.35, inherit.aes = FALSE) +
  # replicate densities (light lines)
  geom_line(data = dens_df, aes(x = x, y = y, group = rep), 
            col = "lightsteelblue2", alpha = 0.3, size = 0.4, inherit.aes = FALSE) +
  # average density (bold line)
  geom_line(data = dens_mean, aes(x = x, y = y), 
            col = "steelblue4", size = 1.2, inherit.aes = FALSE) +
  facet_grid(scenario ~ model) +
  labs(title = "Simulated Values From The Posteriors",
       x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 13))

p_all











# Function to tidy simulations (no y_obs here) 
get_sim_df <- function(tmb_obj, model_name, n_scenarios = 4, n_reps = 100) {
  df_list <- list()
  
  for (k in 1:n_scenarios) {
    mat_sim <- replicate(n_reps, tmb_obj[[k]][[1]]$simulate()$y_sim)
    df_tmp <- as.data.frame(mat_sim)
    
    df_tmp <- df_tmp %>%
      pivot_longer(everything(), names_to = "rep", values_to = "sim") %>%
      mutate(scenario = paste0("Sce. ", k),
             model = model_name)
    
    df_list[[k]] <- df_tmp
  }
  
  bind_rows(df_list)
}

# --- observed data per scenario (take from SPDE for consistency) ---
y_obs_list <- map(1:4, ~{
  tibble(
    y_obs = tmb_spde[[.x]][[1]]$simulate()$y_sim,
    scenario = paste0("Sce. ", .x)
  )
})
y_obs_df <- bind_rows(y_obs_list)

# --- expand y_obs across both models ---
y_obs_df <- y_obs_df %>%
  crossing(model = c("SPDE", "regTPS-KLE"))

# --- Build SPDE + TPS simulation data ---
df_spde <- get_sim_df(tmb_spde, "SPDE")
df_tps  <- get_sim_df(tmb_tps,  "regTPS-KLE")

df_all <- bind_rows(df_spde, df_tps)

# --- set factor levels so SPDE = first column, regTPS-KLE = second ---
model_levels <- c("SPDE", "regTPS-KLE")
scenario_levels <- paste0("Sce. ", 1:4)

df_all$model    <- factor(df_all$model, levels = model_levels)
df_all$scenario <- factor(df_all$scenario, levels = scenario_levels)

y_obs_df$model    <- factor(y_obs_df$model, levels = model_levels)
y_obs_df$scenario <- factor(y_obs_df$scenario, levels = scenario_levels)

# --- Compute densities per replicate ---
dens_list <- df_all %>%
  group_split(model, scenario, rep) %>%
  map(~{
    d <- density(.x$sim, from = min(df_all$sim), to = max(df_all$sim))
    data.frame(
      x = d$x,
      y = d$y,
      model = .x$model[1],
      scenario = .x$scenario[1],
      rep = .x$rep[1]
    )
  })

dens_df <- bind_rows(dens_list)

# Average density across replicates 
dens_mean <- dens_df %>%
  group_by(model, scenario, x) %>%
  summarise(y = mean(y), .groups = "drop")

# Plot 
plot7 <- ggplot(df_all, aes(x = sim, group = rep)) +
  # observed histogram (identical across models within scenario)
  geom_histogram(data = y_obs_df, 
                 aes(x = y_obs, y = after_stat(density)),
                 fill = "grey80", color = "black", binwidth = 0.35,
                 inherit.aes = FALSE) +
  # replicate densities (light lines)
  geom_line(data = dens_df, aes(x = x, y = y, group = rep), 
            col = "lightsteelblue2", alpha = 0.3, size = 0.4, inherit.aes = FALSE) +
  # average density (bold line)
  geom_line(data = dens_mean, aes(x = x, y = y), 
            col = "steelblue4", size = 1.2, inherit.aes = FALSE) +
  facet_grid(scenario ~ model) +
  labs(title = "Simulated Values From The Posteriors",
       x = "Simulated Values", y = "Density") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 13))

plot7







# Save as high-quality PDF
ggsave(filename = "C:/Users/jcavi/OneDrive/Escritorio/KLE/plots/plot7.pdf",
       plot = plot7,        # Replace with your ggplot object name
       device = cairo_pdf,    # Good for embedding text as text
       width = 10,             # Width in inches
       height = 10,            # Height in inches
       dpi = 300              # Only affects raster elements, safe to keep high
)





#=================================
#        RMSE, R2 and MAE
#=================================
extract_spde_field <- function(stan_fit, A_grid) {
  post <- rstan::extract(stan_fit)
  u_draws       <- post$u_tilde
  rho_draws     <- post$rho
  sigma_u_draws <- post$sigma_u
  
  n_iter <- nrow(u_draws)
  n_grid <- nrow(A_grid)
  
  field_draws <- matrix(NA, n_iter, n_grid)
  
  for (iter in seq_len(n_iter)) {
    kappa_iter <- sqrt(8) / rho_draws[iter]
    tau_iter   <- 1.0 / (kappa_iter * sigma_u_draws[iter])
    u_iter     <- u_draws[iter, ] / tau_iter
    field_draws[iter, ] <- as.vector(A_grid %*% u_iter)
  }
  
  field_draws   # n_iter x n_grid matrix, one row per posterior draw
}

#===================================================
# Extracting per-draw field reconstructions (regTPS-KLE)
#===================================================
extract_tps_field <- function(stan_fit, Phi_kle_grid, S_diag_truncated, M_P_null_space) {
  z_tilde_draws <- rstan::extract(stan_fit)$z_tilde
  alpha_draws   <- exp(rstan::extract(stan_fit)$logalpha)
  
  n_iter <- nrow(z_tilde_draws)
  M      <- ncol(z_tilde_draws)
  n_grid <- nrow(Phi_kle_grid)
  
  field_draws <- matrix(NA, n_iter, n_grid)
  
  for (iter in seq_len(n_iter)) {
    alpha_iter <- alpha_draws[iter]
    
    scale_iter <- numeric(M)
    scale_iter[1:M_P_null_space] <- 1.0
    scale_iter[(M_P_null_space + 1):M] <-
      1.0 / sqrt(1.0 + alpha_iter * S_diag_truncated[(M_P_null_space + 1):M])
    
    Z_iter <- scale_iter * z_tilde_draws[iter, ]
    field_draws[iter, ] <- as.vector(Phi_kle_grid %*% Z_iter)
  }
  
  field_draws   # n_iter x n_grid matrix, one row per posterior draw
}

#===================================================
# Compute RMSE for each posterior draw, then average
#===================================================
compute_metrics <- function(field_draws, u_true_grid) {
  rmse_post <- apply(field_draws, 1, function(r) sqrt(mean((r - u_true_grid)^2)))
  r2_post   <- apply(field_draws, 1, function(r) cor(r, u_true_grid)^2)
  mae_post  <- apply(field_draws, 1, function(r) mean(abs(r - u_true_grid)))
  
  list(
    rmse_mean = mean(rmse_post), rmse_q025 = as.numeric(quantile(rmse_post, .025)),
    rmse_q975 = as.numeric(quantile(rmse_post, .975)),
    r2_mean   = mean(r2_post),
    mae_mean  = mean(mae_post),
    rmse_post = rmse_post, r2_post = r2_post, mae_post = mae_post   # keep raw draws too
  )
}

#===================================================
# Apply to all 4 scenarios
#===================================================
spde_field_draws <- vector("list", 4)
tps_field_draws  <- vector("list", 4)
metrics_per_draw <- vector("list", 4)

for (i in 1:4) {
  spde_field_draws[[i]] <- extract_spde_field(stan_fit = get(paste0("mcmc_spde", i)),
                                              A_grid   = tmb_spde[[i]][[4]]$A_grid)
  
  tps_field_draws[[i]] <- extract_tps_field(stan_fit = get(paste0("mcmc_tps", i)),
                                            Phi_kle_grid     = tmb_tps[[i]][[4]]$Phi_kle_grid,
                                            S_diag_truncated = tmb_tps[[i]][[4]]$S_diag_truncated,
                                            M_P_null_space   = tmb_tps[[i]][[4]]$M_P_null_space)
  
  spde_metrics <- compute_metrics(spde_field_draws[[i]], tmb_spde[[i]][[9]])
  tps_metrics  <- compute_metrics(tps_field_draws[[i]],  tmb_tps[[i]][[12]])
  
  metrics_per_draw[[i]] <- data.frame(Metric = c("RMSE", "R2", "MAE"),
                                      SPDE   = c(spde_metrics$rmse_mean, spde_metrics$r2_mean, spde_metrics$mae_mean),
                                      TPS    = c(tps_metrics$rmse_mean,  tps_metrics$r2_mean,  tps_metrics$mae_mean))
  
  cat("\n=== Scenario", i, "(mean-of-per-draw metrics) ===\n")
  print(metrics_per_draw[[i]], row.names = FALSE)
}

# Individually accessible, matching the naming style of your original script:
metrics_scenario1_perdraw <- metrics_per_draw[[1]]
metrics_scenario2_perdraw <- metrics_per_draw[[2]]
metrics_scenario3_perdraw <- metrics_per_draw[[3]]
metrics_scenario4_perdraw <- metrics_per_draw[[4]]





#=================================
# Create LaTeX table from metrics_per_draw
#=================================
library(dplyr)
library(knitr)
library(kableExtra)

#=================================
# Build data frame for the table
#=================================

table_df <- bind_rows(lapply(seq_along(metrics_per_draw), function(i){
  df <- metrics_per_draw[[i]]
  data.frame(Scenario = c(paste0("Sce. ", i), ""),
    Method   = c("SPDE", "regTPS-KLE"),
    RMSE     = c(df$SPDE[df$Metric=="RMSE"], df$TPS[df$Metric=="RMSE"]),
    R2       = c(df$SPDE[df$Metric=="R2"], df$TPS[df$Metric=="R2"]),
    MAE      = c(df$SPDE[df$Metric=="MAE"], df$TPS[df$Metric=="MAE"]))}))

table_df





latex_table <- kbl(table_df, format = "latex", booktabs = TRUE, digits = 3,
    align = c("l","l","r","r","r"),
    escape = FALSE,
    caption = "Comparison of SPDE and regTPS-KLE models with the true GRF using a Matérn covariance function across scenarios.",
    label = "table1") |>
  collapse_rows(columns = 1, valign = "top") |>
  kable_styling(latex_options = c("hold_position"))

cat(latex_table)





