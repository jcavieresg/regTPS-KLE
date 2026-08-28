rm(list = ls())

options(scipen = 999)


library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, rstan, Matrix, fields, sf,
               rnaturalearth, gridExtra, tidyr, patchwork)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

# Reading the outputs
# TMB models
spde_tmb <- readRDS('spde_tmb.RDS')
spde_mcmc <- readRDS('spde_mcmc.RDS')

regTPS_KLE_tmb <- readRDS('regTPS_KLE_tmb.RDS')
regTPS_KLE_mcmc <- readRDS('regTPS_KLE_mcmc.RDS')


#========================================================
# SPDE
#========================================================
post_spde <- rstan::extract(spde_mcmc)
cat("Available SPDE posterior parameters:", paste(names(post_spde), collapse = ", "), "\n")
u_tilde_post <- post_spde$u_tilde   # iterations x n_nodes

# Detect which parameterization was actually fit
if ("rho" %in% names(post_spde) && "sigma_u" %in% names(post_spde)) {
  cat("Detected NATURAL-SCALE parameterization (rho, sigma_u) -- no exp() needed.\n")
  rho_post     <- post_spde$rho
  sigma_u_post <- post_spde$sigma_u
} else if ("logrange" %in% names(post_spde) && "logsigma_u" %in% names(post_spde)) {
  cat("Detected LOG-SCALE (PC-prior) parameterization (logrange, logsigma_u).\n")
  rho_post     <- exp(post_spde$logrange)
  sigma_u_post <- exp(post_spde$logsigma_u)
} else if ("logrho" %in% names(post_spde) && "logsigma_u" %in% names(post_spde)) {
  cat("Detected logrho/logsigma_u parameterization.\n")
  rho_post     <- exp(post_spde$logrho)
  sigma_u_post <- exp(post_spde$logsigma_u)
} else {
  stop("Could not find a recognized (rho, sigma_u) parameterization. ",
       "Available names: ", paste(names(post_spde), collapse = ", "))
}

niter   <- nrow(u_tilde_post)
n_nodes <- ncol(u_tilde_post)
A_grid  <- spde_tmb$A_grid
A_grid_dense <- as.matrix(A_grid)
n_grid <- nrow(A_grid_dense)

y_grid_samples <- matrix(NA_real_, nrow = niter, ncol = n_grid)

for (it in seq_len(niter)) {
  u_tilde_it <- u_tilde_post[it, ]
  kappa_it   <- sqrt(8) / rho_post[it]
  tau_it     <- 1 / (kappa_it * sigma_u_post[it])
  u_it       <- u_tilde_it / tau_it
  y_grid_samples[it, ] <- as.numeric(A_grid_dense %*% u_it)
}

y_grid_samples_spde <- y_grid_samples

y_summary <- t(apply(y_grid_samples_spde, 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(y_summary) <- c("q025", "median", "q975")

df_spde <- data.frame(regTPS_KLE_tmb$grid_total[, 1], regTPS_KLE_tmb$grid_total[, 2],
                      y_summary[, "q025"], y_summary[, "median"], y_summary[, "q975"])
colnames(df_spde) <- c("s1", "s2", "q025", "median", "q975")

#========================================================
# regTPS-KLE
#========================================================
post <- rstan::extract(regTPS_KLE_mcmc)
names(post)

has_z_tilde <- "z_tilde" %in% names(post)
has_z <- "z" %in% names(post)
has_logalpha <- "logalpha" %in% names(post)
has_alpha <- "alpha" %in% names(post)

if (!(has_z_tilde || has_z)) stop("Posterior contains neither 'z_tilde' nor 'z'.")
if (!(has_logalpha || has_alpha)) stop("Posterior does not contain 'logalpha' or 'alpha'.")

if (has_z_tilde) {
  z_tilde_post <- post$z_tilde
} else {
  z_post_direct <- post$z
}

if (has_logalpha) {
  alpha_post <- exp(post$logalpha)
} else {
  alpha_post <- post$alpha
}

if (exists("z_tilde_post")) {
  niter <- nrow(z_tilde_post); M_trunc <- ncol(z_tilde_post)
} else {
  niter <- nrow(z_post_direct); M_trunc <- ncol(z_post_direct)
}

Phi_kle_grid <- regTPS_KLE_tmb$Phi_kle_grid
Phi_mat <- as.matrix(Phi_kle_grid)

if (ncol(Phi_mat) == M_trunc) {
  Phi_rows_are_grid <- TRUE
  n_grid <- nrow(Phi_mat)
} else if (nrow(Phi_mat) == M_trunc) {
  Phi_rows_are_grid <- FALSE
  n_grid <- ncol(Phi_mat)
} else {
  stop("Phi_kle_grid has incompatible dimensions vs M_trunc.")
}

y_grid_samples <- matrix(NA_real_, nrow = niter, ncol = n_grid)
S_vec <- regTPS_KLE_tmb$tmb_data$S_diag_truncated
M_P_null_space <- regTPS_KLE_tmb$tmb_data$M_P_null_space

for (it in seq_len(niter)) {
  if (exists("z_tilde_post")) {
    z_tilde_it <- z_tilde_post[it, ]
    alpha_it <- alpha_post[it]
    prior_sd <- rep(1, M_trunc)
    if (M_P_null_space < M_trunc) {
      k_idx <- (M_P_null_space + 1):M_trunc
      prior_sd[k_idx] <- sqrt(1 / (1 + alpha_it * S_vec[k_idx]))
    }
    z_it <- prior_sd * z_tilde_it
  } else {
    z_it <- z_post_direct[it, ]
  }
  
  if (Phi_rows_are_grid) {
    y_grid_samples[it, ] <- as.numeric(Phi_mat %*% z_it)
  } else {
    y_grid_samples[it, ] <- as.numeric(t(Phi_mat) %*% z_it)
  }
}

y_grid_samples_tps <- y_grid_samples

y_summary <- t(apply(y_grid_samples_tps, 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(y_summary) <- c("q025", "median", "q975")

df_regTPS_KLE <- data.frame(regTPS_KLE_tmb$grid_total[, 1], regTPS_KLE_tmb$grid_total[, 2],
                            y_summary[, "q025"], y_summary[, "median"], y_summary[, "q975"])
colnames(df_regTPS_KLE) <- c("s1", "s2", "q025", "median", "q975")

#========================================================
# spNNGP
#========================================================
grid_matrix <- as.matrix(regTPS_KLE_tmb$grid_total)

pred_spnngp <- predict(spnngp_fits[[1]],
                       X.0 = matrix(1, nrow = nrow(grid_matrix), ncol = 1),
                       coords.0 = grid_matrix,
                       n.omp.threads = no_cores)

# p.y.0 is on the sqrt scale (matches sqrt_transform = TRUE used everywhere else)
y_summary_spnngp <- t(apply(pred_spnngp$p.y.0, 1, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(y_summary_spnngp) <- c("q025", "median", "q975")

df_spnngp <- data.frame(regTPS_KLE_tmb$grid_total[, 1], regTPS_KLE_tmb$grid_total[, 2],
                        y_summary_spnngp[, "q025"], y_summary_spnngp[, "median"], y_summary_spnngp[, "q975"])
colnames(df_spnngp) <- c("s1", "s2", "q025", "median", "q975")




#=====================================================================================
# Held-out comparison: SPDE vs regTPS-KLE vs spNNGP
# 80% train / 20% test split of German stations. All three models fit on the
# same training set; predictive performance evaluated at the SAME held-out
# locations using an identical scoring function.
#=====================================================================================

library(scoringRules)  

set.seed(1234)
n_obs   <- nrow(sp_data)
test_idx  <- sample(seq_len(n_obs), size = round(0.2 * n_obs))
train_idx <- setdiff(seq_len(n_obs), test_idx)

train_data <- sp_data[train_idx, ]
test_data  <- sp_data[test_idx, ]

dim(train_data)
dim(test_data)


y_test_sqrt <- sqrt(test_data$y_obs)     # ground truth, sqrt scale (matches all models' fitting scale)
n_test <- nrow(test_data)

cat("Train n =", nrow(train_data), " | Test n =", n_test, "\n")

#=====================================================================================
# Generic scoring function
#=====================================================================================
logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

score_model <- function(mu_draws, tau2_draws, y_true) {
  n_draws <- nrow(mu_draws)
  n_test  <- ncol(mu_draws)
  stopifnot(length(tau2_draws) == n_draws)
  stopifnot(length(y_true) == n_test)
  
  #---------------------------------------------------
  # Pointwise Log Predictive Density
  #---------------------------------------------------
  lpd_per_point <- numeric(n_test)
  for (i in seq_len(n_test)) {
    logdens_s <- dnorm(y_true[i], mean = mu_draws[, i], sd   = sqrt(tau2_draws), log  = TRUE)
    lpd_per_point[i] <- logSumExp(logdens_s) - log(n_draws)
  }
  
  #---------------------------------------------------
  # Posterior predictive mean
  #---------------------------------------------------
  mu_mean <- colMeans(mu_draws)

  #---------------------------------------------------
  # Posterior predictive samples (for CRPS)
  #---------------------------------------------------
  noise_draws <- sweep(matrix(rnorm(n_draws * n_test), nrow = n_draws, ncol = n_test), 1, sqrt(tau2_draws), "*")
  y_rep <- mu_draws + noise_draws
  crps_per_point <- sapply(seq_len(n_test), function(i) {
    crps_sample(y_true[i], y_rep[, i])
  })
  
  #---------------------------------------------------
  # Prediction metrics
  #---------------------------------------------------
  rmse_val <- sqrt(mean((mu_mean - y_true)^2))
  mae_val <- mean(abs(mu_mean - y_true))
  # R² = squared Pearson correlation
  r2_val <- cor(mu_mean, y_true, use = "complete.obs")^2
  
  #---------------------------------------------------
  # Return metrics
  #---------------------------------------------------
  list(lpd_total = sum(lpd_per_point),
       rmse = rmse_val, mae = mae_val, r2 = r2_val,
       crps_mean = mean(crps_per_point),
       mu_mean = mu_mean, lpd_per_point = lpd_per_point,
       crps_per_point = crps_per_point)
}

#=====================================================================================
# Fit-and-score
#=====================================================================================
fit_and_score_all_three <- function(train_data, test_data, no_cores,
                                    mcmc_iter = 3000, mcmc_warmup = 700, spnngp_samples = 3000) {
  
  y_test_sqrt <- sqrt(test_data$y_obs)
  n_test <- nrow(test_data)
  
  #========================
  # SPDE
  #========================
  t0_spde <- Sys.time()
  spde_train <- run_tmb_spde(train_data, 100)
  spde_mcmc <- tmbstan(spde_train$obj, chains = 3, open_progress = FALSE,
                       control = list(max_treedepth = 12, adapt_delta = 0.9),
                       iter = mcmc_iter, warmup = mcmc_warmup, cores = no_cores,
                       init = 'last.par.best', seed = 12345)
  t1_spde <- Sys.time()
  time_spde_min <- as.numeric(difftime(t1_spde, t0_spde, units = "min"))
  
  post_spde <- rstan::extract(spde_mcmc)
  spde_names <- names(post_spde)
  
  if (all(c("logrho", "logsigma_u", "logsigma_e") %in% spde_names)) {
    rho_post <- exp(post_spde$logrho); sigma_u_post <- exp(post_spde$logsigma_u); sigma_e_post <- exp(post_spde$logsigma_e)
  } else if (all(c("rho", "sigma_u", "sigma") %in% spde_names)) {
    rho_post <- post_spde$rho; sigma_u_post <- post_spde$sigma_u; sigma_e_post <- post_spde$sigma
  } else if (all(c("logrange", "logsigma_u", "logsigma_e") %in% spde_names)) {
    rho_post <- exp(post_spde$logrange); sigma_u_post <- exp(post_spde$logsigma_u); sigma_e_post <- exp(post_spde$logsigma_e)
  } else {
    stop("Unrecognized SPDE parameter names: ", paste(spde_names, collapse = ", "))
  }
  
  u_tilde_post <- post_spde$u_tilde
  A_test_spde <- as.matrix(inla.spde.make.A(mesh = spde_train$mesh, loc = as.matrix(test_data[, c("s1","s2")])))
  
  mu_spde <- matrix(NA, nrow(u_tilde_post), n_test)
  for (d in seq_len(nrow(u_tilde_post))) {
    kappa_d <- sqrt(8) / rho_post[d]
    tau_field_d <- 1 / (kappa_d * sigma_u_post[d])
    mu_spde[d, ] <- as.numeric(A_test_spde %*% (u_tilde_post[d, ] / tau_field_d))
  }
  score_spde <- score_model(mu_spde, sigma_e_post^2, y_test_sqrt)
  
  #========================
  # regTPS-KLE
  #========================
  t0_tps <- Sys.time()
  tps_train <- regTPS_KLE(sp_data = train_data, dim_grid = 100, n_basis_app = 0.95,
                          variance_threshold = 0.95, sqrt_transform = TRUE, expand_grid = 0.05,
                          alpha_ref = 1.0, sigma_e_ref = 1.0)
  if (is.null(tps_train$sm) || is.null(tps_train$evectors)) {
    stop("regTPS_KLE() must return `sm` and `evectors` for out-of-sample prediction.")
  }
  
  tps_mcmc <- tmbstan(tps_train$obj, chains = 3, open_progress = FALSE,
                      control = list(max_treedepth = 12, adapt_delta = 0.9),
                      iter = mcmc_iter, warmup = mcmc_warmup, cores = no_cores,
                      init = 'last.par.best', seed = 12345)
  t1_tps <- Sys.time()
  time_tps_min <- as.numeric(difftime(t1_tps, t0_tps, units = "min"))
  
  post_tps <- rstan::extract(tps_mcmc)
  z_tilde_post <- post_tps$z_tilde
  alpha_post   <- exp(post_tps$logalpha)
  sigma_eps_post_tps <- exp(post_tps$logsigma)
  M_trunc <- ncol(z_tilde_post); M_P_null <- tps_train$M_P_null_space; S_vec <- tps_train$S_diag_truncated
  
  # Basis at TEST locations, using the TRAINING smooth object's knots/basis 
  Phi_test_kle <- mgcv::PredictMat(tps_train$sm, test_data) %*% tps_train$evectors[, 1:M_trunc]
  
  mu_tps <- matrix(NA, nrow(z_tilde_post), n_test)
  for (d in seq_len(nrow(z_tilde_post))) {
    scale_d <- rep(1, M_trunc)
    if (M_P_null < M_trunc) {
      k_idx <- (M_P_null+1):M_trunc
      scale_d[k_idx] <- 1 / sqrt(1 + alpha_post[d] * S_vec[k_idx])
    }
    mu_tps[d, ] <- as.numeric(Phi_test_kle %*% (scale_d * z_tilde_post[d, ]))
  }
  score_tps <- score_model(mu_tps, sigma_eps_post_tps^2, y_test_sqrt)
  
  #========================
  # spNNGP
  #========================
  t0_spnngp <- Sys.time()
  coords_train <- as.matrix(train_data[, c("s1","s2")])
  y_train_sqrt <- sqrt(train_data$y_obs)
  sigma_sq_start <- var(y_train_sqrt)*0.8; tau_sq_start <- var(y_train_sqrt)*0.2
  phi_start <- 3/(0.3*diff(range(coords_train[,1])))
  starting <- list(phi=phi_start, sigma.sq=sigma_sq_start, tau.sq=tau_sq_start, nu=0.5)
  tuning   <- list(phi=0.5, sigma.sq=0.5, tau.sq=0.5, nu=0.05)
  priors   <- list(phi.Unif=c(3/(2*diff(range(coords_train[,1]))), 3/(0.01*diff(range(coords_train[,1])))),
                   sigma.sq.IG=c(2,sigma_sq_start), tau.sq.IG=c(2,tau_sq_start), nu.Unif=c(0.1,2.0))
  
  mu_list <- list(); tau2_list <- list()
  for (chain in 1:3) {
    set.seed(1000+chain)
    fit_c <- spNNGP(y_train_sqrt ~ 1, coords=coords_train, method="latent", family="gaussian",
                    n.neighbors=15, starting=starting, tuning=tuning, priors=priors,
                    cov.model="matern", n.samples=spnngp_samples, n.omp.threads=no_cores, n.report=1000)
    pred_c <- predict(fit_c, X.0=matrix(1,nrow=n_test,ncol=1),
                      coords.0=as.matrix(test_data[,c("s1","s2")]), n.omp.threads=no_cores)
    beta_c <- fit_c$p.beta.samples[,1]
    mu_list[[chain]] <- t(pred_c$p.w.0) + matrix(beta_c, nrow=ncol(pred_c$p.w.0), ncol=n_test)
    theta_c <- fit_c$p.theta.samples
    tau2_list[[chain]] <- theta_c[, grep("tau\\.sq", colnames(theta_c))]
  }
  mu_spnngp <- do.call(rbind, mu_list); tau2_spnngp <- do.call(c, tau2_list)
  t1_spnngp <- Sys.time()
  time_spnngp_min <- as.numeric(difftime(t1_spnngp, t0_spnngp, units = "min"))
  
  score_spnngp <- score_model(mu_spnngp, tau2_spnngp, y_test_sqrt)
  
  data.frame(Model = c("SPDE", "regTPS-KLE", "spNNGP"),
             LogPredDensity = c(score_spde$lpd_total, score_tps$lpd_total, score_spnngp$lpd_total),
             RMSE = c(score_spde$rmse, score_tps$rmse, score_spnngp$rmse),
             MAE = c(score_spde$mae, score_tps$mae, score_spnngp$mae),
             R2 = c(score_spde$r2, score_tps$r2, score_spnngp$r2),
             CRPS = c(score_spde$crps_mean, score_tps$crps_mean, score_spnngp$crps_mean),
             Time_min = c(time_spde_min, time_tps_min, time_spnngp_min))
}

#=====================================================================================
# Repeated CV driver: "repeated" (Monte Carlo, n_reps splits) or "kfold".
#=====================================================================================

run_repeated_cv <- function(sp_data, no_cores, mode = c("repeated", "kfold"),
                            n_reps = 50, K = 10, test_frac = 0.2, seed = 1234,
                            mcmc_iter = 2000, mcmc_warmup = 500, spnngp_samples = 2000) {
  mode <- match.arg(mode)
  n_obs <- nrow(sp_data)
  all_results <- list()
  
  if (mode == "repeated") {
    for (r in seq_len(n_reps)) {
      cat("\n=========================================\n")
      cat("=== Repetition", r, "of", n_reps, "===\n")
      cat("=========================================\n")
      set.seed(seed + r)
      test_idx  <- sample(seq_len(n_obs), size = round(test_frac * n_obs))
      train_idx <- setdiff(seq_len(n_obs), test_idx)
      
      res <- fit_and_score_all_three(sp_data[train_idx, ], sp_data[test_idx, ], no_cores,
                                     mcmc_iter, mcmc_warmup, spnngp_samples)
      res$rep <- r
      all_results[[r]] <- res
    }
  } else {
    set.seed(seed)
    fold_id <- sample(rep(1:K, length.out = n_obs))
    for (k in seq_len(K)) {
      cat("\n=========================================\n")
      cat("=== Fold", k, "of", K, "===\n")
      cat("=========================================\n")
      test_idx  <- which(fold_id == k)
      train_idx <- which(fold_id != k)
      
      res <- fit_and_score_all_three(sp_data[train_idx, ], sp_data[test_idx, ], no_cores,
                                     mcmc_iter, mcmc_warmup, spnngp_samples)
      res$rep <- k
      all_results[[k]] <- res
    }
  }
  bind_rows(all_results)
}

#=====================================================================================
# Run it 
#=====================================================================================
# Option A: 50 repeated random 80/20 splits (trimmed MCMC settings for feasibility)
# cv_results <- run_repeated_cv(sp_data, no_cores, mode = "repeated", n_reps = 50,
#                               mcmc_iter = 1500, mcmc_warmup = 400, spnngp_samples = 1500)

# Option B: 10-fold CV
cv_results <- run_repeated_cv(sp_data, no_cores, mode = "kfold", K = 10,
                              mcmc_iter = 2000, mcmc_warmup = 500, spnngp_samples = 2000)

print(cv_results)

#=====================================================================================
# Aggregate -- now includes R2
#=====================================================================================
cv_summary <- cv_results %>%
  group_by(Model) %>%
  summarise(LogPredDensity_total = sum(LogPredDensity),
            RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
            MAE_mean  = mean(MAE), MAE_sd = sd(MAE),
            R2_mean   = mean(R2), R2_sd = sd(R2),
            CRPS_mean = mean(CRPS), CRPS_sd = sd(CRPS),
            Time_mean_min = mean(Time_min), Time_sd_min   = sd(Time_min),
            .groups = "drop")

