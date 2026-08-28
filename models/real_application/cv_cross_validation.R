setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs")
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
# 1. SPDE: defensive parameter extraction (checks actual names first,
#    mirroring the robust pattern already used for regTPS-KLE below)
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
# 2. regTPS-KLE (unchanged from your version -- already correctly
#    defensive about parameter naming)
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
# 3. spNNGP: predict at the SAME grid, using the already-fitted
#    full-data spnngp_fit (method="latent") from earlier
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

#========================================================
# 4. Shared color scale across ALL THREE methods, then plot together
#========================================================
library(ggplot2)
library(gridExtra)
library(sf)

coords <- data.frame(sp_data[, c(1, 2)])

min_value <- min(min(df_spde$median^2), min(df_regTPS_KLE$median^2), min(df_spnngp$median^2))
max_value <- max(max(df_spde$median^2), max(df_regTPS_KLE$median^2), max(df_spnngp$median^2))


# min_value <- min(min(df_spde$median^2), min(df_regTPS_KLE$median^2))
# max_value <- max(max(df_spde$median^2), max(df_regTPS_KLE$median^2))


make_plot <- function(df, title) {
  ggplot(df) +
    geom_raster(aes(x = s1, y = s2, fill = median^2)) +
    geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.8) +
    scale_fill_viridis_c(option = "C", limits = c(min_value, max_value)) +
    geom_point(data = coords, aes(x = s1, y = s2), color = "black", size = 1.5, alpha = 0.7) +
    coord_sf() +
    labs(title = title, fill = expression(NO[2])) +
    theme_minimal()
}

plot1 <- make_plot(df_spde,        "Posterior Median GRF -- SPDE")
plot2 <- make_plot(df_regTPS_KLE,  "Posterior Median GRF -- regTPS-KLE")
plot3 <- make_plot(df_spnngp,      "Posterior Median GRF -- spNNGP")

grid.arrange(plot1, plot2, plot3, ncol = 3)









#=====================================================================================
# Held-out comparison: SPDE vs regTPS-KLE vs spNNGP
# 80% train / 20% test split of German stations. All three models fit on the
# same training set; predictive performance evaluated at the SAME held-out
# locations using an identical scoring function.
#=====================================================================================

library(loo)
library(scoringRules)   # for CRPS

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
# Generic scoring function -- identical closed-form evaluation for all 3 models,
# given per-draw predictive means (mu_draws: n_draws x n_test) and per-draw
# observation-noise variances (tau2_draws: length n_draws).
#=====================================================================================

logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

score_model <- function(mu_draws, tau2_draws, y_true) {
  
  n_draws <- nrow(mu_draws)
  n_test  <- ncol(mu_draws)
  stopifnot(length(tau2_draws) == n_draws, length(y_true) == n_test)
  
  # -----------------------------------
  # Log predictive density
  # -----------------------------------
  lpd_per_point <- numeric(n_test)
  
  for (i in seq_len(n_test)) {
    logdens_s <- dnorm(y_true[i], mean = mu_draws[, i], sd = sqrt(tau2_draws), log = TRUE)
    
    lpd_per_point[i] <- logSumExp(logdens_s) - log(n_draws)
  }
  lpd_total <- sum(lpd_per_point)
  
  # -----------------------------------
  # RMSE / MAE
  # -----------------------------------
  mu_mean <- colMeans(mu_draws)
  rmse <- sqrt(mean((mu_mean - y_true)^2))
  mae <- mean(abs(mu_mean - y_true))
  
  # -----------------------------------
  # CRPS
  # -----------------------------------
  noise_draws <- matrix(rnorm(n_draws * n_test), nrow = n_draws, ncol = n_test)
  noise_draws <- sweep(noise_draws, MARGIN = 1, STATS = sqrt(tau2_draws), FUN = "*")
  y_rep <- mu_draws + noise_draws
  
  crps_per_point <- sapply(seq_len(n_test),
    function(i)
      crps_sample(y_true[i], y_rep[, i]))
  
  crps_mean <- mean(crps_per_point)
  
  list(lpd_total = lpd_total, lpd_per_point = lpd_per_point,
    rmse = rmse, mae = mae, crps_mean = crps_mean)
}

#=====================================================================================
# --- 1. SPDE: fit on training data, predict at held-out locations ---
#=====================================================================================
cat("\n=== Fitting SPDE on training data ===\n")

dim(train_data)
dim(test_data)


spde_train <- run_tmb_spde(train_data, 100)

spde_mcmc_train <- tmbstan(spde_train$obj, chains = 3, open_progress = FALSE,
                           control = list(max_treedepth = 12, adapt_delta = 0.9),
                           iter = 3000, warmup = 700, cores = no_cores,
                           init = 'last.par.best', seed = 12345)

post_spde <- rstan::extract(spde_mcmc_train)
cat("SPDE posterior parameters:", paste(names(post_spde), collapse = ", "), "\n")

u_tilde_post <- post_spde$u_tilde
n_draws_spde <- nrow(u_tilde_post)

A_test_spde <- as.matrix(inla.spde.make.A(mesh = spde_train$mesh, 
                                          loc = as.matrix(test_data[, c("s1", "s2")])))

dim(A_test_spde)


rho_post <- exp(post_spde$logrho) 
sigma_u_post <- exp(post_spde$logsigma_u) 
sigma_e_post <- exp(post_spde$logsigma_e)

cat("Posterior u_tilde:", dim(post_spde$u_tilde), "\n")
cat("Mesh nodes:", nrow(spde_train$mesh$loc), "\n")
cat("A_test:", dim(A_test_spde), "\n")


mu_spde <- matrix(NA, n_draws_spde, n_test)
for (d in seq_len(n_draws_spde)) {
  kappa_d <- sqrt(8) / rho_post[d]
  tau_field_d <- 1 / (kappa_d * sigma_u_post[d])   # SPDE field-scale tau (NOT observation noise)
  u_d <- u_tilde_post[d, ] / tau_field_d
  mu_spde[d, ] <- as.numeric(A_test_spde %*% u_d)
}


tau2_spde <- sigma_e_post^2   # observation-noise variance draws

score_spde <- score_model(mu_spde, tau2_spde, y_test_sqrt)






#=====================================================================================
# --- 2. regTPS-KLE: fit on training data, predict at held-out locations ---
#=====================================================================================
cat("\n=== Fitting regTPS-KLE on training data ===\n")

tps_train <- regTPS_KLE(sp_data = train_data, dim_grid = 100, n_basis_app = 0.75,
                        variance_threshold = 0.7, 
                        sqrt_transform = TRUE, expand_grid = 0.05,
                        alpha_ref = 1.0,
                        sigma_e_ref = 1.0)


# NOTE: requires `sm` and `evectors` to be included in regTPS_KLE()'s res_list
# (see earlier fix) -- both are needed to predict at NEW, out-of-training locations.
if (is.null(tps_train$sm) || is.null(tps_train$evectors)) {
  stop("regTPS_KLE() must return `sm` and `evectors` in res_list for out-of-sample prediction. ",
       "Add `sm = sm, evectors = evectors` to its return list.")
}

tps_mcmc_train <- tmbstan(tps_train$obj, chains = 3, open_progress = FALSE,
                          control = list(max_treedepth = 12, adapt_delta = 0.9),
                          iter = 3000, warmup = 700, cores = no_cores,
                          init = 'last.par.best', seed = 12345)

post_tps <- rstan::extract(tps_mcmc_train)
z_tilde_post <- post_tps$z_tilde
alpha_post   <- exp(post_tps$logalpha)
sigma_eps_post_tps <- exp(post_tps$logsigma)

M_trunc  <- ncol(z_tilde_post)
M_P_null <- tps_train$M_P_null_space
S_vec    <- tps_train$S_diag_truncated

Phi_test_basis <- mgcv::PredictMat(tps_train$sm, test_data)
Phi_test_kle   <- Phi_test_basis %*% tps_train$evectors[, 1:M_trunc]

n_draws_tps <- nrow(z_tilde_post)
mu_tps <- matrix(NA, n_draws_tps, n_test)
for (d in seq_len(n_draws_tps)) {
  alpha_d <- alpha_post[d]
  scale_d <- rep(1, M_trunc)
  if (M_P_null < M_trunc) {
    k_idx <- (M_P_null + 1):M_trunc
    scale_d[k_idx] <- 1 / sqrt(1 + alpha_d * S_vec[k_idx])
  }
  z_d <- scale_d * z_tilde_post[d, ]
  mu_tps[d, ] <- as.numeric(Phi_test_kle %*% z_d)
}
tau2_tps <- sigma_eps_post_tps^2

score_tps <- score_model(mu_tps, tau2_tps, y_test_sqrt)

#=====================================================================================
# --- 3. spNNGP: fit 3 chains on training data, predict at held-out locations ---
#=====================================================================================
library(spNNGP)
cat("\n=== Fitting spNNGP on training data (3 chains) ===\n")

coords_train <- as.matrix(train_data[, c("s1", "s2")])
y_train_sqrt <- sqrt(train_data$y_obs)

sigma_sq_start <- var(y_train_sqrt) * 0.8
tau_sq_start   <- var(y_train_sqrt) * 0.2
phi_start      <- 3 / (0.3 * diff(range(coords_train[,1])))

starting <- list(phi = phi_start, sigma.sq = sigma_sq_start, tau.sq = tau_sq_start, nu = 0.5)
tuning   <- list(phi = 0.5, sigma.sq = 0.5, tau.sq = 0.5, nu = 0.05)
priors   <- list(phi.Unif    = c(3/(2*diff(range(coords_train[,1]))), 3/(0.01*diff(range(coords_train[,1])))),
                 sigma.sq.IG = c(2, sigma_sq_start),
                 tau.sq.IG   = c(2, tau_sq_start),
                 nu.Unif     = c(0.1, 2.0))

n_chains_spnngp <- 3
spnngp_fits_train <- vector("list", n_chains_spnngp)
for (chain in seq_len(n_chains_spnngp)) {
  set.seed(1000 + chain)
  spnngp_fits_train[[chain]] <- spNNGP(y_train_sqrt ~ 1, coords = coords_train,
                                       method = "latent", family = "gaussian",
                                       n.neighbors = 15, starting = starting, tuning = tuning,
                                       priors = priors, cov.model = "matern",
                                       n.samples = 3000, n.omp.threads = no_cores, n.report = 500)
}

coords_test <- as.matrix(test_data[, c("s1","s2")])
mu_spnngp_list  <- list()
tau2_spnngp_list <- list()

for (chain in seq_len(n_chains_spnngp)) {
  pred_chain <- predict(spnngp_fits_train[[chain]],
                        X.0 = matrix(1, nrow = n_test, ncol = 1),
                        coords.0 = coords_test, n.omp.threads = no_cores)
  
  cat("Chain", chain, "predict() output fields:", paste(names(pred_chain), collapse=", "), "\n")
  if (!("p.w.0" %in% names(pred_chain))) {
    stop("predict.NNGP did not return p.w.0 -- check names(pred_chain) and adapt: ",
         "you need the LATENT process draws at test locations (not just p.y.0) ",
         "to separate the mean prediction from the observation-noise draw, ",
         "for a like-for-like scoring function with SPDE/regTPS-KLE.")
  }
  
  beta_chain <- spnngp_fits_train[[chain]]$p.beta.samples   # n.samples x 1 (intercept)
  w0_chain   <- pred_chain$p.w.0                              # n_test x n.samples
  
  mu_chain <- t(w0_chain) + matrix(beta_chain[,1], nrow = ncol(w0_chain), ncol = n_test)
  mu_spnngp_list[[chain]] <- mu_chain
  
  theta_chain <- spnngp_fits_train[[chain]]$p.theta.samples
  tau_sq_col <- grep("tau\\.sq", colnames(theta_chain))
  tau2_spnngp_list[[chain]] <- theta_chain[, tau_sq_col]
}

mu_spnngp  <- do.call(rbind, mu_spnngp_list)      # combine chains: n_draws_total x n_test
tau2_spnngp <- do.call(c, tau2_spnngp_list)

score_spnngp <- score_model(mu_spnngp, tau2_spnngp, y_test_sqrt)

#=====================================================================================
# --- Final comparison table ---
#=====================================================================================
comparison_table <- data.frame(
  Model = c("SPDE", "regTPS-KLE", "spNNGP"),
  LogPredDensity = c(score_spde$lpd_total, score_tps$lpd_total, score_spnngp$lpd_total),
  RMSE = c(score_spde$rmse, score_tps$rmse, score_spnngp$rmse),
  MAE  = c(score_spde$mae, score_tps$mae, score_spnngp$mae),
  CRPS = c(score_spde$crps_mean, score_tps$crps_mean, score_spnngp$crps_mean)
)

print(comparison_table)


summary(tps_train$S_diag_truncated)

range(tps_train$S_diag_truncated)

quantile(tps_train$S_diag_truncated, probs = c(0, .25, .5, .75, .9, .95, .99, 1))
summary(alpha_post)
quantile(alpha_post, c(.05, .25, .5, .75, .95))









#=====================================================================================
# Spatially blocked cross-validation: SPDE vs regTPS-KLE vs spNNGP
# Instead of a random 80/20 split, partition stations into K spatial blocks
# (via k-means on coordinates) and leave one block out at a time -- a much
# more honest test of spatial extrapolation than random CV.
#=====================================================================================

library(loo)
library(scoringRules)
library(spNNGP)

#=====================================================================================
# 1. Assign spatial blocks via k-means on coordinates
#=====================================================================================

assign_spatial_blocks <- function(coords, n_blocks, seed = 1234) {
  set.seed(seed)
  km <- kmeans(coords, centers = n_blocks, nstart = 25)
  km$cluster
}

n_blocks <- 5   # e.g. 5 spatial folds; adjust based on station density/domain size
coords_all <- as.matrix(sp_data[, c("s1", "s2")])
block_id <- assign_spatial_blocks(coords_all, n_blocks)

cat("Block sizes:\n"); print(table(block_id))

# Quick visual check that blocks are spatially contiguous, not scattered
# plot(coords_all, col = block_id, pch = 19)

#=====================================================================================
# 2. Reusable scoring function (identical to before)
#=====================================================================================


logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

score_model <- function(mu_draws, tau2_draws, y_true) {
  
  n_draws <- nrow(mu_draws)
  n_test  <- ncol(mu_draws)
  stopifnot(length(tau2_draws) == n_draws, length(y_true) == n_test)
  
  # -----------------------------------
  # Log predictive density
  # -----------------------------------
  lpd_per_point <- numeric(n_test)
  
  for (i in seq_len(n_test)) {
    logdens_s <- dnorm(y_true[i], mean = mu_draws[, i], sd = sqrt(tau2_draws), log = TRUE)
    
    lpd_per_point[i] <- logSumExp(logdens_s) - log(n_draws)
  }
  lpd_total <- sum(lpd_per_point)
  
  # -----------------------------------
  # RMSE / MAE
  # -----------------------------------
  mu_mean <- colMeans(mu_draws)
  rmse <- sqrt(mean((mu_mean - y_true)^2))
  mae <- mean(abs(mu_mean - y_true))
  
  # -----------------------------------
  # CRPS
  # -----------------------------------
  noise_draws <- matrix(rnorm(n_draws * n_test), nrow = n_draws, ncol = n_test)
  noise_draws <- sweep(noise_draws, MARGIN = 1, STATS = sqrt(tau2_draws), FUN = "*")
  y_rep <- mu_draws + noise_draws
  
  crps_per_point <- sapply(seq_len(n_test),
                           function(i)
                             crps_sample(y_true[i], y_rep[, i]))
  
  crps_mean <- mean(crps_per_point)
  
  list(lpd_total = lpd_total, lpd_per_point = lpd_per_point,
       rmse = rmse, mae = mae, crps_mean = crps_mean)
}
#=====================================================================================
# 3. Wrap the full three-model fit-and-score pipeline into a reusable function
#=====================================================================================

fit_and_score_all_three <- function(train_data, test_data, no_cores) {
  
  y_test_sqrt <- sqrt(test_data$y_obs)
  n_test <- nrow(test_data)
  
  #========================================
  #--- SPDE ---
  #========================================
  sp_matrix_train <- as.matrix(train_data[, c("s1","s2")])

  spde_train <- run_tmb_spde(train_data, 100)
  spde_mcmc <- tmbstan(spde_train$obj, chains = 3, open_progress = FALSE,
                       control = list(max_treedepth = 12, adapt_delta = 0.9),
                       iter = 3000, warmup = 700, cores = no_cores,
                       init = 'last.par.best', seed = 12345)
  
  post_spde <- rstan::extract(spde_mcmc)
  u_tilde_post <- post_spde$u_tilde
  A_test_spde <- as.matrix(inla.spde.make.A(mesh = spde_train$mesh, 
                                            loc = as.matrix(test_data[, c("s1", "s2")])))
  
  rho_post <- exp(post_spde$logrho) 
  sigma_u_post <- exp(post_spde$logsigma_u) 
  sigma_e_post <- exp(post_spde$logsigma_e)
  
  cat("Posterior u_tilde:", dim(post_spde$u_tilde), "\n")
  cat("Mesh nodes:", nrow(spde_train$mesh$loc), "\n")
  cat("A_test:", dim(A_test_spde), "\n")
  
  n_draws_spde <- nrow(u_tilde_post)
  
  mu_spde <- matrix(NA, n_draws_spde, n_test)
  for (d in seq_len(n_draws_spde)) {
    kappa_d <- sqrt(8) / rho_post[d]
    tau_field_d <- 1 / (kappa_d * sigma_u_post[d])   # SPDE field-scale tau (NOT observation noise)
    u_d <- u_tilde_post[d, ] / tau_field_d
    mu_spde[d, ] <- as.numeric(A_test_spde %*% u_d)
  }
  
  score_spde <- score_model(mu_spde, sigma_e_post^2, y_test_sqrt)
  
  
  #========================================
  #--- regTPS-KLE ---
  #========================================
  tps_train <- regTPS_KLE(sp_data = train_data, dim_grid = 100, n_basis_app = 0.9,
                          variance_threshold = 0.99, sqrt_transform = TRUE, expand_grid = 0.05,
                          alpha_ref = 1.0,
                          sigma_e_ref = 0.5)
  
  tps_mcmc <- tmbstan(tps_train$obj, chains = 3, open_progress = FALSE,
                      control = list(max_treedepth = 12, adapt_delta = 0.9),
                      iter = 3000, warmup = 700, cores = no_cores,
                      init = 'last.par.best', seed = 12345)
  
  post_tps <- rstan::extract(tps_mcmc)
  z_tilde_post <- post_tps$z_tilde
  alpha_post   <- exp(post_tps$logalpha)
  sigma_eps_post_tps <- exp(post_tps$logsigma)
  M_trunc <- ncol(z_tilde_post); M_P_null <- tps_train$M_P_null_space; S_vec <- tps_train$S_diag_truncated
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
  
  #--- spNNGP ---
  coords_train <- as.matrix(train_data[, c("s1","s2")])
  y_train_sqrt <- sqrt(train_data$y_obs)
  sigma_sq_start <- var(y_train_sqrt)*0.8; tau_sq_start <- var(y_train_sqrt)*0.2
  phi_start <- 3/(0.3*diff(range(coords_train[,1])))
  starting <- list(phi = phi_start, sigma.sq = sigma_sq_start, tau.sq   = tau_sq_start,
                   nu = 1.0)
  tuning   <- list(phi=0.5, sigma.sq=0.5, tau.sq=0.5, nu=0)
  priors   <- list(phi.Unif=c(3/(2*diff(range(coords_train[,1]))), 3/(0.01*diff(range(coords_train[,1])))),
                   sigma.sq.IG=c(2,sigma_sq_start), tau.sq.IG=c(2,tau_sq_start), nu.Unif=c(0.1,2.0))
  
  mu_list <- list(); tau2_list <- list()
  for (chain in 1:3) {
    set.seed(1000+chain)
    fit_c <- spNNGP(y_train_sqrt ~ 1, coords=coords_train, method="latent", family="gaussian",
                    n.neighbors=15, starting=starting, tuning=tuning, priors=priors,
                    cov.model="matern", n.samples=2000, n.omp.threads=no_cores, n.report=1000)
    pred_c <- predict(fit_c, X.0=matrix(1,nrow=n_test,ncol=1),
                      coords.0=as.matrix(test_data[,c("s1","s2")]), n.omp.threads=no_cores)
    beta_c <- fit_c$p.beta.samples[,1]
    mu_list[[chain]] <- t(pred_c$p.w.0) + matrix(beta_c, nrow=ncol(pred_c$p.w.0), ncol=n_test)
    theta_c <- fit_c$p.theta.samples
    tau2_list[[chain]] <- theta_c[, grep("tau\\.sq", colnames(theta_c))]
  }
  mu_spnngp <- do.call(rbind, mu_list); tau2_spnngp <- do.call(c, tau2_list)
  score_spnngp <- score_model(mu_spnngp, tau2_spnngp, y_test_sqrt)
  
  data.frame(Model=c("SPDE","regTPS-KLE","spNNGP"),
             LogPredDensity=c(score_spde$lpd_total, score_tps$lpd_total, score_spnngp$lpd_total),
             RMSE=c(score_spde$rmse, score_tps$rmse, score_spnngp$rmse),
             MAE=c(score_spde$mae, score_tps$mae, score_spnngp$mae),
             CRPS=c(score_spde$crps_mean, score_tps$crps_mean, score_spnngp$crps_mean))
}

#=====================================================================================
# 4. Run the spatially blocked CV loop
#=====================================================================================

block_results <- list()

for (b in seq_len(n_blocks)) {
  cat("\n=========================================\n")
  cat("=== Spatial block", b, "of", n_blocks, "held out ===\n")
  cat("=========================================\n")
  
  test_data  <- sp_data[block_id == b, ]
  train_data <- sp_data[block_id != b, ]
  
  block_results[[b]] <- fit_and_score_all_three(train_data, test_data, no_cores)
  block_results[[b]]$block <- b
}

blocked_cv_results <- dplyr::bind_rows(block_results)
print(blocked_cv_results)

#=====================================================================================
# 5. Aggregate across blocks
#=====================================================================================

library(dplyr)
blocked_cv_summary <- blocked_cv_results %>%
  group_by(Model) %>%
  summarise(
    LogPredDensity_total = sum(LogPredDensity),   # log densities: sum across blocks
    RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
    MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
    CRPS_mean = mean(CRPS), CRPS_sd = sd(CRPS),
    .groups = "drop"
  )
print(blocked_cv_summary)













#=====================================================================================
# Scoring function
# Computes:
#   - Log Predictive Density (LPD)
#   - RMSE
#   - MAE
#   - CRPS
#   - R² = squared Pearson correlation between observed and predicted values
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
       
       # Optional: keep vectors for later diagnostics
       mu_mean = mu_mean, lpd_per_point = lpd_per_point,
       crps_per_point = crps_per_point)
}

#=====================================================================================
# Fit-and-score, wrapping your EXACT current settings for all three models.
# Includes an explicit parameter-name check for SPDE, so a naming mismatch
# fails loudly instead of silently. Now reports R2 alongside the other metrics.
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
  
  # Basis at TEST locations, using the TRAINING smooth object's knots/basis --
  # PredictMat does NOT recompute knots, it evaluates the existing basis.
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
# NOTE ON COST: each repetition/fold fits THREE full Bayesian models --
# consider trimming mcmc_iter/spnngp_samples or n_reps if infeasible.
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
# Run it -- pick ONE of the two calls below
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


saveRDS(cv_results, 
        file='C:/Users/jcavi/OneDrive/Escritorio/KLE/real_application/outputs/cv_results.RDS')






#=====================================================================================
# Additional diagnostics: 95% interval coverage, spatial autocorrelation of
# held-out residuals (Moran's I), and a time-vs-accuracy summary plot.
#=====================================================================================

library(ape)      # for Moran.I()
library(ggplot2)

#=====================================================================================
# 1. Extend score_model() to also report 95% predictive interval coverage
#=====================================================================================

score_model <- function(mu_draws, tau2_draws, y_true) {
  n_draws <- nrow(mu_draws); n_test <- ncol(mu_draws)
  stopifnot(length(tau2_draws) == n_draws, length(y_true) == n_test)
  
  lpd_per_point <- numeric(n_test)
  for (i in seq_len(n_test)) {
    logdens_s <- dnorm(y_true[i], mean = mu_draws[, i], sd = sqrt(tau2_draws), log = TRUE)
    lpd_per_point[i] <- logSumExp(logdens_s) - log(n_draws)
  }
  mu_mean <- colMeans(mu_draws)
  noise_draws <- sweep(matrix(rnorm(n_draws * n_test), n_draws, n_test), 1, sqrt(tau2_draws), "*")
  y_rep <- mu_draws + noise_draws
  crps_per_point <- sapply(seq_len(n_test), function(i) crps_sample(y_true[i], y_rep[, i]))
  
  ss_res <- sum((y_true - mu_mean)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  r2_val <- 1 - (ss_res / ss_tot)
  
  # --- 95% predictive interval coverage ---
  q025 <- apply(y_rep, 2, quantile, probs = 0.025)
  q975 <- apply(y_rep, 2, quantile, probs = 0.975)
  covered <- (y_true >= q025) & (y_true <= q975)
  coverage_95 <- mean(covered)
  
  # Residuals for downstream spatial diagnostics (posterior-mean based)
  residuals <- y_true - mu_mean
  
  list(
    lpd_total = sum(lpd_per_point),
    rmse = sqrt(mean((mu_mean - y_true)^2)),
    mae = mean(abs(mu_mean - y_true)),
    crps_mean = mean(crps_per_point),
    r2 = r2_val,
    coverage_95 = coverage_95,
    residuals = residuals   # returned for the Moran's I check below
  )
}

#=====================================================================================
# 2. Moran's I on held-out residuals -- tests for LEFTOVER spatial
#    autocorrelation the model failed to capture. A well-specified spatial
#    model should leave residuals close to spatially uncorrelated (Moran's
#    I near 0, non-significant p-value); a significant positive I indicates
#    the model is missing real spatial structure.
#=====================================================================================

compute_morans_i <- function(residuals, coords) {
  d <- as.matrix(dist(coords))
  # Inverse-distance weights, zero diagonal, avoid division by zero
  w <- 1 / d
  diag(w) <- 0
  w[!is.finite(w)] <- 0
  
  mi <- ape::Moran.I(residuals, w)
  data.frame(observed = mi$observed, expected = mi$expected,
             sd = mi$sd, p_value = mi$p.value)
}

# Example usage after fit_and_score_all_three() has been run for ONE
# train/test split -- residuals are already computed inside score_model(),
# so this just needs the corresponding test coordinates:

# morans_spde <- compute_morans_i(score_spde$residuals, as.matrix(test_data[, c("s1","s2")]))
# morans_tps  <- compute_morans_i(score_tps$residuals,  as.matrix(test_data[, c("s1","s2")]))
# morans_spnngp <- compute_morans_i(score_spnngp$residuals, as.matrix(test_data[, c("s1","s2")]))
#
# morans_table <- dplyr::bind_rows(
#   cbind(Model = "SPDE", morans_spde),
#   cbind(Model = "regTPS-KLE", morans_tps),
#   cbind(Model = "spNNGP", morans_spnngp)
# )
# print(morans_table)

#=====================================================================================
# 3. Time-vs-accuracy trade-off plot, using your existing cv_summary
#=====================================================================================

p_time_vs_rmse <- ggplot(cv_summary, aes(x = Time_mean_min, y = RMSE_mean, color = Model, label = Model)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_sd, ymax = RMSE_mean + RMSE_sd), width = 0) +
  geom_errorbarh(aes(xmin = Time_mean_min - Time_sd_min, xmax = Time_mean_min + Time_sd_min), height = 0) +
  ggrepel::geom_text_repel(show.legend = FALSE) +
  labs(x = "Mean fitting time per fold (min)", y = "Mean RMSE across folds",
       title = "Computational cost vs. predictive accuracy") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

print(p_time_vs_rmse)

#=====================================================================================
# 4. Updated coverage summary across folds (add to your existing cv_results
#    pipeline: extend fit_and_score_all_three()'s output data.frame with a
#    Coverage95 column, analogous to how R2/Time_min were added)
#=====================================================================================

# Inside fit_and_score_all_three()'s final data.frame(...) call, add:
#   Coverage95 = c(score_spde$coverage_95, score_tps$coverage_95, score_spnngp$coverage_95)
#
# Then in cv_summary's summarise(...), add:
#   Coverage95_mean = mean(Coverage95), Coverage95_sd = sd(Coverage95)
#
# A well-calibrated model should show Coverage95_mean close to 0.95.






