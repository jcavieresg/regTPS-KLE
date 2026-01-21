setwd("C:/Users/jcavi/OneDrive/Escritorio/KLE")
rm(list = ls())

options(scipen = 999)

library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA, geigen, fmesher)


# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  


# #==================================
# # Compile the model and load it
compile("spde.cpp", clean = TRUE)
dyn.load(dynlib("spde"))


#==================================
# Compile the model and load it
compile("tps_kle.cpp", clean = TRUE)
dyn.load(dynlib("tps_kle"))


#=====================================================================================
#                               Main Functions
#=====================================================================================

#=========================
# Run SPDE models
#=========================

run_tmb_spde <- function(N_sp, dim_grid, sp_points, mesh, y_obs, u_true, u_grid, Cov_true) {

  # Convert sp_points to matrix
  sp_matrix <- as.matrix(sp_points)

  # Create grid and projection matrices
  A_obs  <- inla.spde.make.A(mesh = mesh, loc = sp_matrix)
  A_grid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(expand.grid(s1 = seq(0, 1, length.out = dim_grid),
                                                                      s2 = seq(0, 1, length.out = dim_grid))))

  # Set up SPDE model
  spde <- inla.spde2.matern(mesh, alpha = 2)
  spde_mat <- spde$param.inla[c("M0", "M1", "M2")]

  
  
  #=================================
  # TMB data
  #=================================
  tmb_data <- list(y = y_obs,
                   A_obs = A_obs,
                   A_grid = A_grid,
                   spde_mat = spde_mat,
                   lambda_rho = -log(0.05)*0.05,
                   lambda_sigma_u = -log(0.05)/5)

  # TMB parameters
  tmb_par <- list(sigma = 0.1,
                  rho   = 0.1,
                  sigma_u = 0.1,
                  u_tilde = rnorm(mesh$n, 0, 1))

  # Make TMB function object and run optimization
  obj_spde <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "spde", random = "u_tilde")
  
  # Run optimization with better settings
  lwr <- c(1e-4, 1e-4, 1e-4)  # Avoid exact zero
  upr <- c(10, 10, 10)        # Add upper bounds to prevent explosion
  
  opt_spde <- try(nlminb(obj_spde$par, obj_spde$fn, obj_spde$gr, 
                         lower = lwr, upper = upr,
                         control = list(eval.max = 2000, iter.max = 1000, 
                                        trace = 0, abs.tol = 1e-6)), 
                  silent = TRUE)
  
  if(inherits(opt_spde, "try-error") || opt_spde$convergence != 0) {
    warning("SPDE optimization had issues: ", 
            if(inherits(opt_spde, "try-error")) "error" else opt_spde$message)
    
    # Try with more robust settings
    opt_spde <- nlminb(obj_spde$par, obj_spde$fn, obj_spde$gr, 
                       lower = lwr, upper = upr,
                       control = list(eval.max = 3000, iter.max = 1500,
                                      step.min = 1e-8, rel.tol = 1e-5))
  }
  
  # sdreport with error handling
  rep_spde <- try(sdreport(obj_spde), silent = TRUE)
  if(inherits(rep_spde, "try-error")) {
    warning("sdreport failed for SPDE model")
    rep_spde <- NULL
  }

  # Return results
  res_list <- list(obj = obj_spde, opt = opt_spde, rep = rep_spde, tmb_data = tmb_data, tmb_par = tmb_par, 
                   sp_points = sp_points,
                   mesh = mesh, spde = spde, u_true_sp = u_true_sp, u_true_grid = u_true_grid, Cov_true = Cov_true)
  return(res_list)
}



#=========================
# Run regTPS-KLE models
#=========================

run_tmb_tps <- function(N_sp, dim_grid, sp_points, mesh, y_obs, u_true_sp, u_true_grid, k_basis, Cov_true,
                        variance_threshold = 0.95) {
  
  n_nodes <- mesh$n
  # Setup Basis and Penalty
  data_smooth <- data.frame(s1 = sp_points$s1, s2 = sp_points$s2, y_obs = y_obs)
  
  sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth, absorb.cons = FALSE)[[1]]
  gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth)
  
  # Get design matrices
  Phi_basis_sp <- predict(gam_fit, newdata = sp_points, type = "lpmatrix")
  Phi_basis_grid <- predict(gam_fit, newdata = expand.grid(
    s1 = seq(0, 1, length.out = dim_grid), 
    s2 = seq(0, 1, length.out = dim_grid)), type = "lpmatrix")
  
  #========================
  # Get penalty matrix S
  S <- sm$S[[1]]
  
  #====================================
  # STANDARD EIGENDECOMPOSITION
  #====================================
  
  cat("  Using standard eigenvalue problem (S ψ = v ψ)\n")
  
  # Standard eigendecomposition: S ψ = v ψ
  S_eig <- eigen(S, symmetric = TRUE)
  S_diag <- S_eig$values
  evectors <- S_eig$vectors
  
  #====================================
  # CLEAN AND ORDER EIGENVALUES
  #====================================
  
  # Clean small eigenvalues
  S_diag[abs(S_diag) < 1e-12] <- 0
  
  # Check for negative eigenvalues (numerical error)
  if(any(S_diag < -1e-10)) {
    warning("Negative eigenvalues detected (max magnitude: ", 
            max(abs(S_diag[S_diag < 0])), ")")
  }
  S_diag[S_diag < 0] <- 0
  
  # Order by INCREASING eigenvalue (smooth to rough)
  # Smallest eigenvalues = null space (v_k ≈ 0)
  # Larger eigenvalues = rougher functions (v_k > 0)
  order_idx <- order(S_diag, decreasing = FALSE)
  S_diag <- S_diag[order_idx]
  evectors <- evectors[, order_idx]
  
  M_P_null_space <- sm$null.space.dim
  
  cat("  Null space dimension:", M_P_null_space, "\n")
  cat("  Total basis functions:", length(S_diag), "\n")
  cat("  Non-zero eigenvalues:", sum(S_diag > 1e-10), "\n")
  
  #====================================
  # VARIANCE-BASED TRUNCATION
  #====================================
  
  cat("  Computing variance-based truncation...\n")
  
  # Estimate alpha from data
  signal_var_est <- max(var(y_obs) - sigma0_error^2, 0.1)
  nonzero_eigs <- S_diag[S_diag > 1e-10]
  
  if(length(nonzero_eigs) == 0) {
    warning("No non-zero eigenvalues found!")
    alpha_est <- 1.0
  } else {
    alpha_est <- signal_var_est / median(nonzero_eigs)
    alpha_est <- max(min(alpha_est, exp(5)), exp(-5))
  }
  
  cat("  Estimated alpha for truncation:", round(alpha_est, 4), "\n")
  
  # Compute KLE eigenvalues: lambda_k = 1/(1 + alpha * v_k)
  lambda_k <- 1 / (1 + alpha_est * S_diag)
  
  # For null space (v_k ≈ 0), lambda_k ≈ 1
  lambda_k[1:M_P_null_space] <- 1.0
  
  # Cumulative variance explained
  total_variance <- sum(lambda_k)
  cumvar <- cumsum(lambda_k) / total_variance
  
  # Find truncation point
  M_truncation <- which(cumvar >= variance_threshold)[1]
  
  if(is.na(M_truncation)) {
    warning("Could not find truncation point at ", variance_threshold*100, "% variance")
    M_truncation <- length(S_diag)
  }
  
  cat("  Initial M_truncation (", variance_threshold*100, "% variance):", M_truncation, "\n")
  
  # Apply constraints
  M_truncation <- max(M_truncation, M_P_null_space + 5)  # At least null space + 5
  M_truncation <- min(M_truncation, k_basis, n_nodes)     # At most available
  
  if (M_truncation < M_P_null_space) {
    warning("M_truncation < null space dimension; increasing to match.")
    M_truncation <- M_P_null_space
  }
  
  # Compute actual variance explained
  var_explained <- sum(lambda_k[1:M_truncation]) / total_variance
  
  cat("  Final M_truncation:", M_truncation, "\n")
  cat("  Variance explained:", round(var_explained * 100, 2), "%\n")
  
  #====================================
  # CREATE TRUNCATED MATRICES
  #====================================
  
  Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
  Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
  S_diag_truncated <- S_diag[1:M_truncation]
  
  #====================================
  # INITIALIZATION
  #====================================
  
  # 1. Sigma initialization
  logsigma_init <- log(sigma0_error)
  
  # 2. Alpha initialization
  signal_var <- max(var(y_obs) - sigma0_error^2, 0.1)
  
  # Use eigenvalues from non-null, non-zero components
  valid_idx <- (M_P_null_space + 1):M_truncation
  valid_eigs <- S_diag_truncated[valid_idx]
  valid_eigs <- valid_eigs[valid_eigs > 1e-10]
  
  if(length(valid_eigs) > 0) {
    # Use 25th percentile for robustness
    target_eig <- quantile(valid_eigs, 0.25)
    # For mass matrix: lambda ≈ 1/(alpha*v) for large alpha*v
    logalpha_init <- log(signal_var / target_eig)
  } else {
    logalpha_init <- 0
  }
  
  # Bound initial values
  logalpha_init <- pmax(pmin(logalpha_init, 3), -3)
  
  cat("  Initial sigma:", round(exp(logsigma_init), 4), "\n")
  cat("  Initial alpha:", round(exp(logalpha_init), 4), "\n")
  
  #====================================
  # TMB DATA
  #====================================
  sigma_prior_s0    <- 0.5   # Upper threshold for sigma
  sigma_prior_alpha <- 0.05  # P(sigma > 0.5) = 0.05
  lambda_sigma <- -log(sigma_prior_alpha) / sigma_prior_s0
  
  tmb_data <- list(
    y_obs = y_obs,
    Phi_kle_sp = Phi_kle_sp,
    Phi_kle_grid = Phi_kle_grid,
    S_diag_truncated = S_diag_truncated,
    M_P_null_space = M_P_null_space,
    lambda_sigma = lambda_sigma
  )
  
  #====================================
  # TMB PARAMETERS
  #====================================
  
  tmb_par <- list(
    z_tilde = rep(0, M_truncation),
    logsigma = logsigma_init,
    logalpha = logalpha_init
  )
  
  #====================================
  # FIT MODEL
  #====================================
  cat("  Fitting TMB model...\n")
  
  obj <- MakeADFun(data = tmb_data, parameters = tmb_par, 
                   DLL = "tps_kle", random = "z_tilde")
  
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 1000, iter.max = 500))
  
  if(opt$convergence != 0) {
    warning("TMB optimization did not converge: ", opt$message)
  } else {
    cat("  Optimization converged successfully\n")
  }
  
  rep_tmb <- try(sdreport(obj), silent = TRUE)
  if(inherits(rep_tmb, "try-error")) {
    warning("sdreport failed")
    rep_tmb <- NULL
  }
  
  #====================================
  # RETURN RESULTS
  #====================================
  
  res_list <- list(
    obj = obj, 
    opt = opt, 
    rep_tmb = rep_tmb, 
    tmb_data = tmb_data, 
    tmb_par = tmb_par, 
    M_truncation = M_truncation, 
    variance_explained = var_explained, 
    k_basis = k_basis, 
    n_nodes = n_nodes,
    sp_points = sp_points,
    u_true_sp = u_true_sp, 
    u_true_grid = u_true_grid, 
    S_diag_full = S_diag, 
    S_diag_truncated = S_diag_truncated, 
    lambda_k = lambda_k,
    evectors = evectors,  
    sm = sm,
    Cov_true = Cov_true
  )
  
  return(res_list)
}


#=====================================================================================
#             Main Script - Fair Comparison Strategy
#=====================================================================================
set.seed(1234)
base_N_sp <- 50
n_scenarios <- 4
dim_grid <- 30
sigma_u <- 1.0
sigma0_error <- 0.3
nu <- 1.5
rho <- 0.3

# Matern covariance function used for simulation
matern_cov <- function(coords, nu, rho, sigma2) {
  D <- as.matrix(dist(coords))
  D[D == 0] <- 1e-10
  scaling <- (sqrt(2 * nu) * D) / rho
  matern_part <- (2^(1 - nu)) / gamma(nu) * scaling^nu * besselK(scaling, nu)
  diag(matern_part) <- 1
  return(sigma2 * matern_part)
}


#====================================
# Creating the lists for the models
#====================================
fits_TMB_spde <- list()
fits_TMB_tps <- list()

for (i in 1:n_scenarios) {
  N_sp <- base_N_sp * i
  
  # Create a common mesh and points for this scenario
  sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))
  sp_matrix <- as.matrix(sp_points)
  bound1 <- fmesher::fm_nonconvex_hull(sp_matrix)
  mesh <- fmesher::fm_rcdt_2d_inla(loc = sp_matrix, boundary = bound1, refine = FALSE, plot.delay = NULL)
  
  # Create a common grid for this scenario
  grid_total <- expand.grid(s1 = seq(0, 1, length.out = dim_grid), s2 = seq(0, 1, length.out = dim_grid))
  
  # Simulate the TRUE latent field on the mesh nodes ONCE
  Cov_true <- matern_cov(mesh$loc, nu = nu, rho = rho, sigma2 = sigma_u^2)
  u_true <- as.numeric(MASS::mvrnorm(1, mu = rep(0, mesh$n), Sigma = Cov_true))
  
  # Create projection matrices to get the true field at other locations
  A_obs_proj <- inla.spde.make.A(mesh = mesh, loc = sp_matrix)
  A_grid_proj <- inla.spde.make.A(mesh = mesh, loc = as.matrix(grid_total))
  
  # Project the true field to observation points and grid points
  u_true_sp <- as.numeric(A_obs_proj %*% u_true)
  u_true_grid <- as.numeric(A_grid_proj %*% u_true)
  
  # Add noise to the observations
  y_obs <- u_true_sp + rnorm(N_sp, 0, sigma0_error)
  
  # Run SPDE model first
  cat("\n--- Running SPDE model ---\n")
  obj_spde <- run_tmb_spde(N_sp, dim_grid, sp_points, mesh, y_obs, u_true_sp, u_true_grid, Cov_true)
  
  # Get number of mesh nodes from SPDE
  n_mesh_nodes <- mesh$n
  
  cat("\n=== Comparison Setup ===\n")
  cat("N observations:", N_sp, "\n")
  cat("SPDE mesh nodes:", n_mesh_nodes, "\n")
  
  #====================================
  # COMPARISON STRATEGY FOR regTPS-KLE
  #====================================
  # Maximum k_basis for TPS (must be < N_sp)
  k_basis_max <- floor(0.99 * N_sp)   # 99% of data points
  k_basis_max <- max(k_basis_max, 10) # At least 10 basis functions
  
  if(k_basis_max < n_mesh_nodes) {
    cat("Note: TPS k_basis (", k_basis_max, ") < SPDE nodes (", n_mesh_nodes, ")\n")
    cat("      This is expected - TPS is constrained by data size\n")
    cat("      Variance-based truncation will select optimal M_truncation\n")
  } else {
    cat("TPS k_basis (", k_basis_max, ") >= SPDE nodes (", n_mesh_nodes, ")\n")
  }
  
  # Run TPS model
  cat("\n--- Running TPS model ---\n")
  obj_tps <- run_tmb_tps(N_sp, dim_grid, sp_points, mesh, y_obs, u_true_sp, u_true_grid, 
                         k_basis = k_basis_max, 
                         Cov_true, 
                         variance_threshold = 0.99)
  
  #====================================
  # COMPARISON SUMMARY
  #====================================
  
  cat("\n=== Final Comparison ===\n")
  cat("SPDE:\n")
  cat("  Basis functions used:", n_mesh_nodes, "\n")
  cat("  (all mesh nodes)\n")
  
  cat("\nTPS:\n")
  cat("  Available k_basis:", k_basis_max, "\n")
  cat("  Selected M_truncation:", obj_tps$M_truncation, "\n")
  cat("  Variance explained:", round(obj_tps$variance_explained * 100, 2), "%\n")
  
  cat("\nEfficiency:\n")
  cat("  TPS uses", obj_tps$M_truncation, "basis functions vs SPDE's", n_mesh_nodes, "\n")
  cat("  Ratio (TPS/SPDE):", round(obj_tps$M_truncation / n_mesh_nodes, 3), "\n")
  
  if(obj_tps$M_truncation < n_mesh_nodes) {
    reduction_pct <- round((1 - obj_tps$M_truncation / n_mesh_nodes) * 100, 1)
    cat("  TPS achieves", reduction_pct, "% reduction in basis functions\n")
    cat("  while maintaining", round(obj_tps$variance_explained * 100, 1), "% variance\n")
  }
  
  fits_TMB_spde[[i]] <- obj_spde
  fits_TMB_tps[[i]] <- obj_tps
}


#====================================
# Saving all the TMB models
#====================================
fits_TMB_spde <- list(fits_TMB_spde[[1]], fits_TMB_spde[[2]], fits_TMB_spde[[3]], fits_TMB_spde[[4]])
saveRDS(fits_TMB_spde, file = 'outputs/fits_TMB_spde.RDS')

fits_TMB_tps <- list(fits_TMB_tps[[1]], fits_TMB_tps[[2]], fits_TMB_tps[[3]], fits_TMB_tps[[4]])
saveRDS(fits_TMB_tps, file = 'outputs/fits_TMB_tps.RDS')




#============================
# SPDE - Stan models
#============================
M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n50"
M[[2]] = list()
M[[2]]$model = "spatial_n100"
M[[3]] = list()
M[[3]]$model = "spatial_n150"
M[[4]] = list()
M[[4]]$model = "spatial_n200"

M[[1]]$formula = fits_TMB_spde[[1]]$obj
M[[2]]$formula = fits_TMB_spde[[2]]$obj
M[[3]]$formula = fits_TMB_spde[[3]]$obj
M[[4]]$formula = fits_TMB_spde[[4]]$obj



#=================
# Run the models
lwr <- c(1e-4, 1e-4, 1e-4)  
upr <- c(10, 10, 10)        

for (i in 1:length(M)){
  startTime <- Sys.time()
  print(paste("Running:  ", M[[i]]$model))
  fit <- tmbstan(M[[i]]$formula,
                 chains= 3, open_progress = FALSE,
                 control = list(max_treedepth= 12,  adapt_delta = 0.9),
                 iter = 3000, warmup= 700, cores=no_cores,
                 lower = lwr, upper = upr, seed = 12345)
                 # init = 'last.par.best', seed = 12345)
  endTime <- Sys.time()
  timeUsed = difftime(endTime, startTime, units='mins')
  print(timeUsed)
  saveRDS(fit, file=paste0('outputs/stan_spde_', i,'.RDS'))
}







#============================
# regTPS-KLE - Stan models
#============================
M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n50"
M[[2]] = list()
M[[2]]$model = "spatial_n100"
M[[3]] = list()
M[[3]]$model = "spatial_n150"
M[[4]] = list()
M[[4]]$model = "spatial_n200"

M[[1]]$formula = fits_TMB_tps[[1]]$obj
M[[2]]$formula = fits_TMB_tps[[2]]$obj
M[[3]]$formula = fits_TMB_tps[[3]]$obj
M[[4]]$formula = fits_TMB_tps[[4]]$obj

#=================
# Run the models

for (i in 1:length(M)){
  startTime <- Sys.time()
  print(paste("Running:  ", M[[i]]$model))
  fit <- tmbstan(M[[i]]$formula,
                 chains= 3, open_progress = FALSE,
                 control = list(max_treedepth= 12,  adapt_delta = 0.9),
                 iter = 3000, warmup= 700, cores=no_cores,
                 init = 'last.par.best', seed = 12345)
  endTime <- Sys.time()
  timeUsed = difftime(endTime, startTime, units='mins')
  print(timeUsed)
  saveRDS(fit, file=paste0('outputs/stan_tps_', i,'.RDS'))
}








