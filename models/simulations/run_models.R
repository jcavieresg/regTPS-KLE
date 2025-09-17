setwd("C:/Users/Usuario/Desktop/KLE")
rm(list = ls())

options(scipen = 999)

library(pacman)
pacman::p_load(tidyverse, dplyr, parallel, ggplot2,
               TMB, tmbstan, mgcv, MASS, INLA)


# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  


#==================================
# Compile the model and load it
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
  # This function now takes the simulated data as input, ensuring consistency.
  
  # Convert sp_points to matrix
  sp_matrix <- as.matrix(sp_points)
  
  # Create grid and projection matrices
  A_obs  <- inla.spde.make.A(mesh = mesh, loc = sp_matrix)
  A_grid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(expand.grid(s1 = seq(0, 1, length.out = dim_grid), 
                                                                      s2 = seq(0, 1, length.out = dim_grid))))
  
  # Set up SPDE model
  spde <- inla.spde2.matern(mesh, alpha = 2)
  spde_mat <- spde$param.inla[c("M0", "M1", "M2")]
  
  # TMB data
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
                  u = rnorm(mesh$n, 0, 1))
  
  # Make TMB function object and run optimization
  obj_spde <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "spde", random = "u")
  lwr <- c(0, 0, 0)
  upr <- c(Inf, Inf, Inf)
  opt_spde = nlminb(obj_spde$par, obj_spde$fn, obj_spde$gr, lower = lwr, upper = upr)
  rep_spde <- sdreport(obj_spde)
  
  # Return results
  res_list <- list(obj = obj_spde, opt = opt_spde, rep = rep_spde, tmb_data = tmb_data, tmb_par = tmb_par, mesh = mesh, 
                   spde = spde,
                   u_true_sp = u_true_sp, 
                   u_true_grid = u_true_grid,
                   Cov_true = Cov_true)
  return(res_list)
}



#=========================
# Run regTPS-KLE models
#=========================
run_tmb_tps <- function(N_sp, dim_grid, sp_points, mesh, y_obs, u_true_sp, u_true_grid, k_basis, Cov_true) {

  n_nodes <- mesh$n
  # Setup Basis and Penalty
  data_smooth <- data.frame(s1 = sp_points$s1, s2 = sp_points$s2, y_obs = y_obs)
  
  sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth, absorb.cons = FALSE)[[1]]
  gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth)
  
  # Get design matrices
  Phi_basis_sp <- predict(gam_fit, newdata = sp_points, type = "lpmatrix")
  Phi_basis_grid <- predict(gam_fit, newdata = expand.grid(s1 = seq(0, 1, length.out = dim_grid), s2 = seq(0, 1, length.out = dim_grid)), type = "lpmatrix")
  
  # Get penalty matrix S and its eigen-decomposition
  S <- sm$S[[1]]
  S_eig <- eigen(S, symmetric = TRUE)
  S_diag <- S_eig$values
  S_diag[abs(S_diag) < 1e-12] <- 0
  evectors <- S_eig$vectors
  order_idx <- order(S_diag, decreasing = TRUE)
  # order_idx <- order(S_diag, decreasing = FALSE)
  S_diag <- S_diag[order_idx]
  evectors <- evectors[, order_idx]
  M_P_null_space <- sm$null.space.dim
  M_truncation <- min(k_basis, n_nodes)
  
  if (M_truncation < M_P_null_space) {
    warning("M_truncation < null space dimension; increasing to match.")
    M_truncation <- M_P_null_space
  }
  
  Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
  Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]
  S_diag_truncated <- S_diag[1:M_truncation]
  
 
  # TMB data and parameters
  tmb_data <- list(y_obs = y_obs,
                   Phi_kle_sp = Phi_kle_sp,
                   Phi_kle_grid = Phi_kle_grid,
                   S_diag_truncated = S_diag_truncated,
                   M_P_null_space = M_P_null_space)
  
  logsigma_init <- log(sd(y_obs))
  logalpha_init <- log(var(y_obs) / mean(S_diag_truncated))
  
  tmb_par <- list(Z = rep(0, M_truncation),
                  logsigma = logsigma_init,
                  logalpha = logalpha_init)
  
  obj <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "tps_kle", random = "Z")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep_tmb <- sdreport(obj)
  
  res_list = list(obj = obj, opt = opt, rep_tmb = rep_tmb, tmb_data = tmb_data, tmb_par = tmb_par, 
                  M_truncation = M_truncation, 
                  n_nodes = n_nodes,
                  u_true_sp = u_true_sp, 
                  u_true_grid = u_true_grid, 
                  S_diag_full = S_diag, 
                  S_diag_truncated = S_diag_truncated, 
                  evectors = evectors, 
                  sm = sm,
                  Cov_true = Cov_true)
  return(res_list)
}

#=====================================================================================
#             Main Script with a unified data simulation
#=====================================================================================

set.seed(1234) # Set the seed once for the entire simulation
base_N_sp <- 50
n_scenarios <- 4
dim_grid <- 30
sigma_u <- 1.0
sigma0_error <- 0.3
nu <- 1.5
rho <- 0.3

# Define a single Matern covariance function to use for simulation
matern_cov <- function(coords, nu, rho, sigma2) {
  D <- as.matrix(dist(coords))
  D[D == 0] <- 1e-10
  scaling <- (sqrt(2 * nu) * D) / rho
  matern_part <- (2^(1 - nu)) / gamma(nu) * scaling^nu * besselK(scaling, nu)
  diag(matern_part) <- 1
  return(sigma2 * matern_part)
}

fits_TMB_spde <- list()
fits_TMB_tps <- list()

for (i in 1:n_scenarios) {
  N_sp <- base_N_sp * i
  
  # Create a common mesh and points for this scenario
  sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))
  sp_matrix <- as.matrix(sp_points)
  bound1 <- inla.nonconvex.hull(sp_matrix)
  mesh <- inla.mesh.create(loc = sp_matrix, boundary = bound1, refine = FALSE, plot.delay = NULL)
  
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
  
  # Run both models with the SAME simulated data
  obj_spde <- run_tmb_spde(N_sp, dim_grid, sp_points, mesh, y_obs, u_true_sp, u_true_grid, Cov_true)
  obj_tps <- run_tmb_tps(N_sp, dim_grid, sp_points, mesh, y_obs, u_true_sp, u_true_grid, k_basis = floor(0.95*N_sp), Cov_true)
  
  fits_TMB_spde[[i]] <- obj_spde
  fits_TMB_tps[[i]] <- obj_tps
}

# Evaluate if the simulated data are the same
y_obs_tps   <- fits_TMB_tps[[1]][[4]]$y_obs
field_tps_sp   <- fits_TMB_tps[[1]][[9]]
field_tps_grid <- fits_TMB_tps[[1]][[8]]

y_obs_spde  <- fits_TMB_spde[[1]][[4]]$y
field_spde_sp  <- fits_TMB_spde[[1]][[9]]
field_spde_grid  <- fits_TMB_spde[[1]][[8]]

cat("TPS vs SPDE  (y_obs):        ", all.equal(y_obs_tps, y_obs_spde), "\n")
cat("TPS vs SPDE  (field_sp):     ", all.equal(field_tps_sp, field_spde_sp), "\n")
cat("TPS vs SPDE  (field_grid):   ", all.equal(field_tps_grid, field_spde_grid), "\n")

# Save results
saveRDS(fits_TMB_spde, file='outputs/fits_TMB_spde.RDS')
saveRDS(fits_TMB_tps, file='outputs/fits_TMB_tps.RDS')






# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
par(mfrow = c(4, 3))
image(matrix(fits_TMB_tps[[1]]$u_true_grid, 30, 30), main = "True GRF 1",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_spde[[1]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (spde=", fits_TMB_spde[[1]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_tps[[1]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (M_KLE=", fits_TMB_tps[[1]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(fits_TMB_tps[[2]]$u_true_grid, 30, 30), main = "True GRF 2",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_spde[[2]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (spde=", fits_TMB_spde[[2]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_tps[[2]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (M_KLE=", fits_TMB_tps[[2]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(fits_TMB_tps[[3]]$u_true_grid, 30, 30), main = "True GRF 3",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_spde[[3]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (spde=", fits_TMB_spde[[3]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_tps[[3]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (M_KLE=", fits_TMB_tps[[3]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(fits_TMB_tps[[4]]$u_true_grid, 30, 30), main = "True GRF 4",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_spde[[4]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (spde=", fits_TMB_spde[[4]][[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(fits_TMB_tps[[4]][[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (M_KLE=", fits_TMB_tps[[4]][[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
