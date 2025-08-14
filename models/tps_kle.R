setwd("C:/Users/cavieresgaet/Desktop/kle")
rm(list = ls())

options(scipen = 999)
library(TMB)
library(tmbstan)
library(mgcv)
library(MASS)
library(INLA)



# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  

#==================================
# Compile the model and load it
compile("kle_model_alpha_M_opt.cpp")
dyn.load(dynlib("kle_model_alpha_M_opt"))




#=====================================================================================
#                Main function Regularized TPS Kernel for KLE
#=====================================================================================

run_tmb <- function(N_sp = N_sp, dim_grid = dim_grid, sp_points = sp_points, mesh = NULL){
set.seed(1234)

# Get triangle count from mesh
# n_triangles <- nrow(mesh$graph$tv)
n_nodes   <- mesh$n 
N_sp <- N_sp
k_basis <- floor(0.95 * N_sp)
sigma0_error <- 0.1

set.seed(1234)
# sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))

matern_cov <- function(coords, nu = 1.5, rho = 0.2, sigma2 = 1.0) {
  D <- as.matrix(dist(coords))
  D[D == 0] <- 1e-10
  scaling <- (sqrt(2 * nu) * D) / rho
  matern_part <- (2^(1 - nu)) / gamma(nu) * scaling^nu * besselK(scaling, nu)
  diag(matern_part) <- 1
  return(sigma2 * matern_part)
}

Cov_true_obs <- matern_cov(sp_points)
true_field_obs <- mvrnorm(1, mu = rep(0, N_sp), Sigma = Cov_true_obs)
y_obs <- true_field_obs + rnorm(N_sp, 0, sigma0_error)

dim_grid <- dim_grid
s1_grid <- seq(0, 1, length.out = dim_grid)
s2_grid <- seq(0, 1, length.out = dim_grid)
grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)


# Simulate true GRF on grid (same kernel as for observations)
Cov_true_grid <- matern_cov(grid_total)
true_field_grid <- MASS::mvrnorm(1, mu = rep(0, nrow(grid_total)), Sigma = Cov_true_grid)


#=====================================================================================
# 3. Setup Basis and Penalty (Optimized for TMB)
#=====================================================================================
data_smooth <- data.frame(s1 = sp_points$s1, s2 = sp_points$s2, y_obs = y_obs)

sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth, absorb.cons = FALSE)[[1]]
gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth)

# Get design matrices
Phi_basis_sp <- predict(gam_fit, newdata = sp_points, type = "lpmatrix")
Phi_basis_grid <- predict(gam_fit, newdata = grid_total, type = "lpmatrix")
n_basis_grid <- ncol(Phi_basis_grid)

# Get penalty matrix S and its eigen-decomposition
S_grid <- sm$S[[1]]
S_eig <- eigen(S_grid, symmetric = TRUE)
S_diag <- S_eig$values
evectors <- S_eig$vectors
order_idx <- order(S_diag, decreasing = TRUE)
S_diag <- S_diag[order_idx]
evectors <- evectors[, order_idx]

M_P_null_space <- sm$null.space.dim

# Choose M_truncation based on mesh
M_truncation <- min(k_basis, n_nodes)
if (M_truncation < M_P_null_space) {
  warning("M_truncation < null space dimension; increasing to match.")
  M_truncation <- M_P_null_space
}

# --- THE CRUCIAL OPTIMIZATION STEP ---
# Pre-compute the KLE basis matrices in R!
Phi_kle_sp <- Phi_basis_sp %*% evectors[, 1:M_truncation]
Phi_kle_grid <- Phi_basis_grid %*% evectors[, 1:M_truncation]

# Also truncate S_diag for the prior
S_diag_truncated <- S_diag[1:M_truncation]

cat(paste0("Pre-computed KLE basis matrices of size ", M_truncation, " columns.\n"))

#=====================================================================================
# 4. Prepare TMB input and run the optimized model
#=====================================================================================
tmb_data <- list(
  y = y_obs,
  Phi_kle_sp = Phi_kle_sp,
  Phi_kle_grid = Phi_kle_grid,
  S_diag_truncated = S_diag_truncated,
  M_P_null_space = M_P_null_space
)

tmb_par <- list(
  Z = rep(0, M_truncation),
  logsigma = log(0.5),
  logalpha = log(1.0)
)

obj <- MakeADFun(data = tmb_data, parameters = tmb_par,
                           DLL = "kle_model_alpha_M_opt", random = "Z")

# The optimization should now be significantly faster
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep_tmb <- sdreport(obj)
res_list = list(obj, opt, rep_tmb, tmb_data, tmb_par, M_truncation, n_nodes, true_field_grid)
return(res_list)
}



# Base parameters
# base_N_sp <- 50
# scenarios <- list()
# 
# # Define custom max.edge values to control number of triangles
# max_edge_list <- list(
#   c(0.16, 0.4),  # Scenario 1 (~100 triangles)
#   c(0.14, 0.35), # Scenario 2 (~130)
#   c(0.12,  0.3),  # Scenario 3 (~160)
#   c(0.11, 0.3)  # Scenario 4 (~200)
# )
# 
# # Define corresponding cutoffs (optional: scale with edge)
# cutoff_list <- c(0.05, 0.04, 0.03, 0.02)
# 
# for (i in 1:4) {
#   N_sp <- base_N_sp * i
#   sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))
#   sp_matrix <- as.matrix(sp_points)
#   
#   # Use pre-tuned mesh parameters
#   max_edge <- max_edge_list[[i]]
#   cutoff <- cutoff_list[i]
#   
#   # Create the mesh
#   mesh <- inla.mesh.2d(loc = sp_matrix, max.edge = max_edge, cutoff = cutoff, offset = c(0.05, 0.1))
#   
#   # Store results
#   scenarios[[paste0("scenario", i)]] <- list(N_sp = N_sp, points = sp_points, mesh = mesh,
#                                         n_triangles = nrow(mesh$graph$tv), n_vertices = mesh$n)}



#=====================================================================================
#         Function to create a list of scenarios using inla.mesh.create
#=====================================================================================

mesh_scenarios <- function(base_N_sp = 50, n_scenarios = 4) {
  scenarios <- list()
  for (i in 1:n_scenarios) {
    set.seed(i)
    N_sp <- base_N_sp * i
    sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))
    sp_matrix <- as.matrix(sp_points)
    bound1 <- inla.nonconvex.hull(sp_matrix)
    mesh <- inla.mesh.create(loc = sp_matrix, boundary = bound1, refine = FALSE, plot.delay = NULL)
    
    # Store the results for this scenario in the list
    scenario_name <- paste0("scenario", i)
    scenarios[[scenario_name]] <- list(N_sp = N_sp, sp_points = sp_points,
      mesh = mesh, n_triangles = nrow(mesh$graph$tv), # Number of triangles in the mesh
      n_nodes = mesh$n # Number of vertices (nodes) in the mesh
    )
    
    # Print the number of vertices for this scenario for a quick check
    cat(paste("Scenario", i, "created with", mesh$n, "vertices.\n"))
  }
  return(scenarios)}


#=====================================================================================
#                    Run the function to generate the scenarios
#=====================================================================================

# Define the base number of points and the number of scenarios
base_N_sp <- 50
n_scenarios <- 4

# Run the function to generate the scenarios
scenarios <- mesh_scenarios(base_N_sp = base_N_sp, n_scenarios = n_scenarios)

plot(scenarios$scenario1$mesh); points(scenarios$scenario1$sp_points, col = "red")
plot(scenarios$scenario2$mesh); points(scenarios$scenario2$sp_points, col = "red")
plot(scenarios$scenario3$mesh); points(scenarios$scenario3$sp_points, col = "red")
plot(scenarios$scenario4$mesh); points(scenarios$scenario4$sp_points, col = "red")


# 
scenarios$scenario1$mesh$n
scenarios$scenario2$mesh$n
scenarios$scenario3$mesh$n
scenarios$scenario4$mesh$n


obj1_tps <- run_tmb(N_sp = scenarios$scenario1$N_sp, dim_grid = 30, sp_points = scenarios$scenario1$sp_points, mesh = scenarios$scenario1$mesh)
obj2_tps <- run_tmb(N_sp = scenarios$scenario2$N_sp, dim_grid = 30, sp_points = scenarios$scenario2$sp_points, mesh = scenarios$scenario2$mesh)
obj3_tps <- run_tmb(N_sp = scenarios$scenario3$N_sp, dim_grid = 30, sp_points = scenarios$scenario3$sp_points, mesh = scenarios$scenario3$mesh)
obj4_tps <- run_tmb(N_sp = scenarios$scenario4$N_sp, dim_grid = 30, sp_points = scenarios$scenario4$sp_points, mesh = scenarios$scenario4$mesh)

obj1[[6]]


M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n50"
M[[2]] = list()
M[[2]]$model = "spatial_n100"
M[[3]] = list()
M[[3]]$model = "spatial_n150"
M[[4]] = list()
M[[4]]$model = "spatial_n200"


M[[1]]$formula = obj1_tps[[1]]
M[[2]]$formula = obj2_tps[[1]]
M[[3]]$formula = obj3_tps[[1]]
M[[4]]$formula = obj4_tps[[1]]




#===========================================
#           Run the models
#===========================================

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
  saveRDS(fit, file=paste0('stan_tps_', i,'.RDS'))
}


#==================================
# Extract and Evaluate Results
#==================================
Z_est <- as.vector(rep_tmb$par.random[names(rep_tmb$par.random) == "Z"])
sigma_est <- exp(rep_tmb$par.fixed["logsigma"])
alpha_est <- exp(rep_tmb$par.fixed["logalpha"])

# Get the approximated GRF on the full grid from TMB's REPORT
tps_field_grid <- as.vector(obj$report()$field_grid)

message(paste("Estimated sigma:", round(sigma_est, 3)))
message(paste("Estimated alpha:", round(alpha_est, 3)))

# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
# par(mfrow = c(4, 2))
# image(matrix(obj1[[8]], 30, 30), main = "True GRF",
#       col = hcl.colors(100, "viridis"), asp = 1)
# image(matrix(obj1[[1]]$report()$field_grid, 30, 30), main = paste0("Estimated GRF (spde=", obj1[[4]]$mesh$n, ")"),
#       col = hcl.colors(100, "viridis"), asp = 1)
# image(matrix(tps_field_grid, dim_grid, dim_grid), main = paste0("Estimated GRF (M_KLE=", M_truncation, ")"),
#       col = hcl.colors(100, "viridis"), asp = 1)




Z_est <- as.vector(obj1_tps[[3]]$par.random[names(obj1_tps[[3]]$par.random) == "Z"])
sigma_est <- exp(obj1_tps[[3]]$par.fixed["logsigma"])
alpha_est <- exp(obj1_tps[[3]]$par.fixed["logalpha"])

# Get the approximated GRF on the full grid from TMB's REPORT
tps_field_grid <- as.vector(obj1_tps[[1]]$report()$field_grid)

message(paste("Estimated sigma:", round(sigma_est, 3)))
message(paste("Estimated alpha:", round(alpha_est, 3)))




mcmc_tps1 <- readRDS('stan_tps_1.RDS')
mcmc_tps2 <- readRDS('stan_tps_2.RDS')
mcmc_tps3 <- readRDS('stan_tps_3.RDS')
mcmc_tps4 <- readRDS('stan_tps_4.RDS')


tps_post_grid1 <- extract(mcmc_tps1)$Z  # Posterior samples of KL coefficients
tps_mean_grid1 <- colMeans(tps_post_grid1)
tps_field_grid1 <- as.vector(obj1_tps[[4]]$Phi_kle_grid %*% tps_mean_grid1)

tps_post_grid2 <- extract(mcmc_tps2)$Z  # Posterior samples of KL coefficients
tps_mean_grid2 <- colMeans(tps_post_grid2)
tps_field_grid2 <- as.vector(obj2_tps[[4]]$Phi_kle_grid %*% tps_mean_grid2)

tps_post_grid3 <- extract(mcmc_tps3)$Z  # Posterior samples of KL coefficients
tps_mean_grid3 <- colMeans(tps_post_grid3)
tps_field_grid3 <- as.vector(obj3_tps[[4]]$Phi_kle_grid %*% tps_mean_grid3)

tps_post_grid4 <- extract(mcmc_tps4)$Z  # Posterior samples of KL coefficients
tps_mean_grid4 <- colMeans(tps_post_grid4)
tps_field_grid4 <- as.vector(obj4_tps[[4]]$Phi_kle_grid %*% tps_mean_grid4)





# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
par(mfrow = c(4, 3))
image(matrix(obj1_tps[[8]], 30, 30), main = "True GRF 1",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid1, 30, 30), main = paste0("Estimated GRF (spde=", obj1_spde[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid1, 30, 30), main = paste0("Estimated GRF (M_KLE=", obj1_tps[[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(obj2_tps[[8]], 30, 30), main = "True GRF 2",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid1, 30, 30), main = paste0("Estimated GRF (spde=", obj2_spde[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid2, 30, 30), main = paste0("Estimated GRF (M_KLE=", obj2_tps[[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(obj3_tps[[8]], 30, 30), main = "True GRF 3",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid3, 30, 30), main = paste0("Estimated GRF (spde=", obj3_spde[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid3, 30, 30), main = paste0("Estimated GRF (M_KLE=", obj3_tps[[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)



image(matrix(obj4_tps[[8]], 30, 30), main = "True GRF 4",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid4, 30, 30), main = paste0("Estimated GRF (spde=", obj4_spde[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(tps_field_grid4, 30, 30), main = paste0("Estimated GRF (M_KLE=", obj4_tps[[6]], ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


# SPDE AND TPS
# OBJ1
rmse_spde <- sqrt(mean((spde_field_grid1 - obj1_spde[[8]])^2));rmse_spde
r2_spde   <- cor(as.vector(spde_field_grid1), obj1_spde[[8]])^2;r2_spde
mae_spde  <- mean(abs(as.vector(spde_field_grid1) - obj1_spde[[8]]));mae_spde

rmse_tps <- sqrt(mean((tps_field_grid1 - obj1_tps[[8]])^2));rmse_tps
r2_tps   <- cor(as.vector(tps_field_grid1), obj1_tps[[8]])^2;r2_tps
mae_tps  <- mean(abs(tps_field_grid1 - obj1_tps[[8]]));mae_tps


# OBJ2
rmse_spde2 <- sqrt(mean((spde_field_grid2 - obj2_spde[[8]])^2));rmse_spde2
r2_spde2   <- cor(as.vector(spde_field_grid2), obj2_spde[[8]])^2;r2_spde2
mae_spde2  <- mean(abs(as.vector(spde_field_grid2) - obj2_spde[[8]]));mae_spde

rmse_tps2 <- sqrt(mean((tps_field_grid2 - obj2_tps[[8]])^2));rmse_tps2
r2_tps2   <- cor(as.vector(tps_field_grid2), obj2_tps[[8]])^2;r2_tps2
mae_tps2  <- mean(abs(tps_field_grid2 - obj2_tps[[8]]));mae_tps2


# OBJ3
rmse_spde3 <- sqrt(mean((spde_field_grid3 - obj3_spde[[8]])^2));rmse_spde3
r2_spde3   <- cor(as.vector(spde_field_grid3), obj3_spde[[8]])^2;r2_spde3
mae_spde3  <- mean(abs(as.vector(spde_field_grid3) - obj3_spde[[8]]));mae_spde3

rmse_tps3 <- sqrt(mean((tps_field_grid3 - obj3_tps[[8]])^2));rmse_tps3
r2_tps3   <- cor(as.vector(tps_field_grid3), obj3_tps[[8]])^2;r2_tps3
mae_tps3  <- mean(abs(tps_field_grid3 - obj3_tps[[8]]));mae_tps3

# OBJ4
rmse_spde4 <- sqrt(mean((spde_field_grid4 - obj4_spde[[8]])^2));rmse_spde4
r2_spde4   <- cor(as.vector(spde_field_grid4), obj4_spde[[8]])^2;r2_spde4
mae_spde4  <- mean(abs(as.vector(spde_field_grid4) - obj4_spde[[8]]));mae_spde4

rmse_tps4 <- sqrt(mean((tps_field_grid4 - obj4_tps[[8]])^2));rmse_tps4
r2_tps4   <- cor(as.vector(tps_field_grid4), obj4_tps[[8]])^2;r2_tps4
mae_tps4  <- mean(abs(tps_field_grid4 - obj4_tps[[8]]));mae_tps4









# Correlation
# correlation <- cor(true_grf_grid, as.numeric(app_grf_grid))
# message(paste("Correlation:", round(correlation, 3)))

# A quick plot to visualize the decay of the eigenvalues.
# These eigenvalues govern the penalty applied to each KLE coefficient.
# Note: A large jump at the beginning is expected for the unpenalized modes (where S_diag=0).
plot(obj1_tps[[4]]$S_diag_truncated, type = "l", log = "y",
     main = "Decay of KLE Eigenvalues (log scale)",
     xlab = "KLE Basis Index (k)",
     ylab = "Eigenvalue (lambda_k)")

# Optionally, you can plot the prior variance of the coefficients,
# which also depends on the regularization parameter 'alpha'.
# We can use the estimated 'alpha' from a previous run or a reasonable guess.
# Here we'll use a guess of alpha=1.
alpha_guess <- 1.0

# Calculate the prior variances based on the formula
prior_variances <- 1 / (1 + alpha_guess * obj1_tps[[4]]$S_diag_truncated)

prior_variances <- 1 / (1 + alpha_est * obj1_tps[[4]]$S_diag_truncated)


# Plot the decay of the prior variances
plot(prior_variances, type = "l", log = "y", col = "blue",
     main = "Decay of Prior Variances of KLE Coefficients",
     xlab = "KLE Basis Index (k)",
     ylab = "Prior Variance")

# Add a point to show where the unpenalized modes end and penalized modes begin.
points(obj1_tps[[4]]$M_P_null_space, prior_variances[obj1_tps[[4]]$M_P_null_space], col = "red", pch = 19)
text(obj1_tps[[4]]$M_P_null_space, prior_variances[obj1_tps[[4]]$M_P_null_space], pos = 4, labels = "End of null space")





n_replicates <- 30

rmse_vec <- numeric(n_replicates)
r2_vec   <- numeric(n_replicates)
mae_vec  <- numeric(n_replicates)

set.seed(2025)

for (rep in 1:n_replicates) {
  
  #=== SIMULATE NEW DATA ===#
  true_field_obs <- MASS::mvrnorm(1, mu = rep(0, N_sp), Sigma = Cov_true_obs)
  y_obs <- true_field_obs + rnorm(N_sp, 0, sigma0_error)
  
  # Recreate data frame and GAM fit
  data_smooth <- data.frame(s1 = sp_points$s1, s2 = sp_points$s2, y_obs = y_obs)
  gam_fit <- gam(y_obs ~ s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth)
  
  Phi_basis_sp   <- predict(gam_fit, newdata = sp_points, type = "lpmatrix")
  Phi_basis_grid <- predict(gam_fit, newdata = grid_total, type = "lpmatrix")
  
  # Same S matrix and decomposition
  sm <- smoothCon(s(s1, s2, k = k_basis, bs = "tp"), data = data_smooth, absorb.cons = FALSE)[[1]]
  S_grid <- sm$S[[1]]
  S_eig  <- eigen(S_grid, symmetric = TRUE)
  S_diag <- S_eig$values
  evectors <- S_eig$vectors
  
  # Reorder
  order_idx <- order(S_diag, decreasing = TRUE)
  S_diag <- S_diag[order_idx]
  evectors <- evectors[, order_idx]
  
  # Null space and M_truncation
  M_P_null_space <- sm$null.space.dim
  evalues <- S_diag
  cumvar <- cumsum(evalues) / sum(evalues)
  M0 <- which(cumvar >= target_var)[1]
  M_truncation <- if (M0 <= M_P_null_space) 0 else min(M0 - M_P_null_space, length(S_diag) - M_P_null_space)
  
  #== Prepare TMB data ==#
  tmb_data <- list(
    Phi_basis_sp = as.matrix(Phi_basis_sp),
    Phi_basis_grid = as.matrix(Phi_basis_grid),
    S_diag = as.vector(S_diag),
    evectors = as.matrix(evectors),
    y = y_obs,
    M_P_null_space = M_P_null_space,
    M_truncation = M_truncation
  )
  
  tmb_pars <- list(
    Z = rep(0, M_truncation),
    logsigma = log(0.5),
    logalpha = log(1.0)
  )
  
  obj <- MakeADFun(data = tmb_data, parameters = tmb_pars, DLL = "kle_model_alpha_M", random = "Z", silent = TRUE)
  opt <- try(nlminb(obj$par, obj$fn, obj$gr), silent = TRUE)
  
  if (inherits(opt, "try-error")) next  # Skip failed replicate
  
  field_grid_est <- try(as.vector(obj$report()$field_grid), silent = TRUE)
  if (inherits(field_grid_est, "try-error")) next
  
  #== Compare with true field on grid ==#
  rmse_vec[rep] <- sqrt(mean((field_grid_est - true_field_grid)^2))
  r2_vec[rep]   <- cor(field_grid_est, true_field_grid)^2
  mae_vec[rep]  <- mean(abs(field_grid_est - true_field_grid))
}

#== Summarize ==#
mean_rmse <- mean(rmse_vec, na.rm = TRUE)
sd_rmse   <- sd(rmse_vec, na.rm = TRUE)

mean_r2   <- mean(r2_vec, na.rm = TRUE)
sd_r2     <- sd(r2_vec, na.rm = TRUE)

mean_mae  <- mean(mae_vec, na.rm = TRUE)
sd_mae    <- sd(mae_vec, na.rm = TRUE)

cat("RMSE: ", mean_rmse, "±", sd_rmse, "\n")
cat("R²:   ", mean_r2, "±", sd_r2, "\n")
cat("MAE:  ", mean_mae, "±", sd_mae, "\n")

# Optional histogram
hist(rmse_vec, main = "RMSE across replicates", xlab = "RMSE")


field_grid_stats <- summary(rep_tmb, "report")
field_grid_idx <- which(rownames(field_grid_stats) == "field_grid")
field_grid_mean <- field_grid_stats[field_grid_idx, "Estimate"]
field_grid_sd   <- field_grid_stats[field_grid_idx, "Std. Error"]


df_grid <- data.frame(
  s1 = grid_total$s1,
  s2 = grid_total$s2,
  mean = field_grid_mean,
  sd = field_grid_sd
)

library(ggplot2)

ggplot(df_grid, aes(s1, s2, fill = sd)) +
  geom_tile() +
  scale_fill_viridis_c() +
  ggtitle("Posterior SD of GRF estimate on grid")












#==================================
# How to evaluate different M values:
#==================================
# To evaluate different M values, you would typically wrap the TMB fitting
# and evaluation part in a loop or a function, varying M_truncation_value.

# Example of a loop to evaluate different M values:
# (Uncomment and run this section to see it in action)

M_values_to_test <- c(M_P_null_space, M_P_null_space + 5, M_P_null_space + 10, 
                      M_P_null_space + 20, M_P_null_space + 30, n_basis_grid)

results_M <- data.frame(M = integer(), MSE = numeric(), Correlation = numeric(), Alpha_Est = numeric(), Sigma_Est = numeric())

for (M_val in M_values_to_test) {
  if (M_val > n_basis_grid) {
    M_val <- n_basis_grid
    warning(paste("Adjusted M_val to n_basis_grid (", n_basis_grid, ") for test M =", M_val))
  }
  if (M_val < M_P_null_space) { # Ensure at least polynomial modes are included
    M_val <- M_P_null_space
  }

  cat(paste("\n--- Running TMB for M =", M_val, "---\n"))

  # Update TMB data for the new M_truncation_value
  tmb_data$M_truncation <- M_val

  # Update TMB parameters for the new Z size
  tmb_pars$Z <- rep(0.0, M_val)

  # Re-make ADFun object for the new data/parameter structure
  obj_M <- MakeADFun(data = tmb_data, parameters = tmb_pars, DLL = "kle_model_alpha_M", random = "Z")
  opt_M <- nlminb(obj_M$par, obj_M$fn, obj_M$gr)
  rep_tmb_M <- sdreport(obj_M)

  Z_est_M <- as.vector(rep_tmb_M$par.random[names(rep_tmb_M$par.random) == "Z"])
  sigma_est_M <- exp(rep_tmb_M$par.fixed["logsigma"])
  alpha_est_M <- exp(rep_tmb_M$par.fixed["logalpha"])
  app_grf_grid_M <- obj_M$report()$field_grid

  mse_M <- mean((true_grf_grid - app_grf_grid_M)^2)
  correlation_M <- cor(true_grf_grid, as.numeric(app_grf_grid_M))

  results_M <- rbind(results_M, data.frame(M = M_val, MSE = mse_M, Correlation = correlation_M, Alpha_Est = alpha_est_M, Sigma_Est = sigma_est_M))

  message(paste("M =", M_val, " | MSE =", round(mse_M, 4), " | Correlation =", round(correlation_M, 4), " | Alpha_Est =", round(alpha_est_M, 4)))

  par(mfrow = c(3, 2))
  image(matrix(true_grf_grid, dim_grid, dim_grid), main = "True GRF", col = hcl.colors(100, "viridis"), asp = 1)
  image(matrix(-app_grf_grid_M, dim_grid, dim_grid), main = paste0("Est GRF (M=", M_val, ", alpha=", round(alpha_est_M,2), ")"), col = hcl.colors(100, "viridis"), asp = 1)
}
#
print(results_M)





#=========================================
# MCMC (HMC of Stan)
#=========================================

startTime <- Sys.time()
mcmc_spde <- tmbstan(obj_spde, chains= 3, open_progress = FALSE,
               control = list(max_treedepth= 12,  adapt_delta = 0.9),
               iter = 3000, warmup= 700, cores=no_cores,
               init = 'last.par.best', seed = 12345)
endTime <- Sys.time()
timeUsed = difftime(endTime, startTime, units='mins')
print(timeUsed)



startTime <- Sys.time()
mcmc_tps <- tmbstan(obj, chains= 3, open_progress = FALSE,
                     control = list(max_treedepth= 12,  adapt_delta = 0.9),
                     iter = 3000, warmup= 700, cores=no_cores,
                     init = 'last.par.best', seed = 12345)
endTime <- Sys.time()
timeUsed = difftime(endTime, startTime, units='mins')
print(timeUsed)


# mcmc_spde <- tmbstan(obj_spde, chains=3, warmup=200, cores=no_cores)
# mcmc_tps  <- tmbstan(obj, chains=3, warmup=200, cores=no_cores)

spde_post_grid <- extract(mcmc_spde)$u  # Posterior samples of KL coefficients
spde_mean_grid <- colMeans(spde_post_grid)
spde_field_grid <- as.vector(A_grid %*% spde_mean_grid)


tps_post_grid <- extract(mcmc_tps)$Z  # Posterior samples of KL coefficients
tps_mean_grid <- colMeans(tps_post_grid)
tps_field_grid_ <- as.vector(Phi_kle_grid %*% tps_mean_grid)



par(mfrow = c(3, 1))
image(matrix(true_field_grid, 30, 30), main = "True GRF")
image(matrix(spde_field_grid, 30, 30), main = "Estimated GRF via SPDE")
image(matrix(tps_field_grid, 30, 30), main = "Estimated GRF via KLE/TPS")




# Mean squared error
mse_spde <- mean((true_field_grid - spde_field_grid)^2); mse_spde
mse_tps <- mean((true_field_grid - tps_field_grid)^2); mse_tps


true_field_grid_plot <- matrix(true_field_grid, nrow = dim_grid, byrow = TRUE)

persp(s1_grid, s2_grid, true_field_grid_plot,
      theta = 60, phi = 30,
      expand = 0.5,
      col = "lightblue",
      ticktype = "detailed",
      main = paste0("3D Surface true GRF"))


spde_field_grid_plot <- matrix(spde_field_grid, nrow = dim_grid, byrow = TRUE)

persp(s1_grid, s2_grid, spde_field_grid_plot,
      theta = 60, phi = 30,
      expand = 0.5,
      col = "lightblue",
      ticktype = "detailed",
      main = paste0("3D Surface TMB"))


tps_field_grid_plot <- matrix(tps_field_grid, nrow = dim_grid, byrow = TRUE)

persp(s1_grid, s2_grid, tps_field_grid_plot,
      theta = 60, phi = 30,
      expand = 0.5,
      col = "lightblue",
      ticktype = "detailed",
      main = paste0("3D Surface Stan"))


# Check Stan's convergence diagnostics (Rhat should be near 1, effective_samples > 400)
print(summary(mcmc_tps)$summary[c("logsigma", "logalpha", "Z[1]", "Z[2]"), c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")])

# Convert samples to a data frame (useful for custom plotting)
posterior_samples <- as.data.frame(mcmc_tps)

# =========================================================================
# Obtain Posterior Distributions and Credible Intervals
# =========================================================================

# 1. For logsigma and logalpha
# You can plot their posterior densities:
library(bayesplot)
mcmc_dens(mcmc_tps, pars = c("logsigma", "logalpha"))
mcmc_trace(mcmc_tps, pars = c("logsigma", "logalpha")) # Check chain mixing

# Transform samples to original scale (sigma and alpha)
posterior_samples$sigma <- exp(posterior_samples$logsigma)
posterior_samples$alpha <- exp(posterior_samples$logalpha)

# Plot transformed posteriors
mcmc_dens(posterior_samples, pars = c("sigma", "alpha"))
mcmc_hist(posterior_samples, pars = c("sigma", "alpha"))

# Get credible intervals for sigma and alpha
sigma_ci <- quantile(posterior_samples$sigma, c(0.025, 0.975))
alpha_ci <- quantile(posterior_samples$alpha, c(0.025, 0.975))

cat("\nPosterior Credible Intervals (from MCMC):\n")
cat(paste0("Sigma 95% CI: [", round(sigma_ci[1], 4), ", ", round(sigma_ci[2], 4), "]\n"))
cat(paste0("Alpha 95% CI: [", round(alpha_ci[1], 4), ", ", round(alpha_ci[2], 4), "]\n"))


# 2. For Z coefficients
# You have many Z coefficients. Plotting all of them isn't practical.
# You might select a few interesting ones or analyze their properties.
# Example: posterior for the first few Z coefficients
mcmc_dens(mcmc_tps, pars = c("Z[1]", "Z[2]", "Z[3]"))
mcmc_trace(mcmc_tps, pars = c("Z[1]", "Z[2]", "Z[3]"))

# If you want to obtain the full posterior distribution of 'field_grid'
# (i.e., the spatial field itself across the grid), you'll need to re-evaluate
# the 'field_grid' calculation for each MCMC sample.

# This involves:
# a. Extracting Z, logsigma, logalpha for each posterior sample.
# b. In R, re-calculating lambda_KL and then field_grid using the basis functions.
# This can be computationally intensive if you have many samples and grid points.

# Example for re-evaluating field_grid for a subset of posterior samples:
# (This is a simplified example; for full field, you need the basis functions)
# For this, you would use your Phi_basis_grid and evectors from tmb_data
# and the lambda_KL_full calculation from your TMB code.

# Extract required data from original tmb_data
Phi_basis_grid_mat <- tmb_data$Phi_basis_grid
evectors_mat <- tmb_data$evectors
S_diag_vec <- tmb_data$S_diag
M_truncation_val <- tmb_data$M_truncation

# Truncate eigenfunctions_grid_full as done in TMB
eigenfunctions_grid_full <- Phi_basis_grid_mat %*% evectors_mat
eigenfunctions_grid_truncated <- eigenfunctions_grid_full[, 1:M_truncation_val]


# Number of posterior samples
n_samples <- nrow(posterior_samples)
n_grid_points <- nrow(Phi_basis_grid_mat)

# Store posterior samples of field_grid
field_grid_posterior <- matrix(NA, nrow = n_samples, ncol = n_grid_points)

cat("\nCalculating posterior of field_grid (this may take a moment)...\n")
for (i in 1:n_samples) {
  current_alpha <- posterior_samples$alpha[i]
  current_Z <- as.numeric(posterior_samples[i, grep("^Z\\[", colnames(posterior_samples))]) # Get Z vector
  
  # Calculate lambda_KL_full for this alpha
  lambda_KL_full_current <- 1 / (1 + current_alpha * S_diag_vec)
  # Ensure Z is multiplied by the correct truncated eigenfunctions
  # Z coefficients correspond to the truncated eigenfunctions
  field_grid_posterior[i, ] <- eigenfunctions_grid_truncated %*% current_Z
}

# Now, field_grid_posterior contains samples of the spatial field on the grid.
# You can get the posterior mean field:
mean_field_grid_posterior <- colMeans(field_grid_posterior)

# And visualize a credible interval for the field at each grid point (e.g., 95% CI)
lower_ci_field_grid <- apply(field_grid_posterior, 2, quantile, probs = 0.025)
upper_ci_field_grid <- apply(field_grid_posterior, 2, quantile, probs = 0.975)

image(matrix(lower_ci_field_grid, dim_grid, dim_grid),
      main = "Posterior Mean of GRF (MCMC)",
      col = hcl.colors(100, "viridis"), asp = 1)

image(matrix(mean_field_grid_posterior, dim_grid, dim_grid),
      main = "Posterior Mean of GRF (MCMC)",
      col = hcl.colors(100, "viridis"), asp = 1)

image(matrix(upper_ci_field_grid, dim_grid, dim_grid),
      main = "Posterior Mean of GRF (MCMC)",
      col = hcl.colors(100, "viridis"), asp = 1)


