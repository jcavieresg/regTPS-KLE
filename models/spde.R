setwd("C:/Users/Usuario/Desktop/KLE")
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
compile("spde.cpp")
dyn.load(dynlib("spde"))



set.seed(1234)
run_tmb <- function(N_sp, dim_grid, sp_points = NULL, mesh = NULL){
set.seed(1234)
  
  sigma0_error <- 0.1

  # Simulate GRF using Matern kernel
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
  
  s1_grid <- seq(0, 1, length.out = dim_grid)
  s2_grid <- seq(0, 1, length.out = dim_grid)
  grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
  
  Cov_true_grid <- matern_cov(grid_total)
  true_fiel_grid <- MASS::mvrnorm(1, mu = rep(0, nrow(grid_total)), Sigma = Cov_true_grid)
  
  # if (is.null(mesh)) {
  #   mesh <- inla.mesh.2d(loc = sp_points, max.edge = c(0.05, 0.39), cutoff = 0.115)
  # }
  
  sp_matrix <- as.matrix(sp_points)
  bound1 <- inla.nonconvex.hull(sp_matrix)
  mesh <- inla.mesh.create(loc = sp_matrix, boundary = bound1, refine = FALSE, plot.delay = NULL)
  
  spde <- inla.spde2.matern(mesh, alpha = 2)
  spde_mat <- spde$param.inla[c("M0", "M1", "M2")]
  
  A_obs  <- inla.spde.make.A(mesh = mesh, loc = as.matrix(sp_points))
  A_grid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(grid_total))
  
  #========================================
  #                TMB data
  #========================================
  tmb_data <- list(y = y_obs, 
                   A_obs = A_obs, 
                   A_grid = A_grid, 
                   spde_mat = spde_mat, 
                   lambda_rho = -log(0.05)*0.05,
                   lambda_sigma_u = -log(0.05)/5)
  #========================================
  #                TMB par
  #========================================
  tmb_par <- list(sigma = 0.1,
                  rho   = 0.1,
                  sigma_u = 0.1,
                  u = rnorm(mesh$n, 0, 1))
  
  obj_spde <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "spde_model2", random = "u")
  lwr <- c(0, 0, 0)
  upr <- c(Inf, Inf, Inf)
  opt_spde = nlminb(obj_spde$par, obj_spde$fn, obj_spde$gr, lower = lwr, upper = upr)
  rep_spde <- sdreport(obj_spde)
  
  res_list <- list(obj = obj_spde, opt = opt_spde, rep = rep_spde, tmb_data = tmb_data, tmb_par = tmb_par, mesh = mesh, 
                   spde = spde, true_fiel_grid = true_fiel_grid)
  return(res_list)
}


#=======================
# Base spatial data 
# base_N_sp <- 50
# scenarios <- list()
# 
# # Define custom max.edge values to control number of triangles
# max_edge_list <- list(
#   c(0.20, 0.4),  # Scenario 1 (~100 triangles)
#   c(0.18, 0.35), # Scenario 2 (~130)
#   c(0.16,  0.3),  # Scenario 3 (~160)
#   c(0.14, 0.3)  # Scenario 4 (~200)
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
#   scenarios[[paste0("scenario", i)]] <- list(N_sp = N_sp, sp_points = sp_points, mesh = mesh,
#                                              n_triangles = nrow(mesh$graph$tv), n_vertices = mesh$n)}



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
    cat(paste("Scenario", i, "created with", mesh$n, "nodes.\n"))
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



# coords = as.matrix(scenarios$scenario4$sp_points)
# bound1 <- inla.nonconvex.hull(coords)
# mesh = inla.mesh.create(coords, plot.delay=NULL, refine=TRUE, boundary = bound1)
# plot(mesh);points(coords, col = "red")


scenarios$scenario1$mesh$n
scenarios$scenario2$mesh$n
scenarios$scenario3$mesh$n
scenarios$scenario4$mesh$n



obj1_spde <- run_tmb(N_sp = scenarios$scenario1$N_sp, dim_grid = 30, sp_points = scenarios$scenario1$sp_points, mesh = scenarios$scenario1$mesh)
obj2_spde <- run_tmb(N_sp = scenarios$scenario2$N_sp, dim_grid = 30, sp_points = scenarios$scenario2$sp_points, mesh = scenarios$scenario2$mesh)
obj3_spde <- run_tmb(N_sp = scenarios$scenario3$N_sp, dim_grid = 30, sp_points = scenarios$scenario3$sp_points, mesh = scenarios$scenario3$mesh)
obj4_spde <- run_tmb(N_sp = scenarios$scenario4$N_sp, dim_grid = 30, sp_points = scenarios$scenario4$sp_points, mesh = scenarios$scenario4$mesh)


u_est <- as.vector(rep_spde$par.random[names(rep_spde$par.random) == "u"])
spde_field_grid <- as.vector(obj_spde$report()$field_grid)

tau_est <- exp(rep_spde$par.fixed["logtau"])
kappa_est <- exp(rep_spde$par.fixed["logkappa"])
sigma_est <- exp(rep_spde$par.fixed["logsigma"])

# SPDE
rmse_spde <- sqrt(mean((spde_field_grid - true_field_grid)^2))
r2_spde   <- cor(spde_field_grid, true_field_grid)^2
mae_spde  <- mean(abs(spde_field_grid - true_field_grid))






M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n50"
M[[2]] = list()
M[[2]]$model = "spatial_n100"
M[[3]] = list()
M[[3]]$model = "spatial_n150"
M[[4]] = list()
M[[4]]$model = "spatial_n200"


M[[1]]$formula = obj1_spde[[1]]
M[[2]]$formula = obj2_spde[[1]]
M[[3]]$formula = obj3_spde[[1]]
M[[4]]$formula = obj4_spde[[1]]




#===========================================
#           Run the models
#===========================================

lwr <- c(0, 0, 0)
upr <- c(Inf, Inf, Inf)

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
  saveRDS(fit, file=paste0('stan_spde_', i,'.RDS'))
}




u_est <- as.vector(obj1[[3]]$par.random[names(obj1[[3]]$par.random) == "u"])
sigma_est <- obj1[[3]]$par.fixed["sigma"]
rho_est <- obj1[[3]]$par.fixed["rho"]
sigmau_est <- obj1[[3]]$par.fixed["sigma_u"]


# Get the approximated GRF on the full grid from TMB's REPORT
spde_field_grid <- as.vector(obj1[[1]]$report()$field_grid)

message(paste("Estimated sigma:", round(sigma_est, 3)))
message(paste("Estimated rho:", round(rho_est, 3)))
message(paste("Estimated sigma_u:", round(sigmau_est, 3)))


mcmc_spde1 <- readRDS('stan_spde_1.RDS')
mcmc_spde2 <- readRDS('stan_spde_2.RDS')
mcmc_spde3 <- readRDS('stan_spde_3.RDS')
mcmc_spde4 <- readRDS('stan_spde_4.RDS')


spde_post_grid1 <- extract(mcmc_spde1)$u  # Posterior samples of KL coefficients
spde_mean_grid1 <- colMeans(spde_post_grid1)
spde_field_grid1 <- as.vector(obj1_spde[[4]]$A_grid %*% spde_mean_grid1)

spde_post_grid2 <- extract(mcmc_spde2)$u  # Posterior samples of KL coefficients
spde_mean_grid2 <- colMeans(spde_post_grid2)
spde_field_grid2 <- as.vector(obj2_spde[[4]]$A_grid %*% spde_mean_grid2)

spde_post_grid3 <- extract(mcmc_spde3)$u  # Posterior samples of KL coefficients
spde_mean_grid3 <- colMeans(spde_post_grid3)
spde_field_grid3 <- as.vector(obj3_spde[[4]]$A_grid %*% spde_mean_grid3)

spde_post_grid4 <- extract(mcmc_spde4)$u  # Posterior samples of KL coefficients
spde_mean_grid4 <- colMeans(spde_post_grid4)
spde_field_grid4 <- as.vector(obj4_spde[[4]]$A_grid %*% spde_mean_grid4)





# Plotting (assuming grid_total coordinates are still dim_grid x dim_grid)
par(mfrow = c(4, 2))
image(matrix(obj1[[8]], 30, 30), main = "True GRF 1",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid1, 30, 30), main = paste0("Estimated GRF (spde=", obj1[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(obj2[[8]], 30, 30), main = "True GRF 2",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid1, 30, 30), main = paste0("Estimated GRF (spde=", obj2[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


image(matrix(obj3[[8]], 30, 30), main = "True GRF 3",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid1, 30, 30), main = paste0("Estimated GRF (spde=", obj3[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)



image(matrix(obj4[[8]], 30, 30), main = "True GRF 4",
      col = hcl.colors(100, "viridis"), asp = 1)
image(matrix(spde_field_grid1, 30, 30), main = paste0("Estimated GRF (spde=", obj4[[6]]$n, ")"),
      col = hcl.colors(100, "viridis"), asp = 1)


# SPDE
rmse_spde <- sqrt(mean((spde_field_grid1 - obj1[[8]])^2))
r2_spde   <- cor(as.vector(spde_field_grid1), obj1[[8]])^2
mae_spde  <- mean(abs(as.vector(spde_field_grid1) - obj1[[8]]))

rmse_spde
r2_spde
mae_spde



rmse_spde <- sqrt(mean((as.vector(obj4[[1]]$report()$field_grid) - obj4[[8]])^2))
r2_spde   <- cor(as.vector(obj4[[1]]$report()$field_grid), obj4[[8]])^2
mae_spde  <- mean(abs(as.vector(obj4[[1]]$report()$field_grid) - obj4[[8]]))

rmse_tps
r2_tps
mae_tps