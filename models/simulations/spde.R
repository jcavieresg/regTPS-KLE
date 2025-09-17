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
compile("spde.cpp")
dyn.load(dynlib("spde"))



run_tmb_spde <- function(N_sp, dim_grid, sp_points = NULL, mesh = NULL, 
                    nu = 1.5, rho = 0.2, sigma_u = 1.0, sigma0_error = 0.1){
set.seed(1234)
  
  # Convert sp_points to matrix
  sp_matrix <- as.matrix(sp_points)
  if (is.null(mesh)) {
    bound1 <- inla.nonconvex.hull(sp_matrix)
    mesh <- inla.mesh.create(loc = sp_matrix, boundary = bound1, refine = FALSE, plot.delay = NULL)
  }
  
  # Matern covariance function
  matern_cov <- function(coords, nu = 1.5, rho = 0.2, sigma2 = 1.0) {
    D <- as.matrix(dist(coords))
    D[D == 0] <- 1e-10
    scaling <- (sqrt(2 * nu) * D) / rho
    matern_part <- (2^(1 - nu)) / gamma(nu) * scaling^nu * besselK(scaling, nu)
    diag(matern_part) <- 1
    return(sigma2 * matern_part)
  }
  
  # Simulate latent field on mesh nodes using matern_cov
  Cov_true <- matern_cov(mesh$loc, nu = nu, rho = rho, sigma2 = sigma_u^2)
  u_true <- as.numeric(MASS::mvrnorm(1, mu = rep(0, mesh$n), Sigma = Cov_true))
  
  # Create grid and projection matrices
  s1_grid <- seq(0, 1, length.out = dim_grid)
  s2_grid <- seq(0, 1, length.out = dim_grid)
  grid_total <- expand.grid(s1 = s1_grid, s2 = s2_grid)
  
  A_obs  <- inla.spde.make.A(mesh = mesh, loc = sp_matrix)
  A_grid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(grid_total))
  
  # Observations and grid truth
  y_obs      <- as.numeric(A_obs %*% u_true + rnorm(N_sp, 0, sigma0_error))
  u_grid <- as.numeric(A_grid %*% u_true)
  
  spde <- inla.spde2.matern(mesh, alpha = 2)
  spde_mat <- spde$param.inla[c("M0", "M1", "M2")]
  

  
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
  
  obj_spde <- MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "spde", random = "u")
  lwr <- c(0, 0, 0)
  upr <- c(Inf, Inf, Inf)
  opt_spde = nlminb(obj_spde$par, obj_spde$fn, obj_spde$gr, lower = lwr, upper = upr)
  rep_spde <- sdreport(obj_spde)
  
  res_list <- list(obj = obj_spde, opt = opt_spde, rep = rep_spde, tmb_data = tmb_data, tmb_par = tmb_par, mesh = mesh, 
                   spde = spde, true_fiel_grid = u_grid, true_field_obs = u_true)
  return(res_list)
}


#=====================================================================================
#         Function to create a list of scenarios using inla.mesh.create
#=====================================================================================

mesh_scenarios <- function(base_N_sp = 50, n_scenarios = 4) {
  set.seed(1234)  
  scenarios <- list()
  for (i in 1:n_scenarios) {
    N_sp <- base_N_sp * i
    # Always consume RNG in sequence → identical results if re-run
    sp_points <- data.frame(s1 = runif(N_sp), s2 = runif(N_sp))
    sp_matrix <- as.matrix(sp_points)
    bound1 <- inla.nonconvex.hull(sp_matrix)
    mesh <- inla.mesh.create(loc = sp_matrix, boundary = bound1,
                             refine = FALSE, plot.delay = NULL)
    scenario_name <- paste0("scenario", i)
    scenarios[[scenario_name]] <- list(N_sp = N_sp,
                                       sp_points = sp_points,
                                       mesh = mesh,
                                       n_triangles = nrow(mesh$graph$tv),
                                       n_nodes = mesh$n)
    
    cat(paste("Scenario", i, "created with", mesh$n, "vertices.\n"))}
  return(scenarios)}


#=====================================================================================
#                    Run the function to generate the scenarios
#=====================================================================================

# Define the base number of points and the number of scenarios
base_N_sp <- 50
n_scenarios <- 4

# Run the function to generate the scenarios
scenarios <- mesh_scenarios(base_N_sp = base_N_sp, n_scenarios = n_scenarios)


obj1_spde <- run_tmb_spde(N_sp = scenarios$scenario1$N_sp, dim_grid = 30, sp_points = scenarios$scenario1$sp_points, mesh = scenarios$scenario1$mesh)
obj2_spde <- run_tmb_spde(N_sp = scenarios$scenario2$N_sp, dim_grid = 30, sp_points = scenarios$scenario2$sp_points, mesh = scenarios$scenario2$mesh)
obj3_spde <- run_tmb_spde(N_sp = scenarios$scenario3$N_sp, dim_grid = 30, sp_points = scenarios$scenario3$sp_points, mesh = scenarios$scenario3$mesh)
obj4_spde <- run_tmb_spde(N_sp = scenarios$scenario4$N_sp, dim_grid = 30, sp_points = scenarios$scenario4$sp_points, mesh = scenarios$scenario4$mesh)


fits_TMB_spde <- list(obj1_spde, obj2_spde, obj3_spde, obj4_spde)

saveRDS(fits_TMB_spde, file='outputs/fits_TMB_spde.RDS')

obj_spde <- obj1_spde[[1]]
rep_spde <- obj1_spde[[3]]

u_est <- as.vector(rep_spde$par.random[names(rep_spde$par.random) == "u"])
spde_field_grid <- as.vector(obj_spde$report()$field_grid)

tau_est <- exp(rep_spde$par.fixed["logtau"])
kappa_est <- exp(rep_spde$par.fixed["logkappa"])
sigma_est <- exp(rep_spde$par.fixed["logsigma"])

# SPDE
spde_field_grid <- obj1_spde[[1]]$report()$field_grid
true_field_grid <- obj1_spde[[8]]

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