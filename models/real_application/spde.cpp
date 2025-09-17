// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;
using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()


// Can also choose which likelihood to use.
// Lognormal density
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// Inverse gamma
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}
// dcauchy for hyperparameters
// template<class Type>
// Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
//   Type logres = 0.0;
//   logres-= log(M_PI);
//   logres-= log(shape);
//   // Note, this is unstable and should switch to log1p formulation
//   logres-= log(1 + pow( (x-mean)/shape ,2));
//   if(give_log) return logres; else return exp(logres);
// }

template<class Type>
Type log1p_custom(Type x) {
  // Accurate for small x
  if (CppAD::abs(x) < 1e-8) {
    return x - x*x/2; // Taylor expansion
  } else {
    return log(Type(1.0) + x);
  }
}


template<class Type>
Type dcauchy_stable(Type x, Type mean, Type scale, int give_log=0){
  Type z = (x - mean) / scale;
  Type logres = -log(M_PI) - log(scale) - log1p_custom(z * z);
  if(give_log) return logres; else return exp(logres);
}


template<class Type>
Type ldhalfnorm(Type x, Type var){
  return 0.5*log(2)-0.5*log(var*M_PI)+pow(x,2)/(2*var);
}


//=====================================================
// Initialize the TMB model
//=====================================================
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  //=========================
  //      DATA SECTION
  //=========================
  DATA_VECTOR(y);                       // Observed data
  DATA_SPARSE_MATRIX(A_obs);            // Projection matrix for observations
  // DATA_SPARSE_MATRIX(A_grid);           // Projection matrix for grid prediction
  DATA_STRUCT(spde_mat, spde_t);        // Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  
  DATA_SCALAR(rho0);        // user-specified, e.g. 10 (units of distance)
  DATA_SCALAR(alpha_rho);   // e.g. 0.05 for P(rho < rho0)=alpha_rho
  DATA_SCALAR(s0_u);        // e.g. 1.0
  DATA_SCALAR(alpha_s_u);   // e.g. 0.05
  
  DATA_SCALAR(cauchy_scale_e); // e.g. 5.0 (unused if using PC exponential obs prior)
  
  // DATA_SCALAR(lambda_rho);
  // DATA_SCALAR(lambda_sigma_u);
  
  //=========================
  //   PARAMETER SECTION
  //=========================
  // PARAMETER(sigma_e);
  // PARAMETER(rho);
  // PARAMETER(sigma_u);
  
  PARAMETER(logsigma_e);
  PARAMETER(logrho);
  PARAMETER(logsigma_u);

  Type sigma_e = exp(logsigma_e);
  Type rho = exp(logrho);
  Type sigma_u = exp(logsigma_u);
  
  // Parameter vector for spatial field
  //PARAMETER_VECTOR(u);
  PARAMETER_VECTOR(u_raw);

// For PC priors
  Type lambda_rho_inv = -rho0 * log(alpha_rho);    // for Exp on rho^{-1}
  Type lambda_sigma_u = -log(alpha_s_u) / s0_u;    // for Exp on sigma_u
  
  // Type rho_inv = Type(1.0)/rho;
  // Type lambda_rho_inv = Type(1.0)/lambda_rho;
  

  // ===================================
  //               Priors
  // ===================================
  Type nlp = 0.0;
  
  // // Prior on sigma (Half-normal)
  // nlp -= dnorm(sigma_e,   Type(0.0), Type(1.0), true);
  // 
  // // Prior for rho and Jacobian adjustment
  // nlp -= dweibull(rho_inv,  lambda_rho_inv, Type(1.0), true);   //
  // nlp += 2*log(rho);
  // 
  // // Prior for sigma_u
  // nlp -= dexp(sigma_u,      lambda_sigma_u,     true);


  // range prior via exponential on rho^{-1} (implemented as logpdf of X=1/rho)
  // 1) Prior on observation sd sigma_e: use PC prior (exponential) on sigma_e
  //    We parameterize via log_sigma_e. Need to include Jacobian term: log|d sigma_e / d log_sigma_e| = log_sigma_e
  //nlp -= dexp(sigma_e, lambda_sigma_e, true);   // subtract log p(sigma_e)
  //nlp -= log_sigma_e;                           // subtract log jacobian (-> add negative log jacobian)

  // If you prefer half-Cauchy for sigma_e instead, comment the above two lines and uncomment:
  nlp -= dcauchy_stable(sigma_e, Type(0.0), cauchy_scale_e, true);
  nlp -= logsigma_e; // Jacobian for log-parameterization
  
  // 2) Prior for rho expressed as exponential on rho^{-1}
  //    X = rho^{-1} = exp(-log_rho)
  Type rho_inv = Type(1.0) / rho; // = exp(-log_rho)
  // log prior for X:
  nlp -= dweibull(rho_inv, lambda_rho_inv, Type(1.0), true); // log p(X)
  // Jacobian term for parameterising by log_rho:
  // X = exp(-log_rho) => dX/d(log_rho) = -exp(-log_rho) => |dX/dt| = exp(-log_rho)
  // log |dX/dt| = -log_rho. For negative log prior we add - log |dX/dt| = + log_rho
  nlp += logrho;
  
  // 3) Prior for sigma_u (marginal sd of spatial field): PC exponential on sigma_u
  nlp -= dexp(sigma_u, lambda_sigma_u, true); // subtract log p(sigma_u)
  nlp -= logsigma_u;                         // subtract log jacobian (log|d sigma_u / d log_sigma_u|)
  
  // Optionally: weak prior on log_rho/log_sigmas to stabilise (not necessary if above are used)
  // e.g. nlp += 0.5 * pow(log_sigma_u / 10.0, 2.0);  

  // Derived spatial quantities
  Type kappa = sqrt(8)/rho;
  Type tau   = 1/(kappa*sigma_u);
  
  SparseMatrix<Type> Q = Q_spde(spde_mat, kappa);
  
  

  //=============================================================================================================
  // Objective function is sum of negative log likelihood components
  int n = y.size();	                   // number of observations 
  Type nll_u = 0.0;		                 // likelihood for the spatial effect
  //nll_u += SCALE(GMRF(Q), 1/ tau)(u);  // returns negative already
  
  nll_u += GMRF(Q)(u_raw); // u_raw has prior independent of tau
  // Build the actual spatial field used in the linear predictor
  vector<Type> u = u_raw / tau; // non-centered transform
  
  // Linear predictor
  vector<Type> mu(n);
  mu = A_obs * u;
  
  //=======================================================
  // Likelihood
  vector<Type> log_lik(n);
  for( int i = 0; i< n; i++){
    log_lik(i) = dnorm(y(i), mu(i), sigma_e, true);
  }
  Type nll = -log_lik.sum(); // total NLL
  

// Calculate joint negative log likelihood
  Type jnll = nll + nll_u + nlp;
  

  //============================================
  // Simulated data from the mu
  //============================================
  vector<Type> y_sim(n);
  SIMULATE {
    for(int i = 0; i < n; i++){
      Type y_tmp = rnorm(mu(i), sigma_e);  // simulate as usual
      y_sim(i) = pow(y_tmp, 2);                // then square it
      // y_sim(i) = rnorm(mu(i), sigma_e);
    } 
    REPORT(y_sim);
  }

  
  //=============================================
  // Derived quantities
  //=============================================
  // Spatial field in A_obs
  vector<Type> field_sp = A_obs * u;

  // REPORT 
  REPORT(sigma_e);
  REPORT(tau);
  REPORT(kappa);
  REPORT(sigma_u);
  REPORT(rho);
  REPORT(field_sp);
  REPORT(Q);
  REPORT(log_lik);
  REPORT(u);
  //REPORT(u_raw);

  
  // ADREPORT
  ADREPORT(field_sp);
  ADREPORT(sigma_e);
  ADREPORT(tau);
  ADREPORT(kappa);
  ADREPORT(u);
  //ADREPORT(u_raw);

  return jnll;
}