// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;

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


// Improved dcauchy with better numerical stability
// template<class Type>
// Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
//   Type logres = Type(0.0);
//   logres -= log(M_PI);
//   logres -= log(shape);
//   Type u = (x - mean) / shape;
//   // Use log1p for even better stability when u^2 is small
//   logres -= log1p(u * u);
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


// Log-normal prior for positive parameters
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = Type(0.0);
  logres -= log(x);
  logres -= log(sdlog);
  logres -= Type(0.5) * log(2.0 * M_PI);
  logres -= Type(0.5) * pow((log(x) - meanlog) / sdlog, 2);
  if(give_log) return logres; else return exp(logres);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(y_obs);
  DATA_MATRIX(Phi_kle_sp); // The pre-computed KLE basis at observation points
  // DATA_MATRIX(Phi_kle_grid); // The pre-computed KLE basis at grid points
  DATA_VECTOR(S_diag_truncated); // Eigenvalues of S (truncated)
  DATA_INTEGER(M_P_null_space); // Number of polynomial modes
  DATA_SCALAR(sigma_prior_s0);     // s0 for P(sigma > s0) = alpha_s (PC prior)
  DATA_SCALAR(sigma_prior_alpha);  // alpha_s for sigma PC prior (e.g. 0.05)
  DATA_SCALAR(logalpha_prior_mean); // prior mean for logalpha (default 0)
  DATA_SCALAR(logalpha_prior_sd);   // prior sd for logalpha (default 1)
  
  
  // PARAMETERS
  //PARAMETER_VECTOR(z); // KLE coefficients (length M_truncation)
  PARAMETER_VECTOR(z_raw);
  PARAMETER(logsigma); // Log of observation noise SD
  PARAMETER(logalpha); // Log of regularization parameter

  // MODEL SETUP
  Type sigma = exp(logsigma);
  Type alpha = exp(logalpha);
  
  //==================================================
  // PRIORS
  //==================================================
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  // Prior on sigma
  // //nlp -= dlognorm(sigma,   Type(0.0),   Type(1.0), true);
  // //nlp += logsigma; // Jacobian
  // nlp -= dcauchy_stable(sigma_e, Type(0.0), cauchy_scale_e, true);
  // nlp -= logsigma_e; // Jacobian for log-parameterization
  // 
  // // Prior on alpha 
  // nlp -= dlognorm(alpha, Type(0.0), Type(1.0), true);
  // nlp += logalpha; // Jacobian

  
  // Prior for sigma: PC prior (exponential) on sigma
  // lambda = -log(alpha_s) / s0  where P(sigma > s0) = alpha_s
  Type lambda_sigma = -log(sigma_prior_alpha) / sigma_prior_s0;
  // prior on sigma: p(sigma) = lambda * exp(-lambda*sigma) for sigma>0
  // we parameterize via logsigma -> need jacobian term (log sigma)
  nlp -= dexp(sigma, lambda_sigma, true); // subtract log p(sigma)
  nlp -= logsigma; // subtract log|d sigma / d logsigma|  (we subtract because nlp is -log prior)
  
  // Prior for alpha: put a Normal prior on logalpha (i.e., log-normal prior on alpha)
  nlp -= dnorm(logalpha, logalpha_prior_mean, logalpha_prior_sd, true);
  // NO jacobian when priors are on logalpha itself
  
  
  int M_trunc = z_raw.size();
  
  // ----- Prior for z_raw: standard normals for all modes (non-centered) -----
  for(int k = 0; k < M_trunc; ++k){
    nlp -= dnorm(z_raw(k), Type(0.0), Type(1.0), true);
  }
  
  // -------------------------
  // Build actual z coefficients (apply prior scaling for penalized modes)
  // -------------------------
  vector<Type> z(M_trunc);
  for(int k = 0; k < M_trunc; ++k){
    if(k < M_P_null_space){
      // unpenalized / null-space modes: prior variance = 1 -> prior_sd = 1
      z(k) = z_raw(k);
    } else {
      Type S_k = S_diag_truncated(k);
      // prior variance: lambda_k = 1 / (1 + alpha * S_k)
      Type prior_var = Type(1.0) / (Type(1.0) + alpha * S_k);
      Type prior_sd = sqrt(prior_var);
      z(k) = prior_sd * z_raw(k);
    }
  }
  
  
  // z has M_truncation elements
  vector<Type> field_sp = Phi_kle_sp * z;

  //====================================================
  // Likelihood
  Type nll = Type(0.0);
  
  int n = y_obs.size();
  vector<Type> log_lik(n);
  
  for(int i = 0; i < y_obs.size(); i++){
    log_lik(i) = dnorm(y_obs(i), field_sp(i), sigma, true);
    nll -= log_lik(i);
  }
  
  // Total negative log-like
  Type jnll = nll + nlp;
  
  
  // Simule data from the mu 
  vector<Type> y_sim(n);
  SIMULATE {
    for(int i = 0; i < n; i++){
      Type y_tmp = rnorm(field_sp(i), sigma);  // simulate as usual
      y_sim(i) = pow(y_tmp, 2);                // then square it
      // y_sim(i) = rnorm(field_sp(i), sigma);
    }
    REPORT(y_sim); // Inside SIMULATE block
  }

  //===========================
  // REPORT
  //===========================
  REPORT(field_sp);
  REPORT(sigma);
  REPORT(alpha);
  REPORT(z_raw);
  REPORT(z);
  REPORT(log_lik);
  REPORT(y_sim);
  REPORT(nll);
  REPORT(nlp);
  
  //====================================================
  // ADREPORT for uncertainty quantification
  //====================================================
  ADREPORT(sigma);
  ADREPORT(alpha);
  ADREPORT(z_raw);
  ADREPORT(z);
  ADREPORT(field_sp); // This will give you uncertainty on predictions
  
  
  return jnll;
}