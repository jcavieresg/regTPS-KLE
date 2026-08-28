// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;


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

  
//====================================================
// DATA
//====================================================
  DATA_VECTOR(y_obs);
  DATA_MATRIX(Phi_kle_sp);            // The pre-computed KLE basis at observation points
  DATA_VECTOR(S_diag_truncated);      // Eigenvalues of S (truncated)
  DATA_INTEGER(M_P_null_space);       // Number of polynomial modes

  DATA_SCALAR(lambda_sigma_e);     
  DATA_SCALAR(mean_logalpha);
  DATA_SCALAR(sigma_logalpha);
  
  
//====================================================
// PARAMETERS
//====================================================
  PARAMETER_VECTOR(z_tilde);
  PARAMETER(logsigma_e); // Log of sigma_e
  PARAMETER(logalpha);   // Log of regularization parameter

  // MODEL SETUP
  Type sigma_e = exp(logsigma_e);
  Type alpha   = exp(logalpha);
  
  //==================================================
  // PRIORS
  //==================================================
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  // Prior for sigma: PC prior (exponential)
  nlp -= dexp(sigma_e, lambda_sigma_e, true); // subtract log p(sigma)
  nlp -= logsigma_e; 
  
  // Prior for alpha: put a Normal prior on logalpha (i.e., log-normal prior on alpha)
  // nlp -= dnorm(logalpha, Type(1.0), lambda_alpha, true);
  nlp -= dnorm(logalpha, mean_logalpha, sigma_logalpha, true);
  
  
  // Standard normal prior on whitened coefficients
  int M = z_tilde.size();
  for(int k = 0; k < M; k++){
    nlp -= dnorm(z_tilde(k), Type(0.0), Type(1.0), true);
  }
  
  
  //====================================================
  // KLE SCALING (non-centered parameterization)
  //====================================================
  vector<Type> scale(M);
  
  for(int k = 0; k < M; k++){
    if(k < M_P_null_space){
      // Null space: unpenalized
      scale(k) = Type(1.0);
    } else {
      // Mass matrix formulation: λ_k = 1/(1 + α v_k)
      // scale_k = sqrt(λ_k) = 1/sqrt(1 + α v_k)
      scale(k) = Type(1.0) / sqrt(Type(1.0) + alpha * S_diag_truncated(k) + Type(1e-10));
    }
  } 
  
  vector<Type> z = scale * z_tilde;
  
  
  //====================================================
  // FIELD CONSTRUCTION
  //====================================================
  vector<Type> field_sp = Phi_kle_sp * z;


  //====================================================
  // Likelihood
  //====================================================
  Type nll = Type(0.0);
  
  int n = y_obs.size();
  vector<Type> log_lik(n);
  
  for(int i = 0; i < n; i++){
    log_lik(i) = dnorm(y_obs(i), field_sp(i), sigma_e, true);
    nll -= log_lik(i);
  }
  

// Total negative log-like
  Type jnll = nll + nlp;
  

//======================================================  
// Simule data from the mu 
//======================================================
  vector<Type> y_sim(n);
  SIMULATE {
    for(int i = 0; i < n; i++){
      Type y_tmp = rnorm(field_sp(i), sigma_e);  // simulate as usual
      y_sim(i) = pow(y_tmp, 2);                // then square it
    }
    REPORT(y_sim); // Inside SIMULATE block
  }

//===========================
// REPORT
//===========================
REPORT(field_sp);
REPORT(sigma_e);
REPORT(alpha);
REPORT(z_tilde);
REPORT(z);
REPORT(log_lik);
REPORT(y_sim);
REPORT(nll);
REPORT(nlp);
  
//====================================================
// ADREPORT for uncertainty quantification
//====================================================
ADREPORT(sigma_e);
ADREPORT(alpha);
ADREPORT(z_tilde);
ADREPORT(z);
ADREPORT(field_sp); // This will give you uncertainty on predictions
return jnll;
}
