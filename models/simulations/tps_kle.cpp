// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;

template<class Type>
Type log1p_custom(Type x) {
  if (CppAD::abs(x) < 1e-8) {
    return x - x*x/2;
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
  DATA_MATRIX(Phi_kle_sp);
  DATA_MATRIX(Phi_kle_grid);
  DATA_VECTOR(S_diag_truncated);
  DATA_INTEGER(M_P_null_space);
  DATA_SCALAR(lambda_sigma);
  
  //====================================================
  // PARAMETERS
  //====================================================
  PARAMETER_VECTOR(z_tilde);
  PARAMETER(logsigma);
  PARAMETER(logalpha);
  
  //====================================================
  // TRANSFORMED PARAMETERS
  //====================================================
  Type sigma = exp(logsigma);
  Type alpha = exp(logalpha);
  
  //====================================================
  // PRIORS
  //====================================================
  Type nlp = Type(0.0);
  
  // Prior for sigma: PC prior (exponential)
  nlp -= dexp(sigma, lambda_sigma, true);
  nlp -= logsigma;
  
  // Prior for alpha: Log-normal
  nlp -= dnorm(logalpha, Type(0.0), Type(3.0), true);

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
      // CRITICAL FIX: Null space should be UNPENALIZED
      scale(k) = Type(1.0);  // No smoothing for polynomial trend
    } else {
      // Penalized components: λ_k = 1/(1 + α v_k)
      scale(k) = Type(1.0) / sqrt(Type(1.0) + alpha * S_diag_truncated(k) + Type(1e-10));
    }
  }
  
  vector<Type> z = scale * z_tilde;
  
  //====================================================
  // FIELD CONSTRUCTION
  //====================================================
  vector<Type> field_sp = Phi_kle_sp * z;
  vector<Type> field_grid = Phi_kle_grid * z;
  
  //====================================================
  // LIKELIHOOD
  //====================================================
  int n_obs = y_obs.size();
  Type nll = Type(0.0);
  
  for(int i = 0; i < n_obs; i++){
    nll -= dnorm(y_obs(i), field_sp(i), sigma, true);
  }
  
  
  //=========================
  // SIMULATION
  //=========================
  vector<Type> y_sim(y_obs.size());
  for( int i=0; i<y_obs.size(); i++){
    SIMULATE {
      y_sim(i) = rnorm(field_sp(i), sigma);
    }
    REPORT(y_sim);
  } 
  
  //====================================================
  // REPORT
  //====================================================
  REPORT(field_sp);
  REPORT(field_grid);
  REPORT(z);
  REPORT(z_tilde);
  REPORT(scale);
  REPORT(alpha);
  REPORT(sigma);
  
  ADREPORT(z);
  ADREPORT(alpha);
  ADREPORT(sigma);
  
  return nll + nlp;
} 