#include <TMB.hpp>

// dcauchy for hyperparameters
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(y);
  DATA_MATRIX(Phi_kle_sp); // The pre-computed KLE basis at observation points
  DATA_MATRIX(Phi_kle_grid); // The pre-computed KLE basis at grid points
  DATA_VECTOR(S_diag_truncated); // Eigenvalues of S (truncated)
  DATA_INTEGER(M_P_null_space); // Number of polynomial modes
 
  // PARAMETERS
  PARAMETER_VECTOR(Z); // KLE coefficients (length M_truncation)
  PARAMETER(logsigma); // Log of observation noise SD
  PARAMETER(logalpha); // Log of regularization parameter

  // MODEL SETUP
  Type sigma = exp(logsigma);
  Type alpha = exp(logalpha);
  
  //==================================================
  // PRIORS
  //==================================================
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  nlp -= dcauchy(sigma,   Type(0.0),   Type(1.0));
  nlp += logsigma; // Jacobian
  
  // Prior on alpha 
  nlp -= dcauchy(alpha, Type(0.0), Type(1.0), true);
  nlp += logalpha; // Jacobian
  
  // Prior for Z_k: Z_k ~ N(0, lambda_k)
  // lambda_k = 1 / (1 + alpha * S_diag(k))
  // For the unpenalized modes (polynomial), S_diag(k) = 0, so lambda_k = 1
  for(int k=0; k < M_P_null_space; ++k){
    nlp -= dnorm(Z(k), Type(0.0), Type(1.0), true);
  }
  
  // Prior for penalized modes
  Type prior_sd;
  for(int k=M_P_null_space; k < S_diag_truncated.size(); ++k){
    Type S_diag_k = S_diag_truncated(k);
    prior_sd = sqrt(Type(1.0) / (Type(1.0) + alpha * S_diag_k));
    nlp -= dnorm(Z(k), Type(0.0), prior_sd, true);
  }

  
  // Z has M_truncation elements
  vector<Type> field_sp = Phi_kle_sp * Z;
  vector<Type> field_grid = Phi_kle_grid * Z; // For reporting
  
  
  //====================================================
  // Likelihood
  vector<Type> log_lik(y.size());
  
  for( int i = 0; i<y.size(); i++){
    log_lik(i) = dnorm(y(i), field_sp(i), sigma, true);
  }
  // Here in nll the likelihood is stored
  Type nll = -log_lik.sum(); // total NLL
  
  // Total negative log-like
  Type jnll = nll + nlp;
  
  
  // Simule data from the mu 
  vector<Type> y_sim(y.size());
  for( int i=0; i<y.size(); i++){
    SIMULATE {
      y_sim(i) = rnorm(field_sp(i), sigma);
      REPORT(y_sim);
    }
  }
  // RECONSTRUCTION

  // REPORTING
  REPORT(field_grid);
  REPORT(sigma);
  REPORT(alpha);
  REPORT(prior_sd);
  return jnll;
}