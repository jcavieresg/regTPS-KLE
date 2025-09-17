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
  DATA_SPARSE_MATRIX(A_grid);           // Projection matrix for grid prediction
  DATA_STRUCT(spde_mat, spde_t);        // Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  
  DATA_SCALAR(lambda_rho);
  DATA_SCALAR(lambda_sigma_u);
  
  //=========================
  //   PARAMETER SECTION
  //=========================
  PARAMETER(sigma);
  PARAMETER(rho);
  PARAMETER(sigma_u);
  
  PARAMETER_VECTOR(u);	          // spatial effects

// Transformed parameters
  Type rho_inv = pow(rho, -1);
  Type lambda_rho_inv = pow(lambda_rho, -1);

  // ===================================
  //               Priors
  // ===================================
  Type nlp = 0.0;
  
  // Prior on sigma
  nlp -= dcauchy(sigma,   Type(0.0), Type(5.0), true);

  nlp -= dweibull(rho_inv,  lambda_rho_inv, Type(1.0), true);   //
  //nlp -= 2.0 * log(rho); // Jacobian correction

  //nlp -= dexp(rho,          lambda_rho,            true);
  nlp -= dexp(sigma_u,      lambda_sigma_u,        true);
  
  Type kappa = sqrt(8)/rho;
  Type tau   = 1/(kappa*sigma_u);
  SparseMatrix<Type> Q = Q_spde(spde_mat, kappa);
  
  

  // Negative log-likelihood
  Type nll = 0.0;
  // Prior on spatial field (GMRF)
  // scaled precision: Q * tau^2
  nll += SCALE(GMRF(Q), 1/ tau)(u); // returns negative already
  
  // Likelihood
  int n = y.size();	                 // number of observations 
  vector<Type> field_sp = A_obs * u;
  for(int i = 0; i < n; i++){
    nll -= dnorm(y(i), field_sp(i), sigma, true);
  }

  // Project to prediction grid
  vector<Type> field_grid = A_grid * u;
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  
  
  // Simulate data from field_sp
  vector<Type> y_sim(n);
  for( int i=0; i<n; i++){
    SIMULATE {
      y_sim(i) = rnorm(field_sp(i), sigma);
      };
    REPORT(y_sim);
  }
  

  // REPORT 
  REPORT(sigma);
  REPORT(tau);
  REPORT(kappa);
  REPORT(sigma_u);
  REPORT(rho);
  REPORT(field_sp);
  REPORT(field_grid);
  REPORT(Q);
  
  // ADREPORT
  ADREPORT(field_sp);
  ADREPORT(field_grid);  // Enable posterior SD for uncertainty
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(kappa);

  return jnll;
}