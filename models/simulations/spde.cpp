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
  // DATA
  //=========================
  DATA_VECTOR(y);
  DATA_SPARSE_MATRIX(A_obs);
  DATA_SPARSE_MATRIX(A_grid);
  DATA_STRUCT(spde_mat, spde_t);
   
  DATA_SCALAR(lambda_rho);
  DATA_SCALAR(lambda_sigma_u);
   
  //=========================
  // PARAMETERS
  //=========================
  PARAMETER(sigma);
  PARAMETER(rho);
  PARAMETER(sigma_u);
   
  // Non-centered latent field
  PARAMETER_VECTOR(u_tilde);   // u_tilde ~ N(0, Q^{-1})
   
  //=========================
  // TRANSFORMS
  //=========================
  Type kappa = sqrt(Type(8.0)) / rho;
  Type tau   = Type(1.0) / (kappa * sigma_u);
   
  SparseMatrix<Type> Q = Q_spde(spde_mat, kappa);
   
  // Transform to centered field
  vector<Type> u = u_tilde / tau;
   
  //=========================
  // PRIORS
  //=========================
  Type nlp = 0.0;
   
  // Prior on observation SD
  nlp -= dcauchy(sigma, Type(0.0), Type(5.0), true);
   
  // Prior on range (via inverse rho)
  Type rho_inv = pow(rho, -1);
  Type lambda_rho_inv = pow(lambda_rho, -1);
  nlp -= dweibull(rho_inv, lambda_rho_inv, Type(1.0), true);
  // Jacobian already implicit in transformation
   
  // Prior on marginal SD
  nlp -= dexp(sigma_u, lambda_sigma_u, true);
   
  //=========================
  // LATENT FIELD PRIOR
  //=========================
  Type nll = 0.0;
   
  // u_tilde ~ N(0, Q^{-1})
  nll += GMRF(Q)(u_tilde);
   
  //=========================
  // LIKELIHOOD
  //=========================
  vector<Type> field_sp = A_obs * u;
   
  for(int i = 0; i < y.size(); i++){
    nll -= dnorm(y(i), field_sp(i), sigma, true);
  }
   
  vector<Type> field_grid = A_grid * u;
   
  //=========================
  // JOINT NLL
  //=========================
  Type jnll = nll + nlp;
   
  //=========================
  // SIMULATION
  //=========================
  vector<Type> y_sim(y.size());
  SIMULATE {
    for(int i = 0; i < y.size(); i++){
      y_sim(i) = rnorm(field_sp(i), sigma);
    } 
    REPORT(y_sim);
  }
   
  //=========================
  // REPORT
  //=========================
  REPORT(sigma);
  REPORT(rho);
  REPORT(sigma_u);
  REPORT(kappa);
  REPORT(tau);
  REPORT(u);
  REPORT(u_tilde);
  REPORT(field_sp);
  REPORT(field_grid);
   
  ADREPORT(sigma);
  ADREPORT(rho);
  ADREPORT(sigma_u);
  ADREPORT(kappa);
  ADREPORT(tau);
  ADREPORT(field_grid);
  ADREPORT(u);
  ADREPORT(u_tilde);
   
  return jnll;
}