#include "log_bf_computer.h"

//this is just a placeholder for now
double log_bf_computer(const model_t * mod){
  return(0.5*(double)(data.num_obs - mod->rank)*log((double)data.num_obs+1.00) +
         -0.5*(double)(data.num_obs - 1)*log(((double)data.num_obs)*(1.00-mod->Rsq)+1.00));
}

double coef_shrinkage_factor_computer(const model_t * mod){
  return((double)data.num_obs/((double)data.num_obs+1.00));
}

double coef_vcov_shrinkage_factor_computer(const model_t * mod){
  return((double)data.num_obs/((double)data.num_obs+1.00));
}