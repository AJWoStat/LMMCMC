#ifndef SHRINKAGE_FACTOR_INTEGRATOR_H
#define SHRINKAGE_FACTOR_INTEGRATOR_H

#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<string.h>
#include"R.h"
#include"Rmath.h"
#include"R_ext/Applic.h"
#include"../cubic_root_finder/cubic_root_finder.h"
#include"../coef_prior_struct/coef_prior_struct.h"
#include"../model_struct/model_struct.h"
#include"../data_store_struct/data_store_struct.h"
#include "../log_bf_integrator/log_bf_integrator.h"


struct shrinkage_factor_integrator_struct{
  //all bayes factors are relative to base model
  //all bayes factors are mixtures of g-over-n priors where mixing is done over w=1/g
  
  //pointer to self
  struct shrinkage_factor_integrator_struct * self_ptr;
  
  double power;
  
  // pull everything else from log_bf_integrator when calling the function
  // just need to compute w_a(t) appropriately
  struct log_bf_integrator_struct * log_bf_integrator_ptr;
  
  //function for w_a
  integr_fn * get_log_w_a_fn;
  
  //integrand function
  integr_fn * shrinkage_factor_integrand_ptr;
  
};

void shrinkage_factor_integrator_struct_constructor();

void shrinkage_factor_integrator_struct_destructor();

double get_shrinkage_factor(struct model_struct * model, double power);

//computations done on log scale and then exponentiated in function
void shrinkage_factor_integrand(double * t, int n, void * ex);

void get_log_w_a_fn_return_scale(double * t, int n, void * ex);
void get_log_w_a_fn_zero_inf(double * t, int n, void * ex);
void get_log_w_a_fn_zero_scale(double * t, int n, void * ex);

#endif