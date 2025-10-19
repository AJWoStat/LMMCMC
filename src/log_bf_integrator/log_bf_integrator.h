#ifndef LOG_BF_INTEGRATOR_H
#define LOG_BF_INTEGRATOR_H

#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<string.h>
#include"R.h"
#include"Rmath.h"
#include"R_ext/Applic.h"
#include"../data_store_struct/data_store_struct.h"
#include"../cubic_root_finder/cubic_root_finder.h"
#include"../coef_prior_struct/coef_prior_struct.h"
#include"../model_struct/model_struct.h"
#include"../data_store_struct/data_store_struct.h"


typedef double get_log_bf_prior_specific_fn_t(struct model_struct * model);

//integrand function
integr_fn * log_bf_integrand; //computations done on log scale and then exponentiated in function, return is log-scale
//function that does the calling
get_log_bf_prior_specific_fn_t * get_log_bf_prior_specific;

struct log_bf_integrator_struct{
  //all bayes factors are relative to base model
  //all bayes factors are mixtures of g-over-n priors where mixing is done over w=1/g
  
  //pointer to self
  struct log_bf_integrator_struct * self_ptr;
  
  //coefficient prior
  //prior paramaters
  double scale; //prior scale of w = 1/g
  double shape_0; //prior shape at 0 for w=1/g
  double shape_1; //prior shape at infty or scale for w=1/g (where applicable)
  //type of bayes factor
  char bf_type[13];
    //seven types supported here, three are just special cases
    //g-prior : no prior on w=1/g, g is equal to parameter scale
    //zellner-siow : special case of gamma prior with shape_0=0.5 and scale = 1
    //hyper-g : special case of beta-prime with shape_0 = (par-2)/2, shape_1 = 1, and scale = 1
    //intrinsic : special case of scaled-beta with shape_0 = shape_1 = 1/2 and scale = 1
    //gamma : prior on w=1/g is const*w^(shape_0-1)*exp(-w/scale) where w>0 (shape_1 is not involved here)
    //beta-prime : prior on w=1/g is const* w^(shape_0-1)/(scale+w)^(shape_0+shape_1) where w>0
    //scaled-beta : prior on w=1/g is const w^(shape_0-1)(scale-w)^(shape_1-1) where 0<w<scale
  
  //data
  //common data
  double n; //number of observations
  double p_0; //base number of covariates
  double r_0; //base model rank
  int eff; //scale g by effective number of parameters? 1 yes, 0 no. 
    //if 1, then normal prior precision is w*p/n*X^T*X where w=1/g
    //if 0, then normal prior precision is w/n*X^T*X where w=1/g
  //model specific data
  double p; //model number of covariates, includes base count
  double r; //model rank, includes base rank
  double q; //scaling of g effective number of covariates
    //q = eff==1 ? r : 1;
  double Rsq; //model Rsq
  
  //integration parameters
  //integral is always from 0 to 1
  //transformation is w(t) = f(g(t)) where f is bf_type specific
  //g(t) = a*t^gamma_0/(a*t^gamma_0+(1-a)*(1-t)^gamma_0)
  
  //transformation parameters for integral
  double gamma_0;
  double gamma_1;
  double t_0; //should set to 0.5
  double wat0; //will always be >0
  double a; //will always be >0
  double max_log_integrand; //subtracted from log integrand before exponentiation and integration
  double gamma_integrand_epsilon; 
  int log_eval;
  
  //coefficient parameters and roots of cubic polynomial
  double A[3];
  double B[3];
  double C[2];
  double D[2];
  double poly_coef[4];
  double poly_root[6];
  double poly_root_epsilon;
  
  //variables for Rdqags
  void * ex;
  double lower;
  double upper;
  double epsabs;
  double epsrel;
  double result;
  double abserr;
  int neval;
  int ier;
  int limit;
  int lenw;
  int last;
  double * lower_ptr;
  double * upper_ptr;
  double * epsabs_ptr;
  double * epsrel_ptr;
  double * result_ptr;
  double * abserr_ptr;
  int * neval_ptr;
  int * ier_ptr;
  int * limit_ptr;
  int * lenw_ptr;
  int * last_ptr;
  int * iwork;
  double * work;
};

extern struct log_bf_integrator_struct log_bf_integrator;

void log_bf_integrator_struct_constructor(struct data_store_struct * data, struct coef_prior_struct * coef_prior);

void log_bf_integrator_struct_destructor();

double get_log_bf(struct model_struct * model);

double get_log_bf_g_prior(struct model_struct * model);
double get_log_bf_gamma(struct model_struct * model);
double get_log_bf_beta_prime(struct model_struct * model);
double get_log_bf_scaled_beta(struct model_struct * model);

void log_bf_g_prior_integrand(double * t, int n, void * ex);
void log_bf_gamma_integrand(double * t, int n, void * ex);
void log_bf_beta_prime_integrand(double * t, int n, void * ex);
void log_bf_scaled_beta_integrand(double * t, int n, void * ex);

//g-prior : no prior on w=1/g, g is equal to parameter scale
//zellner-siow : special case of gamma prior with shape_0=0.5 and scale = 1
//hyper-g : special case of beta-prime with shape_0 = (par-2)/2, shape_1 = 1, and scale = 1
//intrinsic : special case of scaled-beta with shape_0 = shape_1 = 1/2 and scale = 1
//gamma : prior on w=1/g is const*w^(shape_0-1)*exp(-w/scale) where w>0 (shape_1 is not involved here)
//beta-prime : prior on w=1/g is const* w^(shape_0-1)/(scale+w)^(shape_0+shape_1) where w>0
//scaled-beta : prior on w=1/g is const w^(shape_0-1)(scale-w)^(shape_1-1) where 0<w<scale
#endif