#include"shrinkage_factor_integrator.h"

static struct shrinkage_factor_integrator_struct shrinkage_factor_integrator;

void shrinkage_factor_integrator_struct_constructor()
{
  //pointer to self
  shrinkage_factor_integrator.self_ptr = &shrinkage_factor_integrator;
  
  shrinkage_factor_integrator.power = 1.00;
  
  shrinkage_factor_integrator.log_bf_integrator_ptr = &log_bf_integrator;
  
  if(strcmp(shrinkage_factor_integrator.log_bf_integrator_ptr->bf_type, "gamma")==0 || 
     strcmp(shrinkage_factor_integrator.log_bf_integrator_ptr->bf_type, "zellner-siow") ||
     strcmp(shrinkage_factor_integrator.log_bf_integrator_ptr->bf_type, "beta-prime")==0 || 
     strcmp(shrinkage_factor_integrator.log_bf_integrator_ptr->bf_type, "hyper-g")==0)
  {
    shrinkage_factor_integrator.get_log_w_a_fn = &get_log_w_a_fn_zero_inf;
  }else if(strcmp(shrinkage_factor_integrator.log_bf_integrator_ptr->bf_type, "scaled-beta")==0 ||
    strcmp(shrinkage_factor_integrator.log_bf_integrator_ptr->bf_type, "intrinsic")==0)
  {
    shrinkage_factor_integrator.get_log_w_a_fn = &get_log_w_a_fn_zero_scale;
  }else
  {
    shrinkage_factor_integrator.get_log_w_a_fn = &get_log_w_a_fn_return_scale;
  }
  
  shrinkage_factor_integrator.shrinkage_factor_integrand_ptr = &shrinkage_factor_integrand;
  
  return;
}

void shrinkage_factor_integrator_struct_destructor(){
  return;
}

double get_shrinkage_factor(struct model_struct * model, double power){
  //make switch for g-prior to not do integral
  double result, out;
  
  struct log_bf_integrator_struct * par = shrinkage_factor_integrator.log_bf_integrator_ptr;
  
  //switch for g-prior bf_type
  if(strcmp(log_bf_integrator.bf_type, "g-prior")!=0){
    //do integration here
    Rdqags(*shrinkage_factor_integrator.shrinkage_factor_integrand_ptr, (void*) shrinkage_factor_integrator.self_ptr, 
           &(par->lower), &(par->upper),
           &(par->epsabs), &(par->epsrel), 
           &result, 
           &(par->abserr), &(par->neval), &(par->ier), 
           &(par->limit), &(par->lenw),
           &(par->last), par->iwork, par->work);
    //we do no error catching in the integration. It would be shocking if the integral did not work.
    //we should probably do a check on ier here

    out = log(result)+par->max_log_integrand;
          
  }else{
    //nothing else to set
    (*shrinkage_factor_integrator.shrinkage_factor_integrand_ptr)(&out, 1, (void*) shrinkage_factor_integrator.self_ptr);
  }        
  return(exp(out - model->log_BF0));
}

//computations done on log scale and then exponentiated in function
void shrinkage_factor_integrand(double * t, int n, void * ex){
  struct shrinkage_factor_integrator_struct * par = (struct shrinkage_factor_integrator_struct*) ex;
  double * t_for_log_w_a = malloc(n*sizeof(double));
  for(int i=0; i<n; i++) t_for_log_w_a[i] = t[i];
  (*(par->get_log_w_a_fn))(t_for_log_w_a, n, ex);
  double n_obs = par->log_bf_integrator_ptr->n;
  double q = par->log_bf_integrator_ptr->q;
  double log_q_over_n = log(q) - log(n_obs);
  double power = par->power;
  for(int i=0; i<n; i++) t_for_log_w_a[i] = -power*log1pexp(log_q_over_n+t_for_log_w_a[i]);
  par->log_bf_integrator_ptr->log_eval = 1;
  (*(par->log_bf_integrator_ptr->log_bf_integrand))(t, n, (void*) par->log_bf_integrator_ptr);
  for(int i=0; i<n; i++){
    t[i] += t_for_log_w_a[i];
    t[i] = exp(t[i]);
  }
  free(t_for_log_w_a); t_for_log_w_a = NULL;
  return;
}

void get_log_w_a_fn_return_scale(double * t, int n, void * ex){
  struct log_bf_integrator_struct * par = ((struct shrinkage_factor_integrator_struct*) ex)->log_bf_integrator_ptr;
  double log_scale = log(par->scale);
  for(int i=0; i<n; i++) t[i] = log_scale;
  return;
}

void get_log_w_a_fn_zero_inf(double * t, int n, void * ex){
  struct log_bf_integrator_struct * par = ((struct shrinkage_factor_integrator_struct*) ex)->log_bf_integrator_ptr;
  double g_0 = par->gamma_0;
  double g_1 = par->gamma_1;
  double log_a = log(par->a);
  double log_scale = log(par->scale);
  double log_scale_plus_a = log_scale + log_a;
  for(int i=0; i<n; i++) t[i] = log_scale_plus_a + g_0*log(t[i]) - g_1*log1p(-t[i]);
  return;
}

void get_log_w_a_fn_zero_scale(double * t, int n, void * ex){
  struct log_bf_integrator_struct * par = ((struct shrinkage_factor_integrator_struct*) ex)->log_bf_integrator_ptr;
  double g_0 = par->gamma_0;
  double g_1 = par->gamma_1;
  double log_a = log(par->a);
  double log_scale = log(par->scale);
  double log_scale_plus_a = log_scale + log_a;
  for(int i=0; i<n; i++){
    t[i] = -log1pexp(log_scale_plus_a + log_a + g_0*log(t[i]) - g_1*log1p(-t[i]));
  }
  return;
}

