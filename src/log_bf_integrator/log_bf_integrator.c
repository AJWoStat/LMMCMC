#include"log_bf_integrator.h"

// static struct log_bf_integrator_struct log_bf_integrator;
struct log_bf_integrator_struct log_bf_integrator;

void log_bf_integrator_struct_constructor(struct data_store_struct * data, struct coef_prior_struct * coef_prior)
{
  //pointer to self
  log_bf_integrator.self_ptr = &log_bf_integrator;
  
  //coefficient prior
  //prior paramaters
  log_bf_integrator.scale = coef_prior->scale; //prior scale of w = 1/g
  log_bf_integrator.shape_0 = coef_prior->shape_0; //prior shape at 0 for w=1/g
  log_bf_integrator.shape_1 = coef_prior->shape_1; //prior shape at infty or 1 for w=1/g (where applicable)
  //type of bayes factor
  strcpy(log_bf_integrator.bf_type, coef_prior->bf_type);
  
  //integrand and get_log_bf functions
  if(strcmp(log_bf_integrator.bf_type, "gamma")==0 || strcmp(log_bf_integrator.bf_type, "zellner-siow")==0){
    log_bf_integrand = &log_bf_gamma_integrand; 
    get_log_bf_prior_specific = &get_log_bf_gamma;
  }else if(strcmp(log_bf_integrator.bf_type, "beta-prime")==0 || strcmp(log_bf_integrator.bf_type, "hyper-g")==0){
    log_bf_integrand = &log_bf_beta_prime_integrand; 
    get_log_bf_prior_specific = &get_log_bf_beta_prime;
  }else if(strcmp(log_bf_integrator.bf_type, "scaled-beta")==0 || strcmp(log_bf_integrator.bf_type, "intrinsic")==0){
    log_bf_integrand = &log_bf_scaled_beta_integrand;
    get_log_bf_prior_specific = &get_log_bf_scaled_beta;
  }else{
    log_bf_integrand = &log_bf_g_prior_integrand;
    get_log_bf_prior_specific = &get_log_bf_g_prior;
  }
  
  //data
  //common data
  log_bf_integrator.n = (double)(data->num_obs); //number of observations
  log_bf_integrator.p_0 = (double)(data->base_model_size); //base number of covariates
  log_bf_integrator.r_0 = (double)(data->base_model_rank); //base model rank
  log_bf_integrator.eff = coef_prior->eff; //scale g by effective number of parameters? 1 yes, 0 no. 
  //if 1, then normal prior precision is w*p/n*X^T*X where w=1/g
  //if 0, then normal prior precision is w/n*X^T*X where w=1/g
  //model specific data
  log_bf_integrator.p = (double)(log_bf_integrator.p_0); //model number of covariates, includes base count
  log_bf_integrator.r = (double)(log_bf_integrator.r_0); //model rank, includes base rank
  log_bf_integrator.q = log_bf_integrator.eff==1 ? log_bf_integrator.r : 1.00; //scaling of g effective number of covariates
  //q = eff==1 ? r : 1;
  log_bf_integrator.Rsq = 0.00; //model Rsq
  
  //integration parameters
  //integral is always from 0 to 1
  //transformation is w(t) = f(g(t)) where f is bf_type specific
  //g(t) = a*t^gamma_0/(a*t^gamma_0+(1-a)*(1-t)^gamma_0)
  
  //transformation parameters for integral
  log_bf_integrator.gamma_0 = 1.00;
  log_bf_integrator.gamma_1 = 1.00;
  log_bf_integrator.t_0 = 0.5; //good default value
  log_bf_integrator.wat0 = 1.00; //will always be >0
  log_bf_integrator.a = 1.00; //will always be >0
  log_bf_integrator.max_log_integrand = 1.00; //subtracted from log integrand before exponentiation and integration
  log_bf_integrator.gamma_integrand_epsilon = pow(DBL_EPSILON, 0.5);
  log_bf_integrator.log_eval = 0;
  
  //coefficient parameters of cubic polynomial
  log_bf_integrator.A[0] = 1.00;
  log_bf_integrator.A[1] = 1.00;
  log_bf_integrator.A[2] = 1.00;
  log_bf_integrator.B[0] = 1.00;
  log_bf_integrator.B[1] = 1.00;
  log_bf_integrator.B[2] = 1.00;
  log_bf_integrator.C[0] = 1.00;
  log_bf_integrator.C[1] = 1.00;
  log_bf_integrator.D[0] = 1.00;
  log_bf_integrator.D[1] = 1.00;
  log_bf_integrator.poly_coef[0] = 1.00;
  log_bf_integrator.poly_coef[1] = 1.00;
  log_bf_integrator.poly_coef[2] = 1.00;
  log_bf_integrator.poly_coef[3] = 1.00;
  log_bf_integrator.poly_root[0] = 1.00;
  log_bf_integrator.poly_root[1] = 1.00;
  log_bf_integrator.poly_root[2] = 1.00;
  log_bf_integrator.poly_root[3] = 1.00;
  log_bf_integrator.poly_root[4] = 1.00;
  log_bf_integrator.poly_root[5] = 1.00;
  log_bf_integrator.poly_root_epsilon = pow(DBL_EPSILON, 0.75);
  //variables for Rdqags
  log_bf_integrator.ex = log_bf_integrator.self_ptr;
  log_bf_integrator.lower = 0.00;
  log_bf_integrator.upper = 1.00;
  log_bf_integrator.epsabs = pow(DBL_EPSILON, 0.5);
  log_bf_integrator.epsrel = pow(DBL_EPSILON, 0.5);
  log_bf_integrator.result = 0.00;
  log_bf_integrator.abserr = 0.00;
  log_bf_integrator.neval = 0;
  log_bf_integrator.ier = 0;
  log_bf_integrator.limit = 1000;
  log_bf_integrator.lenw = 4*(log_bf_integrator.limit);
  log_bf_integrator.last = 0;
  log_bf_integrator.lower_ptr = &log_bf_integrator.lower;
  log_bf_integrator.upper_ptr = &log_bf_integrator.upper;
  log_bf_integrator.epsabs_ptr = &log_bf_integrator.epsabs;
  log_bf_integrator.epsrel_ptr = &log_bf_integrator.epsrel;
  log_bf_integrator.result_ptr = &log_bf_integrator.result;
  log_bf_integrator.abserr_ptr = &log_bf_integrator.abserr;
  log_bf_integrator.neval_ptr = &log_bf_integrator.neval;
  log_bf_integrator.ier_ptr = &log_bf_integrator.ier;
  log_bf_integrator.limit_ptr = &log_bf_integrator.limit;
  log_bf_integrator.lenw_ptr = &log_bf_integrator.lenw;
  log_bf_integrator.last_ptr = &log_bf_integrator.last;
  log_bf_integrator.iwork = (int*)calloc(log_bf_integrator.limit, sizeof(int));
  log_bf_integrator.work = (double*)calloc(log_bf_integrator.lenw, sizeof(double));
  return;
}

void log_bf_integrator_struct_destructor(){
  free(log_bf_integrator.iwork); log_bf_integrator.iwork = NULL;
  free(log_bf_integrator.work); log_bf_integrator.work = NULL;
  return;
}

double get_log_bf(struct model_struct * model){
  //note use of r and not p for q
  log_bf_integrator.p = (double)(model->size);
  log_bf_integrator.r = (double)(model->rank);
  log_bf_integrator.q = log_bf_integrator.eff == 1 ? log_bf_integrator.r : 1.00;
  //HERE is a problem?
  log_bf_integrator.Rsq = model->Rsq;//1.00-(1.00 - model->Rsq*data.intercept_only_model_SSE)/data.base_model_SSE;
  // log_bf_integrator.Rsq = fmin(log_bf_integrator.Rsq, 1.00);
  // log_bf_integrator.Rsq = fmax(log_bf_integrator.Rsq, 0.00);
  double out = (*get_log_bf_prior_specific)(model);
  return(out);
}

double get_log_bf_g_prior(struct model_struct * model){
  double out = 0.00;
  (*log_bf_integrand)(&out, 1, log_bf_integrator.ex);
  return(out);
}

void log_bf_g_prior_integrand(double * t, int n, void * ex)//actually unused
{ 
  //note use of r and not p
  struct log_bf_integrator_struct * par = (struct log_bf_integrator_struct*) ex;
  double scale_times_q_over_n = par->scale * par->q / par->n ;
  t[0] = 0.5*(par->n - par->r)*log1p(scale_times_q_over_n) +
    -0.5*(par->n - par->r_0)*log1p(scale_times_q_over_n - par->Rsq) +
    0.5*(par->r - par->r_0)*log(scale_times_q_over_n);
  return;
}


double get_log_bf_gamma(struct model_struct * model){
  log_bf_integrator.gamma_0 = 2.00/(2.00*log_bf_integrator.shape_0 + log_bf_integrator.r-log_bf_integrator.r_0) + 1.00;//log_bf_integrator.shape_0>0.50 ? 1.00 : 2.00/log_bf_integrator.shape_0;//4.00/(2.00*log_bf_integrator.shape_0 + log_bf_integrator.r-log_bf_integrator.r_0);//
  log_bf_integrator.gamma_1 = 1.00;
  //note using ranks instead of sizes
  double rat = (log_bf_integrator.gamma_0*(1.00-log_bf_integrator.t_0) + log_bf_integrator.gamma_1*log_bf_integrator.t_0)/(log_bf_integrator.t_0*(1.00-log_bf_integrator.t_0));
  log_bf_integrator.A[2] = 0.00;
  log_bf_integrator.A[1] = 0.5*rat*log_bf_integrator.q*log_bf_integrator.n*((log_bf_integrator.r-log_bf_integrator.r_0)-log_bf_integrator.Rsq*(log_bf_integrator.n-log_bf_integrator.r_0));
  log_bf_integrator.A[0] = 0.5*rat*log_bf_integrator.n*log_bf_integrator.n*((log_bf_integrator.r-log_bf_integrator.r_0)*(1.00-log_bf_integrator.Rsq));
  log_bf_integrator.B[2] = log_bf_integrator.q*log_bf_integrator.q;
  log_bf_integrator.B[1] = log_bf_integrator.q*log_bf_integrator.n*(2.00-log_bf_integrator.Rsq);
  log_bf_integrator.B[0] = log_bf_integrator.n*log_bf_integrator.n*(1.00-log_bf_integrator.Rsq);
  log_bf_integrator.C[1] = -rat/log_bf_integrator.scale;
  log_bf_integrator.C[0] = log_bf_integrator.shape_0*rat + (log_bf_integrator.gamma_1-log_bf_integrator.gamma_0)/(log_bf_integrator.gamma_0*(1.00-log_bf_integrator.t_0) + log_bf_integrator.gamma_1*log_bf_integrator.t_0) - 1.00/log_bf_integrator.t_0 + 1.00/(1.00-log_bf_integrator.t_0);
  log_bf_integrator.D[1] = 0.00;
  log_bf_integrator.D[0] = 1.00;
  log_bf_integrator.poly_coef[0] = log_bf_integrator.A[0]*log_bf_integrator.D[0]+log_bf_integrator.B[0]*log_bf_integrator.C[0];
  log_bf_integrator.poly_coef[1] = log_bf_integrator.A[1]*log_bf_integrator.D[0]+log_bf_integrator.B[1]*log_bf_integrator.C[0]+log_bf_integrator.A[0]*log_bf_integrator.D[1]+log_bf_integrator.B[0]*log_bf_integrator.C[1];
  log_bf_integrator.poly_coef[2] = log_bf_integrator.A[2]*log_bf_integrator.D[0]+log_bf_integrator.B[2]*log_bf_integrator.C[0]+log_bf_integrator.A[1]*log_bf_integrator.D[1]+log_bf_integrator.B[1]*log_bf_integrator.C[1];
  log_bf_integrator.poly_coef[3] = log_bf_integrator.A[2]*log_bf_integrator.D[1]+log_bf_integrator.B[2]*log_bf_integrator.C[1];
  cubic_root_finder(log_bf_integrator.poly_coef, log_bf_integrator.poly_root);
  // we do no error catching in the cubic root solving. It would be shocking if there was an issue.
  // we should probably do a check on INFO here

  for(int i=0; i<3; i++){
    if(log_bf_integrator.poly_root[i]>0 && 
       fabs(log_bf_integrator.poly_root[i+3])<log_bf_integrator.poly_root_epsilon){
      log_bf_integrator.wat0 = log_bf_integrator.poly_root[i];
      break;
    }
  }
  log_bf_integrator.a = log_bf_integrator.wat0/log_bf_integrator.scale*exp(logspace_add(log(1.00-log_bf_integrator.t_0)*log_bf_integrator.gamma_1, log(log_bf_integrator.t_0)*log_bf_integrator.gamma_0));
  // log_bf_integrator.a = log_bf_integrator.wat0*pow(1.00-log_bf_integrator.t_0, log_bf_integrator.gamma_1)/(pow(log_bf_integrator.t_0, log_bf_integrator.gamma_0)*log_bf_integrator.scale);
  log_bf_integrator.max_log_integrand = 0.00;
  double initial_val = log_bf_integrator.t_0;
  log_bf_integrator.log_eval = 1;
  (*log_bf_integrand)(&initial_val, 1, log_bf_integrator.ex);
  log_bf_integrator.max_log_integrand = initial_val;
  log_bf_integrator.log_eval = 0;
  Rdqags(*log_bf_integrand, log_bf_integrator.ex, &log_bf_integrator.lower, &log_bf_integrator.upper,
         &log_bf_integrator.epsabs, &log_bf_integrator.epsrel, &log_bf_integrator.result, &log_bf_integrator.abserr,
         &log_bf_integrator.neval, &log_bf_integrator.ier, &log_bf_integrator.limit, &log_bf_integrator.lenw,
         &log_bf_integrator.last, log_bf_integrator.iwork, log_bf_integrator.work);
  //we do no error catching in the integration. It would be shocking if the integral did not work.
  //we should probably do a check on ier here
  
  return(log(log_bf_integrator.result)+log_bf_integrator.max_log_integrand);
}

void log_bf_gamma_integrand(double * t, int n, void * ex)
{
  struct log_bf_integrator_struct * par = (struct log_bf_integrator_struct*) ex;
  double n_obs = par->n;
  double log_n = log(n_obs);
  double nmrover2 = 0.5*(n_obs-par->r);
  double nmr0over2 = 0.5*(n_obs-par->r_0);
  double g_0 = par->gamma_0;
  double g_1 = par->gamma_1;
  double a = par->a;
  double alpha = par->shape_0;
  double scale = par->scale;
  double saq = scale*a*par->q;
  double log_saq = log(saq);
  double Rsq = par->Rsq;
  double n1mRsq = n_obs*(1.00-Rsq);
  double max_log = par->max_log_integrand;
  double log_const = 0.50*(par->r-par->r_0)*log(saq)+alpha*log(a) - lgammafn(alpha) - max_log;
  double power_t = g_0*(0.50*(par->r-par->r_0)+alpha)-1.00;
  double power_1mt = g_1*alpha+1.00;
  double log_t, log_1mt, g_0log_t, g_1log_1mt;
  double log_gamma_integrand_epsilon = log(par->gamma_integrand_epsilon);
  double C = a*fmax(1.00, g_0)-(alpha+1.00)*log(a) -1.00+0.5*alpha*alpha;
  double delta = 1.00-a*exp(-alpha-sqrt(2.00*(C-log_gamma_integrand_epsilon)));
  for(int i=0; i<n; i++){
    if(t[i]<delta){
      log_t = log(t[i]);
      log_1mt = log1p(-t[i]);
      g_0log_t = g_0*log_t;
      g_1log_1mt = g_1*log_1mt;
      t[i] = log_const+
        nmrover2*logspace_add(log_n + g_1log_1mt, log_saq + g_0log_t)+
        -nmr0over2*logspace_add(log_n + g_1log_1mt + log1p(-Rsq), log_saq + g_0log_t)+
        log(g_0*(1.00-t[i])+g_1*t[i])+
        power_t*log_t+
        -power_1mt*log_1mt+
        -scale*a*exp(g_0log_t - g_1log_1mt);
    }else{
      t[i] = log_gamma_integrand_epsilon;
    }
  }
  if(par->log_eval==0){
    for(int i=0; i<n; i++) t[i] = exp(t[i]);
  }
  return;
}


double get_log_bf_beta_prime(struct model_struct * model){
  log_bf_integrator.gamma_0 = 2.00/(2.00*log_bf_integrator.shape_0 + log_bf_integrator.r-log_bf_integrator.r_0) + 1.00;//log_bf_integrator.shape_0>0.50 ? 1.00 : 2.00/log_bf_integrator.shape_0;//4.00/(2.00*log_bf_integrator.shape_0 + log_bf_integrator.r-log_bf_integrator.r_0);//
  log_bf_integrator.gamma_1 = 1.00/log_bf_integrator.shape_1 + 1.00;//log_bf_integrator.shape_1>1.00 ? 1.00 : 2.00/log_bf_integrator.shape_1;//2.00/log_bf_integrator.shape_1;//
  //note using ranks instead of sizes
  double rat = (log_bf_integrator.gamma_0*(1.00-log_bf_integrator.t_0) + log_bf_integrator.gamma_1*log_bf_integrator.t_0)/(log_bf_integrator.t_0*(1.00-log_bf_integrator.t_0));
  log_bf_integrator.A[2] = 0.00;
  log_bf_integrator.A[1] = 0.5*rat*log_bf_integrator.q*log_bf_integrator.n*((log_bf_integrator.r-log_bf_integrator.r_0)-log_bf_integrator.Rsq*(log_bf_integrator.n-log_bf_integrator.r_0));
  log_bf_integrator.A[0] = 0.5*rat*log_bf_integrator.n*log_bf_integrator.n*((log_bf_integrator.r-log_bf_integrator.r_0)*(1.00-log_bf_integrator.Rsq));
  log_bf_integrator.B[2] = log_bf_integrator.q*log_bf_integrator.q;
  log_bf_integrator.B[1] = log_bf_integrator.q*log_bf_integrator.n*(2.00-log_bf_integrator.Rsq);
  log_bf_integrator.B[0] = log_bf_integrator.n*log_bf_integrator.n*(1.00-log_bf_integrator.Rsq);
  double num_add = (log_bf_integrator.gamma_1-log_bf_integrator.gamma_0)/(log_bf_integrator.gamma_0*(1.00-log_bf_integrator.t_0) + log_bf_integrator.gamma_1*log_bf_integrator.t_0) - 1.00/log_bf_integrator.t_0 + 1.00/(1.00-log_bf_integrator.t_0);
  log_bf_integrator.C[1] = -rat*log_bf_integrator.shape_1 + num_add;
  log_bf_integrator.C[0] = log_bf_integrator.scale*log_bf_integrator.shape_0*rat + log_bf_integrator.scale*num_add;
  log_bf_integrator.D[1] = 1.00;
  log_bf_integrator.D[0] = log_bf_integrator.scale;
  log_bf_integrator.poly_coef[0] = log_bf_integrator.A[0]*log_bf_integrator.D[0]+log_bf_integrator.B[0]*log_bf_integrator.C[0];
  log_bf_integrator.poly_coef[1] = log_bf_integrator.A[1]*log_bf_integrator.D[0]+log_bf_integrator.B[1]*log_bf_integrator.C[0]+log_bf_integrator.A[0]*log_bf_integrator.D[1]+log_bf_integrator.B[0]*log_bf_integrator.C[1];
  log_bf_integrator.poly_coef[2] = log_bf_integrator.A[2]*log_bf_integrator.D[0]+log_bf_integrator.B[2]*log_bf_integrator.C[0]+log_bf_integrator.A[1]*log_bf_integrator.D[1]+log_bf_integrator.B[1]*log_bf_integrator.C[1];
  log_bf_integrator.poly_coef[3] = log_bf_integrator.A[2]*log_bf_integrator.D[1]+log_bf_integrator.B[2]*log_bf_integrator.C[1];
  cubic_root_finder(log_bf_integrator.poly_coef, log_bf_integrator.poly_root);
  //we do no error catching in the cubic root solving. It would be shocking if there was an issue.
  //we should probably do a check on INFO here
  
  for(int i=0; i<3; i++){
    if(log_bf_integrator.poly_root[i]>0 && 
       fabs(log_bf_integrator.poly_root[i+3])<log_bf_integrator.poly_root_epsilon){
      log_bf_integrator.wat0 = log_bf_integrator.poly_root[i];
      break;
    }
  }
  log_bf_integrator.a = log_bf_integrator.wat0/log_bf_integrator.scale*exp(logspace_add(log(1.00-log_bf_integrator.t_0)*log_bf_integrator.gamma_1, log(log_bf_integrator.t_0)*log_bf_integrator.gamma_0));
  log_bf_integrator.max_log_integrand = 0.00;
  double initial_val = log_bf_integrator.t_0;
  log_bf_integrator.log_eval = 1;
  (*log_bf_integrand)(&initial_val, 1, log_bf_integrator.ex);
  log_bf_integrator.max_log_integrand = initial_val;
  log_bf_integrator.log_eval = 0;
  Rdqags(*log_bf_integrand, log_bf_integrator.ex, &log_bf_integrator.lower, &log_bf_integrator.upper,
         &log_bf_integrator.epsabs, &log_bf_integrator.epsrel, &log_bf_integrator.result, &log_bf_integrator.abserr,
         &log_bf_integrator.neval, &log_bf_integrator.ier, &log_bf_integrator.limit, &log_bf_integrator.lenw,
         &log_bf_integrator.last, log_bf_integrator.iwork, log_bf_integrator.work);
         //we do no error catching in the integration. It would be shocking if the integral did not work.
         //we should probably do a check on ier here
         
         return(log(log_bf_integrator.result)+log_bf_integrator.max_log_integrand);
}

void log_bf_beta_prime_integrand(double * t, int n, void * ex)
{
  struct log_bf_integrator_struct * par = (struct log_bf_integrator_struct*) ex;
  double n_obs = par->n;
  double log_n = log(n_obs);
  double nmrover2 = 0.5*(n_obs-par->r);
  double nmr0over2 = 0.5*(n_obs-par->r_0);
  double g_0 = par->gamma_0;
  double g_1 = par->gamma_1;
  double a = par->a;
  double log_a = log(a);
  double alpha = par->shape_0;
  double beta = par->shape_1;
  double scale = par->scale;
  double saq = scale*a*par->q;
  double log_saq = log(saq);
  double Rsq = par->Rsq;
  double n1mRsq = n_obs*(1.00-Rsq);
  double max_log = par->max_log_integrand;
  double log_const = 0.5*(par->r-par->r_0)*log_saq + alpha*log_a - lbeta(alpha, beta) - max_log;
  double power_t = g_0*(0.5*(par->r-par->r_0)+alpha)-1.00;
  double power_1mt = g_1*beta-1.00;
  double log_t, log_1mt, g_0log_t, g_1log_1mt;

  for(int i=0; i<n; i++){
    log_t = log(t[i]);
    log_1mt = log1p(-t[i]);
    g_0log_t = g_0*log_t;
    g_1log_1mt = g_1*log_1mt;
    t[i] = log_const+
      nmrover2*logspace_add(log_n + g_1log_1mt, log_saq + g_0log_t)+
      -nmr0over2*logspace_add(log_n + g_1log_1mt + log1p(-Rsq), log_saq + g_0log_t)+
      log(g_0*(1.00-t[i])+g_1*t[i])+
      -(alpha+beta)*logspace_add(log_a+g_0log_t, g_1log_1mt)+
      power_t*log_t+
      power_1mt*log_1mt;
  }
  if(par->log_eval==0){
    for(int i=0; i<n; i++) t[i] = exp(t[i]);
  }
  return;
}

double get_log_bf_scaled_beta(struct model_struct * model){
  double out = 0.00;

  log_bf_integrator.gamma_0 = 2.00/(2.00*log_bf_integrator.shape_0 + log_bf_integrator.r-log_bf_integrator.r_0) + 1.00;//log_bf_integrator.shape_0>0.50 ? 1.00 : 2.00/log_bf_integrator.shape_0;//4.00/(2.00*log_bf_integrator.shape_0 + log_bf_integrator.r-log_bf_integrator.r_0);//
  log_bf_integrator.gamma_1 = 1.00/log_bf_integrator.shape_1 + 1.00;//log_bf_integrator.shape_1>1.00 ? 1.00 : 2.00/log_bf_integrator.shape_1;//2.00/log_bf_integrator.shape_1;//
  //note using ranks instead of sizes
  double rat = (log_bf_integrator.gamma_0*(1.00-log_bf_integrator.t_0) + log_bf_integrator.gamma_1*log_bf_integrator.t_0)/(log_bf_integrator.t_0*(1.00-log_bf_integrator.t_0));
  log_bf_integrator.A[2] = 0.5*rat/log_bf_integrator.scale*log_bf_integrator.q*log_bf_integrator.n*((log_bf_integrator.r-log_bf_integrator.r_0)-log_bf_integrator.Rsq*(log_bf_integrator.n-log_bf_integrator.r_0));
  log_bf_integrator.A[0] = 0.5*rat*log_bf_integrator.n*log_bf_integrator.n*((log_bf_integrator.r-log_bf_integrator.r_0)*(1.00-log_bf_integrator.Rsq));
  log_bf_integrator.A[1] = -log_bf_integrator.A[0]/log_bf_integrator.scale + log_bf_integrator.A[2]*log_bf_integrator.scale;
  log_bf_integrator.B[2] = log_bf_integrator.q*log_bf_integrator.q;
  log_bf_integrator.B[1] = log_bf_integrator.q*log_bf_integrator.n*(2.00-log_bf_integrator.Rsq);
  log_bf_integrator.B[0] = log_bf_integrator.n*log_bf_integrator.n*(1.00-log_bf_integrator.Rsq);
  double num_add = (log_bf_integrator.gamma_1-log_bf_integrator.gamma_0)/(log_bf_integrator.gamma_0*(1.00-log_bf_integrator.t_0) + log_bf_integrator.gamma_1*log_bf_integrator.t_0) - 1.00/log_bf_integrator.t_0 + 1.00/(1.00-log_bf_integrator.t_0);
  log_bf_integrator.C[1] = -rat/log_bf_integrator.scale*(log_bf_integrator.shape_0+log_bf_integrator.shape_1);
  log_bf_integrator.C[0] = log_bf_integrator.shape_0*rat + num_add;
  log_bf_integrator.D[1] = 0.00;
  log_bf_integrator.D[0] = 1.00;
  log_bf_integrator.poly_coef[0] = log_bf_integrator.A[0]*log_bf_integrator.D[0]+log_bf_integrator.B[0]*log_bf_integrator.C[0];
  log_bf_integrator.poly_coef[1] = log_bf_integrator.A[1]*log_bf_integrator.D[0]+log_bf_integrator.B[1]*log_bf_integrator.C[0]+log_bf_integrator.A[0]*log_bf_integrator.D[1]+log_bf_integrator.B[0]*log_bf_integrator.C[1];
  log_bf_integrator.poly_coef[2] = log_bf_integrator.A[2]*log_bf_integrator.D[0]+log_bf_integrator.B[2]*log_bf_integrator.C[0]+log_bf_integrator.A[1]*log_bf_integrator.D[1]+log_bf_integrator.B[1]*log_bf_integrator.C[1];
  log_bf_integrator.poly_coef[3] = log_bf_integrator.A[2]*log_bf_integrator.D[1]+log_bf_integrator.B[2]*log_bf_integrator.C[1];
  cubic_root_finder(log_bf_integrator.poly_coef, log_bf_integrator.poly_root);
  //we do no error catching in the cubic root solving. It would be shocking if there was an issue.
  //we should probably do a check on INFO here
  
  for(int i=0; i<3; i++){
    if(log_bf_integrator.poly_root[i]>0 && 
       fabs(log_bf_integrator.poly_root[i+3])<log_bf_integrator.poly_root_epsilon){
      log_bf_integrator.wat0 = log_bf_integrator.poly_root[i];
      break;
    }
  }

  double wos = log_bf_integrator.wat0/log_bf_integrator.scale;
  // atg0/(atg0+1mtg1) = wos;
  // atg0*(1-wos) = 1mtg1*wos;
  // a = 1mtg1/tg0*wos/(1-wos);
  log_bf_integrator.a = wos/(1.00-wos)*exp(logspace_add(log1p(-log_bf_integrator.t_0)*log_bf_integrator.gamma_1, log(log_bf_integrator.t_0)*log_bf_integrator.gamma_0));
  
  log_bf_integrator.max_log_integrand = 0.00;
  double initial_val = log_bf_integrator.t_0;
  log_bf_integrator.log_eval = 1;
  (*log_bf_integrand)(&initial_val, 1, log_bf_integrator.ex);
  log_bf_integrator.max_log_integrand = initial_val;
  log_bf_integrator.log_eval = 0;
  Rdqags(*log_bf_integrand, log_bf_integrator.ex, &log_bf_integrator.lower, &log_bf_integrator.upper,
         &log_bf_integrator.epsabs, &log_bf_integrator.epsrel, &log_bf_integrator.result, &log_bf_integrator.abserr,
         &log_bf_integrator.neval, &log_bf_integrator.ier, &log_bf_integrator.limit, &log_bf_integrator.lenw,
         &log_bf_integrator.last, log_bf_integrator.iwork, log_bf_integrator.work);
         //we do no error catching in the integration. It would be shocking if the integral did not work.
         //we should probably do a check on ier here
         
         return(log(log_bf_integrator.result)+log_bf_integrator.max_log_integrand);
  return(out);
}

void log_bf_scaled_beta_integrand(double * t, int n, void * ex)
{
  struct log_bf_integrator_struct * par = (struct log_bf_integrator_struct*) ex;
  double n_obs = par->n;
  double log_n = log(n_obs);
  double nmrover2 = 0.5*(n_obs-par->r);
  double nmr0over2 = 0.5*(n_obs-par->r_0);
  double g_0 = par->gamma_0;
  double g_1 = par->gamma_1;
  double a = par->a;
  double log_a = log(a);
  double alpha = par->shape_0;
  double beta = par->shape_1;
  double scale = par->scale;
  double saq = scale*a*par->q;
  double log_saq = log(saq);
  double log_saq_plus_an = logspace_add(log_saq, log_a+log_n);
  double Rsq = par->Rsq;
  double n1mRsq = n_obs*(1.00-Rsq);
  double log_saq_plus_an1mRsq = logspace_add(log_saq, log_a+log_n+log1p(-Rsq));
  double max_log = par->max_log_integrand;
  double log_const = 0.5*(par->r-par->r_0)*log_saq + alpha*log_a - lbeta(alpha, beta) - max_log;
  double power_t = g_0*(0.5*(par->r-par->r_0)+alpha)-1.00;
  double power_1mt = g_1*beta-1.00;
  double log_t, log_1mt, g_0log_t, g_1log_1mt;
  
  for(int i=0; i<n; i++){
    log_t = log(t[i]);
    log_1mt = log1p(-t[i]);
    g_0log_t = g_0*log_t;
    g_1log_1mt = g_1*log_1mt;
    t[i] = log_const+
      nmrover2*logspace_add(log_n + g_1log_1mt, log_saq_plus_an + g_0log_t)+
      -nmr0over2*logspace_add(log_n + g_1log_1mt + log1p(-Rsq), log_saq_plus_an1mRsq + g_0log_t)+
      log(g_0*(1.00-t[i])+g_1*t[i])+
      -(alpha+beta)*logspace_add(log_a+g_0log_t, g_1log_1mt)+
      power_t*log_t+
      power_1mt*log_1mt;
  }
  if(par->log_eval==0){
    for(int i=0; i<n; i++) t[i] = exp(t[i]);
  }
  return;
}

