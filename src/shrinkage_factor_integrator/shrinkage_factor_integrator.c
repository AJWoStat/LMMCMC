#include"shrinkage_factor_integrator.h"

static struct shrinkage_factor_integrator_struct shrinkage_factor_integrator;

void shrinkage_factor_integrator_struct_constructor(struct data_store_struct * data, struct coef_prior_struct * coef_prior)
{
  //pointer to self
  shrinkage_factor_integrator.self_ptr = &shrinkage_factor_integrator;
  
  //coefficient prior
  //prior paramaters
  shrinkage_factor_integrator.scale = coef_prior->scale; //prior scale of w = 1/g
  shrinkage_factor_integrator.shape_0 = coef_prior->shape_0; //prior shape at 0 for w=1/g
  shrinkage_factor_integrator.shape_1 = coef_prior->shape_1; //prior shape at infty or 1 for w=1/g (where applicable)
  //type of bayes factor
  strcpy(shrinkage_factor_integrator.bf_type, coef_prior->bf_type);
  
  //integrand and get_shrinkage_factor functions
  if(strcmp(shrinkage_factor_integrator.bf_type, "gamma")==0 || strcmp(shrinkage_factor_integrator.bf_type, "zellner-siow")==0){
    integrand = &shrinkage_factor_gamma_integrand; 
    get_shrinkage_factor_prior_specific = &get_shrinkage_factor_gamma;
  }else if(strcmp(shrinkage_factor_integrator.bf_type, "beta-prime")==0 || strcmp(shrinkage_factor_integrator.bf_type, "hyper-g")==0){
    integrand = &shrinkage_factor_beta_prime_integrand; 
    get_shrinkage_factor_prior_specific = &get_shrinkage_factor_beta_prime;
  }else if(strcmp(shrinkage_factor_integrator.bf_type, "scaled-beta")==0 || strcmp(shrinkage_factor_integrator.bf_type, "intrinsic")==0){
    integrand = &shrinkage_factor_scaled_beta_integrand;
    get_shrinkage_factor_prior_specific = &get_shrinkage_factor_scaled_beta;
  }else{
    integrand = &shrinkage_factor_g_prior_integrand;
    get_shrinkage_factor_prior_specific = &get_shrinkage_factor_g_prior;
  }
  
  //data
  //common data
  shrinkage_factor_integrator.n = (double)(data->num_obs); //number of observations
  shrinkage_factor_integrator.p_0 = (double)(data->base_model_size); //base number of covariates
  shrinkage_factor_integrator.r_0 = (double)(data->base_model_rank); //base model rank
  shrinkage_factor_integrator.eff = coef_prior->eff; //scale g by effective number of parameters? 1 yes, 0 no. 
  //if 1, then normal prior precision is w*p/n*X^T*X where w=1/g
  //if 0, then normal prior precision is w/n*X^T*X where w=1/g
  //model specific data
  shrinkage_factor_integrator.p = (double)(shrinkage_factor_integrator.p_0); //model number of covariates, includes base count
  shrinkage_factor_integrator.r = (double)(shrinkage_factor_integrator.r_0); //model rank, includes base rank
  shrinkage_factor_integrator.q = shrinkage_factor_integrator.eff==1 ? shrinkage_factor_integrator.r : 1.00; //scaling of g effective number of covariates
  //q = eff==1 ? r : 1;
  shrinkage_factor_integrator.Rsq = 0.00; //model Rsq
  
  //integration parameters
  //integral is always from 0 to 1
  //transformation is w(t) = f(g(t)) where f is bf_type specific
  //g(t) = a*t^gamma_0/(a*t^gamma_0+(1-a)*(1-t)^gamma_0)
  
  //transformation parameters for integral
  shrinkage_factor_integrator.gamma_0 = 1.00;
  shrinkage_factor_integrator.gamma_1 = 1.00;
  shrinkage_factor_integrator.t_0 = 0.5; //good default value
  shrinkage_factor_integrator.wat0 = 1.00; //will always be >0
  shrinkage_factor_integrator.a = 1.00; //will always be >0
  shrinkage_factor_integrator.max_log_integrand = 1.00; //subtracted from log integrand before exponentiation and integration
  shrinkage_factor_integrator.gamma_integrand_epsilon = pow(DBL_EPSILON, 0.5);
  shrinkage_factor_integrator.log_eval = 0;
  
  //coefficient parameters of cubic polynomial
  shrinkage_factor_integrator.A[0] = 1.00;
  shrinkage_factor_integrator.A[1] = 1.00;
  shrinkage_factor_integrator.A[2] = 1.00;
  shrinkage_factor_integrator.B[0] = 1.00;
  shrinkage_factor_integrator.B[1] = 1.00;
  shrinkage_factor_integrator.B[2] = 1.00;
  shrinkage_factor_integrator.C[0] = 1.00;
  shrinkage_factor_integrator.C[1] = 1.00;
  shrinkage_factor_integrator.D[0] = 1.00;
  shrinkage_factor_integrator.D[1] = 1.00;
  shrinkage_factor_integrator.poly_coef[0] = 1.00;
  shrinkage_factor_integrator.poly_coef[1] = 1.00;
  shrinkage_factor_integrator.poly_coef[2] = 1.00;
  shrinkage_factor_integrator.poly_coef[3] = 1.00;
  shrinkage_factor_integrator.poly_root[0] = 1.00;
  shrinkage_factor_integrator.poly_root[1] = 1.00;
  shrinkage_factor_integrator.poly_root[2] = 1.00;
  shrinkage_factor_integrator.poly_root[3] = 1.00;
  shrinkage_factor_integrator.poly_root[4] = 1.00;
  shrinkage_factor_integrator.poly_root[5] = 1.00;
  shrinkage_factor_integrator.poly_root_epsilon = pow(DBL_EPSILON, 0.75);
  //variables for Rdqags
  shrinkage_factor_integrator.ex = shrinkage_factor_integrator.self_ptr;
  shrinkage_factor_integrator.lower = 0.00;
  shrinkage_factor_integrator.upper = 1.00;
  shrinkage_factor_integrator.epsabs = pow(DBL_EPSILON, 0.5);
  shrinkage_factor_integrator.epsrel = pow(DBL_EPSILON, 0.5);
  shrinkage_factor_integrator.result = 0.00;
  shrinkage_factor_integrator.abserr = 0.00;
  shrinkage_factor_integrator.neval = 0;
  shrinkage_factor_integrator.ier = 0;
  shrinkage_factor_integrator.limit = 1000;
  shrinkage_factor_integrator.lenw = 4*(shrinkage_factor_integrator.limit);
  shrinkage_factor_integrator.last = 0;
  shrinkage_factor_integrator.lower_ptr = &shrinkage_factor_integrator.lower;
  shrinkage_factor_integrator.upper_ptr = &shrinkage_factor_integrator.upper;
  shrinkage_factor_integrator.epsabs_ptr = &shrinkage_factor_integrator.epsabs;
  shrinkage_factor_integrator.epsrel_ptr = &shrinkage_factor_integrator.epsrel;
  shrinkage_factor_integrator.result_ptr = &shrinkage_factor_integrator.result;
  shrinkage_factor_integrator.abserr_ptr = &shrinkage_factor_integrator.abserr;
  shrinkage_factor_integrator.neval_ptr = &shrinkage_factor_integrator.neval;
  shrinkage_factor_integrator.ier_ptr = &shrinkage_factor_integrator.ier;
  shrinkage_factor_integrator.limit_ptr = &shrinkage_factor_integrator.limit;
  shrinkage_factor_integrator.lenw_ptr = &shrinkage_factor_integrator.lenw;
  shrinkage_factor_integrator.last_ptr = &shrinkage_factor_integrator.last;
  shrinkage_factor_integrator.iwork = (int*)calloc(shrinkage_factor_integrator.limit, sizeof(int));
  shrinkage_factor_integrator.work = (double*)calloc(shrinkage_factor_integrator.lenw, sizeof(double));
  return;
}

void shrinkage_factor_integrator_struct_destructor(){
  free(shrinkage_factor_integrator.iwork); shrinkage_factor_integrator.iwork = NULL;
  free(shrinkage_factor_integrator.work); shrinkage_factor_integrator.work = NULL;
  return;
}

double get_shrinkage_factor(struct model_struct * model){
  //note use of r and not p for q
  shrinkage_factor_integrator.p = (double)(model->size);
  shrinkage_factor_integrator.r = (double)(model->rank);
  shrinkage_factor_integrator.q = shrinkage_factor_integrator.eff == 1 ? shrinkage_factor_integrator.r : 1.00;
  //HERE is a problem?
  shrinkage_factor_integrator.Rsq = 1.00-(1.00 - model->Rsq*data.intercept_only_model_SSE)/data.base_model_SSE;
  shrinkage_factor_integrator.Rsq = fmin(shrinkage_factor_integrator.Rsq, 1.00);
  shrinkage_factor_integrator.Rsq = fmax(shrinkage_factor_integrator.Rsq, 0.00);
  //everything done on log scale and only exponetiated at the end
  double out = (*get_shrinkage_factor_prior_specific)(model);
  return(exp(out-model->log_BF0));
}

double get_shrinkage_factor_g_prior(struct model_struct * model){
  double out = 0.00;
  (*integrand)(&out, 1, shrinkage_factor_integrator.ex);
  return(out);
}

void shrinkage_factor_g_prior_integrand(double * t, int n, void * ex)//actually unused
{ 
  //note use of r and not p
  struct shrinkage_factor_integrator_struct * par = (struct shrinkage_factor_integrator_struct*) ex;
  double scale_times_q_over_n = par->scale * par->q / par->n ;
  t[0] = -log1p(scale_times_q_over_n)+ 0.5*(par->n - par->r)*log1p(scale_times_q_over_n) +
    -0.5*(par->n - par->r_0)*log1p(scale_times_q_over_n - par->Rsq) +
    0.5*(par->r - par->r_0)*log(scale_times_q_over_n);
  return;
}


double get_shrinkage_factor_gamma(struct model_struct * model){
  shrinkage_factor_integrator.gamma_0 = 2.00/(2.00*shrinkage_factor_integrator.shape_0 + shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0) + 1.00;//shrinkage_factor_integrator.shape_0>0.50 ? 1.00 : 2.00/shrinkage_factor_integrator.shape_0;//4.00/(2.00*shrinkage_factor_integrator.shape_0 + shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0);//
  shrinkage_factor_integrator.gamma_1 = 1.00;
  //note using ranks instead of sizes
  double rat = (shrinkage_factor_integrator.gamma_0*(1.00-shrinkage_factor_integrator.t_0) + shrinkage_factor_integrator.gamma_1*shrinkage_factor_integrator.t_0)/(shrinkage_factor_integrator.t_0*(1.00-shrinkage_factor_integrator.t_0));
  shrinkage_factor_integrator.A[2] = 0.00;
  shrinkage_factor_integrator.A[1] = 0.5*rat*shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*((shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0)-shrinkage_factor_integrator.Rsq*(shrinkage_factor_integrator.n-shrinkage_factor_integrator.r_0));
  shrinkage_factor_integrator.A[0] = 0.5*rat*shrinkage_factor_integrator.n*shrinkage_factor_integrator.n*((shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0)*(1.00-shrinkage_factor_integrator.Rsq));
  shrinkage_factor_integrator.B[2] = shrinkage_factor_integrator.q*shrinkage_factor_integrator.q;
  shrinkage_factor_integrator.B[1] = shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*(2.00-shrinkage_factor_integrator.Rsq);
  shrinkage_factor_integrator.B[0] = shrinkage_factor_integrator.n*shrinkage_factor_integrator.n*(1.00-shrinkage_factor_integrator.Rsq);
  shrinkage_factor_integrator.C[1] = -rat/shrinkage_factor_integrator.scale;
  shrinkage_factor_integrator.C[0] = shrinkage_factor_integrator.shape_0*rat + (shrinkage_factor_integrator.gamma_1-shrinkage_factor_integrator.gamma_0)/(shrinkage_factor_integrator.gamma_0*(1.00-shrinkage_factor_integrator.t_0) + shrinkage_factor_integrator.gamma_1*shrinkage_factor_integrator.t_0) - 1.00/shrinkage_factor_integrator.t_0 + 1.00/(1.00-shrinkage_factor_integrator.t_0);
  shrinkage_factor_integrator.D[1] = 0.00;
  shrinkage_factor_integrator.D[0] = 1.00;
  shrinkage_factor_integrator.A[2] -= rat*shrinkage_factor_integrator.q*shrinkage_factor_integrator.q;
  shrinkage_factor_integrator.A[1] -= rat*shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*(1.00-shrinkage_factor_integrator.Rsq);
  shrinkage_factor_integrator.poly_coef[0] = shrinkage_factor_integrator.A[0]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[0]*shrinkage_factor_integrator.C[0];
  shrinkage_factor_integrator.poly_coef[1] = shrinkage_factor_integrator.A[1]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[1]*shrinkage_factor_integrator.C[0]+shrinkage_factor_integrator.A[0]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[0]*shrinkage_factor_integrator.C[1];
  shrinkage_factor_integrator.poly_coef[2] = shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[2]*shrinkage_factor_integrator.C[0]+shrinkage_factor_integrator.A[1]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[1]*shrinkage_factor_integrator.C[1];
  shrinkage_factor_integrator.poly_coef[3] = shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[2]*shrinkage_factor_integrator.C[1];
  cubic_root_finder(shrinkage_factor_integrator.poly_coef, shrinkage_factor_integrator.poly_root);
  // we do no error catching in the cubic root solving. It would be shocking if there was an issue.
  // we should probably do a check on INFO here

  for(int i=0; i<3; i++){
    if(shrinkage_factor_integrator.poly_root[i]>0 && 
       fabs(shrinkage_factor_integrator.poly_root[i+3])<shrinkage_factor_integrator.poly_root_epsilon){
      shrinkage_factor_integrator.wat0 = shrinkage_factor_integrator.poly_root[i];
      break;
    }
  }
  shrinkage_factor_integrator.a = shrinkage_factor_integrator.wat0/shrinkage_factor_integrator.scale*exp(logspace_add(log(1.00-shrinkage_factor_integrator.t_0)*shrinkage_factor_integrator.gamma_1, log(shrinkage_factor_integrator.t_0)*shrinkage_factor_integrator.gamma_0));
  // shrinkage_factor_integrator.a = shrinkage_factor_integrator.wat0*pow(1.00-shrinkage_factor_integrator.t_0, shrinkage_factor_integrator.gamma_1)/(pow(shrinkage_factor_integrator.t_0, shrinkage_factor_integrator.gamma_0)*shrinkage_factor_integrator.scale);
  shrinkage_factor_integrator.max_log_integrand = 0.00;
  double initial_val = shrinkage_factor_integrator.t_0;
  shrinkage_factor_integrator.log_eval = 1;
  (*integrand)(&initial_val, 1, shrinkage_factor_integrator.ex);
  shrinkage_factor_integrator.max_log_integrand = initial_val;
  shrinkage_factor_integrator.log_eval = 0;
  Rdqags(*integrand, shrinkage_factor_integrator.ex, &shrinkage_factor_integrator.lower, &shrinkage_factor_integrator.upper,
         &shrinkage_factor_integrator.epsabs, &shrinkage_factor_integrator.epsrel, &shrinkage_factor_integrator.result, &shrinkage_factor_integrator.abserr,
         &shrinkage_factor_integrator.neval, &shrinkage_factor_integrator.ier, &shrinkage_factor_integrator.limit, &shrinkage_factor_integrator.lenw,
         &shrinkage_factor_integrator.last, shrinkage_factor_integrator.iwork, shrinkage_factor_integrator.work);
  //we do no error catching in the integration. It would be shocking if the integral did not work.
  //we should probably do a check on ier here
  
  return(log(shrinkage_factor_integrator.result)+shrinkage_factor_integrator.max_log_integrand);
}

void shrinkage_factor_gamma_integrand(double * t, int n, void * ex)
{
  struct shrinkage_factor_integrator_struct * par = (struct shrinkage_factor_integrator_struct*) ex;
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
  double log_a = log(a);
  double log_saq = log(saq);
  double Rsq = par->Rsq;
  double n1mRsq = n_obs*(1.00-Rsq);
  double max_log = par->max_log_integrand;
  double log_const = 0.50*(par->r-par->r_0)*log_saq+alpha*log(a) - lgammafn(alpha) - max_log;
  double power_t = g_0*(0.50*(par->r-par->r_0)+alpha)-1.00;
  double power_1mt = g_1*alpha+1.00;
  double log_t, log_1mt, g_0log_t, g_1log_1mt;
  double log_gamma_integrand_epsilon = log(par->gamma_integrand_epsilon);
  double C = a*fmax(1.00, g_0)-(alpha+1.00)*log_a -1.00+0.5*alpha*alpha;
  double delta = 1.00-a*exp(-alpha-sqrt(2.00*(C-log_gamma_integrand_epsilon)));
  for(int i=0; i<n; i++){
    if(t[i]<delta){
      log_t = log(t[i]);
      log_1mt = log1p(-t[i]);
      g_0log_t = g_0*log_t;
      g_1log_1mt = g_1*log_1mt;
      t[i] = log_const+
        g_1log_1mt + log_n - logspace_add(g_1log_1mt+log_n, log_saq + g_0log_t) +
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


double get_shrinkage_factor_beta_prime(struct model_struct * model){
  shrinkage_factor_integrator.gamma_0 = 2.00/(2.00*shrinkage_factor_integrator.shape_0 + shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0) + 1.00;//shrinkage_factor_integrator.shape_0>0.50 ? 1.00 : 2.00/shrinkage_factor_integrator.shape_0;//4.00/(2.00*shrinkage_factor_integrator.shape_0 + shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0);//
  shrinkage_factor_integrator.gamma_1 = 1.00/shrinkage_factor_integrator.shape_1 + 1.00;//shrinkage_factor_integrator.shape_1>1.00 ? 1.00 : 2.00/shrinkage_factor_integrator.shape_1;//2.00/shrinkage_factor_integrator.shape_1;//
  //note using ranks instead of sizes
  double rat = (shrinkage_factor_integrator.gamma_0*(1.00-shrinkage_factor_integrator.t_0) + shrinkage_factor_integrator.gamma_1*shrinkage_factor_integrator.t_0)/(shrinkage_factor_integrator.t_0*(1.00-shrinkage_factor_integrator.t_0));
  shrinkage_factor_integrator.A[2] = 0.00;
  shrinkage_factor_integrator.A[1] = 0.5*rat*shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*((shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0)-shrinkage_factor_integrator.Rsq*(shrinkage_factor_integrator.n-shrinkage_factor_integrator.r_0));
  shrinkage_factor_integrator.A[0] = 0.5*rat*shrinkage_factor_integrator.n*shrinkage_factor_integrator.n*((shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0)*(1.00-shrinkage_factor_integrator.Rsq));
  shrinkage_factor_integrator.B[2] = shrinkage_factor_integrator.q*shrinkage_factor_integrator.q;
  shrinkage_factor_integrator.B[1] = shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*(2.00-shrinkage_factor_integrator.Rsq);
  shrinkage_factor_integrator.B[0] = shrinkage_factor_integrator.n*shrinkage_factor_integrator.n*(1.00-shrinkage_factor_integrator.Rsq);
  double num_add = (shrinkage_factor_integrator.gamma_1-shrinkage_factor_integrator.gamma_0)/(shrinkage_factor_integrator.gamma_0*(1.00-shrinkage_factor_integrator.t_0) + shrinkage_factor_integrator.gamma_1*shrinkage_factor_integrator.t_0) - 1.00/shrinkage_factor_integrator.t_0 + 1.00/(1.00-shrinkage_factor_integrator.t_0);
  shrinkage_factor_integrator.C[1] = -rat*shrinkage_factor_integrator.shape_1 + num_add;
  shrinkage_factor_integrator.C[0] = shrinkage_factor_integrator.scale*shrinkage_factor_integrator.shape_0*rat + shrinkage_factor_integrator.scale*num_add;
  shrinkage_factor_integrator.D[1] = 1.00;
  shrinkage_factor_integrator.D[0] = shrinkage_factor_integrator.scale;
  shrinkage_factor_integrator.A[2] -= rat*shrinkage_factor_integrator.q*shrinkage_factor_integrator.q;
  shrinkage_factor_integrator.A[1] -= rat*shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*(1.00-shrinkage_factor_integrator.Rsq);
  shrinkage_factor_integrator.poly_coef[0] = shrinkage_factor_integrator.A[0]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[0]*shrinkage_factor_integrator.C[0];
  shrinkage_factor_integrator.poly_coef[1] = shrinkage_factor_integrator.A[1]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[1]*shrinkage_factor_integrator.C[0]+shrinkage_factor_integrator.A[0]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[0]*shrinkage_factor_integrator.C[1];
  shrinkage_factor_integrator.poly_coef[2] = shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[2]*shrinkage_factor_integrator.C[0]+shrinkage_factor_integrator.A[1]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[1]*shrinkage_factor_integrator.C[1];
  shrinkage_factor_integrator.poly_coef[3] = shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[2]*shrinkage_factor_integrator.C[1];
  cubic_root_finder(shrinkage_factor_integrator.poly_coef, shrinkage_factor_integrator.poly_root);
  //we do no error catching in the cubic root solving. It would be shocking if there was an issue.
  //we should probably do a check on INFO here
  
  for(int i=0; i<3; i++){
    if(shrinkage_factor_integrator.poly_root[i]>0 && 
       fabs(shrinkage_factor_integrator.poly_root[i+3])<shrinkage_factor_integrator.poly_root_epsilon){
      shrinkage_factor_integrator.wat0 = shrinkage_factor_integrator.poly_root[i];
      break;
    }
  }
  shrinkage_factor_integrator.a = shrinkage_factor_integrator.wat0/shrinkage_factor_integrator.scale*exp(logspace_add(log(1.00-shrinkage_factor_integrator.t_0)*shrinkage_factor_integrator.gamma_1, log(shrinkage_factor_integrator.t_0)*shrinkage_factor_integrator.gamma_0));
  shrinkage_factor_integrator.max_log_integrand = 0.00;
  double initial_val = shrinkage_factor_integrator.t_0;
  shrinkage_factor_integrator.log_eval = 1;
  (*integrand)(&initial_val, 1, shrinkage_factor_integrator.ex);
  shrinkage_factor_integrator.max_log_integrand = initial_val;
  shrinkage_factor_integrator.log_eval = 0;
  Rdqags(*integrand, shrinkage_factor_integrator.ex, &shrinkage_factor_integrator.lower, &shrinkage_factor_integrator.upper,
         &shrinkage_factor_integrator.epsabs, &shrinkage_factor_integrator.epsrel, &shrinkage_factor_integrator.result, &shrinkage_factor_integrator.abserr,
         &shrinkage_factor_integrator.neval, &shrinkage_factor_integrator.ier, &shrinkage_factor_integrator.limit, &shrinkage_factor_integrator.lenw,
         &shrinkage_factor_integrator.last, shrinkage_factor_integrator.iwork, shrinkage_factor_integrator.work);
         //we do no error catching in the integration. It would be shocking if the integral did not work.
         //we should probably do a check on ier here
         
         return(log(shrinkage_factor_integrator.result)+shrinkage_factor_integrator.max_log_integrand);
}

void shrinkage_factor_beta_prime_integrand(double * t, int n, void * ex)
{
  struct shrinkage_factor_integrator_struct * par = (struct shrinkage_factor_integrator_struct*) ex;
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
      g_1log_1mt + log_n - logspace_add(g_1log_1mt+log_n, log_saq + g_0log_t) +
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

double get_shrinkage_factor_scaled_beta(struct model_struct * model){
  double out = 0.00;

  shrinkage_factor_integrator.gamma_0 = 2.00/(2.00*shrinkage_factor_integrator.shape_0 + shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0) + 1.00;//shrinkage_factor_integrator.shape_0>0.50 ? 1.00 : 2.00/shrinkage_factor_integrator.shape_0;//4.00/(2.00*shrinkage_factor_integrator.shape_0 + shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0);//
  shrinkage_factor_integrator.gamma_1 = 1.00/shrinkage_factor_integrator.shape_1 + 1.00;//shrinkage_factor_integrator.shape_1>1.00 ? 1.00 : 2.00/shrinkage_factor_integrator.shape_1;//2.00/shrinkage_factor_integrator.shape_1;//
  //note using ranks instead of sizes
  double rat = (shrinkage_factor_integrator.gamma_0*(1.00-shrinkage_factor_integrator.t_0) + shrinkage_factor_integrator.gamma_1*shrinkage_factor_integrator.t_0)/(shrinkage_factor_integrator.t_0*(1.00-shrinkage_factor_integrator.t_0));
  shrinkage_factor_integrator.A[2] = 0.5*rat/shrinkage_factor_integrator.scale*shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*((shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0)-shrinkage_factor_integrator.Rsq*(shrinkage_factor_integrator.n-shrinkage_factor_integrator.r_0));
  shrinkage_factor_integrator.A[0] = 0.5*rat*shrinkage_factor_integrator.n*shrinkage_factor_integrator.n*((shrinkage_factor_integrator.r-shrinkage_factor_integrator.r_0)*(1.00-shrinkage_factor_integrator.Rsq));
  shrinkage_factor_integrator.A[1] = -shrinkage_factor_integrator.A[0]/shrinkage_factor_integrator.scale + shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.scale;
  shrinkage_factor_integrator.B[2] = shrinkage_factor_integrator.q*shrinkage_factor_integrator.q;
  shrinkage_factor_integrator.B[1] = shrinkage_factor_integrator.q*shrinkage_factor_integrator.n*(2.00-shrinkage_factor_integrator.Rsq);
  shrinkage_factor_integrator.B[0] = shrinkage_factor_integrator.n*shrinkage_factor_integrator.n*(1.00-shrinkage_factor_integrator.Rsq);
  double num_add = (shrinkage_factor_integrator.gamma_1-shrinkage_factor_integrator.gamma_0)/(shrinkage_factor_integrator.gamma_0*(1.00-shrinkage_factor_integrator.t_0) + shrinkage_factor_integrator.gamma_1*shrinkage_factor_integrator.t_0) - 1.00/shrinkage_factor_integrator.t_0 + 1.00/(1.00-shrinkage_factor_integrator.t_0);
  shrinkage_factor_integrator.C[1] = -rat/shrinkage_factor_integrator.scale*(shrinkage_factor_integrator.shape_0+shrinkage_factor_integrator.shape_1);
  shrinkage_factor_integrator.C[0] = shrinkage_factor_integrator.shape_0*rat + num_add;
  shrinkage_factor_integrator.D[1] = 0.00;
  shrinkage_factor_integrator.D[0] = 1.00;
  shrinkage_factor_integrator.A[2] -= rat*shrinkage_factor_integrator.q*(shrinkage_factor_integrator.scale*shrinkage_factor_integrator.q+shrinkage_factor_integrator.n)/shrinkage_factor_integrator.scale;
  shrinkage_factor_integrator.A[1] -= rat*shrinkage_factor_integrator.n*(1.00-shrinkage_factor_integrator.Rsq)*(shrinkage_factor_integrator.scale*shrinkage_factor_integrator.q+shrinkage_factor_integrator.n)/shrinkage_factor_integrator.scale;
  shrinkage_factor_integrator.C[1] += rat/shrinkage_factor_integrator.scale;
  shrinkage_factor_integrator.poly_coef[0] = shrinkage_factor_integrator.A[0]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[0]*shrinkage_factor_integrator.C[0];
  shrinkage_factor_integrator.poly_coef[1] = shrinkage_factor_integrator.A[1]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[1]*shrinkage_factor_integrator.C[0]+shrinkage_factor_integrator.A[0]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[0]*shrinkage_factor_integrator.C[1];
  shrinkage_factor_integrator.poly_coef[2] = shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.D[0]+shrinkage_factor_integrator.B[2]*shrinkage_factor_integrator.C[0]+shrinkage_factor_integrator.A[1]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[1]*shrinkage_factor_integrator.C[1];
  shrinkage_factor_integrator.poly_coef[3] = shrinkage_factor_integrator.A[2]*shrinkage_factor_integrator.D[1]+shrinkage_factor_integrator.B[2]*shrinkage_factor_integrator.C[1];
  cubic_root_finder(shrinkage_factor_integrator.poly_coef, shrinkage_factor_integrator.poly_root);
  //we do no error catching in the cubic root solving. It would be shocking if there was an issue.
  //we should probably do a check on INFO here
  
  for(int i=0; i<3; i++){
    if(shrinkage_factor_integrator.poly_root[i]>0 && 
       fabs(shrinkage_factor_integrator.poly_root[i+3])<shrinkage_factor_integrator.poly_root_epsilon){
      shrinkage_factor_integrator.wat0 = shrinkage_factor_integrator.poly_root[i];
      break;
    }
  }

  double wos = shrinkage_factor_integrator.wat0/shrinkage_factor_integrator.scale;
  // atg0/(atg0+1mtg1) = wos;
  // atg0*(1-wos) = 1mtg1*wos;
  // a = 1mtg1/tg0*wos/(1-wos);
  shrinkage_factor_integrator.a = wos/(1.00-wos)*exp(logspace_add(log1p(-shrinkage_factor_integrator.t_0)*shrinkage_factor_integrator.gamma_1, log(shrinkage_factor_integrator.t_0)*shrinkage_factor_integrator.gamma_0));
  
  shrinkage_factor_integrator.max_log_integrand = 0.00;
  double initial_val = shrinkage_factor_integrator.t_0;
  shrinkage_factor_integrator.log_eval = 1;
  (*integrand)(&initial_val, 1, shrinkage_factor_integrator.ex);
  shrinkage_factor_integrator.max_log_integrand = initial_val;
  shrinkage_factor_integrator.log_eval = 0;
  Rdqags(*integrand, shrinkage_factor_integrator.ex, &shrinkage_factor_integrator.lower, &shrinkage_factor_integrator.upper,
         &shrinkage_factor_integrator.epsabs, &shrinkage_factor_integrator.epsrel, &shrinkage_factor_integrator.result, &shrinkage_factor_integrator.abserr,
         &shrinkage_factor_integrator.neval, &shrinkage_factor_integrator.ier, &shrinkage_factor_integrator.limit, &shrinkage_factor_integrator.lenw,
         &shrinkage_factor_integrator.last, shrinkage_factor_integrator.iwork, shrinkage_factor_integrator.work);
         //we do no error catching in the integration. It would be shocking if the integral did not work.
         //we should probably do a check on ier here
         
         return(log(shrinkage_factor_integrator.result)+shrinkage_factor_integrator.max_log_integrand);
  return(out);
}

void shrinkage_factor_scaled_beta_integrand(double * t, int n, void * ex)
{
  struct shrinkage_factor_integrator_struct * par = (struct shrinkage_factor_integrator_struct*) ex;
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
  double log_q = log(par->q);
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
  double funny_const;
  
  for(int i=0; i<n; i++){
    log_t = log(t[i]);
    log_1mt = log1p(-t[i]);
    g_0log_t = g_0*log_t;
    g_1log_1mt = g_1*log_1mt;
    funny_const = logspace_add(log_a+g_0log_t, g_1log_1mt);
    t[i] = log_const+
      funny_const + log_n - logspace_add(funny_const+log_n, log_saq + g_0log_t) +
      nmrover2*logspace_add(log_n + g_1log_1mt, log_saq_plus_an + g_0log_t)+
      -nmr0over2*logspace_add(log_n + g_1log_1mt + log1p(-Rsq), log_saq_plus_an1mRsq + g_0log_t)+
      log(g_0*(1.00-t[i])+g_1*t[i])+
      -(alpha+beta)*funny_const+
      power_t*log_t+
      power_1mt*log_1mt;
  }
  if(par->log_eval==0){
    for(int i=0; i<n; i++) t[i] = exp(t[i]);
  }
  return;
}

