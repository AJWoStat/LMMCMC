#include "model_space_prior.h"

struct model_space_prior_struct model_space_prior;

void model_space_prior_constructor(const char * model_space_family_in, const double * model_space_parameters_in, const int model_space_parameters_length_in, const int trunc_in){
  model_space_prior.family = malloc(50*sizeof(char));
  strcpy(model_space_prior.family, model_space_family_in);
  model_space_prior.parameters_length = model_space_parameters_length_in;
  model_space_prior.parameters = (double*)calloc(model_space_prior.parameters_length, sizeof(double));
  for(int i=0; i<model_space_prior.parameters_length; i++) model_space_prior.parameters[i] = model_space_parameters_in[i];
  model_space_prior.trunc = trunc_in>data.num_test_var || trunc_in<=0 ? data.num_test_var : trunc_in; //this is tested on the R-side now, but no reason not to keep this redundancy
  model_space_prior.log_prior_length = data.num_test_var+1;
  model_space_prior.log_prior_complexity = (double*)calloc(model_space_prior.log_prior_length, sizeof(double));
  model_space_prior.log_prior_models = (double*)calloc(model_space_prior.log_prior_length, sizeof(double));
  
  double p = (double)data.num_test_var;
  double log_nc;
  if(model_space_prior.trunc<data.num_test_var) for(int i=model_space_prior.trunc+1; i<=data.num_test_var; i++) model_space_prior.log_prior_complexity[i] = -INFINITY;
  //specifics for the sequences up to trunc for complexity prior
  if (strcmp(model_space_prior.family, "binomial") == 0){
    for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = dbinom((double)i, p, model_space_prior.parameters[0], 1);
    if(model_space_prior.trunc<data.num_test_var){
      log_nc = pbinom((double)model_space_prior.trunc, p, model_space_prior.parameters[0], 1, 1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if (strcmp(model_space_prior.family, "binomial-complexity") == 0){
    double prob = 1.00/(model_space_prior.parameters[0]+model_space_prior.parameters[1]*pow(p, model_space_prior.parameters[2]));
    for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = dbinom((double)i, p, prob, 1);
    if(model_space_prior.trunc<data.num_test_var){
      log_nc = pbinom((double)model_space_prior.trunc, p, prob, 1, 1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if (strcmp(model_space_prior.family, "beta-binomial") == 0){
    double shape_0 = model_space_prior.parameters[0];
    double shape_1 = model_space_prior.parameters[1];
    for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = lchoose(p, (double)i) + lbeta(shape_0+(double)i, shape_1+p-(double)i) - lbeta(shape_0, shape_1);
    if(model_space_prior.trunc<data.num_test_var){
      log_nc = logspace_sum(model_space_prior.log_prior_complexity, model_space_prior.trunc+1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if (strcmp(model_space_prior.family, "beta-binomial-complexity") == 0){
    double shape_0 = model_space_prior.parameters[0];
    double shape_1 = (model_space_prior.parameters[1]+model_space_prior.parameters[2]*pow(p, model_space_prior.parameters[3]));
    for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = lchoose(p, (double)i) + lbeta(shape_0+(double)i, shape_1+p-(double)i) - lbeta(shape_0, shape_1);
    if(model_space_prior.trunc<data.num_test_var){
      log_nc = logspace_sum(model_space_prior.log_prior_complexity, model_space_prior.trunc+1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if (strcmp(model_space_prior.family, "negative-binomial") == 0){
    for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = dnbinom((double)i, model_space_prior.parameters[0], model_space_prior.parameters[1], 1);
    if(model_space_prior.trunc<data.num_test_var){
      log_nc = pnbinom((double)model_space_prior.trunc, model_space_prior.parameters[0], model_space_prior.parameters[1], 1, 1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if (strcmp(model_space_prior.family, "negative-binomial-complexity") == 0){
    double prob = 1.00 - 1.00/(model_space_prior.parameters[1]+model_space_prior.parameters[2]*pow(p, model_space_prior.parameters[3]));
    for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = dnbinom((double)i, model_space_prior.parameters[0], prob, 1);
    if(model_space_prior.trunc<data.num_test_var){
      log_nc = pnbinom((double)model_space_prior.trunc, model_space_prior.parameters[0], prob, 1, 1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if (strcmp(model_space_prior.family, "matryoshka-doll") == 0){
    if(model_space_prior.parameters[1] < 0.5){ // type = 0 for children, type = 1 for descendants
      double rate = 1.00/model_space_prior.parameters[0];
      log_nc = ppois((double)model_space_prior.trunc, rate, 1, 1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = dpois((double)i, rate, 1) - log_nc;
    }else{
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = -INFINITY;
      model_space_prior.log_prior_complexity[model_space_prior.trunc] = 0.00;
      double log_desc = -INFINITY;
      double log_eta = log(model_space_prior.parameters[0]);
      for(int i=model_space_prior.trunc-1; i>=0; i--){
        log_desc = -INFINITY;
        for(int j=i+1; j<=model_space_prior.trunc; j++) log_desc = logspace_add(log_desc, model_space_prior.log_prior_complexity[j] + lchoose((double)j,(double)i));
        model_space_prior.log_prior_complexity[i] = log_eta + log_desc;
        log_nc = logspace_sum(model_space_prior.log_prior_complexity, model_space_prior.trunc+1);
        for(int j=i; j<=model_space_prior.trunc; j++) model_space_prior.log_prior_complexity[j] -= log_nc;
      }
    }
  }else{ 
    //here we are at custom and must trust the user a bit and just put in what they want with some minimal checks
    //we just blindly copy the parameter vector up to trunc and renormalize
    if(model_space_parameters_length_in <= model_space_prior.trunc){
      for(int i=model_space_parameters_length_in; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = -INFINITY;
      for(int i=0; i<model_space_parameters_length_in; i++) model_space_prior.log_prior_complexity[i] = model_space_prior.parameters[i];
      log_nc = logspace_sum(model_space_prior.log_prior_complexity, model_space_prior.trunc+1);
      for(int i=0; i<model_space_parameters_length_in; i++) model_space_prior.log_prior_complexity[i] -= log_nc; 
    }else{
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] = model_space_prior.parameters[i];
      log_nc = logspace_sum(model_space_prior.log_prior_complexity, model_space_prior.trunc+1);
      for(int i=0; i<=model_space_prior.trunc; i++) model_space_prior.log_prior_complexity[i] -= log_nc; 
    }
  }
  
  //compute prior on models of given complexities
  for(int i=0; i<=data.num_test_var; i++){
    model_space_prior.log_prior_models[i] = model_space_prior.log_prior_complexity[i] - lchoose(p, (double)i);
  }
  return;
}

void model_space_prior_destructor(){
  free(model_space_prior.family); model_space_prior.family = NULL;
  free(model_space_prior.parameters); model_space_prior.parameters = NULL;
  free(model_space_prior.log_prior_complexity); model_space_prior.log_prior_complexity = NULL;
  free(model_space_prior.log_prior_models); model_space_prior.log_prior_models = NULL;
  return;
}
