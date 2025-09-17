#include "model_space_prior.h"

struct model_space_prior_struct model_space_prior;

void model_space_prior_constructor(const char * model_space_family_in, const double * model_space_parameters_in, const int model_space_parameters_length_in, const int trunc_in){
  model_space_prior.family = malloc(50*sizeof(char));
  strcpy(model_space_prior.family, model_space_family_in);
  model_space_prior.parameters_length = model_space_parameters_length_in;
  model_space_prior.parameters = (double*)calloc(model_space_prior.parameters_length, sizeof(double));
  for(int i=0; i<model_space_prior.parameters_length; i++) model_space_prior.parameters[i] = model_space_parameters_in[i];
  model_space_prior.trunc = trunc_in>data.num_test_var || trunc_in<=0 ? data.num_test_var : trunc_in;
  model_space_prior.log_prior_length = data.num_test_var+1;
  model_space_prior.log_prior_complexity = (double*)calloc(model_space_prior.log_prior_length, sizeof(double));
  model_space_prior.log_prior_models = (double*)calloc(model_space_prior.log_prior_length, sizeof(double));

  if (strcmp(model_space_prior.family, "Uniform") == 0){
    double log_prob = -(double)data.num_test_var*M_LN2;
    for(int i=0; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_models[i] = log_prob;
      model_space_prior.log_prior_complexity[i] = model_space_prior.log_prior_models[i] + lchoose((double)data.num_test_var, (double)i);
    }
    model_space_prior.trunc = -1;
  }else if (strcmp(model_space_prior.family, "Bernoulli") == 0){
    double log_prob = log(model_space_prior.parameters[0]);
    double log_1mprob = log1p(-model_space_prior.parameters[0]);
    for(int i=0; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_models[i] = (double)i*log_prob + (double)(data.num_test_var-i)*log_1mprob;
      model_space_prior.log_prior_complexity[i] = model_space_prior.log_prior_models[i] + lchoose((double)data.num_test_var, (double)i);
    }
    model_space_prior.trunc = -1;
  }else if (strcmp(model_space_prior.family, "Beta-Binomial") == 0){
    double a = model_space_prior.parameters[0];
    double b = model_space_prior.parameters[1];
    double log_nc = lbeta(a, b);
    for(int i=0; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_models[i] = lbeta(a + (double)i, b + (double)(data.num_test_var-i)) - log_nc;
      model_space_prior.log_prior_complexity[i] = model_space_prior.log_prior_models[i] + lchoose((double)data.num_test_var, (double)i);
    }
    model_space_prior.trunc = -1;
  }else if  (strcmp(model_space_prior.family, "Trunc-Beta-Binomial") == 0){
    double a = model_space_prior.parameters[0];
    double b = model_space_prior.parameters[1];
    double log_nc = lbeta(a, b);
    for(int i=0; i<=model_space_prior.trunc; i++){
      model_space_prior.log_prior_models[i] = lbeta(a + (double)i, b + (double)(data.num_test_var-i)) - log_nc;
      model_space_prior.log_prior_complexity[i] = model_space_prior.log_prior_models[i] + lchoose((double)data.num_test_var, (double)i);
    }
    for(int i=model_space_prior.trunc+1; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_models[i] =  -INFINITY;
      model_space_prior.log_prior_complexity[i] = -INFINITY;
    }
    //renormalizing from truncation
    double nc = 0;
    for(int i=0; i<=model_space_prior.trunc; i++) nc+=exp(model_space_prior.log_prior_complexity[i]);
    log_nc = log(nc);
    for(int i=0; i<=model_space_prior.trunc; i++){
      model_space_prior.log_prior_models[i] -= log_nc;
      model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if  (strcmp(model_space_prior.family, "Trunc-Poisson") == 0){
    double lambda = model_space_prior.parameters[0];
    double log_nc = lambda;
    double log_lambda = log(lambda);
    for(int i=0; i<=model_space_prior.trunc; i++){
      model_space_prior.log_prior_models[i] = (double)i*log_lambda - lgammafn((double)i) - log_nc;
      model_space_prior.log_prior_complexity[i] = model_space_prior.log_prior_models[i] + lchoose((double)data.num_test_var, (double)i);
    }
    for(int i=model_space_prior.trunc+1; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_models[i] =  -INFINITY;
      model_space_prior.log_prior_complexity[i] = -INFINITY;
    }
    //renormalizing from truncation
    double nc = 0;
    for(int i=0; i<=model_space_prior.trunc; i++) nc+=exp(model_space_prior.log_prior_complexity[i]);
    log_nc = log(nc);
    for(int i=0; i<=model_space_prior.trunc; i++){
      model_space_prior.log_prior_models[i] -= log_nc;
      model_space_prior.log_prior_complexity[i] -= log_nc;
    }
  }else if  (strcmp(model_space_prior.family, "Trunc-Power-Prior") == 0){
    double u = model_space_prior.parameters[0];
    double log_rate_pos = u*log((double)data.num_test_var); //geometric rate is exp(-log_rate_pos)
    double log_nc = 0;
    log_nc += ((double)(model_space_prior.trunc+1)*log_rate_pos)<M_LN2 ? log(-expm1(-(double)(model_space_prior.trunc+1)*log_rate_pos)) : log1p(-exp(-(double)(model_space_prior.trunc+1)*log_rate_pos));
    log_nc -= log_rate_pos<M_LN2 ? log(-expm1(-log_rate_pos)) : log1p(-exp(-log_rate_pos));
    for(int i=0; i<=model_space_prior.trunc; i++){
      model_space_prior.log_prior_complexity[i] = -(double)i*log_rate_pos - log_nc;
      model_space_prior.log_prior_models[i] =  model_space_prior.log_prior_complexity[i] - lchoose((double)data.num_test_var, (double)i);
    }
    for(int i=model_space_prior.trunc+1; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_models[i] =  -INFINITY;
      model_space_prior.log_prior_complexity[i] = -INFINITY;
    }
  }else if (strcmp(model_space_prior.family, "MatryoshkaDoll") == 0){
    double eta = model_space_prior.parameters[0];
    double nc = 0;
    double log_nc = 0;
    double max_val = 0;
    model_space_prior.log_prior_complexity[model_space_prior.log_prior_length-1] = 0.00;
    for(int i=model_space_prior.log_prior_length-2; i>=0; i--){
      model_space_prior.log_prior_complexity[i] = 0;
      for(int j=i+1; j<model_space_prior.log_prior_length; j++){
        model_space_prior.log_prior_complexity[i] += exp(model_space_prior.log_prior_complexity[j]+lchoose((double)j,(double)i));
      }
      model_space_prior.log_prior_complexity[i] *= eta;
      model_space_prior.log_prior_complexity[i] = log(model_space_prior.log_prior_complexity[i]);
      nc = 0;
      max_val = model_space_prior.log_prior_complexity[i];
      for(int j=(i+1); j<model_space_prior.log_prior_length; j++) max_val = fmax(max_val, model_space_prior.log_prior_complexity[j]);
      for(int j=i; j<model_space_prior.log_prior_length; j++){
        model_space_prior.log_prior_complexity[j] -= max_val;
        nc += exp(model_space_prior.log_prior_complexity[j]);
      }
      log_nc = log(nc);
      for(int j=i; j<model_space_prior.log_prior_length; j++) model_space_prior.log_prior_complexity[j] -= log_nc;
    }
    for(int i=0; i<=data.num_test_var; i++) model_space_prior.log_prior_models[i] = model_space_prior.log_prior_complexity[i] - lchoose((double)data.num_test_var, (double)i);
    model_space_prior.trunc = -1;
  }else{ //here we are at custom and must trust the user a bit and just put in what they want with some minimal checks
    //we just blindly copy the parameter vector, which is to be the log_prior_complexity and compute log_prior_models
    for(int i=0; i<model_space_prior.log_prior_length; i++){
      model_space_prior.log_prior_complexity[i] = model_space_prior.parameters[i];
      model_space_prior.log_prior_models[i] = model_space_prior.log_prior_complexity[i] - lchoose((double)data.num_test_var, (double)i);
    }
  }
  //think about adding truncated geometric (generalizing the power prior thing) and truncated negative binomial (generalizing further)
  return;
}
void model_space_prior_destructor(){
  free(model_space_prior.family); model_space_prior.family = NULL;
  free(model_space_prior.parameters); model_space_prior.parameters = NULL;
  free(model_space_prior.log_prior_complexity); model_space_prior.log_prior_complexity = NULL;
  free(model_space_prior.log_prior_models); model_space_prior.log_prior_models = NULL;
  return;
}
