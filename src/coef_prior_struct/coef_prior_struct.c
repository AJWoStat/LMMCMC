#include "coef_prior_struct.h"

struct coef_prior_struct coef_prior;

void coef_prior_struct_constructor(int eff, double scale, double shape[], const char bf_type[]){
  coef_prior.eff = eff;
  coef_prior.scale = 1.00/scale;
  if(strcmp(bf_type, "zellner-siow")==0){
    strcpy(coef_prior.bf_type, bf_type);
    coef_prior.scale *= 2.00;
    coef_prior.shape_0 = 0.50;
    coef_prior.shape_1 = 1.00;
  }else if(strcmp(bf_type, "hyper-g")==0){
    strcpy(coef_prior.bf_type, bf_type);
    coef_prior.shape_0 = shape[0];
    coef_prior.shape_1 = 1.00;
  }else if(strcmp(bf_type, "intrinsic")==0){
    strcpy(coef_prior.bf_type, bf_type);
    coef_prior.shape_0 = 0.50;
    coef_prior.shape_1 = 0.50;
  }else if(strcmp(bf_type, "inverse-gamma")==0){
    strcpy(coef_prior.bf_type, "gamma");
    coef_prior.shape_0 = shape[0];
    coef_prior.shape_1 = 1.00;
  }else if(strcmp(bf_type, "beta-prime")==0){
    strcpy(coef_prior.bf_type, bf_type);
    coef_prior.shape_0 = shape[1];
    coef_prior.shape_1 = shape[0];
  }else if(strcmp(bf_type, "inverse-beta")==0){
    strcpy(coef_prior.bf_type, "scaled-beta");
    coef_prior.shape_0 = shape[0];
    coef_prior.shape_1 = shape[1];
  }else{
    strcpy(coef_prior.bf_type, "g-prior");
    coef_prior.shape_0 = 1.00;
    coef_prior.shape_1 = 1.00;
  }
  
  return;
}
