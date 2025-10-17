#ifndef COEF_PRIOR_STRUCT_H
#define COEF_PRIOR_STRUCT_H

#include<string.h>

struct coef_prior_struct{
  int eff;
  double scale, shape_0, shape_1;
  char bf_type[50];
};

extern struct coef_prior_struct coef_prior;

void coef_prior_struct_constructor(int eff, double scale, double shape[], const char bf_type[]);

#endif