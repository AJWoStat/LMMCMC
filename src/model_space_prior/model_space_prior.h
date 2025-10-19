#ifndef MODEL_SPACE_PRIOR_H
#define MODEL_SPACE_PRIOR_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Rmath.h"
#include "../data_store_struct/data_store_struct.h"

struct model_space_prior_struct{
  char * family;
  int parameters_length;
  double * parameters;
  int trunc;
  int log_prior_length;
  double * log_prior_complexity;
  double * log_prior_models;
};

extern struct model_space_prior_struct model_space_prior;

void model_space_prior_constructor(const char * model_space_family, const double * model_space_parameters, const int model_space_parameters_length, const int trunc);
void model_space_prior_destructor();

#endif
