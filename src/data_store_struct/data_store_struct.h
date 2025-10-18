#ifndef DATA_STORE_STRUCT_H
#define DATA_STORE_STRUCT_H

#include <stdlib.h>
#include <math.h>
#include "../int_bit_converters/int_bit_converters.h"
// #include "int_bit_converters.h"
#include <R_ext/Lapack.h>
// #include "lapack.h"
#include <R_ext/BLAS.h>
// #include "cblas.h"

struct data_store_struct{
  //column of 1s added as first column of X and 0 always included in base model model indices
  //all columns of X (but not the 0th column) have been mean centered and standardized so that their sum of squares is num_obs
  //it might be worth making columns of X not in base model orthogonal to the base model
  //y has been mean centered
  //I think that X is column majorized in BAS to match R
  int num_obs; //number of observations
  int num_var; //includes the base model predictors
  int base_model_size; //number of indices in base model
  int base_model_rank; //rank of base model
  int num_test_var;
  int * base_model_indices; //indices of variables that must be included in base model
  int * test_var_indices; //indices of variables to be considered for inclusion
  int bitrep_length;
  uint32_t * base_model_bitrep;
  double intercept_only_model_SSE;
  double base_model_SSE;
  double * XtX; //weights already in calculation
  double * Xty; //weights already in calculation
  double * weights;
  double sqrt_sum_sq_weights;

  //there could be more things stored here about the variables themselves, especially if we want to force heredity
};

extern struct data_store_struct data;

void data_store_constructor(double * X_in, double * y_in, double * weights_in, int num_obs_in, int num_var_in, int * base_model_indices_in, int base_model_size_in);
void data_store_destructor();

void vec_hadamard_product(const int n, double * out, const double * x, const double * y);

void fit_base_model();

#endif
