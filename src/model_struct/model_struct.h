#ifndef MODEL_STRUCT_H
#define MODEL_STRUCT_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "../int_bit_converters/int_bit_converters.h"
#include "../double_is_equal/double_is_equal.h"
#include "../hash_key_computer/hash_key_computer.h"
#include "R_ext/Print.h"


struct model_struct{
  uint32_t hash_key; // need to compute in bitrep or intarray constructors or when bitrep_replacement or bit_flips is called
  int size, rank; // needed, rank from fit function
  double Rsq, log_prior, log_BF0, residual_sd; //residual_sd computed as root mean square of residuals where the denominator is the degrees of freedom data->num_obs-mod->rank. For models with rank=num_obs the residual sd is reported as nan
  double *coef, *coef_sd, *coef_vcov_upper;

  int bitrep_length; // all functions assume that this is constant over the set of instances
  uint32_t bitrep[];
};

typedef struct model_struct model_t;
typedef model_t* model_ptr_t;

void model_alloc_coef_estimates_and_vcov(model_t * mod, const int include_coef_estimates, const int include_coef_vcov);
void model_realloc_coef_estimates_and_vcov(model_t * mod, const int include_coef_estimates, const int include_coef_vcov);
model_t * model_constructor(const int bitrep_length, const int include_coef_estimates, const int include_coef_vcov);
void model_copy(model_t * mod, const model_t * mod_in);
model_t * model_copy_constructor(const model_t * mod_in);
model_t * model_bitrep_constructor(const uint32_t * bitrep_in, const int bitrep_length_in, const int include_coef_estimates, const int include_coef_vcov);
model_t * model_intarray_constructor(const int * intarray_in, const int intarray_length_in, const int bitrep_length_in, const int include_coef_estimates, const int include_coef_vcov);
void model_destructor(model_t ** mod);
int model_is_equal(const model_t * mod_1, const model_t * mod_2, int only_bitrep);
void model_bitrep_replacement(model_t * mod, const uint32_t * bitrep, const int bitrep_length);
void model_bit_flips(model_t * mod, const int * bits, const int bits_length);
void model_print(const model_t * mod);

#endif
