#ifndef LM_MCMC_FUNCTION_H
#define LM_MCMC_FUNCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "R.h"
#include "Rmath.h"
#include "Rinterface.h"
#include "Rinternals.h"
#include "Rdefines.h"
#include "data_store_struct/data_store_struct.h"
#include "model_struct/model_struct.h"
#include "mcmc_struct/mcmc_struct.h"
#include "hash_table_struct/hash_table_struct.h"
#include "hash_key_computer/hash_key_computer.h"
#include "int_bit_converters/int_bit_converters.h"
#include "model_fit_computer/model_fit_computer.h"
#include "model_space_prior/model_space_prior.h"
#include "log_bf_computer/log_bf_computer.h"
// #include "data_store_struct.h"
// #include "model_struct.h"
// #include "mcmc_struct.h"
// #include "hash_table_struct.h"
// #include "hash_key_computer.h"
// #include "int_bit_converters.h"
// #include "model_fit_computer.h"
// #include "model_space_prior.h"
// #include "log_bf_computer.h"

extern SEXP lm_mcmc_function(SEXP X_in, SEXP y_in, SEXP weights_in, SEXP base_model_indices, //data info
              // SEXP bf_type, SEXP bf_params, SEXP scale_by_rank,  //coefficient prior info - for later
              SEXP msp_family, SEXP msp_params, SEXP trunc, //model space prior info, value of trunc is ignored unless a truncated family is specified
              SEXP include_coef, SEXP include_vcov, //return info (must be 0 for now)
              SEXP start_model_indices, //start model
              SEXP n_burnin, SEXP n_thin, SEXP n_draws, SEXP proposal_probs, //mcmc info
              SEXP hash_table_length, SEXP linked_list_length, SEXP max_load_factor //hash table info
                );


#endif
